#version 440 core

#include <renderer.glsl.h>
#include <random.glsl.h>

layout(std140) uniform Camera
{
	CameraConstants camera;
};

layout(std140, binding = 1) uniform Light
{
	LightConstants light;
};

#define PI 3.1415926535897932384626433832795

#ifdef VERTEX_OUT
out VertexData
{
	vec3 worldPos;
	vec3 worldNormal;
} vout VERTEX_OUT_ARRAY;
#endif

#ifdef VERTEX_IN
in VertexData
{
	vec3 worldPos;
	vec3 worldNormal;
} vin VERTEX_IN_ARRAY;
#endif

#ifdef IN_VS

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;
out vec3[3] cPs;
out float baseWidth;

void main()
{
	// Todo: Sliders for (some of these) input values
	// Patch
	float patchSizeX = 40;
	float patchSizeZ = 40;

	// Grass blade
	float minHeight = 0.2;
	float maxHeight = 1.5;

	// 1.0 = full 360 degree rotation possible
	float rotLimitFactorX = 0.01;
	float rotLimitFactorY = 1.0;
	float rotLimitFactorZ = 0.01;

	// For now only X, Y
	float maxInfluencePointVariationX = 0.5;
	float maxInfluencePointVariationY = 0.1;

	// Top of the blade X variation
	float maxPointVariationX = 0.5;

	float maxBaseWidth = 0.4;
	float minBaseWidth = 0.1;



	vec3 position;
	RandState rng = rand_init(gl_VertexID, 11); 

	// Randomize position
	position.x = patchSizeX * rand_next(rng);
	position.z = patchSizeZ * rand_next(rng);
	position.y = 0.0;

	gl_Position = vec4(position, 1.0);
	vout.worldPos = position;

	// Randomize height, shape
	float height = minHeight + rand_next(rng) * (maxHeight - minHeight);
	float influenceX = 2 * (rand_next(rng) - 0.5) * maxInfluencePointVariationX;
	float influenceY = (height / 2) + 2 * (rand_next(rng) - 0.5) * maxInfluencePointVariationY;

	float pointX = 2 * (rand_next(rng) - 0.5) * maxPointVariationX;

	// Control Points for a curve charaterizing a grass blade
	const vec3 controlPoints[3] = vec3[3]
		(
		vec3(0.0, 0.0, 0.0),
		vec3(influenceX, influenceY, 0.0),
		vec3(pointX, height, 0.0)
		);

	// Randomize rotation
	float angleX = 2 * rotLimitFactorX * rand_next(rng) * PI;
	mat3 rotationX = mat3(1.0, 0.0,		    0.0,
						  0.0, cos(angleX), -sin(angleX),
						  0.0, sin(angleX), cos(angleX));

	float angleY = 2 * rotLimitFactorY * rand_next(rng) * PI;
	mat3 rotationY = mat3(cos(angleY),  0.0, sin(angleY),
					      0.0,			1.0, 0.0,
		                  -sin(angleY), 0.0, cos(angleY));

	float angleZ = 2 * rotLimitFactorZ * rand_next(rng) * PI;
	mat3 rotationZ = mat3(cos(angleZ),  -sin(angleZ), 0.0,
						  sin(angleZ),	cos(angleZ), 0.0,
						  0.0,			0.0,		1.0);

	// Randomize (base) width
	baseWidth = minBaseWidth + rand_next(rng) * (maxBaseWidth - minBaseWidth);

	// "real" position of control points
	for (int i = 0; i < 3; i++) {
		cPs[i] = position + rotationX * rotationY * rotationZ * controlPoints[i];
	}

	// Normal in relation to control points
	vout.worldNormal = normalize(cross((cPs[1] - cPs[0]), (cPs[2] - cPs[0])));
}

#endif

// Tessellation Control Shader
#ifdef IN_TCS

layout(vertices = 3) out;
in vec3[][3] cPs;
in float baseWidth[];
patch out vec3 influencePoint;

#define ID gl_InvocationID

void main() {
	if (ID == 0) {

		influencePoint = cPs[0][1];

		// LOD
		float maxDist = 10;
		float minLod = 2;
		float maxLod = 10;
		float lod = minLod + (1 - length(camera.CamPos - cPs[0][0])/maxDist) * (maxLod - minLod);
		if(lod < minLod) {
			lod = minLod;
		}
		gl_TessLevelInner[0] = 1;
		gl_TessLevelOuter[0] = lod;
		gl_TessLevelOuter[1] = lod;
		gl_TessLevelOuter[2] = 1;

		// Place vertices in in such a way that it has the basic shape of the grass blade characterized by the curve 
		gl_out[ID].gl_Position = vec4(cPs[0][0],0.0) + vec4(-baseWidth[0]/2, 0.0, 0.0, 0.0);
	} else if (ID == 1) {
		gl_out[ID].gl_Position = vec4(cPs[0][0],0.0) + vec4(baseWidth[0]/2, 0.0, 0.0, 0.0);
	} else {
		gl_out[ID].gl_Position = vec4(cPs[0][2],0.0);
	}
}

#endif

// Tessellation Evaluation Shader
#ifdef IN_TES

layout(triangles, equal_spacing, cw) in;
patch in vec3 influencePoint;
out vec3 normal;

void main() {
	vec3 pos = vec3(
		gl_TessCoord[0] * gl_in[0].gl_Position +
		gl_TessCoord[1] * gl_in[1].gl_Position +
		gl_TessCoord[2] * gl_in[2].gl_Position);

	// Modify position with the curve
	vec3 iP_Pos = influencePoint - pos;
	float dist = length(iP_Pos);
	float maxDist = min(min(
		length(vec3(gl_in[0].gl_Position) - influencePoint), 
		length(vec3(gl_in[1].gl_Position) - influencePoint)),
		length(vec3(gl_in[2].gl_Position) - influencePoint));
	float factor = (1 - dist / maxDist);
	vec3 curvePos = pos;
	if(factor > 0)	{
		curvePos = pos + factor * iP_Pos;
	}

	gl_Position = camera.ViewProj * vec4(curvePos, 1.0);
	normal = normalize(cross(vec3(gl_in[1].gl_Position - gl_in[0].gl_Position), vec3(gl_in[2].gl_Position - gl_in[0].gl_Position)));
}

#endif


#ifdef IN_FS

layout(location = 0) out vec4 color0;
in vec3 normal;

void main()
{

	vec4 green = vec4(0.1, 0.4, 0.1, 1.0);

	// (very) Simple Shading
	float factor = max(-dot(normalize(light.Direction), normal), 0.0);
	color0 = vec4((factor * (0.9 * green + 0.1 * light.Color)));
	
	
}

#endif