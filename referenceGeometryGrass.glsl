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

layout(std140, binding = 2) uniform GrassPatch
{
	GrassPatchConstants grassPatch;
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
out vec3[4] quadVertices;
out vec3 vInfluencePoint;

void main()
{
	// For each input vertex create a grass blade, later refined via tessellation
	// Input vertex position is discarded if it exists, it will get a random position in the grass patch

	// Todo: Sliders for (some of these) input values
	// Patch
	float patchSizeX = grassPatch.Size;
	float patchSizeZ = grassPatch.Size;

	// Grass blade
	float minHeight = 0.2;
	float maxHeight = grassPatch.MaxHeight;

	// 1.0 = full 360 degree rotation possible
	float rotLimitFactorX = 0.01;
	float rotLimitFactorY = 1.0;
	float rotLimitFactorZ = 0.01;

	// For now only X, Y
	float maxInfluencePointVariationX = 0.5;
	float maxInfluencePointVariationY = 0.1;

	// Top of the blade X variation
	float maxPointVariationX = 0.5;

	float maxBaseWidth = 0.1;
	float minBaseWidth = 0.01;
	float maxBaseWidthDeviationFactor = 1;



	RandState rng = rand_init(gl_VertexID, 11); 

	// Randomize position
	vec3 bladePosition;
	bladePosition.x = patchSizeX * rand_next(rng);
	bladePosition.z = patchSizeZ * rand_next(rng);
	bladePosition.y = 0.0;

	bladePosition = grassPatch.Position + bladePosition;

	gl_Position = vec4(bladePosition, 1.0);
	vout.worldPos = bladePosition;

	/*mat4 translation = mat4(1.0,  0.0,  0.0, bladePosition.x,
							0.0,  1.0,  0.0, bladePosition.y,
							0.0,  0.0,	1.0, bladePosition.z,
							0.0,  0.0,  0.0, 1.0); */

	// Randomize height, shape
	float height = minHeight + rand_next(rng) * (maxHeight - minHeight);

	float pointX = rand_next(rng) * maxPointVariationX;

	float influenceX = rand_next(rng) * maxInfluencePointVariationX;
	float influenceY = (height / 2) + 2 * (rand_next(rng) - 0.5) * maxInfluencePointVariationY;

	

	// Control Points for a curve charaterizing a grass blade
	vec3[] controlPoints = vec3[3]
		(
		vec3(0.0, 0.0, 0.0),
		vec3(influenceX, influenceY, 0.0),
		vec3(pointX, height, 0.0)
		);

	// Randomize rotation
	float angleX = 2 * rotLimitFactorX * rand_next(rng) * PI;
	mat4 rotationX = mat4(1.0, 0.0,		    0.0,		 0.0,
						  0.0, cos(angleX), -sin(angleX),0.0,
						  0.0, sin(angleX), cos(angleX), 0.0,
						  0.0, 0.0,         0.0,         1.0);

	float angleY = 2 * rotLimitFactorY * rand_next(rng) * PI;
	mat4 rotationY = mat4(cos(angleY),  0.0, sin(angleY), 0.0,
					      0.0,			1.0, 0.0,         0.0,
						  -sin(angleY), 0.0, cos(angleY), 0.0,
						  0.0,          0.0, 0.0,         1.0);

	float angleZ = 2 * rotLimitFactorZ * rand_next(rng) * PI;
	mat4 rotationZ = mat4(cos(angleZ),  -sin(angleZ), 0.0,0.0,
						  sin(angleZ),	cos(angleZ),  0.0,0.0,
						  0.0,			0.0,		  1.0,0.0,
						  0.0,          0.0,          0.0,1.0);


	// Transformation matrix for the blade
	mat4 perBladeTransformation = rotationX * rotationY * rotationZ;

	// Randomize (base) width
	float baseWidth = minBaseWidth + rand_next(rng) * (maxBaseWidth - minBaseWidth);
	float baseWidthTopDeviation = baseWidth * (rand_next(rng) * maxBaseWidthDeviationFactor);

	// Calculate the position of the vertices of the quad to be tessellated representing the grass blade
	// Place vertices in in such a way that it has the basic shape of the grass blade characterized by the curve
	quadVertices[0] = bladePosition + vec3(perBladeTransformation * vec4((controlPoints[0] + vec3(-baseWidth/2, 0.0, 0.0)), 1.0));
	quadVertices[1] = bladePosition + vec3(perBladeTransformation * vec4((controlPoints[2] + vec3(-baseWidth/2 + baseWidthTopDeviation, 0.0, 0.0)), 1.0));
	quadVertices[2] = bladePosition + vec3(perBladeTransformation * vec4((controlPoints[2] + vec3(baseWidth/2 - baseWidthTopDeviation, 0.0, 0.0)), 1.0)); // Todo: A bit higher than other top vertex
	quadVertices[3] = bladePosition + vec3(perBladeTransformation * vec4((controlPoints[0] + vec3(baseWidth/2, 0.0, 0.0)), 1.0));

	vInfluencePoint = bladePosition + (perBladeTransformation * vec4(controlPoints[1], 1.0)).xyz;
}

#endif

// Tessellation Control Shader
#ifdef IN_TCS

layout(vertices = 4) out;
in vec3[] vInfluencePoint;
in vec3[][4] quadVertices;
patch out vec3 influencePoint;

#define ID gl_InvocationID

void main() {
	if (ID == 0) {
		influencePoint = vInfluencePoint[0];

		// LOD
		float maxDist = 5;
		float minLod = 2;
		float maxLod = 10;
		float lod = minLod + (1 - length(camera.CamPos - quadVertices[0][0])/maxDist) * (maxLod - minLod);
		if(lod < minLod) {
			lod = minLod;
		}

		gl_TessLevelInner[0] = 1;
		gl_TessLevelInner[0] = 1;

		gl_TessLevelOuter[0] = lod;
		gl_TessLevelOuter[1] = 1;
		gl_TessLevelOuter[2] = lod;
		gl_TessLevelOuter[3] = 1;

		// Todo: Curve grass blades
	}

	gl_out[ID].gl_Position = vec4(quadVertices[0][ID], 1.0);
}

#endif

// Tessellation Evaluation Shader
#ifdef IN_TES

layout(triangles, equal_spacing, ccw) in;
patch in vec3 influencePoint;
out vec3 normal;

void main() {
	// Curves for left and right edge of the blade (quadratic bezier)
	vec3 r = pow(1 - gl_TessCoord.y, 2) * gl_in[3].gl_Position.xyz + 2 * (1 - gl_TessCoord.y) * gl_TessCoord.y * influencePoint + pow(gl_TessCoord.y, 2) * gl_in[2].gl_Position.xyz;
	vec3 l = r + vec3(gl_in[0].gl_Position - gl_in[3].gl_Position);
	vec3 pos = mix(l, r, gl_TessCoord.x);

	gl_Position = camera.ViewProj * vec4(pos, 1.0);
	normal = (transpose(camera.ViewProjInv) * vec4(normalize(cross(vec3(gl_in[1].gl_Position - gl_in[0].gl_Position), vec3(gl_in[2].gl_Position - gl_in[0].gl_Position))), 0.0)).xyz;
	// Todo: Better normals
}

#endif


#ifdef IN_FS

layout(location = 0) out vec4 color0;
in vec3 normal;

void main()
{

	vec4 green = vec4(0.1, 0.4, 0.1, 1.0);

	// (very) Simple Shading
	float factor = min(max(-dot(normalize(light.Direction), normal), 0.1), 1.0);
	color0 = vec4(factor * green);

	// Todo: Better shader
	
}

#endif