#version 440 core

#include <renderer.glsl.h>

layout(std140) uniform Camera
{
	CameraConstants camera;
};

layout(std140, binding = 1) uniform Light
{
	LightConstants light;
};


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

void main()
{
	gl_Position = camera.ViewProj * vec4(pos, 1.0);
	vout.worldPos = pos;
	vout.worldNormal = normal;
}

#endif

// Tessellation Control Shader
#ifdef IN_HS

void main() {

}

#endif

// Tessellation Evaluation Shader
#ifdef IN_DS

void main() {
	
}

#endif


#ifdef IN_FS

layout(location = 0) out vec4 color0;

void main()
{

	// (very) Simple Shading
	vec3 normal = normalize(vin.worldNormal);
	float factor = max(-dot(normalize(light.Direction), normal), 0.0);
	color0 = vec4((factor * light.Color));
}

#endif