#version 440 core

#include <renderer.glsl.h>

layout(std140) uniform Camera
{
	CameraConstants camera;
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
#ifdef IN_FS

layout(location = 0) out vec4 color0;

void main()
{
	color0 = vec4(0.7, 0.1, 0.1, 1.0);
}

#endif
