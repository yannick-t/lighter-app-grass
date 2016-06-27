#version 450 core

#include <renderer.glsl.h>
#include <random.glsl.h>

void drawWorldPos(vec3 pos, vec4 colour);
void drawNDC(vec3 ndc, vec4 colour);

layout(std140, binding = 1) uniform Camera
{
	CameraConstants camera;
};
layout(std140, binding = 3) uniform Light
{
	LightConstants light;
};
layout(std140, binding = 2) uniform Grass
{
	CSGrassConstants grassConsts;
};

#ifdef IN_CS
layout(rgba32f, binding = 0) uniform image2D result;
layout(local_size_x = 32) in;

void main()
{
	// Todo: Maybe only calc this once per Workgroup?
	vec3 frustumPoints[8];
	vec3 minFrustum;
	vec3 maxFrustum;
	vec3 frustumPlaneNormals[6]
	vec4 temp = camera.ViewProj * vec4(gl_GlocalInvocationID.x / 32, gl_GlobalInvocationID.y / 32, 0.0, 1.0);
	frustumPoints[0] = temp.xyz / temp.w;
	temp = camera.ViewProj * vec4((gl_GlocalInvocationID.x + 1) / 32, (gl_GlobalInvocationID.y) / 32, 0.0, 1.0);
	frustumPoints[1] = ntr.xyz / ntr.w;
	temp = camera.ViewProj * vec4((gl_GlocalInvocationID.x + 1) / 32, (gl_GlobalInvocationID.y + 1) / 32, 0.0, 1.0);
	frustumPoints[2] = ntr.xyz / ntr.w;
	temp = camera.ViewProj * vec4((gl_GlocalInvocationID.x) / 32, (gl_GlobalInvocationID.y + 1) / 32, 0.0, 1.0);
	frustumPoints[3] = ntr.xyz / ntr.w;
	temp = camera.ViewProj * vec4(gl_GlocalInvocationID.x / 32, gl_GlobalInvocationID.y / 32, 1.0, 1.0);
	frustumPoints[4] = temp.xyz / temp.w;
	temp = camera.ViewProj * vec4((gl_GlocalInvocationID.x + 1) / 32, (gl_GlobalInvocationID.y) / 32, 1.0, 1.0);
	frustumPoints[5] = ntr.xyz / ntr.w;
	temp = camera.ViewProj * vec4((gl_GlocalInvocationID.x + 1) / 32, (gl_GlobalInvocationID.y + 1) / 32, 1.0, 1.0);
	frustumPoints[6] = ntr.xyz / ntr.w;
	temp = camera.ViewProj * vec4((gl_GlocalInvocationID.x) / 32, (gl_GlobalInvocationID.y + 1) / 32, 1.0, 1.0);
	frustumPoints[7] = ntr.xyz / ntr.w;

	for(int i = 0; i < 8; i++) {
		
	}
}

void drawWorldPos(vec3 pos, vec4 colour) {
		vec4 p1  = camera.ViewProj * vec4(pos, 1.0);
		vec3 ndc = p1.xyz / p1.w;
		drawNDC(ndc, colour);
}

void drawNDC(vec3 ndc, vec4 colour) {
	vec2 screenCoords = (ndc.xy + vec2(1.0, 1.0)) * 1/2;
	vec2 texCoords = vec2(screenCoords.x, 1.0 - screenCoords.y);
	if(0.0 <= texCoords.x && texCoords.x <= 1.0 &&
		0.0 <= texCoords.y && texCoords.y <= 1.0 &&
		ndc.z >= 0.0 && ndc.z <= 1.0) {
		imageStore(result, ivec2(texCoords * imageSize(result)), colour);
	}
}

#endif
