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
uniform vec3 frontCells[32];

#ifdef IN_CS
layout(rgba32f, binding = 0) uniform image2D result;
layout(local_size_x = 32) in;

void main()
{
	vec3 pos = frontCells[gl_LocalInvocationID.x];

	
	for(float i = 0; i < 32; i++) {
		drawWorldPos(pos + i * grassConsts.FtBDirection, vec4(i / 31, 1.0, i / 31, 1.0));
	 }

	 drawWorldPos(vec3(0.0), vec4(1.0, 1.0, 1.0, 1.0));
	 drawWorldPos(vec3(0.0, 0.0, 32.0), vec4(0.0, 0.0, 1.0, 1.0));
	 drawWorldPos(vec3(32, 0.0, 0.0), vec4(1.0, 0.0, 0.0, 1.0));
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
