#version 450 core

#include <renderer.glsl.h>

layout(std140, binding = 1) uniform Camera
{
	CameraConstants camera;
};
layout(std140, binding = 2) uniform Light
{
	LightConstants light;
};
layout(std140, binding = 3) uniform Grass
{
	CSGrassConstants grassConsts;
};

#ifdef IN_CS
layout(rgba32f, binding = 0) uniform image2D result;
layout(local_size_x = 1) in;

void main()
{
	imageStore(result, ivec2(gl_GlobalInvocationID.xy), vec4(0.1, 0.8, 0.2, 1.0));
	/*
	vec3 pos = vec3(grassConsts.GridStart.x + gl_GlobalInvocationID.x * grassConsts.Step, grassConsts.GridStart.yz);;
	for(float z = 0; z < grassConsts.GridEnd.z; z += grassConsts.Step) {
		vec4 p1  = camera.ViewProj * vec4(pos.xy, z, 1.0);
		vec3 ndc = p1.xyz / p1.w;
		vec2 screenCoords = (ndc.xy + vec2(1.0, 1.0)) * 1/2;
		if(0.0 <= screenCoords.x && screenCoords.x <= 1.0 &&
			0.0 <= screenCoords.y && screenCoords.y <= 1.0) {
			imageStore(result, ivec2(screenCoords * grassConsts.ScreenDim), vec4(0.1, 0.8, 0.2, 1.0));
		}
	 }*/
}

#endif
