#version 330 core

#include <renderer.glsl.h>

layout(std140) uniform Camera
{
	CameraConstants camera;
};

#ifdef IN_VS

out vec3 camDirUnnrm;

void main()
{
	vec2 windowCoords = vec2(
		(gl_VertexID == 2) ? 3.0f : -1.0f,
		(gl_VertexID == 1) ? -3.0f : 1.0f );

	gl_Position = vec4(windowCoords, 0.0, 1.0);

	vec4 rayFromW = camera.ViewProjInv * vec4(windowCoords, -1.25, 1.0);
	vec4 rayOrigW = camera.ViewProjInv * vec4(windowCoords, -1.0, 1.0);
	camDirUnnrm = rayOrigW.xyz / rayOrigW.w - rayFromW.xyz / rayFromW.w;
}

#endif

#ifdef IN_FS

in vec3 camDirUnnrm;
out vec4 color0;

void main()
{
	vec3 camDir = normalize(camDirUnnrm);
	color0 = vec4(0.8 + 0.2 * camDir, 1.0);
//	color0 = vec4(gl_FragCoord.xy * camera.PixelWidth, 0.0, 1.0);
}

#endif