#version 450 core

#include <renderer.glsl.h>

layout(std140, binding = 1) uniform Camera
{
	CameraConstants camera;
};
layout(binding = 0) uniform sampler2D result;

#ifdef IN_VS

out vec2 texCoord;

void main()
{
	vec2 windowCoords = vec2(-1.0, -1.0);
	texCoord = vec2(0.0, 1.0);
	if(gl_VertexID == 1) {
		windowCoords.x = 1.0;
		texCoord.s = 1.0;
	} else if(gl_VertexID == 2) {
		windowCoords.y = 1.0;
		texCoord.t = 0.0;
	} else if(gl_VertexID == 3) {
		windowCoords.x = 1.0;
		windowCoords.y = 1.0;
		texCoord.s = 1.0;
		texCoord.t = 0.0;
	}

	gl_Position = vec4(windowCoords, 0.0, 1.0);
}

#endif

#ifdef IN_FS

in vec2 texCoord;
out vec4 color0;

void main()
{
	color0 = texture(result, texCoord);
}

#endif 
