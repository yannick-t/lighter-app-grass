#version 330 core

#ifdef IN_VS

void main()
{
	vec2 windowCoords = vec2(
		(gl_VertexID == 2) ? 3.0f : -1.0f,
		(gl_VertexID == 1) ? -3.0f : 1.0f );

	gl_Position = vec4(windowCoords, 0.0, 1.0);
}

#endif

#ifdef IN_FS

uniform sampler2D hdrColorIn;
layout(location = 0) out vec4 ldrColorOut;

void main()
{
	ldrColorOut = texelFetch(hdrColorIn, ivec2(gl_FragCoord.xy), 0);
//	ldrColorOut = vec4(1.0f, 0.0f, 0.0f, 1.0f);
//	ldrColorOut.xyz /= (1.0f + dot(ldrColorOut.xyz, vec3(0.3f, 0.5f, 0.2f)));
}

#endif