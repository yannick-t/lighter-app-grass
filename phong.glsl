#version 330 core

#include <renderer.glsl.h>

layout(std140) uniform Camera
{
	CameraConstants camera;
};

#ifdef IN_VS

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 normal;

out vec3 worldPos;
out vec3 normalUnnrm;

void main()
{
	gl_Position = camera.ViewProj * vec4(pos, 1.0);
	worldPos = pos;
	normalUnnrm = normal;
}

#endif
/*
mat3 tangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx(p), dp2 = dFdy(p);
	vec2 duv1 = dFdx(uv), duv2 = dFdy(uv);

	// solve the linear system
	vec3 B = duv1.x * dp2 - duv2.x * dp1;
	vec3 T = duv2.y * dp1 - duv1.y * dp2;

	float D = duv1.x * duv2.y - duv2.x * duv1.y;

	// attention: T, B non-normalized, not orthogonal to N
	return mat3(T * D, B * D, N);
}

mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx(p), dp2 = dFdy(p);
	vec2 duv1 = dFdx(uv), duv2 = dFdy(uv);

	// solve the linear system
	vec3 dp2perp = cross(dp2, N);
	vec3 dp1perp = cross(N, dp1);

	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;

//	float D = 1.0f / (dot(dp1, dp1) * dot(dp2, dp2) - pow(dot(dp1, dp2), 2)); // ||dp1 x dp2||²
	float D = inversesqrt(max(dot(T, T), dot(B, B)));
	// attention: T or B normalized, orthogonal to N
	return mat3(T * D, B * D, N);
}

mat3 cotangent_frame(vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx(p), dp2 = dFdy(p);
	vec2 duv1 = dFdx(uv), duv2 = dFdy(uv);

	// solve the linear system
	vec3 N = cross(dp1, dp2);
	vec3 dp2perp = cross(dp2, N);
	vec3 dp1perp = cross(N, dp1);

	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;

	float DN = inversesqrt(dot(N, N));
	float D = 1.0f / dot(N, N);
	// T, B precisely scaled, N orthonormal
	return mat3(T * D, B * D, N * DN);
}
*/
#ifdef IN_FS

layout(location = 0) out vec4 color0;

in vec3 worldPos;
in vec3 normalUnnrm;

void main()
{
	vec3 normal = normalize(normalUnnrm);
	vec3 color = 0.5 + 0.5 * normal;
	color *= clamp(0.2f + 0.8 * normal.y, 0.f, 1.f);
	color0 = vec4(color, 1.0);
}

#endif