#version 440 core

#include <renderer.glsl.h>

const float Epsilon = 0.001;

layout(std140) uniform Camera
{
	CameraConstants camera;
};

layout(std140, binding = 1) uniform Light
{
	LightConstants light;
};

layout(std140, binding = 2) uniform Ground
{
	GroundConstants ground;
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

Ray[12] calcFrustumRays(vec3[8] frustumPoints);
vec3[4] calcFrustumIntersections(vec3 planePoint, vec3 planeNormal, Ray[12] frustumRays, out int intersectionCount);
float intersectPlane(Ray ray, vec3 point, vec3 normal);
vec3 worldPosToNDC(vec3 pos);
vec3 nDCToWorldPos(vec3 ndc);

// Quick and dirty implementation of a ground beneath the grass blades

void main()
{
	vec3[8] frustumPoints;
	vec3[8] frustumPointsNDC;
	int v = 0;
	for(int i = 0; i < 8; i++) {
		frustumPointsNDC[i] = vec3(int(v / 2) * 2 - 1, -(((mod(v + 1, 2))) * 2 - 1), int(i / 4));
		frustumPoints[i] = nDCToWorldPos(frustumPointsNDC[i]);
		v = int(mod(v + 1, 4));
	}
	Ray[12] frustumRays = calcFrustumRays(frustumPoints);
	int intersectionCount;
	vec3[4] frustumIntersections = calcFrustumIntersections(vec3(0), vec3(0, 1, 0), frustumRays, intersectionCount);

	vec3 position = frustumIntersections[gl_VertexID];

	gl_Position = camera.ViewProj * vec4(position, 1.0);
	vout.worldPos = position;
	vout.worldNormal = vec3(0, 1, 0);
}

Ray[12] calcFrustumRays(vec3[8] frustumPoints) {
	Ray[12] rays;
	rays[0] = Ray(frustumPoints[0], (frustumPoints[1] - frustumPoints[0]));
	rays[1] = Ray(frustumPoints[0], (frustumPoints[2] - frustumPoints[0]));
	rays[2] = Ray(frustumPoints[1], (frustumPoints[3] - frustumPoints[1]));
	rays[3] = Ray(frustumPoints[2], (frustumPoints[3] - frustumPoints[2]));

	rays[4] = Ray(frustumPoints[0], (frustumPoints[4] - frustumPoints[0]));
	rays[5] = Ray(frustumPoints[1], (frustumPoints[5] - frustumPoints[1]));
	rays[6] = Ray(frustumPoints[2], (frustumPoints[6] - frustumPoints[2]));
	rays[7] = Ray(frustumPoints[3], (frustumPoints[7] - frustumPoints[3]));

	rays[8] = Ray(frustumPoints[4], (frustumPoints[5] - frustumPoints[4]));
	rays[9] = Ray(frustumPoints[4], (frustumPoints[6] - frustumPoints[4]));
	rays[10] = Ray(frustumPoints[5], (frustumPoints[7] - frustumPoints[5]));
	rays[11] = Ray(frustumPoints[6], (frustumPoints[7] - frustumPoints[6]));

	return rays;
}

vec3[4] calcFrustumIntersections(vec3 planePoint, vec3 planeNormal, Ray[12] frustumRays, out int intersectionCount) {
	vec3[4] intersections;
	intersectionCount = 0;
	for(int i = 0; i < 12; i++) {
		float t = intersectPlane(frustumRays[i], planePoint, planeNormal);
		if(t >= 0 && t < 1) {
			intersections[intersectionCount] = frustumRays[i].Start + t * frustumRays[i].Dir;
			intersectionCount++;
		}
	}

	return intersections;
}

// methods for intersection
float intersectPlane(Ray ray, vec3 point, vec3 normal) {
	float denom = dot(normal, ray.Dir);
	float t = 0;
	if (abs(denom) > Epsilon) {
		t = dot((point - ray.Start), normal) / denom;
	} else {
		t = 1.0 / 0.0;
	}
	return t;
}

vec3 worldPosToNDC(vec3 pos) {
	vec4 p1 = camera.ViewProj * vec4(pos, 1.0);
	return p1.xyz / p1.w;
}

vec3 nDCToWorldPos(vec3 ndc) {
	vec4 p1 = camera.ViewProjInv * vec4(ndc, 1.0);
	return p1.xyz / p1.w;
}

#endif

#ifdef IN_FS

layout(location = 0) out vec4 color0;

void main()
{
	vec3 normal = vin.worldNormal;
	float factor = max(-dot(light.Direction, normal), 0.1);
	color0 = 0.4 * factor * vec4(0.2, 0.05, 0, 1);

	float dist = distance(vin.worldPos, camera.CamPos);


	float f1 = clamp((dist - ground.MinDist - 0.1) / 0.1, 0, 1);
	float f2 = clamp((ground.MaxDist - dist - 0.1) / 0.1, 0, 1);

	color0 = mix(color0, color0 * 0.4, min(f1, f2));
	
}

#endif
