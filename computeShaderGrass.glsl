#version 450 core

#include <renderer.glsl.h>
#include <random.glsl.h>

bool quadInFrustum(vec3 points[4], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]);
void drawWorldPos(vec3 pos, vec4 colour);
void drawNDC(vec3 ndc, vec4 colour);
vec3 maxVec(vec3 a, vec3 b);
vec3 minVec(vec3 a, vec3 b);

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
layout(local_size_x = 1) in;

void main()
{
	vec2 tileCount = ivec2(imageSize(result) / grassConsts.TileDivisor);

	// Calculate frustum of this workgroup
	// Todo: Maybe only calc this once per Workgroup?

	vec3 frustumPoints[8]; // nbl, ntl, nbr, ntr, fbl, ftl, fbr, ftr
	vec3 minFrustum = vec3((1.0 / 0.0));
	vec3 maxFrustum = vec3(-(1.0 / 0.0));
	vec3 frustumPlaneNormals[6]; // right, left, top, bottom, near, far
	vec4 temp;

	int v = 0;
	for(int i = 0; i < 8; i++) {
		temp = camera.ViewProjInv * vec4(((gl_GlobalInvocationID.x + int(v / 2)) / tileCount.x) * 2 - 1, ((gl_GlobalInvocationID.y + mod(v, 2)) / tileCount.y) * 2 - 1, int(i / 4), 1.0);
		frustumPoints[i] = temp.xyz / temp.w;

		minFrustum = minVec(frustumPoints[i], minFrustum);
		maxFrustum = maxVec(frustumPoints[i], maxFrustum);

		drawWorldPos(frustumPoints[i], vec4(1.0, 1.0, 1.0, 1.0));
		v = int(mod(v + 1, 4));
	}
	
	frustumPlaneNormals[0] = normalize(cross((frustumPoints[3] - frustumPoints[7]), (frustumPoints[6] - frustumPoints[7])));
	frustumPlaneNormals[1] = normalize(cross((frustumPoints[4] - frustumPoints[5]), (frustumPoints[1] - frustumPoints[5])));
	frustumPlaneNormals[2] = normalize(cross((frustumPoints[1] - frustumPoints[5]), (frustumPoints[7] - frustumPoints[5])));
	frustumPlaneNormals[3] = normalize(cross((frustumPoints[6] - frustumPoints[4]), (frustumPoints[0] - frustumPoints[4])));
	frustumPlaneNormals[4] = normalize(cross((frustumPoints[1] - frustumPoints[3]), (frustumPoints[2] - frustumPoints[3])));
	frustumPlaneNormals[5] = normalize(cross((frustumPoints[6] - frustumPoints[7]), (frustumPoints[5] - frustumPoints[7])));

	// Determine where to start / end the iteration based on the regular grid direction
	vec2 start;
	vec2 end; 
	if(dot(normalize(grassConsts.FtBDirection), normalize(maxFrustum - minFrustum)) >= 0.0) {
		start = minFrustum.xz;
		end = maxFrustum.xz;
	} else {
		start = maxFrustum.xz;
		end = minFrustum.xz;
	}

	if(gl_GlobalInvocationID.x == 0 && gl_GlobalInvocationID.y == 0) {
		drawWorldPos(vec3(end.x, 0.0, end.y), vec4(0.0, 0.0, 1.0, 1.0));
		drawWorldPos(vec3(start.x, 0.0, start.y), vec4(0.0, 1.0, 0.0, 1.0));
	}

	// Round down to next cell position
	start = floor(start / grassConsts.Step) * vec2(grassConsts.Step);
	// Round up to next cell position
	end = ceil(end / grassConsts.Step) * vec2(grassConsts.Step);

	
	for (float x = start.x; x <= end.x; x += grassConsts.Step) {
		for (float z = start.y; z <= end.y; z += grassConsts.Step) {
			vec3 cellPos = vec3(x, 0.0, z);
			vec3 quad[4] = {cellPos, cellPos + vec3(grassConsts.Step, 0.0, 0.0), 
						cellPos + vec3(0.0, 0.0, grassConsts.Step), cellPos + vec3(grassConsts.Step, 0.0, grassConsts.Step)};
			if(quadInFrustum(quad, frustumPlaneNormals, frustumPoints)) {
				drawWorldPos(cellPos + vec3(grassConsts.Step / 2, 0.0, grassConsts.Step / 2), vec4(1.0, 0.0, 0.0, 1.0));
			}
		}
	}


	/*
	// Find Line of cells to start with -> vector with right angle to ftb vector at the start vector roughly pointing to end vector
	mat4 r = mat4(0.0,  0.0, 1.0, 0.0,									// Rotate 90 degrees on the around the Y axis
					0.0,	1.0, 0.0, 0.0,
					-1.0, 0.0, 0.0, 0.0,
					0.0,  0.0, 0.0, 1.0);
	vec3 startLine = vec3(r * vec4(grassConsts.FtBDirection, 1.0));
	if(dot(normalize(startLine), normalize(vec3(end.x, 0.0, end.y) - vec3(start.x, 0.0, start.y))) < 0.0) {		// If necessary mirror vector so it points roughly to end
		startLine = -startLine;
	}

	// Iterate through the cells
	int maxDepth = 32;
	bool wasInFrustum = false;

	vec3 cellPos = vec3(start.x, 0.0, start.y) + gl_LocalInvocationID.x * startLine;
	for(int i = 0; i < maxDepth; i++) {
		vec3 quad[4] = {cellPos, cellPos + vec3(grassConsts.Step, 0.0, 0.0), 
						cellPos + vec3(0.0, 0.0, grassConsts.Step), cellPos + vec3(grassConsts.Step, 0.0, grassConsts.Step)};
		bool inFrustum = quadInFrustum(quad, frustumPlaneNormals, frustumPoints);
		if(!inFrustum && wasInFrustum) {
			break;
		}
		if(inFrustum) {
			drawWorldPos(cellPos + vec3(grassConsts.Step / 2, 0.0, grassConsts.Step / 2), vec4(1 - float(i) / maxDepth, float(i) / maxDepth, 0.0, 1.0));
		}
		wasInFrustum = inFrustum;

		cellPos += grassConsts.FtBDirection;
	}*/
}

bool quadInFrustum(vec3 points[4], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]) {
	// Find out if points of the quad are all outside of one of the frustums planes

	if (pointsOutsideOfPlane(camera.CamPos, frustumPlaneNormals[0], points) || pointsOutsideOfPlane(camera.CamPos, frustumPlaneNormals[1], points) || pointsOutsideOfPlane(camera.CamPos, frustumPlaneNormals[2], points) ||
		pointsOutsideOfPlane(camera.CamPos, frustumPlaneNormals[3], points) || pointsOutsideOfPlane(frustumPoints[0], frustumPlaneNormals[4], points) || pointsOutsideOfPlane(frustumPoints[4], frustumPlaneNormals[5], points)) {

		return false;
	}
	return true;
}

bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]) {
	// outside == on the side of the normal
	for (int i = 0; i < 4; i++) {
		if (dot(planeNormal, points[i] - planePoint) < 0) {
			return false;
		}
	}
	return true;
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
		ndc.z >= -0.000001 && ndc.z <= 1.00001) {
		imageStore(result, ivec2(texCoords * imageSize(result)), colour);
	}
}

vec3 maxVec(vec3 a, vec3 b) {
	vec3 r = b;
	if(a.x > b.x) {
		r.x = a.x;
	}
	if(a.y > b.y) {
		r.y = a.y;
	}
	if(a.z > b.z) {
		r.z = a.z;
	}

	return r;
}

vec3 minVec(vec3 a, vec3 b) {
	vec3 r = b;
	if(a.x < b.x) {
		r.x = a.x;
	}
	if(a.y < b.y) {
		r.y = a.y;
	}
	if(a.z < b.z) {
		r.z = a.z;
	}

	return r;
}

#endif
