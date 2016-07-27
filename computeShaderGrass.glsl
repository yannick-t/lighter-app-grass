#version 450 core
#extension GL_NV_shader_thread_group : require


#include <renderer.glsl.h>
#include <random.glsl.h>

const float Epsilon = 0.001;

shared vec3 lastPosition;

bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir);
vec3 floorGridPointByDir(vec3 point, vec3 dir);
bool[6] negate(bool array[6]);
float intersectPlane(Ray ray, vec3 point, vec3 normal);
vec3 intersectFrustum(Ray ray, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
vec3 getPointOnFrustumPlane(int index, vec3 frustumPoints[8]);
bool gridCellOutsideFrustum(vec3 center, float size, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
bool quadOutsideFrustum(vec3 points[4], bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]);
bool pointOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 point, float epsilon);
void drawWorldPos(vec3 pos, vec4 colour);
void drawNDC(vec3 ndc, vec4 colour);

vec4 normalPointPlaneToNormalDistPlane(vec3 normal, vec3 point);

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
layout(local_size_x = 32) in;
layout(rgba32f, binding = 0) uniform image2D result;

void main()
{
	vec2 tileCount = ivec2(imageSize(result) / grassConsts.TileDivisor);

	// Calculate frustum of this workgroup
	// Todo: Maybe only calc this once per Workgroup?

	vec3 frustumPoints[8]; // nbl, ntl, nbr, ntr, fbl, ftl, fbr, ftr
	vec3 minFrustum = vec3((1.0 / 0.0));
	vec3 maxFrustum = vec3(-(1.0 / 0.0));
	vec3 frustumPlaneNormals[6]; // right, left, top, bottom, near, far
	Ray frustumRays[12];

	vec4 temp;

	int v = 0;
	for(int i = 0; i < 8; i++) {
		temp = camera.ViewProjInv * vec4(((gl_WorkGroupID.x + int(v / 2)) / tileCount.x) * 2 - 1, ((gl_WorkGroupID.y + mod(v, 2)) / tileCount.y) * 2 - 1, int(i / 4), 1.0);
		frustumPoints[i] = temp.xyz / temp.w;

		// drawWorldPos(frustumPoints[i], vec4(1.0, 1.0, 1.0, 1.0));
		v = int(mod(v + 1, 4));
	}
	
	frustumPlaneNormals[0] = normalize(cross((frustumPoints[3] - frustumPoints[7]), (frustumPoints[6] - frustumPoints[7])));
	frustumPlaneNormals[1] = normalize(cross((frustumPoints[4] - frustumPoints[5]), (frustumPoints[1] - frustumPoints[5])));
	frustumPlaneNormals[2] = normalize(cross((frustumPoints[1] - frustumPoints[5]), (frustumPoints[7] - frustumPoints[5])));
	frustumPlaneNormals[3] = normalize(cross((frustumPoints[6] - frustumPoints[4]), (frustumPoints[0] - frustumPoints[4])));
	frustumPlaneNormals[4] = normalize(cross((frustumPoints[1] - frustumPoints[3]), (frustumPoints[2] - frustumPoints[3])));
	frustumPlaneNormals[5] = normalize(cross((frustumPoints[6] - frustumPoints[7]), (frustumPoints[5] - frustumPoints[7])));


	frustumRays[0] = Ray(frustumPoints[0], (frustumPoints[1] - frustumPoints[0]));
	frustumRays[1] = Ray(frustumPoints[0], (frustumPoints[2] - frustumPoints[0]));
	frustumRays[2] = Ray(frustumPoints[1], (frustumPoints[3] - frustumPoints[1]));
	frustumRays[3] = Ray(frustumPoints[2], (frustumPoints[3] - frustumPoints[2]));

	frustumRays[4] = Ray(frustumPoints[0], (frustumPoints[4] - frustumPoints[0]));
	frustumRays[5] = Ray(frustumPoints[1], (frustumPoints[5] - frustumPoints[1]));
	frustumRays[6] = Ray(frustumPoints[2], (frustumPoints[6] - frustumPoints[2]));
	frustumRays[7] = Ray(frustumPoints[3], (frustumPoints[7] - frustumPoints[3]));

	frustumRays[8] = Ray(frustumPoints[4], (frustumPoints[5] - frustumPoints[4]));
	frustumRays[9] = Ray(frustumPoints[4], (frustumPoints[6] - frustumPoints[4]));
	frustumRays[10] = Ray(frustumPoints[5], (frustumPoints[7] - frustumPoints[5]));
	frustumRays[11] = Ray(frustumPoints[6], (frustumPoints[7] - frustumPoints[6]));

	// intersect frustum rays with plane on which grass should be drawn
	vec3 planePoint = vec3(0.0);
	vec3 planeNormal = vec3(0.0, 1.0, 0.0);

	vec3 intersections[4];
	int intersectionCount = 0;
	for(int i = 0; i < 12; i++) {
		float t = intersectPlane(frustumRays[i], planePoint, planeNormal);
		if(t >= 0 && t < 1) {
			intersections[intersectionCount] = frustumRays[i].Start + t * frustumRays[i].Dir;
			intersectionCount++;
		}
	}

	if(intersectionCount > 3/* && gl_WorkGroupID.x == 2 && gl_WorkGroupID.y == 1*/) {

		// Find closest point to camera
		float minDist = 1.0 / 0.0;
		int index = 0;
		for(int i = 0; i < intersectionCount; i++) {

			float d = length(camera.CamPos - intersections[i]);
			if(d < minDist) {
				minDist = d;
				index = i;
			}
		}
		vec3 start = intersections[index];

		// Round to next cell position 
		start = floorGridPointByDir(start, grassConsts.FtBDirection);

		// Find line of cells to start with -> vector perpendicular to ftb vector at the start vector roughly pointing to other intersection points
		vec3 startLineDir = grassConsts.PerpFtBDir;
		float dotSum = 0;													// If necessary mirror vector so it points roughly the right direction
		for(int i = 0; i < intersectionCount; i++) {
			if(i != index) {
				dotSum += dot(startLineDir, intersections[i] - intersections[index]);
			}
		}
		if(dotSum < 0.0) {		
			startLineDir = -startLineDir;
		}

		// Find planes of the frustum that are on in the general direction of the startLineDir vector
		// Check these during iteration to find out if we're outside of the frustum
		bool frustumPlanesToCheck[6];
		for(int i = 0; i < 6; i++) {
			frustumPlanesToCheck[i] = dot(startLineDir, frustumPlaneNormals[i]) >= 0;
			// drawNDC(vec3(0.5 + float(i) / 24, 0.5, 0.5), vec4(frustumPlanesToCheck[i] ? 0.0 : 1.0, frustumPlanesToCheck[i] ? 1.0 : 0.0, 0.0, 1.0));
		}

		float maxLength = 0;
		for(int i = 0; i < (intersectionCount - 1); i++) {
			for(int j = i + 1; j < intersectionCount; j++) {
				maxLength = max(maxLength, distance(intersections[i], intersections[j]));
			}
		}
		float maxCellLength = ceil(maxLength / grassConsts.Step);


		// debug
		for(int i = 0; i < intersectionCount; i++) {
			drawWorldPos(intersections[i], vec4(0.0, 1.0, 1.0, 1.0)); // Draw intersection points
		}


		// Iterate through cells
		int lineNumber;
		vec3 line = start;
		vec3 firstLine = start;
		vec3 secondLine;
		vec3 lineStart;
		vec3 currentPos;
		bool outOfFrustum;
		bool done;
		int prevLineIn;
		int threadsIn;
		uint ballot;
		for(int i = 0; i < 1; i++) {
			// Find suitable position for thread to draw at by checking if the currentPos is in the frustum and sharing that information with all threads
			outOfFrustum = true;
			done = false;
			prevLineIn = 0;
			threadsIn = 0;
			lineNumber = 0;
			ballot = 0;

			while(outOfFrustum) {
				// Intersect with frustum to find next line
				vec3 sL = line + grassConsts.FtBDirection;
				sL = intersectFrustum(Ray(sL, -startLineDir), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);
				secondLine = floorGridPointByDir(sL, startLineDir);


				// Check if the intersection of the next line indicates that there are more cells of the line in the frustum
				if(lessThanByDir(firstLine, secondLine, startLineDir)) { // Todo: check if neccesary (more than one cell difference (diagonal?))
					lineStart = firstLine;
				} else {
					lineStart = floorGridPointByDir(secondLine - grassConsts.FtBDirection, startLineDir);
				}
				lineStart = firstLine;


				// Check if position is in frustum and share with other threads
				if(!any(isinf(lineStart))) {
					// Find thread position
					currentPos = lineStart + (gl_LocalInvocationID.x - prevLineIn) * startLineDir;
					// Check if in frustum
					outOfFrustum = gridCellOutsideFrustum(currentPos, grassConsts.Step, frustumPlanesToCheck, frustumPlaneNormals, frustumPoints);
					ballot = ballotThreadNV(outOfFrustum);
				} else if(lineNumber > 1) {
					done = true;
					break;
				}

			
				// thread not in frustum, go to next line
				if(outOfFrustum) {
					// update thread count in frustum
					int firstOut = findLSB(ballot);
					prevLineIn = firstOut;

					// next line
					line += grassConsts.FtBDirection;
					firstLine = secondLine;
					lineNumber++;
				} else if (gl_LocalInvocationID.x == 31) {
					// last thread found position to draw at
					lastPosition = currentPos + startLineDir;
				}
			}

			if(done) {
				break;
			}

			// found position -> draw
			drawWorldPos(currentPos + vec3(grassConsts.Step / 2, 0.0, grassConsts.Step / 2), vec4(1.0 - float(lineNumber) / maxCellLength, float(lineNumber) / maxCellLength, gl_LocalInvocationID.x / 32, 1.0));


			// get position of last thread to find the next cells to draw at
			memoryBarrierShared();
			firstLine = lastPosition;
			line = lastPosition;
		}
	}
}

bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir) {
	bool result = false;
	if(any(isinf(fst)) || any(isinf(snd))) {
		result = true;
	} else {
		if(dot(snd - fst, dir) > 0) {
			result = true;
		} 
	}
	return result;
}

vec3 floorGridPointByDir(vec3 point, vec3 dir) {
	vec3 result = vec3(0.0);
	if(any(isinf(point))) {
		result = vec3(1.0 / 0.0);
	} else {
		// floors point to grid point in the opposite direction of dir 
		if(dir.x > 0) {
			result.x = floor(point.x / grassConsts.Step) * grassConsts.Step;
		} else {
			result.x = ceil(point.x / grassConsts.Step) * grassConsts.Step;
		}
		if(dir.z > 0) {
			result.z = floor(point.z / grassConsts.Step) * grassConsts.Step;
		} else {
			result.z = ceil(point.z / grassConsts.Step) * grassConsts.Step;
		}
	}

	return result;
}

bool[6] negate(bool array[6]) {
	bool r[6];
	for(int i = 0; i < 6; i++) {
		r[i] = !array[i];
	}
	return r;
}

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

vec3 intersectFrustum(Ray ray, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]) {
	vec3 intersection = vec3(1.0 / 0.0);
	bool inFrustum = true;
	for(int i = 0; i < 6; i++) {
		if(!whichPlanes[i]) {
			continue;
		}

		float t = intersectPlane(ray, getPointOnFrustumPlane(i, frustumPoints), frustumPlaneNormals[i]);
		intersection = ray.Start + t * ray.Dir;

		// test if intersection is in the frustum
		inFrustum = true;
		int dontCheck = i + ((i % 2 == 1) ? -1 : 1);// plane index opposite of current plane
		for(int j = 0; j < 6; j++) {
			if(j == i || j == dontCheck) {
				continue;
			}

			if(pointOutsideOfPlane(getPointOnFrustumPlane(j, frustumPoints), frustumPlaneNormals[j], intersection, Epsilon)) {
				inFrustum = false;
				break;
			}
		}

		if(inFrustum) {
			break;
		}
	}

	if(inFrustum) {
		return intersection;
	} else {
		return vec3(1.0/0.0);
	}
}

vec3 getPointOnFrustumPlane(int index, vec3 frustumPoints[8]) {
	vec3 point = camera.CamPos;
	if(index == 4) {
		point = frustumPoints[0];
	} else if(index == 5) {
		point = frustumPoints[4];
	}
	return point;
}

bool gridCellOutsideFrustum(vec3 center, float size, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]) {
	vec3 quad[4] = {center + vec3(-size, 0.0, -size) / 2, center + vec3(-size, 0.0, size) / 2,
					center + vec3(size, 0.0, -size) / 2, center + vec3(size, 0.0, size) / 2};
	return quadOutsideFrustum(quad, whichPlanes, frustumPlaneNormals, frustumPoints);
}

bool quadOutsideFrustum(vec3 points[4], bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]) {
	bool outside = false;
	for(int i = 0; i < 6; i++) {
		if(!whichPlanes[i]) {
			continue;
		}
		outside = pointsOutsideOfPlane(getPointOnFrustumPlane(i, frustumPoints), frustumPlaneNormals[i], points);
		if(outside) {
			break;
		}
	}

	return outside;
}

bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]) {
	// outside == on the side of the normal
	for (int i = 0; i < 4; i++) {
		if (!pointOutsideOfPlane(planePoint, planeNormal, points[i], 0.0)) {
			return false;
		}
	}
	return true;
}

bool pointOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 point, float epsilon) {
	return dot(planeNormal, point - (planePoint + epsilon * planeNormal)) >= 0;
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
		ndc.z >= -0.0 && ndc.z <= 1.0) {

		ivec2 imagePos = ivec2(texCoords * imageSize(result));
		imageStore(result, imagePos, colour);
		//imageStore(result, imagePos + ivec2(0, 1), colour);
		//imageStore(result, imagePos + ivec2(1, 0), colour);
		//imageStore(result, imagePos + ivec2(1, 1), colour);
	}
}


vec4 normalPointPlaneToNormalDistPlane(vec3 normal, vec3 point) {
	return vec4(normal, dot(normal, -point));
}
#endif
