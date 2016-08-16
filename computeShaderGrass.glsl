#version 450 core
#extension GL_NV_shader_thread_group : require


#include <renderer.glsl.h>
#include <random.glsl.h>

const float Epsilon = 0.001;

shared vec3 lastPosition;

float getStepSize(vec3 pos);
float getBlend(vec3 pos, float curentStepSize);
float distPointLine(vec3 a, vec3 d, vec3 point);
bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir);
vec3 floorGridPointByDir(vec3 point, vec3 dir, float stepSize);
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

		float stepSize = grassConsts.Step;

		// Find start point for iteration (minimum by ftb)
		int index = 0;
		for(int i = 1; i < intersectionCount; i++) {

			if(lessThanByDir(intersections[i], intersections[index], grassConsts.FtBDirection)) {
				index = i;
			}
		}
		vec3 start = intersections[index];

		stepSize = getStepSize(start);
		// Round to next cell position 
		start = floorGridPointByDir(start, grassConsts.FtBDirection, stepSize);

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


		// debug
		for(int i = 0; i < intersectionCount; i++) {
			drawWorldPos(intersections[i], vec4(0.0, 1.0, 1.0, 1.0)); // Draw intersection points
		}


		// Iterate through cells
		int i = 0;
		int maxIt = 100;

		vec3 ftb = grassConsts.FtBDirection * stepSize;
		vec3 startLine = startLineDir * stepSize;
		int lineNumber = 0;
		vec3 line = start;
		vec3 firstLine = start;
		vec3 secondLine;
		vec3 lineStart;
		vec3 currentPos;
		bool outOfFrustum;
		bool searching;
		bool done;
		bool maskedOut;
		int prevLineIn;
		int threadsIn;
		uint ballot;
		RandState rng;
		for(;;) {
			// Find suitable position for thread to draw at by checking if the currentPos is in the frustum and sharing that information with all threads
			outOfFrustum = true;
			searching = true;
			done = false;
			maskedOut = false;
			prevLineIn = 0;
			threadsIn = 0;
			// lineNumber = 0;
			ballot = 0;

			while(searching) {
				// Intersect with frustum to find next line
				vec3 sL = line + ftb;
				sL = intersectFrustum(Ray(sL, -startLine), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);
				secondLine = floorGridPointByDir(sL, startLine, stepSize);

				// Calculate step size
				stepSize = getStepSize(firstLine);
				ftb = grassConsts.FtBDirection * stepSize;
				startLine = startLineDir * stepSize;


				// Check if the intersection of the next line indicates that there are more cells of the line in the frustum
				if(lessThanByDir(firstLine, secondLine, startLineDir)) { // Todo: check if neccesary (more than one cell difference (diagonal?))
					lineStart = firstLine;
				} else {
					lineStart = floorGridPointByDir(secondLine - ftb, startLine, stepSize);
				}
				// lineStart = firstLine;


				// Check if position is in frustum and share with other threads
				if(!any(isinf(lineStart))) {
					// Find thread position
					int localCompactId = 0;
					if(gl_LocalInvocationID.x > 0) {
						localCompactId = bitCount(ballotThreadNV(true) >> (32 - gl_LocalInvocationID.x));
					}
					currentPos = lineStart + (localCompactId + prevLineIn) * startLine;
					ivec2 iPos = ivec2(round(currentPos.xz / stepSize));
					ivec2 seedPos = ivec2(round(currentPos.xz / grassConsts.Step));
					// Check if in frustum
					outOfFrustum = gridCellOutsideFrustum(currentPos, stepSize, frustumPlanesToCheck, frustumPlaneNormals, frustumPoints);
					rng = rand_init(seedPos.x, seedPos.y);
					// maskedOut = rand_next(rng) < getBlend(currentPos, stepSize);
					// maskedOut = maskedOut && ((iPos.x | iPos.y) & 1) != 0;
				} else {
					outOfFrustum = true;
				}

				if(prevLineIn == 0 && ballotThreadNV(outOfFrustum) == ballotThreadNV(true)) {
					done = true;
					break;
				}

				prevLineIn += bitCount(ballotThreadNV(true));
				searching = maskedOut || outOfFrustum;
			
				// thread not in frustum, go to next line
				if(ballotThreadNV(outOfFrustum) != 0) {
					// next line
					prevLineIn = 0;
					line += ftb;
					firstLine = secondLine;
					lineNumber++;
				} 
				
				if (gl_LocalInvocationID.x == findMSB(ballotThreadNV(true))) {
					// last thread found position to draw at
					lastPosition = currentPos + startLine;
				}

				i++;
				if(i >= maxIt) {
					done = true;
					break;
				}
			}

			if(done) {
				break;
			}

			// found position -> draw
			drawWorldPos(currentPos, vec4(1.0, 0.0, gl_LocalInvocationID.x / 32, 1.0));


			// get position of last thread to find the next cells to draw at
			memoryBarrierShared();
			firstLine = lastPosition;
			line = lastPosition;
		}
	}
}

// Get step size for regular grid by position representing the line currently worked on
float getStepSize(vec3 pos) {
	float dist = length(pos - camera.CamPos); //distPointLine(pos, grassConsts.PerpFtBDir, camera.CamPos);
	float factor = dist / grassConsts.StepDist; //grassConsts.StepDoubleDist;
	
	if(factor <= 1) {
		factor = 1;
	} else {
		factor = pow(2, floor(log2(factor)));
	}
	return factor * grassConsts.Step;
}

float getBlend(vec3 pos, float currentStepSize) {
	float blend;
	float dist = length(pos - camera.CamPos);
	float factor = dist / grassConsts.StepDist; 

	blend = log2(factor / (currentStepSize / grassConsts.Step));

	return blend;
}

float distPointLine(vec3 a, vec3 d, vec3 point) {
	vec3 pud = ((d * (point - a)) / length(d)) * d; // projection
	return length(point - (a + pud));
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

vec3 floorGridPointByDir(vec3 point, vec3 dir, float stepSize) {
	vec3 result = vec3(0.0);
	if(any(isinf(point))) {
		result = vec3(1.0 / 0.0);
	} else {
		// floors point to grid point in the opposite direction of dir 
		if(dir.x > 0) {
			result.x = floor(point.x / stepSize) * stepSize;
		} else {
			result.x = ceil(point.x / stepSize) * stepSize;
		}
		if(dir.z > 0) {
			result.z = floor(point.z / stepSize) * stepSize;
		} else {
			result.z = ceil(point.z / stepSize) * stepSize;
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
