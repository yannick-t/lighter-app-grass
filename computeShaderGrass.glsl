#version 450 core
#extension GL_NV_shader_thread_group : require

#include <renderer.glsl.h>
#include <random.glsl.h>

const float Epsilon = 0.0000001;

shared vec3 lastPosition;

bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir);
vec3 floorGridPointByDir(vec3 point, vec3 dir);
bool[6] negate(bool array[6]);
float intersectPlane(Ray ray, vec3 point, vec3 normal);
vec3 intersectFrustum(Ray ray, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
vec3 getPointOnFrustumPlane(int index, vec3 frustumPoints[8]);
bool quadOutsideFrustum(vec3 points[4], bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]);
bool pointOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 point);
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
	Ray frustumRays[12];

	vec4 temp;

	int v = 0;
	for(int i = 0; i < 8; i++) {
		temp = camera.ViewProjInv * vec4(((gl_GlobalInvocationID.x + int(v / 2)) / tileCount.x) * 2 - 1, ((gl_GlobalInvocationID.y + mod(v, 2)) / tileCount.y) * 2 - 1, int(i / 4), 1.0);
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

	if(intersectionCount > 3 && gl_GlobalInvocationID.x == 1 && gl_GlobalInvocationID.y == 1) {

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


		// intersect with the frustum to get accurate start position
		vec3 firstLine = intersectFrustum(Ray(start, -startLineDir), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);
		vec3 secondLine = intersectFrustum(Ray(start + grassConsts.FtBDirection, -startLineDir), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);

		// Debug
		drawNDC(vec3(0.0, 0.0, 0.5), vec4(any(isinf(firstLine)) ? 0.0 : 1.0, any(isinf(firstLine)) ? 1.0 : 0.0, 0.0, 1.0));
		drawNDC(vec3(0.1, 0.0, 0.5), vec4(any(isinf(secondLine)) ? 0.0 : 1.0, any(isinf(secondLine)) ? 1.0 : 0.0, 0.0, 1.0));
		
		for(int i = 0; i < 128; i++) {
			vec3 a = intersectFrustum(Ray(start + float(i) * grassConsts.FtBDirection, -startLineDir), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);
			if(!any(isinf(a))) {
				// drawWorldPos(a, vec4(0.5, 1.0, 0.0, 1.0));
			}
		}


		if(gl_GlobalInvocationID.x == 1 && gl_GlobalInvocationID.y == 1) {
			for(int i = 0; i < intersectionCount; i++) {
				drawWorldPos(intersections[i], vec4(0.0, 1.0, 1.0, 1.0)); // Draw intersection points
			}
			drawWorldPos(start, vec4(1.0, 1.0, 1.0, 1.0));
			drawWorldPos(start + grassConsts.FtBDirection, vec4(1.0, 1.0, 1.0, 1.0));
			drawWorldPos(firstLine, vec4(1.0, 0.0, 0.0, 1.0));
			drawWorldPos(secondLine, vec4(0.0, 1.0, 0.0, 1.0));
		}

		if(any(isinf(firstLine))) {
			firstLine = start;
		} else {
			firstLine = floorGridPointByDir(firstLine, startLineDir);
		}
		secondLine = floorGridPointByDir(secondLine, startLineDir);
	
		/*
		vec3 lineStart;
		vec3 currentPos;
		int i = 0;
		int lineNumber = 0;
		bool outOfFrustum = false;
		for(; lineNumber < 256;) {
			if(lessThanByDir(firstLine, secondLine, startLineDir)) {
				lineStart = firstLine;	
			} else {
				lineStart = secondLine - grassConsts.FtBDirection;
			}

			currentPos = lineStart + i * startLineDir;
			i++;
			vec3 quad[4] = {currentPos, currentPos + vec3(grassConsts.Step, 0.0, 0.0), 
							currentPos + vec3(0.0, 0.0, grassConsts.Step), currentPos + vec3(grassConsts.Step, 0.0, grassConsts.Step)};
			outOfFrustum = quadOutsideFrustum(quad, frustumPlanesToCheck, frustumPlaneNormals, frustumPoints);

			if(!outOfFrustum) {
				drawWorldPos(currentPos, vec4(1.0, 0.5, 0.5, 1.0));
			} else if(i == 0 && lineNumber != 0) {
				break;
			} else {
				// Go to next line
				firstLine = secondLine;
				secondLine = currentPos + grassConsts.FtBDirection;
				secondLine = intersectFrustum(Ray(lineStart, -startLineDir), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);
				secondLine = floorGridPointByDir(secondLine, startLineDir);
				i = 0;
				lineNumber++;
			}
		}


		/*

		// Iterate through cells
		vec3 lineStart;
		vec3 currentPos;
		int lineNumber = 0;
		while(true) {
			// Find suitable position for thread to draw at by checking if the currentPos is in the frustum and sharing that information with all threads
			bool outOfFrustum = true;
			bool done = false;
			int prevLineIn = 0;
			int threadsIn = 0;
			while(outOfFrustum) {
				if(lessThanByDir(firstLine, secondLine, startLineDir)) {
					lineStart = firstLine;	
				} else {
					lineStart = secondLine - ftbDirs[0];
				}

				currentPos = lineStart + (gl_LocalInvocationID.x - threadsIn) * startLineDir;
				vec3 quad[4] = {currentPos, currentPos + vec3(grassConsts.Step, 0.0, 0.0), 
								currentPos + vec3(0.0, 0.0, grassConsts.Step), currentPos + vec3(grassConsts.Step, 0.0, grassConsts.Step)};
				outOfFrustum = quadOutsideFrustum(quad, frustumPlanesToCheck, frustumPlaneNormals, frustumPoints);

				uint ballot = ballotThreadNV(outOfFrustum) >> threadsIn;
				
				if(outOfFrustum) {
					int firstOut = findLSB(ballot);
					if(firstOut == 0 && prevLineIn != 0) {
						// last line
						done = true;
						break;
					}
					prevLineIn = firstOut;
					threadsIn += firstOut;


					// Go to next line
					firstLine = secondLine;
					secondLine = currentPos + ftbDirs[lineNumber % ftbDirCount];
					secondLine = intersectFrustum(Ray(lineStart, -startLineDir), negate(frustumPlanesToCheck), frustumPlaneNormals, frustumPoints);
					secondLine = floorGridPointByDir(secondLine, startLineDir);
					lineNumber++;
				}
			}

			if(done) {
				break;
			}

			// drawWorldPos(currentPos + vec3(grassConsts.Step / 2, 0.0, grassConsts.Step / 2), vec4(1.0 - float(lineNumber) / maxCellLength, float(lineNumber) / maxCellLength, gl_LocalInvocationID.x / 32, 1.0));
			drawWorldPos(currentPos + vec3(grassConsts.Step / 2, 0.0, grassConsts.Step / 2), vec4(1.0, 0.5, 0.5, 1.0));

			// Todo: next 32
			break;
		}*/
	}
}

bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir) {
	bool result = false;
	if(dot(snd - fst, dir) > 0) {
		result = true;
	}
	return result;
}

vec3 floorGridPointByDir(vec3 point, vec3 dir) {
	vec3 result = vec3(0.0);
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

			if(pointOutsideOfPlane(frustumPlaneNormals[j], getPointOnFrustumPlane(j, frustumPoints), intersection)) {
				inFrustum = false;
				break;
			}
		}

		drawWorldPos(intersection, vec4(inFrustum ? 0.0 : 1.0, inFrustum ? 1.0 : 0.0, 0.0, 1.0));

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

bool quadOutsideFrustum(vec3 points[4], bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]) {
	bool outside = false;
	for(int i = 0; i < 6; i++) {
		if(!whichPlanes[i]) {
			continue;
		}
		outside = pointsOutsideOfPlane(point, frustumPlaneNormals[i], getPointOnFrustumPlane(i, frustumPoints));
		if(outside) {
			break;
		}
	}

	return outside;
}

bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]) {
	// outside == on the side of the normal
	for (int i = 0; i < 4; i++) {
		if (!pointOutsideOfPlane(planePoint, planeNormal, points[i])) {
			return false;
		}
	}
	return true;
}

bool pointOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 point) {
	return dot(planeNormal, point - planePoint) >= 0;
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


vec4 normalPointPlaneToNormalDistPlane(vec3 normal, vec3 point) {
	return vec4(normal, dot(normal, -point));
}
#endif
