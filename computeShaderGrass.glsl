#version 450 core
#extension GL_NV_shader_thread_group : require


#include <renderer.glsl.h>
#include <random.glsl.h>

const float Epsilon = 0.001;

shared vec3 lastPosition;

void drawGrassBlade(RandState rng, vec3 pos, float cellSize);
float getStepSize(float dist, out float blend);
float distPointRay(vec3 origin, vec3 direction, vec3 point);
bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir);
vec3 floorGridPointByDir(vec3 point, vec3 dir, float stepSize);
vec3 floorGridPointByDirOffset(vec3 point, vec3 dir, float stepSize, vec2 gridOffset);
bool[6] negate(bool array[6]);
float intersectPlane(Ray ray, vec3 point, vec3 normal);
vec3 intersectFrustum(Ray ray, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
vec3 getPointOnFrustumPlane(int index, vec3 frustumPoints[8]);
bool gridCellOutsideFrustum(vec3 pos, float size, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
bool quadOutsideFrustum(vec3 points[4], bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]);
bool pointsOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 points[4]);
bool pointOutsideOfPlane(vec3 planePoint, vec3 planeNormal, vec3 point, float epsilon);
void drawWorldLine(vec3 start, vec3 end, vec4 colour);
void drawImageLine(vec2 imageStart, vec2 imageEnd, vec4 colour);
void drawWorldPos(vec3 pos, vec4 colour);
void drawNDC(vec3 ndc, vec4 colour);
vec2 worldPosToImagePos(vec3 pos);
vec2 ndcToImagePos(vec3 ndc);
vec3 worldPosToNDC(vec3 pos);

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

	if(grassConsts.DrawDebugInfo >= 2) {
		vec2 frustumCorner = vec2((gl_WorkGroupID.x) / tileCount.x * imageSize(result).x, (1 - (gl_WorkGroupID.y) / tileCount.y) * imageSize(result).y);
		vec2 tileSize = imageSize(result) / tileCount;

		drawImageLine(frustumCorner, frustumCorner + vec2(tileSize.x, 0.0), vec4(1.0, 1.0, 1.0, 1.0));
		drawImageLine(frustumCorner + vec2(tileSize.x, 0.0), frustumCorner + vec2(tileSize.x, -tileSize.y), vec4(1.0, 1.0, 1.0, 1.0));
		drawImageLine(frustumCorner + vec2(0.0, -tileSize.y), frustumCorner + vec2(tileSize.x, -tileSize.y), vec4(1.0, 1.0, 1.0, 1.0));
		drawImageLine(frustumCorner + vec2(0.0, -tileSize.y), frustumCorner, vec4(1.0, 1.0, 1.0, 1.0));
	}

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


	if(gl_GlobalInvocationID.x == 0) {
		RandState rng = rand_init(0, 0);
		//drawGrassBlade(rng, vec3(0), 1);
	}

	if(intersectionCount > 3 && gl_WorkGroupID.x == 1 && gl_WorkGroupID.y == 0 /*|| gl_WorkGroupID.x == 4 && gl_WorkGroupID.y == 4*/) {

		// Find line of cells to start with -> vector perpendicular to ftb vector at the start vector roughly pointing away from the camera (for front to back iteration later)
		vec3 lineDirection = grassConsts.PerpFtBDir;
		// If necessary mirror vector so it points roughly the right direction
		if(dot(lineDirection, camera.CamDir) < 0.0) {
			lineDirection = -lineDirection;
		}

		// debug
		vec3 p = vec3(0);
		for (int i = 0; i < 8; i++)
		{
			p += frustumPoints[i];
		}
		// debug

		float stepSize = grassConsts.Step;
		float blend;

		// Find start point for iteration (minimum by ftb)
		int index = 0;
		for(int i = 1; i < intersectionCount; i++) {

			if(lessThanByDir(intersections[i], intersections[index], grassConsts.FtBDirection)) {
				index = i;
			}
		}

		vec3 start = intersections[index];

		stepSize = grassConsts.Step; //getStepSize(1/*distPointRay(start, lineDirection, camera.CamPos)*/, blend);
		// Round to next cell position 
		if(grassConsts.FtBDirection.x != 0.0 && grassConsts.FtBDirection.z != 0) {
			// diagonal workaround
			vec3 s1 = floorGridPointByDirOffset(start, grassConsts.FtBDirection, stepSize, vec2(stepSize / 2));
			vec3 s2 = floorGridPointByDir(start, grassConsts.FtBDirection, stepSize);

			if(distance(start, s1) < distance(start, s2)) {
				start = s1;
			} else {
				start = s2;
			}
		} else {
			start = floorGridPointByDir(start, grassConsts.FtBDirection, stepSize);
		}

		drawWorldPos(start, vec4(0.0,0.0,1.0,1.0));

		// Find planes of the frustum that are on in the general direction of the lineDirection vector
		// Check these during iteration to find out if we're outside of the frustum
		bool frustumPlanesToCheck[6];
		bool frustumPlanesToCheckLine[6];
		bool frustumPlanesToCheckLineNegated[6];
		bool frustumPlanesToCheckFtB[6];
		for(int i = 0; i < 6; i++) {
			frustumPlanesToCheckLine[i] = dot(lineDirection, frustumPlaneNormals[i]) >= 0;
			frustumPlanesToCheckLineNegated[i] = !frustumPlanesToCheckLine[i];
			frustumPlanesToCheckFtB[i] = dot(grassConsts.PerpFtBDir, frustumPlaneNormals[i]) >= 0;
			frustumPlanesToCheck[i] = frustumPlanesToCheckLine[i] || frustumPlanesToCheckFtB[i];
			// drawNDC(vec3(0.5 + float(i) / 24, 0.5, 0.5), vec4(frustumPlanesToCheckLine[i] ? 0.0 : 1.0, frustumPlanesToCheckLine[i] ? 1.0 : 0.0, 0.0, 1.0));
		}


		if(grassConsts.DrawDebugInfo >= 1.0) {
			// debug - draw frustum intersections
			for(int i = 0; i < intersectionCount; i++) {
				if (i == index) {
					drawWorldPos(intersections[i], vec4(1.0, 1.0, 1.0, 1.0));
				} else {
					drawWorldPos(intersections[i], vec4(0.0, 1.0, 1.0, 1.0)); // Draw intersection points
				}
			}
		}

		// Iterate through cells
		int i = 0;
		int maxIt = 5120;
		int searchIt;
		int maxItPerSearch = 32;
		int invalidIntersectionSteps = 0;
		int maxInvalidIntersectionSteps = 5;

		vec3 ftbStep = grassConsts.FtBDirection * stepSize;
		vec3 horizontalStep;
		int lineNumber = 0;
		vec3 lineOneStart = start;
		vec3 lineTwoStart;
		vec3 lineStart;
		vec3 currentPos;
		bool outOfFrustum;
		bool searching;
		bool newLine = true;
		bool newGroup = true;
		bool done = false;
		bool maskedOut;
		int prevThreadsInLine;
		int threadsIn;
		uint activeBallot;
		uint localCompactId;
		float stepSizeLT;
		float localStepSize;
		RandState rng;

		for(int ba = 0;; ba++) {
			// Find suitable position for thread to draw at by checking if the currentPos is in the frustum and sharing that information with all threads
			outOfFrustum = true;
			searching = true;
			maskedOut = false;
			prevThreadsInLine = 0;
			lineNumber = 0;
			searchIt = 0;
			stepSize = getStepSize(distPointRay(lineOneStart, lineDirection, camera.CamPos), blend);

			while(searching) {
				if(newLine || newGroup) {
					// Calculate step size
					
					ftbStep = grassConsts.FtBDirection * stepSize;
					horizontalStep = lineDirection * stepSize;

					vec3 lT = lineOneStart;
					vec3 tempPos = lT;
					vec3 intersectionStep = vec3(0);
					invalidIntersectionSteps = 0;
					do {
						invalidIntersectionSteps++;

						intersectionStep += ftbStep;

						// Intersect with frustum to find next line
						tempPos = lineOneStart + intersectionStep;
						lT = intersectFrustum(Ray(tempPos, -lineDirection), frustumPlanesToCheckLineNegated, frustumPlaneNormals, frustumPoints);
						
						if(invalidIntersectionSteps > maxInvalidIntersectionSteps) {
							break;
						}
					} while(any(isinf(lT))); // search for a valid intersection
					stepSizeLT = getStepSize(distPointRay(lT - (invalidIntersectionSteps - 1) * ftbStep, lineDirection, camera.CamPos), blend);
					lineTwoStart = floorGridPointByDir(lT - (invalidIntersectionSteps - 1) * ftbStep, lineDirection, stepSizeLT);
					
					
					if(newLine) {
						if(lessThanByDir(lT, lineOneStart, horizontalStep)) {
							lineStart = floorGridPointByDir(lT - intersectionStep, lineDirection, stepSize);
						} else {
							lineStart = lineOneStart;
						}
					} else {
						lineStart = lineOneStart;
					}
					
					drawWorldPos(lT, vec4(1.0, 1.0, 1.0, 1.0));

					newLine = false;
					newGroup = false;
				}
				// drawWorldPos(lineTwoStart, vec4(0.0,0.0,1.0,1.0));
				//drawWorldPos(start, vec4(0.0,0.0,1.0,1.0));
				//drawWorldPos(lineStart, vec4(1.0, 1.0, 0.0, 1.0));

				// Check if position is in frustum and share with other threads
				if(!any(isinf(lineStart))) {
					// Find thread position
					uint localCompactId = 0;
					activeBallot = ballotThreadNV(true);
					if(gl_LocalInvocationID.x > 0) {
						localCompactId = bitCount(activeBallot << (32 - gl_LocalInvocationID.x));
					}

					currentPos = lineStart + (localCompactId + prevThreadsInLine) * horizontalStep;
					
					ivec2 iPos = ivec2(round((currentPos.xz / stepSize)));
					ivec2 seedPos = ivec2(round(currentPos.xz / grassConsts.Step));


					// Check if in frustum
					localStepSize = getStepSize(distance(currentPos, camera.CamPos), blend);
					outOfFrustum = gridCellOutsideFrustum(currentPos, stepSize, frustumPlanesToCheckLine, frustumPlaneNormals, frustumPoints);

					
					// drawWorldPos(currentPos, vec4(outOfFrustum ? 1.0 : 0.0, outOfFrustum ? 0.0 : 1.0, 0.0, 1.0));

					
					// mask out cells if the stepSize used is too small for the cell
					// maskedOut = localStepSize > stepSize;
					maskedOut = maskedOut && ((iPos.x | iPos.y) & (int(round(localStepSize / stepSize) - 1))) != 0;

					// create blend between cells of different stepSizes by masking points out depending on how far they are from their optimal step size 
					if(!maskedOut) {
						iPos = ivec2(round((currentPos.xz / localStepSize)));
						rng = rand_init(seedPos.x, seedPos.y);
						// maskedOut = rand_next(rng) < blend;
						maskedOut = maskedOut && ((iPos.x | iPos.y) & 1) != 0;
					}
					
				} else {
					outOfFrustum = true;
				}
			
				
				if(prevThreadsInLine == 0 && ballotThreadNV(outOfFrustum) == ballotThreadNV(true) && lineNumber > 0) {
					done = true;
					break;
				}

				if(!outOfFrustum){
					prevThreadsInLine += bitCount(ballotThreadNV(true));
				}
				searching = maskedOut || outOfFrustum;
			
				// thread not in frustum, go to next line
				if(ballotThreadNV(outOfFrustum) != 0) {
					// next line
					prevThreadsInLine = 0;
					newLine = true;
					lineOneStart = lineTwoStart;
					stepSize = stepSizeLT;
					lineNumber++;
				}
				
				if (gl_LocalInvocationID.x == findMSB(ballotThreadNV(!searching))) {
					// last thread found position to draw at
					lastPosition = currentPos + horizontalStep;
				}

				i++;
				searchIt++;
				if (i >= maxIt/* || searchIt > maxItPerSearch*/) {
					if(grassConsts.DrawDebugInfo >= 3) drawWorldPos(p / 8, vec4(1.0, 1.0, 0.0, 1.0));

					done = true;
					break;
				}
			}

			if(done) {
				break;
			}

			// found position -> draw
			vec3 quad[4] = {currentPos, currentPos + vec3(0.0, 0.0, localStepSize), currentPos + vec3(localStepSize, 0.0, localStepSize),
					currentPos + vec3(localStepSize, 0.0, 0.0)};
			for(int l = 0; l < 4; l++) {
				//drawWorldLine(quad[l], quad[(l + 1) % 4], vec4(1.0 - float(gl_LocalInvocationID.x) / 31, 0.0, float(gl_LocalInvocationID.x) / 31, 1.0));
			}
			// drawWorldPos(currentPos, vec4(1.0 - float(gl_LocalInvocationID.x) / 31, 0.0, float(gl_LocalInvocationID.x) / 31, 1.0));
			// drawWorldPos(currentPos, vec4(/*localStepSize > stepSize ? 0.0 :*/ 1.0, 1.0, 1.0, 1.0));
			// drawGrassBlade(rng, currentPos, localStepSize);

			// get position of last thread to find the next cells to draw at
			memoryBarrierShared();
			lineOneStart = lastPosition;
			drawWorldPos(lineOneStart, vec4(1.0, 1.0, 0.0, 1.0));
			newLine = false;
			newGroup = true;
		}
	}
}


void drawGrassBlade(RandState rng, vec3 pos, float cellSize) {
	int curveSubdivision = 10;

	float minHeight = 0.2;
	float maxHeight = 2.0;
	float maxHorizontalControlPointDerivation = 0.7;
	float maxHorizontalTipDerivation = 0.7;

	// Random generation of grass blade properties
	float cPHeight = minHeight + rand_next(rng) * (maxHeight - minHeight);
	float tipHeight = minHeight + rand_next(rng) * (maxHeight - minHeight);
	

	vec3 root = pos + vec3((rand_next(rng) * cellSize), 0, (rand_next(rng) * cellSize));
	vec3 controlPoint = root + vec3((rand_next(rng) * maxHorizontalControlPointDerivation),
		cPHeight,
		(rand_next(rng) * maxHorizontalControlPointDerivation));
	vec3 tip = root + vec3((rand_next(rng) * maxHorizontalControlPointDerivation),
		tipHeight,
		(rand_next(rng) * maxHorizontalControlPointDerivation));


	vec2 rootProj = worldPosToImagePos(root);
	vec2 cPProj = worldPosToImagePos(controlPoint);
	vec2 tipProj = worldPosToImagePos(tip);

	if(rootProj.x < 0) {
		return;
	}
	// subdivide curve and draw line segments
	vec2 p0 = rootProj;
	vec2 p1;
	
	for(int i = 1; i < curveSubdivision; i++) {
		float u = float(i) / curveSubdivision;
		p1 = pow(1 - u, 2) * rootProj + 2 * (1 - u) * u * cPProj + pow(u, 2) * tipProj;

		drawImageLine(p0, p1, vec4(0.1, 1.0, 0.2, 1.0));
		p0 = p1;
	}
	drawImageLine(p1, tipProj, vec4(0.1, 1.0, 0.2, 1.0));
}


// Get step size for regular grid by distance to the camera
float getStepSize(float dist, out float blend) {
	float factor = dist / grassConsts.StepDist; //grassConsts.StepDoubleDist;
	
	if(factor <= 1) {
		blend = 0;
		factor = 1;
	} else {
		blend = modf(log2(factor), factor);
		factor = pow(2, factor);
	}
	return factor * grassConsts.Step;
}

float distPointRay(vec3 origin, vec3 direction, vec3 point) {
	return length(cross(normalize(direction), point - origin));
}

bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir) {
	bool result = false;
	if(any(isinf(fst)) || any(isinf(snd))) {
		result = false;
	} else {
		if(dot(snd - fst, dir) > 0) {
			result = true;
		} 
	}
	return result;
}

vec3 floorGridPointByDir(vec3 point, vec3 dir, float stepSize) {
	return floorGridPointByDirOffset(point, dir, stepSize, vec2(0,0));
}

vec3 floorGridPointByDirOffset(vec3 point, vec3 dir, float stepSize, vec2 gridOffset) {
	vec3 result = vec3(0.0);
	point = vec3(point.x + gridOffset.x, point.y, point.z + gridOffset.y);

	if(any(isinf(point))) {
		result = vec3(1.0 / 0.0);
	} else {
		// floors point to grid point in the opposite direction of dir 
		if(dir.x > 0) {
			result.x = floor(point.x / stepSize) * stepSize;
		} else if (dir.x == 0) {
			result.x = round(point.x / stepSize) * stepSize;
		} else {
			result.x = ceil(point.x / stepSize) * stepSize;
		}
		if(dir.z > 0) {
			result.z = floor(point.z / stepSize) * stepSize;
		} else if(dir.z == 0) {
			result.z = round(point.z / stepSize) * stepSize;
		} else {
			result.z = ceil(point.z / stepSize) * stepSize;
		}
	}

	return result - vec3(gridOffset.x, 0.0, gridOffset.y);
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
		if(isinf(t)) {
			inFrustum = false;
		} else {
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

bool gridCellOutsideFrustum(vec3 pos, float size, bool whichPlanes[6], vec3 frustumPlaneNormals[6], vec3 frustumPoints[8]) {
	vec3 quad[4] = {pos, pos + vec3(0.0, 0.0, size),
					pos + vec3(size, 0.0, 0.0), pos + vec3(size, 0.0, size)};
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
	vec3 ndc = worldPosToNDC(pos);
	drawNDC(ndc, colour);
}

void drawNDC(vec3 ndc, vec4 colour) {
	ivec2 imagePos = ivec2(ndcToImagePos(ndc));
	if (-1.0 <= ndc.x && ndc.x <= 1.0 &&
		-1.0 <= ndc.y && ndc.y <= 1.0 &&
		ndc.z >= 0.0 && ndc.z <= 1.0) {

		imageStore(result, imagePos, colour);
		// imageStore(result, imagePos + ivec2(0, 1), colour);
		// imageStore(result, imagePos + ivec2(1, 0), colour);
		// imageStore(result, imagePos + ivec2(1, 1), colour);
	}
}

void drawWorldLine(vec3 start, vec3 end, vec4 colour) {
	drawImageLine(worldPosToImagePos(start), worldPosToImagePos(end), colour);
}

void drawImageLine(vec2 imageStart, vec2 imageEnd, vec4 colour) {
	// Digital Differential Algorithm	
	float dx = imageEnd.x - imageStart.x;
	float dy = imageEnd.y - imageStart.y;
	float steps;
	if(abs(dx) > abs(dy)) {
		steps = abs(dx);
	} else {
		steps = abs(dy);
	}
	vec2 increment = vec2(dx / steps, dy / steps);
	vec2 pos = vec2(imageStart);
	for(float v = 0; v < steps; v++) {
		pos += increment;
		imageStore(result, ivec2(pos), colour);
	}
}

vec2 worldPosToImagePos(vec3 pos) {
	vec3 ndc = worldPosToNDC(pos);
	if(ndc.z < 0 || ndc.z > 1) {
		return vec2(-1);
	} else {
		return ndcToImagePos(ndc);
	}
}

vec2 ndcToImagePos(vec3 ndc) {
	vec2 screenCoords = (ndc.xy + vec2(1.0, 1.0)) * 1/2;
	vec2 texCoords = vec2(screenCoords.x, 1.0 - screenCoords.y);
	return vec2(texCoords * imageSize(result));
}

vec3 worldPosToNDC(vec3 pos) {
	vec4 p1 = camera.ViewProj * vec4(pos, 1.0);
	return p1.xyz / p1.w;
}

vec4 normalPointPlaneToNormalDistPlane(vec3 normal, vec3 point) {
	return vec4(normal, dot(normal, -point));
}
#endif
