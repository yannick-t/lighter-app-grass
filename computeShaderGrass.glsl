#version 450 core
#extension GL_NV_shader_thread_group : require


#include <renderer.glsl.h>
#include <random.glsl.h>

#define PI 3.1415926535897932384626433832795

const float Epsilon = 0.001;
const float shininess = 4.0; // Todo: maybe randomize

shared vec3 lastPosition;

shared vec4 pixels[32][32];

void drawGrassBlade(vec3 pos, float stepSize);
bool drawGrassBladePixel(float y, float t);
bool drawGrassBladePixelBlend(int x, int y, vec4 color);
vec4 getGrassBladeShading(vec3 pos, vec3 normal);

void scanlineRasterizeGrassBlade();

void calcWindTranslation(RandState rng);

vec3[6] calcFrustumNormals();
Ray[12] calcFrustumRays();
vec3[4] calcFrustumIntersections(vec3 planePoint, vec3 planeNormal, out int intersectionCount);

float getStepSize(float dist, out float blend);
float distPointRay(vec3 origin, vec3 direction, vec3 point);
float distPointLineSegment(vec3 v, vec3 w, vec3 p);
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

void drawTilePos(vec2 pos, vec4 color);
void drawWorldLine(vec3 start, vec3 end, vec4 color);
void drawImageLine(vec2 imageStart, vec2 imageEnd, vec4 color);
void drawWorldPos(vec3 pos, vec4 color);
void drawNDC(vec3 ndc, vec4 color);

vec2 worldPosToTilePos(vec3 pos);
vec2 worldPosToImagePos(vec3 pos);
vec2 ndcToImagePos(vec3 ndc);
vec3 worldPosToNDC(vec3 pos);
vec3 nDCToWorldPos(vec3 ndc);
vec4 normalPointPlaneToNormalDistPlane(vec3 normal, vec3 point);

bool epsilonEq(float arg1, float arg2);
mat4 rotationMatrix(vec3 axis, float angle);

void plotQuadBezier(int x0, int y0, int x1, int y1, int x2, int y2);
void plotQuadBezierSegAA(int x0, int y0, int x1, int y1, int x2, int y2);
void plotLineAA(int x0, int y0, int x1, int y1);

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
layout(binding = 4) uniform atomic_uint bladeCount;

#ifdef IN_CS
layout(local_size_x = 32) in;
layout(rgba32f, binding = 0) uniform image2D result;

// properties of workgroup tile
ivec2 tileCorner;

vec3 frustumPointsNDC[8]; // nbl, ntl, nbr, ntr, fbl, ftl, fbr, ftr
vec3 frustumPoints[8]; // nbl, ntl, nbr, ntr, fbl, ftl, fbr, ftr
vec3 frustumPlaneNormals[6]; // right, left, top, bottom, near, far
Ray frustumRays[12];
int frustumIntersectionCount = 0;
vec3 frustumIntersections[4];

// properties of grid
float initialStepSize;
float initialStepSizeNDC;

// properties of grass blades
vec3 normal;
vec4 baseColor;

float aoDist;

float width;
float widthPx;
float tipLengthT;

vec3 root;
vec3 cP;
vec3 tip;

vec2 rootProj;
vec2 cPProj;
vec2 tipProj;

mat4 normalPerpRotation;

void main()
{
	// initialize pixel array
	for(int i = 0; i < 32; i++) {
		pixels[gl_LocalInvocationID.x][i] = vec4(0.0);
	}

	vec2 tileCount = vec2(imageSize(result) / float(grassConsts.TileDivisor));
	ivec2 tileSize = ivec2(grassConsts.TileDivisor);
	tileCorner = ivec2(gl_WorkGroupID.x * tileSize.x, gl_WorkGroupID.y * tileSize.y);


	// Calc initial step size
	initialStepSizeNDC = (grassConsts.StepPxAtMinDist / camera.Resolution.x) * 2;
	vec3 referencePos = camera.CamPos + camera.CamDir * grassConsts.MinDist;
	initialStepSize = distance(referencePos, nDCToWorldPos(worldPosToNDC(referencePos) + vec3(1.0, 0.0, 0.0) * initialStepSizeNDC));

	// Calculate frustum of this workgroup
	vec2 minFrustumNDC = vec2((1.0 / 0.0));
	vec2 maxFrustumNDC = vec2(-(1.0 / 0.0));

	int v = 0;
	for(int i = 0; i < 8; i++) {
		frustumPointsNDC[i] = vec3(((gl_WorkGroupID.x + int(v / 2)) / tileCount.x) * 2 - 1, -(((gl_WorkGroupID.y + mod(v + 1, 2)) / tileCount.y) * 2 - 1), int(i / 4));
		minFrustumNDC = min(frustumPointsNDC[i].xy, minFrustumNDC);
		maxFrustumNDC = max(frustumPointsNDC[i].xy, maxFrustumNDC);
		frustumPoints[i] = nDCToWorldPos(frustumPointsNDC[i]);
		v = int(mod(v + 1, 4));
	}

	frustumRays = calcFrustumRays();

	// intersect frustum rays with plane on which grass should be drawn
	vec3 planePoint = vec3(0.0);
	vec3 planeNormal = vec3(0.0, 1.0, 0.0);

	frustumIntersections = calcFrustumIntersections(planePoint, planeNormal, frustumIntersectionCount);

	// only draw tiles containing grass blades further than the min dist
	bool draw = false;
	for(int i = 0; i < frustumIntersectionCount; i++) {
		draw = draw || distance(frustumIntersections[i], camera.CamPos) > grassConsts.MinDist;
	}

	draw = draw && frustumIntersectionCount > 3;
	
	
	if(draw /*&& (gl_WorkGroupID.x == 3 && gl_WorkGroupID.y == 10 ||gl_WorkGroupID.x == 1 && gl_WorkGroupID.y == 12)*/) {

		
		
		// expand frustum to contain grass blades where the cell isn't in the frustum but a part of the blade can be
		// Calculate intersections with the plane of max grass height

		int secondPlaneIntersectionCount = 0;
		vec3[4] secondPlaneIntersections = calcFrustumIntersections(vec3(0, grassConsts.MaxHeight, 0), planeNormal, secondPlaneIntersectionCount);	
		// move intersections to xz-plane and calc min / max in ndc with old intersections
		vec3 ndcMIntersecion;
		for(int i = 0; i < secondPlaneIntersectionCount; i++) {
			ndcMIntersecion = worldPosToNDC(vec3(secondPlaneIntersections[i].x, 0, secondPlaneIntersections[i].z));
			minFrustumNDC = min(ndcMIntersecion.xy, minFrustumNDC);
			maxFrustumNDC = max(ndcMIntersecion.xy, maxFrustumNDC);
		}

		// Include one cell to the left and to the right in case a grass blade grows into this frustum from there ?
		//minFrustumNDC.x -= initialStepSizeNDC;
		//maxFrustumNDC.x += initialStepSizeNDC;

		// Expand frustum with the min / max values (essentially an AABB)
		for(int i = 0; i < 8; i+= 4) {
			frustumPointsNDC[i] = vec3(minFrustumNDC.x, minFrustumNDC.y, frustumPointsNDC[i].z);
			frustumPointsNDC[i + 1] = vec3(minFrustumNDC.x, maxFrustumNDC.y, frustumPointsNDC[i + 1].z);
			frustumPointsNDC[i + 2] = vec3(maxFrustumNDC.x, minFrustumNDC.y, frustumPointsNDC[i + 2].z);
			frustumPointsNDC[i + 3] = vec3(maxFrustumNDC.x, maxFrustumNDC.y, frustumPointsNDC[i + 3].z);
		}
		// Recalc Points
		for(int i = 0; i < 8; i++) {
			frustumPoints[i] = nDCToWorldPos(frustumPointsNDC[i]);
		}


		// Recalc frustum
		frustumRays = calcFrustumRays();
		frustumIntersections = calcFrustumIntersections(planePoint, planeNormal, frustumIntersectionCount);
		frustumPlaneNormals = calcFrustumNormals();


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

		float stepSize = initialStepSize;
		float blend;

		// Find minimum intersection(s) by ftb 
		int index = 0;
		for(int i = 1; i < frustumIntersectionCount; i++) {

			if(lessThanByDir(frustumIntersections[i], frustumIntersections[index], grassConsts.FtBDirection)) {
				index = i;
			}
		}


		// Declare start point for iteration
		vec3 start = frustumIntersections[index];

		// stepSize = getStepSize(distance(start, camera.CamPos), blend);
		stepSize = getStepSize(distPointRay(start, lineDirection, camera.CamPos), blend);
		// Round to next cell position 
		if(!epsilonEq(grassConsts.FtBDirection.x, 0.0) && !epsilonEq(grassConsts.FtBDirection.z, 0.0)) {
			// diagonal workaround
			vec3 s1 = floorGridPointByDirOffset(start, grassConsts.FtBDirection, stepSize, vec2(stepSize / 2));
			vec3 s2 = floorGridPointByDir(start, grassConsts.FtBDirection, stepSize);

			if(distance(start, s1) < distance(start, s2)) {
				start = s1;
				start = floorGridPointByDir(start, lineDirection, stepSize);
			} else {
				start = s2;
			}
		} else {
			start = floorGridPointByDir(start, grassConsts.FtBDirection, stepSize);
		}

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
		}

		// Iterate through cells
		int i = 0;
		int maxIt = 500;
		int searchIt;
		int maxItPerSearch = 32;
		int invalidIntersectionSteps = 0;
		int maxInvalidIntersectionSteps = 2;

		float floorStepSizeLT;
		vec3 ftbStep = grassConsts.FtBDirection * stepSize;
		vec3 horizontalStep;
		int lineNumber = 0;

		vec3 lineOneStart = start;
		vec3 lineTwoStart;
		vec3 lO;
		vec3 lT;
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
		float rand;
		float alpha;
		RandState rng;

		int debugTest = 0;

		for(int ba = 0; ; ba++) {
			// Find suitable position for thread to draw at by checking if the currentPos is in the frustum and sharing that information with all threads
			outOfFrustum = true;
			searching = true;
			maskedOut = false;
			prevThreadsInLine = 0;
			lineNumber = 0;
			searchIt = 0;
			stepSize = getStepSize(distPointRay(lineOneStart, lineDirection, camera.CamPos), blend);
			lO = lineOneStart;
			// vec3 lOE = intersectFrustum(Ray(lineOneStart, lineDirection), frustumPlanesToCheckLine, frustumPlaneNormals, frustumPoints);
			//stepSize = getStepSize(distPointLineSegment(lineOneStart, lOE, camera.CamPos), blend);

			while(searching) {
				alpha = 1.0;

				if(newLine || newGroup) {
					// Calculate steps
					ftbStep = grassConsts.FtBDirection * stepSize;
					horizontalStep = lineDirection * stepSize;

					vec3 tempPos = lineOneStart;
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
					lT = lT - (invalidIntersectionSteps - 1) * ftbStep;
					// vec3 lTE = intersectFrustum(Ray(lT, lineDirection), frustumPlanesToCheckLine, frustumPlaneNormals, frustumPoints);
					// stepSizeLT = getStepSize(distPointLineSegment(lT, lTE, camera.CamPos), blend);
					stepSizeLT = getStepSize(distPointRay(lT, lineDirection, camera.CamPos), blend);
					lineTwoStart = floorGridPointByDir(lT, lineDirection, stepSizeLT);
					
					
					if(newLine) {
						if(lessThanByDir(lT, lO, horizontalStep)) {
							lineOneStart = floorGridPointByDir(lT - intersectionStep, lineDirection, stepSize);
						} else if(lessThanByDir(lO, lT, horizontalStep)) {
							lineTwoStart = floorGridPointByDir(lO + intersectionStep, lineDirection, stepSizeLT);
						}
					}
					lineStart = lineOneStart; //floorGridPointByDir(lineOneStart, lineDirection, stepSize);
					
					
					//drawTilePos(worldPosToTilePos(lineStart) + vec2(10, 0), vec4(1.0, 0.0, 1.0, 1.0));

					//drawWorldPos(lineOneStart, vec4(0.0,0.0,1.0,1.0));
					//drawWorldPos(lineTwoStart, vec4(0.0,1.0,1.0,1.0));

					newLine = false;
					newGroup = false;
				}
				
				// drawWorldPos(lineTwoStart, vec4(0.0,0.0,1.0,1.0));
				// drawWorldPos(start, vec4(0.0,0.0,1.0,1.0));
				// drawWorldPos(lineStart, vec4(stepSize <= initialStepSize ? 1.0 : 0.0, 1.0, 0.0, 1.0));

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
					ivec2 seedPos = ivec2(round(currentPos.xz / initialStepSize));


					// Check if in frustum
					float camDist = distance(currentPos, camera.CamPos);
					localStepSize = getStepSize(camDist, blend);
					outOfFrustum = gridCellOutsideFrustum(currentPos, localStepSize, frustumPlanesToCheckLine, frustumPlaneNormals, frustumPoints);

					
					// drawWorldPos(currentPos, vec4(outOfFrustum ? 1.0 : 0.0, outOfFrustum ? 0.0 : 1.0, 0.0, 1.0));
					

					// begin at MinDist
					maskedOut = camDist < grassConsts.MinDist; // || camDist > grassConsts.MaxDist;
					
					// mask out cells if the stepSize used is too small for the cell
					if(!maskedOut) {
						maskedOut = localStepSize > stepSize;
						maskedOut = maskedOut && ((iPos.x | iPos.y) & (int(round(localStepSize / stepSize) - 1))) != 0;
					}

					
					
					// create blend between cells of different stepSizes by masking points out depending on how far they are from their optimal step size 
					if(!maskedOut) {
						// work with next bigger step size 
						// mask out 3 cells at a time to create one cell of the next step size (expand the 4th one)
						float nextStepSize = localStepSize * 2;

						// decide with the position of the bigger cell if the cell will be masked out
						vec3 relativePos = currentPos / nextStepSize;
						vec3 nextBiggerCellPos = floor(relativePos + vec3(Epsilon)) * nextStepSize;
						seedPos = ivec2(round(nextBiggerCellPos.xz / initialStepSize));

						// get correct blend
						float cornerCellStepSize = getStepSize(distance(nextBiggerCellPos, camera.CamPos), blend);

						// pull random number
						rng = rand_init(seedPos.x, seedPos.y);
						rand = rand_next(rng);

						bool makeBigger = 
							(rand < blend ||									// randomly mask out, the closer to the next step size the more probable to be masked out
							epsilonEq(cornerCellStepSize, nextStepSize))		// If the corner of the bigger cell is already in the next step size mask out (workaround because blend value is off in this case)
							&& cornerCellStepSize >= localStepSize - Epsilon;	// If the corner of the bigger cell is in the next smaller step size don't mask out (workaround because blend value is off in this case)
						if(makeBigger && epsilonEq(currentPos.x, nextBiggerCellPos.x) && epsilonEq(currentPos.z, nextBiggerCellPos.z)) {
							// bottom left corner of the bigger cell
							localStepSize = nextStepSize;
						} else {
							maskedOut = makeBigger;
							// alpha = (rand - blend) / (1 - blend); // Calculate an alpha value representing how close a cell is to being masked out
						}
						
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
					lO = lT;
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
					//drawWorldPos(p / 8, vec4(1.0, 1.0, 0.0, 1.0));

					done = true;
					break;
				}
			}

			if(done) {
				break;
			}

			// found position -> draw
			
			//vec3 quad[4] = {currentPos, currentPos + vec3(0.0, 0.0, localStepSize), currentPos + vec3(localStepSize, 0.0, localStepSize),
					//currentPos + vec3(localStepSize, 0.0, 0.0)};
			//for(int l = 0; l < 4; l++) {
				// drawWorldLine(quad[l], quad[(l + 1) % 4], vec4(1.0 - float(gl_LocalInvocationID.x) / 31, 0.0, float(gl_LocalInvocationID.x) / 31, 1.0));
				// drawWorldLine(quad[l], quad[(l + 1) % 4], vec4(1.0, 1.0, 1.0, 1.0));
			//}
			//drawWorldPos(currentPos, vec4(1.0 - float(gl_LocalInvocationID.x) / 31, 0.0, float(gl_LocalInvocationID.x) / 31, 1.0));
			// drawWorldPos(currentPos + vec3(rand_next(rng) * localStepSize, 0.0, rand_next(rng) * localStepSize) , vec4(1.0, 1.0, 1.0, 1.0));
			//drawWorldPos(vec3(currentPos.x + localStepSize, 0, currentPos.z + localStepSize), vec4(0.9, 0.9, 0.8, 1.0));
			// drawWorldPos(currentPos, vec4(alpha, alpha, alpha, 1));
			drawGrassBlade(currentPos, localStepSize);
			// vec4 c = vec4(0.0,0.1,0.2,1);
			//drawTilePos(worldPosToTilePos(vec3(currentPos.x + localStepSize, 0, currentPos.z + localStepSize)), c);
			// drawTilePos(worldPosToTilePos(vec3(currentPos.x + 0.01 * stepSize, 0, currentPos.z + 0.01 * stepSize)), vec4(1,1,0,1));
			// drawTilePos(worldPosToTilePos(vec3(currentPos.x, 0, currentPos.z)), c);
			// drawTilePos(worldPosToTilePos(currentPos) + vec2(1, 0), c);
			// drawTilePos(worldPosToTilePos(currentPos) + vec2(0, 1), c);
			// drawTilePos(worldPosToTilePos(currentPos) + vec2(1, 1), c);
			// drawTilePos(vec2(0), vec4(0, 0, 0, 1));
			// drawTilePos(vec2(31), vec4(0, 0, 0, 1));


			// get position of last thread to find the next cells to draw at
			memoryBarrierShared();
			lineOneStart = lastPosition;
			newLine = false;
			newGroup = true;
		}
	}


	// Copy workgroup result to framebuffer
	for(int i = 0; i < 32; i++) {
		imageStore(result, tileCorner + ivec2(gl_LocalInvocationID.x, i), pixels[gl_LocalInvocationID.x][i]);
	}

	if(grassConsts.DrawDebugInfo >= 2) {
		drawImageLine(vec2(tileCorner), vec2(tileCorner) + vec2(tileSize.x, 0.0), vec4(1.0, 1.0, 1.0, 1.0));
		drawImageLine(vec2(tileCorner) + vec2(tileSize.x, 0.0), vec2(tileCorner) + vec2(tileSize.x, tileSize.y), vec4(1.0, 1.0, 1.0, 1.0));
		drawImageLine(vec2(tileCorner) + vec2(0.0, tileSize.y), vec2(tileCorner) + vec2(tileSize.x, tileSize.y), vec4(1.0, 1.0, 1.0, 1.0));
		drawImageLine(vec2(tileCorner) + vec2(0.0, tileSize.y), vec2(tileCorner), vec4(1.0, 1.0, 1.0, 1.0));
	}

	if(grassConsts.DrawDebugInfo >= 1.0) {
		// debug - draw frustum intersections
		for(int i = 0; i < frustumIntersectionCount; i++) {
			drawWorldPos(frustumIntersections[i], vec4(0.0, 0.0, 0.0, 1.0)); // Draw intersection points
		}
	}
}

// functions for drawing the grass blade
void drawGrassBlade(vec3 pos, float stepSize) {
	atomicCounterIncrement(bladeCount);
	
	// Calculate a random stable normal sized cell if the cell is bigger than the normal one
	// do this by finding a random one of the smaller cells with half the step size until the stepSize is the initial one
	
	ivec2 seedPos = ivec2(round(pos.xz / initialStepSize));
	RandState rng = rand_init(seedPos.x, seedPos.y); 
	float halfStep;
	while(initialStepSize < stepSize) {
		halfStep = stepSize / 2;
		pos = pos + vec3(floor(rand_next(rng) * stepSize / halfStep) * halfStep, 0.0, floor(rand_next(rng) * stepSize / halfStep) * halfStep);
		seedPos = ivec2(round(pos.xz / initialStepSize));
		rng = rand_init(seedPos.x, seedPos.y);

		stepSize = halfStep;
	}

	float minHeight = grassConsts.MinHeight;
	float maxHeight = grassConsts.MaxHeight;
	float maxHorizontalControlPointDerivation = 1;
	float maxHorizontalTipDerivation = 1;

	// Random generation of grass blade properties
	float cPHeight = (minHeight + rand_next(rng) * (maxHeight - minHeight));
	float tipHeight = cPHeight + rand_next(rng) * (maxHeight - cPHeight);

	aoDist = grassConsts.RelAODist * tipHeight;

	width = grassConsts.MinWidth + rand_next(rng) * (grassConsts.MaxWidth - grassConsts.MinWidth);
	tipLengthT = 0.4 + rand_next(rng) * (1 - 0.4);
	
	root = pos + vec3((rand_next(rng) * initialStepSize), 0, (rand_next(rng) * initialStepSize));
	cP = pos + vec3((rand_next(rng) * maxHorizontalControlPointDerivation * initialStepSize),
		cPHeight,
		(rand_next(rng) * maxHorizontalControlPointDerivation * initialStepSize));
	tip = pos + vec3((rand_next(rng) * maxHorizontalControlPointDerivation * initialStepSize),
		tipHeight,
		(rand_next(rng) * maxHorizontalControlPointDerivation * initialStepSize));

	calcWindTranslation(rng);

	rootProj = worldPosToTilePos(root);
	cPProj = worldPosToTilePos(cP);
	tipProj = worldPosToTilePos(tip);

	//drawTilePos(rootProj, vec4(0.5,0.25,0.2,1.0));
	//drawTilePos(cPProj, vec4(0.85,0.5,0.0,1.0));
	//drawTilePos(tipProj, vec4(0.3,0.85,0.0,1.0));
	// drawWorldPos(root, vec4(1.0,1.0,1.0,1.0));
	//drawImageLine(rootProj, cPProj, vec4(0.5,0.5,0.5,1.0));
	//drawImageLine(cPProj, tipProj, vec4(0.5,0.5,0.5,1.0));
	
	// Generate random color
	baseColor = vec4(rand_next(rng) * 0.2, rand_next(rng) * 1.0, rand_next(rng) * 0.1, 1.0);
	// Generate normal pointing towards the controlPoint
	normal = normalize(cP - (root + tip) / 2);
	vec3 normalPerpRotationAxis = normalize(cross(root - cP, tip - cP));
	normalPerpRotation = rotationMatrix(cross(root - cP, tip - cP), 3 / 2 * PI);
	widthPx = round(distance(worldPosToImagePos(root), worldPosToImagePos(root + normalPerpRotationAxis * width)));

	// scanline rasterization
	scanlineRasterizeGrassBlade();
}

void scanlineRasterizeGrassBlade() {
	float start = 0; // min(min(min(tipProj.y, cPProj.y), rootProj.y), 31); // Todo: maybe better start
	float end = 31; //max(max(max(tipProj.y, cPProj.y), rootProj.y), 0);

	bool cont = true;
	for (float y = start; y <= end; y++) {
		// intersect scanline with curve
		if (epsilonEq(cPProj.y, (rootProj.y + tipProj.y) / 2) && rootProj.y - tipProj.y != 0) {
			float t = (rootProj.y - y) / (rootProj.y - tipProj.y);
			if (t >= 0 && t <= 1) {
				cont = drawGrassBladePixel(y, t);
			}
		} else {
			float a = (2 * cPProj.y - rootProj.y - tipProj.y);
			float det = pow(cPProj.y, 2) - 2 * cPProj.y * y - rootProj.y * tipProj.y + rootProj.y * y + tipProj.y * y;
			if (det >= 0 && a != 0) {
				float t1 = (sqrt(det) + cPProj.y - rootProj.y) / a;
				float t2 = (-sqrt(det) + cPProj.y - rootProj.y) / a;;

				if (t1 >= 0 && t1 <= 1) {
					cont = drawGrassBladePixel(y, t1);
				}

				if (t2 >= 0 && t2 <= 1) {
					cont = drawGrassBladePixel(y, t2);
				}
			}
		}

		if(!cont) break;
	}
}

// function to be called by the scanline rasterization to draw a pixel of a grass blade
// draws grass blade pixels and applies shading and AA
bool drawGrassBladePixel(float y, float t) {
	// Calculate pixel properties
	float x = pow(1 - t, 2) * rootProj.x + 2 * (1 - t) * t * cPProj.x + pow(t, 2) * tipProj.x;
	vec3 pos = pow(1 - t, 2) * root + 2 * (1 - t) * t * cP + pow(t, 2) * tip;
	vec3 tangent = 2 * (-2 * cP * t + cP + root * (t - 1) + t * tip);
	vec3 n = normalize((normalPerpRotation * vec4(tangent, 1.0)).xyz);

	if (x < 0 || x >= 32) {
		return true;
	}

	// get shading
	vec4 color = getGrassBladeShading(pos, n);

	// draw
	int ix;
	int iy = int(y);
	int xStart = int(max(x, 0));
	float widthScaling = 1 - ((t - (1 - tipLengthT)) / tipLengthT);
	float fxEnd = widthScaling > 1 ? x + widthPx : x + widthScaling * widthPx; 
	int xEnd = int(min(fxEnd, 31));

	bool continueDrawing = false;
	for(ix = xStart; ix <= xEnd; ix++) {
		continueDrawing = drawGrassBladePixelBlend(ix, iy, color) || continueDrawing;
	}

	// Antialiasing - in +x and -x because of the width scaling it can be different
	float b = modf(x, x);
	int aaPixel = xStart - 1;
	if(b < 0.5) {
		b = abs((b - 0.5) * 2);

		color.a = b;
		if(aaPixel >= 0 && aaPixel <= 31) {
			continueDrawing = drawGrassBladePixelBlend(aaPixel, iy, color) || continueDrawing;
		}
	}
	b = modf(fxEnd, x);
	aaPixel = ix;
	if(b > 0.5) {
		b = abs((b - 0.5) * 2);
		
		color.a = b;
		if(aaPixel >= 0 && aaPixel <= 31) {
			continueDrawing = drawGrassBladePixelBlend(aaPixel, iy, color) || continueDrawing;
		}
	}

	return continueDrawing;
}

bool drawGrassBladePixelBlend(int x, int y, vec4 color) {
	//imageStore(result, tileCorner + ivec2(x,y), color);
	//pixels[x][y] = color;
	//return true;

	float srcAlpha = pixels[x][y].a;
	if(srcAlpha < 1.0) {
		pixels[x][y].a = srcAlpha + color.a * (1 - srcAlpha);
		pixels[x][y].rgb = (srcAlpha * pixels[x][y].rgb + (1 - srcAlpha) * color.a * color.rgb) / pixels[x][y].a;
		//pixels[x][y] = color;
		return true;
	} else {
		return false;
	}
		
}

// calculates shading and AO for a grass blade pixel
vec4 getGrassBladeShading(vec3 pos, vec3 n) {
	vec4 color = baseColor;

	// fake AO (darker the closer to the ground, curve starts at the bottom)
	float aoFactor = clamp(distance(pos, vec3(pos.x, 0.0, pos.z)) / aoDist, 0.2, 1.0);

	// Shading
	// double sided
	if(dot(n, camera.CamDir) > 0) {
		n = -n;
	}

	float lambertian = dot(n, -light.Direction);
	float specular = 0.0;
	if(lambertian > 0) {
		vec3 halfDir = normalize(-light.Direction - camera.CamDir);
		float specAngle = max(dot(halfDir, n), 0.0);
		specular = 0.1 * pow(specAngle, shininess);

		color.rgb = 0.8 * lambertian * light.Color.rgb * baseColor.rgb + specular * light.Color.rgb;
	} else {
		// back lighting - less bright
		color.rgb = ((- lambertian) / 3) * mix(baseColor.rgb, light.Color.rgb, - lambertian / 2);
	}
	
	color.rgb *= aoFactor;

	return color;
	
}

// Calculate the translation of the characterizing points of the grass blade
void calcWindTranslation(RandState rng) {
	// vary direction a bit
	vec3 randDirVar = vec3(0.8) + vec3(rand_next(rng), rand_next(rng), rand_next(rng)) * vec3(1 - 0.8);
	vec3 dir = randDirVar * grassConsts.WindDirection;

	// vary windspeed a bit over time
	float speedVar = grassConsts.WindSpeed + (0.01 * rand_next(rng) * (sin(grassConsts.TimeStamp)));

	// translation
	float translationFactor = 0.01 * speedVar * (sin(grassConsts.TimeStamp * grassConsts.WindSpeed + rand_next(rng)) + 0.8);

	tip = tip + translationFactor * dir;
	cP = cP + 0.8 * (translationFactor * dir);
}

// Functions to calculate the frustum
vec3[6] calcFrustumNormals() {
	vec3[6] normals;
	normals[0] = normalize(cross((frustumPoints[3] - frustumPoints[7]), (frustumPoints[6] - frustumPoints[7])));
	normals[1] = normalize(cross((frustumPoints[4] - frustumPoints[5]), (frustumPoints[1] - frustumPoints[5])));
	normals[2] = normalize(cross((frustumPoints[1] - frustumPoints[5]), (frustumPoints[7] - frustumPoints[5])));
	normals[3] = normalize(cross((frustumPoints[6] - frustumPoints[4]), (frustumPoints[0] - frustumPoints[4])));
	normals[4] = normalize(cross((frustumPoints[1] - frustumPoints[3]), (frustumPoints[2] - frustumPoints[3])));
	normals[5] = normalize(cross((frustumPoints[6] - frustumPoints[7]), (frustumPoints[5] - frustumPoints[7])));

	return normals;
}

Ray[12] calcFrustumRays() {
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

vec3[4] calcFrustumIntersections(vec3 planePoint, vec3 planeNormal, out int intersectionCount) {
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

// Utility functions to handle points on the regular grid used to go through the grass blades

float distPointRay(vec3 origin, vec3 direction, vec3 point) {
	return length(cross(normalize(direction), point - origin));
}


float distPointLineSegment(vec3 v, vec3 w, vec3 p) {
  const float l2 = pow(distance(v, w), 2);
  if (l2 == 0.0) return distance(p, v); 
  const float t = max(0, min(1, dot(p - v, w - v) / l2));
  const vec3 projection = v + t * (w - v); 
  return distance(p, projection);
}

// Get step size for regular grid by distance to the camera
float getStepSize(float dist, out float blend) {
	vec3 referencePos = camera.CamPos + camera.CamDir * dist;
	float fStepSize = distance(referencePos, nDCToWorldPos(worldPosToNDC(referencePos) + vec3(1.0, 0.0, 0.0) * initialStepSizeNDC));
	float factor = fStepSize / initialStepSize;
	
	if(factor <= 1) {
		blend = 0;
		factor = 1;
	} else {
		blend = modf(log2(factor), factor);
		factor = pow(2, factor);
	}
	return factor * initialStepSize;
}

bool lessThanByDir(vec3 fst, vec3 snd, vec3 dir) {
	bool result = false;
	if(!any(isinf(fst)) && !any(isinf(snd))) {
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
			result.x = floor(point.x / stepSize + Epsilon) * stepSize;
		} else if (dir.x == 0) {
			result.x = round(point.x / stepSize) * stepSize;
		} else {
			result.x = ceil(point.x / stepSize - Epsilon) * stepSize;
		}
		if(dir.z > 0) {
			result.z = floor(point.z / stepSize + Epsilon) * stepSize;
		} else if(dir.z == 0) {
			result.z = round(point.z / stepSize) * stepSize;
		} else {
			result.z = ceil(point.z / stepSize - Epsilon) * stepSize;
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
				if(j == i || !whichPlanes[j] || j == dontCheck) {
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


// outside frustum tests
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


// methods to draw single pixels
void drawTilePos(vec2 pos, vec4 color) {
	if(pos.x >= 0 && pos.x < 32 && pos.y >= 0 && pos.y < 32) {
		pixels[int(pos.x)][int(pos.y)] = color;
	}
}

void drawWorldPos(vec3 pos, vec4 color) {
	vec3 ndc = worldPosToNDC(pos);
	drawNDC(ndc, color);
}

void drawNDC(vec3 ndc, vec4 color) {
	ivec2 imagePos = ivec2(ndcToImagePos(ndc));
	if (-1.0 <= ndc.x && ndc.x <= 1.0 &&
		-1.0 <= ndc.y && ndc.y <= 1.0 &&
		ndc.z >= 0.0 && ndc.z <= 1.0) {

		imageStore(result, imagePos, color);
		// imageStore(result, imagePos + ivec2(0, 1), color);
		// imageStore(result, imagePos + ivec2(1, 0), color);
		// imageStore(result, imagePos + ivec2(1, 1), color);
	}
}


// Simple non-antialiased line drawing
void drawWorldLine(vec3 start, vec3 end, vec4 color) {
	drawImageLine(worldPosToImagePos(start), worldPosToImagePos(end), color);
}

void drawImageLine(vec2 imageStart, vec2 imageEnd, vec4 color) {
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
		imageStore(result, ivec2(pos), color);
	}
}


// coordinate convertion methods
vec2 worldPosToTilePos(vec3 pos) {
	// converts postion to the screen tile of the workgroup
	return worldPosToImagePos(pos) - tileCorner;
}

vec2 worldPosToImagePos(vec3 pos) {
	vec3 ndc = worldPosToNDC(pos);
	/*
	if(ndc.z < 0 || ndc.z > 1) {
		return vec2(-1);
	} else {
		return ndcToImagePos(ndc);
	}*/
	return ndcToImagePos(ndc);
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

vec3 nDCToWorldPos(vec3 ndc) {
	vec4 p1 = camera.ViewProjInv * vec4(ndc, 1.0);
	return p1.xyz / p1.w;
}

vec4 normalPointPlaneToNormalDistPlane(vec3 normal, vec3 point) {
	return vec4(normal, dot(normal, -point));
}

// Utility
bool epsilonEq(float arg1, float arg2) {
	// return arg1 == arg2;
	return (arg1 < arg2 + Epsilon) && arg1 > arg2 - Epsilon;
}

// http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/
mat4 rotationMatrix(vec3 axis, float angle) {
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}
#endif