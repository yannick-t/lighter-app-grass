#ifndef CLSL_FAKE_VEC3
	#define CLSL_RENDERER_VEC3_PADDING(x) x
#else
	#define CLSL_RENDERER_VEC3_PADDING(x)
#endif

struct CameraConstants
{
	vec2 Resolution;
	vec2 PixelWidth;

	uvec2 IntResolution;
	uint FrameIdx;
	uint NumBlocks32X;

	mat4 ViewProj;
	mat4 ViewProjInv;
	
	vec3 CamPos;
	CLSL_RENDERER_VEC3_PADDING(float _pad1;)
	vec3 CamDir;
	CLSL_RENDERER_VEC3_PADDING(float _pad2;)

	float NearPlane;
	float FarPlane;
};

struct LightConstants
{
	vec3 Direction;
	vec4 Color;
};


struct CSGrassConstants 
{
	vec3 FtBDirection;
	int TileDivisor;

	vec3 PerpFtBDir;
	float MinHeight;

	float MaxHeight;
	float RelAODist;
	float MinWidth;
	float MaxWidth;

	float MinDist;
	float MaxDist;
	float StepPxAtMinDist;
	int DrawDebugInfo; // levels of debug info

	vec3 WindDirection;
	float WindSpeed;

	float TimeStamp;
};

struct Ray
{
	vec3 Start;
	vec3 Dir;
};
