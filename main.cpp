#include "pch.h"

#include "stdx" 
#include "appx" 
#include "debug"

#include "mathx"

#include "ogl"
#include "input"

#include "pool"

#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <fstream>

#include "img"
#include "scenex"

#include "text"
#include "ui"

namespace glsl {
	using namespace glm;
#include "renderer.glsl.h"
}

// HELP:
// Movement: W,A,S,D,Space,Q + Mouse
// |- Faster: Left Shift
// |- Slower: Left Ctrl
// Light:
// |- Toogle slow / fast mode: F
// |- Pause / reanimate: P
// Escape: ESC

bool const stdx::is_debugger_present = IsDebuggerPresent() != FALSE;

struct Camera {
	const int maxFootprint = 4;

	glm::vec3 pos;
	glm::mat3 orientation;
	float fov;
	float aspect;
	float nearPlane;
	float farPlane;

	glm::mat4 viewProjection;
	glm::mat4 viewProjectionInverse;

	// frustum defined by points and normals
	glm::vec3 nearPoints[4];
	glm::vec3 farPoints[4];
	glm::vec3 rightNormal;
	glm::vec3 leftNormal;
	glm::vec3 topNormal;
	glm::vec3 bottomNormal;
	glm::vec3 nearNormal;
	glm::vec3 farNormal;

	glm::vec3 max;
	glm::vec3 min;

	bool frustumChanged = true;

	Camera()
		: fov(90.0f)
		  , aspect(1.0f)
		  , nearPlane(0.1f)
		  , farPlane(100.0f) {
	}

	void lookTo(glm::vec3 const& pos, glm::vec3 const& where, glm::vec3 const& up) {
		this->pos = pos;
		orientation[2] = normalize(pos - where);
		reorientate(up);
	}

	void reorientate(glm::vec3 const& up) {
		orientation[0] = normalize(cross(up, orientation[2]));
		orientation[1] = normalize(cross(orientation[2], orientation[0]));
		orientation[2] = normalize(orientation[2]);
		recalcualteFrustum();
	}

	void rePosition(glm::vec3 pos) {
		this->pos = pos;
		recalcualteFrustum();
	}

	void recalcualteFrustum() {
		frustumChanged = true;
		calcViewProj();

		glm::vec4 ntl = (viewProjectionInverse * glm::vec4(-1.0, 1.0, -1.0, 1.0)); // top left
		nearPoints[0] = ntl.xyz / ntl.w;
		glm::vec4 ntr = (viewProjectionInverse * glm::vec4(1.0, 1.0, -1.0, 1.0)); // top right
		nearPoints[1] = ntr.xyz / ntr.w;
		glm::vec4 nbl = (viewProjectionInverse * glm::vec4(-1.0, -1.0, -1.0, 1.0)); // bottom left
		nearPoints[2] = nbl.xyz / nbl.w;
		glm::vec4 nbr = (viewProjectionInverse * glm::vec4(1.0, -1.0, -1.0, 1.0)); // bottom right
		nearPoints[3] = nbr.xyz / nbr.w;


		glm::vec4 ftl = (viewProjectionInverse * glm::vec4(-1.0, 1.0, 1.0, 1.0)); // top left
		farPoints[0] = ftl.xyz / ftl.w;
		glm::vec4 ftr = (viewProjectionInverse * glm::vec4(1.0, 1.0, 1.0, 1.0)); // top right
		farPoints[1] = ftr.xyz / ftr.w;
		glm::vec4 fbl = (viewProjectionInverse * glm::vec4(-1.0, -1.0, 1.0, 1.0)); // bottom left
		farPoints[2] = fbl.xyz / fbl.w;
		glm::vec4 fbr = (viewProjectionInverse * glm::vec4(1.0, -1.0, 1.0, 1.0)); // bottom right
		farPoints[3] = fbr.xyz / fbr.w;

		// Todo: optimize
		rightNormal = glm::normalize(cross(nearPoints[1] - farPoints[1], farPoints[3] - farPoints[1]));
		leftNormal = glm::normalize(cross(farPoints[0] - nearPoints[0], nearPoints[2] - nearPoints[0]));
		topNormal = glm::normalize(cross(farPoints[1] - nearPoints[1], nearPoints[0] - nearPoints[1]));
		bottomNormal = glm::normalize(cross(farPoints[2] - nearPoints[2], nearPoints[3] - nearPoints[2]));
		nearNormal = glm::normalize(cross(nearPoints[0] - nearPoints[1], nearPoints[3] - nearPoints[1]));
		farNormal = glm::normalize(cross(farPoints[1] - farPoints[0], farPoints[2] - farPoints[0]));

		min = glm::vec3(std::numeric_limits<float>::max());
		max = glm::vec3(-std::numeric_limits<float>::max());

		for (int i = 0; i < 4; i++) {
			// printf("i: %d - x: %f, y: %f, z: %f \n", i, nearPoints[i].x, nearPoints[i].y, nearPoints[i].z);
			for (int j = 0; j < 3; j++) {

				if (nearPoints[i][j] < min[j]) {
					min[j] = nearPoints[i][j];
				}
				if (nearPoints[i][j] > max[j]) {
					max[j] = nearPoints[i][j];
				}
			}
		}
		for (int i = 0; i < 4; i++) {
			// printf("i: %d - x: %f, y: %f, z: %f \n", i, farPoints[i].x, farPoints[i].y, farPoints[i].z);
			for (int j = 0; j < 3; j++) {
				if (farPoints[i][j] < min[j]) {
					min[j] = farPoints[i][j];
				}
				if (farPoints[i][j] > max[j]) {
					max[j] = farPoints[i][j];
				}
			}
		}
	}

	float computeFootprintNDCQuad(std::vector<glm::vec3> &points) {
		float f = 0;   
		int j = points.size() - 1; 

		for (int i = 0; i < points.size(); i++)
		{
			f = f + (points[j].x + points[i].x) * (points[j].y - points[i].y);
			j = i;
		}
		f = f / 2;
		for (int i = 0; i < 4; i++) {
			printf("i: %d - x: %f, y: %f, z: %f \n", i, points[i].x, points[i].y, points[i].z);
		}

		printf("footprint: %f \n", abs(f));
		return abs(f);
	}

	bool quadInFrustum(std::vector<glm::vec3> points) {
		// Find out if points of the quad are all outside of one of the frustums planes

		if (pointsOutsideOfPlane(pos, rightNormal, points) || pointsOutsideOfPlane(pos, leftNormal, points) || pointsOutsideOfPlane(pos, topNormal, points) ||
			pointsOutsideOfPlane(pos, bottomNormal, points) || pointsOutsideOfPlane(nearPoints[0], nearNormal, points) || pointsOutsideOfPlane(farPoints[0], farNormal, points)) {

			return false;
		}
		return true;
	}

	bool pointsOutsideOfPlane(glm::vec3 planePoint, glm::vec3 planeNormal, std::vector<glm::vec3> points) {
		// outside == on the side of the normal
		for (int i = 0; i < 4; i++) {
			if (dot(planeNormal, points[i] - planePoint) < 0) {
				return false;
			}
		}
		return true;
	}

	void pointsToNDC(std::vector<glm::vec3> &p, std::vector<glm::vec3> &result) {
		for (int i = 0; i < p.size(); i++) {
			glm::vec4 pc1 = viewProjection * glm::vec4(p[i], 1.0);
			result[i] = (i, pc1.xyz * 1 / pc1.w);
		}
	}

	bool normalizedPointInFrustum(glm::vec3 pointN) {
		return (pointN.x >= -1 && pointN.x <= 1) &&
			(pointN.y >= -1 && pointN.y <= 1) &&
			(pointN.z >= -1 && pointN.z <= 1);
	}

	float distanceRectToCamera(glm::vec3 min, glm::vec3 max) {
		// Rectangle needs to be axis aligned
		float dx = std::max(std::max(min.x - pos.x, 0.0f), pos.x - max.x);
		float dy = std::max(std::max(min.y - pos.y, 0.0f), pos.y - max.y);
		float dz = std::max(std::max(min.z - pos.z, 0.0f), pos.z - max.z);
		return std::sqrtf(dx*dx + dy*dy + dz*dz);
	}

	glm::mat4 view() const {
		return transform(-pos * orientation, glm::vec3(1.0f), transpose(orientation));
	}

	glm::mat4 proj() const {
		return glm::perspective(glm::radians(fov), aspect, nearPlane, farPlane);
	}

	void calcViewProj() {
		viewProjection = proj() * view();
		viewProjectionInverse = glm::inverse(viewProjection);
	}
};

struct GrassPatch {
	glm::vec3 pos;
	float size;
	float density;
};

struct RenderableMesh {
	MOVE_GENERATE(RenderableMesh, MOVE_5, MEMBER, vertices, MEMBER, positions, MEMBER, normals, MEMBER, texCoords, MEMBER, indices)

	ogl::VertexArrays vertices;
	ogl::Buffer positions;
	ogl::Buffer normals;
	ogl::Buffer texCoords;
	ogl::Buffer indices;

	RenderableMesh(nullptr_t = nullptr) {
	}

	template <class V, class N, class T, class I>
	RenderableMesh(V const& vposRange, N const& vnrmRange, T const& vtexRange, I const& idcsRange)
		: vertices(ogl::VertexArrays::create())
		  , positions(ogl::Buffer::init(GL_ARRAY_BUFFER, vposRange))
		  , normals(ogl::Buffer::init(GL_ARRAY_BUFFER, vnrmRange))
		  , texCoords(vtexRange.empty() ? nullptr : ogl::Buffer::init(GL_ARRAY_BUFFER, vtexRange))
		  , indices(ogl::Buffer::init(GL_ELEMENT_ARRAY_BUFFER, idcsRange)) {
		vertices.bind();
		positions.bind(GL_ARRAY_BUFFER);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, sizeof(*vposRange.data()) / sizeof(float), GL_FLOAT, GL_FALSE, 0, nullptr);
		normals.bind(GL_ARRAY_BUFFER);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, sizeof(*vnrmRange.data()) / sizeof(float), GL_FLOAT, GL_FALSE, 0, nullptr);
		if (texCoords) {
			texCoords.bind(GL_ARRAY_BUFFER);
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, sizeof(*vtexRange.data()) / sizeof(float), GL_FLOAT, GL_FALSE, 0, nullptr);
		}
	}

	void bind() {
		vertices.bind();
		indices.bind(GL_ELEMENT_ARRAY_BUFFER);
	}

	void draw(GLenum mode, unsigned indexOffset, unsigned numIndices) {
		glDrawElements(mode, numIndices, GL_UNSIGNED_INT, (void*) indexOffset);
	}
};

namespace main_ex {
}

float baseGrassPatchSize = 200;
float grassPatchMaxHeight = 1.5;
float baseGrassDensity = 0.1; // number of grass blades per unit for largest grass patches
float maxPatchDistanceToSubdivide = 10;
float maxPatchSubdivRecursion = 3;
float subdivDensityFactor = 2;
float maxGrassBladeWidthFactor = 10;

void addAndSubdivide(float x, float z, float size, float density, Camera camera, std::vector<GrassPatch> &patches, int recursion = 0);
void calcPatches(Camera camera, std::vector<GrassPatch> &patches);

void calcPatches(Camera camera, std::vector<GrassPatch> &patches) {
	patches.clear();

	if (camera.min.y < grassPatchMaxHeight && camera.max.y > grassPatchMaxHeight) {
		// Test for patches that are possibly in the frustum if they are
		// Round down to next patch position
		glm::vec2 start = glm::vec2((floorf(camera.min.x / baseGrassPatchSize)) * baseGrassPatchSize,
			(floorf(camera.min.z / baseGrassPatchSize)) * baseGrassPatchSize);
		// Round up to next patch position
		glm::vec2 end = glm::vec2((ceilf(camera.max.x / baseGrassPatchSize)) * baseGrassPatchSize,
			(ceilf(camera.max.z / baseGrassPatchSize)) * baseGrassPatchSize);

		float x = 0;
		float z = 0;
		for (float x = start.x; x <= end.x; x += baseGrassPatchSize) {
			for (float z = start.y; z <= end.y; z += baseGrassPatchSize) {
				addAndSubdivide(x, z, baseGrassPatchSize, baseGrassDensity, camera, patches);
			}
		}
	}
}

void addAndSubdivide(float x, float z, float size, float density, Camera camera, std::vector<GrassPatch> &patches, int recursion) {
	// Check if patch is in view frustum
	std::vector<glm::vec3> points{ glm::vec3(x, grassPatchMaxHeight, z), glm::vec3(x + size, grassPatchMaxHeight, z),
		glm::vec3(x + size, grassPatchMaxHeight, z + size), glm::vec3(x, grassPatchMaxHeight, z + size) };


	if (camera.quadInFrustum(points)) {
		// patch + subpatches
		/*
		if (relativeFootprint < 0) {
			// calc footprint
			
			std::vector<glm::vec3> ndcPoints(4);
			camera.pointsToNDC(points, ndcPoints);

			relativeFootprint = camera.computeFootprintNDCQuad(ndcPoints) / camera.maxFootprint;
		}*/

		GrassPatch patch;
		patch.pos = glm::vec3(x, 0.0, z);
		patch.size = size;
		patch.density = density;
		patches.push_back(patch);


		// Subpatches
		if ((camera.distanceRectToCamera(points[0], points[2]) + (((float)recursion / (float)maxPatchSubdivRecursion) * maxPatchDistanceToSubdivide)) < maxPatchDistanceToSubdivide) { // make sure it doesn't subdivide infinitely
			float subdivSize = size / 2;
			for (int i = 0; i < 4; i++) {
				addAndSubdivide((x + (i % 2) * subdivSize), (z + ((int)(i / 2)) * subdivSize), subdivSize, density * subdivDensityFactor, camera, patches, recursion + 1);
			}
		}
	}
}

int run() {
	ogl::Platform platform(3, 3);
	ogl::Window wnd(1024, 576, "Rise and Shine", nullptr);

	// window & rendering set up
	wnd.makeCurrent();
	ogl::GLEW glew;
	platform.enableVSync(false);
	platform.enableDebugMessages();
	glEnable(GL_FRAMEBUFFER_SRGB);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// input handling
	input::Mouse mouse(wnd);
	input::UiKeyboard keyboard(wnd, ui::InputKeys());

	keyboard.keyEvent[GLFW_KEY_ESCAPE].pressOnce = [&]() {
		glfwSetWindowShouldClose(wnd, true);
	};

	keyboard.keyEvent[GLFW_KEY_ENTER].pressOnce = [&]() {
		mouse.setMouseCapture(wnd, !mouse.mouseCaptured);
	};
	mouse.buttonEvent[GLFW_MOUSE_BUTTON_RIGHT].all = [&](bool pressed) {
		mouse.setMouseCapture(wnd, pressed);
	};

	// Debug options
	// button to stop recalculating the patches
	bool reCalcPatches = true;
	keyboard.keyEvent[GLFW_KEY_C].pressOnce = [&]() {
		reCalcPatches = !reCalcPatches;
	};
	// button to draw the boarders of the patches
	bool drawPatchBorders = false;
	keyboard.keyEvent[GLFW_KEY_B].pressOnce = [&]() {
		drawPatchBorders = !drawPatchBorders;
	};


	// load shaders
	auto preamble = "";
	std::vector<ogl::ProgramWithTime*> shaders;
	ogl::ProgramWithTime backgroundShader("data/background.glsl");
	shaders.push_back(&backgroundShader);
	ogl::ProgramWithTime tonemapShader("data/tonemap.glsl");
	shaders.push_back(&tonemapShader);
	ogl::ProgramWithTime geometryGrassShader("data/referenceGeometryGrass.glsl", "", ogl::ProgramWithTime::HasHS | ogl::ProgramWithTime::HasDS);
	shaders.push_back(&geometryGrassShader);
	ogl::ProgramWithTime computeShaderGrass("data/computeShaderGrass.glsl", "", ogl::ProgramWithTime::HasCS);
	shaders.push_back(&computeShaderGrass);
	ogl::ProgramWithTime csResultShader("data/csResultShader.glsl");
	shaders.push_back(&csResultShader);
	ogl::ProgramWithTime debugShader("data/debug.glsl");
	shaders.push_back(&debugShader);
	ogl::ProgramWithTime textShader("data/text.glsl", "", ogl::ProgramWithTime::HasGS);
	shaders.push_back(&textShader);
	ogl::ProgramWithTime uiShader("data/ui.glsl", "", ogl::ProgramWithTime::HasGS);
	shaders.push_back(&uiShader);

	auto maybeReloadKernels = [&]() {
		bool updated = false;
		for (auto&& s : shaders)
			updated |= s->maybeReload() == 1;
		return updated;
	};
	keyboard.keyEvent[GLFW_KEY_R].pressOnce = [&]() {
		maybeReloadKernels();
	};

	if (false) {
		auto shaderBin = tonemapShader.getBinary();
		std::ofstream shaderBinFile("shaderBin.txt", std::ios::binary);
		shaderBinFile.write(shaderBin.data(), shaderBin.size());
	}

	// camera
	Camera camera;
	camera.lookTo(glm::vec3(-1.0f, 0.0f, 0.0f) * 10.0f, glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	auto camConstBuffer = ogl::Buffer::create(GL_UNIFORM_BUFFER, sizeof(glsl::CameraConstants));

	glm::vec3 lightDirection = normalize(glm::vec3(1.0f, -4.0f, -3.0f));
	glm::vec4 lightColor = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
	auto lightConstBuffer = ogl::Buffer::create(GL_UNIFORM_BUFFER, sizeof(glsl::LightConstants));

	// Grass Patch
	// Todo: more sliders
	auto grassPatchConstBuffer = ogl::Buffer::create(GL_UNIFORM_BUFFER, sizeof(glsl::GrassPatchConstants));
	std::vector<GrassPatch> patches;

	auto csGrassConstBuffer = ogl::Buffer::create(GL_UNIFORM_BUFFER, sizeof(glsl::CSGrassConstants));

	// load environment map
	ogl::Texture envMap = nullptr; {
		//		auto image = img::load_image<float, 3>(stdx::load_binary_file("beach_probe.hdr"));
		//		envMap = ogl::Texture::create2D(GL_TEXTURE_2D, GL_RGBA16F, image.dim.x, image.dim.y, 0, image.pixels.data(), GL_FLOAT, GL_RGB);
	}

	// quad processing necessities
	ogl::NullVertexArray nullVertexArrays;

	// Text
	text::FreeType freeTypeLib;
	auto font = text::Face::create(freeTypeLib, "C:/Windows/Fonts/tahoma.ttf", text::PtSize(10)); // "C:/Windows/Fonts/consola.ttf", "Inconsolata-Regular.ttf", "C:/Windows/Fonts/Andale.ttf"
	auto textUiPtr = ui::UserInterface::create(&freeTypeLib, &font, &textShader, &uiShader);
	auto& textUi = *textUiPtr;

	// framebuffer setup
	glm::uvec2 screenDim;
	appx::resource_pool renderTargetPool;

	wnd.resize = [&](unsigned width, unsigned height) {
		if (width == 0 || height == 0)
			return;

		screenDim.x = width;
		screenDim.y = height;
		glViewport(0, 0, width, height);
		camera.aspect = (float) width / (float) height;

		renderTargetPool.free_unused();
	};
	wnd.initialResize();

	// state
	bool paused = false;
	keyboard.keyEvent[GLFW_KEY_P].pressOnce = [&]() {
		paused = !paused;
	};

	bool fast = false;
	keyboard.keyEvent[GLFW_KEY_F].pressOnce = [&]() {
		fast = !fast;
	};

	bool enableUi = true;
	keyboard.keyEvent[GLFW_KEY_U].pressOnce = [&]() {
		enableUi = !enableUi;
	};

	float camSpeed = 1.0f;

	std::string testString = "test str";

	// test ui
	auto&& tweakUi = [&](ui::UniversalInterface& ui) {
		if (auto uiGroup = ui::Group(ui, nullptr)) {
			ui.addText(nullptr, "Tweak", "", nullptr);
			ui.addSlider(&camSpeed, "cam speed", camSpeed, 10.0f, camSpeed, 2.0f);
			ui.addSlider(&baseGrassDensity, "base grass density for patches", baseGrassDensity, 100, baseGrassDensity, 1);
			ui.addSlider(&maxPatchDistanceToSubdivide, "the maximum distance of a patch to subdivide it", maxPatchDistanceToSubdivide, 100, maxPatchDistanceToSubdivide, 0.5f);
			ui.addSlider(&maxPatchSubdivRecursion, "maximum recursion depth for the subdivision", maxPatchSubdivRecursion, 10, maxPatchSubdivRecursion, 1.0f);
			ui.addSlider(&subdivDensityFactor, "factor applied to the density for subdivided patches", subdivDensityFactor, 20, subdivDensityFactor, 1.0f);
			ui.addSlider(&maxGrassBladeWidthFactor, "maximum scaling factor applied t the grass blade width with distance", maxGrassBladeWidthFactor, 1000, maxGrassBladeWidthFactor);
			//			if (auto uiUnion = ui::Union(ui)) 
			{
				// ui.addButton(3, "test button", nullptr);
				// ui.addInteractiveButton(4, "test button", true, nullptr);
			}
		}
	};

	// Load default preset
	auto defaultIniFile = "default.ini";
	try {
		ui::load_ini_file(defaultIniFile, tweakUi);
	}
	catch (...) {
		appx::print_exception();
	}

	// Always store preset on exit
	struct PresetGuard {
		stdx::fun_ref<void (ui::UniversalInterface&)> ui;

		~PresetGuard() {
			try {
				ui::save_ini_file("lastrun.ini", ui);
			}
			catch (...) {
				appx::print_exception();
			}
		}
	} presetGuard = {tweakUi};

	// main loop
	double lastTime = glfwGetTime();
	float animTime = 0.0f;
	unsigned frameIdx = 0;
	float smoothDt = 1.0f;
	float smoothFDt = 1.0f;

	float fps = 0;

	int bladeCount = 0;
	int patchCount = 0;

	while (!wnd.shouldClose()) {
		glfwPollEvents();

		// frame times
		float dt;
		bool nextSecond;
		bool halfSecond; {
			double thisTime = glfwGetTime();
			dt = (float) (thisTime - lastTime);
			nextSecond = (long long) thisTime > (long long) lastTime;
			halfSecond = thisTime - floor(thisTime) >= 0.5f;
			lastTime = thisTime;
		}
		smoothDt = 0.1f * dt + 0.9f * smoothDt;
		unsigned seed = static_cast<unsigned>(lastTime * 1000.0);

		// reload
		if (nextSecond) {
			maybeReloadKernels();
		}

		// handle camera input
		{
			if (mouse.mouseCaptured) {
				auto rotationDelta = -0.5f * glm::pi<float>() * mouse.relativeMouseDelta(wnd);

				camera.orientation = camera.orientation * (glm::mat3) rotate(rotationDelta.x, glm::vec3(0.0f, 1.0f, 0.0f));
				camera.orientation = camera.orientation * (glm::mat3) rotate(rotationDelta.y, glm::vec3(1.0f, 0.0f, 0.0f));
				camera.reorientate(glm::vec3(0.0f, 1.0f, 0.0f));
			}

			auto moveDelta = dt * 3.0f * camSpeed
				* (1.0f + 3.0f * keyboard.keyState[GLFW_KEY_LEFT_SHIFT])
				* glm::mix(1.0f, 0.3f, keyboard.keyState[GLFW_KEY_LEFT_CONTROL]);

			glm::ivec3 moveInput(
				keyboard.keyState[GLFW_KEY_D] - keyboard.keyState[GLFW_KEY_A]
				, keyboard.keyState[GLFW_KEY_SPACE] - keyboard.keyState[GLFW_KEY_Q]
				, keyboard.keyState[GLFW_KEY_S] - keyboard.keyState[GLFW_KEY_W]
				);

			if (moveInput.x != 0 || moveInput.y != 0 || moveInput.z != 0) {
				camera.rePosition(camera.pos + moveDelta * (camera.orientation * glm::vec3(moveInput))); 
			}

		}

		// animate sun light
		{
			if (!paused)
				animTime += (fast) ? 5.0f * dt : dt;

			float rotSpeed = 0.1f;
			float riseSpeed = 0.05f;

			float height = 0.6f + 0.4f * cos(riseSpeed * animTime);
			float sinH = sqrt(1.0f - 0.99f * height * height);

			//lightDirection = glm::vec3(-0.285239756f, -0.530926347f, -0.797969520);
			lightDirection = glm::vec3(sin(rotSpeed * animTime) * sinH, -height, cos(rotSpeed * animTime) * sinH);
			lightDirection = normalize(lightDirection);
		}

		glsl::CameraConstants camConst;

		// update
		{
			camConst.Resolution = glm::vec2(screenDim);
			camConst.PixelWidth = 1.0f / camConst.Resolution;
			camConst.IntResolution = screenDim;
			camConst.FrameIdx = frameIdx;
			camConst.NumBlocks32X = glm::ceil_div(screenDim.x, 32U);
			camConst.ViewProj = camera.viewProjection;
			camConst.ViewProjInv = inverse(camConst.ViewProj);
			camConst.CamPos = camera.pos;
			camConst.CamDir = -camera.orientation[2];
			camConst.NearPlane = camera.nearPlane;
			camConst.FarPlane = camera.farPlane;
			camConstBuffer.write(GL_UNIFORM_BUFFER, stdx::make_range_n(&camConst, 1));
		}

		glsl::LightConstants lightConst; {
			lightConst.Direction = lightDirection;
			lightConst.Color = lightColor;
			lightConstBuffer.write(GL_UNIFORM_BUFFER, stdx::make_range_n(&lightConst, 1));
		}

		auto hdrTexture = renderTargetPool.acquire(ogl::TextureDesc::make2D(GL_TEXTURE_2D, GL_RGBA16F, screenDim.x, screenDim.y));
		auto hdrDepthBuffer = renderTargetPool.acquire(ogl::RenderBufferDesc::make(GL_DEPTH_COMPONENT, screenDim.x, screenDim.y));
		auto hdrBuffer = renderTargetPool.acquire(ogl::MakeFramebufferDesc::textures().color(hdrTexture).depthBuffer(hdrDepthBuffer));

		// Background
		{
			hdrBuffer.bind(GL_FRAMEBUFFER);
			glClearColor(1.0f, 0.0f, 1.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

			glDisable(GL_DEPTH_TEST);
			camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);
			nullVertexArrays.bind();
			backgroundShader.bind();
			glDrawArrays(GL_TRIANGLES, 0, 3);
		}

		// Reference Geometry Grass
		{

			if (false) {
				bladeCount = 0;

				// Calculate patches if view frustum has changed or O is pressed
				if (camera.frustumChanged && reCalcPatches || keyboard.keyState[GLFW_KEY_O]) {
					calcPatches(camera, patches);
					camera.frustumChanged = false;
				}

				// Draw patches
				for (int i = 0; i < patches.size(); i++) {

					int blades = patches[i].size * patches[i].size * patches[i].density;

					glsl::GrassPatchConstants grassPatchConst;
					{
						grassPatchConst.Position = patches[i].pos;
						grassPatchConst.Size = patches[i].size;
						grassPatchConst.MaxHeight = grassPatchMaxHeight;
						grassPatchConst.MaxWidthFactor = maxGrassBladeWidthFactor;
						grassPatchConst.PatchId = i;
						grassPatchConstBuffer.write(GL_UNIFORM_BUFFER, stdx::make_range_n(&grassPatchConst, 1));
					}

					hdrBuffer.bind(GL_FRAMEBUFFER);

					glEnable(GL_DEPTH_TEST);
					camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);
					lightConstBuffer.bind(GL_UNIFORM_BUFFER, 1);
					grassPatchConstBuffer.bind(GL_UNIFORM_BUFFER, 2);
					nullVertexArrays.bind();
					geometryGrassShader.bind();
					glPatchParameteri(GL_PATCH_VERTICES, 1); // 1 Vertex per (Tessellation) Patch 
					glDrawArrays(GL_PATCHES, 0, blades); // Draw vertices, shader will create grass blades for each one with random positions

					bladeCount += blades;

					// Draw boarders of the patches
					if (drawPatchBorders) {
						float points[] = {
							patches[i].pos.x, grassPatchMaxHeight, patches[i].pos.z,
							patches[i].pos.x, grassPatchMaxHeight, patches[i].pos.z + patches[i].size,
							patches[i].pos.x, grassPatchMaxHeight, patches[i].pos.z + patches[i].size,
							patches[i].pos.x + patches[i].size, grassPatchMaxHeight, patches[i].pos.z + patches[i].size,
							patches[i].pos.x + patches[i].size, grassPatchMaxHeight, patches[i].pos.z + patches[i].size,
							patches[i].pos.x + patches[i].size, grassPatchMaxHeight, patches[i].pos.z,
							patches[i].pos.x + patches[i].size, grassPatchMaxHeight, patches[i].pos.z,
							patches[i].pos.x, grassPatchMaxHeight, patches[i].pos.z
						};

						GLuint vbo = 0;
						glGenBuffers(1, &vbo);
						glBindBuffer(GL_ARRAY_BUFFER, vbo);
						glBufferData(GL_ARRAY_BUFFER, 24 * sizeof(float), points, GL_STATIC_DRAW);

						GLuint vao = 0;
						glGenVertexArrays(1, &vao);
						glBindVertexArray(vao);
						glEnableVertexAttribArray(0);
						glBindBuffer(GL_ARRAY_BUFFER, vbo);
						glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

						debugShader.bind();
						glDrawArrays(GL_LINES, 0, 8);
					}
				}

				patchCount = patches.size();
			}
		}


		// Compute Shader Grass
		auto csResult = renderTargetPool.acquire(ogl::TextureDesc::make2D(GL_TEXTURE_2D, GL_RGBA32F, screenDim.x, screenDim.y));
		{
			glsl::CSGrassConstants csGrassConstants;
			{
				// Todo: sorted array of grid cells for front to back rendering
				csGrassConstants.ScreenDim = screenDim;
				csGrassConstants.GridStart = glm::vec3(0.0);
				csGrassConstants.GridEnd = glm::vec3(40.0, 0.0, 40.0);
				csGrassConstants.Step = 2.0;
				csGrassConstBuffer.write(GL_UNIFORM_BUFFER, stdx::make_range_n(&csGrassConstants, 1));
			}

			camConstBuffer.bind(GL_UNIFORM_BUFFER, 1);
			lightConstBuffer.bind(GL_UNIFORM_BUFFER, 2);
			// csResult.bind(GL_TEXTURE_2D, 2);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, csResult);
			glBindImageTexture(0, csResult, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_R32F);
			csGrassConstBuffer.bind(GL_UNIFORM_BUFFER, 3);
			computeShaderGrass.bind();
			glDispatchCompute(512, 512, 1);

			glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

			// Draw result
			hdrBuffer.bind(GL_FRAMEBUFFER);
			glClearColor(1.0f, 0.0f, 1.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

			glDisable(GL_DEPTH_TEST);
			camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);
			csResult.bind(GL_TEXTURE_2D, 1);
			nullVertexArrays.bind();
			csResultShader.bind();
			glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		}

		// Blit / tonemap
		{
			ogl::Framebuffer::unbind(GL_FRAMEBUFFER);
			glClearColor(1.0f, 1.0f, 0.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

			glDisable(GL_DEPTH_TEST);
			camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);
			hdrTexture.bind(GL_TEXTURE_2D, 0);
			nullVertexArrays.bind();
			tonemapShader.bind();
			glDrawArrays(GL_TRIANGLES, 0, 3);
		}

		// Text
		if (enableUi) {
			glDisable(GL_DEPTH_TEST);
			glEnable(GL_BLEND);

			camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);

			// text test
			{
				float fontBrightness = 0.0f;
				glBlendColor(fontBrightness, fontBrightness, fontBrightness, 1.0f);
				glBlendEquation(GL_FUNC_ADD);
				glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_SRC_COLOR);

				/*				textRenderer.drawText(freeTypeLib, font, glm::ivec2(16, 16),
									"Glyph images are always loaded, transformed, and described in the cartesian coordinate \n"
									"system in FreeType (which means that increasing Y corresponds to upper scanlines), unlike \n"
									"the system typically used for bitmaps (where the topmost scanline has coordinate 0). We \n"
									"must thus convert between the two systems when we define the pen position, and when we \n"
									"compute the topleft position of the bitmap."
									);
								textRenderer.flushText();
				*/
			}

			// ui test
			{
				textUi.state.cursorVisible = halfSecond;
				textUi.setup.rect.min = glm::ivec2(screenDim.x - 220, 300);
				textUi.setup.rect.max = textUi.setup.rect.min + glm::ivec2(200, 500);

				textUi.state.mouse.pos = mouse.lastMousePos;
				textUi.state.mouse.primary = mouse.buttonState[GLFW_MOUSE_BUTTON_LEFT];
				textUi.state.mouse.primaryChanged = mouse.buttonChanged[GLFW_MOUSE_BUTTON_LEFT];
				textUi.state.mouse.secondary = mouse.buttonState[GLFW_MOUSE_BUTTON_RIGHT];
				textUi.state.mouse.secondaryChanged = mouse.buttonChanged[GLFW_MOUSE_BUTTON_RIGHT];

				auto& ui = textUi.reset(textUi.state.mouse, keyboard.inputQueue);
				keyboard.inputQueue.clear();

				glBlendEquation(GL_FUNC_ADD);
				glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

				ui.addSlider(&dt, "dt (ms)", dt * 1000.0f, 500.0f, nullptr);

				char fpsString[20];
				_snprintf(fpsString, 20, "%f FPS", fps);
				ui.addText(nullptr, fpsString, "", nullptr);

				char bladeCountString[30];
				_snprintf(bladeCountString, 30, "%d grass blades", bladeCount);
				ui.addText(nullptr, bladeCountString, "", nullptr);

				char patchCountString[20];
				_snprintf(patchCountString, 20, "%d patches", patchCount);
				ui.addText(nullptr, patchCountString, "", nullptr);

				tweakUi(ui);
				// ui::preset_user_interface(ui, tweakUi, defaultIniFile);

				textUi.flushWidgets();

				float fontBrightness = 0.01f;
				glBlendColor(fontBrightness, fontBrightness, fontBrightness, 1.0f);
				glBlendEquation(GL_FUNC_ADD);
				//				glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_SRC_COLOR);
				glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ONE_MINUS_SRC_COLOR);

				textUi.flushText();
			}

			glDisable(GL_BLEND);
		} {
			wnd.swapBuffers();

			renderTargetPool.free_unused_and_next_frame(5);

			// accept input
			mouse.accept();
			keyboard.accept();
			keyboard.inputQueue.clear();
		}

		++frameIdx;

		if (true) {
			glFinish();
			double frameTime = glfwGetTime();
			float fdt = (float) (frameTime - lastTime);
			smoothFDt = 0.1f * fdt + 0.9f * smoothFDt;
		}
		else
			smoothFDt = smoothDt;

		if (frameIdx % 15 == 0) {
			fps = 1.0f / smoothFDt;
			// std::cout << "Frame time: " << 1000.0f * smoothFDt << " ms; " << 1.0f / smoothFDt << " FPS" << std::endl;
		}

		Sleep(12);
	}

	return 0;
}

int main() {
	try {
		return stdx::dump_on_exception(run);
	}
	catch (std::exception const& excpt) {
		std::cout << "Fatal error: " << excpt.what() << std::endl;
		return -1;
	}
}
