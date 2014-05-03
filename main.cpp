#include "pch.h"

#include "mathx"

#include "ogl"
#include "input"

#include "stdx" 
#include "appx" 

#include <iostream>
#include <string>
#include <vector>

#include <fstream>

#include "img"
#include "obj.h"

#include "text"
#include "ui"

namespace glsl
{
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

struct Camera
{
	glm::vec3 pos;
	glm::mat3 orientation;
	float fov;
	float aspect;
	float nearPlane;
	float farPlane;

	Camera()
		: fov(90.0f)
		, aspect(1.0f)
		, nearPlane(0.1f)
		, farPlane(1000.0f) { }

	void lookTo(glm::vec3 const& pos, glm::vec3 const& where, glm::vec3 const& up)
	{
		this->pos = pos;
		orientation[2] = normalize(pos - where);
		reorientate(up);
	}

	void reorientate(glm::vec3 const& up)
	{
		orientation[0] = normalize( cross(up, orientation[2]) );
		orientation[1] = normalize( cross(orientation[2], orientation[0]) );
		orientation[2] = normalize( orientation[2] );
	}

	glm::mat4 view() const { return translate(glm::mat4(transpose(orientation)), -pos); }
	glm::mat4 proj() const { return glm::perspective(fov, aspect, nearPlane, farPlane); }
	glm::mat4 viewProj() const { return proj() * view(); }
};

struct RenderableMesh
{
	ogl::Buffer positions;
	ogl::Buffer normals;
	ogl::Buffer texCoords;
	ogl::VertexArrays vertices;
	ogl::Buffer indices;

	template <class V, class N, class T, class I>
	RenderableMesh(V const& vposRange, N const& vnrmRange, T const& vtexRange, I const& idcsRange)
		: positions(ogl::Buffer::init(GL_ARRAY_BUFFER, vposRange))
		, normals(ogl::Buffer::init(GL_ARRAY_BUFFER, vnrmRange))
		, texCoords(vtexRange.empty() ? nullptr : ogl::Buffer::init(GL_ARRAY_BUFFER, vtexRange))
		, indices(ogl::Buffer::init(GL_ELEMENT_ARRAY_BUFFER, idcsRange))
	{
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

	void bind()
	{
		vertices.bind();
		indices.bind(GL_ELEMENT_ARRAY_BUFFER);
	}

	void draw(GLenum mode, unsigned indexOffset, unsigned numIndices)
	{
		glDrawElements(mode, numIndices, GL_UNSIGNED_INT, (void*) indexOffset);
	}
};

namespace main_ex
{
}

int main()
{
	try
	{
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

		keyboard.keyEvent[GLFW_KEY_ESCAPE].pressOnce = [&]() { glfwSetWindowShouldClose(wnd, true); };
		
		keyboard.keyEvent[GLFW_KEY_ENTER].pressOnce = [&]() { mouse.setMouseCapture(wnd, !mouse.mouseCaptured); };
		mouse.buttonEvent[GLFW_MOUSE_BUTTON_RIGHT].all = [&](bool pressed) { mouse.setMouseCapture(wnd, pressed); };

		// load shaders
		auto preamble = "";
		ogl::ProgramWithTime backgroundShader("data/background.glsl");
		ogl::ProgramWithTime tonemapShader("data/tonemap.glsl");
		ogl::ProgramWithTime phongShader("data/phong.glsl");
		ogl::ProgramWithTime textShader("data/text.glsl", "", ogl::ProgramWithTime::HasGS);
		ogl::ProgramWithTime uiShader("data/ui.glsl", "", ogl::ProgramWithTime::HasGS);

		auto maybeReloadKernels = [&]()
		{
			backgroundShader.maybeReload();
			tonemapShader.maybeReload();
			phongShader.maybeReload();
			textShader.maybeReload();
			uiShader.maybeReload();
		};
		keyboard.keyEvent[GLFW_KEY_R].pressOnce = [&]() { maybeReloadKernels(); };
		
		if (false)
		{
			auto shaderBin = tonemapShader.getBinary();
			std::ofstream shaderBinFile("shaderBin.txt", std::ios::binary);
			shaderBinFile.write(shaderBin.data(), shaderBin.size());
		}

		// load object
		Obj simpleObj = parse_object(stdx::load_file("data/simple.obj").c_str());
		RenderableMesh simpleMesh(simpleObj.v, simpleObj.n, simpleObj.t, simpleObj.f);
		
		// load environment map
		ogl::Texture envMap = nullptr;
		{
//			auto image = img::load_image<float, 3>(stdx::load_binary_file("beach_probe.hdr"));
//			envMap = ogl::Texture::create2D(GL_TEXTURE_2D, GL_RGBA16F, image.dim.x, image.dim.y, 0, image.pixels.data(), GL_FLOAT, GL_RGB);
		}

		// camera
		Camera camera;
		camera.lookTo(glm::vec3(-0.6f, 0.14f, 0.3f) * 10.0f, glm::vec3(), glm::vec3(0.0f, 1.0f, 0.0f));
		glm::vec3 lightDirection = normalize(glm::vec3(1.0f, -4.0f, -3.0f));
		
		auto camConstBuffer = ogl::Buffer::create(GL_UNIFORM_BUFFER, sizeof(glsl::CameraConstants));

		// quad processing necessities
		ogl::NullVertexArray nullVertexArrays;

		// Text
		text::FreeType freeTypeLib;
		text::Face font(freeTypeLib, "C:/Windows/Fonts/tahoma.ttf", text::PtSize(10)); // "C:/Windows/Fonts/consola.ttf", "Inconsolata-Regular.ttf", "C:/Windows/Fonts/Andale.ttf"
		auto textUiPtr = ui::UserInterface::create(&freeTypeLib, &font, &textShader, &uiShader);
		auto& textUi = *textUiPtr;
		
		// framebuffer setup
		glm::uvec2 screenDim;
		ogl::Texture hdrTexture = nullptr;
		ogl::RenderBuffer depthBuffer = nullptr;
		ogl::Framebuffer hdrBuffer = nullptr;

		wnd.resize = [&](unsigned width, unsigned height)
		{
			if (width == 0 || height == 0)
				return;

			screenDim.x = width;
			screenDim.y = height;
			glViewport(0, 0, width, height);
			camera.aspect = (float) width / (float) height;
			
			auto newHdrTexture = ogl::Texture::create2D(GL_TEXTURE_2D, GL_RGBA16F, width, height);
			auto newDepthBuffer = ogl::RenderBuffer::create(GL_DEPTH_COMPONENT, width, height);
			auto newHdrBuffer = ogl::Framebuffer::create(&newHdrTexture.get(), 1, newDepthBuffer);

			auto paddedWidth = glm::ceil_mul(width, 32U);
			auto paddedHeight = glm::ceil_mul(height, 32U);

			hdrBuffer = std::move(newHdrBuffer);
			hdrTexture = std::move(newHdrTexture);
			depthBuffer = std::move(newDepthBuffer);
		};
		wnd.initialResize();

		std::string testString = "test str";

		// main loop
		bool paused = false;
		keyboard.keyEvent[GLFW_KEY_P].pressOnce = [&]() { paused = !paused; };

		bool fast = false;
		keyboard.keyEvent[GLFW_KEY_F].pressOnce = [&]() { fast = !fast; };

		bool enableUi = true;
		keyboard.keyEvent[GLFW_KEY_U].pressOnce = [&]() { enableUi = !enableUi; };

		double lastTime = glfwGetTime();
		float animTime = 0.0f;
		unsigned frameIdx = 0;
		float smoothDt = 1.0f;
		float smoothFDt = 1.0f;
		float camSpeed = 1.0f;

		while (!wnd.shouldClose())
		{
			glfwPollEvents();
			
			// frame times
			float dt;
			bool nextSecond;
			bool halfSecond;
			{
				double thisTime = glfwGetTime();
				dt = (float) (thisTime - lastTime);
				nextSecond = (long long) thisTime > (long long) lastTime;
				halfSecond = thisTime - floor(thisTime) >= 0.5f;
				lastTime = thisTime;
			}
			smoothDt = 0.1f * dt + 0.9f * smoothDt;
			unsigned seed = static_cast<unsigned>(lastTime * 1000.0);
			
			// reload
			if (nextSecond)
			{
				maybeReloadKernels();
			}

			// handle camera input
			{
				if (mouse.mouseCaptured)
				{
					auto rotationDelta = -90.0f * mouse.relativeMouseDelta(wnd);
				
					camera.orientation = (glm::mat3) glm::rotate(glm::mat4(camera.orientation), rotationDelta.x, glm::vec3(0.0f, 1.0f, 0.0f));
					camera.orientation = (glm::mat3) glm::rotate(glm::mat4(camera.orientation), rotationDelta.y, glm::vec3(1.0f, 0.0f, 0.0f));
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

				camera.pos += moveDelta * (camera.orientation * glm::vec3(moveInput));
			}

			// animate sun light
			{
				if (!paused)
					animTime += (fast) ? 5.0f * dt : dt;

				float rotSpeed = 0.1f;
				float riseSpeed = 0.05f;

				float height = 0.6f + 0.4f * cos(riseSpeed * animTime);
				float sinH = sqrt(1.0f - 0.99f * height * height);

//				lightDirection = glm::vec3(-0.285239756f, -0.530926347f, -0.797969520);
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
				camConst.ViewProj = camera.viewProj();
				camConst.ViewProjInv = inverse(camConst.ViewProj);
				camConst.CamPos = camera.pos;
				camConst.CamDir = -camera.orientation[2];
				camConst.NearPlane = camera.nearPlane;
				camConst.FarPlane = camera.farPlane;
				camConstBuffer.write(GL_UNIFORM_BUFFER, stdx::make_range_n(&camConst, 1));
			}

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

			// Test scene
			{
				hdrBuffer.bind(GL_FRAMEBUFFER);

				glEnable(GL_DEPTH_TEST);
				camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);

				simpleMesh.bind();
				phongShader.bind();
				simpleMesh.draw(GL_TRIANGLES, 0, (unsigned) simpleObj.f.size() * 3);
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
			if (enableUi)
			{
				glDisable(GL_DEPTH_TEST);
				glEnable(GL_BLEND);
				
				camConstBuffer.bind(GL_UNIFORM_BUFFER, 0);

				// text test
				{
					float fontBrightness = 0.0f;
					glBlendColor(fontBrightness, fontBrightness, fontBrightness, 1.0f);
					glBlendEquation(GL_FUNC_ADD);
					glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_SRC_COLOR);

/*					textRenderer.drawText(freeTypeLib, font, glm::ivec2(16, 16),
						"Glyph images are always loaded, transformed, and described in the cartesian coordinate \n"
						"system in FreeType (which means that increasing Y corresponds to upper scanlines), unlike \n"
						"the system typically used for bitmaps (where the topmost scanline has coordinate 0). We \n"
						"must thus convert between the two systems when we define the pen position, and when we \n"
						"compute the topleft position of the bitmap."
						);
					textRenderer.flushText();
*/				}
				
				// ui test
				{
					textUi.state.cursorVisible = halfSecond;
					textUi.setup.rect.min = glm::ivec2(800, 300);
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
					
					if (auto uiGroup = ui::Group(ui, nullptr))
					{
						ui.addText(nullptr, "UI", "", nullptr);

						ui.addText(nullptr, "test str", testString.c_str(), testString);

						if (auto uiGroup = ui::Group(ui, nullptr))
						{
							ui.addText(nullptr, "Tweak", "", nullptr);
							ui.addSlider(&lightDirection, "light dir", 0.1f, 6.0f, nullptr);
							ui.addSlider(&camSpeed, "cam speed", camSpeed, 10.0f, camSpeed, 2.0f);

//							if (auto uiUnion = ui::Union(ui))
							{
								ui.addButton(3, "test button", nullptr);
								ui.addInteractiveButton(4, "test button", true, nullptr);
							}
						}
					}
					
					textUi.flushWidgets();

					float fontBrightness = 0.01f;
					glBlendColor(fontBrightness, fontBrightness, fontBrightness, 1.0f);
					glBlendEquation(GL_FUNC_ADD);
//					glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_SRC_COLOR);
					glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ONE_MINUS_SRC_COLOR);

					textUi.flushText();
				}
			
				glDisable(GL_BLEND);
			}

			{
				wnd.swapBuffers();

				// accept input
				mouse.accept();
				keyboard.accept();
			}
			
			++frameIdx;

			if (true)
			{
				glFinish();
				double frameTime = glfwGetTime();
				float fdt = (float) (frameTime - lastTime);
				smoothFDt = 0.1f * fdt + 0.9f * smoothFDt;
			}
			else
				smoothFDt = smoothDt;

			if (frameIdx % 30 == 0)
			{
				std::cout << "Frame time: " << 1000.0f * smoothFDt << " ms; " << 1.0f / smoothFDt << " FPS" << std::endl;
			}

			Sleep(12);
		}
	}
	catch (std::exception const &excpt)
	{
		std::cout << "Fatal error: " << excpt.what() << std::endl;
		return -1;
	}

	return 0;
}
