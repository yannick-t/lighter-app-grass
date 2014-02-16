#include "pch.h"

#include "mathx"

#include "ogl"

#include "stdx" 
#include "appx" 

#include <iostream>
#include <string>
#include <vector>

#include <fstream>

#include "img"

#include "textx"
#include "uix"

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

struct UiTextRenderer : ui::TextRenderer
{
	text::FreeType& lib;
	text::Face &font;
	text::TextRenderer& renderer;

	UiTextRenderer(text::FreeType& lib, text::Face& font, text::TextRenderer& renderer)
		: lib(lib), font(font), renderer(renderer) { }

	virtual glm::aabb<glm::ivec2> selectChar(glm::ivec2 cursorPos, size_t& charIdx, glm::ivec2 insertPos, char const* text, size_t maxChars, CharType type)
	{
		switch (type)
		{
		case UTF8: default: return renderer.selectChar(lib, font, cursorPos, charIdx, insertPos, text, maxChars);
		case UTF16: return renderer.selectChar(lib, font, cursorPos, charIdx, insertPos, reinterpret_cast<FT_Int16 const*>(text), maxChars);
		case UTF32: return renderer.selectChar(lib, font, cursorPos, charIdx, insertPos, reinterpret_cast<FT_Int32 const*>(text), maxChars);
		}
	}
	virtual glm::aabb<glm::ivec2> boundText(char const* text, size_t maxChars, glm::ivec2 insertPos, CharType type) override
	{
		switch (type)
		{
		case UTF8: default: return renderer.boundText(lib, font, insertPos, text, maxChars);
		case UTF16: return renderer.boundText(lib, font, insertPos, reinterpret_cast<FT_Int16 const*>(text), maxChars);
		case UTF32: return renderer.boundText(lib, font, insertPos, reinterpret_cast<FT_Int32 const*>(text), maxChars);
		}
	}
	virtual glm::aabb<glm::ivec2> drawText(glm::ivec2 insertPos, char const* text, size_t maxChars, CharType type) override
	{
		switch (type)
		{
		case UTF8: default: return renderer.drawText(lib, font, insertPos, text, maxChars);
		case UTF16: return renderer.drawText(lib, font, insertPos, reinterpret_cast<FT_Int16 const*>(text), maxChars);
		case UTF32: return renderer.drawText(lib, font, insertPos, reinterpret_cast<FT_Int32 const*>(text), maxChars);
		}
	}
	using TextRenderer::drawText;
};

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
	ogl::VertexArrays vertices;
	ogl::Buffer indices;

	template <class V, class N, class I>
	RenderableMesh(V const& vposRange, N const& vnrmRange, I const& idcsRange)
		: positions(ogl::Buffer::init(GL_ARRAY_BUFFER, vposRange))
		, normals(ogl::Buffer::init(GL_ARRAY_BUFFER, vnrmRange))
		, indices(ogl::Buffer::init(GL_ELEMENT_ARRAY_BUFFER, idcsRange))
	{
		vertices.bind();
		positions.bind(GL_ARRAY_BUFFER);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, sizeof(*vposRange.data()) / sizeof(float), GL_FLOAT, GL_FALSE, 0, nullptr);
		normals.bind(GL_ARRAY_BUFFER);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, sizeof(*vnrmRange.data()) / sizeof(float), GL_FLOAT, GL_FALSE, 0, nullptr);
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
		ogl::GLFW glfw(3, 3);
		ogl::GLFWWindow wnd(1024, 576, "Molecules + GI", nullptr);

		wnd.makeCurrent();
		ogl::GLEW glew;
		glfw.enableVSync(false);
		glfw.enableDebugMessages();
		glEnable(GL_FRAMEBUFFER_SRGB);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

		// load shaders
		auto preamble = "";
		ogl::ProgramWithTime backgroundShader("data/background.glsl");
		ogl::ProgramWithTime tonemapShader("data/tonemap.glsl");
		ogl::ProgramWithTime textShader("data/text.glsl", "", ogl::ProgramWithTime::HasGS);
		ogl::ProgramWithTime uiShader("data/ui.glsl", "", ogl::ProgramWithTime::HasGS);

		auto maybeReloadKernels = [&]()
		{
			backgroundShader.maybeReload();
			tonemapShader.maybeReload();
			textShader.maybeReload();
			uiShader.maybeReload();
		};
		
		if (false)
		{
			auto shaderBin = tonemapShader.getBinary();
			std::ofstream shaderBinFile("shaderBin.txt", std::ios::binary);
			shaderBinFile.write(shaderBin.data(), shaderBin.size());
		}

		// window & rendering set up
		Camera camera;
		camera.lookTo(glm::vec3(-0.6f, 0.14f, 0.3f) * 10.0f, glm::vec3(), glm::vec3(0.0f, 1.0f, 0.0f));
		glm::vec3 lightDirection = normalize(glm::vec3(1.0f, -4.0f, -3.0f));

		glm::uvec2 screenDim;
		ogl::Texture hdrTexture = nullptr;
		ogl::RenderBuffer depthBuffer = nullptr;
		ogl::Framebuffer hdrBuffer = nullptr;

		// Environment map
		ogl::Texture envMap = nullptr;
		{
//			auto image = img::load_image<float, 3>(stdx::load_binary_file("beach_probe.hdr"));
//			envMap = ogl::Texture::create2D(GL_TEXTURE_2D, GL_RGBA16F, image.dim.x, image.dim.y, 0, image.pixels.data(), GL_FLOAT, GL_RGB);
		}

		// Text
		text::FreeType freeTypeLib;
		text::Face font(freeTypeLib, "C:/Windows/Fonts/tahoma.ttf", text::PtSize(10)); // "C:/Windows/Fonts/consola.ttf", "Inconsolata-Regular.ttf", "C:/Windows/Fonts/Andale.ttf"
		ui::UiRenderer widgetRenderer(&uiShader, 1024);
		text::TextRenderer textRenderer(freeTypeLib, &textShader, 10000);
		auto uiTextRenderer = UiTextRenderer(freeTypeLib, font, textRenderer);
		ui::TextUi textUi(&widgetRenderer, &uiTextRenderer);

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

		// input handling
		bool paused = false;
		bool fast = false;
		bool mouseCaptured = false;

		bool enableUi = true;
		
		glm::ivec2 lastMousePos;
		glm::fvec2 mouseDelta;
		wnd.mouse = [&](int x, int y)
		{
			glm::ivec2 mousePos(x, y);
			mouseDelta = glm::vec2(mousePos - lastMousePos) / glm::vec2(screenDim);
			lastMousePos = mousePos;
		};

		auto setMouseCapture = [&](bool capture)
		{
			if (capture)
				glfwSetInputMode(wnd, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
			else
				glfwSetInputMode(wnd, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			
			double x, y;
			glfwGetCursorPos(wnd, &x, &y);
			lastMousePos = glm::ivec2((int) x, (int) y);
			
			mouseCaptured = capture;
		};
		
		bool buttonState[16] = { };
		bool buttonChanged[16] = { };
		wnd.buttons = [&](int key, bool pressed)
		{
			buttonState[key] = pressed;
			buttonChanged[key] = true;

			if (key == GLFW_MOUSE_BUTTON_RIGHT)
				setMouseCapture(pressed);
		};
		
		std::vector<unsigned> inputQueue;
		wnd.text = [&](unsigned chr)
		{
			inputQueue.push_back(chr);
		};

		bool keyState[GLFW_KEY_LAST] = { };
		bool keyChanged[GLFW_KEY_LAST] = { };
		wnd.keyboard = [&](int key, bool pressed, bool repeat)
		{
			if (pressed) switch (key)
			{
				unsigned key;
				case GLFW_KEY_BACKSPACE: key = ui::InputKeys::DeleteLast; goto INSERT_INPUT_KEY;
				case GLFW_KEY_DELETE: key = ui::InputKeys::DeleteNext; goto INSERT_INPUT_KEY;
				case GLFW_KEY_LEFT: key = ui::InputKeys::MoveLeft; goto INSERT_INPUT_KEY;
				case GLFW_KEY_RIGHT: key = ui::InputKeys::MoveRight; goto INSERT_INPUT_KEY;
				case GLFW_KEY_UP: key = ui::InputKeys::MoveUp; goto INSERT_INPUT_KEY;
				case GLFW_KEY_DOWN: key = ui::InputKeys::MoveDown; goto INSERT_INPUT_KEY;
				INSERT_INPUT_KEY: inputQueue.push_back(key);
			}

			if (!repeat)
			{
				keyState[key] = pressed;
				keyChanged[key] = true;

				if (key == GLFW_KEY_ESCAPE && pressed)
					glfwSetWindowShouldClose(wnd, true);

				if (key == GLFW_KEY_R && pressed)
					maybeReloadKernels();

				if (key == GLFW_KEY_P && pressed)
					paused = !paused;
				if (key == GLFW_KEY_F && pressed)
					fast = !fast;

				if (key == GLFW_KEY_C && pressed)
					enableUi = !enableUi;

				if (key == GLFW_KEY_ENTER && pressed)
					setMouseCapture(!mouseCaptured);
			}
		};

		auto camConstBuffer = ogl::Buffer::create(GL_UNIFORM_BUFFER, sizeof(glsl::CameraConstants));

		auto nullVertexArrays = ogl::VertexArrays();
		nullVertexArrays.bind();
		auto dummyVertexBuffer = ogl::Buffer::create(GL_ARRAY_BUFFER, sizeof(float) * 3);
		glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, nullptr);
		glEnableVertexAttribArray(0);
		nullVertexArrays.unbind();

		std::string testString = "test str";

		// main loop
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
				if (mouseCaptured)
				{
					auto rotationDelta = -90.0f * mouseDelta;
				
					camera.orientation = (glm::mat3) glm::rotate(glm::mat4(camera.orientation), rotationDelta.x, glm::vec3(0.0f, 1.0f, 0.0f));
					camera.orientation = (glm::mat3) glm::rotate(glm::mat4(camera.orientation), rotationDelta.y, glm::vec3(1.0f, 0.0f, 0.0f));
					camera.reorientate(glm::vec3(0.0f, 1.0f, 0.0f));
				}

				auto moveDelta = dt * 3.0f * camSpeed
					* (1.0f + 3.0f * keyState[GLFW_KEY_LEFT_SHIFT])
					* glm::mix(1.0f, 0.3f, keyState[GLFW_KEY_LEFT_CONTROL]);

				glm::ivec3 moveInput(
					  keyState[GLFW_KEY_D] - keyState[GLFW_KEY_A]
					, keyState[GLFW_KEY_SPACE] - keyState[GLFW_KEY_Q]
					, keyState[GLFW_KEY_S] - keyState[GLFW_KEY_W]
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

					textRenderer.drawText(freeTypeLib, font, glm::ivec2(16, 16),
						"Glyph images are always loaded, transformed, and described in the cartesian coordinate \n"
						"system in FreeType (which means that increasing Y corresponds to upper scanlines), unlike \n"
						"the system typically used for bitmaps (where the topmost scanline has coordinate 0). We \n"
						"must thus convert between the two systems when we define the pen position, and when we \n"
						"compute the topleft position of the bitmap."
						);
					textRenderer.flushText();
				}
				
				// ui test
				{
					textUi.state.cursorVisible = halfSecond;
					textUi.setup.rect.min = glm::ivec2(800, 300);
					textUi.setup.rect.max = textUi.setup.rect.min + glm::ivec2(200, 500);

					textUi.state.mouse.pos = lastMousePos;
					textUi.state.mouse.primary = buttonState[GLFW_MOUSE_BUTTON_LEFT];
					textUi.state.mouse.primaryChanged = buttonChanged[GLFW_MOUSE_BUTTON_LEFT];
					textUi.state.mouse.secondary = buttonState[GLFW_MOUSE_BUTTON_RIGHT];
					textUi.state.mouse.secondaryChanged = buttonChanged[GLFW_MOUSE_BUTTON_RIGHT];

					auto& ui = textUi.reset(textUi.state.mouse, inputQueue);
					inputQueue.clear();

					glBlendEquation(GL_FUNC_ADD);
					glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
					
					ui::Group uiGroup(ui, nullptr);
					ui.addText(nullptr, "UI", "", nullptr);

					ui.addText(nullptr, "test str", testString.c_str(), testString);

					{
						ui::Group uiGroup(ui, nullptr);
						ui.addText(nullptr, "Tweak", "", nullptr);
						ui.addSlider(&lightDirection, "light dir", 0.1f, 6.0f, nullptr);
						ui.addSlider(&camSpeed, "cam speed", camSpeed, 10.0f, camSpeed, 2.0f);
					}
					
					widgetRenderer.flushWidgets();

					float fontBrightness = 0.01f;
					glBlendColor(fontBrightness, fontBrightness, fontBrightness, 1.0f);
					glBlendEquation(GL_FUNC_ADD);
//					glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_SRC_COLOR);
					glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ONE_MINUS_SRC_COLOR);

					textRenderer.flushText();
				}
			
				glDisable(GL_BLEND);
			}

			{
				wnd.swapBuffers();

				// reset input
				mouseDelta = glm::vec2();
				memset(keyChanged, 0, sizeof(keyChanged));
				memset(buttonChanged, 0, sizeof(buttonChanged));
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
