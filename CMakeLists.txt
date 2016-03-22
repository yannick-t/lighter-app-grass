cmake_minimum_required(VERSION 3.0.0)
set(LIGHTER_PROJECT riseandshine)
project(${LIGHTER_PROJECT})

add_subdirectory("add/make")

#--------------------------------------------------------------------
# Data
#--------------------------------------------------------------------
set(APPLICATION_DATA
  data/simple.obj
  data/simpler.obj
  )
  
#--------------------------------------------------------------------
# Create target
#--------------------------------------------------------------------
set(APPLICATION_SYSTEM_SRC
  renderer.glsl.h
  )
set(APPLICATION_SHADER_SRC
  background.glsl
  tonemap.glsl
  phong.glsl
  )
set(APPLICATION_CODE
  main.cpp
  obj.cpp
  obj.h
  pch.cpp
  pch.h
  )
set(APPLICATION_SRC
  ${APPLICATION_CODE}
  ${APPLICATION_SHADER_SRC}
  ${APPLICATION_SYSTEM_SRC}
  )
source_group(
  "Source Files" FILES ${APPLICATION_CODE}
  )
source_group(
  "Source Files\\Shaders" FILES ${APPLICATION_SHADER_SRC} ${APPLICATION_SYSTEM_SRC}
  )

set(APPLICATION_DATA_STAMP ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/application_data.stamp)
add_custom_command(OUTPUT ${APPLICATION_DATA_STAMP}
    COMMAND ${CMAKE_COMMAND} -E touch \"${APPLICATION_DATA_STAMP}\"
    DEPENDS ${APPLICATION_SHADER_SRC} ${APPLICATION_DATA} ${LIGHTER_DATA_STAMP})

add_executable(${LIGHTER_PROJECT} ${APPLICATION_SRC} ${APPLICATION_DATA_STAMP})
#add_lighter_dependencies(${LIGHTER_PROJECT})
target_link_libraries(${LIGHTER_PROJECT} libglew_static glfw lighter)
target_precompiled_header(${LIGHTER_PROJECT} pch.h pch.cpp)

configure_file(add/make/cmake/application.vcxproj.user-template ${LIGHTER_PROJECT}.vcxproj.user @ONLY)

#--------------------------------------------------------------------
# Install files other than the application
#--------------------------------------------------------------------
if (LIGHTER_INSTALL)
    install (TARGETS ${LIGHTER_PROJECT} DESTINATION bin)
    install (FILES ${APPLICATION_SYSTEM_SRC} DESTINATION bin)
    install (FILES ${APPLICATION_SHADER_SRC} ${APPLICATION_DATA} DESTINATION ${LIGHTER_INSTALL_DATA_DIRECTORY})
	
    add_custom_command(TARGET ${LIGHTER_PROJECT} POST_BUILD COMMAND ${LIGHTER_INSTALL_COMMAND})
endif()
