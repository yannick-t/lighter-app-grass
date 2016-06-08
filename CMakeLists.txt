cmake_minimum_required(VERSION 3.0.0)
project(riseandshine)

include(add/make/everything.cmake)

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
  random.glsl.h
  )
set(APPLICATION_SHADER_SRC
  background.glsl
  tonemap.glsl
  phong.glsl
  referenceGeometryGrass.glsl
  computeShaderGrass.glsl
  csResultShader.glsl
  debug.glsl
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

add_executable(${PROJECT_NAME} ${APPLICATION_SRC} ${APPLICATION_DATA_STAMP})
#add_lighter_dependencies(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} lighter ${LIGHTER_MAKE_DEFAULT_LIBS})
target_precompiled_header(${PROJECT_NAME} pch.h pch.cpp)

configure_file(add/make/cmake/application.vcxproj.user-template ${PROJECT_NAME}.vcxproj.user @ONLY)

#--------------------------------------------------------------------
# Install files other than the application
#--------------------------------------------------------------------
if (LIGHTER_INSTALL)
    install (TARGETS ${PROJECT_NAME} DESTINATION bin)
    install (FILES ${APPLICATION_SYSTEM_SRC} DESTINATION bin)
    install (FILES ${APPLICATION_SHADER_SRC} ${APPLICATION_DATA} DESTINATION ${LIGHTER_INSTALL_DATA_DIRECTORY})

    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${LIGHTER_INSTALL_COMMAND})
endif()
