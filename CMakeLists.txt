project(${LIGHTER_PROJECT})

#--------------------------------------------------------------------
# Create target
#--------------------------------------------------------------------
set(APPLICATION_SHADER_SRC
  renderer.glsl.h
  background.glsl
  tonemap.glsl
  )
set(APPLICATION_SRC
  main.cpp
  obj.cpp
  obj.h
  pch.cpp
  pch.h
  ${APPLICATION_SHADER_SRC}
  )

add_executable(${LIGHTER_PROJECT} ${APPLICATION_SRC})
add_lighter_dependencies(${LIGHTER_PROJECT})
copy_data_files(${LIGHTER_PROJECT} "${${LIGHTER_PROJECT}_BINARY_DIR}" "${${LIGHTER_PROJECT}_SOURCE_DIR}" "${APPLICATION_SHADER_SRC}")
copy_lighter_data_files(${LIGHTER_PROJECT} "${${LIGHTER_PROJECT}_BINARY_DIR}")
  
#--------------------------------------------------------------------
# Install files other than the application
#--------------------------------------------------------------------
if (LIGHTER_INSTALL)
    install (TARGETS ${LIGHTER_PROJECT} DESTINATION bin)
endif()
