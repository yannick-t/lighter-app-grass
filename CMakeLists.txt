project(${LIGHTER_PROJECT})

#--------------------------------------------------------------------
# Create target
#--------------------------------------------------------------------
set(APPLICATION_SYSTEM_SRC
  renderer.glsl.h
  )
set(APPLICATION_SHADER_SRC
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
  ${APPLICATION_SYSTEM_SRC}
  )

add_executable(${LIGHTER_PROJECT} ${APPLICATION_SRC})
add_lighter_dependencies(${LIGHTER_PROJECT})
copy_data_files(${LIGHTER_PROJECT} "${${LIGHTER_PROJECT}_BINARY_DIR}" "${${LIGHTER_PROJECT}_SOURCE_DIR}" "${APPLICATION_SHADER_SRC}")
copy_lighter_data_files(${LIGHTER_PROJECT} "${${LIGHTER_PROJECT}_BINARY_DIR}")
target_precompiled_header(${LIGHTER_PROJECT} pch.h pch.cpp)
  
#--------------------------------------------------------------------
# Install files other than the application
#--------------------------------------------------------------------
if (LIGHTER_INSTALL)
    install (TARGETS ${LIGHTER_PROJECT} DESTINATION bin)
    install (FILES ${APPLICATION_SYSTEM_SRC} DESTINATION bin)
    install (FILES ${APPLICATION_SHADER_SRC} DESTINATION ${LIGHTER_INSTALL_DATA_DIRECTORY})
endif()
