project(${LIGHTER_PROJECT})
cmake_minimum_required(VERSION 2.8)

#--------------------------------------------------------------------
# Create target
#--------------------------------------------------------------------
set(APPLICATION_SRC
  main.cpp
  )

add_executable(${LIGHTER_PROJECT} ${APPLICATION_SRC})
add_dependencies(${LIGHTER_PROJECT} ${LIGHTER_DEPENDENCIES})
  
#--------------------------------------------------------------------
# Install files other than the application
#--------------------------------------------------------------------
if (LIGHTER_INSTALL)
    install (TARGETS ${LIGHTER_PROJECT} DESTINATION bin)
endif()
