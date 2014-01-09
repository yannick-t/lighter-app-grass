project(${LIGHTER_PROJECT})

#--------------------------------------------------------------------
# Create target
#--------------------------------------------------------------------
set(APPLICATION_SRC
  main.cpp
  )

add_executable(${LIGHTER_PROJECT} ${APPLICATION_SRC})
add_lighter_dependencies(${LIGHTER_PROJECT})
copy_lighter_data_files(${LIGHTER_PROJECT})
  
#--------------------------------------------------------------------
# Install files other than the application
#--------------------------------------------------------------------
if (LIGHTER_INSTALL)
    install (TARGETS ${LIGHTER_PROJECT} DESTINATION bin)
endif()
