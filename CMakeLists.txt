cmake_minimum_required(VERSION 2.8)

PROJECT(ReadPLY)
SET(VTK_DIR /Users/alonzoaltamirano/Repos/BUILD-VTK)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(ReadPLY MACOSX_BUNDLE ReadPLY.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(ReadPLY ${VTK_LIBRARIES})
else()
  target_link_libraries(ReadPLY vtkHybrid vtkWidgets)
endif()
