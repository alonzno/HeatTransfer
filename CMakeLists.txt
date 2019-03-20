cmake_minimum_required(VERSION 2.8)

PROJECT(ReadPLY)
SET(VTK_DIR ~/Repos/BUILD-VTK-OLD)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(ReadPLY ReadPLY.cxx sim.hxx)
target_compile_options(ReadPLY PRIVATE -std=c++11)

if(VTK_LIBRARIES)
  target_link_libraries(ReadPLY ${VTK_LIBRARIES})
else()
  target_link_libraries(ReadPLY vtkHybrid vtkWidgets)
endif()
