#Add the directory containing all header files from pragmatic (note that pragmatic is header only!)
set(PRAGMATIC_INCLUDES "${CMAKE_SOURCE_DIR}/plugins/color_based_refinement/external")
include_directories(${PRAGMATIC_INCLUDES})
message (STATUS "Pragmatic includes: ${PRAGMATIC_INCLUDES}")

#find METIS package
find_package(METIS QUIET)
if (METIS_FOUND)
  message(STATUS "Found METIS")
  message(STATUS "METIS include directories: ${METIS_INCLUDE_DIRS}")
  message(STATUS "METIS libraries: ${METIS_LIBRARIES}")
  include_directories(${METIS_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "METIS not found")
endif()

#Add plugin to ViennaMesh
VIENNAMESH_ADD_PLUGIN(viennamesh-module-color-refinement 	plugin.cpp
                      color_refinement.cpp
                      pragmatic_mesh.cpp)

if (METIS_FOUND)
  target_link_libraries(viennamesh-module-color-refinement ${METIS_LIBRARIES})
endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
target_link_libraries(viennamesh-module-color-refinement ${VTK_LIBRARIES})
