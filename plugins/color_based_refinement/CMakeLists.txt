#Add the directory containing all header files from pragmatic (note that pragmatic is header only!)
#set(PRAGMATIC_INCLUDES "${CMAKE_SOURCE_DIR}/plugins/color_based_refinement/external")
#include_directories(${PRAGMATIC_INCLUDES})
#message (STATUS "Pragmatic includes: ${PRAGMATIC_INCLUDES}")
include_directories("external")

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

find_library(MTMETIS_LIBRARY 
             NAMES mtmetis
             HINTS /usr/local/lib/
             PATH_SUFFIXES lib)
message(STATUS "MTMETIS library: ${MTMETIS_LIBRARY}")

find_library(WILDRIVER_LIBRARY 
             NAMES wildriver
             HINTS /usr/local/lib/
             PATH_SUFFIXES lib)
message(STATUS "WILDRIVER library: ${WILDRIVER_LIBRARY}")

#Set flag if OpenMP has been found on the system or output an error if not, since it is required for the use of pragmatic
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DHAVE_OPENMP)
    message(STATUS "Found OPENMP")
    message(STATUS "${CMAKE_C_FLAGS}")
    message(STATUS "${CMAKE_CXX_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found")
endif()

add_definitions( -DTRILIBRARY -DANSI_DECLARATORS -DEXTERNAL_TEST -DNO_TIMER )

#Add plugin to ViennaMesh

VIENNAMESH_ADD_PLUGIN(viennamesh-module-color-refinement 	plugin.cpp
                      color_refinement.cpp
                      pragmatic_mesh.cpp
                      external/triangle.c
                      triangle_interface.c)

if (METIS_FOUND)
  target_link_libraries(viennamesh-module-color-refinement ${METIS_LIBRARIES})
endif()

target_link_libraries(viennamesh-module-color-refinement ${MTMETIS_LIBRARY} ${WILDRIVER_LIBRARY})


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
target_link_libraries(viennamesh-module-color-refinement ${VTK_LIBRARIES})
