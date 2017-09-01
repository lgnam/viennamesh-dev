#Include Pragmatic and Triangle header files
include_directories(external)

#Add definitions for triangle library
add_definitions( -DTRILIBRARY -DANSI_DECLARATORS -DEXTERNAL_TEST -DNO_TIMER )

#Set flag if OpenMP has been found on the system or output an error if not, since it is required for the use of pragmatic
#find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DHAVE_OPENMP)
    message(STATUS "Found OPENMP")
    message(STATUS "${CMAKE_C_FLAGS}")
    message(STATUS "${CMAKE_CXX_FLAGS}")
#else()
 #   message(FATAL_ERROR "OpenMP not found")
endif()

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

#Add ViennaMesh Plugin
VIENNAMESH_ADD_PLUGIN(viennamesh-module-color-refinement plugin.cpp
                      color_refinement.cpp
                      pragmatic_mesh.cpp)

#Link with VTK library
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
target_link_libraries(viennamesh-module-color-refinement ${VTK_LIBRARIES})

#Link with METIS library
target_link_libraries(viennamesh-module-color-refinement ${METIS_LIBRARIES})
