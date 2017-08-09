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

VIENNAMESH_ADD_PLUGIN(viennamesh-module-viennagrid plugin.cpp
                      affine_transform.cpp
                      extract_boundary.cpp
                      extract_plc_geometry.cpp
                      hyperplane_clip.cpp
                      laplace_smooth.cpp
                      map_regions.cpp
                      merge_meshes.cpp
                      hull_set_regions.cpp
                      split_mesh.cpp
                      make_boundary_box.cpp
                      simplexify.cpp
                      line_coarsening.cpp
                      interpolate_quantities.cpp
                      refine_plc_lines.cpp
                      center_mesh.cpp
                      mark_hull_regions.cpp
                      scale_quantities.cpp
                      stretch_middle.cpp
                      douglas_peucker_line_smoothing.cpp
                      uniform_refine.cpp
                      change_cell_region.cpp
                      chessboard_coloring.cpp
                      parallel_coloring.cpp)