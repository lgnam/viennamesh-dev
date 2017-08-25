/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                ViennaMesh - The Vienna Meshing Framework
                            -----------------

                    http://viennamesh.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include "parallel_coloring.hpp"
#include <set>
#include <numeric>
#include <omp.h>

namespace viennamesh
{
    //Typedefs
    typedef viennagrid::mesh                                                                  MeshType;
    typedef viennagrid::result_of::element<MeshType>::type                                    VertexType;
    typedef viennagrid::result_of::element<MeshType>::type                                    EdgeType;
    typedef viennagrid::result_of::element<MeshType>::type                                    TriangleType;

    //typedef viennagrid::result_of::element_range<MeshType, 0>::type                           VertexRange;
    //typedef viennagrid::result_of::element_range<MeshType, viennagrid::triangle_tag>::type      ElementRange;

    typedef viennagrid::result_of::element_range<MeshType, 2>::type                           TriangleRange; 
    typedef viennagrid::result_of::iterator<TriangleRange>::type                              TriangleIterator; 

    typedef viennagrid::result_of::neighbor_range<MeshType, 1, 2>::type                       TriangleNeighborRangeType;
    typedef viennagrid::result_of::iterator<TriangleNeighborRangeType>::type                  TriangleNeighborIterator;

    typedef viennagrid::result_of::neighbor_range<MeshType, 2, 0>::type                       NeighborTriangleRange;
    typedef viennagrid::result_of::iterator<NeighborTriangleRange>::type                      NeighborTriangleIterator;

    parallel_coloring::parallel_coloring() {}
    std::string parallel_coloring::name() { return "parallel_coloring"; }

    void catalyurek(viennagrid::mesh const & input_mesh, viennagrid::quantity_field & color_field)
    {
        viennamesh::info(1) << "Coloring in parallel with Catalyurek's algorithm" << std::endl;

        //iterate triangles in mesh
        TriangleRange triangles(input_mesh);

        //Iterate over all triangles in the mesh
        viennagrid_element_id * triangle_ids_begin;
        viennagrid_element_id * triangle_ids_end;
        viennagrid_dimension topological_dimension = viennagrid::cell_dimension( input_mesh );	

        viennagrid_element_id * neighbor_begin;
        viennagrid_element_id * neighbor_end;

        viennagrid_dimension connector = 1;

        viennagrid_mesh_elements_get(input_mesh.internal(), topological_dimension, &triangle_ids_begin, &triangle_ids_end);	

        std::vector<int> colors (viennagrid::element_count(input_mesh, 2), -1);
        std::cout << "color size: " << colors.size() << std::endl;

        std::vector<int> work_list(viennagrid::element_count(input_mesh, 2), -1);
        std::iota(work_list.begin(), work_list.end(), 0);

        auto iterations = 0;

        while (work_list.size() > 0)
        {
            std::cout << "RUN " << iterations << std::endl;
            std::cout << "WORK LIST SIZE: " << work_list.size() << std::endl;

            //tentative coloring
            #pragma omp parallel for num_threads(1)
            for (size_t i = 0; i < work_list.size(); ++i)
            {
                std::cerr << omp_get_num_threads() << std::endl;
                viennagrid_element_id *tri = triangle_ids_begin+i;
                int tri_index = viennagrid_index_from_element_id( *tri );

                std::set<int> neigh_colors;

                //std::cout << tri_index << std::endl;

                viennagrid_element_neighbor_elements(input_mesh.internal(), *tri, 0, 2, &neighbor_begin, &neighbor_end);

                for (viennagrid_element_id *n_tri = neighbor_begin; n_tri != neighbor_end; ++n_tri)
                {
                    int n_tri_index = viennagrid_index_from_element_id( *n_tri );

                    //std::cout << "  " << n_tri_index << std::endl;  

                    if (colors[n_tri_index] >= 0 )
                    {
                        //std::cout << "    " << " adding " << colors[n_tri_index] << std::endl;
                        neigh_colors.insert(colors[n_tri_index]);
                    }
                }          

                if (neigh_colors.size() > 0)
                {
                    colors[viennagrid_index_from_element_id( *tri )] = (*neigh_colors.rbegin())+1;
                }

                else    
                {
                    colors[viennagrid_index_from_element_id( *tri )] = 0;
                }
            }//end of tentative coloring

            #pragma omp barrier
            work_list.clear();

            //conflict detection
            /*
            //#pragma omp parallel for
            for (viennagrid_element_id *tri = triangle_ids_begin; tri != triangle_ids_end; ++tri)
            {
                int tri_index = viennagrid_index_from_element_id (*tri);

                viennagrid_element_neighbor_elements(input_mesh.internal(), *tri, 0, 2, &neighbor_begin, &neighbor_end);

                for (viennagrid_element_id *n_tri = neighbor_begin; n_tri != neighbor_end; ++n_tri)
                {
                    int n_tri_index = viennagrid_index_from_element_id( *n_tri );

                    if ( n_tri_index > tri_index && colors[ tri_index ] == colors[n_tri_index] )
                    {
                        work_list.push_back(tri_index);
                    }
                } //end of for loop iterating neighbors
            } //end of conflict detection */

            TriangleRange elements(input_mesh.internal());

            //for (size_t tri = 0; tri < viennagrid::elements<viennagrid::triangle_tag>(input_mesh.internal()).size(); ++tri)
            for (size_t tri = 0; tri < elements.size(); ++tri)
            {
                ;//TriangleNeighborRangeType neighbor_triangles_of_t(input_mesh.internal(), );
            }
/*
            #pragma omp parallel for
            for(size_t i = 0; i < neighbor_triangles_of_t.size(); ++i)
            {
                //EdgeOnVertexRange
            } //end of conflict detection*/
            #pragma omp barrier

            //std::cout << "NEW WORKLIST SIZE: " << work_list.size() << std::endl;
            ++iterations;
        } //end of while loop 

        //DEBUG
    /*    for (size_t i = 0; i < colors.size(); ++i)
        {
            std::cout << i << ": " << colors[i] << std::endl;
        }
        //END OF DEBUG//*/

        //create viennagrid quantity field
        for (size_t i = 0; i < colors.size(); ++i)
        {
            color_field.set(i, colors[i]);
        }
    } //end of catalyurek

    void improved_catalyurek()
    {
        viennamesh::info(1) << "Coloring in parallel with improved Catalyurek's algorithm" << std::endl;
    } //end of improved_catalyurek

    //run function
    bool parallel_coloring::run(viennamesh::algorithm_handle &)
    {
        viennamesh::info(1) << name() << std::endl;

        //get input mesh and create output mesh
        mesh_handle input_mesh = get_required_input<mesh_handle>("mesh");
        mesh_handle output_mesh = make_data<mesh_handle>();
                        
        //get number of vertices and triangles
        int num_vertices = viennagrid::vertex_count(input_mesh());
        int num_triangles = viennagrid::element_count(input_mesh(), 2);
/*
        viennamesh::info(2) << "Number of vertices in mesh : " << num_vertices << std::endl;
        viennamesh::info(2) << "Number of triangles in mesh: " << num_triangles << std::endl;
*/


        viennagrid::quantity_field color_field(2,1);
        color_field.set_name("color");

        catalyurek(input_mesh(), color_field);
        //improved_catalyurek();

        output_mesh = input_mesh;
        quantity_field_handle quantities = make_data<viennagrid::quantity_field>();

        quantities.set(color_field);
   
        set_output("color_field", quantities);
        set_output("mesh", output_mesh());

        return true;
    } //end of bool parallel_coloring::run(viennamesh::algorithm_handle &)

}