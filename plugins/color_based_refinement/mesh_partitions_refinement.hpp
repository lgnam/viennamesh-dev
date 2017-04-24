#ifndef MESH_PARTITIONS_REFINEMENT_HPP
#define MESH_PARTITIONS_REFINEMENT_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "Refine.h"

//viennamesh includes
#include "viennameshpp/plugin.hpp"

//all other includes
#include "mesh_partitions.hpp"

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//class MeshPartitionsRefinement
//
//This class uses MeshPartitions and its coloring to refine the colored partitions in parallel
class MeshPartitionsRefinement{

    public:
        MeshPartitionsRefinement(MeshPartitions& InputMesh);                    //Constructor
        ~MeshPartitionsRefinement();                                            //Destructor

    private:

};

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Helper Functions                                                              //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Implementation                                                                //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//Constructor
//
//Tasks: TODO
MeshPartitionsRefinement::MeshPartitionsRefinement(MeshPartitions& InputMesh, size_t part_id)
{
    viennamesh::info(1) << "Calling MeshPartitionsRefinement Constructor" << std::endl;
    viennamesh::info(5) << "  Refinement with " << InputMesh.get_colors() << " colors." << std::endl;

    //Set Orientation of elements
    for (size_t i = 0; i < InputMesh.pragmatic_partitions[part_id]->get_number_elements(); ++i)
    {
        const int* n = InputMesh.pragmatic_partitions[part_id]->get_element(i);

        //element is marked for deletion
        if (n[0] < 0)
        continue;

        if (dim == 2)
        {
            property = new ElementProperty<double>(InputMesh.pragmatic_partitions[part_id]->get_coords(n[0]), 
            InputMesh.pragmatic_partitions[part_id]->get_coords(n[1]), InputMesh.pragmatic_partitions[part_id]->get_coords(n[2]));
        }

        //TODO: Update for 3D case!!!
    }
}
//end of Constructor

//Destructor
//
//Tasks: TODO
MeshPartitionsRefinement::~MeshPartitionsRefinement()
{
    viennamesh::info(1) << "Calling MeshPartitionsRefinement Destructor" << std::endl;
}
//end of Destructor

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif