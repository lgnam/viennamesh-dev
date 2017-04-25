#ifndef MESH_PARTITIONS_REFINEMENT_HPP
#define MESH_PARTITIONS_REFINEMENT_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "Mesh.h"
#include "Refine.h"
#include "ElementProperty.h"

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
        MeshPartitionsRefinement(Mesh<double>* mesh);           //Constructor
        ~MeshPartitionsRefinement();                            //Destructor

    private:
        Mesh<double>*& input_mesh;

        //TODO: Update those 3 for 3D refinement
        size_t dim = 2;
        size_t nedge = 3;
        size_t nloc = 3;
        size_t msize = 3;
        size_t nthreads = 1;
        size_t splitCnt = 0;
        size_t threadIdx;

        //Refinement Functions
        void refine2D_1(const index_t* newVertex, int eid);
        void refine2D_2(const index_t* newVertex, int eid);
        void refine2D_3(const index_t* newVertex, int eid);

        ElementProperty<double>* property;

        void (MeshPartitionsRefinement::* refineMode2D[3])(const index_t*, int);
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
MeshPartitionsRefinement::MeshPartitionsRefinement(Mesh<double>* mesh) : input_mesh(mesh) 
{
    viennamesh::info(1) << "Calling MeshPartitionsRefinement Constructor" << std::endl;

    //Set Orientation of elements
    for (size_t i = 0; i < input_mesh->get_number_elements(); ++i)
    {
        const int* n = input_mesh->get_element(i);

        //element is marked for deletion
        if (n[0] < 0)
        continue;

        if (dim == 2)
        {
            property = new ElementProperty<double>(input_mesh->get_coords(n[0]), 
            input_mesh->get_coords(n[1]), input_mesh->get_coords(n[2]));
        }

        //TODO: Update for 3D case!!!
    }

    splitCnt = 0;

    //copied from Refine.h Constructor in pragmatic framework
    nthreads = 1;

    // Pre-allocate the maximum size that might be required
    threadIdx=0;

    if(dim==2) 
    {
        //TODO: Update all these functions to work with both partitions and interfaces!!!!
        refineMode2D[0] = &MeshPartitionsRefinement::refine2D_1;
        //refineMode2D[1] = &MeshPartitionsRefinement::refine2D_2;
        //refineMode2D[2] = &MeshPartitionsRefinement::refine2D_3;
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

//MeshPartitions::refine2D_1
void MeshPartitionsRefinement::refine2D_1(const index_t* newVertex, int eid)
{
    //single edge split
    const int* n = input_mesh->get_element(eid);  //works only in test case with partition0
    const int* boundary = input_mesh->get_boundary_ptr(eid*nloc);

    int rotated_ele[3];
    int rotated_boundary[3];
    index_t vertexID = -1;

    for (size_t j = 0; j < 3; ++j)
    {
        if (newVertex[j] >= 0)
        {
            vertexID = newVertex[j];

            rotated_ele[0] = n[j];
            rotated_ele[1] = n[(j+1)%3];
            rotated_ele[2] = n[(j+2)%3];

            rotated_boundary[0] = boundary[j];
            rotated_boundary[1] = boundary[(j+1)%3];
            rotated_boundary[2] = boundary[(j+2)%3];

            break;
        }

        const index_t ele0[] = {rotated_ele[0], rotated_ele[1], vertexID};
        const index_t ele1[] = {rotated_ele[0], vertexID, rotated_ele[2]};

        const index_t ele0_boundary[] = {rotated_boundary[0], 0, rotated_boundary[2]};
        const index_t ele1_boundary[] = {rotated_boundary[0], rotated_boundary[1], 0};

        index_t ele1ID;
        ele1ID = splitCnt;

        // Add rotated_ele[0] to vertexID's NNList
        input_mesh.AddNNList_part(vertexID, rotated_ele[0], part);
        // Add vertexID to rotated_ele[0]'s NNList
        input_mesh.AddNNList_part(rotated_ele[0], vertexID, part);

        // ele1ID is a new ID which isn't correct yet, it has to be
        // updated once each thread has calculated how many new elements
        // it created, so put ele1ID into addNE_fix instead of addNE.
        // Put ele1 in rotated_ele[0]'s NEList
        input_mesh.AddNEList_fix_part(rotated_ele[0], ele1ID, threadIdx, part);

        input_mesh.AddNEList_part(rotated_ele[0], eid, part);
        input_mesh.AddNEList_fix_part(vertexID, ele1ID, threadIdx, part);

        // Replace eid with ele1 in rotated_ele[2]'s NEList
        input_mesh.RemNEList_part(rotated_ele[2], eid, part);
        input_mesh.AddNEList_fix_part(rotated_ele[2], ele1ID, threadIdx,part);

        replace_element(eid, ele0, ele0_boundary, part);
        append_element(ele1, ele1_boundary, part);

        ++splitCnt;
    }
}
//end of MeshPartitions::refine2D_1

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif