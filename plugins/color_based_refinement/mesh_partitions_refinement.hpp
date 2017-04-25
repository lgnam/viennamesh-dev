#ifndef MESH_PARTITIONS_REFINEMENT_HPP
#define MESH_PARTITIONS_REFINEMENT_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "ticker.h"
#include "VTKTools.h"
#include "Mesh.h"
#include "MetricField.h"
#include "Lock.h"
#include "ElementProperty.h"
#include "Edge.h"

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
        MeshPartitionsRefinement(Mesh<double>*& mesh, std::vector<std::set<int>>& nodes_part_ids);           //Constructor
        ~MeshPartitionsRefinement();                            //Destructor

    private:
        Mesh<double>* input_mesh;
        std::vector<std::set<int>> nodes_partition_ids;

        //TODO: Update those 3 for 3D refinement
        size_t dim = 2;
        size_t nedge = 3;
        size_t nloc = 3;
        size_t msize = 3;
        size_t nthreads = 1;
        size_t splitCnt = 0;
        size_t threadIdx;

        //only necessary for comparison with pragmatic
        std::vector<DirectedEdge<index_t>> newVertices;
        std::vector<double> newCoords;
        std::vector<double> newMetric;
        std::vector<index_t> newElements;
        std::vector<int> newBoundaries;
        std::vector<int> newAppearances;
        std::vector<double> newQualities;
        std::vector<index_t> new_vertices_per_element;

        //std::vector<size_t> threadIdx, splitCnt;
        std::vector< DirectedEdge<index_t> > allNewVertices;

        //Refinement Functions
        void refine2D_1(const index_t* newVertex, int eid);
        void refine2D_2(const index_t* newVertex, int eid);
        void refine2D_3(const index_t* newVertex, int eid);

        bool refine(double L_max);
        void refine_edge(index_t n0, index_t n1);
        void refine_element(size_t eid);

        void replace_element(const index_t eid, const index_t *n, const int* boundary);
        void append_element(const index_t* elem, const int* boundary);

        double calculate_quality(const index_t* n);
        double update_quality(index_t element);

        double calc_edge_length(index_t x0, index_t y0);
        double calc_edge_length(index_t x0, index_t y0, double* m);

        size_t edgeNumber(index_t eid, index_t v1, index_t v2);

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
MeshPartitionsRefinement::MeshPartitionsRefinement(Mesh<double>*& mesh, std::vector<std::set<int>>& nodes_part_ids): input_mesh(mesh), 
nodes_partition_ids(nodes_part_ids)

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
        refineMode2D[1] = &MeshPartitionsRefinement::refine2D_2;
        refineMode2D[2] = &MeshPartitionsRefinement::refine2D_3;
    }

    std::cerr << "start refinement" << std::endl;

    refine(sqrt(2.0));

    std::cerr << "refinement done" << std::endl;
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

//MeshPartitionsRefinement::refine
bool MeshPartitionsRefinement::refine(double L_max)
{
    size_t origNNodes = input_mesh->get_number_nodes();
    size_t origNElements = input_mesh->get_number_elements();
    
    size_t splitCnt_wo_vertices = 0;
    splitCnt = 0;

    new_vertices_per_element.resize(nedge*origNElements);
    std::fill(new_vertices_per_element.begin(), new_vertices_per_element.end(), -1);

    /*
    * Average vertex degree in 2D is ~6, so there
    * are approx. (6/2)*NNodes edges in the mesh.
    * In 3D, average vertex degree is ~12.
    */
    size_t reserve_size = nedge*origNNodes;
    newVertices.clear();
    newVertices.reserve(reserve_size);
    newCoords.clear();
    newCoords.reserve(dim*reserve_size);
    newMetric.clear();
    newMetric.reserve(msize*reserve_size);
    newAppearances.clear();
    newAppearances.reserve(reserve_size);

    // Pre-allocate the maximum size that might be required
    allNewVertices.resize(origNElements*nloc);
    
    //Traverse all mesh edges and select them for refinement if length is greater than L_max in metric space
    //i is local index!!!
    for (size_t i = 0; i < origNNodes; ++i)
    {
        //DONT TOUCH BOUNDARY ELEMENTS!!!
        /*if (nodes_partition_ids[i].size() > 1) 
            continue;
*/
        std::vector<index_t> NNList_i = input_mesh->get_nnlist(i);

        for (auto otherVertex : NNList_i) //otherVertex is already a global index
        {
            if (i<otherVertex)
            {
                double length = calc_edge_length(i, otherVertex);

                if (length > L_max)
                {
                    ++splitCnt;
                    refine_edge(i, otherVertex);
                }
            }
        }
    }

    //TODO: do omp atomic capture here!!!
    threadIdx = input_mesh->get_number_nodes();
    input_mesh->set_nnodes(threadIdx + splitCnt); 

    //Resize meshes vectors
    input_mesh->resize_vectors(dim);

    //Append new coords, metric and appearance data to the mesh
    input_mesh->add_coords(threadIdx, newCoords, dim, splitCnt);
    input_mesh->add_metric(threadIdx, newMetric, msize, splitCnt);

    //Fix IDs of new vertices
    for (size_t j = 0; j < splitCnt; ++j)
    {
        newVertices[j].id = threadIdx+j;
    }

    // Accumulate all newVertices in a contiguous array
    memcpy(&allNewVertices[threadIdx-origNNodes], &newVertices[0], newVertices.size() * sizeof(DirectedEdge<index_t>));

    // Mark each element with its new vertices,
    // update NNList for all split edges.
    for (size_t j = 0; j < splitCnt; j++)
    {
        index_t vid = allNewVertices[j].id; 
        index_t firstid = allNewVertices[j].edge.first; 
        index_t secondid = allNewVertices[j].edge.second; 

        //Find which elements share this edge and mark them with their new vertices
        //get NELists for firstid and secondid
        std::set<index_t> NEList_first = input_mesh->get_nelist(firstid);
        std::set<index_t> NEList_second = input_mesh->get_nelist(secondid); 

        std::set<index_t> intersection;
        std::set_intersection(NEList_first.begin(), NEList_first.end(),
                          NEList_second.begin(), NEList_second.end(),
                          std::inserter(intersection, intersection.begin()));

        for(auto element : intersection)
        {
            //index_t eid = element; 
            size_t edgeOffset = edgeNumber(element, firstid, secondid);
            new_vertices_per_element[nedge * element + edgeOffset] = vid;
        }

        //Update NNList for newly created vertices
        input_mesh->add_nnlist(vid, firstid);
        input_mesh->add_nnlist(vid, secondid);

        input_mesh->remove_nnlist(firstid, secondid);
        input_mesh->add_nnlist(firstid, vid);
        input_mesh->remove_nnlist(secondid, firstid);
        input_mesh->add_nnlist(secondid, vid);
    } //end of for loop over splitCnt

    //Start element refinement
    splitCnt = 0;
    newElements.clear();
    newBoundaries.clear();
    newQualities.clear();
    newElements.reserve(dim*dim*origNElements);
    newBoundaries.reserve(dim*dim*origNElements);
    newQualities.reserve(origNElements);

    threadIdx = input_mesh->get_number_elements();
    std::cerr<<"element refinement for " << origNElements << " elements" << std::endl;
    //for loop over originalElements
    for (size_t eid = 0; eid < origNElements; ++eid)
    {
        //if element has been deleted, continue
        const index_t* n = input_mesh->get_element(eid);

        if (n[0] < 0)
            continue;

        for (size_t j = 0; j < nedge; ++j)
        {
            if (new_vertices_per_element[nedge*eid+j] != -1)
            {
                std::cerr << eid << std::endl;
                refine_element(eid);
                break;
            }
        }
    }
    //end of for loop over origalElements

    std::cerr << "TEST" << std::endl;

    //TODO: do omp atomic capture here!!!
    threadIdx = input_mesh->get_number_elements();
    input_mesh->set_nelements(threadIdx + splitCnt);

    //Should we resize???
    if( input_mesh->get_element_node_size() < input_mesh->get_number_elements() * nloc) 
    {
        std::cerr << "RESIZE ENLIST, BOUNDARY AND QUALITY HAS TO BE IMPLEMENTED!!!" << std::endl;
    }

    //Append new elements to the mesh
    input_mesh->add_elements(newElements, splitCnt, threadIdx);
    input_mesh->add_boundaries(newBoundaries, splitCnt, threadIdx);
    input_mesh->add_qualities(newQualities, splitCnt, threadIdx);

    //fix orientation of new elements
    size_t NElements = input_mesh->get_number_elements();

    //for loop over NElements to fix their orientation
    //TODO: do i need this???
    /*
    for (size_t j = 0; j < NElements; ++j)
    {
        const index_t* n = input_mesh->get_element(j);
        std::vector<index_t> _ENList2 = input_mesh->get_element_node();

        //check if element is marked for deletion
        if (n[0] < 0)
        {
            continue;
        }
    }
    //end of for loop over NElements for fixing their orientation
    */

    return true;
}
//end of MeshPartitionsRefinement::refine

//MeshPartitionsRefinement::refine_edge
void MeshPartitionsRefinement::refine_edge(index_t n0, index_t n1)
{  
  //swap indices because lesser id has to be the first
  if(n0 > n1)
  {
    index_t tmp = n0;
    n0 = n1;
    n1 = tmp;
  }

  //newVertices.push_back( DirectedEdge<index_t> (mesh.g2l_vertex_map_partitions[part].at(n0), mesh.g2l_vertex_map_partitions[part].at(n1)) );  //n0 and n1 are global indices    
  newVertices.push_back( DirectedEdge<index_t> (n0, n1));

  //Calculate position of new point
  double x, m;
  double n0_coords[2], n1_coords[2];
  double m0[3], m1[3];
  
  input_mesh->get_coords(n0, n0_coords);
  input_mesh->get_metric(n0, m0);

  input_mesh->get_coords(n1, n1_coords);
  input_mesh->get_metric(n1, m1);

  double weight = 1.0 / (1.0 + sqrt( calc_edge_length(n0, n1, m0) / calc_edge_length(n0, n1, m1) ) );

  //calculate position of new vertex
  for (size_t i = 0; i < dim; ++i)
  {
    x = n0_coords[i] + weight * (n1_coords[i] - n0_coords[i]);
    newCoords.push_back(x);
  }

  //Interpolate new metric
  for (size_t i = 0; i < msize; ++i)
  {
    m = m0[i] + weight * (m1[i] - m0[i]);
    newMetric.push_back(m);    
  }
}
//end of MeshPartitionsRefinement::refine_edge

//MeshPartitionsRefinement::refine_element(size_t eid)
//TODO: only for 2D implemented up to now!!!
void MeshPartitionsRefinement::refine_element(size_t eid)
{
  //works only for test case with partition0

  const int* n = input_mesh->get_element(eid);
  
  // Note the order of the edges - the i'th edge is opposite the i'th node in the element.
  index_t newVertex[3] = {-1, -1, -1};
  newVertex[0] = new_vertices_per_element[nedge*eid];
  newVertex[1] = new_vertices_per_element[nedge*eid+1];
  newVertex[2] = new_vertices_per_element[nedge*eid+2];
/*
  std::cerr << eid << std::endl;
  std::cerr << " " << newVertex[0] << " " << newVertex[1] << " " << newVertex[2] << std::endl;
*/
  int refine_cnt = 0;
  for(size_t i = 0; i < 3; ++i)
  {
    if (newVertex[i] != -1)
    {
      ++refine_cnt;
    }
  }

  if (refine_cnt > 0)
  {
    (this->*refineMode2D[refine_cnt-1])(newVertex, eid);
  }
}
//end of MeshPartitionsRefinement::refine_element(size_t ei)

//MeshPartitions::refine2D_1
void MeshPartitionsRefinement::refine2D_1(const index_t* newVertex, int eid)
{
    //std::cerr << "refine2D_1 " << eid << std::endl;

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
        input_mesh->add_nnlist(vertexID, rotated_ele[0]);
        // Add vertexID to rotated_ele[0]'s NNList
        input_mesh->add_nnlist(rotated_ele[0], vertexID);

        // ele1ID is a new ID which isn't correct yet, it has to be
        // updated once each thread has calculated how many new elements
        // it created, so put ele1ID into addNE_fix instead of addNE.
        // Put ele1 in rotated_ele[0]'s NEList
        input_mesh->add_nelist_fix(rotated_ele[0], ele1ID, threadIdx);
        //input_mesh->add_nelist(rotated_ele[0], ele1ID);
        input_mesh->add_nelist(rotated_ele[0], eid);
        input_mesh->add_nelist_fix(vertexID, ele1ID, threadIdx);
        //input_mesh->add_nelist(vertexID, ele1ID);
 std::cerr << "TEST" << std::endl;
        // Replace eid with ele1 in rotated_ele[2]'s NEList
        input_mesh->remove_nelist(rotated_ele[2], eid);
         std::cerr << "TEST" << std::endl;
        //input_mesh->add_nelist_fix(rotated_ele[2], ele1ID, threadIdx);
        input_mesh->add_nelist(rotated_ele[2], ele1ID);
         std::cerr << "TEST" << std::endl;

        replace_element(eid, ele0, ele0_boundary);
        append_element(ele1, ele1_boundary);

        ++splitCnt;
    }
}
//end of MeshPartitions::refine2D_1

//MeshPartitionsRefinement::refine2D_2
void MeshPartitionsRefinement::refine2D_2(const index_t* newVertex, int eid)
{
    //std::cerr << "refine2D_2 " << eid << std::endl;

    const int* n = input_mesh->get_element(eid);  //works only in test case with partition0
    const int* boundary = input_mesh->get_boundary_ptr(eid*nloc);

    int rotated_ele[3];
    int rotated_boundary[3];
    index_t vertexID[2];
    for(int j=0; j<3; j++) {
        if(newVertex[j] < 0) {
            vertexID[0] = newVertex[(j+1)%3];
            vertexID[1] = newVertex[(j+2)%3];

            rotated_ele[0] = n[j];
            rotated_ele[1] = n[(j+1)%3];
            rotated_ele[2] = n[(j+2)%3];

            rotated_boundary[0] = boundary[j];
            rotated_boundary[1] = boundary[(j+1)%3];
            rotated_boundary[2] = boundary[(j+2)%3];

            break;
        }
    }

    double ldiag0 = calc_edge_length(rotated_ele[1], vertexID[0]);
    double ldiag1 = calc_edge_length(rotated_ele[2], vertexID[1]);

    const int offset = ldiag0 < ldiag1 ? 0 : 1;

    const index_t ele0[] = {rotated_ele[0], vertexID[1], vertexID[0]};
    const index_t ele1[] = {vertexID[offset], rotated_ele[1], rotated_ele[2]};
    const index_t ele2[] = {vertexID[0], vertexID[1], rotated_ele[offset+1]};

    const index_t ele0_boundary[] = {0, rotated_boundary[1], rotated_boundary[2]};
    const index_t ele1_boundary[] = {rotated_boundary[0], (offset==0)?rotated_boundary[1]:0, (offset==0)?0:rotated_boundary[2]};
    const index_t ele2_boundary[] = {(offset==0)?rotated_boundary[2]:0, (offset==0)?0:rotated_boundary[1], 0};

    index_t ele0ID, ele2ID;
    ele0ID = splitCnt;
    ele2ID = ele0ID+1;

    // NNList: Connect vertexID[0] and vertexID[1] with each other
    input_mesh->add_nnlist(vertexID[0], vertexID[1]);
    input_mesh->add_nnlist(vertexID[1], vertexID[0]);

    // vertexID[offset] and rotated_ele[offset+1] are the vertices on the diagonal
    input_mesh->add_nnlist(vertexID[offset], rotated_ele[offset+1]);
    input_mesh->add_nnlist(rotated_ele[offset+1], vertexID[offset]);

    // rotated_ele[offset+1] is the old vertex which is on the diagonal
    // Add ele2 in rotated_ele[offset+1]'s NEList
    input_mesh->add_nelist_fix(rotated_ele[offset+1], ele2ID, threadIdx);

    // Replace eid with ele0 in NEList[rotated_ele[0]]
    input_mesh->remove_nelist(rotated_ele[0], eid);
    input_mesh->add_nelist_fix(rotated_ele[0], ele0ID, threadIdx);

    // Put ele0, ele1 and ele2 in vertexID[offset]'s NEList
    input_mesh->add_nelist(vertexID[offset], eid);
    input_mesh->add_nelist_fix(vertexID[offset], ele0ID, threadIdx);
    input_mesh->add_nelist_fix(vertexID[offset], ele2ID, threadIdx);

    // vertexID[(offset+1)%2] is the new vertex which is not on the diagonal
    // Put ele0 and ele2 in vertexID[(offset+1)%2]'s NEList
    input_mesh->add_nelist_fix(vertexID[(offset+1)%2], ele0ID, threadIdx);
    input_mesh->add_nelist_fix(vertexID[(offset+1)%2], ele2ID, threadIdx);

    replace_element(eid, ele1, ele1_boundary);
    append_element(ele0, ele0_boundary);
    append_element(ele2, ele2_boundary);

    splitCnt += 2;
}
//end of MeshPartitionsRefinement::refine2D_2

//MeshPartitionsRefinement::refine2D_3
void MeshPartitionsRefinement::refine2D_3(const index_t* newVertex, int eid)
{
    //std::cerr << "refine2D_3 " << eid << std::endl;

    //TODO: change to work with any partition!!!
    const int *n=input_mesh->get_element(eid);
    const int* boundary = input_mesh->get_boundary_ptr(eid*nloc);

    const index_t ele0[] = {n[0], newVertex[2], newVertex[1]};
    const index_t ele1[] = {n[1], newVertex[0], newVertex[2]};
    const index_t ele2[] = {n[2], newVertex[1], newVertex[0]};
    const index_t ele3[] = {newVertex[0], newVertex[1], newVertex[2]};

    const int ele0_boundary[] = {0, boundary[1], boundary[2]};
    const int ele1_boundary[] = {0, boundary[2], boundary[0]};
    const int ele2_boundary[] = {0, boundary[0], boundary[1]};
    const int ele3_boundary[] = {0, 0, 0};

    index_t ele1ID, ele2ID, ele3ID;
    ele1ID = splitCnt;
    ele2ID = ele1ID+1;
    ele3ID = ele1ID+2;

    //UpdateNNList
    input_mesh->add_nnlist(newVertex[0], newVertex[1]);
    input_mesh->add_nnlist(newVertex[0], newVertex[2]);
    input_mesh->add_nnlist(newVertex[1], newVertex[0]);
    input_mesh->add_nnlist(newVertex[1], newVertex[2]);
    input_mesh->add_nnlist(newVertex[2], newVertex[0]);
    input_mesh->add_nnlist(newVertex[2], newVertex[1]);

    //UpdateNEList  
    input_mesh->remove_nelist(n[1], eid);
    input_mesh->add_nelist_fix(n[1], ele1ID, threadIdx);
    input_mesh->remove_nelist(n[2], eid);
    input_mesh->add_nelist_fix(n[2], ele2ID, threadIdx);

    input_mesh->add_nelist_fix(newVertex[0], ele1ID, threadIdx);
    input_mesh->add_nelist_fix(newVertex[0], ele2ID, threadIdx);
    input_mesh->add_nelist_fix(newVertex[0], ele3ID, threadIdx);

    input_mesh->add_nelist(newVertex[1], eid);
    input_mesh->add_nelist_fix(newVertex[1], ele2ID, threadIdx);
    input_mesh->add_nelist_fix(newVertex[1], ele3ID, threadIdx);

    input_mesh->add_nelist(newVertex[2], eid);
    input_mesh->add_nelist_fix(newVertex[2], ele1ID, threadIdx);
    input_mesh->add_nelist_fix(newVertex[2], ele3ID, threadIdx);

    replace_element(eid, ele0, ele0_boundary);
    append_element(ele1, ele1_boundary);
    append_element(ele2, ele2_boundary);
    append_element(ele3, ele3_boundary);
}
//end of MeshPartitionsRefinement::refine2D_3

//MeshPartitionsRefinement::replace_element_part(index_t eid, int* boundary)
void MeshPartitionsRefinement::replace_element(const index_t eid, const index_t *n, const int* boundary)
{
  for (size_t i = 0; i < nloc; ++i)
  {
    input_mesh->set_element(eid, n[i], boundary[i], i);
  }

  update_quality(eid);
}
//end of MeshPartitionsRefinement::replace_element_part()

//MeshPartitionsRefinement::append_element(index_t eid, index_t *n, int* boundary)
void MeshPartitionsRefinement::append_element(const index_t* elem, const int* boundary)
{
  for (size_t i = 0; i < nloc; ++i)
  {
    newElements.push_back(elem[i]);
    newBoundaries.push_back(boundary[i]);
  }

  double q = calculate_quality(elem);
  
  newQualities.push_back(q);

  //TODO:
  //update l2g and g2l element mappings
  /*
  size_t origNElements = input_mesh->get_number_elements();
  size_t globalNElements = mesh.num_elements;

  mesh.l2g_element_map_partitions[part].insert( std::make_pair(origNElements+splitCnt, globalNElements+splitCnt) );
  mesh.g2l_element_map_partitions[part].insert( std::make_pair(globalNElements+splitCnt, origNElements+splitCnt) );
 */
  ++splitCnt;
}
//end of MeshPartitionsRefinement::append_element_part()

//MeshPartitionsRefinement::calculate_quality_part(const index_t* n, int part)
double MeshPartitionsRefinement::calculate_quality(const index_t* n)
{

  //works only for test case with partition0
  const double *x0 = input_mesh->get_coords(n[0]);
  const double *x1 = input_mesh->get_coords(n[1]);
  const double *x2 = input_mesh->get_coords(n[2]);

  const double *m0 = input_mesh->get_metric(n[0]);
  const double *m1 = input_mesh->get_metric(n[1]);
  const double *m2 = input_mesh->get_metric(n[2]);

  ElementProperty<double>* prop;
  input_mesh->get_property(prop);

  return prop->lipnikov(x0, x1, x2, m0, m1, m2);
}
//end of MeshPartitionsRefinement::calculate_quality_part(const index_t* n, int part)

//MeshPartitionsRefinement::update_quality(const index_t* element)
double MeshPartitionsRefinement::update_quality(index_t element)
{
  const index_t* n = input_mesh->get_element(element);

  input_mesh->set_quality(element, calculate_quality(n));
}
//end of MeshPartitionsRefinement::update_quality(const index_t* n)

//MeshPartitions::calc_edge_length
double MeshPartitionsRefinement::calc_edge_length(index_t x0, index_t y0)
{
  double length = -1.0;

  //TODO: adapt for 3D space!
  double m[3];
  double m_x[3], m_y[3];

  input_mesh->get_metric(x0, m_x);
  input_mesh->get_metric(y0, m_y);

  m[0] = (m_x[0] + m_y[0])*0.5;
  m[1] = (m_x[1] + m_y[1])*0.5;
  m[2] = (m_x[2] + m_y[2])*0.5;

  //calculate length
  //taken from pragmatic
  /*! Length of an edge as measured in metric space.
   *
   * @param x0 coordinate at start of line segment.
   * @param x1 coordinate at finish of line segment.
   * @param m metric tensor for first point.
   */
  double x_coords[2] {-1.0, -1.0}, y_coords[2] {-1.0, -1.0};

  input_mesh->get_coords(x0, x_coords);
  input_mesh->get_coords(y0, y_coords);

  double x = x_coords[0] - y_coords[0];
  double y = x_coords[1] - y_coords[1];

  // The l-2 norm can fail for anisotropic metrics. In such cases use the l-inf norm.
  double l2 = (m[1]*x + m[2]*y)*y + (m[0]*x + m[1]*y)*x;

  if(l2>0) 
  {
      length = sqrt(l2);
  }
 
  else 
  {
      double linf = std::max((m[2]*y)*y, (m[0]*x)*x);
      length = sqrt(linf);
  }

  return length;
}
//end of MeshPartitions::calc_edge_length

//MeshPartitionsRefinement::calc_edge_length(index x0, index y0, double* m)
//calculates edge length in metric space indices and metric array
double MeshPartitionsRefinement::calc_edge_length(index_t x0, index_t y0, double* m)
{
  double length = -1.0;

  //TODO: adapt for 3D space!

  //calculate length
  //taken from pragmatic
  /*! Length of an edge as measured in metric space.
   *
   * @param x0 coordinate at start of line segment.
   * @param x1 coordinate at finish of line segment.
   * @param m metric tensor for first point.
   */
  double x_coords[2] {-1.0, -1.0}, y_coords[2] {-1.0, -1.0};

  input_mesh->get_coords(x0, x_coords);
  input_mesh->get_coords(y0, y_coords);

  double x = x_coords[0] - y_coords[0];
  double y = x_coords[1] - y_coords[1];

  // The l-2 norm can fail for anisotropic metrics. In such cases use the l-inf norm.
  double l2 = (m[1]*x + m[2]*y)*y + (m[0]*x + m[1]*y)*x;

  if(l2>0) 
  {
      length = sqrt(l2);
  }
 
  else 
  {
      double linf = std::max((m[2]*y)*y, (m[0]*x)*x);
      length = sqrt(linf);
  }

  return length;
}
//end of MeshPartitionsRefinement::calc_edge_length(index x0, index y0, double* m)

//MeshPartitionsRefinement::edgeNumber
size_t MeshPartitionsRefinement::edgeNumber(index_t eid, index_t v1, index_t v2)
{
  const int *n = input_mesh->get_element(eid);

  //if(dim==2) 
  {
      /* In 2D:
        * Edge 0 is the edge (n[1],n[2]).
        * Edge 1 is the edge (n[0],n[2]).
        * Edge 2 is the edge (n[0],n[1]).
        */
      if(n[1]==v1 || n[1]==v2) 
      {
        if(n[2]==v1 || n[2]==v2)
          return 0;
        else
          return 2;
      } 
      else
        return 1;
  }
}
//end of MeshPartitionsRefinement::edgeNumber
//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif