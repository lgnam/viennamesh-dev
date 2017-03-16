#ifndef GROUPED_PARTITIONS_REFINEMENT_HPP
#define GROUPED_PARTITIONS_REFINEMENT_HPP

//pragmatic basic includes for mesh data structure
#include "ticker.h"
#include "VTKTools.h"
#include "Mesh.h"
#include "MetricField.h"
#include "Lock.h"
#include "ElementProperty.h"
#include "Edge.h"

//other includes

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//
class GroupedPartitionsRefinement
{
  public:
    GroupedPartitionsRefinement(GroupedPartitions &GP, int parti);              //Constructor
    ~GroupedPartitionsRefinement();                                  //Destructor

  private:
    GroupedPartitions &mesh;

    void CreateSegmentList();

    double CreateEncroachedSegmentsList();
    double GetTriangleQuality();
    bool CheckDelaunay();

    void PragmaticRefinement(double L_max);
    void GetNNList_assigned(index_t i, std::vector<index_t>& NNList);
    void GetNEList_assigned(index_t i, std::set<index_t>& NEList);
    void RefinementKernel(int part1, int part2, int inter, double L_max);
    void RefinementKernel_2(int part1, int part2, int inter, double L_max);
    void RefinementKernel_2_part(int part, double L_max);

    //struct for using with segments set
    struct _edges : public std::pair<int, int>
    { 
      int a,b;
      _edges(int a, int b) : std::pair<int, int>(a < b ? a : b, a < b ? b : a) 
      {}
    }; //end of struct

    std::set<_edges> segments;
    std::set<_edges> encroached_segments;
    std::vector<std::set<int>> debug_NEList;

    //only necessary for comparison with pragmatic
    std::vector<DirectedEdge<index_t>> newVertices;
    std::vector<double> newCoords;
    std::vector<double> newMetric;
    std::vector<index_t> newElements;
    std::vector<int> newBoundaries;
    std::vector<int> newAppearances;
    std::vector<double> newQualities;
    std::vector<index_t> new_vertices_per_element;

//    std::vector<size_t> threadIdx, splitCnt;
    std::vector< DirectedEdge<index_t> > allNewVertices;

    size_t threadIdx;
    //size_t splitCnt;

    double calc_edge_length(index_t x, index_t y);
    double calc_edge_length(index_t x, index_t y, double* m);
    void refine_edge_part(index_t n0, index_t n1, int part);
    void refine_element_part(size_t eid, int part);
    size_t edgeNumber_part(index_t eid, index_t v1, index_t v2, int part);

    //Refinement Functions
    void refine2D_1_part(const index_t* newVertex, int eid, int part);
    void refine2D_2_part(const index_t* newVertex, int eid, int part);
    void refine2D_3_part(const index_t* newVertex, int eid, int part);

    void replace_element_part(const index_t eid, const index_t *n, const int* boundary, int part);
    void append_element_part(const index_t* elem, const int* boundary, int part);
  
    //TODO: Update those 3 for 3D refinement
    size_t dim = 2;
    size_t nedge = 3;
    size_t nloc = 3;
    size_t msize = 3;
    size_t nthreads = 1;
    size_t splitCnt = 0;

    ElementProperty<double>* property;

    void (GroupedPartitionsRefinement::* refineMode2D[3])(const index_t*, int, int);

    size_t rounds = 0;
}; //end of class

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Helper Functions                                                              //
//----------------------------------------------------------------------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Implementation                                                                //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//Constructor
GroupedPartitionsRefinement::GroupedPartitionsRefinement(GroupedPartitions &GP, int parti) : mesh(GP)
{
  std::cout << std::endl << "Create GroupedPartitionsRefinement Object for patition " << parti << std::endl;
  debug_NEList.resize(mesh.num_nodes);

  //TODO update to work with both partitions and interfaces
  size_t NElements = mesh.pragmatic_partitions[parti]->get_number_elements();
  //std::cerr << "NElements " << NElements << std::endl;

  property = NULL;

  //Set orientation of elements
  for (size_t i =0; i < NElements; ++i)
  {
    //TODO update to work with both partitions and interfaces
    const int* n = mesh.pragmatic_partitions[parti]->get_element(i);

    //element is marked for deletion
    if (n[0] < 0)
      continue;

    if (dim == 2)
    {
      //TODO update to work with both partitions and interfaces
      property = new ElementProperty<double>(mesh.pragmatic_partitions[parti]->get_coords(n[0]), mesh.pragmatic_partitions[parti]->get_coords(n[1]), mesh.pragmatic_partitions[parti]->get_coords(n[2]));
    }
  }

  splitCnt = 0;

  //copied from Refine.h Constructor in pragmatic framework
  //nthreads = pragmatic_nthreads();
  nthreads = 1;
  //newVertices.resize(nthreads);
  //newElements.resize(nthreads);
  //newBoundaries.resize(nthreads);
 // newQualities.resize(nthreads);
 // newCoords.resize(nthreads);
  //newMetric.resize(nthreads);

  // Pre-allocate the maximum size that might be required
  //allNewVertices.resize(mesh._ENList.size());
  threadIdx=0;

  if(dim==2) 
  {
    //TODO: Update all these functions to work with both partitions and interfaces!!!!
    refineMode2D[0] = &GroupedPartitionsRefinement::refine2D_1_part;
    refineMode2D[1] = &GroupedPartitionsRefinement::refine2D_2_part;
    refineMode2D[2] = &GroupedPartitionsRefinement::refine2D_3_part;
  }
/*
  threadIdx.resize(nthreads);
  splitCnt.resize(nthreads);
*/
  //PragmaticRefinement(sqrt(2));
/*
  std::vector<double> _metric = mesh.pragmatic_partitions[0]->get_metric_vector();

  std::cout << "metric size " << _metric.size() << std::endl;

  for (auto it : _metric)
  {
      std::cout << it << std::endl;
  }
  std::cout << std::endl;
*/

 /* //debug
  ofstream nelist_debug_before;
  nelist_debug_before.open("nelist_debug_before.txt");

  std::vector<std::set<int>> NEL_before = mesh.pragmatic_partitions[0]->get_node_element();
  for(size_t i = 0; i < mesh.pragmatic_partitions[0]->get_number_nodes(); ++i)
  {
    nelist_debug_before << i << std::endl;
    for (auto it : NEL_before[i])
    {
      nelist_debug_before << " " << it << std::endl;
    }
  }
  nelist_debug_before.close();
  //end of debug*/
  auto wall_tic = std::chrono::system_clock::now();
  for (int iters = 0; iters < 3; ++iters, ++rounds)
  {
    std::cout <<"ROUND " << iters+1 << std::endl;
    RefinementKernel_2(parti, 1, 0, sqrt(2));
  }
  std::chrono::duration<double> wall_clock_duration = std::chrono::system_clock::now() - wall_tic;
  std::cout << "RefinementKernel_2: " << wall_clock_duration.count() << std::endl;
 
 /* //debug
  std::cerr << "afterwards " << mesh.pragmatic_partitions[0]->get_number_nodes() << " " << mesh.pragmatic_partitions[0]->get_number_elements() << std::endl;

  ofstream nelist_debug;
  nelist_debug.open("nelist_debug.txt");
*//*
  std::vector<std::set<int>> NEL = mesh.pragmatic_partitions[0]->get_node_element();
  for(size_t i = 0; i < mesh.pragmatic_partitions[0]->get_number_nodes(); ++i)
  {
    std::cerr << "NEList for " << i << std::endl;
    for (auto it : NEL[i])
    {
      std::cerr << " " << it << std::endl;
    }
  }
*/
  //nelist_debug.close();
/*
  ofstream nnlist_debug;
  nnlist_debug.open("nnlist_debug.txt");

  std::vector<std::vector<index_t>> nnl = mesh.pragmatic_partitions[0]->copy_nnlist();
  for (size_t i = 0; i < mesh.pragmatic_partitions[0]->get_number_nodes(); ++i)
  {
    nnlist_debug << i << std::endl;
    for (auto it : nnl[i])
    {
      nnlist_debug << " " << it << std::endl;
    }
  }
  nnlist_debug.close();
  //end of debug*/

  //defrag pragmatic data structures
  mesh.pragmatic_partitions[parti]->defragment();
  //mesh.pragmatic_partitions[1]->defragment();

  std::string export_file = "examples/data/pragmatic_metis_partitioning/test_refine_2d_part";
  export_file += std::to_string(parti);

  //export to files
  std::cout << "writing data to: " << export_file << std::endl;
  //std::cout << "nodes:    " << mesh.pragmatic_partitions[parti]->get_number_nodes() << std::endl;
  //std::cout << "elements: " << mesh.pragmatic_partitions[parti]->get_number_elements() << std::endl;
  VTKTools<double>::export_vtu(export_file.c_str(), mesh.pragmatic_partitions[parti]);
  //VTKTools<double>::export_vtu("examples/data/pragmatic_metis_partitioning/test_refine_2d_part1", mesh.pragmatic_partitions[1]);

  //CreateEncroachedSegmentsList();
/*  
  CreateSegmentList();  
  CheckDelaunay();
*/
}
//end of Constructor

//Destructor
GroupedPartitionsRefinement::~GroupedPartitionsRefinement()
{
  //debug
  ofstream nelist_out;
  nelist_out.open("NEList_fetched.txt");
  int ne_elements = 0;
  for (size_t i = 0; i < debug_NEList.size(); ++i)
  {
    nelist_out << i << ": ";

    for (auto vertex : debug_NEList[i])
    {
      nelist_out << " " << vertex << " ";
      ++ne_elements;
    }
    nelist_out << std::endl;
  }
  nelist_out.close();
  //end of debug

 // std::cout << "NEElements_fetched: " << ne_elements << std::endl;
  std::cout << "Delete GroupedPartitionsRefinement Object" << std::endl << std::endl;
}
//end of Destructor

//GroupedPartitionsRefinement::CreateSegmentList()
void GroupedPartitionsRefinement::CreateSegmentList()
{
  std::cout << "  CreateSegmentList" << std::endl;

 // std::set<_edges> segments;

  //double goodAngle = 0.36;

  //generate segments list
  for (size_t i = 0; i < mesh.num_elements; ++i)
  {
    //get vertices of element (use element-node-list implemented in pragmatic)
    int x0 = mesh._ENList[3*i]; 
    int x1 = mesh._ENList[3*i+1];
    int x2 = mesh._ENList[3*i+2];
/*
    double coords_x0[2] {0.0, 0.0};
    double coords_x1[2] {0.0, 0.0};
    double coords_x2[2] {0.0, 0.0};

    mesh.GetCoords(x0, coords_x0);
    mesh.GetCoords(x1, coords_x1);
    mesh.GetCoords(x2, coords_x2);
*/
    segments.insert(_edges(x0, x1));
    segments.insert(_edges(x0, x2));
    segments.insert(_edges(x1, x2));
   
    //check for encroachment 
    //  
    //check all neighboring vertices if they encroach the two vertices of a segment
/*
    //TODO: use InCircleTest() to test all neighboring vertices for encroachment
    std::vector<index_t> _NNList_x0;
    std::vector<index_t> _NNList_x1;
    std::vector<index_t> _NNList_x2;
    
      
    //end of check for encroachment

/*
    //check if edge is a boundary edge, if yes omit refinement
    //TODO: allow also refinement of boundary edges
    if (!mesh.boundary_nodes_mesh[x0] || !mesh.boundary_nodes_mesh[x1])
    {
      //check for encroachment
      //
      //use the third vertex for this approach
      //use dotproduct to check for encroachment (see code in triangle!!!)
      double dotproduct = (coords_x0[0] - coords_x2[0]) * (coords_x1[0] - coords_x2[0]) + (coords_x0[1] - coords_x2[1]) * (coords_x1[1] - coords_x2[2]);
    //( pow( (coords_x0[0] - coords_x2[0]), 2) + pow ( (coords_x0[1] - coords_x2[1]), 2) ) * ( pow( (coords_x1[0] - coords_x2[0]), 2) + pow ( (coords_x1[1] - coords_x2[1]), 2) );    
      if (dotproduct < 0.0)
      {
        if (dotproduct * dotproduct >= pow((2* goodAngle -1.0), 2) * ( pow ((coords_x0[0] - coords_x2[0]), 2) + pow((coords_x0[1] - coords_x2[1]), 2) ) * ( pow((coords_x1[0] - coords_x2[0]) , 2) + 
            pow((coords_x1[1] - coords_x2[2]), 2) ) )
        {
          std::cout << x0 << " " << x1 << " " << x2 << std::endl;
          segments.insert(_edges(x0, x1));
        }
      }
    }
  
    if (!mesh.boundary_nodes_mesh[x0] || !mesh.boundary_nodes_mesh[x2])
    {
      //check for encroachment
      //
      //use the third vertex for this approach
      //use dotproduct to check for encroachment (see code in triangle!!!)
      double dotproduct = (coords_x0[0] - coords_x2[0]) * (coords_x1[0] - coords_x2[0]) + (coords_x0[1] - coords_x2[1]) * (coords_x1[1] - coords_x2[2]);

      //segments.insert(_edges(x0, x2));
    }

    if (!mesh.boundary_nodes_mesh[x1] || !mesh.boundary_nodes_mesh[x2])
    {
      //check for encroachment
      //
      //use the third vertex for this approach
      //use dotproduct to check for encroachment (see code in triangle!!!)
      double dotproduct = (coords_x0[0] - coords_x2[0]) * (coords_x1[0] - coords_x2[0]) + (coords_x0[1] - coords_x2[1]) * (coords_x1[1] - coords_x2[2]);

      //segments.insert(_edges(x1, x2));
    }
    //end of boundary check
/*
    std::cout << " " << i << std::endl;
    std::cout << "  " << std::get<0>(_edges(x0, x1)) << " " << std::get<1>(_edges(x0, x1)) << std::endl;
    std::cout << "  " << std::get<0>(_edges(x0, x2)) << " " << std::get<1>(_edges(x0, x2)) << std::endl;
    std::cout << "  " << std::get<0>(_edges(x1, x2)) << " " << std::get<1>(_edges(x1, x2)) << std::endl;
*/
  } //end of generate segments list
/*
  //debug output
  int ctr = 0;
  for (auto it : segments)
  {
    std::cout << " " << ctr << std::endl;
    std::cout << "  " << std::get<0>(it) << " " << std::get<1>(it) << std::endl;
    //std::cout << "  " << it.first << " " << it.second << std::endl;
    ++ctr;
  }
*/
}
//end of GroupedPartitionsRefinement::CreateSegmentList()

//GroupedPartitionsRefinement::CreateEncroachedSegmentsList()
//TODO: creates presently only a list of segments!!!
double GroupedPartitionsRefinement::CreateEncroachedSegmentsList()
{
  std::cout << "  CreateEncroachedSegmentsList()" << std::endl;

  //iterate over all triangles to get the segments
  for (size_t i = 0; i < mesh.num_elements; ++i)
  {
    //get vertices of element (use element-node-list implemented in pragmatic)
    int x0 = mesh._ENList[3*i]; 
    int x1 = mesh._ENList[3*i+1];
    int x2 = mesh._ENList[3*i+2];
/*
    //get the coordinates of each vertex
    double coords_x0[2] {0.0, 0.0};
    double coords_x1[2] {0.0, 0.0};
    double coords_x2[2] {0.0, 0.0};

    mesh.GetCoords(x0, coords_x0);
    mesh.GetCoords(x1, coords_x1);
    mesh.GetCoords(x2, coords_x2);
*/
    //add edge to segments set
    segments.insert(_edges(x0, x1));
    segments.insert(_edges(x0, x2));
    segments.insert(_edges(x1, x2));
  } // end of iterate over alle triangles to get the segments

  return 0.0;
}
//end of GroupedPartitionsRefinement::CreateEncroachedSegmentsList()

//GroupedPartitionsRefinement::GetTriangleQuality()
double GroupedPartitionsRefinement::GetTriangleQuality()
{
  std::cout << "  GetTriangleQuality" << std::endl;
  return 0.0;
}
//end of GroupedPartitionsRefinement::GetTriangleQuality()

//GroupedPartitionsRefinement::CheckDelaunay
//uses global indices
bool GroupedPartitionsRefinement::CheckDelaunay()
{
  ofstream nelist_output;
  nelist_output.open("NEList.txt", std::ofstream::app);

  std::cout << "  CheckDelaunay" << std::endl;
  //iterate over all segments and check if they are delaunay
  int tmp_ctr = 0;
  for (auto it : segments)
  {
    std::cout << tmp_ctr << "/" << segments.size() << std::endl;
    //get the indices of the segment's vertices
    index_t x0 = std::get<0>(it);
    index_t x1 = std::get<1>(it);
    //end of get the indices of the segment's vertices

    //get the coordinates of each vertex    
    double coords_x0[2] {0.0, 0.0};
    double coords_x1[2] {0.0, 0.0};

    mesh.GetCoords(x0, coords_x0);
    mesh.GetCoords(x1, coords_x1);
    //end of get the coordinates of each vector

    //get NELists of vertices
    std::set<index_t> _NEList_x0;
    std::set<index_t> _NEList_x1;

    //get NEList of x0
    if (mesh.vertex_appearances[x0] == 1)
    {
      int partition {-1};
      int interface {-1};

      mesh.GetNEList(x0, _NEList_x0, &partition, &interface); //returns NEList with local indices

      //debug
      for (auto ne_it : _NEList_x0)
      {
        if (partition >= 0)
        { 
          //std::cout << x0 << " is in partitiong " << partition << " " << mesh.l2g_element_map_partitions[partition].at(ne_it) << std::endl;
          debug_NEList[x0].insert(mesh.l2g_element_map_partitions[partition].at(ne_it)); 
        }
        if (interface >= 0)
        { 
          debug_NEList[x0].insert(mesh.l2g_element_map_interfaces[interface].at(ne_it)); 
        }
      }
      //end of debug
/*
      if (partition >= 0)
      {      
        std::cout << "1 NEList of " << x0 << ": " << mesh.g2l_vertex_map_partitions[partition].at(x0) << " has " << _NEList_x0.size() << " elements" << std::endl;   
      }

      else if (interface >= 0)
      {
        std::cout << "1 NEList of " << x0 << ": " << mesh.g2l_vertex_map_interfaces[interface].at(x0) << " has " << _NEList_x0.size() << " elements" << std::endl;         
      }    
  
      for (auto it : _NEList_x0)
      {
        std::cout << "  " << it << std::endl;
      }
*/
    }
    //check if x0 is on any boundary
    else if (mesh.vertex_appearances[x0] == 2)
    {
      //check if it is on the boundary between one of the assigned partitions and the asigned interface
      if (!mesh.GetNEList_adv(x0, _NEList_x0))  //returns NEList with global indices, since 2 different pragmatic meshes are involved
      {
        //GetNEList_adv returns true if vertex is in one of the
        //assigned partitions and the interface or false otherwise
        //if false is returned, skip this segment!!!

        //however, try to get them for debugging purposes
        std::vector<bool> partitions(mesh.pragmatic_partitions.size());
        std::vector<bool> interfaces(mesh.pragmatic_interfaces.size());

        std::fill(partitions.begin(), partitions.end(), false);
        std::fill(interfaces.begin(), interfaces.end(), false);

        mesh.GetAllVertexPartitionsAndInterfaces(x0, partitions, interfaces);

        //debug: add all elements from all partitions and interfaces to debug NEList
        //loop over partitions
        for(size_t i = 0; i < partitions.size(); ++i)
        {
          if (partitions[i])
          {
            mesh.GetNEList_part(x0, _NEList_x0, i);

            for (auto ne_it : _NEList_x0)
            {
              debug_NEList[x0].insert(mesh.l2g_element_map_partitions[i].at(ne_it));
            }
          }
        } //end of loop over partitions

        //loop over interfaces
        for (size_t i = 0; i < interfaces.size(); ++i)
        {
          if (interfaces[i])
          {
            mesh.GetNEList_inter(x0, _NEList_x0, i);

            for (auto ne_it : _NEList_x0)
            {
              debug_NEList[x0].insert(mesh.l2g_element_map_interfaces[i].at(ne_it));
            }
          }
        } //end of loop over interfaces
        //end of debugging        
  
        continue;
      }
      //debug
      for (auto ne_it : _NEList_x0)
      {
        debug_NEList[x0].insert(ne_it);
      }
      //end of debug
/*              
      for (auto it : _NEList_x0)
      {
        std::cout << "  " << it << std::endl;
      }
*/
    }

    else 
    {
      std::vector<bool> partitions(mesh.pragmatic_partitions.size());
      std::vector<bool> interfaces(mesh.pragmatic_interfaces.size());

      std::fill(partitions.begin(), partitions.end(), false);
      std::fill(interfaces.begin(), interfaces.end(), false);

      mesh.GetAllVertexPartitionsAndInterfaces(x0, partitions, interfaces);

      //debug: add all elements from all partitions and interfaces to debug NEList
      //loop over partitions
      for(size_t i = 0; i < partitions.size(); ++i)
      {
        if (partitions[i])
        {
          mesh.GetNEList_part(x0, _NEList_x0, i);

          for (auto ne_it : _NEList_x0)
          {
            debug_NEList[x0].insert(mesh.l2g_element_map_partitions[i].at(ne_it));
          }
        }
      } //end of loop over partitions

      //loop over interfaces
      for (size_t i = 0; i < interfaces.size(); ++i)
      {
        if (interfaces[i])
        {
          mesh.GetNEList_inter(x0, _NEList_x0, i);

          for (auto ne_it : _NEList_x0)
          {
            debug_NEList[x0].insert(mesh.l2g_element_map_interfaces[i].at(ne_it));
          }
        }
      } //end of loop over interfaces
    } //end of else
    //end of get NEList of x0

    //get NEList of x1
    if (mesh.vertex_appearances[x1] == 1)
    {
      int partition {-1};
      int interface {-1};

      mesh.GetNEList(x1, _NEList_x1, &partition, &interface); //returns NEList with local indices

      //debug
      for (auto ne_it : _NEList_x1)
      {
        if (partition >= 0)
        { 
          debug_NEList[x1].insert(mesh.l2g_element_map_partitions[partition].at(ne_it)); 
        }
        if (interface >= 0)
        { 
          debug_NEList[x1].insert(mesh.l2g_element_map_interfaces[interface].at(ne_it)); 
        }
      }
      //end of debug
    }

    //check if x1 is on any boundary
    else if (mesh.vertex_appearances[x1] == 2)
    {
      //check if it is on the boundary between one of the assigned partitions and the asigned interface
      if (!mesh.GetNEList_adv(x1, _NEList_x1))  //returns NEList with global indices, since 2 different pragmatic meshes are involved
      {
        //GetNEList_adv returns true if vertex is in one of the
        //assigned partitions and the interface or false otherwise
        //if false is returned, skip this segment!!!

        //however try to get them for debugging purposes
        std::vector<bool> partitions(mesh.pragmatic_partitions.size());
        std::vector<bool> interfaces(mesh.pragmatic_interfaces.size());

        std::fill(partitions.begin(), partitions.end(), false);
        std::fill(interfaces.begin(), interfaces.end(), false);

        mesh.GetAllVertexPartitionsAndInterfaces(x1, partitions, interfaces);

        //loop over partitions
        for(size_t i = 0; i < partitions.size(); ++i)
        {
          if (partitions[i])
          {
            mesh.GetNEList_part(x1, _NEList_x1, i);

            for (auto ne_it : _NEList_x1)
            {
              debug_NEList[x1].insert(mesh.l2g_element_map_partitions[i].at(ne_it));
            }
          }
        } //for over partitions

        //loop over interfaces
        for (size_t i = 0; i < interfaces.size(); ++i)
        {       
          if (interfaces[i])
          {
            mesh.GetNEList_inter(x1, _NEList_x1, i);

            for (auto ne_it : _NEList_x1)
            {
              debug_NEList[x1].insert(mesh.l2g_element_map_interfaces[i].at(ne_it));
            }
          }
        } //end of for over interfaces
        //end of debugging
        continue;
      }

      //debug
      for (auto ne_it : _NEList_x1)
      {
        debug_NEList[x1].insert(ne_it);
      }
      //end of debug
    }
    //end of get NEList of x1

    //end of get NELists of vertices

/*
    //skip the edge if any of its vertices appears on more than 2 partitions or interfaces
    //TODO: check if every segment is locally delaunay, if yes the mesh is delaunay (according to the Delaunay Lemma)
    if (mesh.vertex_appearances[x0] >= 2 && mesh.vertex_appearances[x1] >= 2)
    {
      continue;
    }

    //if vertex appears in more than 1 partition or interface, check if it is on the boundary between partition and interface
    if (mesh.vertex_appearances[x0] == 2)
    {      
      //GetNEList_adv returns false, if vertex is in 2 different partitions
      //and true if it is in the assigned partition and interface
      if (!mesh.GetNEList_adv(x0, _NEList_x0)) 
      {
        continue;
      }
    }

    if (mesh.vertex_appearances[x1] == 2)
    {  
      //GetNEList_adv returns false, if vertex is in 2 different partitions
      //and true if it is in the assigned partition and interface
      if (!mesh.GetNEList_adv(x1, _NEList_x1)) 
      {
        continue;
      }
    }
/*
    if (mesh.vertex_appearances[x0] == 1 && mesh.vertex_appearances[x1] ==1)
    {
      std::cout << "Both vertices appear only once!" << std::endl;
    }
*/

    else 
    {
      std::vector<bool> partitions(mesh.pragmatic_partitions.size());
      std::vector<bool> interfaces(mesh.pragmatic_interfaces.size());

      std::fill(partitions.begin(), partitions.end(), false);
      std::fill(interfaces.begin(), interfaces.end(), false);

      mesh.GetAllVertexPartitionsAndInterfaces(x1, partitions, interfaces);

      //loop over partitions
      for(size_t i = 0; i < partitions.size(); ++i)
      {
        if (partitions[i])
        {
          mesh.GetNEList_part(x1, _NEList_x1, i);

          for (auto ne_it : _NEList_x1)
          {
            debug_NEList[x1].insert(mesh.l2g_element_map_partitions[i].at(ne_it));
          }
        }
      } //for over partitions

      //loop over interfaces
      for (size_t i = 0; i < interfaces.size(); ++i)
      {       
        if (interfaces[i])
        {
          mesh.GetNEList_inter(x1, _NEList_x1, i);

          for (auto ne_it : _NEList_x1)
          {
            debug_NEList[x1].insert(mesh.l2g_element_map_interfaces[i].at(ne_it));
          }
        }
      } //end of for over interfaces
    } //end of else
  ++tmp_ctr;
  } //end of iterate over all segments

  return true;
}
//end of GroupedPartitionsRefinement::CheckDelaunay

//GroupedPartitionsRefinement::PragmaticRefinement()
void GroupedPartitionsRefinement::PragmaticRefinement(double L_max)
{/*
  std::cout << "PragmaticRefinement(double L_Max)" << std::endl;

  //TODO: write refinement kernel similar to laplace kernel to do all these computations on the partitions and interfaces!!!
  size_t origNElements = mesh.num_elements;
  size_t origNNodes = mesh.num_nodes;
  size_t edgeSplitCnt = 0;
  splitCnt2= 0 ;

  //Reserve memory

  //taken from pragmatic Refine.h constructor
  newVertices.resize(nthreads);
  //newElements.resize(nthreads);
  //newBoundaries.resize(nthreads);
  //newQualities.resize(nthreads);
  newCoords.resize(nthreads);
  newMetric.resize(nthreads);

  // Pre-allocate the maximum size that might be required
  allNewVertices.resize(mesh._ENList.size());
  //end of constructor part

  size_t reserve_size = nedge*origNNodes/nthreads;
  std::cout << "reserve size " << reserve_size << " " << dim * reserve_size << std::endl;
  newVertices[0].clear();
  newVertices[0].reserve(reserve_size);
  newCoords[0].clear();
  newCoords[0].reserve(dim*reserve_size);
  newMetric[0].clear();
  newMetric[0].reserve(msize*reserve_size);

  std::cout << "origNNodes " << origNNodes << std::endl;

  //traverse all mesh edges
  //select them for refinement if edge length is greater than L_max (in transformed space (metric space))
  for (size_t i = 0; i < origNNodes; ++i)
  {
    std::cerr << i << std::endl;

    //skip node if it is in more than 2 partitions/interfaces
    //TODO: adapt to > 2 and adapt the GetNNList function for the case ==2
    if (mesh.vertex_appearances[i] >= 2)
    {
      continue;
    }

    std::vector<index_t> NNList_i;
    
    GetNNList_assigned(i, NNList_i);  //returns local indices

    for (size_t it = 0; it < NNList_i.size(); ++it)
    {
      index_t otherVertex = NNList_i[it];
      
      if (i < otherVertex)
      {
        double length = calc_edge_length(i, otherVertex);
        //std::cerr << i << " & " << otherVertex << ": " << length << std::endl;
        if (length > L_max)
        {
          //std::cerr << "node " << i << std::endl;
          ++splitCnt2;
          refine_edge(i, otherVertex, 0);
        }
      }
    } //end of loop over NNList_i.size()
  } //end of loop over origNNodes
  std::cerr << "looped over origNNodes " << origNNodes << std::endl;

  edgeSplitCnt = splitCnt2;

  threadIdx2 = mesh.mesh->get_number_nodes() + splitCnt2;

  //Resize meshes vectors
  mesh.mesh->resize_vectors(dim);

  //Append new coords and metric to the mesh
  mesh.AddCoords(threadIdx2, newCoords, dim, splitCnt2);
  mesh.AddMetric(threadIdx2, newMetric, msize, splitCnt2);

  //Fix IDs of new vertices
  for (size_t i = 0; i < splitCnt2; ++i)
  {
    newVertices[0][i].id = threadIdx2+i;
  }
  std::cerr << "fixed IDs" << std::endl;

  //Accumulate all new vertices in a contiguous array
  std::cerr<<'try memcpy' <<std::endl;
  memcpy(&allNewVertices[threadIdx2-origNNodes], &newVertices[0][0], newVertices[0].size()*sizeof(DirectedEdge<index_t>));
  std::cerr<<'memcpy successful' << std::endl;
  // Mark each element with its new vertices,
  // update NNList for all split edges.
  for (size_t i = 0; i < edgeSplitCnt; ++i)
  {
    index_t vid = allNewVertices[i].id;
    index_t firstid = allNewVertices[i].edge.first;
    index_t secondid = allNewVertices[i].edge.second;

    //TODO: Adapt this for the Kernel approach, to get the correct NELists for the assigned partitions and interface
    if(mesh.vertex_appearances[firstid] >= 2 || mesh.vertex_appearances[secondid] >= 2)
    {
      //std::cout << firstid << " " << secondid << std::endl;
      continue;
    }

    //Find which elements share this edge and mark them with their new vertices
    std::set<index_t> intersection;
    std::set<index_t> NEList1, NEList2;

    mesh.GetNEList_adv(firstid, NEList1);
    mesh.GetNEList_adv(secondid, NEList2);

    std::set_intersection( NEList1.begin(), NEList1.end(),
                           NEList2.begin(), NEList2.end(),
                          std::inserter(intersection, intersection.begin()) );

  } //end of for loop over edgeSplitCnt

  std::cout << origNElements << " " << origNNodes << " " << splitCnt2<< std::endl;*/
}
//end of GroupedPartitionsRefinement::PragmaticRefinement()

//GroupedPartitionsRefinement::RefinementKernel(int part1, int part2, int inter)
void GroupedPartitionsRefinement::RefinementKernel(int part1, int part2, int inter, double L_max)
{
  std::cout << "Refinement Kernel" << std::endl;
/*
  //Refinement for partition1 
  std::cout << "NNodes: " << mesh.pragmatic_partitions[part1]->get_number_nodes() << std::endl;
/*
  //for loop over part1
  for (auto it : mesh.g2l_vertex_map_partitions[part1])
  {
    //skip vertices which appear on boundary with an interface!!!
    if (mesh.vertex_appearances[it.first] >= 2)
    {
        continue;
    }
*//*
    size_t origNElements = mesh.pragmatic_partitions[part1]->get_number_elements();
    size_t origNNodes = mesh.pragmatic_partitions[part1]->get_number_nodes();
    size_t edgeSplitCnt = 0;

    new_vertices_per_element.resize(nedge*origNElements);
    std::fill(new_vertices_per_element.begin(), new_vertices_per_element.end(), -1);

    int tid = 0;
    splitCnt[tid] = 0;
    int ctr=0;
    /*
    * Average vertex degree in 2D is ~6, so there
    * are approx. (6/2)*NNodes edges in the mesh.
    * In 3D, average vertex degree is ~12.
    *//*
    size_t reserve_size = nedge*origNNodes/nthreads;
    newVertices[tid].clear();
    newVertices[tid].reserve(reserve_size);
    newCoords[tid].clear();
    newCoords[tid].reserve(dim*reserve_size);
    newMetric[tid].clear();
    newMetric[tid].reserve(msize*reserve_size);
    
    /* Loop through all edges and select them for refinement if
       its length is greater than L_max in transformed space. *//*
    for (size_t i = 0; i < origNNodes; ++i)
    {
      //get NNList for vertex i
      std::vector<index_t> NNList_i = mesh.pragmatic_partitions[part1]->get_reference_nnlist(i);

      for (size_t it = 0; it < NNList_i.size(); ++it)
      {
        index_t otherVertex = NNList_i[it];

        /* Conditional statement ensures that the edge length is only calculated once.
          * By ordering the vertices according to their gnn, we ensure that all processes
          * calculate the same edge length when they fall on the halo.
          *//*
        if(mesh.l2g_vertex_map_partitions[part1].at(i) < mesh.l2g_vertex_map_partitions[part1].at(otherVertex))
        {
          double length = calc_edge_length(mesh.l2g_vertex_map_partitions[part1].at(i), mesh.l2g_vertex_map_partitions[part1].at(otherVertex));
          std::cerr << i << " " << otherVertex << ": " << length<< std::endl;
          if (length > L_max)
          {
            ++ctr;
            std::cerr << "in " << i << " " << otherVertex << std::endl;
            ++splitCnt[tid];
            //refine_edge(i, otherVertex, tid);
          }

            std::cerr << "splitCnt " << splitCnt[tid] << std::endl;
        }
      }
    }//end of for loop over origNNodes

    //TODO: do omp atomic capture here!!!
    std::cerr << "OMP atomic capture missing!" << std::endl;
    std::cerr << origNNodes << " " << threadIdx[tid] << " " << mesh.pragmatic_partitions[part1]->get_number_nodes() << " " << splitCnt[tid] << " " << ctr << std::endl;
    threadIdx[tid] = mesh.pragmatic_partitions[part1]->get_number_nodes();
    mesh.pragmatic_partitions[part1]->set_nnodes(threadIdx[tid] + splitCnt[tid]); 
    std::cerr << threadIdx[tid] << " " << mesh.pragmatic_partitions[part1]->get_number_nodes() << std::endl;

  //} //end of for loop over part1
  //end of refinement for partition 1
*/
}
//end of GroupedPartitionsRefinement::RefinementKernel(int part1, int part2, int inter)

//GroupedPartitionsRefinement::GetNNList_assigned()
//function returns global indices!!!
void GroupedPartitionsRefinement::GetNNList_assigned(index_t i, std::vector<index_t>& NNList)
{
  //if vertex is on the interior of any of the assigned partitions or the interface
  //it is simple to process the NNList
  if (mesh.vertex_appearances[i] == 1)
  {
    int part {-1}, inter {-1};
    std::vector<index_t>* ptr_NNList = mesh.GetNNList(i, &part, &inter);
    
    if (part >= 0)
    {
      for (auto it : *ptr_NNList)
      {
        NNList.push_back(mesh.l2g_vertex_map_partitions[part].at(it));
      }
    }

    else
    {
      for (auto it : *ptr_NNList)
      {
        NNList.push_back(mesh.l2g_vertex_map_interfaces[inter].at(it));
      }
    }   
  }

  //if vertex is on any boundary, check if it is the one between assigned partition and interface
  else if (mesh.vertex_appearances[i] == 2)
  {
    int partition {-1};
    int interface {-1};

    mesh.GetVertexPartitionAndInterface(i, &partition, &interface);

    if (partition == -1 || interface == -1)
    {
      return;
    }
    
    //get NNLists with local indices
    std::vector<index_t>* NNList_partition = mesh.pragmatic_partitions[partition]->get_nnlist(mesh.g2l_vertex_map_partitions[partition].at(i));
    std::vector<index_t>* NNList_interface = mesh.pragmatic_interfaces[interface]->get_nnlist(mesh.g2l_vertex_map_interfaces[interface].at(i));

    //create NNLists with global indices
    std::vector<index_t> global_nnl_part( NNList_partition->size() );
    std::vector<index_t> global_nnl_inter( NNList_interface->size() );

    for (size_t j = 0; j < NNList_partition->size(); ++j)
    {
      global_nnl_part[j] = mesh.l2g_vertex_map_partitions[partition].at(NNList_partition->at(j));
    }

    for (size_t j = 0; j < NNList_interface->size(); ++j)
    {
      //global_nnl_inter[j] = mesh.l2g_vertex_map_interfaces[interface].at(NNList_interface->at(j));
    }
    
    //create union of two vectors with global indices and return the union vector
    std::set_union( global_nnl_part.begin(), global_nnl_part.end(),
                    global_nnl_inter.begin(), global_nnl_inter.end(),
                    inserter(NNList, NNList.begin()));
  }

  //if vertex appears more than 2 times, skip it!
  else
  {
    std::cerr << "no vertex appearance data available!" << std::endl;
  }
}
//end of GroupedPartitionsRefinement::GetNNList_assigned()

//GroupedPartitionsRefinement::GetNEList_assigned()
//function returns global indices!!!
void GroupedPartitionsRefinement::GetNEList_assigned(index_t i, std::set<index_t>& NEList)
{
  //if vertex is on the interior of any of the assigned partitions or the interface
  //it is simple to process the NEList
  if (mesh.vertex_appearances[i] == 1)
  {
    //TODO: adapt to get NEList from any partition or interface
    //std::cerr << "  NEList for " << i << " now only from partition 0! " << std::endl;
    int part {-1}, inter {-1};
    std::set<index_t> NEList_local_indices;
    mesh.GetNEList(i, NEList_local_indices, &part, &inter);
    //std::cerr << "now iterate over " << NEList_local_indices.size() << " elements from partition " << part << std::endl;

    for (auto it : NEList_local_indices)
    {
      //std::cerr << "    " << it << " " << mesh.l2g_element_map_partitions[part].at(it) << std::endl;
      NEList.insert(mesh.l2g_element_map_partitions[part].at(it));
    }
  }

  //if vertex is on any boundary, check if it is the one between assigned partition and interface
  else if (mesh.vertex_appearances[i] == 2)
  {
    //mesh.GetNEList_adv(i, NEList);
    //std::cerr << "NEList for " << i  << " from part0 and inter0" << std::endl;
    //this here works only for the test case with partition0
    int part {-1}, inter {-1};
    std::set<index_t> NEList_local_indices;
    mesh.GetNEList(i, NEList_local_indices, &part, &inter);
    for (auto it : NEList_local_indices)
    {
      NEList.insert(mesh.l2g_element_map_partitions[part].at(it));
    }
  }

  //if vertex appears more than 2 times, skip it!
  else
  {
    ;
  }
}
//end of GroupedPartitionsRefinement::GetNEList_assigned()

//GroupedPartitionsRefinement::RefinementKernel_2(int part1, int part2, int inter, double L_max)
void GroupedPartitionsRefinement::RefinementKernel_2(int part1, int part2, int inter, double L_max)
{
  RefinementKernel_2_part(part1, L_max);
  //RefinementKernel_2_part(part2, L_max);
  //RefinementKernel_2_part(part2, L_max);
}
//end of GroupedPartitionsRefinement::RefinementKernel_2(int part1, int part2, int inter, double L_max)

//GroupedPartitionsRefinement::RefinementKernel_2_part(int part, double L_max)
void GroupedPartitionsRefinement::RefinementKernel_2_part(int part, double L_max)
{

  size_t origNNodes = mesh.pragmatic_partitions[part]->get_number_nodes();
  size_t origNElements = mesh.pragmatic_partitions[part]->get_number_elements();
 
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
    std::vector<index_t> NNList_i;

    if (mesh.vertex_appearances[mesh.l2g_vertex_map_partitions[part].at(i)] == 2)
    {
      //Get NNList for vertex with local index i -> convert to global index
      GetNNList_assigned(mesh.l2g_vertex_map_partitions[part].at(i), NNList_i);
    }

    else 
    {
      //Get NNList for vertex with local index i -> convert to global index
      GetNNList_assigned(mesh.l2g_vertex_map_partitions[part].at(i), NNList_i);
    }

    for (auto otherVertex : NNList_i) //otherVertex is already a global index
    {
      //check if one of the vertices does not reside in partition 0, this is only for test case!
      std::unordered_map<index_t, index_t>::iterator position = mesh.g2l_vertex_map_partitions[part].find(otherVertex);

      if (position == mesh.g2l_vertex_map_partitions[part].end())
      {
        continue;
      }
   
      //order vertices according to their global index
      if( mesh.l2g_vertex_map_partitions[part].at(i) < otherVertex )
      {
        double length = calc_edge_length(mesh.l2g_vertex_map_partitions[part].at(i), otherVertex);

        if (length > L_max)
        {
          ++splitCnt;
          refine_edge_part(mesh.l2g_vertex_map_partitions[part].at(i), otherVertex, part);
        }
      }
    }
  }

  //TODO: do omp atomic capture here!!!
  threadIdx = mesh.pragmatic_partitions[part]->get_number_nodes();
  mesh.pragmatic_partitions[part]->set_nnodes(threadIdx + splitCnt); 

  //Resize meshes vectors
  mesh.mesh->resize_vectors(dim);

  //Append new coords, metric and appearance data to the mesh
  mesh.AddCoords_part(threadIdx, newCoords, dim, splitCnt, part);
  mesh.AddMetric_part(threadIdx, newMetric, msize, splitCnt, part);

  //std::cout << "lnnodes old " << threadIdx << " lnnodes new " << mesh.pragmatic_partitions[part]->get_number_nodes() << " gnnodes old " << mesh.num_nodes << std::endl;
  //Fix IDs of new vertices
  for (size_t j = 0; j < splitCnt; ++j)
  {
    //std::cerr << j << "/" << splitCnt << std::endl;
    newVertices[j].id = threadIdx+j;

    //std::cerr << "  l2g " <<  threadIdx+j << " " << mesh.num_nodes << std::endl;
    //std::cerr << "  g2l " << mesh.num_nodes << " " << threadIdx+j << std::endl;
    //update l2g and g2l vertex maps!!!
    mesh.l2g_vertex_map_partitions[part].insert( std::make_pair(threadIdx+j, mesh.num_nodes) );
    mesh.g2l_vertex_map_partitions[part].insert( std::make_pair(mesh.num_nodes, threadIdx+j) );
    //std::cerr << "updated mappings" << std::endl;
    //TODO atomic increment??
    ++mesh.num_nodes;

    //update vertex appearance data
    mesh.vertex_appearances.push_back(1); //TODO update for general case TODO
  }
  
  // Accumulate all newVertices in a contiguous array
  //TODO Do I need this memcpy here??? //TODO
  memcpy(&allNewVertices[threadIdx-origNNodes], &newVertices[0], newVertices.size() * sizeof(DirectedEdge<index_t>));

  // Mark each element with its new vertices,
  // update NNList for all split edges.

  for (size_t j = 0; j < splitCnt; j++)
  {
    index_t vid = allNewVertices[j].id; //in this test case a local index, but has to be changed to be a global index when using more partitions!!!
    index_t firstid = allNewVertices[j].edge.first; //in this test case a local index, but has to be changed to be a global index when using more partitions!!!
    index_t secondid = allNewVertices[j].edge.second; //in this test case a local index, but has to be changed to be a global index when using more partitions!!!

    //Find which elements share this edge and mark them with their new vertices
    //get NELists for firstid and secondid
    std::set<index_t> NEList_first;
    std::set<index_t> NEList_second;

    GetNEList_assigned(mesh.l2g_vertex_map_partitions[part].at(firstid), NEList_first);
    GetNEList_assigned(mesh.l2g_vertex_map_partitions[part].at(secondid), NEList_second); 

    std::set<index_t> intersection;
    std::set_intersection(NEList_first.begin(), NEList_first.end(),
                          NEList_second.begin(), NEList_second.end(),
                          std::inserter(intersection, intersection.begin()));
   
    for(auto element : intersection)
    {
      index_t eid = mesh.g2l_element_map_partitions[part].at(element); //now it is local (for test case)
      size_t edgeOffset = edgeNumber_part(eid, firstid, secondid, part);
      new_vertices_per_element[nedge * eid + edgeOffset] = vid;
    }
   
    //Update NNList for newly created vertices
    mesh.AddNNList_part(vid, firstid, part);
    mesh.AddNNList_part(vid, secondid, part);

    mesh.RemNNList_part(firstid, secondid, part);
    mesh.AddNNList_part(firstid, vid, part);
    mesh.RemNNList_part(secondid, firstid, part);
    mesh.AddNNList_part(secondid, vid, part);

  } //end of for loop over splitCnt

  //Start element refinement
  splitCnt = 0;
  newElements.clear();
  newBoundaries.clear();
  newQualities.clear();
  newElements.reserve(dim*dim*origNElements);
  newBoundaries.reserve(dim*dim*origNElements);
  newQualities.reserve(origNElements);

  threadIdx = mesh.pragmatic_partitions[part]->get_number_elements();

  //for loop over originalElements
  for (size_t eid = 0; eid < origNElements; ++eid)
  {
    //if element has been deleted, continue
    const index_t* n = mesh.pragmatic_partitions[part]->get_element(eid);

    if (n[0] < 0)
      continue;

    for (size_t j = 0; j < nedge; ++j)
    {
      if (new_vertices_per_element[nedge*eid+j] != -1)
      {
        refine_element_part(eid, part);
        break;
      }
    }
  }
  //end of for loop over origalElements

  //TODO: do omp atomic capture here!!!
  threadIdx = mesh.pragmatic_partitions[part]->get_number_elements();
  mesh.pragmatic_partitions[part]->set_nelements(threadIdx + splitCnt); 

  mesh.num_elements += splitCnt;//increase global node numbering, TODO: atomic increment!!!

  //Should we resize???
  if( mesh.pragmatic_partitions[part]->get_element_node_size() < mesh.pragmatic_partitions[part]->get_number_elements() * nloc) 
  {
    std::cerr << "RESIZE ENLIST, BOUNDARY AND QUALITY HAS TO BE IMPLEMENTED!!!" << std::endl;
  }

  //Append new elements to the mesh
  mesh.pragmatic_partitions[part]->add_elements(newElements, splitCnt, threadIdx);
  mesh.pragmatic_partitions[part]->add_boundaries(newBoundaries, splitCnt, threadIdx);
  mesh.pragmatic_partitions[part]->add_qualities(newQualities, splitCnt, threadIdx);

  //fix orientation of new elements
  size_t NElements = mesh.pragmatic_partitions[part]->get_number_elements();
  
  //for loop over NElements to fix their orientation
  //TODO: do i need this???
  for (size_t j = 0; j < NElements; ++j)
  {
    const index_t* n = mesh.pragmatic_partitions[part]->get_element(j);
    std::vector<index_t> _ENList2 = mesh.pragmatic_partitions[part]->get_element_node();
    //check if element is marked for deletion
    if (n[0] < 0)
    {
      continue;
    }
  }
  //end of for loop over NElements for fixing their orientation
}
//end of GroupedPartitionsRefinement::RefinementKernel_2_part(int part, double L_max)

//GroupedPartitionsRefinement::refine_element(size_t eid, int part)
//TODO: only for 2D implemented up to now!!!
void GroupedPartitionsRefinement::refine_element_part(size_t eid, int part)
{
  //works only for test case with partition0

  const int* n = mesh.pragmatic_partitions[part]->get_element(eid);
  
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
    (this->*refineMode2D[refine_cnt-1])(newVertex, eid, part);
  }
}
//end of GroupedPartitionsRefinement::refine_element(size_t eid, int part)

//edgeNumber_part
size_t GroupedPartitionsRefinement::edgeNumber_part(index_t eid, index_t v1, index_t v2, int part)
{
  const int *n=mesh.pragmatic_partitions[part]->get_element(eid);

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
//end of edgeNumber_part

//calc_edge_length(index_t x0, index_t y0)
//calculates edge length in metric space using global indices
double GroupedPartitionsRefinement::calc_edge_length(index_t x0, index_t y0)
{
  double length = -1.0;

  //TODO: adapt for 3D space!
  double m[3];
  double m_x[3], m_y[3];

  mesh.GetMetric(x0, m_x);
  mesh.GetMetric(y0, m_y);

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

  mesh.GetCoords(x0, x_coords);
  mesh.GetCoords(y0, y_coords);

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

//calc_edge_length(index_t x0, index_t y0, double* m)
//calculates edge length in metric space using global indices and metric array
double GroupedPartitionsRefinement::calc_edge_length(index_t x0, index_t y0, double* m)
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

  mesh.GetCoords(x0, x_coords);
  mesh.GetCoords(y0, y_coords);

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
//end of GroupedPartitionsRefinement::calc_edge_length(index x0, index y0, double* m)

//void GroupedPartitionsRefinement::refine_edge_part(index_t x, index_t y, int part)
void GroupedPartitionsRefinement::refine_edge_part(index_t n0, index_t n1, int part)
{  
  //swap indices because lesser id has to be the first
  if(n0 > n1)
  {
    index_t tmp = n0;
    n0 = n1;
    n1 = tmp;
  }

  newVertices.push_back( DirectedEdge<index_t> (mesh.g2l_vertex_map_partitions[part].at(n0), mesh.g2l_vertex_map_partitions[part].at(n1)) );  //n0 and n1 are global indices    
  
  //Calculate position of new point
  double x, m;
  double n0_coords[2], n1_coords[2];
  double m0[3], m1[3];
  
  mesh.GetCoords(n0, n0_coords);
  mesh.GetMetric(n0, m0);

  mesh.GetCoords(n1, n1_coords);
  mesh.GetMetric(n1, m1);

  double weight = 1.0 / (1.0 + sqrt( calc_edge_length(n0, n1, m0) / calc_edge_length(n0, n1, m1) ) );

  //calculate position of new vertex
  for (size_t i = 0; i < dim; ++i)
  {
    x = n0_coords[i] + weight * (n1_coords[i] - n0_coords[i]);
    newCoords.push_back(x);
    
    //std::cerr << "pushed" <<std::endl;
  }

  //Interpolate new metric
  for (size_t i = 0; i < msize; ++i)
  {
    m = m0[i] + weight * (m1[i] - m0[i]);
    newMetric.push_back(m);    
  }
}
//end of void GroupedPartitionsRefinement::refine_edge_part(index_t x, index_t y, int part)

//GroupedPartitionsRefinement::refine2D_1_part()
void GroupedPartitionsRefinement::refine2D_1_part(const index_t* newVertex, int eid, int part)
{
  //single edge split
  const int* n = mesh.pragmatic_partitions[part]->get_element(eid);  //works only in test case with partition0
  const int* boundary = mesh.pragmatic_partitions[part]->get_boundary_ptr(eid*nloc);

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
    mesh.AddNNList_part(vertexID, rotated_ele[0], part);
    // Add vertexID to rotated_ele[0]'s NNList
    mesh.AddNNList_part(rotated_ele[0], vertexID, part);

    // ele1ID is a new ID which isn't correct yet, it has to be
    // updated once each thread has calculated how many new elements
    // it created, so put ele1ID into addNE_fix instead of addNE.
    // Put ele1 in rotated_ele[0]'s NEList
    mesh.AddNEList_fix_part(rotated_ele[0], ele1ID, threadIdx, part);

    mesh.AddNEList_part(rotated_ele[0], eid, part);
    mesh.AddNEList_fix_part(vertexID, ele1ID, threadIdx, part);

    // Replace eid with ele1 in rotated_ele[2]'s NEList
    mesh.RemNEList_part(rotated_ele[2], eid, part);
    mesh.AddNEList_fix_part(rotated_ele[2], ele1ID, threadIdx,part);

    replace_element_part(eid, ele0, ele0_boundary, part);
    append_element_part(ele1, ele1_boundary, part);

    ++splitCnt;
  }
}
//end of GroupedPartitionsRefinement::refine2D_1()

//GroupedPartitionsRefinement::refine2D_2_part()
void GroupedPartitionsRefinement::refine2D_2_part(const index_t* newVertex, int eid, int part)
{
  const int* n = mesh.pragmatic_partitions[part]->get_element(eid);  //works only in test case with partition0
  const int* boundary = mesh.pragmatic_partitions[part]->get_boundary_ptr(eid*nloc);

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
  mesh.AddNNList_part(vertexID[0], vertexID[1], part);
  mesh.AddNNList_part(vertexID[1], vertexID[0], part);

  // vertexID[offset] and rotated_ele[offset+1] are the vertices on the diagonal
  mesh.AddNNList_part(vertexID[offset], rotated_ele[offset+1], part);
  mesh.AddNNList_part(rotated_ele[offset+1], vertexID[offset], part);

  // rotated_ele[offset+1] is the old vertex which is on the diagonal
  // Add ele2 in rotated_ele[offset+1]'s NEList
  mesh.AddNEList_fix_part(rotated_ele[offset+1], ele2ID, threadIdx, part);

  // Replace eid with ele0 in NEList[rotated_ele[0]]
  mesh.RemNEList_part(rotated_ele[0], eid, part);
  mesh.AddNEList_fix_part(rotated_ele[0], ele0ID, threadIdx, part);

  // Put ele0, ele1 and ele2 in vertexID[offset]'s NEList
  mesh.AddNEList_part(vertexID[offset], eid, part);
  mesh.AddNEList_fix_part(vertexID[offset], ele0ID, threadIdx, part);
  mesh.AddNEList_fix_part(vertexID[offset], ele2ID, threadIdx, part);

  // vertexID[(offset+1)%2] is the new vertex which is not on the diagonal
  // Put ele0 and ele2 in vertexID[(offset+1)%2]'s NEList
  mesh.AddNEList_fix_part(vertexID[(offset+1)%2], ele0ID, threadIdx, part);
  mesh.AddNEList_fix_part(vertexID[(offset+1)%2], ele2ID, threadIdx, part);

  replace_element_part(eid, ele1, ele1_boundary, part);
  append_element_part(ele0, ele0_boundary, part);
  append_element_part(ele2, ele2_boundary, part);
  splitCnt += 2;
}
//end of GroupedPartitionsRefinement::refine2D_2part()

//GroupedPartitionsRefinement::refine2D_3part()
void GroupedPartitionsRefinement::refine2D_3_part(const index_t* newVertex, int eid, int part)
{
  //TODO: change to work with any partition!!!
  const int *n=mesh.pragmatic_partitions[part]->get_element(eid);
  const int* boundary = mesh.pragmatic_partitions[part]->get_boundary_ptr(eid*nloc);

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
  mesh.AddNNList_part(newVertex[0], newVertex[1], part);
  mesh.AddNNList_part(newVertex[0], newVertex[2], part);
  mesh.AddNNList_part(newVertex[1], newVertex[0], part);
  mesh.AddNNList_part(newVertex[1], newVertex[2], part);
  mesh.AddNNList_part(newVertex[2], newVertex[0], part);
  mesh.AddNNList_part(newVertex[2], newVertex[1], part);

  //UpdateNEList  
  mesh.RemNEList_part(n[1], eid, part);
  mesh.AddNEList_fix_part(n[1], ele1ID, threadIdx, part);
  mesh.RemNEList_part(n[2], eid, part);
  mesh.AddNEList_fix_part(n[2], ele2ID, threadIdx, part);

  mesh.AddNEList_fix_part(newVertex[0], ele1ID, threadIdx, part);
  mesh.AddNEList_fix_part(newVertex[0], ele2ID, threadIdx, part);
  mesh.AddNEList_fix_part(newVertex[0], ele3ID, threadIdx, part);

  mesh.AddNEList_part(newVertex[1], eid, part);
  mesh.AddNEList_fix_part(newVertex[1], ele2ID, threadIdx, part);
  mesh.AddNEList_fix_part(newVertex[1], ele3ID, threadIdx, part);

  mesh.AddNEList_part(newVertex[2], eid, part);
  mesh.AddNEList_fix_part(newVertex[2], ele1ID, threadIdx, part);
  mesh.AddNEList_fix_part(newVertex[2], ele3ID, threadIdx, part);

  replace_element_part(eid, ele0, ele0_boundary, part);
  append_element_part(ele1, ele1_boundary, part);
  append_element_part(ele2, ele2_boundary, part);
  append_element_part(ele3, ele3_boundary, part);
}
//end of GroupedPartitionsRefinement::refine2D_3_part()

//GroupedPartitionsRefinement::replace_element_part(index_t eid, int* boundary)
void GroupedPartitionsRefinement::replace_element_part(const index_t eid, const index_t *n, const int* boundary, int part)
{
  for (size_t i = 0; i < nloc; ++i)
  {
    mesh.pragmatic_partitions[part]->set_element(eid, n[i], boundary[i], i);
  }

  mesh.update_quality(eid);
}
//end of GroupedPartitionsRefinement::replace_element_part()

//GroupedPartitionsRefinement::append_element(index_t eid, index_t *n, int* boundary)
void GroupedPartitionsRefinement::append_element_part(const index_t* elem, const int* boundary, int part)
{
  for (size_t i = 0; i < nloc; ++i)
  {
    newElements.push_back(elem[i]);
    newBoundaries.push_back(boundary[i]);
  }

  double q = mesh.calculate_quality(elem);
  
  newQualities.push_back(q);

  //update l2g and g2l element mappings
  size_t origNElements = mesh.pragmatic_partitions[part]->get_number_elements();
  size_t globalNElements = mesh.num_elements;

  mesh.l2g_element_map_partitions[part].insert( std::make_pair(origNElements+splitCnt, globalNElements+splitCnt) );
  mesh.g2l_element_map_partitions[part].insert( std::make_pair(globalNElements+splitCnt, origNElements+splitCnt) );

  ++splitCnt;
}
//end of GroupedPartitionsRefinement::append_element_part()

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif