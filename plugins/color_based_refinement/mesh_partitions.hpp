#ifndef MESH_PARTITIONS_HPP
#define MESH_PARTITIONS_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "Mesh.h"
#include "MetricField.h"
#include "Refine.h"
#include "ElementProperty.h"
#include "Edge.h"
#include "Swapping.h"

//TODO: DEBUG
#include "VTKTools.h"
//END OF DEBUG

//viennamesh includes
#include "viennameshpp/plugin.hpp"

//all other includes
#include "metis.h"
#include <unordered_map>
#include <map>
#include <numeric>  
#include <chrono>
#include <boost/container/flat_map.hpp>

#include "mesh_partitions_refinement.hpp"

#include "outbox.hpp"

extern "C"
{
  #include "triangle_interface.h"
}

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//class MeshPartitions
//
//This class partitions an input mesh and stores its partitions in different data structurs for further processing
class MeshPartitions
{
    public:
        MeshPartitions(Mesh<double> * const original_mesh, int num_regions,
         std::string filename, int thread);                                                   //Constructor
        ~MeshPartitions();                                                                    //Destructor

        //void convert_pra_to_tri(Mesh<double>* partition, struct triangulateio& tri_mesh);
        //void convert_tri_to_pra();

        //REPLACE THESE TWO VECTORS WITH A TEMPLATE VECTOR!!!
        std::vector<Mesh<double>*> pragmatic_partitions;                                      //Vector containing pointers to the pragmatic partitions
        std::vector<triangulateio> triangle_partitions;                                       //Vector containing the triangle data structurs

        bool MetisPartitioning();                                                             //Partition mesh using metis
        bool CreatePragmaticDataStructures_ser();                                             //Create Pragmatic Meshes storing the mesh partitions in serial
        bool CreatePragmaticDataStructures_par(std::vector<double>& threads_log, std::vector<double>& refine_times, std::vector<double>& l2g_build, 
                                               std::vector<double>& l2g_access, std::vector<double>& g2l_build, std::vector<double>& g2l_access,
                                               std::string algorithm, std::string options, std::vector<double>& triangulate_log,
                                               std::vector<double>& ref_detail_log, std::vector<double>& build_tri_ds);   //Create Pragmatic Meshes storing the mesh partitions in parallel
        bool CreateNeighborhoodInformation();                                                 //Create neighborhood information for vertices and partitions
        bool ColorPartitions();                                                               //Color the partitions
        bool WritePartitions();                                                               //ONLY FOR DEBUGGING!
        bool RefineInterior();                                                                //Refinement without refining boundary elements
        bool WriteMergedMesh(std::string filename);                                           //Merges partitions into a single mesh and writes it
        bool RefinementKernel(int part, double L_max);

        int get_colors(){return colors;};
        int get_max(){return max;};
        std::vector<std::vector<int>>& get_color_partitions() {return color_partitions;};
        void GetRefinementStats(int* nodes, int* elements, std::string algorithm);

        //mesh healing functions
        void heal_1(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int& splitCnt, int& threadIdx);
        void heal_2(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int& splitCnt, int& threadIdx);

    private:
/*
        bool MetisPartitioning(Mesh<double>* original_mesh, int num_regions);                 //Partition mesh using metis
        bool CreatePragmaticDataStructures_ser(Mesh<double>* original_mesh, int num_regions); //Create Pragmatic Meshes storing the mesh partitions in serial
        bool CreatePragmaticDataStructures_par(Mesh<double>* original_mesh, int num_regions); //Create Pragmatic Meshes storing the mesh partitions in parallel
        bool CreateNeighborhoodInformation(Mesh<double>* original_mesh, int num_regions);     //Create neighborhood information for vertices and partitions
        bool ColorPartitions(int num_regions);                                                               //Color the partitions
        bool WritePartitions();                                                               //ONLY FOR DEBUGGING!
*/

        template<typename _real_t, int _dim> friend class Refine;

        Mesh<double>* original_mesh;
        int num_regions;
        int nthreads;
        std::string file;

        //Variables for Metis   
        std::vector<idx_t> eptr;
        std::vector<idx_t> eind;
        idx_t num_nodes;
        idx_t num_elements;
        idx_t ncommon;
        //idx_t num_parts;
        idx_t result;
        std::vector<idx_t> epart;
        std::vector<idx_t> npart;

        //Variables used for Pragmatic data structures
        std::vector<std::set<index_t>> nodes_per_partition;
        std::vector<index_t> _ENList;
        std::vector<std::set<index_t>> elements_per_partition;

        //index mappings for the partitions
        std::vector<std::unordered_map<index_t, index_t>> g2l_vertex;
        std::vector<std::vector<int>> l2g_vertex;
        std::vector<std::unordered_map<index_t, index_t>> g2l_element;
        std::vector<std::vector<int>> l2g_element;

        //Neighborhood Information containers
        std::vector<std::set<int>> nodes_partition_ids;
        std::vector<std::set<int>> partition_adjcy;                                           //Stores the IDs of all neighboring partitions

        std::set<int>& get_nodes_partition_ids(int n){return nodes_partition_ids[n];};

        //Color information
        size_t colors;                                                                        //Stores the number of colors used
        std::vector<int> partition_colors;                                                    //Contains the color assigned to each partition
        std::vector<std::vector<int>> color_partitions;                                       //Contains the partition ids assigned to each color

        double calc_edge_length(int part_id, index_t x0, index_t y0);
        double calculate_quality(const index_t* n, int part_id);
        double update_quality(index_t element, int part_id);

        int edgeNumber(Mesh<double>*& partition, index_t eid, index_t v1, index_t v2);

        //Outboxes
        std::vector<Outbox> outboxes;

        //DEBUG
        int max=0;
        //END OF DEBUG

}; //end of class MeshPartitions

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Helper Functions                                                              //
//----------------------------------------------------------------------------------------------------------------------------------------------//

void convert_pra_to_tri(Mesh<double>* partition, triangulateio& tri_mesh)
{
    std::cout << "convert pragmatic to triangle data structure" << std::endl;

    //pointlist   
  /*  if (tri_mesh.pointlist) 
    {
        free(tri_mesh.pointlist);
    }*/
    tri_mesh.numberofpoints = partition->get_number_nodes();
    tri_mesh.pointlist = (REAL*)malloc( sizeof(REAL) * 2 * tri_mesh.numberofpoints);
    tri_mesh.numberofpointattributes = 0;

    std::cout << "point info done " << tri_mesh.numberofpoints << " " << partition->get_number_nodes() << " " << sizeof(tri_mesh.pointlist)  << std::endl;
/*
    if (out_mesh.pointlist) free(out_mesh.pointlist);
    out_mesh.pointlist = (REAL*)malloc( sizeof(REAL) * 2 * partition->get_number_nodes());
    out_mesh.numberofpoints = partition->get_number_nodes();
*/

    //std::vector<int> pragmatic_indices_to_triangle_indices(partition->get_number_nodes());

    int index = 0;
    for (size_t i = 0; i < partition->get_number_nodes(); ++i, ++index)
    {
        tri_mesh.pointlist[2*i] = partition->get_coords(i)[0];
        tri_mesh.pointlist[2*i+1] = partition->get_coords(i)[1];

        //pragmatic_indices_to_triangle_indices[i] = index;
    }

    std::cout << "points done " << std::endl;
    
    //trianglelist
    //if (tri_mesh.trianglelist) free(tri_mesh.trianglelist);
    tri_mesh.numberoftriangles = partition->get_number_elements();
    tri_mesh.trianglelist = (int*)malloc( sizeof(int) * 3 * tri_mesh.numberoftriangles );

    std::cout << "triangle info done " << tri_mesh.numberoftriangles << " " << partition->get_number_elements() << " " << sizeof(tri_mesh.trianglelist) << std::endl;

/*
    if (out_mesh.trianglelist) free(out_mesh.trianglelist);
    out_mesh.trianglelist = (int*)malloc( sizeof(int) * 3 * partition->get_number_elements() );
    out_mesh.numberoftriangles = partition->get_number_elements();
*/
    for (size_t i = 0; i < partition->get_number_elements(); ++i)
    {
        //std::cout << "  " << i << std::endl;
        const int *element_ptr = nullptr;
        element_ptr = partition->get_element(i);
       // std::cout << *(element_ptr++) << " " << *(element_ptr++) << " " << *(element_ptr++) << std::endl;
        
        tri_mesh.trianglelist[3*i] = *(element_ptr++);
        tri_mesh.trianglelist[3*i+1] = *(element_ptr++);
        tri_mesh.trianglelist[3*i+2] = *(element_ptr++);
    }

    std::cout << "triangles done" << std::endl;

    tri_mesh.numberofcorners = 3;
    tri_mesh.numberoftriangleattributes = 0;

    std::cout << "conversion done" << std::endl;

//end of create triangle mesh
} //end of convert_pra_to_tri*/
/*
void convert_tri_to_pra(triangulateio& tri_mesh, Mesh<double>* partition)
{
    std::cout << "convert triangle to pragmatic data structure" << std::endl;
} //end of convert_tri_to_pra*/

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Implementation                                                                //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//Constructor
//
//Tasks: TODO
//MeshPartitions::MeshPartitions(Mesh<double>* original_mesh, int num_regions, std::string filename)
MeshPartitions::MeshPartitions(Mesh<double> * const orig_mesh, int nregions, std::string filename, int threads)
{   /*
    ofstream memout;
    memout.open("memout.csv", ios::app);
    memout << filename << ", ";
    memout.close();
    
    //Create metric field
    MetricField<double,2> metric_field(*orig_mesh);

    double eta = 0.0001;
    std::vector<double> psi(orig_mesh->get_number_nodes());

    for (size_t i = 0; i < orig_mesh->get_number_nodes(); ++i)
    {
        //std::cerr << i << "/" << orig_mesh->get_number_nodes() << std::endl;
        
        double x = 2*orig_mesh->get_coords(i)[0];
        double y = 2*orig_mesh->get_coords(i)[1];
        psi[i] = 0.100000000000000*sin(50*x) + atan2(-0.100000000000000, (double)(2*x - sin(5*y)));
        /*
        double x = partition->get_coords(i)[0];
        double y = partition->get_coords(i)[1];
        psi[i] = x*y;
        
        psi[i]=10;
        //*/
/*    }
    std::cerr << "add_field" << std::endl;
    metric_field.add_field(&(psi[0]), eta, 1);
    std::cerr << "update_mesh" << std::endl;
    metric_field.update_mesh();
    std::cerr << "wow" << std::endl;
/*
    MetricField<double,2> metric_field(*orig_mesh);

    size_t NNodes = orig_mesh->get_number_nodes();
    for(size_t i=0; i<NNodes; i++) {
        //double psi = orig_mesh->get_coords(i)[0] * orig_mesh->get_coords(i)[1];
        double x = orig_mesh->get_coords(i)[0];
        double y = orig_mesh->get_coords(i)[1];
        //std::cerr << x << " " << y << " " << x*y << std::endl;
        double m[] = {0.0000000000001, 0.0, 0.0000000000001};
        metric_field.set_metric(m, i);
    }
    metric_field.update_mesh();
    */
    original_mesh = orig_mesh;
    num_regions = nregions;
    nthreads = threads;
    file = filename;

    /* 
    auto overall_tic = std::chrono::system_clock::now();

    auto wall_tic = std::chrono::system_clock::now();
        MetisPartitioning(original_mesh, num_regions);
    std::chrono::duration<double> partitioning_duration = std::chrono::system_clock::now() - wall_tic;
    //viennamesh::info(1) << "  Partitioning time " << wall_clock_duration.count() << std::endl;

    wall_tic = std::chrono::system_clock::now();
        CreateNeighborhoodInformation(original_mesh, num_regions);
    std::chrono::duration<double> adjacency_duration = std::chrono::system_clock::now() - wall_tic;
    //viennamesh::info(1) << "  Creating adjacency information time " << wall_clock_duration.count() << std::endl;

    wall_tic = std::chrono::system_clock::now();
        ColorPartitions(num_regions);
    std::chrono::duration<double> coloring_duration = std::chrono::system_clock::now() - wall_tic;
    //viennamesh::info(1) << "  Coloring time " << wall_clock_duration.count() << std::endl;

    std::chrono::duration<double> overall_duration = std::chrono::system_clock::now() - overall_tic;
    //viennamesh::info(1) << "  Overall time inside " << overall_duration.count() << std::endl;

    CreatePragmaticDataStructures_ser(original_mesh, num_regions); //Think about where I create the actual data structures!!!

    ofstream csv;
    csv.open("times.csv", ios::app);

    //csv << "File, Vertices, Elements, Partitions, Colors, Partitioning [s], Adjacency Info [s], Coloring [s], Total [s]" << std::endl;
    csv << filename << ", " << original_mesh->get_number_nodes() << ", " << original_mesh->get_number_elements() << ", ";
    csv << num_regions << ", " << max+1 << ", " << colors << ", ";
    csv <<  partitioning_duration.count() << ", ";
    csv <<  adjacency_duration.count() << ", ";
    csv <<  coloring_duration.count() << ", ";
    csv <<  overall_duration.count() << std::endl;
    csv.close();
  
    //WritePartitions();
*/
} //end of Constructor

//Destructor
//
//Tasks: TODO
MeshPartitions::~MeshPartitions()
{/*
    std::cout << "g2l mappings" << std::endl;
    std::cout << "found " << num_nodes << " vertices in the refined mesh" << std::endl;
    std::cout << "found " << num_elements << " elements in the refined mesh" << std::endl;

    std::vector<int> vertex_appearances(num_nodes);
    std::vector<int> element_appearances(num_elements);

    //WRITE LOOP THAT FILLS ABOVE VECTOR WITH THE PARTITION IDS OF THE PARTITIONS WHICH INCLUDE THE VERTEX
    for (size_t part = 0; part < pragmatic_partitions.size(); ++part)
    {
        //std::cout << "Partition " << part << " has color " << partition_colors[part] << std::endl;
        for (size_t local_vertex = 0; local_vertex < pragmatic_partitions[part]->get_number_nodes(); ++local_vertex)
        {
            std::cout << " " << local_vertex << " " << l2g_vertex[part].at(local_vertex) << std::endl;
            vertex_appearances[ l2g_vertex[part].at(local_vertex) ] = part;
        }

        for (size_t local_element = 0; local_element < pragmatic_partitions[part]->get_number_elements(); ++local_element)
        {
            element_appearances[ l2g_element[part].at(local_element) ] = part; 
            //std::cout << local_element << ": " << l2g_element[part].at(local_element) << std::endl;
            //element_appearances[ l2g_element[part].at(local_element) ].push_back(part);;
        }
    }

    //get coords
    std::vector<double> x_coords(num_nodes), y_coords(num_nodes);

    for (size_t vert = 0; vert < num_nodes; ++vert)
    {
        double p[2];
        //pragmatic_partitions[ vertex_appearances[vert].at(0) ]->get_coords( g2l_vertex[ vertex_appearances[vert].at(0) ], p);
        //std::cout << "V: " << vert << " in " << vertex_appearances[vert] << std::endl;
        //std::cout << "g2l: " << g2l_vertex[2].at(3) << std::endl;
        pragmatic_partitions[vertex_appearances[vert]]->get_coords( g2l_vertex[ vertex_appearances[vert] ].at(vert), p );

        x_coords[vert] = p[0];
        y_coords[vert] = p[1];
    }


    //get ENList
    //std::cout << "globalNElement: " << num_elements << std::endl;
    std::vector<int> ENList(num_elements*3);

    for (size_t ele = 0; ele < num_elements; ++ele)
    {
        //std::cout << "E: " << ele << " in " << element_appearances[ele] << " is " << g2l_element[element_appearances[ele]].at(ele) << std::endl;
        const index_t* element_ptr = nullptr;
        element_ptr = pragmatic_partitions[element_appearances[ele]]->get_element( g2l_element[ element_appearances[ele]].at(ele)  );
        
        ENList[3*ele+0] = l2g_vertex[element_appearances[ele]].at(*(element_ptr++));
        ENList[3*ele+1] = l2g_vertex[element_appearances[ele]].at(*(element_ptr++));
        ENList[3*ele+2] = l2g_vertex[element_appearances[ele]].at(*(element_ptr++));
    }

    std::cout << "create new mesh" << std::endl;
    Mesh<double>* merged_output = new Mesh<double> (num_nodes, num_elements, &(ENList[0]), &(x_coords[0]), &(y_coords[0]));
    merged_output->create_boundary();
    //new Mesh<double>(num_points_part, num_elements_part, &(ENList_part[0]), &(x_coords[0]), &(y_coords[0]));

    VTKTools<double>::export_vtu("adapted_mesh", merged_output);

    std::cout << "merged_output: " << merged_output->get_number_nodes() << " " << merged_output->get_number_elements() << std::endl;

    //delete merged_output;
*/
    //free used memory
    for (size_t i = 0; i < pragmatic_partitions.size(); ++i)
    {
        delete pragmatic_partitions[i];
    }

    //delete merged_output;
    std::cout << "freed memory" << std::endl;

    viennamesh::info(5) << "Called Destructor of MeshPartitions" << std::endl;

} //end of Destructor

//MetisPartitioning
//
//Tasks: Partitions the input mesh (in pragmatic data structure) into the specified number of partitions
//bool MeshPartitions::MetisPartitioning(Mesh<double>* const mesh, int num_regions)
bool MeshPartitions::MetisPartitioning()
{
    //get basic mesh information
    num_elements = original_mesh->get_number_elements();
    num_nodes = original_mesh->get_number_nodes();
    ncommon = original_mesh->get_number_dimensions();
   // std::vector<idx_t> bdry = original_mesh->copy_boundary_vector();
   // size_t num_bdry_nodes = std::accumulate(bdry.begin(), bdry.end(), 0);
   // idx_t numflag = 0;     //0...C-style numbering is assumed that starts from 0; 1...Fortran-style numbering is assumed that starts from 1

    //reserve memory
    epart.reserve(num_elements);
    npart.reserve(num_nodes);

 /*   std::vector<idx_t> xadj;
    xadj.resize(num_nodes+1);
    std::vector<idx_t> adjncy;
    adjncy.resize(2*(3*num_nodes - 3 - num_bdry_nodes)); //see: http://math.stackexchange.com/questions/1541125/total-number-of-edges-in-a-triangle-mesh-with-n-vertices
*//*
    int num_edges = 2*(3*num_nodes - 3 - num_bdry_nodes);
    idx_t *xadj;
    xadj = new idx_t[num_nodes+1];
    idx_t *adjncy;
    adjncy= new idx_t[num_edges];

 //   std::cout << adjncy.size() << std::endl;
*/
    //fill eptr and eind, as done in viennamesh plugin "metis", file "mesh_partitionig.cpp"  
    eptr.push_back(0);
    
    for (size_t i = 0; i < original_mesh->get_number_elements(); ++i)
    {
        const index_t* element_ptr = nullptr;
        element_ptr = original_mesh->get_element(i);

        for (size_t j = 0; j < (original_mesh->get_number_dimensions() + 1); ++j)
        {
            eind.push_back( *(element_ptr+j) );
        }               

        eptr.push_back( eind.size() );    
    }  

   /* //DEBUG
    ofstream outfile;
    outfile.open("box300x300.metis");

    std::vector<int> _ENList_dbg = mesh->get_enlist();
    
    outfile << mesh->get_number_elements() << std::endl;

    for (size_t i = 0; i < mesh->get_number_elements(); ++i)
    {
      outfile << _ENList_dbg[3*i]+1 << " " << _ENList_dbg[3*i+1]+1 << " " << _ENList_dbg[3*i+2]+1 << std::endl;
    }
    outfile.close();
    //END OF DEBUG*/
    //Call Metis Partitioning Function (see metis manual for details on the parameters and on the use of the metis API)
    /*METIS_PartMeshDual(&num_elements,
                         &num_nodes,
                         eptr.data(),
                         eind.data(),
                         NULL,
                         NULL,
                         &ncommon,
                         &num_regions,
                         NULL,
                         NULL,
                         &result,
                         epart.data(),
                         npart.data());
                         //*/

   METIS_PartMeshNodal (&num_elements,
                        &num_nodes,
                        eptr.data(),
                        eind.data(),
                        NULL,
                        NULL,
                        &num_regions,
                        NULL,
                        NULL,
                        &result,
                        epart.data(),
                        npart.data());

    /*
    METIS_MeshToDual(&num_elements,
                     &num_nodes,
                     eptr.data(),
                     eind.data(),
                     &ncommon,
                     &numflag,
                     &xadj,
                     &adjncy);
    */

    viennamesh::info(5) << "Created " << num_regions << " mesh partitions using Metis" << std::endl;

    /*//DEBUG
    ofstream epart_stream;
    epart_stream.open("epart.8");
  
    for (size_t i = 0; i < num_elements; ++i)
    {
      epart_stream << epart[i] << std::endl;
    }
    //epart_stream.close();
    //END OF DEBUG*/

/*
    std::cout << "xadj: " << std::endl;
    for (size_t i = 0; i < num_nodes+1; ++i)
        std::cout << " " << i << ": " << xadj[i] << std::endl;
/*
    std::cout << "adjncy: " << std::endl;
    for (size_t i = 0; i < num_edges; ++i)
        std::cout << " " << i << ": " << adjncy[i] << std::endl;
*/

    return true;
}//end of MetisPartitioning

//CreateNeighborhoodInformation
//
//Tasks: Populate Vertex partition container and create adjacency lists for each partition
//bool MeshPartitions::CreateNeighborhoodInformation(Mesh<double>* original_mesh, int num_regions)
bool MeshPartitions::CreateNeighborhoodInformation()
{  
  //prepare a partition id container for vertices and partitions
  //this will be populated with all partitions a vertex is part of
  nodes_partition_ids.resize(num_nodes);
  partition_adjcy.resize(num_regions);

  //populate the container
  for(size_t i = 0; i < num_elements; ++i)
  {
    //get vertices of element
    const index_t *element_ptr = nullptr;
    element_ptr = original_mesh->get_element(i);

    size_t ndims = original_mesh->get_number_dimensions();

    //iterate element vertices and add partition id
    for (size_t j = 0; j < (ndims+1); ++j)
    {
      nodes_partition_ids[*(element_ptr++)].insert(epart[i]);
    }

    //DEBUG
    if (epart[i] > max)
      max = epart[i];
    //END OF DEBUG

  }

  //create partition adjacency information
  for (size_t i = 0; i < nodes_partition_ids.size(); ++i)
  {
      if (nodes_partition_ids[i].size() > 1)
      {
        for (auto set_iter : nodes_partition_ids[i])
        {
            for (auto set_iter2 : nodes_partition_ids[i])
            {
                if (set_iter == set_iter2)
                    continue;

                partition_adjcy[set_iter].insert(set_iter2);
            }
        }
      }
  }

  
  //DEBUG
  for (size_t i = 0; i < num_regions; ++i)
  {
      std::cout << "Partition " << i << " has the following neighbors: " << std::endl;

      for (auto iter : partition_adjcy[i])
      {
          std::cout << "  " << iter << std::endl;
      }
  }
  //END OF DEBUG*/  

  return true;
}
//end of CreateNeighborhoodInformation

//ColorPartitions
//
//Tasks: Color the partitions such that independent sets are created
//bool MeshPartitions::ColorPartitions(int num_regions)
bool MeshPartitions::ColorPartitions()
{
    viennamesh::info(1) << "Coloring partitions" << std::endl;
    //resize vector
    partition_colors.resize(partition_adjcy.size());

    colors = 1;                 //number of used colors
    partition_colors[0] = 0;    //assign first partition color 0

    //visit every partition and assign the smallest color available (not already assigned to on of its neighbors)
    for (size_t i = 1; i < partition_colors.size(); ++i)
    {
        
        //int tmp_color = partition_colors[*(partition_adjcy[i].begin())] + 1;   //assign next color
        int tmp_color = 0; //start with smallest color 
        bool next_color = false;

        do
        {
            //check if assigned color in tmp_color is already assigned to a neighbor
            //since we assign colors to partitions in ascending ID order, check only
            //neighbors with smaller partition ID
            for (auto iter : partition_adjcy[i])
            {
                //if chosen color is already assigned to neighbor, try next color
                if ( i > iter && partition_colors[iter] == tmp_color) 
                {
                    ++tmp_color;
                    next_color = true;
                    break;
                }

                //if chosen color is ok exit loop
                else
                    next_color=false;
            }
        } while(next_color);

       /* for (size_t j = 1; j < partition_adjcy[i].size(); ++j)
        {
            //check if assigned color in tmp_color is already assigned to a neighbor
            if (partition_colors[partition_adjcy[i][j]] >= tmp_color)
            {
                tmp_color = partition_colors[partition_adjcy[i][j]] + 1;
            }
        }*/

        partition_colors[i] = tmp_color;

        if ( (tmp_color + 1) > colors )
        {
            colors = tmp_color + 1;
        }
    }

    //create a vector containing the color information for each partition
    //each vector element is one color and contains the partitions with this color
    color_partitions.resize(colors);

    for (size_t i = 0; i < partition_colors.size(); ++i)
    {
        color_partitions[ partition_colors[i] ].push_back(i);
    }

    //DEBUG
    //std::cout << "Number of used colors: " << colors << std::endl;
    std::cout << "  Partition | Color " << std::endl;
 
    for (size_t i = 0; i < partition_colors.size(); ++i)
    {
        std::cout << "          " << i << " | " << partition_colors[i] << std::endl;
    }
    //*/
/*
    std::cout << std::endl << "      Color | #Partitions " << std::endl;
 
    for (size_t i = 0; i < color_partitions.size(); ++i)
    {
        std::cout << "          " << i << " | " << color_partitions[i].size() << std::endl;
    
   /*     std::cout << "          " << i << " | ";
        for (auto it : color_partitions[i])
        {
           std::cout << it << " ";
        }

        std::cout << std::endl;
    }
    //END OF DEBUG*/

    viennamesh::info(1) << "   Partitions param = " << num_regions << std::endl;
    viennamesh::info(1) << "   Partitions count = " << max+1 << std::endl;
    viennamesh::info(1) << "   Number of colors = " << colors << std::endl;

    return true;
}
//end of ColorPartitions

//CreatePragmaticDataStructures_ser
//
//Tasks: Get and order data needed to create a pragmatic data structure for each partition
//Runs only serial
//bool MeshPartitions::CreatePragmaticDataStructures_ser(Mesh<double>* const original_mesh, int num_regions)
bool MeshPartitions::CreatePragmaticDataStructures_ser()
{
    //reserve memory
    nodes_per_partition.resize(num_regions);
    elements_per_partition.resize(num_regions);
    g2l_vertex.resize(num_regions);
    l2g_vertex.resize(num_regions);

    //get ENList
    _ENList = original_mesh->get_enlist();

    //get the nodes for each partition
    for (int i = 0; i < num_elements; ++i)
    {
        //add nodes
        //TODO: the following hast to be updated for the 3D case!!
        nodes_per_partition[ epart[i] ].insert(_ENList[i*3]);
        nodes_per_partition[ epart[i] ].insert(_ENList[i*3+1]);
        nodes_per_partition[ epart[i] ].insert(_ENList[i*3+2]);

        //add element
        elements_per_partition[ epart[i] ].insert(i);
    } 
    //end of get nodes per partition

    //Create Partitions

    //vectors storing the coordinate
    std::vector< std::vector<double>> x_coords(num_regions);
    std::vector< std::vector<double>> y_coords(num_regions);
    std::vector< std::vector<double>> z_coords(num_regions);

    //vector storing the ENLists of each region
    std::vector<std::vector<index_t>> ENLists_partitions(num_regions);

     //loop over all partitions
    for (size_t i = 0; i < num_regions; ++i)
    {
        //get number of vertices and elements
        int num_points = nodes_per_partition[i].size();
        int num_cells = elements_per_partition[i].size();

        //get the vertex-to-index-mapping between old and new indices
        //and additionally the index-to-vertex-mapping
        //TODO: use unordered_map instead, to speed up the code

        std::unordered_map <index_t, index_t> global_to_local_index_map;
        std::unordered_map <index_t, index_t> local_to_global_index_map;

        index_t new_vertex_id = 0;
        for (auto it : nodes_per_partition[i])
        {
            global_to_local_index_map.insert( std::make_pair(it, new_vertex_id++) );
            //++vertex_appearances[it];
        }

        //global_to_local_index_mappings_partitions[i] = global_to_local_index_map;

        //and get also the index-to-vertex mapping (opposite direction than vertex to index mapping)
        for (auto it : global_to_local_index_map)
        {
            local_to_global_index_map[it.second] = it.first;
        }

        //local_to_global_index_mappings_partitions[i] = local_to_global_index_map;

        //pre-allocate memory
        x_coords[i].reserve(num_points);
        y_coords[i].reserve(num_points);
        ENLists_partitions[i].resize(3*num_cells);

        //get coordinates of each vertex
        int counter = 0;
        for (auto it : nodes_per_partition[i])
        {
            double p[2];
            original_mesh->get_coords( it, p);
            x_coords[i][counter] = p[0];
            y_coords[i][counter] = p[1];
            ++counter;
        }

        //create ENList with respect to the new vertex indices
        counter=0;        
        for (auto it : elements_per_partition[i])
        {      
            const index_t *element_ptr = nullptr;
            element_ptr = original_mesh->get_element(it);
            
            ENLists_partitions[i][counter++] = global_to_local_index_map[*(element_ptr++)];
            ENLists_partitions[i][counter++] = global_to_local_index_map[*(element_ptr++)];
            ENLists_partitions[i][counter++] = global_to_local_index_map[*(element_ptr++)];  
            ///++element_appearances[it];       
        }

        //create pragmatic mesh 
        Mesh<double> *partition_mesh = nullptr;

        //TODO: change for 3D refinement
        //mesh = new Mesh<double> ( num_points, num_cells, &(ENLists_regions[region.id()][0]) ,&(x_coords[region.id()][0]), &(y_coords[region.id()][0]), &(z_coords[region.id()][0]) );        
        partition_mesh = new Mesh<double> ( num_points, num_cells, &(ENLists_partitions[i][0]), &(x_coords[i][0]), &(y_coords[i][0]) );
        partition_mesh->create_boundary();
        pragmatic_partitions.push_back(partition_mesh);

        //Create metric field
        MetricField<double,2> metric_field(*partition_mesh);

        double eta = 0.0001;
        std::vector<double> psi(num_points);

        for (size_t i = 0; i < num_points; ++i)
        {
            //double x = 2*partition->get_coords(i)[0]-1;
            //double y = 2*partition->get_coords(i)[1]-1;
            double x = partition_mesh->get_coords(i)[0];
            double y = partition_mesh->get_coords(i)[1];
    
            psi[i] = 0.100000000000000*sin(50*x) + atan2(-0.100000000000000, (double)(2*x - sin(5*y)));
            //psi[i] = 800000000000;
        }
        
        metric_field.add_field(&(psi[0]), eta, 1);
        metric_field.update_mesh();

        //Create boundary information for refinement algorithm 
        std::vector<int> bdry = partition_mesh->get_boundaryRef();
        std::vector<int> boundary_nodes(partition_mesh->get_number_nodes(), 0);

        for (size_t eid = 0; eid < partition_mesh->get_number_elements(); ++eid)
        {
            const int *n = partition_mesh->get_element(eid);
            const int *boundary=&(bdry[eid*3]);

            //-1 means element is marked for deletion
            if(n[0]==-1)
                continue;

            if( partition_mesh->get_number_dimensions() == 2 ) 
            {
                for(int j=0; j<3; j++) 
                {
                    boundary_nodes[n[(j+1)%3]] = std::max(boundary_nodes[n[(j+1)%3]], bdry[eid*3+j]);
                    boundary_nodes[n[(j+2)%3]] = std::max(boundary_nodes[n[(j+2)%3]], bdry[eid*3+j]);
                }
            }
        }

        Refine<double,2> refiner(*partition_mesh);
        //refiner.refine(0.5, nodes_partition_ids, local_to_global_index_map);

        //delete partial_mesh;        //TODO: creates segfault if comments are removed
    } 
    //end of loop over all partitions

    viennamesh::info(5) << "Created " << pragmatic_partitions.size() << " pragmatic mesh data structures" << std::endl;
/*
    //DEBUG
      //write partitions
      for (size_t i = 0; i < pragmatic_partitions.size(); ++i)
      {
        std::cout << "Writing partition " << i << std::endl;
        std::cout << "  Vertex count = " << pragmatic_partitions[i]->get_number_nodes() << std::endl;
        std::cout << "  Cell count = " << pragmatic_partitions[i]->get_number_elements() << std::endl;
        
        std::string filename;
        filename += "examples/data/color_refinement/";
        filename += "output_";
        filename += "_partition_";
        filename += std::to_string( i );
      
        VTKTools<double>::export_vtu(filename.c_str(), pragmatic_partitions[i]);
      }
    //END OF DEBUG
*/
    return true;
}
//end of CreatePragmaticDataStructures_ser

//CreatePragmaticDataStructures_par
//
//Tasks: Get and order data needed to create a pragmatic data structure for each partition
//Runs in parallel
//bool MeshPartitions::CreatePragmaticDataStructures_par(Mesh<double>* const original_mesh, int num_regions)
bool MeshPartitions::CreatePragmaticDataStructures_par(std::vector<double>& threads_log, std::vector<double>& refine_times, std::vector<double>& l2g_build, 
                                                       std::vector<double>& l2g_access, std::vector<double>& g2l_build, std::vector<double>& g2l_access,
                                                       std::string algorithm, std::string options, std::vector<double>& triangulate_log,
                                                       std::vector<double>& ref_detail_log, std::vector<double>& build_tri_ds)
{    
    viennamesh::info(1) << "Starting mesh adaptation" << std::endl;

    threads_log.resize(nthreads);
    refine_times.resize(nthreads);
    triangulate_log.resize(nthreads);
    ref_detail_log.resize(4, 0.0);
    build_tri_ds.resize(nthreads);
    outboxes.resize(num_regions, Outbox(false));
    
    std::fill(threads_log.begin(), threads_log.end(), 0.0);
    std::fill(refine_times.begin(), refine_times.end(), 0.0);
    std::fill(triangulate_log.begin(), triangulate_log.end(), 0.0);
    std::fill(build_tri_ds.begin(), build_tri_ds.end(), 0.0);
/*
    //debug 
    l2g_build.resize(nthreads);
    l2g_access.resize(nthreads);
    g2l_build.resize(nthreads);
    g2l_access.resize(nthreads);

    std::fill(l2g_build.begin(), l2g_build.end(), 0.0);
    std::fill(l2g_access.begin(), l2g_access.end(), 0.0);
    std::fill(g2l_build.begin(), g2l_build.end(), 0.0);
    std::fill(g2l_access.begin(), g2l_access.end(), 0.0);
    //end of debug
*/
  //  auto prep_tic = std::chrono::system_clock::now();

    std::vector<int> _ENList_orig = original_mesh->get_enlist();
    std::vector<std::set<int>> elements_part(num_regions);
    int dim = original_mesh->get_number_dimensions();

    nodes_per_partition.resize(num_regions);


    //REPLACE THESE TWO WITH TEMPLATE COMMAND
    pragmatic_partitions.resize(num_regions);
    triangle_partitions.resize(num_regions);

  /*  //vectors storing the index mapping information for all partitions
    std::vector<std::unordered_map<int,int>> g2l_vertices(num_regions);
    std::vector<std::vector<int>> l2g_vertices(num_regions);*/

    l2g_vertex.resize(num_regions);
    g2l_vertex.resize(num_regions);
    l2g_element.resize(num_regions);
    g2l_element.resize(num_regions);

    //get vertices and elements from original mesh
    auto nodes_part_tic = std::chrono::system_clock::now();

    for (size_t ele_id = 0; ele_id < original_mesh->get_number_elements(); ++ele_id)
    {
        if (dim == 2)
        {
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id] );
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id+1] );
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id+2] );
        }

        else
        {
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id] );
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id+1] );
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id+2] ); 
            nodes_per_partition[ epart[ele_id] ].insert( _ENList_orig[3*ele_id+3] ); 
        }
        elements_part[ epart[ele_id] ].insert(ele_id);
    }
/*
    std::chrono::duration<double> nodes_part_time = std::chrono::system_clock::now() - nodes_part_tic;

    std::chrono::duration<double> prep_time = std::chrono::system_clock::now() - prep_tic;
/*
    times[0] += prep_time.count();
    times[1] += nodes_part_time.count();
*/

    //iterate colors
    for (size_t color = 0; color < colors; color++)
    {
        std::cout << std::endl << "actual color / # of colors" << std::endl;
        std::cout << color << " / " << colors << std::endl;

        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t part_iter = 0; part_iter < color_partitions[color].size(); ++part_iter)
        {
            auto omp_tic = omp_get_wtime();
   
            size_t part_id = color_partitions[color][part_iter];
            std::cerr << " working on partition " << part_id << std::endl;

            Outbox outbox_data;

            auto msize {0};
            auto nedge {0};  
            auto nloc {0};

            if (dim == 2)
            {
                msize = 3;
                nedge = 3;
                nloc = dim+1;
            }

            else
            {
                msize = 6;
                nedge = 6;
                nloc = dim+1;
            }

            std::vector<int> new_vertices_per_element;

            //get number of vertices and elements of the local partitions
            //int num_points_part = nodes_per_partition[part_id].size();
            int num_points_part = nodes_per_partition[part_id].size() + outbox_data.num_verts();
            int num_elements_part = elements_part[part_id].size();
            
            //create coordinate vectors, g2l- and l2g-index-mappings for the vertices
            //std::unordered_map<int, int> g2l_tmp, l2g_tmp;
            std::unordered_map<int, int> g2l_vertices_tmp(num_points_part);
            std::unordered_map<int, int> g2l_elements_tmp(num_elements_part);
            //std::map<int,int> l2g_tmp, g2l_tmp;
            //std::map<int,int> g2l_tmp;
            //std::unordered_map<int,int> l2g_tmp;

            //boost::container::flat_map<int,int, std::less<int>> g2l_tmp;

            //test: use vector instead of unordered maps!
            std::vector<int> l2g_vertices_tmp(num_points_part);
            std::vector<int> l2g_elements_tmp(num_elements_part);

            int new_vertex_id = 0;    //new vertex id for local partition
            
          //  auto coords_tic = std::chrono::system_clock::now();

            //coords vector
            std::vector<double> x_coords(num_points_part);
            std::vector<double> y_coords(num_points_part); 
            std::vector<double> z_coords;

            if (dim == 3)
                z_coords.reserve(num_points_part);

/*
            double l2g_time {0.0};
            double g2l_time {0.0};
*/
            for (auto it : nodes_per_partition[part_id])
            {
                if (dim == 2)
                {
                    double p[2];
                    original_mesh->get_coords( it, p);
                    
                    x_coords[new_vertex_id] = p[0];
                    y_coords[new_vertex_id] = p[1];
                }

                else
                {
                    double p[3];
                    original_mesh->get_coords( it, p);
                    
                    x_coords[new_vertex_id] = p[0];
                    y_coords[new_vertex_id] = p[1]; 
                    z_coords[new_vertex_id] = p[2];
                }
                //TODO: Update for 3d case!

   //             auto g2l_tic = std::chrono::system_clock::now();
                g2l_vertices_tmp.insert( std::make_pair(it, new_vertex_id) );
   /*             std::chrono::duration<double> g2l_dur = std::chrono::system_clock::now() - g2l_tic;
                g2l_time += g2l_dur.count();             
*/
  //              auto l2g_tic = std::chrono::system_clock::now();
                l2g_vertices_tmp[new_vertex_id++] = it;
    /*            std::chrono::duration<double> l2g_dur = std::chrono::system_clock::now() - l2g_tic;
                l2g_time += l2g_dur.count();
        */        
            } //end of for over nodes_per_partition[part_id]
/*
            if (color > 0)
            {
                for (size_t it = 0; it < outbox_data.num_verts(); ++it)
                {
                    double p[2];
                    pragmatic_partitions[0]->get_coords( it, p);
                    
                    x_coords[new_vertex_id] = p[0];
                    y_coords[new_vertex_id] = p[1];
                }

                //DEBUG
                for (size_t i = 0; i < outbox_data.num_verts(); ++i)
                {
                    std::cout << "Vertex " << outbox_data[3*i+2] << " has been inserted between " << outbox_data[3*i] << " and " << outbox_data[3*i+1] << std::endl;
                }
                //END OF DEBUG
    
                std::cout << "WRITE MESH HEALING ALGORITHM!!!" << " " << x_coords.size() << std::endl;

                std::cout << "FIRST UPDATE NNLISTS!" << std::endl;
                for (size_t i = 0; i < outbox_data.num_verts(); ++i)
                {
                    auto v1 = outbox_data[3*i];
                    auto v2 = outbox_data[3*i+1];
                    auto new_vert = outbox_data[3*i+2];

                    //Find local element ids in actual partition
                    std::cout << v1 << " is " << l2g_vertex[0][v1] << " and " << v2 << " is " << l2g_vertex[0][v2] << std::endl;

                    auto glob_v1 = l2g_vertex[0][v1];
                    auto glob_v2 = l2g_vertex[0][v2];

                    std::cout << glob_v1 << " is here " << g2l_vertices_tmp.at(glob_v1) << " and " << glob_v2 << " is here " << g2l_vertices_tmp.at(glob_v2) << std::endl;

                    /*  // Find which elements share this edge and mark them with their new vertices.
                    std::set<index_t> intersection;
                    std::set<int> NEList_v1 = pragmatic_partitions[part_id]->get_nelist( g2l_vertices_tmp.at(glob_v1) );
                    std::set<int> NEList_v2 = pragmatic_partitions[part_id]->get_nelist( g2l_vertices_tmp.at(glob_v2) );

                    std::set_intersection(NEList_v1.begin(), NEList_v1.end(), NEList_v2.begin(), NEList_v2.end(), std::inserter(intersection, intersection.begin()));
                    /*std::set_intersection(pragmatic_partitions[part_id]->NEList[v1].begin(), pragmatic_partitions[part_id]->NEList[v1].end(),
                                            pragmatic_partitions[part_id]->NEList[v2].begin(), pragmatic_partitions[part_id]->NEList[v2].end(),
                                            std::inserter(intersection, intersection.begin()));*/
                /* std::vector<index_t> NNList_v1 = pragmatic_partitions[part_id]->get_nnlist(g2l_vertices_tmp.at(glob_v1));
                    std::vector<index_t> NNList_v2 = pragmatic_partitions[part_id]->get_nnlist(g2l_vertices_tmp.at(glob_v2));
                }

                std::cout << "AFFECTED ELEMENTS" << std::endl;
                for (size_t i = 0; i < outbox_data.num_verts(); ++i)
                {
                    ;
                }

                }

            }//end of if color > 0 */
/*
            std::chrono::duration<double> coords_time = std::chrono::system_clock::now() - coords_tic;
       
            times[2] += g2l_time;
            times[3] += l2g_time;
            times[4] += coords_time.count();
        */
/*
            g2l_build[omp_get_thread_num()] += g2l_time;
            l2g_build[omp_get_thread_num()] += l2g_time;
*/
            //Create ENLists for local partition and the element-index-mappings
  //          auto enlist_tic = std::chrono::system_clock::now();

            auto ctr {0};          

            std::vector<int>ENList_part;

            if (dim == 2)
            {
                ENList_part.reserve(3*num_elements_part);
            }

            else
            {
                ENList_part.reserve(4*num_elements_part);
            }

            auto counter {0};
            for (auto it : elements_part[part_id])
            {               
                //update l2g and g2l element mappings
                l2g_elements_tmp[ctr/3]=it;
                g2l_elements_tmp.insert( std::make_pair(it, ctr/3) );

                const int *element_ptr = nullptr;
                element_ptr = original_mesh->get_element(it);
                
                ENList_part[ctr++] = g2l_vertices_tmp[*(element_ptr++)];
                ENList_part[ctr++] = g2l_vertices_tmp[*(element_ptr++)];
                ENList_part[ctr++] = g2l_vertices_tmp[*(element_ptr++)]; //three times for triangles

                if (dim == 3)
                    ENList_part[ctr++] = g2l_vertices_tmp[*(element_ptr++)];
                //TODO: Update for 3D case!!!
            }

 /*           std::chrono::duration<double> enlist_time = std::chrono::system_clock::now() - enlist_tic;
            times[5] += enlist_time.count();
            //g2l_access[omp_get_thread_num()] += enlist_time.count();
*/
            //create pragamtic data structure, the partition boundary and put it into the partition vector
//            auto mesh_tic = std::chrono::system_clock::now();
            Mesh<double>* partition = nullptr;

            if (dim == 2)
                partition = new Mesh<double>(num_points_part, num_elements_part, &(ENList_part[0]), &(x_coords[0]), &(y_coords[0]));

            else
                partition = new Mesh<double>(num_points_part, num_elements_part, &(ENList_part[0]), &(x_coords[0]), &(y_coords[0]), &(z_coords[0]));

            std::cout << "  start mesh healing process for partition " << part_id << " with " << num_points_part << " vertices" << std::endl;
/*
            std::cout << "g2l_vertices_tmp" << std::endl;

            for (auto it : g2l_vertices_tmp)
                std::cout << it.first << " " << it.second << std::endl;
*/
            //Heal mesh if the partition has data in its outbox

            auto orig_NNodes {partition->get_number_nodes()};

            if (color > 0)
            {
                auto origNNodes = partition->get_number_nodes();

                for (auto it : partition_adjcy[part_id])
                {
                    //std::cout << std::endl << "   check if color of partition " << it << " is smaller than own color" << std::endl;
                    std::cout << "   check outbox of partition " << it << std::endl;
                    //first check if color of neighbor is smaller than own color, otherwise there is no data in the neighbor's outbox!!!
                    if (partition_colors[it] < color && outboxes[it].num_verts() > 0)
                    {   
                        //std::cout << std::endl << "   yes it is, start with processing outbox data" << std::endl;
                        std::cout << "   process outbox of partition " << it << std::endl;
                        orig_NNodes = partition->get_number_nodes();
                        auto orig_NElements = partition->get_number_elements();

                        //std::cout << " outbox from partition " << it << " has " << outboxes[it].num_verts() << " vertices and " << outboxes[it].size() << " entries" << std::endl;
                        /*for (size_t i = 0; i < outboxes[it].num_verts(); ++i)
                        {
                            std::cout << outboxes[it][4*i] << " " << outboxes[it][4*i+1] << " " << outboxes[it][4*i+2] << " "  << outboxes[it][4*i+3] << std::endl;
                        }*/

                        //std::cout << " get coords from outbox vertices" << std::endl;

                        for (size_t i = 0; i < outboxes[it].num_verts(); ++i)
                        {
                            double p[2] {0.0, 0.0};
                            original_mesh->get_coords(outboxes[it][4*i+1], p); 
                      //      std::cout << outboxes[it][4*i+1] << ": " << p[0] << " " << p[1] << std::endl;

                            original_mesh->get_coords(outboxes[it][4*i+2], p); 
                        //    std::cout << outboxes[it][4*i+2] << ": " << p[0] << " " << p[1] << std::endl;

                            pragmatic_partitions[it]->get_coords(outboxes[it][4*i+3], p); 
                          //  std::cout << outboxes[it][4*i+3] << ": " << p[0] << " " << p[1] << std::endl << std::endl;
                        }

                        auto verts_in_part = outboxes[it].verts_in_part(part_id);
                        //std::cout << " vertices in that outbox for this partition: " << verts_in_part << std::endl;
                        //auto verts_in_outbox_for_my_partition = std::count(outboxes[it].begin(), outboxes[it].end());
                        std::vector<int> outbox_mapping(verts_in_part, -1);
                        //std::vector<int> outbox_mapping(verts_in_outbox_for_my_partition);

                        //resize vector
                        new_vertices_per_element.resize(nedge*num_elements_part);
                        std::fill(new_vertices_per_element.begin(), new_vertices_per_element.end(), -1);

                        //Resize vectors
                        size_t reserve = verts_in_part + partition->NNodes;
                        //size_t reserve = verts_in_outbox_for_my_partition + partition->NNodes;

                        if(partition->_coords.size()<reserve*dim) 
                        {
                            //std::cout << " Resizing! " << reserve << std::endl;
                            partition->_coords.resize(reserve*dim);
                            partition->metric.resize(reserve*msize); //CHANGE 3 to 6 FOR 3D-CASE!!!!
                            partition->NNList.resize(reserve);
                            partition->NEList.resize(reserve);
                            partition->node_owner.resize(reserve);
                            partition->lnn2gnn.resize(reserve);
                        }
                        
                        auto edgeSplitCnt = partition->NNodes - origNNodes;

                        //std::cout << " Append new coords and new metrics to the partition" << std::endl;

                        //Append new coords and new metrics to the partition
                        for (size_t i = 0, j = 0; i < outboxes[it].num_verts(); ++i)
                        {
                            if (outboxes[it][4*i] != part_id)
                            {
                                //std::cout << i << " is not in this partition" << std::endl;
                                continue;
                            }
                            //get 3rd element due to construction of outbox vector (3rd ele is local id of vertex in partition it)
                            //std::cout << "outb " << i << std::endl;
                            double p[2] {0.0, 0.0};
                            pragmatic_partitions[it]->get_coords(outboxes[it][4*i+3], p);

                            double m[3] {0.0, 0.0, 0.0};
                            pragmatic_partitions[it]->get_metric(outboxes[it][4*i+3], m);

                            partition->append_vertex(p, m);
//                            std::cout << partition->get_number_nodes() << std::endl;

                            outbox_mapping[j] = partition->get_number_nodes()-1; 
                            ++j;
                        }

                        //std::cout << " outbox mapping has " << outbox_mapping.size() << " entries" << std::endl;
/*
                        for (size_t i = 0; i < outbox_mapping.size(); ++i)
                        {
                            std::cout << i << ": " << outbox_mapping[i] << std::endl;
                        }
*/
                        std::set<int> elements_to_heal;

                        //std::cout << "mark each element with its new vertices" << std::endl;

                        // Mark each element with its new vertices,
                        // update NNList for all split edges.
                        for (size_t i = 0, j = 0; i < outboxes[it].num_verts(); ++i)
                        {
                            if (outboxes[it][4*i] != part_id)
                            {
                                //std::cout << " skip " << i << std::endl;
                                continue;
                            }

                            //std::cout << " process entry " << i  << std::endl;

                            auto vid = outboxes[it][4*i+3];
                            auto glob_firstid = outboxes[it][4*i+1];
                            auto glob_secondid = outboxes[it][4*i+2];

                            
                            //std::cout << vid << " is " << local_vid << std::endl;
                            //std::cout << glob_firstid << " " << glob_secondid << " " << vid << std::endl;
                            double p[2] {0.0, 0.0};
                            original_mesh->get_coords(glob_firstid, p);
                         //   std::cout << p[0] << " " << p[1] << std::endl;
                            original_mesh->get_coords(glob_secondid, p);
                         //   std::cout << p[0] << " " << p[1] << std::endl;

                            auto firstid = g2l_vertices_tmp.at(glob_firstid);
                          //  std::cout << "firstid ok" << std::endl;
                            auto secondid = g2l_vertices_tmp.at(glob_secondid);
                         //   std::cout << "secondid ok" << std::endl;
                            auto local_vid = outbox_mapping[j];

                         //   std::cout << firstid << " " << secondid << " " << local_vid << std::endl;

                            // Find which elements share this edge and mark them with their new vertices.
                            std::set<int> NEList_firstid = partition->get_nelist(firstid);
                            std::set<int> NEList_secondid = partition->get_nelist(secondid);

                            std::set<index_t> intersection;
                            std::set_intersection(NEList_firstid.begin(), NEList_firstid.end(),
                                                  NEList_secondid.begin(), NEList_secondid.end(),
                                                  std::inserter(intersection, intersection.begin()));

                            for (auto it : intersection)
                            {
                                size_t edgeOffset = edgeNumber(partition, it, firstid, secondid);
                                new_vertices_per_element[nedge*it+edgeOffset] = local_vid;
                                elements_to_heal.insert(it);
                         }

                            /*std::cout << "  Update NNList for newly created vertices " << orig_NNodes << " " << j << std::endl;/*
                            std::cout << "NNList.size: " << partition->NNList.size() << std::endl;
                            std::cout << "NNodes: " << partition->get_number_nodes() << std::endl;
*/
                            //Update NNList for newly created vertices.                        
                            partition->NNList[orig_NNodes+j].push_back(firstid);
                            partition->NNList[orig_NNodes+j].push_back(secondid);

                            partition->remove_nnlist(firstid, secondid);
                            partition->add_nnlist(firstid, local_vid);
                            partition->remove_nnlist(secondid, firstid);
                            partition->add_nnlist(secondid, local_vid);

                            ++j;
                            //*/
                        } //end of update NNList for all split edges

                        /*
                        //DEBUG
                        for (size_t i = 0; i < outboxes[it].num_verts(); ++i)
                        {
                            double p[2] {0.0, 0.0};
                            partition->get_coords(outbox_mapping[i], p);

                            std::cout << "NNList for " << outbox_mapping[i] << " at " << p[0] << " " << p[1] << std::endl;

                            std::vector<int> nnlist = partition->get_nnlist(outbox_mapping[i]);

                            for (auto nnlist_it : nnlist)
                            {
                                std::cout << "  " << nnlist_it << std::endl;
                            }
                        }
                        //END OF DEBUG*/

                        /*
                        //DEBUG
                        for (size_t i = 0; i < new_vertices_per_element.size(); ++i)
                        {
                            if (i % 3 == 0)
                                std::cout << i / 3 << ": " << std::endl;

                            std::cout << " " << new_vertices_per_element[i] << std::endl;
                        }
                        //END OF DEBUG*/                        

                        // Start element healing
                        //std::cout << "start element healing" << std::endl;
                        auto splitCnt {0};

                        int origNElements = partition->get_number_elements();

                        for (auto ele_id : elements_to_heal)
                        {
                            for(size_t j=0; j<nedge; ++j)
                            {
                                if(new_vertices_per_element[nedge*ele_id+j] != -1) 
                                {
                                    //refine_element(eid, tid, l2g_elements, g2l_elements, glob_NElements);
                                    const int *n=partition->get_element(ele_id);

                                     // Note the order of the edges - the i'th edge is opposite the i'th node in the element.
                                    index_t newVertex[3] = {-1, -1, -1};
                                    newVertex[0] = new_vertices_per_element[nedge*ele_id];
                                    newVertex[1] = new_vertices_per_element[nedge*ele_id+1];
                                    newVertex[2] = new_vertices_per_element[nedge*ele_id+2];


                                    int heal_cnt=0;
                                    for(size_t i=0; i<3; ++i)
                                    {
                                        if(newVertex[i]!=-1)
                                        {
                                            ++heal_cnt;
                                        }
                                    }

                                    if (heal_cnt == 1)
                                    {
                                        heal_1(partition, newVertex, ele_id, nloc, splitCnt, origNElements);
                                    }

                                    else if (heal_cnt == 2)
                                    {
                                        heal_2(partition, newVertex, ele_id, nloc, splitCnt, origNElements);
                                    }

                                    break;
                                }
                            }
                        } //end of element healing*/
                    
                        /*
                        //DEBUG
                        for (size_t i = 0; i < new_vertices_per_element.size(); ++i)
                        {
                            if (i % 3 == 0)
                                std::cout << i / 3 << ": " << std::endl;

                            std::cout << " " << new_vertices_per_element[i] << std::endl;
                        } //END OF DEBUG*/
                    } //end of if (partition_colors[it] < color)*/
                }

                //std::cout << "NNodes: " << partition->get_number_nodes() << " " << outbox_data.num_verts() << std::endl;
            } //end of //Heal mesh if the partition has data in its outbox

  //          auto boundary_tic = std::chrono::system_clock::now();
            partition->create_boundary();
    //        std::chrono::duration<double> boundary_time = std::chrono::system_clock::now() - boundary_tic;

            pragmatic_partitions[part_id] = partition;
            l2g_vertex[part_id] = l2g_vertices_tmp;
            g2l_vertex[part_id] = g2l_vertices_tmp;
            l2g_element[part_id] = l2g_elements_tmp;
            g2l_element[part_id] = g2l_elements_tmp;
     //       std::chrono::duration<double> mesh_time = std::chrono::system_clock::now() - mesh_tic;
       /*  
            times[6] += mesh_time.count();
            times[7] += boundary_time.count();
        
            //Create metric field
            auto metric_tic = std::chrono::system_clock::now();
*/
            if (algorithm == "pragmatic")
            {
                //std::cout << "pragmatic" << std::endl;
                if (dim == 2)
                {
                    MetricField<double,2> metric_field(*partition);
    /*
                double eta = 0.0001;
                std::vector<double> psi(num_points_part);

                for (size_t i = 0; i < num_points_part; ++i)
                {
                    //double x = 2*partition->get_coords(i)[0]-1;
                    //double y = 2*partition->get_coords(i)[1]-1;
                    double x = partition->get_coords(i)[0];
                    double y = partition->get_coords(i)[1];
            
                    psi[i] = 0.100000000000000*sin(50*x) + atan2(-0.100000000000000, (double)(2*x - sin(5*y)));
                    //psi[i] = 800000000000;
                }
                
                metric_field.add_field(&(psi[0]), eta, 1);
    */

                    for (auto i {0}; i < partition->get_number_nodes(); ++i)
                    {
                        double m[] = {1.0, 1.0, 0.0};
                        metric_field.set_metric(m, i);
                    }
                    
                    //           auto metric_update_tic = std::chrono::system_clock::now();
                    metric_field.update_mesh();
    /*           std::chrono::duration<double> metric_update_time = std::chrono::system_clock::now() - metric_update_tic;*/
                }

                else
                {
                    MetricField<double,3> metric_field(*partition);

                    for (auto i {0}; i < num_points_part; ++i)
                    {
                        double m[] = {1.0, 1.0, 1.0, 0.0};
                        metric_field.set_metric(m, i);
                    }

                    //auto metric_update_tic = std::chrono::system_clock::now();
                    metric_field.update_mesh();
            /*      std::chrono::duration<double> metric_update_time = std::chrono::system_clock::now() - metric_update_tic;*/
                }
        //     std::chrono::duration<double> metric_time = std::chrono::system_clock::now() - metric_tic;
    /*
                times[8] += metric_time.count();
                times[9] += metric_update_time.count();   
                */
            }//end if algorithm == pragmatic for metric assignment

            //std::cout << "metric field updated" << std::endl;
 
            double int_check_time {0.0};
            double triangulate_time {0.0};
            double tri_ds_time{0.0};
          // auto refine_tic = std::chrono::system_clock::now();
            auto refine_tic = omp_get_wtime();
            if (algorithm == "pragmatic")
            {
                if (dim == 2)
                {
                    Refine<double,2> refiner(*partition);
                    auto triangulate_tic = omp_get_wtime();

                    std::cout << "refine partition " << part_id << std::endl;
                    //std::cout << orig_NNodes << " " << l2g_vertices_tmp.size() << std::endl;

                    //if (color == 0)
                    {
                        refiner.refine(0.0005, nodes_partition_ids, l2g_vertices_tmp, g2l_vertices_tmp, l2g_elements_tmp, g2l_elements_tmp,
                                   &ref_detail_log[0], num_nodes, num_elements, part_id, outbox_data, outboxes, partition_colors,
                                   partition_adjcy[part_id]); //*/
                    }

                    triangulate_time = omp_get_wtime() - triangulate_tic;
                    l2g_vertex[part_id] = l2g_vertices_tmp;
                    g2l_vertex[part_id] = g2l_vertices_tmp;
                    l2g_element[part_id] = l2g_elements_tmp;
                    g2l_element[part_id] = g2l_elements_tmp;
                    outboxes[part_id]=outbox_data;
                }

                else
                {
                // Refine<double,3> refiner(*partition);
                // refiner.refine(0.5, nodes_partition_ids, l2g_tmp, &int_check_time); 
                } 

                //if (color > 0)
               /* {
                    std::cout << std::endl << "Element Quality for Partition " << part_id << std::endl;
                    partition->print_quality();

                    std::cout << std::endl << "Worst Element Quality: " << partition->get_qmin() << std::endl;
                   // Swapping<double, 2> swapper(*partition);
                   // swapper.swap(0.7);
                }//*/
            }

            else if (algorithm == "triangle")
            {       
                struct triangulateio tri_partition, tri_out;
         auto tri_ds_tic = std::chrono::system_clock::now();
                //init tri_partition
                tri_partition.numberofpoints = partition->get_number_nodes();
                tri_partition.numberofpointattributes = 0;
                tri_partition.pointmarkerlist = NULL;
                //tri_partition.pointlist = (REAL *) malloc(tri_partition.numberofpoints * 2 * sizeof(REAL) );
                tri_partition.pointlist = partition->get_coords_pointer();
/*
                for (size_t i = 0; i < tri_partition.numberofpoints; ++i)
                {
                    tri_partition.pointlist[2*i] = partition->get_coords(i)[0];
                    tri_partition.pointlist[2*i+1] = partition->get_coords(i)[1];
                }
  */           
                tri_partition.numberoftriangles = partition->get_number_elements();
                tri_partition.numberofcorners = 3;
                tri_partition.numberoftriangleattributes = 0;
               // tri_partition.trianglelist = (int*) malloc ( tri_partition.numberoftriangles * 3 * sizeof(int) );
                tri_partition.trianglelist = partition->get_enlist_pointer();
/*
                for (size_t i = 0; i < tri_partition.numberoftriangles; ++i)
                {
                    const int *element_ptr = nullptr;
                    element_ptr = partition->get_element(i);
                    
                    tri_partition.trianglelist[3*i+0] = *(element_ptr++);
                    tri_partition.trianglelist[3*i+1] = *(element_ptr++);
                    tri_partition.trianglelist[3*i+2] = *(element_ptr++);
                }
                //end of init tri_partition
*/
                //init tri_out
                tri_out.pointlist = (REAL *) NULL;
                tri_out.pointmarkerlist = (int *) NULL;
                tri_out.pointattributelist = (REAL *) NULL;
                tri_out.trianglelist = (int *) NULL;
                tri_out.numberofpoints = 0;
                tri_out.numberofpointattributes = 0;
/*
                tri_out.triangleattributelist = (REAL *) NULL;
                tri_out.neighborlist = (int *) NULL;
                tri_out.segmentlist = (int *) NULL;
                tri_out.segmentmarkerlist = (int *) NULL;
                tri_out.edgelist = (int *) NULL;
                tri_out.edgemarkerlist = (int *) NULL;

                tri_out.trianglearealist = NULL;
                tri_out.numberoftriangles = 0;
                tri_out.numberofcorners = 0;
                tri_out.numberoftriangleattributes = 0;
                tri_out.numberofsegments = 0;

                tri_out.holelist = NULL;
                tri_out.numberofholes = 0;

                tri_out.regionlist = NULL;
                tri_out.numberofregions = 0;

                tri_out.edgelist = NULL;
                tri_out.edgemarkerlist = NULL;
                tri_out.normlist = NULL;
                tri_out.numberofedges = 0;
                //end of init tri_out*/
std::chrono::duration<double> tri_ds_dur = std::chrono::system_clock::now() - tri_ds_tic;
            tri_ds_time = tri_ds_dur.count();
                //copy options string
                char * options_buffer = new char[options.length()+1];
                std::strcpy(options_buffer, options.c_str());

                //triangulate
                auto triangulate_tic = omp_get_wtime();
                //viennamesh::info(1) << "Making mesh with options " << options << std::endl;
                triangulate (options_buffer, &tri_partition, &tri_out, NULL);
                triangulate_time = omp_get_wtime() - triangulate_tic;
                //end of triangulate*/

                //free all memory (ONLY FREE IF MEMORY HAS BEEN ALLOCATED, NOT IF POINTERS ARE USED)
             //   free(tri_partition.pointlist);
             //   free(tri_partition.trianglelist);
             //   free(tri_partition.pointmarkerlist);

                triangle_partitions[part_id] = tri_out;

               // free(tri_out.pointlist);
               // free(tri_out.trianglelist);

                delete[] options_buffer;
                //end of free all memory                
            }

            auto refine_toc = omp_get_wtime();
         
            auto omp_toc = omp_get_wtime();

            threads_log[omp_get_thread_num()]+= omp_toc - omp_tic;
            refine_times[omp_get_thread_num()]+= refine_toc - refine_tic;
            triangulate_log[omp_get_thread_num()] += triangulate_time;
            build_tri_ds[omp_get_thread_num()] += tri_ds_time;
            //int_check_log[omp_get_thread_num()] += int_check_time;
        }//end parallel for loop
    } //end for loop colors
    viennamesh::info(1) << "Successfully adapted the mesh" << std::endl;
    std::cout << "REWRITE triunsuitable, SEE TRIANGLE.H AND TRIANGLE.C FOR DETAILS!!!" << std::endl;
   /*
    for(size_t i =0; i < pragmatic_partitions.size(); ++i)
    {
        std::cerr << "partition " << i << ": " << std::endl;
        std::cerr << "  " << pragmatic_partitions[i]->get_number_nodes() << std::endl;
        std::cerr << "  " << pragmatic_partitions[i]->get_number_elements() << std::endl;
    }
*/
    return true;
}
//end of CreatePragmaticDataStructures_par

//WritePartitions
//
//Tasks: Writes Partitions
bool MeshPartitions::WritePartitions()
{
    viennamesh::info(1) << "Write " << pragmatic_partitions.size() << " partitions" << std::endl;

    //write partitions
    #pragma omp parallel for num_threads(8)
    for (size_t i = 0; i < pragmatic_partitions.size(); ++i)
    {
        std::cout << "Writing partition " << i << std::endl;
        std::cout << "  Vertex count = " << pragmatic_partitions[i]->get_number_nodes() << std::endl;
        std::cout << "  Cell count = " << pragmatic_partitions[i]->get_number_elements() << std::endl;
        
        std::string filename;
        filename += "examples/data/color_refinement/output/";
        filename += "output_";
        filename += "partition_";
        filename += std::to_string( i );
        
        VTKTools<double>::export_vtu(filename.c_str(), pragmatic_partitions[i]);
    }

  return true;
} 
//end of WritePartitions

//WriteMergedMesh
//
//Tasks: Merges all mesh partitions and writes a single mesh file onto disk
bool MeshPartitions::WriteMergedMesh(std::string filename)
{
    std::cerr << "WriteMergedMesh()" << std::endl;
    return true;
}
//end of WriteMergedMesh

//GetRefinementStats
//
//Tasks: Returns the number of vertices and elements for a partitioned mesh
void MeshPartitions::GetRefinementStats(int* nodes, int* elements, std::string algorithm)
{
    int tmp_nodes {0};
    int tmp_elements {0};

    if (algorithm == "pragmatic")
    {
        for (size_t i =0; i < pragmatic_partitions.size(); ++i)
        {
            tmp_nodes += pragmatic_partitions[i]->get_number_nodes();
            tmp_elements += pragmatic_partitions[i]->get_number_elements();
        }
    }

    else if (algorithm == "triangle")
    {
        for (size_t i =0; i < triangle_partitions.size(); ++i)
        {
            tmp_nodes += triangle_partitions[i].numberofpoints;
            tmp_elements += triangle_partitions[i].numberoftriangles;
        }
    }

    *nodes = tmp_nodes;
    *elements = tmp_elements;
}
//end of GetRefinementStats

//MeshPartitions::edgeNumber
//
//Tasks: Returns edge number
int MeshPartitions::edgeNumber(Mesh<double>*& partition, index_t eid, index_t v1, index_t v2)
{
    const int *n=partition->get_element(eid);

    auto dim = partition->get_number_dimensions();

    if(dim==2) {
        /* In 2D:
            * Edge 0 is the edge (n[1],n[2]).
            * Edge 1 is the edge (n[0],n[2]).
            * Edge 2 is the edge (n[0],n[1]).
            */
        if(n[1]==v1 || n[1]==v2) {
            if(n[2]==v1 || n[2]==v2)
                return 0;
            else
                return 2;
        } else
            return 1;
    } else { //if(dim=3)
        /*
            * In 3D:
            * Edge 0 is the edge (n[0],n[1]).
            * Edge 1 is the edge (n[0],n[2]).
            * Edge 2 is the edge (n[0],n[3]).
            * Edge 3 is the edge (n[1],n[2]).
            * Edge 4 is the edge (n[1],n[3]).
            * Edge 5 is the edge (n[2],n[3]).
            */
        if(n[0]==v1 || n[0]==v2) {
            if(n[1]==v1 || n[1]==v2)
                return 0;
            else if(n[2]==v1 || n[2]==v2)
                return 1;
            else
                return 2;
        } else if(n[1]==v1 || n[1]==v2) {
            if(n[2]==v1 || n[2]==v2)
                return 3;
            else
                return 4;
        } else
            return 5;
    }
} //end of MeshPartitions::edgeNumber

//MeshPartitions::heal_1(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int& splitCnt, int& threadIdx)
//
//Task: Heals mesh after neighboring partition has altered interface
void MeshPartitions::heal_1(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int& splitCnt, int& threadIdx)
{
    // Single edge split.

    const int *n=partition->get_element(eid);
    const int *boundary=&(partition->boundary[eid*nloc]);

    int rotated_ele[3];
    int rotated_boundary[3];
    index_t vertexID = -1;
    for(int j=0; j<3; j++)
        if(newVertex[j] >= 0) {
            vertexID = newVertex[j];

            rotated_ele[0] = n[j];
            rotated_ele[1] = n[(j+1)%3];
            rotated_ele[2] = n[(j+2)%3];

            rotated_boundary[0] = boundary[j];
            rotated_boundary[1] = boundary[(j+1)%3];
            rotated_boundary[2] = boundary[(j+2)%3];

            break;
        }
    assert(vertexID!=-1);

    const index_t ele0[] = {rotated_ele[0], rotated_ele[1], vertexID};
    const index_t ele1[] = {rotated_ele[0], vertexID, rotated_ele[2]};

    const index_t ele0_boundary[] = {rotated_boundary[0], 0, rotated_boundary[2]};
    const index_t ele1_boundary[] = {rotated_boundary[0], rotated_boundary[1], 0};

    index_t ele1ID;
    ele1ID = splitCnt;
    
    partition->add_nnlist(vertexID, rotated_ele[0]);
    partition->add_nnlist(rotated_ele[0], vertexID);

    partition->add_nelist_fix(rotated_ele[0], ele1ID, threadIdx);

    partition->add_nelist(vertexID, eid);
    partition->add_nelist(vertexID, ele1ID);

    partition->remove_nelist(rotated_ele[2], eid);
    partition->add_nelist_fix(rotated_ele[2], ele1ID, threadIdx);

    partition->replace_element(eid, ele0);
    partition->append_element(ele1);

    splitCnt+=1;
}
//end of MeshPartitions::heal_1(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int& splitCnt, int& threadIdx)

//MeshPartitions::heal_2(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int splitCnt, int& threadIdx)
//
//Task: Heals mesh after neighboring partition has altered interface
void MeshPartitions::heal_2(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int& splitCnt, int& threadIdx)
{
    const int *n=partition->get_element(eid);
    const int *boundary=&(partition->boundary[eid*nloc]);

    int rotated_ele[3]; 
    int rotated_boundary[3];
    index_t vertexID[2];

    for(int j=0; j<3; j++) 
    {
        if(newVertex[j] < 0) 
        {
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

    real_t ldiag0 = partition->calc_edge_length(rotated_ele[1], vertexID[0]);
    real_t ldiag1 = partition->calc_edge_length(rotated_ele[2], vertexID[1]);

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
    partition->add_nnlist(vertexID[0], vertexID[1]);
    partition->add_nnlist(vertexID[1], vertexID[0]);

    // vertexID[offset] and rotated_ele[offset+1] are the vertices on the diagonal
    partition->add_nnlist(vertexID[offset], rotated_ele[offset+1]);
    partition->add_nnlist(rotated_ele[offset+1], vertexID[offset]);

    // rotated_ele[offset+1] is the old vertex which is on the diagonal
    // Add ele2 in rotated_ele[offset+1]'s NEList
    partition->add_nelist_fix(rotated_ele[offset+1], ele2ID, threadIdx);

    // Replace eid with ele0 in NEList[rotated_ele[0]]
    partition->remove_nelist(rotated_ele[0], eid);
    partition->add_nelist_fix(rotated_ele[0], ele0ID, threadIdx);

    // Put ele0, ele1 and ele2 in vertexID[offset]'s NEList
    partition->add_nelist(vertexID[offset], eid);
    partition->add_nelist_fix(vertexID[offset], ele0ID, threadIdx);
    partition->add_nelist_fix(vertexID[offset], ele2ID, threadIdx);

    // vertexID[(offset+1)%2] is the new vertex which is not on the diagonal
    // Put ele0 and ele2 in vertexID[(offset+1)%2]'s NEList
    partition->add_nelist_fix(vertexID[(offset+1)%2], ele0ID, threadIdx);
    partition->add_nelist_fix(vertexID[(offset+1)%2], ele2ID, threadIdx);

    partition->replace_element(eid, ele1);
    partition->append_element(ele0);
    partition->append_element(ele2);

    splitCnt += 2;
}
//end of MeshPartitions::heal_2(Mesh<double>*& partition, const index_t *newVertex, int eid, int nloc, int splitCnt, int& threadIdx)
//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif
