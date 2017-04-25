#ifndef MESH_PARTITIONS_HPP
#define MESH_PARTITIONS_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "Mesh.h"
#include "MetricField.h"
#include "Refine.h"
#include "ElementProperty.h"
#include "Edge.h"

//TODO: DEBUG
#include "VTKTools.h"
//END OF DEBUG

//viennamesh includes
#include "viennameshpp/plugin.hpp"

//all other includes
#include "metis.h"
#include <unordered_map>
#include <numeric>  
#include <chrono>

#include "mesh_partitions_refinement.hpp"

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
         std::string filename, int thread);                                                  //Constructor
        ~MeshPartitions();                                                                    //Destructor

        std::vector<Mesh<double>*> pragmatic_partitions;                                      //Vector containing pointers to the pragmatic partitions

        bool MetisPartitioning();                                                             //Partition mesh using metis
        bool CreatePragmaticDataStructures_ser();                                             //Create Pragmatic Meshes storing the mesh partitions in serial
        bool CreatePragmaticDataStructures_par();                                             //Create Pragmatic Meshes storing the mesh partitions in parallel
        bool CreateNeighborhoodInformation();                                                 //Create neighborhood information for vertices and partitions
        bool ColorPartitions();                                                               //Color the partitions
        bool WritePartitions();                                                               //ONLY FOR DEBUGGING!
        bool RefineInterior();                                                                //Refinement without refining boundary elements
        bool WriteMergedMesh(std::string filename);                                           //Merges partitions into a single mesh and writes it
        bool RefinementKernel(int part, double L_max);

        int get_colors(){return colors;};
        int get_max(){return max;};
        std::vector<std::vector<int>>& get_color_partitions() {return color_partitions;};

    private:
/*
        bool MetisPartitioning(Mesh<double>* original_mesh, int num_regions);                 //Partition mesh using metis
        bool CreatePragmaticDataStructures_ser(Mesh<double>* original_mesh, int num_regions); //Create Pragmatic Meshes storing the mesh partitions in serial
        bool CreatePragmaticDataStructures_par(Mesh<double>* original_mesh, int num_regions); //Create Pragmatic Meshes storing the mesh partitions in parallel
        bool CreateNeighborhoodInformation(Mesh<double>* original_mesh, int num_regions);     //Create neighborhood information for vertices and partitions
        bool ColorPartitions(int num_regions);                                                               //Color the partitions
        bool WritePartitions();                                                               //ONLY FOR DEBUGGING!
*/

        Mesh<double>* original_mesh;
        int num_regions;
        int nthreads;

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
        std::vector<std::unordered_map<index_t, index_t>> global_to_local_index_mappings;
        std::vector<std::unordered_map<index_t, index_t>> local_to_global_index_mappings;
        std::vector<std::unordered_map<index_t, index_t>> global_to_local_element_mappings;
        std::vector<std::unordered_map<index_t, index_t>> local_to_global_element_mappings;

        //Neighborhood Information containers
        std::vector<std::set<int>> nodes_partition_ids;
        std::vector<std::set<int>> partition_adjcy;

        std::set<int>& get_nodes_partition_ids(int n){return nodes_partition_ids[n];};

        //Color information
        size_t colors;                                                                        //Stores the number of colors used
        std::vector<int> partition_colors;                                                    //Contains the color assigned to each partition
        std::vector<std::vector<int>> color_partitions;                                       //Contains the partition ids assigned to each color

        double calc_edge_length(int part_id, index_t x0, index_t y0);
        double calculate_quality(const index_t* n, int part_id);
        double update_quality(index_t element, int part_id);

        //DEBUG
        int max=0;
        //END OF DEBUG

}; //end of class MeshPartitions

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Helper Functions                                                              //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Implementation                                                                //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//Constructor
//
//Tasks: TODO
//MeshPartitions::MeshPartitions(Mesh<double>* original_mesh, int num_regions, std::string filename)
MeshPartitions::MeshPartitions(Mesh<double> * const orig_mesh, int nregions, std::string filename, int threads)
{    
    original_mesh = orig_mesh;
    num_regions = nregions;
    nthreads = threads;
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
{
    viennamesh::info(5) << "Calling Destructor of MeshPartitions" << std::endl;
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
    //ncommon = mesh->get_number_dimensions();
   // std::vector<idx_t> bdry = mesh->copy_boundary_vector();
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
    /*METIS_PartMeshDual (&num_elements,
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
                        npart.data());*/

   METIS_PartMeshNodal(&num_elements,
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
/*
    //DEBUG
    if (epart[i] > max)
      max = epart[i];
    //END OF DEBUG
*/
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

/*
  //DEBUG
  for (size_t i = 0; i < num_regions; ++i)
  {
      std::cout << "Partition " << i << " has the following neighbors: " << std::endl;

      for (auto iter : partition_adjcy[i])
      {
          std::cout << "  " << iter << std::endl;
      }
  }
  //END OF DEBUG
*/

  return true;
}
//end of CreateNeighborhoodInformation

//ColorParittions
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
  /*  std::cout << "  Partition | Color " << std::endl;
 
    for (size_t i = 0; i < partition_colors.size(); ++i)
    {
        std::cout << "          " << i << " | " << partition_colors[i] << std::endl;
    }
*//*
    std::cout << std::endl << "      Color | #Partitions " << std::endl;
 
    for (size_t i = 0; i < color_partitions.size(); ++i)
    {
        //std::cout << "          " << i << " | " << color_partitions[i].size() << std::endl;
    
        std::cout << "          " << i << " | ";
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
    global_to_local_index_mappings.resize(num_regions);
    local_to_global_index_mappings.resize(num_regions);

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
bool MeshPartitions::CreatePragmaticDataStructures_par()
{
    std::vector<int> _ENList_orig = original_mesh->get_enlist();
    std::vector<std::set<int>> nodes_part(num_regions);
    std::vector<std::set<int>> elements_part(num_regions);

    global_to_local_index_mappings.resize(num_regions);
    local_to_global_index_mappings.resize(num_regions);
    global_to_local_element_mappings.resize(num_regions);
    local_to_global_element_mappings.resize(num_regions);
    pragmatic_partitions.resize(num_regions);

    //DEBUG
    std::vector<std::vector<int>> threads(2);
    //END OF DEBUG

    //iterate colors
    for (size_t color = 0; color < colors; color++)
    {
        //std::cerr << "color " << color << " has " << color_partitions[color].size() << " partitions" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(1)
        for (size_t part_iter = 0; part_iter < color_partitions[color].size(); ++part_iter)
        {
            size_t part_id = color_partitions[color][part_iter];
            
            //get vertices and elements from original mesh
            for(size_t ele_id = 0; ele_id < original_mesh->get_number_elements(); ele_id++)
            {
                if (epart[ele_id] == part_id)
                {
                    //TODO: Update for 3D case!!!
                    nodes_part[part_id].insert(_ENList_orig[3*ele_id]);
                    nodes_part[part_id].insert(_ENList_orig[3*ele_id+1]);
                    nodes_part[part_id].insert(_ENList_orig[3*ele_id+2]);
                    elements_part[part_id].insert(ele_id);
                }
            }
            
            //get number of vertices and elements of the local partitions
            int num_points = nodes_part[part_id].size();
            int num_elements = elements_part[part_id].size();
            
            //create coordinate vectors, g2l- and l2g-index-mappings for the vertices
            std::unordered_map<int, int> g2l_tmp, l2g_tmp;
            int new_vertex_id = 0;    //new vertex id for local partition
            int ctr = 0;
            
            //coords vector
            std::vector<double> x_coords(num_points);
            std::vector<double> y_coords(num_points); 
            
            for (auto it : nodes_part[part_id])
            {
                double p[2];
                original_mesh->get_coords( it, p);
                
                x_coords[ctr] = p[0];
                y_coords[ctr] = p[1];
                //TODO: Update for 3d case!

                g2l_tmp.insert( std::make_pair(it, new_vertex_id) );
                l2g_tmp.insert( std::make_pair(new_vertex_id++, it) );

                ++ctr;
            }
            
            //put the index mappings into the global vectors
            global_to_local_index_mappings[part_id] = g2l_tmp;
            local_to_global_index_mappings[part_id] = l2g_tmp;

            //Create ENLists for local partition and the element-index-mappings
            std::unordered_map<int, int> g2l_ele_tmp, l2g_ele_tmp;
            ctr = 0;

            std::vector<int> ENList_part(3*num_elements);

            for (auto it : elements_part[part_id])
            {
                g2l_ele_tmp.insert( std::make_pair(it, ctr/3) );
                l2g_ele_tmp.insert( std::make_pair(ctr/3, it) );

                const index_t *element_ptr = nullptr;
                element_ptr = original_mesh->get_element(it);

                ENList_part[ctr++] = g2l_tmp[*(element_ptr++)];
                ENList_part[ctr++] = g2l_tmp[*(element_ptr++)];
                ENList_part[ctr++] = g2l_tmp[*(element_ptr++)]; //three times for triangles
                //TODO: Update for 3D case!!!
            }

            //put the element mappings into the global vectors
            global_to_local_element_mappings[part_id] = g2l_ele_tmp;
            local_to_global_element_mappings[part_id] = l2g_ele_tmp;

            //create pragamtic data structure, the partition boundary and put it into the partition vector
            Mesh<double>* partition = nullptr;
            partition = new Mesh<double>( num_points, num_elements, &(ENList_part[0]), &(x_coords[0]), &(y_coords[0]));
            partition->create_boundary();

            pragmatic_partitions[part_id] = partition;
            
            //TODO: GHOSTING!!!
            //std::cerr << "Ghosting is still waiting for implementation!" << std::endl;

            //Create metric field
            MetricField<double,2> metric_field(*partition);

            double eta = 0.0001;
            std::vector<double> psi(num_points);

            for (size_t i = 0; i < num_points; ++i)
            {
                //double x = 2*partition->get_coords(i)[0]-1;
                //double y = 2*partition->get_coords(i)[1]-1;
        
                //psi[i] = 0.100000000000000*sin(50*x) + atan2(-0.100000000000000, (double)(2*x - sin(5*y)));
                psi[i] = 0.0001;
            }

            metric_field.add_field(&(psi[0]), eta, 1);
            metric_field.update_mesh();

            std::cerr << partition->get_number_elements() << " " << partition->get_number_nodes() << std::endl;

            MeshPartitionsRefinement Refiner(pragmatic_partitions[part_id], nodes_partition_ids);
        

            //TODO: REFINEMENT!!!*/
            //Refine<double, 2> adapt(*partition);
/*
            if (RefinementKernel(part_id, sqrt(2.0)))
                std::cerr << "Refinement OK" << std::endl;

            else 
                std::cerr << "Refinement Error" << std::endl;
                */
        }
        //end parallel for loop
    } //end for loop colors

   /* //DEBUG
    for (size_t i =0; i < threads.size(); ++i)
    {
        std::cerr << "Thread " << i << ": " << std::endl;

        for (auto it : threads[i])
        {
            std::cerr << "  " << it << std::endl;
        }
    }
    //END OF DEBUG*/
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

//RefineInterior
//
//Tasks: Refines mesh partition but leaves boundary elements untouched
//Notes: Refinement runs serially!!!!
bool MeshPartitions::RefineInterior()
{
    std::cout << "RefineInterior()" << std::endl;
}
//end of RefineInterior

//WriteMergedMesh
//
//Tasks: Merges all mesh partitions and writes a single mesh file onto disk
bool MeshPartitions::WriteMergedMesh(std::string filename)
{
    std::cerr << "WriteMergedMesh()" << std::endl;
/*
    std::cerr << "get boundary info" << std::endl; 

    int nfacets;
    //std::vector<int> facets(original_mesh->get_number_elements() * (original_mesh->get_number_dimensions()+1));
    //std::vector<int> ids(original_mesh->get_number_elements() * (original_mesh->get_number_dimensions()+1));
    const int* facets = new int [original_mesh->get_number_elements() * (original_mesh->get_number_dimensions()+1)];
    const int* ids = new int [original_mesh->get_number_elements() * (original_mesh->get_number_dimensions()+1)];

    pragmatic_partitions[0]->get_boundary(&nfacets, &(facets), &(ids));

    for (size_t i = 0; i < nfacets; ++i)
    {
        std::cerr << i << ": " << facets[2*i] << " " << facets[2*i+1] << std::endl;
    }
*/
    /*
    //create merged ENList
    std::vector<index_t> merged_ENList(original_mesh->num_elements*3);

    //iterate over partitions
    int global_element_counter = 0;
    for (size_t i = 0; i < pragmatic_partitions.size(); ++i)  
    {
        for (size_t j = 0; j < pragmatic_partitions[i]->get_number_elements(); ++j)
        {
            const index_t *element_ptr = nullptr;
            element_ptr = pragmatic_partitions[i]->get_element(j);

            merged_ENList[3*global_element_counter  ] = local_to_global_index_mappings[i].at( *(element_ptr++) );
            merged_ENList[3*global_element_counter+1] = local_to_global_index_mappings[i].at( *(element_ptr++) );
            merged_ENList[3*global_element_counter+2] = local_to_global_index_mappings[i].at( *(element_ptr++) );  
            ++global_element_counter;
        }
    }
    //end of create merged ENList

    //ofstream object
    ofstream writer;
    writer.open(filename.c_str(), ios::out);

    //write header
    writer << "<?xml version=\"1.0\"?>" << std::endl;
    writer << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl;
    writer << " <UnstructuredGrid>" << std::endl;
    //end of write header
    */
    return true;
}
//end of WriteMergedMesh

bool MeshPartitions::RefinementKernel(int part_id, double L_max)
{
    ElementProperty<double>* property;
    std::vector< DirectedEdge<index_t>> allNewVertices;

    std::vector<DirectedEdge<index_t>> newVertices;
    std::vector<double> newCoords;
    std::vector<double> newMetric;
    std::vector<int> newAppearances;

    //TODO: Update those 3 for 3D refinement
    size_t nedge = 3;
    size_t dim = pragmatic_partitions[part_id]->get_number_dimensions();
    size_t msize = 3;
    size_t nloc = 3;

    //Set orientation of elements
    for (size_t i =0; i < pragmatic_partitions[part_id]->get_number_elements(); ++i)
    {
        //TODO update to work with both partitions and interfaces
        const int* n = pragmatic_partitions[part_id]->get_element(i);

        //element is marked for deletion
        if (n[0] < 0)
            continue;

        if (pragmatic_partitions[0]->get_number_dimensions() == 2)
        {
            //TODO update to work with both partitions and interfaces
            property = new ElementProperty<double>(pragmatic_partitions[part_id]->get_coords(n[0]), pragmatic_partitions[part_id]->get_coords(n[1]),
            pragmatic_partitions[part_id]->get_coords(n[2]));
        }
    }

    size_t origNNodes = pragmatic_partitions[part_id]->get_number_nodes();
    size_t origNElements = pragmatic_partitions[part_id]->get_number_elements();
    
    size_t splitCnt_wo_vertices = 0;
    size_t splitCnt = 0;

    std::vector<index_t> new_vertices_per_element(nedge*origNElements);
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

        //if vertex is on interface, continue
        if (nodes_partition_ids[i].size() > 1)
        {
            continue;
        }

        else 
        {
            NNList_i = pragmatic_partitions[part_id]->get_nnlist(i);
        }

        for (auto otherVertex : NNList_i) //otherVertex is already a global index
        {
            //check if one of the vertices does not reside in partition 0, this is only for test case!
            std::unordered_map<index_t, index_t>::iterator position = global_to_local_index_mappings[part_id].find(otherVertex);

            if (position == global_to_local_index_mappings[part_id].end())
            {
                continue;
            }
        
            //order vertices according to their global index
            if( local_to_global_index_mappings[part_id].at(i) < otherVertex )
            {
                double length = calc_edge_length(part_id, local_to_global_index_mappings[part_id].at(i), otherVertex);

                if (length > L_max)
                {
                    ++splitCnt;
                    //refine_edge(local_to_global_index_mappings[part_id].at(i), otherVertex, part_id);
                }
            }
        }
    }

    return true;
}

//calc_edge_length(index_t x0, index_t y0)
//calculates edge length in metric space using global indices
double MeshPartitions::calc_edge_length(int part_id, index_t x0, index_t y0)
{
  double length = -1.0;

  //TODO: adapt for 3D space!
  double m[3];
  double m_x[3], m_y[3];

  pragmatic_partitions[part_id]->get_metric(x0, m_x);
  pragmatic_partitions[part_id]->get_metric(y0, m_y);

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

  pragmatic_partitions[part_id]->get_coords(x0, x_coords);
  pragmatic_partitions[part_id]->get_coords(y0, y_coords);

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
/*
//void GroupedPartitionsRefinement::refine_edge_part(index_t x, index_t y, int part)
void MeshPartitions::refine_edge_part(index_t n0, index_t n1, int part_id)
{  
  //swap indices because lesser id has to be the first
  if(n0 > n1)
  {
    index_t tmp = n0;
    n0 = n1;
    n1 = tmp;
  }

  newVertices.push_back( DirectedEdge<index_t> (global_to_local_index_mappings[part_id].at(n0),
  global_to_local_index_mappings[part_id].at(n1)) );  //n0 and n1 are global indices    
  
  //Calculate position of new point
  double x, m;
  double n0_coords[2], n1_coords[2];
  double m0[3], m1[3];
  
  pragmatic_partitions[part_id]->get_coords(n0, n0_coords);
  pragmatic_partitions[part_id]->get_metric(n0, m0);

  pragmatic_partitions[part_id]->get_coords(n1, n1_coords);
  pragmatic_partitions[part_id]->get_metric(n1, m1);

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
*/
//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif