#ifndef MESH_PARTITIONS_HPP
#define MESH_PARTITIONS_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "Mesh.h"

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

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//class MeshPartitions
//
//This class partitions an input mesh and stores its partitions in different data structurs for further processing
class MeshPartitions
{
    public:
        MeshPartitions(Mesh<double> * const original_mesh, int num_regions, std::string filename);   //Constructor
        ~MeshPartitions();                                                                    //Destructor

        std::vector<Mesh<double>*> pragmatic_partitions;                                      //Vector containing pointers to the pragmatic partitions

        bool MetisPartitioning();                                                             //Partition mesh using metis
        bool CreatePragmaticDataStructures_ser();                                             //Create Pragmatic Meshes storing the mesh partitions in serial
        bool CreatePragmaticDataStructures_par();                                             //Create Pragmatic Meshes storing the mesh partitions in parallel
        bool CreateNeighborhoodInformation();                                                 //Create neighborhood information for vertices and partitions
        bool ColorPartitions();                                                               //Color the partitions
        bool WritePartitions();                                                               //ONLY FOR DEBUGGING!
        bool RefineInterior();                                                                //Refinement without refining boundary elements

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

        //Neighborhood Information containers
        std::vector<std::set<int>> nodes_partition_ids;
        std::vector<std::set<int>> partition_adjcy;

        //Color information
        size_t colors;                                                                        //Stores the number of colors used
        std::vector<int> partition_colors;                                                    //Contains the color assigned to each partition
        std::vector<std::vector<int>> color_partitions;                                       //Contains the partition ids assigned to each color


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
MeshPartitions::MeshPartitions(Mesh<double> * const orig_mesh, int nregions, std::string filename)
{    
    original_mesh = orig_mesh;
    num_regions = nregions;
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

    std::cout << std::endl << "      Color | #Partitions " << std::endl;
 *//*
    for (size_t i = 0; i < color_partitions.size(); ++i)
    {
        std::cout << "          " << i << " | " << color_partitions[i].size() << std::endl;
/*      
        std::cout << "          " << i << " | ";
        for (auto it : color_partitions[i])
        {
           std::cout << it << " ";
        }

        std::cout << std::endl;*/
    //}
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

    //iterate colors
    for (size_t color = 0; color < colors; color++)
    {
        #pragma omp parallel for schedule(guided) num_threads(2)
        for (size_t part_id = 0; part_id < color_partitions[color].size(); ++part_id)
        {
            //get vertices and elements from original mesh
            for(size_t ele_id = 0; ele_id < original_mesh->get_number_elements(); ele_id++)
            {
                ;
            }
        }
        //end parallel for loop
    } //end for loop colors
}
//end of CreatePragmaticDataStructures_par

//WritePartitions
//
//Tasks: Writes Partitions
bool MeshPartitions::WritePartitions()
{
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
#endif
