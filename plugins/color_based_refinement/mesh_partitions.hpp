#ifndef MESH_PARTITIONS_HPP
#define MESH_PARTITIONS_HPP

//pragmatic basic includes for data structure
//TODO: implement own data structure
#include "Mesh.h"

//viennamesh includes
#include "viennameshpp/plugin.hpp"

#include "metis.h"
#include <unordered_map>
#include <numeric>  

//all other includes

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//class MeshPartitions
//
//This class partitions an input mesh and stores its partitions in different data structurs for further processing
class MeshPartitions
{
    public:
        MeshPartitions(Mesh<double>* original_mesh, int num_regions);                   //Constructor
        ~MeshPartitions();                                                              //Destructor

        std::vector<Mesh<double>*> pragmatic_partitions;                                //Vector containing pointers to the pragmatic partitions

    private:
        bool MetisPartitioning(Mesh<double>* original_mesh, int num_regions);           //Partition mesh using metis
        bool CreatePragmaticDataStructures(Mesh<double>* original_mesh, int num_regions);                //Create Pragmatic Meshes storing the mesh partitions

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
MeshPartitions::MeshPartitions(Mesh<double>* original_mesh, int num_regions)
{
    MetisPartitioning(original_mesh, num_regions);

    //CreatePragmaticDataStructures(original_mesh, num_regions);
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
bool MeshPartitions::MetisPartitioning(Mesh<double>* const mesh, int num_regions)
{
    //get basic mesh information
    num_elements = mesh->get_number_elements();
    num_nodes = mesh->get_number_nodes();
    ncommon = mesh->get_number_dimensions();
    std::vector<idx_t> bdry = mesh->copy_boundary_vector();
    size_t num_bdry_nodes = std::accumulate(bdry.begin(), bdry.end(), 0);
    idx_t numflag = 0;     //0...C-style numbering is assumed that starts from 0; 1...Fortran-style numbering is assumed that starts from 1

    std::cout << num_bdry_nodes << std::endl;

    //reserve memory
    epart.reserve(num_elements);
    npart.reserve(num_nodes);

 /*   std::vector<idx_t> xadj;
    xadj.resize(num_nodes+1);
    std::vector<idx_t> adjncy;
    adjncy.resize(2*(3*num_nodes - 3 - num_bdry_nodes)); //see: http://math.stackexchange.com/questions/1541125/total-number-of-edges-in-a-triangle-mesh-with-n-vertices
*/
    int num_edges = 2*(3*num_nodes - 3 - num_bdry_nodes);
    idx_t *xadj;
    xadj = new idx_t[num_nodes+1];
    idx_t *adjncy;
    adjncy= new idx_t[num_edges];

 //   std::cout << adjncy.size() << std::endl;

    //fill eptr and eind, as done in viennamesh plugin "metis", file "mesh_partitionig.cpp"  
    eptr.push_back(0);
    
    for (size_t i = 0; i < mesh->get_number_elements(); ++i)
    {
        const index_t* element_ptr = nullptr;
        element_ptr = mesh->get_element(i);

        for (size_t j = 0; j < (mesh->get_number_dimensions() + 1); ++j)
        {
            eind.push_back( *(element_ptr+j) );
        }               

        eptr.push_back( eind.size() );    
    }  

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

 /*   METIS_PartMeshNodal(&num_elements,
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
*/
    METIS_MeshToDual(&num_elements,
                     &num_nodes,
                     eptr.data(),
                     eind.data(),
                     &ncommon,
                     &numflag,
                     &xadj,
                     &adjncy);
    
    viennamesh::info(5) << "Created " << num_regions << " mesh partitions using Metis" << std::endl;

    std::cout << "xadj: " << std::endl;
    for (size_t i = 0; i < num_nodes+1; ++i)
        std::cout << " " << i << ": " << xadj[i] << std::endl;

    std::cout << "adjncy: " << std::endl;
    for (size_t i = 0; i < num_edges; ++i)
        std::cout << " " << i << ": " << adjncy[i] << std::endl;


    return true;
}//end of MetisPartitioning

//CreatePragmaticDataStructures
//
bool MeshPartitions::CreatePragmaticDataStructures(Mesh<double>* const original_mesh, int num_regions)
{
    //reserve memory
    nodes_per_partition.resize(num_regions);
    elements_per_partition.resize(num_regions);
    global_to_local_index_mappings.resize(num_regions);
    local_to_global_index_mappings.resize(num_regions);

    //get ENList
    _ENList = original_mesh->get_element_node();

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

    viennamesh::info(5) << "Created " << pragmatic_partitions.size() << " pragmatic meshes" << std::endl;

    return true;
}
//Tasks: Get and order data needed to create a pragmatic data structure for each partition
//end of CreatePragmaticDataStructures

#endif