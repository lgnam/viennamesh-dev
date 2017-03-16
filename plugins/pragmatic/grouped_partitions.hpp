#ifndef GROUPED_PARTITIONS_HPP
#define GROUPED_PARTITIONS_HPP

//pragmatic basic includes for mesh data structure
#include "ticker.h"
#include "VTKTools.h"
#include "Mesh.h"
#include "MetricField.h"
#include "Lock.h"
#include "ElementProperty.h"

//all other includes
#include <time.h>
//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//
class GroupedPartitions
{
  //list of friend classes
  friend class GroupedPartitionsSmooth;
  friend class GroupedPartitionsRefinement;

  public:
    GroupedPartitions(Mesh<double>* input_mesh, int region_count);              //Constructor
    ~GroupedPartitions();                                                       //Destructor
   
    void WriteMergedMesh(std::string filename = "examples/data/pragmatic_metis_partitioning/merged_mesh.vtu");
    void PrintQuality();

  private:
    Mesh<double>* mesh = nullptr;                                               //input mesh (in pragmatic data structure)

    int num_partitions;
    int num_interfaces;

    std::vector<std::set<int>> _NEList;                                         //TODO: replace this, since its unnecessarily copying data!
    std::vector<index_t> _ENList;                                               //TODO: replace this, since its unnecessarily copying data!
    std::vector<std::set<index_t>> nodes_per_partition;
    std::vector<std::set<index_t>> initial_nodes_per_partition;
    std::vector<std::set<index_t>> elements_per_partition;
    std::vector<std::vector<std::set<index_t>>> interface_nodes;
    std::vector<std::vector<std::vector<index_t>>> interface_sets;
    std::vector<std::vector<std::set<index_t>>> interface_elements_sets;  

    std::vector<size_t> element_counter_interfaces;

    //pragmatic mesh containers
    std::vector<Mesh<double> * > pragmatic_partitions;
    std::vector<Mesh<double> * > pragmatic_interfaces;

    //index mappings for the partitions
    std::vector<std::unordered_map<index_t, index_t>> g2l_vertex_map_partitions;
    std::vector<std::unordered_map<index_t, index_t>> l2g_vertex_map_partitions;

    std::vector<std::unordered_map<index_t, index_t>> g2l_element_map_partitions;
    std::vector<std::unordered_map<index_t, index_t>> l2g_element_map_partitions;

    //vectors storing the mapping information for the interfaces
    std::vector<std::unordered_map<index_t, index_t>> g2l_vertex_map_interfaces;
    std::vector<std::unordered_map<index_t, index_t>> l2g_vertex_map_interfaces;

    std::vector<std::unordered_map<index_t, index_t>> g2l_element_map_interfaces;
    std::vector<std::unordered_map<index_t, index_t>> l2g_element_map_interfaces;

    //vectors storing the number of appearances of the global indices in the local meshes
    std::vector<index_t> element_appearances;
    std::vector<index_t> vertex_appearances;

    //vectors storing the boundaries of the whole mesh, the partitions and the interfaces
    std::vector<int> boundary_nodes_mesh;
    std::vector<std::vector<index_t>> boundary_nodes_partitions;
    std::vector<std::vector<index_t>> boundary_nodes_interfaces;

    //TODO: vectors storing sets of neighboring information about partitions and interfaces 
    //std::vector<std::set<index_t>> partition_neighbors;
    //std::vector<std::set<index_t>> interface_neighbors; //stores the neighboring partitions!
    std::vector<std::vector<index_t>> partition_neighbors;
    std::vector<std::vector<index_t>> interface_neighbors;
  
    //variables used by metis
    std::vector<idx_t> eptr;    //vectors storing the mesh as described in the Metis manual
    std::vector<idx_t> eind;    //vectors storing the mesh as described in the Metis manual
    idx_t num_nodes;
    idx_t num_elements;
    idx_t ncommon;  
    idx_t nparts;
    idx_t result;
    std::vector<idx_t> epart;
    std::vector<idx_t> npart;            //nparts not npart!!!

    //helper functions
    void CreateMetisPartitioning();
    void CreatePartitionsAndInterfaces();
    void CreatePragmaticDataStructures();
    void WritePartitionsAndInterfaces(std::string name = "pragmatic");
    void GetPartitionAndInterfaceBoundaries();
    int GetVertexPartitionOrInterface(index_t n, bool* partition, bool *interface);
    void GetVertexPartitionAndInterface(index_t n, int *partition, int *interface);
    int GetElementPartitionOrInterface(index_t n);
    void GetElementPartitionAndInterface(index_t n, int *partition, int *interface);
    void GetAllVertexPartitionsAndInterfaces(index_t n, std::vector<bool>& partitions, std::vector<bool>& interfaces);

    void GetCoords(index_t n, double *vertex);
    void SetCoords(index_t n, double *vertex);
    std::vector<index_t>* GetNNList(index_t n, int *partitions, int *interfaces);
    bool GetNNList_adv(index_t n, std::vector<index_t> &_NNList);
    bool GetNEList(index_t n, std::set<index_t> &_NEList, int *partitions, int *interfaces); 
    //std::set<index_t> GetNEList(index_t n, std::set<index_t> &_NEList, int *partitions, int *interfaces);
    bool GetNEList_adv(index_t n, std::set<index_t> &_NEList);
    void GetNEList_part(index_t n, std::set<int>& _NEList, int partition);
    void GetNEList_inter(index_t n, std::set<int>& _NEList, int interface);
    void GetCoords_inter(index_t n, double *vertex, int interface);
    void GetCoords_part(index_t n, double *vertex, int partition);
    void SetCoords_part(index_t n, double *vertex, int partition);
    void SetCoords_inter(index_t n, double *vertex, int interface);
    void AddCoords(size_t threadIdx, std::vector<double> coords, size_t dim, size_t splitCnt);
    void AddCoords_part(size_t threadIdx, std::vector<double> coords, size_t dim, size_t splitCnt, size_t part);
    void AddMetric(size_t threadIdx, std::vector<double> metric, size_t msize, size_t splitCnt);
    void AddMetric_part(size_t threadIdx, std::vector<double> metric, size_t msize, size_t splitCnt, size_t part);
    void GetMetric(index_t n, double *metric);
    void AddNNList(size_t id, size_t add_id);
    void AddNNList_part(size_t id, size_t add_id, int part);
    void RemNNList(size_t id, size_t rem_id);
    void RemNNList_part(size_t id, size_t rem_id, int part);
    void AddNEList(size_t id, size_t add_id);
    void AddNEList_part(size_t id, size_t add_id, int part);
    void AddNEList_fix(size_t id, size_t add_id, size_t threadIdx);
    void AddNEList_fix_part(size_t id, size_t add_id, size_t threadIdx, int part);
    void RemNEList(size_t id, size_t rem_id);
    void RemNEList_part(size_t id, size_t rem_id, int part);

    double calculate_quality(const index_t* n);
    double calculate_quality_part(const index_t* n, int part);
    double update_quality(index_t element);
    double update_quality_part(index_t element, int part);

}; //end of class

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Helper Functions                                                              //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//binomial_coefficient(int n, int k)
//TODO: implement a more efficient dynamic version (e.g., http://www.geeksforgeeks.org/dynamic-programming-set-9-binomial-coefficient/)
unsigned int binomial_coefficient (int n, int k)
{ 
  //check inputs
  if (k == 0 || k == n)
  {
    return 1;
  }

  //use recursion for computation
  return (binomial_coefficient(n-1, k-1) + binomial_coefficient(n-1, k));
} //end of binomial_coefficient(int n, int k)

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Implementation                                                                //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//TODO:Constructor
//TODO:conversion from ViennaMesh data structure into pragmatic data structure can be done here!?!?!?!?!?!?
//TODO: REPLACE _NEList and _ENList function, since it's copying data unnecessarily!!!
//TODO: use element initializer list!!!
GroupedPartitions::GroupedPartitions(Mesh<double>* input_mesh, int region_count) : num_nodes(input_mesh->get_number_nodes()), num_elements(input_mesh->get_number_elements()), ncommon(input_mesh->get_number_dimensions()), nparts(region_count), epart(input_mesh->get_number_elements()), npart(input_mesh->get_number_nodes()), nodes_per_partition(region_count), elements_per_partition(region_count), interface_nodes(region_count-1), interface_sets(region_count-1), _NEList(input_mesh->get_node_element()), _ENList(input_mesh->get_element_node()), interface_elements_sets(region_count-1), element_counter_interfaces(input_mesh->get_number_elements(), 0), g2l_vertex_map_partitions(region_count), l2g_vertex_map_partitions(region_count), g2l_element_map_partitions(region_count), l2g_element_map_partitions(region_count), element_appearances(input_mesh->get_number_elements(), 0), vertex_appearances(input_mesh->get_number_nodes(), 0), g2l_vertex_map_interfaces(binomial_coefficient(region_count, 2)), l2g_vertex_map_interfaces(binomial_coefficient(region_count, 2)), g2l_element_map_interfaces(binomial_coefficient(region_count, 2)), l2g_element_map_interfaces(binomial_coefficient(region_count, 2)), boundary_nodes_mesh(input_mesh->get_number_nodes(), 0), boundary_nodes_partitions(region_count), boundary_nodes_interfaces(binomial_coefficient(region_count, 2)), num_partitions(region_count), num_interfaces(binomial_coefficient(region_count, 2)), initial_nodes_per_partition(region_count), partition_neighbors(region_count), interface_neighbors(binomial_coefficient(region_count, 2))
{
  std::cout << "Grouped Partitions Object created" << std::endl;
  mesh = input_mesh;

  clock_t tic = clock();
  CreateMetisPartitioning();
  clock_t toc = clock();
  std::cout << "CreateMetisPartitioning(): " <<  static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;

  tic = clock();
  CreatePartitionsAndInterfaces();
  toc = clock();
  std::cout << "CreatePartitionsAndInterfaces(): " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;

  tic = clock();
  CreatePragmaticDataStructures();
  toc = clock();
  std::cout << "CreatePragmaticDataStructures(): " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;

  tic = clock();
  GetPartitionAndInterfaceBoundaries();
  toc = clock();
  std::cout << "GetPartitionAndInterfaceBoundaries(): " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;

  tic = clock();
  WritePartitionsAndInterfaces();
  toc = clock();
  std::cout << "WritePartitionsAndInterfaces(): " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;
/*
  int counter=0;
  int doubles = 0;
  int triples = 0;
  int quadruples = 0;
  for (size_t i = 0; i < num_nodes; ++i)
  {
    if (vertex_appearances[i] > 1)
      ++counter;

    if (vertex_appearances[i] == 2)
    {
      GetNNList_adv(i);
      ++doubles;
    }

    if (vertex_appearances[i] == 3)
      ++triples;

    if (vertex_appearances[i] == 4)
      ++quadruples;

    if (vertex_appearances[i] > 4)
      std::cout << "oho" << std::endl;
  }

  std::cout << counter << std::endl;
  std::cout << std::accumulate(vertex_appearances.begin(), vertex_appearances.end(), 0) - num_nodes - doubles - 2*triples - 3*quadruples<< std::endl;
  std::cout << "  doubles: " << doubles << "  triples: " << triples << " and quadruples: " << quadruples << std::endl;
*/
}
//end of constructor

//TODO:Destructor
GroupedPartitions::~GroupedPartitions()
{
  //TODO: count if not working as it is supposed to!!!
  int count = std::count_if (vertex_appearances.begin(), vertex_appearances.end(), [&](int n){return ((vertex_appearances[n] > 2) ? true : false);});
  std::cout << "Vertex Appearances > 2: " << count << std::endl; 

  std::cout << "Grouped Partitions Object deleted" << std::endl;
}
//end of destructor

//Create Metis Partitioning
void GroupedPartitions::CreateMetisPartitioning()
{
  std::cout << "Calling METISPartMeshNodal" << std::endl;

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
  METIS_PartMeshDual (&num_elements,
                      &num_nodes,
                      eptr.data(),
                      eind.data(),
                      NULL,
                      NULL,
                      &ncommon,
                      &nparts,
                      NULL,
                      NULL,
                      &result,
                      epart.data(),
                      npart.data());
//*/ 
/*    //This Method leads sometimes to "non-connected areas of partitions", therfore it is better to use METIS_PartMeshDual 
  METIS_PartMeshNodal(&num_elements,
                      &num_nodes,
                      eptr.data(),
                      eind.data(),
                      NULL,
                      NULL,
                      &nparts,
                      NULL,
                      NULL,
                      &result,
                      epart.data(),
                      npart.data());
//*/
}
//end of GroupedPartitions::CreateMetisPartitioning()

void GroupedPartitions::CreatePartitionsAndInterfaces()
{
  //TODO: bugfix interface creation, since dangling elements can occur! (can be seen when using a metric field!!!)
  std::cout << "Creating Partitions and Interfaces" << std::endl;

  //get the nodes per partition
  //this for-loop is not parallelizable, since set.insert() is not thread-safe
  for (size_t i = 0; i < num_elements; ++i)
  {
    //add nodes
    //TODO: the following hast to be updated for the 3D case!!
    initial_nodes_per_partition[ epart[i] ].insert(_ENList[i*3]);
    initial_nodes_per_partition[ epart[i] ].insert(_ENList[i*3+1]);
    initial_nodes_per_partition[ epart[i] ].insert(_ENList[i*3+2]);

    //add element
    elements_per_partition[ epart[i] ].insert(i);
  } 
  //end of get nodes per partition

  //get the interface nodes for all partition boundaries
  int interface_counter = 0;
  for (size_t i = 0; i < nparts; ++i)     
  {        
    for (size_t j = i+1; j < nparts; ++j, ++interface_counter)
    {   
      std::vector<index_t> tmp;

      set_intersection( initial_nodes_per_partition[i].begin(), initial_nodes_per_partition[i].end(),
                        initial_nodes_per_partition[j].begin(), initial_nodes_per_partition[j].end(),
                        inserter(tmp, tmp.begin()) );
      
      interface_sets[i].push_back(tmp);

      //get neighboring information
      if (tmp.size() > 0)
      {
        if ( std::find(partition_neighbors[i].begin(), partition_neighbors[i].end(), j) == partition_neighbors[i].end() )
        {
          partition_neighbors[i].push_back(j);
        }

        if ( std::find(partition_neighbors[j].begin(), partition_neighbors[j].end(), i) == partition_neighbors[j].end() )
        {
          partition_neighbors[j].push_back(i);
        }
        //partition_neighbors[i].insert(j);
        //partition_neighbors[j].insert(i);
        if ( std::find(interface_neighbors[interface_counter].begin(), interface_neighbors[interface_counter].end(), i) == interface_neighbors[interface_counter].end() )
        {
          interface_neighbors[interface_counter].push_back(i);
        }

        if ( std::find(interface_neighbors[interface_counter].begin(), interface_neighbors[interface_counter].end(), j) == interface_neighbors[interface_counter].end() )
        {
          interface_neighbors[interface_counter].push_back(j);
        }
        //interface_neighbors[interface_counter].insert(i);
        //interface_neighbors[interface_counter].insert(j);
      }
    }      
  } 
  //end of get the interface nodes for all partition boundaries

  //remove interface nodes from initial_nodes_per_partition
  for (auto interface_nodes : interface_sets)
  {
    for (size_t i = 0; i < interface_nodes.size(); ++i)
    {
      for (size_t j = 0; j < initial_nodes_per_partition.size(); ++j)
      {
        for (auto it : interface_nodes[i])
        {
          initial_nodes_per_partition[j].erase(it);
        }
      }
    }
  } 
  //end of remove interface nodes from initial_nodes_per_partition

  //get partition interface elements
  for (size_t i = 0; i < interface_sets.size(); ++i)
  {
    for (size_t j = 0; j < interface_sets[i].size(); ++j)
    {
      std::set<index_t> tmp;
      for (size_t k = 0; k < interface_sets[i][j].size(); ++k)
      {            
        for (auto NE_it : _NEList[interface_sets[i][j][k]])
        {  
          //element must not appear in more than one interface mesh!
          if (element_counter_interfaces[NE_it] == 0)
          {          
            tmp.insert(NE_it);
            ++element_counter_interfaces[NE_it];
          }
        }
      }
      interface_elements_sets[i].push_back(tmp);
      tmp.clear();
    }
  }
  //end of get partition interface elements

  //add the vertices of the partition interface elements to interface_sets
  //TODO: this seems to be overload, since a set in interface_sets storing the interface nodes is already available
  for (size_t i = 0; i < interface_elements_sets.size(); ++i)
  { 
    for (size_t j = 0; j < interface_elements_sets[i].size(); ++j)
    {
      std::set<index_t> tmp;

      for (auto interface_element : interface_elements_sets[i][j])
      {
        tmp.insert(_ENList[interface_element*3]);
        tmp.insert(_ENList[interface_element*3+1]);
        tmp.insert(_ENList[interface_element*3+2]);
      }
      interface_nodes[i].push_back(tmp);
    }
  }
  //end of add the vertices of the partition interface elements to interface_sets

  //remove partition interface elements from the actual partitions
  for (auto interface_elements : interface_elements_sets)
  {
    for (size_t i = 0; i < interface_elements.size(); ++i)
    {
      for (size_t j = 0; j < elements_per_partition.size(); ++j)
      {
        for (auto it : interface_elements[i])
        {
          elements_per_partition[j].erase(it);
        }
      }
    }
  } 
  //end of remove partition interface elements from the actual partitions

  //now get vertices from global ENList
  for (size_t i = 0; i < nparts; ++i)
  {
    for (auto element : elements_per_partition[i])
    {
      nodes_per_partition[i].insert( _ENList[3*element] );
      nodes_per_partition[i].insert( _ENList[3*element+1]);
      nodes_per_partition[i].insert( _ENList[3*element+2]);
    }
  }
  //end of get vertices from global ENList
}
//end of void GroupedPartitions::CreatePartitionsAndInterfaces()

void GroupedPartitions::CreatePragmaticDataStructures()
{
  std::cout << "Creating Pragmatic Data Structures" << std::endl;

  //Create Partitions

  //vectors storing the coordinate
  std::vector< std::vector<double>> x_coords(nparts);
  std::vector< std::vector<double>> y_coords(nparts);
  std::vector< std::vector<double>> z_coords(nparts);

  //vector storing the ENLists of each region
  std::vector<std::vector<index_t>> ENLists_partitions(nparts);

  //loop over all partitions
  for (size_t i = 0; i < nparts; ++i)
  {
    //get number of vertices and elements
    int num_points = nodes_per_partition[i].size();
    int num_cells = elements_per_partition[i].size();

    //get the vertex-to-index-mapping between old and new indices
    //and additionally the index-to-vertex-mapping
    std::unordered_map <index_t, index_t> global_to_local_index_map_vertices;
    std::unordered_map <index_t, index_t> local_to_global_index_map_vertices;

    index_t new_vertex_id = 0;
    for (auto it : nodes_per_partition[i])
    {
      global_to_local_index_map_vertices.insert( std::make_pair(it, new_vertex_id++) );
      ++vertex_appearances[it];
    }
    
    g2l_vertex_map_partitions[i] = global_to_local_index_map_vertices;

    //and get also the index-to-vertex mapping (opposite direction than vertex to index mapping)
    //TODO: put this in the for loop above
    for (auto it : global_to_local_index_map_vertices)
    {
      local_to_global_index_map_vertices[it.second] = it.first;
    }

    l2g_vertex_map_partitions[i] = local_to_global_index_map_vertices;

    //pre-allocate memory
    x_coords[i].reserve(num_points);
    y_coords[i].reserve(num_points);
    ENLists_partitions[i].resize(3*num_cells);

    //get coordinates of each vertex
    int counter = 0;
    for (auto it : nodes_per_partition[i])
    {
      double p[2];
      mesh->get_coords( it, p);
      x_coords[i][counter] = p[0];
      y_coords[i][counter] = p[1];
      ++counter;
    }

    //create ENList with respect to the new vertex indices and get also to mapping of indices for the elements
    std::unordered_map <index_t, index_t> global_to_local_index_map_elements;
    std::unordered_map <index_t, index_t> local_to_global_index_map_elements;

    counter=0;      
    for (auto it : elements_per_partition[i])
    {      
      global_to_local_index_map_elements.insert(std::make_pair(it, counter/3));
      local_to_global_index_map_elements.insert(std::make_pair(counter/3, it));

      const index_t *element_ptr = nullptr;
      element_ptr = mesh->get_element(it);
       
      ENLists_partitions[i][counter++] = global_to_local_index_map_vertices[*(element_ptr++)];
      ENLists_partitions[i][counter++] = global_to_local_index_map_vertices[*(element_ptr++)];
      ENLists_partitions[i][counter++] = global_to_local_index_map_vertices[*(element_ptr++)];  
      ++element_appearances[it];       
    }

    l2g_element_map_partitions[i] = local_to_global_index_map_elements;
    g2l_element_map_partitions[i] = global_to_local_index_map_elements;

    //create pragmatic mesh 
    Mesh<double> *partition_mesh = nullptr;

    //TODO: change for 3D refinement
    //mesh = new Mesh<double> ( num_points, num_cells, &(ENLists_regions[region.id()][0]) ,&(x_coords[region.id()][0]), &(y_coords[region.id()][0]), &(z_coords[region.id()][0]) );        
    partition_mesh = new Mesh<double> ( num_points, num_cells, &(ENLists_partitions[i][0]), &(x_coords[i][0]), &(y_coords[i][0]) );
    partition_mesh->create_boundary();

    //create metric
    MetricField<double,2> metric_field(*partition_mesh);
    std::vector<double> metric_vector;
    double eta=0.0001;

    //read data
    ifstream metric_data;
    std::string metric_file = "metric_data_part";
    metric_file += std::to_string(i);
    metric_file += ".txt";
    std::cout << metric_file <<std::endl;
    metric_data.open(metric_file.c_str());

    double metric_value;
    int ctr = 0;
    //for (size_t i = 0; i < num_points; ++i)
    while(metric_data >> metric_value)
    {
        metric_vector.push_back(metric_value);
        ++ctr;
    }

    //for (auto it : psi)
     // std::cerr << it << std::endl;

    for(size_t j = 0; j < num_points; j++) 
    {
    /*
      //std::cerr << i << ": " << j << std::endl;
      double x = 2*mesh->get_coords(j)[0]-1;
      double y = 2*mesh->get_coords(j)[1]-1;
      
      psi[j] = 0.100000000000000*sin(50*x) + atan2(-0.100000000000000, (double)(2*x - sin(5*y)));
      /*
      double x = partition_mesh->get_coords(j)[0];
      double y = partition_mesh->get_coords(j)[1];
      //std::cout << x*y << std::endl;
      //std::cerr << "partition calc psi " << j << std::endl;
      psi[j]= 1;
      //std::cerr << psi[j] << std::endl;*/
      //psi[j] = 10;
      double m[3];
      m[0] = metric_vector[3*j];
      m[1] = metric_vector[3*j+1];
      m[2] = metric_vector[3*j+2];
      //std::cerr << "setting field" << std::endl;
      metric_field.set_metric(m, j);
    }
    //TODO: check pragmatic function Mesh.h --> get_node_patch(nid, min_patch_size)
    //it does not work with test_box_2000x2000.vtu with 4 partitions!!!
    metric_data.close();
    //metric_field.add_field(&(psi[0]), 1);
    metric_field.update_mesh();
    
    //end of create metric

    pragmatic_partitions.push_back(partition_mesh);
  
    //delete partial_mesh;        //TODO: creates segfault if comments are removed
  } 
  //end of loop over all partitions

  //Create Interfaces
  int num_interfaces = binomial_coefficient(nparts, 2);
  
  //vectors storing the coordinate
  std::vector< std::vector<double>> x_coords_interfaces(num_interfaces);
  std::vector< std::vector<double>> y_coords_interfaces(num_interfaces);
  std::vector< std::vector<double>> z_coords_interfaces(num_interfaces);

  //vector storing the ENLists of each region
  std::vector<std::vector<index_t>> ENLists_interfaces(num_interfaces);

  //loop over all interfaces
  int interface_counter = 0 ;
  for(size_t i = 0; i < interface_elements_sets.size(); ++i)
  {
    for (size_t j = 0; j < interface_elements_sets[i].size(); ++j, interface_counter++)
    {
      //get number of vertices and elements
      int num_points = interface_nodes[i][j].size();
      int num_cells = interface_elements_sets[i][j].size();

      //get the vertex-to-index-mapping between old and new indices
      std::unordered_map <index_t, index_t> interface_global_to_local_index_map_vertices;
      std::unordered_map <index_t, index_t> interface_local_to_global_index_map_vertices;

      if (num_points == 0 || num_cells == 0)
      {
         pragmatic_interfaces.push_back(nullptr);
      
         //if interface i does not exist, insert negative numbers into the mapping vectors to mark them as empty
         interface_global_to_local_index_map_vertices.insert( std::make_pair(-1, -1) );
         g2l_vertex_map_interfaces[interface_counter] = interface_global_to_local_index_map_vertices;
         l2g_vertex_map_interfaces[interface_counter] = interface_global_to_local_index_map_vertices;

        continue;
      }

      index_t new_vertex_id = 0;  
      for (auto it : interface_nodes[i][j])
      {
        interface_global_to_local_index_map_vertices.insert( std::make_pair(it, new_vertex_id++) );
        ++vertex_appearances[it]; //TODO: atomic increment for all appearances incrementations for parallelization???
      }

      g2l_vertex_map_interfaces[interface_counter] = interface_global_to_local_index_map_vertices;

      //and get also the index-to-vertex mapping (opposite direction than vertex to index mapping)
      //TODO: put this in the for loop above
      for (auto it : interface_global_to_local_index_map_vertices)
      {
        interface_local_to_global_index_map_vertices[it.second] = it.first;
      }

      l2g_vertex_map_interfaces[interface_counter] = interface_local_to_global_index_map_vertices;

      //pre-allocate memory
      x_coords_interfaces[interface_counter].reserve(num_points);
      y_coords_interfaces[interface_counter].reserve(num_points);
      ENLists_interfaces[interface_counter].resize(3*num_cells);

      //get coordinates of each vertex
      //TODO: combine this loop with the next one!
      size_t counter = 0;
      for (auto it : interface_nodes[i][j])
      {
        double p[2];
        mesh->get_coords( it, p);
        x_coords_interfaces[interface_counter][counter] = p[0];
        y_coords_interfaces[interface_counter][counter] = p[1];           
        ++counter;
      }

      //create ENList with respect to the new vertex indices and get also to mapping of indices for the elements      
      std::unordered_map <index_t, index_t> global_to_local_index_map_elements;
      std::unordered_map <index_t, index_t> local_to_global_index_map_elements;

      counter=0;        
      for (auto it : interface_elements_sets[i][j])
      {         
        global_to_local_index_map_elements.insert(std::make_pair(it, counter/3));
        local_to_global_index_map_elements.insert(std::make_pair(counter/3, it));

        const index_t *element_ptr = nullptr;
        element_ptr = mesh->get_element(it);
          
        ENLists_interfaces[interface_counter][counter++] = interface_global_to_local_index_map_vertices[*(element_ptr++)];
        ENLists_interfaces[interface_counter][counter++] = interface_global_to_local_index_map_vertices[*(element_ptr++)];
        ENLists_interfaces[interface_counter][counter++] = interface_global_to_local_index_map_vertices[*(element_ptr++)];      

        ++element_appearances[it];    
      }

      l2g_element_map_interfaces[interface_counter] = local_to_global_index_map_elements;
      g2l_element_map_interfaces[interface_counter] = global_to_local_index_map_elements;
          
      //create pragmatic mesh 
      Mesh<double> *interface_mesh = nullptr;

      //TODO: change for 3D refinement
      //mesh = new Mesh<double> ( num_points, num_cells, &(ENLists_regions[region.id()][0]) ,&(x_coords[region.id()][0]), &(y_coords[region.id()][0]), &(z_coords[region.id()][0]) );        
      interface_mesh = new Mesh<double> ( num_points, num_cells, &(ENLists_interfaces[interface_counter][0]), &(x_coords_interfaces[interface_counter][0]), &(y_coords_interfaces[interface_counter][0]) );
      interface_mesh->create_boundary();
      
      //create metric
      MetricField<double,2> metric_field(*interface_mesh);
      std::vector<double> metric_vector;
      double eta=0.0001;

      //read data
      ifstream metric_data;
      std::string metric_file = "metric_data_inter";
      metric_file += std::to_string(i);
      metric_file += ".txt";
      std::cout << metric_file <<std::endl;
      metric_data.open(metric_file.c_str());

      double metric_value;
      int ctr = 0;
      //for (size_t i = 0; i < num_points; ++i)
      while(metric_data >> metric_value)
      {
          metric_vector.push_back(metric_value);
          ++ctr;
      }

      for(size_t j = 0; j < num_points; j++) 
      {/*
        double x = 2*mesh->get_coords(j)[0]-1;
        double y = 2*mesh->get_coords(j)[1]-1;
        
        psi[j] = 0.100000000000000*sin(50*x) + atan2(-0.100000000000000, (double)(2*x - sin(5*y)));*/
        /*
        double x = interface_mesh->get_coords(i)[0];
        double y = interface_mesh->get_coords(i)[1];
        //std::cerr << "interface calc psi " << i << std::endl;
        psi[j]= 1;
        //std::cerr << psi[i] << std::endl;*/

        double m[3];
        m[0] = metric_vector[3*j];
        m[1] = metric_vector[3*j+1];
        m[2] = metric_vector[3*j+2];
        //std::cerr << "setting field" << std::endl;
        metric_field.set_metric(m, j);
      }
     
      //metric_field.add_field(&(psi[0]), eta, 1);
      metric_field.update_mesh();
      
      pragmatic_interfaces.push_back(interface_mesh);
    }
  }
  //end of loop over all interfaces
}
//end of void GroupedPartitions::CreatePragmaticDataStructures()

//GroupedPartitions::WritePartitionsAndInterfaces()
void GroupedPartitions::WritePartitionsAndInterfaces(std::string name)
{
  std::cout << "Using filename \033[1;36m" << name << "\033[0m" << std::endl;

  //write partitions
  for (size_t i = 0; i < pragmatic_partitions.size(); ++i)
  {
    std::cout << "Writing partition " << i << std::endl;
    std::cout << "  Vertex count = " << pragmatic_partitions[i]->get_number_nodes() << std::endl;
    std::cout << "  Cell count = " << pragmatic_partitions[i]->get_number_elements() << std::endl;
    
    std::string filename;
    filename += "examples/data/pragmatic_metis_partitioning/";
    filename += name;
    filename += "_partition_";
    filename += std::to_string( i );
  
    VTKTools<double>::export_vtu(filename.c_str(), pragmatic_partitions[i]);
  }
  //end of write partitions

  //write interfaces
  for (size_t i = 0; i < pragmatic_interfaces.size(); ++i)
  {
    if (pragmatic_interfaces[i] == nullptr)
    {
      continue;
    }

    std::cout << "Writing interface " << i << std::endl;
    std::cout << "  Vertex count = " << pragmatic_interfaces[i]->get_number_nodes() << std::endl;
    std::cout << "  Cell count = " << pragmatic_interfaces[i]->get_number_elements() << std::endl;

    std::string filename;
    filename += "examples/data/pragmatic_metis_partitioning/";
    filename += name;
    filename += "_interface_";
    filename += std::to_string( i );
  
    VTKTools<double>::export_vtu(filename.c_str(), pragmatic_interfaces[i]);
  }        
  //end of write interfaces
}
//end of GroupedPartitions::WritePartitionsAndInterfaces()

//GroupedPartitions::GetPartitionAndInterfaceBoundaries()
//TODO: replace function calls with lookups of class variables!
void GroupedPartitions::GetPartitionAndInterfaceBoundaries()
{
  std::cout << "Creating Boundary Information" << std::endl;
  //get global boundary nodes
  std::vector<int> global_boundary = mesh->copy_boundary_vector();

  for (size_t i = 0; i < mesh->get_number_elements(); ++i)
  {
    const int *n = mesh->get_element(i);

    for (size_t j = 0; j < 3; ++j)
    {
      boundary_nodes_mesh[n[(j+1)%3]] = std::max(boundary_nodes_mesh[n[(j+1)%3]], global_boundary[i*3+j]);
      boundary_nodes_mesh[n[(j+2)%3]] = std::max(boundary_nodes_mesh[n[(j+2)%3]], global_boundary[i*3+j]);
    }
  }
  //end of get global boundary nodes

  //get partition boundary nodes
  for (size_t mesh_number = 0; mesh_number < pragmatic_partitions.size(); mesh_number++)
  {
    std::vector<int> boundary_nodes_tmp(pragmatic_partitions[mesh_number]->get_number_nodes(), 0);
    std::vector<int> boundary_mesh = pragmatic_partitions[mesh_number]->copy_boundary_vector();

    for (size_t i = 0; i < pragmatic_partitions[mesh_number]->get_number_elements(); ++i)
    {
      const int *n = pragmatic_partitions[mesh_number]->get_element(i);

      for (size_t j = 0; j < 3; ++j)
      {
        boundary_nodes_tmp[n[(j+1)%3]] = std::max(boundary_nodes_tmp[n[(j+1)%3]], boundary_mesh[i*3+j]);
        boundary_nodes_tmp[n[(j+2)%3]] = std::max(boundary_nodes_tmp[n[(j+2)%3]], boundary_mesh[i*3+j]);
      }
    }
    boundary_nodes_partitions[mesh_number]=boundary_nodes_tmp;
  }
  //end of get partition boundary nodes

  //get interface boundary nodes
  for (size_t mesh_number = 0; mesh_number < pragmatic_interfaces.size(); mesh_number++)
  {
    if  (pragmatic_interfaces[mesh_number] == nullptr)
    {  
      continue;
    }
    std::vector<int> boundary_nodes_tmp(pragmatic_interfaces[mesh_number]->get_number_nodes(), 0);
    std::vector<int> boundary_mesh = pragmatic_interfaces[mesh_number]->copy_boundary_vector();

    for (size_t i = 0; i < pragmatic_interfaces[mesh_number]->get_number_elements(); ++i)
    {
      const int *n = pragmatic_interfaces[mesh_number]->get_element(i);

      for (size_t j = 0; j < 3; ++j)
      {
        boundary_nodes_tmp[n[(j+1)%3]] = std::max(boundary_nodes_tmp[n[(j+1)%3]], boundary_mesh[i*3+j]);
        boundary_nodes_tmp[n[(j+2)%3]] = std::max(boundary_nodes_tmp[n[(j+2)%3]], boundary_mesh[i*3+j]);
      }
    }
    boundary_nodes_interfaces[mesh_number]=boundary_nodes_tmp;
  }
  //end of get partition boundary nodes
} 
//end of GroupedPartitions::GetPartitionAndInterfaceBoundaries()

//GroupedPartitions::GetCoords()
void GroupedPartitions::GetCoords(index_t n, double *vertex)
{
  bool found = false;

  //search partitions
 
  //do a lookup in the g2l_vertex_map_partitions and then return the corresponding coordinates!
  for (size_t i = 0; i < g2l_vertex_map_partitions.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);
    
    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }

    //std::cout << "found vertex " << n << " in partition " << i << std::endl;
    pragmatic_partitions[i]->get_coords(position->second, vertex);
    found = true;
    break;
  }
  //end of seach partitions

  //check interface
  if (found)
    return; 

  for (size_t i = 0; i < g2l_vertex_map_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);

    if(position == g2l_vertex_map_interfaces[i].end())
    {        
        continue;
    }

    //std::cout << "found vertex " << n << " in interface " << i << std::endl;
    pragmatic_interfaces[i]->get_coords(position->second, vertex);
    break;
  }
  //end of search interface
}
//end of GroupedPartitions::GetCoords()

//GroupedPartitions::SetCoords()
//TODO: implementation, and validation check, if the moved vertex is on a valid location!!!
//TODO: implement for 3D
void GroupedPartitions::SetCoords(index_t n, double *vertex)
{
  int to_do = vertex_appearances[n];
  int comitted = 0;

  //update partitions
  for(size_t i = 0; i < g2l_vertex_map_partitions.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);

    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }

      pragmatic_partitions[i]->set_coords(position->second, vertex);
      ++comitted;

      if (to_do == comitted)
        break;
    
  }
  //end of update partitions

  if (to_do == comitted)
    return;

  //update interfaces
  for(size_t i = 0; i < g2l_vertex_map_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);
    
    if(position == g2l_vertex_map_interfaces[i].end())
    {        
        continue;
    }

    else
    {
      pragmatic_interfaces[i]->set_coords(position->second, vertex);
      ++comitted;

      if (to_do == comitted)
        break;
    }
  }
  
  //end of update interfaces
}
//end of GroupedPartitions:SetCoords()

//GroupedPartitions::GetNNList(index_t n, int *partition, int *interfaces)
//Uses global indices
std::vector<index_t>* GroupedPartitions::GetNNList(index_t n, int *partitions, int *interfaces)
{
  bool partition = false;
  bool interface = false;
  int mapping = GetVertexPartitionOrInterface(n, &partition, &interface);

  if (partition)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[mapping].find(n);

    *partitions = mapping;
    return pragmatic_partitions[mapping]->get_nnlist( position->second );
  }

  else if (interface)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[mapping].find(n);

    *interfaces = mapping;
    return pragmatic_interfaces[mapping]->get_nnlist( position->second );
  }
  /*
  //check partitions
  for (size_t i = 0; i < nparts; ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);

    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }

    *partitions = i;
    return pragmatic_partitions[i]->get_nnlist( position->second );
  }
  //end of check partitions

  //check interfaces
  for (size_t i = 0; i < g2l_vertex_map_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);

    if (position == g2l_vertex_map_interfaces[i].end())
    {
      continue;
    }

    *interfaces = i;
    return pragmatic_interfaces[i]->get_nnlist( position->second );
  }  
  //end of check interfaces
  */
  else
    std::cout << "\033[1;31mVertex " << n << " not found\033[0m" << std::endl;
}
//end of GroupedPartitions::GetNNList(index_t n, int *partitions, int *interfaces)

//GroupedPartitions::GetNNList_adv(index_t n, int *partition, int *interfaces)
//this is the GetNNList-function which returns the correct NNList if vertex is on boundary
bool GroupedPartitions::GetNNList_adv(index_t n, std::vector<index_t> &_NNList)
{
  int partition = -1;
  int interface = -1;
  if (vertex_appearances[n] == 2)
  {
    GetVertexPartitionAndInterface(n, &partition, &interface);
  }   

  if (partition == -1 || interface == -1)
    return false;

  //get NNLists with local indices
  std::vector<index_t> *NNList_partition = pragmatic_partitions[partition]->get_nnlist( g2l_vertex_map_partitions[partition].at(n) );
  std::vector<index_t> *NNList_interface = pragmatic_interfaces[interface]->get_nnlist( g2l_vertex_map_interfaces[interface].at(n) );
 
  //create NNLists with global indices
  std::vector<index_t> global_nnl_part( NNList_partition->size() );
  std::vector<index_t> global_nnl_inter( NNList_interface->size() );

  for (size_t i = 0; i < NNList_partition->size(); ++i)
  {
    global_nnl_part[i] = l2g_vertex_map_partitions[partition].at( NNList_partition->at(i) );
  }
  
  for (size_t i = 0; i < NNList_interface->size(); ++i)
  {
    global_nnl_inter[i] = l2g_vertex_map_interfaces[interface].at( NNList_interface->at(i) );
  }

  //create union of the two vectors with global indices and return the union vector
  std::set_union( global_nnl_part.begin(), global_nnl_part.end(),
                  global_nnl_inter.begin(), global_nnl_inter.end(),
                  inserter(_NNList, _NNList.begin()));
  
  //_NNList.resize(it - _NNList.begin());

  //std::cout << n << " in partition: " << partition << " and interface: " << interface << std::endl;
  return true;
}
//end of GroupedPartitions::GetNNList_adv(index_t n, int *partitions, int *interfaces)

//GroupedPartitions::GetNEList(index_t n, int *partition, int *interfaces)
//uses global indices!
//std::set<index_t> GroupedPartitions::GetNEList(index_t n, std::set<index_t> &_NEList, int *partitions, int *interfaces)
bool GroupedPartitions::GetNEList(index_t n, std::set<index_t> &_NEList, int *partitions, int *interfaces)
{
  bool partition;
  bool interface;
  int mapping = GetVertexPartitionOrInterface(n, &partition, &interface);
  //std::cerr << "found mapping" << std::endl;

  if (partition)
  {
    //std::cerr << "search in partition 0 " << std::endl;
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[mapping].find(n);
    //std::cerr << "  " << position->first << " " << position->second << " " << n << std::endl;
    *partitions = mapping;    

    _NEList = pragmatic_partitions[mapping]->get_reference_nelist( position->second );
    //std::cerr << "got nelist" << std::endl;
    //_NEList = pragmatic_partitions[mapping]->get_nelist( g2l_element_map_partitions[mapping].at(n) );

    return true;
  }

  else if (interface)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[mapping].find(n);

    double coords[2];
    mesh->get_coords(n, coords);
  
    std::vector<std::set<int>> NE_global =  mesh->get_node_element();
  
    *interfaces = mapping;
/*
    for (auto it : NE_global[n])
    {
      std::cout << "  " << it << std::endl;
      std::cout << "    " << g2l_element_map_interfaces[mapping].at(it) << std::endl;
    }

    std::cout<<std::endl;
/*
    for (auto it : pragmatic_interfaces[mapping]->get_reference_nelist( position->second ))
    {
      std::cout << "  " << it << std::endl;
    }
*/
    //_NEList = pragmatic_interfaces[mapping]->get_nelist(position->second);
    _NEList = pragmatic_interfaces[mapping]->get_reference_nelist( position->second );
    //return pragmatic_interfaces[mapping]->get_reference_nelist( position->second );
    return true;
  }

  return false;
}
//end of GroupedPartitions::GetNEList(index_t n, int *partitions, int *interfaces)

//GroupedPartitions::GetNEList_adv(index_t n, int *partition, int *interfaces)
//this is the GetNEList-function which returns the correct NEList if vertex is on boundary
bool GroupedPartitions::GetNEList_adv(index_t n, std::set<index_t> &_NEList)
{
  int partition = -1;
  int interface = -1;
    
  if (vertex_appearances[n] == 2)
  {
    GetVertexPartitionAndInterface(n, &partition, &interface);
  }   

  if (partition == -1 || interface == -1)
    return false;

  //get NELists with local indices 
  std::set<index_t> NEList_partition = pragmatic_partitions[partition]->get_reference_nelist( g2l_vertex_map_partitions[partition].at(n) );
  std::set<index_t> NEList_interface = pragmatic_interfaces[interface]->get_reference_nelist( g2l_vertex_map_interfaces[interface].at(n) );

  //std::vector<std::set<int>> NE_global =  pragmatic_interfaces[interface]->get_node_element();

  //create NELists with global indices
  for(auto it : NEList_partition)
  {
    _NEList.insert(l2g_element_map_partitions[partition].at(it));
  }

  for(auto it : NEList_interface)
  {    
    //std::cout << "      " << it << " " << l2g_element_map_interfaces[interface].at(it) << std::endl;
    _NEList.insert(l2g_element_map_interfaces[interface].at(it));
    continue;
  }

  //return true if everything was successful
  return true;
}
//end of GroupedPartitions::GetNEList_adv(index_t n, int *partitions, int *interfaces)

//GroupedPartitions::WriteMergedMesh()
void GroupedPartitions::WriteMergedMesh(std::string filename)
{
  std::cout << "Write merged mesh" << std::endl;
/*
  //boost timer for collecting benchmark data
  boost::timer::cpu_timer t;
*/
  //create merged ENList
  //TODO: Update for 3D case (then num_elements has to be multiplied with 4)
  std::vector<index_t> merged_ENList(num_elements*3);

  //ierate over partitions
  int global_element_counter = 0;
  for (size_t i = 0; i < pragmatic_partitions.size(); ++i)  
  {
    for (size_t j = 0; j < pragmatic_partitions[i]->get_number_elements(); ++j)
    {
      const index_t *element_ptr = nullptr;
      element_ptr = pragmatic_partitions[i]->get_element(j);

      merged_ENList[3*global_element_counter] = l2g_vertex_map_partitions[i].at( *(element_ptr++) );
      merged_ENList[3*global_element_counter+1] = l2g_vertex_map_partitions[i].at( *(element_ptr++) );
      merged_ENList[3*global_element_counter+2] = l2g_vertex_map_partitions[i].at( *(element_ptr++) );  
      ++global_element_counter;
    }
  }

  //iterate over interfaces
  for (size_t i = 0; i < pragmatic_interfaces.size(); ++i)  
  {
    if  (pragmatic_interfaces[i] == nullptr)
    {  
      continue;
    }

    for (size_t j = 0; j < pragmatic_interfaces[i]->get_number_elements(); ++j)
    {
      const index_t *element_ptr = nullptr;
      element_ptr = pragmatic_interfaces[i]->get_element(j);

      merged_ENList[3*global_element_counter] = l2g_vertex_map_interfaces[i].at( *(element_ptr++) );
      merged_ENList[3*global_element_counter+1] = l2g_vertex_map_interfaces[i].at( *(element_ptr++) );
      merged_ENList[3*global_element_counter+2] = l2g_vertex_map_interfaces[i].at( *(element_ptr++) );  
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

  writer << "  <Piece NumberOfPoints=\"" << num_nodes << "\" NumberOfCells=\"" << num_elements << "\">" << std::endl;    

  //write points into file
  writer << "   <Points>" << std::endl;
  writer << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;  

  //iterate over coordinates vector
  for (size_t i = 0; i < num_nodes; ++i)
  {
    //TODO: change for 3D case
    //double x[3];
    double x[2];
    GetCoords(i, x);
    writer << std::fixed << std::setprecision(8) << x[0] << " " << x[1] << " " << "0" << std::endl;      
    //TODO: change for 3D case
    //writer << x[0] << " " << x[1] << " " << x[2] << std::endl;
  }
  //writer << std::endl;
  writer << "    </DataArray>" << std::endl;
  writer << "   </Points> " << std::endl;
  //end of write points into file

  //write cells into file
  writer << "   <Cells> " << std::endl;
  writer << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;

  for (size_t i = 0; i < num_elements; ++i)
  {
    writer << merged_ENList[3*i] << " " << merged_ENList[3*i+1] << " " << merged_ENList[3*i+2] << " ";

    writer << std::endl;
  }

  writer << std::endl;
  writer << "    </DataArray>" << std::endl;
  //end of write cells into file

  //write offset into file
  writer << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    
  for (size_t i = 1; i <= num_elements; ++i)
  {
    //TODO: change for 3D case
    //writer << j*4 << " ";
    writer << i*3 << " ";

    if (i%6 == 0)
    {
      writer << std::endl;
    }
  }

  writer << std::endl;
  writer << "    </DataArray>" << std::endl;
  //end of write offset into file
  
  //write types into file
  writer << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  
  for (size_t j = 0; j < num_elements; ++j) 
  {
    //TODO: change for 3D case
    //writer << 10 << " ";
    writer << 5 << " ";
      
    if (j%6 == 0)
    {
      writer << std::endl;
    }
  }

  writer << std::endl;
  writer << "    </DataArray>" << std::endl;
  //end of write types into file

  writer << "   </Cells>" << std::endl;
  writer << "  </Piece>" << std::endl;

  //write footer
  writer << " </UnstructuredGrid>" << std::endl;
  writer << "</VTKFile>" << std::endl;
  //end of write footer
  
  //close ofstream object
  writer.close();
/*
  boost::timer::cpu_times timer = t.elapsed();

  //runtime benchmark
  ofstream runtime;
  runtime.open("mesh_merging_runtime.txt");
  runtime << "    wall    user    system" << std::endl; 
  runtime << std::fixed << std::setprecision(6) << static_cast<double>(timer.wall) / 1000000000.0 << " " << static_cast<double>(timer.user) / 1000000000.0 << " " << static_cast<double>(timer.system) / 1000000000.0 << std::endl;
  runtime.close();
*/
}
//end of GroupedPartitions::WriteMergedMesh()

//TODO: GroupedPartitions::GetVertexPartitionOrInterface(index_t n, bool *partition, bool *interface)
// this function returns only if a vertex is found in a partition or in an interface
int GroupedPartitions::GetVertexPartitionOrInterface(index_t n, bool *partition, bool *interface)
{
  //check partitions
  for (size_t i = 0; i < nparts; ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);

    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }
    
    *partition = true;
    *interface = false;
    return i;
  }
  //end of check partitions

  //check interfaces
  for (size_t i = 0; i < g2l_vertex_map_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);

    if (position == g2l_vertex_map_interfaces[i].end())
    {
      continue;
    }
    *partition = false;
    *interface = true;
    return i;
  }  
  //end of check interfaces

  //error
  return -1;
}
//end of GroupedPartitions::GetVertexPartitionOrInterface(index_t n, bool *partition, bool *interface)

//GroupedPartitions::GetVertexPartitionAndInterface(index_t n)
// this funtion returns partitions and interfaces for boundary vertex n (global index)
// NOTE: vertex has to be part of any boundary (partition or interface)
void GroupedPartitions::GetVertexPartitionAndInterface(index_t n, int *partition, int *interface)
{
  //check partition
  for (size_t i = 0; i < pragmatic_partitions.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);

    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }
    *partition = i;
    break;
  } //end of check partitions

  //check interfaces
  for (size_t i = 0; i < pragmatic_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);

    if (position == g2l_vertex_map_interfaces[i].end())
    {
      continue;
    }
    
    *interface=i;
    break;
  } //end of check interfaces
}
//end of GroupedPartitions::GetVertexPartitionAndInterface(index_t n)

//TODO: GroupedPartitions::GetElementPartitionOrInterface(index_t n)
int GroupedPartitions::GetElementPartitionOrInterface(index_t n)
{
  return -1;
}
//end of GroupedPartitions::GetElementPartitionOrInterface(index_t n)

void GroupedPartitions::PrintQuality()
{
  std::vector<double> element_area (num_elements);
  std::vector<index_t> tmp_ENList = mesh->get_element_node();

  std::cout << "Print Quality" << std::endl;  

  ofstream quality;
  quality.open("quality.txt");

  //for loop iterating over all elements in the global mesh
  for (size_t i = 0; i < num_elements; ++i)
  {
    //copy paste from Mesh.h - calculate_area()-function
    const double *x1 = mesh->get_coords(tmp_ENList[3*i]);
    const double *x2 = mesh->get_coords(tmp_ENList[3*i+1]);
    const double *x3 = mesh->get_coords(tmp_ENList[3*i+2]);

    // Use Heron's Formula
    long double a;
    {
        long double dx = ((long double)x1[0]-(long double)x2[0]);
        long double dy = ((long double)x1[1]-(long double)x2[1]);
        a = std::sqrt(dx*dx+dy*dy);
    }
    long double b;
    {
        long double dx = ((long double)x1[0]-(long double)x3[0]);
        long double dy = ((long double)x1[1]-(long double)x3[1]);
        b = std::sqrt(dx*dx+dy*dy);
    }
    long double c;
    {
        long double dx = ((long double)x2[0]-(long double)x3[0]);
        long double dy = ((long double)x2[1]-(long double)x3[1]);
        c = std::sqrt(dx*dx+dy*dy);
    }
    long double s = (a+b+c)/2;

    quality << i << ": " << std::sqrt(s*(s-a)*(s-b)*(s-c)) << std::endl;;

  } //end of for loop iterating over all elements in the global mesh
  
  quality.close();
}

//GroupedPartitions::GetAllVertexPartitionsAndInterfaces
void GroupedPartitions::GetAllVertexPartitionsAndInterfaces(index_t n, std::vector<bool>& partitions, std::vector<bool>& interfaces)
{
  //check partitions
  for (size_t i = 0; i < pragmatic_partitions.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);

    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }
    partitions[i] = true;
  } //end of check partitions

  //check interfaces
  for (size_t i = 0; i < pragmatic_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);

    if (position == g2l_vertex_map_interfaces[i].end())
    {
      continue;
    }
    
    interfaces[i]=true;
  } //end of check interfaces
}
//end of GroupedPartitions::GetAllVertexPartitionsAndInterfaces

//GroupedPartitions::GetNEList_part
void GroupedPartitions::GetNEList_part(index_t n, std::set<int>& _NEList, int partition)
{
  std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[partition].find(n); 

  _NEList = pragmatic_partitions[partition]->get_reference_nelist( position->second );
  //_NEList = pragmatic_partitions[mapping]->get_nelist( g2l_element_map_partitions[mapping].at(n) );
}
//end of GroupedPartitions::GetNEList_part

//GroupedPartitions::GetNEList_inter
void GroupedPartitions::GetNEList_inter(index_t n, std::set<int>& _NEList, int interface)
{
  std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[interface].find(n); 

  _NEList = pragmatic_interfaces[interface]->get_reference_nelist( position->second );
}
//end of GroupedPartitions::GetNEList_inter

//GroupedPartitions::GetCoords_part()
inline void GroupedPartitions::GetCoords_part(index_t n, double *vertex, int partition)
{
  pragmatic_partitions[partition]->get_coords(n, vertex);
} //end of GroupedPartitions::GetCoords_part()

//GroupedPartitions::GetCoords_inter()
inline void GroupedPartitions::GetCoords_inter(index_t n, double *vertex, int interface)
{
  pragmatic_interfaces[interface]->get_coords(n, vertex);
} //end of GroupedPartitions::GetCoords_inter

//GroupedPartitions::SetCoords_part
inline void GroupedPartitions::SetCoords_part(index_t n, double *vertex, int partition)
{
  pragmatic_partitions[partition]->set_coords(n, vertex);;
} //end of GroupedPartitions::SetCoords_part

//GroupedPartitions::SetCoords_inter
inline void GroupedPartitions::SetCoords_inter(index_t n, double *vertex, int interface)
{
  pragmatic_interfaces[interface]->set_coords(n, vertex);;
} //end of GroupedPartitions::SetCoords_inter

//GroupedPartitions::GetMetric()
void GroupedPartitions::GetMetric(index_t n, double* metric)
{
  bool found = false;

  //search partitions
 
  //do a lookup in the g2l_vertex_map_partitions and then return the corresponding coordinates!
  for (size_t i = 0; i < g2l_vertex_map_partitions.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_partitions[i].find(n);
    
    if (position == g2l_vertex_map_partitions[i].end())
    {
      continue;
    }

    //std::cout << "found vertex " << n << " in partition " << i << std::endl;
    pragmatic_partitions[i]->get_metric(position->second, metric);
  
    found = true;
    break;
  }
  //end of seach partitions

  //check interface
  if (found)
    return; 

  for (size_t i = 0; i < g2l_vertex_map_interfaces.size(); ++i)
  {
    std::unordered_map<index_t, index_t>::iterator position = g2l_vertex_map_interfaces[i].find(n);

    if(position == g2l_vertex_map_interfaces[i].end())
    {        
        continue;
    }

    //std::cout << "found vertex " << n << " in interface " << i << std::endl;
    pragmatic_interfaces[i]->get_metric(position->second, metric);
    break;
  }
  //end of search interface
}
//end of GroupedPartitions::GetMetric()

//GroupedPartitions::AddCoords()
void GroupedPartitions::AddCoords(size_t threadIdx, std::vector<double> coords, size_t dim, size_t splitCnt)
{
  pragmatic_partitions[0]->add_coords(threadIdx, coords, dim, splitCnt);
}
//end of GroupedPartitions::AddCoords()

//GroupedPartitions::AddCoords_part()
void GroupedPartitions::AddCoords_part(size_t threadIdx, std::vector<double> coords, size_t dim, size_t splitCnt, size_t part)
{
  pragmatic_partitions[part]->add_coords(threadIdx, coords, dim, splitCnt);
}
//end of GroupedPartitions::AddCoords_part()

//GroupedPartitions::AddMetric()
void GroupedPartitions::AddMetric(size_t threadIdx, std::vector<double> metric, size_t msize, size_t splitCnt)
{
  pragmatic_partitions[0]->add_metric(threadIdx, metric, msize, splitCnt);
}
//end of GroupedPartitions::AddMetric()

//GroupedPartitions::AddMetric_part()
void GroupedPartitions::AddMetric_part(size_t threadIdx, std::vector<double> metric, size_t msize, size_t splitCnt, size_t part)
{
  pragmatic_partitions[part]->add_metric(threadIdx, metric, msize, splitCnt);
}
//end of GroupedPartitions::AddMetric_part()

//GroupedPartitions::AddNNList()
void GroupedPartitions::AddNNList(size_t id, size_t add_id)
{
  //std::cerr << "adding " << add_id << " to " << id << std::endl;
  pragmatic_partitions[0]->add_nnlist(id , add_id);
}
//end of GroupedPartitions::AddNNList()

//GroupedPartitions::AddNNList_part()
void GroupedPartitions::AddNNList_part(size_t id, size_t add_id, int part)
{
  //std::cerr << "adding " << add_id << " to " << id << std::endl;
  pragmatic_partitions[part]->add_nnlist(id , add_id);
}
//end of GroupedPartitions::AddNNList_part()

//GroupedPartitions::RemNNList()
void GroupedPartitions::RemNNList(size_t id, size_t rem_id)
{
  //std::cerr << "removing " << rem_id << " from " << id << std::endl;
  pragmatic_partitions[0]->remove_nnlist(id, rem_id);
}
//end of GroupedPartitions::RemNNList()

//GroupedPartitions::RemNNList_part()
void GroupedPartitions::RemNNList_part(size_t id, size_t rem_id, int part)
{
  //std::cerr << "removing " << rem_id << " from " << id << std::endl;
  pragmatic_partitions[part]->remove_nnlist(id, rem_id);
}
//end of GroupedPartitions::RemNNList_part()

//GroupedPartitions::AddNEList()
void GroupedPartitions::AddNEList(size_t id, size_t add_id)
{
  //std::cerr << "adding " << add_id << " to " << id << std::endl;
  pragmatic_partitions[0]->add_nelist(id , add_id);
}
//end of GroupedPartitions::AddNEList()

//GroupedPartitions::AddNEList_part()
void GroupedPartitions::AddNEList_part(size_t id, size_t add_id, int part)
{
  //std::cerr << "adding " << add_id << " to " << id << std::endl;
  pragmatic_partitions[part]->add_nelist(id , add_id);
}
//end of GroupedPartitions::AddNEList_part()

//GroupedPartitions::AddNEList_fix()
void GroupedPartitions::AddNEList_fix(size_t id, size_t add_id, size_t threadIdx)
{
  //std::cerr << "adding " << add_id << " to " << id << std::endl;
  pragmatic_partitions[0]->add_nelist_fix(id , add_id, threadIdx);
}
//end of GroupedPartitions::AddNEList_fix()

//GroupedPartitions::AddNEList_fix_part()
void GroupedPartitions::AddNEList_fix_part(size_t id, size_t add_id, size_t threadIdx, int part)
{
  //std::cerr << "adding " << add_id << " to " << id << std::endl;
  pragmatic_partitions[part]->add_nelist_fix(id , add_id, threadIdx);
}
//end of GroupedPartitions::AddNEList_fix_part()

//GroupedPartitions::RemNEList()
void GroupedPartitions::RemNEList(size_t id, size_t rem_id)
{
  //std::cerr << "removing " << rem_id << " from " << id << std::endl;
  pragmatic_partitions[0]->remove_nelist(id, rem_id);
}
//end of GroupedPartitions::RemNEList()

//GroupedPartitions::RemNEList_part()
void GroupedPartitions::RemNEList_part(size_t id, size_t rem_id, int part)
{
  //std::cerr << "removing " << rem_id << " from " << id << std::endl;
  pragmatic_partitions[part]->remove_nelist(id, rem_id);
}
//end of GroupedPartitions::RemNEList_part()

//GroupedPartitions::calculate_quality(const index_t* n)
double GroupedPartitions::calculate_quality(const index_t* n)
{

  //works only for test case with partition0
  const double *x0 = pragmatic_partitions[0]->get_coords(n[0]);
  const double *x1 = pragmatic_partitions[0]->get_coords(n[1]);
  const double *x2 = pragmatic_partitions[0]->get_coords(n[2]);

  const double *m0 = pragmatic_partitions[0]->get_metric(n[0]);
  const double *m1 = pragmatic_partitions[0]->get_metric(n[1]);
  const double *m2 = pragmatic_partitions[0]->get_metric(n[2]);

  ElementProperty<double>* property;
  pragmatic_partitions[0]->get_property(property);

  return property->lipnikov(x0, x1, x2, m0, m1, m2);
}
//end of GroupedPartitions::calculate_quality(const index_t* n)

//GroupedPartitions::calculate_quality_part(const index_t* n, int part)
double GroupedPartitions::calculate_quality_part(const index_t* n, int part)
{

  //works only for test case with partition0
  const double *x0 = pragmatic_partitions[part]->get_coords(n[0]);
  const double *x1 = pragmatic_partitions[part]->get_coords(n[1]);
  const double *x2 = pragmatic_partitions[part]->get_coords(n[2]);

  const double *m0 = pragmatic_partitions[part]->get_metric(n[0]);
  const double *m1 = pragmatic_partitions[part]->get_metric(n[1]);
  const double *m2 = pragmatic_partitions[part]->get_metric(n[2]);

  ElementProperty<double>* property;
  pragmatic_partitions[part]->get_property(property);

  return property->lipnikov(x0, x1, x2, m0, m1, m2);
}
//end of GroupedPartitions::calculate_quality_part(const index_t* n, int part)

//GroupedPartitions::update_quality(const index_t* element)
double GroupedPartitions::update_quality(index_t element)
{
  const index_t* n = pragmatic_partitions[0]->get_element(element);

  pragmatic_partitions[0]->set_quality(element, calculate_quality(n));
}
//end of GroupedPartitions::update_quality(const index_t* n)

//GroupedPartitions::update_quality_part(const index_t* element, int part)
double GroupedPartitions::update_quality_part(index_t element, int part)
{
  const index_t* n = pragmatic_partitions[part]->get_element(element);

  pragmatic_partitions[part]->set_quality(element, calculate_quality_part(n, part));
}
//end of GroupedPartitions::update_quality_part(const index_t* n, int part)

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif 
