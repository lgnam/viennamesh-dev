#ifndef GROUPED_PARTITIONS_SMOOTH_HPP
#define GROUPED_PARTITIONS_SMOOTH_HPP

#include "grouped_partitions.hpp"

#include <time.h>

#include <omp.h>

#include "dummy.hpp"

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Declaration                                                                   //
//----------------------------------------------------------------------------------------------------------------------------------------------//

class GroupedPartitionsSmooth
{
  public:
    //constructor and destructor
    GroupedPartitionsSmooth(GroupedPartitions &GP);
    ~GroupedPartitionsSmooth();

    //member functions
    bool SimpleLaplace(int iterations);
    bool SimpleLaplaceOnGroups_sections(int iterations);
    bool SimpleLaplaceOnGroups_tasks(int iterations);
    void Evaluate();

  private:
    GroupedPartitions &mesh;

    void LaplaceKernel(int part1, int part2, int inter);
    void AssignPartitionsAndInterface();

    std::vector<std::vector<int>> work_sets;  //vector storing the work sets (value 1 and 2 denote the partitions, value 3 the corresponding interface, parameter is num_threads)
    std::vector<int> num_touches_partitions;
    std::vector<int> num_touches_interfaces;

    std::vector<bool> partition_assigned; //vectors storing if partition i is already assigned to a thread
    std::vector<bool> interface_assigned;

};

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Helper Functions                                                              //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                Implementation                                                                //
//----------------------------------------------------------------------------------------------------------------------------------------------//

//Constructor
GroupedPartitionsSmooth::GroupedPartitionsSmooth(GroupedPartitions &GP) : mesh(GP), num_touches_partitions(mesh.pragmatic_partitions.size(), 0), num_touches_interfaces(mesh.pragmatic_interfaces.size(), 0), work_sets(mesh.nparts/2), partition_assigned(mesh.nparts, false), interface_assigned(mesh.pragmatic_interfaces.size(), false)
{
  int num_interfaces = 0;
  for (size_t i = 0; i < mesh.pragmatic_interfaces.size(); ++i)
  {
    if (mesh.pragmatic_interfaces[i] == nullptr)
    {
      continue;
    }

    ++num_interfaces;
  }

  std::cout << "Create GroupedPartitionsSmooth Object" << std::endl;
  std::cout << "Smooth " << mesh.num_partitions << " partitions and " << num_interfaces << " interfaces" << std::endl;
  std::cout << "Running with " << mesh.nparts/2 << " threads" << std::endl;

  for (size_t i = 0; i < work_sets.size(); ++i)
  {
    work_sets[i].resize(3);
  }
}
//end of Constructor

//Destructor
GroupedPartitionsSmooth::~GroupedPartitionsSmooth()
{
  std::cout << "Delete GroupedPartitionsSmooth Object" << std::endl;
}
//end of Destructor

//Simple Laplace Smoother
bool GroupedPartitionsSmooth::SimpleLaplace(int no_iterations)
{
  std::cout << no_iterations << " iterations of Simple Laplace Smoother" << std::endl;
  //for loop number of iterations
  for (size_t iter = 1; iter <= no_iterations; ++iter)
  {
    volatile double dummy = 0.0;
    std::cout << "iteration " << iter << "/" << no_iterations << std::endl;
    std::cout << mesh.num_nodes << std::endl;
    //for loop over global indices
    for (index_t i = 0; i < mesh.num_nodes; ++i)
    {
      //skip vertex if it is on the global mesh boundary, or if it appears in more than 2 partitions/interfaces
      if ( mesh.boundary_nodes_mesh[i] )
      {  
        continue;
      } 
  
      if (mesh.vertex_appearances[i] > 2)
        continue;

      if ( mesh.vertex_appearances[i] == 2)
      {
        double vertex[2] {0.0, 0.0};
        mesh.GetCoords(i, vertex);
  
        std::vector<index_t> _NNList;
        //GetNNList_adv returns false, if vertex is in 2 different partitions
        //and true if it is in the assigned partition and interface
        if (!mesh.GetNNList_adv(i, _NNList))
        {
          continue;
        }

        int num_neighbors = _NNList.size();

        double q[2] {0.0, 0.0};

        for (size_t j = 0; j < num_neighbors; ++j)
        {
          double q_tmp[2] {0.0, 0.0};

          mesh.GetCoords(_NNList.at(j), q_tmp);

          q[0] += q_tmp[0];
          q[1] += q_tmp[1];
        }

        double h1 = (1.0 / static_cast<double>(_NNList.size()) );

        vertex[0] = h1 * q[0];
        vertex[1] = h1 * q[1]; 

        //dummy cpu expensive operation
        //dummy += atan(pow(3, i)) + sin(50*i) / sqrt (exp(i));
        //dummy = dummy_call(i, vertex[0]);
        
        //set new coordinates of vertex i
        mesh.SetCoords(i, vertex); 

      } //end of if for vertex_appearances == 2

      else if (mesh.vertex_appearances[i] == 1)
      {

        double vertex[2] {0.0, 0.0};

        //partition and interface are set by the GetNNList function to the partition or interface containint the vertex i
        //and to know, from which partition/interface the NNList is coming from, to apply the correct index_mapping
        int partition = -1;
        int interface = -1;
      
        //get vertex coordinates
        mesh.GetCoords(i, vertex);  

        //calculate new position
        //equation from MSc-Thesis "Glättung von Polygonnetzen in medizinischen Visualisierungen - J. Haase, 2005, eq. 2.2"
        std::vector<index_t>* _NNList = mesh.GetNNList(i, &partition, &interface);


        int num_neighbors = _NNList->size();

        double q[2] {0.0, 0.0};
        
        for (size_t j = 0; j < _NNList->size(); ++j)
        {
          double q_tmp[2];

          if (partition >= 0)
          { 

            mesh.GetCoords( mesh.l2g_vertex_map_partitions[partition].at(_NNList->at(j)), q_tmp);
            
            q[0] += q_tmp[0];
            q[1] += q_tmp[1];
          }

          else
          {

            mesh.GetCoords( mesh.l2g_vertex_map_interfaces[interface].at(_NNList->at(j)), q_tmp); 
            
            q[0] += q_tmp[0];
            q[1] += q_tmp[1];
          }
        } 

        double h1 = (1.0 / static_cast<double>(_NNList->size()) );
        
        vertex[0] = h1 * q[0];
        vertex[1] = h1 * q[1];  
        // end of calculate new position

        //dummy cpu expensive operation
        //dummy += atan(pow(3, i)) + sin(50*i) / sqrt (exp(i));
        //dummy = dummy_call(i, vertex[0]);

        //set the new coordinates of vertex i
        mesh.SetCoords(i, vertex);
      } //end of else

      else
      {
        return false;
      }
    } //end of for loop over global indices
  } //end of for loop number of iterations

  return true;
}
//end of Simple Laplace Smoother

//GroupedPartitionsSmooth::SimpleLaplaceOnGroups(int iterations)
//TODO: adapt GetNNList, since it is omitting partition vertex neighbors in an interface!!!
bool GroupedPartitionsSmooth::SimpleLaplaceOnGroups_sections(int iterations)
{
  std::cout << "sections" << std::endl;
  int num_threads = mesh.nparts/2;
  //int num_threads = 1;
  omp_set_dynamic(0);     // explicitly disable dynamic teams
  omp_set_num_threads(num_threads); // use a specified number of threads for all consecutive parallel regions 

/*
  //example work sets 
  work_sets[0].push_back(0);
  work_sets[0].push_back(1); 
  work_sets[0].push_back(0);
  
  work_sets[1].push_back(2);
  work_sets[1].push_back(3); 
  work_sets[1].push_back(5);
*/
  //parallel sections
#pragma omp parallel
{     
  //std::cout << omp_get_num_threads() << std::endl;
  
  #pragma omp single
  {/*
    std::cout << "Running with " << omp_get_num_threads() << " threads" << std::endl;
    std::cout << "iteration 1/2" << std::endl;
*/
    std::cout << "Assigning" << std::endl;
    clock_t tic = clock();
    AssignPartitionsAndInterface();
    clock_t toc = clock();
    std::cerr << "Time for assignement: " <<  static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;
/*    std::cerr << work_sets[0][0] <<  " " << work_sets[0][1] << " " << work_sets[0][2] << std::endl;
    std::cerr << work_sets[1][0] <<  " " << work_sets[1][1] << " " << work_sets[1][2] << std::endl;
   */ 
    for (size_t i = 0; i < work_sets.size(); ++i)
    {
      std::cout << "Thread " << i << ": " << work_sets[i][0] << " " << work_sets[i][1] << " " << work_sets[i][2] << std::endl;
    }
  }  

  #pragma omp sections nowait
  {
    #pragma omp section
    {
      //clock_t tic = clock();
      LaplaceKernel(work_sets[0][0], work_sets[0][1], work_sets[0][2]); 
      //clock_t toc = clock();
      std::cerr << "section 0" << std::endl;
      //LaplaceKernel(work_sets[0][0], work_sets[0][1], work_sets[0][2]); #threads=2
      //LaplaceKernel(0,1,0); //test_box_1000x1000 #threads=2
      //LaplaceKernel(0,1,0);   //test_box_1000x1000, #threads=4
      //LaplaceKernel(0,1,0); //test_box_2000x2000 #threads=2
      //LaplaceKernel(0,1,0); //test_box_2000x2000 #threads=4
    }
    
    #pragma omp section
    {
      //clock_t tic = clock();
      LaplaceKernel(work_sets[1][0], work_sets[1][1], work_sets[1][2]);
      //clock_t toc = clock();
      std::cerr << "section 1 " << std::endl;
      //LaplaceKernel(work_sets[1][0], work_sets[1][1], work_sets[1][2]); #threads=2
      //LaplaceKernel(2,3,5); //test_box_1000x1000, #threads=2
      //LaplaceKernel(2,3,13);  //test_box_1000x1000, #threads=4
      //LaplaceKernel(2,3,5); //test_box_2000x2000 #threads=2
      //LaplaceKernel(2,3,13); //test_box_2000x2000 #threads=4
    }

    //#threads=4

    #pragma omp section
    {
      LaplaceKernel(work_sets[2][0], work_sets[2][1], work_sets[2][2]);
      std::cerr << "section 2" << std::endl;
      //LaplaceKernel(4,5,22); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }

    #pragma omp section
    {
      LaplaceKernel(work_sets[3][0], work_sets[3][1], work_sets[3][2]);
      std::cerr << "section 3" << std::endl;
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }
/*
    //#threads=6

    #pragma omp section
    {
      std::cerr << "section 4" << std::endl;
      LaplaceKernel(work_sets[4][0], work_sets[4][1], work_sets[4][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }

    #pragma omp section
    {
      std::cerr << "section 5" << std::endl;
      LaplaceKernel(work_sets[5][0], work_sets[5][1], work_sets[5][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }
/*
    //#threads=8

    #pragma omp section
    {
      LaplaceKernel(work_sets[6][0], work_sets[6][1], work_sets[6][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }

    #pragma omp section
    {
      LaplaceKernel(work_sets[7][0], work_sets[7][1], work_sets[7][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }
//*/
  } //end of sections

  //#pragma omp barrier   //should not be here, because at the end of pragma omp sections is an implicit barrier
  //now switch assignement of partitions and interfaces
  //parallel sections

  #pragma omp single
  {
    //std::cout << "iteration 2/2" << std::endl;
    std::cout << "Assigning" << std::endl;
    std::fill(partition_assigned.begin(), partition_assigned.end(), false);
    std::fill(interface_assigned.begin(), interface_assigned.end(), false);

    clock_t tic = clock();
    AssignPartitionsAndInterface();
    clock_t toc = clock();
    
    std::cerr << "Time for assignement: " <<  static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;    
/*    std::cout << work_sets[0][0] <<  " " << work_sets[0][1] << " " << work_sets[0][2] << std::endl;
    std::cout << work_sets[1][0] <<  " " << work_sets[1][1] << " " << work_sets[1][2] << std::endl;
  */  

    for (size_t i = 0; i < work_sets.size(); ++i)
    {
      std::cout << "Thread " << i << ": " << work_sets[i][0] << " " << work_sets[i][1] << " " << work_sets[i][2] << std::endl;
    }

  }  

  #pragma omp sections nowait
  {    
    #pragma omp section
    {
      LaplaceKernel(work_sets[0][0], work_sets[0][1], work_sets[0][2]);
      //LaplaceKernel(0, 3, 2); //test_box_300x300, 1000x1000 #threads=2
      //LaplaceKernel(0,6,5);   //test_box_1000x1000, #threads=4
      //LaplaceKernel(0,2,1); //test_box_2000x2000 #threads=2
      //LaplaceKernel(0,2,1); //test_box_2000x2000 #threads=4
    }
    
    #pragma omp section
    {
      LaplaceKernel(work_sets[1][0], work_sets[1][1], work_sets[1][2]);
      //LaplaceKernel(1, 2, 3); //test_box_300x300, 1000x1000 #threads=2
      //LaplaceKernel(2,5,15);   //test_box_1000x1000, #threads=4
      //LaplaceKernel(1,3,4); //test_box_2000x2000 #threads=2
      //LaplaceKernel(4,3,18); //test_box_2000x2000 #threads=4
    }

    //#threads=4

    #pragma omp section
    {
      LaplaceKernel(work_sets[2][0], work_sets[2][1], work_sets[2][2]);
      //LaplaceKernel(4,7,24); //test_box_1000x1000 #threads=4
      //LaplaceKernel(5,7,26); //test_box_2000x2000 #threads=4
    }

    #pragma omp section
    {
      LaplaceKernel(work_sets[3][0], work_sets[3][1], work_sets[3][2]);
      //LaplaceKernel(1,3,8); //test_box_1000x1000 #threads=4
      //LaplaceKernel(6,1,11); //test_box_2000x2000 #threads=4
    }
/*
    //#threads=6

    #pragma omp section
    {
      LaplaceKernel(work_sets[4][0], work_sets[4][1], work_sets[4][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }

    #pragma omp section
    {
      LaplaceKernel(work_sets[5][0], work_sets[5][1], work_sets[5][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }
/*
    //#threads=8

    #pragma omp section
    {
      LaplaceKernel(work_sets[6][0], work_sets[6][1], work_sets[6][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }

    #pragma omp section
    {
      LaplaceKernel(work_sets[7][0], work_sets[7][1], work_sets[7][2]);
      //LaplaceKernel(6,7,27); //test_box_2000x2000 #threads=4, test_box_1000x1000, #threads=4
    }
//*/
  } //end of sections*/
} //end of pragma omp parallel

  //now process the left out interfaces
  for (size_t i = 0; i < mesh.pragmatic_interfaces.size(); ++i)
  {
    if (mesh.pragmatic_interfaces[i] == nullptr)
    {
      continue;
    }

    if (num_touches_interfaces[i] != 2)
    {
      std::cout << "interface " << i << " has to be processed for another " << iterations - num_touches_interfaces[i] << " time(s)" << std::endl;
    }
  }
}
//end of GroupedPartitionsSmooth::SimpleLaplaceOnGroups(int iterations)

//GroupedPartitionsSmooth::SimpleLaplaceOnGroups_tasks(int iterations)
bool GroupedPartitionsSmooth::SimpleLaplaceOnGroups_tasks(int iterations)
{
  int num_threads = mesh.nparts/2;
  omp_set_dynamic(0);               // explicitly disable dynamic teams
  omp_set_num_threads(num_threads); // use a specified number of threads for all consecutive parallel regions 

  #pragma omp parallel
  {
    #pragma omp master
    {
      std::fill(partition_assigned.begin(), partition_assigned.end(), false);
      std::fill(interface_assigned.begin(), interface_assigned.end(), false);
      
      std::cerr << "Assigning" << std::endl;
      AssignPartitionsAndInterface();
    } //end of master region
    #pragma omp barrier

    //this is necessary, otherwise the program segfaults!
    //#pragma omp taskwait

    //now spawn the tasks
    #pragma omp single nowait
    for (size_t i = 0; i < num_threads; ++i)
    {
      std::cerr << "task " << i << std::endl;
      #pragma omp task
      LaplaceKernel(work_sets[i][0], work_sets[i][1], work_sets[i][2]);
    }
    
    #pragma omp taskwait

    //now do a second iteration
    #pragma omp master 
    {
      std::cerr << "Assigning" << std::endl;
      AssignPartitionsAndInterface();
    }

    #pragma omp barrier

    //wait for all tasks to finish, then proceed
    //#pragma omp taskwait

    //now spawn the tasks
    #pragma omp single nowait
    for (size_t i = 0; i < num_threads; ++i)
    {
      std::cerr << "task " << i << std::endl;
      #pragma omp task
      LaplaceKernel(work_sets[i][0], work_sets[i][1], work_sets[i][2]);
    }   

    #pragma omp taskwait 
  } //end of pragma omp parallel
}
//end of GroupedPartitionsSmooth::SimpleLaplaceOnGroups_tasks(int iterations)

//GroupedPartitionsSmooth::LaplaceKernel(int part1, int part2, int inter)
void GroupedPartitionsSmooth::LaplaceKernel(int part1, int part2, int inter)
{
/*
  int loops = 0;
  int bdry_ctr = 0;
  int interior_ctr = 0;
  int bdry_inter_part_ctr = 0;
  int more_parts = 0;
  int failed_inter_check = 0;
*/
  //for loop over part1
  //TODO: idea: iterate over the corresponding global_to_local_index_mappings and update the vertices if necessary in more than 1 pragmatic mesh!    
  for (auto it : mesh.g2l_vertex_map_partitions[part1])  
  {
    volatile double dummy = 0.0;
    //if vertex is on the boundary of the global mesh or appears on more than 2 sub-meshes do not touch it!
    if (mesh.boundary_nodes_mesh[it.first] || mesh.vertex_appearances[it.first] > 2)
    {
      //++bdry_ctr;
      //++loops;
      continue;
    }

    //if vertex is on the interior boundary of the partition check if it is in the assigned interface
    //if yes, update it in both submeshes, if no, do not touch it!
    //TODO: if it is in partition & boundary, the NNLists have to be merged to compute correct results of the Laplace smoothing!!!
    if ( mesh.boundary_nodes_partitions[part1].at(it.second) )
    {
      //check interface
      if ( mesh.g2l_vertex_map_interfaces[inter].find(it.first) != mesh.g2l_vertex_map_interfaces[inter].end() )
      {
        //std::cout << it.second << " is in interface with id " << mesh.g2l_vertex_map_interfaces[inter].at(it.first) << std::endl;

        //get old vertex coordinates
        double vertex[2] {0.0, 0.0};
      
        //mesh.GetCoords(it.first, vertex);
        mesh.GetCoords_inter(mesh.g2l_vertex_map_interfaces[inter].at(it.first), vertex, inter);
        //mesh.GetCoords_part(it.second, vertex, part1);  
  
        //calculate new position
        //equation from MSc-Thesis "Glättung von Polygonnetzen in medizinischen Visualisierungen - J. Haase, 2005, eq. 2.2"

        int partition, interface; //TODO: this is unnecessary with this algorithm!!!
        std::vector<index_t> _NNList;
        if ( !mesh.GetNNList_adv(it.first, _NNList))
        {
          //std::cout << it.first << " is in 2 or more interfaces! " << mesh.g2l_vertex_map_interfaces[inter].at(it.first) << " " << inter << std::endl;
          //++more_parts;
          //++loops;
          continue;
        }
        //std::vector<index_t>* _NNList = mesh.GetNNList_(it.first, &partition, &interface);
        int num_neighbors = _NNList.size();

        double q[2] {0.0, 0.0};

        for (size_t j = 0; j < num_neighbors; ++j)
        {
          double q_tmp[2] {0.0, 0.0};

          //GetNNList_adv already returns global indices, therefore no conversion is needed below!!!
          mesh.GetCoords(_NNList.at(j), q_tmp);

          q[0] += q_tmp[0];
          q[1] += q_tmp[1];
        }

        double h1 = (1.0 / static_cast<double>(num_neighbors) );

        vertex[0] = h1 * q[0];
        vertex[1] = h1 * q[1];  
        // end of calculate new position

        //set the new coordinates of vertex it.first
        //mesh.SetCoords(it.first, vertex);

        //dummy cpu expensive operation
        //dummy += atan(pow(3, it.first)) + sin(50*it.second) / sqrt (exp(it.first));
        //dummy = dummy_call (it.first, vertex[0]);

        mesh.SetCoords_part(it.second, vertex, part1);
        mesh.SetCoords_inter(mesh.g2l_vertex_map_interfaces[inter].at(it.first), vertex, inter);
        
        //++bdry_inter_part_ctr;
      } //end of check interface
/*
      else
        ++failed_inter_check;
*/
    } //end of if mesh.boundary_nodes_partition

    //if vertex is not on any boundary of the assigned partition update its coordinates only in the partition
    else
    {
      //get old vertex coordinates
      double vertex[2] {0.0, 0.0};
      
      mesh.GetCoords_part(it.second, vertex, part1);
      //mesh.GetCoords(it.first, vertex);  
  
      //calculate new position
      //equation from MSc-Thesis "Glättung von Polygonnetzen in medizinischen Visualisierungen - J. Haase, 2005, eq. 2.2"

      int partition, interface; //TODO: this is unnecessary with this algorithm!!!
      // std::vector<index_t> _NNList; 
      // mesh.GetNNList_adv(it.first, _NNList);
      std::vector<index_t>* _NNList = mesh.GetNNList(it.first, &partition, &interface);
      int num_neighbors = _NNList->size();

      double q[2] {0.0, 0.0};

      for (size_t j = 0; j < num_neighbors; ++j)
      {
        double q_tmp[2] {0.0, 0.0};

        mesh.GetCoords(mesh.l2g_vertex_map_partitions[part1].at(_NNList->at(j)), q_tmp);

        q[0] += q_tmp[0];
        q[1] += q_tmp[1];
      }

      double h1 = (1.0 / static_cast<double>(num_neighbors) );

      vertex[0] = h1 * q[0];
      vertex[1] = h1 * q[1];  
      // end of calculate new position

      //dummy cpu expensive operation
      //dummy += atan(pow(3, it.first)) + sin(50*it.second) / sqrt (exp(it.first));
      //dummy = dummy_call (it.first, vertex[0]);

      //set the new coordinates of vertex it.first
      //mesh.SetCoords(it.first, vertex);
      mesh.SetCoords_part(it.second, vertex, part1);

      //++interior_ctr;
    }//end of else

  //++loops;
  } //end of for loop over part1
/*
  std::cerr << "part1 " << part1 << ": " <<  mesh.pragmatic_partitions[part1]->get_number_nodes() << ", loops: " << loops << std::endl;
  std::cerr << part1 << ", bdry_ctr: " << bdry_ctr << ", interior_ctr: " << interior_ctr << ", bdry_inter_part_ctr: "  << bdry_inter_part_ctr << ", more_parts: " << more_parts << ", failed_inter_check: " << failed_inter_check << ", sum: " << bdry_ctr+interior_ctr+bdry_inter_part_ctr+more_parts+failed_inter_check<< std::endl;

  loops = 0;
  bdry_ctr = 0;
  interior_ctr = 0;
  bdry_inter_part_ctr = 0;
  more_parts = 0;
  failed_inter_check = 0;
*/
  //for loop over part2
  for (auto it : mesh.g2l_vertex_map_partitions[part2])  
  {
    volatile double dummy = 0.0;
    //if vertex is on the boundary of the global mesh or appears on more than 2 sub-meshes do not touch it!
    if (mesh.boundary_nodes_mesh[it.first] || mesh.vertex_appearances[it.first] > 2)
    {
      //++bdry_ctr;
      //++loops;
      continue;
    }

    //if mesh is on the interior boundary of the partition check if it is in the assigned interface
    //if yes, update it in both submeshes, if no, do not touch it!
    //TODO: if it is in partition & boundary, the NNLists have to be merged to compute correct results of the Laplace smoothing!!!
    if ( mesh.boundary_nodes_partitions[part2].at(it.second) )
    {
      //check interface
      if ( mesh.g2l_vertex_map_interfaces[inter].find(it.first) != mesh.g2l_vertex_map_interfaces[inter].end() )
      {
        //std::cout << it.second << " is in interface with id " << mesh.g2l_vertex_map_interfaces[inter].at(it.first) << std::endl;

        //get old vertex coordinates
        double vertex[2] {0.0, 0.0};
      
        mesh.GetCoords_inter(mesh.g2l_vertex_map_interfaces[inter].at(it.first), vertex, inter);
        //mesh.GetCoords_part(it.second, vertex, inter);
        //mesh.GetCoords(it.first, vertex);  
  
        //calculate new position
        //equation from MSc-Thesis "Glättung von Polygonnetzen in medizinischen Visualisierungen - J. Haase, 2005, eq. 2.2"

        int partition, interface; //TODO: this is unnecessary with this algorithm!!!
        std::vector<index_t> _NNList;
        if ( !mesh.GetNNList_adv(it.first, _NNList))
        {
          //std::cout << it.first << " is in 2 or more interfaces! " << mesh.g2l_vertex_map_interfaces[inter].at(it.first) << " " << inter << std::endl;
          //++more_parts; 
          //++loops;         
          continue;
        }
        //std::vector<index_t>* _NNList = mesh.GetNNList(it.first, &partition, &interface);
        int num_neighbors = _NNList.size();

        double q[2] {0.0, 0.0};

        for (size_t j = 0; j < num_neighbors; ++j)
        {
          double q_tmp[2] {0.0, 0.0};
        
          //GetNNList_adv already returns global indices, therefore no conversion is needed below!!!
          mesh.GetCoords(_NNList.at(j), q_tmp);

          q[0] += q_tmp[0];
          q[1] += q_tmp[1];
        }

        double h1 = (1.0 / static_cast<double>(num_neighbors) );

        vertex[0] = h1 * q[0];
        vertex[1] = h1 * q[1];  
        // end of calculate new position

        //set the new coordinates of vertex it.first
        //mesh.SetCoords(it.first, vertex);

        //dummy cpu expensive operation
        //dummy += atan(pow(3, it.first)) + sin(50*it.second) / sqrt (exp(it.first));
        //dummy = dummy_call (it.first, vertex[0]);

        mesh.SetCoords_part(it.second, vertex, part2);
        mesh.SetCoords_inter(mesh.g2l_vertex_map_interfaces[inter].at(it.first), vertex, inter);  

        //++bdry_inter_part_ctr;
      } //end of check interface
/*
      else
        ++failed_inter_check;
*/
    } //end of if mesh.boundary_nodes_partition

    //if mesh is not on any boundary of the assigned partition update its coordinates only in the partition
    else
    {
      //get old vertex coordinates
      double vertex[2] {0.0, 0.0};
      
      mesh.GetCoords_part(it.second, vertex, part2);
      //mesh.GetCoords(it.first, vertex);  
  
      //calculate new position
      //equation from MSc-Thesis "Glättung von Polygonnetzen in medizinischen Visualisierungen - J. Haase, 2005, eq. 2.2"

      int partition, interface; //TODO: this is unnecessary with this algorithm!!!
      std::vector<index_t>* _NNList = mesh.GetNNList(it.first, &partition, &interface);
      int num_neighbors = _NNList->size();

      double q[2] {0.0, 0.0};

      for (size_t j = 0; j < num_neighbors; ++j)
      {
        double q_tmp[2] {0.0, 0.0};

        mesh.GetCoords(mesh.l2g_vertex_map_partitions[part2].at(_NNList->at(j)), q_tmp);

        q[0] += q_tmp[0];
        q[1] += q_tmp[1];
      }

      double h1 = (1.0 / static_cast<double>(num_neighbors) );

      vertex[0] = h1 * q[0];
      vertex[1] = h1 * q[1];  
      // end of calculate new position

      //dummy cpu expensive operation
      //dummy += atan(pow(3, it.first)) + sin(50*it.second) / sqrt (exp(it.first));
      //dummy = dummy_call (it.first, vertex[0]);

      //set the new coordinates of vertex it.first
      //mesh.SetCoords(it.first, vertex);
      mesh.SetCoords_part(it.second, vertex, part2);

      //++interior_ctr;
    }//end of else

  //++loops;
  } //end of for loop over part2
/*
  std::cerr << "part2 " << part2 << ": " <<  mesh.pragmatic_partitions[part2]->get_number_nodes()  << ", loops: " << loops << std::endl;
  std::cerr << part2 << ", bdry_ctr: " << bdry_ctr << ", interior_ctr: " << interior_ctr << ", bdry_inter_part_ctr: "  << bdry_inter_part_ctr << ", more_parts: " << more_parts << ", failed_inter_check: " << failed_inter_check << ", sum: " << bdry_ctr+interior_ctr+bdry_inter_part_ctr+more_parts+failed_inter_check<< std::endl;

  loops = 0;
  bdry_ctr = 0;
  interior_ctr = 0;
  bdry_inter_part_ctr = 0;
  more_parts = 0;
*/
  //for loop over interface 
  //adapt only interior vertices, since the boundary vertices have already been updated in the loops over both partitions!!!  
  for (auto it : mesh.g2l_vertex_map_interfaces[inter])  
  {    
    volatile double dummy = 0.0;
    //check if vertex is on any interface boundary (includes already the global mesh boundary)
    //TODO: boundary_nodes_interfaces contains too many elements ==> adapt the algorithm in grouped_partitions.hpp!!!
    if ( mesh.boundary_nodes_interfaces[inter].at(it.second) )
    {
      //++bdry_ctr;
      //++loops;
      continue;
    }
    //get old vertex coordinates
    double vertex[2] {0.0, 0.0};
    
    mesh.GetCoords_inter(it.second, vertex, inter);  
    //mesh.GetCoords(it.first, vertex);  
  
    //calculate new position
    //equation from MSc-Thesis "Glättung von Polygonnetzen in medizinischen Visualisierungen - J. Haase, 2005, eq. 2.2"

    int partition, interface; //TODO: this is unnecessary with this algorithm!!!
    std::vector<index_t>* _NNList = mesh.GetNNList(it.first, &partition, &interface);
    int num_neighbors = _NNList->size();
    double q[2] {0.0, 0.0};

    for (size_t j = 0; j < num_neighbors; ++j)
    {
      double q_tmp[2] {0.0, 0.0};

      mesh.GetCoords(mesh.l2g_vertex_map_interfaces[inter].at(_NNList->at(j)), q_tmp);

      q[0] += q_tmp[0];
      q[1] += q_tmp[1];
    }
    double h1 = (1.0 / static_cast<double>(num_neighbors) );

    vertex[0] = h1 * q[0];
    vertex[1] = h1 * q[1];  
    // end of calculate new position

    //dummy cpu expensive operation
    //dummy += atan(pow(3, it.first)) + sin(50*it.second) / sqrt (exp(it.first));
    //dummy = dummy_call (it.first, vertex[0]);

    //set the new coordinates of vertex i
    //mesh.SetCoords(it.first, vertex);
    mesh.SetCoords_inter(it.second, vertex, inter);
    //++interior_ctr;
    //++loops;
  } //end of for loop over interface
/*  
  std::cerr << "inter " << inter << ": " <<  mesh.pragmatic_interfaces[inter]->get_number_nodes() << ", loops: " << loops << std::endl;
  std::cerr << inter << ", bdry_ctr: " << bdry_ctr << ", interior_ctr: " << interior_ctr << ", bdry_inter_part_ctr: " << bdry_inter_part_ctr << ", sum: " << bdry_ctr+interior_ctr+bdry_inter_part_ctr << std::endl;
*/
}
//end of GroupedPartitionsSmooth::LaplaceKernel(int part1, int part2, int inter)

//GroupedPartitionsSmooth::AssignPartitionsAndInterface()
void GroupedPartitionsSmooth::AssignPartitionsAndInterface()
{
  //std::cout << "Assign Partitions And Interfaces! " << std::endl;
/*
  for (size_t i = 0; i < mesh.partition_neighbors.size(); ++i)
  {
    std::cout << "Neighbors of Partition " << i << std::endl;
    for (auto it : mesh.partition_neighbors[i])
    {
      std::cout << "  " << it << std::endl;
    }
  }
*/
/*
  for (size_t i = 0; i < mesh.interface_neighbors.size(); ++i)
  {
    if (mesh.interface_neighbors[i].size() == 0)
      continue; 

    //std::cout << "Neighbors of Interface " << i << std::endl;
    for (auto it : mesh.interface_neighbors[i])
    {
      //std::cout << "  " << it << std::endl;
    }
  }
*/
  //create work sets
  //IDEA: iterate over num_touches_interfaces to make sure, that every interface and every partition is touched exactly 
  // the number of total iterations (parameter int iterations)

  //assign a unique combination of paritions and interface to each thread  
  //for (size_t thread = 0; thread < work_sets.size(); ++thread)
  //{ 
    int thread = 0;
    int iteration = 1;

    for (size_t i = 0; i < num_touches_interfaces.size(); ++i)
    {
      std::cout << i << ": " << num_touches_interfaces[i] << std::endl;
      //skip if empty interface
      if (mesh.interface_neighbors[i].size() == 0 || interface_assigned[i])
      {
        std::cout << "continuing after " << i << std::endl;
        continue;
      }

      else if (num_touches_interfaces[i] == iteration)
      {
        std::cout << "skip " << i << std::endl;
        continue;
      } 

      //else assign interface and partitions
      else
      {
        //std::set<int>::iterator set_iterator = mesh.interface_neighbors[i].begin();
        //std::cout << interface_assigned[i] << " " << partition_assigned[*(set_iterator++)] << " " << partition_assigned[*(set_iterator)] << std::endl;

        //set_iterator = mesh.interface_neighbors[i].begin();

        //check if interface or its corresponding partitions are already assigned
        //if ( interface_assigned[i] || partition_assigned[*(set_iterator++)] || partition_assigned[*(set_iterator)] ) //this is for interface_assigned and partition_assigned being sets
        if( interface_assigned[i] || partition_assigned[ mesh.interface_neighbors[i][0] ] || partition_assigned[ mesh.interface_neighbors[i][1] ]  )
        {
          //set_iterator = mesh.interface_neighbors[i].begin();
          std::cout << "interface " << i << " can not be assigned!" << std::endl;
          std::cout << i << ": " << interface_assigned[i] << ", " << mesh.interface_neighbors[i][0] << ": " << partition_assigned[ mesh.interface_neighbors[i][0] ] << ", " << mesh.interface_neighbors[i][1] << ": " << partition_assigned[ mesh.interface_neighbors[i][1] ] << std::endl;
          continue;
        }

        //std::cout << "assign " << i << std::endl;
/*
        set_iterator = mesh.interface_neighbors[i].begin();
        std::cout << interface_assigned[i] << " " << partition_assigned[*(set_iterator++)] << " " << partition_assigned[*(set_iterator)] << std::endl;
        
        set_iterator = mesh.interface_neighbors[i].begin();

        work_sets[thread][0] = *(set_iterator);  //partition 1
        ++num_touches_partitions[*(set_iterator)];
        std::cout << "part1: " << *(set_iterator) << std::endl;
        partition_assigned[*(set_iterator++)] = true;

        work_sets[thread][1] = *(set_iterator);  //partition 2
        ++num_touches_partitions[*(set_iterator)];
        std::cout << "part2: " << *(set_iterator) << std::endl;
        partition_assigned[*(set_iterator++)] = true;

        work_sets[thread][2] = i;  //interface
        ++num_touches_interfaces[i];
        interface_assigned[i] = true;
*/

        std::cout << interface_assigned[i] << " " << partition_assigned[ mesh.interface_neighbors[i][0] ] << " " << partition_assigned[ mesh.interface_neighbors[i][1] ] << std::endl;
        
        work_sets[thread][0] = mesh.interface_neighbors[i][0];  //partition 1
        ++num_touches_partitions[ mesh.interface_neighbors[i][0] ];
        std::cout << "part1: " << mesh.interface_neighbors[i][0] << std::endl;
        partition_assigned[ mesh.interface_neighbors[i][0] ] = true;

        work_sets[thread][1] = mesh.interface_neighbors[i][1];  //partition 2
        ++num_touches_partitions[ mesh.interface_neighbors[i][1] ];
        std::cout << "part2: " << mesh.interface_neighbors[i][1] << std::endl;
        partition_assigned[ mesh.interface_neighbors[i][1] ] = true;

        work_sets[thread][2] = i;  //interface
        ++num_touches_interfaces[i];
        interface_assigned[i] = true;

        ++thread;
      }

      if (thread == mesh.nparts/2)
        break;
    } //end of create work sets

    if (thread < mesh.nparts/2)
    {
      //std::cout << "Did not assign work to " << mesh.nparts/2 - thread << " threads" << std::endl;
      //assign the 2 partitions left, together with any available interface
      std::vector<bool>::iterator part;
      
      part = std::find( partition_assigned.begin(), partition_assigned.end(), false );  //find partition
      int index = std::distance(partition_assigned.begin(), part);
      partition_assigned[ index ] = true;
      work_sets[thread][0] = index;
      ++num_touches_partitions[ index ];
      
      part = std::find( partition_assigned.begin(), partition_assigned.end(), false );  //find partition
      index = std::distance(partition_assigned.begin(), part);      
      partition_assigned[ index ]=true;
      work_sets[thread][1] = index;
      ++num_touches_partitions[ index ];

      part = std::find( interface_assigned.begin(), interface_assigned.end(), false );  //find interface
      index = std::distance(interface_assigned.begin(), part); 
      if (num_touches_interfaces[index] == 1)
      {
        int counter=1;    //counter set to 1, otherwise the first try in the do-while loop would be the same interface as found in the first place
        bool interface_found = false;
        do
        {
          part = std::find(part+counter, interface_assigned.end(), false );
          index = std::distance(interface_assigned.begin(), part);
          ++counter;
          if (mesh.pragmatic_interfaces[index] != nullptr && num_touches_interfaces[index] != 1)
          {
            interface_found = true;
          }
        } while (!interface_found);
          
        //std::cout << index << " " << interface_assigned.size() << std::endl;
      }
      interface_assigned[ index ] = true;
      work_sets[thread][2] = index;
      ++num_touches_interfaces[ index ];
      
      ++thread;
      //std::cout << "Assigned leftover partitions together with any available interface" << std::endl;
    }
  //}
} //end of GroupedPartitionsSmooth::AssignPartitionsAndInterface()

//GroupedPartitionsSmooth::Evaluate()
void GroupedPartitionsSmooth::Evaluate()
{
  std::cout << std::endl << "Number of Touches" << std::endl;
  for (size_t i = 0; i < num_touches_partitions.size(); ++i)
  {
    std::cout << "Partition " << i << ": " << num_touches_partitions[i] <<std::endl;
  }
  std::cout << std::endl;
  for (size_t i = 0; i < num_touches_interfaces.size(); ++i)
  {
    if (mesh.pragmatic_interfaces[i] == nullptr)
      continue;
    std::cout << "Interface " << i << ": " << num_touches_interfaces[i] <<std::endl;
  }
}
//end of GroupedPartitionsSmooth::Evaluate()

//----------------------------------------------------------------------------------------------------------------------------------------------//
//                                                                     End                                                                      //
//----------------------------------------------------------------------------------------------------------------------------------------------//

#endif
