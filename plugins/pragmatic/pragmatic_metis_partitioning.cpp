#include "pragmatic_metis_partitioning.hpp"

//standard includes
#include <vector>
#include <string>
#include <cassert>

//Pragmatic includes
#include "pragmatic_mesh.hpp"
#include "Smooth.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "metis.h"

//Defines
#define NUM_THREADS 2

//
// Define the necessary types:
//
typedef viennagrid::result_of::region<MeshType>::type                   RegionType;
typedef viennagrid::result_of::element<MeshType>::type                  VertexType;
typedef viennagrid::result_of::point<MeshType>::type                    PointType;
typedef viennagrid::result_of::element<MeshType>::type      	          CellType;

typedef viennagrid::result_of::cell_range<RegionType>::type             CellRange;
typedef viennagrid::result_of::iterator<CellRange>::type                CellIterator;

typedef viennagrid::result_of::const_cell_range<MeshType>::type         ConstCellRangeType;
typedef viennagrid::result_of::iterator<ConstCellRangeType>::type       ConstCellRangeIterator;

typedef viennagrid::result_of::element<MeshType>::type                  ElementType;
typedef viennagrid::result_of::const_element_range<ElementType>::type   ConstElementRangeType;
typedef viennagrid::result_of::iterator<ConstElementRangeType>::type    ConstElementRangeIterator;


typedef viennagrid::result_of::element<MeshType>::type      	          CellType;
typedef viennagrid::mesh                                                MeshType;

namespace viennamesh
{
		pragmatic_metis_partitioner::pragmatic_metis_partitioner()	{}
		std::string pragmatic_metis_partitioner::name() {return "pragmatic_metis_partitioner";}

    //
    //void laplace_smoother()
    //
    //TODO: implement a laplacian smoothing algorithm
    void laplace_smoother(Mesh<double>*& mesh, size_t max_iterations = 4)
    {
      std::cout << "Laplace smoother" << std::endl;
/*
      std::cout << &mesh << std::endl;
*/
      int NNodes = mesh->get_number_nodes();
      int NElements = mesh->get_number_elements();
      int nloc = mesh->get_nloc();

      std::vector<std::atomic<bool>> is_boundary(NNodes);
      std::vector<std::atomic<bool>> active_vertices(NNodes);

      //find all boundary nodes
      //TODO: find partition boundary nodes and parallalize this section
      
      //initalize is_boundary vector      
      for (size_t i = 0; i < NNodes; ++i)
      {
        is_boundary[i].store(false, std::memory_order_relaxed);
      }

      for(int n=0; n<NNodes; ++n) 
      {
        is_boundary[n].store(false, std::memory_order_relaxed);

        std::vector<index_t>* _NNList = mesh->get_nnlist(n);   
        std::vector<int>::iterator boundary_iterator = mesh->get_boundary();
        std::vector<int>* boundary = mesh->get_boundary_vector();
       // std::vector<std::vector<index_t>>::iterator NNList_iterator = mesh->get_nnlist_iterator();  
/*
        NNList_iterator += n;
       
        std::cout << n << ": " << std::endl;
        
        for (std::vector<index_t>::iterator it = NNList_iterator->begin(); it != NNList_iterator->end(); ++it)
        {
          std::cout << "  " << *it << std::endl;
        }

        //std::cout << n << ": " << (_NNList->size()) << std::endl << std::endl;
        std::cout << n << ": " << std::endl;

        for (size_t i = 0; i < _NNList->size(); ++i)
        {
          std::cout << "  " << _NNList->at(i) << std::endl;
        }
        std::cout << std::endl;
*/ 
/* 
        for (std::vector<index_t>::iterator)
        {
          
        }
*/      
        if(_NNList->empty())
        {
          active_vertices[n].store(false, std::memory_order_relaxed);
        } 

        else 
        {
          active_vertices[n].store(true, std::memory_order_relaxed);
        }
      
        for (size_t i = 0; i < NElements; ++i)
        {
          const int *n = mesh->get_element(i);

          if (n[0] < 0) //if element is marked for deletion (its first node in ENList is set to -1), we skip it
          {
            continue;
          }

          for (size_t j = 0; j < nloc; ++j)
          {
            boundary_iterator = mesh->get_boundary();
            std::advance(boundary_iterator, i*nloc+j);

            if (*boundary_iterator > 0)
            {
              for (size_t k = 1; k < nloc; ++k)
              {
                is_boundary[n[(j+k)%nloc]].store(true, std::memory_order_relaxed);
              }
            }
          }
        }
       }
      //end of find all boundary nodes section

      //output boundary vector for developing and debugging purposes
/*
      std::cout << "IS BOUNDARY" << std::endl;
      
      for (size_t i = 0; i < is_boundary.size(); ++i)
      {
        std::cout << i << ": " << is_boundary[i] << std::endl; 
      }
      std::cout << std::endl;
*/
      //TODO: now implement the laplacian smoothing algorithm here!!!!
      // in a first try, simply smooth the test mesh with all vertices and in a second try use the partition grouping 

      std::vector<double> coords(NNodes*2);
      std::vector<std::atomic<bool>> invalid_area(NNodes);
      //set new position according to the formula from "Glättung von Polygonnetzen in medizinischen Visualisierungen-Jens Haase", p.31
      for(size_t i = 0; i < NNodes; ++i)
      {
        if ( is_boundary[i].load(std::memory_order_relaxed) )
        { 
          continue;
        }
        const int *n = mesh->get_element(i);
        
        ElementProperty<double> *property;
        property = new ElementProperty<double>(mesh->get_coords(n[0]),
                                               mesh->get_coords(n[1]),
                                               mesh->get_coords(n[2]));

        double p[2] = {0.0, 0.0};

        std::vector<index_t>* _NNList = mesh->get_nnlist(i); //get Node-Node-adjacency list
       // const double *coords_ptr = mesh->get_coords(i);
        
        double q[2] = {0.0, 0.0};        

     //  std::cout << i << std::endl;

        for (size_t j = 0; j < _NNList->size()-1; ++j)
        {
          double q_tmp[2] = {0.0, 0.0};
          mesh->get_coords(_NNList->at(j), q_tmp);

          q[0] += q_tmp[0];
          q[1] += q_tmp[1];
        }
          
        double num_neighbors = _NNList->size();

        p[0] = (1/num_neighbors) * q[0];
        p[1] = (1/num_neighbors) * q[1];

        //check for validity
        
        double area;
        for (size_t j = 0; j < _NNList->size(); ++j)
        {
          const int *node = mesh->get_element(j);
          const double *x0 = mesh->get_coords(node[0]);
          const double *x1 = mesh->get_coords(node[1]);
          const double *x2 = mesh->get_coords(node[2]);
        
          if ( node[0] == i )
          {
            area = property->area(p, x1, x2);
          }

          else if ( node[1] == i )
          {
            area = property->area(x0, p, x2);
          }

          else if ( node[2] == i )
          {
            area = property->area(x0, x1, p);
          }

          if ( area < 0)
            {
              invalid_area[i].store(true, std::memory_order_relaxed);
              continue;
            }
          
          invalid_area[i].store(false, std::memory_order_relaxed);
        }

        coords[i*2]   = p[0];
        coords[i*2+1] = p[1];
     
      //  mesh->set_coords(i, p);
     //   std::cout << endl;
        //std::cout << i << ": " << _NNList->size() << std::endl;
      }

      for (size_t n = 0; n < NNodes; ++n)
      {
        double p[2];

        if ( is_boundary[n].load(std::memory_order_relaxed) || invalid_area[n].load(std::memory_order_relaxed))
        {
          mesh->get_coords(n, p);
          std::cout << n << ": " << p[0] << ", " << p[1] << std::endl;
          continue;
        }

        p[0] = coords[n*2];
        p[1] = coords[n*2+1];
    
        std::cout << n << ": " << p[0] << ", " << p[1] << std::endl;
        mesh->set_coords(n, p);
      }

      size_t iter = 0;
      while (iter < max_iterations)
      {
        ++iter;
      } //end of while loop
      //end of laplacian smoothing algorithm
    
      return;
    }

    //
    //binomial_coefficient(int n, int k)
    //
    //Returns the binomial coefficient n! / ( k! (n-k)! )
    //
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

    //
    //create_pragmatic_meshes()
    //
    //create pragmatic meshes for debugging purposes
    void create_pragmatic_meshes()
    {
    } //end of create_pragmatic_meshes()
  
    //
    //pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
    //
    //Uses metis to partition the mesh (represented in pragmatic data structure)
    bool pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
		{
      std::cout << name() << std::endl;

      //create mesh_handle to read input mesh			
		  mesh_handle input_mesh = get_required_input<mesh_handle>("mesh");
            
      //create data_handle for optional inputs
		  data_handle<int> region_count = get_input<int>("region_count");

      //convert viennamesh into pragmatic mesh data structure
      Mesh<double> *mesh = nullptr;

		  mesh = convert(input_mesh(), mesh);
      mesh->create_boundary();
      make_metric(mesh, 2);

      //pair of arrays storing the mesh as described in the metis manual
      std::vector<idx_t> eptr;
      std::vector<idx_t> eind;

      eptr.push_back(0);

      //fill eptr and eind, as done in viennamesh plugin "metis", file "mesh_partitionig.cpp"
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

      idx_t num_nodes = mesh->get_number_nodes();
      idx_t num_elements = mesh->get_number_elements();

      idx_t ncommon = mesh->get_number_dimensions();
      idx_t nparts = region_count();

      std::vector<idx_t> epart(num_elements);
      std::vector<idx_t> npart(num_nodes);                //nparts not npart!!!

      idx_t result;

      std::cout << "Calling METIS_PartMeshNodal" << std::endl;

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
      
      /*METIS_PartMeshNodal(&num_elements,
                            &num_nodes,
                            eptr.data(),
                            eind.data(),
                            NULL,
                            NULL,
                            &nparts,
                            NULL,
                            NULL,
                            NULL,
                            epart.data(),
                            npart.data());
*/

      int nloc = 3;   //nlocal is the number of nodes of an element (3 for a triangle in 2D and 4 for a tetrahedron in 3D)
           
      ofstream output;
      output.open("data.txt");

      //TODO: get boundary elements (boundary is really the boundary of the mesh, not the partition interfaces!)
      //   
      //IDEA: use the is_owned_node function from pragmatic to determine partition boundary and build partition interface
      std::vector<std::set<int>> NEList = mesh->get_node_element(); //TODO: replace this, since its unnecessarily copying data!
      std::vector<index_t> _ENList = mesh->get_element_node();  //TODO: replace this, since its unnecessarily copying data!
      std::vector<int> boundary;

      //get partition interface
      std::vector<std::set<index_t>> nodes_per_partition( nparts );

      for (size_t i = 0; i < num_elements; ++i)
      {
        //TODO: the following hast to be updated for the 3D case!!
        nodes_per_partition[ epart[i] ].insert(_ENList[i*3]);
        nodes_per_partition[ epart[i] ].insert(_ENList[i*3+1]);
        nodes_per_partition[ epart[i] ].insert(_ENList[i*3+2]);
      }

      int counter = 0;
/*
      output << "NODES" << std::endl;
      for (auto it : nodes_per_partition)
      {
        output << "Partition " << counter << std::endl;
        for (auto node : it)
        {
          output << node << std::endl;
        }      
        ++counter;
      }
*/
      //get the interface nodes for all partition boundaries
      std::vector<std::set<index_t>> interface_nodes( binomial_coefficient(nparts, 2) ); //We always intersect 2 partitions ==> k=2
      
      //for (size_t i = 0; i < nparts; ++i)
      for (size_t i = 0; i < binomial_coefficient(nparts, 2); ++i)      
      {        
        for (size_t j = i+1; j < nparts; ++j)
        {   
          set_intersection( nodes_per_partition[i].begin(), nodes_per_partition[i].end(),
                            nodes_per_partition[j].begin(), nodes_per_partition[j].end(),
                            inserter(interface_nodes[counter], interface_nodes[counter].begin()) );
          ++counter;
        }      
      }
 
      output << "INTERFACE NODES" << std::endl;
      counter =0;
      for (auto interfaces : interface_nodes)
      {
        //if (interfaces.size() != 0)
          output << "Interface " << counter << std::endl;

        for (auto it : interfaces)
        {         
          //if (interfaces.size() == 0)
            //continue;

          output << it << std::endl;
        } 
        ++counter;
      }
      output << std::endl;
      //end of get partition interface  

      //get boundary
/*
      boundary.resize(num_elements*nloc);
      std::fill(boundary.begin(), boundary.end(), -2);    

      for(size_t i=0; i<num_elements; i++) 
      {
        for(int j=0; j<nloc; j++) 
        {
          //TODO: the following has to be updated for the 3D case!!
          int n1 = _ENList[i*3+(j+1)%3];
          int n2 = _ENList[i*3+(j+2)%3];

          //change this if-clause to check if the node is in the actual partition //TODO is not really working!!!
          //if(npart[n1] == npart[n2]) 
          //if (true)
          //{
            std::set<int> neighbours;
            set_intersection(NEList[n1].begin(), NEList[n1].end(),
                             NEList[n2].begin(), NEList[n2].end(),
                             inserter(neighbours, neighbours.begin()));

            if(neighbours.size()==2) 
            {
              if(*neighbours.begin()==(int)i)
                boundary[i*3+j] = *neighbours.rbegin();
              else
                boundary[i*3+j] = *neighbours.begin();
            }
           //} 

           else 
           {
            // This is a halo facet.
            boundary[i*3+j] = -1;
           }
         }
      }
*/
      //end TODO get boundary
    
      //TODO: get partition interface elements

      //create vector containing the cell indices of each element
      std::vector< std::set<index_t> > elements_per_region(nparts);
      for(size_t i = 0; i < num_elements; i++)
      {
        elements_per_region[ epart[i] ].insert(i);
      }  
  
      output << "NUMBER OF ELEMENTS PER REGION " << std::endl;
      counter = 0;
      for (auto it : elements_per_region)
      {
        output << "Region " << counter << ": " << it.size() << std::endl;
        ++counter;
      }
      output << std::endl;

      //TODO: get partition boundary elements
      counter = 0;
      std::vector< std::set<index_t> > partition_boundary_elements(interface_nodes.size());
      for(auto interfaces : interface_nodes)
      {
        //output << "Boundary No. " << counter << ":" << std::endl;
        for (auto it : interfaces)
        {            
          //output << it << ": " << std::endl;
          for(auto NE_it : NEList[it])
          {
            //output << "   " << NE_it << std::endl;
            partition_boundary_elements[counter].insert(NE_it);
          }
        } 
        ++counter;
      }

      output << "PARTITION INTERFACE ELEMENTS" << std::endl;  
      counter = 0;
      for (auto it : partition_boundary_elements)
      {
        output << "Interface " << counter << " containing " << it.size() << " elements" << std::endl;
        for (auto it2 : it)
        {
          output << it2 << std::endl;
        }

        ++counter;
      }
      output << std::endl;
  
      output << "DUPLICATED INTERFACE ELEMENTS" << std::endl;
      std::set<index_t> degenerate_elements;
      for (size_t i = 0; i < partition_boundary_elements.size(); ++i)
      {
        std::vector<index_t> temp;
        for (size_t j = i+1; j < partition_boundary_elements.size(); ++j)
        {
          
          set_intersection( partition_boundary_elements[i].begin(), partition_boundary_elements[i].end(),
                            partition_boundary_elements[j].begin(), partition_boundary_elements[j].end(),
                            inserter(temp, temp.begin()) );

          for (auto it : temp)
          {
            degenerate_elements.insert(it);
          }
        }
      }

      for (auto it : degenerate_elements)
      {
        output << it << std::endl;
      }
      output << std::endl;

      //TODO: remove partition interface elements from the actual partitions
      for (size_t i = 0; i < partition_boundary_elements.size(); ++i)
      {
        for (size_t j = 0; j < elements_per_region.size(); ++j)
        {
          for (auto it = partition_boundary_elements[i].begin(); it != partition_boundary_elements[i].end(); ++it)
          { 
            //std::cout << *it << std::endl; 
            elements_per_region[j].erase(*it);
          }
        }
      }

      //for developing and debugging purposes only
      output << "ELEMENTS PER REGION AFTER REMOVING INTERFACE ELEMENTS " << std::endl;
      counter = 0;
      for (auto it : elements_per_region)
      {
        output << "Region " << counter << " now consists of " << it.size() << " elements" << std::endl;
        for (auto it2 : it)
        {
          output << it2 << std::endl;
        }
        ++counter;
      }
      output << std::endl;

      //count element appearances
      
      std::vector<size_t> num_of_appearances(num_elements, 0);

      for (size_t i = 0; i < num_elements; ++i)
      {
        size_t appearances = 0;
        for (auto it : elements_per_region)
        {
          num_of_appearances[i] += count(it.begin(), it.end(), i);
        }    

        for (auto it : partition_boundary_elements)
        {
          num_of_appearances[i] += count(it.begin(), it.end(), i);
        }
      }

      output << "NUMBER OF APPEARANCES" << std::endl;

      for (size_t i = 0; i < num_of_appearances.size(); ++i)
      {
        output << i << ": " << num_of_appearances[i] << std::endl;
      } 
      output << "Accumulate: " << accumulate(num_of_appearances.begin(), num_of_appearances.end(), 0) << std::endl;

      output << std::endl;
      //end!
   
      //output << std::endl;
      //end of get partition boundary elements
      
      //for debugging and test purposes only
      output << "OLD COORDINATES" << std::endl;
      for (size_t i = 0; i<num_nodes; ++i)
      {
        double p[2];
        mesh->get_coords(i, p);
        output << i << ": " << p[0] << ", " << p[1] << std::endl;
      }
      output << std::endl; //end of for debugging and test purposes only

      //laplace smoother
/*
      std::cout << "Plugin" << std::endl;
      std::cout << &mesh << std::endl;
*/
      laplace_smoother(mesh);

      //end of laplace smoother

      //for debugging and test purposes only
      output << "NEW COORDINATES" << std::endl;
      for (size_t i = 0; i < num_nodes; ++i)
      {
        double p[2];
        mesh->get_coords(i, p);
        output << i << ": " << p[0] << ", " << p[1] << std::endl;
      }
      output << std::endl; //end of for debugging and test purposes only
      
      //close output file (used for debugging and developing)
      output.close();

      //write pragmatic data file
      std::string filename;
      filename += "examples/data/pragmatic_metis_partitioning";
  
      VTKTools<double>::export_vtu(filename.c_str(), mesh);

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

      //COPIED FROM MESH_PARTITIONING.CPP INCLUDED IN METIS PLUGIN
      //ONLY HERE FOR DEBUGGING AND DEVELOPMENT PURPOSES
      data_handle<bool> multi_mesh_output = get_input<bool>("multi_mesh_output");
      mesh_handle output_mesh = make_data<mesh_handle>();

      ConstCellRangeType cells( input_mesh() );

      typedef viennagrid::result_of::element_copy_map<>::type ElementCopyMapType;

      if ( multi_mesh_output.valid() && multi_mesh_output() )
      {
        output_mesh.resize( region_count() );
        std::vector<ElementCopyMapType> copy_maps;

        for (int i = 0; i != region_count(); ++i)
        {
          ElementCopyMapType copy_map( output_mesh(i) );
          for (ConstCellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
          {
            int part = epart[(*cit).id().index()];
            if (part == i)
              copy_map( *cit );
          }
        }
      }
      else
      {
        ElementCopyMapType copy_map( output_mesh(), false );

        for (ConstCellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
        {
          ElementType cell = copy_map( *cit );
          viennagrid::add( output_mesh().get_or_create_region(epart[(*cit).id().index()]), cell );
        }
      }

      set_output( "mesh", output_mesh );
      //END OF DEBUGGING AND DEVELOPMENT CODE SNIPPET COPIED FROM METIS PLUGIN
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//


      return true;
    } //end of bool pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
}
