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
#include <unordered_map>
#include <time.h>

#include <chrono>

//TODO: use boost bimap to evade using two unordered maps per mesh for local-to-global- and global-to-local-indexing
#include <boost/bimap.hpp>

#include "grouped_partitions.hpp"
#include "grouped_partitions_smooth.hpp"
#include "grouped_partitions_refinement.hpp"

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
    //pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
    //
    //Uses metis to partition the mesh (represented in pragmatic data structure)
    bool pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
		{
      std::cout << name() << std::endl;

      //create mesh_handle to read input mesh			
		  mesh_handle input_mesh = get_required_input<mesh_handle>("mesh");
      string_handle filename = get_input<string_handle>("filename");
            
      //create data_handle for optional inputs
		  data_handle<int> region_count = get_input<int>("region_count");
/*
      output << "Benchmarking " << filename() << " using " << region_count()/2 << " threads " << std::endl;
      output  << "==================================================" << std::endl;
*/  /*
      //output-file
      ofstream output;     
      size_t pos = filename().find_last_of("/\\");
      size_t pos2 = filename().find(".vtu");
      std::string input_filename = filename().substr(pos+1, pos2-(pos+1));
      //input_filename += "_";
      //input_filename += std::to_string(region_count()/2);
      input_filename += ".txt";           
      //if file is not there yet, write header

      if (!std::ifstream(input_filename.c_str()))
      {
        output.open(input_filename.c_str(), ofstream::app);
        output << "# Threads | Version | Vertex Count | CPU time [s] | Wall-clock time [s] | MinAngleRatioMerged | RadiusRatioRatioMerged | MinAngleRatioOriginal| RadiusRatioRatioOriginal |";
      }

      else
      {        
        output.open(input_filename.c_str(), ofstream::app);
      }     
/*
      //convert viennamesh into pragmatic mesh data structure
      Mesh<double> *mesh = nullptr;

		  mesh = convert(input_mesh(), mesh);
      mesh->create_boundary(); 

      //make_metric(mesh, 2); //it is not necessary to create a metric!
      //output << "SimpleLaplaceOnGroups_sections" << std::endl << "==================================================" << std::endl;
      output << std::endl << region_count()/2 << " sections ";
      output << mesh->get_number_nodes() << " "; 
      GroupedPartitions Mesh1(mesh, region_count());  
      clock_t tic = clock();
      GroupedPartitionsSmooth Smoother1(Mesh1);
      clock_t toc = clock();
      //output << "Create Smoothing Object: " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;

      auto wall_tic = std::chrono::system_clock::now();      
      tic = clock();
      Smoother1.SimpleLaplaceOnGroups_sections(2);
      std::chrono::duration<double> wall_clock_duration = std::chrono::system_clock::now() - wall_tic;
      toc = clock();
      //ProfilerStop();
      //output << "SimpleLaplaceOnGroups: " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;
      output << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << " " << wall_clock_duration.count();

      //Smoother1.Evaluate();

      tic = clock();
      Mesh1.WriteMergedMesh();
      toc = clock();

      delete mesh;

/*
      //output << "WriteMergedMesh: " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl  << "==================================================" << std::endl;

      //convert viennamesh into pragmatic mesh data structure
      Mesh<double> *serial_mesh = nullptr;

		  serial_mesh = convert(input_mesh(), serial_mesh);
      serial_mesh->create_boundary();

      //output << "Serial" << std::endl << "==================================================" << std::endl;
      output << std::endl << "1" << " serial ";      
      output << serial_mesh->get_number_nodes() << "  ";
      GroupedPartitions Mesh2(serial_mesh, region_count());
      clock_t tic = clock();
      GroupedPartitionsSmooth Smoother2(Mesh2);
      clock_t toc = clock();
      //output << "Create Smoothing Object: " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;
      
      auto wall_tic = std::chrono::system_clock::now();
      tic = clock();
      Smoother2.SimpleLaplace(2);
      std::chrono::duration<double> wall_clock_duration = std::chrono::system_clock::now() - wall_tic;
      toc = clock();
      //output << "SimpleLaplace: " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl;
      output << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << " " << wall_clock_duration.count();

      tic = clock();
      //Mesh2.WriteMergedMesh("examples/data/pragmatic_metis_partitioning/SimpleLaplace.vtu");
      Mesh2.WriteMergedMesh();
      toc = clock();
      //output << "WriteMergedMesh: " << static_cast<double>(toc - tic) / CLOCKS_PER_SEC << std::endl << std::endl;

      //delete pointer created at the beginning of pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
      delete serial_mesh;
//*/
      //Mesh1.PrintQuality();

      Mesh<double> *mesh_ref = nullptr;
      mesh_ref = convert(input_mesh(), mesh_ref);
      mesh_ref->create_boundary();

      //debug
      ofstream nelist_original;
      nelist_original.open("NEList_original.txt");
      std::vector<std::set<int>> temp_nelist;
      temp_nelist = mesh_ref->get_node_element();
      int ne_elements2 = 0;
      for (size_t i = 0; i < mesh_ref->get_number_nodes(); ++i)
      {
        nelist_original << i << ": ";
        for (auto it : temp_nelist[i])
        {
          nelist_original << " " << it << " ";
          ++ne_elements2;
        }
        nelist_original << std::endl;
      }
      nelist_original.close();
      
      std::cout << "NEElements_original: " << ne_elements2 << std::endl;
      //end of debug

      GroupedPartitions MeshRefine(mesh_ref, region_count());
      GroupedPartitionsRefinement Refiner(MeshRefine, 1);  

      //delete pointer created at the beginning of pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
      delete mesh_ref;   
//*/

     // output.close();

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

      set_output("mesh", input_mesh);
      
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

      return true;
    } //end of bool pragmatic_metis_partitioner::run(viennamesh::algorithm_handle &)
}
