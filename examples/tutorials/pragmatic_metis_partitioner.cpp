/*
/ PRAGMATIC_METIS_PARTITIONING
/
/ 
/
/ Developer: 	Lukas Gnam
/		Institute for Microelectronics
/		Faculty of Electrical Engineering and Information Engineering
/		Vienna University of Technology
/		gnam@iue.tuwien.ac.at
/		2016
/
/
/
*/

#include <iostream>
#include <string>
#include <vector>
#include <cfloat>
#include <deque>
#include <cassert>
#include <cmath>

//ViennaMesh includes
#include "viennameshpp/core.hpp"

int main(int argc, char **argv)
{
	// Check and read console input 
	int region_count = 0;
	std::string filename;

	if (argc < 3)
	{
		std::cout << "Parameters missing!" << std::endl;
		std::cout << "Correct use of parameters: <input_file> <region_count> " << std::endl;
		return -1;
	}

	else if (!argv[1])
	{
		std::cout << "Parameter <input_file> missing!" << std::endl;
		return -1;
	}

	else if (!atoi(argv[2]))
	{
		std::cout << "Parameter <region_count> missing!" << std::endl;
		return -1;
	}

	else
	{
		filename = argv[1];
		region_count = atoi(argv[2]);
	}

	// Create context handle
	viennamesh::context_handle context;

	// Create algorithm handle for reading the mesh from a file and run it
	viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
	mesh_reader.set_input("filename", filename.c_str());
	mesh_reader.run();

  //create partitions on the pragmatic data structure
  viennamesh::algorithm_handle mesh_partitioner = context.make_algorithm("pragmatic_metis_partitioner");
  mesh_partitioner.set_default_source(mesh_reader);
  mesh_partitioner.set_input("filename", filename.c_str());
  mesh_partitioner.set_input("region_count", region_count);
  mesh_partitioner.set_input("multi_mesh_output", false);
  mesh_partitioner.run();
/*
  //read merged mesh
  viennamesh::algorithm_handle merged_mesh_reader = context.make_algorithm("mesh_reader");
  merged_mesh_reader.set_input("filename", "examples/data/pragmatic_metis_partitioning/merged_mesh.vtu");
  merged_mesh_reader.run();

  //get mesh quality statistics
  viennamesh::algorithm_handle stats = context.make_algorithm("make_statistic");
  stats.set_default_source(merged_mesh_reader);
  //only for benchmarking!!
  std::string input_filename;

  size_t pos = filename.find_last_of("/\\");
  size_t pos2 = filename.find(".vtu");
  input_filename = filename.substr(pos+1, pos2-(pos+1));
  //input_filename += "_";
  //input_filename += std::to_string(region_count/2);
  input_filename += ".txt";      

  stats.set_input("filename", input_filename);   
  //end of only for benchmarking!!           
  stats.set_input("metric_type", "min_angle");
  stats.set_input("good_element_threshold", 0.35);
  stats.run();
  stats.set_input("metric_type", "radius_ratio");
  stats.set_input("good_element_threshold", 2);  
  stats.run();

  viennamesh::algorithm_handle stats_input_mesh = context.make_algorithm("make_statistic");
  stats_input_mesh.set_default_source(mesh_reader);
  stats_input_mesh.set_input("metric_type", "min_angle");
  stats_input_mesh.set_input("good_element_threshold", 0.35);
  //for benchmarks only
  stats_input_mesh.set_input("filename", input_filename);
  //end of for benchmarks only
  stats_input_mesh.run();
  stats_input_mesh.set_input("metric_type", "radius_ratio");
  stats_input_mesh.set_input("good_element_threshold", 2 );
  stats_input_mesh.run();

/*	
	//Write mesh
	viennamesh::algorithm_handle write_merged_mesh = context.make_algorithm("mesh_writer");
	write_merged_mesh.set_default_source(mesh_partitioner);
	write_merged_mesh.set_input("filename", "/home/lgnam/Desktop/software/ViennaMesh/viennamesh-dev/build/examples/data/pragmatic_metis_partitioning/parpartmesh.vtu");
	write_merged_mesh.run();
*/  

  return 0;
}
