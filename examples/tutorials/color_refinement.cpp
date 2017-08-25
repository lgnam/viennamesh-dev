//ViennaMesh includes
#include "viennameshpp/core.hpp"
#include "../../plugins/color_based_refinement/outbox.hpp"

int Outbox::no_outboxes = 0;

int main(int argc, char *argv[])
{
    // Check and read console input 
	int region_count = 0;
	int num_threads = 0;
	std::string filename;
	std::string algorithm;
	std::string options;

	if (argc < 4)
	{
		std::cout << "Parameters missing!" << std::endl;
		std::cout << "Correct use of parameters: <input_file> <region_count> <number_threads> <algorithm>" << std::endl;
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

	else if (!atoi(argv[3]))
	{
		std::cout << "Parameter <num_threads> missing!" << std::endl;
		return -1;
	}

	else
	{
		filename = argv[1];
		region_count = atoi(argv[2]);
		num_threads = atoi(argv[3]);
	}

	if (argv[4])
	{
		algorithm = argv[4];
	}

	if (argv[5])
	{
		options = argv[5];
	}

    // Create context handle
	viennamesh::context_handle context;

	// Create algorithm handle for reading the mesh from a file and run it
	viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
	mesh_reader.set_input("filename", filename.c_str());
	mesh_reader.run();

	//Investigation of memory-leak
	viennamesh::algorithm_handle color = context.make_algorithm("color_refinement");
	color.set_default_source(mesh_reader);
	color.run();
/*
  	// Create algorithm handle for the color-based-refinement
	viennamesh::algorithm_handle color = context.make_algorithm("color_refinement");    
	color.set_default_source(mesh_reader);
	color.set_input("algorithm", algorithm);
	color.set_input("options", options);
	color.set_input("num_partitions", region_count);
	color.set_input("filename", filename.c_str());
	color.set_input("num_threads", num_threads);
	color.set_input("single_mesh_output", false);
	color.run(); //*/
//*/


	//Write output mesh
	viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
	mesh_writer.set_default_source(color);
	
	//construct filename
	
	std::string folder = "test/metis/500x500/dual/partgraphkway/";

	std::string outfilename = filename.substr(filename.find_last_of("/")+1);
	/*//outfilename.replace(outfilename.find(".vtu"), 12, "_initial.vtp");
	outfilename.replace(outfilename.find(".vtu"), 17, "_single_initial_");
	outfilename += std::to_string(num_threads);
	outfilename+= "threads_";
	outfilename+= std::to_string(region_count);
	outfilename+= "parts.vtu"; //*/

	outfilename.replace(outfilename.find(".vtu"), 15, "_partgraphkway_");
	outfilename += std::to_string(num_threads);
	outfilename+= "threads_";
	outfilename+= std::to_string(region_count);
	outfilename+= "parts.vtu";

	folder += outfilename;

	mesh_writer.set_input("filename", folder.c_str());
	//mesh_writer.run();
//*/
    return 0;
}