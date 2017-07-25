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
/*
	viennamesh::algorithm_handle partitioning = context.make_algorithm("metis_mesh_partitioning");
	partitioning.set_default_source(mesh_reader);
	partitioning.set_input("region_count", 	region_count);
	partitioning.set_input("multi_mesh_output", true);
	partitioning.run();
/*
	viennamesh::algorithm_handle boundary = context.make_algorithm("extract_boundary");
	boundary.set_default_source(partitioning);
	boundary.run();

	viennamesh::algorithm_handle geo = context.make_algorithm("extract_plc_geometry");
	geo.set_default_source(boundary);
	geo.set_input("coplanar_tolerance", 1e-4);
	geo.set_input("coplinear_tolerance", 1e-4);
	geo.run();

	viennamesh::algorithm_handle tetgen = context.make_algorithm("tetgen_make_mesh");
	tetgen.set_default_source(geo);
	tetgen.run();

	std::cout << "Tetgen successful" << std::endl;

  */
  	// Create algorithm handle for the color-based-refinement
	viennamesh::algorithm_handle color = context.make_algorithm("color_refinement");    
	color.set_default_source(mesh_reader);
	color.set_input("algorithm", algorithm);
	color.set_input("options", options);
	color.set_input("num_partitions", region_count);
	color.set_input("filename", filename.c_str());
	color.set_input("num_threads", num_threads);
	color.set_input("single_mesh_output", true);
	color.run();
//*/
/*
	viennamesh::algorithm_handle triangle = context.make_algorithm("triangle_make_mesh");
	triangle.set_default_source(color);
	triangle.set_input("colors", color.get_output("colors"));
	triangle.set_input("threads", num_threads);
	//triangle.set_input("option_string", "zp");
	triangle.set_input("cell_size", 0.49);
	triangle.set_input("no_points_on_boundary", true);
	triangle.run();

	std::cout << "Triangle successful" << std::endl;
/*
	viennamesh::algorithm_handle merger = context.make_algorithm("merge_meshes");
	merger.set_default_source(triangle);
	merger.run();
*/
/*
	//Write output mesh
	viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
	mesh_writer.set_default_source(color);
	
	//construct filename
	
	std::string folder = "test/";

	std::string outfilename = filename.substr(filename.find_last_of("/")+1);
	//outfilename.replace(outfilename.find(".vtu"), 12, "_refined.vtu");
	outfilename.replace(outfilename.find(".vtu"), 12, "_refined_");
	outfilename += std::to_string(region_count);
	outfilename += "parts_";
	outfilename += std::to_string(num_threads);
	outfilename += "threads_own_metis.vtu";

	folder += outfilename;

	mesh_writer.set_input("filename", folder.c_str());
	mesh_writer.run();
//*/

/*	viennamesh::algorithm_handle mesh_merger = context.make_algorithm("merge_meshes");
	mesh_merger.set_default_source(color);
	mesh_merger.run();

	viennamesh::algorithm_handle write_merged = context.make_algorithm("mesh_writer");
	write_merged.set_default_source(mesh_merger);
	write_merged.set_input("filename", "test/triangle_merged.vtu");
	write_merged.run();
//*/
    return 0;
}