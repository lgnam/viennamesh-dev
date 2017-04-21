//ViennaMesh includes
#include "viennameshpp/core.hpp"

int main(int argc, char *argv[])
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
  
  // Create algorithm handle for the color-based-refinement
	viennamesh::algorithm_handle color = context.make_algorithm("color_refinement");
//for (region_count =2; region_count <= 8192; region_count*=2)
//{
    
	color.set_default_source(mesh_reader);
	color.set_input("num_partitions", region_count);
	color.set_input("filename", filename.c_str());
	color.run();
//}
	//Write output mesh
	viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
	mesh_writer.set_default_source(color);
	mesh_writer.set_input("filename", "output_mesh.vtu");
	//mesh_writer.run();

    return 0;
}
