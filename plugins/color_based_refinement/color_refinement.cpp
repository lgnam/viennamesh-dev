#include "color_refinement.hpp"

#include "pragmatic_mesh.hpp"
#include "mesh_partitions.hpp"

namespace viennamesh
{
		color_refinement::color_refinement()	{}
		std::string color_refinement::name() {return "color_refinement";}

		bool color_refinement::run(viennamesh::algorithm_handle &)
		{	
			data_handle<pragmatic::pragmatic_mesh> input_mesh = get_required_input<pragmatic::pragmatic_mesh>("mesh");
			data_handle<int> num_partitions = get_required_input<int>("num_partitions");
			
			info(1) << name() << std::endl;

			info(1) << "  Number of vertices: " << input_mesh()->get_number_nodes() << std::endl;

			MeshPartitions InputMesh(input_mesh(), num_partitions());

			set_output("mesh", input_mesh);
			return true;
		} //end run()		
}