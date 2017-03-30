#include "color_refinement.hpp"

#include "pragmatic_mesh.hpp"
#include "mesh_partitions.hpp"

namespace viennamesh
{
		color_refinement::color_refinement()	{}
		std::string color_refinement::name() {return "color_refinement";}

		bool color_refinement::run(viennamesh::algorithm_handle &)
		{	
			data_handle<pragmatic::pragmatic_mesh> output_mesh = get_required_input<pragmatic::pragmatic_mesh>("mesh");
			
			info(1) << name() << std::endl;

			info(1) << "  Number of vertices: " << output_mesh()->get_number_nodes() << std::endl;

			MeshPartitions InputMesh();

			set_output("mesh", output_mesh);
			return true;
		} //end run()		
}