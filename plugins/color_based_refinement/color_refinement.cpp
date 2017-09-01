#include "color_refinement.hpp"
#include "pragmatic_mesh.hpp"

namespace viennamesh
{
		color_refinement::color_refinement()	{}
		std::string color_refinement::name() {return "color_refinement";}

		bool color_refinement::run(viennamesh::algorithm_handle &)
		{
			//mesh_handle input_mesh = get_required_input<mesh_handle>("mesh");
			data_handle<pragmatic_wrapper::mesh> input_mesh = get_required_input<pragmatic_wrapper::mesh>("mesh");
			data_handle<int> num_partitions = get_required_input<int>("num_partitions");
			data_handle<int> num_threads = get_input<int>("num_threads");
			data_handle<bool> single_mesh_output = get_input<bool>("single_mesh_output");
			string_handle input_file = get_input<string_handle>("filename");
			string_handle algorithm = get_input<string_handle>("algorithm");
			
			info(1) << name() << std::endl;

			return VIENNAMESH_SUCCESS;
		}
}
