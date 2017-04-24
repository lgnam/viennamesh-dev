#include "color_refinement.hpp"

#include "pragmatic_mesh.hpp"
#include "mesh_partitions.hpp"
#include "mesh_partitions_refinement.hpp"

//#include <omp.h>

namespace viennamesh
{
		color_refinement::color_refinement()	{}
		std::string color_refinement::name() {return "color_refinement";}

		bool color_refinement::run(viennamesh::algorithm_handle &)
		{	
			data_handle<pragmatic::pragmatic_mesh> input_mesh = get_required_input<pragmatic::pragmatic_mesh>("mesh");
			data_handle<int> num_partitions = get_required_input<int>("num_partitions");
			string_handle input_file = get_input<string_handle>("filename");
			
			info(1) << name() << std::endl;
			
			size_t found = input_file().find_last_of("/");
			info(1) << input_file().substr(found+1) << std::endl;

			info(1) << "  Number of vertices: " << input_mesh()->get_number_nodes() << std::endl;
      		info(1) << "  Dimension: " << input_mesh()->get_number_dimensions() << std::endl;

			MeshPartitions InputMesh(input_mesh(), num_partitions(), input_file().substr(found+1));
		
			//SERIAL PART
			auto overall_tic = std::chrono::system_clock::now();

			auto wall_tic = std::chrono::system_clock::now();
				InputMesh.MetisPartitioning();
			std::chrono::duration<double> partitioning_duration = std::chrono::system_clock::now() - wall_tic;
			//viennamesh::info(1) << "  Partitioning time " << wall_clock_duration.count() << std::endl;

			wall_tic = std::chrono::system_clock::now();
				InputMesh.CreateNeighborhoodInformation();
			std::chrono::duration<double> adjacency_duration = std::chrono::system_clock::now() - wall_tic;
			//viennamesh::info(1) << "  Creating adjacency information time " << wall_clock_duration.count() << std::endl;

			wall_tic = std::chrono::system_clock::now();
				InputMesh.ColorPartitions();
			std::chrono::duration<double> coloring_duration = std::chrono::system_clock::now() - wall_tic;
			//viennamesh::info(1) << "  Coloring time " << wall_clock_duration.count() << std::endl;

			//viennamesh::info(1) << "  Overall time inside " << overall_duration.count() << std::endl;
			wall_tic = std::chrono::system_clock::now();
				InputMesh.CreatePragmaticDataStructures_ser(); //Think about where I create the actual data structures!!!
			std::chrono::duration<double> pragmatic_duration = std::chrono::system_clock::now() - wall_tic;

			/*//PARALLEL PART
			wall_tic = std::chrono::system_clock::now();
				InputMesh.CreatePragmaticDataStructures_par();
			std::chrono::duration<double> pragmatic_duration = std::chrono::system_clock::now() - wall_tic;	
*/

			wall_tic = std::chrono::system_clock::now();
				//InputMesh.RefineInterior();
				MeshPartitionsRefinement RefineMesh(InputMesh);
			std::chrono::duration<double> refinement_duration = std::chrono::system_clock::now() - wall_tic;

			std::chrono::duration<double> overall_duration = std::chrono::system_clock::now() - overall_tic;

			ofstream csv;
			csv.open("times.csv", ios::app);

			//csv << "File, Vertices, Elements, Desired Partitions, Created Partitions, Colors, Partitioning [s], Adjacency Info [s],
			// Coloring [s] ,Create Partitions [s], Refinement [s], Total [s]" << std::endl;
			csv << input_file().substr(found+1) << ", " << input_mesh()->get_number_nodes() << ", " << input_mesh()->get_number_elements() << ", ";
			csv << num_partitions() << ", " << InputMesh.get_max()+1 << ", " << InputMesh.get_colors() << ", ";
			csv << partitioning_duration.count() << ", ";
			csv << adjacency_duration.count() << ", ";
			csv << coloring_duration.count() << ", ";
			csv << pragmatic_duration.count() << ", ";
			csv << refinement_duration.count() << ", ";
			csv << overall_duration.count() << std::endl;
			csv.close();
		
			//InputMesh.WritePartitions();
			/*
			auto overall_tic = std::chrono::system_clock::now();

			auto wall_tic = std::chrono::system_clock::now();
				MetisPartitioning(original_mesh, num_regions);
			std::chrono::duration<double> partitioning_duration = std::chrono::system_clock::now() - wall_tic;
			//viennamesh::info(1) << "  Partitioning time " << wall_clock_duration.count() << std::endl;

			wall_tic = std::chrono::system_clock::now();
				CreateNeighborhoodInformation(original_mesh, num_regions);
			std::chrono::duration<double> adjacency_duration = std::chrono::system_clock::now() - wall_tic;
			//viennamesh::info(1) << "  Creating adjacency information time " << wall_clock_duration.count() << std::endl;

			wall_tic = std::chrono::system_clock::now();
				ColorPartitions(num_regions);
			std::chrono::duration<double> coloring_duration = std::chrono::system_clock::now() - wall_tic;
			//viennamesh::info(1) << "  Coloring time " << wall_clock_duration.count() << std::endl;

			std::chrono::duration<double> overall_duration = std::chrono::system_clock::now() - overall_tic;
			//viennamesh::info(1) << "  Overall time inside " << overall_duration.count() << std::endl;

			CreatePragmaticDataStructures_ser(original_mesh, num_regions); //Think about where I create the actual data structures!!!

			ofstream csv;
			csv.open("times.csv", ios::app);

			//csv << "File, Vertices, Elements, Partitions, Colors, Partitioning [s], Adjacency Info [s], Coloring [s], Total [s]" << std::endl;
			csv << filename << ", " << original_mesh->get_number_nodes() << ", " << original_mesh->get_number_elements() << ", ";
			csv << num_regions << ", " << max+1 << ", " << colors << ", ";
			csv <<  partitioning_duration.count() << ", ";
			csv <<  adjacency_duration.count() << ", ";
			csv <<  coloring_duration.count() << ", ";
			csv <<  overall_duration.count() << std::endl;
			csv.close();
		
			//WritePartitions();
			*/

			set_output("mesh", input_mesh);
			return true;
		} //end run()		
}
