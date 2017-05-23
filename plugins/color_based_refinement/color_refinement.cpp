#include "color_refinement.hpp"

#include "pragmatic_mesh.hpp"
#include "mesh_partitions.hpp"

//#include <omp.h>

namespace viennamesh
{
		color_refinement::color_refinement()	{}
		std::string color_refinement::name() {return "color_refinement";}

		bool color_refinement::run(viennamesh::algorithm_handle &)
		{	
			data_handle<pragmatic::pragmatic_mesh> input_mesh = get_required_input<pragmatic::pragmatic_mesh>("mesh");
			data_handle<int> num_partitions = get_required_input<int>("num_partitions");
			data_handle<int> num_threads = get_input<int>("num_threads");
			string_handle input_file = get_input<string_handle>("filename");
		
			info(1) << name() << std::endl;
			
			size_t found = input_file().find_last_of("/");
			info(1) << input_file().substr(found+1) << std::endl;

			info(1) << "  Number of vertices: " << input_mesh()->get_number_nodes() << std::endl;
      		info(1) << "  Dimension: " << input_mesh()->get_number_dimensions() << std::endl;
			info(1) << "  Threads: " << num_threads() << std::endl;

			MeshPartitions InputMesh(input_mesh(), num_partitions(), input_file().substr(found+1), num_threads());
		
			//SERIAL PART
			auto overall_tic = std::chrono::system_clock::now();

			auto wall_tic = std::chrono::system_clock::now();
				InputMesh.MetisPartitioning();
			std::chrono::duration<double> partitioning_duration = std::chrono::system_clock::now() - wall_tic;
			viennamesh::info(1) << "  Partitioning time " << partitioning_duration.count() << std::endl;

			wall_tic = std::chrono::system_clock::now();
				InputMesh.CreateNeighborhoodInformation();
			std::chrono::duration<double> adjacency_duration = std::chrono::system_clock::now() - wall_tic;
			viennamesh::info(1) << "  Creating adjacency information time " << adjacency_duration.count() << std::endl;

			wall_tic = std::chrono::system_clock::now();
				InputMesh.ColorPartitions();
			std::chrono::duration<double> coloring_duration = std::chrono::system_clock::now() - wall_tic;
			viennamesh::info(1) << "  Coloring time " << coloring_duration.count() << std::endl;
/*
			std::chrono::duration<double> pragmatic_duration;
			double refinement_duration;
			double pragmatic_adjacency;
			double map_duration;
			double metric_duration;
			//if (num_threads() == 1)
			{
				wall_tic = std::chrono::system_clock::now();
					InputMesh.CreatePragmaticDataStructures_ser(); //Think about where I create the actual data structures!!!
				pragmatic_duration = std::chrono::system_clock::now() - wall_tic;
			}
/*
			else
			{*/
			std::chrono::duration<double> pragmatic_duration;
			std::vector<double> timer_log;
			
			std::vector<double> times(13);
			std::fill(times.begin(), times.end(), 0.0);
			
			std::vector<double> l2g_build, l2g_access, g2l_build, g2l_access;

				//PARALLEL PART
				wall_tic = std::chrono::system_clock::now();
					InputMesh.CreatePragmaticDataStructures_par(timer_log, times, l2g_build, l2g_access, g2l_build, g2l_access);
				pragmatic_duration = std::chrono::system_clock::now() - wall_tic;	

			//}

/*
			wall_tic = std::chrono::system_clock::now();
				//InputMesh.RefineInterior();
				MeshPartitionsRefinement RefineMesh(InputMesh, 0);
			std::chrono::duration<double> refinement_duration = std::chrono::system_clock::now() - wall_tic;
*/
			std::chrono::duration<double> overall_duration = std::chrono::system_clock::now() - overall_tic;

			int r_vertices {0};
			int r_elements {0};
			
			InputMesh.GetRefinementStats(&r_vertices, &r_elements);

			ofstream csv;
			csv.open("times.csv", ios::app);

			//csv << "File, Threads, Vertices, Elements, Desired Partitions, Created Partitions, Colors, Metis [s], Adjacency Info [s], 
			//Coloring [s], Parallel DSs [s], Prep [s], Nodes [s], g2l [s], l2g [s], Coords [s], ENList [s], new Mesh [s], Boundary [s], Metric [s],
			// Update Metric [s], Interface Check [s],  Refine [s], Create Refine [s], R-Vertices, R-Elements, Total [s], Thread Times in Color Loop [s]" << std::endl;
			csv << input_file().substr(found+1) << ", " << num_threads() << ", " << input_mesh()->get_number_nodes() << ", ";
			csv << input_mesh()->get_number_elements() << ", "  << num_partitions() << ", " << InputMesh.get_max()+1;
			csv << ", " << InputMesh.get_colors() << ", ";
			csv << std::fixed << std::setprecision(8) << partitioning_duration.count() << ", ";
			csv << adjacency_duration.count() << ", ";
			csv << coloring_duration.count() << ", ";
			/*
			for (size_t i = 0; i < l2g_build.size(); ++i)
			{
				csv << l2g_build[i] << ", ";
			}

			for (size_t i = 0; i < g2l_build.size(); ++i)
			{
				csv << g2l_build[i] << ", ";
			}

			for (size_t i = 0; i < l2g_access.size(); ++i)
			{
				csv << l2g_access[i] << ", ";
			}

			for (size_t i = 0; i < g2l_access.size(); ++i)
			{
				csv << g2l_access[i] << ", ";
			}
			*/
			csv << pragmatic_duration.count() << ", ";
			
			for (size_t i =0; i < times.size(); ++i)
				csv << times[i] << ", ";

			csv << r_vertices << ", ";
			csv << r_elements << ", ";
			csv << overall_duration.count() << ",";
			
			for (size_t i =0; i < timer_log.size(); ++i)
				csv << timer_log[i] << ", ";

			csv << std::endl;
			csv.close();
			
			//InputMesh.WritePartitions();
			//InputMesh.WriteMergedMesh("output.vtu");

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
