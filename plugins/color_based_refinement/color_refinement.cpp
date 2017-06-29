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
			string_handle algorithm = get_input<string_handle>("algorithm");
		
			info(1) << name() << std::endl;
			
			size_t found = input_file().find_last_of("/");
			info(1) << input_file().substr(found+1) << std::endl;

			info(1) << "  Number of vertices: " << input_mesh()->get_number_nodes() << std::endl;
      		info(1) << "  Dimension: " << input_mesh()->get_number_dimensions() << std::endl;
			info(1) << "  Threads: " << num_threads() << std::endl;
			
			std::string algo;
			std::string options;
			if (!algorithm.valid())
			{
				algo = "pragmatic";
			}	

			else if (algorithm() == "pragmatic" || algorithm() == "triangle")
			{
				algo = algorithm();	
				
				if (algo == "triangle")
				{
					string_handle tri_options = get_input<string_handle>("options");
					options = tri_options();
				}
			}	

			else 
			{
				viennamesh::error(1) << "'" << algorithm() << "'" << " is not a valid algorithm!" << std::endl;
				return false;
			}					

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
			std::vector<double> threads_log;
			std::vector<double> refine_times;
			std::vector<double> triangulate_log;
			std::vector<double> int_check_log;
			std::vector<double> build_tri_ds;
			/*
			std::vector<double> times(13);
			std::fill(times.begin(), times.end(), 0.0);
			*/
			std::vector<double> l2g_build, l2g_access, g2l_build, g2l_access;

				//PARALLEL PART
				wall_tic = std::chrono::system_clock::now();
					InputMesh.CreatePragmaticDataStructures_par(threads_log, refine_times, l2g_build, l2g_access, g2l_build, g2l_access, 
																algo, options, triangulate_log, int_check_log, build_tri_ds);
				std::chrono::duration<double> cpds_dur = std::chrono::system_clock::now() - wall_tic;	

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
			
			InputMesh.GetRefinementStats(&r_vertices, &r_elements, algo);

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
			
			csv << cpds_dur.count() << ", ";
			csv << r_vertices << ", ";
	 		csv << r_elements << ", ";
			csv << overall_duration.count() << ",";

			for (size_t i =0; i < refine_times.size(); ++i)
				csv << refine_times[i] << ", ";

			for (size_t i = 0; i < triangulate_log.size(); ++i)
			{
				csv << triangulate_log[i] << ", ";
			}
			
			for (size_t i =0; i < threads_log.size(); ++i)
				csv << threads_log[i] << ", ";

			for (size_t i = 0; i < int_check_log.size(); ++i)
				csv << int_check_log[i] << ",";
/*
			for (size_t i = 0; i < build_tri_ds.size(); ++i)
				csv << build_tri_ds[i] << ",";
/*
			csv << std::endl;
			csv.close();
	*/		
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
			//Convert to ViennaGrid
			auto convert_tic = std::chrono::system_clock::now();

			//ViennaGrid typedefs
			typedef viennagrid::mesh                                                        MeshType;
			typedef viennagrid::result_of::element<MeshType>::type                          VertexType;

			mesh_handle output_mesh = make_data<mesh_handle>();
			int vertices = 0;
			int elements = 0;

			output_mesh.resize(num_partitions());

			//TODO: REPLACE WITH TEMPLATED VERSION!!!
			//convert pragmatic to viennagrid output
			if (algo=="pragmatic")
			{
				std::cout << "Converting Pragmatic to ViennaGrid data structure" << std::endl;

				for (size_t i = 0; i != InputMesh.pragmatic_partitions.size(); ++i)
				{
					//get basic mesh information
					size_t NNodes = InputMesh.pragmatic_partitions[i]->get_number_nodes();
					size_t NElements = InputMesh.pragmatic_partitions[i]->get_number_elements(); 

					//create empty vector of size NNodes containing viennagrid vertices
					std::vector<VertexType> vertex_handles(NNodes);
					
					//iterating all pragmatic vertices and store their coordinates in the viennagrid vertices
					for(size_t j = 0; j < NNodes; ++j)
					{
						vertex_handles[j] = viennagrid::make_vertex( output_mesh(i), viennagrid::make_point(InputMesh.pragmatic_partitions[i]->get_coords(j)[0], 
																	InputMesh.pragmatic_partitions[i]->get_coords(j)[1]) );
					} //end of for loop iterating all pragmatic vertices 

					//iterating all pragmatic elements and createg their corresponding viennagrid triangles
					for (size_t j = 0; j < NElements; ++j)
					{
						index_t const* ENList_ptr = InputMesh.pragmatic_partitions[i]->get_element(j); 
						viennagrid::make_triangle( output_mesh(i), vertex_handles[*(ENList_ptr++)], vertex_handles[*(ENList_ptr++)], vertex_handles[*(ENList_ptr++)] );
					} //end of iterating all pragmatic elements

					//Update total number of vertices and elements
					vertices += InputMesh.pragmatic_partitions[i]->get_number_nodes();
					elements += InputMesh.pragmatic_partitions[i]->get_number_elements();
				}
			} //end of convert pragmatic to viennagrid output

			//convert triangle to viennagrid output
			else
			{
				//std::cout << "Converting Triangle to ViennaGrid data structure" << std::endl;
				
				for (size_t i = 0; i != InputMesh.triangle_partitions.size(); ++i)
				{
					//get basic mesh information
					size_t NNodes = InputMesh.triangle_partitions[i].numberofpoints;
					size_t NElements = InputMesh.triangle_partitions[i].numberoftriangles;

					//create empty vector of size NNodes containing viennagrid vertices
					std::vector<VertexType> vertex_handles(NNodes);

					//iterating all triangle vertices and store their coordinates in the viennagrid vertices
					for(size_t j = 0; j < NNodes; ++j)
					{
						vertex_handles[j] = viennagrid::make_vertex( output_mesh(i), viennagrid::make_point(InputMesh.triangle_partitions[i].pointlist[2*j], 
																	InputMesh.triangle_partitions[i].pointlist[2*j+1]));
					} //end of for loop iterating all pragmatic vertices 

					//iterating all triangle elements and createg their corresponding viennagrid triangles
					for (size_t j = 0; j < NElements; ++j)
					{
						viennagrid::make_triangle( output_mesh(i), 
												   vertex_handles[InputMesh.triangle_partitions[i].trianglelist[3*j]], 
												   vertex_handles[InputMesh.triangle_partitions[i].trianglelist[3*j+1]], 
												   vertex_handles[InputMesh.triangle_partitions[i].trianglelist[3*j+2]] );
					} //end of iterating all triangle elements

					//Update total number of vertices and elements
					vertices += InputMesh.triangle_partitions[i].numberofpoints;
					elements += InputMesh.triangle_partitions[i].numberoftriangles;
				}		

				//free used memory
				for (size_t i = 0; i < output_mesh.size(); ++i)
				{
					free(InputMesh.triangle_partitions[i].pointlist);
					free(InputMesh.triangle_partitions[i].pointattributelist);
					free(InputMesh.triangle_partitions[i].pointmarkerlist);
					//free(InputMesh.triangle_partitions[i].numberofpoints);
					//free(InputMesh.triangle_partitions[i].numberofpointattributes);

					free(InputMesh.triangle_partitions[i].trianglelist);
			/*		//free(InputMesh.triangle_partitions[i].triangleattributelist);
					//free(InputMesh.triangle_partitions[i].trianglearealist);
					//free(InputMesh.triangle_partitions[i].neighborlist);
					//free(InputMesh.triangle_partitions[i].numberoftriangles);
					//free(InputMesh.triangle_partitions[i].numberofcorners);
					//free(InputMesh.triangle_partitions[i].numberoftriangleattributes);

					free(InputMesh.triangle_partitions[i].segmentlist);
					free(InputMesh.triangle_partitions[i].segmentmarkerlist);
					//free(InputMesh.triangle_partitions[i].numberofsegments);

					free(InputMesh.triangle_partitions[i].holelist);
					//free(InputMesh.triangle_partitions[i].numberofholes);

					free(InputMesh.triangle_partitions[i].regionlist);
					//free(InputMesh.triangle_partitions[i].numberofregions);

					free(InputMesh.triangle_partitions[i].edgelist);
					free(InputMesh.triangle_partitions[i].edgemarkerlist);
					free(InputMesh.triangle_partitions[i].normlist);
					//free(InputMesh.triangle_partitions[i].numberofedges);
*/
					//delete InputMesh.pragmatic_partitions[i];
				}
				//end of free used memory*/
			} //end of convert triangle to viennagrid output

			std::chrono::duration<double> convert_dur = std::chrono::system_clock::now() - convert_tic;

			csv << convert_dur.count();
			csv << std::endl;
			csv.close();
/*
			ofstream color_parts;
			color_parts.open("color_parts.txt");

			std::vector<std::vector<int>> col_parts = InputMesh.get_color_partitions();

			for (size_t i = 0; i < InputMesh.get_colors(); ++i)
			{
				for (size_t j = 0; j < col_parts[i].size(); ++j)
				{
					color_parts << col_parts[i][j] << " ";
				}
				color_parts << std::endl;
			}
			
			color_parts.close();
/*
			for (size_t i = 0; i < col_parts.size(); ++i)
			{
				for (size_t j = 0; j < col_parts[i].size(); ++j)
				{
					std::cout << col_parts[i][j] << " ";
				}
				std::cout << std::endl;
			}
*/
			viennamesh::info(1) << "Number of Vertices: " << vertices << std::endl;
			viennamesh::info(1) << "Number of Elements: " << elements << std::endl;
			set_output("mesh", output_mesh);
			set_output("colors", InputMesh.get_colors());

			return true;
		} //end run()		
}
