#include "viennameshpp/core.hpp"
#include <iostream>



int main(int argc, char **argv)
{
	viennamesh::algorithm_handle        coarser;
	viennamesh::algorithm_handle        statistic;
	viennamesh::algorithm_handle        mesh_reader;
	viennamesh::algorithm_handle        mesh_writer;
	viennamesh::context_handle          context;
	
	std::cout << "pre  reader\n";
    mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", "../data/elephant.vtu" );
    mesh_reader.run();
    std::cout << "post reader\n";
    
	std::cout << "pre  coarse\n";
    coarser = context.make_algorithm("cgal_mesh_simplification");
    coarser.set_default_source(mesh_reader);
    coarser.set_input("stop_predicate", "ratio");
    coarser.set_input("ratio",           0.5);

    coarser.set_input("cost_policy",         "lindstrom-turk" );
    coarser.set_input("placement_policy",    "lindstrom-turk");
 
    coarser.run();
    std::cout << "post coarse\n";    
    
    //until here it is basically the same as the cgal mesh simplification plugin
    
    statistic = context.make_algorithm("cgal_statistic");				//normally load the plugin
    statistic.set_input("original mesh",mesh_reader.get_output("mesh"));//set the original mesh
    statistic.set_input("coarse mesh",coarser.get_output("mesh"));		//set the coarsed mesh
    statistic.set_input("default quantities",true);						//use the defautl quantities

    statistic.run();													//run the plugin
    
    
    
    mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(statistic);

    mesh_writer.set_input( "filename", "out.vtu");
    mesh_writer.run();
    return 1;
}
