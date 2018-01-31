#include "viennameshpp/core.hpp"
#include <boost/algorithm/string.hpp>  
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <ctime>
#include <chrono>

struct policy_placement_s
{
    std::string lind="lindstrom-turk";
    std::string mid ="mid-point";
}policy_placement;      //erzeugt structs

struct policy_cost_s
{
    std::string lind="lindstrom-turk";
    std::string edge="edge-length";

}policy_cost;           //erzeugt structs


viennamesh::algorithm_handle        coarser;
viennamesh::algorithm_handle        statistic;
viennamesh::algorithm_handle        mesh_reader;
viennamesh::algorithm_handle        mesh_writer;
viennamesh::context_handle          context;

double ratio                        = 0;
double ratio_step                   = 0.03;
double ratio_start                  = 0.01;
double ratio_stop                   = 0.99;

double lindstrom_volume_weight          = 0;
double lindstrom_volume_weight_step     = 0.3;
double lindstrom_volume_weight_start    = 0;
double lindstrom_volume_weight_stop     = 1;

double lindstrom_boundary_weight        = 0; 
double lindstrom_boundary_weight_step   = 0.3; 
double lindstrom_boundary_weight_start  = 0; 
double lindstrom_boundary_weight_stop   = 1; 

double lindstrom_shape_weight           = 0; 
double lindstrom_shape_weight_step      = 0.3; 
double lindstrom_shape_weight_start     = 0; 
double lindstrom_shape_weight_stop      = 1; 


std::ofstream  file_ofstream;


void coarse(std::string         dest,
            std::string         cost,
            std::string         placement);


int main(int argc, char ** argv)
{
    if(argc != 2 )
    {
        std::cout << "please define input mesh ";
        exit(1);
    }

    std::string file_source         (argv[1]);
    std::string file_destination    =file_source.substr(0,file_source.find_last_of('.')) + "/";
    std::string file_csv            =file_destination + "Time.csv";

   
 
    mkdir(file_source.substr(0,file_source.find_last_of('.')).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    file_ofstream.open(file_csv);


    //mesh reader
    std::cout << "pre  reader\n";
    mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", file_source );
    mesh_reader.run();
    std::cout << "post reader\n";
    //mid edge

    file_ofstream << policy_cost.edge << ";" << policy_placement.mid << "\n";
    file_ofstream << "ratio;volume;boundary;shape;time(milliseconds)\n";
    for(ratio=ratio_start;ratio<=ratio_stop;ratio+=ratio_step)
    {
        coarse(file_destination,policy_cost.edge,policy_placement.mid);
    }

   
    file_ofstream << policy_cost.lind << ";" << policy_placement.mid << "\n";
    file_ofstream << "ratio;volume;boundary;shape;time(milliseconds)\n";
    for(ratio=ratio_start;ratio<=ratio_stop;ratio+=ratio_step)
    {
        for(lindstrom_volume_weight= lindstrom_volume_weight_start;lindstrom_volume_weight<=lindstrom_volume_weight_stop;lindstrom_volume_weight+=lindstrom_volume_weight_step)
            for(lindstrom_boundary_weight=lindstrom_boundary_weight_start; lindstrom_boundary_weight<=lindstrom_boundary_weight_stop;lindstrom_boundary_weight+=lindstrom_boundary_weight_step)
                for(lindstrom_shape_weight=lindstrom_shape_weight_start;lindstrom_shape_weight<=lindstrom_shape_weight_stop;lindstrom_shape_weight+=lindstrom_shape_weight_step)
                    coarse(file_destination,policy_cost.lind,policy_placement.mid);
    }

    
    file_ofstream << policy_cost.edge << ";" << policy_placement.lind << "\n";
    file_ofstream << "ratio;volume;boundary;shape;time(milliseconds)\n";
    for(ratio=ratio_start;ratio<=ratio_stop;ratio+=ratio_step)
    {
        for(lindstrom_volume_weight= lindstrom_volume_weight_start;lindstrom_volume_weight<=lindstrom_volume_weight_stop;lindstrom_volume_weight+=lindstrom_volume_weight_step)
            for(lindstrom_boundary_weight=lindstrom_boundary_weight_start; lindstrom_boundary_weight<=lindstrom_boundary_weight_stop;lindstrom_boundary_weight+=lindstrom_boundary_weight_step)
                for(lindstrom_shape_weight=lindstrom_shape_weight_start;lindstrom_shape_weight<=lindstrom_shape_weight_stop;lindstrom_shape_weight+=lindstrom_shape_weight_step)
                    coarse(file_destination,policy_cost.edge,policy_placement.lind);
    }

    
    file_ofstream << policy_cost.lind << ";" << policy_placement.lind << "\n";
    file_ofstream << "ratio;volume;boundary;shape;time(milliseconds)\n";
    for(ratio=ratio_start;ratio<=ratio_stop;ratio+=ratio_step)
    {
        for(lindstrom_volume_weight= lindstrom_volume_weight_start;lindstrom_volume_weight<=lindstrom_volume_weight_stop;lindstrom_volume_weight+=lindstrom_volume_weight_step)
            for(lindstrom_boundary_weight=lindstrom_boundary_weight_start; lindstrom_boundary_weight<=lindstrom_boundary_weight_stop;lindstrom_boundary_weight+=lindstrom_boundary_weight_step)
                for(lindstrom_shape_weight=lindstrom_shape_weight_start;lindstrom_shape_weight<=lindstrom_shape_weight_stop;lindstrom_shape_weight+=lindstrom_shape_weight_step)
                    coarse(file_destination,policy_cost.lind,policy_placement.lind);
    }



    file_ofstream.close();
    std::cout << "Finished\n";

    return 0;
}

void coarse(std::string         dest,
            std::string         cost,
            std::string         placement)
{
    //time 
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
   
    //coarse
    std::cout << "pre  coarse\n";
    coarser = context.make_algorithm("cgal_mesh_simplification");
    coarser.set_default_source(mesh_reader);
    coarser.set_input("stop_predicate", "ratio");
    coarser.set_input("ratio",           ratio);

    coarser.set_input("cost_policy",         cost );
    coarser.set_input("placement_policy",    placement );

    coarser.set_input("lindstrom_volume_weight",    lindstrom_volume_weight);    
    coarser.set_input("lindstrom_boundary_weight",  lindstrom_boundary_weight);  
    coarser.set_input("lindstrom_shape_weight",     lindstrom_shape_weight); 
    coarser.run();
    std::cout << "post coarse\n";    

    //time
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    file_ofstream << ratio << ";" << lindstrom_volume_weight << ";" <<  lindstrom_boundary_weight << ";" << lindstrom_shape_weight<<";" 
                    << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "\n";

    //statistic    
    statistic = context.make_algorithm("cgal_statistic");
    statistic.set_input("original mesh",mesh_reader.get_output("mesh"));
    statistic.set_input("coarse mesh",coarser.get_output("mesh"));
    statistic.set_input("default quantities",true);

    statistic.run();



    //mesh writer
    mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(statistic);
    std::string file_output=dest + cost      + "_" + placement + "_ratio_" + std::to_string(ratio) + 
                                        "_volume_"      + std::to_string(lindstrom_volume_weight) + 
                                        "_boundary_"    + std::to_string(lindstrom_boundary_weight) +
                                        "_shape_"       + std::to_string(lindstrom_shape_weight) +
                                        ".vtu";

    mesh_writer.set_input( "filename", file_output);
    mesh_writer.run();
}
