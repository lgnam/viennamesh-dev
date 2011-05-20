/* =============================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
                 ViennaMesh - The Vienna Mesh Generation Framework
                             -----------------

   authors:    Josef Weinbub                         weinbub@iue.tuwien.ac.at


   license:    see file LICENSE in the base directory
============================================================================= */

#include <iostream>
#include <vector>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "viennautils/config.hpp"
#include "viennautils/convert.hpp"
#include "viennautils/contio.hpp"
#include "viennautils/file.hpp"
#include "viennautils/io/bnd.hpp"
#include "viennautils/io/hin.hpp"
#ifdef VIENNAMESH_HAVE_GTSIO
#include "viennautils/io/gts.hpp"
#endif

#include "viennagrid/domain.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/io/vtk_reader.hpp"
#include "viennagrid/io/gau_reader.hpp"
#include "viennagrid/algorithm/cell_normals.hpp"

#include "viennamesh/interfaces/cervpt.hpp"
#include "viennamesh/interfaces/netgen.hpp"
#include "viennamesh/interfaces/tetgen.hpp"
#include "viennamesh/wrapper.hpp"
#include "viennamesh/adaptors/orienter.hpp"
#include "viennamesh/adaptors/cell_normals.hpp"
#include "viennamesh/adaptors/hull_quality.hpp"

#include <boost/any.hpp> // removeme
#include <boost/type_traits/remove_pointer.hpp>

namespace viennamesh {

namespace key {
std::string geometry = "geometry";
} // end namespace key

namespace query {

struct input
{
   template<typename ConfigT>
   static std::string type(ConfigT const& config)
   {
      return config.query("/input/type/key/text()");
   }
};


} // end namespace query
} // end namespace viennamesh



int main(int argc, char *argv[])
{
   if(argc != 4)
   {
      std::cerr << "## Error::Parameter - usage: " << argv[0] << " inputfile outputfile configfile" << std::endl;
      std::cerr << "## shutting down .." << std::endl;
      return -1;
   }
   
   std::string inputfile(argv[1]);
   std::string outputfile(argv[2]);
   
   std::string input_extension  = viennautils::file_extension(inputfile);
   std::string output_extension = viennautils::file_extension(outputfile);

   typedef viennautils::config<viennautils::tag::xml>::type    config_type;
   config_type config;
   std::ifstream configstream(argv[3]);
   config.read(configstream);
   configstream.close();


//   if(input_extension == "bnd")
//   {
//      if(viennamesh::query::input::type(config) == viennamesh::key::geometry)
//      {
//         viennautils::io::bnd_reader my_bnd_reader;
//         my_bnd_reader(inputfile); 
// 
//         typedef viennamesh::wrapper<viennamesh::tag::bnd, viennautils::io::bnd_reader>      bnd_wrapper_type;
//         bnd_wrapper_type                 wrapped_data(my_bnd_reader);      

//         typedef viennamesh::result_of::mesh_generator<viennamesh::tag::cervpt>::type        cervpt_mesh_generator_type;
//         cervpt_mesh_generator_type       mesher;       

//         typedef viennamesh::result_of::mesh_adaptor<viennamesh::tag::orienter>::type        orienter_adaptor_type;
//         orienter_adaptor_type            orienter;
//         
//         typedef viennamesh::result_of::mesh_adaptor<viennamesh::tag::cell_normals>::type    cellnormals_adaptor_type;
//         cellnormals_adaptor_type         cell_normals;         

//         typedef cervpt_mesh_generator_type::result_type       result_type;
//         result_type result = cell_normals(mesher(wrapped_data));
//         //result_type result = cell_normals(orienter(mesher(wrapped_data)));

//         viennagrid::io::export_vtk(*result, outputfile);
//      }
//   }
//   else
   if(input_extension == "hin")
   {
      if(viennamesh::query::input::type(config) == viennamesh::key::geometry)
      {
         viennautils::io::hin_reader my_hin_reader;
         my_hin_reader(inputfile);

         typedef viennamesh::wrapper<viennamesh::tag::hin, viennautils::io::hin_reader>      hin_wrapper_type;
         hin_wrapper_type                    wrapped_data(my_hin_reader);      

         typedef viennamesh::result_of::mesh_generator<viennamesh::tag::cervpt>::type        cervpt_hull_mesh_generator_type;
         cervpt_hull_mesh_generator_type     hull_mesher;      

         typedef viennamesh::result_of::mesh_adaptor<viennamesh::tag::orienter>::type        orienter_adaptor_type;
         orienter_adaptor_type               orienter;

         typedef viennamesh::result_of::mesh_adaptor<viennamesh::tag::cell_normals>::type    cell_normals_adaptor_type;
         cell_normals_adaptor_type           cell_normals;         
         
         typedef viennamesh::result_of::mesh_adaptor<viennamesh::tag::hull_quality>::type    hull_quality_adaptor_type;
         hull_quality_adaptor_type           hull_quality;                  

         typedef cervpt_hull_mesh_generator_type::result_type       hull_result_type;

         hull_result_type hull_mesh = hull_quality(cell_normals(orienter(hull_mesher(wrapped_data))));
//         hull_result_type hull_mesh = cell_normals(orienter(hull_mesher(wrapped_data)));

         typedef hull_result_type::value_type                                                      hull_domain_type;
         typedef hull_domain_type::config_type                                                     hull_domain_configuration_type;
         typedef viennagrid::result_of::ncell_type<hull_domain_configuration_type, hull_domain_configuration_type::cell_tag::topology_level>::type     hull_cell_type;
         
         viennagrid::io::vtk_writer<hull_domain_type>  my_hull_vtk_writer;
         
         my_hull_vtk_writer.add_cell_data_normal(
            viennagrid::io::io_data_accessor_segment_based<
               hull_cell_type, viennagrid::seg_cell_normal_tag, viennagrid::seg_cell_normal_data::type
            >(viennagrid::seg_cell_normal_tag()), "cell_normals");
             
         my_hull_vtk_writer.writeDomain(*hull_mesh, "hull_mesh.vtu");

//         typedef viennamesh::result_of::mesh_generator<viennamesh::tag::netgen>::type        netgen_volume_mesh_generator_type;
//         netgen_volume_mesh_generator_type       volume_mesher;      

//         typedef netgen_volume_mesh_generator_type::result_type       volume_result_type;
//         volume_result_type volume_mesh = volume_mesher(hull_mesh);
//         
//         typedef volume_result_type::value_type                                               volume_domain_type;         
//         
//         viennagrid::io::vtk_writer<volume_domain_type>  my_volume_vtk_writer;         
//         my_volume_vtk_writer.writeDomain(*volume_mesh, "volume_mesh.vtu");
      }
   }
//   else
//   if(input_extension == "gau32")
//   { 
//      typedef viennagrid::domain<viennagrid::config::triangular_3d>        domain_type;
//      domain_type domain;
//      
//      viennagrid::io::importGAU(domain, inputfile);      
//      
//      typedef viennamesh::wrapper<viennamesh::tag::viennagrid, domain_type>     gau_wrapper_type;
//      gau_wrapper_type wrapped_data(domain);      
//      
//      typedef viennamesh::result_of::mesh_generator<viennamesh::tag::tetgen>::type   netgen_mesh_generator_type;
//      netgen_mesh_generator_type mesher;      

//      typedef netgen_mesh_generator_type::result_type        netgen_result_type;
//      netgen_result_type result = mesher(wrapped_data);         

//      viennagrid::io::export_vtk(*result, outputfile);
//   }
//   else
//   if(input_extension == "gts")
//   {
//      typedef viennagrid::domain<viennagrid::config::line_2d>                  domain_type;
//      domain_type domain;
//      
//      viennautils::io::gts_reader my_gts_reader;
//      my_gts_reader(domain, inputfile);      
//      
//      typedef viennamesh::wrapper<viennamesh::tag::viennagrid, domain_type>     vgrid_wrapper_type;
//      vgrid_wrapper_type data_in(domain);      
//      
//      typedef viennamesh::result_of::mesh_generator<viennamesh::tag::triangle, vgrid_wrapper_type>::type   mesh_generator_type;
//      mesh_generator_type mesher(data_in);      

//      mesher( boost::fusion::make_map<viennamesh::tag::criteria, viennamesh::tag::size>(viennamesh::tag::conforming_delaunay(), 1.0) );         

//      typedef viennagrid::domain<viennagrid::config::triangular_2d> domain_out_type;
//      domain_out_type domain_out;      
//      
//      typedef viennamesh::transfer<viennamesh::tag::viennagrid>      transfer_type;
//      transfer_type  transfer;
//      transfer(mesher, domain_out);      
//      
//      viennagrid::io::Vtk_writer<domain_out_type> my_vtk_writer;
//      my_vtk_writer.writeDomain(domain_out, outputfile);
//   }
//   else
//   {
//      std::cerr << "## input file format not supported: " << input_extension << std::endl;
//      std::cerr << "## shutting down .. " << std::endl;
//      return -1;   
//   }

   return 0;
}



