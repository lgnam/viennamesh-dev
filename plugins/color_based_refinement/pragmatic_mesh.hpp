#ifndef VIENNAMESH_ALGORITHM_PRAGMATIC_MESH_HPP
#define VIENNAMESH_ALGORITHM_PRAGMATIC_MESH_HPP

//viennamesh includes
#include "viennameshpp/plugin.hpp"

#include "Mesh.h"

namespace viennamesh
{

  namespace pragmatic
  {
    typedef Mesh<double>* pragmatic_mesh;
  }

  //inline viennamesh_error pragmatic_make_mesh(pragmatic::pragmatic_mesh * mesh)

  inline viennamesh_error pragmatic_delete_mesh(pragmatic::pragmatic_mesh mesh)
  {
    std::cerr << std::endl << "pragmatic delete mesh" << std::endl;
    std::cerr << mesh << std::endl;

    std::cerr << mesh->get_number_nodes() << std::endl;

    //delete mesh;

    return VIENNAMESH_SUCCESS;
  }  

  inline viennamesh_error delete_pragmatic_mesh(viennamesh_data data)
  {
    return delete_viennamesh_data<pragmatic::pragmatic_mesh>(data, pragmatic_delete_mesh);
  }

  viennamesh_error convert(viennagrid::mesh const & input, pragmatic::pragmatic_mesh & output);
  viennamesh_error convert(pragmatic::pragmatic_mesh const & input, viennagrid::mesh & output);

  template<>
  viennamesh_error internal_convert<viennagrid_mesh, pragmatic::pragmatic_mesh>(viennagrid_mesh const & input, pragmatic::pragmatic_mesh & output);
  template<>
  viennamesh_error internal_convert<pragmatic::pragmatic_mesh, viennagrid_mesh>(pragmatic::pragmatic_mesh const & input, viennagrid_mesh & output);
 

  namespace result_of
  {
    template<>
    struct data_information<pragmatic::pragmatic_mesh>
    {
      static std::string type_name()
      {
        return "pragmatic_mesh";
      }

      static viennamesh_data_make_function make_function()
      {
        return viennamesh::generic_make<pragmatic::pragmatic_mesh>;
        //return make_pragmatic_mesh;
      }

      static viennamesh_data_delete_function delete_function()
      {
        return viennamesh::generic_delete<pragmatic::pragmatic_mesh>;
        //return delete_pragmatic_mesh;
      }
    };
  }
}
#endif