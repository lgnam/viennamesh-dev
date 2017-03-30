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
      }

      static viennamesh_data_delete_function delete_function()
      {
        return viennamesh::generic_delete<pragmatic::pragmatic_mesh>;
      }
    };
  }
}
#endif