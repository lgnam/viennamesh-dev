#ifndef VIENNAMESH_ALGORITHM_TRIANGLE_GENERATOR_HPP
#define VIENNAMESH_ALGORITHM_TRIANGLE_GENERATOR_HPP

#include "viennamesh/core/algorithm.hpp"
#include "viennamesh/algorithm/triangle/mesh.hpp"

namespace viennamesh
{
  namespace triangle
  {

    class Algorithm : public BaseAlgorithm
    {
    public:

      string name() const { return "Triangle 1.6 mesher"; }

      bool run_impl()
      {
        viennamesh::result_of::const_parameter_handle<triangle::InputMesh>::type input_mesh = inputs.get<triangle::InputMesh>("default");

        if (!input_mesh)
        {
          error(1) << "Input Parameter 'default' (type: mesh) is missing or of non-convertable type" << std::endl;
          return false;
        }

        std::ostringstream options;
        options << "zp";

        ConstDoubleParameterHandle min_angle = inputs.get<double>("min_angle");
        if (min_angle)
          options << "q" << min_angle->get() / M_PI * 180.0;


        ConstDoubleParameterHandle cell_size = inputs.get<double>("cell_size");
        if (cell_size)
          options << "a" << cell_size->get();

        ConstBoolParameterHandle delaunay = inputs.get<bool>("delaunay");
        if (delaunay && delaunay->get())
          options << "D";


        ConstStringParameterHandle algorithm_type = inputs.get<string>("algorithm_type");
        if (algorithm_type)
        {
          if (algorithm_type->get() == "incremental_delaunay")
            options << "i";
          else if (algorithm_type->get() == "sweepline")
            options << "F";
          else if (algorithm_type->get() == "devide_and_conquer")
          {}
          else
          {
            warning(5) << "Algorithm not recognized: '" << algorithm_type->get() << "' supported algorithms:" << std::endl;
            warning(5) << "  'incremental_delaunay'" << std::endl;
            warning(5) << "  'sweepline'" << std::endl;
            warning(5) << "  'devide_and_conquer'" << std::endl;
          }
        }


        triangulateio tmp = input_mesh->get().mesh;

        REAL * tmp_regionlist = NULL;
        REAL * tmp_holelist = NULL;

        typedef viennamesh::result_of::const_parameter_handle<SeedPoint2DContainer>::type ConstSeedPointContainerHandle;
        ConstSeedPointContainerHandle seed_points_handle = inputs.get<SeedPoint2DContainer>("seed_points");
        if (seed_points_handle && !seed_points_handle->get().empty())
        {
          info(5) << "Found seed points" << std::endl;
          SeedPoint2DContainer const & seed_points = seed_points_handle->get();

          tmp_regionlist = (REAL*)malloc( 4*sizeof(REAL)*(tmp.numberofregions+seed_points.size()) );
          memcpy( tmp_regionlist, tmp.regionlist, 4*sizeof(REAL)*tmp.numberofregions );

          for (unsigned int i = 0; i < seed_points.size(); ++i)
          {
            tmp_regionlist[4*(i+tmp.numberofregions)+0] = seed_points[i].first[0];
            tmp_regionlist[4*(i+tmp.numberofregions)+1] = seed_points[i].first[1];
            tmp_regionlist[4*(i+tmp.numberofregions)+2] = REAL(seed_points[i].second);
            tmp_regionlist[4*(i+tmp.numberofregions)+3] = 0;
          }

          tmp.numberofregions += seed_points.size();
          tmp.regionlist = tmp_regionlist;

          options << "A";
        }


        typedef viennamesh::result_of::const_parameter_handle<Point2DContainer>::type ConstPointContainerHandle;
        ConstPointContainerHandle hole_points_handle = inputs.get<Point2DContainer>("hole_points");
        if (hole_points_handle && !hole_points_handle->get().empty())
        {
          info(5) << "Found hole points" << std::endl;
          Point2DContainer const & hole_points = hole_points_handle->get();

          tmp_holelist = (REAL*)malloc( 2*sizeof(REAL)*(tmp.numberofholes+hole_points.size()) );
          memcpy( tmp_holelist, tmp.holelist, 2*sizeof(REAL)*tmp.numberofholes );

          for (std::size_t i = 0; i < hole_points.size(); ++i)
          {
            tmp_holelist[2*(tmp.numberofholes+i)+0] = hole_points[i][0];
            tmp_holelist[2*(tmp.numberofholes+i)+1] = hole_points[i][1];
          }

          tmp.numberofholes += hole_points.size();
          tmp.holelist = tmp_holelist;
        }



        char * buffer = new char[options.str().length()];
        std::strcpy(buffer, options.str().c_str());

        viennamesh::result_of::parameter_handle<triangle::OutputMesh>::type output_mesh = make_parameter<triangle::OutputMesh>();

        viennautils::StdCapture capture;
        capture.start();

        triangulate( buffer, &tmp, &output_mesh->get().mesh, NULL);

        capture.finish();
        info(5) << capture.get() << std::endl;

        delete[] buffer;
        if (seed_points_handle && !seed_points_handle->get().empty())
          free(tmp_regionlist);
        if (hole_points_handle && !hole_points_handle->get().empty())
          free(tmp_holelist);

        outputs.set( "default", output_mesh );

        return true;
      }

    private:
    };

  }


}

#endif
