#ifndef __cgal_statistic
#define __cgal_statistic
#include <CGAL/Polyhedron_3.h>
#define __enable_libigl 1

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#define CGAL_EIGEN3_ENABLED
#include <CGAL/Monge_via_jet_fitting.h>

#include "viennameshpp/plugin.hpp"

#include <vector>

namespace viennamesh
{
    class cgal_statistic : public plugin_algorithm{
        public:
            typedef CGAL::Simple_cartesian<double> Kernel;
            typedef CGAL::Polyhedron_3<Kernel> mesh_t;
            typedef mesh_t::Point_3 Point_3;
            typedef Kernel::Vector_3 Vector_3;
            typedef mesh_t::Face_handle Face;
            typedef mesh_t::Facet Facet_t;
            typedef CGAL::Monge_via_jet_fitting<Kernel>::Monge_form Monge_form;

            typedef CGAL::AABB_face_graph_triangle_primitive<mesh_t> Primitive;
            typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
            typedef CGAL::AABB_tree<Traits> Tree;

            typedef enum quantity_dimention_e
            {
                Vertex=0,
                Edge,
                Facet,
                Point=0,
                Line,
                Area
            } quantity_dimention_t;

            typedef enum quantity_compare_e
            {
                none,
                closest_Point_value,
                closest_Facet_value,
                closest_Points_value,
                closest_Facets_value,
                closest_Facet,
                tree_orig,
                tree_coarse
            } quantity_compare_t;

            typedef union quantity_Funktions_u
            {
                double  (*Funktion_1)(Point_3 a,Point_3 b,Point_3 c);
                double  (*Funktion_2)(Point_3 point,mesh_t& mesh);
                double  (*Funktion_3)(Facet_t& face,std::vector<Facet_t> facets);  

                quantity_Funktions_u(){memset(this,0,sizeof(*this));}
                
                quantity_Funktions_u(double  (*t_Funktion_1)(Point_3 a,Point_3 b,Point_3 c)             ):Funktion_1(t_Funktion_1){}
                quantity_Funktions_u(double  (*t_Funktion_2)(Point_3 point,mesh_t& mesh)                ):Funktion_2(t_Funktion_2){}
                quantity_Funktions_u(double  (*t_Funktion_3)(Facet_t& face,std::vector<Facet_t> facets) ):Funktion_3(t_Funktion_3){}
                
                double  operator()  (Point_3 a,Point_3 b,Point_3 c) {return (*Funktion_1)(a,b,c);}
                double  operator()  (Point_3 point,mesh_t& mesh)    {return (*Funktion_2)(point,mesh);} 
                double  operator()  (Facet_t& face,std::vector<Facet_t> facets){return (*Funktion_3)(face,facets);}  
            } quantity_Funktions_t;

            typedef union compare_Funktions_u
            {
                double (*Funktion_1)(Point_3 point,Tree& tree);
                double (*Funktion_2)(double orig,double corse);
                double (*Funktion_3)(std::vector<double> orig,double corse);

                compare_Funktions_u(){memset(this,0,sizeof(*this));}

                compare_Funktions_u(double  (*t_Funktion_1)(Point_3 point,Tree& tree)):Funktion_1(t_Funktion_1){}
                compare_Funktions_u(double  (*t_Funktion_2)(double orig,double corse)):Funktion_2(t_Funktion_2){}
                compare_Funktions_u(double  (*t_Funktion_3)(std::vector<double> orig,double corse)):Funktion_3(t_Funktion_3){}

                double operator()   (Point_3 point,Tree& tree)  {return (*Funktion_1)(point,tree);}
                double operator()   (double orig,double corse)  {return (*Funktion_2)(orig,corse);}
                double operator()   (std::vector<double> orig,double corse){return (*Funktion_3)(orig,corse);}

            }compare_Funktions_t;

            typedef union point_Funktions_u
            {
                Point_3(*Funktion_1)(Point_3 a,Point_3 b,Point_3 c);

                point_Funktions_u(){memset(this,0,sizeof(*this));}

                point_Funktions_u(Point_3(*t_Funktion_1)(Point_3 a,Point_3 b,Point_3 c)):Funktion_1(t_Funktion_1){}
            
                Point_3 operator()  (Point_3 a,Point_3 b,Point_3 c) {return (*Funktion_1)(a,b,c);}
            }point_Funktions_t;

            typedef struct statistic_data_s
            {
                typedef std::vector<double> data_t;
               
                std::string             name;
                data_t                  data;
                quantity_dimention_t    dimention;
                quantity_compare_t      compare;

                point_Funktions_t point_Funktion;
                compare_Funktions_t compare_Funktion;
                quantity_Funktions_t quantity_Funktion;

                statistic_data_s(   std::string t_name,
                                    quantity_dimention_t    t_dimention,
                                    quantity_Funktions_t    t_quantity_Funktion =quantity_Funktions_t(),
                                    quantity_compare_t      t_compare           =quantity_compare_t::none,
                                    compare_Funktions_t     t_compare_Funktion  =compare_Funktions_t(),
                                    point_Funktions_t       t_point_Funktion    =point_Funktions_t()
                                    ):  name(t_name),dimention(t_dimention),
                                        quantity_Funktion(t_quantity_Funktion),
                                        compare(t_compare),compare_Funktion(t_compare_Funktion),
                                        point_Funktion(t_point_Funktion){}
                
            }statistic_data_t;

            //constructor
            cgal_statistic(mesh_t starting_mesh):old_mesh(starting_mesh){}
            cgal_statistic(){}
            virtual ~cgal_statistic();
        

            //plugin stuff
            static std::string name(){return "cgal_statistic";}
            bool run(viennamesh::algorithm_handle &);
            
            void make_statistic();                                         //ändern
            void make_quantity_fields();
            void default_quantities();
            void add_quantity(statistic_data_t data);

            viennamesh::data_handle<viennagrid_quantity_field> make_data_handle();

            mesh_t& get_old_mesh(){return old_mesh;}
            mesh_t& get_new_mesh(){return new_mesh;}

            //static statisic funktions
            static double quality                       (Point_3 a,Point_3 b,Point_3 c);
            static double quality_diffrence             (double orig,double corse);
            static double relative_quality              (double orig,double corse);
            static double relative_quality_diffrence    (double orig,double corse);
            static double diffrence_of_mean             (std::vector<double> orig,double corse);

            static Point_3 middel_Point                 (Point_3 a,Point_3 b,Point_3 c);
           
            static double squared_distance              (Point_3 point,Tree& tree);
            
            static double mean_curvature                (Point_3 point,mesh_t& mesh);
            static double gaus_curvature                (Point_3 point,mesh_t& mesh);
            static double min_curvature                 (Point_3 point,mesh_t& mesh);
            static double max_curvature                 (Point_3 point,mesh_t& mesh);
            
            static double max_angle                     (Point_3 point,mesh_t& mesh);
            static double min_angle                     (Point_3 point,mesh_t& mesh);
            static double mean_normalized_angle         (Point_3 point,mesh_t& mesh);       //compareable to mean curvature
            static double multiplied_normalized_angle   (Point_3 point,mesh_t& mesh);       //compareable to gaussan curvature


            //needed static funktions for the statistic
            static Monge_form get_monge_form(Point_3 point,mesh_t& mesh);
            static Vector_3 get_surface_normal(mesh_t::Vertex& vertex);
            static mesh_t::Vertex& convert_Point_to_Vertex(Point_3 point,mesh_t& mesh);
            static std::vector<Point_3> closest_points(long amount,Point_3 point,mesh_t& mesh);
            static Vector_3 exproduct(Vector_3 a,Vector_3 b);
            static double norm(Vector_3 e);
            static bool equals(Facet_t& a,Facet_t& b);
            static bool equals(Point_3& a,Point_3& b);

            std::vector<viennagrid_quantity_field> get_quantity_fields(){return quantity_fields;}

            //test curvature igl
#if __enable_libigl
            void curvature_igl(mesh_t& mesh);
#endif //__enable_libigl

        private:
            void make_curvature();
            void statistic();

            long state=0;

            double                  old_area=0;
            double                  new_area=0;

            mesh_t                  old_mesh;
            Tree                    old_tree;

            mesh_t                  new_mesh;
            Tree                    new_tree;

            std::vector<std::vector<Point_3>>   closest_Points;
            std::vector<Point_3>                refrence_points;

            std::vector<std::vector<Facet_t>>   closest_Facets;
            std::vector<Facet_t>                refrence_Facets;

            std::vector<statistic_data_t>   data_Point;     //all quantities that are stored dor Points
            std::vector<statistic_data_t>   data_Edge;      //all quantities that are stored dor Edges
            std::vector<statistic_data_t>   data_Facet;     //all quantities that are stored dor Facets

            std::vector<viennagrid_quantity_field> quantity_fields;
            
            
    };
}

#endif