#ifndef VIENNAMESH_DOMAIN_CGAL_PLC_HPP
#define VIENNAMESH_DOMAIN_CGAL_PLC_HPP

#include "viennagrid/domain/config.hpp"
#include "viennagrid/domain/element_creation.hpp"
#include "viennagrid/algorithm/geometry.hpp"

#include "viennamesh/base/convert.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>




namespace viennamesh
{
    struct cgal_plc_2d_element
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
        typedef CGAL::Triangulation_vertex_base_2<Kernel>               VertexBase;
        typedef CGAL::Delaunay_mesh_face_base_2<Kernel>                      FaceBase;
        typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>        Triangulation_structure;
        typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_structure>  CDT;
        typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
        typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
        
        typedef CDT::Vertex_handle Vertex_handle;
        typedef CDT::Point Point;
      
        
        CDT cdt;
        std::list<Point> cgal_list_of_holes;
    };
    
    struct cgal_plc_2d_domain
    {
        typedef std::vector<cgal_plc_2d_element> cell_container;
        cell_container cells;
    };
    
    
    
    
    
    template<>
    struct convert_impl<viennagrid::config::plc_2d_domain, cgal_plc_2d_domain>
    {
        typedef viennagrid::config::plc_2d_domain vgrid_domain_type;
        typedef viennagrid::config::plc_2d_cell vgrid_element_type;
        
        typedef cgal_plc_2d_domain cgal_domain_type;
        typedef cgal_plc_2d_element cgal_element_type;
        
        
        typedef vgrid_domain_type input_domain_type;
        typedef vgrid_element_type input_element_type;
        
        typedef cgal_domain_type output_domain_type;
        typedef cgal_element_type output_element_type;
        
        static bool convert( vgrid_domain_type const & vgrid_domain, cgal_domain_type & cgal_domain )
        {
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::plc_tag>::type plc_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::plc_tag>::type plc_hook_type;
            
            typedef viennagrid::result_of::const_element_range<vgrid_domain_type, viennagrid::plc_tag>::type plc_range_type;
            typedef viennagrid::result_of::iterator<plc_range_type>::type plc_range_iterator;
            
            typedef viennagrid::result_of::point_type<vgrid_domain_type>::type point_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::vertex_tag>::type vertex_type;
            typedef viennagrid::result_of::const_element_hook<vgrid_domain_type, viennagrid::vertex_tag>::type vertex_const_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::line_tag>::type line_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::line_tag>::type line_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::polygon_tag>::type polygon_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::polygon_tag>::type polygon_hook_type;
            
            
            typedef viennagrid::result_of::const_element_range<vgrid_element_type, viennagrid::vertex_tag>::type vertex_range_type;
            typedef viennagrid::result_of::hook_iterator<vertex_range_type>::type vertex_range_hook_iterator;
            
            typedef viennagrid::result_of::const_element_range<vgrid_element_type, viennagrid::line_tag>::type line_range_type;
            typedef viennagrid::result_of::iterator<line_range_type>::type line_range_iterator;
            
            
            
            
            plc_range_type plcs = viennagrid::elements<viennagrid::plc_tag>( vgrid_domain );
            cgal_domain.cells.resize(plcs.size());
            
            cgal_domain_type::cell_container::iterator jt = cgal_domain.cells.begin();
            
            for (plc_range_iterator it = plcs.begin(); it != plcs.end(); ++it)
            {
                output_element_type & cgal_element = *jt;
                vgrid_element_type const & vgrid_element = *it;
                
                std::map<vertex_const_hook_type, cgal_element_type::Vertex_handle> vertex_handle_map;
                
                vertex_range_type vertices = viennagrid::elements<viennagrid::vertex_tag>(vgrid_element);
                for (vertex_range_hook_iterator it = vertices.hook_begin(); it != vertices.hook_end(); ++it)
                {
                    vertex_const_hook_type const & vtx_hook = *it;
                    vertex_type const & vtx = viennagrid::dereference_hook(vgrid_domain, *it);
                    point_type const & vgrid_point = viennagrid::point(vgrid_domain, vtx );
                    
                    cgal_element_type::Vertex_handle handle = cgal_element.cdt.insert( cgal_element_type::Point(vgrid_point[0], vgrid_point[1]) );
                    
                    vertex_handle_map[vtx_hook] = handle;
                }

                
                line_range_type lines = viennagrid::elements<viennagrid::line_tag>(vgrid_element);
                for (line_range_iterator it = lines.begin(); it != lines.end(); ++it)
                {
                    line_type const & line = *it;
                    
                    vertex_const_hook_type vgrid_v0 = viennagrid::elements<viennagrid::vertex_tag>(line).hook_at(0);
                    vertex_const_hook_type vgrid_v1 = viennagrid::elements<viennagrid::vertex_tag>(line).hook_at(1);
                    
                    cgal_element_type::Vertex_handle cgal_v0 = vertex_handle_map[vgrid_v0];
                    cgal_element_type::Vertex_handle cgal_v1 = vertex_handle_map[vgrid_v1];
                    
                    cgal_element.cdt.insert_constraint(cgal_v0, cgal_v1);
                }
                
                const std::vector<point_type> & vgrid_list_of_holes = viennagrid::hole_points<vgrid_domain_type>(vgrid_element);
                
                for (std::vector<point_type>::const_iterator it = vgrid_list_of_holes.begin(); it != vgrid_list_of_holes.end(); ++it)
                    cgal_element.cgal_list_of_holes.push_back( cgal_element_type::Point( (*it)[0], (*it)[1] ) );
            }
            
            return true;
        }  
    };
    
    
    
    
    
    template<>
    struct convert_impl<cgal_plc_2d_domain, viennagrid::config::triangular_2d_domain>
    {              
        typedef cgal_plc_2d_domain cgal_domain_type;
        typedef cgal_plc_2d_element cgal_element_type;
        
        typedef viennagrid::config::triangular_2d_domain vgrid_domain_type;
        typedef viennagrid::config::triangular_2d_cell vgrid_element_type;
        
        
        typedef cgal_domain_type input_domain_type;
        typedef cgal_element_type input_element_type;
        
        typedef vgrid_domain_type output_domain_type;
        typedef vgrid_element_type output_element_type;

        
        static bool convert( cgal_domain_type const & cgal_domain, vgrid_domain_type & vgrid_domain )
        {
            typedef viennagrid::result_of::point_type<vgrid_domain_type>::type triangle_point_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::vertex_tag>::type triangle_vertex_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::vertex_tag>::type triangle_vertex_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::line_tag>::type triangle_line_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::line_tag>::type triangle_line_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::triangle_tag>::type triangle_triangle_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::triangle_tag>::type triangle_triangle_hook_type;
            
            
            std::map<cgal_element_type::Point, triangle_vertex_hook_type> points;
            
            for (cgal_domain_type::cell_container::const_iterator it = cgal_domain.cells.begin(); it != cgal_domain.cells.end(); ++it)
            {
                cgal_element_type const & cgal_element = *it;
                
                for(cgal_element_type::CDT::Finite_faces_iterator fit = cgal_element.cdt.finite_faces_begin(); fit != cgal_element.cdt.finite_faces_end(); ++fit) 
                {
                    if(fit->is_in_domain())
                    {
                        typedef cgal_element_type::CDT::Triangle Triangle;
                        Triangle tri = cgal_element.cdt.triangle(fit);
                        
                        viennagrid::storage::static_array<triangle_vertex_hook_type, 3> vgrid_vtx;
                        
                        for (int i = 0; i < 3; ++i)
                        {
                            std::map<cgal_element_type::Point, triangle_vertex_hook_type>::iterator pit = points.find( tri[i] );
                            if (pit == points.end())
                            {
                                vgrid_vtx[i] = viennagrid::create_element<triangle_vertex_type>( vgrid_domain, triangle_point_type(tri[i].x(), tri[i].y()) );
                                points[ tri[i] ] = vgrid_vtx[i];
                            }
                            else
                                vgrid_vtx[i] = pit->second;
                        }
                        
                        viennagrid::create_element<triangle_triangle_type>( vgrid_domain, vgrid_vtx.begin(), vgrid_vtx.end() );
                    }
                }
                
            }
            
            return true;
        }
    };
    
    

    
    
    struct cgal_plc_3d_element
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
        typedef CGAL::Triangulation_vertex_base_2<Kernel>               VertexBase;
        typedef CGAL::Delaunay_mesh_face_base_2<Kernel>                      FaceBase;
        typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>        Triangulation_structure;
        typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Triangulation_structure>  CDT;
        typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
        typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
        
        typedef CDT::Vertex_handle Vertex_handle;
        typedef CDT::Point Point;
      
        viennagrid::config::point_type_3d center;
        viennagrid::storage::static_array<viennagrid::config::point_type_3d, 2> projection_matrix;
        
        CDT cdt;
        
        std::list<Point> cgal_list_of_holes;
    };
    
    struct cgal_plc_3d_domain
    {
        typedef std::deque<cgal_plc_3d_element> cell_container;
        cell_container cells;
    };
    
    
    template<>
    struct convert_impl<viennagrid::config::plc_3d_domain, cgal_plc_3d_domain>
    {
        typedef viennagrid::config::plc_3d_domain vgrid_domain_type;
        typedef viennagrid::config::plc_3d_cell vgrid_element_type;
        
        typedef cgal_plc_3d_domain cgal_domain_type;
        typedef cgal_plc_3d_element cgal_element_type;
        
        
        typedef vgrid_domain_type input_domain_type;
        typedef vgrid_element_type input_element_type;
        
        typedef cgal_domain_type output_domain_type;
        typedef cgal_element_type output_element_type;
        
        static bool convert( vgrid_domain_type const & vgrid_domain, cgal_domain_type & cgal_domain )
        {
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::plc_tag>::type plc_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::plc_tag>::type plc_hook_type;
            
            typedef viennagrid::result_of::const_element_range<vgrid_domain_type, viennagrid::plc_tag>::type plc_range_type;
            typedef viennagrid::result_of::iterator<plc_range_type>::type plc_range_iterator;
            
            typedef viennagrid::result_of::point_type<vgrid_domain_type>::type point_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::vertex_tag>::type vertex_type;
            typedef viennagrid::result_of::const_element_hook<vgrid_domain_type, viennagrid::vertex_tag>::type vertex_const_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::line_tag>::type line_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::line_tag>::type line_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::polygon_tag>::type polygon_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::polygon_tag>::type polygon_hook_type;
            
            
            typedef viennagrid::result_of::const_element_range<vgrid_element_type, viennagrid::vertex_tag>::type vertex_range_type;
            typedef viennagrid::result_of::hook_iterator<vertex_range_type>::type vertex_range_hook_iterator;
            
            typedef viennagrid::result_of::const_element_range<vgrid_element_type, viennagrid::line_tag>::type line_range_type;
            typedef viennagrid::result_of::iterator<line_range_type>::type line_range_iterator;
            
            
            
            
            typedef viennagrid::config::point_type_2d point_type_2d;
            
            plc_range_type plcs = viennagrid::elements<viennagrid::plc_tag>( vgrid_domain );
            cgal_domain.cells.resize(plcs.size());
            
            cgal_domain_type::cell_container::iterator jt = cgal_domain.cells.begin();
            
            for (plc_range_iterator it = plcs.begin(); it != plcs.end(); ++it, ++jt)
            {
                output_element_type & cgal_element = *jt;
                vgrid_element_type const & vgrid_element = *it;
                
                std::vector<point_type> plc_points_3d;
                std::vector<point_type_2d> plc_points_2d;
                
                std::map<vertex_const_hook_type, std::size_t> vertex_to_index_map;
                
                vertex_range_type vertices = viennagrid::elements<viennagrid::vertex_tag>(vgrid_element);
                for (vertex_range_hook_iterator it = vertices.hook_begin(); it != vertices.hook_end(); ++it)
                {
                    vertex_const_hook_type const & vtx_hook = *it;
                    vertex_type const & vtx = viennagrid::dereference_hook(vgrid_domain, *it);
                    point_type const & vgrid_point = viennagrid::point(vgrid_domain, vtx);
                    
                    vertex_to_index_map[vtx_hook] = plc_points_3d.size();
                    plc_points_3d.push_back(vgrid_point);
                }

                
                viennagrid::geometry::projection_matrix( plc_points_3d.begin(), plc_points_3d.end(), 1e-6, cgal_element.center, cgal_element.projection_matrix.begin() );
                
//                 std::cout << cgal_element.projection_matrix[0] << " " << cgal_element.projection_matrix[1] << " " << cgal_element.center << std::endl;
                
                plc_points_2d.resize( plc_points_3d.size() );
                viennagrid::geometry::project( plc_points_3d.begin(), plc_points_3d.end(), plc_points_2d.begin(), cgal_element.center, cgal_element.projection_matrix.begin(), cgal_element.projection_matrix.end() );

                
                std::map<vertex_const_hook_type, cgal_element_type::Vertex_handle> vertex_handle_map;
                
                for (vertex_range_hook_iterator it = vertices.hook_begin(); it != vertices.hook_end(); ++it)
                {
                    vertex_const_hook_type const & vtx_hook = *it;
                    vertex_type const & vtx = viennagrid::dereference_hook(vgrid_domain, *it);
                    point_type_2d const & vgrid_point = plc_points_2d[ vertex_to_index_map[vtx_hook] ];
//                     std::cout << "3D Point: " << plc_points_3d[ vertex_to_index_map[vtx_hook] ] << " projected onto " << vgrid_point << std::endl;
                    
                    cgal_element_type::Vertex_handle handle = cgal_element.cdt.insert( cgal_element_type::Point(vgrid_point[0], vgrid_point[1]) );
                    
                    vertex_handle_map[vtx_hook] = handle;
                }

                
                line_range_type lines = viennagrid::elements<viennagrid::line_tag>(vgrid_element);
                for (line_range_iterator it = lines.begin(); it != lines.end(); ++it)
                {
                    line_type const & line = *it;
                    
                    vertex_const_hook_type vgrid_v0 = viennagrid::elements<viennagrid::vertex_tag>(line).hook_at(0);
                    vertex_const_hook_type vgrid_v1 = viennagrid::elements<viennagrid::vertex_tag>(line).hook_at(1);
                    
                    cgal_element_type::Vertex_handle cgal_v0 = vertex_handle_map[vgrid_v0];
                    cgal_element_type::Vertex_handle cgal_v1 = vertex_handle_map[vgrid_v1];
                    
                    cgal_element.cdt.insert_constraint(cgal_v0, cgal_v1);
                }
                
                const std::vector<point_type> & vgrid_list_of_holes = viennagrid::hole_points<vgrid_domain_type>(vgrid_element);
                std::vector<point_type_2d> vgrid_list_of_holes_2d(vgrid_list_of_holes.size());
                
                viennagrid::geometry::project( vgrid_list_of_holes.begin(), vgrid_list_of_holes.end(), vgrid_list_of_holes_2d.begin(), cgal_element.center, cgal_element.projection_matrix.begin(), cgal_element.projection_matrix.end() );
                
                for (std::vector<point_type_2d>::iterator it = vgrid_list_of_holes_2d.begin(); it != vgrid_list_of_holes_2d.end(); ++it)
                    cgal_element.cgal_list_of_holes.push_back( cgal_element_type::Point( (*it)[0], (*it)[1] ) );
            } 
            
            return true;
        }
    };
    
    
    
    

    
    
    
    template<>
    struct convert_impl<cgal_plc_3d_domain, viennagrid::config::triangular_3d_domain>
    {              
        typedef cgal_plc_3d_domain cgal_domain_type;
        typedef cgal_plc_3d_element cgal_element_type;
        
        typedef viennagrid::config::triangular_3d_domain vgrid_domain_type;
        typedef viennagrid::config::triangular_3d_cell vgrid_element_type;
        
        
        typedef cgal_domain_type input_domain_type;
        typedef cgal_element_type input_element_type;
        
        typedef vgrid_domain_type output_domain_type;
        typedef vgrid_element_type output_element_type;

        
        static bool convert( cgal_domain_type const & cgal_domain, vgrid_domain_type & vgrid_domain )
        {
            typedef viennagrid::result_of::point_type<vgrid_domain_type>::type triangle_point_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::vertex_tag>::type triangle_vertex_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::vertex_tag>::type triangle_vertex_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::line_tag>::type triangle_line_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::line_tag>::type triangle_line_hook_type;
            
            typedef viennagrid::result_of::element<vgrid_domain_type, viennagrid::triangle_tag>::type triangle_triangle_type;
            typedef viennagrid::result_of::element_hook<vgrid_domain_type, viennagrid::triangle_tag>::type triangle_triangle_hook_type;
            
            
            typedef std::deque< std::pair<triangle_point_type, triangle_vertex_hook_type> > points_container;
            points_container points;
            
            
//             vgrid_domain_type tmp_domain;
            
            
            for (cgal_domain_type::cell_container::const_iterator it = cgal_domain.cells.begin(); it != cgal_domain.cells.end(); ++it)
            {
//                 std::cout << "Converting back a PLC" << std::endl;
                cgal_element_type const & cgal_element = *it;
                
                for(cgal_element_type::CDT::Finite_faces_iterator fit = cgal_element.cdt.finite_faces_begin(); fit != cgal_element.cdt.finite_faces_end(); ++fit) 
                {
//                     std::cout << " Triangle found!" << std::endl;
                    if(fit->is_in_domain())
                    {
                        typedef cgal_element_type::CDT::Triangle Triangle;
                        Triangle tri = cgal_element.cdt.triangle(fit);
                        
                        viennagrid::storage::static_array<triangle_vertex_hook_type, 3> vgrid_vtx;
                        
                        for (int i = 0; i < 3; ++i)
                        {
                            triangle_point_type point_3d = cgal_element.projection_matrix[0] * tri[i].x() + cgal_element.projection_matrix[1] * tri[i].y() + cgal_element.center;
                            
                            points_container::iterator pit = points.begin();
                            for (; pit != points.end(); ++pit)
                            {
                                if ( viennagrid::norm_2( point_3d - pit->first ) < 1e-6 )
                                {
                                    vgrid_vtx[i] = pit->second;
                                    break;
                                }
                            }
                            
                            if (pit == points.end())
                            {
                                vgrid_vtx[i] = viennagrid::create_element<triangle_vertex_type>( vgrid_domain, point_3d );
                                points.push_back( std::make_pair( point_3d, vgrid_vtx[i] ) );
                            }
                        }
                        
                        viennagrid::create_element<triangle_triangle_type>( vgrid_domain, vgrid_vtx.begin(), vgrid_vtx.end() );
                    }
                }
            }
            
            
            
            return true;
        }
    };
    
    
    
}

#endif