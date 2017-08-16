#include "cgal_statistic.hpp"
#include <CGAL/squared_distance_3.h> //for 3D functions
#include <CGAL/squared_distance_2.h>


// using namespace cgal_statistic
using namespace viennamesh;

//looks if the Facets have the same Points in the same order 
bool cgal_statistic::equals(Facet_t& a,Facet_t& b)
{
    if(a.halfedge()->vertex()->point()==b.halfedge()->vertex()->point())
        if(a.halfedge()->next()->vertex()->point()==b.halfedge()->next()->vertex()->point())
            if(a.halfedge()->next()->next()->vertex()->point()==b.halfedge()->next()->next()->vertex()->point())
                return true;
    return false;
}

//looks if Point a equals Point b
bool cgal_statistic::equals(Point_3& a,Point_3& b)
{
    std::cout << a << ":" << b << "\n";
    return (a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2])?true:false;
}

//exproduct of two vectors
cgal_statistic::Vector_3 cgal_statistic::exproduct(Vector_3 a,Vector_3 b)
{
    return Vector_3(a[1]*b[2]-a[2]*b[1],-a[0]*b[2]+a[2]*b[0],a[0]*b[1]-a[1]*b[0]);
}

//adds all 3 Points and divides them by 3 to get the center of all 3 points
cgal_statistic::Point_3 cgal_statistic::middel_Point(Point_3 a,Point_3 b,Point_3 c)
{
    Point_3 zero(0,0,0);
    return zero+((((a+(b-zero))+(c-zero))-zero)/3);
}

//gets the norm(length) of the vector
double cgal_statistic::norm(Vector_3 e)
{
    return std::sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
}

//return the quality of the facet
double cgal_statistic::quality(Point_3 a,Point_3 b,Point_3 c)
{
    
    double R=0;     //Umkreis
    double r=0;     //Inkreis

    double ab=norm(b-a);    //distanzen
    double bc=norm(c-b);
    double ca=norm(a-c);

    R=(ab*bc*ca)/(4.0*norm(exproduct(b-a,c-b))/2.0);    //https://de.wikipedia.org/wiki/Umkreis
    r=norm(exproduct(b-a,c-b))/(ab+bc+ca);              //https://de.wikipedia.org/wiki/Inkreis

    return 1/2.0*R/r;                                   //quality=1/2*umkreis/inkreis
}

//diffrent quality compare funktions
double cgal_statistic::quality_diffrence(double orig,double coarse)      // the bigger the better
{
    return orig-coarse;
}

double cgal_statistic::relative_quality(double orig,double coarse)       //the bigger the better
{
    return orig/coarse;
}

double cgal_statistic::relative_quality_diffrence(double orig,double coarse)//the bigger the better
{
    return (orig-coarse)/std::max<double>(orig,coarse);
}

//gets the squared distance of the Point to the mesh
double cgal_statistic::squared_distance(Point_3 point,Tree& tree)
{
    return tree.squared_distance(point);
}

// makes the mean of the "array" then subtracts the coarsed
double cgal_statistic::diffrence_of_mean(std::vector<double> orig,double coarse)
{
    double acc=0;
    if(orig.size()==0)return 0;
    //std::cout << std::accumulate(orig.begin(),orig.end(),acc)/orig.size() << "\n";
    return std::accumulate(orig.begin(),orig.end(),acc)/orig.size()-coarse;
}
//retuns the Monge form around a Point of a givn mesh
//From what I understand the Mongeform is a mathematical Plane that is fitted to go through the Points
//http://doc.cgal.org/latest/Jet_fitting_3/Jet_fitting_3_2Mesh_estimation_8cpp-example.html
//https://github.com/CGAL/cgal/tree/master/Jet_fitting_3/examples/Jet_fitting_3
cgal_statistic::Monge_form cgal_statistic::get_monge_form (Point_3 point,mesh_t& mesh)
{
    unsigned int d_fitting = 2;
    unsigned int d_monge = 2;
    unsigned int needed = 6;
    
    Vector_3  normal=get_surface_normal(convert_Point_to_Vertex(point,mesh));

    std::vector<Point_3> points=closest_points(needed,point,mesh);              //<---- These Points
    if(points.size()<needed)
    {
        std::cerr << "not enough points: " << points.size() <<"\n";
        // return 0;
    }

    CGAL::Monge_via_jet_fitting<Kernel> monge_fit;
    CGAL::Monge_via_jet_fitting<Kernel>::Monge_form monge_form 
            = monge_fit(points.begin(),points.end(),d_fitting, d_monge);
               
    monge_form.comply_wrt_given_normal(normal);

    return monge_form;
}

//mean curvature 
double cgal_statistic::mean_curvature(Point_3 point,mesh_t& mesh)
{
    Monge_form monge_form=get_monge_form (point,mesh);
    //monge_form.principal_curvatures(0 is maximum 1 is minimum)
    return (monge_form.principal_curvatures(0)+monge_form.principal_curvatures(1))/2;
}

//gaus curvature
double cgal_statistic::gaus_curvature(Point_3 point,mesh_t& mesh)
{
    Monge_form monge_form=get_monge_form (point,mesh);
    return monge_form.principal_curvatures(0)*monge_form.principal_curvatures(1);
}


//gathers the closest points by collecting the neighbors and there neighbors and so on if nessesary
std::vector<cgal_statistic::Point_3> cgal_statistic::closest_points(long amount,Point_3 point,mesh_t& mesh)
{
    mesh_t::Vertex vertex=convert_Point_to_Vertex(point,mesh);
    std::vector<Point_3> out,tmp;
    std::vector<mesh_t::Vertex> new_vertexes;
    //new_vertexes.reserve(40);
    long id=0,id_end;
   
    //std::cout << "test\n";
    out.push_back(point);
    new_vertexes.push_back(vertex);

    while(out.size()<amount)
    {

        //muss ein long "iterator"
        id_end=new_vertexes.size();
        while(id!=id_end)
        {
             mesh_t::Vertex::Halfedge_around_vertex_circulator at=new_vertexes.at(id).vertex_begin(),end=at;
            do
            {
               
                if(std::find(out.begin(),out.end(),at->opposite()->vertex()->point())==out.end())
                {
                    out.push_back(at->opposite()->vertex()->point());
                    //std::cout << "now adding the vertex\n";
                    new_vertexes.push_back(* (at->opposite()->vertex()));
                }
            
                ++at;
            }while(at!=end);
            ++id;
        }
    }

    //only return the amount of closest points

    // std::vector<double> sort,orig;
    // sort.reserve(out.size());
    // for(std::vector<Point_3>::iterator at=out.begin(),end=out.end();at!=end;++at)
    // {
    //     sort.push_back(norm(*at-point));
    // }
    // orig=sort;
    // std::sort(sort.begin(),sort.end());

    // tmp=out;
    // out.clear();
    // for(long i=0;i<amount;i++)
    // {
    //     out.push_back(tmp.at(std::distance(orig.begin(),std::find(orig.begin(),orig.end(),sort.at(i)))));
    // }
    return out;
}


//gets the vertex that has the same coordinates as the point
cgal_statistic::mesh_t::Vertex& cgal_statistic::convert_Point_to_Vertex(Point_3 point,mesh_t& mesh)
{
    for(mesh_t::Vertex_iterator at=mesh.vertices_begin(),end=mesh.vertices_end();at!=end;++at)
        if(at->point()==point)
            return *at;
    
    std::cerr << "nothing found maby try a diffrent mesh\n";
    static mesh_t::Vertex v(Point_3(0,0,0));
    return v;
}  

//return the facenormal at the vertex by looking at the neighbors
cgal_statistic::Vector_3 cgal_statistic::get_surface_normal(mesh_t::Vertex& vertex)
{
    Vector_3 out(0,0,0),normal;
    mesh_t::Vertex::Halfedge_around_vertex_circulator at=vertex.vertex_begin(),end=at;
    
    do
    {
        if(! at->is_border_edge())
        {
            Vector_3 a=at->opposite()->vertex()->point() - at->vertex()->point();
            Vector_3 b=at->next()->opposite()->vertex()->point() - at->next()->vertex()->point();
            out=out+exproduct(a,b);
        }
        ++at;
    }while(at!=end);

    return out/norm(out);
}       


//test igl curvature
//#include "../statistics/mesh_comparison.hpp"
//#include "../statistics/libigl_convert.hpp"

void convert_cgal_to_igl(cgal_statistic::mesh_t& mesh, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& Vertices, Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& Facets)
{
    typedef cgal_statistic::Point_3 Point_3;
    Vertices.resize(mesh.size_of_vertices(),3);
    Facets.resize(mesh.size_of_facets(),3);
    std::vector<cgal_statistic::Point_3> points_vector(mesh.size_of_vertices());    //needed for faster index search

    //vertixes insert
    typedef cgal_statistic::mesh_t::Point_iterator iterator_vertex;
    long index=0;
    for(iterator_vertex at=mesh.points_begin(),end=mesh.points_end();at!=end;++at,++index)
    {
        Vertices(index,0)= (*at)[0];
        Vertices(index,1)= (*at)[1];
        Vertices(index,2)= (*at)[2];
        points_vector.at(index)=(*at);
    }

    //Facets insert
    typedef cgal_statistic::mesh_t::Facet_iterator iterator_facet;
    index=0;
    for(iterator_facet at=mesh.facets_begin(),end=mesh.facets_end();at!=end;++at,++index)
    {
        //getting the Points of the facet
        Point_3 a=at->halfedge()->vertex()->point();
        Point_3 b=at->halfedge()->next()->vertex()->point();
        Point_3 c=at->halfedge()->next()->next()->vertex()->point();

        //getting the index of the points
        long index_a=std::distance(points_vector.begin(),std::find(points_vector.begin(),points_vector.end(),a));
        long index_b=std::distance(points_vector.begin(),std::find(points_vector.begin(),points_vector.end(),b));
        long index_c=std::distance(points_vector.begin(),std::find(points_vector.begin(),points_vector.end(),c));

        Facets(index,0)=index_a;
        Facets(index,1)=index_b;
        Facets(index,2)=index_c;
    }
}

#if __enable_libigl

#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>


void cgal_statistic::curvature_igl(mesh_t& mesh)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Vertices;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Facets;

    convert_cgal_to_igl(mesh, Vertices, Facets);

    Eigen::VectorXd K;
    // Compute integral of Gaussian curvature
    igl::gaussian_curvature(Vertices, Facets,K);
    
    std::cout << "K size:" << K.size() << "\n";
    

    //extremly ugly only for test purposes
    statistic_data_t data("igl_gaussan_curvvature",quantity_dimention_t::Vertex);
    data.data.resize(K.size());
    Eigen::VectorXd::Map(&data.data[0], K.size()) = K;
    data_Point.push_back(data);
    std::cout << "igl curvature size:" << data.data.size() << "\n";

    //principal curvature
    Eigen::VectorXd PV1,PV2,PD1,PD2;
    igl::principal_curvature(Vertices,Facets,PD1,PD2,PV1,PV2);
    std::vector<double> max(PV1.size()),min(PV2.size());

    Eigen::VectorXd::Map(&max[0], PV1.size()) = PV1;
    Eigen::VectorXd::Map(&min[0], PV2.size()) = PV2;

    std::cout << "PV1.size:" << PV1.size() << "\n";
    statistic_data_t mean("igl_mean_curvvature",quantity_dimention_t::Vertex);
    mean.data.resize(PV1.size());

    for(long index=0;index<PV1.size();++index)
    {
        mean.data.at(index)=(max.at(index)+min.at(index))/2;
    }
    data_Point.push_back(mean);
}  
#endif //__enable_libigl