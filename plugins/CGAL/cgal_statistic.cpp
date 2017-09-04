#include "cgal_statistic.hpp"

#include "cgal_mesh.hpp"

#include <cmath>

using namespace viennamesh;


//plugin run()
bool cgal_statistic::run(algorithm_handle&)
{

    //Get mesh data (REQUIRED)
    std::cout << "original mesh\n";
    data_handle<cgal::polyhedron_surface_mesh> original_mesh    = get_required_input<cgal::polyhedron_surface_mesh>("original mesh");
    std::cout << "coarse mesh\n";
    data_handle<cgal::polyhedron_surface_mesh> coarse_mesh      = get_required_input<cgal::polyhedron_surface_mesh>("coarse mesh"); 

    //optional stuff
    data_handle<bool>   use_default_quantities=get_input<bool>("default quantities");
    //data_handle<std::vector<statistic_data_t>> additional_quantities=get_input<std::vector<statistic_data_t>>("quantities");
    old_mesh=original_mesh();
    new_mesh=coarse_mesh();


    if(use_default_quantities.valid())
    {
        if(use_default_quantities())
        {
            default_quantities();
        }
    }

/*
    if(additional_quantities.valid())
    {
        std::vector<statistic_data_t> tmp=additional_quantities()
        for(std::vector<statistic_data_t>::iterator at=tmp.begin(),end=tmp.end();at!=end;++at)
        {
            add_quantity(*at);
        }
    }
    */


    //run statistic
    operator() ();
    operator() ();

    //set output
    set_output<cgal::polyhedron_surface_mesh>("mesh", new_mesh);
    make_quantity_fields();

    std::vector<viennagrid_quantity_field> quantity_fields=get_quantity_fields();
    std::cout << "quantity fields size: " << quantity_fields.size() << "\n"; 
    if(quantity_fields.size()!=0)
    {
        data_handle<viennagrid_quantity_field> quantity_field_handle=viennamesh::plugin_algorithm::make_data(to_cpp(quantity_fields[0]));            
        if(quantity_fields.size()>1)
            for(std::vector<viennagrid_quantity_field>::iterator at=++ quantity_fields.begin(),end=quantity_fields.end();at!=end;++at)
                quantity_field_handle.push_back(to_cpp(*at));             
        set_output("quantities",quantity_field_handle);
    }
    else 
        std::cout << "no quantities\n";
}

//default quantities
void cgal_statistic::default_quantities()
{
    add_quantity(statistic_data_t(
        "quality",
        quantity_dimention_t::Facet,
        (quantity_Funktions_t)&quality));

    add_quantity(statistic_data_t(
        "quality diffrence",
        quantity_dimention_t::Facet,
        (quantity_Funktions_t)  &quality,
        quantity_compare_t::closest_Facet_value,
        (compare_Funktions_t)   &quality_diffrence,
        (point_Funktions_t)     &middel_Point
    ));

    add_quantity(statistic_data_t(
        "quality diffrence of clostest Facets",
        quantity_dimention_t::Facet,
        (quantity_Funktions_t)  &quality,
        quantity_compare_t::closest_Facets_value,
        (compare_Funktions_t)   &diffrence_of_mean
    ));

    add_quantity(statistic_data_t(
        "realtive quality diffrence",
        quantity_dimention_t::Facet,
        (quantity_Funktions_t)  &quality,
        quantity_compare_t::closest_Facet_value,
        (compare_Funktions_t)   &relative_quality_diffrence,
        (point_Funktions_t)     &middel_Point
    ));

    add_quantity(statistic_data_t(
        "realtive quality",
        quantity_dimention_t::Facet,
        (quantity_Funktions_t)  &quality,
        quantity_compare_t::closest_Facet_value,
        (compare_Funktions_t)   &relative_quality,
        (point_Funktions_t)     &middel_Point
    ));

    add_quantity(statistic_data_t(
        "squared distance",
        quantity_dimention_t::Facet,
        quantity_Funktions_t(),     //hier hätte ich gerne noch was schönes
        quantity_compare_t::tree_orig,
        (compare_Funktions_t)   &squared_distance,
        (point_Funktions_t)     &middel_Point
    ));

    add_quantity(statistic_data_t(
        "mean curvature",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &mean_curvature
    ));

    add_quantity(statistic_data_t(
        "gaus curvature",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &gaus_curvature
    ));

    add_quantity(statistic_data_t(
        "max curvature",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &max_curvature
    ));

    add_quantity(statistic_data_t(
        "min curvature",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &min_curvature
    ));

    add_quantity(statistic_data_t(
        "mean curvature compare to closest points",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &mean_curvature,
        quantity_compare_t::closest_Points_value,
        (compare_Funktions_t)   &diffrence_of_mean
    ));

    add_quantity(statistic_data_t(
        "mean curvature compare to closest point",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &mean_curvature,
        quantity_compare_t::closest_Point_value,
        (compare_Funktions_t)   &quality_diffrence
    ));

    add_quantity(statistic_data_t(
        "max angle",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &max_angle,
        quantity_compare_t::none
    ));


    add_quantity(statistic_data_t(
        "min angle",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &min_angle,
        quantity_compare_t::none
    ));

    add_quantity(statistic_data_t(
        "mean angle",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &mean_normalized_angle,
        quantity_compare_t::none
    ));

    add_quantity(statistic_data_t(
        "gaus angle(gedanklich vergleichbar mit gaus curvature)",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &multiplied_normalized_angle,
        quantity_compare_t::none
    ));
    add_quantity(statistic_data_t(
        "mean angle compared to closest Point",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &mean_normalized_angle,
        quantity_compare_t::closest_Point_value,
        (compare_Funktions_t)   &quality_diffrence
    ));

    add_quantity(statistic_data_t(
        "gaus angle(gedanklich vergleichbar mit gaus curvature) compared to the closest Point",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &multiplied_normalized_angle,
        quantity_compare_t::closest_Point_value,
        (compare_Funktions_t)   &quality_diffrence
    ));

    add_quantity(statistic_data_t(
        "mean angle compared to the mean of the closest points",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &mean_normalized_angle,
        quantity_compare_t::closest_Points_value,
        (compare_Funktions_t)   &diffrence_of_mean
    ));

    add_quantity(statistic_data_t(
        "gaus angle(gedanklich vergleichbar mit gaus curvature) compared to the mean of the closest points",
        quantity_dimention_t::Vertex,
        (quantity_Funktions_t)  &multiplied_normalized_angle,
        quantity_compare_t::closest_Points_value,
        (compare_Funktions_t)   &diffrence_of_mean
    ));
    
}


//add quantity 
//can be called externaly to add a quantity 
void cgal_statistic::add_quantity(statistic_data_t data)
{
    switch(data.dimention)
    {
        case quantity_dimention_t::Vertex :
            std::cout << "adds a Vertex quantity with the name: " << data.name << "\n";
            data_Point.push_back(data);
            break;
        case quantity_dimention_t::Edge :
            std::cout << "adds a Edge   quantity with the name: " << data.name << "\n";
            data_Edge.push_back(data);
            break;
        case quantity_dimention_t::Facet :
            std::cout << "adds a Facet  quantity with the name: " << data.name << "\n";
            data_Facet.push_back(data);
            break;
    }
}


//basically the run function 
//needs to be called twice 
//  once after the initialization or seting the old(original) mesh 
//  once after the new(corsed) mesh was set by get_new_mesh()=coarsed_mesh
//this generates the needed structurs and then runs statisttic()
void cgal_statistic::operator()()
{
    if(state==0)
    {
        old_tree.insert(faces(old_mesh).begin(),faces(old_mesh).end(),old_mesh);
        old_tree.build();
        old_tree.accelerate_distance_queries();

        std::cout << "Mesh size: " << old_mesh.size_of_facets() << "\n";
        std::cout << "Tree size: " << old_tree.size() << "\n";
        ++state;
    }
    else if(state==1)
    {
        if (new_mesh.empty() ) throw "Needs a corased mesh";
        new_tree.insert(faces(new_mesh).begin(),faces(new_mesh).end(),new_mesh);
        new_tree.build();
        new_tree.accelerate_distance_queries();
        statistic();
        std::cout <<"org Area: " << old_area << "corsed Area: " << new_area << "\n";
        ++state;
        curvature_igl(new_mesh);
        // make_curvature();
    }
}


//generates the wanted data (very costly operation)
//linkes every point/facets of the coarse mesh with the closest points/facets of the original mesh
//at the point of documentation 13.8.17 the following combinations exist
//      dimention   |   compare             |   quantity function param | compare function param
// -----------------|-----------------------|---------------------------|------------------------
//    Point/Vertex  |         none          | Point_3,mesh_t new_mesh   |
//                  |  closest_Point_value  | Point_3,mesh_t            |   doible orig,double coarsed
//                  |  closest_Points_value | Point_3,mesh_t            |   std::vector<double> orig,double coarsed
//                  |       default         | not jet implemented error |
//                  |                       |                           |                                       
//   Edge/Line      | not jet implemented error                         |
//                  |                       |                           |
//   Facet/Area     |        none           | Point_3,Point_3,Point_3   |
//                  |                       | all points of the Facet   |
//                  |  closest_Facet_value  | Point_3,Point_3,Point_3   | double orig,double coarsed
//                  |  closest_Facets_value | Point_3,Point_3,Point_3   | std::vector<double> orig,double coarsed
//                  |       tree_orig       | Point_3,Point_3,Point_3   | double coarsed,Tree orig_tree
//                  |       tree_coarse     | Point_3,Point_3,Point_3   | double coarsed,Tree coarsed_tree
//                  |       closest_Facet   | Facet_t,std::vector<Facet_t> closest_Facets
//                  |       default         | not jet implementedd error


void cgal_statistic::statistic()
{
    //closest points
    
    refrence_points.reserve(new_mesh.size_of_vertices());
    closest_Points.resize(new_mesh.size_of_vertices());

    for(mesh_t::Point_iterator at=new_mesh.points_begin(),end=new_mesh.points_end();at!=end;++at)
    {
        refrence_points.push_back(*at);
    }
    for(mesh_t::Point_iterator at=old_mesh.points_begin(),end=old_mesh.points_end();at!=end;++at)
    {
        Point_3 point=new_tree.closest_point(*at); //näheste punkt auf der flache
        long id=0,use;
        double dist=-1;
        for(std::vector<Point_3>::iterator comp=refrence_points.begin(),comp_end=refrence_points.end();comp!=comp_end;++comp,++id)
        {
            if(norm(point-*comp)<dist || dist==-1)
            {
                use=id;
                dist=norm(point-*comp);
            }
        }
        
        
        closest_Points.at(use).push_back(*at);
    }
    
    //closest facets
    refrence_Facets.reserve(new_mesh.size_of_facets());
    closest_Facets.resize(new_mesh.size_of_facets());
    for(mesh_t::Facet_iterator at=new_mesh.facets_begin(),end=new_mesh.facets_end();at!=end;++at)
    {
        refrence_Facets.push_back(*at);
        new_area+=norm(exproduct(at->halfedge()->vertex()->point()-at->halfedge()->next()->vertex()->point(),at->halfedge()->next()->vertex()->point()-at->halfedge()->next()->next()->vertex()->point()))/2;
    } 

    for(mesh_t::Facet_iterator at=old_mesh.facets_begin(),end=old_mesh.facets_end();at!=end;++at)
    {
        Facet_t facet=* new_tree.closest_point_and_primitive(middel_Point(
           at->halfedge()->vertex()->point(),
           at->halfedge()->next()->vertex()->point(),
           at->halfedge()->next()->next()->vertex()->point() 
        )).second;
        old_area+=norm(exproduct(at->halfedge()->vertex()->point()-at->halfedge()->next()->vertex()->point(),at->halfedge()->next()->vertex()->point()-at->halfedge()->next()->next()->vertex()->point()))/2;
    
        // closest_Facets.at(std::distance(refrence_Facets.begin(),std::find(refrence_Facets.begin(),refrence_Facets.end(),facet))).push_back(*at);
        long id=0;
        for(std::vector<Facet_t>::iterator comp=refrence_Facets.begin(),end=refrence_Facets.end();comp!=end;++comp,++id)
            if(equals(*comp,facet))
                break;
        if(id==refrence_Facets.size())
        {
            std::cerr << "compare didn't work\n";
            continue;
        }
        closest_Facets.at(id).push_back(facet);
        // if(closest_Facets.at(id).size()>1)
        //     std::cout << "closest facets " << id << " size " << closest_Facets.at(id).size() << "\n";
    }   


    
    if(!data_Point.empty())
    {
        typedef mesh_t::Point_iterator iterator;
        for(std::vector<statistic_data_t>::iterator at=data_Point.begin(),end=data_Point.end();at!=end;++at)
        {
            at->data.resize(new_mesh.size_of_vertices());
        }
        long id=0;
        for(iterator at=new_mesh.points_begin(),end=new_mesh.points_end();at!=end;++at,++id)
        {
            long id2=0;
            for(std::vector<Point_3>::iterator cmp=refrence_points.begin(),cmp_end=refrence_points.end();cmp!=cmp_end&&*at!=*cmp;++cmp,++id2);
            if(refrence_points.at(id2)!=*at) std::cerr << "wrong point found\n";
            std::vector<Point_3>& current_closest_Points=closest_Points.at(id2);

            for(std::vector<statistic_data_t>::iterator stor_at=data_Point.begin(),stor_end=data_Point.end();stor_at!=stor_end;++stor_at)
            {
                std::vector<double> closest_Values;
                switch(stor_at->compare)
                {
                    case quantity_compare_t::none :
                        stor_at->data.at(id)=stor_at->quantity_Funktion(*at,new_mesh);
                        break;

                    case quantity_compare_t::closest_Point_value :
                        if(current_closest_Points.size()>0)
                        {
                            Point_3 closest_Point;
                            double dist=-1;
                            for(std::vector<Point_3>::iterator cmp=current_closest_Points.begin(),cmp_end=current_closest_Points.end();cmp!=cmp_end;++cmp)
                            {
                               if(norm(*at-*cmp)<dist || dist==-1)
                                {
                                    closest_Point=*cmp;
                                    dist=norm(*at-*cmp);
                                }
                            }
                            stor_at->data.at(id)= stor_at->compare_Funktion(stor_at->quantity_Funktion(closest_Point,old_mesh),stor_at->quantity_Funktion(*at,new_mesh));
                        }
                    break;

                    case quantity_compare_t::closest_Points_value :
                        closest_Values.clear();
                        closest_Values.reserve(current_closest_Points.size());
                        for(std::vector<Point_3>::iterator cmp=current_closest_Points.begin(),cmp_end=current_closest_Points.end();cmp!=cmp_end;++cmp)
                        {
                            double tmp=stor_at->quantity_Funktion(*cmp,old_mesh);
                           // std::cout << tmp << "\n";
                           //if(isnan(tmp))std::cerr << "found nan " << *cmp << "\n";
                            closest_Values.push_back(tmp);
                        }
                        if(current_closest_Points.size()==0) std::cerr << " no closest points\n";
                        stor_at->data.at(id)= stor_at->compare_Funktion(closest_Values,stor_at->quantity_Funktion(*at,new_mesh));
                        break;
                    default:
                        std::cerr << "Not jet implemented\n";
                        throw "not jet implemented";
                }
            }
        }
    }
    else std::cout << "no Point quantities\n";
    if(!data_Edge.empty())
    {
        std::cerr << "Edge quantities not jet implementet\n";
        throw "not jet implemented\n";
    } 
    else std::cout << "no Edge quantities\n";
    if(!data_Facet.empty())
    {
        typedef mesh_t::Facet_iterator iterator;
        for(std::vector<statistic_data_t>::iterator at=data_Facet.begin(),end=data_Facet.end();at!=end;++at)
        {
            at->data.resize(new_mesh.size_of_facets());
        }
        long id=0;
        for(iterator at=new_mesh.facets_begin(),end=new_mesh.facets_end();at!=end;++at,++id)
        {
            Point_3 a=at->facet_begin()->vertex()->point();
            Point_3 b=at->facet_begin()->next()->vertex()->point();
            Point_3 c=at->facet_begin()->next()->next()->vertex()->point();

            long id2=0;
            for(std::vector<Facet_t>::iterator it=refrence_Facets.begin(),tmp=refrence_Facets.end();!equals(*it,*at)&&it!=tmp;++it,++id2);
            
            if(!equals(refrence_Facets.at(id2),*at)) std::cerr << "wrong refernce\n";

            std::vector<Facet_t>& current_closest_Facets=closest_Facets.at(id2);
                
            for(std::vector<statistic_data_t>::iterator stor_at=data_Facet.begin(),stor_end=data_Facet.end();stor_at!=stor_end;++stor_at)
            {
                Face face;
                Point_3 face_a;
                Point_3 face_b;
                Point_3 face_c;
                std::vector<double> closest_Facets_values;
                switch(stor_at->compare)
                {
                    case quantity_compare_t::none :
                        stor_at->data.at(id)=stor_at->quantity_Funktion(a,b,c);
                        break;
                    case quantity_compare_t::closest_Facet_value :
                        face=old_tree.closest_point_and_primitive(stor_at->point_Funktion(a,b,c)).second;
                        face_a=face->halfedge()->vertex()->point();
                        face_b=face->halfedge()->next()->vertex()->point();
                        face_c=face->halfedge()->next()->next()->vertex()->point();
                        stor_at->data.at(id)=stor_at->compare_Funktion(
                                    stor_at->quantity_Funktion(face_a,face_b,face_c),
                                    stor_at->quantity_Funktion(a,b,c)
                                    );
                        break;
                    case quantity_compare_t::tree_orig :
                        stor_at->data.at(id)=stor_at->compare_Funktion(
                            stor_at->point_Funktion(a,b,c),
                            old_tree
                            );
                        break;
                    case quantity_compare_t::tree_coarse :
                        stor_at->data.at(id)=stor_at->compare_Funktion(
                            stor_at->point_Funktion(a,b,c),
                            new_tree
                            );
                        break;
                    case quantity_compare_t::closest_Facets_value :
                        closest_Facets_values.clear();
                        closest_Facets_values.reserve(current_closest_Facets.size());
                        for(std::vector<Facet_t>::iterator qual=current_closest_Facets.begin(),qual_end=current_closest_Facets.end();qual!=qual_end;++qual)
                            {
                                face_a=qual->halfedge()->vertex()->point();
                                face_b=qual->halfedge()->next()->vertex()->point();
                                face_c=qual->halfedge()->next()->next()->vertex()->point();
                                closest_Facets_values.push_back(stor_at->quantity_Funktion(face_a,face_b,face_c));
                            }
                            if(closest_Facets_values.size()==0)std::cerr << "no closest Facet\n";
                        stor_at->data.at(id)=stor_at->compare_Funktion(closest_Facets_values,stor_at->quantity_Funktion(a,b,c));
                        break;
                    case quantity_compare_t::closest_Facet :
                         stor_at->data.at(id)=stor_at->quantity_Funktion(*at,current_closest_Facets);
                        break;
                    case quantity_compare_t::closest_Point_value :
                        throw "this doesnt make sence\"for hole Facet get the closest Point\"";
                        break;
                    default :
                        std::cerr << "Not jet implemented\n";
                        throw "this is not implementet";
                }
            }
        }
    }
    else std::cout << "no Facet quantities\n";

}

//Initializes the quantity fields correctly
void quantity_field_choose(viennagrid_quantity_field quantity_field,cgal_statistic::statistic_data_t data)
{
    if(data.dimention==cgal_statistic::quantity_dimention_t::Vertex)
    {
        viennagrid_quantity_field_init(quantity_field,
                            VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                            VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                            1,                                        // one float per element
                            VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);
    }
    else if(data.dimention==cgal_statistic::quantity_dimention_t::Edge)
    {
        viennagrid_quantity_field_init(quantity_field,
                            VIENNAGRID_ELEMENT_TYPE_EDGE ,                                        // topological dimension of the elements
                            VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                            1,                                        // one float per element
                            VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);
    }
    else if(data.dimention==cgal_statistic::quantity_dimention_t::Facet)
    {
        viennagrid_quantity_field_init(quantity_field,
                            VIENNAGRID_ELEMENT_TYPE_TRIANGLE ,                                        // topological dimension of the elements
                            VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                            1,                                        // one float per element
                            VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);
    }
}

//this dosen't work it has to be made manually in the plugin that uses this class
data_handle<viennagrid_quantity_field> cgal_statistic::make_data_handle()
{
    if(quantity_fields.size()==0) throw "no quantities";
    data_handle<viennagrid_quantity_field> quantity_field_handle=viennamesh::plugin_algorithm::make_data(to_cpp(quantity_fields[0]));            
    if(quantity_fields.size()>1)
        for(std::vector<viennagrid_quantity_field>::iterator at=++ quantity_fields.begin(),end=quantity_fields.end();at!=end;++at)
        {
            quantity_field_handle.push_back(to_cpp(*at));
        }              
            
    return quantity_field_handle;
    
}
//creates  one big quantity vector containing all quantities
void cgal_statistic::make_quantity_fields()
{

    if(data_Point.size()==0 && data_Edge.size()==0 && data_Facet.size()==0) throw "no quantities";
    
    //make one big vector
    std::vector<statistic_data_t> tmp=data_Point;
    tmp.reserve(tmp.size() + distance(data_Edge.begin(),data_Edge.end()));
    tmp.insert(tmp.end(),data_Edge.begin(),data_Edge.end());

    tmp.reserve(tmp.size() + distance(data_Facet.begin(),data_Facet.end()));
    tmp.insert(tmp.end(),data_Facet.begin(),data_Facet.end());
    quantity_fields.reserve(tmp.size());
    for(std::vector<statistic_data_t>::iterator at=tmp.begin(),end=tmp.end();at!=end;++at)
    {
        long id=0;
        viennagrid_quantity_field quantity_field=0;
        viennagrid_quantity_field_create(&quantity_field);
        quantity_field_choose(quantity_field,*at);
        for(statistic_data_t::data_t::iterator data_at=at->data.begin(),data_end=at->data.end();data_at!=data_end;++data_at,++id)
        {
            viennagrid_quantity_field_value_set(quantity_field,id,&(*data_at));
        }
        viennagrid_quantity_field_name_set(quantity_field,at->name.c_str());
        quantity_fields.push_back(quantity_field);
    }
}
