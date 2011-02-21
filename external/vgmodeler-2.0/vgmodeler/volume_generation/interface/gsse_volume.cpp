/* =======================================================================

                     ---------------------------------------
                     ViennaMesh - The Vienna Meshing Library
                     ---------------------------------------
                            
   maintainers:    
               Josef Weinbub                      weinbub@iue.tuwien.ac.at
               Johann Cervenka                   cervenka@iue.tuwien.ac.at

   developers:
               Franz Stimpfl, Rene Heinzl, Philipp Schwaha, Joachim Schoeberl

   license:    
               LGPL, see file LICENSE in the ViennaMesh base directory

   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at

   ======================================================================= */

#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <stlgeom.hpp>
#include <meshing.hpp>


namespace vgmodeler {
#include "nglib.h"
}

using namespace vgmodeler;


namespace vgmodeler 
{

//Destination for messages, errors, ...
void Ng_PrintDest(const char * s)
{
#ifdef DEBUG
  std::cout << s << flush;
#endif
}



}



// constants and types:

namespace vgmodeler
{


// [RH] transition from nginterface
//
void Ng_SetRefinementFlag (Ng_Mesh *mesh, int ei, int flag)
{
  if (((Mesh*)mesh)->GetDimension() == 3)
    {
      ((Mesh*)mesh)->VolumeElement(ei).SetRefinementFlag (flag != 0);
      ((Mesh*)mesh)->VolumeElement(ei).SetStrongRefinementFlag (flag >= 10);
    }
  else
    {
      ((Mesh*)mesh)->SurfaceElement(ei).SetRefinementFlag (flag != 0);
      ((Mesh*)mesh)->SurfaceElement(ei).SetStrongRefinementFlag (flag >= 10);
    }
}

void Ng_set_surface_element (Ng_Mesh *mesh, int ei, int flag)
{
  if (((Mesh*)mesh)->GetDimension() == 3)
    {
      ((Mesh*)mesh)->SurfaceElement(ei).SetRefinementFlag (flag != 0);
      ((Mesh*)mesh)->SurfaceElement(ei).SetStrongRefinementFlag (flag >= 10);
    }
}


void Ng_Refine (Ng_Mesh *mesh, NG_REFINEMENT_TYPE reftype)
{
  BisectionOptions biopt;
  biopt.usemarkedelements = 1;
  biopt.refine_p = 0;
  biopt.refine_hp = 0;
//   if (reftype == NG_REFINE_P)
//     biopt.refine_p = 1;
//   if (reftype == NG_REFINE_HP)
//     biopt.refine_hp = 1;


  Refinement * ref;
  ref =new Refinement();
  // RH tetrahedra bisection
  //
  ref -> Bisect (*((Mesh*)mesh), biopt);

  // RH .. another refinement method
  //
  //ref -> Refine (*((Mesh*)mesh));


  ((Mesh*)mesh) -> UpdateTopology();
  delete ref;
}


void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h)
{
  ((Mesh*)mesh)->RestrictLocalH(Point3d (p[0], p[1], p[2]), h);
  ((Mesh*)mesh)->SetGlobalH( h);

  }




void Ng_GetCellOnVertex(Ng_Mesh *mesh, int vertex_number, int **element_container, int& number_of_elements)
{
  ARRAY<int> element_container_local;    
  ((Mesh*)mesh)->get_elements_on_vertex(vertex_number, element_container_local);

  int *element_c_temp;
  number_of_elements = element_container_local.Size();
  element_c_temp = (int*)realloc (*element_container, number_of_elements* sizeof(int));
  if (element_c_temp !=0)
    *element_container = element_c_temp;
  else
    {
      element_container = NULL;
      return;
    }
  for (int i = 0; i < element_container_local.Size(); ++i)
    {
      //std::cout << "out: "<< element_container_local[i] << std::endl;
      (*element_container)[i] = element_container_local[i];
    }
}




Ng_Mesh * Ng_NewMesh ()
{
  Mesh * mesh = new Mesh;  
  mesh->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
  return (Ng_Mesh*)(void*)mesh;
}


} // namespace end....


