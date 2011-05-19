#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>

#include "stlgeom.hpp"

#include "gsse/domain.hpp"

namespace netgen
{


STLTopology :: STLTopology()
  : trias(), topedges(), points(), ht_topedges(NULL), 
    neighbourtrigs(), trigsperpoint()
{
  ;
}

STLTopology :: ~STLTopology()
{
  ;
}




STLGeometry *  STLTopology :: LoadBinary (istream & ist)
{
  STLGeometry * geom = new STLGeometry();
  ARRAY<STLReadTriangle> readtrigs;

  PrintMessage(1,"Read STL binary file");
  
  if (sizeof(int) != 4 || sizeof(float) != 4) 
    {
      PrintWarning("for stl-binary compatibility only use 32 bit compilation!!!");
    }

  //specific settings for stl-binary format
  const int namelen = 80; //length of name of header in file
  const int nospaces = 2; //number of spaces after a triangle

  //read header: name
  char buf[namelen+1];
  FIOReadStringE(ist,buf,namelen);
  PrintMessage(5,"header = ",buf);

  //Read Number of facets
  int nofacets;
  FIOReadInt(ist,nofacets);
  PrintMessage(5,"NO facets = ",nofacets);

  Point<3> pts[3];
  Vec<3> normal;

  int cntface, j;
  int vertex = 0;
  float f;
  char spaces[nospaces+1];

  for (cntface = 0; cntface < nofacets; cntface++)
    {
      if (cntface % 10000 == 9999) { PrintDot(); } 

      FIOReadFloat(ist,f); normal(0) = f;
      FIOReadFloat(ist,f); normal(1) = f;
      FIOReadFloat(ist,f); normal(2) = f;
      
      for (j = 0; j < 3; j++)
	{
	  FIOReadFloat(ist,f); pts[j](0) = f;
	  FIOReadFloat(ist,f); pts[j](1) = f;
	  FIOReadFloat(ist,f); pts[j](2) = f;	  
	} 

      readtrigs.Append (STLReadTriangle (pts, normal));
      FIOReadString(ist,spaces,nospaces);
    }	    
  

  geom->InitSTLGeometry(readtrigs);

  return geom;
}


void STLTopology :: SaveBinary (const char* filename, const char* aname)
{
  ofstream ost(filename);
  PrintFnStart("Write STL binary file '",filename,"'");

  if (sizeof(int) != 4 || sizeof(float) != 4) 
    {PrintWarning("for stl-binary compatibility only use 32 bit compilation!!!");}

  //specific settings for stl-binary format
  const int namelen = 80; //length of name of header in file
  const int nospaces = 2; //number of spaces after a triangle

  //write header: aname
  int i, j;
  char buf[namelen+1];
  int strend = 0;
  for(i = 0; i <= namelen; i++) 
    {
      if (aname[i] == 0) {strend = 1;}
      if (!strend) {buf[i] = aname[i];}
      else {buf[i] = 0;}
    }

  FIOWriteString(ost,buf,namelen);
  PrintMessage(5,"header = ",buf);

  //RWrite Number of facets
  int nofacets = GetNT();
  FIOWriteInt(ost,nofacets);
  PrintMessage(5,"NO facets = ", nofacets);

  float f;
  char spaces[nospaces+1];
  for (i = 0; i < nospaces; i++) {spaces[i] = ' ';}
  spaces[nospaces] = 0;

  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & t = GetTriangle(i);

      const Vec<3> & n = t.Normal();
      f = n(0); FIOWriteFloat(ost,f);
      f = n(1); FIOWriteFloat(ost,f);
      f = n(2); FIOWriteFloat(ost,f);

      for (j = 1; j <= 3; j++)
	{
	  const Point3d p = GetPoint(t.PNum(j));
	  
	  f = p.X(); FIOWriteFloat(ost,f);
	  f = p.Y(); FIOWriteFloat(ost,f);
	  f = p.Z(); FIOWriteFloat(ost,f);
	}
      FIOWriteString(ost,spaces,nospaces);
    }
  PrintMessage(5,"done");
}


void STLTopology :: SaveSTLE (const char* filename)
{
  ofstream outf (filename);
  int i, j;
  
  outf << GetNT() << endl;
  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & t = GetTriangle(i);
      for (j = 1; j <= 3; j++)
	{
	  const Point3d p = GetPoint(t.PNum(j));
	  outf << p.X() << " " << p.Y() << " " << p.Z() << endl;
	}
    }


  int ned = 0;
  for (i = 1; i <= GetNTE(); i++)
    {
      if (GetTopEdge (i).GetStatus() == ED_CONFIRMED)
	ned++;
    }
  
  outf << ned << endl;

  for (i = 1; i <= GetNTE(); i++)
    {
      const STLTopEdge & edge = GetTopEdge (i);
      if (edge.GetStatus() == ED_CONFIRMED)
	for (j = 1; j <= 2; j++)
	  {
	    const Point3d p = GetPoint(edge.PNum(j));
	    outf << p.X() << " " << p.Y() << " " << p.Z() << endl;
	  }
    }      
}



STLGeometry *  STLTopology :: LoadNaomi (istream & ist)
{
  int i;
  STLGeometry * geom = new STLGeometry();
  ARRAY<STLReadTriangle> readtrigs;

  PrintFnStart("read NAOMI file format");
  
  char buf[100];
  Vec<3> normal;

  int cntface = 0;
  int cntvertex = 0;
  double px, py, pz;
    

  int noface, novertex;
  ARRAY<Point<3> > readpoints;

  ist >> buf;
  if (strcmp (buf, "NODES") == 0)
    {
      ist >> novertex;
      PrintMessage(5,"nuber of vertices = ", novertex);
      for (i = 0; i < novertex; i++)
	{
	  ist >> px;
	  ist >> py;
	  ist >> pz;
	  readpoints.Append(Point<3> (px,py,pz));
	}
    }
  else
    {
      PrintFileError("no node information");
    }


  ist >> buf;
  if (strcmp (buf, "2D_EDGES") == 0)
    {
      ist >> noface;
      PrintMessage(5,"number of faces=",noface);
      int dummy, p1, p2, p3;
      Point<3> pts[3];

      for (i = 0; i < noface; i++)
	{
	  ist >> dummy; //2
	  ist >> dummy; //1
	  ist >> p1;
	  ist >> p2;
	  ist >> p3;
	  ist >> dummy; //0

	  pts[0] = readpoints.Get(p1);
	  pts[1] = readpoints.Get(p2);
	  pts[2] = readpoints.Get(p3);
	  
	  normal = Cross (pts[1]-pts[0], pts[2]-pts[0]) . Normalize();

	  readtrigs.Append (STLReadTriangle (pts, normal));

	}
      PrintMessage(5,"read ", readtrigs.Size(), " triangles");
    }
  else
    {
      PrintMessage(5,"read='",buf,"'\n");
      PrintFileError("ERROR: no Triangle information");
    }

  geom->InitSTLGeometry(readtrigs);

  return geom;
}

void STLTopology :: Save (const char* filename)
{ 
  PrintFnStart("Write stl-file '",filename, "'");

  ofstream fout(filename);
  fout << "solid\n";

  char buf1[50];
  char buf2[50];
  char buf3[50];

  int i, j;
  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & t = GetTriangle(i);

      fout << "facet normal ";
      const Vec3d& n = GetTriangle(i).Normal();

      sprintf(buf1,"%1.9g",n.X());
      sprintf(buf2,"%1.9g",n.Y());
      sprintf(buf3,"%1.9g",n.Z());

      fout << buf1 << " " << buf2 << " " << buf3 << "\n";
      fout << "outer loop\n";

      for (j = 1; j <= 3; j++)
	{
	  const Point3d p = GetPoint(t.PNum(j));
	  
	  sprintf(buf1,"%1.9g",p.X());
	  sprintf(buf2,"%1.9g",p.Y());
	  sprintf(buf3,"%1.9g",p.Z());

	  fout << "vertex " << buf1 << " " << buf2 << " " << buf3 << "\n";
	}

      fout << "endloop\n";
      fout << "endfacet\n"; 
    }
  fout << "endsolid\n";

  
  // write also NETGEN surface mesh:
  ofstream fout2("geom.surf");
  fout2 << "surfacemesh" << endl;
  fout2 << GetNP() << endl;
  for (i = 1; i <= GetNP(); i++)
    {
      for (j = 0; j < 3; j++)
	{
	  fout2.width(8);
	  fout2 << GetPoint(i)(j);
	}

      fout2 << endl;
    }

  fout2 << GetNT() << endl;
  for (i = 1; i <= GetNT(); i++)
    {
      const STLTriangle & t = GetTriangle(i);  
      for (j = 1; j <= 3; j++)
	{
	  fout2.width(8);
	  fout2 << t.PNum(j);
	}
      fout2 << endl;
    }
}


STLGeometry *  STLTopology ::Load (istream & ist)
{
  int i;
  STLGeometry * geom = new STLGeometry();

  ARRAY<STLReadTriangle> readtrigs;

  char buf[100];
  Point<3> pts[3];
  Vec<3> normal;

  int cntface = 0;
  int vertex = 0;
  bool badnormals = 0;

  while (ist.good())
    {
      ist >> buf;

      int n = strlen (buf);
      for (i = 0; i < n; i++)
	buf[i] = tolower (buf[i]);

      if (strcmp (buf, "facet") == 0)
	{
	  cntface++;
	}

      if (strcmp (buf, "normal") == 0)
	{
	  ist >> normal(0)
	      >> normal(1)
	      >> normal(2);
	  normal.Normalize();
	}
      
      if (strcmp (buf, "vertex") == 0)
	{
	  ist >> pts[vertex](0)
	      >> pts[vertex](1)
	      >> pts[vertex](2);

	  vertex++;

	  if (vertex == 3)
	    {
	      if (normal.Length() <= 1e-5)

		{
		  normal = Cross (pts[1]-pts[0], pts[2]-pts[0]);
		  normal.Normalize();
		}

	      else

		{
		  Vec<3> hnormal;
		  hnormal = Cross (pts[1]-pts[0], pts[2]-pts[0]);
		  hnormal.Normalize();

		  if (normal * hnormal < 0.5)
		    {
		      badnormals = 1;
		    }
		}

	      vertex = 0;

	      if ( (Dist2 (pts[0], pts[1]) > 1e-16) &&
		   (Dist2 (pts[0], pts[2]) > 1e-16) &&
		   (Dist2 (pts[1], pts[2]) > 1e-16) )
		
		readtrigs.Append (STLReadTriangle (pts, normal));
	    }
	}
    }
  
  if (badnormals) 
    {
      PrintWarning("File has normal vectors which differ extremly from geometry->correct with stldoctor!!!");
    }

  geom->InitSTLGeometry(readtrigs);
  return geom;
}













void STLTopology :: InitSTLGeometry(const ARRAY<STLReadTriangle> & readtrigs)
{
  int i, j, k;
  
  // const double geometry_tol_fact = 1E6; 
  // distances lower than max_box_size/tol are ignored

  trias.SetSize(0);
  points.SetSize(0);

  PrintMessage(3,"number of triangles = ", readtrigs.Size());

  if (!readtrigs.Size())
    return;
  

  boundingbox.Set (readtrigs[0][0]);
  for (i = 0; i < readtrigs.Size(); i++)
    for (k = 0; k < 3; k++)
      boundingbox.Add (readtrigs[i][k]);
  
  PrintMessage(5,"boundingbox: ", Point3d(boundingbox.PMin()), " - ", 
	       Point3d(boundingbox.PMax()));
  
  int netgen = 0;
  if(netgen)
    {
  // [INFO] .. NETGEN CODE STARTS HERE
  //
  Box<3> bb = boundingbox;
  bb.Increase (1);

  pointtree = new Point3dTree (bb.PMin(), bb.PMax());

  ARRAY<int> pintersect;

  pointtol = boundingbox.Diam() * stldoctor.geom_tol_fact;

  // [FS] .. here the point tolerance decides which two non identical points are treated as one
  //  
  std::cout << "point tolerance = " << pointtol << std::endl;


  for(i = 0; i < readtrigs.Size(); i++)
    {
      const STLReadTriangle & t = readtrigs[i];
      STLTriangle st;
      Vec<3> n = t.Normal();

      std::cout << ".. normal: " << t.Normal() << std::endl;

      st.SetNormal (t.Normal());

      for (k = 0; k < 3; k++)
	{
	  Point<3> p = t[k];

	  Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
	  Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);

	  pointtree->GetIntersecting (pmin, pmax, pintersect);

	  std::cout << ".. after intersection test: " << std::endl;

	  if (pintersect.Size() > 1)
	    PrintError("too many close points");

	  int foundpos = -1;

	  // [FS] .. test: add the interface point a second time	 
	  //

 	  if (pintersect.Size())
 	    foundpos = pintersect[0];

	  std::cout << ".. foundpos: " << foundpos << " :: for point: " << p <<  std::endl;

 	  if (foundpos == -1)
	    {
	      foundpos = AddPoint(p);
	      pointtree->Insert (p, foundpos);
	    }
	  st[k] = foundpos;
	}

      if ( (st[0] == st[1]) ||
	   (st[0] == st[2]) || 
	   (st[1] == st[2]) )
	{
	  PrintError("STL Triangle degenerated");
	}
      else
	{
	  std::cout << ".. triangle added: " << st[0] << " / " << st[1] << " / " << st[2] << std::endl;

	  AddTriangle(st);
	}

    } 
    }
  else
    {

  // [INFO] .. GSSE STARTS HERE
  //
  typedef gsse::detail_topology::unstructured<2>                                unstructured_topology_2t; 
  typedef gsse::get_domain<unstructured_topology_2t, double, double,3>::type    domain_32t;

  typedef gsse::domain_traits<domain_32t>::cell_iterator		  cell_iterator;
  typedef gsse::domain_traits<domain_32t>::vertex_handle              vertex_handle;
  typedef gsse::domain_traits<domain_32t>::cell_2_vertex_mapping    cell_2_vertex_mapping;
  typedef gsse::domain_traits<domain_32t>::segment_iterator	  segment_iterator;
  typedef gsse::domain_traits<domain_32t>::point_t             point_t;
  typedef gsse::domain_traits<domain_32t>::vertex_on_cell_iterator  vertex_on_cell_iterator;
  domain_32t domain;
//   domain.read_file("multiple_surfaces_23.gau32");
  domain.read_file("levelset_test.gau32");


  ARRAY<int> pintersect;
  Box<3> bb = boundingbox;
  bb.Increase (1);
  
  pointtree = new Point3dTree (bb.PMin(), bb.PMax());
  
  typedef std::map<point_t, vertex_handle>::iterator point_map_iterator;
  boost::array<vertex_handle, 3>                       element_container;
  long element_counter = 1;  
  i=0;
  

  segment_iterator segit = domain.segment_begin();
  for( ; segit != domain.segment_end(); ++segit)
    {      
      cell_iterator cit = (*segit).cell_begin();


      std::map<point_t, vertex_handle> point_map;


      for( ; cit != (*segit).cell_end(); ++cit, i++)
	{
	  STLTriangle st;

	  const STLReadTriangle & t = readtrigs[i];
	  st.SetNormal (t.Normal());
	  std::cout << "..normal: " << t.Normal() << std::endl;

	  vertex_on_cell_iterator vocit(*cit);

	  for (int element_index = 0; element_index < 3; element_index++, ++vocit)
	    {
	      point_t point = domain.get_point(*vocit);

	      point_map_iterator pmit = point_map.find(point);
  
	      std::cout << "..checking point: " << point << std::endl;

	      if(pmit == point_map.end())
		{
		  std::cout << ".. adding point: " << point << std::endl;

		  Point<3> p(point[0],point[1],point[2]);		  
		  int foundpos = AddPoint(p);
		  pointtree->Insert (p, foundpos);

//  		  point_map[point] = element_counter;
//  		  element_container[element_index] = element_counter;

 		  point_map[point] = foundpos;
 		  element_container[element_index] = foundpos;
  	  	  
		  element_counter++;
		}
	      else
		{
		   element_container[element_index] = (*pmit).second;
		}
  
// 	      Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
// 	      Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);
  
// 	      pointtree->GetIntersecting (pmin, pmax, pintersect);
  
// 	      std::cout << ".. after intersection test: " << std::endl;
  
// 	      if (pintersect.Size() > 1)
// 		PrintError("too many close points");
  
// 	      int foundpos = -1;
  
// 	      // [FS] .. test: add the interface point a second time	 
// 	      //
  
// 	      if (pintersect.Size())
// 		foundpos = pintersect[0];
  
// 	      std::cout << ".. foundpos: " << foundpos << " :: for point: " << p <<  std::endl;
  
// 	      if (foundpos == -1)
// 		{
// 		  foundpos = AddPoint(p);
// 		  pointtree->Insert (p, foundpos);
// 		}
// 	      st[k] = foundpos;
	    }

	  st[0] = element_container[0];
	  st[1] = element_container[1];
	  st[2] = element_container[2];

	  std::cout << "..st[0]: " << st[0] << std::endl;
	  std::cout << "..st[1]: " << st[1] << std::endl;
	  std::cout << "..st[2]: " << st[2] << std::endl;

	  if ( (st[0] == st[1]) ||
	       (st[0] == st[2]) || 
	       (st[1] == st[2]) )
	    {
	      PrintError("STL Triangle degenerated");
	    }
	  else
	    {
	      AddTriangle(st);
	    }
	}

    }

    } // if netgen


  FindNeighbourTrigs();

  std::cout << "[FS] .. end of InitSTLGeometry: point size: " << GetNP() << std::endl;

}




int STLTopology :: GetPointNum (const Point<3> & p)
{
  Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
  Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);
  
  ARRAY<int> pintersect;

  pointtree->GetIntersecting (pmin, pmax, pintersect);
  if (pintersect.Size() == 1)
    return pintersect[0];
  else 
    return 0;
}



void STLTopology :: FindNeighbourTrigs()
{
  //  if (topedges.Size()) return;

  std::cout << "[FS] .. ENTER Find Neighbour Triangles" << std::endl;

  int i, j, k, l;

  // build up topology tables

  int np = GetNP();
  int nt = GetNT();
  
  std::cout << " .. np: " << np << std::endl;
  std::cout << " .. nt: " << nt << std::endl;


  INDEX_2_HASHTABLE<int> * oldedges = ht_topedges;
  ht_topedges = new INDEX_2_HASHTABLE<int> (GetNP()+1);
  topedges.SetSize(0);
  
  for (i = 1; i <= nt; i++)
    {
      STLTriangle & trig = GetTriangle(i);


      for (j = 1; j <= 3; j++)
	{
	  int pi1 = trig.PNumMod (j+1);
	  int pi2 = trig.PNumMod (j+2);
	  
	  INDEX_2 i2(pi1, pi2);
	  i2.Sort();

	  int enr;
	  int othertn;

	  if (ht_topedges->Used(i2))
	    {
	      enr = ht_topedges->Get(i2);
	      topedges.Elem(enr).TrigNum(2) = i;

	      othertn = topedges.Get(enr).TrigNum(1);
	      STLTriangle & othertrig = GetTriangle(othertn);

	      trig.NBTrigNum(j) = othertn;
	      trig.EdgeNum(j) = enr;
	      for (k = 1; k <= 3; k++)
		if (othertrig.EdgeNum(k) == enr)
		  othertrig.NBTrigNum(k) = i;
	    }
	  else
	    {
	      enr = topedges.Append (STLTopEdge (pi1, pi2, i, 0));
	      ht_topedges->Set (i2, enr);
	      trig.EdgeNum(j) = enr;
	    }
	}
    }

  
  PrintMessage(5,"topology built, checking");

  topology_ok = 1;
  int ne = GetNTE();

  for (i = 1; i <= nt; i++)
    GetTriangle(i).flags.toperror = 0;

  for (i = 1; i <= nt; i++)
    for (j = 1; j <= 3; j++)
      {
	const STLTopEdge & edge = GetTopEdge (GetTriangle(i).EdgeNum(j));
	if (edge.TrigNum(1) != i && edge.TrigNum(2) != i)
	  {
	    topology_ok = 0;
	    GetTriangle(i).flags.toperror = 1;
	  }
      }

  for (i = 1; i <= ne; i++)
    {
      const STLTopEdge & edge = GetTopEdge (i);
      if (!edge.TrigNum(2))
	{
	  topology_ok = 0;
	  GetTriangle(edge.TrigNum(1)).flags.toperror = 1;
	}
    }
 
  if (topology_ok)
    {
      orientation_ok = 1;
      for (i = 1; i <= nt; i++)
	{
	  const STLTriangle & t = GetTriangle (i);
	  for (j = 1; j <= 3; j++)
	    {
	      const STLTriangle & nbt = GetTriangle (t.NBTrigNum(j));
	      if (!t.IsNeighbourFrom (nbt))
		orientation_ok = 0;
	    }
	}
    }
  else
    orientation_ok = 0;
  


  status = STL_GOOD;
  statustext = "";
  if (!topology_ok || !orientation_ok)
    {
      status = STL_ERROR;
      if (!topology_ok)
	statustext = "Topology not ok";
      else
	statustext = "Orientation not ok";
    }


  PrintMessage(3,"topology_ok = ",topology_ok);
  PrintMessage(3,"orientation_ok = ",orientation_ok);
  PrintMessage(3,"topology found");

  // generate point -> trig table

  trigsperpoint.SetSize(GetNP());
  for (i = 1; i <= GetNT(); i++)
    for (j = 1; j <= 3; j++)
      trigsperpoint.Add1(GetTriangle(i).PNum(j),i);


  //check trigs per point:
  /*
  for (i = 1; i <= GetNP(); i++)
    {
      if (trigsperpoint.EntrySize(i) < 3)
	{
	  (*testout) << "ERROR: Point " << i << " has " << trigsperpoint.EntrySize(i) << " triangles!!!" << endl;
	}
    }
  */
  topedgesperpoint.SetSize (GetNP());
  for (i = 1; i <= ne; i++)
    for (j = 1; j <= 2; j++)
      topedgesperpoint.Add1 (GetTopEdge (i).PNum(j), i);

  PrintMessage(5,"point -> trig table generated");



  // transfer edge data:
  // .. to be done
  delete oldedges;



  for (STLTrigIndex ti = 0; ti < GetNT(); ti++)
    {
      STLTriangle & trig = trias[ti];
      for (k = 0; k < 3; k++)
	{
	  STLPointIndex pi = trig[k] - STLBASE;
	  STLPointIndex pi2 = trig[(k+1)%3] - STLBASE;
	  STLPointIndex pi3 = trig[(k+2)%3] - STLBASE;
	  
	  // vector along edge
	  Vec<3> ve = points[pi2] - points[pi];
	  ve.Normalize();

	  // vector along third point
	  Vec<3> vt = points[pi3] - points[pi];
	  vt -= (vt * ve) * ve;
	  vt.Normalize();

	  Vec<3> vn = trig.GeomNormal (points);
	  vn.Normalize();

	  double phimin = 10, phimax = -1; // out of (0, 2 pi)

	  for (j = 0; j < trigsperpoint[pi].Size(); j++)
	    {
	      STLTrigIndex ti2 = trigsperpoint[pi][j] - STLBASE;
	      const STLTriangle & trig2 = trias[ti2];

	      if (ti == ti2) continue;
	      
	      bool hasboth = 0;
	      for (l = 0; l < 3; l++)
		if (trig2[l] - STLBASE == pi2)
		  {
		    hasboth = 1;
		    break;
		  }
	      if (!hasboth) continue;

	      STLPointIndex pi4;
	      for (l = 0; l < 3; l++)
		if (trig2[l] - STLBASE != pi && trig2[l] - STLBASE != pi2)
		  pi4 = trig2[l] - STLBASE;

	      Vec<3> vt2 = points[pi4] - points[pi];
	      
	      double phi = atan2 (vt2 * vn, vt2 * vt);
	      if (phi < 0) phi += 2 * M_PI;
	      
	      if (phi < phimin)
		{
		  phimin = phi;
		  trig.NBTrig (0, (k+2)%3) = ti2 + STLBASE;
		}
	      if (phi > phimax)
		{
		  phimax = phi;
		  trig.NBTrig (1, (k+2)%3) = ti2 + STLBASE;
		}
	    }
	}
    }




  if (status == STL_GOOD)
    {
      // for compatibility:
      neighbourtrigs.SetSize(GetNT());
      for (i = 1; i <= GetNT(); i++)
	for (k = 1; k <= 3; k++)
	  AddNeighbourTrig (i, GetTriangle(i).NBTrigNum(k));
    }
  else
    {
      // assemble neighbourtrigs (should be done only for illegal topology):
      
      neighbourtrigs.SetSize(GetNT());

      int tr, found;
      int wrongneighbourfound = 0;
      for (i = 1; i <= GetNT(); i++)
	{
	  SetThreadPercent((double)i/(double)GetNT()*100.);
	  if (multithread.terminate)
	    {
	      PopStatus();
	      return;
	    }
	  
	  for (k = 1; k <= 3; k++)
	    {
	      for (j = 1; j <= trigsperpoint.EntrySize(GetTriangle(i).PNum(k)); j++)
		{
		  tr = trigsperpoint.Get(GetTriangle(i).PNum(k),j);
		  if (i != tr && (GetTriangle(i).IsNeighbourFrom(GetTriangle(tr))
				  || GetTriangle(i).IsWrongNeighbourFrom(GetTriangle(tr))))
		    {
		      if (GetTriangle(i).IsWrongNeighbourFrom(GetTriangle(tr)))
			{
			  /*(*testout) << "ERROR: triangle " << i << " has a wrong neighbour triangle!!!" << endl;*/
			  wrongneighbourfound ++;
			}
		      
		      found = 0;
		      for (int ii = 1; ii <= NONeighbourTrigs(i); ii++) 
			{if (NeighbourTrig(i,ii) == tr) {found = 1;break;};}
		      if (! found) {AddNeighbourTrig(i,tr);}
		    }
		}
	    }
	  if (NONeighbourTrigs(i) != 3) 
	    {
	      PrintError("TRIG ",i," has ",NONeighbourTrigs(i)," neighbours!!!!");
	      for (int kk=1; kk <= NONeighbourTrigs(i); kk++)
		{
		  PrintMessage(5,"neighbour-trig",kk," = ",NeighbourTrig(i,kk));
		}
	    };
	}
      if (wrongneighbourfound)
	{
	  PrintError("++++++++++++++++++++\n");
	  PrintError(wrongneighbourfound, " wrong oriented neighbourtriangles found!");
	  PrintError("try to correct it (with stldoctor)!");
	  PrintError("++++++++++++++++++++\n");
	  
	  status = STL_ERROR;
	  statustext = "STL Mesh not consistent";

	  multithread.terminate = 1;
#ifdef STAT_STREAM
	  (*statout) << "non-conform stl geometry \\hline" << endl;
#endif
	}
    }

  TopologyChanged();

  PopStatus();

  std::cout << "[FS] .. EXIT Find Neighbour Triangles" << std::endl;
}







void STLTopology :: GetTrianglesInBox (/* 
					  const Point<3> & pmin,
					  const Point<3> & pmax,
				       */
				       const Box<3> & box,
				       ARRAY<int> & trias) const
{
  if (searchtree)

    searchtree -> GetIntersecting (box.PMin(), box.PMax(), trias);
  
  else
    {    
      int i;
      Box<3> box1 = box;
      box1.Increase (1e-4);

      trias.SetSize(0);
   
      int nt = GetNT();
      for (i = 1; i <= nt; i++)
	{
	  if (box1.Intersect (GetTriangle(i).box))
	    {
	      trias.Append (i);
	    }
	}    
    }
}



void STLTopology :: AddTriangle(const STLTriangle& t)
{
  trias.Append(t);
  
  const Point<3> & p1 = GetPoint (t.PNum(1));
  const Point<3> & p2 = GetPoint (t.PNum(2));
  const Point<3> & p3 = GetPoint (t.PNum(3));

  Box<3> box;
  box.Set (p1);
  box.Add (p2);
  box.Add (p3);
  /*
  //  Point<3> pmin(p1), pmax(p1);
  pmin.SetToMin (p2);
  pmin.SetToMin (p3);
  pmax.SetToMax (p2);
  pmax.SetToMax (p3);
  */

  trias.Last().box = box; 
  trias.Last().center = Center (p1, p2, p3);
  double r1 = Dist (p1, trias.Last().center);
  double r2 = Dist (p2, trias.Last().center);
  double r3 = Dist (p3, trias.Last().center);
  trias.Last().rad = max2 (max2 (r1, r2), r3);

  if (geomsearchtreeon)
    {searchtree->Insert (box.PMin(), box.PMax(), trias.Size());}
}




int STLTopology :: GetLeftTrig(int p1, int p2) const
{
  int i;
  for (i = 1; i <= trigsperpoint.EntrySize(p1); i++)
    {
      if (GetTriangle(trigsperpoint.Get(p1,i)).HasEdge(p1,p2)) {return trigsperpoint.Get(p1,i);}
    }
  PrintSysError("ERROR in GetLeftTrig !!!");

  return 0;
}

int STLTopology :: GetRightTrig(int p1, int p2) const
{
  return GetLeftTrig(p2,p1);
}


int STLTopology :: NeighbourTrigSorted(int trig, int edgenum) const
{
  int i, p1, p2;
  int psearch = GetTriangle(trig).PNum(edgenum);

  for (i = 1; i <= 3; i++)
    {
      GetTriangle(trig).GetNeighbourPoints(GetTriangle(NeighbourTrig(trig,i)),p1,p2);
      if (p1 == psearch) {return NeighbourTrig(trig,i);}
    }

  PrintSysError("ERROR in NeighbourTrigSorted");
  return 0;
}






int STLTopology :: GetTopEdgeNum (int pi1, int pi2) const
{
  if (!ht_topedges) return 0;

  INDEX_2 i2(pi1, pi2);
  i2.Sort();

  if (!ht_topedges->Used(i2)) return 0;
  return ht_topedges->Get(i2);
}




void STLTopology :: InvertTrig (int trig)
{
  if (trig >= 1 && trig <= GetNT())
    {
      GetTriangle(trig).ChangeOrientation();
      FindNeighbourTrigs();
    }
  else
    {
      PrintUserError("no triangle selected!");
    }
}




void STLTopology :: DeleteTrig (int trig)
{
  if (trig >= 1 && trig <= GetNT())
    {
      trias.DeleteElement(trig);
      FindNeighbourTrigs();
    }
  else
    {
      PrintUserError("no triangle selected!");
    }
}



void STLTopology :: OrientAfterTrig (int trig)
{
  int starttrig = trig;

  if (starttrig >= 1 && starttrig <= GetNT())
    {

      ARRAY <int> oriented;
      oriented.SetSize(GetNT());
      int i;
      for (i = 1; i <= oriented.Size(); i++)
	{
	  oriented.Elem(i) = 0;
	}
 
      oriented.Elem(starttrig) = 1;
  
      int j = 0,k;
      
      ARRAY <int> list1;
      list1.SetSize(0);
      ARRAY <int> list2;
      list2.SetSize(0);
      list1.Append(starttrig);

      int cnt = 1;
      int end = 0;
      int nt;
      while (!end)
	{
	  end = 1;
	  for (i = 1; i <= list1.Size(); i++)
	    {
	      const STLTriangle& tt = GetTriangle(list1.Get(i));
	      for (k = 1; k <= 3; k++)
		{
		  nt = tt.NBTrigNum (k); // NeighbourTrig(list1.Get(i),k);
		  if (oriented.Get(nt) == 0)
		    {
		      if (tt.IsWrongNeighbourFrom(GetTriangle(nt)))
			{
			  GetTriangle(nt).ChangeOrientation();
			}
		      oriented.Elem(nt) = 1;
		      list2.Append(nt);
		      cnt++;
		      end = 0;
		    }
		}
	    }
	  list1.SetSize(0);
	  for (i = 1; i <= list2.Size(); i++)
	    {
	      list1.Append(list2.Get(i));
	    }
	  list2.SetSize(0);
	}

      PrintMessage(5,"NO corrected triangles = ",cnt);
      if (cnt == GetNT()) 
	{
	  PrintMessage(5,"ALL triangles oriented in same way!");
	}
      else
	{
	  PrintWarning("NOT ALL triangles oriented in same way!");
	}

      //      topedges.SetSize(0);
      FindNeighbourTrigs();
    }
  else
    {
      PrintUserError("no triangle selected!");
    }
}


}
