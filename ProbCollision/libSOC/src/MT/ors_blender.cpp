/*  Copyright 2009 Marc Toussaint
    email: mtoussai@cs.tu-berlin.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a COPYING file of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/> */

#include "ors.h"
//#include "array.cpp"

//#include<fstream>
//<iostream
struct Group{
  String name;
  uintA verts;
};

//template class MT::Array<Group>;
//template class MT::Array<String>;

void readBlender(const char* filename,ors::Mesh& mesh,ors::Graph& bl){
  ifstream is(filename,std::ios::binary);
  CHECK(is.good(),"couldn't open file " <<filename);
  
  arr vertices,normals,frames,tailsHeads;
  uintA faces;
  MT::Array<Group> G;
  uintA graph;
  MT::Array<MT::String> names;
  
  String tag,name;
  char c;
  uint i,j,b;
  arr x0,x1,x2,x3;

  String::readSkipSymbols="\"\n\r\t ";
  String::readEatStopSymbol = true;

  for(;;){
    CHECK(is.good(),"error in scanning the file (previous tag = `" <<tag <<"'");
    String::readStopSymbols="\n\r\t ";  //space terminates tags
    is >>tag;
    String::readStopSymbols="\n\r\t\"";       //quotes terminate names
    if(tag=="vertices"){   vertices.readTagged(is,0);   continue;  }
    if(tag=="normals"){    normals.readTagged(is,0);    continue;  }
    if(tag=="faces"){      faces.readTagged(is,0);      continue;  }
    if(tag=="frames"){     frames.readTagged(is,0);     continue;  }
    if(tag=="tailsHeads"){ tailsHeads.readTagged(is,0); continue;  }
    if(tag=="groups"){
      is >>(const char*)"<" >>j >>(const char*)">";
      G.resize(j);
      for(i=0;i<j;i++){
        is >>G(i).name;
        G(i).verts.readTagged(is,0);
      }
      continue;
    }
    if(tag=="graph"){
      is >>(const char*)"{" >>i >>j >>(const char*)"}";
      graph.resize(j,2);
      c=is.get(); CHECK(c=='\n',"couldn't read newline after ascii tag :-(");
      is.read((char*)graph.p,graph.sizeT*graph.N);
      c=is.get(); CHECK(c=='\n',"couldn't read newline after array buffer :-(");
      continue;
    }
    if(tag=="names"){
      is >>(const char*)"<" >>j >>(const char*)">";
      //fscanf(is.rdbuf()->_File,"<%d>",&j);
      names.resize(j);
      for(i=0;i<j;i++) is >>names(i);
      continue;
    }
    if(tag=="quit"){ CHECK(is.good(),"not perfect import..."); break; }
    HALT("unknown tag `" <<tag <<"'");
  }

  ors::Vector *w;
  ors::Quaternion r; r.setDeg(0,1,0,0); //don't rotate the mesh
  for(i=0;i<vertices.d0;i++){
    w = (ors::Vector*)&vertices(i,0);
    *w = r*(*w);
  }

  mesh.V=vertices;
  mesh.T=faces;
  mesh.Tn=normals;
  mesh.G.resize(vertices.d0); mesh.G=-1;
  for(i=0;i<G.N;i++){
    for(b=0;b<names.N;b++){ if(G(i).name==names(b)) break; }
    if(b==names.N){
      cout <<"unknown body: " <<G(i).name <<endl;
    }else{
      for(j=0;j<G(i).verts.N;j++) mesh.G(G(i).verts(j))=b;
    }
  }

  mesh.computeNormals();

  double v[4],l;
  ors::Body *n,*p;
  ors::Shape *s;
  String parent;
  ors::Graph::edge e;
  ors::Vector h,t;
  ors::Frame f;
  ors::Quaternion ROT; ROT.setDeg(90,1,0,0); //rotate the armature

  for(i=0;i<frames.d0;i++){
    n=new ors::Body(bl.bodies);
    s=new ors::Shape(bl.shapes,n); //always create a shape for a body...
    MT::skip(is);
    n->name=names(i);
    f.p.set(&frames(i,3,0)); f.p=ROT*f.p;
    f.r.setMatrix(frames[i].sub(0,2,0,2).p);
    f.r.invert();
    f.addRelativeRotationDeg(90,1,0,0);   f.r=ROT*f.r;
    t.set(&tailsHeads(i,0,0)); t=ROT*t;
    h.set(&tailsHeads(i,1,0)); h=ROT*h;

#if 0
    n->X.p = f.p;
    n->cog = f.r/(t-h)/2.;
#else
    n->X.p = f.p + (t-h)/2.;
#endif
    n->X.r = f.r;

    s->type=0;
    l=(t-h).length();
    v[0]=v[1]=v[3]=l/20.; v[2]=l;
    memmove(s->size,v,4*sizeof(double));
  }
  for(i=0;i<graph.d0;i++){
    p=bl.bodies(graph(i,0));
    n=bl.bodies(graph(i,1));
    //e=new_edge(p,n,bl.bodies,bl.joints);
    e=new ors::Joint;
    e->index=bl.joints.N;
    bl.joints.append(e);
    e->ifrom=p->index; e->ito=n->index;
    f.p.set(&frames(graph(i,1),3,0));  f.p=ROT*f.p;
    f.r.setMatrix(frames[graph(i,1)].sub(0,2,0,2).p);
    f.r.invert();
    f.addRelativeRotationDeg(90,1,0,0);  f.r=ROT*f.r;

    e->A.setDifference(p->X,f);
    e->B.setDifference(f,n->X); //p=(h-t)/2;
  }
  graphMakeLists(bl.bodies,bl.joints);
  mesh.GF.resize(bl.bodies.N);
  for(i=0;i<bl.bodies.N;i++){
    mesh.GF(i) = &(bl.bodies(i)->X);
  }
  mesh.makeVerticesRelativeToGroup();
  mesh.collectTriGroups();
  //bl.calcNodeFramesFromEdges();

  /*
  bl.topsort();
  bl.orderAsIndexed();
  bl.indexAllAsOrdered();
  */
  MT_MSG("warning - structure is not sorted!");

  String::readEatStopSymbol = false;
}

