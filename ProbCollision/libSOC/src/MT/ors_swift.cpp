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

#ifdef MT_SWIFT

#ifdef MT_SINGLE
#  define SWIFT_USE_FLOAT
#endif

#include <SWIFT.h>
#undef min
#undef max

// internal use
void exportStateToSwift(const ors::Graph& C,SwiftModule& swift);
void importProxiesFromSwift(ors::Graph& C,SwiftModule& swift,bool dumpReport=false);
void swiftQueryExactDistance(SwiftModule& swift);

SwiftModule::~SwiftModule(){
  if(scene) delete scene;
  cout <<" -- SwiftModule closed" <<endl;
}

void SwiftModule::clone(const SwiftModule& S,const ors::Graph& G){
  init(G,S.cutoff);
}

void SwiftModule::init(const ors::Graph& C,double _cutoff){
  ors::Shape *s;
  uint k;
  bool r;

  cutoff=_cutoff;
  scene=new SWIFT_Scene(true, false);
  
  ids.resize(C.shapes.N);

  cout <<" -- SwiftModule init";
  for_list(k,s,C.shapes){
    cout <<'.' <<flush;
    r=false;
    switch(s->type){
    case BNONE: HALT("shapes should have a type - somehow wrong initialization..."); break;
    case BBOX:
      s->mesh.setBox();
      s->mesh.scale(s->size[0],s->size[1],s->size[2]);
      break;
    case BSPHERE:
      s->mesh.setSphere();
      s->mesh.scale(s->size[3],s->size[3],s->size[3]);
      break;
    case BCYLINDER:
      s->mesh.setCylinder(s->size[3],s->size[2]);
      break;
    case BCCYLINDER:
      s->mesh.setCappedCylinder(s->size[3],s->size[2]);
      break;
    case BMESH:
      //check if there is a specific swiftfile!
      MT::String *filename;
      filename=anyListGet<MT::String>(s->ats,"swiftfile",1);
      if(!filename)
        filename=anyListGet<MT::String>(s->body->ats,"swiftfile",1);
      if(filename){
        r=scene->Add_General_Object(*filename, ids(s->index), false);
        if(!r) HALT("--failed!");
      }
      break;
    default:
      NIY;
    }
    if(!r){
      r=scene->Add_Convex_Object(
        s->mesh.V.p, (int*)s->mesh.T.p,
        s->mesh.V.d0, s->mesh.T.d0, ids(s->index), false,
        DEFAULT_ORIENTATION, DEFAULT_TRANSLATION, DEFAULT_SCALE,
        DEFAULT_BOX_SETTING, DEFAULT_BOX_ENLARGE_REL, cutoff);
    }
    if(!r) HALT("--failed!");
  }

  initActivations(C);
  
  exportStateToSwift(C,*this);
  cout <<"...done" <<endl;
}

void SwiftModule::initActivations(const ors::Graph& C){
  ors::Shape *s,*s2;
  ors::Body *b,*b2;
  ors::Joint *e;
  uint j,k,k2;
  
  /* deactivate some collision pairs:
    -- no `cont' -> no collisions with this object at all
    -- no collisions between shapes of same object
    -- no collisions between linked objects
    -- no collisions between objects liked via the tree via 3 links
  */

  cout <<"collision active shapes: ";
  for_list(k,s,C.shapes) if(s->cont) cout <<s->body->name <<' ';
  
  for_list(k,s,C.shapes){
    if(!s->cont) scene->Deactivate( ids(s->index) );
    else         scene->Activate(   ids(s->index) );
  }
  //shapes within a body
  for_list(j,b,C.bodies){
    for_list(k,s,b->shapes) for_list(k2,s2,b->shapes){
      if(k2>k) scene->Deactivate( ids(s->index), ids(s2->index) );
    }
  }
  //deactivate along edges...
  for_list(j,e,C.joints){
    //cout <<"deactivating edge pair"; listWriteNames(TUPLE(e->from,e->to),cout); cout <<endl;
    deactivate(TUPLE(e->from,e->to));
  }
  //deactivate along trees...
  for_list(k,b,C.bodies){
    MT::Array<ors::Body*> group,children;
    group.append(b);
    for(uint l=0;l<3;l++){
      //listWriteNames(group,cout);
      children.clear();
      for_list(k2,b2,group){
        for_list(j,e,b2->outLinks){
          children.setAppend(e->to);
          //listWriteNames(children,cout);
        }
      }
      group.setAppend(children);
      //listWriteNames(group,cout);
    }
    //cout <<"deactivating group "; listWriteNames(group,cout); cout <<endl;
    deactivate(group);
  }
}

void SwiftModule::deactivate(const MT::Array<ors::Body*>& bodies){
  uint i1,i2,k1,k2;
  ors::Shape *s1,*s2;
  ors::Body *b1,*b2;
  for(i1=0;i1<bodies.N;i1++) for(i2=i1+1;i2<bodies.N;i2++){
    b1=bodies(i1);
    b2=bodies(i2);
    for_list(k1,s1,b1->shapes) for_list(k2,s2,b2->shapes){
      scene->Deactivate( ids(s1->index), ids(s2->index) );
    }
  }
}

void SwiftModule::deactivate(const ors::Shape *s1,const ors::Shape *s2){
  scene->Deactivate( ids(s1->index), ids(s2->index) );
}

void exportStateToSwift(const ors::Graph& C,SwiftModule& swift){
  ors::Shape *s;
  uint k;
  arr rot(3,3);
  for_list(k,s,C.shapes) {
    s->X.r.getMatrix(rot.p);
    swift.scene->Set_Object_Transformation( swift.ids(s->index), rot.p, s->X.p.v );
    if(!s->cont) swift.scene->Deactivate( swift.ids(s->index) );
    //else         swift.scene->Activate( swift.ids(s->index) );
  }
}

void importProxiesFromSwift(ors::Graph& C,SwiftModule& swift,bool dumpReport){
  int i,j,k,np;
  int *oids,*num_contacts;
  SWIFT_Real *dists, *nearest_pts, *normals;

  try{
    swift.scene->Query_Contact_Determination(
      false, swift.cutoff, np, 
      &oids, &num_contacts,
      &dists,
      &nearest_pts,
      &normals);
  }catch(const char *msg){
    cout <<"... catching error '"<<msg <<"' -- SWIFT failed! .. no proxies for this posture!!..." <<endl;
    return;
  }
  
  if(dumpReport){
    cout << "contacts: np=" <<np <<endl;
    for(k=0,i=0;i<np;i++){
      cout <<"* Object " <<C.shapes(oids[i<<1])->body->name <<" vs. Object " <<C.shapes(oids[(i<<1)+1])->body->name <<endl;
      cout <<"  #contacts = " <<num_contacts[i] <<endl;
      for(j=0;j<num_contacts[i];j++,k++){
        cout <<"  - contact " <<j <<endl;
        cout <<"    distance= " <<dists[k] << endl;
        cout <<"    points  = " <<nearest_pts[6*k+0] <<' '<<nearest_pts[6*k+1] <<' '<<nearest_pts[6*k+2] <<' '<<nearest_pts[6*k+3] <<' '<<nearest_pts[6*k+4] <<' '<<nearest_pts[6*k+5] << endl;
        cout <<"    normals = " <<normals[3*k+0] <<' ' <<normals[3*k+1] <<' ' <<normals[3*k+2] << endl;
      }
    }
  }

  //count total number of new proxies:
  for(k=0,i=0;i<np;i++){
    if(num_contacts[i]>=0) k+=num_contacts[i];
    if(num_contacts[i]==-1) k++;
  }

  //resize old list
  uint Nold=C.proxies.N;
  for(i=0;i<(int)Nold;i++) C.proxies(i)->age++;
  C.proxies.memMove=true;
  C.proxies.resizeCopy(Nold+k);
  for(i=0;i<k;i++) C.proxies(Nold+i) = new ors::Proxy;
  ors::Proxy *proxy;

  //add contacts to list
  int a,b;
  //  ors::Vector d,p;  d.setZero();
  for(k=0,i=0;i<np;i++){
    a=oids[i<<1];
    b=oids[(i<<1)+1];
    CHECK(swift.ids(a)==a && swift.ids(b)==b, "shape index does not coincide with swift index");
    //non-penetrating pair of objects
    if(num_contacts[i]>=0) for(j=0;j<num_contacts[i];j++,k++){
      proxy=C.proxies(Nold+k);
      proxy->a=a;
      proxy->b=b;
      proxy->age=0;
      proxy->d = dists[k];
      proxy->normal.set(&normals[3*k+0]);
      proxy->normal.normalize();
      //swift returns nearest points in the local frame -> transform them
      proxy->posA.set(&nearest_pts[6*k+0]);  proxy->posA = C.shapes(a)->X * proxy->posA;
      proxy->posB.set(&nearest_pts[6*k+3]);  proxy->posB = C.shapes(b)->X * proxy->posB;
      
      //if(a!=-1) proxy->velA=C.bodies(a).X.v + (C.bodies(a).X.w^(p+d-(C.bodies(a).X.p))); else proxy->velA.setZero();
      //if(b!=-1) proxy->velB=C.bodies(b).X.v + (C.bodies(b).X.w^(p-d-(C.bodies(b).X.p))); else proxy->velB.setZero();
      proxy->velA.setZero();
      proxy->velB.setZero();
      if(a!=-1 && b!=-1) proxy->rel.setDifference(C.shapes(a)->X,C.shapes(b)->X);
      else if(a!=-1) proxy->rel.setInverse(C.shapes(a)->X);
      else if(b!=-1) proxy->rel = C.shapes(b)->X;
      else           proxy->rel.setZero();
    }
    //penetrating pair of objects
    if(num_contacts[i]==-1){
      proxy=C.proxies(Nold+k);
      proxy->a=a;
      proxy->b=b;
      proxy->age=0;
      proxy->d = -.0;
      proxy->posA = C.shapes(a)->X.p;
      proxy->posB = C.shapes(b)->X.p;
      proxy->normal = proxy->posA - proxy->posB; //normal always points from b to a
      proxy->normal.normalize();
      //!! IN PENETRATION we measure d as -1+(distance between object centers) - that gives a well-defined (though non-smooth) gradient!
      //proxy->d += (proxy->posA-proxy->posB).length();
      //proxy->posA -= .5*proxy->normal;
      //proxy->posB += .5*proxy->normal;
      //CHECK(fabs(fabs(proxy->d+.1)-(proxy->posA-proxy->posB).length())<1e-10,"")
      k++;
//       cout << ".";
      //MT_MSG("WARNING - swift penetration!!!");
    }
  }
  CHECK(k+Nold == C.proxies.N,"");

  C.sortProxies();
}

void SwiftModule::computeProxies(ors::Graph& C,bool dumpReport){
  exportStateToSwift(C,*this);
  importProxiesFromSwift(C,*this,dumpReport);
}

void swiftQueryExactDistance(SwiftModule& swift){
  int i,np;
  int *oids;
  SWIFT_Real *dists;
  
  swift.scene->Query_Exact_Distance( false, SWIFT_INFINITY, np, &oids, &dists );
    
  cout << "exact distances: np=" <<np <<endl;
  for(i=0;i<np;i++){
    cout << "    Object " << oids[i<<1] << " vs. Object "
        << oids[(i<<1)+1] << " = " << dists[i] << endl;
  }
}

#else
#include "util.h"
//#warning "MT_SWIFT undefined - using HALT implementations"
void SwiftModule::init(const ors::Graph& C,double _cutoff){ MT_MSG("WARNING - creating dummy SwiftModule"); }
void SwiftModule::deactivate(const MT::Array<ors::Body*>& bodies){}
void SwiftModule::computeProxies(ors::Graph& C,bool dumpReport){}
void SwiftModule::clone(const SwiftModule& S,const ors::Graph& G){}
void swiftQueryExactDistance(SwiftModule& swift){}
SwiftModule::~SwiftModule(){}
#endif
