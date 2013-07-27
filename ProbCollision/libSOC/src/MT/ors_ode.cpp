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

#ifdef MT_ODE

#ifndef dDOUBLE
#  define dDOUBLE
#endif

#  include <ode/ode.h>
#  include <ode/../../ode/src/objects.h>
#  include <ode/../../ode/src/joints/joints.h>
#  include <ode/../../ode/src/collision_kernel.h>
#  include <ode/../../ode/src/collision_transform.h>
#  ifdef MT_MSVC
#    undef HAVE_UNISTD_H
#    undef HAVE_SYS_TIME_H
#  endif

#define OUTs(x) x <<' ' 
#define OUTv(x) '(' <<OUTs(x[0]) <<OUTs(x[1]) <<OUTs(x[2]) <<')'
#define OUTq(x) OUTs(x[0]) <<OUTs(x[1]) <<OUTs(x[2]) <<OUTs(x[3])
#define CP4(x,y) memmove(x,y,4*sizeof(double));
#define CP3(x,y) memmove(x,y,3*sizeof(double));

static bool ODEinitialized=false;

//===========================================================================
//
// Ode implementations
//

OdeModule::OdeModule(){
  time=0.;
  

  noGravity=noContactJoints=false;
  
  ERP= 0.2;     //in [0,1]: rate of error correction (makes more brittle) [default=.2]  
  CFM=1e-3;  // >0: softness (makes more robust) [default=1e-10]
  
	coll_bounce = .0;
	coll_ERP = 0.2;   //usually .2!! stiffness (time-scale of contact reaction)
	coll_CFM = 1e+1;  //softness
	friction = 0.1;   //alternative: dInfinity;
	
	world=0;
  space=0;
  contactgroup=0;
  
  if(!ODEinitialized){  dInitODE();  ODEinitialized=true; }
  clear();
}

OdeModule::~OdeModule(){
  if(contactgroup) dJointGroupDestroy(contactgroup);
  if(space) dSpaceDestroy(space);
  if(world) dWorldDestroy(world);
}

void OdeModule::clear(){
  if(contactgroup) dJointGroupDestroy(contactgroup);
  if(space) dSpaceDestroy(space);
  if(world) dWorldDestroy(world);
  
  world=dWorldCreate();
  space=dSimpleSpaceCreate(0);
  contactgroup=dJointGroupCreate(0);
  //std::cout <<"default ERP="<<dWorldGetERP(world) <<"default CFM=" <<dWorldGetCFM(world) <<std::endl;
  if(noGravity){
    dWorldSetGravity (world,0,0,0);
  }else{
    dWorldSetGravity (world,0,0,-9.81);
  }
  dWorldSetCFM(world,CFM);
  dWorldSetERP(world,ERP);
  dWorldSetContactSurfaceLayer(world,.0); //necessary penetration depth to trigger forces..
  double DAMPING_LINEAR_SCALE = 0.03;
  double DAMPING_ANGULAR_SCALE = 0.03;
  dWorldSetDamping(world, DAMPING_LINEAR_SCALE, DAMPING_ANGULAR_SCALE);

  plane0=dCreatePlane(space,0,0,1,0);
 
  // removed vertical planes,  5. Mar 06 (hh)
  //planex1=dCreatePlane(space, 1,0,0,-10);
  //planex2=dCreatePlane(space,-1,0,0,-10);
  //planey1=dCreatePlane(space,0, 1,0,-10);
  //planey2=dCreatePlane(space,0,-1,0,-10);
}

void OdeModule::staticCallback(void *classP, dGeomID g1, dGeomID g2){
  static dContact contacts[100]; // array with maximum number of contacts
  uint i,n;
  dJointID c;

  dBodyID b1 = dGeomGetBody(g1);
  dBodyID b2 = dGeomGetBody(g2);

  ors::Body *db1 = b1?(ors::Body*)b1->userdata:0;
  ors::Body *db2 = b2?(ors::Body*)b2->userdata:0;

  //-- rule out irrelevant contacts
  // exit without doing anything if the geoms stem from the same body
  if(b1==b2) return;

  // exit without doing anything if the two bodies are connected by a joint  
  if(b1 && b2 && dAreConnectedExcluding(b1,b2,dJointTypeSlider)) return;  

  // exit if fixed body intersects with earth,  4. Mar 06 (hh)
  if(b1==0 && db2->fixed==true) return;
  if(b2==0 && db1->fixed==true) return;

  // exit if we have two fixed bodies,  6. Mar 06 (hh)
  if(db1 && db2 && db1->fixed==true && db2->fixed==true) return;

  // exit if none of the bodies have cont enabled (mt)
  //if(db1 && !db1->cont && db2 && !db2->cont) return;
  if(db1 && db2 && (!db1->shapes(0)->cont || !db2->shapes(0)->cont)) return;

  //-- ok, now compute the contact exactly
  // get contacts:
  n=dCollide(g1,g2,100,&(contacts[0].geom),sizeof(dContact));

  if(!n) return;

  for(i=0;i<n;i++) ((OdeModule*)classP)->conts.append(&contacts[i].geom);

  if(((OdeModule*)classP)->noContactJoints) return;

  for(i=0;i<n;i++){
    contacts[i].surface.mode = dContactSoftCFM | dContactSoftERP | dContactBounce;
    //if(b1 && b2)
      contacts[i].surface.mu	 = ((OdeModule*)classP)->friction; //dInfinity; //friction
    //else contacts[i].surface.mu	 = .01; //dInfinity; //friction
	contacts[i].surface.bounce   = ((OdeModule*)classP)->coll_bounce;
	contacts[i].surface.soft_erp = ((OdeModule*)classP)->coll_ERP;   //usually .2!! stiffness (time-scale of contact reaction)
	contacts[i].surface.soft_cfm = ((OdeModule*)classP)->coll_CFM;  //softness

  if(!db1 || !db2 || (db1->shapes(0)->cont && db2->shapes(0)->cont)){
      //normal contact:
      c=dJointCreateContact(((OdeModule*)classP)->world,((OdeModule*)classP)->contactgroup,&contacts[i]);
      dJointAttach(c,b1,b2);
    }else{
      //special contact:
      //one of them ignores the contact
      //hence we consider it (kinematically) coupled to the background
      //.. we simulated relative motion for (parallel) friction
      ors::Vector no; no.set(contacts[i].geom.normal);
      ors::Vector d1; d1=db1->X.v; d1.makeNormal(no); //bug: we should have sorted b1 and b2 before such that b1!=0
      double v=d1.length();
      if(v>1e-4){
        d1 /= v;
        CP3(contacts[i].fdir1,d1.v);
        contacts[i].surface.mode = contacts[i].surface.mode | dContactMotion1 | dContactFDir1;
	contacts[i].surface.motion1 = v;
      }
      c=dJointCreateContact(((OdeModule*)classP)->world,((OdeModule*)classP)->contactgroup,&contacts[i]);
      if(!db1->shapes(0)->cont){
        dJointAttach(c,0,b2);
      }
      if(!db2->shapes(0)->cont){
        dJointAttach(c,b1,0);
      }
    }
  }
}

void OdeModule::setForceFree(bool free){
  noGravity=free;
  if(free){
    dWorldSetGravity (world,0,0,0);
    dWorldSetCFM(world,CFM);
    dWorldSetERP(world,ERP);
  }else{
    dWorldSetGravity (world,0,0,-9.81);
    dWorldSetCFM(world,CFM);
    dWorldSetERP(world,ERP);
  }
}

void OdeModule::step(double dtime){
  conts.clear();
  dSpaceCollide(space,this,&staticCallback);
  //if(noContactJoints) contactForces();
  //dWorldStep(world,dtime);
  //dWorldQuickStep(world,dtime); 
  dWorldStep(world,dtime); 
  dJointGroupEmpty(contactgroup);
  time+=dtime;
}

void OdeModule::reportContacts(){
  std::cout <<"contacts: " <<conts.N <<std::endl;
  dContactGeom *c;
  for(uint i=0;i<conts.N;i++){
    c = conts(i);
    dBodyID b1 = dGeomGetBody(c->g1);
    dBodyID b2 = dGeomGetBody(c->g2);
    ors::Body *db1 = b1?(ors::Body*)b1->userdata:0;
    ors::Body *db2 = b2?(ors::Body*)b2->userdata:0;
    std::cout
      <<i <<": " <<(b1?db1->name.p:"GROUND") <<'-' <<(b2?db2->name.p:"GROUND")
      <<"\nposition=" <<OUTv(c->pos)
      <<"\nnormal=" <<OUTv(c->normal)
      <<"\npenetration depth=" <<OUTs(c->depth)
      <<std::endl;
  }
}

// - compute depth of penetration and normal for
//   body with touch sensor, 6. Mar 06 (hh).
// - fixed bug (I did not assume that each contact is
//   just represented once), 12. Mar 06 (hh).
// - copy penetration data to Body, 22. May 06 (hh)
void OdeModule::penetration(ors::Vector &p){
  NIY;
#if 0
  dContactGeom *c;
  double d;
  dBodyID b1,b2;

  for(uint i=0;i<conts.N;i++){
    c = conts(i);
    b1 = dGeomGetBody(c->g1);  
    b2 = dGeomGetBody(c->g2);  
    //get body attributes
    double *touch1=0,*touch2=0;
    if(b1) touch1=anyListGet<double>(((ors::Body*)b1->userdata)->ats,"touchsensor",0);
    if(b2) touch2=anyListGet<double>(((ors::Body*)b2->userdata)->ats,"touchsensor",0);

    if(touch1 || touch2){
      d = c->depth;
      p.set(c->normal);
      p *= d;
      if (touch1) ((ors::Body*)b1->userdata)->touch = p;
      else ((ors::Body*)b2->userdata)->touch = -p;
      //printf("%lf\n",d);
      std::cout << "contact "<<i<<": penetration: " << OUTv(p.v) << "\n";
    }
  }
#endif
}


void OdeModule::contactForces(){
  uint i;
  dContactGeom *c;
  dBodyID b1,b2;
  ors::Vector pos,normal,v1,v2,vrel;
  dMass mass;
  double force,d,m1,m2;
  for(i=0;i<conts.N;i++){
    c = conts(i);
    pos.set(c->pos);
    d=c->depth*20.;
    normal.set(c->normal);
    normal *=-1.; //vector from b1 -> b2
    b1 = dGeomGetBody(c->g1);
    b2 = dGeomGetBody(c->g2);
    if(b1){ v1.set(b1->lvel); dBodyGetMass(b1,&mass); m1=mass.mass; }else{ v1.setZero(); m1=0.; }
    if(b2){ v2.set(b1->lvel); dBodyGetMass(b2,&mass); m2=mass.mass; }else{ v2.setZero(); m2=0.; }
    
    // ** normal force:
    force=0.;
    vrel=v2-v1; //inaccurate: should also include rotational velocities at contact point!
    vrel.makeColinear(normal);
    //spring force
    force = .03*d*d*d;
    if(b2) dBodyAddForceAtPos(b2,force*normal(0),force*normal(1),force*normal(2),pos(0),pos(1),pos(2));
    force *= -1.; //invert on other body
    if(b1) dBodyAddForceAtPos(b1,force*normal(0),force*normal(1),force*normal(2),pos(0),pos(1),pos(2));
    
    //inward viscosity
    if(vrel * normal <= 0){
      force += 10.*vrel.length();
      if(b2) dBodyAddForceAtPos(b2,force*m2*normal(0),force*m2*normal(1),force*m2*normal(2),pos(0),pos(1),pos(2));
      force *= -1.; //invert on other body
      if(b1) dBodyAddForceAtPos(b1,force*m1*normal(0),force*m1*normal(1),force*m1*normal(2),pos(0),pos(1),pos(2));
    }
    
    if(!b2) std::cout <<"bodyForce = " <<force*normal <<std::endl;
    
#if 1
    // ** parallel (slip) force:
    force=0.;
    vrel=v2-v1; //inaccurate: should also include rotational velocities at contact point!
    vrel.makeNormal(normal);
    //viscosity
    force += 10.*vrel.length();
    
    vrel.normalize();
    force *= -1.;
    if(b2) dBodyAddForceAtPos(b2,force*m2*vrel(0),force*m2*vrel(1),force*m2*vrel(2),pos(0),pos(1),pos(2));
    force *= -1.;
    if(b1) dBodyAddForceAtPos(b1,force*m1*vrel(0),force*m1*vrel(1),force*m1*vrel(2),pos(0),pos(1),pos(2));
    
    if(!b2) std::cout <<"bodyForce = " <<force*normal <<std::endl;
    std::cout <<"slip force " <<force <<std::endl;
#endif
  }
}



//===========================================================================
//
// graph
//

void OdeModule::exportStateToOde(ors::Graph &C){
  ors::Graph::node n;
  ors::Graph::edge e;
  uint j,i;
  dBodyID b;
  
  C.calcNodeFramesFromEdges();
  
  for_list(i,n,C.bodies){
    CHECK(n->X.r.isNormalized(),"quaternion is not normalized!");
    b = bodies(n->index);
    CP3(b->posr.pos,n->X.p.v);                // dxBody changed in ode-0.6 ! 14. Jun 06 (hh)
    CP4(b->q,n->X.r.q); dQtoR(b->q,b->posr.R);
    CP3(b->lvel,n->X.v.v);
    CP3(b->avel,n->X.w.v);
    //n->copyFrameToOde();
#ifndef MT_ode_nojoints
    for_list(j,e,n->inLinks) if(e->type!=JGLUE){
      dxJointHinge* hj=(dxJointHinge*)joints(e->index);
      dxJointUniversal* uj=(dxJointUniversal*)joints(e->index);
      dxJointSlider* sj=(dxJointSlider*)joints(e->index);
      switch (e->type){ //16. Mar 06 (hh)
      case JHINGE: 
      default:
	CP4(hj->qrel,(e->A.r*e->B.r).q);
	CP3(hj->anchor1,(e->A.p).v);
	CP3(hj->anchor2,(-(e->B.r/e->B.p)).v);
	CP3(hj->axis1,(e->A.r*ors::Vector(1,0,0)).v);
	CP3(hj->axis2,(e->B.r/ors::Vector(1,0,0)).v);
	break;
      case JFIXED: 
	break;
      case JUNIVERSAL: 
	CP4(uj->qrel1,(e->A.r).q);
	CP4(uj->qrel2,(e->B.r).q);
	CP3(uj->anchor1,(e->A.p).v);
	CP3(uj->anchor2,(-(e->B.r/e->B.p)).v);
	CP3(uj->axis1,(e->A.r*ors::Vector(1,0,0)).v);
	CP3(uj->axis2,(e->B.r/ors::Vector(1,0,0)).v);
	break;	
      case JSLIDER:
	CP4(sj->qrel,(e->A.r*e->B.r).q);
	CP3(sj->offset,(e->Q.p).v );
	//printf("offset: (%lf; %lf; %lf)\n",(e->X.p)[0],(e->X.p)[1],(e->X.p)[2]);
	CP3(sj->axis1,(e->A.r*ors::Vector(1,0,0)).v);  
	break;
      }
    }
#endif
  }
}

void OdeModule::exportForcesToOde(ors::Graph &C){
  ors::Graph::node n;
  uint j;
  dBodyID b;
  for_list(j,n,C.bodies){
    b=bodies(n->index);
    CP3(b->facc,n->force.v);
    CP3(b->tacc,n->torque.v);
  }
}

void OdeModule::importStateFromOde(ors::Graph &C){
  ors::Graph::node n;
  uint j;
  dBodyID b;
  for_list(j,n,C.bodies){
    //n->getFrameFromOde();
    b=bodies(n->index);
    CP3(n->X.p.v,b->posr.pos);
    CP4(n->X.r.q,b->q);
    CP3(n->X.v.v,b->lvel);
    CP3(n->X.w.v,b->avel);
    CHECK(n->X.r.isNormalized(),"quaternion is not normalized!");
  }
  C.calcEdgesFromNodeFrames();
}


void OdeModule::setJointForce(ors::Graph& C,ors::Graph::edge e, double f1, double f2){
  NIY; //make internal!
/*  switch (e->type){ //  3. Apr 06 (hh) 
    default:
    case JHINGE:
      dJointAddHingeTorque(joints(e->index),-f1);
      break;
    case JFIXED: // no torque 
      break;
    case JUNIVERSAL:
      dJointAddUniversalTorques(joints(e->index), -f1, -f2);
      break;
    case JSLIDER:
      dJointAddSliderForce(joints(e->index), -f1);
      break;
    }
    */
}

void OdeModule::setJointForce(ors::Graph& C,doubleA& x){
  NIY; //make internal!
/*  edge e;
  uint i=0,n=0;
  
  for_edges(i,e){ // loop over edges, 16. Mar 06 (hh)
    switch (e->type){ //  3. Apr 06 (hh) 
      default:
      case JHINGE:
	dJointAddHingeTorque(joints(e->index),-x(n));
	n++;
	break;
      case JFIXED: // no torque 
	break;
      case JUNIVERSAL:
	dJointAddUniversalTorques(joints(e->index), -x(n), -x(n+1));
	n+=2;
	break;
      case JSLIDER:
	dJointAddSliderForce(joints(e->index), -x(n));
	n++;
	break;
      }
  }
  CHECK(n==x.N,"wrong dimensionality");
  */
}

uint OdeModule::getJointMotorDimension(ors::Graph& C){
  ors::Graph::node n; ors::Graph::edge e;
  uint i=0,j,k;
  for_list(k,n,C.bodies){
    for_list(j,e,n->inLinks) if(!e->fixed && e->motor){ // 16. Mar 06 (hh)
      if (e->type==JUNIVERSAL) i+=2;
      else i++;
    }
  }
  return i;
}

void OdeModule::setJointMotorPos(ors::Graph &C,ors::Graph::edge e,double x0,double maxF,double tau){
  double x=-dJointGetHingeAngle(joints(e->index));
  setJointMotorVel(C,e,.1*(x0-x)/tau,maxF);
}

void OdeModule::setJointMotorVel(ors::Graph &C,ors::Graph::edge e,double v0,double maxF){
  dJointSetAMotorParam(motors(e->index),dParamVel,-v0);
  dJointSetAMotorParam(motors(e->index),dParamFMax,maxF);
}

void OdeModule::setJointMotorPos(ors::Graph &C,doubleA& x,double maxF,double tau){
  //CHECK(x.N==E,"given joint state has wrong dimension, "<<x.N<<"="<<E);
  ors::Graph::node n;
  ors::Graph::edge e;
  uint i=0,j,k;
  for_list(k,n,C.bodies){
    for_list(j,e,n->inLinks) if(e->motor) {
      setJointMotorPos(C,e,x(i),maxF,tau);
      i++;
      break;
    }
  }
  CHECK(x.N==i,"joint motor array had wrong dimension " <<x.N<<"!="<<i);
}

void OdeModule::setJointMotorVel(ors::Graph &C,doubleA& x,double maxF){
  //CHECK(x.N==E,"given joint state has wrong dimension, "<<x.N<<"="<<E);
  ors::Graph::node n;
  ors::Graph::edge e;
  uint i=0,j,k;
  for_list(k,n,C.bodies){
    for_list(j,e,n->inLinks) if(e->motor) {
      setJointMotorVel(C,e,x(i),maxF);
      i++;
      break;
    }
  }
  CHECK(x.N==i,"joint motor array had wrong dimension " <<x.N<<"!="<<i);
}

void OdeModule::unsetJointMotors(ors::Graph &C){
  ors::Graph::node n;
  ors::Graph::edge e;
  uint i,j;
  for_list(j,n,C.bodies) for_list(i,e,n->inLinks) if(e->motor){
    unsetJointMotor(C,e);
    break;
  }
}

void OdeModule::unsetJointMotor(ors::Graph &C,ors::Graph::edge e){
  dJointSetAMotorParam(motors(e->index),dParamFMax,0.);
}

void OdeModule::getJointMotorForce(ors::Graph &C,ors::Graph::edge e,double& f){
  dJointFeedback* fb=dJointGetFeedback(motors(e->index));
  CHECK(fb,"no feedback buffer set for this joint");
  //std::cout <<OUTv(fb->f1) <<' ' <<OUTv(fb->t1) <<' ' <<OUTv(fb->f2) <<' ' <<OUTv(fb->t2) <<std::endl;
  ors::Vector t; t(0)=fb->t1[0]; t(1)=fb->t1[1]; t(2)=fb->t1[2];
  f=t * (e->from->X.r*(e->A.r*ors::Vector(1,0,0)));
  f=-f;
}

void OdeModule::getJointMotorForce(ors::Graph &C,doubleA& f){
  ors::Graph::node n;
  ors::Graph::edge e;
  uint i=0,j,k;
  for_list(k,n,C.bodies){
    for_list(j,e,n->inLinks) if(!e->fixed && e->motor) {
      getJointMotorForce(C,e,f(i));
      i++;
      break;
    }
  }
  CHECK(f.N==i,"joint motor array had wrong dimension " <<f.N<<"!="<<i);
}

void OdeModule::pidJointPos(ors::Graph &C,ors::Graph::edge e,double x0,double v0,double xGain,double vGain,double iGain,double* eInt){
  double x,v,f;
  //double a,b,c,dAcc;

  x=-dJointGetHingeAngle(joints(e->index));
  v=-dJointGetHingeAngleRate(joints(e->index));
  f = xGain*(x0-x) + vGain*(v0-v);
  if(eInt){
    f+=iGain* (*eInt);
    (*eInt) = .99*(*eInt) + .01 *(x0-x);
  }
  std::cout <<"PID:" <<x0 <<' ' <<x <<' ' <<v <<" -> " <<f <<std::endl;
  dJointSetAMotorParam(motors(e->index),dParamFMax,0);
  dJointAddHingeTorque(joints(e->index),-f);
}

void OdeModule::pidJointVel(ors::Graph &C,ors::Graph::edge e,double v0,double vGain){
  double v,f;

  v=-dJointGetHingeAngleRate(joints(e->index));
  f = vGain*(v0-v);
  if(fabs(f)>vGain) f=MT::sign(f)*vGain;
  std::cout <<"PIDv:" <<v0 <<' ' <<v <<" -> " <<f <<std::endl;
  dJointSetAMotorParam(motors(e->index),dParamFMax,0);
  dJointAddHingeTorque(joints(e->index),-f);
}


//===========================================================================
//
// higher level C
//

//static dGeomID lastGeomHack;

void OdeModule::createOde(ors::Graph &C){
  ors::Graph::node n;
  ors::Graph::edge e;
  dBodyID b;
  dMass odeMass;
  dGeomID geom,trans;
  dSpaceID space;
  dxJointHinge* jointH=0;
  dxJointSlider* jointS=0;
  dxJointUniversal* jointU=0;
  dxJointFixed* jointF=0;
  dxJointAMotor* jointM=0;
  dJointFeedback* jointFB=0;
  ors::Vector a;
  //double *mass,*shape,*type,*fixed,*cont,typeD=BCCYLINDER;
  //, *inertiaTensor, *realMass, *centerOfMass; 
  uint i,j;
  //trimeshPhysics triPhys;
  dTriMeshDataID TriData;

  clear();
  C.calcNodeFramesFromEdges();

  bodies.resize(C.bodies.N); bodies=0;
  geoms .resize(C.shapes.N); geoms =0;
  joints.resize(C.joints.N); joints=0;
  motors.resize(C.joints.N); motors=0;
  
  for_list(j,n,C.bodies){
    b=dBodyCreate(world);

    bodies(n->index)=b;
    b->userdata=n;

    //n->copyFrameToOde();
    CHECK(n->X.r.isNormalized(),"quaternion is not normalized!");
    CP3(b->posr.pos,n->X.p.v);                // dxBody changed in ode-0.6 ! 14. Jun 06 (hh)
    CP4(b->q,n->X.r.q); dQtoR(b->q,b->posr.R);
    CP3(b->lvel,n->X.v.v);
    CP3(b->avel,n->X.w.v);

    // need to fix: "mass" is not a mass but a density (in dMassSetBox)
    ors::Shape *s;
    for_list(i,s,n->shapes){
      if(!(s->rel.r.isZero()) || !(s->rel.p.isZero())){ //we need a relative transformation
        trans = dCreateGeomTransform(space);
        space = NULL;
      }else{
        trans = NULL;
        space = space;
      }

      switch(s->type){ 
        default: case BBOX: // box
          dMassSetBox(&odeMass, n->mass, s->size[0], s->size[1], s->size[2]);
          dBodySetMass(b,&odeMass);
          geom=dCreateBox(space, s->size[0], s->size[1], s->size[2]);
          break;
          case BSPHERE: // sphere
            dMassSetSphere(&odeMass, n->mass, s->size[3]);
            dBodySetMass(b,&odeMass);
            geom=dCreateSphere(space,s->size[3]);
            break;
            case BCYLINDER: // capped cylinder, 6. Mar 06 (hh)
              dMassSetCylinder(&odeMass, n->mass, 3, s->size[3], s->size[2]);
              dBodySetMass(b,&odeMass);
              geom=dCreateCylinder(space, s->size[3], s->size[2]);
              break;
              case BCCYLINDER: // capped cylinder, 6. Mar 06 (hh)
                dMassSetCylinder(&odeMass, n->mass, 3, s->size[3], s->size[2]);
                MT_MSG("ODE: setting Cylinder instead of capped cylinder mass");
                dBodySetMass(b,&odeMass);
                geom=dCreateCCylinder(space, s->size[3], s->size[2]);
                break;
                case BMESH: // trimesh loaded from file, 12. Jun 06 (hh & dm)
                  NIY;
                  s->mesh.computeNormals();
                  // get inertia tensor and REAL mass (careful no density)
                  //n->ats.get("I",inertiaTensor,9);
                  //n->ats.get("w",realMass,1);
                  //n->ats.get("X",centerOfMass,3);
            
                  // transform the mesh to ODE trimesh format;
                  //i=0; j=0;
  
                  /*
                  triPhys.reset(s->mesh.numtriangles);
                  if (inertiaTensor && realMass && centerOfMass) // are all important params set in the dcg file
                  {
                  triPhys._mass = *realMass;
                  triPhys.r[0]=centerOfMass[0];
                  triPhys.r[1]=centerOfMass[1];
                  triPhys.r[2]=centerOfMass[2];
                  triPhys.J[0][0]= inertiaTensor[0];
                  triPhys.J[1][1]= inertiaTensor[4];
                  triPhys.J[2][2]= inertiaTensor[8];
                  triPhys.J[0][1]= inertiaTensor[1];
                  triPhys.J[0][2]= inertiaTensor[2];
                  triPhys.J[1][2]= inertiaTensor[5];
                  }
    
                  else // not all parametrs specified in dcg file....need to calculate them
                  {
                  triPhys.calculateODEparams(&s->mesh, s->mass); // note: 2nd parameter is the density
                  }
  
                  dMassSetParameters (&odeMass, triPhys._mass,
                  triPhys.r[0], triPhys.r[1], triPhys.r[2],
                  triPhys.J[0][0], triPhys.J[1][1], triPhys.J[2][2],
                  triPhys.J[0][1], triPhys.J[0][2], triPhys.J[1][2]);  
    
                  dBodySetMass(b,&odeMass);
                  */
                  dMassSetBox(&odeMass, n->mass, s->size[0], s->size[1], s->size[2]);
                  dBodySetMass(b,&odeMass);
  
  
                  if(true){
                    TriData = dGeomTriMeshDataCreate();
                    dGeomTriMeshDataBuildDouble(TriData,
                        s->mesh.V.p, 3*sizeof(double), s->mesh.V.d0,
                                              s->mesh.T.p, s->mesh.T.d0,3*sizeof(uint));
                    //dGeomTriMeshDataPreprocess(TriData);
  
        
                    geom = dCreateTriMesh(space, TriData, 0, 0, 0);
        
                    dGeomTriMeshClearTCCache(geom);
                    //lastGeomHack = geom;
                  }
                  /*else{
                    geom = lastGeomHack;
                    geom->gflags = 28;
                  }*/
                  break; //end of mesh
      }

      geoms(s->index) = geom;
      if(trans){
        //geoms(s->index) = trans;
        dGeomTransformSetGeom(trans,geom);
        dGeomSetPosition(geom,s->rel.p(0),s->rel.p(1),s->rel.p(2));
        dGeomSetQuaternion(geom,s->rel.r.q);
        dGeomSetBody(trans,b); //attaches the geom to the body
      }else{
        dGeomSetBody(geom,b); //attaches the geom to the body
      }
    }//loop through shapes
    
    if(n->fixed){
      jointF=(dxJointFixed*)dJointCreateFixed(world,0);
      dJointAttach(jointF,b,0);
      dJointSetFixed(jointF);
    }
  }
#ifndef MT_ode_nojoints
  for_list(j,n,C.bodies){
    for_list(i,e,n->inLinks){
      switch (e->type){
	default:
	case JHINGE:
	  jointH=(dxJointHinge*)dJointCreateHinge(world,0);
	  if(e->p[1]!=e->p[0]){
	    dJointSetHingeParam(jointH,dParamLoStop,e->p[0]);
	    dJointSetHingeParam(jointH,dParamHiStop,e->p[1]);

	    //dJointSetHingeParam(jointH,dParamCFM,CFM);
	  }
	  dJointAttach(jointH,bodies(e->from->index),bodies(e->to->index));
	  joints(e->index)=jointH;
	  //e->copyFramesToOdeHinge();
	  break;
	case JUNIVERSAL:
	  jointU=(dxJointUniversal*)dJointCreateUniversal(world,0);
	  dJointAttach(jointU,bodies(e->from->index),bodies(e->to->index));
	  joints(e->index)=jointU;
	  //e->copyFramesToOdeUniversal();
	  break;
	case JFIXED:
	  jointF=(dxJointFixed*)dJointCreateFixed(world,0);
	  dJointAttach(jointF,bodies(e->from->index),bodies(e->to->index));
	  dJointSetFixed(jointF);
	  joints(e->index)=jointF;
	  e->fixed=true;
	  break;
	case JSLIDER:
	  jointS=(dxJointSlider*)dJointCreateSlider(world,0);
	  dJointAttach(jointS,bodies(e->from->index),bodies(e->to->index));
	  joints(e->index)=jointS;
	  //e->copyFramesToOdeSlider();
	  break;
	}
      if (!e->fixed && e->motor) {
	jointM = (dxJointAMotor*)dJointCreateAMotor(world,0);
	dJointSetAMotorNumAxes(jointM,1);
	dJointAttach(jointM,bodies(e->from->index),bodies(e->to->index));
	motors(e->index)=jointM;
	a=e->from->X.r*e->A.r*ors::Vector(1,0,0);
	dJointSetAMotorAxis(jointM,0,1,a(0),a(1),a(2));
	//dJointSetAMotorParam(jointM,dParamFMax,1.);
	jointFB = new dJointFeedback;
	dJointSetFeedback(jointM,jointFB);
      }
    }
  }
#endif

  exportStateToOde(C);
}

void OdeModule::step(ors::Graph &C,doubleA& in,doubleA& force,doubleA& out,uint steps){
  C.setFullState(in);
  exportStateToOde(C);
  setJointForce(C,force);
  for(;steps--;) step();
  importStateFromOde(C);
  C.getFullState(out);
}

void OdeModule::step(ors::Graph &C,doubleA& force,doubleA& out,uint steps){
  setJointForce(C,force);
  for(;steps--;) step();
  importStateFromOde(C);
  C.getFullState(out);
}

void OdeModule::step(ors::Graph &C,uint steps,double tau){
  for(;steps--;) step(tau);
  importStateFromOde(C);
}

void OdeModule::getGroundContact(ors::Graph &C,boolA& cts){
  cts.resize(C.bodies.N);
  cts=false;
  //reportContacts();
  uint i;
  dContactGeom *c;
  for(i=0;i<conts.N;i++){
    c = conts(i);
    dBodyID b1,b2;
    b1 = dGeomGetBody(c->g1); 
    b2 = dGeomGetBody(c->g2); 

    if(!b1 && b2) cts(((ors::Graph::node)b2->userdata)->index)=true;
    if(b1 && !b2) cts(((ors::Graph::node)b1->userdata)->index)=true;
  }
}

/*struct dContactGeom {
  dVector3 pos;       // contact position
  dVector3 normal;    // normal vector
  dReal depth;        // penetration depth
  dGeomID g1,g2;      // the colliding geoms
  };*/


void OdeModule::importProxiesFromOde(ors::Graph &C){
  uint i,Nold=C.proxies.N;
  //Nold=0;
  for(i=0;i<Nold;i++) C.proxies(i)->age++;
  C.proxies.memMove=true;
  C.proxies.resizeCopy(Nold+conts.N);
  for(i=0;i<conts.N;i++) C.proxies(Nold+i) = new ors::Proxy;
  dContactGeom *c;
  int a,b;
  ors::Vector d,p;
  for(i=0;i<conts.N;i++){
    c = conts(i);
    dBodyID b1,b2;
    b1 = dGeomGetBody(c->g1); 
    b2 = dGeomGetBody(c->g2); 

    a = b1 ? (((ors::Graph::node)b1->userdata)->index) : (uint)-1;
    b = b2 ? (((ors::Graph::node)b2->userdata)->index) : (uint)-1;

    d.set(c->normal);
    d *= -c->depth/2.;
    p.set(c->pos);

    C.proxies(Nold+i)->a = a;
    C.proxies(Nold+i)->b = b;
    C.proxies(Nold+i)->d = -c->depth;
    C.proxies(Nold+i)->normal.set(c->normal);
    C.proxies(Nold+i)->posA=p+d;
    C.proxies(Nold+i)->posB=p-d;
    if(a!=-1) C.proxies(Nold+i)->velA=C.bodies(a)->X.v + (C.bodies(a)->X.w^(p+d-(C.bodies(a)->X.p))); else C.proxies(Nold+i)->velA.setZero();
    if(b!=-1) C.proxies(Nold+i)->velB=C.bodies(b)->X.v + (C.bodies(b)->X.w^(p-d-(C.bodies(b)->X.p))); else C.proxies(Nold+i)->velB.setZero();
    C.proxies(Nold+i)->posB=p-d;
    if(a!=-1 && b!=-1) C.proxies(Nold+i)->rel.setDifference(C.bodies(a)->X,C.bodies(b)->X);
    else if(a!=-1) C.proxies(Nold+i)->rel.setInverse(C.bodies(a)->X);
    else if(b!=-1) C.proxies(Nold+i)->rel = C.bodies(b)->X;
    else           C.proxies(Nold+i)->rel.setZero();
    C.proxies(Nold+i)->age=0;
  }

  C.sortProxies();
}

void OdeModule::reportContacts2(){
  uint i;
  dBodyID b1,b2;
  ors::Body* b;
  dContactGeom* c;
  ors::Vector x;
  std::cout <<"contacts=" <<conts.N <<std::endl;
  for(i=0;i<conts.N;i++){
    c = conts(i);
    std::cout <<i <<' ';
    b1=dGeomGetBody(c->g1);
    b2=dGeomGetBody(c->g2);
    if(b1){
      b=(ors::Body*)b1->userdata;
      std::cout <<b->name <<' ';
    }else std::cout <<"NIL ";
    if(b2){
      b=(ors::Body*)b2->userdata;
      std::cout <<b->name <<' ';
    }else std::cout <<"NIL ";
    x.set(c->pos);    std::cout <<"pos=" <<x <<' ';
    x.set(c->normal); std::cout <<"normal=" <<x <<' ';
    std::cout <<"depth=" <<c->depth <<std::endl;
  }
}

struct Bound{ ors::Vector p,n; };
bool OdeModule::inFloorContacts(ors::Vector& x){
  uint i,j,k;
  dBodyID b1,b2;
  dContactGeom* c;
  ors::Vector y;

  x(2)=0.;

  MT::Array<ors::Vector> v;
  //collect list of floor contacts
  for(i=0;i<conts.N;i++){
    c = conts(i);
    b1=dGeomGetBody(c->g1);
    b2=dGeomGetBody(c->g2);
    if((!b1 && b2) || (b1 &&!b2)){
      y.set(c->pos);  y(2)=0.;
      v.append(y);
    }
  }

  std::cout <<"\nfloor points: ";
  for(i=0;i<v.N;i++) std::cout <<v(i) <<'\n';

  if(!v.N) return false;

  //construct boundaries
  Bound b;
  MT::Array<Bound> bounds;
  double s;
  bool f;
  for(i=0;i<v.N;i++) for(j=0;j<i;j++){
    b.p=v(i); b.n=(v(j)-v(i)) ^ ors::Vector(0,0,1);
    f=false;
    for(k=0;k<v.N;k++) if(k!=j && k!=i){
      s=(v(k)-b.p) * b.n ;
      if(s<0){
	if(f) break;
	b.n=-b.n;
      }
      if(s!=0.) f=true;
    }
    if(k==v.N) bounds.append(b);
  }

  std::cout <<"\nbounds: ";
  for(i=0;i<bounds.N;i++) std::cout <<bounds(i).p <<' ' <<bounds(i).n <<'\n';

  std::cout <<"\nquery: " <<x <<std::endl;
  //check for internal
  for(i=0;i<bounds.N;i++){
    if((x-bounds(i).p) * bounds(i).n < 0) return false;
  }
  return true;
}

void OdeModule::slGetProxies(ors::Graph &C){
  //C.setJointState(x);
  //dJointGroupEmpty(contactgroup);
  //exportStateToOde(C,ode);
  createOde(C);
  conts.clear();
  noContactJoints=true;
  dSpaceCollide(space,this,OdeModule::staticCallback);
  importProxiesFromOde(C);
}

/*void OdeModule::slGetProxyGradient(arr &dx,const arr &x,ors::Graph &C){
  if(C.proxies.N){
    arr dp,J;
    C.getContactGradient(dp);
    C.jacobianAll(J,x);
    dx = 2.* (~J) * dp;
  }else{
    dx.resize(x.N);
    dx.setZero();
  }
}

void OdeModule::slGetProxyGradient(arr &dx,const arr &x,ors::Graph &C,OdeModule &ode){
  slGetProxies(x,C,ode);
  slGetProxyGradient(dx,x,C);
}*/

#undef OUTs
#undef OUTv
#undef OUTq
#undef CP4
#undef CP3


//===========================================================================
//
// communication with Ode
//

#if 0 //documentation...
/*!\brief copy all frame variables (positions, velocities, etc)
into the ODE engine (createOde had to be called before) */
void OdeModule::exportStateToOde(ors::Graph &C);

/*!\brief copy the current state of the ODE engine into
the frame variables (positions, velocities, etc) */
void OdeModule::importStateFromOde(ors::Graph &C);

/*!\brief copy all frame variables (positions, velocities, etc)
into the ODE engine (createOde had to be called before) */
void OdeModule::exportForcesToOde(ors::Graph &C);

void OdeModule::setJointForce(ors::Graph &C,ors::Graph::edge e,double f1,double f2);

/*!\brief applies forces on the configuration as given by the force vector x */
void OdeModule::setJointForce(ors::Graph &C,arr& f);

/*!\brief returns the number of motors attached to joints */
uint getJointMotorDimension(ors::Graph &C);

/*!\brief set the desired positions of all motor joints */
void OdeModule::setJointMotorPos(ors::Graph &C,arr& x,double maxF,double tau);
/*!\brief set the desired position of a specific motor joint */
void OdeModule::setJointMotorPos(ors::Graph &C,ors::Graph::edge e,double x0,double maxF,double tau);

/*!\brief set the desired velocities of all motor joints */
void OdeModule::setJointMotorVel(ors::Graph &C,arr& v,double maxF=1.);
/*!\brief set the desired velocity of a specific motor joint */
void OdeModule::setJointMotorVel(ors::Graph &C,ors::Graph::edge e,double v0,double maxF=1.);

/*!\brief disable motors by setting maxForce parameter to zero */
void OdeModule::unsetJointMotors(ors::Graph &C,OdeModule& ode);
/*!\brief disable motor by setting maxForce parameter to zero */
void OdeModule::unsetJointMotor(ors::Graph &C,ors::Graph::edge e);

/*!\brief get the forces induced by all motor motor joints */
void OdeModule::getJointMotorForce(ors::Graph &C,arr& f);
/*!\brief get the force induced by a specific motor joint */
void OdeModule::getJointMotorForce(ors::Graph &C,ors::Graph::edge e,double& f);

void OdeModule::pidJointPos(ors::Graph &C,ors::Graph::edge e,double x0,double v0,double xGain,double vGain,double iGain=0,double* eInt=0);

void OdeModule::pidJointVel(ors::Graph &C,ors::Graph::edge e,double v0,double vGain);

//! [obsolete] \ingroup sl
void OdeModule::getGroundContact(ors::Graph &C,boolA& cts);

//! import the information from ODE's contact list into the proximity list \ingroup sl
void OdeModule::importProxiesFromOde(ors::Graph &C,OdeModule& ode);

/*!\brief simulates one time step with the ODE engine, starting from state
    vector in, returns state vector out \ingroup sl */
void OdeModule::step(ors::Graph &C,arr& in,arr& force,arr& out,uint steps=1);

/*!\brief simulates one time step with the ODE engine, starting from state
    vector in, returns state vector out \ingroup sl */
void OdeModule::step(ors::Graph &C,arr& force,arr& out,uint steps=1);

//! \ingroup sl
void OdeModule::step(ors::Graph &C,uint steps=1,double tau=.01);

/*!\brief instantiate the configuration in an ODE engine (first clears the
    ODE engine from all other objects) \ingroup sl */
void OdeModule::createOde(ors::Graph &C,OdeModule& ode);

//! \ingroup sl
void OdeModule::slGetProxies(ors::Graph &C,OdeModule &ode);

//! \ingroup sl
//void OdeModule::slGetProxyGradient(arr &dx,const arr &x,ors::Graph &C,OdeModule &ode);

//! \ingroup sl
void OdeModule::reportContacts(OdeModule& ode);
//! \ingroup sl
bool inFloorContacts(ors::Vector& x);
#endif

#else
OdeModule::OdeModule(){ MT_MSG("WARNING - creating dummy OdeModule"); }
OdeModule::~OdeModule(){}
void OdeModule::OdeModule::step(double dtime){}
void OdeModule::createOde(ors::Graph &C){}
void OdeModule::slGetProxies(ors::Graph &C){}
void OdeModule::unsetJointMotors(ors::Graph &C){}
void OdeModule::exportStateToOde(ors::Graph &C){}
void OdeModule::importStateFromOde(ors::Graph &C){}
void OdeModule::importProxiesFromOde(ors::Graph &C){}
#endif
