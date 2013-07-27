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

#include "soc.h"
#include "opengl.h"
#include "ors.h"

//===========================================================================
//
// SocSystem_Toy
//

// mass on a spring
struct soc::SocSystem_Toy_Workspace{
  real m, d; //mass and spring constant
  real x,v;  //current position and velocity
  real x0,v0;//initial position and velocity

  arr W,H,Q;
  uint T;
  real tau;
};

soc::SocSystem_Toy::SocSystem_Toy(){
  WS = new SocSystem_Toy_Workspace;
  WS->m=1.;
  WS->d=1.3;
  WS->x0=1.;
  WS->v0=0.;

  WS->W.setId(1);
  static MT::Parameter<real> hc("Hcost");
  static MT::Parameter<real> qn("Qnoise");
  WS->H.setDiag(hc,1);
  WS->Q.setDiag(qn,2);
  WS->Q.setZero();
}

soc::SocSystem_Toy::~SocSystem_Toy(){
  delete WS;
}

uint soc::SocSystem_Toy::nTime(){  return WS->T; }
uint soc::SocSystem_Toy::nTasks(){ return 1; }
uint soc::SocSystem_Toy::qDim(){   return 1; }
uint soc::SocSystem_Toy::uDim(){   return 1; }
uint soc::SocSystem_Toy::yDim(uint i){ return 1; }
real soc::SocSystem_Toy::getTau(bool scaled){  return WS->tau; }

void soc::SocSystem_Toy::getq0 (arr& q) { q.resize(1); q(0)=WS->x0; }
void soc::SocSystem_Toy::getqv0(arr& q_){ q_.resize(2); q_(0)=WS->x0; q_(1)=WS->v0; }
void soc::SocSystem_Toy::getqv0(arr& q,arr& qd){ q.resize(1); q(0)=WS->x0; qd.resize(1); qd(0)=WS->v0; }
void soc::SocSystem_Toy::getW  (arr& W) { W=WS->W; }
void soc::SocSystem_Toy::getH  (arr& H) { H=WS->H; }
void soc::SocSystem_Toy::getQ  (arr& Q) { Q=WS->Q; }

void soc::SocSystem_Toy::setq(const arr& q){
  CHECK(q.N==1,"");
  WS->x=q(0);
  WS->v=0.;
}

void soc::SocSystem_Toy::setqv(const arr& q_){
  CHECK(q_.N==2,"");
  WS->x=q_(0);
  WS->v=q_(1);
}

void soc::SocSystem_Toy::setqv(const arr& q,const arr& qd){
  CHECK(q.N==1 && qd.N==1,"");
  WS->x=q(0);
  WS->v=qd(0);
}

void soc::SocSystem_Toy::setq0AsCurrent(){
  WS->x0=WS->x;
  WS->v0=WS->v;
}

void soc::SocSystem_Toy::getMF(arr& M,arr& F){
  M.resize(1,1);
  M(0,0)=WS->m;
  F.resize(1);
  F(0)=-WS->d*WS->x;
}

void soc::SocSystem_Toy::getMinvF(arr& Minv,arr& F){
  Minv.resize(1,1);
  Minv(0,0)=1./WS->m;
  F.resize(1);
  F(0)=-WS->d*WS->x;
}

bool soc::SocSystem_Toy::isConditioned(uint i,uint t){
  CHECK(i==0,"");
  return true;
}

void soc::SocSystem_Toy::getPhi(arr& phiq_i,uint i){
  CHECK(i==0,"");
  phiq_i.resize(1);
  phiq_i(0)=WS->x;
}

void soc::SocSystem_Toy::getJqd(arr& jqd_i,uint i){
  CHECK(i==0,"");
  jqd_i.resize(1);
  jqd_i(0)=WS->v;
}

void soc::SocSystem_Toy::getJtJ(arr& J_i,arr& tJ_i,uint i){
  CHECK(i==0,"");
  J_i .resize(1,1); J_i (0,0)=1.;
  tJ_i.resize(1,1); tJ_i(0,0)=1.;
}

void soc::SocSystem_Toy::getTarget(arr& y_i,uint i,uint t){
  y_i.resize(1);
  y_i=0; //-1.;
}

void soc::SocSystem_Toy::getTargetV(arr& v_i,uint i,uint t){
  v_i.resize(1);
  v_i=0; //3.;
}

void soc::SocSystem_Toy::getPrecision(real& prec,uint i,uint t){
  static MT::Parameter<real> ep("endPrec");
  if(t==WS->T-1) prec=ep;
  else prec=0.;
  //prec=0.;
}

void soc::SocSystem_Toy::getPrecisionV(real& prec,uint i,uint t){
  static MT::Parameter<real> ep("endPrec");
  if(t==WS->T-1) prec=ep;
  else prec=0.;
  //prec=0.;
}



//===========================================================================
//
// problem implementations
//

void toyDrawEnv(void *p){
#ifdef MT_GL
  soc::SocSystem_Toy_Workspace *WS = (soc::SocSystem_Toy_Workspace*)p;
  glStandardLight();
  glDrawAxes(1.);
  glTranslatef(0,0,WS->x);
  glDrawSphere(.1);
#else
  NIY;
#endif
}

void soc::setupOpenGL(SocSystem_Toy &soci){
  if(!soci.gl) soci.gl=new OpenGL;
  soci.gl->add(toyDrawEnv,soci.WS);
  soci.gl->camera.focus(0,0,.8);
}


void soc::createDynamicProblem(SocSystem_Toy &soci,
                          const char *ors_file,
                          real trajectory_time,
                          uint trajectory_steps){
  MT_MSG("*** TOY problem");
  soci.WS->T=trajectory_steps;
  soci.WS->tau=trajectory_time/trajectory_steps;
}
