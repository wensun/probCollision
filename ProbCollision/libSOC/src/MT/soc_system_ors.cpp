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
#include "opengl.h"
#include "plot.h"
#include "soc.h"

//===========================================================================
//
// SocSystem_Ors
//

struct soc::SocSystem_Ors_Workspace{
  arr Qlin,Qoff,Qinv;
  arr q0,v0,W,H,Q,v_act;
  
  uint T;
  real tau;
  bool pseudoDynamic;
  bool newedOrs;
};

soc::SocSystem_Ors::SocSystem_Ors(){
  WS = new SocSystem_Ors_Workspace;
  ors   = NULL;
  swift = NULL;
  vars = NULL;
  WS->pseudoDynamic=false;
  WS->newedOrs=false;
}

soc::SocSystem_Ors::~SocSystem_Ors(){
  if(WS->newedOrs){
    listDelete(vars);
    delete swift;
    delete gl;
    delete ors;
  }
  delete WS;
}

void soc::SocSystem_Ors::clone(const SocSystem_Ors& sys){
  os = sys.os;
  dynamic = sys.dynamic;
  
  if(!sys.WS->newedOrs){
    gl = sys.gl;
    ors   = sys.ors;
    swift = sys.swift;
  }else{
    gl = new OpenGL;            gl->clone(*sys.gl);
    ors=new ::ors::Graph;     ors->clone(*sys.ors);
    swift=new SwiftModule;  swift->clone(*sys.swift,*ors);
    gl->remove(ors::glDrawGraph,sys.ors);
    gl->add   (ors::glDrawGraph,ors);
  }
  WS->newedOrs = sys.WS->newedOrs;

  listCopy(vars, sys.vars);
  uint i;  TaskVariable *y;
  for_list(i,y,vars) y->ors=ors;
  
  WS->Qlin = sys.WS->Qlin;
  WS->Qoff = sys.WS->Qoff;
  WS->Qinv = sys.WS->Qinv;
  WS->q0 = sys.WS->q0;
  WS->v0 = sys.WS->v0;
  WS->W = sys.WS->W;
  WS->H = sys.WS->H;
  WS->Q = sys.WS->Q;
  WS->v_act = sys.WS->v_act;
  WS->T = sys.WS->T;
  WS->tau = sys.WS->tau;
  WS->pseudoDynamic = sys.WS->pseudoDynamic=false;
}

void soc::SocSystem_Ors::initKinematic(ors::Graph *_ors,SwiftModule *_swift,OpenGL *_gl,uint trajectory_length, arr*W){
  ors   = _ors;
  swift = _swift;
  gl        = _gl;
  WS->T=trajectory_length;
  WS->tau=.01;
  ors->getJointState(WS->q0,WS->v0);
  if(W){
    if(W->nd==1){ CHECK(W->N==WS->q0.N,""); WS->W.setDiag(*W); }
    else NIY;
  }else{
    ors->computeNaturalQmetric(WS->W);
    cout <<"automatic W initialization =" <<diag(WS->W) <<endl;
//     graphWriteDirected(cout,ors->bodies,ors->joints);
  }
  static MT::Parameter<real> wc("Wcost");
  WS->W *= wc();
  WS->H = WS->W;
  dynamic=false;
}

void soc::SocSystem_Ors::initPseudoDynamic(ors::Graph *_ors,SwiftModule *_swift,OpenGL *_gl, real trajectory_time,  uint trajectory_steps, arr *W){
  ors   = _ors;
  swift = _swift;
  gl        = _gl;
  setTimeInterval(trajectory_time, trajectory_steps);
  setq0AsCurrent();
  if(W){
    if(W->nd==1){ CHECK(W->N==WS->q0.N,""); WS->W.setDiag(*W); }
    else NIY;
  }else{
    ors->computeNaturalQmetric(WS->W);
  }
  static MT::Parameter<real> hc("Hcost");
  static MT::Parameter<real> qn("Qnoise",1e-10);
  uint n=WS->q0.N;
  WS->H = hc()*WS->W;     //u-metric for torque control
  WS->Q.setDiag(qn,2*n);  //covariance \dot q-update
  dynamic=true;
  WS->pseudoDynamic=true;
}

void soc::SocSystem_Ors::initStandardReachProblem(uint rand_seed,uint T){
  MT::String orsFile = MT::getParameter<MT::String>("orsFile");
  if(!T) T = MT::getParameter<uint>("trajectoryLength");
  dynamic = MT::getParameter<bool>("isDynamic");
  real margin = MT::getParameter<real>("margin");
  MT::String endeffName= MT::getParameter<MT::String>("endeffName");
  MT::String endeffRel = MT::getParameter<MT::String>("endeffRel");
  bool standOnFoot = MT::getParameter<bool>("standOnFoot");
  bool useTruncation = MT::getParameter<bool>("useTruncation");

  
  //-- setup the simulator, swift, opengl
  WS->newedOrs=true;
  ors=new ors::Graph;
  ors->init(orsFile);
  if(standOnFoot){
    ors->reconfigureRoot(ors->getName("rfoot"));
    ors->calcNodeFramesFromEdges();
  }
  if(rand_seed>0){
    rnd.seed(rand_seed);
    ors::Body &t=*ors->getName("target");
    t.X.p(0) += .05*rnd.gauss();
    t.X.p(1) += .05*rnd.gauss();
    t.X.p(2) += .05*rnd.gauss();
  }
  
  swift=new SwiftModule;
  swift->init(*ors,margin);
  
  gl = new OpenGL;
  gl->add(glStandardScene);
  gl->add(ors::glDrawGraph,ors);
  gl->camera.setPosition(5,-10,10);
  gl->camera.focus(0,0,1);
  //gl->watch("loaded configuration - press ENTER");

  if(dynamic){
    initPseudoDynamic(ors,swift,gl,3.,T);
  }else{
    initKinematic(ors,swift,gl,T);
  }
  os=&std::cout;
  
  real endPrec=MT::getParameter<real>("endPrec");
  real midPrec=MT::getParameter<real>("midPrec");
  real colPrec=MT::getParameter<real>("colPrec");
  real balPrec=MT::getParameter<real>("balPrec");
   //-- setup the control variables (problem definition)
  TaskVariable *pos = new TaskVariable("position" , *ors, posTVT, endeffName, endeffRel,0,0, ARR());
  TaskVariable *col;
  if(!useTruncation) col = new TaskVariable("collision", *ors, collTVT,0,0,0,0,ARR(margin));
  else               col = new TaskVariable("collision", *ors, colConTVT,0,0,0,0,ARR(margin));
  TaskVariable *com = new TaskVariable("balance", *ors, comTVT,0,0,0,0,ARR());
  setTaskVariables(TUPLE(pos,col,com));
      
  pos->y_target = arr(ors->getName("target")->X.p.v,3);
  pos->setInterpolatedTargetsEndPrecisions(T,midPrec,endPrec,0.,10*endPrec);
  if(col->type==collTVT){
    col->y        = ARR(0.);
    col->y_target = ARR(0.);
    col->setInterpolatedTargetsConstPrecisions(T,colPrec,0.);
  }else col->active=true;
  if(balPrec){
    com->y_target = com->y;
    com->setInterpolatedTargetsConstPrecisions(T,balPrec,0.);
  }else com->active=false;
}

void soc::SocSystem_Ors::initStandardBenchmark(uint rand_seed){
  uint K = MT::getParameter<uint>("segments");
  uint T = MT::getParameter<uint>("trajectoryLength");
  dynamic = MT::getParameter<bool>("isDynamic");
  real margin = MT::getParameter<real>("margin");
  bool useTruncation = MT::getParameter<bool>("useTruncation");

  //generate the configuration
  ors::Body *b,*target,*endeff;  ors::Shape *s;  ors::Joint *j;
  MT::String str;
  ors=new ors::Graph;
  //the links
  for(uint k=0;k<=K;k++){
    b=new ors::Body(ors->bodies);
    b->name = STRING("body"<<k);
    if(!k) b->fixed=true;
    s=new ors::Shape(ors->shapes,b);
    s->type=2;
    s->size[0]=.0; s->size[1]=.0; s->size[2]=1./K; s->size[3]=.2/K;
    s->rel.setText(STRING("<t(0 0 "<<.5/K<<")>"));
    if(k&1){ s->color[0]=.5; s->color[1]=.2; s->color[2]=.2; }
    else   { s->color[0]=.2; s->color[1]=.2; s->color[2]=.2; }
    s->cont=true;
    if(k){
      j=new ors::Joint(ors->joints,ors->bodies(k-1),ors->bodies(k));
      j->type = JHINGE;
      j->Q.setText("<d(45 1 0 0)>");
      if(k&1){ //odd -> rotation around z
        j->A.setText(STRING("<t(0 0 "<<1./K<<") d(-90 0 1 0)>"));
        j->B.setText("<d(90 0 1 0)>");
      }else{  //even -> rotation around x
        j->A.setText(STRING("<t(0 0 "<<1./K<<")>"));
      }
    }
  }
  endeff=b;
  //the target
  b=new ors::Body(ors->bodies);
  b->name = "target";
  b->X.setText("<t(.2 0 0)>");
  s=new ors::Shape(ors->shapes,b);
  s->read(STREAM("type=1 size=[.0 .0 .1 .02] color=[0 0 1]"));
  s=new ors::Shape(ors->shapes,b);
  s->read(STREAM("type=0 rel=<t(0 -.1 0)> size=[.2 .01 .3 .0] color=[0 0 0] contact"));
  s=new ors::Shape(ors->shapes,b);
  s->read(STREAM("type=0 rel=<t(0 .1 0)> size=[.2 .01 .3 .0] color=[0 0 0] contact"));
  s=new ors::Shape(ors->shapes,b);
  s->read(STREAM("type=0 rel=<t(.1 0 0)> size=[.01 .2 .3 .0] color=[0 0 0] contact"));
  graphMakeLists(ors->bodies,ors->joints);
  ors->calcNodeFramesFromEdges();
  target=b;
  
  if(rand_seed>0){
    rnd.seed(rand_seed);
    target->X.p(0) += .05*rnd.gauss();
    target->X.p(1) += .05*rnd.gauss();
    //target->X.p(2) += .05*rnd.gauss();
  }
  
  swift=new SwiftModule;
  swift->init(*ors,margin);
  
  gl = new OpenGL;
  gl->add(glStandardScene);
  gl->add(ors::glDrawGraph,ors);
  gl->camera.setPosition(0,-5,5);
  gl->camera.focus(0,0,.5);
  gl->watch("loaded configuration - press ENTER");

  if(dynamic){
    initPseudoDynamic(ors,swift,gl,3.,T);
  }else{
    initKinematic(ors,swift,gl,T);
  }
  os=&std::cout;
  
  real endPrec=MT::getParameter<real>("endPrec");
  real midPrec=MT::getParameter<real>("midPrec");
  real colPrec=MT::getParameter<real>("colPrec");
   //-- setup the control variables (problem definition)
  TaskVariable *pos = new TaskVariable("position" , *ors, posTVT, endeff->name, STRING("<t(0 0 "<<.5/K<<")>"),0,0, ARR());
  TaskVariable *col;
  if(!useTruncation) col = new TaskVariable("collision", *ors, collTVT,0,0,0,0,ARR(margin));
  else               col = new TaskVariable("collision", *ors, colConTVT,0,0,0,0,ARR(margin));
  setTaskVariables(TUPLE(pos,col));
      
  pos->y_target = arr(ors->getName("target")->X.p.v,3);
  pos->setInterpolatedTargetsEndPrecisions(T,midPrec,endPrec,0.,10*endPrec);
  if(col->type==collTVT){
    col->y        = ARR(0.);
    col->y_target = ARR(0.);
    col->setInterpolatedTargetsConstPrecisions(T,colPrec,0.);
  }else col->active=true;
}

/* OLD VERSION!!
void soc::createEndeffectorReachProblem(SocSystem_Ors &sys,
                                        const char *ors_file,
                                        uint trajectory_length,
                                        int rand_seed)
{

  //setup the workspace
  MT::load(*sys.ors,ors_file);
  if(rand_seed>0){
    rnd.seed(rand_seed);
    ors::Body &t=*sys.ors->getName("target");
    t.X.p(0) += .05*rnd.gauss();
    t.X.p(1) += .05*rnd.gauss();
    t.X.p(2) += .05*rnd.gauss();
  }
  sys.ors->calcNodeFramesFromEdges();
  sys.ors->reconfigureRoot(sys.ors->getName("rfoot"));
  sys.ors->getJointState(sys.WS->q0,sys.WS->v0);
  sys.swift->init(*sys.ors);


  sys.WS->T=trajectory_length;
  sys.WS->tau=.01;
  arr Wdiag(sys.WS->q0.N);
  Wdiag=1.;
  MT_MSG("Warning - need to change this");
  Wdiag <<"[20 20 20 10 10 10 10 1 1 1 1 10 10 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 20 20 10 10 10 10 10 10 ]";
  //Wdiag <<"[20 20 20 10 10 10 10 1 1 1 1 10 10 1 1 1 20 20 10 10 10 10 10 10 ]";
  CHECK(Wdiag.N==sys.WS->q0.N,"wrong W matrix!");
  sys.WS->W.setDiag(Wdiag);
  sys.WS->H = sys.WS->W;

  //set task variables
  TaskVariable *x0,*x1,*x2;
  x0 = new TaskVariable("finger-tip",*sys.ors,posTVT ,"effector",0,0,0,0);
  x1 = new TaskVariable("COM",       *sys.ors,comTVT ,0,0,0,0,0);
  x2 = new TaskVariable("collision", *sys.ors,collTVT,0,0,0,0,ARR(.02));
  sys.vars = TUPLE(x0,x1,x2);

  updateState(globalSpace);
  //reportAll(globalSpace,cout);


  real midPrec,endPrec,balPrec,colPrec;
  MT::getParameter(midPrec,"midPrec");
  MT::getParameter(endPrec,"endPrec");
  MT::getParameter(balPrec,"balPrec");
  MT::getParameter(colPrec,"colPrec");

  x0->y_target.setCarray(sys.ors->getName("target")->X.p.v,3);
  x0->setInterpolatedTargetTrajectory(sys.WS->T);
  x0->setPrecisionTrajectoryFinal(sys.WS->T,midPrec,endPrec);

  x1->y_target(0)=0.;
  x1->y_target(1)=-.02;
  x1->setInterpolatedTargetTrajectory(sys.WS->T);
  x1->setPrecisionTrajectoryConstant(sys.WS->T,balPrec);
  if(!balPrec) x1->active=false;

  x2->y_target = 0.;
  x2->setInterpolatedTargetTrajectory(sys.WS->T);
  x2->setPrecisionTrajectoryConstant(sys.WS->T,colPrec);
  if(!colPrec) x2->active=false;

  sys.dynamic=false;
}
*/
    
void soc::SocSystem_Ors::setTimeInterval(real trajectory_time,  uint trajectory_steps){
  WS->T=trajectory_steps;
  WS->tau=trajectory_time/trajectory_steps;
}

void soc::SocSystem_Ors::setTaskVariables(const TaskVariableList& _CVlist){
  vars=_CVlist;
}

//! report on some 
void soc::SocSystem_Ors::reportOnState(ostream& os){
  os <<"OrsSocImplementat - state report:\n";
  os <<"** control variables:" <<endl;
  reportAll(vars,os);
  os <<"** proxies:" <<endl;
  ors->reportProxies(&os);
}

//overload the display method to include variances
void soc::SocSystem_Ors::displayState(const arr& q,const arr *invQ,const char *text){
  setq(q);
  if(text) gl->text.clr() <<text;
  if(invQ){
    arr Q;
    inverse_SymPosDef(Q,*invQ);
    if(gl->drawers.last().classP!=&plotModule) gl->add(glDrawPlot,&plotModule);
    plotClear();
    TaskVariable *v;
    uint i;
    for_list(i,v,vars){
      if(v->type==posTVT){
        arr X(3,3);
        X = v->J * Q * v->tJ;
        plotCovariance(v->y,X);
        break;
      }
    }
  }
  gl->update();
  //gl->watch();
  //gl->timedupdate(getTau()*(T-1)/(display-1));
  //if(invQ) gl->drawers.popLast();
}


uint soc::SocSystem_Ors::nTime(){  return WS->T>>scalePower; }
uint soc::SocSystem_Ors::nTasks(){ return vars.N; }
uint soc::SocSystem_Ors::qDim(){   if(!WS->Qlin.N) return WS->q0.N; else return WS->Qlin.d1; }
uint soc::SocSystem_Ors::uDim(){   return qDim(); }
uint soc::SocSystem_Ors::yDim(uint i){ return vars(i)->y.N; }
real soc::SocSystem_Ors::getTau(bool scaled){
  real tau=WS->tau;
  if(scaled) for(uint i=0;i<scalePower;i++) tau *=2.;
  return tau;
}

void soc::SocSystem_Ors::getq0 (arr& q) { if(!WS->Qlin.N) q=WS->q0; else q=WS->Qinv*(WS->q0-WS->Qoff); }
void soc::SocSystem_Ors::getv0 (arr& v) { if(!WS->Qlin.N) v=WS->v0; else v=WS->Qinv*WS->v0; }
void soc::SocSystem_Ors::getqv0(arr& q_){
   if(!WS->Qlin.N) q_.setBlockVector(WS->q0,WS->v0);
   else{
     arr q,qd;
     getqv0(q,qd);
     q_.setBlockVector(q,qd);
   }
}
void soc::SocSystem_Ors::getqv0(arr& q,arr& qd){ getq0(q); getv0(qd); }
void soc::SocSystem_Ors::getW  (arr& W) {
  if(!WS->Qlin.N) W=WS->W; else W = ~WS->Qlin*WS->W*WS->Qlin;
  if(scalePower){
    arr invW;
    inverse_SymPosDef(invW,W);
    for(uint i=0;i<scalePower;i++) invW = invW+invW;
    inverse_SymPosDef(W,invW);
  }
}
void soc::SocSystem_Ors::getWinv(arr& Winv) {
  arr W;
  if(!WS->Qlin.N) W=WS->W; else W = ~WS->Qlin*WS->W*WS->Qlin;
  inverse_SymPosDef(Winv,W);
  for(uint i=0;i<scalePower;i++) Winv = Winv+Winv;
}
void soc::SocSystem_Ors::getH  (arr& H) {
   if(!WS->Qlin.N) H=WS->H; else H = ~WS->Qlin*WS->H*WS->Qlin;
   for(uint i=0;i<scalePower;i++) H = H+H;
}
void soc::SocSystem_Ors::getQ  (arr& Q) {
  if(!WS->Qlin.N) Q=WS->Q; else{
    arr Qbig;
    Qbig.resize(2*WS->Qlin.d0,2*WS->Qlin.d1); Qbig.setZero();
    Qbig.setMatrixBlock(WS->Qlin,0,0); Qbig.setMatrixBlock(WS->Qlin,WS->Qlin.d0,WS->Qlin.d1);
    //cout <<Qbig <<endl;
    Q = ~Qbig*WS->Q*Qbig;
  }
}

void soc::SocSystem_Ors::setq(const arr& _q){
  arr q;
  if(!WS->Qlin.N) q.referTo(_q);
  else q = WS->Qlin*_q + WS->Qoff;
  ors->setJointState(q);
  ors->calcNodeFramesFromEdges();
  swift->computeProxies(*ors,false);
  WS->v_act.resizeAs(q);
  WS->v_act.setZero();
  uint i;
  TaskVariable *v;
  for_list(i,v,vars) if(v->active){
    v->updateState();
    v->updateJacobian();
  }
}

void soc::SocSystem_Ors::setqv(const arr& _q,const arr& _qd){
  arr q,qd;
  if(!WS->Qlin.N){ q.referTo(_q); qd.referTo(_qd); }
  else{ q = WS->Qlin*_q + WS->Qoff;  qd = WS->Qlin*_qd; }
  ors->setJointState(q,qd);
  ors->calcNodeFramesFromEdges();
  swift->computeProxies(*ors,false);
  WS->v_act=qd;
  uint i;
  TaskVariable *v;
  for_list(i,v,vars) if(v->active){
    v->updateState();
    v->updateJacobian();
  }
}

void soc::SocSystem_Ors::setqv(const arr& q_){
  uint n=q_.N/2;
  CHECK(q_.N==2*n,"");
  arr q,v;
  q.referToSubRange(q_,0,n-1);
  v.referToSubRange(q_,n,2*n-1);
  setqv(q,v);
}
void soc::SocSystem_Ors::setq0AsCurrent(){
  ors->getJointState(WS->q0,WS->v0);
}


void soc::SocSystem_Ors::getMF(arr& M,arr& F){
  if(!WS->pseudoDynamic){
    if(WS->Qlin.N) NIY;
    ors->clearForces();
    ors->gravityToForces();
    //ors->frictionToForces(1.1);
    ors->equationOfMotion(M,F,WS->v_act);
    //M.setId();  F = .1;
    //Minv *= .2;//1e-1;
  }else{
    uint n=qDim();
    M.setId(n);
    F.resize(n); F.setZero();
  }
}

void soc::SocSystem_Ors::getMinvF(arr& Minv,arr& F){
  arr M;
  getMF(M,F);
  inverse(Minv,M);
}

bool soc::SocSystem_Ors::isConditioned(uint i,uint t){
  if(vars(i)->type==colConTVT) return false;
  return vars(i)->active;
}

bool soc::SocSystem_Ors::isConstrained(uint i,uint t){
  return vars(i)->active && vars(i)->type==colConTVT;
}

const char* soc::SocSystem_Ors::taskName(uint i){
  return vars(i)->name.p;
}

void soc::SocSystem_Ors::getPhi(arr& phiq_i,uint i){
  //vars(i)->updateState();
  phiq_i=vars(i)->y;
}

void soc::SocSystem_Ors::getJqd(arr& jqd_i,uint i){
  arr q,qd;
  ors->getJointState(q,qd);
  //vars(i)->updateJacobian();
  jqd_i = vars(i)->J * qd;
}

void soc::SocSystem_Ors::getJtJ(arr& J_i,arr& tJ_i,uint i){
  //vars(i)->updateJacobian();
  if(!WS->Qlin.N){
    J_i=vars(i)->J;
    tJ_i=vars(i)->tJ;
  }else{
    J_i=vars(i)->J*WS->Qlin;
    transpose(tJ_i,J_i);
  }
}

void soc::SocSystem_Ors::getHessian(arr& H_i,uint i){
  if(WS->Qlin.N) NIY;
  vars(i)->getHessian(H_i);
}

void soc::SocSystem_Ors::getTarget(arr& y_i,real& y_prec,uint i,uint t){
  if(!t && vars(i)->targetType!=trajectoryTT){
    TaskVariable *v=vars(i);
    v->updateChange(-1);
    y_i  = v->y+v->y_change;
    y_prec = vars(i)->y_prec;
    return;
  }
  y_i    = vars(i)->y_trajectory[t];
  y_prec = vars(i)->y_prec_trajectory(t);
}

void soc::SocSystem_Ors::getTargetV(arr& v_i,real& y_prec,uint i,uint t){
  if(!t && vars(i)->targetType!=trajectoryTT){
    v_i    = vars(i)->v_target;
    y_prec = vars(i)->v_prec;
    return;
  }
  v_i    = vars(i)->v_trajectory[t];
  y_prec = vars(i)->v_prec_trajectory(t);
}


//===========================================================================
//
// problem implementations
//

void drawOrsSocEnv(void*){
  glStandardLight();
  //glDrawFloor(1.,.4,.4,.4);
  //DrawAxes(1.);
}

void soc::setupOpenGL(SocSystem_Ors &sys){
  if(!sys.gl) sys.gl=new OpenGL();
  sys.gl->add(drawOrsSocEnv,0);
  sys.gl->add(ors::glDrawGraph,sys.ors);
  //sys.gl->add(plotDrawOpenGL,&plotData);
  sys.gl->camera.focus(0,0,.8);
}


void soc::createDynamicProblem(SocSystem_Ors &sys,
                          const char *ors_file,
                          real trajectory_time,
                          uint trajectory_steps){

  //setup the workspace
  MT::load(*sys.ors,ors_file);
  sys.ors->calcNodeFramesFromEdges();
  sys.ors->getJointState(sys.WS->q0,sys.WS->v0);
  sys.swift->init(*sys.ors);

  uint n=sys.WS->q0.N;
  sys.WS->T=trajectory_steps;
  sys.WS->tau=trajectory_time/trajectory_steps;
  sys.WS->W.setDiag(1e-6,n);  //q-metric for inverse kinematics (initialization)
  static MT::Parameter<real> hc("Hcost");
  static MT::Parameter<real> qn("Qnoise");
  sys.WS->H.setDiag(hc,n);  //u-metric for torque control
  sys.WS->Q.setDiag(qn,2*n);  //covariance \dot q-update

  //set task variables
  TaskVariable *x0;
  x0 = new TaskVariable("finger-tip",*sys.ors,posTVT ,"eff","t(0 0 .15)",0,0,0);
  sys.vars = TUPLE(x0);

  updateState(sys.vars);
  //reportAll(globalSpace,cout);

  real midPrec,endPrec;
  MT::getParameter(midPrec,"midPrec");
  MT::getParameter(endPrec,"endPrec");

  x0->y_target.setCarray(sys.ors->getName("target")->X.p.v,3);
  x0->v_target <<"[2 0 0]";
  x0->setInterpolatedTargetTrajectory(sys.WS->T);
  x0->setPrecisionTrajectoryFinal (sys.WS->T,midPrec,endPrec);
  x0->setPrecisionVTrajectoryFinal(sys.WS->T,0.,endPrec);
  

  sys.dynamic=true;
}

void createNikolayReachProblem(soc::SocSystem_Ors& sys,
                                   ors::Graph& _ors,
                                   SwiftModule& _swift,
                                   uint trajectory_length,
                                   const arr& endeffector_target,
                                   const char* endeffector_name,
                                   const arr& W){
  static soc::SocSystem_Ors_Workspace WS;

  //setup the workspace
  //WS.vars = globalSpace;
  sys.ors=&_ors;
  sys.ors->getJointState(WS.q0,WS.v0);
  WS.T=trajectory_length;
  WS.W=W;
  sys.swift=&_swift;

  //set task variables
  TaskVariable *x0,*x1;
  x0 = new TaskVariable("finger-tip",*sys.ors,posTVT ,endeffector_name,"",0,0,0);
  x1 = new TaskVariable("collision", *sys.ors,collTVT,0,0,0,0,0);
  sys.vars = TUPLE(x0,x1);

  updateState(sys.vars);
  //reportAll(globalSpace,cout);

  x0->y_prec=1e3;  x0->setGainsAsAttractor(30.);  x0->y_target = endeffector_target;
  x0->setTrajectory(WS.T,0,1e10);
  x0->setPrecisionTrajectoryFinal(WS.T,1e1,1e4);

  x1->y_prec=1e6;  x1->setGainsAsAttractor(30.);  x1->y_target = 0.;
  x1->setTrajectory(WS.T,0);
  x1->setPrecisionTrajectoryConstant(WS.T,1e5);

  sys.WS=&WS;
  sys.dynamic=false;
}
