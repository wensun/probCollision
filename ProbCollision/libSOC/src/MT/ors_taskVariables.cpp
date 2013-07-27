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

TaskVariable::TaskVariable(){
  active=false;
  type=noneTVT;
  targetType=noneTT;
  y_prec=0.; v_prec=0.; a=b=0.; state=-1; state_tol=.05; err=derr=0.;
}

TaskVariable::TaskVariable(
    const char* _name,
    ors::Graph& _sl,
    TVtype _type,
    const char *iname,const char *iframe,
    const char *jname,const char *jframe,
    const arr& _params){
  active=false;
  type=noneTVT;
  targetType=noneTT;
  y_prec=0.; v_prec=0.; a=b=0.; state=-1; state_tol=.05;
  set(
    _name,_sl,_type,
    iname  ? (int)_sl.getName(iname)->index      : -1,
    iframe ? ors::Frame().setText(iframe) : ors::Frame(),
    jname  ? (int)_sl.getName(jname)->index      : -1,
    jframe ? ors::Frame().setText(jframe) : ors::Frame(),
    _params);
}

TaskVariable::~TaskVariable(){
}

void TaskVariable::set(
    const char* _name,
    ors::Graph &_ors,
    TVtype _type,
    int _i,const ors::Frame& _irel,
    int _j,const ors::Frame& _jrel,
    const arr& _params){
  type=_type;
  name=_name;
  ors=&_ors;
  i=_i;
  irel=_irel;
  j=_j;
  jrel=_jrel;
  params=_params;
  updateState();
  updateJacobian();
  y_target=y;
  v_target=v;
}

/*void TaskVariable::set(const char* _name,ors::Graph& _sl,TVtype _type,const char *iname,const char *jname,const char *reltext){
  set(
    _name,_sl,_type,
    _sl.getName(iname)->index,
    _sl.getName(jname)->index,
    ors::Frame().setText(reltext));
}*/

void TaskVariable::setGains(real _a,real _b){
  targetType=gainsTT;
  active=true;
  a=_a;
  b=_b;
  if(!y_prec) y_prec=100.;
}

void TaskVariable::setGainsAsAttractor(real decaySteps,real oscillations){
  targetType=attractorTT;
  active=true;
  b=2./decaySteps;
  a=b*b*(MT_PI*MT_PI*oscillations*oscillations + .25);
  if(!y_prec) y_prec=100.;
}

//compute an y_trajectory and y_prec_trajectory which connects y with y_target and 0 with y_prec
void TaskVariable::setTrajectory(uint T,real funnelsdv,real funnelvsdv){
  OPS;
  targetType=trajectoryTT;
  active=true;
  uint t;
  real a;
  y_trajectory.resize(T,y.N);
  y_prec_trajectory.resize(T);
  v_trajectory.resize(T,y.N);
  v_prec_trajectory.resize(T);
  for(t=0;t<T;t++){
    a = (real)t/(T-1);
    y_trajectory[t]()  = ((real)1.-a)*y + a*y_target;
    y_prec_trajectory(t) = (real)1./MT::sqr(sqrt((real)1./y_prec) + ((real)1.-a)*funnelsdv);

    v_trajectory[t]()  = ((real)1.-a)*v + a*v_target;
    v_prec_trajectory(t) = (real)1./MT::sqr(sqrt((real)1./v_prec) + ((real)1.-a)*funnelvsdv);
  }
}

//compute an y_trajectory and y_prec_trajectory which connects y with y_target and 0 with y_prec
void TaskVariable::setConstantTargetTrajectory(uint T){
  OPS;
  targetType=trajectoryTT;
  active=true;
  uint t;
  y_trajectory.resize(T,y.N);
  v_trajectory.resize(T,y.N);
  for(t=0;t<T;t++){
    y_trajectory[t]()  = y_target;
    v_trajectory[t]()  = v_target;
  }
}

//compute an y_trajectory and y_prec_trajectory which connects y with y_target and 0 with y_prec
void TaskVariable::setInterpolatedTargetTrajectory(uint T){
  OPS;
  targetType=trajectoryTT;
  active=true;
  uint t;
  real a;
  y_trajectory.resize(T,y.N);
  v_trajectory.resize(T,y.N);
  for(t=0;t<T;t++){
    a = (real)t/(T-1);
    y_trajectory[t]()  = ((real)1.-a)*y + a*y_target;
    v_trajectory[t]()  = ((real)1.-a)*v + a*v_target;
  }
}

void TaskVariable::setInterpolatedTargetsEndPrecisions(uint T,real inter_y_prec,real end_y_prec,real inter_v_prec,real end_v_prec){
  targetType=trajectoryTT;
  active=true;
  uint t;
  real a;
  y_trajectory.resize(T+1,y.N);   y_prec_trajectory.resize(T+1);
  v_trajectory.resize(T+1,y.N);   v_prec_trajectory.resize(T+1);
  for(t=0;t<=T;t++){
    a = (real)t/T;
    y_trajectory[t]()  = ((real)1.-a)*y + a*y_target;
    v_trajectory[t]()  = ((real)1.-a)*v + a*v_target;
  }
  for(t=0;t<T;t++){
     y_prec_trajectory(t) = inter_y_prec;
     v_prec_trajectory(t) = inter_v_prec;
  }
  y_prec_trajectory(T) = end_y_prec;
  v_prec_trajectory(T) = end_v_prec;
}

void TaskVariable::setInterpolatedTargetsConstPrecisions(uint T,real y_prec,real v_prec){
  targetType=trajectoryTT;
  active=true;
  uint t;
  real a;
  y_trajectory.resize(T+1,y.N);    y_prec_trajectory.resize(T+1);
  v_trajectory.resize(T+1,y.N);   v_prec_trajectory.resize(T+1);
  for(t=0;t<=T;t++){
    a = (real)t/T;
    y_trajectory[t]()  = ((real)1.-a)*y + a*y_target;
    v_trajectory[t]()  = ((real)1.-a)*v + a*v_target;
  }
  for(t=0;t<=T;t++){
    y_prec_trajectory(t) = y_prec;
    v_prec_trajectory(t) = v_prec;
  }
}

//compute an y_trajectory and y_prec_trajectory which connects y with y_target and 0 with y_prec
void TaskVariable::setPrecisionTrajectoryFinal(uint T,real intermediate_prec,real final_prec){
  OPS;
  active=true;
  uint t;
  y_prec_trajectory.resize(T);
  for(t=0;t<T-1;t++) y_prec_trajectory(t) = intermediate_prec;
  y_prec_trajectory(T-1) = final_prec;
}

//compute an y_trajectory and y_prec_trajectory which connects y with y_target and 0 with y_prec
void TaskVariable::setPrecisionTrajectoryConstant(uint T,real const_prec){
  OPS;
  active=true;
  y_prec_trajectory.resize(T);
  y_prec_trajectory = const_prec;
}

void TaskVariable::setPrecisionVTrajectoryFinal(uint T,real intermediate_v_prec,real final_v_prec){
  OPS;
  active=true;
  uint t;
  v_prec_trajectory.resize(T);
  for(t=0;t<T-1;t++) v_prec_trajectory(t) = intermediate_v_prec;
  v_prec_trajectory(T-1) = final_v_prec;
}

//compute an y_trajectory and y_prec_trajectory which connects y with y_target and 0 with y_prec
void TaskVariable::setPrecisionVTrajectoryConstant(uint T,real const_prec){
  OPS;
  active=true;
  v_prec_trajectory.resize(T);
  v_prec_trajectory = const_prec;
}

void TaskVariable::updateState(){
  arr p;
  arr q,qv;
  ors::Vector pi,pj,c;
  arr zi,zj;
  ors::Frame f,fi,fj;

  v_old=v;
  y_old=y;
  
  //get state
  switch(type){
  case posTVT:
    if(j==-1){ ors->kinematics(y,i,&irel); break; }
    pi = ors->bodies(i)->X.p + ors->bodies(i)->X.r * irel.p;
    pj = ors->bodies(j)->X.p + ors->bodies(j)->X.r * jrel.p;
    c = ors->bodies(j)->X.r / (pi-pj);
    y.resize(3); y.setCarray(c.v,3);
    break;
  case zoriTVT:
    if(j==-1){ ors->kinematicsZ(y,i,&irel); break; }
    //relative
    MT_MSG("warning - don't have a correct Jacobian for this TVType yet");
    fi = ors->bodies(i)->X; fi.addRelativeFrame(irel);
    fj = ors->bodies(j)->X; fj.addRelativeFrame(jrel);
    f.setDifference(fi,fj);
    f.r.getZ(c);
    y.setCarray(c.v,3);
    break;
  case rotTVT:       y.resize(3); y.setZero(); break; //the _STATE_ of rot is always zero... the Jacobian not... (hack)
  case contactTVT:   ors->getPenetrationState(p); y.resize(1);  y(0) = p(i);  break;
  case gripTVT:      ors->getGripState(y,i);      break;
  case qLinearTVT:   ors->getJointState(q,qv);    y = params * q;   break;
  case qSquaredTVT:  ors->getJointState(q,qv);    y.resize(1);  y(0) = scalarProduct(params,q,q);  break;
  case qSingleTVT:   ors->getJointState(q,qv);    y.resize(1);  y(0)=q(-i);  break;
  case qLimitsTVT:   ors->getLimitsMeasure(y,params);    break;
  case comTVT:       ors->getCenterOfMass(y);     y.resizeCopy(2);  break;
  case collTVT:      ors->getContactMeasure(y,params(0));   break;
  case colConTVT:    ors->getContactConstraints(y);  break;
  case skinTVT:
    y.resize(params.N);
    y.setZero();
    break;
  case zalignTVT:
    ors->kinematicsZ(zi,i,&irel);
    if(j==-1){
      ors::Vector world_z;
      if(params.N==3) world_z.set(params.p); else world_z=VEC_z;
      zj.setCarray((jrel*world_z).v,3);
    }
    else ors->kinematicsZ(zj,j,&jrel);
    y.resize(1);
    y(0) = scalarProduct(zi,zj);
    break;
  default:  HALT("no such TVT");
  }

  if(y_old.N!=y.N){
    y_old=y;
    v.resizeAs(y); v.setZero();
    v_old=v;
  }

  //v = .5*v + .5*(y - y_old);
  v = y - y_old;

  if(y_target.N==y.N){
    err=norm(y - y_target);
    derr=err - norm(y_old - y_target);
    if(err < state_tol && fabs(derr) < state_tol){
      state = 1;
    }else{
      state = 0;
    }
  }
}
  
void TaskVariable::updateJacobian(){
  arr q,qv;
  ors::Vector normal,d,vi,vj,r,jk,pi,pj;
  arr zi,zj,Ji,Jj,JRj;
  uint k;

  switch(type){
  case posTVT:
    if(j==-1){ ors->jacobian(J,i,&irel); break; }
    pi = ors->bodies(i)->X.p + ors->bodies(i)->X.r * irel.p;
    pj = ors->bodies(j)->X.p + ors->bodies(j)->X.r * jrel.p;
    ors->jacobian(Ji,i,&irel);
    ors->jacobian(Jj,j,&jrel);
    ors->jacobianR(JRj,j);
    J.resize(3,Jj.d1);
    for(k=0;k<Jj.d1;k++){
      vi.set(Ji (0,k),Ji (1,k),Ji (2,k));
      vj.set(Jj (0,k),Jj (1,k),Jj (2,k));
      r .set(JRj(0,k),JRj(1,k),JRj(2,k));
      jk =  ors->bodies(j)->X.r / (vi - vj);
      jk += ors->bodies(j)->X.r / ((-r) ^ (pi - pj));
      J(0,k)=jk(0); J(1,k)=jk(1); J(2,k)=jk(2); 
    }
    break;
  case zoriTVT:  ors->jacobianZ(J,i,&irel);   break;
  case rotTVT:   ors->jacobianR(J,i);   break;
  case contactTVT:  NIY;  break;
  case gripTVT:  NIY;  break;
  case qLinearTVT:  J = params;   break;
  case qSquaredTVT:
    ors->getJointState(q,qv);
    J = params * q;
    J *= (real)2.;
    J.reshape(1,q.N);
    break;
  case qSingleTVT:
    J.resize(1,ors->getJointStateDimension());
    J.setZero();
    J(0,-i) = 1.;
    break;
  case qLimitsTVT:  ors->getLimitsGradient(J,params);  break;
  case comTVT:      ors->getComGradient(J);  J.resizeCopy(2,J.d1);  break;
  case collTVT:     ors->getContactGradient(J,params(0));  break;
  case colConTVT:   ors->getContactConstraintsGradient(J);  break;
  case skinTVT:
    J.clear();
    uint k,l;
    for(k=0;k<params.N;k++){
      l=(uint)params(k);
      ors->jacobian(Ji,l,NULL);
      ors->bodies(l)->X.r.getY(vi);
      vi *= -1.;
      zi.setCarray(vi.v,3);
      J.append(zi*Ji);
    }
    J.reshape(params.N,J.N/params.N);
    break;
  case zalignTVT:
    ors->kinematicsZ(zi,i,&irel);
    ors->jacobianZ(Ji,i,&irel);
    if(j==-1){
      ors::Vector world_z;
      if(params.N==3) world_z.set(params.p); else world_z=VEC_z;
      zj.setCarray((jrel*world_z).v,3);
      Jj.resizeAs(Ji);
      Jj.setZero();
    }else{
      ors->kinematicsZ(zj,j,&jrel);
      ors->jacobianZ(Jj,j,&jrel);
    }
    J = ~zj * Ji + ~zi * Jj;
    J.reshape(1,ors->jd);
    break;
  default:  NIY;
  }
  transpose(tJ,J);
}

void TaskVariable::getHessian(arr& H){
  switch(type){
  case posTVT:
    if(j==-1){ ors->hessian(H,i,&irel); break; }
  default:  NIY;
  }
}

void TaskVariable::updateChange(int t){
  CHECK(y.N,"variable needs to be updated before!");
  arr yt,vt;
  if(t!=-1){
    yt.referToSubDim(y_trajectory,t);
    vt.referToSubDim(v_trajectory,t);
    //y_prec     = y_prec_trajectory(t);
  }else{
    yt.referTo(y_target);
    vt.referTo(v_target);
  }
  CHECK(yt.N==y.N,"targets need to be set before");
  CHECK(vt.N==v.N,"targets need to be set before");
  switch(targetType){
  case directTT:{
    y_change = yt - y;
    break;
  }
  case gainsTT:{
    y_change = a*(yt - y) + b*(vt - v);
    break;
  }
  case attractorTT:{
    y_change = a*(yt - y) + b*(vt - v);
    y_change = v + y_change;
    break;
  }
  case trajectoryTT:{
    y_change = yt - y;
    break;
  }
  default:
    HALT("needs a target type! set targets before!");
  }
}

/*  switch(type){
  case zoriTVT:
    dx.resize(3);
    normal.set(y.p);
    d.set(y_change.p);
    d.makeNormal(normal);
    dx.setCarray(d.v,3);
    break;
  default:
    dx = y_change;
  default:
    break;
  }
    */

void TaskVariable::write(ostream &os) const{
  os <<"CV `" <<name;
  switch(type){
  case posTVT:     os <<"  (pos " <<ors->bodies(i)->name <<")"; break;
  //case relPosTVT:  os <<"  (relPos " <<ors->bodies(i)->name <<'-' <<ors->bodies(j)->name <<")"; break;
  case zoriTVT:    os <<"  (zori " <<ors->bodies(i)->name <<")"; break;
  case rotTVT:     os <<"  (rot " <<ors->bodies(i)->name <<")"; break;
  case contactTVT: os <<"  (contact " <<ors->bodies(i)->name <<' '<<params(0) <<")"; break;
  case gripTVT:    os <<"  (grip " <<ors->bodies(i)->name <<")"; break;
  case qLinearTVT: os <<"  (qLinear " <<sum(params) <<")"; break;
  case qSquaredTVT:os <<"  (qSquared " <<sum(params) <<")"; break;
  case qSingleTVT: os <<"  (qSingle " <<ors->joints(-i)->from->name <<'-' <<ors->joints(-i)->to->name <<")"; break;
  case qLimitsTVT: os <<"  (qLimitsTVT " <<sum(params) <<")"; break;
  case comTVT:     os <<"  (COM)"; break;
  case collTVT:    os <<"  (COLL)"; break;
  case colConTVT:  os <<"  (colCon)"; break;
  case zalignTVT:  os <<"  (zalign " <<ors->bodies(i)->name <<'-' <<ors->bodies(j)->name <<")"; break;
  default: HALT("CV::write - no such TVT");
  }
  os
    <<"\n  y=" <<y
    <<"\n  v=" <<v
    <<"\n  y_target=" <<y_target
    <<"\n  v_target=" <<v_target
    <<"\n  change=" <<y_change
    <<"\n  y_prec=" <<y_prec
    <<"\n  a=" <<a <<"  b=" <<b
    <<"\n  state=" <<state
    <<endl;
}

//===========================================================================
//
// TaskVariableList functions
//






void reportAll(TaskVariableList& CS,ostream& os,bool onlyActives){
  for(uint i=0;i<CS.N;i++) if(!onlyActives || CS(i)->active){
    os <<'[' <<i <<"] " <<*CS(i);
  }
}

void reportNames (TaskVariableList& CS,ostream& os,bool onlyActives){
  uint i,j,n=1;
  os <<"CVnames = {";
  for(i=0;i<CS.N;i++) if(!onlyActives || CS(i)->active){
    for(j=0;j<CS(i)->y.N;j++){
      os <<"'" <<n <<'-' <<CS(i)->name <<j <<"' ";
      n++;
    }
  }
  os <<"};" <<endl;
}

void reportState (TaskVariableList& CS,ostream& os,bool onlyActives){;
  uint i;
  MT::IOraw=true;
  for(i=0;i<CS.N;i++) if(!onlyActives || CS(i)->active){
    os <<CS(i)->y;
  }
  os <<endl;
}

void reportErrors(TaskVariableList& CS,ostream& os,bool onlyActives,int t){
  uint i;
  real e,E=0.;
  for(i=0;i<CS.N;i++) if(!onlyActives || CS(i)->active){
    if(t!=-1)
      if(t) e=norm(CS(i)->y - CS(i)->y_trajectory[t-1]);
      else  e=0.;
    else
      e=norm(CS(i)->y - CS(i)->y_target);
    os <<e <<' ';
    E += e; //*CS(i)->y_prec;
  }
  os <<E <<endl;
}

void activateAll(TaskVariableList& CS,bool active){
  for(uint i=0;i<CS.N;i++) CS(i)->active=active;
}

void updateState(TaskVariableList& CS){
  for(uint i=0;i<CS.N;i++){
    CS(i)->updateState();
  }
}

void updateJacobian(TaskVariableList& CS){
  for(uint i=0;i<CS.N;i++){
    CS(i)->updateJacobian();
  }
}

void updateChanges(TaskVariableList& CS,int t){
  for(uint i=0;i<CS.N;i++) if(CS(i)->active){
    CS(i)->updateChange(t);
  }
}

void getJointJacobian(TaskVariableList& CS,arr& J){
  uint i,n=0;
  J.clear();
  for(i=0;i<CS.N;i++) if(CS(i)->active){
    CS(i)->updateJacobian();
    J.append(CS(i)->J);
    n=CS(i)->J.d1;
  }
  J.reshape(J.N/n,n);
}

void getJointXchange(TaskVariableList& CS,arr& y_change){
  uint i;
  y_change.clear();
  for(i=0;i<CS.N;i++) if(CS(i)->active){
    y_change.append(CS(i)->y_change);
  }
  y_change.reshape(y_change.N);
}

real getCost_obsolete(TaskVariableList& CS,const arr& W,int t){
  uint i;
  real e,C=0.;
  for(i=0;i<CS.N;i++) if(CS(i)->active){
    if(t!=-1){
      e=sumOfSqr(CS(i)->y - CS(i)->y_trajectory[t]);
    }else{
      e=sumOfSqr(CS(i)->y - CS(i)->y_target);
    }
    C += e*CS(i)->y_prec;
    //cout <<"cost(" <<CS(i)->name <<") = " <<e <<"," <<e*CS(i)->y_prec <<endl;
  }
  return C;
}

void getCostGradient_obsolete(TaskVariableList& CS,arr& dCdq,const arr& W,int t){
  uint i,n=W.d0;
  dCdq.resize(n);
  dCdq.setZero();
  arr e,J,tJ,dx;
  for(i=0;i<CS.N;i++) if(CS(i)->active){
    CS(i)->updateJacobian();
    if(t!=-1){
      e=CS(i)->y - CS(i)->y_trajectory[t];
    }else{
      e=CS(i)->y - CS(i)->y_target;
    }
    dCdq += CS(i)->tJ * (e*((real)2.*CS(i)->y_prec));
  }
}

void hierarchicalControl_obsolete(TaskVariableList& CS,arr& dq,const arr& W){
  uint i,n=W.d0;
  dq.resize(n);
  dq.setZero();
  arr hatJ,invHatJ,N;
  N.setId(n);
  arr invW;
  inverse_SymPosDef(invW,W);
  for(i=0;i<CS.N;i++) if(CS(i)->active){
    CS(i)->updateJacobian();

    hatJ = CS(i)->J * N;
    pseudoInverse(invHatJ,hatJ,invW,1e-5);
    dq += invHatJ * (CS(i)->y_change - CS(i)->J * dq);
    N  -= invHatJ * hatJ;
  }
}

void bayesianControl_obsolete(TaskVariableList& CS,arr& dq,const arr& W){
  uint n=W.d0;
  dq.resize(n);
  dq.setZero();
  uint i;
  arr a(n),A(n,n),invA(n,n);
  //arr Q,JQ;
  A=W;
  a.setZero();
  arr w(3);
  for(i=0;i<CS.N;i++) if(CS(i)->active){
    CS(i)->updateJacobian();
    a += CS(i)->y_prec * CS(i)->tJ * CS(i)->y_change;
    A += CS(i)->y_prec * CS(i)->tJ * CS(i)->J;
  }
  inverse_SymPosDef(invA,A);
  dq = invA * a;
}

/*void bayesianIterateControl_obsolete(TaskVariableList& CS,
                            arr& qt,const arr& qt_1,const arr& W,real eps,uint maxIter){
  uint j;
  qt=qt_1;
  arr dq;
  for(j=0;j<maxIter;j++){
    top.setq(qt);
    bayesianIKControl(top,dq,W);
    if(j<3) qt+=dq;
    //else if(j<10) qt+=.8*dq;
    else qt+=.8*dq;
    if(dq.absMax()<eps) break;
  }
  if(j==maxIter) HALT("warning: IK didn't converge (|last step|="<<dq.absMax()<<")");
  else cout <<"IK converged after steps=" <<j <<endl;
}*/

void additiveControl_obsolete(TaskVariableList& CS,arr& dq,const arr& W){
  dq.resize(W.d0);
  dq.setZero();
  uint i,n=0;
  arr invJ;
  arr invW;
  inverse_SymPosDef(invW,W);
  for(i=CS.N;i--;) if(CS(i)->active){
    CS(i)->updateJacobian();
    pseudoInverse(invJ,CS(i)->J,invW,0.);
    dq += invJ * CS(i)->y_change;
    n++;
  }
  dq/=(real)n;
}

/*OLD
void bayesianPlanner_obsolete(ors::Graph *ors,TaskVariableList& CS,SwiftModule *swift,OpenGL *gl,
                     arr& q,uint T,const arr& W,uint iterations,
                     std::ostream* os,int display,bool repeat){
  //FOR THE OLD VERSION, SEE SMAC.CPP IN THE DEPOSIT
  uint n=W.d0,i;
  arr J,tJ,phiHatQ,dx;
  arr q0,qv0,invW;
  inverse_SymPosDef(invW,W);

  byteA img(300,500,3);

  ors->getJointState(q0,qv0);

  arr tmp,
    a(T,n),invA(T,n,n),
    z(T,n),invZ(T,n,n),
    b(T,n),B(T,n,n),invB(T,n,n),
    r(T,n),R(T,n,n),
    hatq(T,n);
  a[0]=q0;
  invA[0].setDiag(1e10);
  b[0]=q0;
  B[0].setDiag(1e-10);
  invB[0].setDiag(1e10);
  r[0]=0.;
  R[0]=0.;
  z.setZero();
  invZ.setZero();

  MT::timerStart();
  
  uint k,t,dt,t0;
  for(k=0;k<iterations;k++){
    if(!(k&1)){ dt=1; t0=1; }else{ dt=(uint)-1; t0=T-1; }
    for(t=t0;t<T && t>0;t+=dt){
      //compute (a,A)
      inverse_SymPosDef(tmp,invA[t-1] + R[t-1]);
      a[t] = tmp * (invA[t-1]*a[t-1] + r[t-1]);
      inverse_SymPosDef(invA[t](), invW + tmp);
      
      //cout <<"a\n" <<a[t] <<endl <<invA[t] <<endl;

      //compute (z,Z)
      if(k && t<T-1){
        inverse_SymPosDef(tmp,invZ[t+1] + R[t+1]);
        z[t] = tmp * (invZ[t+1]*z[t+1] + r[t+1]);
        inverse_SymPosDef(invZ[t](), invW + tmp);
      }
      if(k && t==T-1){
        z[t] = b[t];
        invZ[t].setDiag(1e10); //fixes the end posture!, use 1e-5 otherwise
      }

      //cout <<"z\n" <<z[t] <<endl <<invZ[t] <<endl;
      
      //compute (r,R)
      //if(k) hatq[t]()=.2*b[t]+.8*hatq[t]; else hatq[t]()=a[t];
      if(k) hatq[t]()=b[t]; else hatq[t]()=a[t];
      ors->setJointState(hatq[t]);
      ors->calcNodeFramesFromEdges();
      computeProxiesUsingSwift(*ors,*swift,false);
      //slGetProxies(*ors,*ode);
      r[t].setZero();
      R[t].setZero();
      for(i=0;i<CS.N;i++) if(CS(i)->active){
        CS(i)->updateState();    phiHatQ=CS(i)->y;
        CS(i)->updateJacobian(); J=CS(i)->J; tJ=CS(i)->tJ;
        dx = CS(i)->y_trajectory[t] - phiHatQ;
        r[t]() += CS(i)->y_prec_trajectory(t) * tJ * dx;
        R[t]() += CS(i)->y_prec_trajectory(t) * tJ * J;
      }
      r[t]() += R[t] * hatq[t];

      //cout <<"r\n" <<r[t] <<endl <<R[t] <<endl;

      //compute (b,B);
      invB[t] = invA[t] + invZ[t] + R[t];
      //cout <<"invB\n" <<invB[t] <<endl;
      inverse_SymPosDef(B[t](), invB[t]);
      b[t] = B[t] * (invA[t]*a[t] + invZ[t]*z[t] + r[t]);

      //cout <<"b\n" <<b[t] <<endl <<B[t] <<endl;

      //display
      if(display>0){
        ors->setJointState(b[t]);
        ors->calcNodeFramesFromEdges();
        //if(t==1 || !(t%display)){ gl->text.clr() <<k <<':' <<t; gl->update(); }
        //glGrabImage(img); write_ppm(img,STRING("imgs/plan_"<<std::setfill('0') <<std::setw(3) <<k<<std::setfill('0') <<std::setw(3) <<((k&1)?T-t:t)<<".ppm"),true);
      }

      if(repeat){
        //meassure offset
        real off=sqrDistance(W,b[t],hatq[t]);
        //cout <<"off = " <<off <<endl;
        if(false && k>0 && off>.05){
          //cout <<t <<" REPEAT: off = " <<off <<endl;
          t-=dt;
        }
      }
    }


    //evaluate trajectory
    //cout <<"variances over time = ";
    //for(t=0;t<T;t++) cout <<' ' <<trace(B[t]);
    real cost_t,cost1=.0,cost2=.0,length=0.;
    for(t=0;t<T;t++){
      ors->setJointState(b[t]);
      ors->calcNodeFramesFromEdges();
      computeProxiesUsingSwift(*ors,*swift,false);
      //slGetProxies(*ors,*ode);
      updateState(CS);
      if(t>0) cost2 += sqrDistance(W,b[t-1],b[t]);
      if(t>0) length += metricDistance(W,b[t-1],b[t]);
      cost1 += cost_t = getCost(CS,W,t);  //cout <<"cost = " <<cost_t <<endl;
    }
    *os <<std::setw(3) <<k
        <<"  time " <<MT::timerRead(false)
        <<"  cost1 " <<cost1
        <<"  cost2 " <<cost2
        <<"  length " <<length
        <<"  total-cost " <<cost1+cost2 <<endl;
  }

  q = b;
}
*/

#if 0
void SMAC::readCVdef(std::istream& is){
  char c;
  TaskVariable *cv;
  MT::String name,ref1,ref2;
  ors::Frame f;
  uint i,j,k;
  arr mat;
  MT::String::readSkipSymbols=" \n\r\t";
  MT::String::readStopSymbols=" \n\r\t";
  for(;;){
    MT::skip(is);
    is.get(c);
    if(!is.good()) break;
    cv=&CVs.append();
    switch(c){
    case 'p':
      is >>name >>ref1 >>"<" >>f >>">";
      cv->initPos(name,*ors,ors->getName(ref1)->index,f);
      break;
    case 'j':
      is >>name >>ref1;
      cv->initQSingle(name,*ors,ors->getName(ref1)->firstIn->index);
      break;
    case 'l':
      is >>name >>ref1 >>ref2;
      i=ors->getName(ref1)->firstIn->index;
      j=ors->getName(ref2)->firstIn->index;
      mat.resize(j-i+1,ors->getJointStateDimension());
      mat.setZero();
      for(k=0;k<=j-i;k++) mat(k,i+k)=1.;
      cv->initQLinear(name,*ors,mat);
      break;
    case 's':
      is >>name >>ref1 >>ref2;
      i=ors->getName(ref1)->firstIn->index;
      j=ors->getName(ref2)->firstIn->index;
      mat.resize(j-i+1,ors->getJointStateDimension());
      mat.setZero();
      for(k=0;k<=j-i;k++) mat(k,i+k)=1.;
      cv->initQSquared(name,*ors,mat);
      break;
    case 'c':
      is >>name >>ref1;
      cv->initContact(name,*ors,ors->getName(ref1)->index);
      break;
    case 'X':
      is >>name >>ref1;
      cv->initGrip(name,*ors,ors->getName(ref1)->index);
      break;
    case 'o':
      is >>name >>ref1 >>ref2 >>"<" >>f >>">";
      if(ref2=="*")
        cv->initOri(name,*ors,ors->getName(ref1)->index,-1,f);
      else
	cv->initOri(name,*ors,ors->getName(ref1)->index,ors->getName(ref2)->index,f);
      break;
    case 'M':
      is >>name;
      cv->initCom(name,*ors);
      break;
    case 'C':
      is >>name;
      cv->initCollision(name,*ors);
      break;
    default:
      NIY;
    }
    if(is.fail()) HALT("error reading `"<<c<<"' variable in smac");
  }
  is.clear();
}

void SMAC::CVclear(){
  uint i;
  for(i=0;i<CS.N;i++)
    CS(i)->y.clear();
}

void SMAC::plotCVs(){
  plotData.points.resize(2*CS.N);
  uint i;
  for(i=0;i<CS.N;i++){
    plotData.points(2*i) = CS(i)->y;
    plotData.points(2*i+1) = CS(i)->y_target;
  }
}

#endif
