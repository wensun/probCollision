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
#include "plot.h"

uint countMsg=0,countSetq=0;

//===========================================================================
//
// general documentation
//

/*! \brief libSOC -- Stochastic Optimal Control library

    This is the core namespace of libSOC.  See the <a
    href="../guide.pdf">guide</a> for an introduction.

    Please see also the header <a href="soc_8h-source.html">MT/soc.h</a> . */
namespace soc{};


//===========================================================================
//
// trivial helpers
//

//! \brief get the velocity vt of a trajectory q at time t
void soc::getVelocity(arr& vt,const arr& q,uint t,real tau){
  if(!t) vt = (q[0]-q[0])/tau;
  else   vt = (q[t]-q[t-1])/tau;
}

  //! compute the full (q,v) trajectory from a trajectory q
void soc::getPhaseTrajectory(arr& _q,const arr& q,real tau){
  uint T=q.d0,n=q.d1,t;
  arr vt;
  _q.resize(T,2,n);
  for(t=0;t<T;t++){
    getVelocity(vt,q,t,tau);
    _q.subDim(t,0)=q[t];
    if(t) _q.subDim(t,1)=(q[t]-q[t-1])/tau;
    else  _q.subDim(t,1)=0.;
  }
  _q.reshape(T,2*n);
}

//! simply get the q-trajectory from a (q,v)-trajectory
void soc::getPositionTrajectory(arr& q,const arr& _q){
  uint T=_q.d0,n=_q.d1/2,i,t;
  CHECK(2*n==_q.d1,"")
  q.resize(T,n);
  for(t=0;t<T;t++) for(i=0;i<n;i++) q(t,i)=_q(t,i);
}

void soc::interpolateTrajectory(arr& qNew,const arr& q,real step){
  uint t,T=floor(q.d0/step);
  qNew.resize(T,q.d1);
  real tref=0.,mod;
  for(t=0;t<T;t++){
    uint a=floor(tref),b=ceil(tref);
    mod = tref - a;
    qNew[t] = ((real)1.-mod)*q[a] + mod*q[b];
    tref+=step;
  }
}

/*! \brief use regularized Inverse Kinematics to compute a joint
    trajectory from a given task trajectory x for the 0-th task variable */
void soc::getJointFromTaskTrajectory(SocSystemAbstraction& soci,arr& q,const arr& x){
  uint t,T=x.d0,n=soci.qDim();
  arr phiq,J,tJ,invJ,W,invW;
  q.resize(T,n);
  soci.getq0(q[0]());
  soci.getW(W);
  inverse_SymPosDef(invW,W);
  for(t=1;t<T;t++){
    soci.setq(q[t-1]);
    soci.getPhi(phiq,0);
    soci.getJtJ(J,tJ,0);
    pseudoInverse(invJ,J,invW,1e-5);
    q[t]() = q[t-1] + invJ*(x[t]-phiq);
  }
}

/*! \brief use regularized Inverse Kinematics to compute a joint
    trajectory from the task trajectory previously specifies for the
    taskid-th task variable */
void soc::straightTaskTrajectory(SocSystemAbstraction& soci,arr& q,uint taskid){
  uint t,T=soci.nTime(),n=soci.qDim();
  arr phiq,xt,J,tJ,invJ,W,invW;
  real prec;
  q.resize(T+1,n);
  soci.getq0(q[0]());
  soci.getW(W);
  inverse_SymPosDef(invW,W);
  for(t=1;t<=T;t++){
    soci.setq(q[t-1]);
    soci.getPhi(phiq,taskid);
    soci.getJtJ(J,tJ,taskid);
    soci.getTarget(xt,prec,taskid,t);
    pseudoInverse(invJ,J,invW,1e-5);
    q[t]() = q[t-1] + invJ*(xt-phiq);
  }
  
}

//! not-yet-implemented
void soc::partialJointFromTaskTrajectory(SocSystemAbstraction& soci,arr& dCdx,const arr& delCdelq,const arr& q,const arr& x){
  NIY;
}

//===========================================================================
//
// SocSystemAbstraction defines routines that optimizers want access to
//

soc::SocSystemAbstraction::SocSystemAbstraction(){
  gl=NULL;
  os=NULL;
  scalePower=0;
}

soc::SocSystemAbstraction::~SocSystemAbstraction(){
}

uint soc::SocSystemAbstraction::uDim(){ NIY; }
void soc::SocSystemAbstraction::getqv0(arr& q_){ NIY; }
void soc::SocSystemAbstraction::getqv0(arr& q,arr& qd){ NIY; }
real soc::SocSystemAbstraction::getTau(bool scaled){ NIY; }
void soc::SocSystemAbstraction::setqv (const arr& q_){ NIY; }
void soc::SocSystemAbstraction::setqv (const arr& q,const arr& qd){ NIY; }
void soc::SocSystemAbstraction::getH  (arr& H){ NIY; }
void soc::SocSystemAbstraction::getQ  (arr& Q){ NIY; }
void soc::SocSystemAbstraction::getJqd(arr& Jqd_i,uint i){ NIY; }
void soc::SocSystemAbstraction::getHessian(arr& H_i,uint i){ NIY; }
void soc::SocSystemAbstraction::getTargetV(arr& v_i,real& prec,uint i,uint t){ NIY; }
//void soc::SocSystemAbstraction::getLinearConstraint(arr& c,real& coff,uint i,uint t){ NIY; }
bool soc::SocSystemAbstraction::isConstrained(uint i,uint t){ NIY; return false; }
void soc::SocSystemAbstraction::getMF(arr& M,arr& F){ NIY; }
void soc::SocSystemAbstraction::getMinvF(arr& Minv,arr& F){ NIY; }

void soc::SocSystemAbstraction::getProcess(arr& A,arr& a,arr& B){
  uint n=qDim();
  if(!dynamic){
    A.setId(n);
    B.setId(n);
    a.resize(n);
    a.setZero();
  }else{
    real tau=getTau(false);
    arr I,Z,Minv,F;
    I.setId(n);
    Z.resize(n,n); Z.setZero();
    
    getMinvF(Minv,F);
    
    A.setBlockMatrix(I,tau*I,Z,I);
    //real alpha = .1; //with fricion
    //A.setBlockMatrix(I,tau*I-tau*alpha*Minv,Z,I-tau*alpha*Minv);
    
    B.resize(2*n,n);
    B.setZero();
    B.setMatrixBlock(tau*tau*Minv,0,0);
    B.setMatrixBlock(tau*Minv,n,0);

    a.resize(2*n);
    a.setZero();
    a.setVectorBlock(tau*tau*Minv*F,0);
    a.setVectorBlock(tau*Minv*F,n);
  }
  for(uint i=1;i<scalePower;i++){
    a = A*a + a;
    B = A*B + B;
    A = A*A;
  }
}

void soc::SocSystemAbstraction::getProcess(arr& A,arr& tA,arr& invA,arr& invtA,arr& a,arr& B,arr& tB){
  getProcess(A,a,B);
  //A^{-1} is A transpose and the lower-left matrix negative.. BLOCKMATRIX(Id, -2^scale*tau*Id, 0, Id)
  invA=A;
  uint n=qDim();
  for(uint i=0;i<n;i++) invA(i,n+i) *= -1.;
  transpose(tA,A);
  transpose(tB,B);
  transpose(invtA,invA);
}

real soc::SocSystemAbstraction::getCosts(arr& R,arr& r,uint t,const arr& qt){
  uint i,m=nTasks(),n=qDim();
  real C=0.;
  if(!dynamic){ //kinematic
    arr phiHatQ,J,tJ,x;
    real prec;
    R.resize(n,n); R.setZero();
    r.resize(n);   r.setZero();
    for(i=0;i<m;i++) if(isConditioned(i,t<<scalePower)){
      getPhi      (phiHatQ,i);
      getTarget   (x,prec,i,t<<scalePower);
      C += prec*sqrDistance(x,phiHatQ);
      getJtJ      (J,tJ,i);
      R += prec*tJ*J;
      r += prec*tJ*(x - phiHatQ + J*qt);
    }
  }else{
    uint n2=2*n;
    arr phiHatQ,Jqd,J,tJ,x,v,Ri,ri;
    real prec,precv;
    R.resize(n2,n2); R.setZero();
    r.resize(n2);    r.setZero();
    Ri.resize(n2,n2); Ri.setZero();
    ri.resize(n2);    ri.setZero();
    for(i=0;i<m;i++) if(isConditioned(i,t<<scalePower)){
      getPhi       (phiHatQ,i);
      getTarget    (x,prec,i,t<<scalePower);
      C += prec*sqrDistance(x,phiHatQ);
      getJqd       (Jqd,i);
      getTargetV   (v,precv,i,t<<scalePower);
      C += precv*sqrDistance(v,Jqd);
      getJtJ       (J,tJ,i);
      Ri.setMatrixBlock(prec*tJ*J,0,0);
      Ri.setMatrixBlock(precv*tJ*J,n,n);
      ri.setVectorBlock(prec*tJ*(x - phiHatQ + J*qt),0);
      ri.setVectorBlock(precv*tJ*v,n);
      R += Ri;
      r += ri;
    }
  }
  /*for(uint i=1;i<scalePower;i++){
    r = 2.*r;
    R = 2.*R;
  }*/
  return C;
}

void soc::SocSystemAbstraction::getConstraints(arr& cdir,arr& coff,uint t,const arr& qt){
  uint i,j,m=nTasks(),con=0,n=qDim();
  arr phiHatQ,J,tJ;
  cdir.clear();
  coff.clear();
  for(i=0;i<m;i++) if(isConstrained(i,t)){
    getPhi(phiHatQ,i);
    getJtJ(J,tJ,i);
    for(j=0;j<phiHatQ.N;j++){ //loop through all constraints in the constraint vector
      //if(phiHatQ(j)>.5) continue; //that's good enough -> don't add the constraint
      J *= (real).5;
      coff.append(-phiHatQ(j) + scalarProduct(J[j],qt));
      cdir.append(J[j]);
      con++;
    }
    /*
    CHECK(phiHatQ.N==1,"so far, constraints work only on 1D task variables!");
    if(phiHatQ(0)>0.){ //potential violation, else discard
      getJtJ(J,tJ,i);
      cdir.append(-J);
      *if(phiHatQ(0)>1.){
        cout <<"constraint violated: "<<phiHatQ <<" -> making it more graceful.." <<endl;
        phiHatQ=1.;
      }*
      coff.append(phiHatQ-J*qt-5.);
      //cout <<"qt= " <<qt <<"\nJ*qt="  <<J*qt <<"\nphiHatQ= " <<phiHatQ <<endl;
      con++;
    }
    */
  }
  cdir.reshape(con,n);
  coff.reshape(con);
}

void soc::SocSystemAbstraction::constantTrajectory(arr& q){
  CHECK(!dynamic,"");
  uint t,T=nTime(),n=qDim();
  q.resize(T,n);
  getq0(q[0]());
  for(t=1;t<T;t++) q[t]() = q[0];
}

void soc::SocSystemAbstraction::passiveDynamicTrajectory(arr& q){
  uint t,T=nTime(),n=qDim();
  real tau=getTau();
  arr qd,Minv,F;
  qd.resize(T,n);
  q .resize(T,n);
  getqv0(q[0](),qd[0]());
  for(t=0;t<T-1;t++){
    setqv(q[t],qd[t]);
    getMinvF(Minv,F);
    qd[t+1]() = qd[t] + tau*Minv*F;
    q [t+1]() = q [t] + tau*qd[t+1];
  }
}

void soc::SocSystemAbstraction::controlledDynamicTrajectory(arr& q,const arr& u){
  uint t,T=nTime(),n=qDim();
  if(!dynamic){
    q.resize(T,n);
    getq0(q[0]());
    for(t=1;t<T;t++) q[t] = q[t-1]+u[t-1];
    return;
  }
  real tau=getTau();
  arr qd,Minv,F;
  qd.resize(T,n);
  q .resize(T,n);
  getqv0(q[0](),qd[0]());
  for(t=0;t<T-1;t++){
    setqv(q[t],qd[t]);
    getMinvF(Minv,F);
    qd[t+1]() = qd[t] + tau*Minv*(F+u[t]);
    q [t+1]() = q [t] + tau*qd[t+1];  //note: this is qd[t+1] on the RHS
#if 0
    arr xx,x,A,a,B;
    x.setBlockVector(q[t],qd[t]);
    getProcessDynamic(A,a,B);
    xx = A*x + a + B*u[t];
    cout <<"xx = " <<xx <<"\nq=" <<q[t+1] <<qd[t+1] <<endl;
#endif
  }
}

void soc::SocSystemAbstraction::getControlFromTrajectory(arr& u,const arr& q){
  uint t,T=nTime(),n=qDim();
  if(!dynamic){
    u.resize(T,n);
    u[T-1]().setZero();
    for(t=0;t<T-1;t++) u[t] = q[t+1]-q[t];
    return;
  }
  real tau,tau_1,tau_2;
  tau=getTau();
  tau_1 = 1./tau;
  tau_2 = tau_1*tau_1;
  arr M,F;
  u.resize(T,n);
  u[T-1]().setZero();
  for(t=0;t<T-1;t++){
    if(!t) setq(q[t]);
    else   setqv(q[t],tau_1*(q[t]-q[t-1]));
    getMF(M,F);
    if(t<T-1 && t>0)
      u[t]() = tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t]) - F;
    if(t==0)
      u[t]() = tau_2*M*(q[t+1]-q[t]) - F;
  }
}

real soc::getTaskCost(soc::SocSystemAbstraction& soci,int t,int i){
  real C=0.;
  if(i==-1){ //sum over all tasks
    uint m=soci.nTasks();
    for(i=0;i<(int)m;i++) if(soci.isConditioned(i,t)) C += getTaskCost(soci,t,i);
  }else{     //one specific task
    arr phiHatQ,Jqd,x,v;
    real prec,precv;
    soci.getPhi      (phiHatQ,i);
    soci.getTarget   (x,prec,i,t);
    C += prec*sqrDistance(x,phiHatQ);
    if(soci.dynamic){
      soci.getJqd       (Jqd,i);
      soci.getTargetV   (v,precv,i,t);
      C += precv*sqrDistance(v,Jqd);
    }
  }
  return C;
}

void soc::getTaskCostGradient(soc::SocSystemAbstraction& soci,arr& dCdq,int t){
  uint n=soci.qDim(),m=soci.nTasks();
  uint i;
  dCdq.resize(n);
  dCdq.setZero();
  arr phiHatQ,J,tJ,x,v,Jqd;
  real prec,precv;
  for(i=0;i<m;i++) if(soci.isConditioned(i,t)){
    soci.getPhi      (phiHatQ,i);
    soci.getTarget   (x,prec,i,t);
    soci.getJtJ      (J,tJ,i);
    dCdq += tJ * ((phiHatQ-x)*((real)2.*prec));
    if(soci.dynamic){
      soci.getJqd       (Jqd,i);
      soci.getTargetV   (v,precv,i,t);
      dCdq += tJ * ((v - Jqd)*((real)2.*precv));
    }
  }
}

real soc::SocSystemAbstraction::computeTotalCost(const arr& q,bool plot){
  //if(dynamic){ return getTotalDynamicCost(*this,q,os); }
  uint t,T=nTime();
  CHECK(q.nd==2 && q.d0==T+1 && q.d1==qDim(),"q has wrong dimension: "<<q.getDim());
  real taskCost=0., ctrlCost=0.;
  arr W,H,M,F;
  real tau,tau_1=1.,tau_2=1.;
  if(!dynamic){
    getW(W);
  }else{
    getH(H);
    tau=getTau();
    tau_1 = 1./tau;
    tau_2 = tau_1*tau_1;
  }
  arr taskC(T+1);  taskC.setZero();
  arr ctrlC(T+1);  ctrlC.setZero();
  for(t=0;t<=T;t++){
    countSetq++;
    if(!dynamic || !t) setq(q[t]);
    else setqv(q[t],tau_1*(q[t]-q[t-1]));
    taskCost += taskC(t) = getTaskCost(*this,t);  //cout <<"cost = " <<cost_t <<endl;
    if(!dynamic){
      if(t>0) ctrlC(t) = sqrDistance(W,q[t-1],q[t]);
    }else{
      getMF(M,F);
      if(t<T && t>0)
        ctrlC(t) = sqrDistance(H,tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t]),F);
      if(t==0)
        ctrlC(t) = sqrDistance(H,tau_2*M*(q[t+1]-q[t]),F);
    }
    ctrlCost += ctrlC(t);
  }
#ifdef NIKOLAY
 // plot = false;
#endif
  if(plot){
    std::ofstream fil;
    MT::open(fil,"z.trana");
    for(t=0;t<T;t++){
      fil
          <<"time " <<t
          <<"  ctrlC " <<ctrlC(t)
          <<"  taskC " <<taskC(t)
          <<"  totC "  <<ctrlC(t)+taskC(t)
          <<"  q " <<q[t]
          <<endl;
    }
    gnuplot("plot 'z.trana' us 0:4 title 'ctrl costs', 'z.trana' us 0:6 title 'task costs', 'z.trana' us 0:8 title 'tot costs'");
  }
  
#ifdef NIKOLAY
      if(os) *os
        <<" " <<taskCost
        <<" " <<ctrlCost
        <<" " <<taskCost+ctrlCost <<endl;
#else
       if(os) *os
        <<"  task-cost " <<taskCost
        <<"  control-cost " <<ctrlCost
        <<"  total-cost " <<taskCost+ctrlCost <<endl;
#endif
 
  return taskCost+ctrlCost;
}

void getTotalDynamicCostGradient(soc::SocSystemAbstraction& soci,const arr& q,arr& dCdq);

real soc::SocSystemAbstraction::computeTotalCostGradient(arr& dCdq,const arr& q){
  if(dynamic){ getTotalDynamicCostGradient(*this,q,dCdq); return 0.; }
  uint t,T=nTime();
  CHECK(q.nd==2 && q.d0==T+1 && q.d1==qDim(),"");
  dCdq.resizeAs(q);
  dCdq.setZero();
  real length=.0, cost1=0., cost2=0.;
  real cost_t;
  arr W;
  getW(W);
  for(t=0;t<=T;t++){
    countSetq++;
    setq(q[t]);
    cost1 += cost_t = getTaskCost(*this,t);  //cout <<"cost = " <<cost_t <<endl;
    getTaskCostGradient(*this,dCdq[t](),t);
    if(t>0) cost2 += sqrDistance(W,q[t-1],q[t]);
    if(t>0) length += metricDistance(W,q[t-1],q[t]);
    if(t>0)   dCdq[t]() += (real)2.*W*(q[t]-q[t-1]);
    if(t+1<T) dCdq[t]() -= (real)2.*W*(q[t+1]-q[t]);
  }
  return cost1+cost2;
}

real getTotalDynamicCost_obsolete(soc::SocSystemAbstraction& soci,const arr& q,std::ostream* os){
  uint t,T=soci.nTime(),i,m=soci.nTasks();
  CHECK(q.nd==2 && q.d0==T && q.d1==soci.qDim(),"q has wrong dimension: "<<q.getDim());
  real taskCost=0., controlCost=0.;
  arr H;
  soci.getH(H);
  real tau,tau_1,tau_2;
  tau=soci.getTau();
  tau_1 = 1./tau;
  tau_2 = tau_1*tau_1;
  arr phiHatQ,x,v,J,tJ,M,F;
  real prec,precv;
  for(t=0;t<T;t++){
    if(!t) soci.setq(q[t]);
    else   soci.setqv(q[t],tau_1*(q[t]-q[t-1]));

    for(i=0;i<m;i++) if(soci.isConditioned(i,t)){
      soci.getPhi      (phiHatQ,i);
      soci.getTarget   (x,prec,i,t);
      taskCost += prec * sqrDistance(phiHatQ,x);

      soci.getJtJ      (J,tJ,i);
      soci.getTargetV  (v,precv,i,t);
      if(t>0)
        taskCost += precv * sqrDistance(J*(tau_1*(q[t]-q[t-1])),v);
    }

    soci.getMF(M,F);
    if(t<T-1 && t>0)
      controlCost += sqrDistance(H,tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t]),F);
    if(t==0)
      controlCost += sqrDistance(H,tau_2*M*(q[t+1]-q[t]),F);
  }
#if 1
  if(os) *os
    <<"  task-cost " <<taskCost
    <<"  control-cost " <<controlCost
    <<"  total-cost " <<taskCost+controlCost <<endl;
#else
  if(os) *os
    <<' ' <<taskCost
    <<' ' <<controlCost
    <<' ' <<taskCost+controlCost <<endl;
#endif
  return taskCost+controlCost;
}

void getTotalDynamicCostGradient(soc::SocSystemAbstraction& soci,const arr& q,arr& dCdq){
  uint t,T=soci.nTime(),i,m=soci.nTasks();
  CHECK(q.nd==2 && q.d0==T+1 && q.d1==soci.qDim(),"");
  dCdq.resizeAs(q);
  dCdq.setZero();
  arr H;
  soci.getH(H);
  real tau,tau_1,tau_2;
  tau=soci.getTau();
  tau_1 = 1./tau;
  tau_2 = tau_1*tau_1;
  arr phiHatQ,x,v,J,tJ,M,Mt,F;
  real prec,precv;
  for(t=0;t<=T;t++){
    countSetq++;
    if(!t) soci.setq(q[t]);
    else   soci.setqv(q[t],tau_1*(q[t]-q[t-1]));

    for(i=0;i<m;i++) if(soci.isConditioned(i,t)){
      soci.getPhi      (phiHatQ,i);
      soci.getTarget   (x,prec,i,t);
      soci.getJtJ      (J,tJ,i);
      dCdq[t]() += ((real)2.*prec)*(tJ * (phiHatQ-x));

      soci.getTargetV   (v,precv,i,t);
      if(t>0){
        dCdq[t  ]() += ((real)2.*precv*tau_1)*(tJ * (J*(tau_1*(q[t]-q[t-1])) - v));
        dCdq[t-1]() -= ((real)2.*precv*tau_1)*(tJ * (J*(tau_1*(q[t]-q[t-1])) - v));
      }
    }

    soci.getMF(M,F);
    transpose(Mt,M);
    if(t<T && t>0){
      dCdq[t  ]() -= ((real)4.*tau_2)*Mt*H*(tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t])-F);
      dCdq[t+1]() += ((real)2.*tau_2)*Mt*H*(tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t])-F);
      dCdq[t-1]() += ((real)2.*tau_2)*Mt*H*(tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t])-F);
    }
    if(t==0){
      dCdq[t  ]() -= ((real)2.*tau_2)*Mt*H*(tau_2*M*(q[t+1]-q[t])-F);
      dCdq[t+1]() += ((real)2.*tau_2)*Mt*H*(tau_2*M*(q[t+1]-q[t])-F);
    }
  }
}

//! play the trajectory using OpenGL
void soc::SocSystemAbstraction::displayState(const arr& q,const arr *invQ,const char *text){
  if(gl){
    setq(q);
    if(text) gl->text.clr() <<text;
    gl->update();
    //gl->timedupdate(getTau()*(T-1)/(display-1));
  }else{
  }
}

void soc::SocSystemAbstraction::displayTrajectory(const arr& q,const arr *invQ,int steps,const char *tag){
  uint k,t,T=nTime();
  if(!gl || !steps) return;
  uint num;
  if(steps==1 || steps==-1) num=T; else num=steps;
  for(k=0;k<=(uint)num;k++){
    t = k*T/num;
    if(invQ) displayState(q[t],&(*invQ)[t](),STRING(tag <<" (time " <<std::setw(3) <<t <<'/' <<T <<')'));
    else     displayState(q[t],NULL         ,STRING(tag <<" (time " <<std::setw(3) <<t <<'/' <<T <<')'));
    if(steps==-1) gl->watch();
  }
  if(steps==1) gl->watch();
}

//! computes separate costs for each ctrl variable
real soc::SocSystemAbstraction::analyzeTrajectory(const arr& q,bool plot){
  uint t,T=nTime(),i,m=nTasks();
  CHECK(q.nd==2 && q.d0==T+1 && q.d1==qDim(),"");
  real taskCost=0., controlCost=0.;
  arr W,H,M,F;
  real tau=1.,tau_1=1.,tau_2=1.;
  if(!dynamic){
    getW(W);
  }else{
    getH(H);
    tau=getTau();
    tau_1 = 1./tau;
    tau_2 = tau_1*tau_1;
  }

  arr phiHatQ,x,v,Jqd;
  real dx,dv,prec,precv;

  arr taskC(T+1);  taskC.setZero();
  arr ctrlC(T+1);    ctrlC.setZero();
  arr taskCi(T+1,m);  taskC.setZero();
  arr taskDx(T+1,m); taskDx.setZero();
  arr taskDv(T+1,m); taskDv.setZero();
  for(t=0;t<=T;t++){
    if(!dynamic || !t) setq(q[t]);
    else setqv(q[t],tau_1*(q[t]-q[t-1]));

    for(i=0;i<m;i++){
      if(isConditioned(i,t<<scalePower)){
        getPhi      (phiHatQ,i);
        getTarget   (x,prec,i,t<<scalePower);
        dx = sqrDistance(x,phiHatQ);
        taskCost += prec*dx*(1<<scalePower);
        taskC(t) += prec*dx;
        taskCi(t,i)=prec*dx;
        taskDx(t,i)=sqrt(dx);

        if(dynamic){
          getJqd       (Jqd,i);
          getTargetV   (v,precv,i,t<<scalePower);
          dv = sqrDistance(v,Jqd);
          taskCost += precv*dv*(1<<scalePower);
          taskC(t) += precv*dv;
          taskCi(t,i)+=precv*dv;
          taskDv(t,i)=sqrt(dv);
        }
      }
      if(isConstrained(i,t<<scalePower)){
        getPhi      (phiHatQ,i);
        //taskCost += 1.-phiHatQ;
        //taskC(t) += prec*dx;
        taskCi(t,i)=sum(phiHatQ);
        taskDx(t,i)=sum(phiHatQ);
      }
    }

    if(!dynamic){
      if(t>0) ctrlC(t) = sqrDistance(W,q[t-1],q[t]);
    }else{
      getMF(M,F);
      if(t<T && t>0)
        ctrlC(t) = sqrDistance(H,tau_2*M*(q[t+1]+q[t-1]-(real)2.*q[t]),F);
      if(t==0)
        ctrlC(t) = sqrDistance(H,tau_2*M*(q[t+1]-q[t]),F);
    }
    controlCost += ctrlC(t);
  }
#ifdef NIKOLAY
  plot = false;
#endif
  if(plot){
    std::ofstream fil;
    MT::open(fil,"z.trana");
    for(t=0;t<=T;t++){
      fil
	<<"time " <<t*tau
	<<"  ctrlC " <<ctrlC(t)
	<<"  taskC " <<taskC(t);
      fil <<"  taskCi "; taskCi[t].writeRaw(fil);
      fil <<"  taskDx "; taskDx[t].writeRaw(fil);
      fil <<"  taskDv "; taskDv[t].writeRaw(fil);
      fil <<"  q "; q[t].writeRaw(fil);
      fil <<endl;
    }
    MT::String cmd;
    cmd <<"plot 'z.trana' us 0:4 title 'ctrlC', 'z.trana' us 0:6 title 'taskC'";
    for(i=0;i<m;i++) if(isConditioned(i,0)||isConstrained(i,0)) cmd <<", 'z.trana' us 0:"<<8+i<<" title '"<<taskName(i)<<"'";
    gnuplot(cmd);
  }
  /*plotClear();
  plotFunctions(q);
  plot(false);*/
#ifdef NIKOLAY
  if(os) *os
    <<" " <<taskCost
    <<" " <<controlCost
    <<" " <<taskCost+controlCost <<endl;
#else
  if(os) *os
    <<"  task-cost " <<taskCost
    <<"  control-cost " <<controlCost
    <<"  total-cost " <<taskCost+controlCost <<endl;
#endif
  return taskCost+controlCost;
}


void soc::getQuadraticTaskCost(soc::SocSystemAbstraction& soci,int t,arr& R,arr& r,const arr& qt){
  uint m=soci.nTasks(),n=soci.qDim();
  uint i;
  arr phiHatQ,J,tJ,x;
  real prec;
  R.resize(n,n); R.setZero();
  r.resize(n);   r.setZero();
  for(i=0;i<m;i++) if(soci.isConditioned(i,t<<soci.scalePower)){
    soci.getPhi      (phiHatQ,i);
    soci.getTarget   (x,prec,i,t<<soci.scalePower);
    soci.getJtJ      (J,tJ,i);
    R += prec*tJ*J;
    r -= (real)2.*prec*tJ*(x - phiHatQ + J*qt);
  }
}



real getTaskLogLikelihood(soc::SocSystemAbstraction& soci,int t){
  uint m=soci.nTasks();
  uint i;
  real L=0.;
  arr phiHatQ,Jqd,x,v;
  real prec;
  for(i=0;i<m;i++) if(soci.isConditioned(i,t)){
    soci.getPhi(phiHatQ,i);
    soci.getTarget (x,prec,i,t);
    NIY;
    //L += logNNprec(x,phiHatQ,prec);
    if(soci.dynamic){
      soci.getJqd(Jqd,i);
      soci.getTargetV(v,prec,i,t);
      //L += logNNprec(v,Jqd,prec);
    }
  }
  return L;
}

real getFilterCostMeassure(soc::SocSystemAbstraction& soci,arr& q,real& cost1,real& cost2,std::ostream *os){
  //evaluate trajectory
  real cost_t,length=0.;
  cost1=.0;
  cost2=.0;
  real tau=soci.getTau(),tau2=tau*tau;
  arr v;
  arr W,H,invH,Q,invQ,Minv,F,tmp1,tmp2,tmp3,tmp4,d;
  soci.getW(W);
  soci.getH(H);
  soci.getQ(Q);
  inverse_SymPosDef(invQ,Q);
  inverse_SymPosDef(invH,H);
  uint t,T=q.d0;
  for(t=0;t<T;t++){
    soc::getVelocity(v,q,t,tau);
    soci.setqv(q[t],v);
    if(t>1){
      soci.getMinvF(Minv,F);
      tmp1 = (real)2.*q[t-1]+q[t-2]+tau2*Minv*F - q[t];
      cost2 += scalarProduct(invQ,tmp1,tmp1);
    }
    if(t==1){
      tmp1 = q[t-1] - q[t];
      cost2 += scalarProduct(invQ,tmp1,tmp1);
    }
    if(t>0){ d=q[t]-q[t-1]; length += ::sqrt(scalarProduct(W,d,d)); }
    CHECK(soci.dynamic==false,"");
    cost1 += cost_t = getTaskCost(soci,t);
  }
  if(os){
    *os <<std::setw(3) <<0
        <<"  time " <<MT::timerRead(false)
        <<"  cost1 " <<cost1
        <<"  cost2 " <<cost2
        <<"  length " <<length
        <<"  total-cost " <<cost1+cost2 <<endl;
  }
  return cost1+cost2;
}

