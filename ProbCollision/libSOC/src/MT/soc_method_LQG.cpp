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
#include "util.h"

/*! \brief iterated LQG (linear quadratic Gaussian) applied to the
    kinematic trajectory optimization case. This means that we assume
    A=B=1 (simple additive control) with control cost H=W
    (corresponding directly to the q-space metric). The quadratic cost
    terms are computed from the tast constraints */
real soc::iLQG::stepKinematic(){
  uint t,T=sys->nTime(),n=sys->qDim();

  //arr V(T+1,n,n),v(T+1,n);
  //arr R(T+1,n,n),r(T+1,n);
  arr HVinv(T+1,n,n),VHVinv;
  //we assume A=B=\id for now
  arr H;

  sys->getW(H);

  //remember the old trajectory
  arr q_old(q);

  if(q.nd==2 && q.d0==(T>>1)+1){ //upscale...
    q.resize(T+1,q.d1);
    for(t=0;t<=T;t+=2){
      q[t] = q_old[t>>1];  if(t<T){ q[t+1] = (real).5*(q_old[t>>1] + q_old[(t>>1)+1]); }
    }
    q_old=q;
  }
     
  //trajectory needs to be initialized!
  CHECK(q.nd==2 && q.d0==T+1 && q.d1==n,"please initialize trajectory!");
  //possible initialization routines:
  //  straightTaskTrajectory(soci,q,0);
  //  bayesianIKTrajectory(soci,q);
  //  sys->passiveDynamicTrajectory(q);


  /*MT::timerStart();
  if(sys->os){
#ifdef NIKOLAY
    *sys->os <<std::setw(3) <<-1 << " " <<MT::timerRead(false);
    sys->computeTotalCost(q,false);
#else
    *sys->os <<"iLQG " <<std::setw(3) <<-1 <<"  time " <<MT::timerRead(false);
    sys->analyzeTrajectory(q,display>0);
    //sys->computeTotalCost(q);
#endif
  }
  if(sys->gl){
    sys->displayTrajectory(q,NULL,display,STRING("iLQG optimization -- iteration "<<-1));
  }*/

  V.resize(T+1,n,n);  v.resize(T+1,n);
  R.resize(T+1,n,n);  r.resize(T+1,n);

  //linearize around current trajectory
  for(t=0;t<=T;t++){
    countSetq++;
    sys->setq(q[t]);
    getQuadraticTaskCost(*sys,t,R[t](),r[t](),q[t]);
  }

  //bwd Ricatti equations
  V[T]() = R[T];
  v[T]() = r[T];
  for(t=T;t--;){
    inverse_SymPosDef(HVinv[t+1](), H+V[t+1]);
    VHVinv = V[t+1]*HVinv[t+1];
    V[t]() = R[t] + V[t+1] - VHVinv*V[t+1];
    v[t]() = r[t] + v[t+1] - VHVinv*v[t+1];
  }

  //fwd with optimal control
  sys->getq0(q[0]());
  for(t=1;t<=T;t++){
    if(convergenceRate==1.){
      q[t]() = q[t-1] - HVinv[t]*((real).5*v[t] + V[t]*q[t-1]);
    }else{
      q[t]() = ((real)1.-convergenceRate)*q[t]
	+ convergenceRate*(q[t-1] - HVinv[t]*((real).5*v[t] + V[t]*q[t-1]));
    }
  }

  double diff = -1.;
  if(q_old.N==q.N) diff=maxDiff(q_old,q);

  //display or evaluate
  MT::timerPause();
  if(sys->os){
    *sys->os <<"iLQGk("<<sys->scalePower<<") " <<std::setw(3) << " time " <<MT::timerRead(false) <<" diff " <<diff;
    sys->analyzeTrajectory(q,display>0);
    //sys->computeTotalCost(q);
  }
  if(display){
    sys->displayTrajectory(q,NULL,display,STRING("iLQG optimization -- iteration "));
  }
  MT::timerResume();

  return diff;
}

//===========================================================================

/*! \brief iterated LQG (linear quadratic Gaussian) for the general
    Stochastic Optimal Control case */
real soc::iLQG::stepDynamic(){
  CHECK(sys->dynamic,"assumed dynamic SOC abstraction");
  uint t,T=sys->nTime(),n=sys->qDim(),n2;

  n2=2*n;
  arr A(T+1,n2,n2),a(T+1,n2),B(T+1,n2,n);
  arr At,Bt,VA,HVinv,BIG,u; //helpers
  arr H;

  sys->getH(H);
  real tau=sys->getTau();

  //remember the old trajectory
  arr q_old(q);

  //trajectory needs to be initialized!
  CHECK(q.nd==2 && q.d0==T+1 && q.d1==n,"please initialize trajectory!");

  //get velocity profile
  getPhaseTrajectory(_q,q,tau);

  V.resize(T+1,n2,n2);  v.resize(T+1,n2);
  R.resize(T+1,n2,n2);  r.resize(T+1,n2);
  
  //linearize around current trajectory
  for(t=0;t<=T;t++){
    sys->setqv(_q[t]);
    sys->getCosts  (R[t](),r[t](),t,q[t]);
    sys->getProcess(A[t](),a[t](),B[t]());
    //cout <<"t=" <<t <<" A=" <<A[t] <<" a=" <<a[t] <<" B=" <<B[t] <<endl;
    //cout <<"t=" <<t <<" R=" <<R[t] <<" r=" <<r[t] <<endl;
  }

  //bwd Ricatti equations
  V[T]() = R[T];
  v[T]() = r[T];
  for(t=T;t--;){
    transpose(At,A[t]);
    transpose(Bt,B[t]);
    innerProduct(VA,V[t+1],A[t]);
    inverse_SymPosDef(HVinv, H+Bt*V[t+1]*B[t]);
    BIG = (~VA)*B[t]*HVinv*Bt;
    V[t]() = R[t] + (At-BIG)*VA;
    v[t]() = r[t] + (At-BIG)*(v[t+1]-V[t+1]*a[t]);
    //if(t==0) cout <<"\nt=" <<t <<"\nV=" <<V[t] <<"\nv=" <<v[t] <<endl;
  }

  //fwd with optimal control
  //we assume that _q[0] is fine and fixed!    (from getPhaseTrajectory)
  real ctrlC=0;
  for(t=0;t<T;t++){
    sys->setqv(_q[t]);
    sys->getProcess(A[t](),a[t](),B[t]());
    //if(t>T-10) cout <<"t=" <<t <<"\nhatq=" <<_q[t] <<endl;

    transpose(Bt,B[t]);
    inverse_SymPosDef(HVinv, H+Bt*V[t+1]*B[t]);
    u = - HVinv*Bt*(V[t+1]*(A[t]*_q[t]+a[t]) - v[t+1]);
    ctrlC += scalarProduct(H,u,u);
    //cout <<"ilqg time "<<t<<" u=" <<u <<" q=" <<_q[t] <<endl;
#if 0
    _q[t+1]() = A[t]*_q[t] + a[t] + B[t]*u;
#else
    //if(t>T-10) cout <<"t=" <<t+1 <<"\nb=" <<A[t]*_q[t] + a[t] + B[t]*u <<endl;
      
    if(convergenceRate==1.){
      _q[t+1]() = A[t]*_q[t] + a[t] + B[t]*u;
    }else{
      _q.reshape(T+1,2,n);
      _q.subDim(t+1,1)
	= ((real)1.-convergenceRate)*_q.subDim(t+1,1)
	+ convergenceRate*(_q.subDim(t,1) + a[t].sub(n,-1) + B[t].sub(n,-1,0,-1)*u);
      _q.subDim(t+1,0) = _q.subDim(t,0) + tau*_q.subDim(t+1,1);
      _q.reshape(T+1,2*n);
    }
#endif
  }

  _q.reshape(T+1,2,n);
#if 0
  if(convergenceRate==1.){
    for(t=0;t<T;t++) q[t]=_q.subDim(t,0);
  }else{
    for(t=0;t<T;t++) q[t]=(1.-convergenceRate)*q[t] + convergenceRate*_q.subDim(t,0);
  }
#else
  for(t=0;t<=T;t++) q[t]=_q.subDim(t,0);
#endif

  _q.reshape(T+1,2*n);
  double diff = -1.;
  getPositionTrajectory(q,_q);
  if(q_old.N==q.N) diff=maxDiff(q_old,q);

  //display or evaluate
  MT::timerPause();
  if(sys->os){
    *sys->os <<"iLQGd("<<sys->scalePower<<") " <<std::setw(3) <<" time " <<MT::timerRead(false) <<" diff " <<diff;
    sys->analyzeTrajectory(q,display>0);
  }
  if(sys->gl){
    sys->displayTrajectory(q,NULL,display,STRING("iLQG optimization -- iteration "));
  }
  MT::timerResume();

  return diff;
}
