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
#include "truncatedGaussian.h"

//#define TightMode

//===========================================================================
//
// helpers - fwd declarations
//

void setBlockVec(arr& X, real a0,real a1, const arr& x);
void setBlockVec(arr& X,const arr& x1,const arr& x2);
void setBlockDiagMat(arr& X,const arr& A,const arr& B);
void setBlockMat(arr& X, real a00,real a01,real a10,real a11, const arr& x);
void blockTimesVec(arr& Y, real a00,real a01,real a10,real a11, const arr& X);
void bl_ockTimesMat(arr& Y, real a00,real a01,real a10,real a11, const arr& X);
void blockTransMat(arr& Y, real a,real b,real c,real d, const arr& X);
real getDynamicCostMeassure(soc::SocSystemAbstraction& soci,arr& q,real& cost1,real& cost2,std::ostream *os);
real getTransitionLogLikelihood(soc::SocSystemAbstraction& soci,arr& b0,arr& b1);

void computeEPmessage(arr& a,arr &invA,const arr& b_from,const arr& invB_from,const arr& b_to,const arr& invB_to){
  invA = invB_to - invB_from;
  lapack_Ainv_b_sym(a, invA, invB_to*b_to - invB_from*b_from);
}

//===========================================================================
//
// methods
//

/*! \brief compute a single control step from current state to state
    at time t (taking into account all targets and precisions
    associated with all declared task variables) */
void soc::bayesianIKControl(SocSystemAbstraction& soci,
                            arr& dq,uint t){
  uint n=soci.qDim(),m=soci.nTasks();
  dq.resize(n);
  dq.setZero();
  arr a(n),A(n,n),invA(n,n);
  arr target,actual,J,tJ;
  real prec;
  soci.getW(A);
  a.setZero();
  //arr w(3);
  for(uint i=0;i<m;i++) if(soci.isConditioned(i,t)){
    soci.getPhi   (actual,i);
    soci.getJtJ   (J,tJ,i);
    soci.getTarget(target,prec,i,t);
    a += prec * tJ * (target-actual);
    A += prec * tJ * J;
  }
  inverse_SymPosDef(invA,A);
  dq = invA * a;
}

void soc::bayesianIKControl2(SocSystemAbstraction& soci,
                             arr& q,const arr& q_1,uint t,arr *v, arr *invV){
  CHECK(!soci.dynamic,"assumed non-dynamic SOC abstraction");
  uint n=soci.qDim();
  
  //-- access necessary information
  arr W;
  soci.getW(W);
  
  //fwd message
  arr s(n),invS(n,n);
  s = q_1;
  invS = W;
  
  //task message
  arr R,r;
  soci.getCosts(R,r,t,q_1);
  
  //v,invV are optional bwd messages!

  //belief
  arr invB,b;
  if(!v){
    invB = invS + R;
    lapack_Ainv_b_sym(b, invB, invS*s + r);
  }else{
    invB = invS + (*invV) + R;
    lapack_Ainv_b_sym(b, invB, invS*s + (*invV)*(*v) + r);
  }
  
  //constraints
  arr cdir,coff;
  soci.getConstraints(cdir,coff,t,q_1);
  if(cdir.d0){
    arr _b,_B,_invB;
    inverse_SymPosDef(_B,invB);
    _b=b;
    //plotClear();  plotCovariance(_b,_B);
    for(uint i=0;i<cdir.d0;i++){ //one-by-one truncate the constraint from the belief
      TruncateGaussian(_b,_B,cdir[i],coff(i));
      //plotTruncationLine(cdir[i],coff[i]);  plotCovariance(_b,_B);  plot();
    }
    //compute the EP message and 'add' it to the task message
    inverse_SymPosDef(_invB,_B);
    R += _invB - invB;
    r += _invB * _b - invB*b;

    //recompute (b,B);
    invB = invS + R;
    lapack_Ainv_b_sym(b, invB, invS*s + r);
  }



  q=b;
}
                            
/*! \brief standard IK -- follows the first (active) task variable
    exactly, but fails if more than one task variable is
    active. regularization=make it singularity robust*/
void soc::pseudoIKControl(SocSystemAbstraction& soci, arr& dq, uint t,real regularization){
  uint n=soci.qDim(),m=soci.nTasks();
  dq.resize(n);
  dq.setZero();
  arr invJ,W,invW;
  arr target,actual,J,tJ;
  real prec;
  soci.getW(W);
  inverse_SymPosDef(invW,W);
  uint k=0;
  for(uint i=0;i<m;i++) if(soci.isConditioned(i,t)){
    CHECK(!k,"pseudo IK only works for a single task variable");
    soci.getPhi   (actual,i);
    soci.getJtJ   (J,tJ,i);
    soci.getTarget(target,prec,i,t);
    pseudoInverse(invJ,J,invW,regularization);
    dq = invJ * (target-actual);
    k++;
  }
}

/*! \brief hierarchical IK: follows the 1st task variable exactly, the
    2nd exactly up to the 1st, etc.., might be come
    brittle. regularization=make it singularity robust */
void soc::hierarchicalIKControl(SocSystemAbstraction& soci,arr& dq,uint t,real regularization){
  uint n=soci.qDim(),m=soci.nTasks();
  dq.resize(n);
  dq.setZero();
  arr hatJ,invHatJ,N,W,invW;
  arr target,actual,J,tJ;
  real prec;
  N.setId(n);
  soci.getW(W);
  inverse_SymPosDef(invW,W);
  for(uint i=0;i<m;i++) if(soci.isConditioned(i,t)){
    soci.getPhi   (actual,i);
    soci.getJtJ   (J,tJ,i);
    soci.getTarget(target,prec,i,t);
    hatJ = J * N;
    pseudoInverse(invHatJ,hatJ,invW,regularization);
    dq += invHatJ * ((target-actual) - J * dq);
    N  -= invHatJ * hatJ;
  }
}

//===========================================================================

/*! \brief compute a single control step from state at time t-1 to
    state at time t. If eps=0, this is equivalen to bayesianIKControl;
    for eps>0 the IK step is repeated until convergence up to
    tolerance eps. qt=output, qt_1=state at time t-1. */
void soc::bayesianIterateIKControl(SocSystemAbstraction& soci,
                                   arr& qt,const arr& qt_1,uint t,real eps,uint maxIter){
  uint j;
  if(&qt!=&qt_1) qt=qt_1;
  arr dq;
  for(j=0;j<maxIter;j++){
    soci.setq(qt);
    bayesianIKControl(soci,dq,t);
    if(j<3) qt+=dq;
    //else if(j<20) qt+=.8*dq;
    else qt+=(real).8*dq;
    if(dq.absMax()<eps) break;
    //cout <<"IK iteration " <<j <<" dq=" <<dq <<endl;
  }
  if(j==maxIter) HALT("warning: IK didn't converge (|last step|="<<dq.absMax()<<")")
  else cout <<"IK converged after steps=" <<j <<endl;
}

//===========================================================================


/*! \brief compute a single control step from current state to target of time t.
    qv=output, qv_1=state at time t-1 */
void soc::bayesianDynamicControl(SocSystemAbstraction& soci, arr& qv, const arr& qv_1, uint t, arr *v, arr *invV){
  CHECK(soci.dynamic,"assumed dynamic SOC abstraction");
  uint n=soci.qDim();
  
  //-- access necessary information
  arr A,a,B,tB;
  soci.getProcess(A,a,B);
  transpose(tB,B);
  
  arr Q,H,invH;
  soci.getQ(Q);
  soci.getH(H);
  inverse_SymPosDef(invH,H);
  
  //fwd message
  arr s(2*n),S(2*n,2*n),invS(2*n,2*n);
  S = Q;
  S += B*invH*tB;
  s = a + A*qv_1;
  inverse_SymPosDef(invS,S);

  //task message
  arr R,r,q_1;
  q_1.referToSubRange(qv_1,0,n-1);
  soci.getCosts(R,r,t,q_1);

  //v,invV are optional bwd messages!
  
  //belief
  arr invB,b;
  if(!v){
    invB = invS + R;
    lapack_Ainv_b_sym(b, invB, invS*s + r);
  }else{
    invB = invS + (*invV) + R;
    lapack_Ainv_b_sym(b, invB, invS*s + (*invV)*(*v) + r);
  }
  
  qv=b;
}
                           

//===========================================================================

/*! \brief compute a trajectory using inverse kinematics (iterating
    bayesianIK forward). If eps=0, no IK step is repeated, for eps>0
    IK steps are repeated until they converge up to tolerance eps (see
    \ref bayesianIterateIKControl) */
void soc::bayesianIKTrajectory(SocSystemAbstraction& soci,arr& q,real eps){
  uint t,T=soci.nTime(),n=soci.qDim();
  q.resize(T+1,n);
  arr dq;
  soci.getq0(q[0]());
  for(t=1;t<=T;t++){
    if(eps<0){
      soci.setq(q[t-1]);
      //bayesianIKControl(soci,dq,t);   q[t]() = q[t-1] + dq;
      bayesianIKControl2(soci,q[t](),q[t-1],t);
    }else{
      bayesianIterateIKControl(soci,q[t](),q[t-1],t,eps,20);
    }
  }
}

/*
void bayesianIterateIKTrajectory(SocSystemAbstraction& soci,arr& q,real eps,uint maxIter){
  uint t,T=soci.nTime(),n=soci.qDim();
  q.resize(T,n);
  arr dq;
  soci.getq0(q[0]());
  for(t=1;t<T;t++){
    bayesianIterateIKControl(soci,q[t](),q[t-1],t,eps,maxIter);
}
}
*/


//===========================================================================

//#define CANON

void lapack_A_Binv_A_sym(arr& X,const arr& A,const arr& B){
  static arr Binv,tmp;
  inverse_SymPosDef(Binv,B);
  blas_MM(tmp,Binv,A);
  blas_MM(X,A,tmp);
}

void soc::AICO::initMessages(){
  uint n=sys->qDim();
  uint T=sys->nTime();
  arr q0;
  if(!sys->dynamic){
    sys->getq0(q0);
  }else{
    sys->getqv0(q0);
    n*=2;
  }
  s.resize(T+1,n);  invS.resize(T+1,n,n);  s[0]=q0;      invS[0].setDiag(1e10);
  v.resize(T+1,n);  invV.resize(T+1,n,n);  v.setZero();  invV.setZero();
  b.resize(T+1,n);  invB.resize(T+1,n,n);  b[0]=q0;      invB[0].setDiag(1e10);
  r.resize(T+1,n);  R.resize(T+1,n,n);     r[0]=0.;      R[0]=0.;
  hatq.resize(T+1,n);                      hatq[0]=q0;

  //locks.resize(T+1);

#ifdef CANON
//CHECK OLD VERSIONS:
   //I once tried to implement the algo in the canonical representation - didn't work better
  s[0]=invS[0]*q0;
#endif
}

#if 0
void soc::AICO::getMultiScaleMessages(arr& s_,arr& S_,arr& v_,arr& V_,uint t,real upMixS,real selfMixS,real dnMixS,real upMixV,real selfMixV,real dnMixV){
  CHECK(multiScales(scale)==this,"");
  arr s,S,v,V;
  s.resize(b.d1);      s.setZero();
  S.resize(b.d1,b.d1); S.setZero();
  v.resize(b.d1);      v.setZero();
  V.resize(b.d1,b.d1); V.setZero();
  AICO *A;
  real pS,pSsum=0.,pV,pVsum=0.;;
  uint i,t_;
  for_list(i,A,multiScales){
    CHECK(i==A->scale,"");
    if(!A->s.N) continue; //assume this wasn't computed yet..
    if(i <scale){
      pS=dnMixS;   pV=dnMixV;   t_=t<<(scale-i);
    }
    if(i==scale){ pS=selfMixS; pV=selfMixV; t_=t; }
    if(i >scale){
      pS=upMixS;   pV=upMixV;   t_=t>>(i-scale);
      if(t<<(i-scale)!=t) continue; //don't down-share messages intermediate steps
    }
    s += pS*A->s[t_];   S += pS*(A->S[t_] + (A->s[t_]^A->s[t_]));
    v += pV*A->v[t_];   V += pV*(A->V[t_] + (A->v[t_]^A->v[t_]));
    pSsum += pS;
    pVsum += pV;
  }
  CHECK(pSsum>1e-10 && pVsum>1e-10,"");
  s/=pSsum;  S/=pSsum;
  v/=pVsum;  V/=pVsum;
  S -= s^s;
  V -= v^v;
  s_=s;  S_=S;
  v_=v;  V_=V;
}

void soc::AICO::getMultiScaleMessages(arr& s_,arr& S_,arr& v_,arr& V_,uint t,real upMixS,real selfMixS,real dnMixS,real upMixV,real selfMixV,real dnMixV){
  NIY;
#if 0
  //get parallel messages in multi-scale case
  if(multiScales.N){
    real mix=::pow(.8,sweep);
    if(sweep<2) mix=1.; else{ if(sweep<4) mix=.5; else mix=0.; }
    if(scale+1<multiScales.N){
      AICO *A = multiScales(scale+1);
      arr s_left,s_mid,invS_left,invS_mid,v_mid,v_right,invV_mid,invV_right,hatq_left,hatq_mid,hatq_right;
      if(t&1){ //t-1 and t+1 are even
	uint tLeft =(t-1)>>1;
	uint tRight=(t+1)>>1;
	{ A->locks(tLeft) .readLock(); s_left =A->s[tLeft ]; invS_left =A->invS[tLeft ]; hatq_left =A->hatq[tLeft ]; A->locks(tLeft ).unlock(); }
	{ A->locks(tRight).readLock(); v_right=A->v[tRight]; invV_right=A->invV[tRight]; hatq_right=A->hatq[tRight]; A->locks(tRight).unlock(); }
	/*if(t>0 && dt==-1){
	  locks(t-1).writeLock(STRING("-1 scale="<<scale <<" sweep=" <<sweep <<" t=" <<t));
	  s   [t-1] = (1.-mix)*s   [t-1] + mix*   s_left;
	  invS[t-1] = (1.-mix)*invS[t-1] + mix*invS_left;
	  locks(t-1).unlock();
	  }
	  if(t<T && dt==+1){
	  locks(t+1).writeLock(STRING("+1 scale="<<scale <<" sweep=" <<sweep <<" t=" <<t));
	  v   [t+1] = (1.-mix)*v   [t+1] + mix*   v_right;
	  invV[t+1] = (1.-mix)*invV[t+1] + mix*invV_right;
	  locks(t+1).unlock();
	  }*/
	hatq[t] = (1.-mix)*hatq[t] + mix*.5*(hatq_left + hatq_right);
      }else{ //t is even
	uint tLeft =(t>>1)-1;
	uint tMid  =(t>>1);
	uint tRight=(t>>1)+1;
	if(t>1 && dt==-1){ A->locks(tLeft) .readLock(); s_left =A->s[tLeft ]; invS_left =A->invS[tLeft ]; hatq_left =A->hatq[tLeft ]; A->locks(tLeft ).unlock(); }
	A->locks(tMid)  .readLock(); s_mid=A->s[tMid]; invS_mid=A->invS[tMid]; v_mid=A->v[tMid]; invV_mid=A->invV[tMid]; hatq_mid=A->hatq[tMid]; A->locks(tMid).unlock();
	if(t<T && dt==+1){ A->locks(tRight).readLock(); v_right=A->v[tRight]; invV_right=A->invV[tRight]; hatq_right=A->hatq[tRight]; A->locks(tRight).unlock(); }
	/*if(t>1 && dt==-1){
	  locks(t-1).writeLock(STRING("-1 scale="<<scale <<" sweep=" <<sweep <<" t=" <<t));
	  s   [t-1] = (1.-mix)*s   [t-1] + mix*.5*(s_mid+s_left);
	  invS[t-1] = (1.-mix)*invS[t-1] + mix*.5*(invS_mid+invS_left);
	  locks(t-1).unlock();
	  }
	  if(t<T && dt==+1){
	  locks(t+1).writeLock(STRING("+1 scale="<<scale <<" sweep=" <<sweep <<" t=" <<t));
	  v   [t+1] = (1.-mix)*v   [t+1] + mix*.5*(v_mid+v_right);
	  invV[t+1] = (1.-mix)*invV[t+1] + mix*.5*(invV_mid+invV_right);
	  locks(t+1).unlock();
	  }*/
	hatq[t] = (1.-mix)*hatq[t] + mix*hatq_mid;
      }
      useFwdMessageAsInitialHatq=false;
    }
  }
#endif
}
#endif

soc::AICO* soc::AICO_solver(SocSystemAbstraction& sys,
                      arr& q,double tolerance,
                      real convergenceRate,real repeatThreshold, real recomputeTaskThreshold,
                      uint display){
  soc::AICO *aico=new soc::AICO;
  aico->init(sys,convergenceRate,repeatThreshold,recomputeTaskThreshold,display,0,LIST<AICO>());
  aico->q=q;
  for(uint k=0;k<100;k++){
    double d;
    if(!sys.dynamic) d=aico->stepKinematic();
    else             d=aico->stepDynamic();
    if(k && d<tolerance) break;
  }
  q=aico->q;
  return aico;
}

soc::AICO* soc::AICO_multiScaleSolver(SocSystemAbstraction& sys,
                                arr& q,
                                double tolerance,
                                real convergenceRate,real repeatThreshold, real recomputeTaskThreshold,
                                uint display,
                                uint scalePowers){
  MT::Array<soc::AICO*> aicos(scalePowers);
  for(uint i=0;i<aicos.N;i++) aicos(i) = new soc::AICO;
  for(uint i=0;i<aicos.N;i++){
    sys.scalePower=i;
    aicos(i)->init(sys,convergenceRate,repeatThreshold,recomputeTaskThreshold,display,i,aicos);
  }
  for(uint i=aicos.N;i--;) for(int k=0;k<100;k++){
    sys.scalePower=i;
    double d;
    if(!sys.dynamic) d=aicos(i)->stepKinematic();
    else             d=aicos(i)->stepDynamic();
    if(k && d<tolerance) break;
  }
  q=aicos(0)->q;
  return aicos(0);
}

/*void soc::AICO_kinematicTol(SocSystemAbstraction& soci,
                            arr& q,real tolerance,
                            real convergenceRate,real repeatThreshold,
                            uint display,
                            uint scalePower){
  soc::AICO aico;
  aico.scale=scalePower;
  aico.q=q;
  aico.soci=&soci;
  for(uint k=0;;k++){
    real d=aico.stepKinematic(convergenceRate,repeatThreshold,display);
    q=aico.q;
    if(k && d<tolerance) break;
  }
}*/

//! Approximate Inference Control (AICO) in the kinematic case
real soc::AICO::stepKinematic(){
  CHECK(!sys->dynamic,"assumed dynamic SOC abstraction");
  MT_DEBUG(uint n=sys->qDim());
  uint T=sys->nTime();
  uint t,t0=0;
  int dt;

  //variables for the dynamics
  arr invW;
  sys->getWinv(invW);
  
  //temporary variables
  arr Vt,St,barS,barV,K,K2;

  //initializations (initial q or multiscale)
  bool useFwdMessageAsInitialHatq=true;
  if(!sweep){
    /*if(multiScales.N){
      if(scale+1<multiScales.N){
        AICO *A = multiScales(scale+1);
        arr q_sub = A->q;
        CHECK(q_sub.nd==2 && q_sub.d0==(T>>1)+1 && q_sub.d1==n,"initial trajectory was wrong dimensionality");
        q.resize(T+1,n);
        for(t=0;t<=T;t+=2){
          q[t] = q_sub[t>>1];  if(t<T){ q[t+1] = (real).5*(q_sub[t>>1] + q_sub[(t>>1)+1]); }
        }
      }
    }*/

    if(q.N){
      CHECK(q.nd==2 && q.d0==T+1 && q.d1==n,"initial trajectory was wrong dimensionality");
      useFwdMessageAsInitialHatq=false;
      hatq=q;
      b=hatq;
      v=hatq;  for(uint t=0;t<=T;t++){ invV[t].setDiag(1e6);  }
      if(sys->os){//type initial value
        *sys->os <<"AICO " <<std::setw(3) <<-1 <<" time " <<MT::timerRead(false);
        sys->computeTotalCost(q);
      }
      if(sys->gl){
        sys->displayTrajectory(q,NULL,display,STRING("AICO_kinematic - iteration -1"));
      }
    }

    //-- initialize messages from lower scale
    if(multiScales.N){
      if(scale+1<multiScales.N){
        AICO *A = multiScales(scale+1);
        for(t=0;t<=T;t+=2){
          //s[t] = A->s[t>>1]; invS[t] = A->invS[t>>1];   if(t<T){ s[t+1]=s[t]; invS[t+1]=invS[t]; }  //don't need to copy fwd messages
          v[t] = A->v[t>>1]; invV[t] = A->invV[t>>1];  if(t<T){ v[t+1]=A->v[(t>>1)+1]; invV[t+1]=A->invV[(t>>1)+1]; }
//           if(t<T){ v[t+1]=.5*(A->v[t>>1] + A->v[(t>>1)+1]); invV[t+1]=.5*(A->invV[t>>1] + A->invV[(t>>1)+1]); }
          hatq[t] = A->hatq[t>>1];  if(t<T){ hatq[t+1] = (real).5*(A->hatq[t>>1] + A->hatq[(t>>1)+1]); }
          //b   [t] = A->b   [t>>1];   if(t<T) b   [t+1] = .5*(A->b   [t>>1] + A->b   [(t>>1)+1]);    //don't need to copy the belief
        }
        useFwdMessageAsInitialHatq=false;
      }
    }
  }

  //remember the old trajectory and old hatq
  arr q_old(q),hatq_old(hatq);
  
  //MT::timerStart();
  uint repeatCount=0;

  for(dt=1;dt>=-1;dt-=2){ //fwd & bwd
    if(dt==1)  t0=1;
    if(dt==-1) t0=T;
    for(t=t0;t<=T && t>0;t+=dt){
      //if(multiScales.N) locks(t).writeLock(STRING("ME scale="<<scale <<" sweep=" <<sweep <<" t=" <<t));
      //if(multiScales.N) getMultiScaleMessages(arr& s_,arr& S_,arr& v_,arr& V_,uint t,real upMixS,real selfMixS,real dnMixS,real upMixV,real selfMixV,real dnMixV){
      
      //compute (s,S)
      if(dt==1 && !repeatCount){ //only on fwd pass and non-repeats
        countMsg++;
#ifndef TightMode
#if 1 //natural
        inverse_SymPosDef(barS,invS[t-1] + R[t-1]);
        s[t] = barS * (invS[t-1]*s[t-1] + r[t-1]);
        St = invW + barS;
        inverse_SymPosDef(invS[t](), St);
#else //canonical (matrix multiplications become slow...)
        barS = invS[t-1] + R[t-1];
        lapack_Ainv_b_sym(s[t](), barS, invS[t-1]*s[t-1] + r[t-1]);
        inverse_SymPosDef(K,barS + W);
        //invS[t] = barS - barS*K*barS;
        blas_MsymMsym(K2,K,barS); blas_MsymMsym(K,barS,K2); invS[t] = barS - K;
#endif
        //cout <<"s\n" <<s[t] <<endl <<invS[t] <<endl;
#else
        s[t] = hatq[t-1];
        St = invW;
        inverse_SymPosDef(invS[t](),St);
#endif
      }

      //compute (v,V)
      if(dt==-1 && !repeatCount){ //only on bwd pass and non-repeats
        countMsg++;
        if(t<T){
#if 1 //natural
          inverse_SymPosDef(barV,invV[t+1] + R[t+1]);   //eq (*)
          v[t] = barV * (invV[t+1]*v[t+1] + r[t+1]);
          Vt = invW + barV;
          inverse_SymPosDef(invV[t](), Vt);
#else //canonical
          barV = invV[t+1] + R[t+1];
          lapack_Ainv_b_sym(v[t](), barV, invV[t+1]*v[t+1] + r[t+1]);
          inverse_SymPosDef(K,barV + W);
          //invV[t] = barV - barV*K*barV;
          blas_MsymMsym(K2,K,barV); blas_MsymMsym(K,barV,K2); invV[t] = barV - K;
#endif
        }
        if(t==T){ //last time slice
          v[t] = hatq[t]; //alternatives: hatq or b
#ifndef TightMode
          invV[t].setDiag(1e-0); //regularization, makes eq (*) above robust
#else
          invV[t].setDiag(1e-1); //regularization, makes eq (*) above robust
#endif
        }
        //cout <<"v\n" <<v[t] <<endl <<invV[t] <<endl;
      }
        
      //first sweep: set hatq equal to fwd message
      if(!sweep && useFwdMessageAsInitialHatq) hatq[t]()=s[t];

      //compute (r,R)
      if(!sweep || maxDiff(hatq[t],hatq_old[t])>recomputeTaskThreshold){ //recompute only when significant change of state
        countSetq++;
        sys->setq(hatq[t]);
        arr Rt,rt;
        sys->getCosts(Rt,rt,t,hatq[t]);
        if(!sweep){
          R[t] = Rt; r[t] = rt;
        }else{
          double eps=5./sweep;
          if(eps>1.) eps=1.;
          R[t] = (1.-eps)*R[t] + eps*Rt;
          r[t] = (1.-eps)*r[t] + eps*rt;
        }
        hatq_old[t]() = hatq[t];
        //cout <<"r\n" <<r[t] <<endl <<R[t] <<endl;
      }
      //else cout <<"skip." <<flush;
        
      //compute (b,B);
      invB[t] = invS[t] + invV[t] + R[t];
      lapack_Ainv_b_sym(b[t](), invB[t], invS[t]*s[t] + invV[t]*v[t] + r[t]);
      //cout <<"b\n" <<b[t] <<endl <<B[t] <<endl;

#if USE_TRUNCATION //PRELIMINARY - hard constraints handled with truncating Gaussians
      //sys->displayState(b[t],&invB[t],STRING("AICO kinematic (online) t="<<t));
      //account for constraints:
      arr cdir,coff;
      sys->setq(b[t]);
      //sys->gl->watch(STRING("time "<<t));
      sys->getConstraints(cdir,coff,t<<scale,b[t]);
      if(cdir.d0){
        //cout <<"t=" <<t <<' ';
        arr _b,_B,_invB;
        inverse_SymPosDef(_B,invB[t]);
        _b=b[t];
        //plotClear();  plotCovariance(_b,_B);
        for(uint i=0;i<cdir.d0;i++){ //one-by-one truncate the constraint from the belief
          TruncateGaussian(_b,_B,cdir[i],coff(i));
          //plotTruncationLine(cdir[i],coff[i]);  plotCovariance(_b,_B);  plot();
        }
        //compute the EP message and 'add' it to the task message
        inverse_SymPosDef(_invB,_B);
        R[t]() += _invB - invB[t];
        r[t]() += _invB * _b - invB[t]*b[t];

        //recompute (b,B);
        invB[t] = invS[t] + invV[t] + R[t];
        lapack_Ainv_b_sym(b[t](), invB[t], invS[t]*s[t] + invV[t]*v[t] + r[t]);
        //cout <<"b\n" <<b[t] <<endl <<B[t] <<endl;

        //sys->displayState(b[t],&invB[t],STRING("AICO kinematic (after truncation) t="<<t));
      }
#endif

#ifndef TightMode
      //decide on \hat q
      if(sweep || !useFwdMessageAsInitialHatq){
        if(convergenceRate) hatq[t]()=((real)1.-convergenceRate)*hatq[t] + convergenceRate*b[t];
        else hatq[t]()=b[t];
      }
#else
      //update hatq
      if(dt==1){
        if(convergenceRate) hatq[t]()=((real)1.-convergenceRate)*hatq[t] + convergenceRate*b[t];
        else hatq[t]()=b[t];
      }
#endif

      //if(multiScales.N) locks(t).unlock();

      //decide whether to repeat this time slice
      if(sweep && repeatThreshold && t!=T){
        real off=sqrDistance(b[t],hatq[t]);  //sqrDistance(W,b[t],hatq[t]);
        if(off>repeatThreshold){
          //cout <<t <<" REPEAT: off=" <<off <<" (repeatCount=" <<repeatCount <<")" <<endl;
          if(repeatCount<20){
            t-=dt;
            repeatCount++;
          }else{
            cout <<" ** no convergence! ** (skipping repeat) at t=" <<t <<endl;
            repeatCount=0;
          }
        }else{
          repeatCount=0;
        }
      }
    } //loop t
    sweep++;
#ifdef TightMode
    b=hatq;
#endif

  }//loop over dt in {-1,1}
    
#if 1
  q = b;
#else
  getControlledTrajectory(q,*this);
#endif
  double diff = -1.;
  if(q_old.N==q.N) diff=maxDiff(q_old,q);
  
  //display or evaluate
  MT::timerPause();
  if(sys->os){
#ifdef NIKOLAY
    *sys->os <<std::setw(3) <<sweep <<"  " <<MT::timerRead(false);
#else
    *sys->os <<"AICOk("<<scale<<") " <<std::setw(3) <<sweep <<" time " <<MT::timerRead(false) <<" diff " <<diff;
#endif
    sys->analyzeTrajectory(q,display>0);
  }
  
  if(display){
    sys->displayTrajectory(q,NULL,display,STRING("AICO_kinematic - sweep "<<sweep)); //&invB
  }
  MT::timerResume();

  return diff;
}

//! Approximate Inference Control (AICO) in the general (e.g. dynamic) case
real soc::AICO::stepDynamic(){
  CHECK(sys->dynamic,"assumed dynamic SOC abstraction");
  uint n=sys->qDim();
  uint T=sys->nTime();
  uint n2=2*n;
  uint t,t0=0;
  int dt;

  //variables for the dynamics
  arr q0;
  A.resize(T+1,n2,n2);  tA.resize(T+1,n2,n2);  invA.resize(T+1,n2,n2);  invtA.resize(T+1,n2,n2);  a.resize(T+1,n2);  B.resize(T+1,n2,n);  tB.resize(T+1,n,n2); //fwd dynamics
  arr Q,H,invH;
  sys->getqv0(q0);
  sys->getQ(Q);
  sys->getH(H);
  inverse_SymPosDef(invH,H);
  CHECK(q0.N==n2,"");
  sys->setqv(q0);
  sys->getProcess(A[0](),tA[0](),invA[0](),invtA[0](),a[0](),B[0](),tB[0]());

  //temporary variables
  arr St,Vt,barS,barV,invR;
  
  //initializations (initial q or multiscale)
  bool useFwdMessageAsInitialHatq=true;
  if(!sweep){
    if(q.N){
      CHECK(q.nd==2 && q.d0==T+1 && q.d1==n,"initial trajectory was wrong dimensionality");
      useFwdMessageAsInitialHatq=false;
      getPhaseTrajectory(hatq,q,sys->getTau());
      b=hatq;
      v=hatq;  for(uint t=0;t<=T;t++){ invV[t].setDiag(1e6);  }
      if(sys->os){//type initial value
        *sys->os <<"AICO " <<std::setw(3) <<-1 <<" time " <<MT::timerRead(false);
        sys->computeTotalCost(q);
      }
      if(sys->gl){
        sys->displayTrajectory(q,NULL,display,STRING("AICO_dynamic - iteration -1"));
      }
    }

    //-- initialize messages from lower scale
    if(multiScales.N){
      if(scale+1<multiScales.N){
        AICO *A = multiScales(scale+1);
        for(t=0;t<=T;t+=2){
          //s[t] = A->s[t>>1]; invS[t] = A->invS[t>>1];   if(t<T){ s[t+1]=s[t]; invS[t+1]=invS[t]; }  //don't need to copy fwd messages
          v[t] = A->v[t>>1]; invV[t] = A->invV[t>>1];  if(t<T){ v[t+1]=A->v[(t>>1)+1]; invV[t+1]=A->invV[(t>>1)+1]; }
//           if(t<T){ v[t+1]=.5*(A->v[t>>1] + A->v[(t>>1)+1]); invV[t+1]=.5*(A->invV[t>>1] + A->invV[(t>>1)+1]); }
          hatq[t] = A->hatq[t>>1];  if(t<T){ hatq[t+1] = (real).5*(A->hatq[t>>1] + A->hatq[(t>>1)+1]); }
          //b   [t] = A->b   [t>>1];   if(t<T) b   [t+1] = .5*(A->b   [t>>1] + A->b   [(t>>1)+1]);    //don't need to copy the belief
        }
        useFwdMessageAsInitialHatq=false;
      }
    }
  }

  //remember the old trajectory and old hatq
  double diff = -1.;
  arr q_old(q),hatq_old(hatq);

  //MT::timerStart();
  uint repeatCount=0;

  //for(uint k=0;k<2;k++){
  for(dt=1;dt>=-1;dt-=2){ //fwd & bwd
    if(dt==1)  t0=1;
    if(dt==-1) t0=T;
    for(t=t0;t<=T && t>0;t+=dt){
      //compute (s,S)
      if(dt==1 && !repeatCount){ //only on fwd pass and non-repeats
#ifndef TightMode
	inverse_SymPosDef(barS,invS[t-1] + R[t-1]);
	St = Q;
	St += B[t-1]*invH*tB[t-1];
	St += A[t-1]*barS*tA[t-1];
	s[t] = a[t-1] + A[t-1]*(barS*(invS[t-1]*s[t-1] + r[t-1]));
	inverse_SymPosDef(invS[t](),St);
	//cout <<"s\n" <<s[t] <<endl <<invS[t] <<endl;
#else
	St = Q;
	St += B[t-1]*invH*tB[t-1];
	s[t] = a[t-1] + A[t-1]*hatq[t-1];
	inverse_SymPosDef(invS[t](),St);
#endif
      }

      //compute (v,V)
      if(dt==-1 && !repeatCount){ //only on bwd pass and non-repeats
	if(t<T){
	  inverse_SymPosDef(barV,invV[t+1] + R[t+1]);
	  //cout <<"R[t+1]=" <<R[t+1] <<"invV[t+1]=" <<invV[t+1] <<"barV=" <<barV <<endl;
	  Vt = Q;
	  Vt += B[t]*invH*tB[t];
	  Vt += barV;
	  Vt = invA[t]*Vt*invtA[t];
	  v[t] = invA[t]*(-a[t] + barV*(invV[t+1]*v[t+1] + r[t+1]));
	  inverse_SymPosDef(invV[t](), Vt);
	}
	if(t==T){  //last time slice
	  v[t] = b[t]; //alternative: hatq
#ifndef TightMode
          invV[t].setDiag(1e-0); //regularization, makes eq (*) above robust
#else
          invV[t].setDiag(1e-1); //regularization, makes eq (*) above robust
#endif
        }
	//cout <<"v\n" <<v[t] <<endl <<invV[t] <<endl;
      }

      if(!sweep && useFwdMessageAsInitialHatq) hatq[t]()=s[t];
        
      //compute (r,R)
      if(!sweep || maxDiff(hatq[t],hatq_old[t])>recomputeTaskThreshold){ //recompute only when significant change of state
        countSetq++;
        sys->setqv(hatq[t]);
        arr Rt,rt;
        sys->getCosts(Rt,rt,t,hatq[t].sub(0,n-1));
        if(!sweep){
          R[t] = Rt; r[t] = rt;
        }else{
          double eps=5./sweep;
          if(eps>1.) eps=1.;
          R[t] = (1.-eps)*R[t] + eps*Rt;
          r[t] = (1.-eps)*r[t] + eps*rt;
        }
	sys->getProcess(A[t](),tA[t](),invA[t](),invtA[t](),a[t](),B[t](),tB[t]());
        hatq_old[t]() = hatq[t];
	//cout <<"r\n" <<r[t] <<endl <<R[t] <<endl;
      }
      //else cout <<"skip." <<flush;
          
      //compute (b,B);
      invB[t] = invS[t] + invV[t] + R[t];
      lapack_Ainv_b_sym(b[t](), invB[t], invS[t]*s[t] + invV[t]*v[t] + r[t]);
      //cout <<"b\n" <<b[t] <<endl <<B[t] <<endl;

#ifndef TightMode
      //decide on \hat q
      if(sweep || !useFwdMessageAsInitialHatq){
	if(convergenceRate) hatq[t]()=((real)1.-convergenceRate)*hatq[t] + convergenceRate*b[t];
	else hatq[t]()=b[t];
      }
#else
      //update hatq
      if(dt==1){
	if(!convergenceRate){
	  hatq[t]()=b[t];
	}else{
	  hatq.reshape(T+1,2,n);
	  hatq.subDim(t,1)
	    = ((real)1.-convergenceRate)*hatq.subDim(t,1) + convergenceRate*(b[t].sub(n,-1));
	  hatq.subDim(t,0) = hatq.subDim(t-1,0) + sys->getTau()*hatq.subDim(t,1);
	  hatq.reshape(T+1,2*n);
	}
      }
#endif
        
#if 0 //don't check for task discounting
      arr invD,d;
      real like;
      do{
        //cout <<"b\n" <<b[t] <<endl <<invB[t] <<endl;
	
        //compute likelihood
	invD = invS[t] + invV[t];
	lapack_Ainv_b_sym(d, invD, invS[t]*s[t] + invV[t]*v[t]);
	like=metricDistance(invD,d,b[t]);
        /*cout <<t
          <<" - fwd like=" <<metricDistance(invS[t],s[t],b[t])
          <<" - bwd like=" <<metricDistance(invV[t],v[t],b[t])
          <<" - r like=" <<like
          <<endl;*/
	if(like>100.){
	  cout <<" REDUCING task precision" <<endl;
	  R[t]() *= .1;
	  r[t]() *= .1;
	}
	//compute (b,B);
	invB[t] = invS[t] + invV[t] + R[t];
	lapack_Ainv_b_sym(b[t](), invB[t], invS[t]*s[t] + invV[t]*v[t] + r[t]);
      }while(like>100.);
#endif
      
      //decide whether to repeat this time slice
      if(sweep && repeatThreshold && t!=T){
	//real off=sqrDistance(Q,b[t],hatq[t]);
	real off=sqrDistance(b[t],hatq[t]);
	if(off>repeatThreshold){
	  //cout <<t <<" REPEAT: off=" <<off <<" (repeatCount=" <<repeatCount <<")" <<endl;
	  if(repeatCount<20){
	    t-=dt;
	    repeatCount++;
	  }else{
	    cout <<" ** no convergence! ** (skipping repeat) at t=" <<t <<endl;
	    repeatCount=0;
	  }
	}else{
	  repeatCount=0;
	}
      }
    } //loop t
    sweep++;
#ifdef TightMode
    b=hatq;
#endif
  }//loop over dt in {-1,1}
  
#if 1
  getPositionTrajectory(q,b);
#else
  getControlledTrajectory(q,*this);
#endif
  if(q_old.N==q.N) diff=maxDiff(q_old,q);
  
  //display or evaluate
  MT::timerPause();
  if(sys->os){
    *sys->os <<"AICOd("<<scale<<") " <<std::setw(3) <<sweep <<" time " <<MT::timerRead(false) <<" diff " <<diff;
    getPositionTrajectory(q,b);
    cost = sys->analyzeTrajectory(q,display>0);
  }
  if(sys->gl){
    if(!sys->os) getPositionTrajectory(q,b);
    sys->displayTrajectory(q,NULL,display,STRING("AICO_dynamic - iteration "<<sweep));
  }
  MT::timerResume();

  return diff;
}

