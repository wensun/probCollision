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

/** \file soc.h
    \brief Stochastic Optimal Control library */

#ifndef MT_soc_h
#define MT_soc_h

#include "array.h"

//-- fwd declarations
class OpenGL;
struct SwiftModule;
namespace ors{ struct Graph; }
struct TaskVariable;
typedef MT::Array<TaskVariable*> TaskVariableList;


namespace soc{

/// gradient method options
enum { ConjGrad=0, LevMar=1, Rprop=2, RpropConjGrad=3, SQP=4, Attractor=5 };


//===========================================================================
//
// SocSystemAbstraction
//

/** \brief defines an abstraction of stochastic optimal control
    problems which interfaces between solution methods and system simulators
    -- see section 3.2 of the <a href="../guide.pdf">guide</a> */
struct SocSystemAbstraction{

  ///@name data fields
  std::ostream *os; ///< if non-NULL, some routines might give output
  OpenGL *gl;       ///< if non-NULL, some routines might give output
  bool dynamic;     ///< determines whether this problem is dynamic or not
  uint scalePower;  ///< if non-zero, all routines assume an horizon T=T/2^scalePower
  
  ///@name initialization
  SocSystemAbstraction();
  virtual ~SocSystemAbstraction();
  //virtual void clone(const SocSystemAbstraction& sys); ///< deeply clones sys, i.e., creates copies of simulators etc
    
  ///@name low level access routines: need to be implemented by the simulator

  // general problem information
  virtual uint nTime() = 0;            ///< total time of the trajectory
  virtual uint nTasks() = 0;           ///< number of task variables
  virtual uint qDim() = 0;             ///< dimensionality of q-space
  virtual uint uDim();                 ///< dimensionality of control
  virtual uint yDim(uint i) = 0;       ///< dimensionality of the i-th task
  virtual void getq0 (arr& q) = 0;     ///< start joint configuration
  virtual void getv0 (arr& v) = 0;     ///< start joint velocity
  virtual void getqv0(arr& q_);        ///< start joint configuration and velocity
  virtual void getqv0(arr& q,arr& qd); ///< start joint configuration and velocity
  virtual real getTau(bool scaled=true);    ///< time step size (for dynamic problems)
  void getx0(arr& x){ if(dynamic) getqv0(x); else getq0(x); }

  // set x-state (following calls to getPhi and getJ are w.r.t. this x)
  virtual void setq  (const arr& q) = 0;
  virtual void setqv (const arr& q_);
  virtual void setqv (const arr& q,const arr& qd);
  void setx(const arr& x){ if(dynamic) setqv(x); else setq(x); }
  virtual void setq0AsCurrent() = 0;

  //motion prior, or control cost
  //virtual void geth  (arr& h) = 0; NIY...
  virtual void getW   (arr& W) = 0;     ///< kinematic step cost metric: cost = dq^T W dq
  virtual void getWinv(arr& Winv){ throw("NIY"); } ///< kinematic step cost metric: cost = dq^T W dq
  virtual void getH   (arr& H);         ///< dynamic control cost metric: cost = u^T H u
  virtual void getQ   (arr& Q);         ///< process stochasticity or integration noise Q (e.g., setDiag(1e-10,qDim()) )

  // task coupling
  virtual bool isConditioned(uint i,uint t) = 0;
  virtual bool isConstrained(uint i,uint t);
  virtual const char* taskName(uint i){ return NULL; };
  virtual void getPhi(arr& phiq_i,uint i){ throw("NIY"); }
  virtual void getJtJ(arr& J_i,arr& tJ_i,uint i){ throw("NIY"); }
  virtual void getJqd(arr& Jqd_i,uint i);
  virtual void getHessian(arr& H_i,uint i);
  virtual void getTarget (arr& y_i,real& prec,uint i,uint t){ throw("NIY"); }
  virtual void getTargetV(arr& v_i,real& prec,uint i,uint t);
  //virtual void getLinearConstraint(arr& c,real& coff,uint i,uint t); ///< defines a cost 1 iff [[c^T y_i + coff>0]]

  // dynamic model
  virtual void getMF(arr& M,arr& F);
  virtual void getMinvF(arr& Minv,arr& F);


  ///@name high level methods: they are being accessed by the solvers

  // abstract SOC interface
  virtual void getProcess(arr& A,arr& a,arr& B);
  virtual void getProcess(arr& A,arr& tA,arr& invA,arr& invtA,arr& a,arr& B,arr& tB);
  virtual real getCosts(arr& R,arr& r,uint t,const arr& qt);
  virtual void getConstraints(arr& c,arr& coff,uint t,const arr& qt);
      
  // cost info
  real computeTotalCost(const arr& q,bool plot=false);
  real computeTotalCostGradient(arr& dCdq,const arr& q);
  
  virtual void displayState(const arr& q,const arr *invQ,const char *text);
  virtual void displayTrajectory(const arr& q,const arr *invQ,int steps,const char *tag);

  //-- convenience (prelim...)
  real analyzeTrajectory(const arr& q,bool plot);
  void constantTrajectory(arr& q);
  void passiveDynamicTrajectory(arr& q);
  void controlledDynamicTrajectory(arr& q,const arr& u);
  void getControlFromTrajectory(arr& u,const arr& q);
};


//===========================================================================
//
///@name     trivial helpers
// @{

void getVelocity(arr& vt,const arr& q,uint t,real tau);
void getPhaseTrajectory(arr& _q,const arr& q,real tau);
void getPositionTrajectory(arr& q,const arr& _q);
void interpolateTrajectory(arr& qNew,const arr& _q,real step);

//only for the first task so far!
void getJointFromTaskTrajectory(SocSystemAbstraction& soci,arr& q,const arr& x);
void partialJointFromTaskTrajectory(SocSystemAbstraction& soci,arr& dx,const arr& dq,const arr& q,const arr& x);
void straightTaskTrajectory(SocSystemAbstraction& soci,arr& q,uint taskid);


//===========================================================================
// @}
///@name     inverse kinematics control
// @{

void bayesianIKControl(SocSystemAbstraction& soci, arr& dq, uint t);
void pseudoIKControl(SocSystemAbstraction& soci, arr& dq, uint t,real regularization=1e-8);
void hierarchicalIKControl(SocSystemAbstraction& soci, arr& dq, uint t,real regularization=1e-8);
void bayesianIterateIKControl(SocSystemAbstraction& soci,
                              arr& qt,const arr& qt_1,uint t,real eps,uint maxIter);
void bayesianIKTrajectory  (SocSystemAbstraction& soci, arr& q, real eps=-1);
void bayesianDynamicControl(SocSystemAbstraction& soci, arr& qv, const arr& qv_1, uint t, arr *v=NULL, arr *invV=NULL);
void bayesianIKControl2    (SocSystemAbstraction& soci, arr& q , const arr& q_1 , uint t, arr *v=NULL, arr *invV=NULL);



//===========================================================================
// @}
///@name     gradient optimization
// @{

void gradientOptimization(SocSystemAbstraction& soci,
                                arr& q,
//                                int gradient_method,
                                uint maxIterations,
                                uint spline_points,
                                uint spline_degree,
                                real stoppingTolerance,
                                bool checkGradient,
                                uint display);

void gradientTaskOptimization(SocSystemAbstraction& soci,
                                arr& q,
                                uint spline_points,
                                uint spline_degree,
                                int gradient_method,
                                uint maxIterations,
                                real stoppingTolerance,
                                bool checkGradient,
                                uint display);

void gradientAttractorTaskOptimization(SocSystemAbstraction& soci,
                                arr& q,
                                uint spline_points,
                                uint spline_degree,
                                int gradient_method,
                                uint maxIterations,
                                real stoppingTolerance,
                                bool checkGradient,
                                uint display);

void SQPOptimization(SocSystemAbstraction& soci,
                           arr& q,uint iterations,
                           uint spline_points,
                           uint spline_degree,
                           uint display);

//===========================================================================
// @}
///@name     iLQG methods
// @{

struct iLQG{
  //parameters
  SocSystemAbstraction *sys;
  real convergenceRate;
  uint display;
  //messages & state info
  arr q,_q;
  arr v,V,r,R;      //!< fwd, bwd, and task messages
  real cost;                      //!< cost of MAP trajectory
  uint sweep;                     //!< #sweeps so far
  uint scale;                     //!< scale of this iLQG in a multi-scale approach
  
  iLQG(){ sweep=0; scale=0; }
  
  void init(SocSystemAbstraction& _sys,
            double _convergenceRate, uint _display,
            uint _scale)
  {
    sys = &_sys;  //sys.clone(_sys);
    convergenceRate=_convergenceRate;
    display=_display;
    scale=_scale;
  }

  real stepDynamic  ();
  real stepKinematic();
};

inline iLQG* iLQG_solve(SocSystemAbstraction& sys,
                       arr& q,double tolerance,
                       real convergenceRate,
                       uint display)
{
  iLQG *ilqg=new iLQG;
  ilqg->init(sys,convergenceRate,display,0);
  ilqg->q=q;
  if(!sys.dynamic){ for(uint k=0;k<100;k++) if(ilqg->stepKinematic()<tolerance) break; }
  else            { for(uint k=0;k<100;k++) if(ilqg->stepDynamic  ()<tolerance) break; }
  q=ilqg->q;
  return ilqg;
}

inline iLQG* iLQG_multiScaleSolver(SocSystemAbstraction& sys,
                                  arr& q,
                                  double tolerance,
                                  real convergenceRate,
                                  uint display,
                                  uint scalePowers)
{
  iLQG *ilqg=new iLQG;
  ilqg->init(sys,convergenceRate,display,scalePowers-1);
  ilqg->q=q;
  for(uint i=scalePowers;i--;) for(int k=0;;k++){
    sys.scalePower=i;
    double d=ilqg->stepKinematic();
    if(k && d<tolerance) break;
  }
  q=ilqg->q;
  return ilqg;
}

//===========================================================================
// @}
///@name     bayesian inference methods
// @{

/** \brief Apprioximate Inference Control */
struct AICO{
  //parameters
  SocSystemAbstraction *sys;
  real convergenceRate,repeatThreshold,recomputeTaskThreshold;
  uint display;

  //messages & state info
  arr s,invS,v,invV,r,R;      //!< fwd, bwd, and task messages
  arr b,invB;                     //!< beliefs
  arr q,hatq;                     //!< trajectory (MAP), and point of linearization
  arr A,tA,invA,invtA,a,B,tB; //!< processes...
  real cost;                      //!< cost of MAP trajectory
  uint sweep;                     //!< #sweeps so far
  
  uint scale;                     //!< scale of this AICO in a multi-scale approach
  //MT::Array<Lock> locks;
  MT::Array<AICO*> multiScales;   //!< list of all AICOs in a multi-scale approach
  AICO(){ sweep=0; scale=0; }
  
  void init(SocSystemAbstraction& _sys,
            double _convergenceRate,double _repeatThreshold, double _recomputeTaskThreshold, uint _display,
            uint _scale,const MT::Array<soc::AICO*>& _multiScales){
    sys = &_sys;  //sys.clone(_sys);
    convergenceRate=_convergenceRate;
    repeatThreshold=_repeatThreshold;
    recomputeTaskThreshold=_recomputeTaskThreshold;
    display=_display;
    scale=_scale;
    multiScales=_multiScales;
    initMessages();
  }

  void initMessages();
  real stepDynamic  ();
  real stepKinematic();
  real stepDynamicIlqg();

  //--internal use
  void getMultiScaleMessages(arr& s_,arr& S_,arr& v_,arr& V_,uint t,real upMixS,real selfMixS,real dnMixS,real upMixV,real selfMixV,real dnMixV);
};

soc::AICO* AICO_solver(SocSystemAbstraction& soci,
                       arr& q,double tolerance,
                       real convergenceRate,real repeatThreshold, real recomputeTaskThreshold,
                       uint display);

soc::AICO* AICO_multiScaleSolver(SocSystemAbstraction& sys,
                           arr& q,
                           double tolerance,
                           real convergenceRate,real repeatThreshold, real recomputeTaskThreshold,
                           uint display,
                           uint scalePowers);

void BayesianDynamicMotionPlanning_obsolete(SocSystemAbstraction& soci,
                                     arr& q,uint iterations,
                                     real stepSmooth,real repeatThreshold,
                                     uint display);

void BayesianDynamicFilter_obsolete(SocSystemAbstraction& soci,
				    arr& q,uint iterations,
				    real stepSmooth,real repeatThreshold,
				    uint display);
                                    
inline void getController(arr& G,arr& g,const soc::AICO& aico){
  //we can only compute a controller for time steps 0 to T-1 (based on V_{t+1})
  uint T=aico.s.d0-1;
  uint n=aico.s.d1;
  if(!aico.sys->dynamic){
    G.resize(T,n,n);
    g.resize(T,n);
  }else{
    G.resize(T,n/2,n);
    g.resize(T,n/2);
  }
  arr H;
  aico.sys->getH(H);
  for(uint t=0;t<T;t++){
    arr Vstar,barv,VstarH;
    if(!aico.sys->dynamic){
      //controller model u_mean = G*x+g
      Vstar = aico.invV[t+1] + aico.R[t+1];
      lapack_Ainv_b_sym(barv, Vstar, aico.invV[t+1]*aico.v[t+1] + aico.r[t+1]);
      inverse_SymPosDef(VstarH,Vstar + H);
      G[t] = - VstarH * Vstar; // * aico.A[t];
      g[t] = VstarH * Vstar * (barv); // - aico.a[t]);
    }else{
      Vstar = aico.invV[t+1] + aico.R[t+1];
      lapack_Ainv_b_sym(barv, Vstar, aico.invV[t+1]*aico.v[t+1] + aico.r[t+1]);
      inverse_SymPosDef(VstarH,aico.tB[t]*Vstar*aico.B[t] + H);
      G[t] = - VstarH * aico.tB[t] * Vstar * aico.A[t];
      g[t] = VstarH * aico.tB[t] * Vstar * (barv - aico.a[t]);
    }
  }
}

inline void forwardSimulateTrajectory(arr& q,const arr& G,const arr& g,soc::SocSystemAbstraction& sys,const soc::AICO& aico){
  uint t,T=sys.nTime(),n=sys.qDim();
  if(!aico.sys->dynamic){
    q.resize(T+1,n);
    sys.getq0(q[0]());
    for(t=0;t<T;t++) q[t+1]() = q[t] + (G[t]*q[t] + g[t]); //A=1, B=1
  }else{
    q.resize(T+1,2*n);
    sys.getqv0(q[0]());
    for(t=0;t<T;t++) q[t+1]() = aico.A[t]*q[t] + aico.B[t]*(G[t]*q[t] + g[t]) + aico.a[t];
    arr q_sub;
    getPositionTrajectory(q_sub,q);
    q=q_sub;
  }
}

inline void getControlledTrajectory(arr& q,const soc::AICO& aico){
  arr G,g;
  getController(G,g,aico);
  forwardSimulateTrajectory(q,G,g,*aico.sys,aico);
}

//===========================================================================
// @}
///@name preliminary or obsolete
// @{

#if 1
class SocSystem_Ors;
class SocSystem_Toy;

real getDynamicCostMeassure(SocSystemAbstraction& soci,arr& q,real& cost1,real& cost2,std::ostream *os=0);

real getFilterCostMeassure(SocSystemAbstraction& soci,arr& q,real& cost1,real& cost2,std::ostream *os=0);

real getTaskCost(soc::SocSystemAbstraction& soci,int t,int i=-1);
void getTaskCostGradient(soc::SocSystemAbstraction& soci,arr& dCdq,int t);
void getQuadraticTaskCost(soc::SocSystemAbstraction& soci,int t,arr& R,arr& r,const arr& qt);

void setupOpenGL(SocSystem_Ors &soci);

/*void createNikolayReachProblem(SocSystem_Ors &soci,
                               slGraph &ors,
                               SwiftModule& swift,
                               uint trajectory_length,
                               const arr& endeffector_target,
                               const char* endeffector_name,
                               const arr& W);*/

void createDynamicProblem(SocSystem_Ors &soci,
                          const char *ors_file,
                          real trajectory_time,
                          uint trajectory_steps);

void setupOpenGL(SocSystem_Toy &soci);

void createEndeffectorReachProblem(SocSystem_Toy &soci,
                                   const char *ors_file,
                                   uint trajectory_length,
                                   int rand_seed);

/*void createNikolayReachProblem(SocSystem_Toy &soci,
                               slGraph &ors,
                               SwiftModule& swift,
                               uint trajectory_length,
                               const arr& endeffector_target,
                               const char* endeffector_name,
                               const arr& W);*/

void createDynamicProblem(SocSystem_Toy &soci,
                          const char *ors_file,
                          real trajectory_time,
                          uint trajectory_steps);
#endif

//===========================================================================
// @}
// ORS simulator implementation of the SocAbstration
//

struct SocSystem_Ors_Workspace;

/** \brief an implementation of the SocSystemAbstraction using the \ref ors
    simulator */
struct SocSystem_Ors: public virtual SocSystemAbstraction{
  ors::Graph *ors;
  SwiftModule *swift;
  MT::Array<TaskVariable*> vars;
  SocSystem_Ors_Workspace *WS;

  SocSystem_Ors();
  virtual ~SocSystem_Ors();
  void clone(const SocSystem_Ors& _sys);
  
  //initialization methods
  void initKinematic(ors::Graph *ors,SwiftModule *swift,OpenGL *gl,uint trajectory_length, arr*W=NULL);
  void initPseudoDynamic(ors::Graph *ors,SwiftModule *swift,OpenGL *gl,
                         real trajectory_time,uint trajectory_steps, arr*W=NULL);
  void setTimeInterval(real trajectory_time,uint trajectory_steps);
  void setTaskVariables(const TaskVariableList& CVlist);

  //--exemplary problem setups: read specifications from MT.cfg
  void initStandardReachProblem(uint rand_seed=0,uint T=0);
  void initStandardBenchmark(uint rand_seed=0);
  
  //info
  void reportOnState(std::ostream& os);
  void displayState(const arr& q,const arr *Q,const char *text=NULL);

  //implementations of virtual methods
  uint nTime();
  uint nTasks();
  uint qDim();
  uint uDim();
  uint yDim(uint i);
  void getq0 (arr& q);
  void getv0 (arr& v);
  void getqv0(arr& q_);
  void getqv0(arr& q,arr& qd);
  bool isDynamic();
  void setq  (const arr& q);
  void setqv (const arr& q_);
  void setqv (const arr& q,const arr& qd);
  void setq0AsCurrent();
  void geth  (arr& h);
  void getW  (arr& W);
  void getWinv(arr& Winv);
  void getH  (arr& H);
  void getQ  (arr& Q);
  bool isConditioned(uint i,uint t);
  bool isConstrained(uint i,uint t);
  const char* taskName(uint i);
  void getPhi(arr& phiq_i,uint i);
  void getJtJ(arr& J_i,arr& tJ_i,uint i);
  void getHessian(arr& H_i,uint i);
  void getJqd(arr& Jqd_i,uint i);
  void getTarget (arr& y_i,real& prec,uint i,uint t);
  void getTargetV(arr& v_i,real& prec,uint i,uint t);
  //void getC   (arr& C_i,uint i,uint t);
  //void getCV  (arr& D_i,uint i,uint t);
  real getTau(bool scaled=true);
  void getMF(arr& M,arr& F);
  void getMinvF(arr& Minv,arr& F);
};

//========= untidy stuff ...



//===========================================================================
//
// toy implementation of the SocAbstration
//

struct SocSystem_Toy_Workspace;

/** \brief an implementation of the SocSystemAbstraction that simulates a
    single 1D point mass on a spring */
struct SocSystem_Toy: public virtual SocSystemAbstraction{
  SocSystem_Toy_Workspace *WS;

  SocSystem_Toy();
  virtual ~SocSystem_Toy();

  //implementations of virtual methods
  uint nTime();
  uint nTasks();
  uint qDim();
  uint uDim();
  uint yDim(uint i);
  void getq0 (arr& q);
  void getv0 (arr& v){throw("NIY");}
  void getqv0(arr& q_);
  void getqv0(arr& q,arr& qd);
  bool isDynamic();
  void setq  (const arr& q);
  void setqv (const arr& q_);
  void setqv (const arr& q,const arr& qd);
  void setq0AsCurrent();
  void geth  (arr& h);
  void getW  (arr& W);
  void getWinv(arr& Winv){ throw("NIY"); };
  void getH  (arr& H);
  void getQ  (arr& Q);
  bool isConditioned(uint i,uint t);
  void getPhi(arr& phiq_i,uint i);
  void getJtJ(arr& J_i,arr& tJ_i,uint i);
  void getJqd(arr& Jqd_i,uint i);
  void getTarget (arr& y_i,uint i,uint t);
  void getTargetV(arr& v_i,uint i,uint t);
  void getC   (arr& C_i,uint i,uint t);
  void getCV  (arr& D_i,uint i,uint t);
  void getPrecision (real& prec,uint i,uint t);
  void getPrecisionV(real& prec,uint i,uint t);
  real getTau(bool scaled=true);
  void getMF(arr& M,arr& F);
  void getMinvF(arr& Minv,arr& F);
};

}


//===========================================================================
//
// implementations
//

extern uint countMsg,countSetq;

#ifdef MT_IMPLEMENTATION
#  include "soc.cpp"
#  include "soc_method_AICO.cpp"
#  include "soc_method_LQG.cpp"
#  include "soc_method_gradient.cpp"
//#  include "soc_method_attractor.cpp"
#  include "soc_system_ors.cpp"
#  include "soc_system_analytical.cpp"
#  include "soc_system_toy.cpp"
#endif

#endif
