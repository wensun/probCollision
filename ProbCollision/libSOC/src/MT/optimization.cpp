#include "optimization.h"

#define CHECK_EPS 1e-8

void checkGradient(OptimizationProblem &p,
                   const arr& x, real tolerance)
{
  arr J,dx,JJ;
  real y,dy;
  y=p.f(&J,x);
  
  JJ.resize(x.N);
  real eps=CHECK_EPS;
  uint i;
  for(i=0;i<x.N;i++){
    dx=x;
    dx.elem(i) += eps;
    dy = p.f(NULL,dx);
    dy = (dy-y)/eps;
    JJ(i)=dy;
  }
  JJ.reshapeAs(J);
  real md=maxDiff(J,JJ,0);
//   MT::save(J,"z.J");
//   MT::save(JJ,"z.JJ");
  if(md>tolerance){
    MT_MSG("checkGradient -- FAILURE -- max diff=" <<md <<" (stored in files z.J and z.JJ)");
    MT::save(J,"z.J");
    MT::save(JJ,"z.JJ");
    //cout <<"\nmeasured grad=" <<JJ <<"\ncomputed grad=" <<J <<endl;
    //HALT("");
  }else{
    cout <<"checkGradient -- SUCCESS (max diff error=" <<md <<")" <<endl;
  }
}

void checkGradient_vec(OptimizationProblem &p,
                       const arr& x,real tolerance)
{
  arr y,J,dx,dy,JJ;
  p.F(y,&J,x);
  
  JJ.resize(y.N,x.N);
  real eps=CHECK_EPS;
  uint i,k;
  for(i=0;i<x.N;i++){
    dx=x;
    dx.elem(i) += eps;
    p.F(dy,NULL,dx);
    dy = (dy-y)/eps;
    for(k=0;k<y.N;k++) JJ(k,i)=dy.elem(k);
  }
  JJ.reshapeAs(J);
  real md=maxDiff(J,JJ,&i);
//   MT::save(J,"z.J");
//   MT::save(JJ,"z.JJ");
  if(md>tolerance){
    MT_MSG("checkGradient -- FAILURE -- max diff=" <<md <<" (stored in files z.J and z.JJ)");
    MT::save(J,"z.J");
    MT::save(JJ,"z.JJ");
    //cout <<"\nmeasured grad=" <<JJ <<"\ncomputed grad=" <<J <<endl;
    //HALT("");
  }else{
    cout <<"checkGradient -- SUCCESS (max diff error=" <<md <<")" <<endl;
  }
}


//===========================================================================
//
// Rprop
//

int _sgn(real x){ if (x > 0) return 1; if (x < 0) return -1; return 0; }
real _mymin(real x,real y){ return x < y ? x : y; }
real _mymax(real x,real y){ return x > y ? x : y; }


Rprop::Rprop(){
  incr   = 1.2;
  decr   = .33;
  dMax = 50;
  dMin = 1e-6;
  rMax = 0;
  delta0 = 1.;
}

void Rprop::init(real _delta0){
  stepSize.resize(0);
  lastGrad.resize(0);
  delta0 = _delta0;
}

bool Rprop::done(){
  real maxStep = stepSize(stepSize.maxIndex());
  return maxStep < incr*dMin;
}

void Rprop::step(real& w,const real& grad){
  static arr W,GRAD;
  W.referTo(&w,1); GRAD.referTo(&grad,1);
  step(W,GRAD);
}

void Rprop::step(arr& w,const arr& grad,uint *singleI){
  if(!stepSize.N){ //initialize
    stepSize.resize(w.N);
    lastGrad.resize(w.N);
    lastGrad.setZero();
    stepSize = delta0;
  }
  CHECK(grad.N==stepSize.N,"Rprop: gradient dimensionality changed!");
  CHECK(w.N==stepSize.N   ,"Rprop: parameter dimensionality changed!");

  uint i=0,I=w.N;
  if(singleI){ i=*(singleI); I=i+1; }
  for(;i<I;i++){
    if(grad.elem(i) * lastGrad(i) > 0){        //same direction as last time
      if(rMax) dMax=fabs(rMax*w.elem(i));
      stepSize(i) = _mymin(dMax, incr * stepSize(i)); //increase step size
      w.elem(i) += stepSize(i) * -_sgn(grad.elem(i)); //step in right direction
      lastGrad(i) = grad.elem(i);                    //memorize gradient
    }else if(grad.elem(i) * lastGrad(i) < 0){  //change of direction
      stepSize(i) = _mymax(dMin, decr * stepSize(i)); //decrease step size
      w.elem(i) += stepSize(i) * -_sgn(grad.elem(i)); //step in right direction
      lastGrad(i) = 0;                               //memorize to continue below next time
    }else{                                     //after change of direcion
      w.elem(i) += stepSize(i) * -_sgn(grad.elem(i)); //step in right direction
      lastGrad(i) = grad.elem(i);                    //memorize gradient
    }
  }
}

void Rprop::step(arr& x,OptimizationProblem& p){
   arr grad;
   p.f(&grad,x);
   step(x,grad);
}

//----- the rprop wrapped with stopping criteria
int Rprop::loop(arr& _x,
                OptimizationProblem& p,
                real *fmin_return,
                real stoppingTolerance,
                uint maxIterations)
{
  arr x,J(_x.N),xmin;
  real y,ymin=0;
  uint lost_steps=0,small_steps=0;
  x=_x;
  
  uint i;
  for(i=0;i<maxIterations;i++){
    cout <<"RPROP iter= "<<i <<flush;
    //checkGradient(p,x,stoppingTolerance);
    //compute value and gradient at x
    y = p.f(&J,x);
    //update best-so-far
    if(!i){ ymin=y; xmin=x; }
    if(y<ymin){
      ymin=y; xmin=x;
      lost_steps=0;
    }else{
      lost_steps++;
      if(lost_steps>10){
        stepSize*=(real).1;
        lastGrad=(real)0.;
        x=xmin;
        lost_steps=0;
      }
    }
    //update x
    step(x,J);
    //check stopping criterion based on step-length in x
    double diff=maxDiff(x,xmin);
    cout <<"  x-diff= " <<diff <<"  f= " <<y <<endl;
    if(diff<stoppingTolerance){ small_steps++; }else{ small_steps=0; }
    if(small_steps>10)  break;
  }
  if(fmin_return) *fmin_return=ymin;
  _x=xmin;
  return i;
}


#ifdef MT_GSL
#include <gsl/gsl_cdf.h>
bool DecideSign::step(double x)
{
  N++;
  sumX+=x;
  sumXX+=x*x;
  if(N<=10) return false;
  if(!sumX) return true;
  double m=sumX/N;
  double s=sqrt((sumXX-sumX*m)/(N-1));
  double t=sqrt(N)*fabs(m)/s;
  double T=gsl_cdf_tdist_Pinv(1.-1e-6,N-1); //decide with error-prob 1e-4 if sign is significant
  if(t>T) return true;
  return false;
}
#else
             bool DecideSign::step(double x){ NIY; }
#endif
