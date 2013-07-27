// ProbabilityCollision.cpp : Defines the entry point for the console application.
//

#define M_1_SQRT2PI 0.39894228040143267793994605993438

#include "stdafx.h"
#include "igammaf.h"
#include "matrix.h"

// uniform random sample in [0,1]
inline double uniform() {
  return (double) rand() / RAND_MAX;
}

// random sample from N(0,1)
inline double normal() {
  double u, v, s(0);

  while (s == 0 || s >= 1) {
    u = 2.0*uniform()-1.0;
    v = 2.0*uniform()-1.0;
    s = u*u + v*v;
  }
  return u * sqrt(-2.0*log(s)/s);
}

// random sample from N(mean,stddev^2)
inline double normal(double mean, double stddev) {
  return stddev * normal() + mean;
}

template <size_t dim>
inline Matrix<dim> sampleGaussian() {
  Matrix<dim> sample;
  for (int j = 0; j < dim; ++j) {
    sample[j] = normal();
  }
  return sample;
}

template <size_t dim>
inline Matrix<dim> sampleGaussian(const Matrix<dim>& mean, const Matrix<dim, dim>& var) {
  Matrix<dim> sample = sampleGaussian<dim>();
  Matrix<dim, dim> SVec, SVal;
  jacobi(var, SVec, SVal);
  for (int i = 0; i < dim; ++i) {
    if (SVal(i,i) < 0) {
      SVal(i,i) = 0;
    } else {
      SVal(i,i) = sqrt(SVal(i,i));
    }
  }
  return SVec * SVal * sample + mean;
}

inline double pdf(double x) {
  return exp(-0.5*x*x) * M_1_SQRT2PI;
}

inline double cdf(double x) {
  if (x < 0) {
    return 0.5 - 0.5*incompletegamma(0.5, 0.5*x*x);
  } else {
    return 0.5 + 0.5*incompletegamma(0.5, 0.5*x*x);
  }
}

inline void truncate(double x, double mean, double var, double& newMean, double& newVar) {
  double stddev = sqrt(var);
  double y = (x - mean) / stddev;
  double z = pdf(y) / cdf(y);
  newMean = mean - z*stddev;
  newVar = var*(1.0 - y*z - z*z);
}


int _tmain(int argc, _TCHAR* argv[])
{
  
  time_t t;
  srand(time(&t));

  
  Matrix<2,2> xVar;
  Matrix<2> xMean;
  Matrix<1,2> C1;
  Matrix<1,2> C2;
  double c1, c2;
  
  xMean[0] = 2.0;
  xMean[1] = 1.0;

  xVar(0,0) = 2.0; xVar(0,1) = 0.5;
  xVar(1,0) = 0.5; xVar(1,1) = 1.0;

  //C1(0,0) = -3.0; C1(0,1) = 2.0;
  C1(0,0) = -1.0; C1(0,1) = 0.0;
  c1 = 0.0;

  C2(0,0) = 1.0; C2(0,1) = 0.0;
  c2 = 3.0;

  double y1Mean = tr(C1*xMean);
  double y1Var = tr(C1*xVar*~C1);

  double y1NewMean, y1NewVar;
  truncate(c1, y1Mean, y1Var, y1NewMean, y1NewVar);
  
  double y2Mean = tr(C2*xMean);
  double y2Var = tr(C2*xVar*~C2);

  double y2NewMean, y2NewVar;
  truncate(c2, y2Mean, y2Var, y2NewMean, y2NewVar);

  Matrix<2,1> xy1Covar = xVar*~C1;
  Matrix<2,1> xy2Covar = xVar*~C2;

  Matrix<2> xNewMean;
  Matrix<2,2> xNewVar;

  Matrix<2,1> L1, L2;
  L1 = xy1Covar / y1Var; L2 = xy2Covar / y2Var;
  xNewMean = xMean + L1*(y1NewMean - y1Mean) + L2*(y2NewMean - y2Mean);
  xNewVar = xVar + L1*(y1NewVar - y1Var)*~L1 + L2*(y2NewVar - y2Var)*~L2;

  /*Matrix<2,1> L;

  Matrix<2> xNewMean1, xNewMean2;
  Matrix<2,2> xNewVar1, xNewVar2;

  L = xy1Covar / y1Var;
  xNewMean1 = xMean + L*(y1NewMean - y1Mean);
  xNewVar1 = xVar + L*(y1NewVar - y1Var)*~L;

  L = xy2Covar / y2Var;
  xNewMean2 = xMean + L*(y2NewMean - y2Mean);
  xNewVar2 = xVar + L*(y2NewVar - y2Var)*~L;

  xNewMean = (!xNewVar1 + !xNewVar2)%(!xNewVar1*xNewMean1 + !xNewVar2*xNewMean2);
  xNewVar = !(!xNewVar1 + !xNewVar2);*/

  std::vector<Matrix<2> > samples;

  int numSamples = 10000000;
  for (int i = 0; i < numSamples; ++i) {
    samples.push_back(sampleGaussian(xMean, xVar));
  }

  Matrix<2> sampleMean = zeros<2>();
  Matrix<2,2> sampleVar = zeros<2,2>();
  int numTruncSamples = 0;
  for (int i = 0; i < numSamples; ++i) {
    if (tr(C1*samples[i]) <= c1 && tr(C2*samples[i]) <= c2) {
      sampleMean += samples[i];
      ++numTruncSamples;
    }
  }
  sampleMean /= numTruncSamples;
  for (int i = 0; i < numSamples; ++i) {
    if (tr(C1*samples[i]) <= c1 && tr(C2*samples[i]) <= c2) {
      sampleVar += (samples[i] - sampleMean)*~(samples[i] - sampleMean);
    }
  }
  sampleVar /= (numTruncSamples - 1);
  

  

  std::cout 
    << "Mean x: " << std::endl << ~xNewMean << std::endl
    << "Var x:  " << std::endl << xNewVar  << std::endl
    << "Sample Mean x: " << std::endl << ~sampleMean << std::endl
    << "Sample Var x:  " << std::endl << sampleVar  << std::endl;

  int k;
  std::cin >> k;
	return 0;
}

