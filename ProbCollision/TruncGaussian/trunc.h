#ifndef __TRUNC_GAUSSIAN_H__
#define __TRUNC_GAUSSIAN_H__

#define _CRT_RAND_S
#include <cstdlib>

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <algorithm>

#include "matrix.h"

#include <float.h>

#include "callisto.h"
#include "glut.h"
#include "gammafunc.h"

// constant definitions
#define M_SQRT_PI_2 sqrt(M_PI/2.0)
#define M_1_SQRT_2PI 1.0/sqrt(M_PI*2.0)
#define M_DEG2RAD M_PI/180.0

// var definitions
#define SDIM 2

// Cumulative probability for the Standard Normal Distribution
inline double erf(double x)
{
	double t, z, retval;
	z = fabs( x );
	t = 1.0 / ( 1.0 + 0.5 * z );
	retval = t * exp( -z * z - 1.26551223 + t *
		( 1.00002368 + t *
		( 0.37409196 + t *
		( 0.09678418 + t *
		( -0.18628806 + t *
		( 0.27886807 + t *
		( -1.13520398 + t *
		( 1.48851587 + t *
		( -0.82215223 + t *
		0.1708727 ) ) ) ) ) ) ) ) );
	if( x < 0.0 ) return retval - 1.0;
	return 1.0 - retval;
}
  
// Integral of N(0,1) from -infty to x
inline double gaussInt(double x) {
	return 0.5 * (1.0 + erf(x*M_SQRT1_2));
}

#define AXETS 1280
#define AXETR 10.0
// Approximate exp (sets up a static value table)
inline double approxExp(double x)
{
	static bool initialized=false;
	static double ExpTable[AXETS]; //table ranges from x=-10 to x=10
	int i;
	if(!initialized) {
		for(i=0;i<AXETS;++i) ExpTable[i]=exp(AXETR*(2*i-AXETS)/AXETS);
		initialized=true;
	}
	x *= 0.5*AXETS/AXETR;
	i = (int)x;
	x -= (double)i; //x = residual
	i += AXETS/2; //zero offset
	if(i >= AXETS-1) return ExpTable[AXETS-1];
	if(i <= 0) return 0.0; //ExpTable[0];
	return (1.0 - x)*ExpTable[i] + x*ExpTable[i+1];
}

// Expectation \f$\int_x^\infty {\cal N}(x) x dx\f$ when integrated from -infty to x
inline double gaussIntExpectation(double x) 
{
	double norm = gaussInt(x)*M_1_SQRT_2PI;
	return (-norm * approxExp(-0.5*x*x));
}

inline double normcdf(double x) {
	return 0.5*(1.0 + erf(x*M_SQRT1_2));
}

inline double maxDiff(const Matrix<SDIM>& v, const Matrix<SDIM>& w) 
{
	double d, t = 0.0;
	for(int i = 0; i < SDIM; ++i){ 
		d = fabs(v[i] - w[i]);
		if( d > t ) t = d; 
	}
	return t;
}

// compute a rotation matrix that rotates a onto v in arbitrary dimensions
inline void rotationFromAtoB(Matrix<SDIM, SDIM>& R, const Matrix<SDIM>& a, const Matrix<SDIM>& v)
{
	if( maxDiff(a,v) < 1e-10 ) { //nothing to rotate!
		R = identity<SDIM>();
		return; 
	} 
	
	// compute b orthogonal to a such that (a,b) span the rotation plane
	Matrix<SDIM> b;
	b = v - a*tr(~a*v);
	double norm_b = tr(~b*b);
	if (norm_b == 0) {
		R = -identity<SDIM>();
		return;
	} else {
		b /= sqrt(tr(~b*b));
	}

	// compute rotation coefficients within the (a,b) plane, namely, R_2D = (v_a  -v_b ;  v_b  v_a)
	double v_a = tr(~v*a);     //component along a
	double v_b = tr(~v*b);     //component along b
	
	// compute columns of R:
	Matrix<SDIM> x, x_res;
	double x_a, x_b;
	for(int i = 0; i < SDIM; ++i){
		x.reset(); 
		x[i] = 1.0;      //x=i-th unit vector
		x_a = tr(~x*a);  //component along a
		x_b = tr(~x*b);  //component along b
		x_res = x - x_a*a - x_b*b;  //residual (rest) of the vector
		// rotated vector = residual + rotated a-component + rotated b-component
		x = x_res + (v_a*x_a - v_b*x_b)*a + (v_b*x_a + v_a*x_b)*b;
		for(int j = 0; j < SDIM; ++j) {
			R(j,i) = x[j]; //store as column of the final rotation
		}
	}
}

inline double uniform() {
	unsigned int i;
	rand_s(&i);
	return ((double) i) / UINT_MAX;
}

inline double uniform(double low, double high) { return low + uniform()*(high-low); }

inline double normal() {
	double u_1 = 0;
	while (u_1 == 0) {
		u_1 = uniform();
	}
	double u_2 = 0;
	while (u_2 == 0) {
		u_2 = uniform();
	}
	return sqrt(-2*log(u_1)) * sin(2*M_PI*u_2);
}

inline double normal2() {
	double w,v,rsq,fac;
	do{
		v = 2.0 * uniform() - 1.0;
		w = 2.0 * uniform() - 1.0;
		rsq = v*v + w*w;
	} while( rsq >= 1.0 || rsq == 0 );
	fac  = sqrt( -2.0*log(rsq)/rsq);
	return v*fac;
}

template <size_t _numRows, size_t _numColumns>
inline void rndMatrix(Matrix<_numRows, _numColumns>& q, double low, double high){
	for(size_t i = 0; i < _numRows; ++i) {
		for(size_t j = 0; j < _numColumns; ++j) {
			q(i,j) = uniform(low,high);
		}
	}
}

template <size_t size>
inline Matrix<size> sampleGaussian(const Matrix<size>& mean, const Matrix<size, size>& var) {
	Matrix<size> sample;
	for (int j = 0; j < size; ++j) {
		sample[j] = normal();
	}
	Matrix<size, size> SVec, SVal;
	jacobi(var, SVec, SVal);
	for (int i = 0; i < size; ++i) {
		SVal(i,i) = sqrt(SVal(i,i));
	}
	return SVec * SVal * sample + mean;
}


inline double pdf(double x) {
	return exp(-0.5*x*x) * M_1_SQRT_2PI;
}

inline double cdf(double x) {
	if (x < 0) {
		return 0.5 - 0.5*incompletegamma(0.5, 0.5*x*x);
	} else {
		return 0.5 + 0.5*incompletegamma(0.5, 0.5*x*x);
	}
}

inline double truncate(double x, double mean, double var, double& newMean, double& newVar) {
	double stddev = sqrt(var);
	double y = (x - mean) / stddev;
	double z = pdf(y) / cdf(y);
	newMean = mean - z*stddev;
	newVar = var*(1.0 - y*z - z*z);
	return (1.0 - y*z - z*z);
}

// Global variables - Callisto
int cal_environment;
int cal_bbox;
int cal_samples;
int cal_constraint;
int cal_ellipse;

inline double mod2pi(double x) 
{
	// Returns a value 0 <= x < 2*PI
	double result = fmod(x, 2*M_PI);
	if (result < 0) { result += (2*M_PI); }
	return result;
}

void drawUnitCircle(void *userdef, float time)
{
	glLineWidth(3.0f);
	//glDisable(GL_LIGHTING);
	glColor3f(0.0f, 0.0f, 0.0f);
	
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 36; i++) {
		float theta = (float)(i*10.0*M_DEG2RAD);
		glVertex3f(cosf(theta), sinf(theta), 0.0f);
	}
	glEnd();
	//glEnable(GL_LIGHTING);
	glLineWidth(2.0f);
}

inline void drawCovariance(const Matrix<SDIM>& mean, const Matrix<SDIM, SDIM>& var)
{
	int id;
	CAL_CreateUserDrawn(cal_ellipse, drawUnitCircle, NULL, (float)mean[0], (float)mean[1], 0.001f, &id);
	CAL_SetObjectScaling(id, 0.025f, 0.025f, 1.0f);

	CAL_CreateUserDrawn(cal_ellipse, drawUnitCircle, NULL, 0.0f, 0.0f, 0.0f, &id);
	Matrix<SDIM,SDIM> V, E;
	jacobi(var, V, E);
	CAL_SetObjectPosition(id, (float)mean[0], (float)mean[1], 0.001f);
	CAL_SetObjectOrientation(id, 0.0f, 0.0f, (float)mod2pi(atan2(V(1,0), V(0,0))) );
	CAL_SetObjectScaling(id, (float)(sqrt(E(0,0))), (float)(sqrt(E(1,1))), 1.0f);
}

// Environment initialization
inline void initEnvironment()
{
	CAL_Initialisation(true, true, true);

	// set view params
	CAL_SetViewParams(0, 1.0f, 1.5f, 10, 1.0f, 1.5f, 0.0f, 0.0f, 1.0f, 0.0f);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");

	// create bounding box
	CAL_CreateGroup(&cal_bbox, 0, false, "Bbox");
	CAL_SetGroupColor(cal_bbox, 0.0f, 0.0f, 0.0f);

	/*
	int np[1] = {2};
	float side1[6] = {-2.0f, -2.0f, 0.0f, 4.0f, -2.0f, 0.0f};
	float side2[6] = {4.0f, -2.0f, 0.0f, 4.0f, 5.0f, 0.0f};
	float side3[6] = {4.0f, 5.0f, 0.0f, -2.0f, 5.0f, 0.0f};
	float side4[6] = {-2.0f, 5.0f, 0.0f, -2.0f, -2.0f, 0.0f};
	CAL_CreatePolyline(cal_bbox, 1, np, side1);
	CAL_CreatePolyline(cal_bbox, 1, np, side2);
	CAL_CreatePolyline(cal_bbox, 1, np, side3);
	CAL_CreatePolyline(cal_bbox, 1, np, side4);
	*/

	// constraint line
	CAL_CreateGroup(&cal_constraint, 0, false, "Constraint");
	CAL_SetGroupColor(cal_constraint, 1.0f, 0.0f, 0.0f);

	// ellipses
	CAL_CreateGroup(&cal_ellipse, 0, false, "Covariances");
	CAL_SetGroupColor(cal_ellipse, 0.0f, 0.0f, 0.0f);

	// samples
	CAL_CreateGroup(&cal_samples, 0, false, "Samples");
	CAL_SetGroupColor(cal_samples, 0.0f, 0.0f, 0.0f);
}

struct Constraint
{
	Matrix<SDIM> a;
	double b;

	Constraint(const Matrix<SDIM>& cdir, const Matrix<SDIM>& coff) {
		a = cdir;
		b = tr(~cdir*coff);
	}

	Constraint(const Matrix<SDIM>& ca, double cb) {
		a = ca; 
		b = cb;
	}

	bool satisfy(const Matrix<SDIM>& x) {
		return (tr(~a*x) - b >= 0);
	}

	void drawConstraint()
	{
		double p1x, p1y, p2x, p2y;
		//p1x = -2.0; p2x = 4.0;
		p1x = -5.0; p2x = 5.0;
		if (a[1] != 0) {
			p1y = (b - a[0]*p1x)/a[1];
			p2y = (b - a[0]*p2x)/a[1];
		} else {
			p1y = p2y = 0;
		}
		int np[1] = {2};
		float side[6] = {(float)p1x, (float)p1y, 0.001f, (float)p2x, (float)p2y, 0.001f};
		int id;
		CAL_CreatePolyline(cal_constraint, 1, np, side, &id);
	}

	double probFailure(const Matrix<SDIM>& m, const Matrix<SDIM, SDIM>& S)
	{
		std::cout << "b: " << b << std::endl;
		std::cout << "tr(~a*m): " << tr(~a*m) << std::endl;
		std::cout << (b - tr(~a*m)) << std::endl;
		double ps = normcdf((b - tr(~a*m))/(sqrt(tr(~a*S*a))));
		return ps;
	}
};

#endif