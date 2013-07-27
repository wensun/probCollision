#ifndef __COMMON_H__
#define __COMMON_H__

#define _CRT_RAND_S

// constant definitions
#define M_SQRT_PI_2 sqrt(M_PI/2.0)
#define M_1_SQRT2PI 1.0/sqrt(M_PI*2.0)
#define M_DEG2RAD M_PI/180.0

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <set>
#include <deque>

#include <windows.h>

#include <math.h>
#include <time.h>
#include <float.h>

//#define _SECURE_SCL 0

#include "callisto.h"
#include "glut.h"
#include "matrix.h"
#include "gammafunc.h"

//*****************************************************************************
// utility routines
inline double mod2pi(double x) {
	double result = fmod(x, (2*M_PI));
	if (result < 0) { result += (2*M_PI); }
	return result;
}

inline double sqr(double x) {
	return x*x;
}

inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) 
{
	double x = R(2,1) - R(1,2);
	double y = R(0,2) - R(2,0);
	double z = R(1,0) - R(0,1);
	double r = sqrt(x*x+y*y+z*z);
	double t = R(0,0) + R(1,1) + R(2,2);
	double angle = atan2(r,t-1);
	if (angle != 0) {
		x /= r;
		y /= r;
		z /= r;
	} else {
		x = 0;
		y = 0;
		z = 0;
	}
	Matrix<4,1> q;
	q(0,0) = sin(angle/2)*x;
	q(1,0) = sin(angle/2)*y;
	q(2,0) = sin(angle/2)*z;
	q(3,0) = cos(angle/2);

	return q;
}

// uniform random sample in [0,1]
inline double uniform() {
	return (double) rand() / RAND_MAX;
}


inline double uniform(double low, double high) { 
	return low + ((double)rand()/(double)RAND_MAX)*(high-low); 
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

template <size_t size>
inline Matrix<size> sampleUniform(const Matrix<size>& mean, const Matrix<size, size>& var) 
{
	Matrix<size> sample;
	for (int j = 0; j < size; ++j) {
		sample[j] = 2*(uniform()-0.5) * sqrt(3.0 * var(j,j));
	}
	return sample + mean;
}

template <size_t _numRows, size_t _numColumns>
inline void rndMatrix(Matrix<_numRows, _numColumns>& q, double low, double high){
	for(size_t i = 0; i < _numRows; ++i) {
		for(size_t j = 0; j < _numColumns; ++j) {
			q(i,j) = uniform(low,high);
		}
	}
}

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

inline double pdf(double x) {
	return exp(-0.5*x*x) * M_1_SQRT2PI;
}

inline double cdf(double x) {
	if (x < 0) {
		return 0.5 - 0.5*incompletegamma(0.5, 0.5*x*x);
	} else {
		return 0.5 + 0.5*incompletegamma(0.5, 0.5*x*x);
	}

	//return 0.5*(1.0 + erf(x*M_SQRT1_2));
}

inline void truncate(double x, double mean, double var, double& newMean, double& newVar) {
	double stddev = sqrt(var);
	double y = (x - mean) / stddev;
	double z = pdf(y) / cdf(y);
	//std::cout << "y: " << y << " z: " << z << " (y+z): " << (y+z) << " pdf(y): " << pdf(y) << " cdf(y): " << cdf(y) << std::endl;
	newMean = mean - z*stddev;
	//newVar = var*(1.0 - y*z - z*z);
	newVar = var*(1.0 - z*(y + z));
	//std::cout << "newVar: " << newVar << std::endl;
}

inline void drawUnitCircle(void *userdef, float time)
{
	glDisable(GL_LIGHTING);
	glLineWidth(4.0f);
	//glColor4f(0.0f, 1.0f, 0.0f, 0.5f);
	//glColor3f(0.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 0.0f);

	glBegin(GL_LINE_LOOP);
	//glBegin(GL_TRIANGLE_FAN);
	for (int i = 0; i < 72; i++) {
		float theta = (float)(i*5.0*M_DEG2RAD);
		glVertex3f(cosf(theta), sinf(theta), 0.0f);
	}
	glEnd();
	glLineWidth(1.0f);
	glEnable(GL_LIGHTING);
}

//*****************************************************************************

#include "bsp2D.h"

#include "car2D.h"
//#include "point2D.h"

std::vector<Matrix<XDIM, XDIM> > A; // process matrix
std::vector<Matrix<XDIM, UDIM> > B; // input matrix
std::vector<Matrix<ZDIM, XDIM> > H; // measurement matrix

std::vector<Matrix<XDIM, UDIM> > V; // process noise matrix
std::vector<Matrix<ZDIM, ZDIM> > W; // measurement noise matrix

Matrix<XDIM, XDIM> P0; // initial state covariance
Matrix<UDIM, UDIM> M; // process noise covariance
Matrix<ZDIM, ZDIM> N; // measurement noise covariance

Matrix<XDIM, XDIM> C; // LQR state deviation cost
Matrix<UDIM, UDIM> D; // LQR input deviation cost

std::vector<Matrix<UDIM, XDIM> > L; // feedback matrix
std::vector<Matrix<XDIM, ZDIM> > K; // Kalman-gain matrix

std::vector<Matrix<2*XDIM> > y; // a priori means
std::vector<Matrix<2*XDIM, 2*XDIM> > R; // a priori beliefs

std::vector<double> alpha;

struct RRTNode {
	Matrix<XDIM> x;
	Matrix<UDIM> u;
	int bp;
};

struct PathNode {
	Matrix<XDIM> x;
	Matrix<UDIM> u;

	void serialize(std::ofstream & fptr) {
		serializeMatrix(fptr, x); 
		serializeMatrix(fptr, u);
	}

	void deserialize(std::ifstream & fptr) {
		deserializeMatrix(fptr, x); 
		deserializeMatrix(fptr, u);
	}
};

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

void truncateGaussian(Matrix<SDIM>& m, Matrix<SDIM, SDIM>& S, const Matrix<SDIM>& a, double b)
{
	// linear transform to make (a,A) to a standard Gaussian
	Matrix<SDIM, SDIM> M;
	cholesky(S, M);

	// transform and rescale the constraint 
	double z = tr(~a*m) - b;

	Matrix<SDIM> v = M*a;
	double norm_v = sqrt(tr(~v*v));
	if (norm_v < 1e-10) {
		std::cerr << "No gradient for Gaussian trunctation!" << std::endl;
		std::exit(-1);
	}
	z /= norm_v;
	v /= norm_v;

	// build rotation matrix for constraint to be along the x-axis
	Matrix<SDIM> e_1;
	e_1.reset(); e_1[0] = 1.0;
	Matrix<SDIM, SDIM> R;
	rotationFromAtoB(R, e_1, v);

	// get mean and variance along x-axis
	double mean, var;

	double norm = M_SQRT_PI_2 * (1.0 - erf(-z*M_SQRT1_2));
	if(norm < 1e-2) { std::cout << "Likelihood of that truncation is very low: " << norm << std::endl; }
	mean = exp(z*-z*0.5)/norm;
	var  = 1.0 - z*mean - mean*mean;

	Matrix<SDIM, SDIM> T;
	Matrix<SDIM> t;
	t.reset(); t[0] = mean;
	T = identity<SDIM>(); T(0,0) = var;

	// transform back
	S = ~M*R*T*~R*M;
	m = ~M*R*t + m;

	//checkNaN(m);
	//checkNaN(S);
}

inline double computeConfidence(const Matrix<SDIM>& pos, const Matrix<SDIM,SDIM>& R) 
{
	Matrix<SDIM,SDIM> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~EVec);
	Matrix<4,1> q = quatFromRot(Temp);

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

	Matrix<SDIM,SDIM> invScale = zeros<2,2>();
	invScale(0,0) = 1/(float)sqrt(EVal(0,0));
	invScale(1,1) = 1/(float)sqrt(EVal(1,1));
	Matrix<SDIM> transPos =  invScale * ~EVec * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);

	int num_pairs;
	CAL_GetClosestPairs (cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete[] results;

	return distance;
}

inline void drawEllipse2d(const Matrix<SDIM>& x, const Matrix<SDIM, SDIM>& S, int groupid, bool contour = false)
{
	Matrix<SDIM, SDIM> V, E;
	jacobi(S, V, E);
	Matrix<3,3> R = identity<3>();
	R.insert(0,0, V);
	Matrix<3,3> T = identity<3>();
	if (!contour) {
		T(1,1) = T(2,2) = 0.0;
		T(1,2) = -1.0; T(2,1) = 1.0;
	}
	Matrix<4,1> q = quatFromRot(R*T);

	int obj;
	if (contour) {
		CAL_CreateUserDrawn(groupid, drawUnitCircle, NULL, 0.0f, 0.0f, 0.0f, &obj);
		CAL_SetObjectPosition(obj, (float) x[0], (float) x[1], 0.05f);
	} else {
		CAL_CreateCylinder(groupid, 1.0f, 0.02f, 0.0f, 0.0f, 0.0f, &obj);
		CAL_SetObjectPosition(obj, (float) x[0], (float) x[1], 0.01f);
	}
	CAL_SetObjectQuaternion(obj, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	if (contour) {
		CAL_SetObjectScaling(obj, (float) (3*sqrt(E(0,0))), (float) (3*sqrt(E(1,1))), 1.0f);
	} else {
		CAL_SetObjectScaling(obj, (float) (3*sqrt(E(0,0))), 1.0f, (float) (3*sqrt(E(1,1))));
	}
}

inline void drawPath2d(const std::vector<PathNode>& path, int groupid)
{
	//int l = (int)path.size();
	//int np[1] = {l};
	//float* p = new float[3*np[0]];
	//for (int i = 0; i < l; ++i) {
	//	p[3*i] = (float) path[i].x[0]; p[3*i+1] = (float) path[i].x[1]; p[3*i+2] = 0;
	//}
	//CAL_CreatePolyline(groupid, 1, np, p);
	//delete[] p;

	
	int l = (int)path.size();
	Vector2d p1, p2, mid, dv;
	double len;
	p1.x = path[0].x[0]; p1.y = path[0].x[1];

	for(int i = 1; i < l; ++i) 
	{
		p2.x = path[i].x[0]; p2.y = path[i].x[1];
		dv = Vector2d(p2.x - p1.x, p2.y - p1.y);
		len = norm(dv);
		dv = dv/len;
		mid = Vector2d(0.5*(p1.x + p2.x), 0.5*(p1.y + p2.y));

		int obj;
		CAL_CreateBox(groupid, (float)(len + 0.1), 0.1f, 0.01f, (float)mid.x, (float)mid.y, 0.02f, &obj);
		CAL_SetObjectOrientation(obj, 0.0f, 0.0f, (float)(mod2pi(atan2(dv.y, dv.x))));

		p1 = p2;
	}


}

class Timer
{
public:
	Timer() {
		this->init();
	}

	~Timer(){};

	bool start() {

		bool bSuccess = false;

		if(!m_bTimerRunning && m_bTimerSupported)
		{
			m_startCount = 0;
			m_stopCount = 0;
			m_interval = 0;

			if(QueryPerformanceCounter((LARGE_INTEGER*)&m_startCount))
			{
				m_bTimerRunning = true;
				bSuccess = true;
			}
		}

		return bSuccess;
	}

	
	bool stop() {

		bool bSuccess = false;

		if(m_bTimerRunning && m_bTimerSupported)
		{
			if(QueryPerformanceCounter((LARGE_INTEGER*)&m_stopCount))
			{
				m_bTimerRunning = false;
				bSuccess = true;
			}
		}

		return bSuccess;
	}

	void reset() {
		this->init();
	}

	double interval_S() {
		return ((double)(m_stopCount - m_startCount) - m_adjustCount) / (double)m_frequency;
	}

	double interval_mS() {
		return (((m_stopCount - m_startCount) - m_adjustCount) * 1000.0) / (double)m_frequency;
	}

	double interval_uS() {
		return (((m_stopCount - m_startCount) - m_adjustCount) * 1000000.0) / (double)m_frequency;
	}

	double resolution_S() {
		return 1.0 / (double)m_frequency;
	}

	double resolution_mS() {
		return 1000.0 / (double)m_frequency;
	}

	double resolution_uS() {
		return 1000000.0 / (double)m_frequency;
	}

	double correction_uS() {
		return (m_adjustCount * 1000000.0) / (double)m_frequency;
	}

	void init() {

		m_frequency = 0;
		m_adjustCount = 0;
		m_bTimerSupported = false;
		m_bTimerRunning = false;

		if(QueryPerformanceFrequency((LARGE_INTEGER*)&m_frequency))
		{
			m_bTimerSupported = true;

			// Measure the 'Stop' function call overhead
			const int iNumSamples = 10;
			__int64 samples[iNumSamples];
			__int64 countTot = 0;
			double dAvCount = 0.0;
			double dAvDeviance = 0.0;

			for(int i = 0; i < iNumSamples; i++)
			{
				this->start();
				this->stop();

				samples[i] = m_stopCount - m_startCount;
				countTot += samples[i];
			}

			dAvCount = (double)countTot / (double)iNumSamples;

			// Get the average deviance
			for(int i = 0; i < iNumSamples; i++)
			{
				dAvDeviance += fabs(((double)samples[i]) - dAvCount);
			}

			// Average deviance only required for debug
			dAvDeviance /= iNumSamples;
			m_adjustCount = (__int64)dAvCount;
		}
	}

	void set_start_time(__int64 & c) {
		c = 0;
		QueryPerformanceCounter((LARGE_INTEGER*) &c);
	}

	double elapsed_time(__int64 & c) {
		__int64 end = 0;
		QueryPerformanceCounter((LARGE_INTEGER*)&end);
		return (double)((end - c) - m_adjustCount) / (double)m_frequency;
	}

private:
	__int64 m_interval;
	__int64 m_frequency;
	__int64 m_startCount;
	__int64 m_stopCount;
	__int64 m_adjustCount;
	bool m_bTimerSupported;
	bool m_bTimerRunning;
};

#endif