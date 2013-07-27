#ifndef __NEEDLE3D_H__
#define __NEEDLE3D_H__

#define _CRT_RAND_S

// constant definitions
#define M_SQRT_PI_2 sqrt(M_PI/2.0)
#define M_1_SQRT2PI 1.0/sqrt(M_PI*2.0)
#define M_DEG2RAD M_PI/180.0
#define INFTY 9e9

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
#include <stack>

#include <math.h>
#include <time.h>
#include <float.h>

//#define _SECURE_SCL 0

#include "callisto.h"
#include "glut.h"
#include "matrix.h"
#include "gammafunc.h"

#define ZDIM 3

#define SAVEFILE "needlePathX.txt"

int cal_environment;
int cal_obstacles;
int cal_box;
int cal_goal;
int cal_rrt;
int cal_paths, cal_paths_trunc;
int cal_ellipse, cal_ellipse_trunc;
int cal_point, cal_plane;

int skip = 4;

double Xrange, Yrange, Zrange;
double speed;
double needleSteeringRadius;
double maxCurvature;
double maxTwist;

double dt; 
double factor;

Matrix<3> goal;

std::vector<Matrix<6,6> > A; // process matrix
std::vector<Matrix<6,3> > B; // input matrix
std::vector<Matrix<6,6> > N; // state noise covariance

std::vector<Matrix<ZDIM,6> > H; // measurement matrix
std::vector<Matrix<3,6> > L; // feedback matrix
std::vector<Matrix<6,ZDIM> > K; // Kalman-gain matrix

Matrix<6,6> P0; // initial state covariance
Matrix<ZDIM, ZDIM> Q; // measurement noise covariance
Matrix<6,6> M; // process noise covariance
Matrix<6,6> C; // LQR state deviation cost
Matrix<3,3> D; // LQR input deviation cost

std::vector<Matrix<12> > y; // a priori means
std::vector<Matrix<12,12> > R; // a priori beliefs

struct TreeNode {
	Matrix<4,4> T;
	double v, w, k;

	TreeNode() {
		T = identity<4>();
		v = 0; w = 0; k = maxCurvature;
	}
	
	int bp;
	double depth;
	
	std::vector<int> children;
};

struct PathNode {
	Matrix<4,4> T;
	double v, w, k;

	PathNode() {
		T = identity<4>();
		v = 0; w = 0; k = maxCurvature;
	}

	void serialize(std::ofstream & fptr) {
		serializeMatrix(fptr, T);
		fptr << " " << v << " " << w << " " << k;
	}

	void deserialize(std::ifstream & fptr) {
		deserializeMatrix(fptr, T);
		fptr >> v >> w >> k;
	}
};

struct Obstacle {
	Matrix<3> a;
	double b;
	Matrix<3> intpt;

	Obstacle(Matrix<3>& avec, double bval, Matrix<3>& pt): a(avec), b(bval), intpt(pt) {};
};

std::vector<Obstacle> obs;

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

inline Matrix<4,1> quatFromAA(const Matrix<3>& axis, double angle)
{
	Matrix<4,1> q;
	double sa = sin(-angle*0.5);
	double ca = cos(angle*0.5);
	q(0,0) = sa*axis[0];
	q(1,0) = sa*axis[1];
	q(2,0) = sa*axis[2];
	q(3,0) = ca;

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

inline Matrix<3,3> cpMatrix(const Matrix<3>& a) {
	Matrix<3,3> A;
	A(0,0) = 0;       A(0,1) = -a(2,0); A(0,2) = a(1,0);
	A(1,0) = a(2,0);  A(1,1) = 0;       A(1,2) = -a(0,0);
	A(2,0) = -a(1,0); A(2,1) = a(0,0);  A(2,2) = 0;
	return A;
}

inline Matrix<3,3> rotFromErr(const Matrix<3,1>& q) {
	double rr = q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0);
	if (rr == 0) {
		return identity<3>();
	} else {
		double r = sqrt(rr);
		//printf_s("%24.24g ", cos(r) - 1);
		return cpMatrix(q * (sin(r) / r)) + identity<3>() * cos(r) + (q*~q) * ((1 - cos(r)) / rr);
		//return cpMatrix(q * (sin(r) / r)) + identity<3>() * cos(r) + (q*~q)/rr - (q*~q)*(cos(r)/rr);
	}
}

inline Matrix<3,1> errFromRot(const Matrix<3,3>& R) {
	Matrix<3,1> q;
	q(0,0) = R(2,1) - R(1,2);
	q(1,0) = R(0,2) - R(2,0);
	q(2,0) = R(1,0) - R(0,1);

	double r = sqrt(q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0));
	double t = R(0,0) + R(1,1) + R(2,2) - 1;

	if (r == 0) {
		return zeros<3,1>();
	} else {
		return q * (atan2(r, t) / r);
	}
}

inline Matrix<4,4> transFromErr(const Matrix<6,1>& x) {
	Matrix<4,4> X;
	X = identity<4>();
	X.insert(0,0, rotFromErr(x.subMatrix<3,1>(3,0)));
	X.insert(0,3, x.subMatrix<3,1>(0,0));
	return X;
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

void drawUnitSphere(void *userdef, float time)
{
	glDisable(GL_LIGHTING);
	glLineWidth(2.0f);
	glColor3f(0.0f, 0.0f, 0.0f);
	
	glutWireSphere(1.0, 15, 15);

	glLineWidth(1.0f);
	glEnable(GL_LIGHTING);
}

inline void drawEllipse3d(const Matrix<3>& x, const Matrix<3,3>& S, int groupid, bool contour = false)
{
	Matrix<3,3> V, E;
	jacobi(S, V, E);
	Matrix<3,3> R = identity<3>();
	R.insert(0,0, V);
	Matrix<4,1> q = quatFromRot(R);

	int obj;
	if (contour) {
		CAL_CreateUserDrawn(groupid, drawUnitSphere, NULL, 0.0f, 0.0f, 0.0f, &obj);
	} else {
		CAL_CreateSphere(groupid, 1.0f, 0.0f, 0.0f, 0.0f, &obj);
	}
	CAL_SetObjectPosition(obj, (float) x[0], (float) x[1], (float) x[2]);
	CAL_SetObjectQuaternion(obj, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	
	double scale = 3.0;
	CAL_SetObjectScaling(obj, (float) (scale*sqrt(E(0,0))), (float) (scale*sqrt(E(1,1))), (float) (scale*sqrt(E(2,2))));
}

inline void drawPath3d(const std::vector<PathNode>& path, int groupid)
{
	int l = (int)path.size();
	
	//int np[1] = {l};
	//float* p = new float[3*np[0]];
	//for (int i = 0; i < l; ++i) {
	//	p[3*i] = (float)path[i].T(0,3); p[3*i+1] = (float)path[i].T(1,3); p[3*i+2] = (float)path[i].T(2,3);
	//}
	//CAL_CreatePolyline(groupid, 1, np, p);
	//delete[] p;

	Matrix<3> p1, p2, m, t, c;
	Matrix<4,1> q;
	float len1, len2;
	int obj;

	p1[0] = (float)path[0].T(0,3); p1[1] = (float)path[0].T(1,3); p1[2] = (float)path[0].T(2,3);
	for (int i = 1; i < l; ++i) {
		p2[0] = (float)path[i].T(0,3); p2[1] = (float)path[i].T(1,3); p2[2] = (float)path[i].T(2,3);
		m = (p1 + p2)*0.5;
		len1 = sqrt(tr((p2 - p1)*~(p2 - p1)));
		t = (p2 - p1)/len1;
		c[0] = -t[2]; c[1] = 0.0; c[2] = t[0];
		len2 = sqrt(t[0]*t[0] + t[2]*t[2]);
		c = c/len2;
		CAL_CreateCylinder(groupid, 0.04f, len1, 0.0f, 0.0f, 0.0f, &obj);
		q = quatFromAA(c, acos(t[1]));
		CAL_SetObjectQuaternion(obj, (float) q[0], (float) q[1], (float) q[2], (float) q[3]);
		CAL_SetObjectPosition(obj, (float)m[0], (float)m[1], (float)m[2]);
		p1 = p2;
	}
}

inline void drawPolyline3d(int pathlen, float* pts, int groupid)
{
	Matrix<3> p1, p2, m, t, c;
	Matrix<4,1> q;
	float len1, len2;
	int obj;

	p1[0] = (float)pts[0]; p1[1] = (float)pts[1]; p1[2] = (float)pts[2];
	for (int i = 1; i < pathlen; ++i) {
		p2[0] = (float)pts[3*i]; p2[1] = (float)pts[3*i+1]; p2[2] = (float)pts[3*i+2];
		m = (p1 + p2)*0.5;
		len1 = sqrt(tr((p2 - p1)*~(p2 - p1)));
		t = (p2 - p1)/len1;
		c[0] = -t[2]; c[1] = 0.0; c[2] = t[0];
		len2 = sqrt(t[0]*t[0] + t[2]*t[2]);
		c = c/len2;
		CAL_CreateCylinder(groupid, 0.04f, len1, 0.0f, 0.0f, 0.0f, &obj);
		q = quatFromAA(c, acos(t[1]));
		CAL_SetObjectQuaternion(obj, (float) q[0], (float) q[1], (float) q[2], (float) q[3]);
		CAL_SetObjectPosition(obj, (float)m[0], (float)m[1], (float)m[2]);
		p1 = p2;
	}
}

#endif