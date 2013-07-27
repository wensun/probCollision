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

#undef abs
#include <algorithm>
#include "ors.h"

#define SL_DEBUG_LEVEL 1
#define SL_DEBUG(l,x) if(l<=SL_DEBUG_LEVEL) x;


#define Qstate

const ors::Vector VEC_x(1,0,0);
const ors::Vector VEC_y(0,1,0);
const ors::Vector VEC_z(0,0,1);

//===========================================================================
//
// general documentation
//

/*! \brief Open Robot Simulation Toolkit

    This namespace defines some core data structures for robot
    simulation and linking to external simulation engines. In
    particular, using ors we can implement a soc::SocSystemAbstraction. See
    the <a href="../guide.pdf">guide</a> for an introduction.

    Please see also the header <a
    href="ors_8h-source.html">MT/ors.h</a> -- it gives a good
    overview. */
namespace ors{};


namespace ors{

  //! copy operator
  /*
  Vector& Vector::operator=(const Vector& b){
    v[0]=b.v[0]; v[1]=b.v[1]; v[2]=b.v[2]; return *this;
  }*/

  //! copy operator
  /*Vector& Vector::operator=(const real* b){
    v[0]=b[0]; v[1]=b[1]; v[2]=b[2]; return *this;
  }*/

  //! set all entries to same value
  /*Vector& Vector::operator=(real b){
    v[0]=v[1]=v[2]=b; return *this;
  }*/


//{ access
  //! lhs reference
  real& Vector::operator()(int i){ CHECK(i>=0 && i<3,"ors::Vector access - out of range"); return v[i]; }
  const real& Vector::operator()(int i) const{ CHECK(i>=0 && i<3,"ors::Vector access - out of range"); return v[i]; }
  //real& Vector::operator[](int i){ CHECK(i>=0 && i<3,"ors::Vector access - out of range"); return v[i]; }

#ifdef MT_MSVC
  //! real-pointer access
  //Vector::operator real*(){ return v; }
#endif

  //! real-pointer access
  //Vector::operator const real*() const{ return v; }

  //! set the vector
  void Vector::set(real x,real y,real z){ v[0]=x; v[1]=y; v[2]=z; }

  //! set the vector
  void Vector::set(real* x){ v[0]=x[0]; v[1]=x[1]; v[2]=x[2]; }

  //! set the vector
  void Vector::setZero(){ v[0]=v[1]=v[2]=0.; }

  //! a random vector in [-1,1]^3
  void Vector::setRandom(){ v[0]=rnd.uni(-1,1); v[1]=rnd.uni(-1,1); v[2]=rnd.uni(-1,1); }

//{ vector operations

  /*
  void Vector::multiply(real c){
    v[0]*=c; v[1]*=c; v[2]*=c;
  }

  void Vector::divide(real c){
    v[0]/=c; v[1]/=c; v[2]/=c;
  }*/

  //! this=this+b
  void Vector::add(real x,real y,real z){ v[0]+=x; v[1]+=y; v[2]+=z; }

  //! this=this-b
  void Vector::subtract(real x,real y,real z){ v[0]-=x; v[1]-=y; v[2]-=z; }

  //! this=this/length(this)
  void Vector::normalize(){ (*this)/=length(); }

  //! this=this*l/length(this)
  void Vector::setLength(real l){
    if(isZero()) MT_MSG("can't change length of null vector");
    (*this)*=l/length();
  }

  //! this=component of this normal to \c b, (unnormalized!)
  void Vector::makeNormal(const Vector& b){
    real l=b.length(),s=v[0]*b.v[0]+v[1]*b.v[1]+v[2]*b.v[2];
    s/=l*l;
    v[0]-=s*b.v[0]; v[1]-=s*b.v[1]; v[2]-=s*b.v[2];
  }

  //! this=component of this colinear to \c b, (unnormalized!)
  void Vector::makeColinear(const Vector& b){
    // *this = ((*this)*b)/b.length()) * (*this);
    real l=b.length(),s=v[0]*b.v[0]+v[1]*b.v[1]+v[2]*b.v[2];
    s/=l*l;
    v[0]=s*b.v[0]; v[1]=s*b.v[1]; v[2]=s*b.v[2];
  }

//{ measuring the vector

  //! is zero?
  bool Vector::isZero() const{ return (v[0]==0. && v[1]==0. && v[2]==0.); }

  //! is it normalized?
  bool Vector::isNormalized() const{ return fabs(lengthSqr()-1.)<1e-6; }

  //! returns the length of this
  real Vector::length() const{ return ::sqrt(lengthSqr()); }

  //! returns the square of length |a|^2
  real Vector::lengthSqr() const{ return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }

  //! angle in [0..pi] between this and b
  real Vector::angle(const Vector& b) const{
    real x=((*this)*b)/(length()*b.length());
    if(x<-1.) x=-1.;
    if(x>1.) x=1.;
    return ::acos(x);
  }

  /*!\brief if \c this and \c b are colinear, it returns the factor c such
      that this=c*b; otherwise it returns zero */
  real Vector::isColinear(const Vector& b) const{
    real c=v[0]/b.v[0];
    if(v[1]==c*b.v[1] && v[2]==c*b.v[2]) return c;
    return 0.;
  }


//{ sphere coordinates

  //! the radius in the x/y-plane
  real Vector::radius() const{ return ::sqrt(v[0]*v[0]+v[1]*v[1]); }
  //! the angle in the x/y-plane in [-pi, pi]
  real Vector::phi() const{
    real p;
    if(v[0]==0. || ::fabs(v[0])<1e-10) p=MT_PI/2.; else p=::atan(v[1]/v[0]);
    if(v[0]<0.){ if(v[1]<0.) p-=MT_PI; else p+=MT_PI; }
    return p;
  }
  //! the angle from the x/y-plane
  real Vector::theta() const{ return ::atan(v[2]/radius())+MT_PI/2.; }


//{ I/O
  void Vector::write(std::ostream& os) const{
    if(!MT::IOraw) os <<'(' <<v[0] <<' ' <<v[1] <<' ' <<v[2] <<')';
    else os <<' ' <<v[0] <<' ' <<v[1] <<' ' <<v[2];
  }
  void Vector::read(std::istream& is){
    if(!MT::IOraw) is >>"(" >>v[0] >>v[1] >>v[2] >>")";
    else is >>v[0] >>v[1] >>v[2];
  }
  //}

  //! scalar product (inner product)
  real operator*(const Vector& a,const Vector& b){
    return a.v[0]*b.v[0]+a.v[1]*b.v[1]+a.v[2]*b.v[2];
  }

  //! cross product (corresponds to antisymmetric exterior product)
  Vector operator^(const Vector& b,const Vector& c){
    Vector a;
    a.v[0]=b.v[1]*c.v[2]-b.v[2]*c.v[1];
    a.v[1]=b.v[2]*c.v[0]-b.v[0]*c.v[2];
    a.v[2]=b.v[0]*c.v[1]-b.v[1]*c.v[0];
    return a;
  }

  //! sum of two vectors 
  Vector operator+(const Vector& b,const Vector& c){
    Vector a;
    a.v[0]=b.v[0]+c.v[0];
    a.v[1]=b.v[1]+c.v[1];
    a.v[2]=b.v[2]+c.v[2];
    return a;
  }

  //! difference between two vectors
  Vector operator-(const Vector& b,const Vector& c){
    Vector a;
    a.v[0]=b.v[0]-c.v[0];
    a.v[1]=b.v[1]-c.v[1];
    a.v[2]=b.v[2]-c.v[2];
    return a;
  }

  //! multiplication with a scalar
  Vector operator*(real b,const Vector& c){
    Vector a;
    a.v[0]=b*c.v[0];
    a.v[1]=b*c.v[1];
    a.v[2]=b*c.v[2];
    return a;
  }

  //! multiplication with a scalar
  Vector operator*(const Vector& b,real c){ return c*b; }

  //! division by a scalar
  Vector operator/(const Vector& b,real c){ return (1./c)*b; }

  //! multiplication with a scalar
  Vector& operator*=(Vector& a,real c){
    a.v[0]*=c; a.v[1]*=c; a.v[2]*=c;
    return a;
  }
  //! divide by a scalar
  Vector& operator/=(Vector& a,real c){
    a.v[0]/=c; a.v[1]/=c; a.v[2]/=c;
    return a;
  }
  //! add a vector
  Vector& operator+=(Vector& a,const Vector& b){
    a.v[0]+=b.v[0]; a.v[1]+=b.v[1]; a.v[2]+=b.v[2];
    return a;
  }
  //! subtract a vector
  Vector& operator-=(Vector& a,const Vector& b){
    a.v[0]-=b.v[0]; a.v[1]-=b.v[1]; a.v[2]-=b.v[2];
    return a;
  }
  //! return the negative of a vector
  Vector operator-(const Vector& b){
    Vector a;
    a.v[0]=-b.v[0]; a.v[1]=-b.v[1]; a.v[2]=-b.v[2];
    return a;
  }

  //! const access via two row and column indices
  const real& Matrix::operator()(int i,int j) const{ return m[i*3+j]; }
  //! LHS access via two row and column indices
  real& Matrix::operator()(int i,int j){ return m[i*3+j]; }

  //! copy operator
  /*Matrix& Matrix::operator=(const Matrix& b){
    memmove(m,b.m,9*sizeof(real)); return *this;
  }*/

  //! reset to zero
  void Matrix::setZero(){
    m[0]=m[4]=m[8]=0.;
    m[1]=m[2]=m[3]=m[5]=m[6]=m[7]=0.;
  }

  //! reset to identity
  void Matrix::setId(){
    m[0]=m[4]=m[8]=1.;
    m[1]=m[2]=m[3]=m[5]=m[6]=m[7]=0.;
  }

  //! assign the matrix to the transformation from unit frame to given XYZ frame
  void Matrix::setFrame(Vector& X,Vector& Y,Vector& Z){
    m[0]=X(0); m[1]=Y(0); m[2]=Z(0);
    m[3]=X(1); m[4]=Y(1); m[5]=Z(1);
    m[6]=X(2); m[7]=Y(2); m[8]=Z(2);
  }

  //! assign the matrix to the transformation from the ORTHOGONAL XYZ frame to the unit frame
  void Matrix::setInvFrame(Vector& X,Vector& Y,Vector& Z){
    m[0]=X(0); m[1]=X(1); m[2]=X(2);
    m[3]=Y(0); m[4]=Y(1); m[5]=Y(2);
    m[6]=Z(0); m[7]=Z(1); m[8]=Z(2);
  }

  //! assign the matrix to a rotation around the X-axis with angle a (in rad units)
  void Matrix::setXrot(real a){
    m[0]=1.; m[1]=0.;     m[2]=0.;
    m[3]=0.; m[4]=cos(a); m[5]=-sin(a);
    m[6]=0.; m[7]=sin(a); m[8]= cos(a);
  }

  void Matrix::setSkew(const Vector& a){
    m[0]=   0.; m[1]=-a(2); m[2]= a(1);
    m[3]= a(2); m[4]=   0.; m[5]=-a(0);
    m[6]=-a(1); m[7]= a(0); m[8]=   0.;
  }

  void Matrix::setExponential(const Vector& a){
    Matrix S;
    real phi=a.length();
    if(phi<1e-10){ setId(); return; }
    S.setSkew(a/phi);
    *this = sin(phi)*S + (1.-cos(phi))*S*S;
    m[0]+=1.; m[4]+=1.; m[8]+=1.;
  }

  void Matrix::setOdeMatrix(real* o){
    m[0]=o[0]; m[1]=o[1]; m[2]=o[2];
    m[3]=o[4]; m[4]=o[5]; m[5]=o[6];
    m[6]=o[8]; m[7]=o[9]; m[8]=o[10];
  }

  void Matrix::setTensorProduct(const Vector& b,const Vector& c){
    m[0]=b.v[0]*c.v[0]; m[1]=b.v[0]*c.v[1]; m[2]=b.v[0]*c.v[2];
    m[3]=b.v[1]*c.v[0]; m[4]=b.v[1]*c.v[1]; m[5]=b.v[1]*c.v[2];
    m[6]=b.v[2]*c.v[0]; m[7]=b.v[2]*c.v[1]; m[8]=b.v[2]*c.v[2];
  }

  void Matrix::write(std::ostream& os) const{
    os <<"\n " <<m[0] <<' ' <<m[1] <<' ' <<m[2];
    os <<"\n " <<m[3] <<' ' <<m[4] <<' ' <<m[5];
    os <<"\n " <<m[6] <<' ' <<m[7] <<' ' <<m[8];
    os <<endl;
  }
  void Matrix::read(std::istream& is){
    NIY;
  }
  //}

    //! multiplication of two matrices
  Matrix operator*(const Matrix& b,const Matrix& c){
    Matrix a;
    a(0,0)=b(0,0)*c(0,0)+b(0,1)*c(1,0)+b(0,2)*c(2,0);
    a(0,1)=b(0,0)*c(0,1)+b(0,1)*c(1,1)+b(0,2)*c(2,1);
    a(0,2)=b(0,0)*c(0,2)+b(0,1)*c(1,2)+b(0,2)*c(2,2);

    a(1,0)=b(1,0)*c(0,0)+b(1,1)*c(1,0)+b(1,2)*c(2,0);
    a(1,1)=b(1,0)*c(0,1)+b(1,1)*c(1,1)+b(1,2)*c(2,1);
    a(1,2)=b(1,0)*c(0,2)+b(1,1)*c(1,2)+b(1,2)*c(2,2);

    a(2,0)=b(2,0)*c(0,0)+b(2,1)*c(1,0)+b(2,2)*c(2,0);
    a(2,1)=b(2,0)*c(0,1)+b(2,1)*c(1,1)+b(2,2)*c(2,1);
    a(2,2)=b(2,0)*c(0,2)+b(2,1)*c(1,2)+b(2,2)*c(2,2);
    return a;
  }
  //! sum of two matrices
  Matrix operator+(const Matrix& b,const Matrix& c){
    Matrix a;
    a.m[0]=b.m[0]+c.m[0]; a.m[1]=b.m[1]+c.m[1]; a.m[2]=b.m[2]+c.m[2];
    a.m[3]=b.m[3]+c.m[3]; a.m[4]=b.m[4]+c.m[4]; a.m[5]=b.m[5]+c.m[5];
    a.m[6]=b.m[6]+c.m[6]; a.m[7]=b.m[7]+c.m[7]; a.m[8]=b.m[8]+c.m[8];
    return a;
  }
  //! transformation of a vector
  Vector operator*(const Matrix& b,const Vector& c){
    Vector a;
    a.v[0]=b.m[0]*c.v[0]+b.m[1]*c.v[1]+b.m[2]*c.v[2];
    a.v[1]=b.m[3]*c.v[0]+b.m[4]*c.v[1]+b.m[5]*c.v[2];
    a.v[2]=b.m[6]*c.v[0]+b.m[7]*c.v[1]+b.m[8]*c.v[2];
    return a;
  }
  //! multiplication with a scalar
  Matrix& operator*=(Matrix& a,real c){
    a.m[0]*=c; a.m[1]*=c; a.m[2]*=c;
    a.m[3]*=c; a.m[4]*=c; a.m[5]*=c;
    a.m[6]*=c; a.m[7]*=c; a.m[8]*=c;
    return a;
  }
  //! multiplication with scalar
  Matrix operator*(real b,const Matrix& c){
    Matrix a;
    a=c;
    a*=b;
    return a;
  }
  //! sum of two matrices
  Matrix& operator+=(Matrix& a,const Matrix& b){
    a.m[0]+=b.m[0]; a.m[1]+=b.m[1]; a.m[2]+=b.m[2];
    a.m[3]+=b.m[3]; a.m[4]+=b.m[4]; a.m[5]+=b.m[5];
    a.m[6]+=b.m[6]; a.m[7]+=b.m[7]; a.m[8]+=b.m[8];
    return a;
  }

    //! initializes to identity
  Quaternion::Quaternion(){
    setZero();
  }

  //! inverts the current rotation
  void Quaternion::invert(){ q[0]=-q[0]; }

  //! multiplies the rotation by a factor f (i.e., makes f-times the rotation)
  void Quaternion::multiply(real f){
    if(q[0]==1. || f==1.) return;
    real phi=acos(q[0]);
    phi*=f;
    q[0]=cos(phi);
    f=sin(phi)/sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    q[1]*=f; q[2]*=f; q[3]*=f;
  }

  bool Quaternion::isNormalized() const{
    real n=q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    return fabs(n-1.)<1e-6;
  }

  void Quaternion::normalize(){
    real n=q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    n=sqrt(n);
    q[0]/=n; q[1]/=n; q[2]/=n; q[3]/=n;
  }

  /*!\brief roughly, removes all ``components'' of the rotation that are not
      around the given vector v. More precisely, aligns/projects
      the rotation axis (given by q[1],q[2],q[3] of the quaternion)
      with v and re-normalizes afterwards. */
  void Quaternion::alignWith(const Vector& v){
    real s=q[1]*v(0) + q[2]*v(1) + q[3]*v(2);
    if(!s){ setZero(); return; }// are orthogonal
    s/=v*v;
    q[1]=s*v(0); q[2]=s*v(1); q[3]=s*v(2);
    normalize();
  }


  //! set the quad
  void Quaternion::set(real* x){ q[0]=x[0]; q[1]=x[1]; q[2]=x[2]; q[3]=x[3]; }
  //! set the quad
  void Quaternion::set(real q0,real x,real y,real z){ q[0]=q0; q[1]=x; q[2]=y; q[3]=z; }
  //! reset the rotation to identity
  void Quaternion::setZero(){ q[0]=1; q[1]=q[2]=q[3]=0; }
  //! samples the rotation uniformly from the whole SO(3)
  void Quaternion::setRandom(){
    real s,s1,s2,t1,t2;
    s=rnd.uni();
    s1=sqrt(1-s);
    s2=sqrt(s);
    t1=MT_2PI*rnd.uni();
    t2=MT_2PI*rnd.uni();
    q[0]=cos(t2)*s2;
    q[1]=sin(t1)*s1;
    q[2]=cos(t1)*s1;
    q[3]=sin(t2)*s2;
  }

  //! sets this to a smooth interpolation between two rotations
  void Quaternion::setInterpolate(real t,const Quaternion& a,const Quaternion b){
    real sign=1.;
    if(scalarProduct(a,b)<0) sign=-1.;
    q[0]=a.q[0]+t*(sign*b.q[0]-a.q[0]);
    q[1]=a.q[1]+t*(sign*b.q[1]-a.q[1]);
    q[2]=a.q[2]+t*(sign*b.q[2]-a.q[2]);
    q[3]=a.q[3]+t*(sign*b.q[3]-a.q[3]);
    normalize();
  }

  //! assigns the rotation to \c a DEGREES around the vector (x,y,z)
  void Quaternion::setDeg(real degree,real x,real y,real z){ setRad(degree*MT_PI/180.,x,y,z); }
  void Quaternion::setDeg(real degree,const Vector& vec){ setRad(degree*MT_PI/180.,vec(0),vec(1),vec(2)); }
  //! assigns the rotation to \c a RADIANTS (2*PI-units) around the vector (x,y,z)
  void Quaternion::setRad(real angle,real x,real y,real z){
    real l = x*x + y*y + z*z;
    if(l<1e-15){ setZero(); return; }
    angle/=2.;
    l=sin(angle)/sqrt(l);
    q[0]=cos(angle);
    q[1]=x*l;
    q[2]=y*l;
    q[3]=z*l;
  }
  //! ..
  void Quaternion::setRad(real angle,const Vector &axis){
    setRad(angle,axis(0),axis(1),axis(2));
  }
  //! assigns the rotation to \c a RADIANTS (2*PI-units) around the current axis
  void Quaternion::setRad(real angle){
    real l = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    if(l<1e-15){ setZero(); return; }
    angle/=2.;
    l=sin(angle)/sqrt(l);
    q[0]=cos(angle);
    q[1]*=l;
    q[2]*=l;
    q[3]*=l;
  }
  //! rotation around X-axis by given radiants
  void Quaternion::setRadX(real angle){
    angle/=2.;
    q[0]=cos(angle);
    q[1]=sin(angle);
    q[2]=q[3]=0.;
  }
  //! rotation around X-axis by given radiants
  void Quaternion::setRadY(real angle){
    angle/=2.;
    q[0]=cos(angle);
    q[2]=sin(angle);
    q[1]=q[3]=0.;
  }
  //! rotation around the given vector with angle (in rad) equal to norm of the vector
  void Quaternion::setVec(Vector w){
    real phi=w.length();
    setRad(phi,w(0),w(1),w(2));
  }
  //! rotation that will rotate 'from' to 'to' on direct path
  void Quaternion::setDiff(const Vector& from,const Vector& to){
    real phi=acos(from*to/(from.length()*to.length()));
    if(!phi) return;
    Vector axis(from^to);
    if(axis.isZero()) axis=Vector(0,0,1)^to;
    setRad(phi,axis);
  }

  //! is identical
  bool Quaternion::isZero() const{ return q[0]==1.; }

#ifdef MT_MSVC
  //! real-pointer access
  //Quaternion::operator real*(){ return q; }
#endif

  //! real-pointer access
  //Quaternion::operator const real*() const{ return q; }

  //! gets rotation angle (in rad [0,2pi])
  real Quaternion::getRad() const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)) return 0;
    return 2.*acos(q[0]);
  }

  //! gets rotation angle (in degree [0,360])
  real Quaternion::getDeg() const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)) return 0;
    return 360./MT_PI*acos(q[0]);
  }

  //! gets rotation angle (in degree [0,360]) and vector
  void Quaternion::getDeg(real& degree,Vector& vec) const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)){ degree=0.; vec.set(0.,0.,1.); return; }
    degree=acos(q[0]);
    real s=sin(degree);
    degree*=360./MT_PI;
    vec(0)=q[1]/s; vec(1)=q[2]/s; vec(2)=q[3]/s;
  }

  //! gets rotation angle (in rad [0,2pi]) and vector
  void Quaternion::getRad(real& angle,Vector& vec) const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)){ angle=0.; vec.set(0.,0.,1.); return; }
    angle=acos(q[0]);
    real s=sin(angle);
    angle*=2;
    vec(0)=q[1]/s; vec(1)=q[2]/s; vec(2)=q[3]/s;
    CHECK(angle>=0. && angle<=MT_2PI,"");
  }

  //! gets the axis rotation vector with length equal to the rotation angle in rad
  Vector& Quaternion::getVec(Vector& v) const{
    if(q[0]>=1. || q[0]<=-1. || (q[1]==0. && q[2]==0. && q[3]==0.)){ v.setZero(); return v; }
    real phi=acos(q[0]);
    real s=2.*phi/sin(phi);
    v(0)=s*q[1]; v(1)=s*q[2]; v(2)=s*q[3];
    return v;
  }

  Vector& Quaternion::getX(Vector& Rx) const{
    real q22 = 2.*q[2]*q[2];
    real q33 = 2.*q[3]*q[3];
    real q12 = 2.*q[1]*q[2];
    real q13 = 2.*q[1]*q[3];
    real q02 = 2.*q[0]*q[2];
    real q03 = 2.*q[0]*q[3];
    Rx.v[0]=1-q22-q33;
    Rx.v[1]=q12+q03;
    Rx.v[2]=q13-q02;
    return Rx;
  }
  Vector& Quaternion::getY(Vector& Ry) const{ Ry = (*this)*Vector(0,1,0);  return Ry; }
  Vector& Quaternion::getZ(Vector& Rz) const{ Rz = (*this)*Vector(0,0,1);  return Rz; }

  void Quaternion::setMatrix(real* m){
    q[0]=sqrt(1.-(3.-(m[0]+m[4]+m[8]))/4.);
    q[3] = (m[3]-m[1])/(4.*q[0]);
    q[2] = (m[2]-m[6])/(4.*q[0]);
    q[1] = (m[7]-m[5])/(4.*q[0]);
    normalize();
    //CHECK(normalized(),"failed :-(");
  }

  //! exports the rotation to a real[9] matrix, row-by-row
  real* Quaternion::getMatrix(real* m) const{
    real q11 = 2.*q[1]*q[1];
    real q22 = 2.*q[2]*q[2];
    real q33 = 2.*q[3]*q[3];
    real q12 = 2.*q[1]*q[2];
    real q13 = 2.*q[1]*q[3];
    real q23 = 2.*q[2]*q[3];
    real q01 = 2.*q[0]*q[1];
    real q02 = 2.*q[0]*q[2];
    real q03 = 2.*q[0]*q[3];
    m[0]=1-q22-q33; m[1]=q12-q03;    m[2]=q13+q02;
    m[3]=q12+q03;   m[4]=1-q11-q33;  m[5]=q23-q01;
    m[6]=q13-q02;   m[7]=q23+q01;    m[8]=1-q11-q22;
    return m;
  }

  //! exports the rotation to an ODE format matrix of type real[12]
  real* Quaternion::getMatrixOde(real* m) const{
    real q11 = 2.*q[1]*q[1];
    real q22 = 2.*q[2]*q[2];
    real q33 = 2.*q[3]*q[3];
    real q12 = 2.*q[1]*q[2];
    real q13 = 2.*q[1]*q[3];
    real q23 = 2.*q[2]*q[3];
    real q01 = 2.*q[0]*q[1];
    real q02 = 2.*q[0]*q[2];
    real q03 = 2.*q[0]*q[3];
    m[0]=1-q22-q33; m[1]=q12-q03;   m[2] =q13+q02;
    m[4]=q12+q03;   m[5]=1-q11-q33; m[6] =q23-q01;
    m[8]=q13-q02;   m[9]=q23+q01;   m[10]=1-q11-q22;
    m[3]=m[7]=m[11]=0.;
    return m;
  }
  //! exports the rotation to an OpenGL format matrix of type real[16]
  real* Quaternion::getMatrixGL(real* m) const{
    real q11 = 2.*q[1]*q[1];
    real q22 = 2.*q[2]*q[2];
    real q33 = 2.*q[3]*q[3];
    real q12 = 2.*q[1]*q[2];
    real q13 = 2.*q[1]*q[3];
    real q23 = 2.*q[2]*q[3];
    real q01 = 2.*q[0]*q[1];
    real q02 = 2.*q[0]*q[2];
    real q03 = 2.*q[0]*q[3];
    m[0]=1-q22-q33; m[4]=q12-q03;  m[8] =q13+q02;
    m[1]=q12+q03;   m[5]=1-q11-q33; m[9] =q23-q01;
    m[2]=q13-q02;   m[6]=q23+q01;  m[10]=1-q11-q22;
    m[3]=m[7]=m[11]=m[12]=m[13]=m[14]=0.;
    m[15]=1.;
    return m;
  }

  void Quaternion::writeNice(std::ostream& os) const{ Vector v; os <<"Quaternion: " <<getDeg() <<" around " <<getVec(v) <<"\n"; }
  void Quaternion::write(std::ostream& os) const{
    if(!MT::IOraw) os <<'(' <<q[0] <<' ' <<q[1] <<' ' <<q[2] <<' ' <<q[3] <<')';
    else os <<' ' <<q[0] <<' ' <<q[1] <<' ' <<q[2] <<' ' <<q[3];
  }
  void Quaternion::read(std::istream& is){ is >>"(" >>q[0] >>q[1] >>q[2]  >>q[3] >>")"; }
  //}

  //! compound of two rotations (A=B*C)
  Quaternion operator*(const Quaternion& b,const Quaternion& c){
    Quaternion a;
    a.q[0] = b.q[0]*c.q[0] - b.q[1]*c.q[1] - b.q[2]*c.q[2] - b.q[3]*c.q[3];
    a.q[1] = b.q[0]*c.q[1] + b.q[1]*c.q[0] + b.q[2]*c.q[3] - b.q[3]*c.q[2];
    a.q[2] = b.q[0]*c.q[2] + b.q[2]*c.q[0] + b.q[3]*c.q[1] - b.q[1]*c.q[3];
    a.q[3] = b.q[0]*c.q[3] + b.q[3]*c.q[0] + b.q[1]*c.q[2] - b.q[2]*c.q[1];
    return a;
  }

  //! A=B*C^{-1}
  Quaternion operator/(const Quaternion& b,const Quaternion& c){
    Quaternion a;
    a.q[0] =-b.q[0]*c.q[0] - b.q[1]*c.q[1] - b.q[2]*c.q[2] - b.q[3]*c.q[3];
    a.q[1] = b.q[0]*c.q[1] - b.q[1]*c.q[0] + b.q[2]*c.q[3] - b.q[3]*c.q[2];
    a.q[2] = b.q[0]*c.q[2] - b.q[2]*c.q[0] + b.q[3]*c.q[1] - b.q[1]*c.q[3];
    a.q[3] = b.q[0]*c.q[3] - b.q[3]*c.q[0] + b.q[1]*c.q[2] - b.q[2]*c.q[1];
    return a;
  }

  //! transform of a vector by a rotation
  Vector operator*(const Quaternion& b,const Vector& c){
    real m[9];
    b.getMatrix(m);
    Vector a;
    a.v[0]=m[0]*c.v[0]+m[1]*c.v[1]+m[2]*c.v[2];
    a.v[1]=m[3]*c.v[0]+m[4]*c.v[1]+m[5]*c.v[2];
    a.v[2]=m[6]*c.v[0]+m[7]*c.v[1]+m[8]*c.v[2];
    return a;
  }
  
  //! inverse transform of a vector by a rotation
  Vector operator/(const Quaternion& b,const Vector& c){
    real m[9];
    b.getMatrix(m);
    Vector a;
    a.v[0]=m[0]*c.v[0]+m[3]*c.v[1]+m[6]*c.v[2];
    a.v[1]=m[1]*c.v[0]+m[4]*c.v[1]+m[7]*c.v[2];
    a.v[2]=m[2]*c.v[0]+m[5]*c.v[1]+m[8]*c.v[2];
    return a;
  }

  //! transform of a vector by a frame
  Vector operator*(const Frame& X,const Vector& c){
    Vector a;
    a = X.r * c;
    a += X.p;
    return a;
  }
  
  //! inverse transform of a vector by a frame
  Vector operator/(const Frame& X,const Vector& c){
    Vector a(c);
    a -= X.p;
    a = X.r / a;
    return a;
  }

  //! constructor
  Frame::Frame(){ setZero(); }

  //! copy operator
  //Frame& Frame::operator=(const Frame& f){
  //  p=f.p; v=f.v; r=f.r; w=f.w; /*a=f.a; b=f.b; s=f.s;*/ return *this; }

  //! copy operator
  //Frame& Frame::set(const Frame& f){
    //p=f.p; v=f.v; r=f.r; w=f.w; /*a=f.a; b=f.b; s=f.s;*/ return *this; }

  //! initialize by reading from the string
  Frame& Frame::setText(const char* txt){ read(MT::String(txt)()); return *this; }

  //! initialize by reading from the stream
  //Frame& Frame::set(istream &is){ setZero(); read(is); return *this; }

  //! resets the position to origin, rotation to identity, velocities to zero, scale to unit
  void Frame::setZero(){
    p.setZero(); v.setZero(); r.setZero(); w.setZero(); /*a.setZero(); b.setZero(); s=1.;*/ }

  //! randomize the frame
  void Frame::setRandom(){
    p.setRandom(); v.setRandom(); r.setRandom(); w.setRandom(); /*a.setRandom(); b.setRandom(); s=rnd.uni();*/ }

  /*!\brief moves the frame according to the current velocities \c v and \c w
      and the time span \c time (if time is given in seconds, v has
      dimension units/sec, and w has dimension rad/sec) */
  /*void Frame::step(real time){
    Quaternion W;
    W.setVec(w);
    W.multiply(time);
    p+=v;
    r=W*r;
  }*/

  //! multiply the current scale by f
  /*void Frame::scale(real f){
    s*=f;
  }*/
  //! move the turtle by the vector (x,z,y) WITH RESPECT TO the current orientation/scale
  void Frame::addRelativeTranslation(real x,real y,real z){
    Vector X(x,y,z);
    //X=r*(s*X); //in global coords
    X=r*X; //in global coords
    p+=X;
    v+=w^X;
  }
  //! add a velocity to the turtle's inertial frame
  void Frame::addRelativeVelocity(real x,real y,real z){
    Vector X(x,y,z);
    //v+=r*(s*X);
    v+=r*X;
  }
  //! add an angular velocity to the turtle inertial frame
  void Frame::addRelativeAngVelocityDeg(real degree,real x,real y,real z){
    Vector W(x,y,z); W.normalize();
    W*=degree*MT_PI/180.;
    w+=r*W;
  }
  //! add an angular velocity to the turtle inertial frame
  void Frame::addRelativeAngVelocityRad(real rad,real x,real y,real z){
    Vector W(x,y,z); W.normalize();
    W*=rad;
    w+=r*W;
  }
  //! add an angular velocity to the turtle inertial frame
  void Frame::addRelativeAngVelocityRad(real wx,real wy,real wz){
    Vector W(wx,wy,wz);
    w+=r*W;
  }
  //! rotate the turtle orientation by an angle (given in DEGREE) around the vector (x,y,z) (given relative to the current orientation)
  void Frame::addRelativeRotationDeg(real degree,real x,real y,real z){
    Quaternion R;
    R.setDeg(degree,x,y,z);
    r=r*R;
  }
  //! rotate the turtle orientation by an angle (given in radiants) around the vector (x,y,z) (given relative to the current orientation)
  void Frame::addRelativeRotationRad(real rad,real x,real y,real z){
    Quaternion R;
    R.setRad(rad,x,y,z);
    r=r*R;
  }
  //! rotate the turtle orientation as given by a quaternion
  void Frame::addRelativeRotationQuat(real s,real x,real y,real z){
    Quaternion R;
    R.q[0]=s; R.q[1]=x; R.q[2]=y; R.q[3]=z;
    r=r*R;
  }
  /*!\brief transform the turtle into the frame f,
      which is interpreted RELATIVE to the current frame
      (new = f * old) */
  void Frame::addRelativeFrame(const Frame& f){
#ifdef ORS_NO_DYNAMICS_IN_FRAMES
      p += r*f.p;
      r = r*f.r;
#else
      //Vector P(r*(s*f.p)); //relative offset in global coords
      //Vector V(r*(s*f.v)); //relative vel in global coords
      Matrix R;
      r.getMatrix(R.m);
      Vector P(R*f.p); //relative offset in global coords
      Vector V(R*f.v); //relative vel in global coords
      Vector W(R*f.w); //relative ang vel in global coords
      p += P;
      v += w^P;
      v += V;
      //a += b^P;
      //a += w^((w^P) + 2.*V);
      //a += r*(s*f.a);
      //b += w^W;
      //b += r*f.b;
      w += W;
      r = r*f.r;
      //s*=f.s;
#endif
  }
  //! inverse transform (new = f^{-1} * old) or (old = f * new) 
  void Frame::subRelativeFrame(const Frame& f){
    //s/=f.s;
    r=r/f.r;
    w-=r*f.w;
    v-=r*f.v;
    //v-=r*(s*f.v);
    //Vector P(r*(s*f.p)); //frame offset in global coords
    Vector P(r*f.p); //frame offset in global coords
    v-=w^P;
    p-=P;
  }
  //! new = old * f
  void Frame::makeRelative(const Frame& f){
    Frame newf;
    newf=f;
    newf.addRelativeFrame(*this);
    *this=newf;
  }
  //! new = old * f^{-1}
  void Frame::makeRelativeToInv(const Frame& f){
    Frame newf;
    newf.setInverse(f);
    newf.addRelativeFrame(*this);
    *this=newf;
  }
  //! this = f^{-1}
  void Frame::setInverse(const Frame& f){ setZero(); subRelativeFrame(f); }

  //!  to = new * from
  void Frame::setDifference(const Frame& from,const Frame& to){
    //s=to.s/from.s;
    r=Quaternion()/from.r *to.r;
    w=from.r/(to.w-from.w);
    /*v=(1./from.s) * (from.r/(to.v-from.v));
    v-=(1./from.s) * (from.r/(from.w^(to.p-from.p)));
    p=(1./from.s) * (from.r/(to.p-from.p));*/
    v = from.r/(to.v-from.v);
    v-= from.r/(from.w^(to.p-from.p));
    p = from.r/(to.p-from.p);
  }

  //! get the current position/orientation/scale in an OpenGL format matrix (of type real[16])
  real* Frame::getAffineMatrix(real *m) const{
    real M[9]; r.getMatrix(M);
    m[0] =M[0]; m[1] =M[1]; m[2] =M[2]; m[3] =p(0);
    m[4] =M[3]; m[5] =M[4]; m[6] =M[5]; m[7] =p(1);
    m[8] =M[6]; m[9] =M[7]; m[10]=M[8]; m[11]=p(2);
    m[12]=0.;   m[13]=0.;   m[14]=0.;   m[15]=1.;
    return m;
  }

  //! get inverse OpenGL matrix for this frame (of type real[16])
  real* Frame::getInverseAffineMatrix(real *m) const{
    real M[9]; r.getMatrix(M);
    Vector pinv; pinv=r/p;
    m[0] =M[0]; m[1] =M[3]; m[2] =M[6]; m[3] =-pinv(0);
    m[4] =M[1]; m[5] =M[4]; m[6] =M[7]; m[7] =-pinv(1);
    m[8] =M[2]; m[9] =M[5]; m[10]=M[8]; m[11]=-pinv(2);
    m[12]=0.;   m[13]=0.;   m[14]=0.;   m[15]=1.;
    return m;
  }

  //! get the current position/orientation/scale in an OpenGL format matrix (of type real[16])
  real* Frame::getAffineMatrixGL(real *m) const{
    real M[9]; r.getMatrix(M);
    m[0]=M[0]; m[4]=M[1]; m[8] =M[2]; m[12]=p(0);
    m[1]=M[3]; m[5]=M[4]; m[9] =M[5]; m[13]=p(1);
    m[2]=M[6]; m[6]=M[7]; m[10]=M[8]; m[14]=p(2);
    m[3]=0.;   m[7]=0.;   m[11]=0.;   m[15]=1.;
    return m;
  }

  //! get inverse OpenGL matrix for this frame (of type real[16]) */
  real* Frame::getInverseAffineMatrixGL(real *m) const{
    real M[9]; r.getMatrix(M);
    Vector pinv; pinv=r/p;
    m[0]=M[0]; m[4]=M[3]; m[8] =M[6]; m[12]=-pinv(0);
    m[1]=M[1]; m[5]=M[4]; m[9] =M[7]; m[13]=-pinv(1);
    m[2]=M[2]; m[6]=M[5]; m[10]=M[8]; m[14]=-pinv(2);
    m[3]=0.;   m[7]=0.;   m[11]=0.;   m[15]=1.;
    return m;
  }

  //! operator<<
  void Frame::write(std::ostream& os) const{
    bool space=false;
    os <<"<";
    if(!p.isZero()){ os <<"t" <<p;  space=true; }
    if(!v.isZero()){ if(space) os <<' ';  os <<"v" <<v;  space=true; }
    if(!w.isZero()){ if(space) os <<' ';  os <<"w" <<w;  space=true; }
    if(!r.isZero()){ if(space) os <<' ';  os <<"q" <<r;  space=true; }
    //if(s!=1.) os <<" s(" <<s <<") ";
    os <<">";
  }
  //! operator>>
  void Frame::read(std::istream& is){
    setZero();
    char c;
    real x[4];
    MT::parse(is,"<");
    for(;;){
      is >>c;
      if(is.fail()) return; //EOF I guess
      //if(c==';') break;
      //if(c==',') is >>c;
      switch(c){
      //case '<': break; //do nothing -- assume this is an opening tag
      case 't': is>>"(">>x[0]>>x[1]>>x[2]>>")";	      addRelativeTranslation(x[0],x[1],x[2]); break;
      case 'q': is>>"(">>x[0]>>x[1]>>x[2]>>x[3]>>")"; addRelativeRotationQuat(x[0],x[1],x[2],x[3]); break;
      case 'r': is>>"(">>x[0]>>x[1]>>x[2]>>x[3]>>")"; addRelativeRotationRad(x[0],x[1],x[2],x[3]); break;
      case 'd': is>>"(">>x[0]>>x[1]>>x[2]>>x[3]>>")"; addRelativeRotationDeg(x[0],x[1],x[2],x[3]); break;
      case 'v': is>>"(">>x[0]>>x[1]>>x[2]>>")";	      addRelativeVelocity(x[0],x[1],x[2]); break;
      case 'w': is>>"(">>x[0]>>x[1]>>x[2]>>")";       addRelativeAngVelocityRad(x[0],x[1],x[2]); break;
      //case 's': is>>"(">>x[0]>>")";                   scale(x[0]); break;
      case '|': is.putback('<'); return;
      case '>': return; //those symbols finish the reading without error
      default: MT_MSG("unknown Frame read tag: " <<c <<"abort reading this frame"); is.putback(c); return;
      }
      if(is.fail()) HALT("error reading '"<<c<<"' parameters in frame");
    }
    if(is.fail()) HALT ("could not read Frame struct");
  }

  Mesh::Mesh(){ }


  void Mesh::makeVerticesRelativeToGroup(){
    uint i;
    int g;
    Vector *v;
    for(i=0;i<V.d0;i++) if((g=G(i))!=-1){
      v = (Vector*)&V(i,0);
      *v = GF(g)->r/((*v) - GF(g)->p);
      v = (Vector*)&Vn(i,0);
      *v = GF(g)->r/(*v);
    }
  }

  void Mesh::collectTriGroups(){
    uint i;
    int g;
    GT.resize(GF.N+1);
    for(i=0;i<T.d0;i++){
      g=G(T(i,0));
      if(g!=-1 && g==G(T(i,1)) && g==G(T(i,2))){
        GT(g).append(i);
      }else{
        GT(GF.N).append(i);
      }
    }
  }
}

//! use as similarity measure (distance = 1 - |scalarprod|)
real scalarProduct(const ors::Quaternion& a,const ors::Quaternion& b){
  return a.q[0]*b.q[0]+a.q[1]*b.q[1]+a.q[2]*b.q[2]+a.q[3]*b.q[3];
}

std::istream& operator>>(std::istream& is,ors::Vector& x)    { x.read(is); return is; }
std::istream& operator>>(std::istream& is,ors::Matrix& x)    { x.read(is); return is; }
std::istream& operator>>(std::istream& is,ors::Quaternion& x){ x.read(is); return is; }
std::istream& operator>>(std::istream& is,ors::Frame& x)     { x.read(is); return is; }
std::ostream& operator<<(std::ostream& os,const ors::Vector& x)    { x.write(os); return os; }
std::ostream& operator<<(std::ostream& os,const ors::Matrix& x)    { x.write(os); return os; }
std::ostream& operator<<(std::ostream& os,const ors::Quaternion& x){ x.write(os); return os; }
std::ostream& operator<<(std::ostream& os,const ors::Frame& x)     { x.write(os); return os; }


//================================================================================
//
// Mesh code
//

void ors::Mesh::clear(){
  V.clear(); Vn.clear(); T.clear(); Tn.clear(); C.clear(); //strips.clear();
}

void ors::Mesh::setBox(){
  real verts[24] = {
    -.5, -.5, -.5 ,
    +.5, -.5, -.5 ,
    +.5, +.5, -.5 ,
    -.5, +.5, -.5 ,
    -.5, -.5, +.5 ,
    +.5, -.5, +.5 ,
    +.5, +.5, +.5 ,
    -.5, +.5, +.5 };
  uint   tris [36] = {
    0, 3, 2, 2, 1, 0,
    4, 5, 6, 6, 7, 4,
    1, 5, 4, 4, 0, 1,
    2, 6, 5, 5, 1, 2,
    3, 7, 6, 6, 2, 3,
    0, 4, 7, 7, 3, 0 };
  V.setCarray(verts,24);
  T.setCarray(tris ,36);
  V.reshape(8,3);
  T.reshape(12,3);
  //cout <<V <<endl;  for(uint i=0;i<4;i++) cout <<norm(V[i]) <<endl;
}

void ors::Mesh::setTetrahedron(){
  real s2=MT_SQRT2/3.,s6=sqrt(6.)/3.;
  real verts[12] = { 0., 0., 1. , 2.*s2, 0., -1./3., -s2, s6, -1./3., -s2, -s6, -1./3. };
  uint   tris [12] = { 0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2 };
  V.setCarray(verts,12);
  T.setCarray(tris ,12);
  V.reshape(4,3);
  T.reshape(4,3);
  //cout <<V <<endl;  for(uint i=0;i<4;i++) cout <<norm(V[i]) <<endl;
}

void ors::Mesh::setOctahedron(){
  real verts[18] = {
    1, 0, 0,
   -1, 0, 0,
    0, 1, 0,
    0,-1, 0,
    0, 0, 1,
    0, 0,-1
  };
  uint   tris [24] = {
    4, 0, 2,  4, 2, 1,  4, 1, 3,  4, 3, 0,
    5, 2, 0,  5, 1, 2,  5, 3, 1,  5, 0, 3
  };
  V.setCarray(verts,18);
  T.setCarray(tris ,24);
  V.reshape(6,3);
  T.reshape(8,3);
  //cout <<V <<endl;  for(uint i=0;i<4;i++) cout <<norm(V[i]) <<endl;
}

void ors::Mesh::setDodecahedron(){
  real a = 1/sqrt(3.), b = sqrt((3.-sqrt(5.))/6.), c=sqrt((3.+sqrt(5.))/6.);
  real verts[60] = {
    a, a, a,
      a, a,-a,
      a,-a, a,
      a,-a,-a,
      -a, a, a,
      -a, a,-a,
      -a,-a, a,
      -a,-a,-a,
      b, c, 0,
      -b, c, 0,
      b,-c, 0,
      -b,-c, 0,
      c, 0, b,
      c, 0,-b, 
      -c, 0, b, 
      -c, 0,-b,
      0, b, c,
      0,-b, c,
      0, b,-c,
      0,-b,-c
  };
  uint tris [108] = {
    0, 8, 9, 0, 9, 4, 0, 4, 16, 0, 12, 13, 0, 13, 1, 0, 1, 8,
      0, 16, 17, 0, 17, 2, 0, 2, 12, 8, 1, 18, 8, 18, 5, 8, 5, 9,
      12, 2, 10, 12, 10, 3, 12, 3, 13, 16, 4, 14, 16, 14, 6, 16, 6, 17,
      9, 5, 15, 9, 15, 14, 9, 14, 4, 6, 11, 10, 6, 10, 2, 6, 2, 17,
      3, 19, 18, 3, 18, 1, 3, 1, 13, 7, 15, 5, 7, 5, 18, 7, 18, 19,
      7, 11, 6, 7, 6, 14, 7, 14, 15, 7, 19, 3, 7, 3, 10, 7, 10, 11 
  };
  V.setCarray(verts,60);
  T.setCarray(tris ,108);
  V.reshape(20,3);
  T.reshape(36,3);
}

void ors::Mesh::setSphere(uint fineness){
  setOctahedron();
  for(uint k=0;k<fineness;k++){
    subDevide();
    for(uint i=0;i<V.d0;i++) V[i]() /= norm(V[i]);
  }
}

void ors::Mesh::setHalfSphere(uint fineness){
  setOctahedron();
  V.resizeCopy(5,3);
  T.resizeCopy(4,3);
  for(uint k=0;k<fineness;k++){
    subDevide();
    for(uint i=0;i<V.d0;i++) V[i]() /= norm(V[i]);
  }
}

void ors::Mesh::setCylinder(real r,real l,uint fineness){
  uint div = 4 * (1<<fineness);
  V.resize(2*div+2,3);
  T.resize(4*div,3);
  uint i,j;
  real phi;
  for(i=0;i<div;i++){ //vertices
    phi=MT_2PI*i/div;
    V(i,0)=r*::cos(phi);
    V(i,1)=r*::sin(phi);
    V(i,2)=.5*l;
    V(i+div,0)=V(i,0);
    V(i+div,1)=V(i,1);
    V(i+div,2)=-.5*l;
  }
  V(2*div+0,0)=V(2*div+0,1)=.0;  V(2*div+0,2)=+.5*l; //upper center
  V(2*div+1,0)=V(2*div+1,1)=.0;  V(2*div+1,2)=-.5*l; //lower center
  for(i=0;i<div;i++){ //triangles
    j=(i+1)%div;
    T(4*i  ,0)=i;
    T(4*i  ,1)=j+div;
    T(4*i  ,2)=j;
    
    T(4*i+2,0)=i;
    T(4*i+2,1)=j;
    T(4*i+2,2)=2*div+0;
    
    T(4*i+1,0)=i;
    T(4*i+1,1)=i+div;
    T(4*i+1,2)=j+div;

    T(4*i+3,0)=j+div;
    T(4*i+3,1)=i+div;
    T(4*i+3,2)=2*div+1;
  }
}

void ors::Mesh::setCappedCylinder(real r,real l,uint fineness){
#if 0
  setCylinder(m,r,l,fineness);
  ors::Mesh cap;
  setHalfSphere(cap,fineness);
  scale(cap,-r,r,r);
  MeshTranslate(cap,0,0,.5*l);
  MeshAddMesh(m,cap);
  scale(cap,-1,1,-1);
  MeshAddMesh(m,cap);
  fuseNearVertices(1e-5);
#else
  uint i;
  setSphere(fineness);
  scale(r,r,r);
  for(i=0;i<V.d0;i++) if(V(i,2)>0.) V(i,2)+=l;
  translate(0,0,-.5*l);
#endif
}

/*!\brief add triangles according to the given grid; grid has to be a 2D
  Array, the elements of which are indices referring to vertices in
  the vertex list (V) */
void ors::Mesh::setGrid(uint X,uint Y){
  CHECK(X>1 && Y>1,"grid has to be at least 2x2");
  CHECK(V.d0==X*Y,"don't have X*Y mesh-vertices to create grid faces");
  uint i,j,k=T.d0;
  T.resizeCopy(k+(Y-1)*2*(X-1),3);
  for(j=0;j<Y-1;j++){
    for(i=0;i<X-1;i++){
      T(k,0)=j*X+i; T(k,1)=(j+1)*X+i; T(k,2)=(j+1)*X+(i+1);
      k++;
      T(k,0)=j*X+i; T(k,1)=(j+1)*X+(i+1); T(k,2)=j*X+(i+1);
      k++;
    }
  }
}

#ifndef MT_Lewiner
void ors::Mesh::setImplicitSurface(real(*fct)(real,real,real),real lo,real hi,uint res){
  NIY;
}
#else
#  include "ors_lewiner.cpp"
#endif

void ors::Mesh::subDevide(){
  uint v=V.d0,t=T.d0;
  V.resizeCopy(v+3*t,3);
  uintA newT(4*t,3);
  uint a,b,c,i,k,l;
  for(i=0,k=v,l=0;i<t;i++){
    a=T(i,0); b=T(i,1); c=T(i,2);
    V[k+0]() = (real).5*(V[a] + V[b]);
    V[k+1]() = (real).5*(V[b] + V[c]);
    V[k+2]() = (real).5*(V[c] + V[a]);
    newT(l,0)=a;   newT(l,1)=k+0; newT(l,2)=k+2; l++;
    newT(l,0)=k+0; newT(l,1)=b;   newT(l,2)=k+1; l++;
    newT(l,0)=k+0; newT(l,1)=k+1; newT(l,2)=k+2; l++;
    newT(l,0)=k+2; newT(l,1)=k+1; newT(l,2)=c;   l++;
    k+=3;
  }
  T = newT;
}

void ors::Mesh::scale(real f){  V *= f; }

void ors::Mesh::scale(real sx,real sy,real sz){
  uint i;
  for(i=0;i<V.d0;i++){  V(i,0)*=sx;  V(i,1)*=sy;  V(i,2)*=sz;  }
}

void ors::Mesh::translate(real dx,real dy,real dz){
  uint i;
  for(i=0;i<V.d0;i++){  V(i,0)+=dx;  V(i,1)+=dy;  V(i,2)+=dz;  }
}

void ors::Mesh::center(){
  arr mean(3);
  mean.setZero();
  uint i;
  for(i=0;i<V.d0;i++) mean += V[i];
  mean /= (real)V.d0;
  for(i=0;i<V.d0;i++) V[i]() -= mean;
}

void ors::Mesh::box(){
  real x,X,y,Y,z,Z,m;
  x=X=V(0,0);
  y=Y=V(0,1);
  z=Z=V(0,2);
  for(uint i=0;i<V.d0;i++){
    if(V(i,0)<x) x=V(i,0);
    if(V(i,0)>X) X=V(i,0);
    if(V(i,1)<y) y=V(i,1);
    if(V(i,1)>Y) Y=V(i,1);
    if(V(i,2)<z) z=V(i,2);
    if(V(i,2)>Z) Z=V(i,2);
  }
  translate(-.5*(x+X),-.5*(y+Y),-.5*(z+Z));
  m=X-x;
  if(Y-y>m) m=Y-y;
  if(Z-z>m) m=Z-z;
  scale(1./m);
}

void ors::Mesh::addMesh(const ors::Mesh& mesh2){
  uint n=V.d0,t=T.d0;
  V.append(mesh2.V);
  T.append(mesh2.T);
  for(;t<T.d0;t++){  T(t,0)+=n;  T(t,1)+=n;  T(t,2)+=n;  }
}

/*!\brief calculate the normals of all triangles (Tn) and the average
  normals of the vertices (N); average normals are averaged over
  all adjacent triangles that are in the triangle list or member of
  a strip */
void ors::Mesh::computeNormals(){
  uint i;
  Vector a,b,c;
  Tn.resize(T.d0,3);
  Tn.setZero();
  Vn.resize(V.d0,3);
  Vn.setZero();
  //triangle normals and contributions
  for(i=0;i<T.d0;i++){      
    a.set(&V(T(i,0),0)); b.set(&V(T(i,1),0)); c.set(&V(T(i,2),0));
    b-=a; c-=a; a=b^c; a.normalize();
    Tn(i,0)=a(0);  Tn(i,1)=a(1);  Tn(i,2)=a(2);
    Vn(T(i,0),0)+=a(0);  Vn(T(i,0),1)+=a(1);  Vn(T(i,0),2)+=a(2);
    Vn(T(i,1),0)+=a(0);  Vn(T(i,1),1)+=a(1);  Vn(T(i,1),2)+=a(2);
    Vn(T(i,2),0)+=a(0);  Vn(T(i,2),1)+=a(1);  Vn(T(i,2),2)+=a(2);
  }
  Vector *d;
  for(i=0;i<Vn.d0;i++){ d=(Vector*)&Vn(i,0); d->normalize(); }
}
  
/*!\brief add triangles according to the given grid; grid has to be a 2D
  Array, the elements of which are indices referring to vertices in
  the vertex list (V) */
/*void ors::Mesh::gridToTriangles(const uintA &grid){
  uint i,j,k=T.d0;
  T.resizeCopy(T.d0+2*(grid.d0-1)*(grid.d1-1),3);
  for(i=0;i<grid.d0-1;i++) for(j=0;j<grid.d1-1;j++){
    if((i+j)&1){
      T(k,0)=grid(i+1,j  );
      T(k,1)=grid(i  ,j  );
      T(k,2)=grid(i  ,j+1);
      k++;
      T(k,0)=grid(i+1,j  );
      T(k,1)=grid(i  ,j+1);
      T(k,2)=grid(i+1,j+1);
      k++;
    }else{
      T(k,0)=grid(i+1,j  );
      T(k,1)=grid(i  ,j  );
      T(k,2)=grid(i+1,j+1);
      k++;
      T(k,0)=grid(i+1,j+1);
      T(k,1)=grid(i  ,j  );
      T(k,2)=grid(i  ,j+1);
      k++;
    }
  }
}*/
  
/*!\brief add strips according to the given grid (sliced in strips along
  the x-axis (the first index)); grid has to be a 2D Array, the
  elements of which are indices referring to vertices in the vertex
  list (V) */
/*void ors::Mesh::gridToStrips(const uintA& grid){
  CHECK(grid.d0>1 && grid.d1>1,"grid has to be at least 2x2");
  uint i,j,k=strips.N,l;
  strips.resizeCopy(strips.N+grid.d0-1);
  for(i=0;i<grid.d0-1;i++){
    strips(k).resize(2*grid.d1);
    l=0;
    for(j=0;j<grid.d1;j++){
      strips(k)(l)=grid(i+1,j); l++;
      strips(k)(l)=grid(i  ,j); l++;
    }
#if 0 //code to make it less symmetric
      //}else{
    strips(k)(l)=grid(i,0); l++;
    for(j=0;j<grid.d1;j++){
      strips(k)(l)=grid(i  ,j); l++;
      strips(k)(l)=grid(i+1,j); l++;
    }
#endif
    k++;
  }
}*/
  
/*!\brief add strips according to the given grid (sliced in strips along
  the x-axis (the first index)); it is assumed that the vertices in
  the list V linearly correspond to points in the XxY grid */
/*void ors::Mesh::gridToStrips(uint X,uint Y){
  CHECK(X>1 && Y>1,"grid has to be at least 2x2");
  uint i,j,k=strips.N,l;
  strips.resizeCopy(strips.N+Y-1);
  for(j=0;j<Y-1;j++){
    strips(k).resize(2*X);
    l=0;
    for(i=0;i<X;i++){
      strips(k)(l)=(j+1)*X+i;
      l++;
      strips(k)(l)=    j*X+i;
      l++;
    }
    k++;
  }
}*/
  
void deleteZeroTriangles(ors::Mesh& m){
  uintA newT;
  newT.resizeAs(m.T);
  uint i,j;
  for(i=0,j=0;i<m.T.d0;i++){
    if(m.T(i,0)!=m.T(i,1) && m.T(i,0)!=m.T(i,2) && m.T(i,1)!=m.T(i,2))
      memmove(&newT(j++,0), &m.T(i,0), 3*newT.sizeT);
  }
  newT.resizeCopy(j,3);
  m.T=newT;
}
  
void permuteVertices(ors::Mesh& m,uintA& p){
  CHECK(p.N==m.V.d0,"");
  uint i;
  arr x(p.N,3);
  for(i=0;i<p.N;i++){ x(i,0)=m.V(p(i),0); x(i,1)=m.V(p(i),1); x(i,2)=m.V(p(i),2); }
  m.V=x;
  if(m.Vn.N){
    for(i=0;i<p.N;i++){ x(i,0)=m.Vn(p(i),0); x(i,1)=m.Vn(p(i),1); x(i,2)=m.Vn(p(i),2); }
    m.Vn=x;
  }
  if(m.C.N){
    for(i=0;i<p.N;i++){ x(i,0)=m.C(p(i),0); x(i,1)=m.C(p(i),1); x(i,2)=m.C(p(i),2); }
    m.C=x;
  }
  uintA y(m.T.d0,3);
  uintA p2(p.N); //inverse permutation
  for(i=0;i<p.N;i++) p2(p(i))=i;
  for(i=0;i<m.T.d0;i++){ y(i,0)=p2(m.T(i,0)); y(i,1)=p2(m.T(i,1)); y(i,2)=p2(m.T(i,2)); }
  m.T=y;
}

/*!\brief delete all void triangles (with vertex indices (0,0,0)) and void
  vertices (not used for triangles or strips) */
void ors::Mesh::deleteUnusedVertices(){
  if(!V.N) return;
  uintA p;
  uintA u;
  uint i,Nused;

  deleteZeroTriangles(*this);

  //count vertex usage
  u.resize(V.d0);
  u.setZero();
  for(i=0;i<T.d0;i++){ u(T(i,0))++; u(T(i,1))++; u(T(i,2))++; }
  //for(i=0;i<strips.N;i++) for(j=0;j<strips(i).N;j++) u(strips(i)(j))=true;

  //find proper permutation of vertex list
  p.setStraightPerm(V.d0);
  Nused=p.N;
  for(i=0;i<Nused;i++) if(!u(i)){ Nused--; p.permute(i,Nused); u.permute(i,Nused); i--; }

  permuteVertices(*this,p);
  V.resizeCopy(Nused,3);
}

arr *COMP_V;
bool COMP(uint i,uint j){
  bool r=(*COMP_V)[i]<(*COMP_V)[j];
  return r;
}

/*!\brief delete all void triangles (with vertex indices (0,0,0)) and void
  vertices (not used for triangles or strips) */
void ors::Mesh::fuseNearVertices(real tol){
  if(!V.N) return;
  uintA p;
  uint i,j;

  cout <<"fusing vertices: #V=" <<V.d0 <<", sorting.." <<std::flush;
  //cout <<V <<endl;
  //sort vertices lexically
  p.setStraightPerm(V.d0);
  COMP_V=&V;
  std::sort(p.p,p.pstop,COMP);
  permuteVertices(*this,p);

  cout <<"permuting.." <<std::flush;
  //cout <<V <<endl;
  p.setStraightPerm(V.d0);
  for(i=0;i<V.d0;i++){
    if(p(i)!=i) continue; //i has already been fused with p(i), and p(i) has already been checked...
    for(j=i+1;j<V.d0;j++){
      if(V(j,0)-V(i,0)>tol) break;
      if(MT::sqr(V(j,0)-V(i,0))+MT::sqr(V(j,1)-V(i,1))+MT::sqr(V(j,2)-V(i,2))<tol*tol){
	//cout <<"fusing " <<i <<" " <<j <<" " <<V[i] <<" " <<V[j] <<endl;
	p(j)=i;
      }
    }
  }

  uintA y(T.d0,3);
  for(i=0;i<T.d0;i++){ y(i,0)=p(T(i,0)); y(i,1)=p(T(i,1)); y(i,2)=p(T(i,2)); }
  T=y;

  cout <<"deleting tris.." <<std::flush;
  deleteZeroTriangles(*this);

  cout <<"deleting verts.." <<std::flush;
  deleteUnusedVertices();

  cout <<"#V=" <<V.d0 <<", done" <<endl;
}

void getVertexNeighorsList(const ors::Mesh& m,intA& Vt,intA& VT){
  uint i,j;
  Vt.resize(m.V.d0);  Vt.setZero();
  VT.resize(m.V.d0,100);
  for(i=0;i<m.T.d0;i++){
    j=m.T(i,0);  VT(j,Vt(j))=i;  Vt(j)++;
    j=m.T(i,1);  VT(j,Vt(j))=i;  Vt(j)++;
    j=m.T(i,2);  VT(j,Vt(j))=i;  Vt(j)++;
  }
}

void getTriNormals(const ors::Mesh& m,arr& Tn){
  uint i;
  ors::Vector a,b,c;
  Tn.resize(m.T.d0,3); //tri normals
  for(i=0;i<m.T.d0;i++){
    a.set(&m.V(m.T(i,0),0)); b.set(&m.V(m.T(i,1),0)); c.set(&m.V(m.T(i,2),0));
    b-=a; c-=a; a=b^c; a.normalize();
    Tn(i,0)=a(0);  Tn(i,1)=a(1);  Tn(i,2)=a(2);
  }
}

void ors::Mesh::flipFaces(){
  uint i,a;
  for(i=0;i<T.d0;i++){
    a=T(i,0);
    T(i,0)=T(i,1);
    T(i,1)=a;
  }
}

void ors::Mesh::clean(){
  uint i,j,idist=0;
  Vector a,b,c,m;
  real mdist=0.;
  arr Tc(T.d0,3); //tri centers
  arr Tn(T.d0,3); //tri normals
  uintA Vt(V.d0);
  intA VT(V.d0,100); //tri-neighbors to a vertex
  Vt.setZero(); VT=-1;

  for(i=0;i<T.d0;i++){
    a.set(&V(T(i,0),0)); b.set(&V(T(i,1),0)); c.set(&V(T(i,2),0));

    //tri center
    m=(a+b+c)/3.;
    Tc(i,0)=m(0);  Tc(i,1)=m(1);  Tc(i,2)=m(2);

    //farthest tri
    if(m.length()>mdist){ mdist=m.length(); idist=i; }

    //tri normal
    b-=a; c-=a; a=b^c; a.normalize();
    Tn(i,0)=a(0);  Tn(i,1)=a(1);  Tn(i,2)=a(2);

    //vertex neighbor count
    j=T(i,0);  VT(j,Vt(j))=i;  Vt(j)++;
    j=T(i,1);  VT(j,Vt(j))=i;  Vt(j)++;
    j=T(i,2);  VT(j,Vt(j))=i;  Vt(j)++;
  }

  //step through tri list and flip them if necessary
  boolA Tisok(T.d0); Tisok=false;
  uintA Tok; //contains the list of all tris that are ok oriented
  uintA Tnew(T.d0,T.d1);
  Tok.append(idist);
  Tisok(idist)=true;
  int A=0,B=0,C=0,D;
  uint r,k,l;
  intA neighbors;
  for(k=0;k<Tok.N;k++){
    i=Tok(k);
    Tnew(k,0)=T(i,0); Tnew(k,1)=T(i,1); Tnew(k,2)=T(i,2); 

    for(r=0;r<3;r++){
      if(r==0){ A=T(i,0);  B=T(i,1);  C=T(i,2); }
      if(r==1){ A=T(i,1);  B=T(i,2);  C=T(i,0); }
      if(r==2){ A=T(i,2);  B=T(i,0);  C=T(i,1); }

      //check all triangles that share A & B
      setSection(neighbors,VT[A],VT[B]);
      neighbors.removeAllValues(-1);
      if(neighbors.N>2) MT_MSG("edge shared by more than 2 triangles " <<neighbors);
      neighbors.removeValue(i);
      //if(!neighbors.N) cout <<"mesh.clean warning: edge has only one triangle that shares it" <<endl;

      //orient them correctly
      for(l=0;l<neighbors.N;l++){
        j=neighbors(l); //j is a neighboring triangle sharing A & B
        D=-1;
        //align the neighboring triangle and let D be its 3rd vertex
        if((int)T(j,0)==A && (int)T(j,1)==B) D=T(j,2);
        if((int)T(j,0)==A && (int)T(j,2)==B) D=T(j,1);
        if((int)T(j,1)==A && (int)T(j,2)==B) D=T(j,0);
        if((int)T(j,1)==A && (int)T(j,0)==B) D=T(j,2);
        if((int)T(j,2)==A && (int)T(j,0)==B) D=T(j,1);
        if((int)T(j,2)==A && (int)T(j,1)==B) D=T(j,0);
        if(D==-1) HALT("dammit");
      	//determine orientation
        if(!Tisok(j)){
          T(j,0)=B;  T(j,1)=A;  T(j,2)=D;
          Tok.append(j);
          Tisok(j)=true;
        }else{
	        //check if consistent!
        }
      }

#if 0
      //compute their rotation
      if(neighbors.N>1){
        real phi,phimax;
        int jmax=-1;
        Vector ni, nj;
        for(l=0;l<neighbors.N;l++){
          j=neighbors(l); //j is a neighboring triangle sharing A & B

          a.set(&V(T(i,0),0)); b.set(&V(T(i,1),0)); c.set(&V(T(i,2),0));
          b-=a; c-=a; a=b^c; a.normalize();
          ni = a;

          a.set(&V(T(j,0),0)); b.set(&V(T(j,1),0)); c.set(&V(T(j,2),0));
          b-=a; c-=a; a=b^c; a.normalize();
          nj = a;

          Quaternion q;
          q.setDiff(ni,-nj);
          q.getDeg(phi,c);
          a.set(&V(A,0)); b.set(&V(B,0));
          if(c*(a-b) < 0.) phi+=180.;

          if(jmax==-1 || phi>phimax){ jmax=j; phimax=phi; }
        }
        if(!Tisok(jmax)){
          Tok.append(jmax);
          Tisok(jmax)=true;
        }
      }else{
        j = neighbors(0);
        if(!Tisok(j)){
          Tok.append(j);
          Tisok(j)=true;
        }
      }
#endif
    }
  }
  if(k<T.d0){
    cout <<"mesh.clean warning: not all triangles connected: "<<k <<"<" <<T.d0 <<endl;
    cout <<"WARNING: cutting of all non-connected triangles!!" <<endl;
    Tnew.resizeCopy(k,3);
    T=Tnew;
    deleteUnusedVertices();
  }
  computeNormals();
}

void getEdgeNeighborsList(const ors::Mesh& m,uintA& EV,uintA& Et,intA& ET){
  intA Vt,VT;
  getVertexNeighorsList(m,Vt,VT);

  uint A=0,B=0,C=0,t,tt,i,r,k;
  //build edge list
  EV.resize(m.T.d0*3,2);   EV=0;     //edge vert neighbors
  ET.resize(m.T.d0*3,10);  ET=-1;    //edge tri neighbors
  Et.resize(m.T.d0*3); Et.setZero(); //#edge tri neighbors
  boolA done(m.T.d0); done=false;
  for(t=0,k=0;t<m.T.d0;t++){
    for(r=0;r<3;r++){
      if(r==0){ A=m.T(t,0);  B=m.T(t,1);  C=m.T(t,2); }
      if(r==1){ A=m.T(t,1);  B=m.T(t,2);  C=m.T(t,0); }
      if(r==2){ A=m.T(t,2);  B=m.T(t,0);  C=m.T(t,1); }

      //has AB already been taken care of?
      bool yes=false;
      for(i=0;i<(uint)Vt(A);i++){
        tt=VT(A,i);
	if(m.T(tt,0)==B || m.T(tt,1)==B || m.T(tt,2)==B){
	  if(done(tt)) yes=true;
	}
      }
      if(yes) continue;

      //if not, then do it
      EV(k,0)=A;
      EV(k,1)=B;
      for(i=0;i<(uint)Vt(A);i++){
        tt=VT(A,i);
	if(m.T(tt,0)==B || m.T(tt,1)==B || m.T(tt,2)==B){
	  ET(k,Et(k))=tt;
	  Et(k)++;
	}
      }
      k++;
    }
    done(t)=true;
  }

  EV.resizeCopy(k,2);
  ET.resizeCopy(k,10);
  Et.resizeCopy(k);

  cout <<"\n#edges=" <<k
    <<"\nedge=\n" <<EV
    <<"\n@neighs=\n" <<Et
    <<"\nneighs=\n" <<ET <<endl;
}

void getTriNeighborsList(const ors::Mesh& m,uintA& Tt,intA& TT){
  intA Vt,VT;
  getVertexNeighorsList(m,Vt,VT);

  uint A=0,B=0,C=0,t,tt,r,i;
  Tt.resize(m.T.d0,3);     Tt.setZero();
  TT.resize(m.T.d0,3,100); TT=-1;
  for(t=0;t<m.T.d0;t++){
    for(r=0;r<3;r++){
      if(r==0){ A=m.T(t,0);  B=m.T(t,1);  C=m.T(t,2); }
      if(r==1){ A=m.T(t,1);  B=m.T(t,2);  C=m.T(t,0); }
      if(r==2){ A=m.T(t,2);  B=m.T(t,0);  C=m.T(t,1); }

      for(i=0;i<(uint)Vt(A);i++){
        tt=VT(A,i);
	if(tt!=t && (m.T(tt,0)==B || m.T(tt,1)==B || m.T(tt,2)==B)){
  	  TT(t,r,Tt(t,r))=tt;
	  Tt(t,r)++;
	}
      }
    }
  }

  //cout <<Tt <<TT <<endl;
}

void ors::Mesh::skin(uint start){
  intA TT;
  uintA Tt;
  getTriNeighborsList(*this,Tt,TT);
  arr Tn;
  getTriNormals(*this,Tn);

  uintA goodTris;
  boolA added(T.d0);
  goodTris.append(start);
  added=false;
  added(start)=true;
  uint t,tt,r,i,k;
  int m;
  real p,mp=0;
  for(k=0;k<goodTris.N;k++){
    t=goodTris(k);
    for(r=0;r<3;r++){
      //select from all neighbors the one most parallel
      m=-1;
      for(i=0;i<Tt(t,r);i++){
        tt=TT(t,r,i);
        p=scalarProduct(Tn[t],Tn[tt]);
        if(m==-1 || p>mp){ m=tt; mp=p; }
      }
      if(m!=-1 && !added(m)){ goodTris.append(m); added(m)=true; }
    }
  }

  uintA Tnew(k,3);
  for(k=0;k<goodTris.N;k++){
    t=goodTris(k);
    Tnew(k,0)=T(t,0); Tnew(k,1)=T(t,1); Tnew(k,2)=T(t,2); 
  }
  T=Tnew;
  cout <<T <<endl;
}

void ors::Mesh::readFile(const char* filename){
  bool loaded=false;
  const char *type = filename+(strlen(filename)-3);
  //cout <<"reading mesh file '" <<filename <<"' of type '" <<type<< "'" <<endl;
  if(!strcmp(type,"obj")){ readObjFile(filename); loaded=true; }
  if(!strcmp(type,"off")){ readOffFile(filename); loaded=true; }
  if(!strcmp(type,"ply")){ readPlyFile(filename); loaded=true; }
  if(!strcmp(type,"tri")){ readTriFile(filename); loaded=true; }
  if(!strcmp(type,"stl")){ readStlFile(filename); loaded=true; }
  if(!loaded) HALT("can't read file type '"<<type<<"'");
}

void ors::Mesh::writeTriFile(const char* filename){
  bool tmp=MT::IOraw;
  MT::IOraw=true;
  ofstream os;
  MT::open(os,filename);
  os <<"TRI" <<endl <<endl
     <<V.d0 <<endl
     <<T.d0 <<endl <<endl
     <<V <<endl <<endl
     <<T <<endl;
  MT::IOraw=tmp;
}

void ors::Mesh::readTriFile(const char* filename){
  ifstream is;
  MT::open(is,filename);
  uint i,nV,nT;
  is >>(const char*)"TRI" >>nV >>nT;
  V.resize(nV,3);
  T.resize(nT,3);
  for(i=0;i<V.N;i++) is >>V.elem(i);
  for(i=0;i<T.N;i++) is >>T.elem(i);
}

void ors::Mesh::writeOffFile(const char* filename){
  ofstream os;
  MT::open(os,filename);
  uint i;
  os <<"OFF\n" <<V.d0 <<' ' <<T.d0 <<' ' <<0 <<endl;
  for(i=0;i<V.d0;i++) os <<V(i,0) <<' ' <<V(i,1) <<' ' <<V(i,2) <<endl;
  for(i=0;i<T.d0;i++) os <<3 <<' ' <<T(i,0) <<' ' <<T(i,1) <<' ' <<T(i,2) <<endl;
}

void ors::Mesh::readOffFile(const char* filename){
  ifstream is;
  MT::open(is,filename);
  uint i,k,nVertices,nFaces,nEdges;
  is >>(const char*)"OFF" >>nVertices >>nFaces >>nEdges;
  CHECK(!nEdges,"can't read edges in off file");
  V.resize(nVertices,3);
  T.resize(nFaces   ,3);
  for(i=0;i<V.N;i++) is >>V.elem(i);
  for(i=0;i<T.d0;i++){
    is >>k;
    CHECK(k==3,"can only read triangles from OFF");
    is >>T(i,0) >>T(i,1) >>T(i,2);
  }
}

void ors::Mesh::readPlyFile(const char* filename){
  ifstream is;
  MT::open(is,filename);
  uint i,k,nVertices,nFaces;
  is >>"ply" >>"format ascii 1.0";
  is >>"element vertex" >>nVertices;
  is >>"property float32 x" >>"property float32 y" >>"property float32 z";
  is >>"property float32 nx" >>"property float32 ny" >>"property float32 nz";
  is >>"element face" >>nFaces;
  is >>"property list uint8 int32 vertex_indices" >>"end_header";
  V.resize(nVertices,3);
  T.resize(nFaces   ,3);
  real nx,ny,nz;
  for(i=0;i<V.d0;i++){
    is >>V(i,0) >>V(i,1) >>V(i,2) >>nx >>ny >>nz;
  }
  for(i=0;i<T.d0;i++){
    is >>k >>T(i,0) >>T(i,1) >>T(i,2);
    CHECK(k==3,"can only read triangles from ply");
  }
}

void ors::Mesh::readStlFile(const char* filename){
  ifstream is;
  MT::open(is,filename);
  MT::String name;
  is >>"solid" >>name;
  uint i,k=0,k0;
  real x,y,z;
  cout <<"reading STL file '"<<filename<<"' object name '"<<name <<"'..."<<endl;
  V.resize(10000);
  //1st pass
  for(i=0,k=0;;i++){
    k0=k;
    if(k>V.N-10) V.resizeCopy(2*V.N);
    if(!(i%100)) cout <<"\r" <<i <<' ' <<i*7;
    if(MT::peerNextChar(is)!='f') break;
    is >>(const char*)"facet";
    is >>(const char*)"normal" >>x >>y >>z;  MT::skip(is);
    is >>(const char*)"outer" >>(const char*)"loop";      MT::skip(is);
    is >>(const char*)"vertex">>V(k++); is>>V(k++); is>>V(k++);   MT::skip(is);
    is >>(const char*)"vertex">>V(k++); is>>V(k++); is>>V(k++);   MT::skip(is);
    is >>(const char*)"vertex">>V(k++); is>>V(k++); is>>V(k++);   MT::skip(is);
    is >>(const char*)"endloop";             MT::skip(is);
    is >>(const char*)"endfacet";            MT::skip(is);
    if(!is.good()){
      MT_MSG("reading error - skipping facet "<<i<<" (line "<<i*7+2<<")");
      is.clear();
      cout <<1 <<endl;
      MT::skipUntil(is,"endfacet");
      cout <<2 <<endl;
      k=k0;
    }
  }
  is >>"endsolid";
  if(!is.good()) MT_MSG("couldn't read STL end tag (line" <<i*7+2);
  cout <<"... STL file read: #tris=" <<i <<" #lines=" <<i*7+2 <<endl;
  CHECK(!(k%9),"not mod 9..");
  V.resizeCopy(k/3,3);
  T.resize(k/9,3);
  for(i=0;i<T.N;i++){ T.elem(i)=i; }
}

/*void ors::Mesh::getOBJ(char* filename){
  if (!glm){
  glm = glmReadOBJ(filename);
  glmReverseWinding(glm);
  }
  
  ////glmUnitize(glm);
  glmFacetNormals(glm);
  glmVertexNormals(glm, 90.0);
  
  // creates a display list for the OBJ
  ////	g._pmodel_displaylist = glmList(glm, GLM_SMOOTH | GLM_MATERIAL);
  }*/

uint& Tni(uint,uint){ static uint dummy; return dummy; } //normal index

uint& Tti(uint,uint){ static uint dummy; return dummy; } //texture index



// dm 07.06.2006
/*!\ initialises the ascii-obj file "filename"*/
void ors::Mesh::readObjFile(const char* filename){
  FILE* file;
  
  // open the file 
  file = fopen(filename, "r");
  if(!file) HALT("readObjFile() failed: can't open data file "<<filename);
  
  // make a first pass through the file to get a count of the number
  // of vertices, normals, texcoords & triangles 
  uint nV;        
  uint nN;        
  uint nTex;       
  uint nT;       
  int v, n, t;
  char buf[128];
  
  nV = nN = nTex = nT = 0;

  while(fscanf(file, "%s", buf) != EOF) {
    switch(buf[0]) {
    case '#':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;  // comment
    case 'v':
      switch(buf[1]) {
      case '\0': nV++;    CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;  // vertex
      case 'n':  nN++;    CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;  // normal
      case 't':  nTex++;  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;  // texcoord
      default: HALT("firstPass(): Unknown token '"<<buf<<"'");	break;
      }
      break;
    //case 'm':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  sscanf(buf, "%s %s", buf, buf);  break;
      //mtllibname = strdup(buf);  glmReadMTL(model, buf);
    //case 'u':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;
    //case 'g':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  sscanf(buf, "%s", buf);  break;
    case 'f':               // face 
      v = n = t = 0;
      CHECK(fscanf(file, "%s", buf),"fscan failed");
      // can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d 
      if (strstr(buf, "//")) {
	// v//n 
	sscanf(buf, "%d//%d", &v, &n);
	CHECK(fscanf(file, "%d//%d", &v, &n),"fscan failed");
	CHECK(fscanf(file, "%d//%d", &v, &n),"fscan failed");
	nT++;
	//// group->numtriangles++;
	while(fscanf(file, "%d//%d", &v, &n) > 0) nT++;
      } else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3) {
	// v/t/n 
	CHECK(fscanf(file, "%d/%d/%d", &v, &t, &n),"fscan failed");
	CHECK(fscanf(file, "%d/%d/%d", &v, &t, &n),"fscan failed");
	nT++;
	//// group->numtriangles++;
	while(fscanf(file, "%d/%d/%d", &v, &t, &n) > 0) nT++;
      } else if (sscanf(buf, "%d/%d", &v, &t) == 2) {
	// v/t 
	CHECK(fscanf(file, "%d/%d", &v, &t),"fscan failed");
	CHECK(fscanf(file, "%d/%d", &v, &t),"fscan failed");
	nT++;
	////group->numtriangles++;
	while(fscanf(file, "%d/%d", &v, &t) > 0) nT++;
      } else {
	// v 
	CHECK(fscanf(file, "%d", &v),"fscan failed");
	CHECK(fscanf(file, "%d", &v),"fscan failed");
	nT++;
	while(fscanf(file, "%d", &v) > 0) nT++;
      }
      break;
	
    default:  MT_MSG("unsupported .obj file tag '"<<buf[0]<<"'");  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;
    }
  }

  //allocate memory
  V.resize(nV,3);
  Vn.resize(nN,3);
  T.resize(nT,3);
  Tn.resize(nT,3);
  //if(nVN) N.resize(nVN,3);
  //if(nTex) Tex.tesize(nTex,2);
  
 
  // rewind to beginning of file and read in the data this pass 
  rewind(file);

  /* on the second pass through the file, read all the data into the
     allocated arrays */
  nV = nN = nTex = nT = 0;
  ////_material = 0;

  while(fscanf(file, "%s", buf) != EOF) {
    switch(buf[0]) {
    case '#':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;  //comment
    case 'v':               // v, vn, vt 
      switch(buf[1]) {
      case '\0': CHECK(fscanf(file, "%lf %lf %lf", &V(nV,0),&V(nV,1),&V(nV,2)),"fscan failed");  nV++;  break;  //vertex
      case 'n':  CHECK(fscanf(file, "%lf %lf %lf", &Vn(nN,0),&Vn(nN,1),&Vn(nN,2)),"fscan failed");  nN++;  break;  //normal
      case 't':  /*CHECK(fscanf(file, "%f %f", &Tex(nTex,0), &Tex(nTex,1)),"fscan failed");  nTex++;*/  break;  //texcoord
      }
      break;
    //case 'u':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  sscanf(buf, "%s %s", buf, buf);  break;
      //group->material = material = glmFindMaterial(model, buf);*/
    //case 'g':  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  sscanf(buf, "%s", buf);  break;
      //  group = glmFindGroup(model, buf);
      //  group->material = material;
    case 'f':               // face 
      v = n = t = 0;
      CHECK(fscanf(file, "%s", buf),"fscan failed");
      //can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d 
      if (strstr(buf, "//")) {
	// v//n 
	sscanf(buf, "%d//%d", &v, &n);

	T(nT,0) = v < 0 ? v + nV : v;
	Tni(nT,0) = n < 0 ? n + nN : n;
	CHECK(fscanf(file, "%d//%d", &v, &n),"fscan failed");
	T(nT,1) = v < 0 ? v + nV : v;
	Tni(nT,1) = n < 0 ? n + nN : n;
	CHECK(fscanf(file, "%d//%d", &v, &n),"fscan failed");
	T(nT,2) = v < 0 ? v + nV : v;
	Tni(nT,2) = n < 0 ? n + nN : n;
	//// group->triangles[group->nT++] = nT;
	nT++;
	while(fscanf(file, "%d//%d", &v, &n) > 0) {
	  T(nT,0) = T(nT-1,0);
	  Tni(nT,0) = Tni(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  Tni(nT,1) = Tni(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  Tni(nT,2) = n < 0 ? n + nN : n;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      } else if (sscanf(buf, "%d/%d/%d", &v, &t, &n) == 3) {
	// v/t/n 
	T(nT,0) = v < 0 ? v + nV : v;
	Tti(nT,0) = t < 0 ? t + nTex : t;
	Tni(nT,0) = n < 0 ? n + nN : n;
	CHECK(fscanf(file, "%d/%d/%d", &v, &t, &n),"fscan failed");
	T(nT,1) = v < 0 ? v + nV : v;
	Tti(nT,1) = t < 0 ? t + nTex : t;
	Tni(nT,1) = n < 0 ? n + nN : n;
	CHECK(fscanf(file, "%d/%d/%d", &v, &t, &n),"fscan failed");
	T(nT,2) = v < 0 ? v + nV : v;
	Tti(nT,2) = t < 0 ? t + nTex : t;
	Tni(nT,2) = n < 0 ? n + nN : n;
	//// group->triangles[group->numtriangles++] = numtriangles;
	nT++;
	while(fscanf(file, "%d/%d/%d", &v, &t, &n) > 0) {
	  T(nT,0) = T(nT-1,0);
	  Tti(nT,0) = Tti(nT-1,0);
	  Tni(nT,0) = Tni(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  Tti(nT,1) = Tti(nT-1,2);
	  Tni(nT,1) = Tni(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  Tti(nT,2) = t < 0 ? t + nTex : t;
	  Tni(nT,2) = n < 0 ? n + nN : n;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      } else if (sscanf(buf, "%d/%d", &v, &t) == 2) {
	// v/t 
	      
	T(nT,0) = v < 0 ? v + nV : v;
	Tti(nT,0) = t < 0 ? t + nTex : t;
	CHECK(fscanf(file, "%d/%d", &v, &t),"fscan failed");
	T(nT,1) = v < 0 ? v + nV : v;
	Tti(nT,1) = t < 0 ? t + nTex : t;
	CHECK(fscanf(file, "%d/%d", &v, &t),"fscan failed");
	T(nT,2) = v < 0 ? v + nV : v;
	Tti(nT,2) = t < 0 ? t + nTex : t;
	//// group->triangles[group->numtriangles++] = numtriangles;
	nT++;
	while(fscanf(file, "%d/%d", &v, &t) > 0) {
	  T(nT,0) = T(nT-1,0);
	  Tti(nT,0) = Tti(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  Tti(nT,1) = Tti(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  Tti(nT,2) = t < 0 ? t + nTex : t;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      } else {
	// v 
	sscanf(buf, "%d", &v);
	T(nT,0) = v < 0 ? v + nV : v;
	CHECK(fscanf(file, "%d", &v),"fscan failed");
	T(nT,1) = v < 0 ? v + nV : v;
	CHECK(fscanf(file, "%d", &v),"fscan failed");
	T(nT,2) = v < 0 ? v + nV : v;
	//// group->triangles[group->numtriangles++] = nT;
	nT++; // fixed bug, 21. Jun 06 (hh)
	while(fscanf(file, "%d", &v) > 0) {
	  T(nT,0) = T(nT-1,0);
	  T(nT,1) = T(nT-1,2);
	  T(nT,2) = v < 0 ? v + nV : v;
	  //// group->triangles[group->numtriangles++] = numtriangles;
	  nT++;
	}
      }
      break;
	    
    default:  CHECK(fgets(buf, sizeof(buf), file),"fgets failed");  break;
    }
  }

  //CONVENTION!: start counting vertex indices from 0!!
  T -= (uint)1;
  
  // close the file 
  fclose(file);
}



void inertiaSphere(real *I, real& mass, real density, real radius){
  real r2=radius*radius;
  if(density) mass=density*4./3.*MT_PI*r2*radius;
  I[1]=I[2]=I[3]=I[5]=I[6]=I[7]=0.;
  I[0]=.4*mass*r2;
  I[4]=.4*mass*r2;
  I[8]=.4*mass*r2;
}

void inertiaBox(real *I, real& mass, real density, real dx,real dy,real dz){
  if(density) mass=density*dx*dy*dz;
  I[1]=I[2]=I[3]=I[5]=I[6]=I[7]=0.;
  real x2=dx*dx,y2=dy*dy,z2=dz*dz;
  I[0]=mass/12.*(y2+z2);
  I[4]=mass/12.*(x2+z2);
  I[8]=mass/12.*(x2+y2);
}

void inertiaCylinder(real *I, real& mass, real density, real height,real radius){
  real r2=radius*radius, h2=height*height;
  if(density) mass=density*MT_PI*r2*height;
  I[1]=I[2]=I[3]=I[5]=I[6]=I[7]=0.;
  I[0]=mass/12.*(3.*r2+h2);
  I[4]=mass/12.*(3.*r2+h2);
  I[8]=mass/2.*r2;
}

//================================================================================
//
// Spline
//

void ors::Spline::plotBasis(){
#ifdef MT_plot_h
  plotClear();
  arr b_sum(T+1);
  tensorMarginal(b_sum,basis_trans,TUP(1));
  plotFunction(b_sum,-1,1);
  for(uint i=0;i<=K;i++) plotFunction(basis_trans[i],-1,1);
  plot();
#else
  NIY;
#endif
}

void ors::Spline::setBasis(){
  uint i,t,p;
  real time,x,y;
  CHECK(times.N-1==K+1+degree,"wrong number of time knots");
  arr b(K+1,T+1),b_0(K+1,T+1);
  for(p=0;p<=degree;p++){
    if(p>0) b_0=b;
    for(i=0;i<=K;i++) for(t=0;t<=T;t++){
      time = (real)t/(real)T;
      if(!p){
	b(i,t) = 0.;
	if(times(i)<=time && time<times(i+1)) b(i,t)=1.;
	if(t==T && i==K && time==times(i+1)) b(i,t)=1.;
      }else{
	x=MT::DIV(time-times(i), times(i+p)-times(i), true);
	b(i,t) = x * b_0(i,t);
	if(i<K){
	  y=MT::DIV(times(i+p+1)-time, times(i+p+1)-times(i+1), true);
	  b(i,t) += y * b_0(i+1,t);
	}
      }
    }
  }
  basis_trans=b;
  transpose(basis,b);
}

void ors::Spline::setBasisAndTimeGradient(){
  uint i,j,t,p,m=times.N-1;
  real time,x,xx,y,yy;
  CHECK(m==K+1+degree,"wrong number of time knots");
  arr b(K+1,T+1),b_0(K+1,T+1),dbt(m+1,K+1,T+1),dbt_0(m+1,K+1,T+1);
  for(p=0;p<=degree;p++){
    if(p>0){ b_0=b; dbt_0=dbt; }
    for(i=0;i<=K;i++) for(t=0;t<=T;t++){
      time = (real)t/(real)T;
      if(!p){
	b(i,t) = 0.;
	if(times(i)<=time && time<times(i+1)) b(i,t)=1.;
	if(t==T && i==K && time==times(i+1)) b(i,t)=1.;
	for(j=0;j<=m;j++) dbt(j,i,t)=0.;
      }else{
	xx=times(i+p)-times(i);
	x=MT::DIV(time-times(i),xx,true);
	if(i<K){
	  yy=times(i+p+1)-times(i+1);
	  y=MT::DIV(times(i+p+1)-time,yy,true);
	}else{
	  yy=1.;
	  y=0.;
	}
	b(i,t) = x * b_0(i,t);
	if(i<K) b(i,t) += y * b_0(i+1,t);
	for(j=0;j<=m;j++){
	  dbt(j,i,t) = x * dbt_0(j,i,t);
	  if(i<K) dbt(j,i,t) += y * dbt_0(j,i+1,t);
	  if(j==i)            dbt(j,i,t) += MT::DIV((x-1),xx,true) * b_0(i,t);
	  if(j==i+p)          dbt(j,i,t) -= MT::DIV(   x ,xx,true) * b_0(i,t);
	  if(i<K && j==i+1)   dbt(j,i,t) += MT::DIV(   y ,yy,true) * b_0(i+1,t);
	  if(i<K && j==i+p+1) dbt(j,i,t) -= MT::DIV((y-1),yy,true) * b_0(i+1,t);
	}
      }
    }
  }
  basis_trans=b;
  transpose(basis,b);
  basis_timeGradient=dbt;
}

void ors::Spline::setUniformNonperiodicBasis(uint _T,uint _K,uint _degree){
  T=_T; K=_K; degree=_degree;
  uint i,m;
  m=K+1+degree;
  times.resize(m+1);
  for(i=0;i<=m;i++){
    if(i<=degree) times(i)=.0;
    else if(i>=m-degree) times(i)=1.;
    else times(i) = real(i-degree)/real(m-2*degree);
  }
  //setBasis(T,K,degree);
  setBasisAndTimeGradient();
}

void ors::Spline::evalF(arr& f_t,uint t) const{ f_t = basis[t]*points; };
void ors::Spline::evalF(arr& f) const{ f = basis*points; };

void ors::Spline::partial(arr& dCdx,const arr& dCdf) const{
  CHECK(dCdf.d0==T+1 && dCdf.d1==points.d1,"");
  dCdx = basis_trans * dCdf;
}

void ors::Spline::partial(arr& dCdx,arr& dCdt,const arr& dCdf,bool constrain) const{
  CHECK(dCdf.d0==T+1 && dCdf.d1==points.d1,"");
  CHECK(basis_timeGradient.N,"");
  uint n=dCdf.d1, m=K+1+degree, j;
  dCdx = basis_trans * dCdf;
  arr X;
  X.referTo(points);
  X.reshape((K+1)*n);
  arr B;
  B.referTo(basis_timeGradient);
  B.reshape((m+1)*(K+1),T+1);
  arr Z = B * dCdf; Z.reshape(m+1,(K+1)*n);
  dCdt = Z*X;
  if(constrain){
    for(j=0;j<=degree;j++) dCdt(j)=0.;
    for(j=m-degree;j<=m;j++) dCdt(j)=0.;
  }
  dCdt(0)=dCdt(m)=0.;
}


#define LEN .2


//===========================================================================
//
// Body implementations
//


/* dm 15.06.2006--these had to be changed with the switch to ODE 0.6
   void Body::copyFrameToOde(){
   CHECK(X.r.isNormalized(),"quaternion is not normalized!");
   CP3(b->posr.pos,X.p);                // dxBody changed in ode-0.6 ! 14. Jun 06 (hh)
   CP4(b->q,X.r); dQtoR(b->q,b->posr.R);
   CP3(b->lvel,X.v);
   CP3(b->avel,X.w);
   }
   void Body::getFrameFromOde(){
   CP3(X.p.v,b->posr.pos);
   CP4(X.r.q,b->q);
   CP3(X.v.v,b->lvel);
   CP3(X.w.v,b->avel);
   CHECK(X.r.isNormalized(),"quaternion is not normalized!");
   }
*/

ors::Body::~Body(){ reset(); }

void ors::Body::reset(){
   listDelete(ats);
   X.setZero();
   fixed=false;
   shapes.memMove=true;
   com.setZero();
   mass=1.;
   inertia.setId();
   inertia*=.2;
}

void ors::Body::write(std::ostream& os) const{
  os <<"X=" <<X <<" ";
  listWrite(ats,os);
}

#define RERR(x) { HALT("ORS FILE ERROR: " <<x); is.clear(); return; }

void ors::Body::read(std::istream& is){
  reset();

  //old convention: the first frame is X
  if(MT::peerNextChar(is)=='<'){
    is >>X;
    if(!is.good()) RERR("READING ABORT - body '"<<name<<"' read error: could not read Frame tag properly");
  }

  anyListRead(ats,is);
  if(!is.good()) HALT("body '"<<name<<"' read error: in ");
  
  //interpret some of the attributes
  real *dval;
  MT::String *sval;
  sval=anyListGet<MT::String>(ats,"X",1);    if(sval) X.setText(*sval);

  //shape declared in body attributes..
  dval=anyListGet<real>(ats,"type",1);     if(dval){
    if(!shapes.N) shapes.append(new Shape());
    shapes(0)->body=this;
    shapes(0)->ibody=index;
    shapes(0)->type=(int)(*dval);
  }
  dval=anyListGet<real>(ats,"size",4);     if(dval) memmove(shapes(0)->size,dval,4*sizeof(real));
  dval=anyListGet<real>(ats,"color",3);    if(dval) memmove(shapes(0)->color,dval,3*sizeof(real));
  sval=anyListGet<MT::String>(ats,"rel",1);  if(sval) shapes(0)->rel.setText(*sval);
  sval=anyListGet<MT::String>(ats,"mesh",1); if(sval) shapes(0)->mesh.readFile(*sval);
  if(anyListGet<real>(ats,"contact",0))    shapes(0)->cont=true;

  //mass properties
  dval=anyListGet<real>(ats,"mass",1);     if(dval){
    mass=*dval;
#if 1
    inertia.setId();
    inertia *= .2*(*dval);
#else
    switch(shapes(0)->type){
    case BSPHERE:   inertiaSphere  ( inertia.m, mass, 1000., shapes(0)->size[3]);  break;
    case BBOX:      inertiaBox     ( inertia.m, mass, 1000., shapes(0)->size[0], shapes(0)->size[1], shapes(0)->size[2]);  break;
    case BCCYLINDER:
    case BCYLINDER: inertiaCylinder( inertia.m, mass, 1000., shapes(0)->size[2], shapes(0)->size[3]);  break;
    case BNONE:
    default: HALT("can't set inertia tensor for this type");
    }
#endif
  }


  if(anyListGet<real>(ats,"fixed",0))  fixed=true; else fixed=false;

}


//===========================================================================
//
// Shape implementations
//

void ors::Shape::read(std::istream& is){
  reset();
  anyListRead(ats,is);
  if(!is.good()) HALT("shape read error");
  //listWrite(ats,cout); cout <<endl;
  
  real *dval;
  MT::String *sval;
  sval=anyListGet<MT::String>(ats,"rel",1);  if(sval) rel.setText(*sval);
  dval=anyListGet<real>(ats,"color",3);    if(dval) memmove(color,dval,3*sizeof(real));
  dval=anyListGet<real>(ats,"size",4);     if(dval) memmove(size,dval,4*sizeof(real));
  dval=anyListGet<real>(ats,"type",1);     if(dval) type=(int)(*dval);

  sval=anyListGet<MT::String>(ats,"mesh",1); if(sval) mesh.readFile(*sval);
  if(anyListGet<real>(ats,"contact",0))    cont=true;
}

void ors::Shape::write(std::ostream& os) const{
  listWrite(ats,os);
}

void ors::Shape::reset(){
  type=BNONE;
  size[0]=size[1]=size[2]=size[3]=1.;
  color[0]=color[1]=color[2]=.8;
  contactOrientation.setZero();
  listDelete(ats);
  rel.setZero();
  mesh.V.clear();
  cont=false;
}


//===========================================================================
//
// Joint implementations
//

void ors::Joint::write(std::ostream& os) const{
  os <<"A=" <<A <<" Q=" <<Q <<" B=" <<B <<' ';
  listWrite(ats,os);
}

void ors::Joint::read(std::istream& is){
  reset();
  
  //old convention: the first frame is A | Q | B
  if(MT::peerNextChar(is)=='<'){
    is >>A >>Q >>B;
    if(!is.good()) RERR("READING ABORT - joint read error: could not read Frame tag properly");
  }
  
  //read all generic attributes
  anyListRead(ats,is);
  if(!is.good()) HALT("joint ("<<from->name<<' '<<to->name<<") read read error");

  //interpret some of the attributes
  real *dval;
  MT::String *sval;
  sval=anyListGet<MT::String>(ats,"A",1);  if(sval) A.setText(*sval);
  sval=anyListGet<MT::String>(ats,"B",1);  if(sval) B.setText(*sval);
  sval=anyListGet<MT::String>(ats,"Q",1);  if(sval) Q.setText(*sval);
  dval=anyListGet<real>(ats,"type",1);   if(dval) type=(int)(*dval); else type=JHINGE;
}


//===========================================================================
//
// Graph implementations
//

void ors::Graph::init(const char* filename){
  MT::load(*this,filename,true);
  calcNodeFramesFromEdges();
}

void ors::Graph::clear(){
  sd=jd=td=0;
  listDelete(proxies);
  listDelete(joints);
  listDelete(shapes);
  listDelete(bodies);
}

void ors::Graph::clone(const ors::Graph& G){
  sd=G.sd;  jd=G.jd;  td=G.td;
  listCopy(proxies,G.proxies);
  listCopy(joints, G.joints);
  listCopy(shapes, G.shapes);
  listCopy(bodies, G.bodies);
  graphMakeLists(bodies,joints);
  uint i;  Shape *s;  Body *b;
  for_list(i,s,shapes){
    b=bodies(s->ibody);
    s->body=b;
    b->shapes.append(s);
  }
}

/*!\brief transforms (e.g., translates or rotates) the joints coordinate system):
  `adds' the transformation f to A and its inverse to B */
void ors::Graph::transformEdge(ors::Graph::edge e,const ors::Frame &f){
  e->A.addRelativeFrame(f);
  e->B.makeRelativeToInv(f);
}
  
//! [prelim] some kind of gyroscope
void ors::Graph::getGyroscope(ors::Vector& up){
  up.set(0,0,1);
  up=bodies(0)->X.r*up;
}

/*!\brief KINEMATICS: given the (absolute) frames of root nodes and the relative frames
    on the edges, this calculates the absolute frames of all other nodes (propagating forward
    through trees and testing consistency of loops). */
void ors::Graph::calcNodeFramesFromEdges(){
  node n;
  edge e;
  Shape *s;
  uint i,j;
  ors::Frame f;
  for_list(j,n,bodies){
    for_list(i,e,n->inLinks){
      f = e->from->X;
      f.addRelativeFrame(e->A);
      e->Xworld=f;
      f.addRelativeFrame(e->Q);
      f.addRelativeFrame(e->B);
      if(e==n->inLinks(0))
        n->X=f;
      else{
        MT_MSG("loopy geometry - code missing: check if n->X and f are consistent!");
      }
    }
    for_list(i,s,n->shapes){
      s->X = n->X;
      s->X.addRelativeFrame(s->rel);
    }
  }
}

/*!\brief given the absolute frames of all nodes and the two rigid (relative)
    frames A & B of each edge, this calculates the dynamic (relative) joint
    frame X for each edge (which includes joint transformation and errors) */
void ors::Graph::calcEdgesFromNodeFrames(){
  edge e;
  uint i;
  for_list(i,e,joints){
    ors::Frame A(e->from->X),B(e->to->X);
    A.addRelativeFrame(e->A);
    B.subRelativeFrame(e->B);
    e->Q.setDifference(A,B);
  }
}

/*!\brief in all edge frames: remove any displacements,velocities and non-x rotations.
    After this, edges and nodes are not coherent anymore. You might want to call
    calcNodeFramesFromEdges() */
void ors::Graph::clearEdgeErrors(){
  edge e;
  ors::Vector xaxis(1,0,0);
  uint i;
  for_list(i,e,joints){
    e->Q.p.setZero();
    e->Q.v.setZero();
    e->Q.r.alignWith(xaxis);
    e->Q.w.makeColinear(xaxis);
  }
}
  
  /*!\brief invert all velocity variables of all frames */
void ors::Graph::invertTime(){
  node n;
  edge e;
  uint i,j;
  for_list(j,n,bodies){
    n->X.v*=-1.;
    n->X.w*=-1.;
    for_list(i,e,n->inLinks){
      e->Q.v*=-1.;
      e->Q.w*=-1.;
    }
  }
}

void ors::Graph::computeNaturalQmetric(arr& W){
  //compute generic q-metric depending on tree depth
  uint i,j;
  arr BM(bodies.N);
  BM=1.;
  for(i=BM.N;i--;){
    for(j=0;j<bodies(i)->outLinks.N;j++){
      BM(i) += BM(bodies(i)->outLinks(j)->to->index);
    }
  }
  if(!jd) jd = getJointStateDimension();
  arr Wdiag(jd);
  for(i=0;i<jd;i++) Wdiag(i)=::pow(BM(joints(i)->to->index),1.);
  W.setDiag(Wdiag);
}

  /*!\brief revert the topological orientation of a joint (edge),
     e.g., when choosing another body as root of a tree */
void ors::Graph::revertEdge(ors::Graph::edge e){
  cout <<"reverting edge (" <<e->from->name <<' ' <<e->to->name <<")" <<endl;
  //revert
  uint i=e->ifrom;  e->ifrom=e->ito;  e->ito=i;
  graphMakeLists(bodies,joints);
  ors::Frame f;     //!< transformation from parent body to joint (attachment, usually static)
  f=e->A;
  e->A.setInverse(e->B);
  e->B.setInverse(f);
  f=e->Q;
  e->Q.setInverse(f); 
}

  /*!\brief re-orient all joints (edges) such that n becomes
    the root of the configuration */
void ors::Graph::reconfigureRoot(const ors::Graph::node& n){
  MT::Array<node> list,list2;
  node *m;
  edge e;
  list.append(n);
  uintA level(bodies.N);
  level=0;
  int i=0;
  uint j;
  
  while(list.N>0){
    i++;
    list2.clear();
    for(m=list.p;m!=list.pstop;m++){
      level((*m)->index)=i;
      for_list(j,e,(*m)->inLinks){
      //for_in_edges_save(e,es,(*m))
        if(!level(e->from->index)){ revertEdge(e); j--; }
      }
      for_list(j,e,(*m)->outLinks) list2.append(e->to);
    }
    list=list2;
  }
  
  graphTopsort(bodies,joints);
}

  /*!\brief returns the dimensionality of the full dynamic state vector (0
     DOFs for fixed bodies, 13 DOFs for free bodies (including a
     quaternion), 2 DOFs for jointed bodies) */
uint ors::Graph::getFullStateDimension(){
  node n;
  uint i=0,j;
  for_list(j,n,bodies){
    if(!n->inLinks.N && !n->fixed) i+=6;
    else if(n->inLinks.N){
      switch(n->inLinks(0)->type){
      case 0: i+=1; break;
#ifndef Qstate
      case 4: i+=3; break;
#else
      case 4: i+=4; break;
#endif
      }
    }
  }
  sd=i;
  return i;
}

/*!\brief returns the joint (actuator) dimensionality */
uint ors::Graph::getJointStateDimension() const{
  edge e;
  uint n=0,i;
  
  for_list(i,e,joints){
    switch(e->type) //  3. Apr 06 (hh)
      {
      case JHINGE: case JSLIDER: n++;  break;
      case JGLUE:  case JFIXED:        break;
      case JUNIVERSAL:           n+=2; break;
      default: NIY; break;
      }
  }
  ((Graph*)this)->jd=n;
  return n;
}

  /*!\brief returns the full state vector */
void ors::Graph::getFullState(arr& x){
  HALT("outdated");
  node n;
  edge e;
  uint i=0,j;

  if(!sd) sd=getFullStateDimension();
  x.resize(sd);
  
  ors::Vector rot;
  for_list(j,n,bodies){
    if(!n->inLinks.N && !n->fixed){
      memmove(&x(i),&(n->X),13*sizeof(real));
      i+=13;
    }
    e=n->inLinks(0);
    if(e){
      switch(e->type){
      case 0:
        e->Q.r.getRad(x(i),rot);
        if(x(i)>MT_PI) x(i)-=MT_2PI;
        if(rot*VEC_x<0.) x(i)=-x(i); //MT_2PI-x(i);
        i+=1;
        x(i)=e->Q.w.length();
        if(e->Q.w*VEC_x<0.) x(i)=-x(i);
        i+=1;
        break;
      case 4: NIY; break;
      }
    }
  }
  CHECK(i==sd,"");
}

  /*!\brief returns the full state vector */
void ors::Graph::getFullState(arr& x,arr& v){
  node n;
  edge e;
  uint i=0,j;
  
  if(!sd) sd=getFullStateDimension();
  x.resize(sd);
  v.resize(sd);
  
  ors::Vector rot;
  ors::Quaternion q;
  for_list(j,n,bodies) if(!n->fixed){
    if(!n->inLinks.N){
      n->X.r.getVec(rot);
      memmove(&x(i)  ,&(n->X.p),3*x.sizeT);
      memmove(&x(i+3),&(rot)   ,3*x.sizeT);
      memmove(&v(i)  ,&(n->X.v),3*x.sizeT);
      memmove(&v(i+3),&(n->X.w),3*x.sizeT);
      i+=6;
    }else{
      e=n->inLinks(0);
      switch(e->type){
      case 0:
        e->Q.r.getRad(x(i),rot);
        if(x(i)>MT_PI) x(i)-=MT_2PI;
        if(rot*VEC_x<0.) x(i)=-x(i); //MT_2PI-x(i);
        v(i)=e->Q.w.length();
        if(e->Q.w*VEC_x<0.) v(i)=-v(i);
        i+=1;
        break;
#ifndef Qstate
      case 4:
        e->Q.r.getVec(rot);
        memmove(&x(i),&(rot)   ,3*x.sizeT);
        memmove(&v(i),&(e->Q.w),3*x.sizeT);
        i+=3;
        break;
#else
      case 4:
        memmove(&x(i),&(e->Q.r),4*x.sizeT);
        q.q[0]=0.; q.q[1]=.5*e->Q.w(0); q.q[2]=.5*e->Q.w(1); q.q[3]=.5*e->Q.w(2);
        q = e->Q.r*q;
        v(i)=q.q[0]; v(i+1)=q.q[1]; v(i+2)=q.q[2]; v(i+3)=q.q[3];
        i+=4;
	break;
#endif
      }
    }
  }
  CHECK(i==sd,"");
}

  /*!\brief sets the state of the cofiguration as given by the state vector x;
      if clearJointErrors is true, the joints are assigned zero translational
      (and translational velocity) errors */
void ors::Graph::setFullState(const arr& x,bool clearJointErrors){
  HALT("outdated");
  node n;
  edge e;
  uint i=0,j;

  if(!sd) sd=getFullStateDimension();
  CHECK(x.N==sd,"state doesn't have right dimension ("<<x.N<<"!="<<sd<<")");
  
  for_list(j,n,bodies){
    if(!n->inLinks.N && !n->fixed){
      memmove(&(n->X),&x(i),13*sizeof(real));
      if(!n->X.r.isNormalized()){
	MT_MSG("normalizing quaternion of input state");
	n->X.r.normalize();
	memmove((void*)&x(i),&(n->X),13*sizeof(real));
      }
      i+=13;
    }
    e=n->inLinks(0);
    if(e){
      switch(e->type){
      case 0:
        e->Q.r.setRadX(x(i)); 
        i+=1;
        if(e->Q.w.isZero()) e->Q.w=VEC_x;
        if(e->Q.w*VEC_x<0.) e->Q.w.setLength(-x(i)); else e->Q.w.setLength(x(i));
        if(clearJointErrors){
          e->Q.p.setZero();
          e->Q.v.setZero();
	        //truely, also the rotations X.r and X.w should be made orthgonal to the x-axis
        }
        i+=1;
        break;
      case 4:
        NIY;
        i+=1;
        break;
      }
    }
  }
}

  /*!\brief sets the state of the cofiguration as given by the state vector x;
      if clearJointErrors is true, the joints are assigned zero translational
      (and translational velocity) errors */
void ors::Graph::setFullState(const arr& x,const arr& v,bool clearJointErrors){
  node n;
  edge e;
  uint i=0,j;

  if(!sd) sd=getFullStateDimension();
  CHECK(x.N==sd,"state doesn't have right dimension ("<<x.N<<"!="<<sd<<")");

  ors::Vector rot;
  ors::Quaternion q;
  for_list(j,n,bodies) if(!n->fixed){
    if(!n->inLinks.N){
      memmove(&(n->X.p),&x(i)  ,3*x.sizeT);
      memmove(&(rot)   ,&x(i+3),3*x.sizeT);
      memmove(&(n->X.v),&v(i)  ,3*x.sizeT);
      memmove(&(n->X.w),&v(i+3),3*x.sizeT);
      n->X.r.setVec(rot);
      i+=6;
    }
    e=n->inLinks(0);
    if(e){
      switch(e->type){
      case 0:
	e->Q.r.setRadX(x(i)); 
	if(e->Q.w.isZero()) e->Q.w=VEC_x;
	if(e->Q.w*VEC_x<0.) e->Q.w.setLength(-v(i)); else e->Q.w.setLength(v(i));
	if(clearJointErrors){
	  e->Q.p.setZero();
	  e->Q.v.setZero();
	  //truely, also the rotations X.r and X.w should be made orthgonal to the x-axis
	}
	i+=1;
	break;
#ifndef Qstate
      case 4:
	memmove(&(rot)   ,&x(i),3*x.sizeT);
	memmove(&(e->Q.w),&v(i),3*x.sizeT);
	e->Q.r.setVec(rot);
	i+=3;
	break;
#else
      case 4:
	memmove(&(e->Q.r),&x(i),4*x.sizeT);
	e->Q.r.normalize();
	memmove(q.q,&v(i),4*x.sizeT);
	//e->Q.r.invert();
	q = e->Q.r/q;
	//e->Q.r.invert();
	e->Q.w(0)=q.q[1]; e->Q.w(1)=q.q[2]; e->Q.w(2)=q.q[3];
	e->Q.w *=2.;
	i+=4;
	break;
#endif
      }
    }
  }
  CHECK(i==sd,"");
}


void ors::Graph::zeroGaugeEdges(){
  node n;
  edge e;
  ors::Vector w;
  uint j;
  for_list(j,n,bodies) if(!n->fixed){
    e=n->inLinks(0);
    if(e){
      w=e->Q.r / e->Q.w; e->Q.w.setZero();
      e->A.addRelativeFrame(e->Q);
      e->Q.setZero();
      e->Q.w=w;
    }
  }
}

/*!\brief returns the joint state vectors separated in positions and
  velocities */
void ors::Graph::getJointState(arr& x,arr& v){
  edge e;
  uint n=0,i;
  ors::Vector rotv;
  ors::Quaternion rot;

  if(!jd) jd = getJointStateDimension();
  x.resize(jd); v.resize(jd); // 16. Mar 06 (hh)

  for_list(i,e,joints){ // 16. Mar 06 (hh)
    switch(e->type){
    case JHINGE: 
      //angle
      e->Q.r.getRad(x(n),rotv);
      if(x(n)>MT_PI) x(n)-=MT_2PI;
      if(rotv*VEC_x<0.) x(n)=-x(n); //MT_2PI-x(i);
      //velocity
      v(n)=e->Q.w.length();
      if(e->Q.w*VEC_x<0.) v(n)=-v(n);
      n++;
      break;   
    case JUNIVERSAL:
      //angle
      rot = e->Q.r;
      if(fabs(rot.q[0])>1e-15) {
        x(n) = 2.0 * atan(rot.q[1]/rot.q[0]); 
        x(n+1) = 2.0 * atan(rot.q[2]/rot.q[0]); 
      }else{
        x(n) = MT_PI;
        x(n+1) = MT_PI;
      }

      // velocity: need to fix
	
      n+=2; 
      break;
    case JSLIDER: // 26. May 06 (hh)
      x(n)=(e->Q.p)(0); 
      v(n)=(e->Q.v)(0);
      n++;
      break;
    case JGLUE:
    case JFIXED:
      break;
      default: NIY; break;
    }
  }
}

  /*!\brief returns the joint positions only */
void ors::Graph::getJointState(arr& x){ arr v; getJointState(x,v); }

  /*!\brief sets the joint state vectors separated in positions and
    velocities */
void ors::Graph::setJointState(const arr& x,const arr& v,bool clearJointErrors){
  edge e;
  uint n=0,i;
  real tempAngleDeg;	// temporarily stores the joint angle to check with boundaries
  real tempDeflection;// temporarily stores the joint deflection in a slider joint
  ors::Quaternion rot1,rot2;
  
  if(!jd) jd = getJointStateDimension();
  CHECK(x.N==jd && v.N==jd,"wrong joint state dimensionalities"); // 16. Mar 06 (hh)
  
  for_list(i,e,joints){
    switch(e->type){
    case JHINGE: 
      //angle
      e->Q.r.setRadX(x(n)); 
      
      // check boundaries
      if(e->p[0] < e->p[1]){
	tempAngleDeg = x(n); //dm *180.0/MT_PI; 
	if(tempAngleDeg <= e->p[0]){ // joint angle is smaller than lower bound (e->p[0])--> clip it
	  e->Q.r.setDeg(e->p[0],VEC_x); 
	  //	cout<<"lower clipping "<< tempAngleDeg <<endl;
	}else if(tempAngleDeg >= e->p[1]){ // joint angle is larger than upper bound (e->p[1])--> clip it
	  e->Q.r.setDeg(e->p[1],VEC_x);
	  //	cout<<"upper clipping "<< tempAngleDeg <<endl;
	}
      }

      //velocity
      e->Q.w.set(v(n),0.,0.);
      //if(e->Q.w.isZero()) e->Q.w=VEC_x;
      //if(e->Q.w*VEC_x<0.) e->Q.w.setLength(-v(n)); else e->Q.w.setLength(v(n));

      if(clearJointErrors){
	e->Q.p.setZero();
	e->Q.v.setZero();
	//truely, also the rotations X.r and X.w should be made orthogonal to the x-axis
      }
      n++;
      break;
      
    case JUNIVERSAL:
      //angle
      rot1.setRadX(x(n)); 
      rot2.setRadY(x(n+1)); 
      
      //check boundaries
      if(e->p[0] < e->p[1]){
	// TODO: both angles are restricted to the same boundaries. Could be enhanced 
	//		 in order to be able to restrict the two angles differently
	tempAngleDeg = x(n); //dm *180.0/MT_PI; 
	if(tempAngleDeg <= e->p[0]){ // joint angle is smaller than lower bound (e->p[0])--> clip it
	  rot1.setRadX(e->p[0]); 
	  rot2.setRadY(e->p[0]);
	  //	cout<<"lower clipping "<< tempAngleDeg <<endl;
	}else if(tempAngleDeg >= e->p[1]){ // joint angle is larger than upper bound (e->p[1])--> clip it
	  rot1.setRadX(e->p[1]); 
	  rot2.setRadY(e->p[1]);
	  //	cout<<"upper clipping "<< tempAngleDeg <<endl;
	}
      }
      
      e->Q.r = rot1*rot2;
      //velocity
      // need to fix
      
      if(clearJointErrors){
        e->Q.p.setZero();
        e->Q.v.setZero();	  
      }
      n+=2; 
      break;
    case JSLIDER: // 26. May 06 (hh)	
      e->Q.p = x(n)*VEC_x;
      
      // check boundaries
      if(e->p[0] < e->p[1]){
	tempDeflection = x(n); 
	if(tempDeflection <= e->p[0]){ // joint angle is smaller than lower bound (e->p[0])--> clip it
	  e->Q.p = e->p[0]*VEC_x;
	  cout<<"lower clipping "<< tempDeflection <<endl;
	}else if(tempDeflection >= e->p[1]){ // joint angle is larger than upper bound (e->p[1])--> clip it
	  e->Q.p = e->p[1]*VEC_x;
	  cout<<"upper clipping "<< tempDeflection <<endl;
	}
      }
      
      //velocity
      e->Q.v = v(n)*VEC_x;
      e->Q.r.setZero();
      e->Q.w.setZero();
      n++;
      break;	
    case JGLUE:
    case JFIXED:
      e->Q.setZero();
      break;
    default:  NIY;  break;
    }
  }
}

  /*!\brief sets the joint angles with velocity zero - e.g. for kinematic
    simulation only */
void ors::Graph::setJointState(const arr& x,bool clearJointErrors){

  //uint n=0;
  arr v;
  
  if(!jd) jd = getJointStateDimension();
  /*
  if(x.N==2*jd){ //assume it contains x and v
    x.reshape(2,jd);
    setJointState(x[0](),x[1](),clearJointErrors); // need to fix, works only for one joint (do we need this part at all?)
    x.reshape(x.N);
    return;
  }*/
  v.resize(jd); v = 0.0;
  setJointState(x,v,clearJointErrors); // 26. May 06 (hh)
}

//===========================================================================
//===========================================================================
//===========================================================================
//
// core: kinematics and dynamics
//
// essential papers: 
// David Baraff: "Linear-Time Dynamics using Lagrange Multipliers"
// Roy Featherstone, David Orin: "Robot Dynamics: Equations and Algorithms"

  /*!\brief return the position \f$x = \phi_i(q)\f$ of the i-th body (3 vector) */
void ors::Graph::kinematics(arr& y,uint a,ors::Frame *rel){
  ors::Frame f=bodies(a)->X;
  if(rel) f.addRelativeFrame(*rel);
  y.setCarray(f.p.v,3);
}

  /*!\brief return the jacobian \f$J = \frac{\partial\phi_i(q)}{\partial q}\f$ of the position
    of the i-th body (3 x n tensor)*/
void ors::Graph::jacobian(arr& J,uint a,ors::Frame *rel){
  uint i;
  ors::Frame Xa,Xi;
  edge ei;
  ors::Vector tmp,ti;

  if(!jd) jd = getJointStateDimension();

  //initialize Jacobian
  J.resize(3,jd);
  J.setZero();

  //get reference frame
  Xa = bodies(a)->X;
  if(rel) Xa.addRelativeFrame(*rel);

  if(!bodies(a)->inLinks.N) return;
  ei=bodies(a)->inLinks(0);
  while(ei){
    i=ei->index;
    if(ei->index>=jd){
      CHECK(ei->type==JGLUE || ei->type==JFIXED,"");
      if(!ei->from->inLinks.N) break;
      ei=ei->from->inLinks(0);
      continue;
    }
    CHECK(ei->type!=JGLUE && ei->type!=JFIXED,"resort joints so that fixed and glued are last");

#if 0
    Xi = ei->from->X;
    Xi.addRelativeFrame(ei->A);
#else
    Xi = ei->Xworld;
#endif
    Xi.r.getX(ti);

    tmp = ti ^ (Xa.p-Xi.p);

    J(0,i) = tmp.v[0];
    J(1,i) = tmp.v[1];
    J(2,i) = tmp.v[2];

    if(!ei->from->inLinks.N) break;
    ei=ei->from->inLinks(0);
  }
}

  /*!\brief return the Hessian \f$H = \frac{\partial^2\phi_i(q)}{\partial q\partial q}\f$ of the position
    of the i-th body (3 x n x n tensor) */
void ors::Graph::hessian(arr& H,uint a,ors::Frame *rel){
  uint i,j;
  ors::Frame Xa,Xi,Xj;
  edge ei,ej;
  ors::Vector r,ti,tj;

  if(!jd) jd = getJointStateDimension();

  //initialize Jacobian
  H.resize(3,jd,jd);
  H.setZero();

  //get reference frame
  Xa = bodies(a)->X;
  if(rel) Xa.addRelativeFrame(*rel);

  if(!bodies(a)->inLinks.N) return;
  ei=bodies(a)->inLinks(0);
  while(ei){
    i=ei->index;

    Xi = ei->from->X;
    Xi.addRelativeFrame(ei->A);
    Xi.r.getX(ti);

    ej=ei;
    while(ej){
      j=ej->index;

      Xj = ej->from->X;
      Xj.addRelativeFrame(ej->A);
      Xj.r.getX(tj);

      r = tj ^ (ti ^ (Xa.p-Xi.p));

      H(0,i,j) = H(0,j,i) = r.v[0];
      H(1,i,j) = H(1,j,i) = r.v[1];
      H(2,i,j) = H(2,j,i) = r.v[2];

      if(!ej->from->inLinks.N) break;
      ej=ej->from->inLinks(0);
    }
    if(!ei->from->inLinks.N) break;
    ei=ei->from->inLinks(0);
  }
}

  /*!\brief return the configuration's inertia tensor $M$ (n x n tensor)*/
void ors::Graph::inertia(arr& M){
  uint a,i,j;
  ors::Frame Xa,Xi,Xj;
  edge ei,ej;
  ors::Vector vi,vj,ti,tj;
  real tmp;

  if(!jd) jd = getJointStateDimension();

  //initialize Jacobian
  M.resize(jd,jd);
  M.setZero();

  for(a=0;a<bodies.N;a++){
    //get reference frame
    Xa = bodies(a)->X;
    
    ei=bodies(a)->inLinks(0);
    while(ei){
      i=ei->index;
      
      Xi = ei->from->X;
      Xi.addRelativeFrame(ei->A);
      Xi.r.getX(ti);
      
      vi = ti ^ (Xa.p-Xi.p);

      ej=ei;
      while(ej){
        j=ej->index;
	
        Xj = ej->from->X;
        Xj.addRelativeFrame(ej->A);
        Xj.r.getX(tj);
	
        vj = tj ^ (Xa.p-Xj.p);
	
        tmp = bodies(a)->mass * (vi*vj);
        //tmp += scalarProduct(bodies(a)->a.inertia,ti,tj);

        M(i,j) += tmp;
	
        if(!ej->from->inLinks.N) break;
        ej=ej->from->inLinks(0);
      }
      if(!ei->from->inLinks.N) break;
      ei=ei->from->inLinks(0);
    }
  }
  //symmetric: fill in other half
  for(i=0;i<jd;i++) for(j=0;j<i;j++) M(j,i) = M(i,j);
}

void ors::Graph::equationOfMotion(arr& M,arr& F,const arr& qd){
  static ors::LinkTree tree;
  if(!tree.N) GraphToTree(tree,*this);
  else updateGraphToTree(tree,*this);
  ors::equationOfMotion(M,F,tree,qd);
}

  /*!\brief return the joint accelerations \f$\ddot q\f$ given the
    joint torques \f$\tau\f$ (computed via Featherstone's Articulated Body Algorithm in O(n)) */
void ors::Graph::dynamics(arr& qdd,const arr& qd,const arr& tau){
#if 0
  Featherstone::Robot r;
  r.C = this;
  arr dummy;
  Featherstone::fwdDynamics_old(qdd,r,qd,tau,dummy);
#else
  static ors::LinkTree tree;
  if(!tree.N) GraphToTree(tree,*this);
  else updateGraphToTree(tree,*this);
  //cout <<tree <<endl; 
  //ors::fwdDynamics_aba_1D(qdd,tree,qd,tau);
  //ors::fwdDynamics_aba_nD(qdd,tree,qd,tau);
  ors::fwdDynamics_MF(qdd,tree,qd,tau);
#endif
}

  /*!\brief return the necessary joint torques \f$\tau\f$ to achieve joint accelerations
    \f$\ddot q\f$ (computed via the Recursive Newton-Euler Algorithm in O(n)) */
void ors::Graph::inverseDynamics(arr& tau,const arr& qd,const arr& qdd){
#if 0
  Featherstone::Robot r;
  r.C = this;
  arr dummy;
  Featherstone::invdyn_old(tau,r,qd,qdd,dummy);
#else
  static ors::LinkTree tree;
  if(!tree.N) GraphToTree(tree,*this);
  else updateGraphToTree(tree,*this);
  ors::invDynamics(tau,tree,qd,qdd);
#endif
}

/*void ors::Graph::impulsePropagation(arr& qd1,const arr& qd0){
  static MT::Array<Featherstone::Link> tree;
  if(!tree.N) GraphToTree(tree,*this);
  else updateGraphToTree(tree,*this);
  mimickImpulsePropagation(tree);
  Featherstone::RF_abd(qdd,tree,qd,tau);
}*/

  //! kinematis of the i-th body's z-orientation vector
void ors::Graph::kinematicsZ(arr& y,uint a,ors::Frame *rel){
  ors::Frame f=bodies(a)->X;
  if(rel) f.addRelativeFrame(*rel);
  ors::Vector v;
  f.r.getZ(v);
  y.setCarray(v.v,3);
}

/* takes the joint state x and returns the jacobian dz of
   the position of the ith body (w.r.t. all joints) -> 2D array */
  //! Jacobian of the i-th body's z-orientation vector
void ors::Graph::jacobianZ(arr& J,uint a,ors::Frame *rel){
  uint i;
  ors::Frame Xa,Xi;
  edge ei;
  ors::Vector r,ta,ti;

  if(!jd) jd = getJointStateDimension();

  //initialize Jacobian
  J.resize(3,jd);
  J.setZero();

  //get reference frame
  Xa = bodies(a)->X;
  if(rel) Xa.addRelativeFrame(*rel);
  Xa.r.getZ(ta);

  if(!bodies(a)->inLinks.N) return;
  ei=bodies(a)->inLinks(0);
  while(ei){
    i=ei->index;
    if(ei->index>=jd){
      CHECK(ei->type==JGLUE || ei->type==JFIXED,"");
      if(!ei->from->inLinks.N) break;
      ei=ei->from->inLinks(0);
      continue;
    }
    CHECK(ei->type!=JGLUE && ei->type!=JFIXED,"resort joints so that fixed and glued are last");

    Xi = ei->Xworld;
    Xi.r.getX(ti);

    r = ti ^ ta;

    J(0,i) = r.v[0];
    J(1,i) = r.v[1];
    J(2,i) = r.v[2];

    if(!ei->from->inLinks.N) break;
    ei=ei->from->inLinks(0);
  }
}

/* takes the joint state x and returns the jacobian dz of
   the position of the ith body (w.r.t. all joints) -> 2D array */
void ors::Graph::jacobianR(arr& J,uint a){
  uint i;
  ors::Frame Xi;
  edge ei;
  ors::Vector r,ti;

  if(!jd) jd = getJointStateDimension();

  //initialize Jacobian
  J.resize(3,jd);
  J.setZero();

  //get reference frame -- in this case we always take
  //the Z and X-axis of the world system as references
  // -> don't need to compute explicit reference for object a
  //  object a is relevant in the sense that only the tree-down
  //  joints contribute to this rotation

  if(!bodies(a)->inLinks.N) return;
  ei=bodies(a)->inLinks(0);
  while(ei){
    i=ei->index;
    if(ei->index>=jd){
      CHECK(ei->type==JGLUE || ei->type==JFIXED,"");
      if(!ei->from->inLinks.N) break;
      ei=ei->from->inLinks(0);
      continue;
    }
    CHECK(ei->type!=JGLUE && ei->type!=JFIXED,"resort joints so that fixed and glued are last");

    Xi = ei->Xworld;
    Xi.r.getX(ti);

    //we measure how much the rotation around joint i contributes
    // to a displacement of the world Z- and X-axis
    r = ti ^ VEC_z;
    J(0,i) = r.v[0];
    J(1,i) = r.v[1];

    r = ti ^ VEC_x;
    J(2,i) = r.v[1];

    if(!ei->from->inLinks.N) break;
    ei=ei->from->inLinks(0);
  }
}


/*void ors::Graph::kinematicsOri2(arr& z,uint a,const ors::Vector &rel){
  CHECK(rel.isNormalized(),"I take only normalized orientations");
  z.resize(2);
  ors::Frame *Xa=&(bodies(a)->X);
  ors::Vector v,ey;
  ors::Quaternion up;
  up.setDiff(rel,ors::Vector(0,0,1));
  Xa->r.getY(ey);
  z(0) = 1. - ey * rel; //rel.angle(ey);
  v=up*ey;
  z(1) = .01* MT::phi(v(0),v(1));
}
void ors::Graph::jacobianOri2(arr& J,uint a,const ors::Vector &rel){
  CHECK(rel.isNormalized(),"I take only normalized orientations");
  uint i;
  ors::Frame Xa,Xi;
  edge ei;
  ors::Vector r,ti;

  if(!jd) jd = getJointStateDimension();

  //initialize Jacobian
  J.resize(2,jd);
  J.setZero();

  //get reference frame
  Xa = bodies(a)->X;

  //get features of this frame
  ors::Vector ey,v,dv;
  ors::Quaternion up;
  up.setDiff(rel,ors::Vector(0,0,1));
  Xa.r.getY(ey);

  if(!bodies(a)->inLinks.N) return;
  ei=bodies(a)->inLinks(0);
  while(ei){
    i=ei->index;
    if(ei->index>=jd){
      CHECK(ei->fixed,"");
      ei=ei->from->inLinks(0);
      continue;
    }

    Xi = ei->from->X;
    Xi.addRelativeFrame(ei->A);
    Xi.r.getX(ti); //ti is the rotation axis of the i-th joint

    r = ti ^ (Xa.p-Xi.p);

    J(0,i) = - (ti ^ ey) * rel;
    v=up*ey; dv=up*(ti^ey);
    J(1,i) = .01* MT::dphi(v(0),v(1),dv(0),dv(1));//(up * (ti ^ ey))(0);

    if(!ei->from->inLinks.N) break;
    ei=ei->from->inLinks(0);
  }
}*/

  //! Jacobian of _all_ body positions ( (3k) x n tensor, where k is the number of bodies)
/* takes the joint state x and returns the jacobian dz of
   the position of the ith body (w.r.t. all joints) -> 2D array */
void ors::Graph::jacobianAll(arr& J){
  arr x,v;
  getJointState(x,v);
  uint i,j,n=x.N;
  J.resize(bodies.N,3,n);
  for(j=0;j<n;j++){
    v.setZero();
    v(j)=1.;
    setJointState(x,v);
    calcNodeFramesFromEdges(); // transform each body's frame relative to the world frame
    for(i=0;i<bodies.N;i++){
      J(i,0,j) = bodies(i)->X.v(0);
      J(i,1,j) = bodies(i)->X.v(1);
      J(i,2,j) = bodies(i)->X.v(2);
    }
  }
  J.reshape(bodies.N*3,n);
}

  //! [prelim] some heuristic measure for the joint errors
real ors::Graph::getJointErrors(){
  edge e;
  real err=0.0;
  uint i;

  for_list(i,e,joints) err+=e->Q.p.lengthSqr();

  return ::sqrt(err);
}

  /*!\brief checks if all names of the bodies are disjoint */
bool ors::Graph::checkUniqueNames() const{
  node n,m;
  uint i,j;
  for_list(i,n,bodies) for_list(j,m,bodies){
    if(n==m) break;
    if(n->name==m->name) return false;
  }
  return true;
}

  //! find body with specific name
ors::Graph::node ors::Graph::getBodyByName(const char* name){
  node n;
  uint j;
  for_list(j,n,bodies){
    if(strcmp(n->name,name)==0) return n; // 16. Mar 06 (hh)
  }
  MT_MSG("cannot find Body named '"<<name<<"' in Graph");
  return 0;
}

  //! find body with specific name
ors::Shape* ors::Graph::getShapeByName(const char* name){
  Shape *s;
  uint j;
  for_list(j,s,shapes){
    if(strcmp(s->name,name)==0) return s; // 16. Mar 06 (hh)
  }
  MT_MSG("cannot find Shape named '"<<name<<"' in Graph");
  return NULL;
}

  //! find joint connecting two bodies with specific names
ors::Graph::edge ors::Graph::getName(const char* from,const char* to){
  node f=NULL,t=NULL;
  uint j;
  for_list(j,f,bodies) if(strcmp(f->name,from)==0) break;
  for_list(j,t,bodies) if(strcmp(t->name,to  )==0) break;
  if(!f || !t) return 0;
  return graphGetEdge<Body,Joint>(f,t);
}

  /*!\brief creates uniques names by prefixing the node-index-number to each name */
void ors::Graph::prefixNames(){
  node n;
  uint j;
  for_list(j,n,bodies) n->name=n->name.prepend(n->index);
}

  /*!\brief prototype for \c operator<< */
void ors::Graph::write(std::ostream& os) const{
  node n;
  edge e;
  Shape *s;
  uint i,j;
  for_list(j,n,bodies){
    os <<"body " <<n->name <<" { ";
    n->write(os);  os <<" }\n";
  }
  os <<std::endl;
  for_list(i,s,shapes){
    os <<"shape (" <<s->body->name <<"){ ";
    s->write(os);  os <<" }\n";
  }
  os <<std::endl;
  for_list(i,e,joints){
    os <<"joint (" <<e->from->name <<' ' <<e->to->name <<"){ ";
    e->write(os);  os <<" }\n";
  }
  os <<"</slGraph>" <<std::endl;
}

#define DEBUG(x) //x

  /*!\brief prototype for \c operator>> */
void ors::Graph::read(std::istream& is){
  node n=NULL,f=NULL,t=NULL; edge e;
  Shape *s;
  String tag,name,node1,node2;
  char c;
  uint j;
  c=MT::peerNextChar(is);
  if(c=='<'){
    MT_MSG("this file format is deprecated! Please convert using the sed filter!");
    is >>"<slGraph>";
    clear();
    for(;;){
      MT::skip(is);
      is.get(c);
      if(is.fail()) RERR("failed to read tag (no proper end of file?)");
      if(c=='<') break;
      switch(c){
        case '%': case '#': MT::skipLine(is); break;
        case 'B':
          n=new Body(bodies);
          s=new Shape(shapes,n); //always create a shape for a body...
          MT::skip(is);
          n->name.read(is," \n\r\t"," \n\r\t");
          n->read(is);
          break;
        case 'J':
          MT::skip(is);
          name.read(is," \n\r\t"," -\n\r\t",false);
          for_list(j,n,bodies) if(n->name==name){ f=n; break; }
          if(!n) RERR("reading edge: don't know from-name "<<name);
          is >>"-";
          MT::skip(is);
          name.read(is,""," \n\r\t");
          for_list(j,n,bodies) if(n->name==name){ t=n; break; }
          if(!n) RERR("reading edge: don't know to-name "<<name);
          e=new Joint(joints,f,t);
          e->read(is);
          break;
        default:
          RERR("READING ABORT -- don't know tag "<<c<<" when reading slGraph");
      }
    }
    is >>"/slGraph>";
    CHECK(is.good(),"could not read slGraph stop tag");
  }else{
    clear();
    for(;;){
      tag.read(is," \t\n\r"," \t\n\r({",false);
      if(!tag.N()) break; //end of file
      DEBUG(cout <<"tag=" <<tag <<endl);
      if(tag=="body"){ //node
        name.read(is," \t\n\r"," \t\n\r({",false);
        DEBUG(cout <<"name=" <<name <<endl);
        n=new Body(bodies);
        n->name = name;
        MT::parse(is,"{");
        if(is.fail()) RERR("can't read opening brace for body ("<<n->name <<")");
        try{ n->read(is); }catch(...) RERR("error in parsing body ("<<n->name <<")");
        MT::parse(is,"}");
        if(is.fail()) RERR("can't read closing brace for body ("<<n->name <<")");
        if(n->shapes.N==1){ //parsing has implicitly added a shape...
          s=n->shapes(0);
          s->index=shapes.N;
          shapes.append(s);
        }
        continue;
      }
      if(tag=="joint"){ //edge
        name.read(is," \t\n\r"," \t\n\r({",false); //potential name - not used
        DEBUG(cout <<"name=" <<name <<endl);
        t=f=NULL;
        MT::parse(is,"(");
        node1.read(is," "," ,)",true);
        DEBUG(cout <<"node1=" <<node1 <<endl);
        for_list(j,n,bodies) if(n->name==node1){ f=n; break; }
        if(!f) RERR("reading edge: don't know from-name "<<node1);
        node2.read(is," "," ,)",true);
        DEBUG(cout <<"node2=" <<node2<<endl);
        for_list(j,n,bodies) if(n->name==node2){ t=n; break; }
        if(!t) RERR("reading edge: don't know to-name "<<node2);
        e=new Joint(joints,f,t);
        MT::parse(is,"{");
        if(is.fail()) RERR("can't read opening brace for edge ("<<e->from->name<<' '<<e->to->name<<")");
        try{ e->read(is); }catch(...) RERR("error in parsing edge ("<<e->from->name<<' '<<e->to->name<<")");
        MT::parse(is,"}");
        if(is.fail()) RERR("can't read closing brace for edge ("<<e->from->name<<' '<<e->to->name<<")");
        continue;
      }
      if(tag=="shape"){ //edge
        name.read(is," \t\n\r"," \t\n\r({",false); //potential name - not used
        DEBUG(cout <<"name=" <<name <<endl);
        f=NULL;
        MT::parse(is,"(");
        node1.read(is," "," )",true);
        DEBUG(cout <<"node1=" <<node1 <<endl);
        for_list(j,n,bodies) if(n->name==node1){ f=n; break; }
        if(!f) RERR("reading shape: don't know from-name "<<node1);
        s=new Shape(shapes,f);
        s->name=name;
        MT::parse(is,"{");
        if(is.fail()) RERR("can't read opening brace for shape ("<<f->name<<'-'<<name<<")");
        try{ s->read(is); }catch(...) RERR("error in parsing shape ("<<f->name<<'-'<<name<<")");
        MT::parse(is,"}");
        if(is.fail()) RERR("can't read closing brace for shape ("<<f->name<<'-'<<name<<")");
        continue;
      }
      RERR("can't parse element of type '"<<tag<<"'");
    }
    is.clear();
  }
  graphMakeLists(bodies,joints);
}


  //! dump the list of current proximities on the screen
void ors::Graph::reportProxies(std::ostream *os){
  uint i;
  int a,b;
  (*os) <<"Proximity report: #" <<proxies.N <<endl;
  for(i=0;i<proxies.N;i++){
    a=proxies(i)->a;
    b=proxies(i)->b;
    (*os)
      <<i <<" ("
      <<a <<':' <<(a!=-1?shapes(a)->body->name.p:"earth") <<")-("
      <<b <<':' <<(b!=-1?shapes(b)->body->name.p:"earth")
      <<") [" <<proxies(i)->age
      <<"] d=" <<proxies(i)->d
      //<<" posA=" <<proxies(i)->posA
      //<<" posB=" <<proxies(i)->posB
      //<<" norm=" <<proxies(i)->posB-proxies(i)->posA
      <<endl;
  }
}

bool ProxySortComp(const ors::Proxy *a,const ors::Proxy *b){
  return (a->a < b->a) || (a->a==b->a && a->b<b->b) || (a->a==b->a && a->b==b->b && a->d < b->d);
}

void ors::Graph::sortProxies(bool deleteMultiple,bool deleteOld){
  uint i;
  if(deleteOld){
    for(i=0;i<proxies.N;i++) if(proxies(i)->age){
      delete proxies(i);
      proxies.remove(i);
      i--;
    }
  }

  std::sort(proxies.p,proxies.pstop,ProxySortComp);

  for(i=0;i<proxies.N;i++) if(proxies(i)->age){
    if(
      (i+1==proxies.N) || //this is the last one
      (i && proxies(i)->a==proxies(i-1)->a && proxies(i)->b==proxies(i-1)->b) || //the previous is older
      (proxies(i)->a!=proxies(i+1)->a || proxies(i)->b!=proxies(i+1)->b) || //the next one is between different objects
      (proxies(i+1)->d>=0.) //the next one is non-colliding
      ){
      delete proxies(i);
      proxies.remove(i);
      i--;
    }
  }

  if(deleteMultiple){
    for(i=0;i<proxies.N;i++) if(!proxies(i)->age){
      if(i && !proxies(i-1)->age && proxies(i)->a==proxies(i-1)->a && proxies(i)->b==proxies(i-1)->b){
        delete proxies(i);
	proxies.remove(i);
	i--;
      }
    }
  }
}

  /*!\brief dump a list body pairs for which the upper conditions hold */
void ors::Graph::reportGlue(std::ostream *os){
  uint i,A,B;
  node a,b;
  bool ag,bg;
  (*os) <<"Glue report: " <<endl;
  for(i=0;i<proxies.N;i++){
    A=proxies(i)->a; a=(A==(uint)-1?NULL:bodies(A));
    B=proxies(i)->b; b=(B==(uint)-1?NULL:bodies(B));
    if(!a || !b) continue;
    ag=anyListGet<real>(a->ats,"glue",0);
    bg=anyListGet<real>(b->ats,"glue",0);

    if(ag || bg){
      (*os)
          <<i <<' '
          <<a->index <<',' <<a->name <<'-'
          <<b->index <<',' <<b->name
          <<" d=" <<proxies(i)->d
	        //<<" posA=" <<proxies(i)->posA
	        //<<" posB=" <<proxies(i)->posB
          <<" norm=" <<proxies(i)->posB-proxies(i)->posA
          <<endl;
    }
  }
}

void ors::Graph::glueBodies(node f,node t){
  Joint *e;
  e=newEdge(f->index,t->index,joints);
  graphMakeLists(bodies,joints);
  e->A.setDifference(f->X,t->X);
  e->A.v.setZero();
  e->A.w.setZero();
  e->type=JGLUE;
  e->Q.setZero();
  e->B.setZero();
}

  /*!\brief if two bodies touch, the are not yet connected, and one of them has
    the `glue' attribute, add a new edge of FIXED type between them */
void ors::Graph::glueTouchingBodies(){
  uint i,A,B;
  node a,b;//,c;
  bool ag,bg;
  for(i=0;i<proxies.N;i++){
    A=proxies(i)->a; a=(A==(uint)-1?NULL:bodies(A));
    B=proxies(i)->b; b=(B==(uint)-1?NULL:bodies(B));
    if(!a || !b) continue;
    ag=anyListGet<real>(a->ats,"glue",0);
    bg=anyListGet<real>(b->ats,"glue",0);

    if(ag || bg){
      //if(a->index > b->index){ c=a; a=b; b=c; } //order them topolgically
      if(graphGetEdge<Body,Joint>(a,b)) continue; //they are already connected
      glueBodies(a,b);
      //a->cont=b->cont=false;
    }
  }
}

  //! clear all forces currently stored at bodies
void ors::Graph::clearForces(){
  node n;
  uint j;
  for_list(j,n,bodies){
    n->force.setZero();
    n->torque.setZero();
  }
}

  //! apply a force on body n at position pos (in world coordinates)
void ors::Graph::addForce(ors::Vector force,node n,ors::Vector pos){
  n->force += force;
  NIY;
  //n->torque += (pos - n->X.p) ^ force;
}

void ors::Graph::frictionToForces(real coeff){
  HALT("never do this: add it directly in the equations...");
  edge e;
  ors::Vector a;
  ors::Frame X;
  real v;
  uint i;
  for_list(i,e,joints){
    X = e->from->X;
    X.addRelativeFrame(e->A);
    X.r.getX(a);//rotation axis

    v=e->Q.w.length();
    if(e->Q.w*VEC_x<0.) v=-v;

    e->from->torque -= (coeff*v)*a;
    e->to->torque   += (coeff*v)*a;
  }
}

void ors::Graph::gravityToForces(){
  node n;
  uint j;
  ors::Vector g(0,0,-9.81);
  for_list(j,n,bodies) n->force += n->mass * g;
}

  //! compute forces from the current contacts
void ors::Graph::contactsToForces(real hook,real damp){
  ors::Vector trans,transvel,force;
  real d;
  uint i;
  int a,b;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<0.){
    a=proxies(i)->a; b=proxies(i)->b;
    
    //if(!i || proxies(i-1).a!=a || proxies(i-1).b!=b) continue; //no old reference sticking-frame
    //trans = proxies(i)->rel.p - proxies(i-1).rel.p; //translation relative to sticking-frame
    trans    = proxies(i)->posB-proxies(i)->posA;
    transvel = proxies(i)->velB-proxies(i)->velA;
    d=trans.length();

    force.setZero();
    force += (hook) * trans; //*(1.+ hook*hook*d*d)
    force += damp * transvel;
    SL_DEBUG(1,cout <<"applying force: [" <<a <<':' <<b <<"] " <<force <<endl);

    if(a!=-1) addForce( force,shapes(a)->body,proxies(i)->posA);
    if(b!=-1) addForce(-force,shapes(b)->body,proxies(i)->posB);
  }
}

//! measure (=scalar kinematics) for the contact cost summed over all bodies
void ors::Graph::getContactMeasure(arr &x,real margin,bool linear){
  x.resize(1);
  x=0.;
  uint i;
  Shape *a,*b;
  real d,discount;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<margin){
    a=shapes(proxies(i)->a); b=shapes(proxies(i)->b);
    d=1.-proxies(i)->d/margin;
    //NORMALS ALWAYS GO FROM b TO a !!
    discount = 1.;
    if(!a->contactOrientation.isZero()){ //object has an 'allowed contact orientation'
      CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
      CHECK(a->contactOrientation.isNormalized(),"contact orientation is not normalized");
      //real theta = ::acos( proxies(i)->normal*a->contactOrientation);
      real theta = .5*( proxies(i)->normal*a->contactOrientation-1.);
      discount *= theta*theta;
    }
    if(!b->contactOrientation.isZero()){
      CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
      CHECK(a->contactOrientation.isNormalized(),"contact orientation is not normalized");
      //real theta = ::acos(-proxies(i)->normal*b->contactOrientation);
      real theta = .5*(-proxies(i)->normal*a->contactOrientation-1.);
      discount *= theta*theta;
    }
    if(!linear) x(0) += discount*d*d;
    else        x(0) += discount*d;
  }
}

//! gradient (=scalar Jacobian) of this contact cost
real ors::Graph::getContactGradient(arr &grad,real margin,bool linear){
  ors::Vector normal;
  uint i;
  Shape *a,*b;
  real d,discount;
  real cost=0.;
  arr J,dnormal;
  grad.resize(1,jd);
  grad.setZero();
  ors::Frame arel,brel;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<margin){
    a=shapes(proxies(i)->a); b=shapes(proxies(i)->b);
    d=1.-proxies(i)->d/margin;
    discount = 1.;
    if(!a->contactOrientation.isZero()){ //object has an 'allowed contact orientation'
      CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
      CHECK(a->contactOrientation.isNormalized(),"contact orientation is not normalized");
      //real theta = ::acos( proxies(i)->normal*a->contactOrientation);
      real theta = .5*( proxies(i)->normal*a->contactOrientation-1.);
      discount *= theta*theta;
    }
    if(!b->contactOrientation.isZero()){
      CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
      CHECK(a->contactOrientation.isNormalized(),"contact orientation is not normalized");
      //real theta = ::acos(-proxies(i)->normal*b->contactOrientation);
      real theta = .5*(-proxies(i)->normal*a->contactOrientation-1.);
      discount *= theta*theta;
    }
    if(!linear) cost += discount*d*d;
    else        cost += discount*d;
    
    arel.setZero();  arel.p=a->X.r/(proxies(i)->posA-a->X.p);
    brel.setZero();  brel.p=b->X.r/(proxies(i)->posB-b->X.p);
    
    CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
    dnormal.referTo(proxies(i)->normal.v,3); dnormal.reshape(1,3);
    if(!linear){
      jacobian(J,a->body->index,&arel); grad -= ((real)2.*discount*d)/margin*(dnormal*J);
      jacobian(J,b->body->index,&brel); grad += ((real)2.*discount*d)/margin*(dnormal*J);
    }else{
      jacobian(J,a->body->index,&arel); grad -= discount/margin*(dnormal*J);
      jacobian(J,b->body->index,&brel); grad += discount/margin*(dnormal*J);
    }
  }

  return cost;
}

//! measure (=scalar kinematics) for the contact cost summed over all bodies
void ors::Graph::getContactConstraints(arr& y){
  y.clear();
  uint i;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age){
    y.append(proxies(i)->d);
  }
}

//! gradient (=scalar Jacobian) of this contact cost
void ors::Graph::getContactConstraintsGradient(arr &dydq){
  dydq.clear();
  ors::Vector normal;
  uint i,con=0;
  Shape *a,*b;
  arr J,dnormal,grad(1,jd);
  ors::Frame arel,brel;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age){
    a=shapes(proxies(i)->a); b=shapes(proxies(i)->b);
    
    arel.setZero();  arel.p=a->X.r/(proxies(i)->posA-a->X.p);
    brel.setZero();  brel.p=b->X.r/(proxies(i)->posB-b->X.p);
    
    CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
    dnormal.referTo(proxies(i)->normal.v,3); dnormal.reshape(1,3);
    grad.setZero();
    jacobian(J,a->body->index,&arel); grad += dnormal*J; //moving a long normal b->a increases distance
    jacobian(J,b->body->index,&brel); grad -= dnormal*J; //moving b long normal b->a decreases distance
    dydq.append(grad);
    con++;
  }
  dydq.reshape(con,jd);
}


#if 0 //alternative implementation : cost=1 -> contact, other discounting...
real ors::Graph::getContactGradient(arr &grad,real margin){
  ors::Vector normal;
  uint i;
  Shape *a,*b;
  real d,discount;
  real cost=0.;
  arr J,dnormal;
  grad.resize(1,jd);
  grad.setZero();
  ors::Frame arel,brel;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<margin){
    a=shapes(proxies(i)->a); b=shapes(proxies(i)->b);
    discount = 1.;
    if(!a->contactOrientation.isZero()){ //object has an 'allowed contact orientation'
      CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
      CHECK(a->contactOrientation.isNormalized(),"contact orientation is not normalized");
      //real theta = ::acos( proxies(i)->normal*a->contactOrientation);
      real theta = .5*( proxies(i)->normal*a->contactOrientation-1.);
      discount *= theta*theta;
    }
    if(!b->contactOrientation.isZero()){
      CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
      CHECK(a->contactOrientation.isNormalized(),"contact orientation is not normalized");
      //real theta = ::acos(-proxies(i)->normal*b->contactOrientation);
      real theta = .5*(-proxies(i)->normal*a->contactOrientation-1.);
      discount *= theta*theta;
    }
    real marg=(discount+.1)*margin;
    d=1.-proxies(i)->d/marg;
    if(d<0.) continue;
    cost += d*d;
    
    arel.setZero();  arel.p=a->X.r/(proxies(i)->posA-a->X.p);
    brel.setZero();  brel.p=b->X.r/(proxies(i)->posB-b->X.p);
    
    CHECK(proxies(i)->normal.isNormalized(),"proxy normal is not normalized");
    dnormal.referTo(proxies(i)->normal.v,3); dnormal.reshape(1,3);
    jacobian(J,a->body->index,&arel); grad -= (2.*d/marg)*(dnormal*J);
    jacobian(J,b->body->index,&brel); grad += (2.*d/marg)*(dnormal*J);
  }

  return cost;
}
#endif
                                                  
void ors::Graph::getLimitsMeasure(arr &x,const arr& limits,real margin){
  CHECK(limits.d0==jd && limits.d1==2,"joint limits parameter size mismatch");
  x.resize(1);
  x=0.;
  arr q;
  getJointState(q);
  uint i;
  real d;
  for(i=0;i<q.N;i++){
    d = q(i) - limits(i,0);
    if(d<margin){  d-=margin;  x(0) += d*d;  }
    d = limits(i,1) - q(i);
    if(d<margin){  d-=margin;  x(0) += d*d;  }
  }
}

real ors::Graph::getLimitsGradient(arr &grad,const arr& limits,real margin){
  CHECK(limits.d0==jd && limits.d1==2,"");
  uint i;
  real d;
  real cost=0.;
  arr J;
  grad.resize(1,jd);
  grad.setZero();
  arr q;
  getJointState(q);
  for(i=0;i<q.N;i++){
    d = q(i) - limits(i,0);
    if(d<margin){  d-=margin;  grad(0,i) += 2.*d;  }
    d = limits(i,1) - q(i);
    if(d<margin){  d-=margin;  grad(0,i) -= 2.*d;  }
  }
  return cost;
}

//! center of mass of the whole configuration (3 vector)
real ors::Graph::getCenterOfMass(arr& x_){
  real M=0.;
  node n;
  uint j;
  ors::Vector x;
  x.setZero();
  for_list(j,n,bodies){
    M+=n->mass;
    x+=n->mass*n->X.p;
  }
  x/=M;
  x_.setCarray(x.v,3);
  return M;
}

//! gradient (Jacobian) of the COM w.r.t. q (3 x n tensor)
void ors::Graph::getComGradient(arr &grad){
  real M=0.;
  node n;
  uint j;
  arr J(3,getJointStateDimension());
  grad.resizeAs(J); grad.setZero();
  for_list(j,n,bodies){
    M += n->mass;
    jacobian(J,n->index);
    grad += n->mass * J;
  }
  grad/=M;
}

  /*!\brief returns a k-dim vector containing the penetration depths of all bodies */
void ors::Graph::getPenetrationState(arr &vec){
  vec.resize(bodies.N);
  vec.setZero();
  ors::Vector d;
  uint i;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<0.){
    d=proxies(i)->posB - proxies(i)->posA;

    if(proxies(i)->a!=-1) vec(proxies(i)->a) += d.length();
    if(proxies(i)->b!=-1) vec(proxies(i)->b) += d.length();
  }
}

ors::Proxy* ors::Graph::getContact(uint a,uint b){
  uint i;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<0.){
    if(proxies(i)->a==(int)a && proxies(i)->b==(int)b) return proxies(i);
    if(proxies(i)->a==(int)b && proxies(i)->b==(int)a) return proxies(i);
  }
  return NULL;
}

  /*!\brief a vector describing the incoming forces (penetrations) on one object */
void ors::Graph::getGripState(arr& grip,uint j){
  ors::Vector d,p;
  ors::Vector sumOfD; sumOfD.setZero();
  ors::Vector torque; torque.setZero();
  real sumOfAbsD = 0.;
  real varOfD = 0.;

  p.setZero();
  uint i,n=0;
  for(i=0;i<proxies.N;i++) if(!proxies(i)->age && proxies(i)->d<0.){
    if(proxies(i)->a!=(int)j && proxies(i)->b!=(int)j) continue;

    n++;

    if(proxies(i)->a==(int)j){
      d=proxies(i)->posB - proxies(i)->posA;
      p=proxies(i)->posA;
    }
    if(proxies(i)->b==(int)j){
      d=proxies(i)->posA - proxies(i)->posB;
      p=proxies(i)->posB;
    }

    sumOfAbsD += d.length();
    sumOfD    += d;
    varOfD    += d.lengthSqr();
    torque    += (p - bodies(j)->X.p) ^ d;

  }
  if(n){ varOfD = (varOfD - sumOfD*sumOfD) / n; }
  
  grip.resize(8);
  grip(0)=sumOfAbsD;
  grip(1)=varOfD;
  memmove(grip.p+2,sumOfD.v,3*sizeof(real));
  memmove(grip.p+5,torque.v,3*sizeof(real));
}

  //! returns the number of touch-sensors
uint ors::Graph::getTouchDimension(){
  node n;
  uint i=0,j;

  // count touchsensors
  for_list(j,n,bodies) if(anyListGet<real>(n->ats,"touchsensor",0)) i++;

  td=i;
  return i;
}

  //! returns the touch vector (penetrations) of all touch-sensors
void ors::Graph::getTouchState(arr& touch){
  if(!td) td=getTouchDimension();
  arr pen;
  getPenetrationState(pen);
  node n;
  uint i=0,j;
  for_list(j,n,bodies){
    if(anyListGet<real>(n->ats,"touchsensor",0)){
      touch(i)=pen(n->index);
      i++;
    }
  }
  CHECK(i==td,"");
}

  /*!\brief */
real ors::Graph::getEnergy(){
  node n;
  uint j;
  real m,v,E;
  ors::Matrix I;
  ors::Vector w;

  E=0.;
  for_list(j,n,bodies){
    m=n->mass;
    I=n->inertia;
    v=n->X.v.length();
    w=n->X.w;
    //I.setZero(); I(0,0)=I(1,1)=I(2,2)=.1*m;
    E += .5*m*v*v;
    E += 9.81 * m * n->X.p(2);
    E += .5*(w*(I*w));
  }

  return E;
}

// ------------------ end slGraph ---------------------


//===========================================================================
//
// helper routines -- in a classical C interface
//



/*!\brief get the center of mass, total velocity, and total angular momemtum */ 
void ors::Graph::getTotals(ors::Vector& c,ors::Vector& v,ors::Vector& l,ors::Quaternion& ori){
  node n;
  uint j;
  real m,M;
  
  //dMass mass;
  ors::Matrix ID;
  //ors::Matrix TP;
  ors::Vector r,o;

  ID.setId();
  c.setZero();
  v.setZero();
  l.setZero();
  o.setZero();
  //Iall.setZero();
  M=0.;
  for_list(j,n,bodies){
    l+=n->inertia*n->X.w;
    //TP.setTensorProduct(n->X.p,n->X.p);
    //Iall+=m*((n->X.p*n->X.p)*ID + TP);
    
    m=n->mass;
    l+=m*(n->X.p ^ n->X.v);
    o+=m*n->X.r.getVec(r);
    
    M+=m;
    c+=m*n->X.p;
    v+=m*n->X.v;
  }
  c/=M;
  v/=M;
  o/=M;
  ori.setVec(o);
}



#undef LEN

