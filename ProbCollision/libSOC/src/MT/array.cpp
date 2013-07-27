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

#include "array.h"
#include "util.h"

static real MT_SIGN_SVD(real a,real b){ return b>0 ? ::fabs(a) : -::fabs(a); }
#define MT_max_SVD(a,b) ( (a)>(b) ? (a) : (b) )
#define MT_SVD_MINVALUE .0 //1e-10
#ifndef MT_NOCHECK
//#  define MT_CHECK_INVERSE 1e-5
//#  define MT_CHECK_SVD 1e-5
#endif

namespace MT{
  bool useLapack=true;
  uint64 memoryTotal=0, memoryBound=10ull<<30; //this is 1GB
  bool memoryStrict=false;
}
//int ARRAYOLDREAD=0;

/* LAPACK notes
Use the documentation at
  http://www.netlib.org/lapack/double/
  http://www.netlib.org/lapack/individualroutines.html
to find the right function! Also use the man tools with Debian package lapack-doc installed.

I've put the clapack.h directly into the MT directory - one only has to link to the Fortran lib
*/


//===========================================================================
//
//!@name cblas and Lapack support
//

//===========================================================================
//
//!@name matrix operations
//

//! returns the identity matrix
arr Identity(uint dim){
  arr z;
  z.setId(dim);
  return z;
}

//! make symmetric \f$A=(A+A^T)/2\f$
void makeSymmetric(arr& A){
  CHECK(A.nd==2 && A.d0==A.d1,"not symmetric");
  uint n=A.d0,i,j;
  for(i=1;i<n;i++) for(j=0;j<i;j++) A(j,i) = A(i,j) = .5 * (A(i,j) + A(j,i));
}

//! make its transpose \f$A \gets A^T\f$
void transpose(arr& A){
  CHECK(A.nd==2 && A.d0==A.d1,"not symmetric");
  uint n=A.d0,i,j;
  real z;
  for(i=1;i<n;i++) for(j=0;j<i;j++){ z=A(j,i); A(j,i)=A(i,j); A(i,j)=z; }
}

namespace MT{
  //! use this to turn on Lapack routines [default true if MT_LAPACK is defined]
  extern bool useLapack;
}



//===========================================================================
//
//!@name SVD etc
//

//! called from svd if MT_LAPACK is not defined
uint own_SVD(const arr& A,
			    arr& U,
			    arr& w,
			    arr& V,
			    bool sort);

/*!\brief Singular Value Decomposition (from Numerical Recipes);
  computes \f$U, D, V\f$ with \f$A = U D V^T\f$ from \f$A\f$ such that
  \f$U\f$ and \f$V\f$ are orthogonal and \f$D\f$ diagonal (the
  returned array d is 1-dimensional) -- uses LAPACK if MT_LAPACK is
  defined */
uint svd(const arr& A,arr& U,arr& d,arr& V,bool sort){
  uint r;
#ifdef MT_LAPACK
  if(MT::useLapack){
    r=lapack_SVD(A,U,d,V);
    V=~V;
  }else{
    r=own_SVD(A,U,d,V,sort);
  }
#else
  r=own_SVD(A,U,d,V,sort);
#endif
  
#ifdef MT_CHECK_SVD
  bool uselapack=MT::useLapack;
  MT::useLapack=false;
  real err;
  arr dD,I;
  setDiagonal(dD,d);
  //cout <<U <<dD <<Vt;
  //Atmp = V * D * U;
  arr Atmp;
  Atmp = U * dD * ~V;
  //cout <<"\nA=" <<A <<"\nAtmp=" <<Atmp <<"U=" <<U <<"W=" <<dD <<"~V=" <<~V <<endl;
  std::cout <<"SVD is correct:  " <<(err=maxDiff(Atmp,A)) <<' ' <<endl;    CHECK(err<MT_CHECK_SVD,"");
  if(A.d0<=A.d1){
    I.setId(U.d0);
    std::cout <<"U is orthogonal: " <<(err=maxDiff(U * ~U,I)) <<' ' <<endl;  CHECK(err<MT_CHECK_SVD,"");
    I.setId(V.d1);
    std::cout <<"V is orthogonal: " <<(err=maxDiff(~V * V,I)) <<endl;        CHECK(err<MT_CHECK_SVD,"");
  }else{
    I.setId(U.d1);
    std::cout <<"U is orthogonal: " <<(err=maxDiff(~U * U,I)) <<' ' <<endl;  CHECK(err<MT_CHECK_SVD,"");
    I.setId(V.d0);
    std::cout <<"V is orthogonal: " <<(err=sqrDistance(V * ~V,I)) <<endl;        CHECK(err<1e-5,"");
  }
  MT::useLapack=uselapack;
#endif
  
  return r;
}

//! gives a decomposition \f$A = U V^T\f$
void svd(const arr& A,arr& U,arr& V){
  arr d,D;
  ::svd(A,U,d,V);
  D.resize(d.N,d.N); D=0.;
  for(uint i=0;i<d.N;i++) D(i,i)=::sqrt(d(i));
  U=U*D;
  V=V*D;
  //CHECK(maxDiff(A,U*~V) <1e-4,"");
}

void check_inverse(const arr& invA,const arr& A){
#ifdef MT_CHECK_INVERSE
  arr D,_D; D.setId(A.d0);
  uint me;
  _D=A*invA;
  real err=maxDiff(_D,D,&me);
  cout <<"inverse is correct: " <<err <<endl;
  if(A.d0<10){
    CHECK(err<MT_CHECK_INVERSE ,"inverting failed, error="<<err <<" " <<_D.elem(me) <<"!=" <<D.elem(me) <<"\nA="<<A <<"\ninvA=" <<invA <<"\nA*invA=" <<_D);
  }else{
    CHECK(err<MT_CHECK_INVERSE ,"inverting failed, error="<<err <<" " <<_D.elem(me) <<"!=" <<D.elem(me) );
  }
#endif
}

uint inverse(arr& invA,const arr& A){
  uint r=inverse_SVD(invA,A);
  //MT::inverse_LU(inverse,A); return A.d0;
  return r;
}

//! calls inverse(B,A) and returns B
arr inverse(const arr& A){ arr B; inverse(B,A); return B; }

//! Pseudo Inverse based on SVD; computes \f$B\f$ such that \f$ABA = A\f$
uint inverse_SVD(arr& invA,const arr& A){
  unsigned i,j,k,m=A.d0,n=A.d1,r;
  arr U,V,w,winv;
  invA.resize(n,m);
  if(m==0 || n==0) return 0;
  if(m==n && m==1){ invA(0,0)=1./A(0,0); return 0; }
  if(m==n && m==2){ inverse2d(invA,A); return 0; }

  r=svd(A,U,w,V,true);
  
  //arr W;
  //setDiagonal(W,w);
  //CHECK(fabs(maxDiff(A,U*W*~V))<1e-10,"");
  
  winv.resizeAs(w);
  for(i=0;i<r;i++){
    if(w(i)>1e-10) winv(i) = 1./w(i); else winv(i) = 1e10;
  }
  for(   ;i<w.N;i++) winv(i) = 0.;

#if 0
  //arr W;
  setDiagonal(W,winv);
  invA = V * W * ~U;
#else
  real *invAij=&invA(0,0);
  for(i=0;i<n;i++) for(j=0;j<m;j++){
    real* vi = &V(i,0);
    real* uj = &U(j,0);
    real  t  = 0.;
    for(k=0;k<w.N;k++) t += vi[k] * winv.p[k] * uj[k];
    *invAij = t;
    invAij++;
  }
#endif

#ifdef MT_CHECK_INVERSE
  check_inverse(invA,A);
#endif
  return r;
}

void inverse_LU(arr& invX,const arr& X){
  NIY;
#if 0
  CHECK(X.nd==2 && X.d0==X.d1,"");
  uint n=X.d0,i,j;
  invX.resize(n,n);
  if(n==0) return;
  if(n==n && n==1){ invX(0,0)=1./X(0,0); return; }
  if(n==n && n==2){ inverse2d(invX,X); return; }
  arr LU,piv;
  lapackLU(X,LU,piv);
  arr col(n);
  for(j=0;j<n;j++){
    col.setZero();
    col(j)=1.0;
    lubksb(LU.pp,n,idx,col.p);
    for(i=0;i<n;i++) invX(i,j)=col(i);
  }

  delete[] idx;
  delete[] d;
  
#ifdef MT_CHECK_INVERSE
  check_inverse(invX,X);
#endif
#endif
}

void inverse_SymPosDef(arr& invA,const arr& A){
  CHECK(A.d0==A.d1,"");
#ifdef MT_LAPACK
  lapack_inverseSymPosDef(invA,A);
#else
  inverse_SVD(invA,A);
#endif
#ifdef MT_CHECK_INVERSE
  check_inverse(invA,A);
#endif
}

void pseudoInverse(arr& invA,const arr& A,const arr& invW,real eps){
  arr tA,E,invAWA;
  transpose(tA,A);
  if(!eps){
    inverse_SymPosDef(invAWA,A*invW*tA);
  }else{
    E.setDiag(eps,A.d0);
    inverse_SymPosDef(invAWA,E+A*invW*tA);
  }
  invA = invW * tA * invAWA;
}

//! the determinant of a 2D squared matrix
real determinant(const arr& A);

/*!\brief the cofactor is the determinant of a 2D squared matrix after removing
  the ith row and the jth column */
real cofactor(const arr& A,uint i,uint j);

void gaussFromData(arr& a,arr& A,const arr& X){
  CHECK(X.nd==2,"");
  uint N=X.d0,n=X.d1;
  arr ones(N); ones=1.;
  a = ones*X/(real)N; a.reshape(n);
  A = (~X*X)/(real)N - (a^a);
}

/* compute a rotation matrix that rotates a onto v in arbitrary dimensions */
void rotationFromAtoB(arr& R,const arr& a,const arr& v){
  CHECK(a.N==v.N,"");
  CHECK(fabs(1.-norm(a))<1e-10 && fabs(1.-norm(v))<1e-10,"");
  uint n=a.N,i,j;
  if(maxDiff(a,v)<1e-10){ R.setId(n); return; } //nothing to rotate!!
  R.resize(n,n);
  //-- compute b orthogonal to a such that (a,b) span the rotation plane
  arr b;
  b = v - a*scalarProduct(a,v);
  b /= norm(b);
  //-- compute rotation coefficients within the (a,b) plane, namely, R_2D=(v_a  -v_b ;  v_b  v_a)
  real v_a,v_b;
  v_a=scalarProduct(v,a);     //component along a
  v_b=scalarProduct(v,b);     //component along b
  //-- compute columns of R:
  arr x(n),x_res;
  real x_a,x_b;
  for(i=0;i<n;i++){
    x.setZero(); x(i)=1.;       //x=i-th unit vector
    x_a=scalarProduct(x,a);     //component along a
    x_b=scalarProduct(x,b);     //component along b
    x_res = x - x_a*a - x_b*b;  //residual (rest) of the vector
    //rotated vector = residual + rotated a-component + rotated b-component
    x = x_res + (v_a*x_a-v_b*x_b)*a + (v_b*x_a+v_a*x_b)*b;
    for(j=0;j<n;j++) R(j,i)=x(j); //store as column of the final rotation
  }
}

//===========================================================================
//
//!@name gnuplot fun
//

//! calls gnuplot to display the (n,2) or (n,3) array (n=number of points of line or surface)
void gnuplot(const arr& X);

//! write 2 arrays in parallel columns in a file
void write(const arr& X,const arr& Y,const char* name);

//! write 3 arrays in parallel columns in a file
void write(const arr& X,const arr& Y,const arr& Z,const char* name);


//===========================================================================
//
//!@name simple image formats
//

/*! save data as ppm or pgm. Images are (height,width,[0,2,3,4])-dim
  byte arrays, where the 3rd dimension determines whether it's a grey
  (0), grey-alpha (2), RGB (3), or RGBA (4) image */
void write_ppm(const byteA &img,const char *file_name,bool swap_rows);

/*! read data from an ppm or pgm file */
void read_ppm(byteA &img,const char *file_name,bool swap_rows);

//! add an alpha channel to an image array
void add_alpha_channel(byteA &img,byte alpha);

//! make grey scale image
void make_grey(byteA &img);

//! make a grey image and RGA image
void make_RGB(byteA &img);


//===========================================================================
//
// Array class
//



uint own_SVD(const arr& A,
			    arr& U,
			    arr& w,
			    arr& V,
			    bool sort){
  //MT::Array<real*> Apointers,Upointers,Vpointers;
  unsigned m = A.d0; /* rows */
  unsigned n = A.d1; /* cols */
  U.resize(m,n);
  V.resize(n,n);
  w.resize(n);
  real **a = A.getCarray(); //Pointers(Apointers); /* input matrix */
  real **u = U.getCarray(); //Pointers(Upointers); /* left vectors */
  real **v = V.getCarray(); //Pointers(Vpointers); /* right vectors */

  int flag;
  unsigned i,its,j,jj,k,l,nm(0),r;
  real anorm,c,f,g,h,s,scale,x,y,z,t;

  arr rv1(n);

  /* copy A to U */
  for(i=0;i<m;i++) for(j=0;j<n;j++) u[i][j] = a[i][j];

  /* householder reduction to pickBiagonal form */
  g = scale = anorm = 0.0;

  for(i=0;i<n;i++){
    l = i + 1;
    rv1(i) = scale * g;
    g = s = scale = 0.0;

    if(i<m){
      for(k=i;k<m;k++) scale += fabs(u[k][i]);

      if(scale!=0.0){
	for(k=i;k<m;k++){
	  u[k][i] /= scale;
	  s += u[k][i] * u[k][i];
	}

	f = u[i][i];
	g = -MT_SIGN_SVD(sqrt(s),f);
	h = f * g - s;
	u[i][i] = f - g;

	for(j=l;j<n;j++){
	  s = 0.0;
	  for(k=i;k<m;k++) s += u[k][i] * u[k][j];

	  f = s / h;
	  for(k=i;k<m;k++) u[k][j] += f * u[k][i];
	}

	for(k=i;k<m;k++) u[k][i] *= scale;
      }
    }

    w(i) = scale * g;
    g = s = scale = 0.0;

    if(i<m && i!=n-1){
      for(k=l;k<n;k++)scale += fabs(u[i][k]);

      if(scale!=0.0){
	for(k=l;k<n;k++){
	  u[i][k] /= scale;
	  s += u[i][k] * u[i][k];
	}

	f = u[i][l];
	g = -MT_SIGN_SVD(sqrt(s),f);
	h = f * g - s;
	u[i][l] = f - g;

	for(k=l;k<n;k++) rv1(k) = u[i][k] / h;

	for(j=l;j<m;j++){
	  s = 0.0;
	  for(k=l;k<n;k++) s += u[j][k] * u[i][k];

	  for(k=l;k<n;k++) u[j][k] += s * rv1(k);
	}

	for(k=l;k<n;k++) u[i][k] *= scale;
      }
    }

    anorm = MT_max_SVD(anorm,fabs(w(i)) + fabs(rv1(i)));
  }

  /* accumulation of right-hand transformations */
  for(l=i=n;i--;l--){
    if(l<n){
      if(g!=0.0){
	/* double division avoids possible underflow */
	for(j=l;j<n;j++) v[j][i] = (u[i][j] / u[i][l]) / g;

	for(j=l;j<n;j++){
	  s = 0.0;
	  for(k=l;k<n;k++) s += u[i][k] * v[k][j];

	  for(k=l;k<n;k++) v[k][j] += s * v[k][i];
	}
      }

      for(j=l;j<n;j++) v[i][j] = v[j][i] = 0.0;
    }

    v[i][i] = 1.0;
    g = rv1(i);
  }

  /* accumulation of left-hand transformations */
  for(l=i=(m<n?m:n);i--;l--){
    g = w(i);

    for(j=l;j<n;j++) u[i][j] = 0.0;

    if(g!=0.0){
      g = 1.0 / g;

      for(j=l;j<n;j++){
	s = 0.0;
	for(k=l;k<m;k++) s += u[k][i] * u[k][j];

	/* double division avoids possible underflow */
	f = (s / u[i][i]) * g;

	for(k=i;k<m;k++) u[k][j] += f * u[k][i];
      }

      for(j=i;j<m;j++) u[j][i] *= g;
    }else{
      for(j=i;j<m;j++) u[j][i] = 0.0;
    }

    u[i][i]++;
  }

  /* diagonalization of the pickBiagonal form */
  for(k=n;k--;){
    for(its=1;its<=30;its++){
      flag = 1;

      /* test for splitting */
      for(l = k + 1;l--;){
	/* rv1 [0] is always zero, so there is no exit */
	nm = l - 1;

	if(fabs(rv1(l)) + anorm == anorm){
	  flag = 0;
	  break;
	}

	//if(!l) break; //(mt 07-01-16)
	if(fabs(w(nm)) + anorm == anorm) break;
      }

      if(flag){
	/* cancellation of rv1 [l] if l greater than 0 */
	c = 0.0;
	s = 1.0;

	for(i=l;i<=k;i++){
	  f = s * rv1(i);
	  rv1(i) *= c;

	  if(fabs(f) + anorm == anorm) break;

	  g = w(i);
	  h = hypot(f,g);
	  w(i) = h;
	  h = 1.0 / h;
	  c = g * h;
	  s = -f * h;

	  for(j=0;j<m;j++){
	    y = u[j][nm];
	    z = u[j][i];
	    u[j][nm] = y * c + z * s;
	    u[j][i] = z * c - y * s;
	  }
	}
      }

      /* test for convergence */
      z = w(k);

      if(l==k){
	if(z<0.0){
	  w(k) = -z;
	  for(j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }

      if(its==50) HALT("svd failed");
      //if(its==30) throw k;

      /* shift from bottom 2 by 2 minor */
      x = w(l);
      nm = k - 1;
      y = w(nm);
      g = rv1(nm);
      h = rv1(k);
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = hypot(f,1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + MT_SIGN_SVD(g,f))) - h)) / x;

      /* next qr transformation */
      c = s = 1.0;

      for(j=l;j<k;j++){
	i = j + 1;
	g = rv1(i);
	y = w(i);
	h = s * g;
	g *= c;
	z = hypot(f,h);
	rv1(j) = z;
	c = f / z;
	s = h / z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;

	for(jj=0;jj<n;jj++){
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x * c + z * s;
	  v[jj][i] = z * c - x * s;
	}

	z = hypot(f,h);
	w(j) = z;

	/* rotation can be arbitrary if z is zero */
	if(z!=0.0){
	  z = 1.0 / z;
	  c = f * z;
	  s = h * z;
	}

	f = c * g + s * y;
	x = c * y - s * g;

	for(jj=0;jj<m;jj++){
	  y = u[jj][j];
	  z = u[jj][i];
	  u[jj][j] = y * c + z * s;
	  u[jj][i] = z * c - y * s;
	}
      }

      rv1(l) = 0.0;
      rv1(k) = f;
      w(k) = x;
    }
  }

  //sorting:
  if(sort){
    unsigned i,j,k;
    real   p;

    for(i=0;i<n-1;i++){
      p = w(k=i);

      for(j=i+1;j<n;j++) if(w(j)>=p) p = w(k=j);

      if(k!=i){
	w(k) = w(i);
	w(i) = p;

	for(j=0;j<n;j++){
	  p       = v[j][i];
	  v[j][i] = v[j][k];
	  v[j][k] = p;
	}

	for(j=0;j<m;j++){
	  p       = u[j][i];
	  u[j][i] = u[j][k];
	  u[j][k] = p;
	}
      }
    }
  }

  //rank analysis

  for(r=0;r<n && w(r)>MT_SVD_MINVALUE;r++){};

  t = r < n ? fabs(w(n-1)) : 0.0;
  r = 0;
  s = 0.0;
  while(r<n && w(r)>t && w(r)+s>s) s += w(r++);

  return r;
}

real _determinant(real **A,uint n){
  if(n==1) return A[0][0];
  if(n==2) return A[0][0]*A[1][1]-A[0][1]*A[1][0];
  uint i,j;
  real d=0;
  real **B=new real*[n-1];
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(j<i) B[j]=&A[j][1];
      if(j>i) B[j-1]=&A[j][1];
    }
    d+=((i&1)?-1.:1.) * A[i][0] * _determinant(B,n-1);
  }
  delete[] B;
  return d;
}

real determinant(const arr& A){
  CHECK(A.nd==2 && A.d0==A.d1,"determinants require a squared 2D matrix");
  //MT::Array<real*> B;
  A.getCarray(); //Pointers(B);
  return _determinant(A.pp,A.d0);
}

real cofactor(const arr& A,uint i,uint j){
  CHECK(A.nd==2 && A.d0==A.d1,"determinants require a squared 2D matrix");
  arr B=A;
  B.delRow(i);
  B.delColumns(j,1);
  return ((i&1)^(j&1)?-1.:1) * determinant(B);
}

void SUS(const arr& p,uint n,uintA& s){
  //following T. Baeck "EA in Theo. and Prac." p120
  s.resize(n);
  real sum=0,ptr=MT::rnd.uni();
  uint i,j=0;
  for(i=0;i<p.N;i++){
    sum+=p(i)*n;
    while(sum>ptr){ s(j)=i; j++; ptr+=1.; }
  }
  //now, 'sum' should = 'n' and 'ptr' has been 'n'-times increased -> 'j=n'
  CHECK(j==n,"error in rnd::SUS(p,n,s) -> p not normalized?");
}

uint SUS(const arr& p){
  real sum=0,ptr=MT::rnd.uni();
  uint i;
  for(i=0;i<p.N;i++){
    sum+=p(i);
    if(sum>ptr) return i;
  }
  HALT("error in rnd::SUS(p) -> p not normalized? " <<p);
  return 0;
}

void gnuplot(const arr& X){
  if(X.nd==2 && X.d1!=2){ //assume array -> splot
    MT::IOraw=true;
    MT::save(X,"z.pltX");
    gnuplot("splot 'z.pltX' matrix with pm3d,'z.pltX' matrix with lines");
    MT::IOraw=false;
    return;
  }
  if(X.nd==2 && X.d1==2){ //assume curve -> plot
    MT::IOraw=true;
    MT::save(X,"z.pltX");
    gnuplot("plot 'z.pltX' us 1:2");
    MT::IOraw=false;
    return;
  }
  if(X.nd==1){ //assume curve -> plot
    MT::IOraw=true;
    arr Y;
    Y.referTo(X);
    Y.resize(Y.N,1);
    MT::save(Y,"z.pltX");
    gnuplot("plot 'z.pltX' us 1");
    MT::IOraw=false;
    return;
  }
}

void write(const arr& X,const arr& Y,const char* name){
  std::ofstream os;
  MT::open(os,name);
  MT::IOraw=true;
  uint i,j;
  if(X.nd==1){
    for(i=0;i<X.N;i++) os <<X(i) <<' ' <<Y(i) <<std::endl;
  }
  if(X.nd==2){
    for(i=0;i<X.d0;i++){
      for(j=0;j<X[i].N;j++) os <<X[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Y[i].elem(j) <<' ';
      os <<std::endl;
    }
  }
}

void write(const arr& X,const arr& Y,const arr& Z,const char* name){
  std::ofstream os;
  MT::open(os,name);
  MT::IOraw=true;
  uint i,j;
  if(X.nd==1){
    for(i=0;i<X.N;i++) os <<X(i) <<' ' <<Y(i) <<' ' <<Z(i) <<std::endl;
  }
  if(X.nd==2){
    for(i=0;i<X.d0;i++){
      for(j=0;j<X[i].N;j++) os <<X[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Y[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Z[i].elem(j) <<' ';
      os <<std::endl;
    }
  }
}

void write(const arr& X,const arr& Y,const arr& Z,const arr& A,const char* name){
  std::ofstream os;
  MT::open(os,name);
  MT::IOraw=true;
  uint i,j;
  if(X.nd==1){
    for(i=0;i<X.N;i++) os <<X(i) <<' ' <<Y(i) <<' ' <<Z(i) <<' ' <<A(i) <<std::endl;
  }
  if(X.nd==2){
    for(i=0;i<X.d0;i++){
      for(j=0;j<X[i].N;j++) os <<X[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Y[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Z[i].elem(j) <<' ';
      for(j=0;j<A[i].N;j++) os <<A[i].elem(j) <<' ';
      os <<std::endl;
    }
  }
}

void write(const arr& X,const arr& Y,const arr& Z,const arr& A,const arr& B,const char* name){
  std::ofstream os;
  MT::open(os,name);
  MT::IOraw=true;
  uint i,j;
  if(X.nd==1){
    for(i=0;i<X.N;i++) os <<X(i) <<' ' <<Y(i) <<' ' <<Z(i) <<' ' <<A(i) <<' ' <<B(i) <<std::endl;
  }
  if(X.nd==2){
    for(i=0;i<X.d0;i++){
      for(j=0;j<X[i].N;j++) os <<X[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Y[i].elem(j) <<' ';
      for(j=0;j<Y[i].N;j++) os <<Z[i].elem(j) <<' ';
      for(j=0;j<A[i].N;j++) os <<A[i].elem(j) <<' ';
      for(j=0;j<B[i].N;j++) os <<B[i].elem(j) <<' ';
      os <<std::endl;
    }
  }
}

void write_ppm(const byteA &img,const char *file_name,bool swap_rows){
  CHECK(img.nd==2 || (img.nd==3 && img.d2==3),"only rgb or gray images to ppm");
  ofstream os;
  os.open(file_name, std::ios::out | std::ios::binary);
  if(!os.good()) HALT("could not open file `" <<file_name <<"' for output");
  switch(img.d2){
  case 0:  os <<"P5 "<<img.d1 <<' ' <<img.d0 <<" 255\n";  break; //PGM
  case 3:  os <<"P6 "<<img.d1 <<' ' <<img.d0 <<" 255\n";  break; //PPM
  }
  if(!swap_rows){
    os.write((char*)img.p,img.N);
  }else{
    for(uint i=img.d0;i--;) os.write((char*)&img(i,0,0),img.d1*img.d2);
  }
}

void read_ppm(byteA &img,const char *file_name,bool swap_rows){
  uint mode, width, height, max;
  ifstream is;
  is.open(file_name, std::ios::in | std::ios::binary);
  if(!is.good()) HALT("could not open file `" <<file_name <<"' for input");
  if(is.get()!='P') HALT("NO PPM FILE:" <<file_name);
  is >>mode;
  if(MT::peerNextChar(is)=='#') MT::skipLine(is);
  is >>width >>height >>max; 
  is.get(); //MUST be a white character if everything went ok
  switch(mode){
  case 5:  img.resize(height,width);    break; //PGM
  case 6:  img.resize(height,width,3);  break; //PPM
  }
  if(!swap_rows){
    is.read((char*)img.p,img.N);
  }else{
    for(uint i=img.d0;i--;) is.read((char*)&img(i,0,0),img.d1*img.d2);
  }
}

void add_alpha_channel(byteA &img,byte alpha){
  uint w=img.d1,h=img.d0;
  img.reshape(h*w,3);
  img.insColumns(3,1);
  for(uint i=0;i<img.d0;i++) img(i,3)=alpha;
  img.reshape(h,w,4);
}

void flip_image(byteA &img){
  uint h=img.d0,n=img.N/img.d0;
  byteA line(n);
  byte *a,*b,*c;
  for(uint i=0;i<h/2;i++){
    a=img.p+i*n;
    b=img.p+(h-1-i)*n;
    c=line.p;
    memmove(c,a,n);
    memmove(a,b,n);
    memmove(b,c,n);
  }
}

void make_grey(byteA &img){
  CHECK(img.nd==3 && (img.d2==3 || img.d1==4),"makeGray requires color image as input");
  byteA tmp;
  tmp.resize(img.d0,img.d1);
  for(uint i=0;i<img.d0;i++) for(uint j=0;j<img.d1;j++){
    tmp(i,j) = ((uint)img(i,j,0) + img(i,j,1) + img(i,j,2))/3; 
  }
  img=tmp;
}

void make_RGB(byteA &img){
  CHECK(img.nd==2,"make_RGB requires grey image as input");
  byteA tmp;
  tmp.resize(img.d0,img.d1,3);
  for(uint i=0;i<img.d0;i++) for(uint j=0;j<img.d1;j++){
    tmp(i,j,0) = img(i,j);
    tmp(i,j,1) = img(i,j);
    tmp(i,j,2) = img(i,j);
  }
  img=tmp;
}

#ifdef MT_EXPRESSIONS
void assign(arr& x){
  CHECK(x.ex,"self-assignment only if it is an expression");
  MT::Ex *e=x.ex;
  x.init();
  x.ex=e;
  assign(x,x);
  delete x.ex;
  x.ex=0;
}

void assign(arr& x,const arr& a){
  if(!a.ex){ x=a; return; }
  MT::Ex &e=*a.ex;
  if(e.op==MT::UNI){
    arr *A=(arr*)e.A;
    if(A->ex) assign(*A);
    if(!e.trans && e.mul==1 && e.add==0 ){ x=*A; return; }
    if(!e.trans && e.mul==1 ){ scalarPlus(x,*A,*((real*)&e.add)); return; }
    if(!e.trans && e.add==0 ){ scalarMultiplication(x,*A,*((real*)&e.mul)); return; }
    if(e.mul==1 && e.add==0 ){ transpose(x,*A); return; }
    HALT("");
  }else{
    arr *A=(arr*)e.A,*B=(arr*)e.B;
    if(A->ex) assign(*A);
    if(B->ex) assign(*B);
    //bool at,bt;
    //real ac,bc,ap,bp;
    switch(e.op){
    case MT::PROD:
      if(!A->ex && !B->ex){ innerProduct(x,*A,*B); return; }
      HALT("prod");
      break;
    case MT::MUL:
      if(!A->ex && !B->ex){ mult(x,*A,*B); return; }
      HALT("mult");
      break;
    case MT::Div:
      if(!A->ex && !B->ex){ div(x,*A,*B); return; }
      HALT("mult");
      break;
    case MT::OUT:
      if(!A->ex && !B->ex){ outerProduct(x,*A,*B); return; }
      HALT("out");
      break;
    case MT::PLUS:
      if(!A->ex && !B->ex){ plus(x,*A,*B); return; }
      //if(A->ex){ ap=A->ex->add; ac=A->ex->mul; at=A->ex->trans; A=(arr*)A->ex->A; }else{ ap=0; ac=1; at=false; }
      //if(B->ex){ bp=B->ex->add; bc=B->ex->mul; bt=B->ex->trans; B=(arr*)B->ex->A; }else{ bp=0; bc=1; bt=false; }
      //if(!at && !bt && !ap && !bp){ plus(x,ac,*A,bc,*B); return; }
      //if(!at && !bt && !B){ scalarPlus(x,*A,bc); return; }
      HALT("plus");
      break;
    case MT::MINUS:
      if(!A->ex && !B->ex){ minus(x,*A,*B); return; }
      //if(A->ex){ ap=A->ex->add; ac=A->ex->mul; at=A->ex->trans; A=(arr*)A->ex->A; }else{ ap=0; ac=1; at=false; }
      //if(B->ex){ bp=B->ex->add; bc=B->ex->mul; bt=B->ex->trans; B=(arr*)B->ex->A; }else{ bp=0; bc=1; bt=false; }
      //if(!at && !bt && !ap && !bp){ plus(x,ac,*A,-bc,*B); return; }
      //if(!at && !bt && !B){ scalarPlus(x,*A,bc); return; }
      HALT("minus");
      break;
    case MT::UNI:
      HALT("shouldn't be here!");
      break;
    }
    HALT("yet undefined expression");
  }
}
#endif



void getIndexTuple(uintA &I,uint i,const uintA &d){
  uint j;
  CHECK(i<product(d),"out of range");
  I.resize(d.N);
  I.setZero();
  for(j=d.N;j--;){
    I.p[j] = i%d.p[j];
    i -= I.p[j];
    i /= d.p[j];
  }
}

void lognormScale(arr& P,real& logP,bool force){
#ifdef MT_NoLognormScale
  return;
#endif
  real Z=0.;
  for(uint i=0;i<P.N;i++) Z += fabs(P.elem(i));
  if(!force && Z>1e-3 && Z<1e3) return;
  if(fabs(Z-1.)<1e-10) return;
  if(Z>1e-100){
    logP+=::log(Z);
    P/=Z;
  }else{
    logP+=::log(Z);
    P=1.;
    MT_MSG("ill-conditioned table factor for norm scaling");
  }
}

void sparseProduct(arr& y,arr& A,const arr& x){
  CHECK(x.nd==1 && A.nd==2 && x.d0==A.d1,"not a proper matrix multiplication");
  if(!A.sparse && !x.sparse){
    innerProduct(y,A,x);
    return;
  }
  if(A.sparse && !x.sparse){
    uint i,j,*k,*kstop;
    y.resize(A.d0); y.setZero();
    real *Ap=A.p;
    uintA *elems=A.sparse;
    for(k=elems->p, kstop=elems->pstop; k!=kstop; Ap++){
      i=*k; k++;
      j=*k; k++;
      y.p[i] += (*Ap) * x.p[j];
    }
    return;
  }
  if(A.sparse && x.sparse){
    uint i,j,n,*k,*kstop,*l,*lstop;
    y.clear(); y.nd=1; y.d0=A.d0; y.sparse=new uintA [2]; y.sparse[1].resize(y.d0); y.sparse[1]=(uint)-1;
    real *xp=x.p;
    uintA *elems,*col;
    elems=x.sparse;
    uint *slot;
    for(k=elems->p, kstop=elems->pstop; k!=kstop; xp++){
      j=*k; k++;
      col=A.sparse+(1+j);
      for(l=col->p, lstop=col->pstop;l!=lstop;){
	i =*l; l++;
	n =*l; l++;
	slot=&y.sparse[1](i);
	if(*slot==(uint)-1){
	  *slot=y.N;
	  y.resizeMEM(y.N+1,true); y(y.N-1)=0.;
	  y.sparse[0].append(i);
	  CHECK(y.sparse[0].N==y.N,"");
	}
	i=*slot;
	y(i) += A.elem(n) * (*xp);
      }
    }
    return;
  }
  if(!A.sparse && x.sparse){
    uint i,j,*k,*kstop,d1=A.d1;
    y.resize(A.d0); y.setZero();
    real *xp=x.p;
    uintA *elems;
    elems=x.sparse;
    for(k=elems->p, kstop=elems->pstop; k!=kstop; xp++){
      j=*k; k++;
      for(i=0; i<A.d0; i++){
        y.p[i] += A.p[i*d1+j] * (*xp);
      }
    }
    return;
  }
}

void scanArrFile(const char* name){
  ifstream is(name,std::ios::binary);
  CHECK(is.good(),"couldn't open file " <<name);
  arr x;
  String tag;
  for(;;){
    tag.read(is," \n\r\t"," \n\r\t");
    if(!is.good() || tag.N()==0) return;
    x.readTagged(is,NULL);
    x.writeTagged(cout,tag);  cout <<endl;
    if(!is.good()) return;
  }
}

//===========================================================================
//
// lists
//

void anyListRead(AnyList& ats,std::istream& is){
  char c,delim;
  MT::String tag,str;
  real d;
  arr reals;
  MT::Array<MT::String> strings;
  
  //read all generic attributes
  for(;;){
    tag.read(is," \t\r\n"," \t=},;([\n\r",false);
    if(!tag.N()){
      MT::skip(is," \t\r\n;");
      is.clear();
      break;
    }
    MT::skip(is);  is.get(c);
    if(c=='='){ MT::skip(is); is.get(c); }
    switch(c){
      case '(':{  //vector of strings
        delim=c;
        strings.clear();
        for(;;){
          MT::skip(is);
          is.get(c);
          if(c==')') break; else is.putback(c);
          str.read(is,"","), \t\r\n",false);
          strings.append(str);
        }
        if(strings.N==1){ //not nice - one should clearly distinguish between a vector and scalar...
          ats.append(anyNew<MT::String>(tag,strings(0)));
        }else{
          ats.append(anyNew<MT::String>(tag,strings.p,strings.N,delim));
        }
      }break;
      case '[':{  //vector of reals
        delim=c;
        reals.clear();
        for(;;){
          MT::skip(is);
          is.get(c);
          if(c==']' || c==')') break; else is.putback(c);
          is >>d;
          reals.append(d);
          if(!is.good()) HALT("ERROR");
        }
        if(reals.N==1){ //not nice - one should clearly distinguish between a vector and scalar...
          ats.append(anyNew<real>(tag,reals(0)));
        }else{
          ats.append(anyNew<real>(tag,reals.p,reals.N,delim));
        }
      }break;
      case '\'':{  //strings
        str.read(is,"","\'",true);
        ats.append(anyNew<MT::String>(tag,str));
      }break;
      case '\"':{
        str.read(is,"","\"",true);
        ats.append(anyNew<MT::String>(tag,str));
      }break;
      case '<':{
        is.putback(c);
        str.read(is,"",">",true);
        str <<'>';
        ats.append(anyNew<MT::String>(tag,str));
      }break;
      default:{
        is.putback(c);
        if(MT::contains(".0123456789",c)){ //single real
          is >>d;
          ats.append(anyNew<real>(tag,d));
        }else{ //bool
          ats.append(anyNew<real>(tag,(real*)0,0,0));
        }
      }break;
    }
    MT::skip(is," \n\t,");
  }
}

