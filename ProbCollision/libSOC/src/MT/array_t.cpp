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

#ifndef MT_array_t_cpp
#define MT_array_t_cpp

#include "array.h"
#include "util.h"
#include <math.h>
#include <algorithm>

#define maxRank 30


//===========================================================================
//
// Array class
//

template<class T> char MT::Array<T>::memMoveInit=-1;
template<class T> int MT::Array<T>::sizeT=-1;

/*!\brief Simple array container to store arbitrary-dimensional arrays
  (tensors); can buffer more memory than necessary for faster
  resize; enables non-const reference of subarrays; enables fast
  memove for elementary types; implements many standard
  array/matrix/tensor operations. Please see the fully public attributes at the
  bottom of this page -- everthing is meant to be perfectly
  transparent. Interfacing with ordinary C-buffers is simple,
  e.g. via \c Array::referTo (Cbuffer,size) and \c Array::p and \c
  Array::pp. Please see also the reference for the \ref array.h
  header, which contains lots of functions that can be applied on
  Arrays. */

template<class T> void MT::Array<T>::init(){
  reference=false;
  memMove=false;
  if(sizeT==-1) sizeT=sizeof(T);
  if(memMoveInit==-1){
    memMoveInit=0;
    if(typeid(T)==typeid(bool) ||
       typeid(T)==typeid(char) ||
       typeid(T)==typeid(unsigned char) ||
       typeid(T)==typeid(int) ||
       typeid(T)==typeid(unsigned int) ||
       typeid(T)==typeid(short) ||
       typeid(T)==typeid(unsigned short) ||
       typeid(T)==typeid(long) ||
       typeid(T)==typeid(unsigned long) ||
       typeid(T)==typeid(float) ||
       typeid(T)==typeid(double)) memMoveInit=1;
  }
  memMove=(memMoveInit==1);
  flexiMem=true;
  p=pstop=NULL;
  M=N=nd=d0=d1=d2=0;
  pp=NULL;
  sparse=NULL;
  p_device=NULL;
  d=&d0;
  mtype=full;
}


//***** constructors

//! standard constructor -- this becomes an empty array
template<class T> MT::Array<T>::Array(){ init(); }

//! copy constructor
template<class T> MT::Array<T>::Array(const MT::Array<T>& a){ init(); operator=(a); }

//! constructor with resize
template<class T> MT::Array<T>::Array(uint i){ init(); resize(i); }

//! constructor with resize
template<class T> MT::Array<T>::Array(uint i,uint j){ init(); resize(i,j); }

//! constructor with resize
template<class T> MT::Array<T>::Array(uint i,uint j,uint k){ init(); resize(i,j,k); }

//! constructor with resize
template<class T> MT::Array<T>::Array(uint i,uint j,uint k,uint l){ init(); resize(i,j,k,l); }

//! this becomes a reference on a subdimension of \c a
template<class T> MT::Array<T>::Array(const MT::Array<T>& a,uint i){ init(); referToSubDim(a,i); }

//! this becomes a reference on a subdimension of \c a
template<class T> MT::Array<T>::Array(const MT::Array<T>& a,uint i,uint j){ init(); referToSubDim(a,i,j); }

//! this becomes a reference on the C-array \c p
template<class T> MT::Array<T>::Array(const T* p,uint size){ init(); referTo(p,size); }

template<class T> MT::Array<T>::~Array(){
  freeMEM();
}
  

//***** resize

//! frees all memory; this becomes an empty array
template<class T> void MT::Array<T>::clear(){ freeMEM(); }

//! resize 1D array, discard the previous contents
template<class T> MT::Array<T>& MT::Array<T>::resize(uint D0){ nd=1; d0=D0; resetD(); resizeMEM(d0,false); return *this; }

//! resize but copy the previous contents
template<class T> MT::Array<T>& MT::Array<T>::resizeCopy(uint D0){ nd=1; d0=D0; resetD(); resizeMEM(d0,true); return *this; }

//! reshape the dimensionality (e.g. from 2D to 1D); throw an error if this actually requires to resize the memory
template<class T> MT::Array<T>& MT::Array<T>::reshape(uint D0){ CHECK(N==D0,"reshape must preserve total memory size"); nd=1; d0=D0; d1=d2=0; resetD(); return *this; }

//! same for 2D ...
template<class T> MT::Array<T>& MT::Array<T>::resize(uint D0,uint D1){ nd=2; d0=D0; d1=D1; resetD(); resizeMEM(d0*d1,false); return *this; }

//! ...
template<class T> MT::Array<T>& MT::Array<T>::resizeCopy(uint D0,uint D1){ nd=2; d0=D0; d1=D1; resetD(); resizeMEM(d0*d1,true); return *this; }

//! ...
template<class T> MT::Array<T>& MT::Array<T>::reshape(uint D0,uint D1){ CHECK(N==D0*D1,"reshape must preserve total memory size"); nd=2; d0=D0; d1=D1; d2=0; resetD(); return *this; }

//! same for 3D ...
template<class T> MT::Array<T>& MT::Array<T>::resize(uint D0,uint D1,uint D2){ nd=3; d0=D0; d1=D1; d2=D2; resetD(); resizeMEM(d0*d1*d2,false); return *this; }

//! ...
template<class T> MT::Array<T>& MT::Array<T>::resizeCopy(uint D0,uint D1,uint D2){ nd=3; d0=D0; d1=D1; d2=D2; resetD(); resizeMEM(d0*d1*d2,true); return *this; }

//! ...
template<class T> MT::Array<T>& MT::Array<T>::reshape(uint D0,uint D1,uint D2){ CHECK(N==D0*D1*D2,"reshape must preserve total memory size"); nd=3; d0=D0; d1=D1; d2=D2; resetD(); return *this; }

//! same for 4D ...
template<class T> MT::Array<T>& MT::Array<T>::resize(uint D0,uint D1,uint D2,uint D3){ uint dim[4]; dim[0]=D0;dim[1]=D1;dim[2]=D2;dim[3]=D3; resize(4,dim); return *this; }

//! resize to multi-dimensional tensor
template<class T> MT::Array<T>& MT::Array<T>::resize(uint ND,uint *dim){
  nd=ND; d0=d1=d2=0; resetD();
  uint j;
  for(j=0;j<nd && j<3;j++){ (&d0)[j]=dim[j]; }
  if(nd>3){ d=new uint[nd];  memmove(d,dim,nd*sizeof(uint)); }
  uint64 S;
  for(S=1,j=0;j<nd;j++) S*=dim[j];
  if(S>=(1ull<<32)) HALT("Array #elements "<<(S>>30)<<"G is >= 2^32");
  resizeMEM(S,false);
  return *this;
}

//! resize to multi-dimensional tensor
template<class T> MT::Array<T>& MT::Array<T>::resizeCopy(uint ND,uint *dim){
  nd=ND; d0=d1=d2=0; resetD();
  uint j;
  for(j=0;j<nd && j<3;j++){ (&d0)[j]=dim[j]; }
  if(nd>3){ d=new uint[nd];  memmove(d,dim,nd*sizeof(uint)); }
  uint64 S;
  for(S=1,j=0;j<nd;j++) S*=dim[j];
  if(S>=(1ull<<32)) HALT("Array #elements "<<(S>>30)<<"G is >= 2^32");
  resizeMEM(S,true);
  return *this;
}

//! resize to multi-dimensional tensor
template<class T> MT::Array<T>& MT::Array<T>::reshape(uint ND,uint *dim){
  nd=ND; d0=d1=d2=0; resetD();
  uint j,S;
  for(j=0;j<nd && j<3;j++){ (&d0)[j]=dim[j]; }
  if(nd>3){ d=new uint[nd];  memmove(d,dim,nd*sizeof(uint)); }
  for(S=1,j=0;j<nd;j++) S*=dim[j];
  CHECK(N==S,"reshape must preserve total memory size"); 
  return *this;
}


//! resize to multi-dimensional tensor
template<class T> MT::Array<T>& MT::Array<T>::resize(const Array<uint> &newD){ resize(newD.N,newD.p); return *this; }

//! resize to multi-dimensional tensor
template<class T> MT::Array<T>& MT::Array<T>::resizeCopy(const Array<uint> &newD){ resizeCopy(newD.N,newD.p); return *this; }

//! resize to multi-dimensional tensor
template<class T> MT::Array<T>& MT::Array<T>::reshape(const Array<uint> &newD){ reshape(newD.N,newD.p); return *this; }


template<class T> MT::Array<T>& MT::Array<T>::resizeAs(const MT::Array<T>& a){
  CHECK(this!=&a,"never do this!!!");
  nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2;
  resetD();
  if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
  resizeMEM(a.N,false);
  return *this;
}

//! make it the same size as \c a and copy previous content
template<class T> MT::Array<T>& MT::Array<T>::resizeCopyAs(const MT::Array<T>& a){
  CHECK(this!=&a,"never do this!!!");
  nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2; resetD();
  if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
  resizeMEM(a.N,true);
  return *this;
}

template<class T> MT::Array<T>& MT::Array<T>::reshapeAs(const MT::Array<T>& a){
  CHECK(this!=&a,"never do this!!!");
  CHECK(N==a.N,"reshape must preserve total memory size");
  nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2; resetD();
  if(nd>3){ d=new uint[nd];  memmove(d,a.d,nd*sizeof(uint)); }
  return *this;
}

//! return the size of memory allocated in bytes
template<class T> uint MT::Array<T>::getMemsize() const{ return M*sizeof(T); }

//! multi-dimensional (tensor) reference
template<class T> T& MT::Array<T>::operator()(const Array<uint> &I) const{
  CHECK(I.N == nd, "wrong dimensions");
  uint i,j;
  i=0;
  for(j=0;j<nd;j++) i = i*dim(j) + I(j);
  return p[i];
}

//! I becomes the index tuple for the absolute index i
template<class T> void MT::Array<T>::getIndexTuple(Array<uint> &I,uint i) const{
  uint j;
  CHECK(i<N,"out of range");
  I.resize(nd);
  I.setZero();
  for(j=nd;j--;){
    I.p[j] = i % d[j];
    i -= I.p[j];
    i /= d[j];
  }
}

//! return the k-th dimensionality
template<class T> uint MT::Array<T>::dim(uint k) const{
  CHECK(k<nd,"dimensionality range check error: " <<k <<"!<" <<nd);
  if(!d && k<3) return (&d0)[k]; else return d[k];
}


//***** sparse arrays

//! return fraction of non-zeros in the array
template<class T> real MT::Array<T>::sparsity(){
  uint i,m=0;
  for(i=0;i<N;i++) if(elem(i)) m++;
  return ((real)m)/N;
}

//! make sparse: create the \ref sparse index
template<class T> void MT::Array<T>::makeSparse(){
  CHECK(!sparse,"only once yet");
  uint n=0;
  if(nd==1){
    uint i;
    sparse=new Array<uint> [2];
    sparse[1].resize(d0); sparse[1]=-1;
    for(i=0;i<d0;i++) if(p[i]){
      sparse[0].append(i); //list of entries (maps n->i)
      sparse[1](i)=n;      //index list to entries (maps i->n)
      permute(i,n);
      n++;
    }
    N=n; resizeMEM(n,true);
    return;
  }
  if(nd==2){
    uint i,j;
    Array<uint> pair(2);
    sparse=new Array<uint> [1+d1+d0];
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) if(p[i*d1+j]){
      pair(0)=i; pair(1)=j; sparse[0].append(pair);   sparse[0].reshape(n+1,2);
      permute(i*d1+j,n);
      //register entry in columns an row indices
      pair(0)=i; pair(1)=n; sparse[1+j]   .append(pair); sparse[1+j]   .reshape(sparse[1+j]   .N/2,2);
      pair(0)=j; pair(1)=n; sparse[1+d1+i].append(pair); sparse[1+d1+j].reshape(sparse[1+d1+j].N/2,2);
      n++;
    }
    N=n; resizeMEM(n,true);
    return;
  }
}


//***** internal memory routines (probably not for external use)

//! allocate memory (maybe using \ref flexiMem)
template<class T> void MT::Array<T>::resizeMEM(uint n,bool copy){
  if(n==N) return;
  CHECK(!reference,"real resize of subarray is not allowed! (only a resize without changing memory size)");
  uint i;
  T *pold=p;
  uint Mold=M;
  //determine a new M (number of allocated items)
  if(!flexiMem){
    M=n;
  }else{
    if(n>0 && Mold==0){
      M=n;      //first time: exact allocation
    }else if(n>M || 10+2*n<M/2){
      M=10+2*n; //big resize: allocate with some extra space
    }else{
      M=Mold;   //small resize: don't really resize memory
    }
  }
  //if M changed, allocate the memory
  if(M!=Mold){
    uint64 memoryNew = ((uint64)M)*sizeT;
    memoryTotal -=  ((uint64)Mold)*sizeT;
    memoryTotal += memoryNew;
    if(M){
      if(memoryTotal>memoryBound){
        if(memoryStrict){
	        MT_MSG("allocating "<<(memoryNew>>20)<<"MB (total="<<(memoryTotal>>20)<<"M, bound="<<(memoryBound>>20)<<"M)...");
          //first undo changes...
          memoryTotal += ((uint64)Mold)*sizeT;
          memoryTotal -= memoryNew;
          p=pold; M=Mold;
          //...then throw an error...
          HALT("...throwing exception");
        }else{
          if(memoryNew>>20 || memoryTotal-memoryNew<=memoryBound){
	          MT_MSG("allocating "<<(memoryNew>>20)<<"MB (total="<<(memoryTotal>>20)<<"M, bound="<<(memoryBound>>20)<<"M)...");
	        }
          //MT_MSG("...are you sure?");
          //MT::wait();
        }
      }
      p=new T [M];      //p=(T*)malloc(M*sizeT);
      if(!p) HALT("MT::Array failed memory allocation of "<<M*sizeT<<"bytes");
      if(copy && !memMove) for(i=N<n?N:n;i--;) p[i]=pold[i];
      if(copy && memMove) memmove(p,pold,sizeT*(N<n?N:n));
    }else{
      p=0;
    }
    CHECK((pold && Mold) || (!pold && !Mold),"");
    if(Mold) delete[] pold;      //if(Mold) free(pold);
  }
  N=n;
  pstop=p+N;
  if(pp){ delete[] pp; pp=0; }
}

//! free all memory and reset all pointers and sizes
template<class T> void MT::Array<T>::freeMEM(){
  memoryTotal -= M*sizeT;
  if(M) delete[] p;
  //if(M) free(p);
  if(pp) delete[] pp;
  if(sparse) delete[] sparse;
  if(d && d!=&d0) delete[] d;
  p=pstop=NULL;
  M=N=nd=d0=d1=d2=0;
  pp=NULL;
  sparse=NULL;
  d=&d0;
}

//! reset the dimensionality pointer d to point to &d0
template<class T> void MT::Array<T>::resetD(){
  if(d && d!=&d0) delete[] d;
  d=&d0;
}



//***** append, insert & remove

//! append an (uninitialized) element to the array and return its reference -- the array becomes 1D!
template<class T> T& MT::Array<T>::append(){
  if(nd==2 && d1==1)
    resizeCopy(d0+1,d1);
  else
    resizeCopy(N+1);
  return p[N-1];
}

//! append an element to the array -- the array becomes 1D!
template<class T> void MT::Array<T>::append(const T& x){ append()=x; }

//! append another array to the array (by copying it) -- the array might become 1D!
template<class T> void MT::Array<T>::append(const MT::Array<T>& x){
  uint oldN=N,i;
  if(nd==2 && x.nd==1 && d1==x.N) 
    resizeCopy(d0+1,d1);
  else if(nd==2 && x.nd==2 && d1==x.d1)
    resizeCopy(d0+x.d0,d1);
  else
    resizeCopy(N+x.N);
  if(!memMove){
    for(i=0;i<x.N;i++) p[oldN+i]=x.p[i];
  }else memmove(p+oldN,x.p,sizeT*x.N);
}

//! append a C array to the array (by copying it) -- the array might become 1D!
template<class T> void MT::Array<T>::append(const T *q,uint n){
  uint oldN=N,i;
  if(nd==2 && d1==n) 
    resizeCopy(d0+1,d1);
  else
    resizeCopy(N+n);
  if(!memMove){
    for(i=0;i<n;i++) p[n+i]=q[i];
  }else memmove(p+oldN,q,sizeT*n);
}

//! append an element to the array if it is not included yet -- the array becomes 1D! [TL]
template<class T> void MT::Array<T>::setAppend(const T& x){ 
  if (findValue(x) < 0)
    append(x);
}

//! append elements of another array to the array which are not included yet -- the array might become 1D! [TL]
template<class T> void MT::Array<T>::setAppend(const MT::Array<T>& x){
  uint i;
  FOR1D(x, i) {
    setAppend(x(i));
  }
}



//! the array contains `copies' copies of the old one
template<class T> void MT::Array<T>::replicate(uint copies){
  if(copies<2) return;
  uint i,oldN=N;
  resizeCopy(copies*N);
  if(!memMove){
    NIY;
  }else{
    for(i=0;i<copies;i++) memmove(p+i*oldN,p,sizeT*oldN);
  }
}

//! inserts x at the position i -- the array becomes 1D! [only with memMove!]
template<class T> void MT::Array<T>::insert(uint i,const T& x){
  CHECK(memMove,"only with memMove");
  uint Nold=N;
  resizeCopy(Nold+1);
  if(i<Nold) memmove(p+i+1,p+i,sizeT*(Nold-i));
  p[i]=x;
}

//! remove (delete) a subsequence of the array -- the array becomes 1D!  [only with memMove!]
template<class T> void MT::Array<T>::remove(uint i,uint n){
  if(i==N-n){ resizeCopy(N-n); return; }
  CHECK(memMove,"only with memMove");
  if(N>i+n) memmove(p+i,p+i+n,sizeT*(N-i-n));
  resizeCopy(N-n);
}

//! remove some element by permuting the last element in its place -- the array becomes 1D!
template<class T> void MT::Array<T>::removePerm(uint i){
  p[i]=p[N-1];
  resizeCopy(N-1);
}

//! remove (delete) a subsequence of the array -- the array becomes 1D!  [only with memMove!]
template<class T> void MT::Array<T>::removeValue(const T& x){
  CHECK(memMove,"only with memMove");
  uint i;
  for(i=0;i<N;i++) if(p[i]==x) break;
  CHECK(i<N,"value to remove not found");
  remove(i,1);
}

//! remove (delete) a subsequence of the array -- the array becomes 1D!  [only with memMove!]
template<class T> void MT::Array<T>::removeAllValues(const T& x){
  CHECK(memMove,"only with memMove");
  uint i;
  for(i=0;i<N;i++) if(p[i]==x){ remove(i,1); i--; }
}
  
//! remove (delete) a subsequence of the array -- the array becomes 1D!  [only with memMove!]
//! Returns true if value was found and deleted.
//! Returns false if value was not found.
template<class T> bool MT::Array<T>::removeValueSafe(const T& x){
  CHECK(memMove,"only with memMove");
  uint i;
  for(i=0;i<N;i++) if(p[i]==x) break;
  if (i >= N)
    return false;
  else {
    remove(i,1);
    return true;
  }
}

//! replace n elements at pos i by the sequence x -- the array becomes 1D!  [only with memMove!]
template<class T> void MT::Array<T>::replace(uint i,uint n,const MT::Array<T>& x){
  CHECK(memMove,"only with memMove");
  uint Nold=N;
  if(n==x.N){
    memmove(p+i,x.p,sizeT*(x.N));
  }else if(n>x.N){
    memmove(p+i+x.N,p+i+n,sizeT*(Nold-i-n));
    if(i+n<Nold) memmove(p+i,x.p,sizeT*(x.N));
    resizeCopy(Nold-n+x.N);
  }else{
    resizeCopy(Nold+x.N-n);
    if(i+n<Nold) memmove(p+i+x.N,p+i+n,sizeT*(Nold-i-n));
    memmove(p+i,x.p,sizeT*(x.N));
  }
}

//! deletes the i-th row [must be 2D]
template<class T> void MT::Array<T>::delRow(uint i){
  CHECK(memMove,"only with memMove");
  CHECK(nd==2,"only for matricies");
  CHECK(i<d0,"range check error");
  uint n=d1;
  if(i+1<d0) memmove(p+i*n,p+(i+1)*n,sizeT*(d0-i-1)*n);
  resizeCopy(d0-1,n);
}

//! inserts k rows at the i-th row [must be 2D]
template<class T> void MT::Array<T>::insRows(uint i,uint k){
  CHECK(memMove,"only with memMove");
  CHECK(nd==2,"only for matricies");
  CHECK(i<=d0,"range check error");
  uint n=d0;
  resizeCopy(d0+k,d1);
  memmove(p+(i+k)*d1,p+i*d1,sizeT*d1*(n-i));
  memset (p+ i   *d1,0     ,sizeT*d1*k);
}

//! deletes k columns starting from the i-th (i==d1 -> deletes the last k columns)
template<class T> void MT::Array<T>::delColumns(uint i,uint k){
  CHECK(memMove,"only with memMove");
  CHECK(k>0,"");
  CHECK(nd==2,"only for matricies");
  if(i==d1) i=d1-k; 
  CHECK(i+k<=d1,"range check error");
  uint n=d1;
  for(uint j=0;j<d0;j++){
    memmove(p+j*(n-k)  ,p+j*n      ,sizeT*i);
    memmove(p+j*(n-k)+i,p+j*n+(i+k),sizeT*(n-i-k));
  }
  resizeCopy(d0,n-k);
}

//! inserts k columns at the i-th column [must be 2D]
template<class T> void MT::Array<T>::insColumns(uint i,uint k){
  CHECK(memMove,"only with memMove");
  CHECK(k>0,"");
  CHECK(nd==2,"only for matricies");
  CHECK(i<=d1,"range check error");
  uint n=d1;
  resizeCopy(d0,n+k);
  for(uint j=d0;j--;){
    memmove(p+j*d1+(i+k),p+j*n+i,sizeT*(n-i));
    memset (p+j*d1+i    ,0      ,sizeT*k);
    memmove(p+j*d1      ,p+j*n  ,sizeT*i);
  }
}

//! changes the range of one dimension (generalization of ins/delColumn to arbitrary tensors)
template<class T> void MT::Array<T>::resizeDim(uint k,uint dk){
  if(dim(k)==dk) return;
  uint i;
  uint L=1,R=1,subR=1;
  uintA oldDim(nd);
  for(i=0;i<nd;i++){
    oldDim(i)=dim(i);
    if(i<k) L*=dim(i);
    if(i>=k) R*=dim(i);
    if(i>k) subR*=dim(i);
  }
  //cout <<oldDim <<" L=" <<L <<" R=" <<R <<" subR=" <<subR <<endl;
  reshape(L,R);
  if(dk>oldDim(k)){ //need to increase dim
    i=dk-oldDim(k);
    insColumns(R,i*subR);
    oldDim(k)+=i;
  }else{
    i=oldDim(k)-dk;
    delColumns(R,i*subR);
    oldDim(k)-=i;
  }
  reshape(oldDim);
}

//***** access operations

//! the \c ith element
template<class T> T& MT::Array<T>::elem(uint i) const{
  CHECK(i<N,"range error ("<<i<<">="<<N<<")");
  return p[i];
}

/* scalar reference (legal iff N==1) */
/*operator T&() const{ 
  CHECK(N==1,"scalar reference ("<<N<<"!=1)");
  return *p;
  }*/

//! scalar reference
template<class T> T& MT::Array<T>::scalar() const{
  CHECK(nd==0 && N==1,"scalar range error (N="<<N<<")");
  return *p;
}

template<class T> T& MT::Array<T>::last() const{
  CHECK(N,"can't take last from empty");
  return p[N-1];
}

//! 1D reference
template<class T> T& MT::Array<T>::operator()(uint i) const{ 
  CHECK(nd==1 && i<d0,
	"1D range error ("<<nd<<"=1,"<<i<<"<"<<d0<<")");
  return p[i];
}

//! 2D reference
template<class T> T& MT::Array<T>::operator()(uint i,uint j) const{
  CHECK(nd==2 && i<d0 && j<d1 && !sparse,
	"2D range error ("<<nd<<"=2,"<<i<<"<"<<d0<<","<<j<<"<"<<d1<<")");
  return p[i*d1+j];
}

//! 3D reference
template<class T> T& MT::Array<T>::operator()(uint i,uint j,uint k) const{
  CHECK(nd==3 && i<d0 && j<d1 && k<d2 && !sparse,
	"3D range error ("<<nd<<"=3,"<<i<<"<"<<d0<<","<<j<<"<"<<d1<<","<<k<<"<"<<d2<<")");
  return p[(i*d1+j)*d2+k];
}

//! 4D reference
template<class T> T& MT::Array<T>::operator()(uint i,uint j,uint k,uint l) const{
  CHECK(nd==4 && i<d0 && j<d1 && k<d2 && l<d[3] && !sparse,
        "4D range error ("<<nd<<"=3,"<<i<<"<"<<d0<<","<<j<<"<"<<d1<<","<<k<<"<"<<d2<<","<<l<<"<"<<d[3]<<")");
  return p[((i*d1+j)*d2+k)*d[3]+l];
}

//! get a subarray; use in conjuction with operator()() to get a reference
template<class T> MT::Array<T> MT::Array<T>::operator[](uint i) const{ return Array(*this,i); }

//! get a subarray; use in conjuction with operator()() to get a reference
template<class T> MT::Array<T> MT::Array<T>::subDim(uint i,uint j) const{ return Array(*this,i,j); }

//! convert a subarray into a reference (e.g. a[3]()+=.123)
template<class T> MT::Array<T>& MT::Array<T>::operator()(){ return (*this); }

//! reference to the max entry
template<class T> T& MT::Array<T>::max() const{ CHECK(N,""); uint i,m=0; for(i=1;i<N;i++) if(p[i]>p[m]) m=i; return p[m]; }

//! reference to the min entry
template<class T> T& MT::Array<T>::min() const{ CHECK(N,""); uint i,m=0; for(i=1;i<N;i++) if(p[i]<p[m]) m=i; return p[m]; }

//! gets the min and max
template<class T> void MT::Array<T>::minmax(T& minVal, T& maxVal) const{
  CHECK(N,""); 
  uint i;
  minVal=maxVal=p[0];
  for(i=1;i<N;i++){
    if(p[i]<minVal) minVal=p[i];
    else if(p[i]>maxVal) maxVal=p[i];
  }
}

//! get absolute maximum (using fabs)
template<class T> T MT::Array<T>::absMax() const{
  CHECK(N,"");
  uint i;
  T t((T)::fabs(p[0]));
  for(i=1;i<N;i++) if(fabs(p[i])>t) t=(T)::fabs(p[i]);
  return t;
}


/*
//! also returns the index (argmin) \c ind
template<class T> T minA(const MT::Array<T>& v, uint & ind, uint start, uint end){
  CHECK(v.N>0,"");
  CHECK(v.N>start,"");	
  CHECK(v.N>=end,"");
  CHECK(end>=start,"");
  T t(v(start));
  ind=start;
  if (end==0)	end=v.N;
  for(uint i=start; i<end; ++i) if(v.p[i]<t){
    t  =v.p[i];
    ind=i;
  }
  return t;
}

//! also returns the index (argmax) \c ind
template<class T> T maxA(const MT::Array<T>& v, uint & ind, uint start, uint end){
  CHECK(v.N>0,"");
  CHECK(v.N>start,"");
  CHECK(v.N>=end,"");
  CHECK(end>=start,"");
  T t(v(start));
  ind=start;
  if(end==0){
    end=v.N;
  }
  for(uint i=start; i<end; ++i) if(v.p[i]>t){
    t  =v.p[i];
    ind=long(i);
  }
  return t;
}
*/


/*!\brief the index of the maxium; precondition: the comparision operator
  > exists for type T */ 
template<class T> uint MT::Array<T>::maxIndex() const{ uint i,m=0; for(i=0;i<N;i++) if(p[i]>p[m]) m=i; return m; }

/*!\brief the index of the maxium; precondition: the comparision operator
  > exists for type T */ 
template<class T> void MT::Array<T>::maxIndex(uint& i,uint& j) const{ CHECK(nd==2,"needs 2D array"); j=maxIndex(); i=j/d1; j=j%d1; }

/*!\brief the index of the maxium; precondition: the comparision operator
  > exists for type T */ 
template<class T> void MT::Array<T>::maxIndex(uint& i,uint& j,uint& k) const{ CHECK(nd==3,"needs 3D array"); k=maxIndex(); i=k/(d1*d2); k=k%(d1*d2); j=k/d2; k=k%d2; }

template<class T> void MT::Array<T>::maxIndeces(uint& m1,uint& m2) const{
  uint i;
  if(p[0]>p[1]){ m1=0; m2=1; }else{ m1=1; m2=0; }
  for(i=2;i<N;i++){
    if(p[i]>p[m2]){ //greater than m2
      if(p[i]>p[m1]){ //greater than m1
        m2=m1;
        m1=i;
      }else{
        m2=i;
      }
    }
  }
}

//! ...
template<class T> uint MT::Array<T>::minIndex() const{ uint i,m=0; for(i=0;i<N;i++) if(p[i]<p[m]) m=i; return m; }

template<class T> int MT::Array<T>::findValue(const T& x) const{ uint i; for(i=0;i<N;i++) if(p[i]==x) return i; return -1; }

template<class T> void MT::Array<T>::findValues(MT::Array<uint>& indices, const T& x) const {
  indices.clear();
  uint i;
  for(i=0;i<N;i++)
    if(p[i]==x) indices.append(i);
}

// TL 17.07.08
/*!\brief whether at least one object is contained more than once
   */
template<class T> bool MT::Array<T>::containsDoubles() const{ 
  if (N < 2)
    return false;
  uint i, j;
  for(i=0;i<N-1;i++) {
      for (j=i+1; j<N; j++) {
          if (p[i]==p[j])
              return true;
      }
  }
  return false;
}


/*!\brief a sub array of a 1D Array (corresponds to matlab [i:I]); when
  the upper limit I is -1, it is replaced by the max limit (like
  [i:]) */
template<class T> MT::Array<T> MT::Array<T>::sub(uint i,int I) const{
  CHECK(nd==1,"1D range error ");
  MT::Array<T> x;
  if(I==-1) I=d0-1;
  if((uint)I+1==i){ x.resize(0); return x; }
  CHECK(i<=(uint)I,"lower limit higher than upper!");
  x.resize(I-i+1);
  uint k;
  for(k=i;k<=(uint)I;k++) x(k-i)=operator()(k);
  return x;
}

/*!\brief copies a sub array of a 2D Array (corresponds to matlab [i:I,j:J]);
  when the upper limits I or J are -1, they are replaced by the
  max limit (like [i:,j:]) */
template<class T> MT::Array<T> MT::Array<T>::sub(uint i,int I,uint j,int J) const{
  CHECK(nd==2,"2D range error ");
  MT::Array<T> x;
  if(I==-1) I=d0-1;
  if(J==-1) J=d1-1;
  CHECK(i<=(uint)I && j<=(uint)J,"lower limit must be higher than upper!");
  x.resize(I-i+1,J-j+1);
  uint k,l;
  for(k=i;k<=(uint)I;k++) for(l=j;l<=(uint)J;l++) x(k-i,l-j)=operator()(k,l);
  return x;
}

/*!\brief copies a sub array of a 2D Array (corresponds to matlab [i:I,j:J]);
  when the upper limits I or J are -1, they are replaced by the
  max limit (like [i:,j:]) */
template<class T> MT::Array<T> MT::Array<T>::sub(uint i,int I,uint j,int J,uint k,int K) const{
  CHECK(nd==3,"3D range error ");
  MT::Array<T> x;
  if(I==-1) I=d0-1;
  if(J==-1) J=d1-1;
  if(K==-1) K=d2-1;
  CHECK(i<=(uint)I && j<=(uint)J && k<=(uint)K,"lower limit must be higher than upper!");
  x.resize(I-i+1,J-j+1,K-k+1);
  uint ii,jj,kk;
  for(ii=i;ii<=(uint)I;ii++) for(jj=j;jj<=(uint)J;jj++)  for(kk=k;kk<=(uint)K;kk++) x(ii-i,jj-j,kk-k)=operator()(ii,jj,kk);
  return x;
}

/*!\brief copies a selection of columns from a 2D array, the first index (rows)
  runs from i to I (as explained above) while the second index runs
  over the columns explicitly referred to by cols. (col doesn't have
  to be ordered or could also contain some columns multiply) */
template<class T> MT::Array<T> MT::Array<T>::sub(uint i,int I,Array<uint> cols) const{
  CHECK(nd==2,"2D range error ");
  MT::Array<T> x;
  if(I==-1) I=d0-1;
  CHECK(i<=(uint)I,"lower limit higher than upper!");
  x.resize(I-i+1,cols.N);
  uint k,l;
  for(k=i;k<=(uint)I;k++) for(l=0;l<cols.N;l++) x(k-i,l)=operator()(k,cols(l));
  return x;
}


//***** C-array interfacing

//! allocates, sets and returns the \c Array::pp pointer (of type \c T**)
template<class T> T** MT::Array<T>::getCarray() const{
  CHECK(nd>=2,"only 2D or higher-D arrays gives C-array of type T**");
  if(pp) return pp;
  ((MT::Array<T>*)this)->pp=new T* [d0];
  uint skip;
  if(nd==2) skip=d1; else skip=N/d0;
  for(uint i=0;i<d0;i++) pp[i]=p+i*skip;
  return pp;
}

//! makes this array a reference to the C buffer
template<class T> void MT::Array<T>::referTo(const T *buffer,uint n){
  freeMEM();
  reference=true;
  nd=1; d0=n; d1=d2=0; N=n;
  p=(T*)buffer;
  pstop=p+N;
}

#if 0
/*!\brief returns an ordinary 2-dimensional C-pointer to the Array content.
  Requires the Array<T*> as buffer. */
template<class T> T** MT::Array<T>::getPointers(Array<T*>& array2d) const{
  CHECK(nd==2,"only 2D array gives C-array of type T**");
  array2d.resize(d0);
  for(uint i=0;i<d0;i++) array2d(i)=p+i*d1;
  return array2d.p;
}

//! returns an ordinary 3-dimensional C-pointer-array
template<class T> T*** MT::Array<T>::getPointers(Array<T**>& array3d,Array<T*>& array2d) const{
  CHECK(nd==3,"only 3D array gives C-array of type T*** ");
  array2d.resize(d0,d1);
  for(uint i=0;i<d0;i++){
    for(uint j=0;j<d1;j++) array2d(i,j)=&operator()(i,j,0);
    array3d(i)=&array2d(i,0);
  }
  return array3d.p;
}
#endif


//***** assignments

//! set all elements to value \c v 
template<class T> MT::Array<T>& MT::Array<T>::operator=(const T& v){
  uint i;
  //if(memMove && typeid(T)==typeid(T)) memset(p,*((int*)&v),N); else
  for(i=0;i<N;i++) p[i]=v;
  return *this;
}

//! copy operator
template<class T> MT::Array<T>& MT::Array<T>::operator=(const MT::Array<T>& a){
  //if(a.temp){ takeOver(*((MT::Array<T>*)&a)); return *this; }
  resizeAs(a);
  uint i;
  if(memMove) memmove(p,a.p,sizeT*N);
  else for(i=0;i<N;i++) p[i]=a.p[i];
  return *this;
}

//! reads it from the C string (by using std::istringstream)
template<class T> void MT::Array<T>::setText(const char* str){
  std::istringstream s(str);
  read(s);
}

/*!\brief same as memset(p,0,sizeT*N); precondition: memMove is
  true! */
template<class T> void MT::Array<T>::setZero(byte zero){
  CHECK(memMove,"can set array's memory to zero only if memMove option is true");
  memset(p,zero,sizeT*N);
}

//! set all entries to same value x [default: don't change dimension]
template<class T> void MT::Array<T>::setUni(const T& x,int d){
  if(d!=-1) resize(d);
  uint i;
  for(i=0;i<N;i++) elem(i)=x;
}

/*!\brief becomes the n-dim identity matrix [default:
  don't change dimension (already has to be squared matrix)] */
template<class T> void MT::Array<T>::setId(int d){
  CHECK(d!=-1 || (nd==2 && d0==d1),"need squared matrix to set to identity");
  if(d!=-1) resize(d,d);
  setZero();
  for(uint i=0;i<d0;i++) operator()(i,i)=(T)1;
}

template<class T> void MT::Array<T>::setDiag(const T& x,int d){
  CHECK(d!=-1 || nd==2,"need squared matrix to set to diagonal");
  if(d!=-1) resize(d,d);
  if(d==-1) d=MT::MIN(d0,d1);
  setZero();
  uint i;
  for(i=0;i<(uint)d;i++) operator()(i,i)=x;
}

//! sets x to be the diagonal matrix with diagonal v
template<class T> void MT::Array<T>::setDiag(const MT::Array<T>& v){
  CHECK(v.nd==1,"can only give diagonal of 1D array");
  resize(v.d0,v.d0);
  setZero();
  uint i;
  for(i=0;i<v.d0;i++) operator()(i,i)=v(i);
}


/*
//! constructs the block matrix X=[A,B ; C,D] where D is a scalar
template<class T> void blockMatrix(MT::Array<T>& X,const MT::Array<T>& A,const MT::Array<T>& B,const MT::Array<T>& C,const T& D){
  CHECK(A.nd==2 && B.nd==1 && C.nd==1,"");
  CHECK(A.d0==B.d0 && A.d1==C.d0,"");
  uint i,j,a=A.d0,b=A.d1;
  X.resize(A.d0+1,A.d1+1);
  for(i=0;i<a;i++) for(j=0;j<b;j++) X(i,j)=A(i,j);
  for(i=0;i<a;i++) X(i,b)=B(i);
  for(j=0;j<b;j++) X(a,j)=C(j);
  X(a,b)=D;
}
*/

//! constructs the block matrix X=[A,B ; C,D]
template<class T> void MT::Array<T>::setBlockMatrix(const MT::Array<T>& A,const MT::Array<T>& B,const MT::Array<T>& C,const MT::Array<T>& D){
  CHECK(A.nd==2 && B.nd==2 && C.nd==2 && D.nd==2,"");
  CHECK(A.d0==B.d0 && A.d1==C.d1 && B.d1==D.d1 && C.d0==D.d0,"");
  uint i,j,a=A.d0,b=A.d1;
  resize(A.d0+C.d0,A.d1+B.d1);
  for(i=0;i<A.d0;i++) for(j=0;j<A.d1;j++) operator()(i  ,j  )=A(i,j);
  for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) operator()(i  ,j+b)=B(i,j);
  for(i=0;i<C.d0;i++) for(j=0;j<C.d1;j++) operator()(i+a,j  )=C(i,j);
  for(i=0;i<D.d0;i++) for(j=0;j<D.d1;j++) operator()(i+a,j+b)=D(i,j);
}

//! constructs a vector x=[a,b]
template<class T> void MT::Array<T>::setBlockVector(const MT::Array<T>& a,const MT::Array<T>& b){
  CHECK(a.nd==1 && b.nd==1,"");
  uint i;
  resize(a.N+b.N);
  for(i=0;i<a.N;i++) operator()(i    )=a(i);
  for(i=0;i<b.N;i++) operator()(i+a.N)=b(i);
}

template<class T> void MT::Array<T>::setMatrixBlock(const MT::Array<T>& B,uint lo0,uint lo1){
  CHECK(nd==2 && B.nd==2 && lo0+B.d0<=d0 && lo1+B.d1<=d1,"");
  uint i,j;
  for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) operator()(lo0+i,lo1+j)=B(i,j);
}

template<class T> void MT::Array<T>::getMatrixBlock(MT::Array<T>& B,uint lo0,uint lo1) const{
  CHECK(nd==2 && B.nd==2 && lo0+B.d0<=d0 && lo1+B.d1<=d1,"");
  uint i,j;
  for(i=0;i<B.d0;i++) for(j=0;j<B.d1;j++) B(i,j)=operator()(lo0+i,lo1+j);
}

template<class T> void MT::Array<T>::setVectorBlock(const MT::Array<T>& B,uint lo){
  CHECK(nd==1 && B.nd==1 && lo+B.N<=N,"");
  uint i;
  for(i=0;i<B.N;i++) operator()(lo+i)=B(i);
}

template<class T> void MT::Array<T>::getVectorBlock(MT::Array<T>& B,uint lo) const{
  CHECK(nd==1 && B.nd==1 && lo+B.N<=N,"");
  uint i;
  for(i=0;i<B.N;i++) B(i)=operator()(lo+i);
}

//! sorted permutation of length \c n
template<class T> void MT::Array<T>::setStraightPerm(int n){
  if(n!=-1) resize(n);
  for(uint i=0;i<N;i++) elem(i)=(T)i;
}

//! reverse sorted permutation of lenth \c N
template<class T> void MT::Array<T>::setReversePerm(int n){
  if(n!=-1) resize(n);
  for(uint i=0;i<N;i++) elem(N-1-i)=(T)i;
}

//! permute all elements randomly
template<class T> void MT::Array<T>::setRandomPerm(int n){
  setStraightPerm(n);
  int j,r;
  for(j=N-1;j>=1;j--){
    r=rnd(j+1);
    permute(r,j);
  }
}

template<class T> void MT::Array<T>::setCarray(T *buffer,uint D0){
  resize(D0);
  uint i;
  if(memMove && typeid(T)==typeid(T))
    memmove(p,buffer,sizeT*d0);
  else for(i=0;i<d0;i++) operator()(i)=(T)buffer[i];
}

template<class T> void MT::Array<T>::setCarray(T **buffer,uint D0,uint D1){
  resize(D0,D1);
  uint i,j;
  for(i=0;i<d0;i++){
    if(memMove && typeid(T)==typeid(T))
      memmove(p+i*d1,buffer[i],sizeT*d1);
    else for(j=0;j<d1;j++) operator()(i,j)=(T)buffer[i][j];
  }
}

template<class T> void MT::Array<T>::copyInto(T *buffer) const{
  CHECK(nd==1,"can only copy 1D Array into 1D C-array");
  uint i;
  if(memMove && typeid(T)==typeid(T)) memmove(buffer,p,sizeT*d0);
  else for(i=0;i<d0;i++) buffer[i]=(T)operator()(i);
}

template<class T> void MT::Array<T>::copyInto2D(T **buffer) const{
  CHECK(nd==2,"can only copy 2D Array into 2D C-array");
  uint i,j;
  for(i=0;i<d0;i++){
    if(memMove && typeid(T)==typeid(T)) memmove(buffer[i],p+i*d1,sizeT*d1);
    else for(j=0;j<d1;j++) buffer[i][j]=(T)operator()(i,j);
  }
}

//! make this array a reference to the array \c a
template<class T> void MT::Array<T>::referTo(const MT::Array<T>& a){
  freeMEM();
  reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
  N=a.N; nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2;
  p=a.p;
  pstop=a.pstop;
}

//! make this array a subarray reference to \c a
template<class T> void MT::Array<T>::referToSubRange(const MT::Array<T>& a,uint i,int I){
  CHECK(a.nd<=3,"not implemented yet");
  freeMEM();
  resetD();
  reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
  if(I==-1) I=a.d0-1;
  CHECK(i<a.d0 && (uint)I<a.d0,"SubRange range error ("<<i<<"<"<<a.d0<<","<<I<<"<"<<a.d0<<")");
  if(a.nd==1){
    nd=1;  d0=I+1-i; d1=0; d2=0;  N=d0;
    p=a.p+i;
  }
  if(a.nd==2){
    nd=2;  d0=I+1-i; d1=a.d1; d2=0;  N=d0*d1;
    p=a.p+i*d1;
  }
  if(a.nd==3){
    nd=3;  d0=I+1-i; d1=a.d1; d2=a.d2;  N=d0*d1*d2;
    p=a.p+i*d1*d2;
  }
  pstop=p+N;
}

//! make this array a subarray reference to \c a
template<class T> void MT::Array<T>::referToSubDim(const MT::Array<T>& a,uint dim){
  CHECK(a.nd>1,"can't create subarray of array less than 2 dimensions");
  CHECK(dim<a.d0,"SubDim range error ("<<dim<<"<"<<a.d0<<")");
  freeMEM();
  reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
  if(a.nd==2){
    nd=1; d0=a.d1; d1=d2=0; N=d0;
  }
  if(a.nd==3){
    nd=2; d0=a.d1; d1=a.d2; d2=0; N=d0*d1;
  }
  if(a.nd>3){
    nd=a.nd-1; d0=a.d1; d1=a.d2; d2=a.d[3]; N=a.N/a.d0;
    resetD();
    if(nd>3){ d=new uint[nd];  memmove(d,a.d+1,nd*sizeof(uint)); }
  }
  p=a.p+dim*N;
  pstop=p+N;
}

//! make this array a subarray reference to \c a
template<class T> void MT::Array<T>::referToSubDim(const MT::Array<T>& a,uint i,uint j){
  CHECK(a.nd>2,"can't create subsubarray of array less than 3 dimensions");
  CHECK(i<a.d0 && j<a.d1,"SubDim range error ("<<i<<"<"<<a.d0<<","<<j<<"<"<<a.d1<<")");
  freeMEM();
  reference=true; memMove=a.memMove; flexiMem=a.flexiMem;
  if(a.nd==3){
    nd=1; d0=a.d2; d1=0; d2=0; N=d0;
    p=&a(i,j,0);
    pstop=p+N;
  }else{
    NIY;
  }
}

/*
//! redirects a reference (previously constructed) to a new sub dim of a
template<class T> void MT::Array<T>::redirect(const MT::Array<T>& a,uint i){
  CHECK(reference && a.nd==2,"can only do that hack with references");
  p=a.p+i*a.d1;
}
*/

/*! takes over the memory buffer from a; afterwards, this is a
  proper array with own memory and a is only a reference on the
  memory */
template<class T> void MT::Array<T>::takeOver(MT::Array<T>& a){
  freeMEM();
  memMove=a.memMove; flexiMem=a.flexiMem;
  N=a.N; nd=a.nd; d0=a.d0; d1=a.d1; d2=a.d2;
  p=a.p;
  pstop=a.pstop;
  M=a.M;
  a.reference=true;
  a.M=0;
}

/*!\brief return a `dim'-dimensional grid with `steps' intervals
  filling the range [lo,hi] in each dimension */
template<class T> void MT::Array<T>::setGrid(uint dim,T lo,T hi,uint steps){
  uint i,j;
  if(dim==1){
    resize(steps+1,1);
    for(i=0;i<d0;i++) operator()(i,0)=lo+(hi-lo)*i/steps;
    return;
  }
  if(dim==2){
    resize(steps+1,steps+1,2);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++){
      operator()(i,j,0)=lo+(hi-lo)*j/steps;
      operator()(i,j,1)=lo+(hi-lo)*i/steps;
    }
    reshape(d0*d1,2);
    return;
  }
  HALT("not implemented yet");
}



//----- sorting etc
template<class T> uint MT::Array<T>::findInSorted(const T& x,ElemCompare comp){
  if(!N) return 0;
  if(comp(x,elem(0))) return 0;
  if(!comp(x,elem(N-1))) return N;
  uint lo=0,up=N-1,mi;
  for(;;){
    if(lo+1==up) return up;
    mi=lo+(up-lo)/2;
    if(comp(x,elem(mi))) up=mi; else lo=mi;
  }
  HALT("you shouldn't be here!");
  return 0;
}
template<class T> uint MT::Array<T>::insertSorted(const T& x,ElemCompare comp){
  CHECK(memMove,"");
  if(!N){ append(x); return 1; }
  uint i=findInSorted(x,comp);
  insert(i,x);
  return i;
}
template<class T> void MT::Array<T>::removeValueInSorted(const T& x,ElemCompare comp){
  uint i=findInSorted(x,comp);
  CHECK(elem(i)==x,"value not found");
  remove(i);
}


//***** permutations

//! permute the elements \c i and \c j
template<class T> void MT::Array<T>::permute(uint i,uint j){ T x=p[i]; p[i]=p[j]; p[j]=x; }

//! permutes the entries according to the current this
template<class T> void MT::Array<T>::permute(const uintA& permutation){
  CHECK(permutation.N<=N,"array smaller than permutation ("<<N<<"<"<<permutation.N<<")");
  MT::Array<T> b=(*this);
  for(uint i=0;i<N;i++) elem(i)=b.elem(permutation(i));
}

//! permutes the entries of a according to the inverse of the current this
template<class T> void MT::Array<T>::permuteInv(const uintA& permutation){
  CHECK(permutation.N<=N,"array smaller than permutation ("<<N<<"<"<<permutation.N<<")");
  MT::Array<T> b=(*this);
  for(uint i=0;i<N;i++) elem(permutation(i))=b.elem(i);
}

//! randomly permute all entries of a
template<class T> void MT::Array<T>::permuteRandomly(){
  uintA perm;
  perm.setRandomPerm(N);
  permute(perm);
}

//! push all elements backward
template<class T> void MT::Array<T>::shiftPerm(int offset){
  static MT::Array<T> tmp;
  CHECK(memMove,"pushing only works with memMove enabled");
  uint m=offset>0?offset:-offset;
  tmp.resize(m);
  if(offset>0){
    memmove(tmp.p,p+N-m,sizeT*m);
    memmove(p+m,p,sizeT*(N-m));
    memmove(p,tmp.p,sizeT*m);
  }
  if(offset<0){
    memmove(tmp.p,p,sizeT*m);
    memmove(p,p+m,sizeT*(N-m));
    memmove(p+(N-m),tmp.p,sizeT*m);
  }
}


//***** I/O

//! prototype for operator<<
template<class T> void MT::Array<T>::write(std::ostream& os,const char *ELEMSEP,const char *LINESEP,bool dimTag,bool binary) const{
  CHECK(!binary || memMove,"binary write works only for memMoveable data");
  uint i,j,k;
  
  if(binary){
    writeDim(os);
    os <<std::endl;
    os.put(0);
    os.write((char*)p,sizeT*N);
    os.put(0);
    os <<std::endl;
  }else{
    if(dimTag){ writeDim(os); if(nd>=2) os <<'\n'; else os <<' '; }
    if(!IOraw) os <<'[';
    if(nd==0 && N==0){ os <<']'; return; }
    if(nd==0 && N==1){ os <<scalar()<<']'; return; }
    if(nd==1){
      //os <<' ';
      for(i=0;i<N;i++) os <<ELEMSEP  <<operator()(i);
    }
    if(nd==2) for(j=0;j<d0;j++){
      if(j) os <<LINESEP;
      for(i=0;i<d1;i++) os <<ELEMSEP <<operator()(j,i);
    }
    if(nd==3) for(k=0;k<d0;k++){
      if(k) os <<LINESEP;
      for(j=0;j<d1;j++){
        for(i=0;i<d2;i++) os <<ELEMSEP <<operator()(k,j,i);
        os <<LINESEP;
      }
    }
    if(nd>3){
      CHECK(d && d!=&d0,"");
      //Array<uint> I;
      for(i=0;i<N;i++){
        if(i && !(i%d[nd-1])) os <<LINESEP;
        if(nd>1 && !(i%(d[nd-2]*d[nd-1]))){
          /*getIndexTuple(I,i);
          os <<LINESEP <<'<' <<I(0);
          for(j=1;j<nd;j++) os <<' ' <<I(j);
          os <<':' <<i <<'>' <<LINESEP;
          */
          os <<LINESEP;
        }
        os <<ELEMSEP <<elem(i);
      }
    }
    if(!IOraw) os <<ELEMSEP <<']';
  }
}

//! prototype for operator>>
template<class T> void MT::Array<T>::read(std::istream& is){
  uint d,i;
  char c;
  T x;
  c=MT::peerNextChar(is);
  switch(c){
  case '<':
    readDim(is);
    c=MT::peerNextChar(is);
    if(c=='['){ //fast ascii read
      is >>"[";
      for(i=0;i<N;i++){
        if(is.fail()) HALT("could not read "<<i<<"-th element of an array");
        is >>p[i];
      }
      is >>"]";
      if(is.fail()) HALT ("could not read array end tag");
    }else if(c==0){//binary read
      c=is.get(); CHECK(c==0,"couldn't read newline before binary data block :-(");
      is.read((char*)p,sizeT*N);
      if(is.fail()) HALT("could not binary data");
      c=is.get(); CHECK(c==0,"couldn't read newline after binary data block :-(");
    }else{ //just directly read numbers
      for(i=0;i<N;i++){
        if(is.fail()) HALT("could not read "<<i<<"-th element of an array");
        is >>p[i];
      }
    }
    break;
  case '[': //slow read
    is >>"[";
  default:
    uint i=0;
    d=0;
    for(;;){
      MT::skip(is," \r\t");
      is.get(c);
      if(c==']' || !is.good()){ is.clear(); break; }
      if(c==';' || c=='\n'){ //set an array width
        if(!d) d=i; else CHECK(!(i%d),"Array::read: mis-structured array in row "<<i/d);
        continue;
      }
      if(c!=',') is.putback(c);
      is >>x;
      if(!is.good()){ is.clear(); break; }
      if(i>=N) resizeCopy(i+1000);
      elem(i)=x;
      i++;
    }
    resizeCopy(i);
    if(d){
      CHECK(!(N%d),"Array::read: mis-structured array in last row");
      reshape(N/d,d);
    }
    break;
  }
}

/*
template<class T> void MT::Array<T>::readOld(std::istream& is){
  uint d,i,j,k;
  char c;
  MT::skip(is);
  is.get(c);
  switch(c){
  case '[':
    is >>d;
    if(is.fail()) HALT ("could not read array tag");
    if(d==0){ is >>">"; return; }
    if(d==1){
      is >>":" >>i;
      if(is.fail()) HALT ("could not read array's dimensions");
      resize(i);
    }
    if(d==2){
      is >>":" >>i >>"," >>j;
      if(is.fail()) HALT ("could not read array's dimensions");
      resize(i,j);
    }
    if(d==3){
      is >>":" >>i >>"," >>j >>"," >>k;
      if(is.fail()) HALT ("could not read array's dimensions");
      resize(i,j,k);
    }
    is >>"]";
    if(is.fail()) HALT ("could not read array end tag");
    for(i=0;i<N;i++){
      if(is.fail()) HALT("could not read "<<i<<"-th element of an array");
      is >>p[i];
    }
    break;
  default:
    NIY;
  }
}
*/

template<class T> void MT::Array<T>::writeDim(std::ostream& os) const{
  uint i;
  os <<'<';
  if(nd) os <<dim(0); else os <<0;
  for(i=1;i<nd;i++) os <<' ' <<dim(i);
  os <<'>';
}

template<class T> void MT::Array<T>::readDim(std::istream& is){
  char c;
  uint ND,dim[10];
  is >>"<";
  for(ND=0;;ND++){
    is >>dim[ND];
    is.get(c);
    if(c=='>') break;
    CHECK(c==' ',"error in reading dimensionality");
  }
  resize(ND+1,dim);
}

//! IO-manip: arrays will be streamed as raw (without tags)
template<class T> const MT::Array<T>& MT::Array<T>::ioraw() const{ IOraw=true; return *this; }

//! IO-manip: arrays will be streamed non-raw (with tags)
//template<class T> const MT::Array<T>& MT::Array<T>::ionoraw() const{ IOraw=false; return *this; }

//! array must have correct size! simply streams in all elements sequentially
template<class T> void MT::Array<T>::readRaw(std::istream& is){
  uint i;
  for(i=0;i<N;i++){
    is >>p[i];
    if(is.fail()) HALT("could not read "<<i<<"-th element of an array");
  }
}

template<class T> void MT::Array<T>::writeRaw(std::ostream& os) const{
  bool tmp=IOraw;
  IOraw=true;
  write(os);
  IOraw=tmp; 
}

// TL 15.07.08
template<class T> void MT::Array<T>::writeWithIndex(std::ostream& os) const{
    uint i;
    FOR1D((*this), i) {
        os<<i<<" "<<elem(i)<<endl;
    }
}

//***** generic data files of real arrays

template<class T> void MT::Array<T>::writeTagged(std::ostream& os,const char* tag,bool binary) const{
  os <<tag <<' ';
  write(os," ","\n ",true,binary);
}

template<class T> bool MT::Array<T>::readTagged(std::istream& is,const char *tag){
  if(tag){
    String read_tag;
    read_tag.read(is," \t\n\r"," \t\n\r");
    if(!is.good() || read_tag.N()==0) return false;
    CHECK(read_tag==tag,"read `"<<read_tag<<"' instead of `"<<tag<<"' in arr file");
  };
  read(is);
  return true;
}



//! x = y^T
template<class T> void transpose(MT::Array<T>& x,const MT::Array<T>& y){
  CHECK(y.nd<=3,"can only transpose up to 3D arrays");
  if(y.nd==3){
    uint i,j,k,d0=y.d2,d1=y.d1,d2=y.d0;
    x.resize(d0,d1,d2);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) for(k=0;k<d2;k++)
      x(i,j,k) = y(k,j,i);
    //x.p[(i*d1+j)*d2+k]=y.p[(k*d1+j)*d0+i];
    return;
  }
  if(y.nd==2){
    uint i,j,d0=y.d1,d1=y.d0;
    x.resize(d0,d1);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) x.p[i*d1+j]=y.p[j*d0+i];
    return;
  }
  if(y.nd==1){
    x=y;
    x.reshape(1,y.N);
    return;
  }
  HALT("transpose not implemented for this dims");
}

//! returns the diagonal x = diag(y) (the diagonal-vector of the symmetric 2D matrix y)
template<class T> void getDiag(MT::Array<T>& x,const MT::Array<T>& y){
  CHECK(y.nd==2 && y.d0==y.d1,"can only give diagonal of symmetric 2D matrix");
  x.resize(y.d0);
  uint i;
  for(i=0;i<x.d0;i++) x(i)=y(i,i);
}

//! inverse of a 2d matrix
template<class T> void inverse2d(MT::Array<T>& invA,const MT::Array<T>& A){
  invA.resize(2,2);
  invA(0,0)=A(1,1); invA(1,1)=A(0,0); invA(0,1)=-A(0,1); invA(1,0)=-A(1,0);
  invA/=(T)(A(0,0)*A(1,1)-A(0,1)*A(1,0));
}



//===========================================================================
//
//!@name probability distribution operations
//

/*!\brief entropy \f$H_i = - \sum_x p(X_i=x) \ln p(X_i=x)\f$, where \f$x\f$ ranges
  from \c 0 to \c range-1, of the \c ith variable */
template<class T> T entropy(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--;) if(v.p[i]) t-=(T)(v.p[i]*::log(v.p[i]));
  return (T)(t/MT_LN2);
}

//! v = v / sum(v)
template<class T> T normalizeDist(MT::Array<T>& v){
  T Z=sum(v);
  if(Z>1e-100) v/=Z; else v=(T)1./(T)v.N;
  return Z;
}

//! v = v / sum(v)
template<class T> void makeConditional(MT::Array<T>& P){
  MT_MSG("makeConditional: don't use this anymore because it normalizes over the second index!!!, rather use tensorCondNormalize and condition on _later_ indices");
  CHECK(P.nd==2,"");
  uint i,j;
  T pi;
  for(i=0;i<P.d0;i++){
    pi=(T)0;
    for(j=0;j<P.d1;j++) pi+=P(i,j);
    for(j=0;j<P.d1;j++) P(i,j) /= pi;
  }
}

//inline uintA TUP(uint i,uint j,uint k,uint l){                      uintA z(4); z(0)=i; z(1)=j; z(2)=k; z(3)=l; return z; }

//! check whether this is a distribution over the first index w.r.t. the later indices
template<class T> void checkNormalization(MT::Array<T>& v,real tol){
  real p;
  uint i,j,k,l;
  switch(v.nd){
  case 1:
    for(p=0,i=0;i<v.d0;i++) p+=v(i);
    CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v);
    break;
  case 2:
    for(j=0;j<v.d1;j++){
      for(p=0,i=0;i<v.d0;i++) p+=v(i,j);
      CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v);
    }
    break;
  case 3:
    for(j=0;j<v.d1;j++) for(k=0;k<v.d2;k++){
      for(p=0,i=0;i<v.d0;i++) p+=v(i,j,k);
      CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v);
    }
    break;
  case 4:
    for(j=0;j<v.d1;j++) for(k=0;k<v.d2;k++){
      for(l=0;l<v.N/(v.d0*v.d1*v.d2);l++){
        for(p=0,i=0;i<v.d0;i++) p+=v(TUP(i,j,k,l));
	CHECK(fabs(1.-p)<tol,"distribution is not normalized: "<<v << " " << p);
      }
    }
    break;
  default: 
    CHECK(false, "Not implemented yet");
    break;
  }
}

template<class T> void eliminate(MT::Array<T>& x,const MT::Array<T>& y,uint d){
  CHECK(y.nd==2,"only implemented for 2D yet");
  uint i,j;
  if(d==1){
    x.resize(y.d0); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) x(i)+=y(i,j);
  }
  if(d==0){
    x.resize(y.d1); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) x(j)+=y(i,j);
  }
}

template<class T> void eliminate(MT::Array<T>& x,const MT::Array<T>& y,uint d,uint e){
  CHECK(y.nd==3,"only implemented for 3D yet");
  uint i,j,k;
  if(d==1 && e==2){
    x.resize(y.d0); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(i)+=y(i,j,k);
  }
  if(d==0 && e==2){
    x.resize(y.d1); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(j)+=y(i,j,k);
  }
  if(d==0 && e==1){
    x.resize(y.d2); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(k)+=y(i,j,k);
  }
}

// Eliminates one-dimension, d, from a 3D-tensor, y, and puts the result in x.
template<class T> void eliminatePartial(MT::Array<T>& x,const MT::Array<T>& y,uint d){
  CHECK(y.nd==3,"only implemented for 3D yet");
  uint i,j,k;
  if(d==2){
    x.resize(y.d0,y.d1); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(i,j)+=y(i,j,k);
  }
  if(d==1){
    x.resize(y.d0,y.d2); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(i,k)+=y(i,j,k);
  }
  if(d==0){
    x.resize(y.d1,y.d2); x=(T)0;
    for(i=0;i<y.d0;i++) for(j=0;j<y.d1;j++) for(k=0;k<y.d2;k++) x(j,k)+=y(i,j,k);
  }
}

/*! stochastic universal sampling: given a probability vector p
  (that is positive in each component and sums to one), the returned
  array s is resized to length n and contains the states of n
  samples. */
void SUS(const arr& p,uint n,uintA& s);

/*! as above, but returns only the index of one sample */
uint SUS(const arr& p);


//===========================================================================
//
//!@name distances
//

//! \f$\sum_i (v^i-w^i)^2\f$
template<class T> 
T sqrDistance(const MT::Array<T>& v, const MT::Array<T>& w){
  CHECK(v.N==w.N, 
	"sqrDistance on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0);
  for(uint i=v.N;i--;){ d=v.p[i]-w.p[i]; t+=d*d; }
  return t;
}

template<class T> T maxDiff(const MT::Array<T>& v,const MT::Array<T>& w,uint *im){
  CHECK(v.N==w.N, 
	"maxDiff on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0);
  if(!im)
      for(uint i=v.N;i--;){ 
        d=(T)::fabs(v.p[i]-w.p[i]);
        if(d>t) t=d; }
  else{
    *im=0;
    for(uint i=v.N;i--;){ d=(T)::fabs(v.p[i]-w.p[i]); if(d>t){ t=d; *im=i; } }
  }
  return t;
}

template<class T> T maxRelDiff(const MT::Array<T>& v,const MT::Array<T>& w,T tol){
  CHECK(v.N==w.N, 
	"maxDiff on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0),a,b,c;
  for(uint i=v.N;i--;){
    a=(T)::fabs(v.p[i]) + tol;
    b=(T)::fabs(w.p[i]) + tol;
    if(a<b){ c=a; a=b; b=c; }
    d=a/b-(T)1;
    if(d>t) t=d;
  }
  return t;
}

//! \f$\sum_{i|{\rm mask}_i={\rm true}} (v^i-w^i)^2\f$
/*template<class T> 
  T sqrDistance(const MT::Array<T>& v, const MT::Array<T>& w,const MT::Array<bool>& mask){
  CHECK(v.N==w.N, 
  "sqrDistance on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T d,t(0);
  for(uint i=v.N;i--;) if(mask(i)){ d=v.p[i]-w.p[i]; t+=d*d; }
  return t;
  }*/


template<class T> MT::Array<T>& minusA(MT::Array<T>& x,const MT::Array<T>& y,const MT::Array<T>& z);

//! \f$\sqrt{\sum_{ij} g_{ij} (v^i-w^i) (v^j-w^j)}\f$
template<class T> T sqrDistance(const MT::Array<T>& g,const MT::Array<T>& v,const MT::Array<T>& w){
  MT::Array<T> d;
  ::minusA(d,v,w);
  return scalarProduct(g,d,d);
}

//! \f$\sqrt{\sum_i (v^i-w^i)^2}\f$
template<class T> 
T euclideanDistance(const MT::Array<T>& v, const MT::Array<T>& w){
  return (T)::sqrt(sqrDistance(v,w));
}

//! \f$\sqrt{\sum_i (v^i-w^i)^2}\f$
template<class T> 
T metricDistance(const MT::Array<T>& g,const MT::Array<T>& v, const MT::Array<T>& w){
  return (T)::sqrt(sqrDistance(g,v,w));
}


//===========================================================================
//
//!@name running sums
//

//! \f$\sum_i x_i\f$
template<class T> T sum(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--; t+=v.p[i]){};
  return t;
}

//! \f$\sum_i |x_i|\f$
template<class T> T sumOfAbs(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--; t+=(T)::fabs(v.p[i])){};
  return t;
}

//! \f$\sum_i x_i^2\f$
template<class T> T sumOfSqr(const MT::Array<T>& v){
  T t(0);
  for(uint i=v.N; i--; t+=v.p[i]*v.p[i]){};
  return t;
}

//! \f$\sqrt{\sum_i x_i^2}\f$
template<class T> T norm(const MT::Array<T>& v){ return (T)::sqrt(sumOfSqr(v)); }

//! \f$\sqrt{\sum_i x_i^2}\f$
template<class T> T mean(const MT::Array<T>& v){ return sum(v)/v.N; }

template<class T> T var(const MT::Array<T>& v){ T m=mean(v); return sumOfSqr(v)/v.N-m*m; }

//! \f$\sum_i x_{ii}\f$
template<class T> T trace(const MT::Array<T>& v){
  CHECK(v.nd==2 && v.d0==v.d1,"only for squared matrix");
  T t(0);
  for(uint i=0;i<v.d0;i++) t+=v(i,i);
  return t;
}

template<class T> T minDiag(const MT::Array<T>& v){
  CHECK(v.nd==2 && v.d0==v.d1,"only for squared matrix");
  T t=v(0,0);
  for(uint i=1;i<v.d0;i++) if(v(i,i)<t) t=v(i,i);
  return t;
}



//===========================================================================
//
//!@name products
//

//! \f$\prod_i x_i\f$
template<class T> T product(const MT::Array<T>& v){
  T t(1);
  for(uint i=v.N; i--; t *= v.p[i]);
  return t;
}

/*!\brief inner product (also ordinary matrix or scalar product):
  \f$\forall_{ik}:~ x_{ik} = \sum_j v_{ij}\, w_{jk}\f$ but also:
  \f$\forall_{i}:~ x_{i} = \sum_j v_{ij}\, w_{j}\f$*/
template<class T> 
void innerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  /*
    if(y.nd==2 && z.nd==2 && y.N==z.N && y.d1==1 && z.d1==1){  //elem-wise
    HALT("make element-wise multiplication explicite!");
    mult(x,y,z);
    return;
    }
  */
  if(y.nd==2 && z.nd==1){ //matrix x vector -> vector
    CHECK(y.d1==z.d0,"wrong dimensions for inner product");
#ifdef MT_LAPACK
    if(MT::useLapack && typeid(T)==typeid(real)){ blas_Mv(x,y,z); return; }
#endif
    uint i,d0=y.d0,dk=y.d1;
    T *a,*astop,*b,*c;
    x.resize(d0); x.setZero();
    c=x.p;
    for(i=0;i<d0;i++){
      //for(s=0.,k=0;k<dk;k++) s+=y.p[i*dk+k]*z.p[k];
      //this is faster:
      a=y.p+i*dk; astop=a+dk; b=z.p;
      for(;a!=astop; a++, b++) (*c)+=(*a) * (*b);
      c++;
    }
    return;
  }
  if(y.nd==2 && z.nd==2){ //plain matrix multiplication
    CHECK(y.d1==z.d0,"wrong dimensions for inner product");
#ifdef MT_LAPACK
    if(MT::useLapack && typeid(T)==typeid(real)){ blas_MM(x,y,z); return; }
#endif
    uint i,j,d0=y.d0,d1=z.d1,dk=y.d1;
#if 0
    if(y.mtype==MT::Array<T>::diag){
      x.resize(d0,d1);
      for(i=0;i<d0;i++) for(j=0;j<d1;j++) x(i,j) = y(i,i) * z(i,j);
      return;
    }
    if(z.mtype==MT::Array<T>::diag){
      x.resize(d0,d1);
      for(i=0;i<d0;i++) for(j=0;j<d1;j++) x(i,j) *= y(i,j) * z(j,j);
      return;
    }
#endif
    T *a,*astop,*b,*c;
    x.resize(d0,d1); x.setZero();
    c=x.p;
    for(i=0;i<d0;i++) for(j=0;j<d1;j++){
      //for(s=0.,k=0;k<dk;k++) s+=y.p[i*dk+k]*z.p[k*d1+j];
      //this is faster:
      a=y.p+i*dk; astop=a+dk; b=z.p+j;
      for(;a!=astop; a++, b+=d1) (*c)+=(*a) * (*b);
      c++;
    }
    return;
  }
  if(y.nd==1 && z.nd==1 && z.N==1){ //vector multiplied with scalar (disguised as 1D vector)
    uint k,dk=y.N;
    x.resize(y.N);
    for(k=0;k<dk;k++) x.p[k]=y.p[k]*z.p[0];
    return;
  }
  /*if(y.nd==1 && z.nd==2 && z.d0==1){ //vector x vector^T -> matrix (outer product)
    CHECK(y.d0==z.d1,"wrong dimensions for inner product");
    uint i,j,d0=y.d0,d1=z.d1;
    x.resize(d0,d1);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) x(i,j)=y(i)*z(0,j);
    return;
  }*/
  if(y.nd==1 && z.nd==2){ //vector^T x matrix -> vector^T
    CHECK(y.d0==z.d0,"wrong dimensions for inner product");
    uint i,k,d0=z.d1,dk=y.d0;
    x.resize(d0);
    T s;
    for(i=0;i<d0;i++){
      for(s=0,k=0;k<dk;k++) s+=y.p[k]*z.p[k*d0+i];
      x.p[i]=s;
    }
    return;
  }
  if(y.nd==2 && z.nd==3){
    MT::Array<T> zz; zz.referTo(z);
    zz.reshape(z.d0,z.d1*z.d2);
    innerProduct(x,y,zz);
    x.reshape(y.d0,z.d1,z.d2);
    return;
  }
  if(y.nd==3 && z.nd==2){
    MT::Array<T> yy; yy.referTo(y);
    yy.reshape(y.d0*y.d1,y.d2);
    innerProduct(x,yy,z);
    x.reshape(y.d0,y.d1,z.d1);
    return;
  }
  if(y.nd==3 && z.nd==1){
    MT::Array<T> yy; yy.referTo(y);
    yy.reshape(y.d0*y.d1,y.d2);
    innerProduct(x,yy,z);
    x.reshape(y.d0,y.d1);
    return;
  }
  if(y.nd==1 && z.nd==3){
    MT::Array<T> zz; zz.referTo(z);
    zz.reshape(z.d0,z.d1*z.d2);
    innerProduct(x,y,zz);
    x.reshape(z.d1,z.d2);
    return;
  }
  if(y.nd==1 && z.nd==1){ //should be scalar product, but be careful
    HALT("what do you want? scalar product or element wise multiplication?");
    CHECK(y.d0==z.d0,"wrong dimensions for inner product");
    uint k,dk=y.d0;
    x.resize(1);
    T s;
    for(s=0,k=0;k<dk;k++) s+=y.p[k]*z.p[k];
    x.p[0]=s;
    return;
  }
  HALT("inner product - not yet implemented for these dimensions");
}

/*!\brief outer product (also exterior or tensor product): \f$\forall_{ijk}:~
  x_{ijk} = v_{ij}\, w_{k}\f$ */
template<class T> 
void outerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  if(y.nd==1 && z.nd==1){
    uint i,j,d0=y.d0,d1=z.d0;
    x.resize(d0,d1);
    for(i=0;i<d0;i++) for(j=0;j<d1;j++) x.p[i*d1+j]=y.p[i]*z.p[j];
    return;
  }
  HALT("outer product - not yet implemented for these dimensions");
}

//! \f$\sum_i v_i\, w_i\f$, or \f$\sum_{ij} v_{ij}\, w_{ij}\f$, etc.
template<class T> 
T scalarProduct(const MT::Array<T>& v, const MT::Array<T>& w){
  CHECK(v.N==w.N, 
	"scalar product on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T t(0);
  for(uint i=v.N; i--; t+=v.p[i]*w.p[i]);
  return t;
}

//! \f$\sum_{ij} g_{ij}\, v_i\, w_i\f$
template<class T> 
T scalarProduct(const MT::Array<T>& g,const MT::Array<T>& v, const MT::Array<T>& w){
  CHECK(v.N==w.N && g.nd==2 && g.d0==v.N && g.d1==w.N,
	"scalar product on different array dimensions ("<<v.N<<", "<<w.N<<")");
  T t(0);
  uint i,j;
  T *gp=g.p,*vp=v.p;
  for(i=0;i<g.d0;i++){
    for(j=0;j<g.d1;j++){
      t+=(*gp)*(*vp)*w.p[j];
      gp++;
    }
    vp++;
  }
  return t;
}


//===========================================================================
//
//!@name tensor operations
//

template<class T> std::ostream& operator<<(std::ostream& os,const MT::Array<T>& x);
void getIndexTuple(uintA &I,uint i,const uintA &d);

#define DEBUG_TENSOR(x) //x
//!@name tensor

/*! makes X to be a distribution over the left leftmost-indexed
  variables and normalizes it */
template<class T> void tensorCondNormalize(MT::Array<T>& X,int left){
  uint i,j,dl=1,dr;
  T sum;
  if(left>=0){ //normalize over the left variables
    for(j=0;j<(uint)left;j++) dl*=X.dim(j);
    dr=X.N/dl;
    CHECK(dl*dr==X.N,"");
    for(i=0;i<dr;i++){
      sum=(T)0;
      for(j=0;j<dl;j++)  sum += X.p[j*dr + i];
      if(sum) for(j=0;j<dl;j++) X.p[j*dr + i] /= sum;
      else    for(j=0;j<dl;j++) X.p[j*dr + i] = (T)1/dl;
    }
  }else{ //normalize over the right variables
    for(j=0;j<(uint)-left;j++) dl*=X.dim(j);
    dr=X.N/dl;
    CHECK(dl*dr==X.N,"");
    for(i=0;i<dl;i++){
      sum=(T)0;
      for(j=0;j<dr;j++)  sum += X.p[i*dr + j];
      if(sum) for(j=0;j<dr;j++) X.p[i*dr + j] /= sum;
      else    for(j=0;j<dr;j++) X.p[i*dr + j] = (T)1/dr;
    }
  }
}

/*! makes X to be a distribution over the left leftmost-indexed
  variables and normalizes it */
template<class T> void tensorCondMax(MT::Array<T>& X,uint left){
  uint i,j,dl=1,dr,jmax;
  T pmax;
  for(j=0;j<left;j++) dl*=X.dim(j);
  dr=X.N/dl;
  CHECK(dl*dr==X.N,"");
  for(i=0;i<dr;i++){
    jmax=0;
    pmax=X.p[jmax*dr + i];
    X.p[i]=(T)0;
    for(j=1;j<dl;j++){
      if(X.p[j*dr + i]>pmax){  jmax=j;  pmax=X.p[jmax*dr + i];  }
      X.p[j*dr + i]=(T)0;
    }
    X.p[jmax*dr + i]=(T)1;
  }
}

template<class T> void tensorCond11Rule(MT::Array<T>& X,uint left,double rate){
  uint i,j,dl=1,dr,jmax1,jmax2;
  for(j=0;j<left;j++) dl*=X.dim(j);
  dr=X.N/dl;
  CHECK(dl*dr==X.N,"");
  arr X_i(dl);
  double amin=10.;
  for(i=0;i<dr;i++){
    for(j=0;j<dl;j++) X_i(j) = X.p[j*dr + i]; //copy this row into one array
    //find the 1st and 2nd highest in this conditional part
    X_i.maxIndeces(jmax1,jmax2);
    if(X_i(jmax1)==X_i(jmax2)) continue;
    CHECK(X_i(jmax1)>X_i(jmax2),"must be really greater...");
    //compute the log ratio
    double a = ::log(rate)/::log(X_i(jmax1)/X_i(jmax2));
    //rescale everything [no -> globally rescale]
    //for(j=0;j<dl;j++)  X.p[j*dr + i] = pow(X.p[j*dr + i],a);
    if(a<amin) amin=a;
  }
  if(amin>1.) for(i=0;i<X.N;i++) X.elem(i) = pow(X.elem(i),amin);
}

/*! makes X to be a distribution over the left leftmost-indexed
  variables and normalizes it */
template<class T> void tensorCondSoftMax(MT::Array<T>& X,uint left,real beta){
  uint i;
  for(i=0;i<X.N;i++) X.p[i] = (T)::exp(beta*X.p[i]);
  tensorCondNormalize(X,left);
}

/*!\brief checks whether X is a normalized distribution over the left leftmost-indexed
  variables */
template<class T> void tensorCheckCondNormalization(const MT::Array<T>& X,uint left,real tol){
  uint i,j,dl=1,dr;
  real sum;
  for(j=0;j<left;j++) dl*=X.dim(j);
  dr=X.N/dl;
  CHECK(dl*dr==X.N,"");
  for(i=0;i<dr;i++){
    sum=0.;
    for(j=0;j<dl;j++) sum += X.p[j*dr + i]; 
    CHECK(fabs(1.-sum)<tol,"distribution is not normalized: "<<X);
  }
}

template<class T> void tensorCheckCondNormalization_with_logP(const MT::Array<T>& X,uint left, real logP, real tol){
  uint i,j,dl=1,dr;
  real sum, coeff=::exp(logP);
  for(j=0;j<left;j++) dl*=X.dim(j);
  dr=X.N/dl;
  CHECK(dl*dr==X.N,"");
  for(i=0;i<dr;i++){
    sum=0.;
    uintA checkedIds;
    for(j=0;j<dl;j++){ sum += X.p[j*dr + i]*coeff; checkedIds.append(j*dr + i); }
    CHECK(fabs(1.-sum)<tol,"distribution is not normalized for parents-config#"<<i<<endl<<checkedIds<<endl<<" "<<X);
  }
}

/*! X becomes a tensor product (maybe with some variables summed out)
  of A and B. pickA and pickB indicate which slots of A and B are
  associated with which slots of C and the summation sign. More
  precisely, if we have \f$C_{i_0i_1i_2} = \sum_{i_3i_4}
  A_{i_4i_2i_1}~ B_{i_3i_0}\f$ then you should call
  tensor(C,A,TUP(4,2,1),B,TUP(3,0),2); Here, the `2` indicates that
  the last two indices of i_0,..,i_4 are summed over, and C only
  becomes a 3rd rank instead of 5th rank tensor */
template<class T> void tensorEquation(MT::Array<T> &X,const MT::Array<T> &A,const uintA &pickA,const MT::Array<T> &B,const uintA &pickB,uint sum){
  CHECK(&X!=&A && &X!=&B,"output tensor must be different from input tensors");
  CHECK(A.nd==pickA.N && B.nd==pickB.N,"miss-sized tensor references: " <<A.nd <<"!=" <<pickA.N <<" " <<B.nd <<"!=" <<pickB.N);
  
  uint n=1+MT::MAX(pickA.max(),pickB.max());
  uint i,j,r,s,N,res;
  intA a(n),b(n);
  uintA d(n),dx(n-sum),I,Ia(A.nd),Ib(B.nd);
  
  DEBUG_TENSOR(cout <<"pickA=" <<pickA <<" pickB=" <<pickB <<endl;);
  
  // permutation for A
  a=-1;
  for(i=0;i<A.nd;i++) a(pickA(i))=i;
  //j=A.nd;
  //for(i=0;i<n;i++) if(a(i)==-1){ a(i)=j; j++;  }
  DEBUG_TENSOR(cout <<"permutation for A: " <<a <<endl;);
  
  //permutation for B
  b=-1;
  for(i=0;i<B.nd;i++) b(pickB(i))=i;
  //j=B.nd;
  //for(i=0;i<n;i++) if(b(i)==-1){ b(i)=j; j++; }
  DEBUG_TENSOR(cout <<"permutation for B: " <<b <<endl;);
  
  //dimensionalities
  for(i=0;i<n;i++){
    if(a(i)!=-1) r=A.dim(a(i)); else r=0;
    if(b(i)!=-1) s=B.dim(b(i)); else s=0;
    CHECK(!r || !s || r==s,"inconsistent sharing dimensionalities: " <<r <<"!=" <<s);
    d(i)=MT::MAX(r,s);
  }
  DEBUG_TENSOR(cout <<"full dimensionality d=" <<d <<endl;);
  
  //total elements:
  N=product(d);
  if(!sum){
    res=1;
    //X.resizeTensor(d);
    CHECK(d==X.d,"for security, please set size before");
  }else{
    dx.resize(d.N-sum);
    res=1;
    for(j=0;j<dx.N;j++) dx(j)=d(j);
    for(   ;j<d .N;j++) res*=d(j);
    //X.resizeTensor(dx);
    CHECK(dx==X.d,"for security, please set size before");
  }
  CHECK(N==X.N*res,"");
  DEBUG_TENSOR(cout <<"dx=" <<dx <<" res=" <<res <<endl; );
  
  //here the copying and multiplying takes place...
  X.setZero();
  for(i=0;i<N;i++){
    getIndexTuple(I,i,d);
    for(j=0;j<A.nd;j++) Ia(j)=I(pickA(j));
    for(j=0;j<B.nd;j++) Ib(j)=I(pickB(j));
    //DEBUG_TENSOR(cout <<"i="<<i <<" I="<<I <<" i/res="<<i/res <<" Ia="<<Ia <<" Ib="<<Ib <<endl;)
    if(!sum){
      X.elem(i) = A(Ia) * B(Ib);
    }else{
      X.elem(i/res) += A(Ia) * B(Ib);
    }
  }
}

template<class T> void tensorEquation_doesntWorkLikeThat(MT::Array<T> &X,const MT::Array<T> &A,const uintA &pickA,const MT::Array<T> &B,const uintA &pickB,uint sum){
  CHECK(A.nd==pickA.N && B.nd==pickB.N,"miss-sized tensor references: " <<A.nd <<"!=" <<pickA.N <<" " <<B.nd <<"!=" <<pickB.N);
  
  uint n=1+MT::MAX(pickA.max(),pickB.max());
  uint i,j;
  uintA a(n),b(n);
  
  DEBUG_TENSOR(cout <<"pickA=" <<pickA <<" pickB=" <<pickB <<endl;);
  
  // permutations for A & B
  a=-1;  for(i=0;i<A.nd;i++) a(pickA(i))=i;
  b=-1;  for(i=0;i<B.nd;i++) b(pickB(i))=i;
  DEBUG_TENSOR(cout <<"permutation for A: " <<a <<"\npermutation for B: " <<b <<endl;);

  //permute tensors
  arr Aperm,Bperm;
  tensorPermutation(Aperm,A,a);
  tensorPermutation(Bperm,B,b);
  
  //dimensionalities: left-sum-right
  uint ldim=1,sdim=1,rdim=1;
  for(i=0;i<Aperm.nd-sum;i++){ ldim *= Aperm.d[i]; }
  for(i=0;i<sum;i++){ j = Aperm.d[sum+i]; CHECK(j==Bperm.d[i],""); sdim*=j; }
  for(i=0;i<Bperm.nd-sum;i++){ rdim *= Bperm.d[sum+i]; }
  DEBUG_TENSOR(cout <<"ldim=" <<ldim <<" sdim=" <<sdim <<" rdim=" <<rdim <<endl;)
      
  //reshape to matrices
  Aperm.reshape(ldim,sdim);
  Bperm.reshape(sdim,rdim);

  //matrix multiplication
  innerProduct(X,Aperm,Bperm);

  //reshape
}

inline void getStride(uintA& stride,uintA& dim){
  stride.resize(dim.N+1);
  stride(dim.N) = 1;
  for(uint i=dim.N;i--;){ stride(i) = stride(i+1)*dim(i); }
}

//index and limit are w.r.t is the GLOBAL indexing!, j_stride w.r.t. the permuted
inline void multiDimIncrement(uint& Ycount,uint* index,uint* limit,uint* Yinc,uint* Ydec,uint nd){
  uint k;
  for(k=nd;k--;){
    Ycount+=Yinc[k];
    index[k]++;
    if(index[k]<limit[k]) break;  //we need the index only to decide when to overflow -- is there a more efficient way?
    index[k]=0;
    Ycount-=Ydec[k];
  }
}

inline void getMultiDimIncrement(const uintA& Xdim,const uintA &Yid,uint* Ydim,uint* Yinc,uint* Ydec){
  uint i;
  memset(Ydim,0,sizeof(uint)*maxRank);  for(i=0;i<Xdim.N;i++) if(i<Yid.N) Ydim[i]=Xdim(Yid.p[i]);  //dimension of Y
  memset(Yinc,0,sizeof(uint)*maxRank);  Yinc[Yid.p[Yid.N-1]]=1;  for(i=Yid.N-1;i--;) Yinc[Yid.p[i]] = Ydim[i+1] * Yinc[Yid.p[i+1]];  //stride of Y
  for(i=Xdim.N;i--;) Ydec[i] = Xdim(i)*Yinc[i];
  //cout <<"Xdim=" <<Xdim <<"\nYid =" <<Yid <<"\nYdim=" <<uintA(Ydim,Yid.N) <<"\nYinc=" <<uintA(Yinc,Xdim.N) <<"\nYdec=" <<uintA(Ydec,Xdim.N) <<endl;
}

/*! \f$Y_{i_Yid(0),i_Yid(1)} = \sum_{i_1} X_{i_0,i_1,i_2}\f$. Get the marginal Y
  from X, where Y will share the slots `Yid' with X */
template<class T> void tensorMarginal(MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //handle scalar case
  if(!Yid.N){  Y.resize(1);  Y.nd=0;  Y.scalar()=sum(X);  return;  }

  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);
  Y.resize(Yid.N,Ydim);
  Y.setZero();

  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
#if 0 //use this to check looping -- all routines below don't have this check anymore!
    cout <<"current globalIndex=" <<Xcount <<' ' <<Ycount <<' ' <<uintA(I,X.nd) <<endl;
    //check if incrementing Y worked out
    uint k,jj=0;
    for(k=0;k<Yid.N;k++){ jj*=Ydim[k]; jj+=I[Yid(k)]; }
    CHECK(jj==Ycount,"");
    //check if incrementing I worked out
    uintA II;
    getIndexTuple(II,Xcount,uintA(X.d,X.nd));
    CHECK(II==I,"not equal: "<<II <<uintA(I,X.nd));
#endif
    Y.p[Ycount] += X.p[Xcount];
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

/*! \f$Y_{i_Yid(0),i_Yid(1)} = \sum_{i_1} X_{i_0,i_1,i_2}\f$. Get the marginal Y
  from X, where Y will share the slots `Yid' with X */
template<class T> void tensorPermutation(MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N==X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");
  
  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);
  Y.resize(Yid.N,Ydim);
  
  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
    Y.p[Ycount] = X.p[Xcount];
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

                        
/*! \f$Y_{i_2,i_0} = {\rm max}_{i_1} X_{i_0,i_1,i_2}\f$. Get the ``max-marginal'' Y
  from X, where Y will share the slots `Yid' with X (basis of max-product BP) */
template<class T> void tensorMaxMarginal(MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);
  Y.resize(Yid.N,Ydim);
  Y.setZero();
  HALT("WRONG IMPLEMENTATION! - zero doens't guarantee max...");

  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
    if(X.p[Xcount]>Y.p[Ycount]) Y.p[Ycount] = X.p[Xcount];
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

/*!\brief \f$X_{i_0,i_1,i_2} \gets X_{i_0,i_1,i_2}~ Y_{i_Yid(0),i_Yid(1)}\f$. Multiply Y onto X,
  where Y shares the slots `Yid' with X */
template<class T> void tensorAdd_old(MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N==Y.nd,"need to specify "<<Y.nd<<" slots, not "<<Yid.N);
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //handle scalar case
  if(!Yid.N){ CHECK(Y.N==1,"");  X+=Y.scalar();  return; } //Y is only a scalar

  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);
  Y.resize(Yid.N,Ydim);
  Y.setZero();

  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
    X.p[Xcount] += Y.p[Ycount];
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

template<class T> void tensorMarginal_old(MT::Array<T> &y,const MT::Array<T> &x,const uintA &xd,const uintA &ids){
  uint i,j,k,n=product(xd);
  CHECK(x.N==n,"");
  //size y
  uintA yd(ids.N);
  for(i=0;i<ids.N;i++) yd(i)=xd(ids(i));
  //y.resize(yd); y.setZero();
  y.resize(product(yd)); y.setZero();

  uintA xt(xd.N); xt.setZero();
  for(i=0;i<n;i++){
    //compute j
    for(j=0,k=0;k<ids.N;k++){ j*=yd(k); j+=xt(ids(k)); }
    //uintA yt(ids.N); for(k=ids.N;k--;){ yt(k)=xt(ids(k)); }
    //cout <<"i=" <<i <<" j=" <<j <<" yt=" <<yt <<" xt=" <<xt <<endl;
    y(j)+=x.elem(i);
    //increment xt
    for(k=xt.N;k--;){ xt(k)++; if(xt(k)<xd(k)) break; else xt(k)=0; }
  }
}

/*! \f$X_{i_0,i_1,i_2} \gets X_{i_0,i_1,i_2}~ Y_{i_Yid(0),i_Yid(1)}\f$. Multiply Y onto X,
  where Y shares the slots `Yid' with X */
template<class T> void tensorMultiply(MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N==Y.nd,"need to specify "<<Y.nd<<" slots, not "<<Yid.N);
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //handle scalar case
  if(!Yid.N){ CHECK(Y.N==1,"");  X*=Y.scalar();  return; } //Y is only a scalar
  
  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);

  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
    X.p[Xcount] *= Y.p[Ycount];
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

/*! \f$X_{i_0,i_1,i_2} \gets X_{i_0,i_1,i_2}~ Y_{i_Yid(0),i_Yid(1)}\f$. Multiply Y onto X,
  where Y shares the slots `Yid' with X */
// TODO cope with division by 0, in particular 0/0
template<class T> void tensorDivide(MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N==Y.nd,"need to specify "<<Y.nd<<" slots, not "<<Yid.N);
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //handle scalar case
  if(!Yid.N){ CHECK(Y.N==1,"");  X/=Y.scalar();  return; } //Y is only a scalar
  
  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);

  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
    // TODO division by zero??
    X.p[Xcount] = MT::DIV(X.p[Xcount], Y.p[Ycount]);
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

template<class T> void tensorAdd(MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid){
  uint Xcount,Ycount;
  CHECK(Yid.N==Y.nd,"need to specify "<<Y.nd<<" slots, not "<<Yid.N);
  CHECK(Yid.N<=X.nd,"can't take slots "<<Yid<<" from "<<X.nd<<"D tensor");

  //handle scalar case
  if(!Yid.N){ CHECK(Y.N==1,"");  X+=Y.scalar();  return; } //Y is only a scalar
  
  //initialize looping
  uint I[maxRank];     memset(I,0,sizeof(uint)*maxRank);  //index on X
  uint Ydim[maxRank],Yinc[maxRank],Ydec[maxRank];
  getMultiDimIncrement(X.getDim(),Yid,Ydim,Yinc,Ydec);

  //loop
  for(Xcount=0,Ycount=0;Xcount<X.N;Xcount++){
    X.p[Xcount] += Y.p[Ycount];
    multiDimIncrement(Ycount,I,X.d,Yinc,Ydec,X.nd);
  }
}

/*! multiply y onto x, where x has dimensions `d' and y shares the
  dimensions `ids' with x */
template<class T> void tensorMultiply_old(MT::Array<T> &x,const MT::Array<T> &y,const uintA &d,const uintA &ids){
  uint i,j,k,n=x.N;
  CHECK(n==product(d),"");

  uintA yd(ids.N);
  for(i=0;i<ids.N;i++) yd(i)=d(ids(i));
  CHECK(y.N==product(yd),"");

  uintA I(d.N); I.setZero();
  for(i=0;i<n;i++){
    for(j=0,k=0;k<ids.N;k++){ j*=yd(k); j+=I(ids(k)); }
    x.elem(i) *= y.elem(j);
    for(k=I.N;k--;){ I(k)++; if(I(k)<d(k)) break; else I(k)=0; }
  }
}



//===========================================================================
//
//!@name set operations
//

//! x becomes the section of y and z
template<class T> 
void setSection(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  uint i,j;
  x.clear();
  for(i=0;i<y.N;i++){
    for(j=0;j<z.N;j++) if(y(i)==z(j)){ x.append(y(i)); break; }
  }
}

//! x becomes the section of y and z
template<class T> 
void setUnion(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z){
  uint i,j;
  if(&x!=&y) x=y;
  for(i=0;i<z.N;i++){
    for(j=0;j<y.N;j++) if(z(i)==y(j)) break;
    if(j==y.N) x.append(z(i));
  }
}

//! x becomes the section of y and z
template<class T> 
void setMinus(MT::Array<T>& x,const MT::Array<T>& y){
  uint i,j;
  for(i=0;i<x.N;){
    for(j=0;j<y.N;j++){
      if(x(i)==y(j)){
	x.remove(i);
	break;
      }
    }
    if(j==y.N) i++;
  }
}

//! share x and y at least one element?
template<class T> uint numberSharedElements(const MT::Array<T>& x,const MT::Array<T>& y) {
  MT::Array<T> z;
  setSection(z, x, y);
  return z.N;
}


//===========================================================================
//
//!@name randomizations
//

//! Assign all elements of \c a to a uniformly distributed discrete value in {low,..,hi}
template<class T> void rndInteger(MT::Array<T>& a,int low,int high,bool add){
  if(!add) for(uint i=0;i<a.N;i++) a.p[i] =low+(int)rnd.num(1+high-low);
  else     for(uint i=0;i<a.N;i++) a.p[i]+=low+(int)rnd.num(1+high-low);
}

//! Assign all elements of \c a to a uniformly distributed continuous value in [low,hi]
template<class T> void rndUniform(MT::Array<T>& a,double low,double high,bool add){
  if(!add) for(uint i=0;i<a.N;i++) a.p[i] =(T)rnd.uni(low,high);
  else     for(uint i=0;i<a.N;i++) a.p[i]+=(T)rnd.uni(low,high);
}

template<class T> void rndNegLogUniform(MT::Array<T>& a,double low,double high,bool add){
  if(!add) for(uint i=0;i<a.N;i++) a.p[i] =(T)(-::log(rnd.uni(low,high)));
  else     for(uint i=0;i<a.N;i++) a.p[i]+=(T)(-::log(rnd.uni(low,high)));
}

//! Assign all elements of \c a to a Gaussian random variable
//! where _each_ elelemt has the given stdDev
//! -- devide by sqrt(N) if you want the multivariate Gaussian
//! to have a given standard deviation!
template<class T> void rndGauss(MT::Array<T>& a,real stdDev,bool add){
  //MT_MSG("rndGauss(Array): watch out what sdv means - inconsistent...");
  /*if(a.nd==2 && a.d0!=a.d1)  //assumes that it is a batch of data
    stdDev/=::sqrt((real)a.d1);
  else
    stdDev/=::sqrt((real)a.N);
  */
  if(!add) for(uint i=0;i<a.N;i++) a.p[i] =(T)(stdDev*rnd.gauss());
  else     for(uint i=0;i<a.N;i++) a.p[i]+=(T)(stdDev*rnd.gauss());
}


//! a gaussian random vector with Id covariance matrix (sdv = sqrt(dimension))
/*template<class T> void rndGauss(MT::Array<T>& a,bool add){
  if(!add) for(uint i=0;i<a.N;i++) a.p[i]=rnd.gauss();
  else     for(uint i=0;i<a.N;i++) a.p[i]+=rnd.gauss();
  }*/

//! returns an array with \c dim Gaussian noise elements
/*template<class T> MT::Array<T>& rndGauss(real stdDev,uint dim){
  static MT::Array<T> z;
  stdDev/=::sqrt(dim);
  z.resize(dim);
  rndGauss(z,stdDev);
  return z;
  }*/

/*!\brief from a vector of numbers, calculates the softmax distribution
  soft(i) = exp(beta*a(i)), and returns a random sample from
  this distribution (an index in {0,..,a.N-1}) */
template<class T> uint softMax(const MT::Array<T>& a,arr& soft,real beta){
  real norm=0.,r;
  uint i; int sel=-1;
  resizeAs(soft,a);
  for(i=0;i<a.N;i++){
    soft(i)=exp(beta*a(i));
    norm+=soft(i);
  }
  r=rnd.uni();
  for(i=0;i<a.N;i++){
    soft(i)/=norm;
    r-=soft(i);
    if(sel==-1 && r<0.) sel=i;
  }
  return sel;
}


//===========================================================================
//
//!@name certain initializations
//








namespace MT{
  //! transpose
  template<class T> Array<T> operator~(const Array<T>& y){ Array<T> x; transpose(x,y); return x; }
  //! negative
  template<class T> Array<T> operator-(const Array<T>& y){ Array<T> x; negative(x,y);  return x; }
  //! outer product (notation borrowed from the wedge product, though not anti-symmetric)
  template<class T> Array<T> operator^(const Array<T>& y,const Array<T>& z){ Array<T> x; outerProduct(x,y,z); return x; }

  //! inner product
  template<class T> Array<T> operator*(const Array<T>& y,const Array<T>& z){ Array<T> x; innerProduct(x,y,z); return x; }
  //! scalar multiplication
  template<class T> Array<T> operator*(const Array<T>& y,T z){             Array<T> x; multS(x,y,z); return x; }
  //! scalar multiplication
  template<class T> Array<T> operator*(T y,const Array<T>& z){             Array<T> x; multS(x,y,z); return x; }


#define BinaryOperator( op, name)					\
  template<class T> Array<T> operator op(const Array<T>& y,const Array<T>& z){ Array<T> x; name(x,y,z); return x; } \
    template<class T> Array<T> operator op(T y,const Array<T>& z){             Array<T> x; name##S(x,y,z); return x; } \
    template<class T> Array<T> operator op(const Array<T>& y,T z){             Array<T> x; name##S(x,y,z); return x; } \

  BinaryOperator( + , plusA)
    BinaryOperator( - , minusA)
    BinaryOperator( % , mult)
    BinaryOperator( / , div)
#undef BinaryOperator

#define CompoundAssignmentOperator( op )				\
    template<class T> Array<T>& operator op (Array<T>& x,const Array<T>& y){ \
      CHECK(x.N==y.N,							\
	    "binary operator on different array dimensions ("<<x.N<<", "<<y.N<<")"); \
      T *xp=x.p;							\
      const T *yp=y.p;							\
      for(; xp!=x.pstop; xp++,yp++) *xp op *yp;				\
      return x;						\
    }									\
									\
    template<class T> Array<T>& operator op ( Array<T>& x, T y ){	\
      T *xp=x.p;							\
      for(; xp!=x.pstop; xp++) *xp op y;				\
      return x;						\
    }

    CompoundAssignmentOperator( |= )
    CompoundAssignmentOperator( ^= )
      CompoundAssignmentOperator( &= )
      CompoundAssignmentOperator( += )
      CompoundAssignmentOperator( -= )
      CompoundAssignmentOperator( *= )
      CompoundAssignmentOperator( /= )
      CompoundAssignmentOperator( %= )
#undef CompoundAssignmentOperator
}


//! calls MT::Array<T>::read
template<class T> std::istream& operator>>(std::istream& is,MT::Array<T>& x){ x.read(is);return is; }

//! allows a notation such as x <<"[0 1; 2 3]"; to initialize an array x 
template<class T> MT::Array<T>& operator<<(MT::Array<T>& x,const char* str){ std::istringstream ss(str); ss >>x; return x; }

//! calls MT::Array<T>::write
template<class T> std::ostream& operator<<(std::ostream& os,const MT::Array<T>& x){
  x.write(os); return os;
}

//! check for Nans in the array (checks x.elem(i)==x.elem(i) for all elements)
template<class T> void checkNan(const MT::Array<T>& x){
  for(uint i=0;i<x.N;i++){
    //CHECK(x.elem(i)!=NAN,"found a NaN" <<x.elem(i) <<'[' <<i <<']');
    CHECK(x.elem(i)==x.elem(i),"inconsistent number: " <<x.elem(i) <<'[' <<i <<']');
  }
}

//! equal in size and all elements
template<class T> bool operator==(const MT::Array<T>& v,const MT::Array<T>& w){
  if(!samedim(v,w)) return false;
  const T *iv,*iw;
  for(iv=v.p,iw=w.p; iv!=v.pstop; iv++,iw++)
    if (*iv != *iw) return false;
  return true;
}

//! equal in size and all elements
template<class T> bool operator==(const MT::Array<T>& v,const T *w){
  const T *iv,*iw;
  for(iv=v.p,iw=w; iv!=v.pstop; iv++,iw++)
    if (*iv != *iw) return false;
  return true;
}

//! not equal
template<class T> bool operator!=(const MT::Array<T>& v,const MT::Array<T>& w){
  return !(v==w);
}

//! lexical comparison
template<class T> bool operator<(const MT::Array<T>& v,const MT::Array<T>& w){
  if(v.N==w.N){
    for(uint i=0;i<v.N;i++){
      if(v.p[i]>w.p[i]) return false;
      if(v.p[i]<w.p[i]) return true;
    }
    return false; //they are equal
  }
  return v.N<w.N;
}


//===========================================================================
//
//!@name arithmetic operators
//


//---------- unary operators

#define UnaryOperation( name, op )					\
  template<class T>							\
  MT::Array<T>& name (MT::Array<T>& x,const MT::Array<T>& y){	\
    x.resizeAs(y);							\
    T *xp=x.p;								\
    const T *iy=y.p;							\
    for(; iy!=y.pstop; iy++,xp++) *xp= op *iy;				\
    return x;								\
  }

UnaryOperation( negative, - );
#undef UnaryOperator

//---------- binary function

#define BinaryOperation( name, op )					\
  template<class T>							\
  MT::Array<T>& name(MT::Array<T>& x,const MT::Array<T>& y,const MT::Array<T>& z){ \
    CHECK(y.N==z.N,							\
	  "binary operator on different array dimensions ("<<y.N<<", "<<z.N<<")"); \
    if(&x!=&y) x.resizeAs(y);							\
    T *xp=x.p;								\
    const T *zp=z.p,*yp=y.p;						\
    for(; xp!=x.pstop; xp++,yp++,zp++) *xp = *yp op *zp;		\
    return x;								\
  }									\
									\
  template<class T>							\
  MT::Array<T>& name##S(MT::Array<T>& x,const MT::Array<T>& y,T z){ \
    if(&x!=&y) x.resizeAs(y);							\
    T *xp=x.p;								\
    const T *yp=y.p;							\
    for(; xp!=x.pstop; xp++,yp++) *xp = *yp op z;			\
    return x;								\
  }									\
									\
  template<class T>							\
  MT::Array<T>& name##S(MT::Array<T>& x, T y, const MT::Array<T>& z){ \
    if(&x!=&z) x.resizeAs(z);							\
    T *xp=x.p;								\
    const T *zp=z.p;							\
    for(; xp!=x.pstop; xp++,zp++) *xp = y op *zp;			\
    return x;								\
  }

BinaryOperation( plusA , + );
BinaryOperation( minusA , - );
BinaryOperation( mult , * );
//     BinaryOperation( div , / )
#undef BinaryOperation


// To be able to cope with division by 0, we need to take the div out of the preprocessor-templates.
template<class T> MT::Array<T>& div(MT::Array<T>& x,const MT::Array<T>& y,const MT::Array<T>& z){
    CHECK(y.N==z.N,	"binary operator on different array dimensions ("<<y.N<<", "<<z.N<<")");
    if(&x!=&y) x.resizeAs(y);
    T *xp=x.p;
    const T *zp=z.p,*yp=y.p;
    for(; xp!=x.pstop; xp++,yp++,zp++) *xp = MT::DIV(*yp, *zp, true);
    return x;
  }

template<class T> MT::Array<T>& divS(MT::Array<T>& x,const MT::Array<T>& y,T z){
    if(&x!=&y) x.resizeAs(y);
    T *xp=x.p;
    const T *yp=y.p;
    for(; xp!=x.pstop; xp++,yp++) *xp = MT::DIV(*yp, z, true);
    return x;
  }

template<class T> MT::Array<T>& divS(MT::Array<T>& x, T y, const MT::Array<T>& z){
    if(&x!=&z) x.resizeAs(z);
    T *xp=x.p;
    const T *zp=z.p;
    for(; xp!=x.pstop; xp++,zp++) *xp = MT::DIV(y, *zp, true);
    return x;
  }





  //---------- compound assignment operators



  //---------- unary functions

#define UnaryFunction( func )					\
    template<class T>						\
    MT::Array<T> func (const MT::Array<T>& y){		\
      MT::Array<T> x;						\
      if(&x!=&y) x.resizeAs(y);					\
      T *xp=x.p;						\
      const T *yp=y.p;						\
      for(; xp!=x.pstop; xp++,yp++) *xp = ::func( *yp );	\
      return x;					\
    }

// trigonometric functions
UnaryFunction( acos );
UnaryFunction( asin );
UnaryFunction( atan );
UnaryFunction( cos );
UnaryFunction( sin );
UnaryFunction( tan );

// hyperbolic functions
UnaryFunction( cosh );
UnaryFunction( sinh );
UnaryFunction( tanh );
UnaryFunction( acosh );
UnaryFunction( asinh );
UnaryFunction( atanh );

// exponential and logarithmic functions
UnaryFunction( exp );
UnaryFunction( log );
UnaryFunction( log10 );

//roots
UnaryFunction( sqrt );
UnaryFunction( cbrt );

// nearest integer and absolute value
UnaryFunction( ceil );
UnaryFunction( fabs );
UnaryFunction( floor );
#undef UnaryFunction


//---------- binary functions

#define BinaryFunction( func )						\
  template<class T>							\
  MT::Array<T> func(const MT::Array<T>& y,const MT::Array<T>& z){ \
    CHECK(y.N==z.N,							\
	  "binary operator on different array dimensions ("<<y.N<<", "<<z.N<<")"); \
    MT::Array<T> x;							\
    x.resizeAs(y);							\
    for(uint i=x.N;i--; ) x.p[i]= func(y.p[i],z.p[i]);			\
    return x;						\
  }									\
									\
  template<class T>							\
  MT::Array<T> func(const MT::Array<T>& y,T z){			\
    MT::Array<T> x;							\
    x.resizeAs(y);							\
    for(uint i=x.N;i--; ) x.p[i]= func(y.p[i],z);			\
    return x;						\
  }									\
									\
  template<class T>							\
  MT::Array<T> func(T y,const MT::Array<T>& z){			\
    MT::Array<T> x;							\
    x.resizeAs(z);							\
    for(uint i=x.N;i--; ) x.p[i]= func(y,z.p[i]);			\
    return x;						\
  }

BinaryFunction( atan2 );
BinaryFunction( pow );
BinaryFunction( fmod );
#undef BinaryFunction

#ifndef MT_doxy // exclude these macros when generating the documentation
  // (doxygen can't handle them...)
#endif //(doxygen exclusion)


  /*

  //! element-wise linear combination (plus with scalar factors for each array)
  template<class T> void plusSASA(MT::Array<T>& x,T a,const MT::Array<T>& y,T b,const MT::Array<T>& z){
  CHECK(y.N==z.N,"must have same size for adding!");
  uint i,n=y.N;
  x.resizeAs(y);
  for(i=0;i<n;i++) x.p[i]=a*y.p[i]+b*z.p[i];
  }*/



#if 0
  #define IMPLEMENT_Array(x) void implement_Array_##x(){ MT::Array<x> dummy; }
  #define IMPLEMENT_Array_(x,key) void implement_Array_##key(){ MT::Array<x> dummy; }

  template class MT::Array<double>;
  template class MT::Array<int>;
  template class MT::Array<uint>;
  template class MT::Array<float>;
  template class MT::Array<byte>;
  template class MT::Array<char>;
  template class MT::Array<bool>;
  #define T double
  #  include "array_instantiate.cpp"
  #define T int
  #  include "array_instantiate.cpp"
#endif

//===========================================================================
//
// lists
//

template<class T> void listWrite(const MT::Array<T*>& L,std::ostream& os,const char *ELEMSEP=" ",const char *delim=NULL){
  uint i;
  if(delim) os <<delim[0];
  for(i=0;i<L.N;i++){ if(i) os <<ELEMSEP;  if(L.elem(i)) os <<*L.elem(i); else os <<"<NULL>"; }
  if(delim) os <<delim[1] <<std::flush;
}
               
template<class T> void listWriteNames(const MT::Array<T*>& L,std::ostream& os){
  uint i;
  os <<'(';
  for(i=0;i<L.N;i++){ if(i) os <<' ';  if(L.elem(i)) os <<L.elem(i)->name; else os <<"<NULL>"; }
  os <<')' <<std::flush;
}

/*template<class T> void listRead(MT::Array<T*>& L,std::istream& is){
                 CHECK(!L.N,"delete the list before reading!");
                 MT::parse(is,"{");
                 char c;
                 for(;;){
                 c=MT::peerNextChar(is);
                 if(c=='}') break;
                 L.append(new T());
                 is >>*L.last();
}
}*/

template<class T> void listClone(MT::Array<T*>& L,const MT::Array<T*>& M){
  L.resize(M.N);
  uint i;
  for(i=0;i<L.N;i++) L(i)=M(i)->new_clone();
}

template<class T> void listCopy(MT::Array<T*>& L,const MT::Array<T*>& M){
  L.resize(M.N);
  uint i;
  for(i=0;i<L.N;i++) L(i)=new T(*M(i));
}

template<class T> void listDelete(MT::Array<T*>& L){
  uint i;
  for(i=0;i<L.N;i++) delete L.elem(i);
  L.clear();
}

template<class T> T* listFindName(const MT::Array<T*>& L,const char* name){
  uint i;
  T *e;
  for_list(i,e,L) if(!strcmp(e->name,name)) return e;
  std::cerr <<"\n*** name '"<<name<<"' not in this list!" <<std::endl;
  return NULL;
}

template<class T> T* listFindType(const MT::Array<T*>& L,const char* type){
  uint i;
  T *e;
  for_list(i,e,L) if(!strcmp(e->type,type)) return e;
  std::cerr <<"type '"<<type<<"' not in this list!" <<std::endl;
  return NULL;
}

template<class T,class LowerOperator> void listSort(MT::Array<T*>& L, LowerOperator lowerop){
  std::sort(L.p,L.pstop,lowerop);
  for(uint i=0;i<L.N;i++) L(i)->index=i;
}


//===========================================================================
//
// graphs
//

template<class vert,class edge> void graphDelete(MT::Array<vert*>& V,MT::Array<edge*>& E){
  listDelete(E);
  listDelete(V);
}

template<class vert,class edge> edge* graphGetEdge(vert *from, vert *to){
  edge *e;
  uint i;
  for_list(i,e,from->outLinks) if(e->from==to) return e;
  return NULL;
}

template<class vert,class edge> void graphRandomUndirected(MT::Array<vert*>& V,MT::Array<edge*>& E,uint N,real connectivity){
  uint i,j;
  for(i=0;i<N;i++) V.append(new vert);
  for(i=0;i<N;i++) for(j=i+1;j<N;j++){
    if(rnd.uni()<connectivity) newEdge(i,j,E);
  }
}

template<class vert,class edge> void graphRandomLinear(MT::Array<vert*>& V,MT::Array<edge*>& E, uint N){
  uint i;
  for(i=0;i<N;i++) V.append(new vert);
  for(i=1;i<N;i++) newEdge(i-1,i,E);
}

template<class vert,class edge> void graphRandomTree(MT::Array<vert*>& V,MT::Array<edge*>& E, uint N, uint roots){
  uint i;
  for(i=0;i<N;i++) V.append(new vert);
  CHECK(roots>=1,"");
  for(i=roots;i<N;i++)  newEdge(rnd(i),i,E);
}

template<class vert,class edge>
void graphMaximumSpanningTree(MT::Array<vert*>& V,MT::Array<edge*>& E,const arr& W){
  CHECK(W.nd==2 && W.d0==W.d1,"");
  uint i;
  for(i=0;i<W.d0;i++) new_elem(V);

  boolA done(V.N);  done=false;
  uintA addedNodes;
  
  if(!V.N) return;

  i=rnd(V.N);
  addedNodes.append(i); done(i)=true;
  
  uint j,k;
  real Wmax;
  uintA m;
  while(addedNodes.N<V.N){
    m.clear();
    for(i=0;i<addedNodes.N;i++){
      j=addedNodes(i);
      for(k=0;k<V.N;k++) if(!done(k) && (!m.N || W(j,k)>Wmax)){ m=TUP(j,k); Wmax=W(j,k); }
    }
    CHECK(m.N,"graph is not connected!");
    if(done(m(1))) m.permute(0,1);
    newEdge(m(0),m(1),E);
    CHECK(!done(m(1)),"node added twice??");
    done(m(1))=true;
    addedNodes.append(m(1));
  }
  graphMakeLists(V,E);
}

/*template<class vert,class edge> void graphRandomFixedDegree(MT::Array<vert*>& V,MT::Array<edge*>& E,uint N,uint degree){
  uint i;
  for(i=0;i<N;i++) V.append(new vert);
  uintA D(N); D.setZero(); //degrees of vertices
  for(i=0;i<N;i++) for(j=i+1;j<N;j++){
    if(rnd.uni()<connectivity) newEdge(i,j,E);
  }
  for(i=0;i<I.N;i++) I(i)=i/degree; //a list, where each node appears degree-times
  I.permuteRandomly();
  cout <<I <<endl;
  for(i=0;i<I.N;i+=2) newEdge(I(i),I(i+1),E);
  graphConnectUndirected(V,E);
  graphWriteUndirected(cout,V,E);
}*/

template<class vert,class edge> void graphRandomFixedDegree(MT::Array<vert*>& V,MT::Array<edge*>& E,uint N,uint d){
  // --- from Joris' libDAI!!
  // Algorithm 1 in "Generating random regular graphs quickly"
  // by A. Steger and N.C. Wormald
  //
  // Draws a random graph with size N and uniform degree d
  // from an almost uniform probability distribution over these graphs
  // (which becomes uniform in the limit that d is small and N goes
  // to infinity).
    
  CHECK((N*d)%2==0,"");

  uint i,j;
  edge *e;
  for(i=0;i<N;i++) V.append(new vert);

  bool ready = false;
  uint tries = 0;
  while( !ready ) {
    tries++;

    // Start with N*d points {0,1,...,N*d-1} (N*d even) in N groups.
    // Put U = {0,1,...,N*d-1}. (U denotes the set of unpaired points.)
    uintA U;
    U.setStraightPerm(N*d);

    // Repeat the following until no suitable pair can be found: Choose
    // two random points i and j in U, and if they are suitable, pair
    // i with j and delete i and j from U.
    listDelete(E);
    bool finished = false;
    while(!finished){
      U.permuteRandomly();
      uint i1, i2;
      bool suit_pair_found = false;
      for(i1=0; i1<U.N-1 && !suit_pair_found; i1++){
        for(i2=i1+1; i2<U.N && !suit_pair_found; i2++){
          if((U(i1)/d) != (U(i2)/d)){  // they are suitable (refer to different nodes)
            suit_pair_found = true;
            newEdge( U(i1)/d, U(i2)/d, E);
            U.remove(i2);  // first remove largest
            U.remove(i1);  // remove smallest
          }
          if( !suit_pair_found || !U.N )  finished = true;
        }
      }
    }

    if(!U.N){
      // G is a graph with edge from vertex r to vertex s if and only if
      // there is a pair containing points in the r'th and s'th groups.
      // If G is d-regular, output, otherwise return to Step 1.
      uintA degrees(N);
      degrees.setZero();
      for_list(j,e,E){
        degrees(e->ifrom)++;
        degrees(e->ito)  ++;
      }
      ready = true;
      for(uint n=0; n<N; n++){
        CHECK(degrees(n)<=d,"");
        if(degrees(n)!=d){
          ready = false;
          break;
        }
      }
    }else ready=false;
  }

  graphConnectUndirected(V,E);
  //graphWriteUndirected(cout,V,E);
}



template<class vert,class edge> void graphLayered(MT::Array<vert*>& V,MT::Array<edge*>& E,const uintA& layers,bool interConnections){
  uint i,j,a=0,b=0;
  uint l,L=layers.N;
  for(l=0;l<L;l++){
    for(i=0;i<layers(l);i++){
      V.append(new vert);
      if(l) for(j=0;j<layers(l-1);j++) newEdge(b+j,a+i,E);
      if(l && interConnections) for(j=0;j<i;j++) newEdge(a+j,a+i,E);
    }
    a+=layers(l);             //a is the offset of current layer
    if(l) b+=layers(l-1);     //b is the offset of previous layer
  }
}

template<class vert,class edge> void graphMakeLists(MT::Array<vert*>& V,MT::Array<edge*>& E){
  vert *v;
  edge *e;
  uint i;
  for_list(i,v,V){
    v->outLinks.clear();
    v-> inLinks.clear();
  }
  for_list(i,e,E){
    e->from = V(e->ifrom);
    e->to   = V(e->ito);
    e->from->outLinks.append(e);
    e->to  -> inLinks.append(e);
  }
}

template<class vert,class edge> void graphConnectUndirected(MT::Array<vert*>& V,MT::Array<edge*>& E){
  vert *v;
  edge *e;
  uint i;
  for_list(i,v,V) v->edges.clear();
  for_list(i,e,E){
    e->from = V(e->ifrom);
    e->to   = V(e->ito);
    e->from->edges.append(e);
    e->to  ->edges.append(e);
  }
}

template<class vert,class edge> edge *newEdge(vert *a,vert *b,MT::Array<edge*>& E){
  return newEdge(a->index,b->index,E);
}

template<class edge> edge *newEdge(uint a,uint b,MT::Array<edge*>& E){
  edge *e=new edge;
  e->index=E.N;
  E.append(e);
  e->ifrom=a; e->ito=b;
  return e;
}

template<class vert,class edge> edge *new_edge_deprecated(vert *a,vert *b,MT::Array<edge*>& E){
  edge *e=new edge;
  e->index=E.N;
  E.append(e);
  e->a=a->index; e->b=b->index;
  return e;
}

template<class vert,class edge> edge *del_edge(edge *e,MT::Array<vert*>& V,MT::Array<edge*>& E,bool remakeLists){
  //uint i,k=e-E.p;
  //E.remove(k);
  uint i;
  for(i=0;i<E.N;i++) if(E(i)==e) break;
  E.remove(i);
  if(remakeLists) graphMakeLists(V,E);
  return e;
}


template<class vert,class edge> void graphWriteDirected(std::ostream& os,const MT::Array<vert*>& V,const MT::Array<edge*>& E){
  vert *v;
  edge *e;
  uint i,j;
  for_list(j,v,V){
    for_list(i,e,v->inLinks) os <<e->ifrom <<' ';
    os <<"-> ";
    os <<j <<" -> ";
    for_list(i,e,v->outLinks) os <<e->ito <<' ';
    os <<'\n';
  }
  for_list(j,e,E) os <<e->ifrom <<"->"<<e->ito <<'\n';
  //for_list(j,e,E) os <<e->from->name <<"->"<<e->to->name <<'\n';
}

template<class vert,class edge> void graphWriteUndirected(std::ostream& os,const MT::Array<vert*>& V,const MT::Array<edge*>& E){
  vert *v;
  edge *e;
  uint i,j;
  for_list(i,v,V){
    os <<i <<": ";
    for_list(j,e,v->edges) if(e->ifrom==i) os <<e->ito <<' '; else os <<e->ifrom <<' ';
    os <<'\n';
  }
  for_list(j,e,E) os <<e->ifrom <<"-"<<e->ito <<'\n';
}

template<class vert,class edge> bool graphTopsort(MT::Array<vert*>& V,MT::Array<edge*>& E){
  MT::Array<vert*> noInputs;
  noInputs.memMove=true;
  uintA newIndex(V.N);
  intA inputs(V.N);
  
  uint count=0;
  vert *v;
  edge *e;
  uint i;
  
  for_list(i,v,V){
    inputs(v->index)=v->inLinks.N;
    if(!inputs(v->index)) noInputs.append(v);
  }
  
  while(noInputs.N){
    v=noInputs.popFirst();
    newIndex(v->index)=count;
    count++;
    for_list(i,e,v->outLinks){
      inputs(e->to->index)--;
      if(!inputs(e->to->index)) noInputs.append(e->to);
    }
  }
  
  if(count!=V.N) return false;
  
  //success!
  //give each edge and vert new index:
  for_list(i,e,E){
    e->ifrom=newIndex(e->ifrom);
    e->ito  =newIndex(e->ito);
  }
  for_list(i,v,V){
    v->index=newIndex(v->index);
  }
  //permute vertex array:
  V.permuteInv(newIndex);
  graphMakeLists(V,E);
  
  return true;
}
         
template<class vert,class edge>
void graphRevertEdge(MT::Array<vert*>& V,MT::Array<edge*>& E, edge *e){
  uint i=e->ifrom;  e->ifrom=e->ito;  e->ito=i;
  graphMakeLists(V,E);
}

template<class vert,class edge,class CompareOp>
void maximumSpanningTree(MT::Array<vert*>& V,MT::Array<edge*>& E,const CompareOp& cmp){
  uint i;
  vert *n;
  edge *e,*m;  uint ei;
  boolA nodeAdded(V.N);  nodeAdded=false;
  boolA edgeAdded(E.N);  edgeAdded=false;
  uintA addedNodes;
  
  if(!V.N) return;

  n=V(rnd(V.N));
  addedNodes.append(n->index); nodeAdded(n->index)=true;
  while(addedNodes.N<V.N){
    m=0;
    for(i=0;i<addedNodes.N;i++){
      n=V(addedNodes(i));
      for_list(ei,e,n->outLinks) if(!nodeAdded(e->to  ->index) && (!m || cmp(e,m))) m=e;
      for_list(ei,e,n->inLinks ) if(!nodeAdded(e->from->index) && (!m || cmp(e,m))) m=e;
    }
    CHECK(m,"graph is not connected!");
    edgeAdded(m->index)=true;
    if(nodeAdded(m->to->index)) graphRevertEdge(V,E,m);
    CHECK(!nodeAdded(m->to->index),"node added twice??");
    nodeAdded(m->to->index)=true;
    addedNodes.append(m->to->index);
  }
  for_list_rev(ei,e,E) if(!edgeAdded(ei)) del_edge(e,V,E,true);
  graphMakeLists(V,E);
}
    
    
template<class T> MT::Array<T> get(const AnyList& L,const char* tag){
  uint i;
  Any *a;
  for_list(i,a,L) if(!strcmp(a->tag,tag)) break;
  if(i==L.N) HALT("tag '"<<tag<<"' is not in this AnyList");
  if(!a->n) return MT::Array<T>((T*)a->p,1);
  return MT::Array<T>((T*)a->p,a->n);
}


template<class T> T* anyListGet(const AnyList& L,const char *tag,uint n){
  uint i;
  for(i=0;i<L.N;i++){
    if(!strcmp(tag,L(i)->tag)){
      if(strcmp(typeid(T).name(),L(i)->type)){
        HALT("ABORT GETTING ATTRIBUTE -- found tag ("<<tag<<") but with different type ("<<typeid(T).name()<<") than requested ("<<L(i)->type<<")");
        return NULL;
      }
      if(n>1 && n!=L(i)->n){
        HALT("ABORT GETTING ATTRIBUTE -- found tag ("<<tag<<") but with different size ("<<L(i)->n<<") than requested ("<<n<<")");
        return NULL;
      }
      if(n==1 && !L(i)->p){
        HALT("ABORT GETTING ATTRIBUTE -- found tag ("<<tag<<") but as boolean flag instead of size 1)");
        return NULL;
      }
      if(n==0 && L(i)->p){
        HALT("ABORT GETTING ATTRIBUTE -- found tag ("<<tag<<") but with size ("<<L(i)->n<<") instead fo boolean flag)");
        return NULL;
      }
      if(!L(i)->p) return (T*)1; //boolean return - tag found but no data
      return (T*)(L(i)->p);
    }
  }
  //HALT("ABORT GETTING ATTRIBUTE -- couldn't find tag ("<<tag<<") in anylist");
  return NULL;
}



#endif
