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

/*! \file array.h
    \brief defines the MT::Array class */

#ifndef MT_array_h
#define MT_array_h

#include <iostream>

#define FOR1D(x,i)   for(i=0;i<x.d0;i++)
#define FOR1D_DOWN(x,i)   for(i=x.d0;i--;)
#define FOR2D(x,i,j) for(i=0;i<x.d0;i++) for(j=0;j<x.d1;j++)
#define FOR3D(x,i,j,k) for(i=0;i<x.d0;i++) for(j=0;j<x.d1;j++) for(k=0;k<x.d2;k++)
#define FOR_ALL(x,i)   for(i=0;i<x.N;i++)
#define forAll(i,A)  for(i=A.p;i!=A.pstop;i++)

#define for_index(i,X)  for(i=0;i<X.N;i++)
#define for_elem(e,X)   for(e=X.p;e!=X.pstop;e++)
#define for_list(i,e,X) for(i=0;i<X.N && ((e=X(i)) || true);i++)
#define for_list_rev(i,e,X) for(i=X.N;i-- && ((e=X(i)) || true);)

#ifdef MT_SINGLE
typedef float real;
#else
typedef double real;
#endif
typedef unsigned char byte;
typedef unsigned int uint;
typedef unsigned long long int uint64;

//-- global memory information and options
namespace MT{
  extern bool useLapack;
  extern uint64 memoryTotal, memoryBound;
  extern bool memoryStrict;
}

//===========================================================================
//
// Array class
//

namespace MT{
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

  template<class T> class Array{   //template <class T> class Array
  public:
    typedef Array<T> ArrayT; //typedef Array<T> ArrayT
    typedef bool (*ElemCompare)(const T& a,const T& b);

    //!@name data fields
    // @{
    T *p;     //!< the pointer on the linear memory allocated
    uint N;   //!< number of elements
    uint nd;  //!< number of dimensions
    uint d0;  //!< 0th dim
    uint d1;  //!< 1st dim
    uint d2;  //!< 2nd dim
    uint *d;  //!< pointer to dimensions (for nd<3 points to d0)
    T *pstop; //!< end of memory (pstop is already out of the bound)
    uint M;   //!< size of actually allocated memory (may be greater than N)
    bool reference;//!< true if this refers to some external memory
    T **pp;   //!< C-style 2D pointer (only valid if \c Array::getCarray was called)


    //options
    /*!\brief be careful!!!; this option will modify the \c
      Array::resizeCopy : instead of calling the copy operator= for
      each element, the memory is moved with the \c memmove command -
      works only if you know what you do ... [default is true for
      simple types (like double, int, etc) and false for unknown
      ones] */
    bool memMove;

    /*!\brief if flexiMem is true (which is default!) the resize method will
      (1) at the first call allocate the exact amount of memory, (2)
      at further calls of increasing memory allocate twice the memory
      needed or (3) at further calls of decreasing memory only free
      the memory if the new size is smaller than a fourth */
    bool flexiMem;

    //! stores the sizeof(T)
    static int sizeT;

#if 1//? garbage
    Array<uint> *sparse; //!< to allocate element, column and row index lists (1+d0+d1 lists)
    enum matrixType{ full, diag };
    matrixType mtype;
    T *p_device;
#endif
    // @}

    static char memMoveInit;
  private:
    void init();

  public:
    //!@name constructors
    Array();
    Array(const ArrayT& a);
    Array(uint i);
    Array(uint i,uint j);
    Array(uint i,uint j,uint k);
    Array(uint i,uint j,uint k,uint l);
    Array(const ArrayT& a,uint i);
    Array(const ArrayT& a,uint i,uint j);
    Array(const T* p,uint size);
    ~Array();

    ArrayT& operator=(const T& v);
    ArrayT& operator=(const ArrayT& a);

    //!@name resizing
    ArrayT& resize    (uint D0);
    ArrayT& resizeCopy(uint D0);
    ArrayT& reshape   (uint D0);
    ArrayT& resize    (uint D0,uint D1);
    ArrayT& resizeCopy(uint D0,uint D1);
    ArrayT& reshape   (uint D0,uint D1);
    ArrayT& resize    (uint D0,uint D1,uint D2);
    ArrayT& resizeCopy(uint D0,uint D1,uint D2);
    ArrayT& reshape   (uint D0,uint D1,uint D2);
    ArrayT& resize    (uint D0,uint D1,uint D2,uint D3);
    ArrayT& resize    (uint ND,uint *dim);
    ArrayT& resizeCopy(uint ND,uint *dim);
    ArrayT& reshape   (uint ND,uint *dim);
    ArrayT& resize    (const Array<uint> &newD);
    ArrayT& resizeCopy(const Array<uint> &newD);
    ArrayT& reshape   (const Array<uint> &newD);
    ArrayT& resizeAs    (const ArrayT& a);
    ArrayT& resizeCopyAs(const ArrayT& a);
    ArrayT& reshapeAs   (const ArrayT& a);

    //!@name initializing/setting entries
    void clear();
    void setText(const char* str);
    void setZero(byte zero=0);
    void setUni(const T& scalar,int d=-1);
    void setId(int d=-1);
    void setDiag(const T& scalar,int d=-1);
    void setDiag(const ArrayT& vector);
    void setBlockMatrix(const ArrayT& A,const ArrayT& B,const ArrayT& C,const ArrayT& D);
    void setBlockVector(const ArrayT& a,const ArrayT& b);
    void setMatrixBlock(const ArrayT& B,uint lo0,uint lo1);  
    void setVectorBlock(const ArrayT& B,uint lo);            
    void setStraightPerm(int n=-1);
    void setReversePerm(int n=-1);
    void setRandomPerm(int n=-1);
    void setCarray(T *buffer,uint D0);
    void setCarray(T **buffer,uint D0,uint D1);
    void referTo(const T *buffer,uint n);
    void referTo(const ArrayT& a);
    void referToSubRange(const ArrayT& a,uint i,int I);
    void referToSubDim(const ArrayT& a,uint dim);
    void referToSubDim(const ArrayT& a,uint i,uint j);
    //void redirect(const ArrayT& a,uint i); //? remove; hardly more efficient than referToSubDim
    void takeOver(ArrayT& a);                   //a becomes a reference on its previously owned memory!
    void setGrid(uint dim,T lo,T hi,uint steps);

    //!@name access
    uint dim(uint k) const;
    Array<uint> getDim() const{ Array<uint> dd; dd.referTo(d,nd); return dd; }
    T& elem(uint i) const;
    T& scalar() const;
    T& last() const;
    T& operator()(uint i) const;
    T& operator()(uint i,uint j) const;
    T& operator()(uint i,uint j,uint k) const;
    T& operator()(uint i,uint j,uint k,uint l) const;
    T& operator()(const Array<uint> &I) const;
    ArrayT operator[](uint i) const;
    ArrayT subDim(uint i,uint j) const;
    ArrayT& operator()();
    ArrayT sub(uint i,int I) const;
    ArrayT sub(uint i,int I,uint j,int J) const;
    ArrayT sub(uint i,int I,uint j,int J,uint k,int K) const;
    ArrayT sub(uint i,int I,Array<uint> cols) const;
    void getMatrixBlock(ArrayT& B,uint lo0,uint lo1) const;
    void getVectorBlock(ArrayT& B,uint lo) const;
    T** getCarray() const;
    void copyInto  (T *buffer) const;
    void copyInto2D(T **buffer) const;
    T& min() const;
    T& max() const;
    T absMax() const;
    void minmax(T& minVal, T& maxVal) const;
    uint minIndex() const;
    uint maxIndex() const;
    void maxIndeces(uint& m1,uint& m2) const; //best and 2nd best
    void maxIndex(uint& i,uint& j) const;
    void maxIndex(uint& i,uint& j,uint& k) const;
    int findValue(const T& x) const;
    void findValues(MT::Array<uint>& indices, const T& x) const;
    bool containsDoubles() const; // TL 17.07.08; whether at least one object is contained more than once
    uint getMemsize() const;
    void getIndexTuple(Array<uint> &I,uint i) const;

    //!@name appending etc
    T& append();
    void append(const T& x);
    void append(const ArrayT& x);
    void append(const T *p,uint n);
    void replicate(uint copies);
    void insert(uint i,const T& x);
    void replace(uint i,uint n,const ArrayT& x);
    void remove(uint i,uint n=1);
    void removePerm(uint i);          //more efficient for sets, works also for non-memMove arrays
    void removeValue(const T& x);
    bool removeValueSafe(const T& x); //? same as if((i=findValue(x))!=-1) remove[Perm](i);
    void removeAllValues(const T& x);
    void delRow(uint i);
    void delColumns(uint i,uint k=1);
    void insRows(uint i,uint k=1);
    void insColumns(uint i,uint k=1);
    void resizeDim(uint k,uint dk);
    void setAppend(const T& x); //? same as if(findValue(x)==-1) append(x)
    void setAppend(const ArrayT& x);
    T popFirst(){ T x; x=elem(0);   remove(0);   return x; }
    T popLast (){ T x; x=elem(N-1); remove(N-1); return x; }

    //!@name sorting etc
    uint findInSorted(const T& x,ElemCompare comp);
    uint insertSorted(const T& x,ElemCompare comp);
    void removeValueInSorted(const T& x,ElemCompare comp);

    //!@name apply permutations on this array
    void permute(uint i,uint j);
    void permute(const Array<uint>& permutation);
    void permuteInv(const Array<uint>& permutation);
    void permuteRandomly();
    void shiftPerm(int offset);

    //!@name sparse matrices [to be moved outside]
    real sparsity();
    void makeSparse();

    //!@name I/O
    void write(std::ostream& os=std::cout,const char *ELEMSEP=" ",const char *LINESEP="\n ",bool dimTag=false,bool binary=false) const;
    void read(std::istream& is);
    //void readOld(std::istream& is); //? garbage
    void writeDim(std::ostream& os=std::cout) const;
    void readDim(std::istream& is);
    void writeRaw(std::ostream& os) const;
    void readRaw(std::istream& is);
    void writeWithIndex(std::ostream& os=std::cout) const;
    const ArrayT& ioraw() const;
    void writeTagged(std::ostream& os,const char* tag,bool binary=false) const;
    bool readTagged (std::istream& is,const char *tag);
    
    //!@name kind of private
    void resizeMEM(uint n,bool copy);
    void freeMEM();
    void resetD();
  };
}


//===========================================================================
//
//!@name standard types
// @{

typedef MT::Array<real>   arr;
typedef MT::Array<double> doubleA;
typedef MT::Array<float>  floatA;
typedef MT::Array<uint>   uintA;
typedef MT::Array<int>    intA;
typedef MT::Array<char>   charA;
typedef MT::Array<byte>   byteA;
typedef MT::Array<bool>   boolA;
struct Any;
typedef MT::Array<Any*>   AnyList;


//===========================================================================
// @}
//!@name non-template functions //? needs most cleaning
// @{

arr Identity(uint dim);
void makeSymmetric(arr& A);
void transpose(arr& A);
void SUS(const arr& p,uint n,uintA& s);
uint SUS(const arr& p);

namespace MT{
  //! use this to turn on Lapack routines [default true if MT_LAPACK is defined]
  extern bool useLapack;
}

uint svd(const arr& A,arr& U,arr& d,arr& V,bool sort=true);
void svd(const arr& A,arr& U,arr& V);

uint inverse(arr& invA,const arr& A);
arr  inverse(const arr& A);
uint inverse_SVD(arr& inverse,const arr& A);
void inverse_LU(arr& invX,const arr& X);
void inverse_SymPosDef(arr& invA,const arr& A);
void pseudoInverse(arr& invA,const arr& A,const arr& invW,real robustnessEps);
void gaussFromData(arr& a,arr& A,const arr& X);
void rotationFromAtoB(arr& R,const arr& a,const arr& v);

real determinant(const arr& A);
real cofactor(const arr& A,uint i,uint j);

//void getIndexTuple(uintA &I,uint i,const uintA &d);  //? that also exists inside of array!
void lognormScale(arr& P,real& logP,bool force=true);

void gnuplot(const arr& X);
//void save(const MT::Array<const arr*>&, const char *filename, const char *comments);
void write(const arr& X,const arr& Y,const char* name);
void write(const arr& X,const arr& Y,const arr& Z,const char* name);
void write(const arr& X,const arr& Y,const arr& Z,const arr& W,const char* name);

void write_ppm(const byteA &img,const char *file_name,bool swap_rows=false);
void read_ppm(byteA &img,const char *file_name,bool swap_rows=false);
void add_alpha_channel(byteA &img,byte alpha);
void make_grey(byteA &img);
void make_RGB(byteA &img);
void flip_image(byteA &img);

void scanArrFile(const char* name);


//===========================================================================
// @}
//!@name template functions
// @{

template<class T> void transpose(MT::Array<T>& x,const MT::Array<T>& y);
template<class T> void getDiag(MT::Array<T>& x,const MT::Array<T>& y);
template<class T> MT::Array<T> diag(MT::Array<T>& x){  MT::Array<T> y;  getDiag(y,x);  return y;  }
template<class T> void inverse2d(MT::Array<T>& invA,const MT::Array<T>& A);

template<class T> T entropy(const MT::Array<T>& v);
template<class T> T normalizeDist(MT::Array<T>& v);
template<class T> void makeConditional(MT::Array<T>& P);
template<class T> void checkNormalization(MT::Array<T>& v,real tol);
template<class T> void checkNormalization(MT::Array<T>& v){ checkNormalization(v,1e-10); }
template<class T> void eliminate(MT::Array<T>& x,const MT::Array<T>& y,uint d);
template<class T> void eliminate(MT::Array<T>& x,const MT::Array<T>& y,uint d,uint e);
template<class T> void eliminatePartial(MT::Array<T>& x,const MT::Array<T>& y,uint d);

#ifndef SWIG
template<class T> T sqrDistance(const MT::Array<T>& v, const MT::Array<T>& w);
template<class T> T maxDiff(const MT::Array<T>& v,const MT::Array<T>& w,uint *maxi=0);
template<class T> T maxRelDiff(const MT::Array<T>& v,const MT::Array<T>& w,T tol);
//template<class T> T sqrDistance(const MT::Array<T>& v, const MT::Array<T>& w,const MT::Array<bool>& mask);
template<class T> T sqrDistance(const MT::Array<T>& g,const MT::Array<T>& v,const MT::Array<T>& w);
template<class T> T euclideanDistance(const MT::Array<T>& v, const MT::Array<T>& w);
template<class T> T metricDistance(const MT::Array<T>& g,const MT::Array<T>& v, const MT::Array<T>& w);

template<class T> T sum(const MT::Array<T>& v);
template<class T> T sumOfAbs(const MT::Array<T>& v);
template<class T> T sumOfSqr(const MT::Array<T>& v);
template<class T> T norm(const MT::Array<T>& v);
template<class T> T mean(const MT::Array<T>& v);
template<class T> T product(const MT::Array<T>& v);

template<class T> T trace(const MT::Array<T>& v);
template<class T> T var(const MT::Array<T>& v);
template<class T> T minDiag(const MT::Array<T>& v);

template<class T> void innerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z);
template<class T> void outerProduct(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z);
template<class T> T scalarProduct(const MT::Array<T>& v,const MT::Array<T>& w);
template<class T> T scalarProduct(const MT::Array<T>& g,const MT::Array<T>& v, const MT::Array<T>& w);

    //tensor operations [to be moved out of class declaration]
template<class T> void tensorCondNormalize(MT::Array<T> &X,int left);
template<class T> void tensorCondMax      (MT::Array<T> &X,uint left);
template<class T> void tensorCondSoftMax  (MT::Array<T> &X,uint left,real beta);
template<class T> void tensorCheckCondNormalization          (const MT::Array<T> &X,uint left,real tol=1e-10);
template<class T> void tensorCheckCondNormalization_with_logP(const MT::Array<T> &X,uint left, real logP, real tol=1e-10);

template<class T> void tensorEquation    (MT::Array<T> &X,const MT::Array<T> &A,const uintA &pickA,const MT::Array<T> &B,const uintA &pickB,uint sum=0);
template<class T> void tensorPermutation (MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid);
template<class T> void tensorMarginal    (MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid);
template<class T> void tensorMaxMarginal (MT::Array<T> &Y,const MT::Array<T> &X,const uintA &Yid);
template<class T> void tensorMarginal_old(MT::Array<T> &y,const MT::Array<T> &x,const uintA &xd,const uintA &ids);
template<class T> void tensorMultiply    (MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid);
template<class T> void tensorAdd         (MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid);
template<class T> void tensorMultiply_old(MT::Array<T> &x,const MT::Array<T> &y,const uintA &d,const uintA &ids);
template<class T> void tensorDivide      (MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid);
template<class T> void tensorAdd         (MT::Array<T> &X,const MT::Array<T> &Y,const uintA &Yid);

//arrays interpreted as a set (container)
template<class T> MT::Array<T> cat(const MT::Array<T>& y, const MT::Array<T>& z){ MT::Array<T> x; x.append(y); x.append(z); return x; }
template<class T> MT::Array<T> cat(const MT::Array<T>& y, const MT::Array<T>& z, const MT::Array<T>& w){ MT::Array<T> x; x.append(y); x.append(z); x.append(w); return x; }
template<class T> MT::Array<T> cat(const MT::Array<T>& a, const MT::Array<T>& b, const MT::Array<T>& c, const MT::Array<T>& d){ MT::Array<T> x; x.append(a); x.append(b); x.append(c); x.append(d); return x; }
template<class T> MT::Array<T> cat(const MT::Array<T>& a, const MT::Array<T>& b, const MT::Array<T>& c, const MT::Array<T>& d, const MT::Array<T>& e){ MT::Array<T> x; x.append(a); x.append(b); x.append(c); x.append(d); x.append(e); return x; }
template<class T> void setUnion(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z);
template<class T> void setSection(MT::Array<T>& x,const MT::Array<T>& y, const MT::Array<T>& z);
template<class T> MT::Array<T> setUnion(const MT::Array<T>& y, const MT::Array<T>& z){ MT::Array<T> x; setUnion(x,y,z); return x; }
template<class T> MT::Array<T> setSection(const MT::Array<T>& y, const MT::Array<T>& z){ MT::Array<T> x; setSection(x,y,z); return x; }
template<class T> void setMinus(MT::Array<T>& x,const MT::Array<T>& y);
template<class T> uint numberSharedElements(const MT::Array<T>& x,const MT::Array<T>& y);
template<class T> void rndInteger(MT::Array<T>& a,int low=0,int high=1,bool add=false);
template<class T> void rndUniform(MT::Array<T>& a,double low=0.,double high=1.,bool add=false);
template<class T> void rndNegLogUniform(MT::Array<T>& a,double low=0.,double high=1.,bool add=false);
template<class T> void rndGauss(MT::Array<T>& a,double stdDev=1.,bool add=false);
//template<class T> void rndGauss(MT::Array<T>& a,bool add=false);
//template<class T> MT::Array<T>& rndGauss(real stdDev,uint dim);
template<class T> uint softMax(const MT::Array<T>& a,arr& soft,real beta);
template<class T> MT::Array<T> sqr(const MT::Array<T>& y){ MT::Array<T> x; x.resizeAs(y); for(uint i=0;i<x.N;i++) x.elem(i)=y.elem(i)*y.elem(i); return x; }


//===========================================================================
// @}
//!@name macros to specify tuples and lists
// @{

template<class T> MT::Array<T> TUPLE(){                                    MT::Array<T> z(0); return z; }
template<class T> MT::Array<T> TUPLE(const T& i){                                    MT::Array<T> z(1); z(0)=i; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j){                               MT::Array<T> z(2); z(0)=i; z(1)=j; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k){                          MT::Array<T> z(3); z(0)=i; z(1)=j; z(2)=k; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k,const T& l){                     MT::Array<T> z(4); z(0)=i; z(1)=j; z(2)=k; z(3)=l; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k,const T& l,const T& m){                MT::Array<T> z(5); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n){           MT::Array<T> z(6); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n,const T& o){      MT::Array<T> z(7); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; z(6)=o; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n,const T& o,const T& p){ MT::Array<T> z(8); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; z(6)=o; z(7)=p; return z; }
template<class T> MT::Array<T> TUPLE(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n,const T& o,const T& p,const T& q){ MT::Array<T> z(9); z(0)=i; z(1)=j; z(2)=k; z(3)=l; z(4)=m; z(5)=n; z(6)=o; z(7)=p; z(8)=q; return z; }

#define TUP TUPLE<uint>
#define ARR TUPLE<real>

template<class T> MT::Array<T*> LIST(){                                    MT::Array<T*> z(0); return z; }
template<class T> MT::Array<T*> LIST(const T& i){                                    MT::Array<T*> z(1); z(0)=(T*)&i; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j){                               MT::Array<T*> z(2); z(0)=(T*)&i; z(1)=(T*)&j; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j,const T& k){                          MT::Array<T*> z(3); z(0)=(T*)&i; z(1)=(T*)&j; z(2)=(T*)&k; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j,const T& k,const T& l){                     MT::Array<T*> z(4); z(0)=(T*)&i; z(1)=(T*)&j; z(2)=(T*)&k; z(3)=(T*)&l; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j,const T& k,const T& l,const T& m){                MT::Array<T*> z(5); z(0)=(T*)&i; z(1)=(T*)&j; z(2)=(T*)&k; z(3)=(T*)&l; z(4)=(T*)&m; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n){           MT::Array<T*> z(6); z(0)=(T*)&i; z(1)=(T*)&j; z(2)=(T*)&k; z(3)=(T*)&l; z(4)=(T*)&m; z(5)=(T*)&n; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n,const T& o){      MT::Array<T*> z(7); z(0)=(T*)&i; z(1)=(T*)&j; z(2)=(T*)&k; z(3)=(T*)&l; z(4)=(T*)&m; z(5)=(T*)&n; z(6)=(T*)&o; return z; }
template<class T> MT::Array<T*> LIST(const T& i,const T& j,const T& k,const T& l,const T& m,const T& n,const T& o,const T& p){ MT::Array<T*> z(8); z(0)=(T*)&i; z(1)=(T*)&j; z(2)=(T*)&k; z(3)=(T*)&l; z(4)=(T*)&m; z(5)=(T*)&n; z(6)=(T*)&o; z(7)=(T*)&p; return z; }


//===========================================================================
// @}
//!@name basic operators
// @{

namespace MT{
template<class T> Array<T> operator~(const Array<T>& y);
template<class T> Array<T> operator-(const Array<T>& y);
template<class T> Array<T> operator^(const Array<T>& y,const Array<T>& z);
template<class T> Array<T> operator*(const Array<T>& y,const Array<T>& z);
template<class T> Array<T> operator*(const Array<T>& y,T z);
template<class T> Array<T> operator*(T y,const Array<T>& z);
  
#define BinaryOperator( op, name)					\
  template<class T> Array<T> operator op(const Array<T>& y,const Array<T>& z); \
  template<class T> Array<T> operator op(T y,const Array<T>& z);	\
  template<class T> Array<T> operator op(const Array<T>& y,T z)
  BinaryOperator( + , plusA);
  BinaryOperator( - , minusA);
  BinaryOperator( % , mult);
  BinaryOperator( / , div);
#undef BinaryOperator
  
#ifndef SWIG
#define CompoundAssignmentOperator( op )				\
  template<class T> Array<T>& operator op (Array<T>& x,const Array<T>& y); \
  template<class T> Array<T>& operator op ( Array<T>& x, T y )
  CompoundAssignmentOperator( |= );
  CompoundAssignmentOperator( ^= );
  CompoundAssignmentOperator( &= );
  CompoundAssignmentOperator( += );
  CompoundAssignmentOperator( -= );
  CompoundAssignmentOperator( *= );
  CompoundAssignmentOperator( /= );
  CompoundAssignmentOperator( %= );
#undef CompoundAssignmentOperator
#endif
}

#define UnaryOperation( name, op )					\
  template<class T> MT::Array<T>& name (MT::Array<T>& x,const MT::Array<T>& y); 
UnaryOperation( negative, - )
#undef UnaryOperation

#define BinaryOperation( name, op )					\
  template<class T> MT::Array<T>& name(MT::Array<T>& x,const MT::Array<T>& y,const MT::Array<T>& z); \
  template<class T> MT::Array<T>& name##S(MT::Array<T>& x,const MT::Array<T>& y,T z); \
  template<class T> MT::Array<T>& name##S(MT::Array<T>& x, T y, const MT::Array<T>& z)
BinaryOperation( plusA , + );
BinaryOperation( minusA , - );
BinaryOperation( mult , * );
BinaryOperation( div , / );
#undef BinaryOperation

#ifndef SWIG
#define UnaryFunction( func )						\
  template<class T> MT::Array<T> func (const MT::Array<T>& y)
UnaryFunction( acos );
UnaryFunction( asin );
UnaryFunction( atan );
UnaryFunction( cos );
UnaryFunction( sin );
UnaryFunction( tan );
UnaryFunction( cosh );
UnaryFunction( sinh );
UnaryFunction( tanh );
UnaryFunction( acosh );
UnaryFunction( asinh );
UnaryFunction( atanh );
UnaryFunction( exp );
UnaryFunction( log );
UnaryFunction( log10 );
UnaryFunction( sqrt );
UnaryFunction( cbrt );
UnaryFunction( ceil );
UnaryFunction( fabs );
UnaryFunction( floor );
#undef UnaryFunction

#define BinaryFunction( func )						\
  template<class T> MT::Array<T> func(const MT::Array<T>& y,const MT::Array<T>& z); \
  template<class T> MT::Array<T> func(const MT::Array<T>& y,T z);	\
template<class T> MT::Array<T> func(T y,const MT::Array<T>& z)
BinaryFunction( atan2 );
BinaryFunction( pow );
BinaryFunction( fmod );
#undef BinaryFunction

template<class T> std::istream& operator>>(std::istream& is,MT::Array<T>& x);
template<class T> MT::Array<T>& operator<<(MT::Array<T>& x,const char* str);
template<class T> std::ostream& operator<<(std::ostream& os,const MT::Array<T>& x);
template<class T> void checkNan(const MT::Array<T>& x);
template<class T> bool operator==(const MT::Array<T>& v,const MT::Array<T>& w);
template<class T> bool operator==(const MT::Array<T>& v,const T *w);
template<class T> bool operator!=(const MT::Array<T>& v,const MT::Array<T>& w);
template<class T> bool operator<(const MT::Array<T>& v,const MT::Array<T>& w);
#endif

//===========================================================================
// @}
//!@name double template functions
// @{

#ifndef SWIG
template<class T,class S> void resizeAs    (MT::Array<T>& x,const MT::Array<S>& a){
  x.nd=a.nd; x.d0=a.d0; x.d1=a.d1; x.d2=a.d2;
  x.resetD();
  if(x.nd>3){ x.d=new uint[x.nd];  memmove(x.d,a.d,x.nd*sizeof(uint)); }
  x.resizeMEM(a.N,false);
}
template<class T,class S> void resizeCopyAs(MT::Array<T>& x,const MT::Array<S>& a);
template<class T,class S> void reshapeAs   (MT::Array<T>& x,const MT::Array<S>& a);
template<class T,class S> void copy        (MT::Array<T>& x,const MT::Array<S>& a){
  resizeAs(x,a);
  for(uint i=0;i<x.N;i++) x.elem(i)=(T)a.elem(i);
}
//! check whether this and \c a have same dimensions
template<class T,class S>
bool samedim(const MT::Array<T>& a,const MT::Array<S>& b){
  return (b.nd==a.nd && b.d0==a.d0 && b.d1==a.d1 && b.d2==a.d2);
}
#endif


//===========================================================================
// @}
//!@name lapack interfaces
// @{

void blas_Mv(arr& y,const arr& A,const arr& x);
void blas_MM(arr& X,const arr& A,const arr& B);
void blas_MsymMsym(arr& X,const arr& A,const arr& B);
void lapack_cholesky(arr& C,const arr& A);
uint lapack_SVD(const arr& A, arr& U, arr& d, arr& Vt);
void lapack_LU(arr& LU, const arr& A);
void lapack_EigenDecomp(const arr& symmA, arr& Evals, arr& Evecs);
bool lapack_isPositiveSemiDefinite(const arr& symmA);
void lapack_inverseSymPosDef(arr& invA,const arr& A);
real lapack_determinantSymPosDef(const arr& A);
void lapack_Ainv_b_sym(arr& x,const arr& A, const arr& b);


//===========================================================================
// @}
//!@name lists
// @{

template<class T> MT::Array<T*> getList(const MT::Array<T>& A){
   MT::Array<T*> L;
   resizeAs(L,A);
   for(uint i=0;i<A.N;i++) L.elem(i) = &A.elem(i);
   return L;
}
    
template<class T> void listWrite(const MT::Array<T*>& L,std::ostream& os,const char *ELEMSEP=" ",const char *delim=NULL);
template<class T> void listWriteNames(const MT::Array<T*>& L,std::ostream& os);
/*template<class T> void listRead(MT::Array<T*>& L,std::istream& is);*/
template<class T> void listClone(MT::Array<T*>& L,const MT::Array<T*>& M);
template<class T> void listCopy(MT::Array<T*>& L,const MT::Array<T*>& M);
template<class T> void listDelete(MT::Array<T*>& L);
template<class T> T* listFindName(const MT::Array<T*>& L,const char* name);
template<class T> T* listFindType(const MT::Array<T*>& L,const char* type);

template<class T> T* new_elem(MT::Array<T*>& L){ T *e=new T; e->index=L.N; L.append(e); return e; }


// -- AnyLists
void anyListRead(AnyList& ats,std::istream& is);
template<class T> T* anyListGet(const AnyList& L,const char *tag,uint n);
template<class T> MT::Array<T> get(const AnyList& L,const char* tag);

//===========================================================================
// @}
//!@name graphs
// @{

template<class vert,class edge> edge* graphGetEdge(vert *from, vert *to);
template<class vert,class edge> void graphMakeLists(MT::Array<vert*>& V,MT::Array<edge*>& E);
template<class vert,class edge> void graphRandomUndirected(MT::Array<vert*>& V,MT::Array<edge*>& E,uint N,real connectivity);
template<class vert,class edge> void graphRandomFixedDegree(MT::Array<vert*>& V,MT::Array<edge*>& E,uint N,uint degree);
template<class vert,class edge> void graphConnectUndirected(MT::Array<vert*>& V,MT::Array<edge*>& E);
template<class vert,class edge> void graphLayered(MT::Array<vert*>& V,MT::Array<edge*>& E,const uintA& layers,bool interConnections);
template<class vert,class edge> edge *newEdge(vert *a,vert *b,MT::Array<edge*>& E);
template<class edge> edge *newEdge(uint a ,uint b,MT::Array<edge*>& E);
template<class vert,class edge> edge *del_edge(edge *e,MT::Array<vert*>& V,MT::Array<edge*>& E,bool remakeLists);
template<class vert,class edge> void graphWriteDirected(std::ostream& os,const MT::Array<vert*>& V,const MT::Array<edge*>& E);
template<class vert,class edge> void graphWriteUndirected(std::ostream& os,const MT::Array<vert*>& V,const MT::Array<edge*>& E);
template<class vert,class edge> bool graphTopsort(MT::Array<vert*>& V,MT::Array<edge*>& E);
template<class vert,class edge> void graphDelete(MT::Array<vert*>& V,MT::Array<edge*>& E);


//===========================================================================
// @}
// implementations
//

#if defined MT_IMPLEMENT_TEMPLATES | defined MT_IMPLEMENTATION
#  include "array_t.cpp"
#endif

#ifdef MT_IMPLEMENTATION
#  include"array.cpp"
#  include"array_lapack.cpp"
#endif
               
#endif
#endif

// (note: http://www.informit.com/articles/article.aspx?p=31783&seqNum=2)
        
