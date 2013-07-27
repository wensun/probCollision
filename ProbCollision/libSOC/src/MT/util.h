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

/*! \file util.h
    \brief many utilities: contains a bunch of standard includes, macros, operators,
    etc. that are used by most of my code */
#ifndef MT_util_h
#define MT_util_h

//----- if no system flag, I assume Linux
#if !defined MT_MSVC && !defined MT_Cygwin && !defined MT_Linux && !defined MT_MinGW
#  define MT_Linux
#endif

//----- standard includes:
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <math.h>
#include <string.h>
#include <cmath>
#include <stdint.h>
#if defined MT_Linux || defined MT_Cygwin
#  include<stdint.h>
#  include<sys/time.h>
#  include<sys/times.h>
#  include<sys/resource.h>
#endif
#if defined MT_MSVC
#  include<time.h>
#  include<sys/timeb.h> 
#  include<windows.h>
#  undef min
#  undef max
#  define MT_TIMEB
#  define __FUNCTION__ ""
#  ifdef MT_QT
#    undef  NOUNICODE
#    define NOUNICODE
#  endif
#  pragma warning(disable: 4305 4244 4250 4355 4786 4996)
#endif
#if defined MT_MinGW
#  include<unistd.h>
#  include<sys/time.h>
#  include<sys/timeb.h>
#  define MT_TIMEB
#endif


//----- basic defs:
#define MT_PI 3.14159265358979323846
#define MT_LN2 0.69314718055994528622676398299518041312694549560546875
#define MT_2PI 6.283195307179587
#define MT_LnSqrt2Pi -0.9189385332046727417803296 
#define MT_SQRT2 1.414213562373095049
#define MT_SQRTPI 1.772453850905516027
#ifdef MT_SINGLE
typedef float real;
#else
typedef double real;
#endif
typedef unsigned char byte;            //!< byte
typedef unsigned short int uint16;     //!< 2 bytes
typedef unsigned int uint;             //!< unsigned integer
typedef unsigned long long int uint64; //!< 8 bytes


//----- macros to define the standard << and >> operatos for most my classes:
#define stdInPipe(type)\
  inline std::istream& operator>>(std::istream& is,type& x){ x.read(is);return is; }
#define stdOutPipe(type)\
  inline std::ostream& operator<<(std::ostream& os,const type& x){ x.write(os); return os; }
#define stdPipes(type)\
  inline std::istream& operator>>(std::istream& is,type& x){ x.read(is);return is; }\
  inline std::ostream& operator<<(std::ostream& os,const type& x){ x.write(os); return os; }
#define inPipe(type)\
  inline type& operator<<(type& x,const char* str){ std::istringstream ss(str); ss >>x; return x; }



//----- macros for piping doubles EXACTLY (without rounding errors) in hex coding:
#define OUTHEX(y) "0x" <<std::hex <<*((unsigned long*)&y) <<std::dec
#define INHEX(y)  std::hex >>*((unsigned long*)&y) >>std::dec


//===========================================================================
//
// standard little methods in my namespace
//

namespace MT{
  extern int argc;
  extern char** argv;
  extern std::ifstream cfgFile;
  extern bool cfgOpenFlag,cfgLock;
  extern bool IOraw;  //!< stream modifier for some classes (Mem in particular)
  extern bool noLog;  //!< no logfile: default=true, becomes false when MT::init is called
  extern uint lineCount;

  //----- files
  void open(std::ofstream& fs,const char *name,const char *errmsg="");
  void open(std::ifstream& fs,const char *name,const char *errmsg="");
  template<class T> void save(const T& x,const char *filename);
  template<class T> void load(T& x,const char *filename,bool change_directory=false);
  void decomposeFilename(char *&_path,char *&name,const char *filename);
  std::ofstream& log(const char *name="MT.log");
  
  //----- strings and streams
  bool contains(const char *s,char c);
  char skip(std::istream& is,const char *skipchars=" \n\r\t",bool skipCommentLines=true);
  void skipLine(std::istream& is);
  char peerNextChar(std::istream& is,const char *skipchars=" \n\r\t",bool skipCommentLines=true);
  void parse(std::istream& is,const char *str);
  bool skipUntil(std::istream& is,const char *tag);

  //----- functions
  byte bit(byte *str,uint i);
  void flip(byte& b,uint i);
  void flip(int& b,uint i);
  real modMetric(real x,real y,real mod);
  real sign(real x);
  real linsig(real x);
  void   constrain(real& x,real a,real b);
  real phi(real dx,real dy);
  real dphi(real x,real y,real dx,real dy);
  real DIV(real x,real y,bool force=false);
  real sigmoid11(real x);
  real sigmoid(real x);
  real dsigmoid(real x);
  real approxExp(real x);
  real Log(real x);
  uint   Log2(uint n);
  real sqr(real x);
  real sinc(real x);
  real cosc(real x);
  template<class T> T MIN(T a,T b){ return a<b?a:b; }
  template<class T> T MAX(T a,T b){ return a>b?a:b; }
  real erf(real x);
  real gaussInt(real x);
  real gaussIntExpectation(real x);

  //----- time access
  real realTime();
  real cpuTime();
  real sysTime();
  real totalTime();
  char *date();
  void wait(real sec);
  bool wait();

  //----- memory
  long mem();

  //----- timer functions
  void   timerStart(bool useRealTime=false);
  real timerRead(bool reset=false);
  real timerRead(bool reset,real startTime);
  real timerPause();
  void   timerResume();

  //----- command line handling
  void initCmdLine(int _argc, char *_argv[]);
  bool checkCmdLineTag(const char *tag);
  char *getCmdLineArgument(const char *tag);

  //----- config file handling
  void openConfigFile(const char *name=0);

  //----- parameter grabbing from command line, config file, or default value
  template<class T> T getParameter(const char *tag);
  template<class T> T getParameter(const char *tag,const T& Default);
  template<class T> void getParameter(T& x,const char *tag,const T& Default);
  template<class T> void getParameter(T& x,const char *tag);
  template<class T> bool checkParameter(const char *tag);
}

std::istream& operator>>(std::istream& is,const char *str);
std::istream& operator>>(std::istream& is,char *str);


//===========================================================================
//
// String class
//

#define STRING(x) (((MT::String&)(MT::String().stream()<<x)).p)
#define STREAM(x) (((MT::String&)(MT::String().stream()<<x)).stream())

namespace MT{
  typedef std::iostream IOStream;
  
  /*!\brief String implements the functionalities of an ostream and an
  istream, but also can be send to an ostream or read from an
  istream. It is based on a simple streambuf derived from the
  MT::Mem class */
  class String:public IOStream{
    private:
      class StringBuf:public std::streambuf{
        public:
          String *string;
          virtual int overflow(int C = traits_type::eof());
          virtual int sync();
          void setIpos(char *p);
          char *getIpos();
      } buffer;
      void init();
      void append(char x);
      void resize(uint n,bool copy);
    
    public:
      //!@name data fields
      char *p;    //!< pointer to memory
      uint memN;  //!< \# elements (excluding zero)
      uint M;     //!< actual buffer size (in terms of # elements)
      static const char *readSkipSymbols; //!< default argument to read method (also called by operator>>)
      static const char *readStopSymbols; //!< default argument to read method (also called by operator>>)
      static int readEatStopSymbol;       //!< default argument to read method (also called by operator>>)
      void (*flushHandler)(String&);

      //!@name constructors
      String();
      String(const String& s);
      explicit String(const char *s);
      ~String();

      //!@name access
      uint N() const;
      operator char*();
      operator const char*() const;
      char &operator()(uint i) const;
      IOStream& stream();
      String& operator()();

      //!@name setting, appending
      String& operator=(const String& s);
      void operator=(const char *s);
      void set(const char *s,uint n);
      template<class T> String operator+(const T& v) const{ String news(*this); news <<v; return news; }
      template<class T> String prepend(const T& v) const{ String news; news <<v <<*this; return news; }
      
      //!@name resetting
      String& clr();
      String& resetI();

      //!@name equality
      bool operator==(const char *s);
      bool operator==(const String& s);

      //!@name I/O
      void write(std::ostream& os) const;
      void read(std::istream& is,const char* skipSymbols=NULL,const char *stopSymbols=NULL,int eatStopSymbol=-1);
  };
}
stdPipes(MT::String)


//===========================================================================
//
// macros for halting/MSGs etc
//


//----- declare my namespace for the first time:
//! Marc Toussaint namespace
namespace MT{
  extern String errString;

  inline void breakPoint(){
    int i;
    i=5;    //set a break point here, if you want to catch errors directly
  }
}

//----- error handling:
#ifdef MT_MSVC
#  define MT_HERE "@"<<(strrchr(__FILE__,'\\')?strrchr(__FILE__,'\\')+1:__FILE__)<<':'<<__LINE__<<':'<<__FUNCTION__<<": "
#else
#  define MT_HERE "@"<<(strrchr(__FILE__,'/')?strrchr(__FILE__,'/')+1:__FILE__)<<':'<<__LINE__<<':'<<__FUNCTION__<<": "
#endif
#ifndef MT_MSG
#  define MT_MSG(msg){ std::cerr <<MT_HERE <<msg <<std::endl; MT::breakPoint(); }
#endif
#ifndef HALT
#  define HALT(msg)  { MT::errString.clr() <<MT_HERE <<msg <<" --- HALT"; std::cerr <<MT::errString <<std::endl; MT::breakPoint(); throw MT::errString.p; }
#  define NIY HALT("not implemented yet")
#  define OPS HALT("obsolete")
#endif


//----- check macros:
#ifndef MT_NOCHECK
#  define MT_DEBUG(x) x
#  define CHECK(cond,msg) if(!(cond)) HALT("CHECK failed: "<<msg);
#else
#  define MT_DEBUG(x)
#  define CHECK(cond,msg)
#endif


//----- other macros:
#define MEM_COPY_OPERATOR(x) memmove(this,&x,sizeof(this));

//===========================================================================
//
// a shared memory
//

namespace MT{

#define maxShmGuests 20
#define shmMaxBlocks 20
  struct SHMinfoBlock{
    const char pagename[20];
    uint size,used;
    uint guests;
    int guestPids[maxShmGuests];
    char blockNames[shmMaxBlocks][20];
    uint blockOffsets[shmMaxBlocks];
    pthread_rwlock_t blockLocks[shmMaxBlocks];
  };

  struct Lock{
  //documentation: http://cs.pub.ro/~apc/2003/resources/pthreads/uguide/document.htm
    pthread_rwlock_t lock;
    void init(bool process_shared=true){
      pthread_rwlockattr_t   att;
      int r;
      r=pthread_rwlockattr_init(&att); if(r) HALT("");
      if(process_shared){
        r=pthread_rwlockattr_setpshared(&att, PTHREAD_PROCESS_SHARED); if(r) HALT("");
      }
      r=pthread_rwlock_init(&lock, &att); if(r) HALT("");
    }
    void read(){
      int r=pthread_rwlock_rdlock(&lock);
      if(r) HALT("");
    }
  
    void write(){
      int r=pthread_rwlock_wrlock(&lock);
      if(r) HALT("");
    }
  
    void unlock(){
      int r=pthread_rwlock_unlock(&lock);
      if(r) HALT("");
    }
  };

  struct SHM{
    int systemId;
    uint guestId;
    byte *p;
    bool created,opened;
    SHMinfoBlock *info;
    SHM(){ created=opened=false; p=NULL; systemId=0; info=NULL; guestId=-1; }
    void open(const char *pagename,uint size,bool directDestroy=false);
    void close();
    void destroy();
    void readLock(uint i);
    void writeLock(uint i);
    void unlock(uint i);
    void report();
    template<class T> T* newBlock(const char* objName,uint *i=NULL);
  };

  template<class T> struct Shared{
    uint blockId;
    T *p;
    SHM *shm;
    Shared(SHM& _shm,const char* name){ shm=&_shm; p=shm->newBlock<T>(name,&blockId); }
    void readLock() { shm->readLock(blockId); }
    void writeLock(){ shm->writeLock(blockId); }
    void unlock()   { shm->unlock(blockId); }
    void set(const T& x){ writeLock(); *p=x; unlock(); }
    void get(T& x)      { readLock();  x=*p; unlock(); }
    T& operator()(){ return *p; }
  };
}

//===========================================================================
//
// Parameter class - I use it frequently to read parameters from file or cmd line
// 

#if 1
namespace MT{
  /*!\brief A parameter that initializes itself from the command line
    (use \c MT::init), parameter file, or a default value (priority in
    this order).  Initialization is done on the fly the _first_ time
    its value is queried (i.e., referenced by the cast operators).*/
  template<class type>
  class Parameter{
  public:
    const char *typeName;
    type value,Default;
    const char *tag;
    bool initialized,hasDefault;

  public:
    //!@name constructors

    //! Determines the tag to search for in parameter file/command line
    explicit Parameter(const char *_tag){
      typeName=typeid(type).name();
      initialized=false;
      tag=_tag;
      hasDefault=false;
    };

    /*!\brief specifies also a default value -- parameter does not have to but
      can be specified in the parameter file/command line */
    Parameter(const char *_tag,const type& _default){
      typeName=typeid(type).name();
      initialized=false;
      tag=_tag;
      hasDefault=true;
      Default=_default;
    };

    ~Parameter(){}

    //!@name value access

    //! standard type conversion: returns a const of the parameter value
    operator type(){ if(!initialized) initialize(); return value; }

    //! ()-operator: returns an lvalue of the parameter value
    type& operator()(){ if(!initialized) initialize(); return value; }


    //!@name manipulation

    //! assigs a value to the parameter -- no further initialization needed
    type& operator=(const type v){ initialized=true; value=v; return value; }

    //! set the tag (replacing the one from the constructor)
    void setTag(char *_tag){ tag=_tag; }

    /*!\brief enforces that the parameter is reinitialized from the parameter
      file/command line, the next time it is referenced -- even if it
      has been initialized before */
    void reInitialize(){ initialized=false; }


    //!@name explicit grabbing


  private:
    void initialize();
  };

}
#endif


//===========================================================================
//
// random number generator
//

namespace MT{
  /*!\brief A random number generator. An global instantiation \c
    MT::rnd of a \c Rnd object is created. Use this one object to get
    random numbers.*/
  class Rnd{
  private:
    bool ready;
    int32_t rpoint;     /* Feldindex    */
    int32_t rfield[256];   /* Schieberegisterfeld  */


  public:
    //! ...
    Rnd(){ ready=false; };


  public://!@name initialization
    //! initialize with a specific seed
    uint32_t seed(uint32_t n);

    //! use Parameter<uint>("seed") as seed
    uint32_t seed();

    //! uses the internal clock to generate a seed
    uint32_t clockSeed();
    
  public://!@name access
    //! a initeger random number uniformly distributed in [0,?]
    uint32_t num(){ if(!ready) seed(); return (uint32_t)rnd250() >> 5; }
    //! same as \c num()
    uint32_t operator ()(){ return num(); }
    //! a initeger random number uniformly distributed in [0,\c i-1]
    uint32_t num(uint32_t limit){
      CHECK(limit,"zero limit in rnd.num()"); return num() % limit; }
    uint32_t num(int32_t lo, int32_t hi){ return lo+num(hi-lo+1); }
    //! same as \c num(i)
    uint32_t operator ()(uint32_t i){ return num(i); }
    uint32_t operator ()(int32_t lo, int32_t hi){ return num(lo,hi); }
    //! a random variable uniformly distributed in [0,1]
    real uni(){ return ((real)num(1<<22))/(1<<22); }
    //! a random variable uniformly distributed in [\c low, \c high]
    real uni(real low,real high){ return low+uni()*(high-low); }
    //! a gaussian random variable with mean zero
    real gauss();
    /*!\brief a positive integer drawn from a poisson distribution with given
      \c mean; is case \c mean>100, a (positive) gauss number
      \c floor(mean+gauss(sqrt(mean))+.5) is returned */
    uint32_t poisson(real mean);
    // @}


  private:
    int32_t rnd250(){
      rpoint = (rpoint+1) & 255;          // Index erhoehen
      return rfield[rpoint] =  rfield[(rpoint-250) & 255]
	                     ^ rfield[(rpoint-103) & 255];
    }

    void seed250(int32_t seed);
  };

  //! The global Rnd object
  extern Rnd rnd;
}


//===========================================================================
//
// generic any container
//

struct Any{
  char *tag;
  const char *type;
  void *p;
  uint n;
  char delim;
  virtual ~Any(){};
  //virtual void read(std::ostream &os) = 0;
  virtual void write(std::ostream &os) const = 0;
  virtual void free() = 0;
  virtual Any *new_clone() = 0;
  //virtual void operator=(const Any&) = 0;
  //virtual void set(const void *_p,uint _n) = 0;
};
stdOutPipe(Any)

template<class T> Any* anyNew(const char* tag,const T &x);
template<class T> Any* anyNew(const char* tag,const T *x,uint n,char delim);


//===========================================================================
//
// gnuplot calls
//

void gnuplot(const char *command,const char* EPSfile=NULL);
void gnuplotClose();


//===========================================================================
//
// implementations
//

#if defined MT_IMPLEMENTATION | defined MT_IMPLEMENT_TEMPLATES
#  include "util_t.cpp"
#endif

#ifdef MT_IMPLEMENTATION
#  include"util.cpp"
#endif



//===========================================================================
//
// USING's
//

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using MT::rnd;
using MT::String;


#endif
