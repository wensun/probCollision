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

#include "util.h"

#ifdef MT_QT
#  undef scroll
#  undef Unsorted
#  include <Qt/qmetatype.h>
#  include <Qt/qdatastream.h>
#  include <Qt/qapplication.h>
#endif

#ifndef MT_ConfigFileName
#  define MT_ConfigFileName "MT.cfg"
#  define MT_LogFileName "MT.log"
#endif

#ifndef MT_MSVC
#  include <sys/ipc.h>
#  include <sys/shm.h>
#endif

//===========================================================================
//
// Bag container
//

const char *MT::String::readSkipSymbols = " \t";
const char *MT::String::readStopSymbols = "\n\r";
int   MT::String::readEatStopSymbol     = 1;
MT::String MT::errString;


//===========================================================================
//
// utilities in MT namespace
//

namespace MT{
  int argc;
  char** argv;
  std::ifstream cfgFile;
  bool cfgOpenFlag=false;
  bool cfgLock=false;
  bool IOraw=false;
  bool noLog=true;
  uint lineCount=0;
  
#ifndef MT_TIMEB
  timeval startTime;
#else
  _timeb startTime;
#endif
  real timerStartTime=0.;
  real timerPauseTime=-1.;
  bool timerUseRealTime=false;

#ifdef MT_QT
  QApplication *myApp=NULL;
#endif

  class Demon{
  public:
    std::ofstream logFile;
    int logstat;

    Demon(){
      logstat=0;
      timerStartTime=MT::cpuTime();
#ifndef MT_TIMEB
      gettimeofday(&startTime,0);
#else
      _ftime(&startTime);
#endif
    }
    ~Demon(){
      if(logstat){ //don't open a log file anymore in the destructor
        char times[200];
        sprintf(times,"Ellapsed real time:  %.3lfsec\nProcess  user time:   %.3lfsec", realTime(),cpuTime());
        MT::log() <<"Execution stop:      " <<date()
          <<times <<std::endl;
        //"Ellapsed real time:  " <<::dtoa(realTime()) <<"sec\n"
          //<<"Process  user time:   "<<::dtoa(cpuTime()) <<"sec" <<std::endl;
#ifndef MT_TIMEB
        MT::log() <<"Process system time: " <<sysTime() <<"sec" <<std::endl;
#endif
        MT::log().close();
      }
#ifdef MT_QT
      if(myApp){
        myApp->processEvents();
        myApp->quit();
      }
#endif
    }

    std::ofstream& log(const char *name){
      if(!logstat && noLog) logstat=1;
      if(!logstat){
	logFile.open(name);
	if(!logFile.good()) MT_MSG("could not open log-file `"<<name<<"' for output");
	logstat=1;
      }
      return logFile;
    }
  }demon;

  //! access to the log-file
  std::ofstream& log(const char *name){ return demon.log(name); }

  //! open an output-file with name '\c name'
  void open(std::ofstream& fs,const char *name,const char *errmsg){
    fs.clear();
    fs.open(name);
    log() <<"opening output file `" <<name <<"'" <<std::endl;
    if(!fs.good()) MT_MSG("could not open file `" <<name <<"' for output"<<errmsg);
  }

  //! open an input-file with name '\c name'
  void open(std::ifstream& fs,const char *name,const char *errmsg){ 
    fs.clear();
    fs.open(name);
    log() <<"opening input file `" <<name <<"'" <<std::endl;
    if(!fs.good()) HALT("could not open file `" <<name <<"' for input"<<errmsg);
  }

  //! change to the directory of the given filename
  void decomposeFilename(char *&_path,char *&name,const char *filename){
    static char path[128];
    uint i=strlen(filename);
    CHECK(i<128,"");
    memmove(path,filename,i);
    for(;i--;) if(path[i]=='/' || path[i]=='\\') break;
    if(!(i+1)) path[0]=0; else path[i]=0;
    name = (char*)filename+i+1;
    _path = path;
  }


  //! returns true if the (0-terminated) string s contains c 
  bool contains(const char *s,char c){
    for(uint i=0;s[i];i++) if(s[i]==c) return true;
    return false;
  }

  //! skips the chars (typically white characters) when parsing from the istream
  char skip(std::istream& is,const char *skipchars,bool skipCommentLines){
    char c;
    for(;;){
      c=is.get();
      if(skipCommentLines && c=='#'){ skipLine(is); continue; }
      if(c=='\n') lineCount++;
      if(contains(skipchars,c)) continue;
      break;
    }
    is.putback(c);
    return c;
  }

  //! skips a newline character (same as skip(is,"\n");)
  void skipLine(std::istream& is){
    char c;
    do{ c=is.get(); } while(c!='\n');
    lineCount++;
  }

  //! tell you about the next char (after skip()) but puts it back in the stream
  char peerNextChar(std::istream& is,const char *skipchars,bool skipCommentLines){
    char c;
    skip(is,skipchars,skipCommentLines);
    is.get(c);
    is.putback(c);
    return c;
  }

  //! skip throught a stream until the tag is found (which is eaten)
  bool skipUntil(std::istream& is,const char *tag){
    unsigned n=strlen(tag);
    char *buf=new char [n+1];
    memset(buf,0,n+1);
    while(is.good()){
      memmove(buf,buf+1,n);
      buf[n-1]=is.get();
      if(buf[n-1]=='\n') lineCount++;
      buf[n]=0;
      if(!strcmp(tag,buf)){ delete[] buf; return true; }
    };
    delete[] buf;
    return false;
  }

  //! a global operator to scan (parse) strings from a stream
  void parse(std::istream& is,const char *str){
    if(!is.good()){ MT_MSG("bad stream tag when scanning for `" <<str <<"'"); return; }//is.clear(); }
    uint i,n=strlen(str);
    char *buf=new char [n+1]; buf[n]=0;
    MT::skip(is," \n\r\t");
    is.read(buf,n);
    if(!is.good() || strcmp(str,buf)){
      for(i=n;i--;) is.putback(buf[i]);
      is.setstate(std::ios::failbit);
      MT_MSG("scanning of constant string `" <<str
          <<"' failed! (read instead: `" <<buf <<"')");
    }
    delete[] buf;
  }


  //! returns the i-th of str 
  byte bit(byte *str,uint i){ return (str[i>>3] >> (7-(i&7))) & 1; }
  //byte bit(byte b,uint i){ return (b >> (7-(i&7))) & 1; }
  //void set(byte *state,uint i){ state[i>>3] |= 1 << (7-(i&7)); }
  //void del(byte *state,uint i){ state[i>>3] &= (byte)~(1 << (7-(i&7))); }
  //void flip(byte *str,uint i){ str[i>>3] ^= 1 << (7-(i&7)); }

  //! flips the i-th bit of b
  void flip(byte& b,uint i){ b ^= 1 << (7-(i&7)); }

  //! filps the i-th bit of b
  void flip(int& b,uint i){ b ^= 1 << (7-(i&7)); }

  /*!\brief the distance between x and y w.r.t.\ a circular topology
      (e.g. modMetric(1,8,10)=3) */
  real modMetric(real x,real y,real mod){
    real d=fabs(x-y);
    d=fmod(d,mod);
    if(d>mod/2.) d=mod-d;
    return d;
  }

  //! the sign (+/-1) of x (+1 for zero)
  real sign(real x){ if(x<0.) return -1.; return 1.; }

  //! returns 0 for x<0, 1 for x>1, x for 0<x<1
  real linsig(real x){ if(x<0.) return 0.; if(x>1.) return 1.; return x; }

  //! x ends up in the interval [a,b]
  void constrain(real& x,real a,real b){ if(x<a) x=a; if(x>b) x=b; }

  //! the angle of the vector (x,y) in [-pi,pi]
  real phi(real x,real y){
    if(x==0. || ::fabs(x)<1e-10){ if(y>0.) return MT_PI/2.; else return -MT_PI/2.; }
    real p=::atan(y/x);
    if(x<0.){ if(y<0.) p-=MT_PI; else p+=MT_PI; }
    if(p>MT_PI)  p-=2.*MT_PI;
    if(p<-MT_PI) p+=2.*MT_PI;
    return p;
  }

  //! the change of angle of the vector (x,y) when displaced by (dx,dy)
  real dphi(real x,real y,real dx,real dy){
    //return (dy*x - dx*y)/sqrt(x*x+y*y);
    if(x==0. || ::fabs(x)<1e-10){ if(y>0.) return -dx/y; else return dx/y; }
    real f=y/x;
    return 1./(1.+f*f)*(dy/x - f/x*dx);
  }

  /*!\brief save division, checks for division by zero; force=true will return
    zero if y=0 */
  real DIV(real x,real y,bool force){
    if(x==0.) return 0.;
    if(force){ if(y==0.) return 0.; }else CHECK(y!=0,"Division by Zero!");
    return x/y;
  }

  real sigmoid11(real x){
    return x/(1.+::fabs(x));
  }

  real sigmoid(real x){
    return 1./(1.+exp(-x));
  }

  real dsigmoid(real x){
    real y=sigmoid(x);
    return y*(1.-y);
  }

#define AXETS 1280
#define AXETR 10.
  //! approximate exp (sets up a static value table)
  real approxExp(real x){
    static bool initialized=false;
    static real ExpTable [AXETS]; //table ranges from x=-10 to x=10
    int i;
    if(!initialized){
      for(i=0;i<AXETS;i++) ExpTable[i]=::exp(AXETR*(2*i-AXETS)/AXETS);
      initialized=true;
    }
    x*=.5*AXETS/AXETR;
    i=(int)x;
    x-=(real)i; //x = residual
    i+=AXETS/2; //zero offset
    if(i>=AXETS-1) return ExpTable[AXETS-1];
    if(i<=0) return 0.;//ExpTable[0];
    return (1.-x)*ExpTable[i] + x*ExpTable[i+1];
  }

  //! ordinary Log, but cutting off for small values
  real Log(real x){
    if(x<.001) x=.001; return ::log(x);
  }

  //! integer log2
  uint Log2(uint n){
    uint l=0;
    n=n>>1;
    while(n){ l++; n=n>>1; }
    return l;
  }

  //! square of a real
  real sqr(real x){ return x*x; }

  real sinc(real x){
    if(fabs(x)<1e-10) return 1.-.167*x*x;
    return ::sin(x)/x;
  }

  real cosc(real x){
    if(fabs(x)<1e-10) return 1.-.167*x*x;
    return ::cos(x)/x;
  }




  /*!\brief real time since start of the process in floating-point seconds
    (probably in micro second resolution) -- Windows checked! */
  real realTime(){
#ifndef MT_TIMEB
    static timeval t; gettimeofday(&t,0);
    return ((real)(t.tv_sec-startTime.tv_sec-1) +
	    (real)((long)1000000+t.tv_usec-startTime.tv_usec)/1000000.);
#else
    static _timeb t; _ftime(&t);
    return ((real)(t.time-startTime.time-1) +
	    (real)((unsigned short)1000+t.millitm-startTime.millitm)/1000.);
#endif
  }

  /*!\brief user CPU time of this process in floating-point seconds (pure
    processor time) -- Windows checked! */
  real cpuTime(){
#ifndef MT_TIMEB
    static tms t; times(&t);
    return ((real)t.tms_utime)/sysconf(_SC_CLK_TCK);
#else
    return ((real)clock())/CLOCKS_PER_SEC; //MSVC: CLOCKS_PER_SEC=1000
#endif
  }

  /*!\brief system CPU time of this process in floating-point seconds (the
    time spend for file input/output, x-server stuff, etc.)
    -- not implemented for Windows! */
  real sysTime(){
#ifndef MT_TIMEB
    static tms t; times(&t);
    return ((real)(t.tms_stime))/sysconf(_SC_CLK_TCK);
#else
    HALT("sysTime() is not implemented for Windows!");
    return 0.;
#endif
  }

  /*!\brief total CPU time of this process in floating-point seconds (same
    as cpuTime + sysTime) -- not implemented for Windows! */
  real totalTime(){
#ifndef MT_TIMEB
    static tms t; times(&t);
    return ((real)(t.tms_utime+t.tms_stime))/sysconf(_SC_CLK_TCK);
#else
    HALT("totalTime() is not implemented for Windows!");
    return 0.;
#endif
  }

  //! the absolute real time and date as string
  char *date(){ static time_t t; time(&t); return ctime(&t); }

  //void wait(real sec){ real t=realTime(); while(realTime()-t<sec); }

  //! wait real time
  void wait(real sec){
    timespec ts;
    ts.tv_sec = (long)(floor(sec));
    sec -= (real)ts.tv_sec;
    ts.tv_nsec = (long)(floor(1000000000. * sec));
    int rc = clock_nanosleep(CLOCK_MONOTONIC,0,&ts,NULL);
    if(rc){
      if(rc==0){ MT_MSG("clock_nanosleep() interrupted by signal"); }
      else MT_MSG("clock_nanosleep() failed " <<rc);
    }

#if 0    
#ifndef MT_TIMEB
    timeval tv;
    tv.tv_sec = (int)(floor(sec));
    sec -= (real)tv.tv_sec;
    tv.tv_usec = (int)(floor(1000000. * sec));
    int r = select(1, NULL, NULL, NULL, &tv);
    if(r==-1) MT_MSG("select() failed");
    /* r=0 time is up
       r!=0 data in NULL stream available (nonsense)
       */
#else
#  ifdef _WINDOWS_
    Sleep((int)(1000.*sec));
    //MsgWaitForMultipleObjects( 0, NULL, FALSE, (int)(1000.*sec), QS_ALLEVENTS);
#  else
    real t=realTime(); while(realTime()-t<sec);
#  endif
#endif
#endif
  }

  //! wait for an ENTER at the console
  bool wait(){
    char c[10];
    std::cout <<" -- hit a key to continue..." <<std::flush;
    //cbreak(); getch();
    std::cin.getline(c,10);
    std::cout <<"\r" <<std::flush;
    if(c[0]==' ') return true;
    else return false;
    return true;
  }


  //! the integral shared memory size -- not implemented for Windows!
  long mem(){
#ifndef MT_TIMEB
    static rusage r; getrusage(RUSAGE_SELF,&r);
    return r.ru_idrss;
#else
    HALT("MT::mem() is not implemented for Windows!");
    return 0;
#endif
  }

#ifdef MT_MSVC
  struct SHM{
    HANDLE systemId;
    void *p;
    bool created,opened;
    SHM(){ created=opened=false; p=NULL; systemId=NULL; }
  };
  
  void *openSharedMemory(SHM& shm,const char *name,uint size){
    CHECK(!created && !opened,"shm object already in use");
    
    systemId = OpenFileMapping(FILE_MAP_ALL_ACCESS, FALSE, name);
    if(systemId == NULL){
      systemId = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, size, name);
      std::cout <<"** creating new shared memory `" <<name<<"' with systemId " <<systemId <<std::endl;
      if(systemId==NULL)
        HALT("can't create shared memory `"<<name<<"'. CreateFileMapping error code "<<GetLastError());
      created=true;
    }else{
      std::cout <<"** opening existing shared memory `" <<name<<"' with systemId " <<systemId <<std::endl;
    }
    
    p = MapViewOfFile(systemId, FILE_MAP_ALL_ACCESS, 0, 0, size);
    if(p==NULL)
      HALT("Could not map view of file. MapViewOfFile error code "<<GetLastError());
    
    opened=true;
    return p;
  }
  
  void closeSharedMemory(SHM& shm){
    CHECK(opened,"can't close what's not open");
    UnmapViewOfFile(p);
    CloseHandle(systemId);
    systemId=NULL;
    p=NULL;
    created=opened=false;
  }

  void destroySharedMemory(SHM& shm){
    std::cout <<"** destroying shared memory with systemId "<<systemId <<" -- not supported under MSVC" <<std::endl;
  }

#else
#ifdef MT_PTHREAD
  void SHM::open(const char *name,uint size,bool directDestroy){
    key_t key=ftok(name,0);
    systemId=shmget(key, size,0666);
    if(systemId==-1){
      systemId=shmget(key, size, IPC_CREAT|0666);
      std::cout <<"** creating new shared memory `" <<name<<"' with systemId " <<systemId <<std::endl;
      if(systemId==-1)
        HALT("can't create shared memory `"<<name);
      created=true;
    }else{
      std::cout <<"** opening existing shared memory `" <<name<<"' with systemId " <<systemId <<std::endl;
    }
    opened=true;
    
    if(directDestroy){ destroy(); return; }

    p=(byte*)shmat(systemId, NULL, 0);
    if(p==(byte*) -1)
      HALT("can't access memory of shared memory `"<<name);

    info=(SHMinfoBlock*)p;
    uint i;
    if(created){
      if(size<=sizeof(SHMinfoBlock)) HALT("need to allocate larger SHM to hold info block!");
      strcpy((char*)info->pagename,name);
      info->size=size;
      info->used=sizeof(SHMinfoBlock);
      info->guests=1;
      guestId=0;
      for(i=0;i<maxShmGuests;i++) info->guestPids[i]=-1;
      info->guestPids[guestId]=getpid();
      for(i=0;i<shmMaxBlocks;i++) info->blockNames[i][0]=0;
      for(i=0;i<shmMaxBlocks;i++) info->blockOffsets[i]=0;
      //create mutexes for all (potential) blockes:
      pthread_rwlockattr_t   att;
      int r;
      r=pthread_rwlockattr_init(&att); CHECK(!r,"");
      r=pthread_rwlockattr_setpshared(&att, PTHREAD_PROCESS_SHARED); CHECK(!r,"");
      for(i=0;i<shmMaxBlocks;i++){
        r=pthread_rwlock_init(&info->blockLocks[i], &att);
        CHECK(!r,"");
      }
    }else{
      //CHECK(!strcmp(name,info->pagename),"requested shm has different name");
      //CHECK(size==info->size,"requested shm has different size");
      info->guests++;
      for(i=0;i<maxShmGuests;i++) if(info->guestPids[i]==-1) break;
      if(i==maxShmGuests) HALT("exceeded max number of clients "<<maxShmGuests);
      guestId=i;
      info->guestPids[guestId]=getpid();
    }
    opened=true;

    report();
  }

  void SHM::close(){
    CHECK(opened,"can't close what's not open");
    info->guests--;
    info->guestPids[guestId]=-1;
    shmdt(p);
    created=opened=false; p=NULL; systemId=NULL; info=NULL; guestId=-1;
  }
  
  void SHM::destroy(){
    std::cout <<"** destroying shared memory with systemId "<<systemId <<std::endl;
    shmid_ds buf;
    if(shmctl(systemId, IPC_RMID, &buf)<0){
      HALT("failed destroying the shm");
    }
    created=opened=false; p=NULL; systemId=NULL; info=NULL; guestId=-1;
  }
  
  //documentation: http://cs.pub.ro/~apc/2003/resources/pthreads/uguide/document.htm
  void SHM::readLock(uint i){
    int r=pthread_rwlock_rdlock(&info->blockLocks[i]);
    if(r) HALT("");
  }
  
  void SHM::writeLock(uint i){
    int r=pthread_rwlock_wrlock(&info->blockLocks[i]);
    if(r) HALT("");
  }
  
  void SHM::unlock(uint i){
    int r=pthread_rwlock_unlock(&info->blockLocks[i]);
    if(r) HALT("");
  }
  
  void SHM::report(){
    std::cout <<"** SHM report:";
    if(!opened){ std::cout <<" -- not open" <<std::endl; return; }
    uint i;
    std::cout
        <<"\n  pagename=" <<info->pagename
        <<"\n  size=" <<info->size
        <<"\n  used=" <<info->used
        <<"\n  guests=" <<info->guests
        <<"\n  guestPids=";
    for(i=0;i<maxShmGuests;i++) std::cout <<' ' <<info->guestPids[i];
    std::cout
        <<"\n  local guestId=" <<guestId
        <<"\n  local systemId=" <<systemId
        <<"\n  local p=" <<(void*)p
        <<"\n  locally created=" <<created
        <<"\n  allocated blocks:";
    for(i=0;i<shmMaxBlocks;i++) if(info->blockOffsets[i]) std::cout <<"\n    " <<i <<": name=" <<info->blockNames[i]<<" offset=" <<info->blockOffsets[i];
    std::cout <<std::endl;
  }
#endif //PTHREAD
#endif
  
  
  //! start and reset the timer (user CPU time)
  void timerStart(bool useRealTime){
    if(useRealTime) timerUseRealTime=true; else timerUseRealTime=false;
    timerPauseTime=-1.;
    timerStartTime=(timerUseRealTime?realTime():cpuTime());
  }

  //! read the timer and optionally also reset it (user CPU time)
  real timerRead(bool reset){
    real c;
    if(timerPauseTime!=-1.) c=timerPauseTime; else c=(timerUseRealTime?realTime():cpuTime())-timerStartTime;
    if(reset) timerStart(timerUseRealTime);
    return c;
  }

  //! read the timer relative to a given start time (user CPU time)
  real timerRead(bool reset,real startTime){
    real c=(timerUseRealTime?realTime():cpuTime())-startTime;
    if(reset) timerStart(timerUseRealTime);
    return c;
  }

  //! read and pause the timer
  real timerPause(){
    timerPauseTime=(timerUseRealTime?realTime():cpuTime())-timerStartTime;
    return timerPauseTime;
  }

  //! resume the timer after a pause, neglecting the time inbetween
  void timerResume(){
    timerStartTime=(timerUseRealTime?realTime():cpuTime())-timerPauseTime;
    timerPauseTime=-1.;
  }

  //! memorize the command line arguments and open a log file
  void initCmdLine(int _argc, char *_argv[]){
    argc=_argc; argv=_argv;
    time_t t; time(&t);
    const char *name;
    if(checkCmdLineTag("nolog")) noLog=true; else noLog=false;
    name=getCmdLineArgument("log");
    if(!name) name=MT_LogFileName;
    log(name);
    if(!log().good()) MT_MSG(" -- use `-nolog' or `-log' option to specify the log file");
    log() <<"Compiled at:     " <<__DATE__ <<" "<<__TIME__ <<"\n";
    log() <<"Execution start: " <<ctime(&t);
    log() <<"Program call:    '"; for(int i=0;i<argc;i++) log() <<argv[i] <<" ";
    log() <<"\b'" <<std::endl;
  }

  //! returns true if the tag was found on command line
  bool checkCmdLineTag(const char *tag){
    for(int n=1;n<argc;n++) if(argv[n][0]=='-' && !strcmp(tag,argv[n]+1)){
      return true;
    }
    return false;
  }

  //! returns the argument after the cmd-line tag; NULL if the tag is not found
  char *getCmdLineArgument(const char *tag){
    int n;
    for(n=1;n<argc;n++) if(argv[n][0]=='-' && !strcmp(tag,argv[n]+1)){
      if(n+1==argc) return (char*)"1";
      return argv[n+1];
    }
    return NULL;
  }

  /*!\brief Open a (possibly new) config file with name '\c name'.<br> If
    \c name is not specified, it searches for a command line-option
    '-cfg' and, if not found, it assumes \c name=MT.cfg */
  void openConfigFile(const char *name){
    log() <<"opening config file ";
    if(!name) name=getCmdLineArgument("cfg");
    if(!name) name=MT_ConfigFileName;
    if(cfgOpenFlag){ cfgFile.close(); log() <<"(old config file closed) "; }
    log() <<"'" <<name <<"'";
    cfgFile.clear();
    cfgFile.open(name);
    cfgOpenFlag=true;
    if(!cfgFile.good()){
      //MT_MSG("couldn't open config file "<<name);
      log() <<" - failed";
    }
    log() <<std::endl;
  }
}

//! a global operator to scan (parse) strings from a stream
std::istream& operator>>(std::istream& is,const char *str){
  MT::parse(is,str); return is;
}
//! the same global operator for non-const string
std::istream& operator>>(std::istream& is,char *str){
  MT::parse(is,(const char*)str); return is;
}


//===========================================================================
//
// String class
//

          int MT::String::StringBuf::overflow(int C){
            string->append(C);
            return C;
          }
          int MT::String::StringBuf::sync(){
            if(string->flushHandler) string->flushHandler(*string);
            return 0;
          }
          void MT::String::StringBuf::setIpos(char *p){ setg(string->p,p,string->p+string->memN); }
          char *MT::String::StringBuf::getIpos(){ return gptr(); }
    
      //-- direct memory operations
          void MT::String::append(char x){ resize(memN+1,true); operator()(memN-1)=x; }
          void MT::String::resize(uint n,bool copy){
        if(memN==n && M>memN) return;
        char *pold=p;
        uint Mold=M;
        //flexible allocation (more than needed in case of multiple resizes)
        if(M==0){ //first time
          M=n+1;
        }else if(n+1>M || 10+2*n<M/2){
          M=11+2*n;
        }
        if(M!=Mold){
          p=new char [M];
          if(!p) HALT("MT::Mem failed memory allocation of "<<M<<"bytes");
          if(copy) memmove(p,pold,memN<n?memN:n);
          if(Mold) delete[] pold;
        }
        memN=n;
        p[memN]=0;
        resetI();
      }
    
      void MT::String::init(){ p=0; memN=0; M=0; buffer.string=this; flushHandler=NULL; }
      //! standard constructor
      MT::String::String():IOStream(&buffer){ init(); clr(); }
      //! copy constructor
      MT::String::String(const String& s):IOStream(&buffer){ init(); this->operator=(s); }
      //! copy constructor for an ordinary C-string (needs to be 0-terminated)
      MT::String::String(const char *s):IOStream(&buffer){ init(); this->operator=(s); }
      MT::String::~String(){ if(M) delete[] p; }

      
      //! returns a reference to this
      MT::IOStream& MT::String::stream(){ return (MT::IOStream&)(*this); }
      //! returns a reference to this
      MT::String& MT::String::operator()(){ return *this; }
      /*!\brief returns the true memory buffer (C-string) of this class (which is
      always kept 0-terminated) */
      MT::String::operator char*(){ return p; }
      //! as above but const
      MT::String::operator const char*() const{ return p; }
      //! returns the i-th char
      char& MT::String::operator()(uint i) const{ CHECK(i<=memN,"String range error ("<<i<<"<="<<memN<<")"); return p[i]; }

      //! copy operator
      MT::String& MT::String::operator=(const String& s){
        resize(s.memN,false);
        memmove(p,s.p,memN);
        return *this;
      }

      //! copies from the C-string
      void MT::String::operator=(const char *s){ resize(strlen(s),false); memmove(p,s,strlen(s)); }
      void MT::String::set(const char *s,uint n){ resize(n,false); memmove(p,s,n); }
    //! shorthand for the !strcmp command
      bool MT::String::operator==(const char *s){ return !strcmp(p,s); }
    //! shorthand for the !strcmp command
      bool MT::String::operator==(const String& s){ return !strcmp(p,s.p); }

    //! deletes all memory and resets all stream flags
      MT::String& MT::String::clr(){ resize(0,false); return *this; }
    /*!\brief when using this String as an istream (to read other variables
      from it), this method resets the reading-pointer to the beginning
      of the string and also clears all flags of the stream */
      MT::String& MT::String::resetI(){ buffer.setIpos(p); clear(); return *this; }
    //! set the reading-pointer to an arbitrary place
      //void setIpos(char *p){ buffer.setIpos(p); }
    //! get the reading-pointer
      //char *getIpos(){ return buffer.getIpos(); }

    //! returns the length of the string (excluding the terminating 0)
      uint MT::String::N() const{ return memN; }

      //! writes the string into some ostream
      void MT::String::write(std::ostream& os) const{ os <<p; }
      /*!\brief reads the string from some istream: first skip until one of the stopSymbols
      is encountered (default: newline symbols) */
      void MT::String::read(std::istream& is,const char* skipSymbols,const char *stopSymbols,int eatStopSymbol){
        if(!skipSymbols) skipSymbols=readSkipSymbols;
        if(!stopSymbols) stopSymbols=readStopSymbols;
        if(eatStopSymbol==-1) eatStopSymbol=readEatStopSymbol;
        MT::skip(is,skipSymbols);
        clr();
        char c=is.get();
        while(c!=-1 && is.good() && !MT::contains(stopSymbols,c)){
          append(c);
          c=is.get();
        }
        if(c==-1) is.clear();
        if(c!=-1 && !eatStopSymbol) is.putback(c);
        //is.clear();
      }


//===========================================================================
//
// random number generator
//

namespace MT{ Rnd rnd; }

uint32_t MT::Rnd::seed(uint32_t n){
  uint32_t s,c;
  if(n>12345){ s=n; c=n%113; }
  else{ s=12345; c=n; }
  while(c--) s*=65539;
  s=s>>1;
  seed250(s);
  ready=true;
  return n;
}

uint32_t MT::Rnd::seed(){
  return seed(MT::getParameter<uint32_t>("seed",0));
}

uint32_t MT::Rnd::clockSeed(){
  uint32_t s;
#ifndef MT_TIMEB
  timeval t; gettimeofday(&t,0); s=1000000L*t.tv_sec+t.tv_usec;
#else
  _timeb t; _ftime(&t); s=1000L*t.time+t.millitm;
#endif
  log() <<"random clock seed: " <<s <<std::endl;
  return seed(s);
}

real MT::Rnd::gauss(){
  real w,v,rsq,fac;
  do{
    v   = 2 * uni() - 1;
    w   = 2 * uni() - 1;
    rsq = v*v + w*w;
  } while( rsq >= 1 || rsq == 0 );
  fac  = ::sqrt( -2 * ::log( rsq ) / rsq );
  return v*fac;
}

uint32_t MT::Rnd::poisson(real mean){
  if(mean>100){
    uint32_t i=(uint32_t)::floor(mean+::sqrt(mean)*gauss()+.5);
    return (i>0)?(uint32_t)i:0;
  }
  uint32_t count = 0;
  real bound, product;
  if(mean>=0){
    bound=::exp(-mean);
    product=uni();
    while(product>=bound){
      count++;
      product*=uni( );
    }
  }
  return count;
}

void  MT::Rnd::seed250(int32_t seed){
  int32_t      i; 
  int32_t     j, k;
  
  if(seed<=0) seed=1;
    
  for(i=0;i<250;++i){           // Schleife ueber Zufallsfeld
    k = seed / 127773;          // Modulozufallszahlengenerator
    seed = 16807 * ( seed - k*127773 ) - 2836 * k;
    if(seed<0) seed += 0x7FFFFFFF;
    rfield[i] = seed;
  }

  // Masken ueberlagern
  k = 0x7FFFFFFF;       
  j = 0x40000000;
  for(i=1;i<250;i+=8) rfield[i] = (rfield[i] & k) | j;
    
  // rpoint initialisieren
  rpoint = 249;            
    
  // Anfangszahlen verwerfen
  for(i=0;i<4711;++i) rnd250();
}    


//===========================================================================
//
// gnuplot calls
//

static FILE *MT_gp=NULL;
void gnuplotClose(){
  if(MT_gp){ fflush(MT_gp); fclose(MT_gp); }
}
void gnuplot(const char *command,const char *EPSfile){
  if(!MT_gp){
    MT_gp=popen("/usr/bin/gnuplot","w"); // -persist
    CHECK(MT_gp,"could not open gnuplot pipe");
    fprintf(MT_gp,"set data style lines\n");
  }
  MT::String cmd;
  if(EPSfile){
    cmd <<"set terminal postscript 24 dl 5; set key samplen 2\n"
        <<"set output '" <<EPSfile <<"'\n"
        <<"set data style lines\n"
        <<command <<std::endl;
  }else{
    cmd <<"set terminal pop\n"
        <<"set title '(MT/plot.h -> gnuplot pipe)'\n"
        <<command <<std::endl;
  }
  fputs(cmd.p,MT_gp);
  fflush(MT_gp) ;
#if 0 //old!
  std::ofstream os("z.plt");
#ifndef MT_MSVC
  os
    <<"set terminal x11\n" //persist noraise\n"
    <<"set data style lines\n"
    <<command <<std::endl
    <<"pause mouse" <<std::endl;
  os.close();
  if(system("gnuplot z.plt &"))
    MT_MSG("attempt to call gnuplot on via `gnuplot z.plt' failed");
#else
  os
    <<"set data style lines\n"
    <<command <<std::endl
    <<"pause mouse" <<std::endl;
  os.close();
  if(system("wgnuplot z.plt &"))
    MT_MSG("attempt to call gnuplot on windows via `wgnuplot z.plt' failed (probably not installed)");
#endif
#endif
}


//===========================================================================
//
// Cumulative probability for the Standard Normal Distribution 
//

namespace MT{
  real erf(real x){
    real t, z, retval;
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
  
  //! the integral of N(0,1) from -infty to x
  real gaussInt(real x){
    return .5*(1.+erf(x/MT_SQRT2));
  }

  //! expectation \f$\int_x^\infty {\cal N}(x) x dx\f$ when integrated from -infty to x
  real gaussIntExpectation(real x){
    real norm=gaussInt(x) / (::sqrt(MT_2PI));
    return - norm*MT::approxExp(-.5*x*x);
  }
}
