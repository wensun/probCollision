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

#define MT_util_t_cpp

#include "util.h"

namespace MT{
  /*!\brief a standard method to save an object into a file. The same as
    std::ofstream file; MT::open(file,filename); file <<x;
    file.close(); */
  template<class T> void save(const T& x,const char *filename){
    std::ofstream file;
    open(file,filename);
    file <<x;
    file.close();
  }

    /*!\brief a standard method to load object from a file. The same as
    std::ifstream file; MT::open(file,filename); file >>x;
    file.close(); */
  template<class T> void load(T& x,const char *filename,bool change_directory){
    if(!change_directory){
      std::ifstream file;
      open(file,filename);
      file >>x;
      file.close();
    }else{
      char *path,*name,cwd[200];
      MT::decomposeFilename(path,name,filename);
      if(!getcwd(cwd,200)) HALT("couldn't get current dir");
      if(path[0]) if(chdir(path)) HALT("couldn't change to directory "<<path);
      std::ifstream file;
      open(file,name);
      file >>x;
      file.close();
      if(path[0]) if(chdir(cwd)) HALT("couldn't change to directory "<<cwd);
    }
  }

    /*!\brief Search for a command line option \c -tag and, if found, pipe the
     next command line option into \c value by the
     \c operator>>(istream&,type&). Returns false on failure. */
  template<class T>
  bool getFromCmdLine(T& x,const char *tag){
    char *opt=getCmdLineArgument(tag);
    if(!opt) return false;
    std::istringstream s(opt);
    s >>x;
    if(s.fail()) HALT("error when reading parameter from command line: " <<tag);
    return true;
  }

  /*!\brief Search the first occurence of a sequence '\c tag:'
  in the config file (opened automatically) and, if found, pipes
  it in \c value. Returns false if parameter is not found. */
  template<class T>
  bool getFromCfgFile(T& x,const char *tag){
    if(!cfgOpenFlag) openConfigFile();
    CHECK(!cfgLock,"cfg file is locked");
    cfgLock=true;
    cfgFile.clear();
    cfgFile.seekg(std::ios::beg);
    if(!cfgFile.good()){ cfgLock=false; return false; }
    
    unsigned n=strlen(tag);
    char *buf=new char [n+2]; memset(buf,0,n+2);
    while(cfgFile.good()){
      memmove(buf,buf+1,n);
      buf[n]=cfgFile.get();
      if(buf[n]==' ' || buf[n]=='\t' || buf[n]==':' || buf[n]=='='){ buf[n]=0; if(!strcmp(tag,buf)) break; buf[n]=':'; }
    };
    delete[] buf;
    
    if(!cfgFile.good()){ cfgLock=false; return false; }
    
    skip(cfgFile," :=\n\r\t");
    cfgFile >>x;
    
    if(cfgFile.fail()) HALT("error when reading parameter " <<tag);
    cfgLock=false;
    return true;
  }
  
  template<class T>
  bool getParameterBase(T& x,const char *tag,bool hasDefault,const T* Default){
    log() <<std::setw(20) <<tag <<" = " <<std::setw(5);
    log().flush();
    
    if(getFromCmdLine(x,tag)){
      log() <<x <<" [" <<typeid(x).name() <<"] (cmd line!)" <<std::endl;
      return true;
    }
    
    if(getFromCfgFile(x,tag)){
      log() <<x <<" [" <<typeid(x).name() <<"]" <<std::endl;
      return true;
    }
    
    if(hasDefault){
      if(Default){
        x=*Default;
        log() <<x <<" [" <<typeid(x).name() <<"] (default!)"<<std::endl;
      }
      return false;
    }
    
    HALT("could not initialize parameter `" <<tag
	     <<"': parameter has no default;\n     either use command option `-"
	     <<tag <<" ...' or specify `"
	     <<tag <<"= ...' in the config file `MT.cfg'");
  }

  template<class T> T getParameter(const char *tag){
    T x;
    getParameterBase<T>(x,tag,false,(T*)NULL);
    return x;
  }
  template<class T> T getParameter(const char *tag,const T& Default){
    T x;
    getParameterBase<T>(x,tag,true,&Default);
    return x;
  }
  template<class T> void getParameter(T& x,const char *tag,const T& Default){
    getParameterBase<T>(x,tag,true,&Default);
  }
  template<class T> void getParameter(T& x,const char *tag){
    getParameterBase(x,tag,false,(T*)NULL);
  }
  template<class T> bool checkParameter(const char *tag){
    T x;
    return getParameterBase(x,tag,true,(T*)NULL);
  }
  template<class T> void Parameter<T>::initialize(){
    if(!initialized){
      getParameterBase(value,tag,hasDefault,&Default);
      initialized = true;
    }
  }
}


//===========================================================================
//
// generic any container
//
  
template<class T> struct ANY:public Any{
  virtual ~ANY(){ free(); };
  ANY(const char* _tag,const T &x){          tag=NULL; p=NULL; set(_tag,&x,0,0);  }
  ANY(const char* _tag,const T *_p,uint _n,char _delim){ tag=NULL; p=NULL; set(_tag,_p,_n,_delim); }
  virtual void write(std::ostream &os) const{
    if(!p){ os <<tag; return; } //boolean
    os <<tag <<"="; //<<"[" <<type <<"] = ";
    if(!n){
      if(typeid(T)==typeid(const char*) || typeid(T)==typeid(char*) || typeid(T)==typeid(MT::String)) os <<'\'' <<*((T*)p) <<'\'';
      else os <<*((T*)p);
    }else{
      T *t=(T*)p;
      os <<delim <<t[0];
      for(uint i=1;i<n;i++) os <<' ' <<t[i];
      if(delim=='(') os <<')';
      else if(delim=='[') os <<']';
      else if(delim=='{') os <<'}';
      else os <<delim;
    }
  }
  virtual void free(){
    if(!tag){ CHECK(!p,""); return; }
    delete[] tag;
    if(!p) return;
    if(!n) delete ((T*)p);
    else   delete[] ((T*)p);
    p=NULL;
  }
  virtual void set(const char* _tag,const T *_p,uint _n,char _delim){
    free();
    type=typeid(T).name();
    tag=new char[strlen(_tag)+1];
    strcpy(tag,_tag);
    if(!_p){ p=NULL; n=0; delim=0; return; } //assume this is a ``boolean tag'' without data
    n=_n;
    delim=_delim;
    if(!n){
      p = new T(_p[0]);
    }else{
      p = new T[n];
      T *t=(T*)p;
      for(uint i=0;i<n;i++) t[i]=_p[i];
    }
  }
  virtual Any* new_clone(){ return new ANY<T>(tag,(T*)p,n,delim); }
};

template<class T> Any* anyNew(const char* tag,const T &x){        return new ANY<T>(tag,x); }
template<class T> Any* anyNew(const char* tag,const T *x,uint n,char delim){ return new ANY<T>(tag,x,n,delim); }

#if 0 //don't need the bag anymore, use an any list instead!
namespace MT{
  class BagItem{
  public:
    BagItem *next; 
    void *value;
    int type;
    BagItem(){ next=0; value=0; type=-1; }
    void operator=(const BagItem&){ NIY; }
    virtual ~BagItem(){}
    virtual void write(std::ostream &os) = 0;
    virtual void free() = 0;
    virtual void set(void *x) = 0;
    virtual BagItem *newClone() = 0;
    //virtual get(void *x) = 0;
  };
  

  struct BagItemType{ char *tag; const char *typeidname; uint n; };
  extern Mem<BagItemType> BagItemTypes;

  template<class T> class BagTypedItem:public BagItem{
  public:
    virtual void write(std::ostream &os){
      if(!value) return;
      if(!BagItemTypes(type).n){
        os <<*((T*)value);
      }else{
        T *t=(T*)value;
        uint n=BagItemTypes(type).n;
        for(uint i=0;i<n;i++,t++) os <<' ' <<*t;
      }
    }
    virtual void free(){
      if(!value) return;
      if(!BagItemTypes(type).n){
        delete ((T*)value);
      }else{
        delete[] ((T*)value);
      }
      value=0;
    }
    virtual void set(void *_x){
      if(!_x){ value=0; return; }
      uint n=BagItemTypes(type).n;
      T *t=(T*)value;
      T *x=(T*)_x;
      if(value) free();
      if(!n){
        value = new T;
        t=(T*)value;
        *t=*x;
      }else{
        value = new T[n];
        t=(T*)value;
        for(uint i=0;i<n;i++) t[i]=x[i];
      }
    }
    virtual BagItem *newClone(){ return new BagTypedItem; }
  };

#ifndef MT_MSVC
  //! add an arbitrary-type item to the bag
  template<class T> void Bag::add(const char *tag,T *x,uint n){
    BagItem *a=new BagTypedItem<T>;
    a->type=getType(tag,typeid(T).name(),n);
    if(a->type==-1){ delete a; return; }
    a->set(x);
    a->next=first;
    first=a;
  }
  
  //! get the item with tag from the bag
  template<class T> BagItem* Bag::get(const char *tag,T *&x,uint n){
    int i;
    BagItem *a;
    i=getType(tag,typeid(T).name(),n);
    for(a=first;a;a=a->next) if(a->type==i) break;
    if(a) x=(T*)a->value; else x=NULL;
    return a;
  }
#endif
}
#endif

template<class T> T* MT::SHM::newBlock(const char* objName,uint *_i){
  if(!opened) open("shmPageName",1000);
  uint i;
  for(i=0;i<shmMaxBlocks;i++) if(!strcmp(objName,info->blockNames[i])) break;
  if(i==shmMaxBlocks){ //allocate new
    for(i=0;i<shmMaxBlocks;i++) if(!info->blockOffsets[i]) break;
    if(i==shmMaxBlocks) HALT("exceeded number of shm blocks");
    strcpy(info->blockNames[i],objName);
    info->blockOffsets[i] = info->used;
    info->used += sizeof(T);
    if(info->used>info->size) HALT("not enough space in shared memory");
  }
  T *x=(T*)(p+info->blockOffsets[i]);
  if(_i) *_i=i;
  return x;
}
