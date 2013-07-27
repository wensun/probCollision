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

/*! \file opengl.h
    \brief defines the OpenGL interface to freeglut or Qt */

#ifndef MT_opengl_h
#define MT_opengl_h

#undef MT_GL
#if defined MT_QTGLUT || defined MT_FREEGLUT
#  include <GL/gl.h>
#  include <GL/glu.h>
#  define MT_GL
#  ifdef MT_QTGLUT
#    undef MT_FREEGLUT
#  endif
#  ifdef MT_MSVC
#    include<windows.h>
#    undef min //I hate it that windows defines these macros!
#    undef max
#  endif
#endif

#ifdef MT_QTGLUT
//#  if !defined MT_Cygwin 
#    include <GL/glut.h>
//#  endif
#  undef scroll
#  undef border
// #  define QT3_SUPPORT
// #  include <Qt/Qt3Support>
// #  include <Qt/qtimer.h>
#  include <QtGui/QApplication>
#  include <QtOpenGL/QGLWidget>
// #  include <Qt/qobject.h>
// #  include <Qt/qevent.h>
#  include <QtOpenGL/QtOpenGL>
#  if defined MT_Cygwin //|| defined MT_Linux
#    define GLformat QGL::DirectRendering | QGL::DepthBuffer | QGL::Rgba
#    define GLosformat QGL::DirectRendering | QGL::DepthBuffer | QGL::Rgba
#  else
#    define GLformat QGL::DirectRendering | QGL::DepthBuffer | QGL::Rgba
#    define GLosformat QGL::DirectRendering | QGL::DepthBuffer | QGL::Rgba
#  endif
#endif

#ifdef MT_FREEGLUT
#  define FREEGLUT_STATIC
#  include <GL/freeglut.h>
extern "C"{
  void FGAPIENTRY glutMainLoopMT( void );
  void fgDeinitialize( void );
}
#endif

#ifdef MT_GL2PS
#  include<gl2ps.h>
#endif

#include "util.h"
#include "array.h"


//===========================================================================
//
// utility functions
//

void glStandardLight();
void glStandardScene(void*);
void glColor(float r, float g, float b,float a=1.f);
void glColor(int col);
void glDrawText(const char* txt,float x,float y,float z);
//void glShadowTransform();
void glTransform (const double pos[3], const double R[12]);
void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
		float x3,float y3,float z3,float x4,float y4,float z4,
		float r,float g,float b);
void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
		float x3,float y3,float z3,float x4,float y4,float z4);
void glDrawRect(float x,float y,float z,float rad);
void glDrawFloor(float x,float r,float g,float b);
void glDrawBox(float x,float y,float z);
void glDrawDiamond(float dx,float dy,float dz);
void glDrawDiamond(float x,float y,float z,float dx,float dy,float dz);
void glDrawSphere(float radius);
void glDrawDisk(float radius);
void glDrawCylinder(float radius, float length,bool closed=true);
void glDrawCappedCylinder(float radius, float length);
void glDrawAxes(double scale);
void glDrawGridBox(float x);
void glDrawGridBox(float x1,float y1,float z1,float x2,float y2,float z2);
void glDrawKhepera();
void glMakeSquare(int num);
void glMakeStdSimplex(int num);
void glMakeTorus(int num);
void glDrawRobotArm(float a,float b,float c,float d,float e,float f);
uint glImageTexture(const byteA &img);
void glDrawTexQuad(uint texture,
		   float x1,float y1,float z1,float x2,float y2,float z2,
	           float x3,float y3,float z3,float x4,float y4,float z4,
		   float mulX=1.,float mulY=1.);
void glGrabImage(byteA& img);
void glGrabDepth(byteA& depth);
void glGrabDepth(floatA& depth);
void glRasterImage(int x,int y,byteA &img,float zoom=1.);


//===========================================================================
//
// standalone draw routines for larget data structures
//

void glDrawDots(arr& dots);
void glDrawDots(void *dots);


//===========================================================================
//
// basic gui routines - wrapped for compatibility
//

void MTEvents();
void MTenterLoop();
void MTexitLoop();



//===========================================================================
//
// OpenGL class
//

namespace ors{
  struct Frame;
  struct Vector;
  struct Camera{
    Frame *f;
    Vector *foc,*offset;
    bool ownFrame;
    float phi,theta;
    
    float heightAbs;
    float heightAngle;
    float whRatio;
    float zNear,zFar;

    Camera(Frame *frame=0);
    Camera(const Camera& c){ *this=c; }
    ~Camera();
    Camera& operator=(const Camera& c);
    
    void setZero();
    void setHeightAngle(float a);
    void setHeightAbs(float h);
    void setZRange(float znear,float zfar);
    void setWHRatio(float ratio);
    void setPosition(float x,float y,float z);
    void setOffset(float x,float y,float z);
    void focusOrigin();
    void focus(float x,float y,float z);
    void focus(const Vector& v);
    void focus();
    void watchDirection(float x,float y,float z);
    void upright();
    void glSetProjectionMatrix();
    void glConvertToTrueDepth(double &d);
    void glConvertToLinearDepth(double &d);
  };
}

struct OpenGLWorkspace;

#ifndef MT_QTGLUT
class QGLWidget{};
#undef  Q_OBJECT
#define Q_OBJECT
#endif

/*!\brief A class to display and control 3D scenes using OpenGL and Qt.

    Minimal use: call \ref add to add routines or objects to be drawn
    and \ref update or \ref watch to start the display. */
class OpenGL:public QGLWidget{
  Q_OBJECT
public:
  //!@name little structs to store objects and callbacks
  struct GLDrawer   { void *classP; void (*call)(void*); };
  struct GLInitCall { void *classP; bool (*call)(void*,OpenGL*); };
  struct GLHoverCall{ void *classP; bool (*call)(void*,OpenGL*); };
  struct GLClickCall{ void *classP; bool (*call)(void*,OpenGL*); };
  struct GLEvent    { int button,key,x,y; float dx,dy; void set(int b,int k,int _x,int _y,float _dx,float _dy){ button=b; key=k; x=_x; y=_y; dx=_dx; dy=_dy; } };
  struct GLSelect   { int name; double dmin,dmax; };
    
  //!@name data fields
  OpenGLWorkspace *WS;               //!< ..
  static MT::Array<OpenGL*> glwins;  //!< global window list
  int windowID;                      //!< id of this window in the global glwins list
  MT::Array<GLDrawer> drawers;         //!< list of draw routines
  MT::Array<GLInitCall> initCalls;     //!< list of initialization routines
  MT::Array<GLHoverCall> hoverCalls;   //!< list of hover callbacks
  MT::Array<GLClickCall> clickCalls;   //!< list of click callbacks
  ors::Camera camera; //!< the camera used for projection
  String text;        //!< the text to be drawn as title within the opengl frame
  float clearR,clearG,clearB,clearA;  //!< colors of the beackground (called in glClearColor(...))
  bool reportEvents,reportSelects;    //!< flags for verbosity
  bool selectOnHover;
  int pressedkey;           //!< stores the key pressed
  const char *exitkeys;     //!< additional keys to exit watch mode
  int mouse_button;         //!< stores which button was pressed
  int mouseposx,mouseposy;  //!< current x- and y-position of mouse
  bool mouseIsDown;
  MT::Array<GLSelect> selection; //!< list of all selected objects
  GLSelect *topSelection;        //!< top selected object
  bool immediateExitLoop;
  bool drawFocus;
  byteA *img;
  double zoom;
  
  //!@name constructors & destructors
#ifdef MT_FREEGLUT
  OpenGL(const char* title="MT::OpenGL(Freeglut)",int w=400,int h=400,int posx=-1,int posy=-1);
#endif
#ifdef MT_QT
  OpenGL(const char* title="MT::OpenGL(Qt)",int w=400,int h=400,int posx=-1,int posy=-1);
  OpenGL(QWidget *parent,const char* title,int width=400,int height=400,int posx=-1,int posy=-1);
#endif
#ifndef MT_GL
  OpenGL(const char* title="MT::OpenGL(NONE)",int w=400,int h=400,int posx=-1,int posy=-1);
#endif
  
  void clone(const OpenGL& gl);
  ~OpenGL();

  //!@name adding drawing routines and callbacks
  void add(void (*call)(void*),const void* classP=0);
  void remove(void (*call)(void*),const void* classP=0);
  template<class T> void add(const T& x){ add(x.staticDraw,&x); } //!< add a class or struct with a staticDraw routine
  void clear();
  void addHoverCall(bool (*call)(void*,OpenGL*),const void* classP=0);
  void clearHoverCalls();
  void addClickCall(bool (*call)(void*,OpenGL*),const void* classP=0);
  void clearClickCalls();

  //!@name the core draw routines (actually only for internal use)
  void Draw(int w,int h,ors::Camera *cam=NULL);
  void Select();  

  //!@name showing, updating, and watching
  bool update(const char *text=NULL);
  int  watch(const char *text=NULL);
  int  timedupdate(double sec);
  void resize(int w,int h);
  void setClearColors(float r,float g,float b,float a);
  void unproject(double &x,double &y,double &z);

  //!@name info & I/O
  void reportSelection();
  int width();
  int height();
  void saveEPS(const char *filename);
  void about(std::ostream& os=std::cout);

  //!@name to display image data (kind of misuse)
  void watchImage(const byteA &img,bool wait,float zoom);
  void watchImage(const floatA &img,bool wait,float zoom);
  void displayGrey(const arr &x,uint d0,uint d1,bool wait,uint win);
  void displayRedBlue(const arr &x,uint d0,uint d1,bool wait,uint win);

  //!@name capture routines
  void capture(byteA &img,int w,int h,ors::Camera *cam);
  void captureStereo(byteA &imgL,byteA &imgR,int w,int h,ors::Camera *cam,double baseline);

#ifdef MT_QTGLUT
  void createOffscreen(int width,int height);
  void offscreenGrab(byteA& image);
  void offscreenGrab(byteA& image,byteA& depth);
  void offscreenGrabDepth(byteA& depth);
  void offscreenGrabDepth(floatA& depth);
private:
  void setOffscreen(int width,int height);
#endif

  private:
    GLEvent lastEvent;
    static uint selectionBuffer[1000];
#ifdef MT_QTGLUT
    bool quitLoopOnTimer;
    QPixmap *osPixmap;      // the paint device for off-screen rendering
    QGLContext *osContext;  //the GL context for off-screen rendering
#endif

    void init(); //initializes camera etc
  //general callbacks (used by QT & Freeglut)
    void Key(unsigned char key, int x, int y);
    void Mouse(int button, int updown, int x, int y);
    void Motion(int x, int y);
    void PassiveMotion(int x, int y);
    void Close(){ }
    void Reshape(int w,int h);
    void Special(int key, int x, int y);
    void MouseWheel(int wheel, int direction, int x, int y);
#ifdef MT_FREEGLUT
  //hooks for FREEGLUT (static callbacks)
    static void _Void(){ }
    static void _Draw(){ OpenGL *gl=glwins(glutGetWindow()); gl->Draw(gl->width(),gl->height()); glutSwapBuffers(); }
    static void _Key(unsigned char key, int x, int y){ glwins(glutGetWindow())->Key(key,x,y); }
    static void _Mouse(int button, int updown, int x, int y){ glwins(glutGetWindow())->Mouse(button,updown,x,y); }
    static void _Motion(int x, int y){ glwins(glutGetWindow())->Motion(x,y); }
    static void _PassiveMotion(int x, int y){ glwins(glutGetWindow())->PassiveMotion(x,y); }
    static void _Close(){ glwins(glutGetWindow())->Close(); }
    static void _Reshape(int w,int h){ glwins(glutGetWindow())->Reshape(w,h); }
    static void _Special(int key, int x, int y){ glwins(glutGetWindow())->Special(key,x,y); }
    static void _MouseWheel(int wheel, int direction, int x, int y){ glwins(glutGetWindow())->MouseWheel(wheel,direction,x,y); }
#endif
#ifdef MT_QTGLUT
  //hooks for Qt (overloading virtuals)
    void paintGL(){ Draw(width(),height()); }
    void initializeGL(){ }
    void resizeGL(int w,int h){ Reshape(w,h); }
    void keyPressEvent(QKeyEvent *e){ pressedkey=e->text().toAscii()[0]; Key(pressedkey,mouseposx,mouseposy); }
    void timerEvent(QTimerEvent*){ if(quitLoopOnTimer) MTexitLoop(); }
    void mouseMoveEvent(QMouseEvent* e){
      if(!mouseIsDown) PassiveMotion(e->x(),e->y()); else Motion(e->x(),e->y());
    }
    void mousePressEvent(QMouseEvent* e){
      if(e->button()==Qt::LeftButton) { Mouse(0,0,e->x(),e->y()); }
      if(e->button()==Qt::MidButton)  { Mouse(1,0,e->x(),e->y()); }
      if(e->button()==Qt::RightButton){ Mouse(2,0,e->x(),e->y()); }
    }
    void mouseReleaseEvent(QMouseEvent* e){
      if(e->button()==Qt::LeftButton) { Mouse(0,1,e->x(),e->y()); }
      if(e->button()==Qt::MidButton)  { Mouse(1,1,e->x(),e->y()); }
      if(e->button()==Qt::RightButton){ Mouse(2,1,e->x(),e->y()); }
    }
#endif

};


//===========================================================================
//
// simplest UI
//

class glUI{
public:
  int top;
  struct Button{ byteA img1,img2; bool hover; uint x,y,w,h; const char* name; };
  MT::Array<Button> buttons;

  glUI(){ top=-1; }

  void addButton(uint x,uint y,const char *name,const char *img1=0,const char *img2=0);
  void glDraw();
  bool checkMouse(int _x,int _y);
};

void glDrawUI(void *p);
bool glHoverUI(void *p,OpenGL *gl);
bool glClickUI(void *p,OpenGL *gl);
//--------- implementation


#ifdef MT_IMPLEMENTATION
#  include "opengl.cpp"
#endif

#endif
