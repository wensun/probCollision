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

#include "opengl.h"
#include "ors.h"

//===========================================================================
//
// camera class
//


/*!\brief constructor; specify a frame if the camera is to be attached
   to an existing frame. Otherwise the camera creates its own
      frame */
ors::Camera::Camera(Frame *frame){
  if(!frame){ f=new Frame; ownFrame=true; }else{ f=frame; ownFrame=false; }
  offset = new ors::Vector;
  foc = new ors::Vector;
  setZero();
}

ors::Camera::~Camera(){
  if(ownFrame) delete f;
  delete offset;
  delete foc;
}

ors::Camera& ors::Camera::operator=(const ors::Camera& c){
  //MEM_COPY_OPERATOR(c);
  memmove(this,&c,sizeof(this));
  offset = new ors::Vector(*c.offset);
  foc = new ors::Vector(*c.foc);
  if(ownFrame) f = new ors::Frame(*c.f);
  return *this;
}

void ors::Camera::setZero(){
  f->setZero();
  offset->setZero();
  theta=phi=0.;
  heightAngle=90.;
  heightAbs=10.;
  whRatio=1.;
  zNear=.1;
  zFar=1000.;
}

  //! the height angle (in degrees) of the camera perspective; set it 0 for orthogonal projection
      void ors::Camera::setHeightAngle(float a){ heightAngle=a; }
  //! the absolute height of the camera perspective (automatically also sets heightAngle=0)
      void ors::Camera::setHeightAbs(float h){ heightAngle=0; heightAbs=h; }
  //! the z-range (depth range) visible for the camera
      void ors::Camera::setZRange(float znear,float zfar){ zNear=znear; zFar=zfar; }
  //! set the width/height ratio of your viewport to see a non-distorted picture
      void ors::Camera::setWHRatio(float ratio){ whRatio=ratio; }
  //! the frame's position
      void ors::Camera::setPosition(float x,float y,float z){ f->p.set(x,y,z); }
  //! the camera's additional offset within the frame
      void ors::Camera::setOffset(float x,float y,float z){ offset->set(x,y,z); }
  //! rotate the frame to focus the absolute coordinate origin (0,0,0)
      void ors::Camera::focusOrigin(){ foc->setZero(); focus(); }
  //! rotate the frame to focus the point (x,y,z)
      void ors::Camera::focus(float x,float y,float z){ foc->set(x,y,z); focus(); }
  //! rotate the frame to focus the point given by the vector
      void ors::Camera::focus(const Vector& v){ *foc=v; focus(); }
  //! rotate the frame to focus (again) the previously given focus
      void ors::Camera::focus(){ Vector v(*foc-f->p); watchDirection(v(0),v(1),v(2)); } //f->Z=f->p; f->Z-=foc; f->Z.normalize(); upright(); }
  //! rotate the frame to watch in the direction vector D
      void ors::Camera::watchDirection(float x,float y,float z){
#if 1
        Vector D(x,y,z),V(0,0,-1);
        f->r.setZero();
        if(D(0)==0. && D(1)==0.){
          if(D(2)>0) f->r.setDeg(180,1,0,0);
          return;
        }
        D.normalize();
        V=f->r*V;
        Quaternion r;
        r.setDiff(V,D);
        f->r=r*f->r;
        upright();
#else
        f->Z.set(-D[0],-D[1],-D[2]); f->Z.normalize();
        upright();
#endif
      }
  //! rotate the frame to set it upright (i.e. camera's y aligned with 's z)
      void ors::Camera::upright(){
#if 1
    //construct desired X:
        Vector V(0,0,-1),X(1,0,0),dX,up;
        X=f->r*X; //true X
        V=f->r*V;
        if(fabs(V(2))<1.) up.set(0,0,1); else up.set(0,1,0);
        dX=up^V; //desired X
        if(dX*X<=0) dX=-dX;
        Quaternion r;
        r.setDiff(X,dX);
        f->r=r*f->r;
#else
        if(f->Z[2]<1.) f->Y.set(0,0,1); else f->Y.set(0,1,0);
        f->X=f->Y^f->Z; f->X.normalize();
        f->Y=f->Z^f->X; f->Y.normalize();
#endif
      }

  //}

  /*! sets OpenGL's GL_PROJECTION matrix accordingly -- should be
      called in an opengl draw routine */
      void ors::Camera::glSetProjectionMatrix(){
#ifdef MT_GL
        if(heightAngle==0)
          glOrtho(-whRatio*heightAbs/2,whRatio*heightAbs/2,
                   -heightAbs/2,heightAbs/2,zNear,zFar);
        else
          gluPerspective(heightAngle,whRatio,zNear,zFar);
        glRotatef(theta,1,0,0);
        glRotatef(phi,0,1,0);
        glTranslatef(-(*offset)(0),-(*offset)(1),-(*offset)(2));
        double m[16];
        glMultMatrixd(f->getInverseAffineMatrixGL(m));
#else
        NIY;
#endif
      }

  //! convert from gluPerspective's non-linear [0,1] depth to the true [zNear,zFar] depth
      void ors::Camera::glConvertToTrueDepth(double &d){
        d = zNear + (zFar-zNear)*d/(zFar/zNear*(1.-d)+1.);
      }

  //! convert from gluPerspective's non-linear [0,1] depth to the linear [0,1] depth
      void ors::Camera::glConvertToLinearDepth(double &d){
        d = d/(zFar/zNear*(1.-d)+1.);
      }


struct OpenGLWorkspace{
  ors::Vector downVec,downPos,downFoc;
  ors::Quaternion downRot;
};


//===========================================================================
//
// static objects
//

OpenGL *staticgl [10]; //ten pointers to be potentially used as display windows
uint nrWins=0;
uint OpenGL::selectionBuffer[1000];
MT::Array<OpenGL*> OpenGL::glwins;


//===========================================================================
//
// basic gui routines - wrapped for compatibility
//

#ifdef MT_FREEGLUT
void MTprocessEvents(){ glutMainLoopEvent(); }
void MTenterLoop(){     glutMainLoopMT(); }
void MTexitLoop(){      glutLeaveMainLoop(); }
#endif
#ifdef MT_QTGLUT
void MTprocessEvents(){ qApp->processEvents(); }
void MTenterLoop(){     qApp->exec(); }
void MTexitLoop(){      qApp->exit(); }
#endif
#ifndef MT_GL
void MTprocessEvents(){ }
void MTenterLoop(){     }
void MTexitLoop(){      }
#endif


//===========================================================================
//
// utility implementations
//

#ifdef MT_GL
void glStandardLight(){
  glEnable(GL_LIGHTING);
#if 1
  static GLfloat ambient[]   = { .5, .5, .5, 1.0 };
  static GLfloat diffuse[]   = { .2, .2, .2, 1.0 };
  static GLfloat specular[]  = { .3, .3, .3, 1.0 };
  static GLfloat position[]  = { 100.0, -100.0, 100.0, 1.0 };
  static GLfloat direction[] = { -1.0, 1.0, -1.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, direction);
  glLighti (GL_LIGHT0, GL_SPOT_CUTOFF,   90);
  glLighti (GL_LIGHT0, GL_SPOT_EXPONENT, 10);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glEnable(GL_LIGHT0);
#endif
#if 0
  static GLfloat diffuse1[]   = { 0.5, 0.5, 0.5, 1.0 };
  static GLfloat specular1[]  = { 0.1, 0.1, 0.1, 1.0 };
  static GLfloat position1[]  = { -100.0, 20.0, -100.0, 1.0 };
  static GLfloat direction1[] = { 1.0, -.2, 1.0 };
  glLightfv(GL_LIGHT1, GL_POSITION, position1);
  glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
  glLighti (GL_LIGHT1, GL_SPOT_CUTOFF,   90);
  glLighti (GL_LIGHT1, GL_SPOT_EXPONENT, 10);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  diffuse1);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specular1);
  glEnable(GL_LIGHT1);
#endif
}

void glStandardScene(void*){
  glStandardLight();
//   glDrawFloor(10,.8,.8,.8);
	glDrawFloor(10,1.5,0.83,.0);
}

void glColor(int col){
  static const GLfloat colorsTab[6][4] = {
    {0.2, 0.2, 1.0, 1.0}, // blue
    {1.0, 0.8, 0.0, 1.0}, // gold
    {1.0, 0.0, 0.0, 1.0}, // red
    {0.7, 0.7, 0.7, 1.0}, // gray
    {1.0, 1.0, 1.0, 1.0}, // white
    {0.2, 1.0, 0.2, 1.0}}; // green

  if (col<0) col=0; if (col>5) col=5;
  glColor(colorsTab[col][0],colorsTab[col][1],colorsTab[col][2],colorsTab[col][3]);
}

void glColor(float r, float g, float b, float alpha){
  float amb=1.f,diff=1.f,spec=1.f;
  GLfloat ambient[4],diffuse[4],specular[4];
  ambient[0] = r*amb;
  ambient[1] = g*amb;
  ambient[2] = b*amb;
  ambient[3] = alpha;
  diffuse[0] = r*diff;
  diffuse[1] = g*diff;
  diffuse[2] = b*diff;
  diffuse[3] = alpha;
  specular[0] = spec*.5*(1.+r);
  specular[1] = spec*.5*(1.+g);
  specular[2] = spec*.5*(1.+b);
  specular[3] = alpha;
#if 0
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, diffuse);
#else
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
  glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
#endif
  glColor4f(r,g,b,alpha);
}

void glColor(float *rgb){ glColor(rgb[0],rgb[1],rgb[2],1.); }

/* // shadows do not work with a light source;
   // thus, we need to leave this out. 4. Mar 06 (hh)
void glShadowTransform()
{
  GLfloat matrix[16];
  for (int i=0; i<16; i++) matrix[i] = 0;
  matrix[0]=1;
  matrix[5]=1;
  matrix[8]=-1;  //light_x
  matrix[9]=-1;  //light_y
  matrix[14]=.02; //ground offset
  matrix[15]=1;
  glPushMatrix();
  glMultMatrixf (matrix);
}
*/

void glTransform(const double pos[3], const double R[12]){
  GLfloat matrix[16];
  matrix[0]=R[0];
  matrix[1]=R[4];
  matrix[2]=R[8];
  matrix[3]=0;
  matrix[4]=R[1];
  matrix[5]=R[5];
  matrix[6]=R[9];
  matrix[7]=0;
  matrix[8]=R[2];
  matrix[9]=R[6];
  matrix[10]=R[10];
  matrix[11]=0;
  matrix[12]=pos[0];
  matrix[13]=pos[1];
  matrix[14]=pos[2];
  matrix[15]=1;
  glPushMatrix();
  glMultMatrixf(matrix);
}

static GLboolean glLightIsOn = false;
void glPushLightOff(){ glGetBooleanv(GL_LIGHTING,&glLightIsOn); glDisable(GL_LIGHTING); }
void glPopLight(){ if(glLightIsOn) glEnable(GL_LIGHTING); }

void glDrawText(const char* txt,float x,float y,float z){
#if defined __glut_h__ || (defined __FREEGLUT_H__ && !defined MT_Cygwin)
  glPushLightOff();
  glRasterPos3f(x,y,z);
  void *font=GLUT_BITMAP_HELVETICA_12;
  while(*txt){ 
    switch(*txt){
    case '\n':
      y+=15;
      glRasterPos3f(x,y,z);
      break;
    case '\b':
      if(font==GLUT_BITMAP_HELVETICA_12) font=GLUT_BITMAP_HELVETICA_18;
      else font=GLUT_BITMAP_HELVETICA_12;
      break;
    default:
      glutBitmapCharacter(font,*txt); 
    }
    txt++;
  }
  glPopLight();
#endif
}

void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
	        float x3,float y3,float z3,float x4,float y4,float z4){
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_POLYGON);
  glVertex3f(x1,y1,z1);
  glVertex3f(x2,y2,z2);
  glVertex3f(x3,y3,z3);
  glVertex3f(x4,y4,z4);
  glVertex3f(x1,y1,z1);
  glEnd();
}

void glDrawRect(float x1,float y1,float z1,float x2,float y2,float z2,
	        float x3,float y3,float z3,float x4,float y4,float z4,
		float r,float g,float b){
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_POLYGON);
  glColor(r,g,b);
  glVertex3f(x1,y1,z1);
  glVertex3f(x2,y2,z2);
  glVertex3f(x3,y3,z3);
  glVertex3f(x4,y4,z4);
  glVertex3f(x1,y1,z1);
  glEnd();
}

void glDrawRect(float x,float y,float z,float r){
  glDrawRect(x-r,y-r,z,x-r,y+r,z,x+r,y+r,z,x+r,y-r,z);
}

void glDrawFloor(float x,float r,float g,float b){
  x/=2.;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor(r,g,b);
  glBegin(GL_POLYGON);
  glNormal3f(0,0,1);
  glVertex3f(-x,-x,0.);
  glVertex3f( x,-x,0.);
  glVertex3f( x, x,0.);
  glVertex3f(-x, x,0.);
  glVertex3f(-x,-x,0.);
  glEnd();
#if 0
  float grey=0.75;
  glColor(grey,grey,grey);
  for(float i=0.0; i<= 1.1;i= i+0.1){
    for(float j=0.0; j <= 1.1;j=j+0.1){
      glBegin(GL_LINES);
      glVertex3f(i,-1,0.001);
      glVertex3f(i,1,0.001);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-i,-1,0.001);
      glVertex3f(-i,1,0.001);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-1,i,0.001);
      glVertex3f(1,i,0.001);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(-1,-i,0.001);
      glVertex3f(1,-i,0.001);
      glEnd();
    }
  }

  grey=0.25;
  glColor(grey,grey,grey);
  glBegin(GL_LINE_STRIP);
  glVertex3f(-1,-1,0.002);
  glVertex3f(-1,1,0.002);
  glVertex3f(1,1,0.002);
  glVertex3f(1,-1,0.002);
  glVertex3f(-1,-1,0.002);
  glEnd();
#endif
}

void glDrawBox(float x,float y,float z){
  static GLfloat n[6][3] = {
    {-1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, -1.0, 0.0},
    {0.0, 0.0, 1.0},
    {0.0, 0.0, -1.0}};
  static GLint faces[6][4] = {
    {0, 1, 2, 3},
    {3, 2, 6, 7},
    {7, 6, 5, 4},
    {4, 5, 1, 0},
    {5, 6, 2, 1},
    {7, 4, 0, 3}};
  GLfloat v[8][3];
  GLint i;

  v[0][0] = v[1][0] = v[2][0] = v[3][0] = -x / 2;
  v[4][0] = v[5][0] = v[6][0] = v[7][0] =  x / 2;
  v[0][1] = v[1][1] = v[4][1] = v[5][1] =  -y / 2;
  v[2][1] = v[3][1] = v[6][1] = v[7][1] =  y / 2;
  v[0][2] = v[3][2] = v[4][2] = v[7][2] =  -z / 2;
  v[1][2] = v[2][2] = v[5][2] = v[6][2] =  z / 2;

  for (i = 5; i >= 0; i--) {
    glBegin(GL_QUADS);
    glNormal3fv(&n[i][0]);
    glVertex3fv(&v[faces[i][0]][0]);
    glVertex3fv(&v[faces[i][1]][0]);
    glVertex3fv(&v[faces[i][2]][0]);
    glVertex3fv(&v[faces[i][3]][0]);
    glEnd();
  }
}

void glDrawDiamond(float x,float y,float z){
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(.0,.0, z);
  glVertex3f( x,.0,.0);
  glVertex3f(.0, y,.0);
  glVertex3f(-x,.0,.0);
  glVertex3f(.0,-y,.0);
  glVertex3f( x,.0,.0);
  glEnd();
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(.0,.0,-z);
  glVertex3f( x,.0,.0);
  glVertex3f(.0, y,.0);
  glVertex3f(-x,.0,.0);
  glVertex3f(.0,-y,.0);
  glVertex3f( x,.0,.0);
  glEnd();
}

void glDrawDiamond(float x,float y,float z,float dx,float dy,float dz){
  glPushMatrix();
  glTranslated(x,y,z);
  glDrawDiamond(dx,dy,dz);
  glPopMatrix();
}

void glDrawAxes(double scale){
  GLUquadric *style=gluNewQuadric();

  for(uint i=0;i<3;i++){
    glPushMatrix();
    glScalef(scale,scale,scale);
    switch(i){
    case 0:  glColor(1,0,0);  break;
    case 1:  glColor(0,1,0);  glRotatef(90,0,0,1);  break;
    case 2:  glColor(0,0,1);  glRotatef(90,0,-1,0);  break;
    }
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(.95,0,0);
    glEnd();
    glTranslatef(.8,0,0);
    glRotatef(90,0,1,0);
    gluCylinder(style,.08,0,.2,20,1); 
    glPopMatrix();
  }
  
  gluDeleteQuadric(style);
}

void drawCoordinateFrame(){
	// x-axis = green
	glBegin(GL_LINE_STRIP);
  glColor(4.0,0.0,0.0);
  glVertex3f(0.0,0.0,0.001);
  glVertex3f(1.0,0.0,0.001);
  glVertex3f(0.0,0.0,0.002);
  glVertex3f(1.0,0.0,0.002);
  glVertex3f(0.0,0.0,0.003);
  glVertex3f(1.0,0.0,0.003);
  glVertex3f(0.0,0.0,0.004);
  glVertex3f(1.0,0.0,0.004);
  glEnd();

	// y-axis = green
	glBegin(GL_LINE_STRIP);
  glColor(0.0,4.0,0.0);
  glVertex3f(0.0,0.0,0.001);
  glVertex3f(0.0,1.0,0.001);
  glVertex3f(0.0,0.0,0.002);
  glVertex3f(0.0,1.0,0.002);
  glVertex3f(0.0,0.0,0.003);
  glVertex3f(0.0,1.0,0.003);
  glVertex3f(0.0,0.0,0.004);
  glVertex3f(0.0,1.0,0.004);
  glEnd();

	// z-axis = blue
	glBegin(GL_LINE_STRIP);
  glColor(0.0,0.0,4.0);
  glVertex3f(0.001,0.0,0.0);
  glVertex3f(0.001,0.0,1.0);
  glVertex3f(0.002,0.0,0.0);
  glVertex3f(0.002,0.0,1.0);
  glVertex3f(0.003,0.0,0.0);
  glVertex3f(0.003,0.0,1.0);
  glVertex3f(0.004,0.0,0.0);
  glVertex3f(0.004,0.0,1.0);
  glEnd();
}

void glDrawSphere(float radius){
  GLUquadric *style=gluNewQuadric();
  gluSphere(style,radius,10,10); // last two value for detail
  gluDeleteQuadric(style);
}

void glDrawDisk(float radius){
  GLUquadric *style=gluNewQuadric();
  gluDisk(style,0,radius,10,1);
  gluDeleteQuadric(style);
}

void glDrawCylinder(float radius, float length,bool closed){
  GLUquadric *style=gluNewQuadric();
  glTranslatef(0,0,-length/2);
  gluCylinder(style,radius,radius,length,20,1); 
  if(closed){
    gluDisk(style,0,radius,20,1); 
    glTranslatef(0,0,length);
    gluDisk(style,0,radius,20,1); 
    glTranslatef(0,0,-length/2);
  }
  gluDeleteQuadric(style);
}

void glDrawCappedCylinder(float radius, float length){
  GLUquadric *style1=gluNewQuadric();
  GLUquadric *style2=gluNewQuadric();
  GLUquadric *style3=gluNewQuadric();

  glTranslatef(0,0,-length/2);
  gluCylinder(style1,radius,radius,length,20,1);
  glTranslatef(0,0,length);
  gluSphere(style2,radius,10,10);
  glTranslatef(0,0,-length);
  gluSphere(style3,radius,10,10);

  gluDeleteQuadric(style1);
  gluDeleteQuadric(style2);
  gluDeleteQuadric(style3);
}

void glDrawGridBox(float x=10.){
  x/=2.;
  glBegin(GL_LINE_LOOP);
  glVertex3f(-x,-x,-x);
  glVertex3f(-x,-x, x);
  glVertex3f(-x, x, x);
  glVertex3f(-x, x,-x);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f( x,-x,-x);
  glVertex3f( x,-x, x);
  glVertex3f( x, x, x);
  glVertex3f( x, x,-x);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f( x, x, x);
  glVertex3f(-x, x, x);
  glVertex3f( x,-x, x);
  glVertex3f(-x,-x, x);
  glVertex3f( x, x,-x);
  glVertex3f(-x, x,-x);
  glVertex3f( x,-x,-x);
  glVertex3f(-x,-x,-x);
  glEnd();
}

void glDrawGridBox(float x1,float y1,float z1,float x2,float y2,float z2){
  glBegin(GL_LINE_STRIP);
  glVertex3f(x1,y1,z1+0.001);
  glVertex3f(x2,y1,z1+0.001);
  glVertex3f(x2,y2,z1+0.001);
  glVertex3f(x1,y2,z1+0.001);
  glVertex3f(x1,y1,z1+0.001);
  glEnd();
 
  glBegin(GL_LINES);
   glVertex3f( x2, y2,z1 +0.001);
  glVertex3f( x2, y2,z2 +0.001);
  glEnd();
}

void glDrawKhepera(){
  GLUquadric *style=gluNewQuadric();
  glPushMatrix();
  glRotatef(-90,1,0,0);
  glColor3f (.3,.3,.3);
  gluCylinder(style,1.5,1.5,2,20,1);
  glPopMatrix();

  glColor3f (1,0,0);
  glBegin(GL_LINES);
  glVertex3f(0,2,0);
  glVertex3f(0,2,-2.5);
  glEnd();
  gluDeleteQuadric(style);
}

void glMakeSquare(int num){
  glNewList(num,GL_COMPILE);
  glColor3f(1.,0.,0.);
  glBegin( GL_LINE_LOOP );
  glVertex3f( -1 , -1 , 0. );
  glVertex3f( -1 , +1 , 0. );
  glVertex3f( +1 , +1 , 0. );
  glVertex3f( +1 , -1 , 0. );
  glEnd();
  glEndList();
}

void glMakeStdSimplex(int num){
  glNewList(num,GL_COMPILE);
  //glPolygonMode(GL_BACK,GL_FILL);
  glShadeModel(GL_SMOOTH);
  glBegin(GL_TRIANGLE_FAN);
  glColor3f (1.,1.,1.);
  glVertex3f(0.,0.,0.);
  glColor3f (1.,0.,0.);
  glVertex3f(1.,0.,0.);
  glColor3f (0.,1.,0.);
  glVertex3f(0.,1.,0.);
  glColor3f (0.,0.,1.);
  glVertex3f(0.,0.,1.);
  glColor3f (1.,0.,0.);
  glVertex3f(1.,0.,0.);
  glEnd();
  /*
    glColor4f(.5,.5,.5,.9);
    glBegin(GL_POLYGON);
    glVertex3f( 1. , 0. , 0. );
    glVertex3f( 0. , 1. , 0. );
    glVertex3f( 0. , 0. , 1. );
    glEnd();
  */
  glEndList();
}

void glMakeTorus(int num){
  glNewList(num,GL_COMPILE);

  GLint i, j, rings, sides;
  float theta1, phi1, theta2, phi2;
  float v0[03], v1[3], v2[3], v3[3];
  float t0[03], t1[3], t2[3], t3[3];
  float n0[3], n1[3], n2[3], n3[3];
  float innerRadius=0.4;
  float outerRadius=0.8;
  float scalFac;

  rings = 8;
  sides = 10;
  scalFac=1/(outerRadius*2);

  for(i=0;i<rings;i++){
    theta1 = (float)i * 2.0 * MT_PI / rings;
    theta2 = (float)(i + 1) * 2.0 * MT_PI / rings;
    for(j=0;j<sides;j++){
      phi1 = (float)j * 2.0 * MT_PI / sides;
      phi2 = (float)(j + 1) * 2.0 * MT_PI / sides;
	
      v0[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi1));
      v0[1] =-sin(theta1) * (outerRadius + innerRadius * cos(phi1));
      v0[2] = innerRadius * sin(phi1);

      v1[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi1));
      v1[1] =-sin(theta2) * (outerRadius + innerRadius * cos(phi1));
      v1[2] = innerRadius * sin(phi1);

      v2[0] = cos(theta2) * (outerRadius + innerRadius * cos(phi2));
      v2[1] =-sin(theta2) * (outerRadius + innerRadius * cos(phi2));
      v2[2] = innerRadius * sin(phi2);

      v3[0] = cos(theta1) * (outerRadius + innerRadius * cos(phi2));
      v3[1] =-sin(theta1) * (outerRadius + innerRadius * cos(phi2));
      v3[2] = innerRadius * sin(phi2);

      n0[0] = cos(theta1) * (cos(phi1));
      n0[1] =-sin(theta1) * (cos(phi1));
      n0[2] = sin(phi1);

      n1[0] = cos(theta2) * (cos(phi1));
      n1[1] =-sin(theta2) * (cos(phi1));
      n1[2] = sin(phi1);

      n2[0] = cos(theta2) * (cos(phi2));
      n2[1] =-sin(theta2) * (cos(phi2));
      n2[2] = sin(phi2);

      n3[0] = cos(theta1) * (cos(phi2));
      n3[1] =-sin(theta1) * (cos(phi2));
      n3[2] = sin(phi2);

      t0[0] = v0[0]*scalFac + 0.5;
      t0[1] = v0[1]*scalFac + 0.5;
      t0[2] = v0[2]*scalFac + 0.5;

      t1[0] = v1[0]*scalFac + 0.5;
      t1[1] = v1[1]*scalFac + 0.5;
      t1[2] = v1[2]*scalFac + 0.5;

      t2[0] = v2[0]*scalFac + 0.5;
      t2[1] = v2[1]*scalFac + 0.5;
      t2[2] = v2[2]*scalFac + 0.5;

      t3[0] = v3[0]*scalFac + 0.5;
      t3[1] = v3[1]*scalFac + 0.5;
      t3[2] = v3[2]*scalFac + 0.5;

      if((i+j)%2) glColor3f(0.,1.,0.);
      else glColor3f(0.,0.,1.);

      glBegin(GL_POLYGON);
      glNormal3fv(n3); glTexCoord3fv(t3); glVertex3fv(v3);
      glNormal3fv(n2); glTexCoord3fv(t2); glVertex3fv(v2);
      glNormal3fv(n1); glTexCoord3fv(t1); glVertex3fv(v1);
      glNormal3fv(n0); glTexCoord3fv(t0); glVertex3fv(v0);
      glEnd();
    }
  }
  glEndList();
}

uint glImageTexture(const byteA &img){
  GLuint texName;

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glGenTextures(1, &texName);

  glBindTexture(GL_TEXTURE_2D, texName);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  switch(img.d2){
  case 0:
  case 1:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img.p);
    break;
  case 2:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, img.p);
    break;
  case 3:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_RGB, GL_UNSIGNED_BYTE, img.p);
    break;
  case 4:
    glTexImage2D(GL_TEXTURE_2D, 0, 4, img.d1, img.d0, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.p);
    break;
  default:
    HALT("no image fomat");
  }
  return texName;
}

void glDrawTexQuad(uint texture,
                   float x1,float y1,float z1,float x2,float y2,float z2,
                   float x3,float y3,float z3,float x4,float y4,float z4,
                   float mulX,float mulY){
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); //GL_MODULATE);
  glBindTexture(GL_TEXTURE_2D, texture);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0,  0.0);  glVertex3f(x1,y1,z1);
  glTexCoord2f(0.0,  mulY); glVertex3f(x2,y2,z2);
  glTexCoord2f(mulX, mulY); glVertex3f(x3,y3,z3);
  glTexCoord2f(mulX, 0.0);  glVertex3f(x4,y4,z4);
  glEnd();
  glFlush();
  glDisable(GL_TEXTURE_2D);
}

#ifdef MT_FREEGLUT
/*!\brief return the RGBA-image of scenery drawn just before; the image
  buffer has to have either 2 dimensions [width,height] for a
  gray-scale luminance image or 3 dimensions [width,height,4] for an
  RGBA-image. */
void glGrabImage(byteA& image){
  if(!image.N) image.resize(glutGet(GLUT_WINDOW_HEIGHT),glutGet(GLUT_WINDOW_WIDTH),3);
  CHECK(image.nd==2 ||image.nd==3,"not an image format");
  GLint w=image.d1,h=image.d0;
  CHECK(w<=glutGet(GLUT_WINDOW_WIDTH) && h<=glutGet(GLUT_WINDOW_HEIGHT),"grabbing large image from small window");

  //glPixelStorei(GL_PACK_SWAP_BYTES,0); 
  switch(image.d2){
  case 0:
  case 1:
    glPixelTransferf(GL_RED_SCALE,.3333);
    glPixelTransferf(GL_GREEN_SCALE,.3333);
    glPixelTransferf(GL_BLUE_SCALE,.3333);
    glReadPixels(0,0,w,h,GL_LUMINANCE,GL_UNSIGNED_BYTE,image.p);
    glPixelTransferf(GL_RED_SCALE,1.);
    glPixelTransferf(GL_GREEN_SCALE,1.);
    glPixelTransferf(GL_BLUE_SCALE,1.);
    break;
  case 2:
    //glReadPixels(0,0,w,h,GL_GA,GL_UNSIGNED_BYTE,image.p);
    break;
  case 3:
    glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE,image.p);
    break;
  case 4:
#if defined MT_SunOS
    glReadPixels(0,0,w,h,GL_ABGR_EXT,GL_UNSIGNED_BYTE,image.p);
#else
#if defined MT_Cygwin
    glReadPixels(0,0,w,h,GL_RGBA,GL_UNSIGNED_BYTE,image.p);
#else
    //glReadPixels(0,0,w,h,GL_BGRA_EXT,GL_UNSIGNED_BYTE,image.p);
    glReadPixels(0,0,w,h,GL_RGBA,GL_UNSIGNED_BYTE,image.p);
#endif
#endif
    break;
  default: HALT("wrong image format");
  }
}
#else
void glGrabImage(byteA& image){ NIY; }
#endif

/*!\brief return the depth map of the scenery drawn just before; the depth
    buffer has to be a 2-dimensional [width,height] and is filled with
    depth values between 0 and 1. */
void glGrabDepth(byteA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  GLint w=depth.d1,h=depth.d0;
  glReadPixels(0,0,w,h,GL_DEPTH_COMPONENT,GL_UNSIGNED_BYTE,depth.p);
}

/*!\brief return the depth map of the scenery drawn just before; the depth
    buffer has to be a 2-dimensional [width,height] and is filled with
    depth values between 0 and 1. */
void glGrabDepth(floatA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  GLint w=depth.d1,h=depth.d0;
  glReadPixels(0,0,w,h,GL_DEPTH_COMPONENT,GL_FLOAT,depth.p);
}

void glRasterImage(int x,int y,byteA &img,float zoom){
  glRasterPos2i(x,y); //(int)(y+zoom*img.d0)); (the latter was necessary for other pixel/raster coordinates)
  glPixelZoom(zoom,-zoom);
  if(img.d1%4){ //necessary: extend the image to have width mod 4
    uint P=img.d2;
    if(!P) P=1;
    uint add=4-(img.d1%4);
    img.reshape(img.d0,img.d1*P);
    img.insColumns(img.d1,add*P);
    if(P>1) img.reshape(img.d0,img.d1/P,P);
  }

  switch(img.d2){
  case 0:
  case 1:  glDrawPixels(img.d1,img.d0,GL_LUMINANCE,GL_UNSIGNED_BYTE,img.p);        break;
  case 2:  glDrawPixels(img.d1,img.d0,GL_LUMINANCE_ALPHA,GL_UNSIGNED_BYTE,img.p);  break;
  case 3:  glDrawPixels(img.d1,img.d0,GL_RGB,GL_UNSIGNED_BYTE,img.p);              break;
  case 4:  glDrawPixels(img.d1,img.d0,GL_RGBA,GL_UNSIGNED_BYTE,img.p);             break;
  default: HALT("no image format");
  };
}

void OpenGL::watchImage(const floatA &_img,bool wait,float _zoom){
  static byteA img;
  resizeAs(img,_img);
  float x;
  for(uint i=0;i<img.N;i++){
    x=_img.elem(i);
    img.elem(i) = (x<0.)?0:((x>255.)?255:x);
  }
  watchImage(img,wait,_zoom);
}

void OpenGL::watchImage(const byteA &_img,bool wait,float _zoom){
#ifdef MT_FREEGLUT
  img=(byteA*)&_img;
  zoom=_zoom;
  
  int w=(int)(zoom*img->d1+10),h=(int)(zoom*img->d0+10);
  //set projection
  glutSetWindow(windowID);
  glutReshapeWindow(w,h);
  glViewport(0,0,w,h);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glOrtho(-5,5+zoom*img->d1,5+zoom*img->d0,-5,-1.,1.);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glDisable(GL_DEPTH_TEST);
  glClearColor(.5,.5,.5,0.);
  
  //open window
  if(wait) watch(); else update();
  img=NULL;
#else
  NIY
#endif
}

/*void glWatchImage(const floatA &x,bool wait,float zoom){
  double ma=x.max();
  double mi=x.min();
  if(wait) cout <<"watched image min/max = " <<mi <<' ' <<ma <<endl;
  byteA img;
  img.resize(x.d0*x.d1);
  img.setZero();
  for(uint i=0;i<x.N;i++){
    img(i)=(byte)(255.*(x.elem(i)-mi)/(ma-mi));
  }
  img.reshape(x.d0,x.d1);
  glWatchImage(img,wait,20);
}*/

#ifdef MT_FREEGLUT
void OpenGL::displayGrey(const arr &x,uint d0,uint d1,bool wait,uint win){
  if(!d0) d0=x.d0;
  if(!d1) d1=x.d1;
  glutSetWindow(windowID);
  double ma=x.max();
  text.clr() <<"display"<<win<<" max="<<ma<<endl;
  byteA img;
  img.resize(d0*d1);
  img.setZero();
  for(uint i=0;i<x.N;i++){
    if(x.elem(i)>0.) img(i)=(byte)(255.*x.elem(i)/ma);
    if(x.elem(i)<0.) MT_MSG("warning: negative entry");
  }
  img.reshape(d0,d1);
  watchImage(img,wait,20);
}

void OpenGL::displayRedBlue(const arr &x,uint d0,uint d1,bool wait,uint win){
  if(!d0) d0=x.d0;
  if(!d1) d1=x.d1;
  glutSetWindow(windowID);
  double mi=x.min(),ma=x.max();
  text.clr() <<"display"<<win<<" max="<<ma<<"min="<<mi<<endl;
  cout <<"\rdisplay"<<win<<" max="<<ma<<"min="<<mi;
  byteA img;
  img.resize(d0*d1,4);
  img.setZero();
  for(uint i=0;i<x.N;i++){
    if(x.elem(i)>0.) img(i,0)=(byte)(255.*x.elem(i)/ma);
    if(x.elem(i)<0.) img(i,2)=(byte)(255.*x.elem(i)/mi);
  }
  img.reshape(d0,d1,4);
  watchImage(img,wait,20);
}
#endif

void glDrawUI(void *p){
  glPushName(0x10);
  ((glUI*)p)->glDraw();
  glPopName();
}

bool glHoverUI(void *p,OpenGL *gl){
  //bool b=
  ((glUI*)p)->checkMouse(gl->mouseposx,gl->mouseposy);
  //if(b) glutPostRedisplay();
  return true;
}

bool glClickUI(void *p,OpenGL *gl){
  bool b=((glUI*)p)->checkMouse(gl->mouseposx,gl->mouseposy);
  if(b) gl->update();
  int t=((glUI*)p)->top;
  if(t!=-1){
    cout <<"CLICK! on button #" <<t <<endl;
    MTexitLoop(); //glutLeaveMainLoop();
    return false;
  }
  return true;
}

#ifdef MT_FREEGLUT
void glSelectWin(uint win){
  if(!staticgl[win]) staticgl[win]=new OpenGL;
  glutSetWindow(staticgl[win]->windowID);
}
#endif

#else //! MT_GL
void glColor(float,float,float,float){ NIY; }
void glDrawDiamond(float,float,float,float,float,float){ NIY; }
void glStandardLight(){ NIY; }
void glStandardScene(void*){ NIY; };
#endif


//===========================================================================
//
// standalone draw routines for larget data structures
//

#ifdef MT_GL
void glDrawDots(void *dots){ glDrawDots(*(arr*)dots); }

void glDrawDots(arr& dots){
  CHECK(dots.nd==2 && dots.d1==3,"wrong dimension");
  glBegin(GL_POINTS);
  for(uint i=0;i<dots.d0;i++) glVertex3dv(&dots(i,0));
  glEnd();
}
#endif

//===========================================================================
//
// OpenGL implementations
//

#ifdef MT_FREEGLUT
  //! constructor
OpenGL::OpenGL(const char* title,int w,int h,int posx,int posy){
  init();

  if(!nrWins){
    int argc=1;
    char *argv[1]={(char*)"x"};
    glutInit(&argc, argv);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_CONTINUE_EXECUTION);
  }
  nrWins++;

  glutInitWindowSize(w,h);
  glutInitWindowPosition(posx,posy);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

  windowID = glutCreateWindow(title);

  //OpenGL initialization
  //two optional thins:
  glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE); glFrontFace(GL_CW);
  glDisable(GL_CULL_FACE);
  glDepthFunc(GL_LESS);
  glShadeModel(GL_SMOOTH);
  glShadeModel(GL_FLAT);

  if(glwins.N<(uint)windowID+1) glwins.resizeCopy(windowID+1);
  glwins(windowID) = this;

  glutDisplayFunc( _Draw );
  glutKeyboardFunc( _Key );
  glutMouseFunc ( _Mouse ) ;
  glutMotionFunc ( _Motion ) ;
  glutPassiveMotionFunc ( _PassiveMotion ) ;
  glutCloseFunc ( _Close ) ;
  glutReshapeFunc( _Reshape );
  glutSpecialFunc( _Special );
  glutMouseWheelFunc ( _MouseWheel ) ;

  //  glutVisibilityFunc( Visibility );
  //  glutKeyboardUpFunc( KeyUp );
  //  glutSpecialUpFunc( SpecialUp );
  //  glutJoystickFunc( Joystick, 100 );
  //  glutEntryFunc ( Entry ) ;
}
// freeglut destructor
OpenGL::~OpenGL(){
  glutDestroyWindow(windowID);
  glwins(windowID)=0;
  nrWins--;
  if(!nrWins) fgDeinitialize();
  delete WS;
}
#endif


#ifdef MT_QTGLUT
OpenGL::OpenGL(const char* title,int width,int height,int posx,int posy)
  :QGLWidget(QGLFormat(GLformat)){
  QGLWidget::move(posx,posy);
  QGLWidget::resize(width,height);
  QWidget::setMouseTracking(true);
  QWidget::setWindowTitle(title);
  init();
  windowID=(int)winId();
}
  //! Qt constructor when window is parent of another window
OpenGL::OpenGL(QWidget *parent,const char* title,int width,int height,int posx,int posy)
  :QGLWidget(QGLFormat(GLformat),parent){
  QGLWidget::move(posx,posy);
  QGLWidget::resize(width,height);
  QWidget::setMouseTracking(true);
  QWidget::setWindowTitle(title);
  init();
  windowID=(int)winId();
}
  //! destructor
OpenGL::~OpenGL(){
  if(osContext) delete osContext;
  if(osPixmap) delete osPixmap;
  delete WS; 
};
#endif

void OpenGL::clone(const OpenGL& gl){
  *WS = *gl.WS;
  drawers = gl.drawers;
  camera = gl.camera;
    /*float clearR,clearG,clearB,clearA;  //!< colors of the beackground (called in glClearColor(...))
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
  double zoom;*/
}

void OpenGL::init(){
  WS = new OpenGLWorkspace;

  camera.setPosition(0.,0.,10.);
  camera.focus(0,0,0);
  camera.setZRange(.1,1000.);
  camera.setHeightAngle(12.);

  drawFocus=false;
  clearR=clearG=clearB=1.; clearA=0.;
  drawers.memMove=true;
  mouseposx=mouseposy=0;
  mouse_button=0;
  mouseIsDown=false;
  
  reportEvents=false;
  reportSelects=false;
  selectOnHover=false;
  immediateExitLoop=false;
  exitkeys="";

  img=NULL;
  zoom=1;
#ifdef MT_QTGLUT
  osPixmap=0;
  osContext=0;
  quitLoopOnTimer=reportEvents=false;
#endif
};

//! add a draw routine
void OpenGL::add(void (*call)(void*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  GLDrawer d; d.classP=(void*)classP; d.call=call;
  drawers.append(d);
}

//! remove a draw routine
void OpenGL::remove(void (*call)(void*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  uint i;
  for(i=0;i<drawers.N;i++) if(drawers(i).call==call && drawers(i).classP==classP) break;
  CHECK(i<drawers.N,"value to remove not found");
  drawers.remove(i,1);
}

  //! clear the list of all draw routines
void OpenGL::clear(){
  drawers.resize(0);
}

  //! add a hover callback
void OpenGL::addHoverCall(bool (*call)(void*,OpenGL*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  GLHoverCall c; c.classP=(void*)classP; c.call=call;
  hoverCalls.append(c);
}

  //! clear hover callbacks
void OpenGL::clearHoverCalls(){
  hoverCalls.clear();
}

  //! add a click callback
void OpenGL::addClickCall(bool (*call)(void*,OpenGL*),const void* classP){
  CHECK(call!=0,"OpenGL: NULL pointer to drawing routine");
  GLClickCall c; c.classP=(void*)classP; c.call=call;
  clickCalls.append(c);
}

  //! clear click callbacks
void OpenGL::clearClickCalls(){
  clickCalls.clear();
}

#ifdef MT_GL
void OpenGL::Draw(int w,int h,ors::Camera *cam){
  if(img){ //raster an image
    glClear(GL_COLOR_BUFFER_BIT);
    glRasterImage(0,0,*img,zoom);
#ifdef MT_FREEGLUT
    glutSwapBuffers();
#else
    swapBuffers();
#endif
    return;
  }
  //clear bufferer
  glViewport(0,0,w,h);
  glClearColor(clearR,clearG,clearB,clearA);
  if(drawers.N)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  //OpenGL initialization
  //two optional things:
  glEnable(GL_DEPTH_TEST);
  //glEnable(GL_CULL_FACE); glFrontFace(GL_CCW); //CCW is default!
  glDepthFunc(GL_LESS);
  //glShadeModel(GL_SMOOTH);

  //projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if(!cam) camera.glSetProjectionMatrix();
  else     cam->glSetProjectionMatrix();
  //glLineWidth(2);

  //draw focus?
  if(drawFocus){
    glColor(1.,.7,.3);
    double size = .005 * (camera.f->p-*camera.foc).length();
    glDrawDiamond((*camera.foc)(0),(*camera.foc)(1),(*camera.foc)(2),size,size,size);
  }

  //std color: black:
  glColor(.3,.3,.5);

  //draw objects
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if(!drawers.N){ MT_MSG("OpenGL: nothing to be drawn -- add draw routines!"); if(!text.N()) text <<"<nothing to draw>"; }
  for(uint i=0;i<drawers.N;i++) (*drawers(i).call)(drawers(i).classP);

  //check matrix stack
  GLint s;
  glGetIntegerv(GL_MODELVIEW_STACK_DEPTH,&s);
  CHECK(s==1,"OpenGL matrix stack has not depth 1 (pushs>pops)");

  //draw text
  glGetIntegerv(GL_RENDER_MODE,&s);
  if(text.N() && s!=GL_SELECT){
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glOrtho(0.,(double)w,(double)h,.0,-1.,1.);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(clearR+clearG+clearB>1.){
      glColor(0.0,0.0,0.0,1.0); // clear text  5. Mar 06 (hh)
    }else{
      glColor(1.0,1.0,1.0,1.0); // white colored text,  5. Mar 06 (hh)
    }
    glDrawText(text,10,20,0);
  }
}

void OpenGL::Select(){
  uint i,k;
  int j;

  glSelectBuffer(1000,selectionBuffer);
  glRenderMode(GL_SELECT);

  //projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  GLint viewport[4];
  viewport[0]=0; viewport[1]=0; viewport[2]=width(); viewport[3]=height();
  gluPickMatrix((GLdouble) mouseposx, (GLdouble) (height()-mouseposy), 2., 2., viewport);
  camera.glSetProjectionMatrix();

  //draw objects
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glInitNames();
  for(i=0;i<drawers.N;i++){
    glLoadName(i);
    (*drawers(i).call)(drawers(i).classP);
    glGetIntegerv(GL_NAME_STACK_DEPTH,&j);
    CHECK(j==0,"OpenGL name stack has not depth 1 (pushs>pops)");
  }

  GLint n;
  n=glRenderMode(GL_RENDER);
  selection.resize(n);

  GLuint *obj,maxD=(GLuint)(-1);
  topSelection=0;
  for(j=0,i=0;i<(uint)n;i++){
    obj=selectionBuffer+j;
    j+=3+obj[0];

    //get name as superposition of all names
    selection(i).name = 0;
    for(k=0;k<obj[0];k++) selection(i).name |= obj[3+k];

    //get dim and dmax
    selection(i).dmin=(double)obj[1]/maxD;  camera.glConvertToTrueDepth(selection(i).dmin);
    selection(i).dmax=(double)obj[2]/maxD;  camera.glConvertToTrueDepth(selection(i).dmax);

    //get top-most selection
    if(!topSelection || selection(i).dmin < topSelection->dmin) topSelection = &selection(i);
  }
  if(reportSelects) reportSelection();
}
#endif

  //! update the view (in Qt: also starts displaying the window)
bool OpenGL::update(const char *txt){
  pressedkey=0;
  if(txt) text.clr() <<txt;
#ifdef MT_FREEGLUT
  glutSetWindow(windowID);
  glutPostRedisplay();
#endif
#ifdef MT_QTGLUT
  show();
  QGLWidget::update();
#endif
  MTprocessEvents();
  return !pressedkey;
}

  /*!\brief watch in interactive mode and wait for an exiting event
    (key pressed or right mouse) */
int OpenGL::watch(const char *txt){
  update(txt);
  MTenterLoop();
  MTprocessEvents();
  return pressedkey;
}

  //! waits some msecons before updating
int OpenGL::timedupdate(double sec){
  static double lasttime=-1;
  double now;
  now=MT::realTime();
  if(lasttime>0. && now-lasttime<sec) MT::wait(lasttime+sec-now);
  lasttime=now;
  return update();
#if 0//def MT_QTGLUT
  int i;
  quitLoopOnTimer=true;
  i=startTimer(msec);
  MTenterLoop();
  killTimer(i);
  return update();
#endif
}

  //! resize the window
void OpenGL::resize(int w,int h){
#ifdef MT_FREEGLUT
  glutSetWindow(windowID);
  glutReshapeWindow(w,h);
#elif defined MT_QTGLUT
  QGLWidget::resize(w,h);
#endif
  MTprocessEvents();
}

  //! set the four clear colors
void OpenGL::setClearColors(float r,float g,float b,float a){
  clearR=r; clearG=g; clearB=b; clearA=a;
}

  /*!\brief inverse projection: given a 2D+depth coordinates in the
    camera view (e.g. as a result of selection) computes the world 3D
    coordinates */
void OpenGL::unproject(double &x,double &y,double &z){
#ifdef MT_GL
  double _x,_y,_z;
  GLdouble modelMatrix[16],projMatrix[16];
  GLint viewPort[4];
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  camera.glSetProjectionMatrix();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
  glGetIntegerv(GL_VIEWPORT,viewPort);
  gluUnProject(x,y,z, modelMatrix, projMatrix, viewPort, &_x, &_y, &_z);
  x=_x; y=_y; z=_z;
#else
  NIY;
#endif
}


void OpenGL::capture(byteA &img,int w,int h,ors::Camera *cam){
#ifdef MT_FREEGLUT
  glutSetWindow(windowID);
  glutPostRedisplay();
  MTprocessEvents();
  Draw(w,h,cam);
  img.resize(h,w,3);
  glGrabImage(img);
#endif
}

void OpenGL::captureStereo(byteA &imgL,byteA &imgR,int w,int h,ors::Camera *cam,double baseline){
#ifdef MT_FREEGLUT
  glutSetWindow(windowID);
  glutPostRedisplay();
  MTprocessEvents();
  Draw(w,h,cam);
  imgR.resize(h,w,3);
  glGrabImage(imgR);
  double xorg=cam->f->p(0);
  cam->f->p(0) -= baseline;
  Draw(w,h,cam);
  imgL.resize(h,w,3);
  glGrabImage(imgL);
  cam->f->p(0) = xorg;
#endif
}

#ifdef MT_QTGLUT
  /*!\brief creates a off-screen rendering context for future backround
      rendering routines -- the off-screen context cannot be
      resized... */
void OpenGL::createOffscreen(int width,int height){
  if(osContext && (width>osPixmap->width() || height>osPixmap->height())){
    delete osContext;
    delete osPixmap;
    osContext=0;
  }
  if(!osContext){
    osPixmap=new QPixmap(width,height);
    if(!osPixmap) MT_MSG("can't create off-screen Pixmap");
    osContext=new QGLContext(QGLFormat(GLosformat),osPixmap);
    if(!osContext->create()) MT_MSG("can't create off-screen OpenGL context");
  }
}

  /*!\brief return the RGBA-image of the given perspective; rendering is done
      off-screen (on an internal QPixmap) */
void OpenGL::offscreenGrab(byteA& image){
  if(image.nd==3){ CHECK(image.d2==4,"3rd dim of image has to be 4 for RGBA");}
  else{ CHECK(image.nd==2,"image has to be either 2- or 3(for RGBA)-dimensional");}
  setOffscreen(image.d1,image.d0);
  Draw(image.d1,image.d0);
  glGrabImage(image);
}

  /*!\brief return the RGBA-image of the given perspective; rendering
      is done off-screen (on an internal QPixmap) */
void OpenGL::offscreenGrab(byteA& image,byteA& depth){
  if(image.nd==3){ CHECK(image.d2==4,"3rd dim of image has to be 4 for RGBA");}
  else{ CHECK(image.nd==2,"image has to be either 2- or 3(for RGBA)-dimensional");}
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  setOffscreen(image.d1,image.d0);
  Draw(image.d1,image.d0);
  glGrabImage(image);
  glGrabDepth(depth);
}

  /*!\brief return only the depth gray-scale map of given perspective;
      rendering is done off-screen (on an internal QPixmap) */
void OpenGL::offscreenGrabDepth(byteA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  setOffscreen(depth.d1,depth.d0);
  Draw(depth.d1,depth.d0);
  glGrabDepth(depth);
}

  /*!\brief return only the depth gray-scale map of given perspective;
      rendering is done off-screen (on an internal QPixmap) */
void OpenGL::offscreenGrabDepth(floatA& depth){
  CHECK(depth.nd==2,"depth buffer has to be either 2-dimensional");
  setOffscreen(depth.d1,depth.d0);
  Draw(depth.d1,depth.d0);
  glGrabDepth(depth);
}

void OpenGL::setOffscreen(int width,int height){
  createOffscreen(width,height);
  CHECK(width<=osPixmap->width() && height<=osPixmap->height(),
	"width ("<<width<<") or height ("<<height
	<<") too large for the created pixmap - create and set size earlier!");
  osContext->makeCurrent();
  //if(initRoutine) (*initRoutine)();
}
#endif


  //! print some info on the selection buffer
void OpenGL::reportSelection(){
  uint i;
  std::cout <<"selection report: mouse=" <<mouseposx <<" " <<mouseposy <<" -> #selections=" <<selection.N <<std::endl;
  for(i=0;i<selection.N;i++){
    if(topSelection == &selection(i)) std::cout <<"  TOP: "; else std::cout <<"       ";
    std::cout
      <<"name = 0x" <<std::hex <<selection(i).name <<std::dec
      <<" min-depth:" <<selection(i).dmin <<" max-depth:" <<selection(i).dmax
      <<endl;
  }
}

#ifdef MT_FREEGLUT
int OpenGL::width(){  glutSetWindow(windowID); return glutGet(GLUT_WINDOW_WIDTH); }
int OpenGL::height(){ glutSetWindow(windowID); return glutGet(GLUT_WINDOW_HEIGHT); }
#endif
#ifdef MT_QTGLUT
int OpenGL::width(){  return QGLWidget::width(); }
int OpenGL::height(){ return QGLWidget::height(); }
#endif
#ifndef MT_GL
  //! get width
int OpenGL::width(){ return 0; }
  //! get height
int OpenGL::height(){ return 0; }
#endif

#ifdef MT_GL2PS
  /*!\brief generates a ps from the current OpenGL display, using gl2ps */
void OpenGL::saveEPS(const char *filename){
  FILE *fp = fopen(filename, "wb");
  GLint buffsize = 0, state = GL2PS_OVERFLOW;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  while(state==GL2PS_OVERFLOW){
    buffsize+=1024*1024;
    gl2psBeginPage( "Marc Toussaint", "MT", viewport,
		    GL2PS_EPS, GL2PS_BSP_SORT, GL2PS_SILENT |
		    GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
		    GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
		    GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
		    fp, filename );
    Draw(width(),height());
    state = gl2psEndPage();
  }
  fclose(fp);
}
#else
void OpenGL::saveEPS(const char*){
  MT_MSG("WARNING: OpenGL::saveEPS was called without MT_GL2PS configured!");
}
#endif

#ifdef MT_QTGLUT
void OpenGL::about(std::ostream& os){
  os <<"Widget's OpenGL capabilities:\n";
  QGLFormat f=format();
  os <<"direct rendering: " <<f.directRendering() <<"\n"
     <<"double buffering: " <<f.doubleBuffer()  <<"\n"
     <<"depth:            " <<f.depth() <<"\n"
     <<"rgba:             " <<f.rgba() <<"\n"
     <<"alpha:            " <<f.alpha() <<"\n"
     <<"accum:            " <<f.accum() <<"\n"
     <<"stencil:          " <<f.stencil() <<"\n"
     <<"stereo:           " <<f.stereo() <<"\n"
     <<"overlay:          " <<f.hasOverlay() <<"\n"
     <<"plane:            " <<f.plane() <<std::endl;

  if(!osContext){
    os <<"no off-screen context created yet" <<std::endl;
  }else{
    os <<"Off-screen pixmaps's OpenGL capabilities:\n";
    f=osContext->format();
    os <<"direct rendering: " <<f.directRendering() <<"\n"
       <<"double buffering: " <<f.doubleBuffer()  <<"\n"
       <<"depth:            " <<f.depth() <<"\n"
       <<"rgba:             " <<f.rgba() <<"\n"
       <<"alpha:            " <<f.alpha() <<"\n"
       <<"accum:            " <<f.accum() <<"\n"
       <<"stencil:          " <<f.stencil() <<"\n"
       <<"stereo:           " <<f.stereo() <<"\n"
       <<"overlay:          " <<f.hasOverlay() <<"\n"
       <<"plane:            " <<f.plane() <<std::endl;
  }
}
#else
  /*!\brief report on the OpenGL capabilities (the QGLFormat) */
void OpenGL::about(std::ostream& os){ MT_MSG("NIY"); }
#endif


//===========================================================================
//
// modified main loop routine
// only change to free glut: don't deinitialize the glut engine!
//

#if defined MT_FREEGLUT //&& !defined MT_Cygwin
extern "C"{
#include "opengl_freeglut.c"
}
#endif


//===========================================================================
//
// callbacks
//

#if 1
#  define CALLBACK_DEBUG(x) if(reportEvents) x
#else
#  define CALLBACK_DEBUG(x)
#endif

void getSphereVector(ors::Vector& vec,int _x,int _y,int width,int height){
  int minwh = width<height?width:height;
  double x,y;
  x=(double)_x;  x=x-.5*width;  x*= 2./minwh;
  y=(double)_y;  y=y-.5*height; y*=-2./minwh;
  vec(0)=x;
  vec(1)=y;
  vec(2)=.5-(x*x+y*y);
  if(vec(2)<0.) vec(2)=0.;
}

void OpenGL::Reshape(int width, int height){
  CALLBACK_DEBUG(printf("Window %d Reshape Callback:  %d %d\n", windowID, width, height ));
  camera.setWHRatio((double)width/height);
  //update();
}

void OpenGL::Key(unsigned char key, int x, int y){
  CALLBACK_DEBUG(printf("Window %d Keyboard Callback:  %d (`%c') %d %d\n", windowID, key, (char)key, x, y ));
  pressedkey=key;
  if(key==13 || key==32 || key==27) MTexitLoop();
  if(MT::contains(exitkeys,key)) MTexitLoop();
}

void OpenGL::Mouse(int button, int updown, int _x, int _y){
  CALLBACK_DEBUG(printf("Window %d Mouse Click Callback:  %d %d %d %d\n", windowID, button, updown, _x, _y ));
  mouse_button=1+button;
  if(updown) mouse_button=-1-mouse_button;
  mouseposx=_x; mouseposy=_y;
  ors::Vector vec;
  getSphereVector(vec,_x,_y,width(),height());
  lastEvent.set(mouse_button,-1,_x,_y,0.,0.);

  if(!updown){
    CHECK(!mouseIsDown,"I thought the mouse is up...");
    mouseIsDown=true;
    drawFocus=true;
  }else{
    CHECK(mouseIsDown,"I thought the mouse is down...");
    mouseIsDown=false;
    drawFocus=false;
  }
  //store where you've clicked
  WS->downVec=vec;
  WS->downRot=camera.f->r;
  WS->downPos=camera.f->p;
  WS->downFoc=*camera.foc;

  //check object clicked on
  if(!updown){
    if(reportSelects) Select();
  }
  //step through all callbacks
  bool cont=true;
  if(!updown){
    for(uint i=0;i<clickCalls.N;i++) cont=cont && (*clickCalls(i).call)(clickCalls(i).classP,this);
  }
  update();
}

void OpenGL::Motion(int _x, int _y){
#ifdef MT_GL
  CALLBACK_DEBUG(printf("Window %d Mouse Motion Callback:  %d %d\n", windowID, _x, _y ));
  mouseposx=_x; mouseposy=_y;
  ors::Vector vec;
  getSphereVector(vec,_x,_y,width(),height());
  lastEvent.set(mouse_button,-1,_x,_y,vec(0)-WS->downVec(0),vec(1)-WS->downVec(1));
#ifndef MT_Linux
    int modifiers=glutGetModifiers();
#else
    //int modifiers=0;
#endif
  CHECK(mouseIsDown,"I thought the mouse is down...");
  if(mouse_button==1){ // && !(modifiers&GLUT_ACTIVE_SHIFT) && !(modifiers&GLUT_ACTIVE_CTRL)){
    ors::Quaternion rot;
    if(WS->downVec(2)<.1){
      rot.setDiff(vec,WS->downVec);  //consider imagined sphere rotation of mouse-move
    }else{
      rot.setVec((vec-WS->downVec) ^ ors::Vector(0,0,1)); //consider only xy-mouse-move
    }
    camera.f->r = WS->downRot * rot;   //rotate camera's direction
    rot = WS->downRot * rot / WS->downRot; //interpret rotation relative to current viewing
    camera.f->p = rot * WS->downPos;   //rotate camera's position
    //camera.f->r = rot * WS->downRot;   //rotate camera's direction
    //camera.f->p = WS->downFoc + (camera.f->r/WS->downRot)* (WS->downPos - WS->downFoc);   //rotate camera's position
    //camera.focus();
    update();
    if(immediateExitLoop) MTexitLoop();
  }
  if(mouse_button==3){ // || (mouse_button==1 && (modifiers&GLUT_ACTIVE_SHIFT) && !(modifiers&GLUT_ACTIVE_CTRL))){
    ors::Vector trans = WS->downVec - vec;
    trans(2)=0.;
    trans = WS->downRot*trans;
    camera.f->p = WS->downPos + trans;
    update();
  }
  if(mouse_button==2){ // || (mouse_button==1 && !(modifiers&GLUT_ACTIVE_SHIFT) && (modifiers&GLUT_ACTIVE_CTRL))){
    double dy = WS->downVec(1) - vec(1);
    if(dy<-.99) dy = -.99;
    camera.f->p = WS->downPos + WS->downRot*ors::Vector(0,0,1) * dy * WS->downPos.length();
    update();
  }
#else
  NIY;
#endif
}

void OpenGL::PassiveMotion(int x, int y){
  CALLBACK_DEBUG(printf("Window %d Mouse Passive Motion Callback:  %d %d\n", windowID, x, y ));
  static int calls=0;
  if(calls) return;
  calls++;
  mouseposx=x; mouseposy=y;
  if(selectOnHover) Select();
  bool ud=false;
  for(uint i=0;i<hoverCalls.N;i++) ud=ud || (*hoverCalls(i).call)(hoverCalls(i).classP,this);
  if(ud) update();
  calls--;
}

void OpenGL::Special(int key, int x, int y){}
void OpenGL::MouseWheel(int wheel, int direction, int x, int y){}

//===========================================================================
//
// GUI implementation
//

void glUI::addButton(uint x,uint y,const char *name,const char *img1,const char *img2){
  Button &b = buttons.append();
  byteA img;
  b.hover=false;
  b.x=x; b.y=y; b.name=name;
  //read_png(img,tex1);
  if(img1){
    read_ppm(img,img1,true);
  }else{
    img.resize(18,strlen(name)*9+10,3);
    img=255;
  }
  b.w=img.d1; b.h=img.d0;
  b.img1=img;    add_alpha_channel(b.img1,100);
  if(img2){
    read_ppm(img,img1,true);
    CHECK(img.d1==b.w && img.d0==b.h,"mismatched size");
  }
  b.img2=img;    add_alpha_channel(b.img2,200);
}

void glUI::glDraw(){
#ifdef MT_GL
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  GLint viewPort[4];
  glGetIntegerv(GL_VIEWPORT,viewPort);
  glOrtho(0.,viewPort[2],viewPort[3],.0,-1.,1.);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
    
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  Button *b;
  float x,y,w,h;
  for(uint i=0;i<buttons.N;i++){
    b = &buttons(i);
    x=b->x-b->w/2.;
    y=b->y-b->h/2.;
    w=b->w;
    h=b->h;
    glColor(0,0,0,1);
    glDrawText(b->name,x+5,y+h-5,0.);
    if((int)i==top) glRasterImage((int)x,(int)y,b->img2);
    else	    glRasterImage((int)x,(int)y,b->img1);
  }
#else
  NIY;
#endif
}

bool glUI::checkMouse(int _x,int _y){
  float x,y,w,h;
  Button *b;
  int otop=top;
  top=-1;
  for(uint i=0;i<buttons.N;i++){
    b = &buttons(i);
    x=b->x-b->w/2.;
    y=b->y-b->h/2.;
    w=b->w;
    h=b->h;
    if(_x>=x && _x <=x+w && _y>=y && _y<=y+h) top = i;
  }
  if(otop==top) return false;
  //glutPostRedisplay();
  return true;
}

#ifdef MT_QTGLUT
#if   defined MT_MSVC
#  include"opengl_MSVC.moccpp"
#elif defined MT_SunOS
#  include"opengl_SunOS.moccpp"
#elif defined MT_Linux
#  include"opengl_moc.cpp"
#elif defined MT_Cygwin
#  include"opengl_Cygwin.moccpp"
#endif
#endif


#ifndef MT_GL
OpenGL::OpenGL(const char* title,int w,int h,int posx,int posy){
  MT_MSG("WARNING - creating dummy OpenGL");
  init();
}
OpenGL::~OpenGL(){ delete WS; }
void OpenGL::Draw(int w,int h, ors::Camera*){}
void OpenGL::Select(){}
void OpenGL::watchImage(const byteA &_img,bool wait,float _zoom){}
#endif
