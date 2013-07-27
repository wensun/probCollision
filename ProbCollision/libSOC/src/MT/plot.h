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

#ifndef MT_plot_h
#define MT_plot_h

//===========================================================================

class OpenGL;
struct Gaussian;
namespace MT{  template<class T> class Array;  }

typedef unsigned int uint;
typedef double real;
typedef MT::Array<real> arr;
typedef MT::Array<Gaussian> GaussianA;
typedef MT::Array<Gaussian*> GaussianL;


//===========================================================================

typedef enum{ opengl, xfig, gnupl } PlotMode;

struct PlotModuleWorkspace;
struct PlotModule{
  PlotModuleWorkspace *WS;
  PlotMode mode;
  OpenGL *gl;
  bool light,grid,colors,drawBox,drawDots;
  uint thickLines;//display options
  PlotModule();
  ~PlotModule();
};
extern PlotModule plotModule;

//===========================================================================

void plotGnuplot();
void plotOpengl();
void plotOpengl(bool threeD,real xl,real xh,real yl=-1.,real yh=1.,real zl=-1.,real zh=1.);

void plot(bool wait=true);
void plotClear();
void plotFunction(const arr& f,real x0=0.,real x1=0.);
void plotFunctions(const arr& f,real x0=0.,real x1=0.);
void plotFunction(const arr& x,const arr& f);
void plotSurface(const arr& X);
void plotArray(const arr& X);
void plotPoint(real x,real y,real z);
void plotPoint(const arr& x);
void plotPoints(const arr& X);
void plotClearPoints();
void plotLine(const arr& X);
void plotPoints(const arr& X,const arr& Y);
void writeGnuplotFiles();
void plotCovariance(const arr& mean,const arr& cov);
void plotVectorField(const arr& X,const arr& dX);
void plotVectorField(arr& dX);
void plotGaussians(const GaussianA& G);
void plotGaussians(const GaussianL& G);

void glDrawPlot(void *module);

#ifdef MT_IMPLEMENTATION
#  include"plot.cpp"
#endif

#endif

