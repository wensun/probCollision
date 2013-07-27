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

#ifndef MT_IMPLEMENT_TEMPLATES
#  define MT_IMPLEMENT_TEMPLATES
#endif

#include <MT/soc.h>
#include <MT/util.h>
#include <MT/opengl.h>
#include <MT/array.h>
#include <MT/ors.h>

//explicit template instantiations

//-- from util.h
template void MT::getParameter(double&,const char*);
template void MT::getParameter(uint&,const char*);
template void MT::getParameter(int&,const char*);
template void MT::getParameter(int&,const char*,const int&);
template void MT::getParameter(bool&,const char*,const bool&);
template int  MT::getParameter(const char*,const int&);
template uint MT::getParameter(const char*,const uint&);
template bool MT::getParameter(const char*,const bool&);
template MT::String MT::getParameter(const char*);
template MT::String MT::getParameter(const char*,const MT::String&);
template void MT::save<uintA>(const uintA&,const char*);

template void MT::Parameter<MT::String>::initialize();
template void MT::Parameter<bool>::initialize();
template void MT::Parameter<double>::initialize();
template void MT::Parameter<int>::initialize();

//-- from array.h
#define T double
#  include <MT/array_instantiate.cpp>
#undef T
#define NOFLOAT
#define T uint
#  include <MT/array_instantiate.cpp>
#undef T
#define T int
#  include <MT/array_instantiate.cpp>
#undef T
#define T byte
#  include <MT/array_instantiate.cpp>
#undef T
#undef NOFLOAT

//-- from ors.h
template MT::Array<ors::Body*>::Array();
template MT::Array<ors::Body*>::~Array();
template MT::Array<ors::Joint*>::Array();
template MT::Array<ors::Joint*>::~Array();
template MT::Array<ors::Shape*>::Array();
template MT::Array<ors::Shape*>::~Array();
template MT::Array<ors::Proxy*>::Array();
template MT::Array<ors::Proxy*>::~Array();
template MT::Array<ors::Frame*>::Array();
template MT::Array<ors::Frame*>::~Array();
template MT::Array<TaskVariable*>::Array(uint);
template MT::Array<TaskVariable*>::~Array();
template TaskVariable*& MT::Array<TaskVariable*>::operator()(uint) const;

template MT::Array<MT::Array<uint> >::~Array();

template MT::Array<Any*>::Array();
template MT::Array<Any*>::~Array();

