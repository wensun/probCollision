#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include <cfloat>
#include <math.h>
#include <map>
#include <vector>
#include <list>
#include <deque>

#include "matrix.h"

#include "callisto.h"
#include "glut.h"

#define dmax DBL_MAX
#define dmin DBL_MIN
#define deps 1e-08

enum {FRONT, BACK, INTERSECTING, COINCIDENT};

//*****************************************************************************
struct normalDist {
	double x, y;
	Matrix<2,2> M;
};

//*****************************************************************************
struct Vector2d
{
public:  
	double x;
	double y;

public: 
	Vector2d():x(0),y(0){}
	Vector2d(double px,double py):x(px),y(py){};
	~Vector2d(){}
};

bool operator< (const Vector2d& v1, const Vector2d& v2); 
const Vector2d operator + (const Vector2d& v1, const Vector2d& v2);
const Vector2d operator - (const Vector2d& v1, const Vector2d& v2);
const Vector2d operator - (const Vector2d& v);
const double operator * (const Vector2d& v1, const Vector2d& v2);
const Vector2d operator * (const double& d, const Vector2d& v);
const Vector2d operator * (const Vector2d& v, const double& d);
const Vector2d operator / (const Vector2d& v, const double& d);
double det(const Vector2d& v1, const Vector2d& v2);
double norm(const Vector2d &v);
double sqnorm(const Vector2d& v);
Vector2d normal(const Vector2d& v1, const Vector2d& v2);
Vector2d dir(const Vector2d & v1, const Vector2d& v2);
Vector2d bisector(const Vector2d &v1, const Vector2d &v2);
double distSqPointSegment(const Vector2d& p1, const Vector2d& p2, const Vector2d& x);
double signDistance(const Vector2d& p1, const Vector2d& p2, const Vector2d& x);
//*****************************************************************************


//*****************************************************************************
struct bspLine
{
public:
	Vector2d n; //normal to split line
	double d; //nx*x+ny*y = d line equation

public: 
	//signed distance from line
	double signDistance(const Vector2d &p) const;
	//signed intersection of line with an edge
	Vector2d intersectLine(const Vector2d &v1,const Vector2d &v2) const;
};
//*****************************************************************************


class bspTree;
class bspNode;

//*****************************************************************************
class bspEdge
{
private:
	static bspTree* btree;
	int e0, e1;
	Vector2d dv;
	bspLine line;
	Vector2d mid;
	double len;

public: 
	~bspEdge() {};
	bspEdge(bspTree * tree) { btree = tree; };
	bspEdge(bspTree * tree, const int v0, const int v1);
	void addEdge(const int v0, const int v1);
	Vector2d* getVertex(int id) const;
	Vector2d* getDirVec() { return &dv; };
	Vector2d* getMidPt() { return &mid; };
	double getLength() { return len; };
	bspLine* getLine() { return &line; };
	int getVertexId(int id) const;
	int checkAgainstLine(const bspLine &line);
};
//*****************************************************************************

//*****************************************************************************
class bspNode
{
private:
	static bspTree* btree;
	bspNode* front;
	bspNode* back;
	bspLine line;
	int level;
	std::vector<int> evec;
	int sedge;

public: 
	bspNode(bspTree * tree):front(0),back(0),level(0){
		btree = tree;
		evec.clear();
	};
	bspNode(bspNode * parent):front(0),back(0),level(parent->level+1){
		btree = parent->btree;
		evec.clear();
	};

	~bspNode(){ 
		delete front; delete back;
		evec.clear();
	}

	bool isLeaf() { return (front==0 && back==0); };
	bspNode* getFront() { return front; };
	bspNode* getBack () { return back; };
	bspLine* getLine() { return &line; };
	void setSplitEdge(int se) { sedge = se; };
	int getSplitEdge() { return sedge; };
	int getLevel() { return level; };
	int getNumEdges() { return (int)evec.size(); };
	bspEdge* getEdge(int id);
	int getEdgeId(int id) { return evec[id]; };
	Vector2d getBarycenter(){ return computeBarycenter(); }
	bool split(const bspLine &line);
	int addEdge(int id) {
		evec.push_back(id);
		return (int)evec.size()-1;
	}
	Vector2d computeBarycenter();
	void print();

	//computes the intersection of an edge with the line
	void intersect(Vector2d &v1, Vector2d &v2, const bspLine &line, int &i); 
};
//*****************************************************************************


//*****************************************************************************
class bspTree
{
private:
	std::vector<Vector2d> vertices;
	std::vector<bspEdge*> edges;
	bspNode* root;

public: 
	bspTree(); 
	~bspTree();

	void build(bspNode * node);

	int selectOptimalSplit(bspNode* node);
	int addVertex(double x,double y) { 
		vertices.push_back(Vector2d(x,y)); 
		return (int)vertices.size()-1;
	}
	int addVertex(const Vector2d &v){ 
		vertices.push_back(v); 
		return (int)vertices.size()-1;
	}
	Vector2d* getVertex(int id) { return &vertices[id]; };
	
	int addEmptyEdge() { 
		edges.push_back(new bspEdge(this));
		return (int)edges.size()-1;
	} 
	int addEdge(bspEdge* e) {
		edges.push_back(e);
		return (int)edges.size()-1;
	}
	bspEdge* getEdge(int id) { return (edges[id]); }
	bspNode* getRoot() { return root; }
	
	void initEnvironment();

	void traverse();
	void query(const normalDist& nd, std::multimap<double, int>& nbrs);
	void testQuery();

public:
	int cal_environment;
	int cal_obstacles;
	int cal_ellipse;
	
	void initVisualization();
};
//*****************************************************************************