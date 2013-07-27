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
	Matrix<3> p;
	Matrix<3,3> M;
};

//*****************************************************************************
struct Vector3d
{
public:  
	double x, y, z;

public: 
	Vector3d(): x(0),y(0),z(0){}
	Vector3d(double px, double py, double pz): x(px), y(py), z(pz){};
	~Vector3d(){}
};

bool operator< (const Vector3d& v1, const Vector3d& v2); 
const Vector3d operator + (const Vector3d& v1, const Vector3d& v2);
const Vector3d operator - (const Vector3d& v1, const Vector3d& v2);
const Vector3d operator - (const Vector3d& v);
const double operator * (const Vector3d& v1, const Vector3d& v2);
const Vector3d operator * (const double& d, const Vector3d& v);
const Vector3d operator * (const Vector3d& v, const double& d);
const Vector3d operator / (const Vector3d& v, const double& d);
double norm(const Vector3d &v);
double sqnorm(const Vector3d& v);
Vector3d normal(const Vector3d& v1, const Vector3d& v2);
Vector3d dir(const Vector3d & v1, const Vector3d& v2);
Vector3d bisector(const Vector3d &v1, const Vector3d &v2);
double signDistance(const Vector3d& p1, const Vector3d& p2, const Vector3d& x);
//*****************************************************************************


//*****************************************************************************
struct bspPlane
{
public:
	Vector3d n; //normal to split line
	double d; //nx*x+ny*y+nz*z = d line equation

public: 
	//signed distance from line
	double signDistance(const Vector3d &p) const;
	//signed intersection of plane with an edge
	Vector3d intersectLine(const Vector3d &v1,const Vector3d &v2) const;
};
//*****************************************************************************


class bspTree;
class bspNode;

//*****************************************************************************
class bspPolygon
{
private:
	static bspTree* btree;
	std::list<int> mlist;
	
//constru'-.,_,.-'"^`^"'-.,_,.-'"^`^"'-.,_,.-'"^`^"'-.,_,.-'"^`^"
public: 
  //bspPolygon(bspTree * tree):m_tree(tree),m_n_points(0),m_points(0),m_list(0){};
bspPolygon(bspTree * tree):m_list(0){m_tree=tree;};

//destruc'-.,_,.-'"^`^"'-.,_,.-'"^`^"'-.,_,.-'"^`^"'-.,_,.-'"^`^"'
public: 
 ~bspPolygon()
  {
    //delete[] m_points;
    //return m_points[id];
    bsplistnode * node = m_list;
    bsplistnode * lastnode = 0;
    //go to last element
    for(int i=0; true ;++i)
    {
      delete lastnode;
      if (node == 0) break;
      lastnode = node;
      node = node->next;
    }     
  }

//accesso'-.,_,.-'"^`^"'-.,_,.-'"^`^"'-.,_,.-'"^`^"'-.,_,.-'"^`^"' 
public: 
  //append a vertex
  void        add_vertex(const int id) ;
  //get the vertex id
  int         get_vertex_id(int id) const;
  //next vertex 
  int         get_next(int id) const {if ((++id)>=size()) id=0; return id;}
  //previus vertex 
  int         get_prev(int id) const {if ((--id)<0) id=size()-1; return id;}
  //return the vertex geometry
  bspPoint3D* get_vertex_geometry (int id) const;
  //number of points
  int         size() const;;
  //check if a given plane intersect the polygon
  int         check_against_plane( const bspPlane &plane);   


public: 
	~bspPolygon() {};
	bspPolygon(bspTree * tree) { btree = tree; };
	bspPolygon(bspTree * tree, const int v0, const int v1);
	void addEdge(const int v0, const int v1);
	Vector3d* getVertex(int id) const;
	Vector3d* getDirVec() { return &dv; };
	Vector3d* getMidPt() { return &mid; };
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
	bspPolygon* getEdge(int id);
	int getEdgeId(int id) { return evec[id]; };
	Vector3d getBarycenter(){ return computeBarycenter(); }
	bool split(const bspLine &line);
	int addEdge(int id) {
		evec.push_back(id);
		return (int)evec.size()-1;
	}
	Vector3d computeBarycenter();
	void print();

	//computes the intersection of an edge with the line
	void intersect(Vector3d &v1, Vector3d &v2, const bspLine &line, int &i); 
};
//*****************************************************************************


//*****************************************************************************
class bspTree
{
private:
	std::vector<Vector3d> vertices;
	std::vector<bspPolygon*> edges;
	bspNode* root;

public: 
	bspTree(); 
	~bspTree();

	void build(bspNode * node);

	int selectOptimalSplit(bspNode* node);
	int addVertex(double x,double y) { 
		vertices.push_back(Vector3d(x,y)); 
		return (int)vertices.size()-1;
	}
	int addVertex(const Vector3d &v){ 
		vertices.push_back(v); 
		return (int)vertices.size()-1;
	}
	Vector3d* getVertex(int id) { return &vertices[id]; };
	
	int addEmptyEdge() { 
		edges.push_back(new bspPolygon(this));
		return (int)edges.size()-1;
	} 
	int addEdge(bspPolygon* e) {
		edges.push_back(e);
		return (int)edges.size()-1;
	}
	bspPolygon* getEdge(int id) { return (edges[id]); }
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