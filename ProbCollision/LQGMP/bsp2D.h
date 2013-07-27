#ifndef __BSP2D_H__
#define __BSP2D_H__

#define dmax DBL_MAX
#define dmin DBL_MIN
#define deps 1e-08

#define SDIM 2

enum {FRONT, BACK, INTERSECTING, COINCIDENT};

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
struct bspLine2d
{
public:
	Vector2d n; //normal to split line
	double d; //nx*x+ny*y = d line equation

public: 
	double signDistance(const Vector2d &p) const;
	Vector2d intersectLine(const Vector2d &v1,const Vector2d &v2) const;
};
//*****************************************************************************


class bspTree2d;
class bspNode2d;

//*****************************************************************************
class bspEdge2d
{
private:
	static bspTree2d* btree;
	int e0, e1;
	Vector2d dv;
	bspLine2d line;
	Vector2d mid;
	double len;

public: 
	~bspEdge2d() {};
	bspEdge2d(bspTree2d * tree) { btree = tree; };
	bspEdge2d(bspTree2d * tree, const int v0, const int v1);
	void addEdge(const int v0, const int v1);
	Vector2d* getVertex(int id) const;
	Vector2d* getDirVec() { return &dv; };
	Vector2d* getMidPt() { return &mid; };
	double getLength() { return len; };
	bspLine2d* getLine() { return &line; };
	int getVertexId(int id) const;
	int checkAgainstLine(const bspLine2d &line);
};
//*****************************************************************************

//*****************************************************************************
class bspNode2d
{
private:
	static bspTree2d* btree;
	bspNode2d* front;
	bspNode2d* back;
	bspLine2d line;
	int level;
	std::vector<int> evec;
	int sedge;

public: 
	bspNode2d(bspTree2d * tree):front(0),back(0),level(0){
		btree = tree;
		evec.clear();
	};
	bspNode2d(bspNode2d * parent):front(0),back(0),level(parent->level+1){
		btree = parent->btree;
		evec.clear();
	};

	~bspNode2d(){ 
		delete front; delete back;
		evec.clear();
	}

	bool isLeaf() { return (front==0 && back==0); };
	bspNode2d* getFront() { return front; };
	bspNode2d* getBack () { return back; };
	bspLine2d* getLine() { return &line; };
	void setSplitEdge(int se) { sedge = se; };
	int getSplitEdge() { return sedge; };
	int getLevel() { return level; };
	int getNumEdges() { return (int)evec.size(); };
	bspEdge2d* getEdge(int id);
	int getEdgeId(int id) { return evec[id]; };
	Vector2d getBarycenter(){ return computeBarycenter(); }
	bool split(const bspLine2d &line);
	int addEdge(int id) {
		evec.push_back(id);
		return (int)evec.size()-1;
	}
	Vector2d computeBarycenter();
	void print();

	//computes the intersection of an edge with the line
	void intersect(Vector2d &v1, Vector2d &v2, const bspLine2d &line, int &i); 
};
//*****************************************************************************


//*****************************************************************************
class bspTree2d
{
private:
	std::vector<Vector2d> vertices;
	std::vector<bspEdge2d*> edges;
	bspNode2d* root;

public: 
	bspTree2d(); 
	~bspTree2d();

	void build(bspNode2d * node);

	int selectOptimalSplit(bspNode2d* node);
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
		edges.push_back(new bspEdge2d(this));
		return (int)edges.size()-1;
	} 
	int addEdge(bspEdge2d* e) {
		edges.push_back(e);
		return (int)edges.size()-1;
	}
	int getNumEdges() { return (int)edges.size(); };
	bspEdge2d* getEdge(int id) { return (edges[id]); }
	bspNode2d* getRoot() { return root; }

	void traverse();
	void query(const Matrix<SDIM>& x, const Matrix<SDIM, SDIM>& S, std::vector<std::pair<Matrix<SDIM>, double> >& cvx);
	void testQuery();
};
//*****************************************************************************

#endif