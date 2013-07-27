#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <deque>

#include <math.h>
#include <time.h>
#include <float.h>

#include "matrix.h"

#include "bsp2D.h"

//*****************************************************************************
// Vector2D
bool operator< (const Vector2d& v1, const Vector2d& v2 ) 
{
	if(fabs(v1.x - v2.x) < deps) {
		if(fabs(v1.y - v2.y) < deps) { return false ; }
		else { return (v1.y < v2.y); }
	}
	else { return (v1.x < v2.x); }
}
// sum
const Vector2d operator + (const Vector2d& v1, const Vector2d& v2) {
	return Vector2d(v1.x+v2.x,v1.y+v2.y);
}
// difference
const Vector2d operator - (const Vector2d& v1, const Vector2d& v2) {
	return Vector2d(v1.x-v2.x,v1.y-v2.y);
}
// unary minus
const Vector2d operator - (const Vector2d& v) {
	return Vector2d(-v.x,-v.y);
}
// dot product
const double operator * (const Vector2d& v1, const Vector2d& v2) {
	return (v1.x*v2.x + v1.y*v2.y);
}
// scalar product
const Vector2d operator * (const double& d, const Vector2d& v) {
	return Vector2d(v.x*d,v.y*d);
}
// scalar product
const Vector2d operator * (const Vector2d& v, const double& d) {
	return Vector2d(v.x*d,v.y*d);
}
// scalar division
const Vector2d operator / (const Vector2d& v, const double& d) {
	return Vector2d(v.x/d,v.y/d);
}
// cross product
double det(const Vector2d& v1, const Vector2d& v2) {
	return (v1.x * v2.y - v2.x * v1.y); 
}
// norm
double norm(const Vector2d &v) {
	return sqrt(v.x*v.x+v.y*v.y);
}
// squared norm
double sqnorm(const Vector2d& v) {
	return (v.x*v.x+v.y*v.y);
}
// normal vector
Vector2d normal(const Vector2d& v1, const Vector2d& v2) { 
	Vector2d nv = Vector2d(-(v2.y-v1.y), (v2.x-v1.x));
	nv = nv/norm(nv);
	return nv;
}
// direction vector
Vector2d dir(const Vector2d & v1, const Vector2d& v2) {
	Vector2d dv = Vector2d(v2.x-v1.x, v2.y-v1.y);
	dv = dv/norm(dv);
	return dv;
}
// bisector
Vector2d bisector(const Vector2d &p1, const Vector2d &p2)
{
	Vector2d bis = p1+p2;
	bis = bis/norm(bis);
	return bis;
}
// squared distance of point from line segment
double distSqPointSegment(const Vector2d& p1, const Vector2d& p2, const Vector2d& x)
{
	double r = ((x - p1)*(p2 - p1))/sqnorm(p2 - p1);
	if (r < 0) { return sqnorm(x - p1); } 
	else if (r > 1) { return sqnorm(x - p2); } 
	else { return sqnorm(x - (p1 + r*(p2-p1))); }
}

// signed distance from a line segment to a point
double signDistance(const Vector2d& p1, const Vector2d& p2, const Vector2d& x) {
	return det(p1 - x, p2 - p1);
}

//*****************************************************************************


//*****************************************************************************
// bspLine2d
// signed distance from line
double bspLine2d::signDistance( const Vector2d &p ) const {
	return (p.x*n.x+p.y*n.y-d);
}

// signed intersection of line with an edge
Vector2d bspLine2d::intersectLine(const Vector2d &v1,const Vector2d &v2) const {
	// convex coordinate of the intersection
	double t = (d - n.x*v1.x - n.y*v1.y) / (n.x*(v2.x-v1.x) + n.y*(v2.y-v1.y));
	return Vector2d(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y));
}
//*****************************************************************************


//*****************************************************************************
// bspEdge2d
bspTree2d* bspEdge2d::btree;

bspEdge2d::bspEdge2d(bspTree2d * tree, const int v0, const int v1) { 
	btree = tree;
	e0 = v0; e1 = v1;
	Vector2d* p0 = btree->getVertex(e0);
	Vector2d* p1 = btree->getVertex(e1);

	dv = dir(*p0, *p1);
	mid = (*p0 + *p1)*0.5;
	len = norm(*p1 - *p0);
	line.n = normal(*p0, *p1);
	line.d = (*p1*line.n);
};

Vector2d* bspEdge2d::getVertex(int id) const {
	if (id == 0) {
		return btree->getVertex(e0);
	} else if (id == 1) {
		return btree->getVertex(e1);
	} else {
		std::cerr << "Invalid edge vertex id" << std::endl;
		std::exit(-1);
	}
}

int bspEdge2d::checkAgainstLine(const bspLine2d &line)
{
	Vector2d* p0 = getVertex(0);
	double d0 = line.signDistance(*p0);

	Vector2d* p1 = getVertex(1);
	double d1 = line.signDistance(*p1);

	double m = std::min(d0, d1);
	double M = std::max(d0, d1);

	if (m >= -deps && M <= deps) return COINCIDENT;
	else if (m > deps && M > deps) return FRONT;
	else if (m < -deps && M < -deps) return BACK;
	else return INTERSECTING;
}

int bspEdge2d::getVertexId(int id) const {
	if (id == 0) { return e0; } else { return e1; }
}

void bspEdge2d::addEdge(const int v0, const int v1) { 
	e0 = v0; e1 = v1; 
	Vector2d* p0 = btree->getVertex(e0);
	Vector2d* p1 = btree->getVertex(e1);

	dv = dir(*p0, *p1);
	mid = (*p0 + *p1)*0.5;
	len = norm(*p1 - *p0);
	line.n = normal(*p0, *p1);
	line.d = (*p1*line.n);
}
//*****************************************************************************


//*****************************************************************************
// bspNode2d
bspTree2d* bspNode2d::btree;

bspEdge2d* bspNode2d::getEdge(int id) { 
	return btree->getEdge(evec[id]); 
}

Vector2d bspNode2d::computeBarycenter() {
	Vector2d sum(0,0);
	double count = 0;
	int ne = (int)evec.size();
	for(int it = 0; it < ne; ++it)
	{
		bspEdge2d * actual = btree->getEdge(evec[it]);
		sum = sum + *(actual->getVertex(0)) + *(actual->getVertex(1));
		count += 2;
	}
	return sum/count;
}

void bspNode2d::intersect(Vector2d &v1, Vector2d &v2, const bspLine2d &line, int &i) {  
	Vector2d v = line.intersectLine(v1,v2);
	i = btree->addVertex(v);
}

void bspNode2d::print() {
	std::cout << "level: " << level << " line: " << line.n.x << " " << line.n.y << " " << line.d;
	bspEdge2d * e = btree->getEdge(sedge);
	Vector2d* p0 = btree->getVertex(e->getVertexId(0));
	Vector2d* p1 = btree->getVertex(e->getVertexId(1));
	std::cout << " edge: " << sedge << " [ " << p0->x << ", " << p0->y << " ] -> [ " << p1->x << ", " << p1->y << " ] " << std::endl;
}

bool bspNode2d::split(const bspLine2d& sline)
{	   
	if (!isLeaf()) {
		std::cout << "Cannot split a non leaf node! \n";
		return false;
	}

	// Create the two child nodes
	front = new bspNode2d(this);
	back  = new bspNode2d(this);

	std::vector<int> intEdges;
	// separate edges according to front, back, intersection
	int ne = (int)evec.size();
	for(int i = 0; i < ne; ++i)
	{
		if (evec[i] != sedge) {
			bspEdge2d * ei = btree->getEdge(evec[i]);
			int val = ei->checkAgainstLine(sline);    
			switch(val)
			{
			case FRONT:
				front->addEdge(evec[i]);
				break;
			case BACK:
				back->addEdge(evec[i]);
				break;
			case INTERSECTING:
				intEdges.push_back(evec[i]);
				break;
			case COINCIDENT:
				if (*ei->getDirVec() * sline.n < deps) {
					front->addEdge(evec[i]);
				} else {
					back->addEdge(evec[i]);
				}
				break;
			}    
		}
	}

	// Split each intersecting edge and add them to front and back
	for(std::vector<int>::const_iterator it = intEdges.begin(); it != intEdges.end(); ++it)
	{ 
		bspEdge2d *e  = btree->getEdge(*it);

		// Check the sides of the edge intersection
		double d0 = sline.signDistance(*e->getVertex(0));
		double d1 = sline.signDistance(*e->getVertex(1));

		if ( (fabs(d0) < deps && d1 <= -deps) || (fabs(d1) < deps && d0 <= -deps) ) {
			back->addEdge(*it);
		}
		else if ( (fabs(d0) < deps && d1 >= deps) || (fabs(d1) < deps && d0 >= deps) ) {
			front->addEdge(*it);
		}
		else // intersection lies on the edge (not passing through vertices)
		{	
			int intId = -1;
			intersect(*e->getVertex(0),*e->getVertex(1), sline, intId);

			bspEdge2d *ef = new bspEdge2d(btree); 
			bspEdge2d *eb = new bspEdge2d(btree);

			if (d0 >= deps && d1 <= -deps) {
				ef->addEdge(e->getVertexId(0), intId);
				eb->addEdge(intId, e->getVertexId(1));
			} else if (d0 <= -deps && d1 >= deps) {
				eb->addEdge(e->getVertexId(0), intId);
				ef->addEdge(intId, e->getVertexId(1));
			}

			front->addEdge(btree->addEdge(ef));
			back->addEdge(btree->addEdge(eb));
		}	
	}

	line = sline;
	front->computeBarycenter();
	back->computeBarycenter();
	return true;
}
//*****************************************************************************


//*****************************************************************************
// bspTree2d
int bspTree2d::selectOptimalSplit(bspNode2d * node)
{
	int optimalSplit = 0;
	int ne = node->getNumEdges();

	int minFront = ne;
	int minBack = ne;
	bool flag = false;

	bspLine2d* sline;

	for (int i = 0; i < ne; ++i) 
	{
		int frontSize = 0;
		int backSize = 0;

		bspEdge2d * ei = node->getEdge(i);
		sline = ei->getLine();

		for(int j = 0; j < ne; ++j) 
		{
			if (i == j) {
				continue;
			}

			bspEdge2d * ej = node->getEdge(j);
			int val = ej->checkAgainstLine(*sline);    

			if (val == FRONT) { ++frontSize; }
			else if (val == BACK) { ++backSize; }
			else if (val == INTERSECTING) { ++frontSize; ++backSize; }
			else if (val == COINCIDENT) { (*ej->getDirVec() * sline->n < deps)? ++frontSize : ++backSize; }

			flag = (std::make_pair(std::max(frontSize, backSize), std::min(frontSize, backSize)) >= std::make_pair(std::max(minFront, minBack), std::min(minFront, minBack)));
			if (flag) break;
		}

		if (!flag) {
			minFront = frontSize;
			minBack = backSize;
			optimalSplit = i;
		}
	}
	return optimalSplit;
}

void bspTree2d::build(bspNode2d * node)
{
	std::vector<bspNode2d*> Q;
	Q.push_back(node);  

	while (Q.size() > 0)
	{
		bspNode2d * r = Q.back();
		Q.pop_back();

		int id = selectOptimalSplit(r);

		bspEdge2d * se = r->getEdge(id);
		bspLine2d* sline = se->getLine();

		//std::cout << "Split line from edge: " << r->getEdgeId(id) << " " << sline.n.x << " " << sline.n.y << " " << sline.d << std::endl;

		r->setSplitEdge(r->getEdgeId(id));
		r->split(*sline);

		if (r->getFront()->getNumEdges() > 0) 
			Q.push_back(r->getFront());
		if (r->getBack()->getNumEdges() > 0) 
			Q.push_back(r->getBack());
	}
}

bspTree2d::bspTree2d() {
	root = new bspNode2d(this);
}

bspTree2d::~bspTree2d() 
{
	delete root;
	for( std::vector<bspEdge2d*>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
		delete *it;
	}
}

void bspTree2d::traverse() 
{ 
	std::deque<bspNode2d*> Q;
	Q.push_back(root);

	while (!Q.empty())
	{
		bspNode2d* node = Q.front();
		Q.pop_front();

		// print node details
		node->print();

		if(!node->isLeaf()) {
			if(node->getFront()->getNumEdges() > 0)
				Q.push_back(node->getFront());
			if (node->getBack()->getNumEdges() > 0)
				Q.push_back(node->getBack());
		}  
	}
}

void bspTree2d::query(const Matrix<SDIM>& x, const Matrix<SDIM, SDIM>& S, std::vector<std::pair<Matrix<SDIM>, double> >& cvx)
{
	std::deque<bspNode2d*> Q;
	std::multimap<double, int> nbrs;
	Q.push_back(root);

	Matrix<SDIM, SDIM> V, E;
	jacobi(S, V, E);

	double sx = 3.0*sqrt(E(0,0));
	double sy = 3.0*sqrt(E(1,1));
	double phi = std::atan2(V(1,0), V(0,0));

	double sxinv = 1.0/sx, syinv = 1.0/sy;
	double sphi = sin(phi), cphi = cos(phi);
	Vector2d p(x[0], x[1]);
	Vector2d pt(sxinv*(cphi*x[0]+sphi*x[1]), syinv*(-sphi*x[0]+cphi*x[1]));

	double radiusSq = 1.0;

	while (!Q.empty()) 
	{
		bspNode2d* node = Q.front();
		Q.pop_front();

		if (node->isLeaf()) {
			continue;
		}

		bspEdge2d * e = getEdge(node->getSplitEdge());
		Vector2d * p0 = e->getVertex(0);
		Vector2d * p1 = e->getVertex(1);

		// Inverse transform edge
		Vector2d p0t(sxinv*(cphi*p0->x+sphi*p0->y), syinv*(-sphi*p0->x+cphi*p0->y));
		Vector2d p1t(sxinv*(cphi*p1->x+sphi*p1->y), syinv*(-sphi*p1->x+cphi*p1->y));

		double sqdist = distSqPointSegment(p0t, p1t, pt);

		bspLine2d sline;
		sline.n = normal(p0t, p1t);
		sline.d = (p1t*sline.n);

		double sgndist = sline.signDistance(pt);

		//std::cout << std::setprecision(12) << "edge: " << node->getSplitEdge() << " sqdist: " << sqdist << " sgndist: " << sgndist << std::endl;

		if (sqdist <= radiusSq+deps) {
			nbrs.insert(std::make_pair<double, int>(sqdist, node->getSplitEdge()));
		}

		if (sgndist*sgndist > radiusSq) {
			(sgndist > 0)? Q.push_back(node->getFront()) : Q.push_back(node->getBack());
		} else {
			Q.push_back(node->getFront());
			Q.push_back(node->getBack());
		}
	}

	// filter edges in the order of increasing distance
	Matrix<SDIM> a;
	double b;

	std::multimap<double, int>::iterator it1, it2;
	it1 = nbrs.begin();
	while(it1 != nbrs.end()) 
	{
		bspEdge2d* e1 = getEdge(it1->second);
		bspLine2d* l1 = e1->getLine();
		double sgn1 = l1->signDistance(p);

		it2 = it1;
		++it2;
		while(it2 != nbrs.end()) 
		{
			//std::cout << "Processing: " << it1->second << " and " << it2->second << std::endl;
			bspEdge2d* e2 = getEdge(it2->second);
			double sgn2 = l1->signDistance(*e2->getMidPt());

			if ( fabs(it2->first - it1->first) < deps ) 
			{
				// hack!
				Vector2d* cp;
				if (it1->second < it2->second) {
					cp = getEdge(it1->second)->getVertex(1);
				} else {
					cp = getEdge(it1->second)->getVertex(0);
				}
				Vector2d nv = (p - *cp);
				nv = nv/norm(nv);

				a[0] = -nv.x; a[1] = -nv.y;
				b = a[0]*cp->x + a[1]*cp->y;

				cvx.push_back(std::make_pair(a, b));
				//std::cout << "Pushing perpendicular constraint" << std::endl;

				it2 = nbrs.erase(it2);
				it1 = nbrs.erase(it1);
				it1 = it2;
			} else if (sgn1*sgn2 <= 0) {
				it2 = nbrs.erase(it2);
			}
			if (it2 != nbrs.end()) {
				++it2;
			}
		}
		if (it1 != nbrs.end()) {
			++it1;
		}
	}

	for(std::multimap<double, int>::iterator it = nbrs.begin(); it != nbrs.end(); ++it) {
		bspLine2d * l = getEdge(it->second)->getLine();
		a[0] = -l->n.x; a[1] = -l->n.y;
		b = -l->d;
		cvx.push_back(std::make_pair(a, b));
	}

	// filter edges in the order of increasing distance
	/*
	std::multimap<double, int>::iterator it1, it2;
	it1 = nbrs.begin();
	while(it1 != nbrs.end()) 
	{
		bspEdge2d* e1 = getEdge(it1->second);
		bspLine2d* l1 = e1->getLine();
		double sgn1 = l1->signDistance(p);

		it2 = it1;
		++it2;
		while(it2 != nbrs.end()) {
			//std::cout << "Processing: " << it1->second << " and " << it2->second << std::endl;
			bspEdge2d* e2 = getEdge(it2->second);
			double sgn2 = l1->signDistance(*e2->getMidPt());

			if ( fabs(it2->first - it1->first) < deps ) {
				double sgn3 = e2->getLine()->signDistance(p);
				if ( (sgn1*sgn3 < 0) || (sgn1*sgn3 > 0 && fabs(sgn3) > fabs(sgn1))) {
					it1 = nbrs.erase(it1);
					it2 = it1;
					++it2;
				} else if (sgn1*sgn3 > 0) {
					it2 = nbrs.erase(it2);
				}
				continue;
			}
			
			if (sgn1*sgn2 <= 0) {
				it2 = nbrs.erase(it2);	
			} else {
				++it2;
			}
		}
		++it1;
	}

	for(std::multimap<double, int>::iterator it = nbrs.begin(); it != nbrs.end(); ++it) {
		bspLine2d * l = getEdge(it->second)->getLine();
		Matrix<SDIM> a;
		a[0] = -l->n.x; a[1] = -l->n.y;
		double b = -l->d;
		cvx.push_back(std::make_pair(a, b));
	}
	*/
}