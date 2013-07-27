#include "bsp2D.h"

//*****************************************************************************
// Vector2D
bool operator< (const Vector3d& v1, const Vector3d& v2 ) 
{
	if(fabs(v1.x - v2.x) < deps) {
		if(fabs(v1.y - v2.y) < deps) { return false ; }
		else { return (v1.y < v2.y); }
	}
	else { return (v1.x < v2.x); }
}
// sum
const Vector3d operator + (const Vector3d& v1, const Vector3d& v2) {
	return Vector3d(v1.x+v2.x,v1.y+v2.y);
}
// difference
const Vector3d operator - (const Vector3d& v1, const Vector3d& v2) {
	return Vector3d(v1.x-v2.x,v1.y-v2.y);
}
// unary minus
const Vector3d operator - (const Vector3d& v) {
	return Vector3d(-v.x,-v.y);
}
// dot product
const double operator * (const Vector3d& v1, const Vector3d& v2) {
	return (v1.x*v2.x + v1.y*v2.y);
}
// scalar product
const Vector3d operator * (const double& d, const Vector3d& v) {
	return Vector3d(v.x*d,v.y*d);
}
// scalar product
const Vector3d operator * (const Vector3d& v, const double& d) {
	return Vector3d(v.x*d,v.y*d);
}
// scalar division
const Vector3d operator / (const Vector3d& v, const double& d) {
	return Vector3d(v.x/d,v.y/d);
}
// cross product
double det(const Vector3d& v1, const Vector3d& v2) {
	return (v1.x * v2.y - v2.x * v1.y); 
}
// norm
double norm(const Vector3d &v) {
	return sqrt(v.x*v.x+v.y*v.y);
}
// squared norm
double sqnorm(const Vector3d& v) {
	return (v.x*v.x+v.y*v.y);
}
// normal vector
Vector3d normal(const Vector3d& v1, const Vector3d& v2) { 
	Vector3d nv = Vector3d(-(v2.y-v1.y), (v2.x-v1.x));
	nv = nv/norm(nv);
	return nv;
}
// direction vector
Vector3d dir(const Vector3d & v1, const Vector3d& v2) {
	Vector3d dv = Vector3d(v2.x-v1.x, v2.y-v1.y);
	dv = dv/norm(dv);
	return dv;
}
// bisector
Vector3d bisector(const Vector3d &p1, const Vector3d &p2)
{
	Vector3d bis = p1+p2;
	bis = bis/norm(bis);
	return bis;
}
// squared distance of point from line segment
double distSqPointSegment(const Vector3d& p1, const Vector3d& p2, const Vector3d& x)
{
	double r = ((x - p1)*(p2 - p1))/sqnorm(p2 - p1);
	if (r < 0) { return sqnorm(x - p1); } 
	else if (r > 1) { return sqnorm(x - p2); } 
	else { return sqnorm(x - (p1 + r*(p2-p1))); }
}

// signed distance from a line segment to a point
double signDistance(const Vector3d& p1, const Vector3d& p2, const Vector3d& x) {
    return det(p1 - x, p2 - p1);
}

//*****************************************************************************


//*****************************************************************************
// bspLine
// signed distance from line
double bspLine::signDistance( const Vector3d &p ) const {
	return (p.x*n.x+p.y*n.y-d);
}

// signed intersection of line with an edge
Vector3d bspLine::intersectLine(const Vector3d &v1,const Vector3d &v2) const {
	// convex coordinate of the intersection
	double t = (d - n.x*v1.x - n.y*v1.y) / (n.x*(v2.x-v1.x) + n.y*(v2.y-v1.y));
	return Vector3d(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y));
}
//*****************************************************************************


//*****************************************************************************
// bspPolygon
bspTree* bspPolygon::btree;

bspPolygon::bspPolygon(bspTree * tree, const int v0, const int v1) { 
	btree = tree;
	e0 = v0; e1 = v1;
	Vector3d* p0 = btree->getVertex(e0);
	Vector3d* p1 = btree->getVertex(e1);

	dv = dir(*p0, *p1);
	mid = (*p0 + *p1)*0.5;
	len = norm(*p1 - *p0);
	line.n = normal(*p0, *p1);
	line.d = (*p1*line.n);
};

Vector3d* bspPolygon::getVertex(int id) const {
	if (id == 0) {
		return btree->getVertex(e0);
	} else if (id == 1) {
		return btree->getVertex(e1);
	} else {
		std::cerr << "Invalid edge vertex id" << std::endl;
		std::exit(-1);
	}
}

int bspPolygon::checkAgainstLine(const bspLine &line)
{
	Vector3d* p0 = getVertex(0);
	double d0 = line.signDistance(*p0);
	
	Vector3d* p1 = getVertex(1);
	double d1 = line.signDistance(*p1);

	double m = std::min(d0, d1);
	double M = std::max(d0, d1);
	
	if (m >= -deps && M <= deps) return COINCIDENT;
	else if (m > deps && M > deps) return FRONT;
	else if (m < -deps && M < -deps) return BACK;
	else return INTERSECTING;
}

int bspPolygon::getVertexId(int id) const {
	if (id == 0) { return e0; } else { return e1; }
}

void bspPolygon::addEdge(const int v0, const int v1) { 
	e0 = v0; e1 = v1; 
	Vector3d* p0 = btree->getVertex(e0);
	Vector3d* p1 = btree->getVertex(e1);

	dv = dir(*p0, *p1);
	mid = (*p0 + *p1)*0.5;
	len = norm(*p1 - *p0);
	line.n = normal(*p0, *p1);
	line.d = (*p1*line.n);
}
//*****************************************************************************


//*****************************************************************************
// bspNode
bspTree* bspNode::btree;

bspPolygon* bspNode::getEdge(int id) { 
	return btree->getEdge(evec[id]); 
}

Vector3d bspNode::computeBarycenter() {
	Vector3d sum(0,0);
	double count = 0;
	int ne = (int)evec.size();
	for(int it = 0; it < ne; ++it)
	{
		bspPolygon * actual = btree->getEdge(evec[it]);
		sum = sum + *(actual->getVertex(0)) + *(actual->getVertex(1));
		count += 2;
	}
	return sum/count;
}

void bspNode::intersect(Vector3d &v1, Vector3d &v2, const bspLine &line, int &i) {  
	Vector3d v = line.intersectLine(v1,v2);
	i = btree->addVertex(v);
}

void bspNode::print() {
	std::cout << "level: " << level << " line: " << line.n.x << " " << line.n.y << " " << line.d;
	bspPolygon * e = btree->getEdge(sedge);
	Vector3d* p0 = btree->getVertex(e->getVertexId(0));
	Vector3d* p1 = btree->getVertex(e->getVertexId(1));
	std::cout << " edge: " << sedge << " [ " << p0->x << ", " << p0->y << " ] -> [ " << p1->x << ", " << p1->y << " ] " << std::endl;
}

bool bspNode::split(const bspLine& sline)
{	   
	if (!isLeaf()) {
		std::cout << "Cannot split a non leaf node! \n";
		return false;
	}

	// Create the two child nodes
	front = new bspNode(this);
	back  = new bspNode(this);

	std::vector<int> intEdges;
	// separate edges according to front, back, intersection
	int ne = (int)evec.size();
	for(int i = 0; i < ne; ++i)
	{
		if (evec[i] != sedge) {
			bspPolygon * ei = btree->getEdge(evec[i]);
			int val = ei->checkAgainstLine(sline);    
			switch(val)
			{
			case FRONT:
				front->addEdge(evec[i]);
				//std::cout << "Adding edge: " << i << " to FRONT" << std::endl;
				break;
			case BACK:
				back->addEdge(evec[i]);
				//std::cout << "Adding edge: " << i << " to BACK" << std::endl;
				break;
			case INTERSECTING:
				intEdges.push_back(evec[i]);
				//std::cout << "Adding edge: " << i << " to INTERSECTING" << std::endl;
				break;
			case COINCIDENT:
				if (*ei->getDirVec() * sline.n < deps) {
					front->addEdge(evec[i]);
					//std::cout << "Adding edge: " << i << " to COINCIDENT -> FRONT" << std::endl;
				} else {
					back->addEdge(evec[i]);
					//std::cout << "Adding edge: " << i << " to COINCIDENT -> BACK" << std::endl;
				}
				break;
			}    
		}
	}

	// Split each intersecting edge and add them to front and back
	for(std::vector<int>::const_iterator it = intEdges.begin(); it != intEdges.end(); ++it)
	{ 
		bspPolygon *e  = btree->getEdge(*it);
		
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

			bspPolygon *ef = new bspPolygon(btree); 
			bspPolygon *eb = new bspPolygon(btree);
			
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
// bspTree
int bspTree::selectOptimalSplit(bspNode * node)
{
	int optimalSplit = 0;
	int ne = node->getNumEdges();
	
	int minFront = ne;
	int minBack = ne;
	bool flag = false;

	bspLine* sline;

	for (int i = 0; i < ne; ++i) 
	{
		int frontSize = 0;
		int backSize = 0;

		bspPolygon * ei = node->getEdge(i);
		sline = ei->getLine();

		for(int j = 0; j < ne; ++j) 
		{
			if (i == j) {
				continue;
			}
			
			bspPolygon * ej = node->getEdge(j);
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

void bspTree::build(bspNode * node)
{
	std::vector<bspNode*> Q;
	Q.push_back(node);  

	while (Q.size() > 0)
	{
		bspNode * r = Q.back();
		Q.pop_back();

		int id = selectOptimalSplit(r);

		bspPolygon * se = r->getEdge(id);
		bspLine* sline = se->getLine();

		//std::cout << "Split line from edge: " << r->getEdgeId(id) << " " << sline.n.x << " " << sline.n.y << " " << sline.d << std::endl;

		r->setSplitEdge(r->getEdgeId(id));
		r->split(*sline);

		/*
		int nfe = r->getFront()->getNumEdges();
		std::cout << "Front: ";
		for(int i = 0; i < nfe; ++i) {
			std::cout << r->getFront()->getEdgeId(i) << " ";
		}
		std::cout << std::endl;
		int nfb = r->getBack()->getNumEdges();
		std::cout << "Back: ";
		for(int i = 0; i < nfb; ++i) {
			std::cout << r->getBack()->getEdgeId(i) << " ";
		}
		std::cout << std::endl;
		*/

		if (r->getFront()->getNumEdges() > 0) 
			Q.push_back(r->getFront());
		if (r->getBack()->getNumEdges() > 0) 
			Q.push_back(r->getBack());
	}
}

bspTree::bspTree()
{
	root = new bspNode(this);
	
	initEnvironment();
	initVisualization();

	root->computeBarycenter();
}


bspTree::~bspTree() 
{
	delete root;
	for( std::vector<bspPolygon*>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
		delete *it;
	}
}

void bspTree::initEnvironment()
{
	// Environment 1
	int v[24];

	v[0] = addVertex(-9.5, -9.5);
	v[1] = addVertex(-5.5, -9.5);
	v[2] = addVertex(-5.5, -7.5);
	v[3] = addVertex(-0.5, -7.5);
	v[4] = addVertex(-0.5, -9.5);
	v[5] = addVertex(9.5, -9.5);
	v[6]  = addVertex(9.5, -5.5);
	v[7]  = addVertex(7.5, -5.5);
	v[8]  = addVertex(7.5, -0.5);
	v[9] = addVertex(9.5, -0.5);
	v[10] = addVertex(9.5, 9.5);
	v[11] = addVertex(-0.5, 9.5);
	v[12] = addVertex(-0.5, 7.5);
	v[13] = addVertex(-5.5, 7.5);
	v[14] = addVertex(-5.5, 9.5);
	v[15] = addVertex(-9.5, 9.5);
	v[16] = addVertex(-9.5, -0.5);
	v[17] = addVertex(-7.5, -0.5);
	v[18] = addVertex(-7.5, -5.5);
	v[19] = addVertex(-9.5, -5.5);

	v[20] = addVertex(-5.5, -5.5);
	v[21] = addVertex(-5.5, -0.5);
	v[22] = addVertex(-0.5, -0.5);
	v[23] = addVertex(-0.5, -5.5);

	int p = 0;
	bspPolygon * e;
	for(int i = 0; i < 19; ++i) {
		e = new bspPolygon(this, v[i], v[i+1]); p = addEdge(e); root->addEdge(p);
	}
	e = new bspPolygon(this, v[19], v[0]); p = addEdge(e); root->addEdge(p);
	for(int i = 20; i < 23; ++i) {
		e = new bspPolygon(this, v[i], v[i+1]); p = addEdge(e); root->addEdge(p); 
	}
	e = new bspPolygon(this, v[23], v[20]); p = addEdge(e); root->addEdge(p); 

	// Environment 2
	/*
	int v[16];

	v[0] = addVertex(-10.0, -10.0);
	v[1] = addVertex(10.0, -10.0);
	v[2] = addVertex(10.0, 10.0);
	v[3] = addVertex(-10.0, 10.0);

	v[4] = addVertex(-6.0, -2.0);
	v[5] = addVertex(-1.5, -4.5);
	v[6] = addVertex(-4.0, -5.0);
	v[7] = addVertex(-5.0, -8.5);
	v[8] = addVertex(-8.0, -3.5);
	v[9] = addVertex(-5.5, -5.0);

	v[10] = addVertex(6.0, 0.0);
	v[11] = addVertex(1.5, 2.5);
	v[12] = addVertex(4.0, 3.0);
	v[13] = addVertex(5.0, 6.5);
	v[14] = addVertex(8.0, 1.5);
	v[15] = addVertex(5.5, 3.0);

	int p = 0;
	bspPolygon * e;
	for(int i = 0; i < 3; ++i) {
		e = new bspPolygon(this, v[i], v[i+1]); p = addEdge(e); root->addEdge(p);
	}
	e = new bspPolygon(this, v[3], v[0]); p = addEdge(e); root->addEdge(p);
	for(int i = 4; i < 9; ++i) {
		e = new bspPolygon(this, v[i], v[i+1]); p = addEdge(e); root->addEdge(p); 
	}
	e = new bspPolygon(this, v[9], v[4]); p = addEdge(e); root->addEdge(p); 

	for(int i = 10; i < 15; ++i) {
		e = new bspPolygon(this, v[i], v[i+1]); p = addEdge(e); root->addEdge(p); 
	}
	e = new bspPolygon(this, v[15], v[10]); p = addEdge(e); root->addEdge(p); 
	*/
}

inline double mod2pi(double x) {
	double result = fmod(x, (2*M_PI));
	if (result < 0) { result += (2*M_PI); }
	return result;
}

void bspTree::initVisualization() 
{
	CAL_SetViewParams(0, 0.0f, 0.0f, 26.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f); 

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0.7f, 0.0f, 0.0f, 1.0f);

	int ne = (int)edges.size();
	for(int i = 0; i < ne; ++i) 
	{
		bspPolygon* e = edges[i];
		Vector3d* mid = e->getMidPt();
		Vector3d* dv = e->getDirVec();

		int obj;
		CAL_CreateBox(cal_obstacles, (float)(e->getLength()+0.2), 0.2f, 0.2f, (float)mid->x, (float)mid->y, 0.0f, &obj);
		CAL_SetObjectOrientation(obj, 0.0f, 0.0f, (float)(mod2pi(atan2(dv->y, dv->x))));
	}
}

void bspTree::traverse() 
{ 
	std::deque<bspNode*> Q;
	Q.push_back(root);
	
	while (!Q.empty())
	{
		bspNode* node = Q.front();
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

void bspTree::query(const normalDist& nd, std::multimap<double, int>& nbrs)
{
	std::deque<bspNode*> Q;
	Q.push_back(root);
	
	Matrix<2,2> V, E;
	jacobi(nd.M, V, E);

	double sx = 3.0*sqrt(E(0,0));
	double sy = 3.0*sqrt(E(1,1));
	double phi = std::atan2(V(1,0), V(0,0));

	double sxinv = 1.0/sx, syinv = 1.0/sy;
	double sphi = sin(phi), cphi = cos(phi);
	Vector3d p(nd.x, nd.y);
	Vector3d pt(sxinv*(cphi*nd.x+sphi*nd.y), syinv*(-sphi*nd.x+cphi*nd.y));

	double radiusSq = 1.0;

	while (!Q.empty()) 
	{
		bspNode* node = Q.front();
		Q.pop_front();
		
		if (node->isLeaf()) {
			continue;
		}

		bspPolygon * e = getEdge(node->getSplitEdge());
		Vector3d * p0 = e->getVertex(0);
		Vector3d * p1 = e->getVertex(1);

		// Inverse transform edge
		Vector3d p0t(sxinv*(cphi*p0->x+sphi*p0->y), syinv*(-sphi*p0->x+cphi*p0->y));
		Vector3d p1t(sxinv*(cphi*p1->x+sphi*p1->y), syinv*(-sphi*p1->x+cphi*p1->y));
		
		double sqdist = distSqPointSegment(p0t, p1t, pt);
		
		bspLine sline;
		sline.n = normal(p0t, p1t);
		sline.d = (p1t*sline.n);

		double sgndist = sline.signDistance(pt);
	
		std::cout << std::setprecision(12) << "edge: " << node->getSplitEdge() << " sqdist: " << sqdist << " sgndist: " << sgndist << std::endl;
		
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

	for(std::multimap<double, int>::iterator it = nbrs.begin(); it != nbrs.end(); ++it) {
		std::cout << "Neighbor: " << it->second << " dist: " << it->first << std::endl;
	}

	// filter edges in the order of increasing distance
	std::multimap<double, int>::iterator it1, it2;
	it1 = nbrs.begin();
	while(it1 != nbrs.end()) 
	{
		bspPolygon* e1 = getEdge(it1->second);
		bspLine* l1 = e1->getLine();
		double sgn1 = l1->signDistance(p);

		it2 = it1;
		++it2;
		while(it2 != nbrs.end()) {
			std::cout << "Processing: " << it1->second << " and " << it2->second << std::endl;
			bspPolygon* e2 = getEdge(it2->second);
			double sgn2 = l1->signDistance(*e2->getMidPt());

			if ( fabs(it2->first - it1->first) < deps ) {
				double sgn3 = e2->getLine()->signDistance(p);
				if (sgn1 < 0 && sgn3 > 0) {
					it1 = nbrs.erase(it1);
					it2 = it1;
					++it2;
				} else if (sgn1 > 0 && sgn3 > 0) {
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
	
	// debug
	for(std::multimap<double, int>::iterator it = nbrs.begin(); it != nbrs.end(); ++it) {
		//std::cout << "Updated Neighbor: " << it->second << " dist: " << it->first << std::endl;
		bspPolygon * e = getEdge(it->second);
		Vector3d* p0 = getVertex(e->getVertexId(0));
		Vector3d* p1 = getVertex(e->getVertexId(1));
		std::cout << "Edge: " << it->second << " [ " << p0->x << ", " << p0->y << " ] -> [ " << p1->x << ", " << p1->y << " ] " << " dist: " << it->first << std::endl;
	}
}

inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) 
{
	double x = R(2,1) - R(1,2);
	double y = R(0,2) - R(2,0);
	double z = R(1,0) - R(0,1);
	double r = sqrt(x*x+y*y+z*z);
	double t = R(0,0) + R(1,1) + R(2,2);
	double angle = atan2(r,t-1);
	if (angle != 0) {
		x /= r;
		y /= r;
		z /= r;
	} else {
		x = 0;
		y = 0;
		z = 0;
	}
	Matrix<4,1> q;
	q(0,0) = sin(angle/2)*x;
	q(1,0) = sin(angle/2)*y;
	q(2,0) = sin(angle/2)*z;
	q(3,0) = cos(angle/2);

	return q;
}

inline double uniform(double low, double high) { return low + ((double)rand()/(double)RAND_MAX)*(high-low); }

template <size_t _numRows, size_t _numColumns>
inline void rndMatrix(Matrix<_numRows, _numColumns>& q, double low, double high){
	for(size_t i = 0; i < _numRows; ++i) {
		for(size_t j = 0; j < _numColumns; ++j) {
			q(i,j) = uniform(low,high);
		}
	}
}

void bspTree::testQuery() 
{
	std::multimap<double, int> nbrs;

	CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
	CAL_SetGroupScaling(cal_environment, 1, 1, 1);

	normalDist nd;


	nd.x = -5.95528; nd.y = -5.51622;
	nd.M(0,0) = 0.0265071; nd.M(0,1) = 0.000708107;
	nd.M(1,0) = 0.000708107; nd.M(1,1) = 0.00281007;
	
	//srand(time(NULL));
	//for(int i = 0; i < 1000; ++i) {
	//	rand();
	//}

	//nd.x = uniform(-10.0, 10.0); nd.y = uniform(-10.0, 10.0);
	//Matrix<2,2> S;
	//rndMatrix(S, 0.1, 1.0);
	//nd.M = ~S*S;

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.6f);

	Matrix<2,2> V, E;
	jacobi(nd.M, V, E);
	Matrix<3,3> T = identity<3>();
	T.insert(0,0, V);
	Matrix<3,3> R = identity<3>();
	R(1,1) = R(2,2) = 0.0;
	R(1,2) = -1.0; R(2,1) = 1.0;
	Matrix<4,1> q = quatFromRot(T*R);

	int obj;
	CAL_CreateCylinder(cal_ellipse, 1.0f, 0.1f, 0.0f, 0.0f, 0.0f, &obj);
	CAL_SetObjectPosition(obj, (float)nd.x, (float)nd.y, 0.1f);
	CAL_SetObjectQuaternion(obj, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetObjectScaling(obj, (float) (3.0*sqrt(E(0,0))), 1.0f, (float) (3.0*sqrt(E(1,1))));

	CAL_CreateSphere(cal_ellipse, 0.15f, (float)nd.x, (float)nd.y, 0.1f, &obj);
	CAL_SetObjectColor(obj, 0.0f, 0.0f, 0.8f);
	
	query(nd, nbrs);
}

//******************************************************************************

int main()
{
	CAL_Initialisation(true, true, true);

	bspTree bsp;
	bsp.build(bsp.getRoot());

	std::cout << "Done constructing bsp tree" << std::endl;

	//bsp.traverse();
	
	bsp.testQuery();

	int num;
	std::cin >> num;

	CAL_End();
	return 0;
}