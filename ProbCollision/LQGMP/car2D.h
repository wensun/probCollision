#ifndef __CAR2D_H__
#define __CAR2D_H__

#define SDIM 2
#define XDIM 4
#define UDIM 2
#define ZDIM 3

#define SAVEFILE "carPathX.txt"

int cal_environment;
int cal_obstacles;
int cal_goal;
int cal_rrt;
int cal_paths, cal_paths_trunc;
int cal_ellipse, cal_ellipse_trunc;
int cal_point;
int cal_beacons;

int skip = 5;

Matrix<UDIM> u_min, u_max;
Matrix<XDIM> x_min, x_max;
Matrix<XDIM> start;
Matrix<SDIM> goal;
double goalRadius;
double dt;
double car_l;
double b1x, b1y, b2x, b2y;

bspTree2d * btree;

Matrix<XDIM> f(const Matrix<XDIM>& x, const Matrix<UDIM>& u, const Matrix<UDIM>& m) 
{
	Matrix<XDIM> x_new;
	x_new = x;
	x_new[0] += dt*x[3]*cos(x[2]);
	x_new[1] += dt*x[3]*sin(x[2]);
	x_new[2] += dt*x[3]*tan(u[1] + m[1])/car_l;
	x_new[3] += dt*(u[0] + m[0]);
	return x_new;
}

void dfdx(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<XDIM, XDIM>& J) {
	J(0,0) = 1; J(0,1) = 0; J(0,2) = -dt*x[3]*sin(x[2]); J(0,3) = dt*cos(x[2]);       
	J(1,0) = 0; J(1,1) = 1; J(1,2) = dt*x[3]*cos(x[2]);  J(1,3) = dt*sin(x[2]);       
	J(2,0) = 0; J(2,1) = 0; J(2,2) = 1;                  J(2,3) = dt*tan(u[1])/car_l; 
	J(3,0) = 0; J(3,1) = 0; J(3,2) = 0;                  J(3,3) = 1;                           
}

void dfdu(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<XDIM, UDIM>& J) {
	J(0,0) = 0;  J(0,1) = 0;
	J(1,0) = 0;  J(1,1) = 0;
	J(2,0) = 0;  J(2,1) = dt*x[3]*(1+tan(u[1])*tan(u[1]))/car_l;
	J(3,0) = dt; J(3,1) = 0;
}

void dfdm(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<XDIM, UDIM>& J) {
	J(0,0) = 0;  J(0,1) = 0;
	J(1,0) = 0;  J(1,1) = 0;
	J(2,0) = 0;  J(2,1) = dt*x[3]*(1+tan(u[1])*tan(u[1]))/car_l;
	J(3,0) = dt; J(3,1) = 0;
}

void dhdx(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<ZDIM, XDIM>& J) {
	//J(0,0) = 1; J(0,1) = 0; J(0,2) = 0; J(0,3) = 0;
	//J(1,0) = 0; J(1,1) = 1; J(1,2) = 0; J(1,3) = 0;
	
	double b1 = 1.0/(sqr(sqr(x[0] - b1x) + sqr(x[1] - b1y) + 1.0));
	double b2 = 1.0/(sqr(sqr(x[0] - b2x) + sqr(x[1] - b2y) + 1.0));

	J(0,0) = -2.0*(x[0] - b1x)*b1; J(0,1) = -2.0*(x[1] - b1y)*b1; J(0,2) = 0; J(0,3) = 0;
	J(1,0) = -2.0*(x[0] - b2x)*b2; J(1,1) = -2.0*(x[1] - b2y)*b2; J(1,2) = 0; J(1,3) = 0;
	J(2,0) = 0; J(2,1) = 0; J(2,2) = 0; J(2,3) = 1;
}

void dhdn(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<ZDIM, ZDIM>& J) {
	J = identity<ZDIM>();
}

void initMatrices(Matrix<XDIM, XDIM>& C, Matrix<UDIM, UDIM>& D, Matrix<XDIM, XDIM>& P0, Matrix<UDIM, UDIM>& M, Matrix<ZDIM, ZDIM>& N)
{
	// Initialization
	C = identity<XDIM>();
	D = identity<UDIM>();

	//P0 = identity<XDIM>() * sqr(0.075);
	//M = identity<UDIM>() * sqr(0.05);
	//N = identity<ZDIM>() * sqr(0.05);

	P0 = identity<XDIM>() * sqr(0.075);
	M = identity<UDIM>() * sqr(0.015);
	N = identity<ZDIM>() * sqr(0.01);
}

void initEnvironment() 
{
	//start[0] = -7.5; start[1] = -7.5; start[2] = M_PI*0.25; start[3] = 0.1;
	//goal[0] = 7.5; goal[1] = 7.5;

	//start[0] = 5.5; start[1] = 7.5; start[2] = M_PI*1.5; start[3] = 0.1;
	//start[0] = 2.5; start[1] = -3.5; start[2] = M_PI*0.5; start[3] = 0.1;
	
	start[0] = 7.5; start[1] = 3.5; start[2] = M_PI*0.75; start[3] = 0.1;
	goal[0] = -7.5; goal[1] = -7.5;
	
	goalRadius = 1;

	// carPath1.txt
	dt = 0.2;

	// carPath6.txt
	//dt = 0.3;

	car_l = 1;

	//x_min[0] = -10; x_min[1] = -10; x_min[2] = -0.75*M_PI; x_min[3] = 0;
	//x_max[0] = 10;  x_max[1] = 10;  x_max[2] = 1.25*M_PI;  x_max[3] = 2;

	x_min[0] = -10; x_min[1] = -10; x_min[2] = 0.0*M_PI; x_min[3] = 0;
	x_max[0] = 10;  x_max[1] = 10;  x_max[2] = 2.0*M_PI;  x_max[3] = 2;

	u_min[0] = -0.5; u_min[1] = -M_PI/4;
	u_max[0] = 0.5;  u_max[1] = M_PI/4;

	b1x = -7.5; b1y = 4.5;
	b2x = 4.5; b2y = -7.5;

	// Callisto init
	CAL_SetViewParams(0, 0.0f, 0.0f, 26.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f); 
	//CAL_SetViewParams(0, -4.0f, -2.5f, 12.0f, -4.0f, -2.5f, 0.0f, 0.0f, 1.0f, 0.0f); 

	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	CAL_CreateGroup(&cal_paths, 0, false, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 0);

	CAL_CreateGroup(&cal_paths_trunc, 0, false, "Truncated Paths");
	CAL_SetGroupColor(cal_paths_trunc, 1, 0, 1);
	//CAL_SetGroupColor(cal_paths_trunc, 0.5, 0, 0.5);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	//CAL_SetGroupColor(cal_goal, 0.5, 0.0, 0.5, 0.5);
	CAL_SetGroupColor(cal_goal, 0.0, 0.3, 0.7, 1.0f);
	//CAL_CreateCylinder(cal_goal, (float) goalRadius, 0.05f, 0, 0, 0);
	//CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], -0.05f);
	//CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);
	Vector2d p1, p2, mid, dv;
	double len;
	p1.x = goal[0] + 1.0; p1.y = goal[1] + 0.0;
	for(int i = 1; i < 72; ++i) 
	{
		double theta = (i*5.0*M_DEG2RAD);
		p2.x = goal[0] + goalRadius*cos(theta); p2.y = goal[1] + goalRadius*sin(theta);
		dv = Vector2d(p2.x - p1.x, p2.y - p1.y);
		len = norm(dv);
		dv = dv/len;
		mid = Vector2d(0.5*(p1.x + p2.x), 0.5*(p1.y + p2.y));

		int obj;
		CAL_CreateBox(cal_goal, (float)(len + 0.15), 0.15f, 0.01f, (float)mid.x, (float)mid.y, 0.02f, &obj);
		CAL_SetObjectOrientation(obj, 0.0f, 0.0f, (float)(mod2pi(atan2(dv.y, dv.x))));

		p1 = p2;
	}

	CAL_CreateGroup(&cal_beacons, 0, false, "Beacons");
	//CAL_SetGroupColor(cal_beacons, 0, 0.3, 1.0, 1.0);
	CAL_SetGroupColor(cal_beacons, 0.0, 0.6, 0.5, 1.0);
	CAL_CreateBox(cal_beacons, 0.5f, 0.5f, 0.1f, (float)b1x, (float)b1y, 0.0f);
	CAL_CreateBox(cal_beacons, 0.5f, 0.5f, 0.1f, (float)b2x, (float)b2y, 0.0f);
	
	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	//CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.6f);
	//CAL_SetGroupColor(cal_ellipse, 0, 0.6, 0, 1.0f);
	CAL_SetGroupColor(cal_ellipse, 0.2, 0.2, 0.2, 0.5f);
	//CAL_SetGroupColor(cal_ellipse, 0.0, 0.5, 1.0, 0.5f);
	//int obj;
	//CAL_CreateCylinder(cal_ellipse, 1, 0.1f, 0, 0, 0.1f, &obj);
	//CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "EllipseTruncated");
	//CAL_SetGroupColor(cal_ellipse_trunc, 0, 0, 0.7, 1.0f);
	CAL_SetGroupColor(cal_ellipse_trunc, 1.0, 1.0, 1.0, 1.0f);
	
	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	//CAL_SetGroupColor(cal_obstacles, 0.7f, 0.0f, 0.0f, 1.0f);
	//CAL_SetGroupColor(cal_obstacles, 1.0f, 0.62f, 0.25f, 1.0f);
	CAL_SetGroupColor(cal_obstacles, 1.0f, 0.643f, 0.1f, 1.0f);
	
	CAL_CreateBox(cal_obstacles, 21.0f, 1.0f, 0.1f, 0.0f, -10.0f, -0.05f);
	CAL_CreateBox(cal_obstacles, 1.0f, 21.0f, 0.1f, -10.0f, 0.0f, -0.05);
	CAL_CreateBox(cal_obstacles, 1.0f, 21.0f, 0.1f, 10.0f, 0.0f, -0.05);
	CAL_CreateBox(cal_obstacles, 21.0f, 1.0f, 0.1f, 0.0f, 10.0f, -0.05);

	CAL_CreateBox(cal_obstacles, 3.0f, 5.0f, 0.1f, -9.0f, -3.0f, -0.05);
	CAL_CreateBox(cal_obstacles, 6.0f, 5.0f, 0.1f, 7.0f, -3.0f, -0.05);
	CAL_CreateBox(cal_obstacles, 5.0f, 6.0f, 0.1f, -3.0f, 7.0f, -0.05);
	CAL_CreateBox(cal_obstacles, 5.0f, 3.0f, 0.1f, -3.0f, -9.0f, -0.05);

	//CAL_CreateBox(cal_obstacles, 5.0f, 15.0f, 0.1f, -3.0f, 7.0f, -0.05);

	CAL_CreateBox(cal_obstacles, 5.0f, 5.0f, 0.1f, -3.0f, -3.0f, -0.05);

	// init BSP tree
	btree = new bspTree2d();
	bspNode2d * root = btree->getRoot();

	int v[24];

	v[0] = btree->addVertex(-9.5, -9.5);
	v[1] = btree->addVertex(-5.5, -9.5);
	v[2] = btree->addVertex(-5.5, -7.5);
	v[3] = btree->addVertex(-0.5, -7.5);
	v[4] = btree->addVertex(-0.5, -9.5);
	v[5] = btree->addVertex(9.5, -9.5);
	v[6] = btree->addVertex(9.5, -5.5);
	//v[7] = btree->addVertex(7.5, -5.5);
	//v[8] = btree->addVertex(7.5, -0.5);

	v[7] = btree->addVertex(4.0, -5.5);
	v[8] = btree->addVertex(4.0, -0.5);
	
	v[9] = btree->addVertex(9.5, -0.5);
	v[10] = btree->addVertex(9.5, 9.5);
	v[11] = btree->addVertex(-0.5, 9.5);
	//v[12] = btree->addVertex(-0.5, 7.5);
	//v[13] = btree->addVertex(-5.5, 7.5);
	
	v[12] = btree->addVertex(-0.5, 4.0);
	v[13] = btree->addVertex(-5.5, 4.0);
	
	v[14] = btree->addVertex(-5.5, 9.5);
	v[15] = btree->addVertex(-9.5, 9.5);
	v[16] = btree->addVertex(-9.5, -0.5);
	v[17] = btree->addVertex(-7.5, -0.5);
	v[18] = btree->addVertex(-7.5, -5.5);
	v[19] = btree->addVertex(-9.5, -5.5);

	v[20] = btree->addVertex(-5.5, -5.5);
	v[21] = btree->addVertex(-5.5, -0.5);
	v[22] = btree->addVertex(-0.5, -0.5);
	v[23] = btree->addVertex(-0.5, -5.5);

	int p = 0;
	bspEdge2d * e;
	for(int i = 0; i < 19; ++i) {
		e = new bspEdge2d(btree, v[i], v[i+1]); p = btree->addEdge(e); root->addEdge(p);
	}
	e = new bspEdge2d(btree, v[19], v[0]); p = btree->addEdge(e); root->addEdge(p);
	for(int i = 20; i < 23; ++i) {
		e = new bspEdge2d(btree, v[i], v[i+1]); p = btree->addEdge(e); root->addEdge(p); 
	}
	e = new bspEdge2d(btree, v[23], v[20]); p = btree->addEdge(e); root->addEdge(p); 

	/*
	int ne = btree->getNumEdges();
	for(int i = 0; i < ne; ++i) 
	{
		bspEdge2d* e = btree->getEdge(i);
		Vector2d* mid = e->getMidPt();
		Vector2d* dv = e->getDirVec();

		int obj;
		CAL_CreateBox(cal_obstacles, (float)(e->getLength()+0.2), 0.2f, 0.2f, (float)mid->x, (float)mid->y, 0.0f, &obj);
		CAL_SetObjectOrientation(obj, 0.0f, 0.0f, (float)(mod2pi(atan2(dv->y, dv->x))));
	}
	*/

	root->computeBarycenter();
	btree->build(root);

	std::cout << "Done constructing bsp tree" << std::endl;
	//btree->traverse();
}

void initEnvironment2() 
{
	start[0] = -8.25; start[1] = 0.0; start[2] = M_PI*0.0; start[3] = 0.1;
	goal[0] = 8.25; goal[1] = 0.0;
	goalRadius = 1;
	dt = 1.0;

	// car length
	car_l = 1;

	x_min[0] = -10; x_min[1] = -10; x_min[2] = -0.75*M_PI; x_min[3] = 0;
	x_max[0] = 10;  x_max[1] = 10;  x_max[2] = 1.25*M_PI;  x_max[3] = 2;

	u_min[0] = -0.5; u_min[1] = -M_PI/4;
	u_max[0] = 0.5;  u_max[1] = M_PI/4;

	CAL_SetViewParams(0, 0.0f, 0.0f, 26.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f); 

	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	CAL_CreateGroup(&cal_paths, 0, false, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 1);

	CAL_CreateGroup(&cal_paths_trunc, 0, false, "Truncated Paths");
	CAL_SetGroupColor(cal_paths_trunc, 1, 0, 0);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	CAL_SetGroupColor(cal_goal, 0, 1, 1, 0.5);
	CAL_CreateCylinder(cal_goal, (float) goalRadius, 0.02f, 0, 0, 0);
	CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], 0.0f);
	CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.6f);
		
	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "EllipseTruncated");
	CAL_SetGroupColor(cal_ellipse_trunc, 0, 0, 1, 0.6f);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0.7f, 0.0f, 0.0f, 1.0f);

	CAL_CreateBox(cal_obstacles, 21.0f, 1.0f, 0.1f, 0.0f, -10.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 1.0f, 21.0f, 0.1f, -10.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 1.0f, 21.0f, 0.1f, 10.0f, 0.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 21.0f, 1.0f, 0.1f, 0.0f, 10.0f, 0.0f);

	CAL_CreateBox(cal_obstacles, 14.0f, 6.0f, 0.1f, 0.0f, 4.0f, 0.0f);
	CAL_CreateBox(cal_obstacles, 14.0f, 6.0f, 0.1f, 0.0f, -4.0f, 0.0f);

	// init BSP tree
	btree = new bspTree2d();
	bspNode2d * root = btree->getRoot();

	int v[12];

	v[0] = btree->addVertex(-9.5, -9.5);
	v[1] = btree->addVertex(9.5, -9.5);
	v[2] = btree->addVertex(9.5, 9.5);
	v[3] = btree->addVertex(-9.5, 9.5);

	v[4] = btree->addVertex(-7.0, 1.0);
	v[5] = btree->addVertex(-7.0, 7.0);
	v[6] = btree->addVertex(7.0, 7.0);
	v[7] = btree->addVertex(7.0, 1.0);
	
	v[8] = btree->addVertex(-7.0, -1.0);
	v[9] = btree->addVertex(7.0, -1.0);
	v[10] = btree->addVertex(7.0, -7.0);
	v[11] = btree->addVertex(-7.0, -7.0);

	int p = 0;
	bspEdge2d * e;
	for(int i = 0; i < 3; ++i) {
		e = new bspEdge2d(btree, v[i], v[i+1]); p = btree->addEdge(e); root->addEdge(p); 
	}
	e = new bspEdge2d(btree, v[3], v[0]); p = btree->addEdge(e); root->addEdge(p); 
	for(int i = 4; i < 7; ++i) {
		e = new bspEdge2d(btree, v[i], v[i+1]); p = btree->addEdge(e); root->addEdge(p); 
	}
	e = new bspEdge2d(btree, v[7], v[4]); p = btree->addEdge(e); root->addEdge(p); 
	for(int i = 8; i < 11; ++i) {
		e = new bspEdge2d(btree, v[i], v[i+1]); p = btree->addEdge(e); root->addEdge(p); 
	}
	e = new bspEdge2d(btree, v[11], v[8]); p = btree->addEdge(e); root->addEdge(p); 
	
	
	root->computeBarycenter();
	btree->build(root);

	std::cout << "Done constructing bsp tree" << std::endl;
}

void cleanup() {
	delete btree;
	CAL_End();
}

#endif