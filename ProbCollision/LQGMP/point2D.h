#ifndef __POINT2D_H__
#define __POINT2D_H__

#define SDIM 2
#define XDIM 2
#define UDIM 2
#define ZDIM 1

#define SAVEFILE "point2.txt"

int cal_environment;
int cal_obstacles;
int cal_goal;
int cal_rrt;
int cal_paths, cal_paths_trunc;
int cal_ellipse, cal_ellipse_trunc;
int cal_point;

int skip = 1;

Matrix<UDIM> u_min, u_max;
Matrix<XDIM> x_min, x_max;
Matrix<XDIM> start;
Matrix<SDIM> goal;
double goalRadius;
double dt;

bspTree2d * btree;

Matrix<XDIM> f(const Matrix<XDIM>& x, const Matrix<UDIM>& u, const Matrix<UDIM>& m) 
{
	Matrix<XDIM> x_new;
	x_new = x;
	x_new[0] += dt*(u[0] + m[0]);
	x_new[1] += dt*(u[1] + m[1]);
	return x_new;
}

void dfdx(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<XDIM, XDIM>& J) {
	J(0,0) = 1; J(0,1) = 0;
	J(1,0) = 0; J(1,1) = 1;                          
}

void dfdu(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<XDIM, UDIM>& J) {
	J(0,0) = dt;  J(0,1) = 0;
	J(1,0) = 0;  J(1,1) = dt;
}

void dfdm(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<XDIM, UDIM>& J) {
	J(0,0) = dt;  J(0,1) = 0;
	J(1,0) = 0;  J(1,1) = dt;
}

void dhdx(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<ZDIM, XDIM>& J) {
	J(0,0) = 1; J(0,1) = 0;
	//J(1,0) = 0; J(1,1) = 1;
}

void dhdn(const Matrix<XDIM>& x, const Matrix<UDIM>& u, Matrix<ZDIM, ZDIM>& J) {
	J = identity<ZDIM>();
}

void initMatrices(Matrix<XDIM, XDIM>& C, Matrix<UDIM, UDIM>& D, Matrix<XDIM, XDIM>& P0, Matrix<UDIM, UDIM>& M, Matrix<ZDIM, ZDIM>& N)
{
	C = identity<XDIM>();
	D = identity<UDIM>();

	P0 = identity<XDIM>() * sqr(0.1);
	M = identity<UDIM>() * sqr(0.07);
	N = identity<ZDIM>() * sqr(0.07);
}

void initEnvironment() 
{
	//start[0] = -8.25; start[1] = 0.0;
	//goal[0] = 8.25; goal[1] = 0.0;
	
	start[0] = -8.25; start[1] = 0.75;
	goal[0] = 8.25; goal[1] = 0.75;
	
	goalRadius = 1;
	dt = 0.5;

	//x_min[0] = -10; x_min[1] = -10;
	//x_max[0] = 10;  x_max[1] = 10;

	x_min[0] = -10; x_min[1] = 0.4;
	x_max[0] = 10;  x_max[1] = 1.0;

	u_min[0] = -1.0; u_min[1] = -1.0;
	u_max[0] = 1.0;  u_max[1] = 1.0;

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
	CAL_CreateCylinder(cal_goal, (float) goalRadius, 0.1f, 0, 0, 0);
	CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], 0.0f);
	CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.6f);
	//int obj;
	//CAL_CreateCylinder(cal_ellipse, 1, 0.1f, 0, 0, 0.1f, &obj);
	//CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);

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