#include "needle3D.h"
#include "utils.h"

inline Matrix<4,4> stepT(const Matrix<4,4>& T, const Matrix<3>& v, const Matrix<3>& w, double dt) {
	Matrix<4,4> U = zeros<4,4>();
	U.insert(0,0, cpMatrix(w));
	U.insert(0,3, v);
	return T*exp(dt*U);
}

inline void createInputs(Matrix<3>& v, Matrix<3>& w, double twist, double curvature) 
{
	v(0,0) = 0;
	v(1,0) = 0;
	v(2,0) = speed; 

	w(0,0) = speed * curvature;
	w(1,0) = 0;
	w(2,0) = twist;  
}

inline void createABN(Matrix<6,6>& A, Matrix<6,3>& B, Matrix<6,6>& N, const Matrix<3>& v, const Matrix<3>& w, const Matrix<6,6>& M, double dt) {
	Matrix<6,6> F = zeros<6,6>();
	F.insert(0,0, -cpMatrix(w));
	F.insert(0,3, -cpMatrix(v));
	F.insert(3,3, -cpMatrix(w));

	Matrix<6,3> G = zeros<6,3>();
	G(2,0) = 1;
	G(3,2) = 1;
	G(5,1) = 1;

	A = exp(dt*F);
	B = (dt/6.0) * (G + 4*exp(dt*F*0.5)*G + A*G); // trapezoid integration approximation
	N = (dt/6.0) * (M + 4*exp(dt*F*0.5)*M*~exp(dt*F*0.5) + A*M*~A); // trapezoid integration approximation
}

inline double computeConfidence(const Matrix<3>& pos, const Matrix<3,3>& R) 
{
	Matrix<3,3> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<4,1> q = quatFromRot(~EVec);

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1/(float)sqrt(EVal(2,2)));

	Matrix<3,3> invScale = zeros<3,3>();
	invScale(0,0) = 1/sqrt(EVal(0,0));
	invScale(1,1) = 1/sqrt(EVal(1,1));
	invScale(2,2) = 1/sqrt(EVal(2,2));
	Matrix<3> transPos =  invScale * ~EVec * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], (float) transPos[2]);

	int num_pairs;
	CAL_GetClosestPairs (cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete[] results;

	CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
	CAL_SetGroupScaling(cal_environment, 1, 1, 1);

	return distance;
}

inline double dist(const Matrix<4,4>& T, const Matrix<3>& point) 
{
	Matrix<3> proj_point = ~T.subMatrix<3,3>(0,0) * (point - T.subMatrix<3,1>(0,3));
	double y = proj_point(2,0);
	if (y > 0) { 
		return INFTY;
	}
	double x = sqrt(proj_point(0,0)*proj_point(0,0) + proj_point(1,0)*proj_point(1,0));
	if (x == 0) {
		return -y;
	}
	double r = (x*x + y*y) / (2*x);
	if (r < needleSteeringRadius) {
		return INFTY;
	}
	//return sqrt(x*x + y*y);
	return atan2(-y, r-x) * r;
}

int nearestNeighbor(const Matrix<3>& point, const std::vector<TreeNode>& tree) 
{
	int closest = -1;
	double mindist = INFTY;
	int col = 1;
	CAL_CheckLineCollision(cal_obstacles, (float)point(0,0), (float)point(1,0), (float)point(2,0), (float)tree[0].T(0,3), (float)tree[0].T(1,3), (float)tree[0].T(2,3), false, &col);
	if (col == 0) {
		closest = 0;
		mindist = factor*sqrt(tr(~(tree[0].T.subMatrix<3,1>(0,3) - point) * (tree[0].T.subMatrix<3,1>(0,3) - point)));
	} 

	std::stack<int> st;
	for (int i = 0; i < (int) tree[0].children.size(); ++i) {
		st.push(tree[0].children[i]);
	}

	while (!st.empty()) {
		int i = st.top();
		st.pop();
		double d = dist(tree[i].T, point);
		if (tree[i].T(2,3) < 0) { // do not expand nodes outside entry plane
			d = INFTY;
		}
		if (d < INFTY) {
			if (factor*d + tree[i].depth * speed * dt < mindist) {
				col = 1;
				CAL_CheckLineCollision(cal_obstacles, (float)point(0,0), (float)point(1,0), (float)point(2,0), 
					(float)tree[i].T(0,3), (float)tree[i].T(1,3), (float)tree[i].T(2,3), 
					false, &col);
				if (col == 0) {
					closest = i;
					mindist = factor*d + tree[i].depth * speed * dt;
				}
			}
			for (int c = 0; c < (int) tree[i].children.size(); ++c) {
				st.push(tree[i].children[c]);
			}
		}
	}

	return closest;
}

void initObstacles() 
{
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles", false);
	CAL_SetGroupColor(cal_obstacles, 1.0f, 0.643f, 0.1f, 1.0f);

	//CAL_CreateBox(cal_obstacles, (float)Xrange, (float)(Yrange * 0.25), 0.5f, (float)(0.5*Xrange), (float)(0.125*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)Xrange, (float)(Yrange * 0.25), 0.5f, (float)(0.5*Xrange), (float)(0.875*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.25), (float)Yrange, 0.5f, (float)(0.125*Xrange), (float)(0.5*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.25), (float)Yrange, 0.5f, (float)(0.875*Xrange), (float)(0.5*Yrange), (float)(0.75*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.4), (float)(Yrange * 0.4), 0.5f, (float)(0.55*Xrange), (float)(0.55*Yrange), (float)(0.75*Zrange));

	//CAL_CreateBox(cal_obstacles, (float)Xrange, (float)(Yrange * 0.25), 0.5f, (float)(0.5*Xrange), (float)(0.125*Yrange), (float)(0.25*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)Xrange, (float)(Yrange * 0.25), 0.5f, (float)(0.5*Xrange), (float)(0.875*Yrange), (float)(0.25*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.25), (float)Yrange, 0.5f, (float)(0.125*Xrange), (float)(0.5*Yrange), (float)(0.25*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.25), (float)Yrange, 0.5f, (float)(0.875*Xrange), (float)(0.5*Yrange), (float)(0.25*Zrange));
	//CAL_CreateBox(cal_obstacles, (float)(Xrange * 0.4), (float)(Yrange * 0.4), 0.5f, (float)(0.45*Xrange), (float)(0.45*Yrange), (float)(0.25*Zrange));

	//CAL_CreateBox(cal_obstacles, 2.0f, 2.0f, 2.0f, 5.0f, 5.0f, 7.0f);

	// Polygons
	Matrix<3> v1, v2, a, intpt;
	double b, len;
	int npoly, npts;

	std::ifstream fptr("env3.txt", std::ios::in);
	fptr >> npoly;

	for(int i = 0; i < npoly; ++i) {
		fptr >> npts;
		float * p = new float[3*npts];
		for(int j = 0; j < npts; ++j) {
			fptr >> p[3*j] >> p[3*j+1] >> p[3*j+2];
		}
		
		CAL_CreatePolygon(cal_obstacles, npts, p);
		v1[0] = p[6] - p[3]; v1[1] = p[7] - p[4]; v1[2] = p[8] - p[5];
		v2[0] = p[0] - p[3]; v2[1] = p[1] - p[4]; v2[2] = p[2] - p[5];
		a[0] = v1[1]*v2[2] - v2[1]*v1[2]; 
		a[1] = v2[0]*v1[2] - v1[0]*v2[2];
		a[2] = v1[0]*v2[1] - v2[0]*v1[1];
		len = sqrt(tr(a*~a));
		a = -a / len;
		b = a[0]*p[3] + a[1]*p[4] + a[2]*p[5];
		intpt[0] = 0.5*(p[0] + p[6]); intpt[1] = 0.5*(p[1] + p[7]); intpt[2] = 0.5*(p[2] + p[8]);
		obs.push_back(Obstacle(a, b, intpt));
		delete[] p;
	}

	if (fptr) { fptr.close(); }
}

void initEnvironment() 
{
	Xrange = 10.0;
	Yrange = 10.0;
	Zrange = 10.0;

	speed = 1.0;
	needleSteeringRadius = 5.0;
	maxTwist = 2.0*M_PI;
	maxCurvature = 1.0/needleSteeringRadius;

	dt = 0.25;
	factor = 1.3;

	goal[0] = 0.5 * Xrange; goal[1] = 0.5 * Yrange; goal[2] = 1.0 * Zrange;

	// Callisto init
	//CAL_SetViewParams(0, 23.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 0.0f, 1.0f, 0.0f);
	CAL_SetViewParams(0, 23.5f, 6.5f, 11.5f, 5.0f, 5.0f, 5.0f, 0.0f, 1.0f, 0.0f);

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	
	initObstacles();

	// Bounding box
	CAL_CreateGroup (&cal_box, 0, false, "Box");
	CAL_SetGroupColor(cal_box, 0, 0, 0);
	//int np[1] = {2};
	//float side1[6] = {0, 0, 0, (float)Xrange, 0, 0};
	//float side2[6] = {0, 0, 0, 0, (float)Yrange, 0};
	//float side3[6] = {0, 0, 0, 0, 0, (float)Zrange};
	//float side4[6] = {(float)Xrange, (float)Yrange, (float)Zrange, 0, (float)Yrange, (float)Zrange};
	//float side5[6] = {(float)Xrange, (float)Yrange, (float)Zrange, (float)Xrange, 0, (float)Zrange};
	//float side6[6] = {(float)Xrange, (float)Yrange, (float)Zrange, (float)Xrange, (float)Yrange, 0};
	//float side7[6] = {0, (float)Yrange, (float)Zrange, 0, 0, (float)Zrange};
	//float side8[6] = {0, (float)Yrange, (float)Zrange, 0, (float)Yrange, 0};
	//float side9[6] = {(float)Xrange, 0, (float)Zrange, 0, 0, (float)Zrange};
	//float side10[6] = {(float)Xrange, 0, (float)Zrange, (float)Xrange, 0, 0};
	//float side11[6] = {(float)Xrange, (float)Yrange, 0, 0, (float)Yrange, 0};
	//float side12[6] = {(float)Xrange, (float)Yrange, 0, (float)Xrange, 0, 0};
	//CAL_CreatePolyline(cal_box, 1, np, side1);
	//CAL_CreatePolyline(cal_box, 1, np, side2);
	//CAL_CreatePolyline(cal_box, 1, np, side3);
	//CAL_CreatePolyline(cal_box, 1, np, side4);
	//CAL_CreatePolyline(cal_box, 1, np, side5);
	//CAL_CreatePolyline(cal_box, 1, np, side6);
	//CAL_CreatePolyline(cal_box, 1, np, side7);
	//CAL_CreatePolyline(cal_box, 1, np, side8);
	//CAL_CreatePolyline(cal_box, 1, np, side9);
	//CAL_CreatePolyline(cal_box, 1, np, side10);
	//CAL_CreatePolyline(cal_box, 1, np, side11);
	//CAL_CreatePolyline(cal_box, 1, np, side12);
	float w = 0.07f;
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, 0.0f, 0.0f);
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, (float)Yrange, 0.0f);
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, (float)Yrange, (float)Zrange);
	CAL_CreateBox(cal_box, (float)Xrange, w, w, 0.5*(float)Xrange, 0.0f, (float)Zrange);

	CAL_CreateBox(cal_box, w, (float)Yrange, w, 0.0f, 0.5*(float)Yrange, 0.0f);
	CAL_CreateBox(cal_box, w, (float)Yrange, w, (float)Xrange, 0.5*(float)Yrange, 0.0f);
	CAL_CreateBox(cal_box, w, (float)Yrange, w, (float)Xrange, 0.5*(float)Yrange, (float)Zrange);
	CAL_CreateBox(cal_box, w, (float)Yrange, w, 0.0f, 0.5*(float)Yrange, (float)Zrange);

	CAL_CreateBox(cal_box, w, w, (float)Zrange, 0.0f, 0.0f, 0.5*(float)Zrange);
	CAL_CreateBox(cal_box, w, w, (float)Zrange, (float)Xrange, 0.0f, 0.5*(float)Zrange);
	CAL_CreateBox(cal_box, w, w, (float)Zrange, (float)Xrange, (float)Yrange, 0.5*(float)Zrange);
	CAL_CreateBox(cal_box, w, w, (float)Zrange, 0.0f, (float)Yrange, 0.5*(float)Zrange);

	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 0, 1, 0);

	CAL_CreateGroup(&cal_paths, 0, false, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 1);

	CAL_CreateGroup(&cal_paths_trunc, 0, false, "Truncated Paths");
	CAL_SetGroupColor(cal_paths_trunc, 1, 0, 1);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	CAL_SetGroupColor(cal_goal, 0.0, 0.3, 0.7, 0.4f);
	CAL_CreateSphere(cal_goal, 0.5f, (float) goal[0], (float) goal[1], (float) goal[2]);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0.2, 0.2, 0.2, 0.6f);
	//int obj;
	//CAL_CreateCylinder(cal_ellipse, 1, 0.1f, 0, 0, 0.1f, &obj);
	//CAL_SetObjectOrientation(obj, (float) (M_PI*0.5), 0, 0);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "EllipseTruncated");
	CAL_SetGroupColor(cal_ellipse_trunc, 0, 0, 0, 1.0);
}

void rrt(std::ofstream & fptr) 
{
	std::vector<TreeNode> tree;
	//CAL_EmptyGroup(cal_rrt);

	TreeNode n;
	n.T = identity<4>();
	n.T.insert(0,3, goal);
	n.bp = -1;
	n.depth = -1;
	tree.push_back(n);

	bool found = false;
	
	int node;
	Matrix<3> v; 
	Matrix<3> w;
		
	while (!found) 
	{
		Matrix<3> point; 
		do {
			if (uniform() < 0.2) {
				point(0,0) = Xrange * uniform(); point(1,0) = Yrange * uniform(); point(2,0) = -1.0;
			} else {
				point(0,0) = Xrange * uniform(); point(1,0) = Yrange * uniform(); point(2,0) = Zrange * uniform();
			}
			node = nearestNeighbor(point, tree);
		} while (node == -1);
	
		if (node == 0) {
			TreeNode goalnode;
			goalnode.bp = 0;
			goalnode.depth = 0;
			goalnode.T = tree[0].T;
			goalnode.w = 0;

			Matrix<3> z = (tree[0].T.subMatrix<3,1>(0,3) - point); 
			z = z / sqrt(tr(~z*z)); 
			Matrix<3> random_point; 
			do {
				random_point(0,0) = uniform(); random_point(1,0) = uniform(); random_point(2,0) = uniform();
			} while (tr(~random_point * random_point) > 1);
			Matrix<3> y = random_point - z * tr(~random_point * z);
			y = y / sqrt(tr(~y*y));
			Matrix<3> x = cpMatrix(y) * z;

			goalnode.T.insert(0,0, x);
			goalnode.T.insert(0,1, y);
			goalnode.T.insert(0,2, z);

			node = (int) tree.size();
			tree[0].children.push_back(node);
			tree.push_back(goalnode);
		}
		
		//createInputs(v, w, (2.0 * maxTwist * uniform()) - maxTwist, maxCurvature * uniform());
		createInputs(v, w, (2.0 * maxTwist * uniform()) - maxTwist, maxCurvature);
		
		TreeNode newnode;
		newnode.bp = node;
		newnode.depth = tree[node].depth + speed*dt;
		newnode.T = stepT(tree[node].T, v, w, -dt);
		newnode.w = w(2,0);
		newnode.k = w(0,0)/speed;
		newnode.v = speed;

		float line[6] = {(float) newnode.T(0,3), (float) newnode.T(1,3), (float) newnode.T(2,3), (float) tree[newnode.bp].T(0,3), (float) tree[newnode.bp].T(1,3), (float) tree[newnode.bp].T(2,3)};
		int col = 1;
		CAL_CheckLineCollision(cal_obstacles, line[0], line[1], line[2], line[3], line[4], line[5], false, &col);
		if (col != 0) { // collision
			if (tree[newnode.bp].bp == 0) {
				tree.pop_back();
				tree[0].children.pop_back();
			}
		} else {
			int np[1] = {2};
			CAL_CreatePolyline(cal_rrt, 1, np, line);

			tree[node].children.push_back((int) tree.size());
			tree.push_back(newnode);

			if (newnode.T(2,3) < 0.0) {
				found = true;
			}
		}
	}

	//int k;
	//std::cin >> k;

	std::vector<PathNode> path;
	int i = (int) tree.size() - 1;

	while (i != 0) {
		PathNode stage;

		stage.T = tree[i].T;

		stage.w = tree[i].w;
		stage.v = tree[i].v;
		stage.k = tree[i].k;
		
		path.push_back(stage);

		i = tree[i].bp;
	}

	fptr << (int) path.size();
	for (int i = 0; i < (int) path.size(); ++i) {
		path[i].serialize(fptr);
	}
	fptr << std::endl;
	
	//std::cout << "Finished RRT" << std::endl;
}

//-------------------------------------------------------------------------------------------------

/*
void query2(const Matrix<3>& x, const Matrix<3, 3>& S, std::vector<std::pair<Matrix<3>, double> >& cvx)
{
	Matrix<3,3> V, E;
	jacobi(S, V, E);

	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~V);
	Matrix<4,1> q = quatFromRot(Temp);

	double scale = 1.0;
	double sx = scale*sqrt(E(0,0));
	double sy = scale*sqrt(E(1,1));
	double sz = scale*sqrt(E(2,2));
	Matrix<3,3> invScale = zeros<3,3>();
	invScale(0,0) = 1/sx;
	invScale(1,1) = 1/sy;
	invScale(2,2) = 1/sz;
	Matrix<3,3> T = invScale * ~V; // ? == !(V * !invscale)
	Matrix<3> transPos =  T * x;
	Matrix<3,3> invT = !T;

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, (float) invScale(0,0), (float) invScale(1,1), (float) invScale(2,2));

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], (float) transPos[2]);

	int num_pairs;
	int retcode = CAL_GetClosestPairs (cal_point, cal_environment, &num_pairs);

	//std::cout << "Ret code for CAL_GetClosestPairs: " << retcode << std::endl;

	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	for(int it1 = 0; it1 < num_pairs; ++it1) 
	{
		Matrix<3> p0, p1, vec, normal;
		p0[0] = (double)results[it1].vector0[0];
		p0[1] = (double)results[it1].vector0[1];
		p0[2] = (double)results[it1].vector0[2];

		p1[0] = (double)results[it1].vector1[0];
		p1[1] = (double)results[it1].vector1[1];
		p1[2] = (double)results[it1].vector1[2];

		vec = (p1 - p0);
		double vecnorm = tr(~vec*vec);
		vec = vec/vecnorm;

		normal[0] = p1[0] + vec[1]; 
		normal[1] = p1[1] - vec[0];
		normal[2] = p1[2] - vec[0];

		normal = invT*normal; p1 = invT*p1;

		vec = (normal - p1);
		vecnorm = sqrt(tr(~vec*vec));
		vec = vec/vecnorm;

		normal[0] = -vec[1]; normal[1] = vec[0];

		Matrix<3> a = normal;
		double b = tr(~normal*p1);

		cvx.push_back(std::make_pair(a, b));
	}

	delete[] results;

	//std::cout << "Num truncation constraints: " << (int)cvx.size() << std::endl;

	CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
	CAL_SetGroupScaling(cal_environment, 1, 1, 1);
}
*/

void query(const Matrix<3>& pos, const Matrix<3,3>& R, std::vector<std::pair<Matrix<3>, double> >& cvx)
{
	Matrix<3,3> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<4,1> q = quatFromRot(~EVec);

	//std::cout << "pos: " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1/(float)sqrt(EVal(2,2)));

	Matrix<3,3> invScale = zeros<3,3>();
	invScale(0,0) = 1/sqrt(EVal(0,0));
	invScale(1,1) = 1/sqrt(EVal(1,1));
	invScale(2,2) = 1/sqrt(EVal(2,2));
	Matrix<3,3> T = invScale * ~EVec;
	Matrix<3,3> invT = !T;
	Matrix<3> transPos =  T * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], (float) transPos[2]);

	int num_pairs;
	CAL_GetClosestPairs (cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	//for(int it1 = 0; it1 < num_pairs; ++it1) 
	//{
	//	if (results[it1].distance < 3.5) {
	//		std::cout << "Polygon: " << results[it1].objID1 << " dist: " << results[it1].distance << std::endl;
	//	}
	//}
	
	std::set<int> valid, invalid;
	for(int it1 = 0; it1 < num_pairs; ++it1) 
	{
		if (invalid.find(it1) != invalid.end()) 
			continue;

		if (results[it1].distance > 3.5)
			break;

		// hack!
		if (results[it1].distance == results[it1+1].distance) 
		{
			Matrix<3> p0, p1;
			p0[0] = (double)results[it1].vector0[0];
			p0[1] = (double)results[it1].vector0[1];
			p0[2] = (double)results[it1].vector0[2];

			p1[0] = (double)results[it1].vector1[0];
			p1[1] = (double)results[it1].vector1[1];
			p1[2] = (double)results[it1].vector1[2];

			p0 = invT*p0; p1 = invT*p1;

			Matrix<3> n = (p0 - p1);
			n = -n/sqrt(tr(~n*n));

			cvx.push_back(std::make_pair(n, tr(~n*p1)));			
			//std::cout << "n: " << n[0] << " " << n[1] << " " << n[2] << " " << tr(~n*p1) << std::endl;

			invalid.insert(it1);
			invalid.insert(it1+1);
		} // equidistant
		else
		{
			//std::cout << "Polygon: " << results[it1].objID1 << " dist: " << results[it1].distance << std::endl;
			valid.insert(results[it1].objID1);

			Obstacle& o1 = obs[results[it1].objID1];

			for(int it2 = it1 + 1; it2 < num_pairs; ++it2) {
				if (results[it2].distance > 3.5)
					break;
			
				Obstacle& o2 = obs[results[it2].objID1];
				if ( tr(~o1.a*o2.intpt) >= o1.b) {
					invalid.insert(it2);
				}
			}
		} // not equidistant
	}

	delete[] results;

	//std::cout << "Num collision pairs: " << (int)valid.size() << std::endl;

	for(std::set<int>::iterator it = valid.begin(); it != valid.end(); ++it) {
		Obstacle& o = obs[*it];
		cvx.push_back(std::make_pair(o.a, o.b));
		//std::cout << "a: " << o.a[0] << " " << o.a[1] << " " << o.a[2] << " " << o.b << std::endl;
	}

	//std::cout << "Num truncation planes: " << (int)cvx.size() << std::endl;
	
	CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
	CAL_SetGroupScaling(cal_environment, 1, 1, 1);
}

double probSuccessBoole(const Matrix<3>& x, const Matrix<3,3>& S, std::vector<std::pair<Matrix<3>, double> >& cvx) 
{
	double pb = 0.0;
	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; ++i) {
		Matrix<3>& a = cvx[i].first;
		double b = cvx[i].second;
		pb += (1.0 - cdf((b - tr(~a*x))/(sqrt(tr(~a*S*a)))));
	}
	return (1.0 - pb);
}

double probSuccess(const Matrix<3>& x, const Matrix<3,3>& S, std::vector<std::pair<Matrix<3>, double> >& cvx) 
{
	double ps = 1.0;
	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; ++i) {
		Matrix<3>& a = cvx[i].first;
		double b = cvx[i].second;
		ps *= cdf((b - tr(~a*x))/(sqrt(tr(~a*S*a))));
	}
	return ps;
}

double lqgmpTruncation(const Matrix<12>& x, const Matrix<12,12>& S, const Matrix<3>& phat, const Matrix<3,3>& Rhat,
					   Matrix<12>& xTrunc, Matrix<12,12>& STrunc, std::vector<std::pair<Matrix<3>, double> >& cvx)
{
	xTrunc = x;
	STrunc = S;

	double ps = 1.0;

	Matrix<3> xp = phat + Rhat*x.subMatrix<3,1>(0,0);
	Matrix<3,3> Sp = Rhat * S.subMatrix<3,3>(0,0) * ~Rhat;

	Matrix<3> xDelta;
	Matrix<3,3> SDelta;

	xDelta.reset();
	SDelta.reset();

	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; ++i) {
		std::pair<Matrix<3>, double>& c = cvx[i];

		Matrix<3> a;
		a = c.first;

		double b = c.second;

		double yMean, yVar, yNewMean, yNewVar;
		Matrix<3> xyCovar, L;
	
		yMean = tr(~a*xp);
		yVar = tr(~a*Sp*a);

		truncate(b, yMean, yVar, yNewMean, yNewVar);

		xyCovar = Sp*a;
		L = (xyCovar/yVar);

		double ratio1 = (yNewVar/yVar);
		xDelta += L*(yNewMean - yMean);

		SDelta += ((ratio1 - 1.0)*yVar)*(L*~L);

		double ratio2 = cdf((b - tr(~a*xp))/(sqrt(tr(~a*Sp*a))));
		//SDelta += ((ratio2 - 1.0)*yVar)*(L*~L);

		ps *= ratio2;
		//std::cout << "yNewMean/yMean: " << yNewMean/yMean << " yNewVar/yVar: " << yNewVar/yVar << " prob mass: " << cdf((b - tr(~a*xp))/(sqrt(tr(~a*Sp*a)))) << std::endl;
	}
	
	xp += xDelta;
	Sp += SDelta;

	Matrix<3> pbar = ~Rhat*(xp - phat);
	Matrix<3,3> Sbar = ~Rhat*Sp*Rhat;
	
	// Conditional
	Matrix<12,3> Lambda;
	Lambda.reset();
	Lambda.insert<3,3>(0,0,identity<3>());
	
	Matrix<3,3> SigmaY = ~Lambda * S * Lambda;
	Matrix<12,3> SigmaXY = S * Lambda;

	Matrix<12,3> K = SigmaXY * !SigmaY;

	xTrunc += K*(pbar - ~Lambda*x);
	STrunc += K*(Sbar - ~Lambda*S*Lambda)*~K;

	checkNaN(xTrunc);
	checkNaN(STrunc);

	return ps;
}

void preprocess(const std::vector<PathNode>& path) 
{
	int l = (int) path.size();

	A.resize(l); // process matrix
	B.resize(l); // input matrix
	N.resize(l); // process noise matrix
	H.resize(l); // measurement matrix
	L.resize(l); // feedback matrix
	K.resize(l); // Kalman-gain matrix

	P0 = identity<6>() * 0.05;
	Q = identity<ZDIM>() * 0.1;
	M = identity<6>() * 0.001;
	C = identity<6>();
	D = identity<3>();

	// Jacobians
	Matrix<3> v;
	Matrix<3> w;
	l = (int) path.size() - 1;
	for (int k = 0; k <= l; ++k) {
		createInputs(v, w, path[k].w, path[k].k);
		createABN(A[k], B[k], N[k], v, w, M, dt);
	}

	// LQR
	Matrix<6,6> S;
	S = C;
	//L[l - 1] = zeros<3,6>();
	for (int k = l - 1; k != -1; --k) {
		L[k] = -!(~B[k]*S*B[k] + D)*~B[k]*S*A[k];
		S = C + ~A[k]*S*A[k] + ~A[k]*S*B[k]*L[k];
	}

	// Kalman
	Matrix<6,6> P;
	P = P0;
	//K[0] = zeros<6, ZDIM>();
	for (int k = 1; k <= l; ++k) {
		P = A[k-1]*P*~A[k-1] + N[k-1];

		H[k].reset();

		// dependent on measurement model !!!!!!
		H[k].insert(0,0, path[k].T.subMatrix<1,3>(0,0));
		H[k].insert(1,0, path[k].T.subMatrix<1,3>(1,0));
		H[k].insert(2,0, path[k].T.subMatrix<1,3>(2,0));

		K[k] = P*~H[k]*!(H[k]*P*~H[k] + Q);
		P = (identity<6>() - K[k]*H[k])*P;
	}
}

double simulate(const std::vector<PathNode>& path, int num_samples, bool vis = false) 
{
	int l = (int) path.size() -1;

	preprocess(path);

	// Simulate
	Matrix<6> x_est;
	Matrix<4,4> T_true;
	
	Matrix<3> u;
	Matrix<6> m;

	Matrix<ZDIM> n;
	Matrix<ZDIM> z;
	Matrix<6,6> P;

	int fail = 0;

	int* cal_samples;
	if (vis) {
		cal_samples = new int[num_samples];
	}

	Timer t;
	t.start();

	int progress = 0;
	for (int s = 0; s < num_samples; ++s) 
	{
		//if ((100 * s) / num_samples != progress) {
		//	progress = (100 * s) / num_samples;
		//	std::cout << progress << " ";
		//}

		if (vis) {
			CAL_CreateGroup(&(cal_samples[s]), 0, false);
			CAL_SetGroupColor(cal_samples[s], 1, 0, 0);
		}

		P = P0; 
		x_est = zeros<6,1>();
		T_true = path[0].T * transFromErr(sampleGaussian(x_est, P0));

		float pos[3] = {(float)T_true(0,3), (float)T_true(1,3), (float)T_true(2,3)};

		if (vis) {
			int obj;
			CAL_CreateSphere(cal_samples[s], 0.05f, pos[0], pos[1], pos[2], &obj);
		}

		for (int k = 1; k <= l; ++k) 
		{
			Matrix<3,1> v;
			Matrix<3,1> w;
			createInputs(v, w, path[k-1].w, path[k-1].k);
			
			u = L[k-1]*x_est;
			
			v(2,0) += u(0,0);
			w(2,0) += u(1,0);
			w(0,0) += u(2,0);

			m = sampleGaussian(zeros<6,1>(), N[k-1]);
			v += m.subMatrix<3,1>(0,0);
			w += m.subMatrix<3,1>(3,0);

			// Cache earlier values for collision checks
			double xprev = T_true(0,3), yprev = T_true(1,3), zprev = T_true(2,3);

			T_true = stepT(T_true, v, w, dt);
			x_est = A[k-1]*x_est + B[k-1]*u;

			if (vis) {
				float pos[3] = {(float)T_true(0,3), (float)T_true(1,3), (float)T_true(2,3)};

				if ((k)%skip == 0) {
					int obj;
					//CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
					//CAL_SetObjectScaling(obj, 0.075f, 0.075f, 1.0f);
					CAL_CreateSphere(cal_samples[s], 0.05f, pos[0], pos[1], pos[2], &obj);
				}
			}
			int col = 1;
			CAL_CheckLineCollision(cal_environment, (float)xprev, (float)yprev, (float)zprev, (float)T_true(0,3), (float)T_true(1,3), (float)T_true(2,3), false, &col);
			if (col != 0) {
				++fail;
				//if (vis) {
				//	CAL_DestroyGroup(cal_samples[s]);
				//}
				break;
			}

			// Measurement update
			n = sampleGaussian(zeros<ZDIM,1>(), Q);

			// dependent of measurement model !!!!
			z(0,0) = T_true(0,3) - path[k].T(0,3);
			z(1,0) = T_true(1,3) - path[k].T(1,3);
			z(2,0) = T_true(2,3) - path[k].T(2,3);

			z += n;
			x_est = x_est + K[k] * (z - H[k]*x_est);
		}
		if (vis) {
			float pos[3] = {(float)T_true(0,3), (float)T_true(1,3), (float)T_true(2,3)};
			int obj;
			//CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
			//CAL_SetObjectScaling(obj, 0.075f, 0.075f, 1.0f);
			CAL_CreateSphere(cal_samples[s], 0.05f, pos[0], pos[1], pos[2], &obj);
		}
	}
	t.stop();
	//std::cout << "Monte-Carlo 1 run time (ms): " << (double)t.interval_mS()/(double)num_samples << std::endl;

	if (vis) {
		delete[] cal_samples;
	}

	//std::cout << "Num failures: " << fail << std::endl;
	//std::cout << "Num samples: " << num_samples << std::endl;

	return (double)fail/(double)num_samples;
}

void computeQuality(const std::vector<PathNode>& path, bool vis, bool lqgmp, bool boole, double & plqgmp, double & pboole) 
{
	int l = (int) path.size();

	preprocess(path);

	R.resize(l);

	R[0].reset();
	R[0].insert(0,0,P0);

	for (int k = 1; k < l; ++k) {

		Matrix<12,12> Y;
		Matrix<12,6+ZDIM> V;

		Y.insert(0,0, A[k-1]);           
		Y.insert(0,6, B[k-1]*L[k-1]);
		Y.insert(6,0, K[k]*H[k]*A[k-1]);
		Y.insert(6,6, A[k-1] + B[k-1]*L[k-1] - K[k]*H[k]*A[k-1]);

		V.insert(0,0, identity<6>()); 
		V.insert(0,6, zeros<6,3>());
		V.insert(6,0, K[k]*H[k]);     
		V.insert(6,6, K[k]);

		Matrix<6+ZDIM,6+ZDIM> Z = zeros<6+ZDIM,6+ZDIM>();
		Z.insert(0,0, N[k-1]);
		Z.insert(6,6, Q);

		R[k] = Y*R[k-1]*~Y + V*Z*~V;
	}

	if (lqgmp) {
		double quality = 1;
		for (int k = 0; k < l; ++k) {
			double conf = computeConfidence(path[k].T.subMatrix<3,1>(0,3), path[k].T.subMatrix<3,3>(0,0) * R[k].subMatrix<3,3>(0,0) * ~path[k].T.subMatrix<3,3>(0,0));
			double prob = incompletegamma(1.5, 0.5*conf*conf);
			quality *= prob;
		}
		plqgmp = (1 - quality);
	
		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);
	}

	if (vis) {
		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);

		drawPath3d(path, cal_paths);

		for (int i = 0; i < l; i+=skip) {
			drawEllipse3d(path[i].T.subMatrix<3,1>(0,3), path[i].T.subMatrix<3,3>(0,0) * R[i].subMatrix<3,3>(0,0) * ~path[i].T.subMatrix<3,3>(0,0), cal_ellipse);
		}
		drawEllipse3d(path[l-1].T.subMatrix<3,1>(0,3), path[l-1].T.subMatrix<3,3>(0,0) * R[l-1].subMatrix<3,3>(0,0) * ~path[l-1].T.subMatrix<3,3>(0,0), cal_ellipse);
	}

	//int node = 32;
	//CAL_CreateSphere(cal_box, 0.1f, (float)path[node].T(0,3), (float)path[node].T(1,3), (float)path[node].T(2,3));
	//Matrix<3> xk = path[node].T.subMatrix<3,1>(0,3);
	//Matrix<3,3> Sk = path[node].T.subMatrix<3,3>(0,0) * R[node].subMatrix<3,3>(0,0) * ~path[node].T.subMatrix<3,3>(0,0);
	//drawEllipse3d(xk, Sk, cal_ellipse);
	//std::vector<std::pair<Matrix<3>, double> > cvx;
	//query(xk, Sk, cvx);
	//std::cout << "p(boole): " << probSuccessBoole(xk, Sk, cvx) << " p(trunc): " << probSuccess(xk, Sk, cvx) << std::endl;

	if (boole) 
	{
		double pb = 1.0, pt = 1.0;

		// for each time step
		for(int k = 0; k < l; ++k) 
		{
			Matrix<3> xk = path[k].T.subMatrix<3,1>(0,3);
			Matrix<3,3> Sk = path[k].T.subMatrix<3,3>(0,0) * R[k].subMatrix<3,3>(0,0) * ~path[k].T.subMatrix<3,3>(0,0);

			std::vector<std::pair<Matrix<3>, double> > cvx;
			query(xk, Sk, cvx);

			pb *= probSuccessBoole(xk, Sk, cvx);
			//pt *= probSuccess(xk, Sk, cvx);
			//std::cout << "pb: " << pb << " pt: " << pt << std::endl;
		}
		pboole = (1 - pb);
	}
	
	//std::cout << "Quality: " << quality << std::endl;
	//std::cout << "Prob success Boole: " << pb << std::endl;
	//std::cout << "Prob success Trunc: " << pt << std::endl;
}

double computeLQGMPTruncate(const std::vector<PathNode>& path, bool vis = false) 
{
	Timer t;
	t.start();

	int l = (int) path.size();
	//std::cout << "Length of path: " << l << std::endl;

	preprocess(path);

	y.resize(l);
	R.resize(l);
	
	y[0].reset();
	R[0].reset();
	R[0].insert(0,0,P0);

	Matrix<12> yTrunc;
	Matrix<12,12> RTrunc;

	double ps = 1.0;

	// Kalman
	Matrix<6,6> P = P0;
	
	for (int k = 1; k < l; ++k) 
	{
		//std::cout << "Stage: " << k << " ";
		
		Matrix<12,12> Y;
		Matrix<12,6+ZDIM> V;

		// Truncate
		std::vector<std::pair<Matrix<3>, double> > cvx;
		Matrix<3> phat = path[k-1].T.subMatrix<3,1>(0,3);
		Matrix<3,3> Rhat = path[k-1].T.subMatrix<3,3>(0,0);
		query(phat + Rhat*y[k-1].subMatrix<3,1>(0,0), Rhat * R[k-1].subMatrix<3,3>(0,0) * ~Rhat, cvx);
		
		ps *= lqgmpTruncation(y[k-1], R[k-1], phat, Rhat, yTrunc, RTrunc, cvx);
		
		// TODO: Jacobians (constant), recompute other matrices??

		P = A[k-1]*P*~A[k-1] + N[k-1];

		H[k].reset();

		// dependent on measurement model !!!!!!
		Matrix<3,3> Rbar = rotFromErr(y[k-1].subMatrix<3,1>(3,0));
		Matrix<3,3> RM = Rhat*Rbar;
		
		H[k].insert(0,0, path[k].T.subMatrix<1,3>(0,0));
		H[k].insert(1,0, path[k].T.subMatrix<1,3>(1,0));
		H[k].insert(2,0, path[k].T.subMatrix<1,3>(2,0));

		//H[k].insert(0,0, RM.subMatrix<1,3>(0,0));
		//H[k].insert(1,0, RM.subMatrix<1,3>(1,0));
		//H[k].insert(2,0, RM.subMatrix<1,3>(2,0));

		K[k] = P*~H[k]*!(H[k]*P*~H[k] + Q);
		P = (identity<6>() - K[k]*H[k])*P;

		Y.insert(0,0, A[k-1]);           
		Y.insert(0,6, B[k-1]*L[k-1]);
		Y.insert(6,0, K[k]*H[k]*A[k-1]);
		Y.insert(6,6, A[k-1] + B[k-1]*L[k-1] - K[k]*H[k]*A[k-1]);

		V.insert(0,0, identity<6>()); 
		V.insert(0,6, zeros<6,3>());
		V.insert(6,0, K[k]*H[k]);     
		V.insert(6,6, K[k]);
		
		Matrix<6+ZDIM,6+ZDIM> Z = zeros<6+ZDIM,6+ZDIM>();
		Z.insert(0,0, N[k-1]);
		Z.insert(6,6, Q);

		//std::cout << "yTrunc: " << ~yTrunc << std::endl;
		//std::cout << "RTrunc: " << RTrunc << std::endl;

		y[k] = Y*yTrunc;
		R[k] = Y*RTrunc*~Y + V*Z*~V;

		//std::cout << "ps: " << ps << std::endl;
	}

	//std::cout << "LQG truncation prob success: " << ps << std::endl;
	t.stop();
	std::cout << std::setprecision(12) << "LQG truncation time (ms): " << t.interval_mS() << std::endl;

	if (vis) 
	{
		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);

		//drawPath3d(path, cal_paths);

		CAL_EmptyGroup(cal_paths_trunc);
		
		float* pt = new float[3*l];
		for (int i = 0; i < l; ++i) {
			Matrix<3> phat = path[i].T.subMatrix<3,1>(0,3);
			Matrix<3,3> Rhat = path[i].T.subMatrix<3,3>(0,0);
			Matrix<3> pbar = y[i].subMatrix<3,1>(0,0);
			Matrix<3> p = phat + Rhat*pbar;
			pt[3*i] = (float)p[0]; pt[3*i+1] = (float)p[1]; pt[3*i+2] = p[2];
		}
		drawPolyline3d(l, pt, cal_paths_trunc);
		
		for (int i = 0; i < l; i+=skip) {
			Matrix<3,3> Rhat = path[i].T.subMatrix<3,3>(0,0);
			Matrix<3> p;
			p[0] = pt[3*i]; p[1] = pt[3*i+1]; p[2] = pt[3*i+2];
			drawEllipse3d(p, Rhat * R[i].subMatrix<3,3>(0,0) * ~Rhat, cal_ellipse_trunc, true);
		}
		Matrix<3> pf;
		pf[0] = pt[3*(l-1)]; pf[1] = pt[3*(l-1)+1]; pf[2] = pt[3*(l-1)+2];
		drawEllipse3d(pf, path[l-1].T.subMatrix<3,3>(0,0) * R[l-1].subMatrix<3,3>(0,0) * ~path[l-1].T.subMatrix<3,3>(0,0), cal_ellipse_trunc, true);

		delete[] pt;
	}

	return (1.0 - ps);
}

/*
double computePSSampling(const std::vector<PathNode>& path, int num_samples, bool vis = false) 
{
	int l = (int) path.size();
	//std::cout << "Length of path: " << l << std::endl;

	preprocess(path);

	// Simulate
	Matrix<6> x_est;
	Matrix<6> x_true;
	Matrix<3> u;
	Matrix<3> m;
	Matrix<6> x_true_old;
	Matrix<ZDIM> n;
	Matrix<ZDIM> z;
	Matrix<6,6> P;

	std::vector<std::vector<Matrix<3> > > pts(l);

	int fail = 0;

	int* cal_samples;
	if (vis) {
		cal_samples = new int[num_samples];
	}

	int progress = 0;
	for (int s = 0; s < num_samples; ++s) 
	{
		if ((100 * s) / num_samples != progress) {
			progress = (100 * s) / num_samples;
			std::cout << progress << " ";
		}

		if (vis) {
			CAL_CreateGroup(&(cal_samples[s]), 0, false);
			CAL_SetGroupColor(cal_samples[s], 0, 0, 1);
		}

		P = P0;
		x_est = path[0].x;
		x_true = path[0].x + sampleGaussian(zeros<6>(), P0);

		pts[0].push_back(x_true.subMatrix<3>(0,0));

		for (int k = 1; k < l; ++k) 
		{
			u = path[k-1].u + L[k-1]*(x_est - path[k-1].x);
			m = sampleGaussian(zeros<3>(), M);

			x_true_old = x_true;

			x_true = f(x_true, u, m);
			//pts[k].push_back(x_true.subMatrix<3>(0,0));

			int col;
			CAL_CheckLineCollision(cal_environment, (float) x_true_old[0], (float) x_true_old[1], 0, (float) x_true[0], (float) x_true[1], 0, false, &col);
			if (col != 0) {
				++fail;
				//if (vis) {
				//	CAL_DestroyGroup(cal_samples[s]);
				//}
				break;
			} else {
				pts[k].push_back(x_true.subMatrix<3>(0,0));
				if (vis) {
					float pos[3] = {(float)x_true[0], (float)x_true[1], 0.0f};

					if (k%skip == 0) {
						int obj;
						CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
						CAL_SetObjectScaling(obj, 0.05f, 0.05f, 1.0f);
					}
				}
			}

			// Jacobian around current estimate
			dfdx(x_est, u, A[k]);
			dfdm(x_est, u, V[k]);

			// Process update
			x_est = f(x_est, u, zeros<3>());
			P = A[k]*P*~A[k] + V[k]*M*~V[k];

			// Measurement update
			n = sampleGaussian(zeros<ZDIM,1>(), N);
			z = H[k]*x_true + W[k]*n;

			K[k] = P*~H[k]*!(H[k]*P*~H[k] + W[k]*N*~W[k]);
			x_est = x_est + K[k] * (z - H[k]*x_est);
			P = (identity<6>() - K[k]*H[k])*P;
		}
	}
	std::cout << std::endl;

	// compute mean and variance
	std::vector<Matrix<3> > mu(l); // a priori means
	std::vector<Matrix<3,3> > Sigma(l); // a priori beliefs

	for(int i = 0; i < l; ++i) {
		mu[i].reset();
		int nss = (int)pts[i].size();
		for(int j = 0; j < nss; ++j) {
			mu[i] += pts[i][j];
		}
		mu[i] /= nss;
	}

	for(int i = 0; i < l; ++i) {
		Sigma[i].reset();
		int nss = (int)pts[i].size();
		for(int j = 0; j < nss; ++j) {
			Sigma[i] += (pts[i][j] - mu[i])*~(pts[i][j] - mu[i]);
		}
		Sigma[i] /= (nss - 1);
	}
	//std::cout << "Num surviving samples: " << pts[l-1].size() << std::endl;

	double ps = 1.0;

	for (int k = 0; k < l; ++k) 
	{
		std::vector<std::pair<Matrix<3>, double> > cvx;
		btree->query(mu[k], Sigma[k], cvx);

		int ncvx = (int)cvx.size();
		for(int i = 0; i < ncvx; ++i) {
			std::pair<Matrix<3>, double>& c = cvx[i];

			Matrix<3>& a = c.first;
			double b = c.second;

			//std::cout << "Constraint: " << a[0] << " " << a[1] << " " << b << std::endl;

			ps *= cdf((b - tr(~a*mu[k]))/(sqrt(tr(~a*Sigma[k]*a))));
		}
		std::cout << "Processing belief: " << k << " num survivors: " << pts[k].size() << " num constraints: " << ncvx << " ps: " << ps << std::endl;
	}

	for(int i = 0; i < l; ++i) {
		pts[i].clear();
	}
	pts.clear();

	if (vis) 
	{
		delete[] cal_samples;

		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);

		//drawPath2d(path, cal_paths);

		int np[1] = {l};
		float* pt = new float[3*np[0]];
		for (int i = 0; i < l; ++i) {
			pt[3*i] = (float) mu[i][0]; pt[3*i+1] = (float) mu[i][1]; pt[3*i+2] = 0;
			//CAL_CreateSphere(cal_paths, 0.1f, (float) mu[i][0], (float) mu[i][1], 0.0f);
		}
		//CAL_CreatePolyline(cal_paths_trunc, 1, np, pt);
		//CAL_CreatePolyline(cal_paths, 1, np, pt);
		delete[] pt;

		for (int i = 0; i < l; i+=skip) {
			//drawEllipse2d(mu[i], Sigma[i], cal_ellipse, true);
		}
	}

	std::cout << "Probability of success (Truncation): " << ps << std::endl;
	return ps;
}

void showPaths(const std::vector<std::vector<PathNode> >& paths) 
{
	int np = (int)paths.size();
	for (int j = 0; j < np; ++j) {
		drawPath3d(paths[j], cal_paths);
	}
}
*/

void testPlan()
{
	
	std::ifstream fptr("needleTest.txt", std::ios::out);
	int len;
	fptr >> len;
	std::vector<PathNode> path(len);
	for (int i = 0; i < len; ++i) {
		PathNode node;
		node.deserialize(fptr);
		path[i] = node;
	}
	
	double plqgmp, pboole, ptrunc, psampling;

	//computeQuality(path, true, false, false, plqgmp, pboole);
	
	//psampling = simulate(path, 100000, false);
	
	//ptrunc = computeLQGMPTruncate(path, false);

	//std::cout << "psampling: " << psampling << " plqgmp: " << plqgmp << " pboole: " << pboole << " ptrunc: " << ptrunc << std::endl;
	//std::cout << "psampling: " << psampling << " ptrunc: " << ptrunc << std::endl;

	int ns[19] = {50,100,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,2750,3000};
	std::ofstream fptr2("mc-needle.txt",std::ios::out);

	Timer t;
	for(int i = 0; i < 19; ++i) {
		t.start();
		for(int j = 0; j < 100; ++j) {
			psampling = simulate(path, ns[i], false);
			fptr2 << std::setprecision(10) << psampling << " ";
		}
		fptr2 << std::endl;
		t.stop();
		std::cout << "Finished processing " << ns[i] << " samples in " << t.interval_S() << " secs" << std::endl;
	}
	fptr2.close();

	std::cout << "Done" << std::endl;

	//computePSSampling(path, 3000, true);
}

void generateStats()
{
	CAL_SuspendVisualisation();
	std::ofstream results("results.txt", std::ios::out);

	std::ifstream fptr(SAVEFILE, std::ios::in);
	int len;

	clock_t tstart, tstop;
	clock_t ttruncate = 0, tsampling = 0, tlqgmp = 0, tboole = 0;

	double pctruncate = 0.0, pcsampling = 0.0, pclqgmp = 0.0, pcboole = 0.0;

	for(int p = 0; p < 200; ++p) 
	{
		fptr >> len;
		std::vector<PathNode> path(len);
		for (int i = 0; i < len; ++i) {
			PathNode node;
			node.deserialize(fptr);
			path[i] = node;
		}

		tstart = clock();
		pcsampling = simulate(path, 2500, false);
		tstop = clock();
		tsampling += (tstop - tstart);
		results << pcsampling << " ";

		tstart = clock();
		pctruncate = computeLQGMPTruncate(path, false);
		tstop = clock();
		ttruncate += (tstop - tstart);
		results << pctruncate << " ";

		tstart = clock();
		computeQuality(path, false, true, false, pclqgmp, pcboole);
		tstop = clock();
		tlqgmp += (tstop - tstart);
		results << pclqgmp << " ";

		tstart = clock();
		computeQuality(path, false, false, true, pclqgmp, pcboole);
		tstop = clock();
		tboole += (tstop - tstart);
		results << pcboole << std::endl;

		std::cout << p << " ";
		path.clear();
	}
	//----------------------------------------------------
	
	results.close();

	std::cout << std::endl;
	std::cout << "Clocks per sec: " << (double) (CLOCKS_PER_SEC) << std::endl;
	std::cout << "Truncation time: " << ttruncate << std::endl;
	std::cout << "Sampling time: " << tsampling << std::endl;
	std::cout << "lqgmp time: " << tlqgmp << std::endl;
	std::cout << "Boole time: " << tboole << std::endl;

	
	CAL_ResumeVisualisation();
}

int main()
{
	time_t t;
	srand(unsigned(time(&t)));

	// Create Obstacles in Callisto
	CAL_Initialisation (true, true, true);

	initEnvironment();

	//----------------------------------------------------
	// Construct Paths
	/*
	clock_t startTime = clock();
	int numPaths = 10;
	std::cout << numPaths << std::endl;
	for (int i = 0; i < numPaths; ++i) {
	rrt();
	std::cerr << i << " ";
	}
	CAL_End();
	clock_t endTime = clock();
	std::cerr << std::endl << (double) (endTime - startTime) / CLOCKS_PER_SEC << std::endl;
	showPaths(paths);
	return 0;
	*/
	
	//rrt();
	//int num;
	//std::cin >> num;

	std::ofstream frrt(SAVEFILE, std::ios::out);

	clock_t tstart = clock();
	for(int i = 0; i < 200; ++i) {
		rrt(frrt);
		std::cout << i << " ";
	}

	frrt.close();
	std::cout << "\nRRT generation time: " << (clock() - tstart) / (double) (CLOCKS_PER_SEC) << std::endl;

	int num;
	std::cin >> num;


	//----------------------------------------------------
	
	testPlan();
	//generateStats();

	//----------------------------------------------------

	int k;
	std::cin >> k;

	CAL_End();

	return 0;
}