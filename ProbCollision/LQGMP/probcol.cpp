#include "common.h"

double probSuccessBoole(const Matrix<SDIM>& x, const Matrix<SDIM, SDIM>& S, std::vector<std::pair<Matrix<SDIM>, double> >& cvx) 
{
	double pb = 0.0;
	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; ++i) {
		Matrix<SDIM>& a = cvx[i].first;
		double b = cvx[i].second;
		pb += (1.0 - cdf((b - tr(~a*x))/(sqrt(tr(~a*S*a)))));
	}
	return (1.0 - pb);
} 

double probSuccess(const Matrix<SDIM>& x, const Matrix<SDIM, SDIM>& S, std::vector<std::pair<Matrix<SDIM>, double> >& cvx) 
{
	double ps = 1.0;
	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; ++i) {
		Matrix<SDIM>& a = cvx[i].first;
		double b = cvx[i].second;
		ps *= cdf((b - tr(~a*x))/(sqrt(tr(~a*S*a))));
	}
	return ps;
}

double lqgmpTruncation(const Matrix<2*XDIM>& x, const Matrix<2*XDIM, 2*XDIM>& S, Matrix<2*XDIM>& xDelta, 
					   Matrix<2*XDIM, 2*XDIM>& SDelta, std::vector<std::pair<Matrix<SDIM>, double> >& cvx)
{
	double yMean, yVar, yNewMean, yNewVar;
	Matrix<2*XDIM> xyCovar, L;
	xDelta.reset();
	SDelta.reset();

	double ps = 1.0;

	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; ++i) {
		std::pair<Matrix<SDIM>, double>& c = cvx[i];

		Matrix<2*XDIM> A;
		A.reset();
		A.insert<SDIM,1>(0,0,c.first);
		
		double B = c.second;
		
		yMean = tr(~A*x);
		yVar = tr(~A*S*A);
	
		truncate(B, yMean, yVar, yNewMean, yNewVar);
		
		xyCovar = S*A;
		L = (xyCovar/yVar);

		double ratio1 = (yNewVar/yVar);
		double ratio2 = cdf((B - tr(~A*x))/(sqrt(tr(~A*S*A))));
		
		xDelta += L*(yNewMean - yMean);
		
		//SDelta += (ratio1 - 1.0)*(L*~xyCovar);
		//SDelta += (ratio2 - 1.0)*(L*~xyCovar);
		
		SDelta += ((ratio1 - 1.0)*yVar)*(L*~L);
		//SDelta += ((ratio2 - 1.0)*yVar)*(L*~L);
		//SDelta += (((ratio1 + ratio2)*0.5 - 1.0)*yVar)*(L*~L);
		
		//std::cout << "yNewMean/yMean: " << yNewMean/yMean << " yNewVar/yVar: " << yNewVar/yVar << " prob mass: " << cdf((B - tr(~A*x))/(sqrt(tr(~A*S*A)))) << std::endl;

		ps *= ratio2;
	}

	checkNaN(xDelta);
	checkNaN(SDelta);

	return ps;
}

void preprocess(const std::vector<PathNode>& path) 
{
	int l = (int) path.size();

	A.clear(); B.clear(); V.clear(); H.clear(); W.clear(); L.clear(); K.clear(); y.clear(); R.clear();

	A.resize(l); // process matrix
	B.resize(l); // input matrix
	V.resize(l); // process noise matrix
	H.resize(l); // measurement matrix
	W.resize(l); // measurement noise matrix
	L.resize(l); // feedback matrix
	K.resize(l); // Kalman-gain matrix

	initMatrices(C, D, P0, M, N);

	// Jacobians
	for (int k = 1; k < l; ++k) {
		const Matrix<XDIM>& x = path[k-1].x;
		const Matrix<UDIM>& u = path[k-1].u;

		dfdx(x, u, A[k]);
		dfdu(x, u, B[k]);
		dfdm(x, u, V[k]);
		
		dhdx(x, u, H[k]);
		dhdn(x, u, W[k]);
	}

	// LQR
	Matrix<XDIM, XDIM> S;
	S = C;
	L[l - 1] = zeros<UDIM, XDIM>();
	for (int k = l - 2; k >= 0; --k) {
		L[k] = -!(~B[k+1]*S*B[k+1] + D)*~B[k+1]*S*A[k+1];
		S = C + ~A[k+1]*S*A[k+1] + ~A[k+1]*S*B[k+1]*L[k];
	}

	// Kalman
	Matrix<XDIM, XDIM> P;
	P = P0;
	K[0] = zeros<XDIM, ZDIM>();
	for (int k = 1; k < l; ++k) {
		P = A[k]*P*~A[k] + V[k]*M*~V[k];
		K[k] = P*~H[k]*!(H[k]*P*~H[k] + W[k]*N*~W[k]);
		P = (identity<XDIM>() - K[k]*H[k])*P;
	}
}

double simulate(const std::vector<PathNode>& path, int num_samples, bool vis = false) 
{
	int l = (int) path.size();

	preprocess(path);

	// Simulate
	Matrix<XDIM> x_est;
	Matrix<XDIM> x_true;
	Matrix<UDIM> u;
	Matrix<UDIM> m;
	Matrix<XDIM> x_true_old;
	Matrix<ZDIM> n;
	Matrix<ZDIM> z;
	Matrix<XDIM, XDIM> P;

	int fail = 0;

	int* cal_samples;
	if (vis) {
		cal_samples = new int[num_samples];
	}

	for (int s = 0; s < num_samples; ++s) {
		if (vis) {
			CAL_CreateGroup(&(cal_samples[s]), 0, false);
			CAL_SetGroupColor(cal_samples[s], 1, 0, 0);
		}

		P = P0;
		x_est = path[0].x;
		x_true = path[0].x + sampleGaussian(zeros<XDIM,1>(), P0);

		for (int k = 1; k < l; ++k) 
		{
			u = path[k-1].u + L[k-1]*(x_est - path[k-1].x);
			m = sampleGaussian(zeros<UDIM,1>(), M);

			x_true_old = x_true;

			x_true = f(x_true, u, m);

			if (vis) {
				float pos[3] = {(float)x_true_old[0], (float)x_true_old[1], 0.01f};
				
				if ((k-1)%skip == 0) {
					int obj;
					//CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
					//CAL_SetObjectScaling(obj, 0.075f, 0.075f, 1.0f);
					CAL_CreateSphere(cal_samples[s], 0.05f, pos[0], pos[1], pos[2], &obj);
					
				}
			}
			int col;
			CAL_CheckLineCollision(cal_environment, (float) x_true_old[0], (float) x_true_old[1], 0, (float) x_true[0], (float) x_true[1], 0, false, &col);
			if (col != 0) {
				++fail;
				if (vis) {
					CAL_DestroyGroup(cal_samples[s]);
				}
				break;
			}
			
			// Jacobian around current estimate
			dfdx(x_est, u, A[k]);
			dfdm(x_est, u, V[k]);

			// Process update
			x_est = f(x_est, u, zeros<UDIM,1>());
			P = A[k]*P*~A[k] + V[k]*M*~V[k];

			// Measurement update
			n = sampleGaussian(zeros<ZDIM,1>(), N);
			z = H[k]*x_true + W[k]*n;

			K[k] = P*~H[k]*!(H[k]*P*~H[k] + W[k]*N*~W[k]);
			x_est = x_est + K[k] * (z - H[k]*x_est);
			P = (identity<XDIM>() - K[k]*H[k])*P;
		}
		if (vis) {
			float pos[3] = {(float)x_true[0], (float)x_true[1], 0};
			int obj;
			//CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
			//CAL_SetObjectScaling(obj, 0.075f, 0.075f, 1.0f);
			CAL_CreateSphere(cal_samples[s], 0.05f, pos[0], pos[1], pos[2], &obj);
		}
	}

	// Divergence
	/*
	for (int k = 0; k < l; ++k) {
		var[k] /= (num_samples - 1);
		Matrix<XDIM, XDIM> var0 = R[k].subMatrix<XDIM,XDIM>(0,0);

		double div = 0.5 * (tr(!var[k]*var0) - log(det(var0)/det(var[k])) - XDIM) + 0.5 * (tr(!var0*var[k]) - log(det(var[k])/det(var0)) - XDIM);
		std::cout << k << " " << div << std::endl;
	}
	*/

	if (vis) {
		delete[] cal_samples;
	}

	std::cout << "Num failures: " << fail << std::endl;
	std::cout << "Num samples: " << num_samples << std::endl;

	return (1 - (double)fail/(double)num_samples);
}

void computeQuality(const std::vector<PathNode>& path, bool vis, bool lqgmp, bool boole, double & pclqgmp, double& pcboole) 
{
	int l = (int) path.size();

	preprocess(path);

	std::vector<Matrix<2*XDIM, 2*XDIM> > F(l);
	std::vector<Matrix<2*XDIM, UDIM+ZDIM> > G(l);
	
	R.resize(l);

	// Combination of LQR and Kalman
	Matrix<UDIM+ZDIM, UDIM+ZDIM> Q = zeros<UDIM+ZDIM, UDIM+ZDIM>();
	Q.insert(0,0, M); Q.insert(UDIM, UDIM, N);

	R[0].reset();
	R[0].insert(0,0,P0);

	for (int k = 1; k < l; ++k) {
		F[k].insert(0,0,     A[k]);           F[k].insert(0,XDIM,     B[k]*L[k-1]);
		F[k].insert(XDIM,0, K[k]*H[k]*A[k]); F[k].insert(XDIM,XDIM, A[k] + B[k]*L[k-1] - K[k]*H[k]*A[k]);

		G[k].insert(0,0,     V[k]);           G[k].insert(0,UDIM,     zeros<XDIM,ZDIM>());
		G[k].insert(XDIM,0, K[k]*H[k]*V[k]); G[k].insert(XDIM,UDIM, K[k]*W[k]);

		R[k] = F[k]*R[k-1]*~F[k] + G[k]*Q*~G[k];
	}

	if (lqgmp) {
		double quality = 1;
		for (int k = 0; k < l; ++k) {
			double conf = computeConfidence(path[k].x.subMatrix<SDIM,1>(0,0), R[k].subMatrix<SDIM,SDIM>(0,0));
			double prob = incompletegamma(0.5*SDIM, 0.5*conf*conf);
			quality *= prob;
		}
		pclqgmp = 1 - quality;
		//std::cout << "pclqgmp: " << pclqgmp << std::endl;
		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);
	}

	if (vis) {
		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);

		//drawPath2d(path, cal_paths);

		for (int i = 0; i < l; i+=skip) {
			//drawEllipse2d(path[i].x.subMatrix<SDIM,1>(0,0), R[i].subMatrix<SDIM,SDIM>(0,0), cal_ellipse, false);
			drawEllipse2d(path[i].x.subMatrix<SDIM,1>(0,0), R[i].subMatrix<SDIM,SDIM>(0,0), cal_ellipse, false);
		}
		drawEllipse2d(path[l-1].x.subMatrix<SDIM,1>(0,0), R[l-1].subMatrix<SDIM,SDIM>(0,0), cal_ellipse, false);
	}

	if (boole) {
		double pb = 1.0, pt = 1.0;

		// for each time step
		for(int k = 0; k < l; ++k) 
		{
			Matrix<XDIM> xk = path[k].x;
			Matrix<XDIM, XDIM> Sk = R[k].subMatrix<XDIM,XDIM>(0,0);

			std::vector<std::pair<Matrix<SDIM>, double> > cvx;
			btree->query(xk.subMatrix<SDIM,1>(0,0), Sk.subMatrix<SDIM,SDIM>(0,0), cvx);

			pb *= probSuccessBoole(xk.subMatrix<SDIM,1>(0,0), Sk.subMatrix<SDIM,SDIM>(0,0), cvx);

			//Matrix<XDIM> xkTrunc;
			//Matrix<XDIM, XDIM> SkTrunc;
			//pt *= probSuccessTrunc(xk, Sk, xkTrunc, SkTrunc, cvx);
		}
		pcboole = (1.0-pb);
	}

	//std::cout << "Quality: " << quality << std::endl;
	//std::cout << "Prob success Boole: " << pb << std::endl;
	//std::cout << "Prob success Trunc: " << pt << std::endl;
}

void query(const Matrix<SDIM>& x, const Matrix<SDIM, SDIM>& S, std::vector<std::pair<Matrix<SDIM>, double> >& cvx)
{
	Matrix<SDIM,SDIM> V, E;
	jacobi(S, V, E);

	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~V);
	Matrix<4,1> q = quatFromRot(Temp);

	double sx = 3.0*sqrt(E(0,0));
	double sy = 3.0*sqrt(E(1,1));
	Matrix<SDIM,SDIM> invScale = zeros<SDIM,SDIM>();
	invScale(0,0) = 1/sx;
	invScale(1,1) = 1/sy;
	Matrix<SDIM,SDIM> T = invScale * ~V; // ? == !(V * !invscale)
	Matrix<SDIM> transPos =  T * x;
	Matrix<SDIM,SDIM> invT = !T;

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, (float) invScale(0,0), (float) invScale(1,1), 1.0f);

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0.0f);

	int num_pairs;
	int retcode = CAL_GetClosestPairs (cal_point, cal_environment, &num_pairs);

	//std::cout << "Ret code for CAL_GetClosestPairs: " << retcode << std::endl;

	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	for(int it1 = 0; it1 < num_pairs; ++it1) 
	{
		Matrix<SDIM> p0, p1, vec, normal;
		p0[0] = (double)results[it1].vector0[0];
		p0[1] = (double)results[it1].vector0[1];

		p1[0] = (double)results[it1].vector1[0];
		p1[1] = (double)results[it1].vector1[1];

		vec = (p1 - p0);
		double vecnorm = tr(~vec*vec);
		vec = vec/vecnorm;

		normal[0] = p1[0] + vec[1]; normal[1] = p1[1] - vec[0];

		normal = invT*normal; p1 = invT*p1;

		vec = (normal - p1);
		vecnorm = sqrt(tr(~vec*vec));
		vec = vec/vecnorm;

		normal[0] = -vec[1]; normal[1] = vec[0];

		Matrix<SDIM> a = normal;
		double b = tr(~normal*p1);

		cvx.push_back(std::make_pair(a, b));
	}

	delete[] results;

	//std::cout << "Num truncation constraints: " << (int)cvx.size() << std::endl;

	CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
	CAL_SetGroupScaling(cal_environment, 1, 1, 1);
}

double computeLQGMPTruncate(const std::vector<PathNode>& path, bool vis = false) 
{
	Timer t;
	t.start();

	int l = (int) path.size();
	//std::cout << "Length of path: " << l << std::endl;

	preprocess(path);

	std::vector<Matrix<2*XDIM, 2*XDIM> > F(l);
	std::vector<Matrix<2*XDIM, UDIM+ZDIM> > G(l);
	
	y.resize(l);
	R.resize(l);

	// Combination of LQR and Kalman
	Matrix<UDIM+ZDIM, UDIM+ZDIM> Q = zeros<UDIM+ZDIM, UDIM+ZDIM>();
	Q.insert(0,0, M); Q.insert(UDIM, UDIM, N);

	y[0].reset();
	R[0].reset();
	R[0].insert(0,0,P0);

	Matrix<2*XDIM,1> ykm1;
	
	Matrix<2*XDIM> yDelta;
	Matrix<2*XDIM, 2*XDIM> RDelta;
	
	double ps = 1.0;

	for (int k = 1; k < l; ++k) 
	{
		// Truncate
		ykm1.reset();
		ykm1.insert<XDIM,1>(0,0,path[k-1].x);
		ykm1.insert<XDIM,1>(XDIM,0,path[k-1].x);

		std::vector<std::pair<Matrix<SDIM>, double> > cvx;
		
		btree->query((ykm1 + y[k-1]).subMatrix<SDIM,1>(0,0), R[k-1].subMatrix<SDIM, SDIM>(0,0), cvx);
		//query((ykm1 + y[k-1]).subMatrix<SDIM,1>(0,0), R[k-1].subMatrix<SDIM, SDIM>(0,0), cvx);

		ps *= lqgmpTruncation((ykm1 + y[k-1]), R[k-1], yDelta, RDelta, cvx);

		// Recompute Jacobians
		Matrix<XDIM> x = (ykm1 + yDelta).subMatrix<XDIM,1>(0,0);
		Matrix<UDIM> u = path[k-1].u;

		dfdx(x, u, A[k]);
		dfdu(x, u, B[k]);
		dfdm(x, u, V[k]);
		
		dhdx(x, u, H[k]);
		dhdn(x, u, W[k]);
		
		F[k].insert(0,0,     A[k]);           F[k].insert(0,XDIM,     B[k]*L[k-1]);
		F[k].insert(XDIM,0, K[k]*H[k]*A[k]); F[k].insert(XDIM,XDIM, A[k] + B[k]*L[k-1] - K[k]*H[k]*A[k]);

		G[k].insert(0,0,     V[k]);           G[k].insert(0,UDIM, zeros<XDIM,ZDIM>());
		G[k].insert(XDIM,0, K[k]*H[k]*V[k]); G[k].insert(XDIM,UDIM, K[k]*W[k]);
		
		//y[k] = F[k]*y[k-1] + yDelta;
		//R[k] = F[k]*R[k-1]*~F[k] + G[k]*Q*~G[k] + RDelta;
		
		y[k] = F[k]*(y[k-1] + yDelta);
		R[k] = F[k]*(R[k-1] + RDelta)*~F[k] + G[k]*Q*~G[k];
	}
	t.stop();
	std::cout << std::setprecision(12) << "LQG truncation time: " << t.interval_mS() << std::endl;

	//std::cout << "LQG truncation prob success: " << ps << std::endl;
	
	if (vis) 
	{
		CAL_SetGroupQuaternion(cal_obstacles, 0, 0, 0, 1);
		CAL_SetGroupScaling(cal_environment, 1, 1, 1);

		//drawPath2d(path, cal_paths);
		
		//CAL_EmptyGroup(cal_paths_trunc);
		//int np[1] = {l};
		//float* pt = new float[3*np[0]];
		//for (int i = 0; i < l; ++i) {
		//	pt[3*i] = (float)(path[i].x[0] + y[i][0]); pt[3*i+1] = (float)(path[i].x[1] + y[i][1]); pt[3*i+2] = 0;
		//}
		//CAL_CreatePolyline(cal_paths_trunc, 1, np, pt);
		//delete[] pt;

		drawPath2d(path, cal_paths);

		Vector2d p1, p2, mid, dv;
		double len;
		p1.x = path[0].x[0] + y[0][0]; p1.y = path[0].x[1] + y[0][1];

		for(int i = 1; i < l; ++i) 
		{
			p2.x = path[i].x[0] + y[i][0]; p2.y = path[i].x[1] + y[i][1];
			dv = Vector2d(p2.x - p1.x, p2.y - p1.y);
			len = norm(dv);
			dv = dv/len;
			mid = Vector2d(0.5*(p1.x + p2.x), 0.5*(p1.y + p2.y));

			int obj;
			CAL_CreateBox(cal_paths_trunc, (float)(len + 0.1), 0.1f, 0.01f, (float)mid.x, (float)mid.y, 0.02f, &obj);
			CAL_SetObjectOrientation(obj, 0.0f, 0.0f, (float)(mod2pi(atan2(dv.y, dv.x))));

			p1 = p2;
		}

		for (int i = 0; i < l; i+=skip) {
			//drawEllipse2d(path[i].x.subMatrix<SDIM,1>(0,0) + y[i].subMatrix<SDIM,1>(0,0), R[i].subMatrix<SDIM,SDIM>(0,0), cal_ellipse_trunc, false);
			drawEllipse2d(path[i].x.subMatrix<SDIM,1>(0,0) + y[i].subMatrix<SDIM,1>(0,0), R[i].subMatrix<SDIM,SDIM>(0,0), cal_ellipse_trunc, true);
		}
		//drawEllipse2d(path[l-1].x.subMatrix<SDIM,1>(0,0) + y[l-1].subMatrix<SDIM,1>(0,0), R[l-1].subMatrix<SDIM,SDIM>(0,0), cal_ellipse_trunc, false);
		drawEllipse2d(path[l-1].x.subMatrix<SDIM,1>(0,0) + y[l-1].subMatrix<SDIM,1>(0,0), R[l-1].subMatrix<SDIM,SDIM>(0,0), cal_ellipse_trunc, true);
	}

	return (1-ps);
}

double computePSSampling(const std::vector<PathNode>& path, int num_samples, bool vis = false) 
{
	int l = (int) path.size();
	//std::cout << "Length of path: " << l << std::endl;

	preprocess(path);

	// Simulate
	Matrix<XDIM> x_est;
	Matrix<XDIM> x_true;
	Matrix<UDIM> u;
	Matrix<UDIM> m;
	Matrix<XDIM> x_true_old;
	Matrix<ZDIM> n;
	Matrix<ZDIM> z;
	Matrix<XDIM, XDIM> P;

	std::vector<std::vector<Matrix<SDIM> > > pts(l);
	
	int fail = 0;

	int* cal_samples;
	if (vis) {
		cal_samples = new int[num_samples];
	}

	int progress = 0;
	Timer t;
	t.start();
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
		x_est = path[0].x;
		x_true = path[0].x + sampleGaussian(zeros<XDIM,1>(), P0);

		pts[0].push_back(x_true.subMatrix<SDIM,1>(0,0));

		if (vis) {
			float pos[3] = {(float)x_true[0], (float)x_true[1], 0.1f};
			int obj;
			CAL_CreateSphere(cal_samples[s], 0.035f, pos[0], pos[1], pos[2], &obj);
		}

		for (int k = 1; k < l; ++k) 
		{
			u = path[k-1].u + L[k-1]*(x_est - path[k-1].x);
			m = sampleGaussian(zeros<UDIM,1>(), M);

			x_true_old = x_true;

			x_true = f(x_true, u, m);
			//pts[k].push_back(x_true.subMatrix<SDIM,1>(0,0));
			
			int col;
			CAL_CheckLineCollision(cal_environment, (float) x_true_old[0], (float) x_true_old[1], 0, (float) x_true[0], (float) x_true[1], 0, false, &col);
			if (col != 0) {
				++fail;
				//if (vis) {
				//	CAL_DestroyGroup(cal_samples[s]);
				//}
				break;
			} else {
				pts[k].push_back(x_true.subMatrix<SDIM,1>(0,0));
				if (vis) {
					float pos[3] = {(float)x_true[0], (float)x_true[1], 0.1f};
		
					if (k == l-1 || k%skip == 0) {
						int obj;
						CAL_CreateSphere(cal_samples[s], 0.035f, pos[0], pos[1], pos[2], &obj);
						//CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
						//CAL_SetObjectScaling(obj, 0.05f, 0.05f, 1.0f);
					}
				}
			}
			
			/*
			pts[k].push_back(x_true.subMatrix<SDIM,1>(0,0));
			if (vis) {
				float pos[3] = {(float)x_true[0], (float)x_true[1], 0.0f};

				if (k%skip == 0) {
					int obj;
					CAL_CreateSphere(cal_samples[s], 0.05f, pos[0], pos[1], pos[2], &obj);
					//CAL_CreateUserDrawn(cal_samples[s], drawUnitCircle, NULL, pos[0], pos[1], pos[2], &obj);
					//CAL_SetObjectScaling(obj, 0.05f, 0.05f, 1.0f);
				}
			}
			*/

			// Jacobian around current estimate
			dfdx(x_est, u, A[k]);
			dfdm(x_est, u, V[k]);

			// Process update
			x_est = f(x_est, u, zeros<UDIM,1>());
			P = A[k]*P*~A[k] + V[k]*M*~V[k];

			// Measurement update
			n = sampleGaussian(zeros<ZDIM,1>(), N);
			z = H[k]*x_true + W[k]*n;

			K[k] = P*~H[k]*!(H[k]*P*~H[k] + W[k]*N*~W[k]);
			x_est = x_est + K[k] * (z - H[k]*x_est);
			P = (identity<XDIM>() - K[k]*H[k])*P;
		}
	}
	//std::cout << std::endl;
	t.stop();
	//std::cout << std::setprecision(12) << "MC sampling (one iter) time: " << (double)t.interval_mS()/(double)num_samples << std::endl;

	// compute mean and variance
	std::vector<Matrix<SDIM> > mu(l); // a priori means
	std::vector<Matrix<SDIM, SDIM> > Sigma(l); // a priori beliefs

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

	return 1.0 - ((double)pts[l-1].size()/(double)num_samples);

	double ps = 1.0;

	for (int k = 0; k < l; ++k) 
	{
		std::vector<std::pair<Matrix<SDIM>, double> > cvx;
		btree->query(mu[k], Sigma[k], cvx);

		int ncvx = (int)cvx.size();
		for(int i = 0; i < ncvx; ++i) {
			std::pair<Matrix<SDIM>, double>& c = cvx[i];

			Matrix<SDIM>& a = c.first;
			double b = c.second;

			//std::cout << "Constraint: " << a[0] << " " << a[1] << " " << b << std::endl;

			ps *= cdf((b - tr(~a*mu[k]))/(sqrt(tr(~a*Sigma[k]*a))));
		}
		//std::cout << "Processing belief: " << k << " num survivors: " << pts[k].size() << " num constraints: " << ncvx << " ps: " << ps << std::endl;
	}
	
	for(int i = 0; i < l; ++i) {
		pts[i].clear();
	}
	pts.clear();
	
	if (vis) 
	{
		delete[] cal_samples;
	
		/*
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
		*/
	}

	//std::cout << "Probability of success (Truncation): " << ps << std::endl;
	return ps;
}

void showPaths(const std::vector<std::vector<PathNode> >& paths) 
{
	int np = (int)paths.size();
	for (int j = 0; j < np; ++j) {
		drawPath2d(paths[j], cal_paths);
	}
}

void rrt(std::ofstream & fptr) 
{
	std::vector<RRTNode> rrtTree;
	//CAL_EmptyGroup(cal_rrt);

	RRTNode startNode;
	startNode.x = start;

	rrtTree.push_back(startNode);

	int maxIter = 5000;
	int iter = 0;

	while (iter < maxIter && tr(~(rrtTree.back().x.subMatrix<SDIM,1>(0,0) - goal) * (rrtTree.back().x.subMatrix<SDIM,1>(0,0) - goal)) > 0.25* goalRadius * goalRadius) 
	{
		++iter;
		Matrix<XDIM> sample;
		for (int i = 0; i < XDIM; ++i) {
			sample[i] = uniform(x_min[i], x_max[i]);
		}
		int col;

		CAL_CheckPointCollision(cal_environment, (float) sample[0], (float) sample[1], 0, false, &col);
		if (col != 0) {
			continue;
		}

		int node = -1;
		double mindist = 9e99;
		for (int j = 0; j < (int) rrtTree.size(); ++j) {
			double ddist = tr(~(rrtTree[j].x - sample) * (rrtTree[j].x - sample));
			if (ddist < mindist) {
				mindist = ddist;
				node = j;
			}
		}
		if (node == -1) {
			continue;
		}

		Matrix<UDIM> input;
		for (int i = 0; i < UDIM; ++i) {
			input[i] = uniform() * (u_max[i] - u_min[i]) + u_min[i];
		}

		Matrix<XDIM> x_new = f(rrtTree[node].x, input, zeros<UDIM,1>());

		bool valid = true;
		for (int i = 0; i < XDIM; ++i) {
			if (x_new[i] < x_min[i] || x_new[i] > x_max[i]) {
				valid = false;
				break;
			}
		}
		if (!valid) {
			continue;
		}

		CAL_CheckLineCollision(cal_environment, (float) rrtTree[node].x[0], (float) rrtTree[node].x[1], 0, (float) x_new[0], (float) x_new[1], 0, false, &col);
		if (col != 0) {
			continue;
		}

		RRTNode newnode;
		newnode.x = x_new;
		newnode.u = input;
		newnode.bp = node;

		rrtTree.push_back(newnode); 

		//int np[1] = {2};
		//float p[6] = {(float) rrtTree[node].x[0], (float) rrtTree[node].x[1], 0, (float) x_new[0], (float) x_new[1], 0};
		//CAL_CreatePolyline(cal_rrt, 1, np, p);
	}

	if (iter == maxIter) {
		return;
	}

	//int k;
	//std::cin >> k;

	std::deque<PathNode> path;

	int i = (int) rrtTree.size() - 1;
	PathNode node;
	node.x = rrtTree[i].x;
	node.u = zeros<UDIM, 1>();
	path.push_front(node);

	while (i != 0) {
		node.u = rrtTree[i].u;
		i = rrtTree[i].bp;
		node.x = rrtTree[i].x;

		path.push_front(node);
	}

	/*
	std::cout << path.size() << std::endl;
	for (int i = 0; i < (int) path.size(); ++i) {
		std::cout << ~path[i].x;
		std::cout << ~path[i].u;
	}
	*/

	fptr << (int) path.size();
	for (int i = 0; i < (int) path.size(); ++i) {
		path[i].serialize(fptr);
	}
	fptr << std::endl;

	std::cout << "Finished RRT" << std::endl;
}

/*
void selectBestPath() 
{
	std::vector<std::vector<PathNode> > paths;
	
	// Read Paths
	int numPaths;
	std::cin >> numPaths;
	paths.resize(numPaths);
	for (int i = 0; i < numPaths; ++i) {
		int length;
		std::cin >> length;
		std::vector<PathNode> path(length);
		for (int j = 0; j < length; ++j) {
			PathNode node;
			std::cin >> node.x;
			std::cin >> node.u;
			path[j] = node;
		}
		paths[i] = path;
	}
	
	// Compute Optimal Path
	clock_t startTime = clock();

	std::multimap<double, int> bestPaths;
	int bestPath = 0;
	double bestQuality = -DBL_MAX;

	for (int i = 0; i < (int) paths.size(); ++i) {
		double quality = computeQuality(paths[i]);
		//double quality = simulate(paths[i], 10000);
		//std::cout << i << " " << quality << std::endl;
		std::cerr << i << " ";

		if (bestPaths.size() < 4 || quality > bestPaths.begin()->first) {
			bestPaths.insert(std::make_pair(quality, i));
			if (bestPaths.size() > 4) {
				bestPaths.erase(bestPaths.begin());
			}
		}

		if (quality > bestQuality) {
			bestPath = i;
			bestQuality = quality;

			std::cerr << std::endl << "Best path so far: " << bestPath << " " << bestQuality << std::endl;
		}
	}

	clock_t endTime = clock();

	std::cerr << std::endl << "Running time: " << (double) (endTime - startTime) / CLOCKS_PER_SEC << std::endl;
	//std::cout << (double) (endTime - startTime) / CLOCKS_PER_SEC << std::endl;

	//std::cout << paths[bestPath].size() << std::endl;
	//for (int i = 0; i < (int) paths[bestPath].size(); ++i) {
	//	std::cout << ~paths[bestPath][i].x;
	//	std::cout << ~paths[bestPath][i].u;
	//}
	
	//std::cout << std::endl << bestPaths.size() << std::endl;
	//for (std::multimap<double, int>::iterator j = bestPaths.begin(); j != bestPaths.end(); ++j) {
	//  std::cout << paths[j->second].size() << std::endl;
	//  for (int i = 0; i < (int) paths[j->second].size(); ++i) {
	//    std::cout << ~paths[j->second][i].x;
	//    std::cout << ~paths[j->second][i].u;
	//  }
	//}
}
*/

void generateStats()
{
	CAL_SuspendVisualisation();
	std::ofstream results("results.txt", std::ios::out);

	std::ifstream fptr(SAVEFILE, std::ios::in);
	int len;

	clock_t tstart, tstop;
	clock_t ttruncate = 0, tsampling = 0, tlqgmp = 0, tboole = 0;

	double pctruncate = 0.0, pcsampling = 0.0, pclqgmp = 0.0, pcboole = 0.0;

	for(int p = 0; p < 158; ++p) 
	{
		fptr >> len;
		std::vector<PathNode> path(len);
		for (int i = 0; i < len; ++i) {
			PathNode node;
			deserializeMatrix(fptr, node.x);
			deserializeMatrix(fptr, node.u);
			path[i] = node;
		}

		tstart = clock();
		pctruncate = computeLQGMPTruncate(path, false);
		tstop = clock();
		ttruncate += (tstop - tstart);
		results << pctruncate << " ";

		tstart = clock();
		pcsampling = computePSSampling(path, 10000, false);
		tstop = clock();
		tsampling += (tstop - tstart);
		results << pcsampling << " ";

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

void testPlan()
{
	std::ifstream fptr("carPath1.txt", std::ios::in);
	int len;

	fptr >> len;
	std::vector<PathNode> path(len);
	for (int i = 0; i < len; ++i) {
		PathNode node;
		deserializeMatrix(fptr, node.x);
		deserializeMatrix(fptr, node.u);
		path[i] = node;
	}
	
	double pctruncate = 0, pcsampling = 0, pclqgmp = 0, pcboole = 0;

	pctruncate = computeLQGMPTruncate(path, false);
	pcsampling = computePSSampling(path, 100000, false);
	computeQuality(path, true, false, true, pclqgmp, pcboole);
	std::cout << "PC Truncate: " << pctruncate << " PC Sampling: " << pcsampling << " " << " PC boole: " << pcboole << std::endl;
	//std::cout << "PC Sampling: " << pcsampling << std::endl;

	int ns[19] = {50,100,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,2750,3000};
	std::ofstream fptr2("mc-car.txt",std::ios::out);

	for(int i = 0; i < 19; ++i) {
		clock_t startTime = clock();
		for(int j = 0; j < 100; ++j) {
			pcsampling = computePSSampling(path, ns[i], false);
			fptr2 << std::setprecision(10) << pcsampling << " ";
		}
		fptr2 << std::endl;
		clock_t endTime = clock();
		std::cout << "Finished processing " << ns[i] << " samples in " << (double) (endTime - startTime) / CLOCKS_PER_SEC << " secs" << std::endl;
	}
	fptr2.close();

	std::cout << "Done" << std::endl;
	
	//std::cout << "PC Truncate: " << pctruncate << " PC Sampling: " << pcsampling << " " << " PC boole: " << pcboole << std::endl;
}



int main()
{
	time_t t;
	srand(unsigned(time(&t)));

	// Create Obstacles in Callisto
	CAL_Initialisation (true, true, true);
	
	initEnvironment();
	//initEnvironment2();

	//----------------------------------------------------
	// Construct Paths
	//std::ofstream frrt(SAVEFILE, std::ios::out);
	//clock_t startTime = clock();
	//int numPaths = 1;
	//std::cout << numPaths << std::endl;
	//for (int i = 0; i < numPaths; ++i) {
	//	rrt(frrt);
	//	std::cerr << i << " ";
	//}
	//frrt.close();
	//CAL_End();
	//clock_t endTime = clock();
	//std::cerr << std::endl << (double) (endTime - startTime) / CLOCKS_PER_SEC << std::endl;
	////showPaths(paths);
	//return 0;
	
	//std::ofstream frrt(SAVEFILE, std::ios::out);

	//clock_t tstart = clock();
	//int col;
	//for(int i = 0; i < 500; ++i) 
	//{
	//	do {
	//		start[0] = -8.0 + uniform()*16.0; start[1] = 1.5 + uniform()*7.0; start[2] = (0.5+(uniform()*1.5))*M_PI; start[3] = 0.1;
	//		CAL_CheckPointCollision(cal_environment, (float) start[0], (float) start[1], 0, false, &col);
	//	} while (col != 0);

	//	rrt(frrt);
	//	std::cout << i << " ";
	//}

	//frrt.close();
	//std::cout << "\nRRT generation time: " << (clock() - tstart) / (double) (CLOCKS_PER_SEC) << std::endl;

	//int num;
	//std::cin >> num;

	//----------------------------------------------------
	
	testPlan();
	//generateStats();

	int k;
	std::cin >> k;
	
	cleanup();
	return 0;
}