#include "trunc.h"

// input z: determines the heaviside function [[x>z]], 
// output: mean and variances of the remaining probability mass when everything left of z is cut off
void TruncStandardGaussian(double& mean, double& var, double z)
{
	double norm = M_SQRT_PI_2 * (1.0 - erf(z*M_SQRT1_2));
	if(norm < 1e-2) { std::cout << "Likelihood of that truncation is very low: " << norm << std::endl; }
	mean = exp(-z*z*0.5)/norm;
	var  = 1.0 + z*mean - mean*mean;
}

void TruncGaussian(Matrix<SDIM>& m, Matrix<SDIM, SDIM>& S, const Constraint& c)
{
	// linear transform to make (a,A) to a standard Gaussian
	Matrix<SDIM, SDIM> M;
	cholesky(S, M);

	// transform and rescale the constraint 
	double z = tr(~c.a*m) - c.b;
	
	Matrix<SDIM> v = M*c.a;
	double norm_v = sqrt(tr(~v*v));
	if (norm_v < 1e-10) {
		std::cerr << "No gradient for Gaussian trunctation!" << std::endl;
		std::exit(-1);
	}
	z /= norm_v;
	v /= norm_v;
	
	// build rotation matrix for constraint to be along the x-axis
	Matrix<SDIM> e_1;
	e_1.reset(); e_1[0] = 1.0;
	Matrix<SDIM, SDIM> R;
	rotationFromAtoB(R, e_1, v);
	
	//std::cout << "R = " << R << ~R << !R << ~v << std::endl << R*e_1 << std::endl;

	// get mean and variance along x-axis
	double mean, var;
	TruncStandardGaussian(mean, var, -z);

	//std::cout << "z: " << z << " mean: " << mean << " var: " << var << std::endl;
	
	Matrix<SDIM, SDIM> B;
	Matrix<SDIM> b;
	b.reset(); b[0] = mean;
	B = identity<SDIM>(); B(0,0) = var;

	// transform back
	S = ~M*R*B*~R*M;
	m = ~M*R*b + m;

	//checkNaN(m);
	//checkNaN(S);
}

void TruncGaussianBySampling(Matrix<SDIM>& m, Matrix<SDIM, SDIM>& S, std::vector<Constraint>& cons, bool drawsamples = false)
{
	// generate samples from the Gaussian
	Matrix<SDIM, SDIM> M;
	Matrix<SDIM> x;
	std::vector<Matrix<SDIM> > X;
	cholesky(S, M);

	int ns = 10000;
	for(int i = 0; i < ns; ++i)
	{
		for(int j = 0; j < SDIM; ++j) x[j] = normal();
		x = ~M*x;
		x += m;

		//x = sampleGaussian(m, S);
		bool cond = true;
		for(std::vector<Constraint>::iterator it = cons.begin(); it != cons.end() && cond; ++it) {
			cond = cond && it->satisfy(x);
		}

		if(cond) {
			X.push_back(x);
			if (drawsamples) {
				int id;
				CAL_CreateUserDrawn(cal_samples, drawUnitCircle, NULL, (float)x[0], (float)x[1], 0.001f, &id);
				CAL_SetObjectScaling(id, 0.02f, 0.02f, 1.0f);
			}
		}
	}
	if(X.empty()){
		std::cout <<"TruncGaussianBySampling: no samples survived!" << std::endl;
		return;
	}

	std::cout << "Sampling probability of success: " << (double)X.size()/(double)ns << std::endl;

	int n = X.size();
	for(int i = 0; i < n; ++i) {
		m += X[i];
	}
	m = m/n;
	for(int i = 0; i < n; ++i) {
		S += (X[i] - m)*~(X[i] - m);
	}
	S = S/n;

	std::cout << "Sampled mean: " << ~m << std::endl;
	std::cout << "Sampled var: " << S << std::endl;
}

void testTruncGauss()
{
	// random covariance and mean
	Matrix<SDIM, SDIM> S;
	Matrix<SDIM> m, cdir, coff;
	std::vector<Constraint> cons;

	/*
	Matrix<SDIM, SDIM> M;
	rndMatrix(M, 0.1, 1.0);
	S = ~M*M;
	std::cout << "S: " << S << std::endl;
	
	rndMatrix(m, -1.0, 1.0);
	std::cout << "m: " << ~m << std::endl;
	//m[0] = 0.0;  m[1] = 0.0;  
	//S(0,0) = 1.0; S(0,1) = 1.0; S(1,0) = 1.0; S(1,1) = 2.0;

	//random constraint
	rndMatrix(cdir, -1.0, 1.0);
	rndMatrix(coff, -1.0, 1.0);
	//cdir[0] = 1.0; cdir[1] = 1.0; coff[0] = 1.0; coff[1] = 0.0; 
	std::cout << "cdir: " << ~cdir << std::endl;
	std::cout << "coff: " << ~coff << std::endl;
	*/

	/*
	S(0,0) = 1.54641; S(0,1) = 1.226096;
	S(1,0) = 1.226096; S(1,1) = 1.361004;

	m[0] = 1.0; m[1] = 1.5;
	
	// First constraint
	cdir[0] = 0.806037; cdir[1] = 0.389368;
	coff[0] = 0.151189; coff[1] = 0.81633;
	cons.push_back(Constraint(cdir, coff));

	// Second constraint
	cdir[0] = 0.0; cdir[1] = -1.0;
	coff[0] = 2.251189; coff[1] = 2.31633;
	cons.push_back(Constraint(cdir, coff));

	// Third constraint
	cdir[0] = -1.295725; cdir[1] = 0.24875;
	coff[0] = 1.251189; coff[1] = 0.31633;
	cons.push_back(Constraint(cdir, coff));
	*/

	/*
	for(int i = 0; i < 4; ++i) {
		//random constraint
		rndMatrix(cdir, -1.0, 1.0);
		rndMatrix(coff, -1.0, 1.0);
		cons.push_back(Constraint(cdir, coff));
	}
	*/

	S(0,0) = 0.104258; S(0,1) = 0.00860636;
	S(1,0) = 0.00860636; S(1,1) = 0.00646838;

	m[0] = -7.44845; m[1] = -1.14356;
	
	Matrix<SDIM> a;
	double b;
	// First constraint
	a[0] = 1.0; a[1] = 0.0;
	b = -7.5;
	cons.push_back(Constraint(a, b));

	// Second constraint
	a[0] = -1.0; a[1] = 0.0;
	b = 5.5;
	cons.push_back(Constraint(a, b));
	
	int np[1] = {2};
	int id;
	float side1[6] = {(float)-7.5, (float)-5.0, 0.001f, (float)-7.5, (float)5.0, 0.001f};
	CAL_CreatePolyline(cal_constraint, 1, np, side1, &id);

	float side2[6] = {(float)-5.5, (float)-5.0, 0.001f, (float)-5.5, (float)5.0, 0.001f};
	CAL_CreatePolyline(cal_constraint, 1, np, side2, &id);

	//std::random_shuffle(cons.begin(), cons.end());

	// analytical solution
	Matrix<SDIM> mn(m);
	Matrix<SDIM, SDIM> Sn(S);

	drawCovariance(mn, Sn);

	clock_t startClk = clock();
	double ps2 = 1.0;
	for(std::vector<Constraint>::iterator it = cons.begin(); it != cons.end(); ++it) {
		std::cout << "Processing constraint: " << ~it->a << " " << it->b << std::endl;
		//drawCovariance(mn, Sn);
		std::cout << "ps: " << (1.0 - it->probFailure(mn, Sn)) << std::endl;
		ps2 *= (1.0 - it->probFailure(mn, Sn));
		TruncGaussian(mn, Sn, (*it));
	}
	std::cout << "Probability of success (Method 2): " << ps2 << std::endl;
	std::cout << "Analytical: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
	drawCovariance(mn, Sn);
	std::cout << "Mean: " << ~mn << std::endl;
	std::cout << "Var: " << Sn << std::endl;

	// Boole's inequality - prob success
	mn = m; Sn = S;
	double ps1 = 1.0;
	for(std::vector<Constraint>::iterator it = cons.begin(); it != cons.end(); ++it) {
		ps1 -= it->probFailure(mn, Sn);
	}
	std::cout << "Probability of success (Boole): " << ps1 << std::endl;

	//sampling solution
	mn = m; Sn = S;
	startClk = clock();
	TruncGaussianBySampling(mn, Sn, cons, false);
	std::cout << "Sampling: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
	//drawCovariance(mn, Sn);
}

void testTruncGauss2()
{
	// random covariance and mean
	Matrix<SDIM, SDIM> S;
	Matrix<SDIM> m, a;
	double b;
	std::vector<Constraint> cons;

	S(0,0) = 0.104258; S(0,1) = 0.00860636;
	S(1,0) = 0.00860636; S(1,1) = 0.00646838;

	m[0] = -7.44845; m[1] = -1.14356;
	
	int np[1] = {2};
	int id;
	
	// First constraint
	a[0] = 1.0; a[1] = 0.0;
	b = -7.5;
	cons.push_back(Constraint(a, b));
	
	float side1[6] = {(float)3.0, (float)-5.0, 0.001f, (float)3.0, (float)5.0, 0.001f};
	CAL_CreatePolyline(cal_constraint, 1, np, side1, &id);

	// Second constraint
	a[0] = -1.0; a[1] = 0.0;
	b = 5.5;
	cons.push_back(Constraint(a, b));

	cons[1].drawConstraint();
	float side2[6] = {(float)0.0, (float)-5.0, 0.001f, (float)0.0, (float)5.0, 0.001f};
	CAL_CreatePolyline(cal_constraint, 1, np, side2, &id);
	
	//std::random_shuffle(cons.begin(), cons.end());

	// analytical solution
	Matrix<SDIM> mn(m);
	Matrix<SDIM, SDIM> Sn(S); 

	drawCovariance(mn, Sn);

	clock_t startClk = clock();
	double yMean, yVar, yNewMean, yNewVar;
	Matrix<SDIM> xyCovar, L;
	double ps = 1.0, z;
	int nc = (int)cons.size();

	for(int i = 0; i < nc; ++i) 
	{
		Constraint & c = cons[i];
		yMean = tr(~c.a*m);
		yVar = tr(~c.a*S*c.a);

		z = truncate(c.b, yMean, yVar, yNewMean, yNewVar);
		std::cout << "z: " << z << std::endl;
		ps *= (1.0 - z);	

		xyCovar = S*c.a;
		L = xyCovar/yVar;
		mn += L*(yNewMean - yMean);
		Sn += L*(yNewVar - yVar)*~L;
	}
	std::cout << "Probability of success: " << ps << std::endl;
	std::cout << "Analytical: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
	drawCovariance(mn, Sn);
	std::cout << "Mean: " << ~mn << std::endl;
	std::cout << "Var: " << Sn << std::endl;

	/*
	//sampling solution
	mn = m; Sn = S;
	startClk = clock();
	TruncGaussianBySampling(mn, Sn, cons, false);
	std::cout << "Sampling: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
	drawCovariance(mn, Sn);
	*/
}

int main()
{
	srand (unsigned(time(NULL)));

	initEnvironment();
	
	testTruncGauss();
	//testTruncGauss2();
	
	int num;
	std::cin >> num;

	CAL_End ();

	return 0;
}

