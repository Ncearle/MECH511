//=================================================
// Source file for Compressible 1D Sod Problem
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "vecmat_funcs.h"
#include "phys_funcs.h"

void ghost(vector<vector<double> > &U)
{
	// density -- Neumann boundary condition (no flux across boundary)
	U[1][0] = U[2][0];
	U[0][0] = U[1][0];
	U[imax-2][0] = U[imax-3][0];
	U[imax-1][0] = U[imax-2][0];

	// density * velocity -- Neumann boundary condition (assumed a slip wall)
	U[1][1] = U[2][1];
	U[0][1] = U[1][1];
	U[imax-2][1] = U[imax-3][1];
	U[imax-1][1] = U[imax-2][1];

	// Energy -- Neumann boundary condition (assumed a slip wall)
	U[1][2] = U[2][2];
	U[0][2] = U[1][2];
	U[imax-2][2] = U[imax-3][2];
	U[imax-1][2] = U[imax-2][2];
}

void init(vector<vector<double> > &U)
{
	for (int i = 2; i < (imax-4)/2+2; i++)
	{
		double x = (i - 0.5) / (imax - 2);

		// density
		U[i][0] = rhoL;
		// density * velocity
		U[i][1] = uL*rhoL;
		// energy
		U[i][2] = PL/(gam-1) + rhoL*uL*uL/2.;
	}
	for (int i = (imax-4)/2+2; i < imax-2; i++)
	{
		double x = (i - 0.5) / (imax - 2);

		// density
		U[i][0] = rhoR;
		// density * velocity
		U[i][1] = uR*rhoR;
		// energy
		U[i][2] = PR/(gam-1) + rhoR*uR*uR/2.;
	}
	ghost(U);
}

vector<vector<double> > jac(vector<double> &Upm, string inv="")	// Set inv="inv" to return inverse of the jacobian
{
	vector<vector<double> > J(3, vector<double>(3));
	double rho = density(Upm);
	double u = uVel(Upm);

	J[0][0] = 1.0;
	J[1][0] = u;
	J[1][1] = rho;
	J[2][0] = u*u/2.0;
	J[2][1] = rho*u;
	J[2][2] = 1.0/(gam-1.0);

	if (inv == "inv")
	{
		J[0][0] = 1.0;
		J[1][0] = -u/rho;
		J[1][1] = 1.0/rho;
		J[2][0] = (gam-1.0) * u*u/2.0;
		J[2][1] = -(gam-1.0) * u;
		J[2][2] = gam-1.0;
	}
	return J;
}

vector<vector<double> > Xmatrix(vector<double> &Upm, string LR) // Input L or R depending
{
	vector<vector<double> > X(3, vector<double>(3));
	double c = speedofsound(Upm);
	double rho = density(Upm);
	if (LR == "R"){
		X[0][0] = 1.0;
		X[0][1] = rho/c;
		X[0][2] = -X[0][1];
		X[1][1] = 1;
		X[1][2] = 1;
		X[2][1] = rho*c;
		X[2][2] = -X[2][1];
	}
	else if (LR == "L"){
		X[0][0] = 1.0;
		X[0][2] = -1.0/(c*c);
		X[1][1] = 0.5;
		X[1][2] = 1/(2*rho*c);
		X[2][1] = 0.5;
		X[2][2] = -X[1][2];
	}
	else{
		cout << "Xmatrix input error";
		exit;
	}
	return X;
}

vector<double> eig(vector<double> &Uplus, vector<double> &Uminus)
{
	vector<double> lambda(3);
	double uminus = uVel(Uminus);
	double uplus = uVel(Uplus);
	double cminus = speedofsound(Uminus);
	double cplus = speedofsound(Uplus);

	lambda[0] = (uminus + uplus)/2;
	lambda[1] = (uminus + cminus + uplus + cplus)/2;
	lambda[2] = (uminus - cminus + uplus - cplus)/2;
	return lambda;
}

vector<vector<double> > bigLamb(vector<double> &Uplus, vector<double> &Uminus, string PM)		// input i & PM = "+" for right side, i+1 & PM = "-" for left side
{
	vector<vector<double> > bigLamb(3, vector<double>(3));
	double eps = pow(10,-8);
	if (PM == "+"){
		vector<double> littleLamb = eig(Uplus, Uminus);
		double u = uVel(Uplus);
		double c = speedofsound(Uplus);
		if (littleLamb[0] > 0.0){bigLamb[0][0] = 0.5*(littleLamb[0] + sqrt(pow(littleLamb[0],2)+eps));}
		if (littleLamb[1] > 0.0){bigLamb[1][1] = 0.5*(littleLamb[1] + sqrt(pow(littleLamb[1],2)+eps));}
		if (littleLamb[2] > 0.0){bigLamb[2][2] = 0.5*(littleLamb[2] + sqrt(pow(littleLamb[2],2)+eps));}
	}
	else if (PM == "-"){
		vector<double> littleLamb = eig(Uplus, Uminus);
		double u = uVel(Uminus);
		double c = speedofsound(Uminus);
		if (littleLamb[0] < 0.0){bigLamb[0][0] = 0.5*(littleLamb[0] - sqrt(pow(littleLamb[0],2)+eps));}
		if (littleLamb[1] < 0.0){bigLamb[1][1] = 0.5*(littleLamb[1] - sqrt(pow(littleLamb[1],2)+eps));}
		if (littleLamb[2] < 0.0){bigLamb[2][2] = 0.5*(littleLamb[2] - sqrt(pow(littleLamb[2],2)+eps));}
	}
	return bigLamb;
}

vector<vector<double> > absBigLamb(vector<double> &Uplus, vector<double> &Uminus)
{
	vector<vector<double> > bigLamb(3, vector<double>(3));
	double eps = pow(10,-8);
	vector<double> littleLamb = eig(Uplus, Uminus);
	for (int i = 0; i < 3; i++){
		bigLamb[i][i] = sqrt(pow(littleLamb[i],2) + eps);
	}
	return bigLamb;
}

vector<vector<double> > getFlux(vector<vector<double> > &U, string scheme, int order, string limiter=" ")
{
	vector<vector<double> > flux(imax, vector<double>(3));
	// Steger-Warming scheme
	if (scheme == "SW" || scheme == "Steger-Warming"){
		for (int i = 2; i < imax-2; i++)
		{

			vector<double> Fminushalf(3);
			vector<double> Fplushalf(3);
			for (int k = 0; k <= 1; k++)
			{
				// Start with left interface F_i-1/2
				vector<double> Uplus = getPlus(U, i-1+k, 2, limiter);
				vector<double> Uminus = getMinus(U, i+k, 2, limiter);

				//F+
				vector<vector<double> > UoverV = jac(Uplus);
				vector<vector<double> > XR = Xmatrix(Uplus, "R");
				vector<vector<double> > Lambda = bigLamb(Uplus, Uminus, "+");
				vector<vector<double> > XL = Xmatrix(Uplus, "L");
				vector<vector<double> > VoverU = jac(Uplus, "inv");

				vector<vector<double> >Aplus = MM(UoverV, XR);
				Aplus = MM(Aplus, Lambda);
				Aplus = MM(Aplus, XL);
				Aplus = MM(Aplus, VoverU);
				vector<double> Fplus = MVM(Aplus, Uplus);

				// F-
				UoverV = jac(Uminus);
				XR = Xmatrix(Uminus, "R");
				Lambda = bigLamb(Uplus, Uminus, "-");
				XL = Xmatrix(Uminus, "L");
				VoverU = jac(Uminus, "inv");

				vector<vector<double> > Aminus = MM(UoverV, XR);
				Aminus = MM(Aminus, Lambda);
				Aminus = MM(Aminus, XL);
				Aminus = MM(Aminus, VoverU);
				vector<double> Fminus = MVM(Aminus, Uminus);

				if (k == 0){
					Fminushalf = Vadd(Fplus, Fminus);
				}
				else if (k == 1){
					Fplushalf = Vadd(Fplus, Fminus);
				}
			}
			flux[i] = Vsub(Fplushalf, Fminushalf);
			flux[i] = ScaV(-1.0/dx, flux[i]);
		}
}
	// Roe's Flux Differencing Scheme
	else if (scheme == "Roe" || scheme == "roe"){
		vector<vector<double> > F = getF(U);

		for (int i = 2; i < imax-2; i++)
		{
			vector<double> Fminushalf(3);
			vector<double> Fplushalf(3);

			for (int k = 0; k <= 1; k++){
				vector<double> Uplus = getPlus(U, i-1+k, 2, limiter);
				vector<double> Uminus = getMinus(U, i+k, 2, limiter);
				vector<double> Fplus = getPlus(F, i-1+k, 2, limiter);
				vector<double> Fminus = getMinus(F, i+k, 2, limiter);

				// Roe average
				double rightRho = density(Uminus);
				double leftRho = density(Uplus);
				double rightu = uVel(Uminus);
				double leftu = uVel(Uplus);
				double rightP = pressure(Uminus);
				double leftP = pressure(Uplus);
				double rightH = (Uminus[2]+rightP)/rightRho;
				double leftH = (Uplus[2]+leftP)/leftRho;

				double rhoTilde = sqrt(rightRho*leftRho);
				double uTilde = (sqrt(rightRho)*rightu + sqrt(leftRho)*leftu)/(sqrt(rightRho)*sqrt(leftRho));
				double hTilde = (sqrt(rightRho)*rightH + sqrt(leftRho)*leftH)/(sqrt(rightRho)+sqrt(leftRho));
				double PTilde = (hTilde - 0.5*uTilde*uTilde)*(gam-1)*rhoTilde/gam;
				double ETilde = PTilde/(gam-1) + 0.5*rhoTilde*uTilde*uTilde;
				vector<double> UTilde {rhoTilde, rhoTilde*uTilde, ETilde};

				vector<vector<double> > UoverV = jac(UTilde);
				vector<vector<double> > XR = Xmatrix(UTilde, "R");
				vector<vector<double> > Lambda = absBigLamb(Uplus, Uminus);
				vector<vector<double> > XL = Xmatrix(UTilde, "L");
				vector<vector<double> > VoverU = jac(UTilde, "inv");

				vector<vector<double> >ATilde = MM(UoverV, XR);
				ATilde = MM(ATilde, Lambda);
				ATilde = MM(ATilde, XL);
				ATilde = MM(ATilde, VoverU);

				vector<double> term1 = Vadd(Fplus, Fminus);
				term1 = ScaV(0.5, term1);
				vector<double> term2 = Vsub(Uminus, Uplus);
				term2 = MVM(ATilde, term2);
				term2 = ScaV(0.5, term2);

				if (k == 0){
					Fminushalf = Vsub(term1, term2);
				}
				else if (k == 1){
					Fplushalf = Vsub(term1, term2);
				}
			}
			flux[i] = Vsub(Fplushalf, Fminushalf);
			flux[i] = ScaV(-1.0/dx, flux[i]);
		}
	}
	else{
		cout << "Not a recognised scheme, guess again";
		exit;
	}
	return flux;
}

vector<vector<double> > EE(vector<vector<double> > &U, string scheme)
{
	vector<vector<double> > Unew(imax, vector<double>(3));
	int order = 1;
	vector<vector<double> > flux = getFlux(U, scheme, order);
	flux = ScaM(dt, flux);
	Unew = Madd(U, flux);
	ghost(Unew);
	return Unew;
}

vector<vector<double> > RK2(vector<vector<double> > &U, string scheme, string limiter=" ")
{
	vector<vector<double> > Unew(imax, vector<double>(3));
	int order = 2;
	// Intermediate Step
	vector<vector<double> > flux = getFlux(U, scheme, order, limiter);
	flux = ScaM(dt/2.0, flux);
	Unew = Madd(U, flux);
	ghost(Unew);

	// Full Step
	flux = getFlux(Unew, scheme, order, limiter);
	flux = ScaM(dt, flux);
	Unew = Madd(U, flux);
	ghost(Unew);

	return Unew;
}

int main()
{
	vector<vector<double> > U(imax, vector<double>(3));
	init(U);

	double t_fin = 0.15;
	int t_steps = t_fin/dt+1;

	for (int n = 0; n < t_steps; n++)
	{
		cout << n+1 << endl;
		U = RK2(U, "Roe", "minmod");

	}
	cout << "Unew:" << endl;
	U = transpose(U);
	printVec2D(U);
	U = transpose(U);

	vector<double> P(imax-4);
	vector<double> T(imax-4);
	vector<double> rho(imax-4);
	vector<double> u(imax-4);
	for (int i = 0; i < imax-4; i++)
	{
		P[i] = pressure(U[i+2]);
		T[i] = temperature(U[i+2]);
		rho[i] = density(U[i+2]);
		u[i] = uVel(U[i+2]);
	}

	string Pname = "P.dat";
	vec1D2File(Pname, P);
	string Tname = "T.dat";
	vec1D2File(Tname, T);
	string rhoname = "rho.dat";
	vec1D2File(rhoname, rho);
	string uname = "u.dat";
	vec1D2File(uname, u);

	string Dname = "data.dat";
	vec2D2File(Dname, U);

	return 0;
}
