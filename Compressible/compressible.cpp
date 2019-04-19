//=================================================
// Source file for Navier-Stokes Equation
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "vecmat_funcs.h"
#include "phys_funcs.h"
#include "blocktri_vector.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<double> > &U)
{
	// // Top and Bottom Ghost Cells
	// for (int i = 0; i < imax; i++)
	// {
	// 	double x = (i - 0.5) / (imax - 2);
	// 	// pressure -- Neumann boundary condition (no flux across boundary)
	// 	U[0][0][i] = U[0][1][i];
	// 	U[0][jmax-1][i] = U[0][jmax-2][i];
	//
	// 	// density * velocity -- Neumann boundary condition (assumed a slip wall)
	// 	U[1][0][i] = U[1][1][i];
	// 	U[1][jmax-1][i] = U[1][jmax-2][i];
	//
	// 	// energy -- Neumann boundary condition (assumed a slip wall)
	// 	U[2][0][i] = U[2][1][i];
	// 	U[2][jmax-1][i] = U[2][jmax-2][i];
	// }

	// // Left and Right Ghost Cells
	// for (int j = 1; j < jmax-1; j++)
	// {
	// 	double y = (j - 0.5) / (jmax - 2);
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

// Initializes the domain (Pressure and Velocity)
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
		U[i][2] = rhoL*Cv*inittemp(U[i][0], PL) + rhoL*uL*uL/2.;
	}
	for (int i = (imax-4)/2+2; i < imax-2; i++)
	{
		double x = (i - 0.5) / (imax - 2);

		// density
		U[i][0] = rhoR;
		// density * velocity
		U[i][1] = uR*rhoR;
		// energy
		U[i][2] = rhoR*Cv*inittemp(U[i][0], PR) + rhoR*uR*uR/2.;
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
	if (PM == "+"){
		vector<double> littleLamb = eig(Uplus, Uminus);
		double u = uVel(Uplus);
		double c = speedofsound(Uplus);
		if (littleLamb[0] > 0){bigLamb[0][0] = 0.5*(u + abs(u));}
		else if (littleLamb[1] > 0){bigLamb[1][1] = 0.5*(u+c + abs(u+c));}
		else if (littleLamb[2] > 0){bigLamb[2][2] = 0.5*(u-c + abs(u-c));}
	}
	else if (PM == "-"){
		vector<double> littleLamb = eig(Uplus, Uminus);
		double u = uVel(Uminus);
		double c = speedofsound(Uminus);
		if (littleLamb[0] < 0){bigLamb[0][0] = 0.5*(u - abs(u));}
		else if (littleLamb[1] < 0){bigLamb[1][1] = 0.5*(u+c - abs(u+c));}
		else if (littleLamb[2] < 0){bigLamb[2][2] = 0.5*(u-c - abs(u-c));}
	}
	return bigLamb;
}

vector<vector<double> > absBigLamb(vector<double> &Uplus, vector<double> &Uminus)
{
	vector<vector<double> > bigLamb(3, vector<double>(3));
	vector<double> littleLamb = eig(Uplus, Uminus);
	double u = uVel(Uminus);
	double c = speedofsound(Uminus);
	for (int i = 0; i < 3; i++){
		bigLamb[i][i] = abs(littleLamb[i]);
	}
	return bigLamb;
}

vector<vector<double> > getFlux(vector<vector<double> > &U, string scheme)
{
	vector<vector<double> > flux(imax, vector<double>(3));
	// Steger-Warming scheme
	if (scheme == "SW"){
		for (int i = 2; i < imax-2; i++)
		{
			vector<double> Fminushalf(3);
			vector<double> Fplushalf(3);
			for (int k = 0; k <= 1; k++)
			{
				vector<double> Uplus = getPlus(U, i-1+k);
				vector<double> Uminus = getMinus(U, i+k);

				//F+i-1
				vector<vector<double> > UoverV = jac(Uplus);
				vector<vector<double> > XR = Xmatrix(Uplus, "R");
				vector<vector<double> > lamPlus = bigLamb(Uplus, Uminus, "+");
				vector<vector<double> > XL = Xmatrix(Uplus, "L");
				vector<vector<double> > VoverU = jac(Uplus, "inv");

				vector<vector<double> >Aplus = MM(UoverV, XR);
				Aplus = MM(Aplus, lamPlus);
				Aplus = MM(Aplus, XL);
				Aplus = MM(Aplus, VoverU);
				vector<double> Fplus = MVM(Aplus, Uplus);

				// F-i
				UoverV = jac(Uminus);
				XR = Xmatrix(Uminus, "R");
				vector<vector<double> > lamMinus = bigLamb(Uplus, Uminus, "-");
				XL = Xmatrix(Uminus, "L");
				VoverU = jac(Uminus, "inv");

				vector<vector<double> > Aminus = MM(UoverV, XR);
				Aminus = MM(Aminus, lamMinus);
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
				vector<double> Uplus = getPlus(U, i-1+k);
				vector<double> Uminus = getMinus(U, i+k);
				vector<double> Fplus = getPlus(F, i-1+k);
				vector<double> Fminus = getMinus(F, i+k);

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
				double hTilde = (sqrt(rightRho)*rightH + sqrt(leftRho)*leftH)/(sqrt(rightRho)*sqrt(leftRho));
				double PTilde = sqrt(rightP*leftP);
				double ETilde = hTilde*rhoTilde-PTilde;
				vector<double> UTilde = {rhoTilde, rhoTilde*uTilde, ETilde};

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

// Two Stage Runge Kutta time adt_stepsvance
vector<vector<double> > RK2(vector<vector<double> > &U, string scheme)
{
	vector<vector<double> > Unew(imax, vector<double>(3));

	// Intermediate Step
	// cout << "Intermediate step:" << endl;
	vector<vector<double> > flux = getFlux(U, scheme);
	// flux = transpose(flux);
	// cout << "flux:"<< endl;
	// printVec2D(flux);
	// flux = transpose(flux);

	flux = ScaM(dt/2.0, flux);
	// flux = transpose(flux);
	// cout << "flux*dt/2:"<< endl;
	// printVec2D(flux);
	// flux = transpose(flux);

	Unew = Madd(U, flux);
	ghost(Unew);
	// Unew = transpose(Unew);
	// cout << "Unew:"<< endl;
	// printVec2D(Unew);
	// Unew = transpose(Unew);


	// Full Step
	// cout << "Full step:" << endl;
	flux = getFlux(Unew, scheme);
	// flux = transpose(flux);
	// cout << "flux:"<< endl;
	// printVec2D(flux);
	// flux = transpose(flux);

	flux = ScaM(dt, flux);
	// flux = transpose(flux);
	// cout << "flux*dt:"<< endl;
	// printVec2D(flux);
	// flux = transpose(flux);

	Unew = Madd(U, flux);
	ghost(Unew);
	// Unew = transpose(Unew);
	// cout << "Unew:"<< endl;
	// printVec2D(Unew);
	// Unew = transpose(Unew);

	return Unew;
}

int main()
{
	vector<vector<double> > U(imax, vector<double>(3));
	init(U);
	U = transpose(U);
	printVec2D(U);
	U = transpose(U);
	vector<vector<double> > FI = getFlux(U, "Roe");
	// FI = transpose(FI);
	// printVec2D(FI);
	// FI = transpose(FI);

	double t_fin = 0.15;
	int t_steps = t_fin/dt+1;

	for (int n = 0; n < t_steps; n++)
	{
		cout << n+1 << endl;
		U = RK2(U, "SW");
		U = transpose(U);
		printVec2D(U);
		U = transpose(U);
	}
	//
	// string Dname = "data.dat";
	// vec2D2File(Dname, U);

	return 0;
}
