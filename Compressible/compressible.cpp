//=================================================
// Source file for Navier-Stokes Equation
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
// #include "exact.h"
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
		U[0][1] = U[0][2];
		U[0][0] = 3.0*U[0][1];
		U[0][imax-2] = U[0][imax-3];
		U[0][imax-1] = 3.0*U[0][imax-2];

		// density * velocity -- Neumann boundary condition (assumed a slip wall)
		U[1][1] = U[1][2];
		U[1][0] = 3.0*U[1][1];
		U[1][imax-2] = U[1][imax-3];
		U[1][imax-1] = 3.0*U[1][imax-2];

		// Energy -- Neumann boundary condition (assumed a slip wall)
		U[2][1] = U[2][2];
		U[2][0] = 3.0*U[2][1];
		U[2][imax-2] = U[2][imax-3];
		U[2][imax-1] = 3.0*U[2][imax-2];
}

// Initializes the domain (Pressure and Velocity)
void init(vector<vector<double> > &U)
{
	for (int i = 2; i < (imax-4)/2+2; i++)
	{
		double x = (i - 0.5) / (imax - 2);

		// density
		U[0][i] = rhoL;
		// density * velocity
		U[1][i] = uL*rhoL;
		// energy
		U[2][i] = rhoL*Cv*inittemp(U[0][i], PL) + rhoL*uL*uL/2.;
	}
	for (int i = (imax-4)/2+2; i < imax-2; i++)
	{
		double x = (i - 0.5) / (imax - 2);

		// density
		U[0][i] = rhoR;
		// density * velocity
		U[1][i] = uR*rhoR;
		// energy
		U[2][i] = rhoR*Cv*inittemp(U[0][i], PR) + rhoR*uR*uR/2.;
	}
	ghost(U);
}

vector<vector<double> > jac(vector<vector<double> > &U, int i, string PM, string inv="")	// Set inv="inv" to return inverse of the jacobian
{
	vector<vector<double> > J(3, vector<double>(3));
	double rho = density(U, i, PM);
	double u = uVel(U, i, PM);

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

vector<vector<double> > Xmatrix(vector<vector<double> > &U, int i, string PM, string LR) // Input L or R depending
{
	vector<vector<double> > X(3, vector<double>(3));
	double c = speedofsound(U, i, PM);
	double rho = density(U, i, PM);
	if (LR == "R")
	{
		X[0][0] = 1.0;
		X[0][1] = rho/c;
		X[0][2] = -X[0][1];
		X[1][1] = 1;
		X[1][2] = 1;
		X[2][1] = rho*c;
		X[2][2] = -X[2][1];
	}
	else if (LR == "L")
	{
		X[0][0] = 1.0;
		X[0][2] = -1.0/(c*c);
		X[1][1] = 0.5;
		X[1][2] = 1/(2*rho*c);
		X[2][1] = 0.5;
		X[2][2] = -X[1][2];
	}
}

vector<double> eig(vector<vector<double> > &U, int i, string LR)
{
	vector<double> lambda(3);
	double uminus, uplus, cminus, cplus;
	if (LR == "L"){
		uplus = uVel(U, i-1, "+");
		uminus = uVel(U, i, "-");
		cplus = speedofsound(U, i-1, "+");
		cminus = speedofsound(U, i, "-");
	}
	else if (LR == "R"){
		uplus = uVel(U, i, "+");
		uminus = uVel(U, i+1, "-");
		cplus = speedofsound(U, i, "+");
		cminus = speedofsound(U, i+1, "-");
	}
	lambda[0] = (uminus + uplus)/2;
	lambda[1] = (uminus + cminus + uplus + cplus)/2;
	lambda[2] = (uminus - cminus + uplus - cplus)/2;
	return lambda;
}

vector<vector<double> > bigLamb(vector<vector<double> > &U, int i, string PM)		// input i & PM = "+" for right side, i+1 & PM = "-" for left side
{
	vector<vector<double> > bigLamb(3, vector<double>(3));
	if (PM == "+"){
		vector<double> eig = eig(U, i);
		double u = uVel(U, i, "+");
		double c = speedofsound(U, i, "+");
		if (eig[0] > 0){bigLamb[0][0] = 0.5*(u + abs(u));}
		else if (eig[1] > 0){bigLamb[1][1] = 0.5*(u+c + abs(u+c));}
		else if (eig[2] > 0){bigLamb[2][2] = 0.5*(u-c + abs(u-c));}
	}
	else if (PM == "-"){
		vector<double> eig = eig(U, i);
		double u = uVel(U, i, "-");
		double c = speedofsound(U, i, "-");
		if (eig[0] < 0){bigLamb[0][0] = 0.5*(u - abs(u));}
		else if (eig[1] < 0){bigLamb[1][1] = 0.5*(u+c - abs(u+c));}
		else if (eig[2] < 0){bigLamb[2][2] = 0.5*(u-c - abs(u-c));}
	}
	return bigLamb;
}

vector<vector<double> > flux(vector<vector<double> > &U, string scheme)
{
	vector<vector<double> > flux(3, vector<double>(imax));
	vector<vector<double> > Fplushalf(3, vector<double>(imax));
	vector<vector<double> > Fminushalf(3, vector<double>(imax));
	// Steger-Warming scheme
	if (scheme == "SW"){
		// F_i+
		for (int i = 2; i < imax-2; i++)
		{
			// Fi+ for F_i+1/2
			vector<vector<double> > UoverV = jac(U, i);
			vector<vector<double> > XR = Xmatrix(U, i, "R");
			vector<vector<double> > lamPlus = bigLamb(U, i, "+");
			vector<vector<double> > XL = Xmatrix(U, i, "L");
			vector<vector<double> > VoverU = jac(U, i, "inv");

			vector<vector<double> >term1 = MM(UoverV, XR);
			term1 = MM(term1, lamPlus);
			term1 = MM(term1, XL);
			term1 = MM(term1, VoverU);
			vector<double> Uplus = Uplus(U, i);
			vector<double> Fplus = MVM(term1, Uplus);

			// Fi+1- for F_i+1/2
			UoverV = jac(U, i+1);
			XR = Xmatrix(U, i+1, "R");
			vector<vector<double> > lamminus = bigLamb(U, i+1, "-");
			XL = Xmatrix(U, i+1, "L");
			VoverU = jac(U, i+1, "inv");

			term1 = MM(UoverV, XR);
			term1 = MM(term1, lamPlus);
			term1 = MM(term1, XL);
			term1 = MM(term1, VoverU);
			vector<double> Uminus = Uminus(U, i+1);
			vector<double> Fminus = MVM(term1, Uminus);

			Fplushalf[i] = Vadd(Fplus, Fminus);
			// ========================================================================
			// Fi+ for F_i-1/2
			UoverV = jac(U, i-1);
			XR = Xmatrix(U, i-1, "R");
			lamPlus = bigLamb(U, i-1, "+");
			XL = Xmatrix(U, i-1, "L");
			VoverU = jac(U, i-1, "inv");

			term1 = MM(UoverV, XR);
			term1 = MM(term1, lamPlus);
			term1 = MM(term1, XL);
			term1 = MM(term1, VoverU);
			Uplus = Uplus(U, i-1);
			Fplus = MVM(term1, Uplus);

			// Fi+1- for F_i-1/2
			UoverV = jac(U, i);
			XR = Xmatrix(U, i, "R");
			lamminus = bigLamb(U, i, "-");
			XL = Xmatrix(U, i, "L");
			VoverU = jac(U, i, "inv");

			term1 = MM(UoverV, XR);
			term1 = MM(term1, lamPlus);
			term1 = MM(term1, XL);
			term1 = MM(term1, VoverU);
			Uminus = Uminus(U, i);
			Fminus = MVM(term1, Uminus);

			Fminushalf[i] = Vadd(Fplus, Fminus);

		}
		flux = -1/dx * Msub(Fplushalf, Fminushalf);
	}
	return flux;
}

// Two Stage Runge Kutta time advance
// void RK2(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &S)
// {
// 	vector<vector<double> > T0(jmax, vector<double>(imax));
// 	int it = 0;
// 	double delta = 1.0;
// 	while (abs(delta) > tol)
// 	{
// 		it++;
// 		// Intermediate Step
// 		vector<vector<double> > FIint = FI2C(T, u, v, S);
// 		vector<vector<double> > Tint(jmax, vector<double>(imax));
// 		for (int j = 1; j < jmax-1; j++)
// 		{
// 			for (int i = 1; i < imax-1; i++)
// 			{
// 				T0[j][i] = T[j][i];
// 				Tint[j][i] = T[j][i] + dt/2.0 * FIint[j][i];
// 			}
// 		}
// 		ghost(Tint, u, v);
//
// 		// Full Step
// 		vector<vector<double> > FI = FI2C(Tint, u, v, S);
// 		for (int j = 1; j < jmax-1; j++)
// 		{
// 			for (int i = 1; i < imax-1; i++)
// 			{
// 				T[j][i] = T[j][i] + dt * FI[j][i];
// 			}
// 		}
// 		ghost(T, u, v);
//
// 		delta = maxChange(T0, T);
// 	}
// 	printVec2D(T);
// 	cout << setprecision(6) << delta << endl;
// 	vector<vector<double> > ExT = exactTemp();
// 	vector<vector<double> > err = error(T, ExT);
// 	double L2 = L2Norm(err);
// 	cout << setprecision(6) << "L2 norm: " << L2 << endl;
// 	cout << "Iterations: " << it << endl;
// 	cout << "Timestep: " << dt << endl;
// 	cout << "Tolerance: " << tol << endl;
// }

// void Imp(vector<vector<vector<double>>> &U)
// {
// 	init(U);
//
// 	vector<vector<double>> I = Id(3);						// 3 x 3 identity matrix
// 	vector<double> ZV(3, 0.0);								// vector of zeros
// 	vector<vector<double>> ZM(3, vector<double>(3, 0.0));	// 3 x 3 matrix of zeros
// 	vector<vector<double>> L2(0, vector<double>(3));	// L2 norms for each timestep based on previous results
//
// 	int it = 0;
// 	double tol = pow(10, -12);
// 	double maxL2 = 1.0;
// 	double t = dt;
// 	while (maxL2 > tol)
// 	{
// 		it ++;
// 		t = dt * it;
// 		vector<vector<vector<double>>> FI = flux(U);
// 		// cout << "Flux Integral: " << endl;
// 		// printTable(FI);
// 		vector<vector<vector<double>>> U0 = copy3(U);
// 		vector<vector<vector<vector<double>>>> DX(imax, vector<vector<vector<double>>>(3, vector<vector<double>>(3, vector<double>(3))));
// 		vector<vector<double>> FIx(imax, vector<double>(3));
// 		vector<vector<vector<double>>> Utilda(jmax, vector<vector<double>>(imax, vector<double>(3)));
// 		// double A = 0.00001;
//
// 		for (int j = 1; j < jmax - 1; j++)
// 		{
// 			for (int i = 1; i < imax - 1; i++)
// 			{
// 				vector<vector<double>> Ax = jac(U, "F", -1, -1, j, i);
// 				vector<vector<double>> Bxp = jac(U, "F", 1, -1, j, i);
// 				vector<vector<double>> Bxm = jac(U, "F", -1, 1, j, i);
// 				vector<vector<double>> Cx = jac(U, "F", 1, 1, j, i);
// 				vector<vector<double>> Bx = Msub(Bxp, Bxm);
//
// 				Ax = ScaM(-dt, Ax);
// 				Bx = ScaM(dt, Bx);
// 				Cx = ScaM(dt, Cx);
//
// 				DX[i][0] = Ax;
// 				DX[i][1] = Madd(I, Bx);
// 				DX[i][2] = Cx;
//
// 				FIx[i] = ScaV(dt, FI[j][i]);
// 				// printVec2D(FIx);
// 				// double lapP = ((FI[j][i+1][0] - 2 * FI[j][i][0] + FI[j][i-1][0])/pow(dx,2) + (FI[j+1][i][0] - 2 * FI[j][i][0] + FI[j-1][i][0])/pow(dy,2));
// 				// FIx[i][0] += A*dx*dy*lapP;
// 				// printVec2D(FIx);
// 			}
//
// 			DX[0][0] = ZM;
// 			DX[0][1] = I;
// 			DX[0][2] = I;
// 			DX[0][2][0][0] = -1;
//
// 			DX[imax-1][0] = I;
// 			DX[imax-1][1] = I;
// 			DX[imax-1][1][0][0] = -1;
// 			DX[imax-1][2] = ZM;
//
// 			FIx[0] = ZV;
// 			// FIx[0][1] = 2;
// 			FIx[imax-1] = ZV;
//
// 			SolveBlockTri(DX, FIx, imax);
//
// 			for (int i = 0; i < imax; i++)
// 			{
// 				Utilda[j][i] = FIx[i];
// 			}
// 		}
// 		// cout << "Solving for lines of constant J:" << endl;
// 		// printTable(Utilda);
//
// 		vector<vector<vector<vector<double>>>> DY(jmax, vector<vector<vector<double>>>(3, vector<vector<double>>(3, vector<double>(3))));
// 		vector<vector<double>> Uty(jmax, vector<double>(3));
// 		vector<vector<vector<double>>> deltaU(jmax, vector<vector<double>>(imax, vector<double>(3)));
//
//
// 		for (int i = 1; i < imax - 1; i++)
// 		{
// 			for (int j = 1; j < jmax - 1; j++)
// 			{
// 				vector<vector<double>> Ay = jac(U, "G", -1, -1, j, i);
// 				vector<vector<double>> Byp = jac(U, "G", 1, -1, j, i);
// 				vector<vector<double>> Bym = jac(U, "G", -1, 1, j, i);
// 				vector<vector<double>> Cy = jac(U, "G", 1, 1, j, i);
// 				vector<vector<double>> By = Msub(Byp, Bym);
//
// 				Ay = ScaM(-dt, Ay);
// 				By = ScaM(dt, By);
// 				Cy = ScaM(dt, Cy);
//
// 				DY[j][0] = Ay;
// 				DY[j][1] = Madd(I, By);
// 				DY[j][2] = Cy;
//
// 				Uty[j] = Utilda[j][i];
// 				// double lapP = ((FI[j][i+1][0] - 2 * FI[j][i][0] + FI[j][i-1][0])/pow(dx,2) + (FI[j+1][i][0] - 2 * FI[j][i][0] + FI[j-1][i][0])/pow(dy,2));
// 				// Uty[j][0] += A*dx*dy*lapP;
// 			}
//
// 			DY[0][0] = ZM;
// 			DY[0][1] = I;
// 			DY[0][2] = I;
// 			DY[0][2][0][0] = -1;
//
//
// 			DY[jmax-1][0] = I;
// 			DY[jmax-1][1] = I;
// 			DY[jmax-1][1][0][0] = -1;
// 			DY[jmax-1][2] = ZM;
//
// 			Uty[0] = ZV;
// 			// Uty[0][1] = 2;
// 			Uty[jmax-1] = ZV;
//
// 			SolveBlockTri(DY, Uty, jmax);
//
// 			for (int j = 0; j < jmax; j++)
// 			{
// 				deltaU[j][i] = Uty[j];
// 			}
// 		}
// 		// cout << "Solving for lines of constant I: " << endl;
// 		// printTable(deltaU);
//
// 		deltaU = ScaM3(w, deltaU);
// 		U = Madd3D(U, deltaU);
// 		ghost(U);
// 		vector<double> L2it = L2Norm(Msub3D(U, U0));
// 		maxL2 = MaxV(L2it);
// 		L2.push_back (L2it);
// 	}
//
// 	cout << "Iterations: " << it << endl;
// 	// cout << "Time: " << t << endl;
// 	cout << "L2 norm: " << endl;
// 	printVec(L2[it-1]);
//
// 	// printTable(U);
//
// 	// printVec2D(L2);
// 	string L2name = "L2_U1m.dat";
// 	vec2D2File(L2name, L2);
//
// }

int main()
{
	vector<vector<double> > U(3, vector<double>(imax));
	init(U);
	printVec2D(U);
	vector<vector<double> > F = flux(U, "SW");
	printVec2D(F);
	return 0;
}