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
		U[0][0] = 3.0*U[1][0];
		U[imax-2][0] = U[imax-3][0];
		U[imax-1][0] = 3.0*U[imax-2][0];

		// density * velocity -- Neumann boundary condition (assumed a slip wall)
		U[1][1] = U[2][1];
		U[0][1] = 3.0*U[1][1];
		U[imax-2][1] = U[imax-3][1];
		U[imax-1][1] = 3.0*U[imax-2][1];

		// Energy -- Neumann boundary condition (assumed a slip wall)
		U[1][2] = U[2][2];
		U[0][2] = 3.0*U[1][2];
		U[imax-2][2] = U[imax-3][2];
		U[imax-1][2] = 3.0*U[imax-2][2];
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

vector<vector<double> > flux(vector<vector<double> > &U, string scheme)
{
	vector<vector<double> > F(imax, vector<double>(3));
	vector<vector<double> > Fplushalf(imax, vector<double>(3));
	vector<vector<double> > Fminushalf(imax, vector<double>(3));
	// Steger-Warming scheme
	if (scheme == "SW"){
		for (int i = 2; i < imax-2; i++)
		{
			//================
		 	// F_i-1/2
			//================
			vector<double> Up = Uplus(U, i-1);
			vector<double> Um = Uminus(U, i);

			//F+i-1
			vector<vector<double> > UoverV = jac(Up);
			vector<vector<double> > XR = Xmatrix(Up, "R");
			vector<vector<double> > lamPlus = bigLamb(Up, Um, "+");
			vector<vector<double> > XL = Xmatrix(Up, "L");
			vector<vector<double> > VoverU = jac(Up, "inv");

			vector<vector<double> >Aplus = MM(UoverV, XR);
			Aplus = MM(Aplus, lamPlus);
			Aplus = MM(Aplus, XL);
			Aplus = MM(Aplus, VoverU);
			vector<double> Fplus = MVM(Aplus, Up);


			// F-i
			UoverV = jac(Um);
			XR = Xmatrix(Um, "R");
			vector<vector<double> > lamMinus = bigLamb(Up, Um, "-");
			XL = Xmatrix(Um, "L");
			VoverU = jac(Um, "inv");

			vector<vector<double> > Aminus = MM(UoverV, XR);
			Aminus = MM(Aminus, lamMinus);
			Aminus = MM(Aminus, XL);
			Aminus = MM(Aminus, VoverU);
			vector<double> Fminus = MVM(Aminus, Um);

			Fminushalf[i] = Vadd(Fplus, Fminus);

			// ===============
			// F_i+1/2
			// ==============
			Up = Uplus(U, i);
			Um = Uminus(U, i+1);

			//F+i-1
			UoverV = jac(Up);
			XR = Xmatrix(Up, "R");
			lamPlus = bigLamb(Up, Um, "+");
			XL = Xmatrix(Up, "L");
			VoverU = jac(Up, "inv");

			Aplus = MM(UoverV, XR);
			Aplus = MM(Aplus, lamPlus);
			Aplus = MM(Aplus, XL);
			Aplus = MM(Aplus, VoverU);
			Fplus = MVM(Aplus, Up);

			// F-i
			UoverV = jac(Um);
			XR = Xmatrix(Um, "R");
			lamMinus = bigLamb(Up, Um, "-");
			XL = Xmatrix(Um, "L");
			VoverU = jac(Um, "inv");

			Aminus = MM(UoverV, XR);
			Aminus = MM(Aminus, lamPlus);
			Aminus = MM(Aminus, XL);
			Aminus = MM(Aminus, VoverU);
			Fminus = MVM(Aminus, Um);

			Fplushalf[i] = Vadd(Fplus, Fminus);
		// }
		F = Msub(Fplushalf, Fminushalf);
		F = ScaM(-1.0/dx, F);
	}
}
	else{
		cout << "Not a recognised scheme";
		exit;
	}
	return F;
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
//
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
	vector<vector<double> > U(imax, vector<double>(3));
	init(U);
	printVec2D(U);
	vector<vector<double> > F = flux(U, "SW");
	printVec2D(F);
	cout << dx;
	return 0;
}
