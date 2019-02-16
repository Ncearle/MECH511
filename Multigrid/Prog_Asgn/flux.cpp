//=================================================
// Source file for Navier-Stokes Equation
//=================================================
#include "constant.h"
#include "print_fcns.h"
#include "error_fcns.h"
#include "exact.h"
#include "thomas.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<vector<double>>> &U)
{
	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		double x = (i - 0.5) / (imax - 2);
		// u velocity
		U[0][i][1] = -U[1][i][1];
		U[jmax-1][i][1] = -U[jmax-2][i][1];

		//v velocity
		U[0][i][2] = -U[1][i][2];
		U[jmax-1][i][2] = -U[jmax-2][i][2];
	}

	// Left and Right Ghost Cells
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		U[j][0][1] = -U[j][1][1];
		U[j][imax-1][1] = -U[j][imax-2][1];

		U[j][0][2] = -U[j][1][2];
		U[j][imax-1][2] = -U[j][imax-2][2];
	}
}

// Initializes the domain (Pressure and Velocity)
void init(vector<vector<vector<double>>> &U)
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			// Pressure
			U[j][i][0] = P0 * cos(pi*x) * cos(pi * y);

			// u velocity
			U[j][i][1] = u0 * sin(pi*x) * sin(2*pi*y);

			// v velocity
			U[j][i][2] = v0 * sin(2*pi*x) * sin(pi*y);
		}
	}
	ghost(U);
}

vector<vector<vector<double>>> flux(vector<vector<vector<double>>> &U)
{
	vector<vector<vector<double>>> flux(jmax, vector<vector<double>>(imax, vector<double>(3)));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			vector<double> Fp(3);	// F_i+1/2,j
			vector<double> Fm(3);	// F_i-1/2,j
			vector<double> Gp(3);	// G_i,j+1/2
			vector<double> Gm(3);	// G_i,j-1/2

			Fp[0] = (U[j][i+1][1] + U[j][i][1]) / (2*B);
			Fp[1] = pow((U[j][i+1][1] + U[j][i][1])/2, 2) + (U[j][i+1][0] + U[j][i][0])/2 - (U[j][i+1][1] - U[j][i][1])/(Re*dx);
			Fp[2] = ((U[j][i+1][1] + U[j][i][1]) / 2) * ((U[j][i+1][2] + U[j][i][2]) / 2) - (U[j][i+1][2] - U[j][i][2])/(Re*dx);

			Fm[0] = (U[j][i][1] + U[j][i-1][1]) / (2*B);
			Fm[1] = pow((U[j][i][1] + U[j][i-1][1])/2, 2) + (U[j][i][0] + U[j][i-1][0])/2 - (U[j][i][1] - U[j][i-1][1])/(Re*dx);
			Fm[2] = ((U[j][i][1] + U[j][i-1][1]) / 2) * ((U[j][i][2] + U[j][i-1][2]) / 2) - (U[j][i][2] - U[j][i-1][2])/(Re*dx);

			Gp[0] = (U[j+1][i][2] + U[j][i][2]) / (2*B);
			Gp[1] = ((U[j+1][i][1] + U[j][i][1]) / 2) * ((U[j+1][i][2] + U[j][i][2]) / 2) - (U[j+1][i][1] - U[j][i][1])/(Re*dy);
			Gp[2] = pow((U[j+1][i][2] + U[j][i][2])/2, 2) + (U[j+1][i][0] + U[j][i][0])/2 - (U[j+1][i][2] - U[j][i][2])/(Re*dy);

			Gm[0] = (U[j][i][2] + U[j-1][i][2]) / (2*B);
			Gm[1] = ((U[j][i][1] + U[j-1][i][1]) / 2) * ((U[j][i][2] + U[j-1][i][2]) / 2) - (U[j][i][1] - U[j-1][i][1])/(Re*dy);
			Gm[2] = pow((U[j][i][2] + U[j-1][i][2])/2, 2) + (U[j][i][0] + U[j-1][i][0])/2 - (U[j][i][2] - U[j-1][i][2])/(Re*dy);


			for (int k = 0; k < 3; k++)
			{
				flux[j][i][k] = - (Fp[k] - Fm[k])/dx - (Gp[k] - Gm[k])/dy;
			}
		}
	}
	return flux;
}

int main()
{
	vector<vector<vector<double>>> U(jmax, vector<vector<double>>(imax, vector<double>(3)));

	init(U);
	vector<vector<vector<double>>> FI = flux(U);

	vector<vector<vector<double>>> ExFlux = exactFlux();
	

	int k = 2;	// Index for solution: 0 = pressure; 1 = u velocity; 2 = v velocity
	// printVec3D(U, k, p);
	// printVec3D(ExFlux, k, p);
	// printVec3D(FI, k, p);

	vector<vector<vector<double>>> E = error(FI, ExFlux);
	vector<double> L2 = L2Norm(E);

	printVec(L2);


	getchar();
	return 0;
}
