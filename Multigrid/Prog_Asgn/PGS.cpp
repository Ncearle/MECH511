//=============================
// Source File for Multi-grid
//=============================

#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<double>> &U)
{
	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		double x = (i - 0.5) / (imax - 2);
		U[0][i] = U[1][i];
		U[jmax-1][i] = U[jmax-2][i];
	}

	// Left and Right Ghost Cells
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		U[j][0] = U[j][1];		
		U[j][imax-1] = U[j][imax-2];
	}
}

// Initialize the domain
void init(vector<vector<double>> &U)
{
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 1; i < imax-1; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			U[j][i][0] = cos(pi*x) * cos(pi * y);
		}
	}
	ghost(U);
}

// Point Gauss Seidel with over-relaxation
void PGS(vector<vector<double>> &U){
	vector<vector<double>> delta(jmax, vector<double>(imax));
	double maxDelta = 1.0;
	int it = 1;

	while (maxDelta > tol) {
		for (int i = 1; i < imax - 1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				double Uij0 = U[i][j];
				double Uij = pow(dy,2)/(2*(pow(dx,2)+pow(dy,2))) * (domain[i + 1][j] + domain[i - 1][j]) + pow(dx, 2) / (2 * (pow(dx, 2) + pow(dy, 2)))* (domain[i][j + 1] + domain[i][j - 1]) - (pow(dy,2) * pow(dx, 2) / (2 * (pow(dx, 2) + pow(dy, 2)))) * source[i][j] - domain_ij0;

				U[i][j] = Uij0 + w * Uij;
				double Uij1 = U[i][j];
				delta[i][j] = abs(Uij1 - Uij0);
			}
		}
		maxDelta = maxChange(delta);
		vec.push_back(maxDelta);
		//setBoundaryT(domain);
		setBoundaryP(domain);
		it++;
	}
}

int main()
{
	vector<vector<double>> U(jmax, vector<double>(imax));
} 
