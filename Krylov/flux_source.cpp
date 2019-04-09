//====================================================
// Flux and Source term verification for Parts 1 and 2
//====================================================
#include "constant.h"
#include "print_fcns.h"
#include "error_fcns.h"
#include "exact.h"

void init(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			T[j][i] = T0 * cos(pi*x) * sin(pi*y);	
			u[j][i] = u0 * y * sin(pi*x);
			v[j][i] = v0 * x * cos(pi*y);
		}
	}
}

vector<vector<double>> source(vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> S(jmax, vector<double>(imax));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{

			S[j][i] = Ec / Re * (2 * pow((u[j][i + 1] - u[j][i - 1]) / (2 * dx), 2) + 2 * pow((v[j + 1][i] - v[j - 1][i]) / (2 * dy), 2) + pow((v[j][i + 1] - v[j][i - 1]) / (2 * dx) + (u[j + 1][i] - u[j - 1][i]) / (2 * dy), 2));
		}
	}
	return S;
}

vector<vector<double>> FI2C(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> FI(jmax, vector<double>(imax));
	vector<vector<double>> con(jmax, vector<double>(imax));
	vector<vector<double>> dif(jmax, vector<double>(imax));

	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			// Convective term
			con[j][i] = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - v[j - 1][i] * T[j - 1][i]) / (2 * dx);

			// Diffusive term
			dif[j][i] = ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2)) / (Re * Pr);

			FI[j][i] = -(con[j][i] + dif[j][i]);
		}
	}
	return FI;
}

vector<vector<double>> FI2C_2(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> FI(jmax, vector<double>(imax));
	vector<vector<double>> con(jmax, vector<double>(imax));
	vector<vector<double>> dif(jmax, vector<double>(imax));

	vector<vector<double>> S = source(u,v);

	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			// Convective term
			con[j][i] = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - v[j - 1][i] * T[j - 1][i]) / (2 * dx);

			// Diffusive term
			dif[j][i] = ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2)) / (Re * Pr);

			FI[j][i] = -(con[j][i] + dif[j][i]) + S[j][i];
		}
	}
	return FI;
}

int main()
{
	vector<vector<double>> T(jmax, vector<double>(imax));
	vector<vector<double>> u(jmax, vector<double>(imax));
	vector<vector<double>> v(jmax, vector<double>(imax));

	init(T, u, v);
	// printVec2D(T);
	
	vector<vector<double>> S = source(u, v);
	vector<vector<double>> ExS = exactSource();
	vector<vector<double>> FI = FI2C(T, u, v);
	vector<vector<double>> ExFI = exactFlux();
	
	vector<vector<double>> FE = error(FI, ExFI);
	vector<vector<double>> SE = error(S, ExS);

	vector<vector<double>> FI_2 = FI2C_2(T, u, v);	

	vector<vector<double>> TEx(jmax, vector<double>(imax));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			T[j][i] = FI[j][i] + S[j][i];
			TEx[j][i] = ExFI[j][i] + ExS[j][i];
		}
	}
	
	printVec2D(S);

	//printVec2D(FI_2);
	//printVec2D(T);

	//printVec2D(TEx);

	vector<vector<double>> TE = error(T, TEx);

	double L2TE = L2Norm(TE);
	// double L2SE = L2Norm(SE);

	cout << "Flux L2: " << L2TE << endl;
	// cout << "Source L2: " << L2SE << endl;

	getchar();
	return 0;
}
