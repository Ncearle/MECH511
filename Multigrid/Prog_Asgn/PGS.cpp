//=============================
// Source File for Multi-grid
//=============================

#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<double>> &U)
{
	int jmax = U.size();
	int imax = U[0].size();

	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		double x = (i - 0.5) / (imax - 2);
		U[0][i] = -U[1][i];
		U[jmax-1][i] = -U[jmax-2][i];

		U[0][i] = sin(pi*x)*sin(pi*(-0.5/jmax-2));
		U[jmax-1][i] = sin(pi*x)*sin(pi*((jmax-1.5)/jmax-2));
	}

	// Left and Right Ghost Cells
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		U[j][0] = -U[j][1];
		U[j][imax-1] = -U[j][imax-2];

		U[j][0] = sin(pi*y)*sin(pi*(-0.5/imax-2));
		U[j][imax-1] = sin(pi*y)*sin(pi*((imax-1.5)/imax-2));
	}
}

// Initialize the domain
void init(vector<vector<double>> &U)
{
	int jmax = U.size();
	int imax = U[0].size();

	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 1; i < imax-1; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			U[j][i] = sin(pi*x) * sin(pi*y);
		}
	}
	ghost(U);
}

// Fine-to-coarse transfer operator
vector<vector<double>> F2C(vector<vector<double>> &Uf)
{
	vector<vector<double>> Uc(Uf.size()/2+1, vector<double>(Uf[0].size()/2+1));
	for (int j = 1; j < Uc.size()-1; j++)
	{
		for (int i = 1; i < Uc[j].size()-1; i++)
		{
			Uc[j][i] = 0.25 * (Uf[2*j][2*i] + Uf[2*j-1][2*i] + Uf[2*j][2*i-1] + Uf[2*j-1][2*i-1]);
		}
	}
	ghost(Uc);
	return Uc;
}

// // Coarse-to-fine transfer operator by injection
vector<vector<double>> C2F_inj(vector<vector<double>> &Uc)
{
	vector<vector<double>> Uf(Uc.size()*2-2, vector<double>(Uc[0].size()*2-2));
	for (int j = 1; j < Uc.size()-1; j++)
	{
		for (int i = 1; i < Uc[j].size()-1; i++)
		{
			Uf[j*2][i*2] = Uc[j][i];
			Uf[j*2-1][i*2] = Uc[j][i];
			Uf[j*2][i*2-1] = Uc[j][i];
			Uf[j*2-1][i*2-1] = Uc[j][i];
		}
	}
	ghost(Uf);
	return Uf;
}
//
// // Coarse-to-fine transfer operator by interpolation
vector<vector<double>> C2F_int(vector<vector<double>> &Uc)
{
	vector<vector<double>> Uf(Uc.size()*2-2, vector<double>(Uc[0].size()*2-2));
	for (int j = 1; j < Uc.size()-1; j++)
	{
		for (int i = 1; i < Uc[j].size()-1; i++)
		{
			Uf[j*2][i*2] = 0.5625*Uc[j][i] + 0.1875*Uc[j][i+1] + 0.1875*Uc[j+1][i] + 0.0625*Uc[j+1][i+1];
			Uf[j*2-1][i*2] = 0.5625*Uc[j][i] + 0.1875*Uc[j][i+1] + 0.1875*Uc[j-1][i] + 0.0625*Uc[j-1][i+1];
			Uf[j*2][i*2-1] = 0.5625*Uc[j][i] + 0.1875*Uc[j][i-1] + 0.1875*Uc[j+1][i] + 0.0625*Uc[j+1][i-1];
			Uf[j*2-1][i*2-1] = 0.5625*Uc[j][i] + 0.1875*Uc[j][i-1] + 0.1875*Uc[j-1][i] + 0.0625*Uc[j-1][i-1];
		}
	}
	ghost(Uf);
	return Uf;
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
				double Uij = pow(dy,2)/(2*(pow(dx,2)+pow(dy,2))) * (U[i + 1][j] + U[i - 1][j]) + pow(dx, 2) / (2 * (pow(dx, 2) + pow(dy, 2)))* (U[i][j + 1] + U[i][j - 1]) - Uij0;

				U[i][j] = Uij0 + w * Uij;
				double Uij1 = U[i][j];
				delta[i][j] = abs(Uij1 - Uij0);
			}
		}
		ghost(U);
		it++;
	}
}

int main()
{
	vector<vector<double>> U(jmax, vector<double>(imax));
	vector<vector<double>> Ucoarse(jmax/2+1, vector<double>(imax/2+1));
	init(U);
	init(Ucoarse);
	ghost(U);
	ghost(Ucoarse);
	vector<vector<double>> Uc = F2C(U);
	vector<vector<double>> Uf_inj = C2F_inj(Uc);
	vector<vector<double>> Uf_int = C2F_int(Uc);

	vector<vector<double>> E_coarse = error(Uc, Ucoarse);
	vector<vector<double>> E_inject = error(Uf_inj, U);
	vector<vector<double>> E_interp = error(Uf_int, U);

	double L2_coarse = L2Norm(E_coarse);
	double L2_inject = L2Norm(E_inject);
	double L2_interp = L2Norm(E_interp);

	cout << "Transferred Coarse grid error: " << L2_coarse << endl;
	cout << "Transferred Fine Injection grid error: " << L2_inject << endl;
	cout << "Transferred Fine interpolation grid error: " << L2_interp << endl;

	string OG_name = "original.dat";
	vec2D2File(OG_name, U);
	string OGc_name = "OG_coarse.dat";
	vec2D2File(OGc_name, Ucoarse);
	string Cname = "coarse.dat";
	vec2D2File(Cname, Uc);
	string Ftname = "fine_int.dat";
	vec2D2File(Ftname, Uf_int);
	string Fjname = "fine_inj.dat";
	vec2D2File(Fjname, Uf_inj);
	string ECname = "Ecoarse.dat";
	vec2D2File(ECname, E_coarse);
	string Ejname = "Einject.dat";
	vec2D2File(Ejname, E_inject);
	string Etname = "Einterp.dat";
	vec2D2File(Etname, E_interp);
	return 0;
}
