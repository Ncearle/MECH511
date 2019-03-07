//=============================
// Source File for Multi-grid
//=============================

#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "vecmat_funcs.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<double>> &U)
{
	int jmax = U.size();
	int imax = U[0].size();
	double dx = xmax / (imax - 2);
	double dy = ymax / (jmax - 2);

	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		double x = (i - 0.5) / (imax - 2);
		U[0][i] = -U[1][i];
		U[jmax-1][i] = -U[jmax-2][i];

	}

	// Left and Right Ghost Cells
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		U[j][0] = -U[j][1];
		U[j][imax-1] = -U[j][imax-2];

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

void source(vector<vector<double>> &S)
{
	int jmax = S.size();
	int imax = S[0].size();

	for (int j = 1; j < jmax - 1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 1; i < imax - 1; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			// S[j][i] = 0;
			S[j][i] = x*(1-x)*y*(1-y);
		}
	}
}

// Fine-to-coarse transfer operator
// vector<vector<double>> F2C(vector<vector<double>> &Uf)
void F2C(vector<vector<double>> &U)
{
	vector<vector<double>> Uc(U.size()/2+1, vector<double>(U[0].size()/2+1));
	for (int j = 1; j < Uc.size()-1; j++)
	{
		for (int i = 1; i < Uc[j].size()-1; i++)
		{
			Uc[j][i] = 0.25 * (U[2*j][2*i] + U[2*j-1][2*i] + U[2*j][2*i-1] + U[2*j-1][2*i-1]);
		}
	}
	ghost(Uc);
	resizeMat(U, Uc.size(),Uc[0].size());
	AequalB(U,Uc);
	// return Uc;
}

// // Coarse-to-fine transfer operator by injection
// vector<vector<double>> C2F_inject(vector<vector<double>> &Uc)
void C2F_inject(vector<vector<double>> &U)
{
	vector<vector<double>> Uf(U.size()*2-2, vector<double>(U[0].size()*2-2));
	for (int j = 1; j < U.size()-1; j++)
	{
		for (int i = 1; i < U[j].size()-1; i++)
		{
			Uf[j*2][i*2] = U[j][i];
			Uf[j*2-1][i*2] = U[j][i];
			Uf[j*2][i*2-1] = U[j][i];
			Uf[j*2-1][i*2-1] = U[j][i];
		}
	}
	ghost(Uf);
	resizeMat(U, Uf.size(),Uf[0].size());
	AequalB(U,Uf);
	// return Uf;
}
//
// // Coarse-to-fine transfer operator by interpolation
// vector<vector<double>> C2F_interp(vector<vector<double>> &Uc)
void C2F_interp(vector<vector<double>> &U)
{
	vector<vector<double>> Uf(U.size()*2-2, vector<double>(U[0].size()*2-2));
	for (int j = 1; j < U.size()-1; j++)
	{
		for (int i = 1; i < U[j].size()-1; i++)
		{
			Uf[j*2][i*2] = 0.5625*U[j][i] + 0.1875*U[j][i+1] + 0.1875*U[j+1][i] + 0.0625*U[j+1][i+1];
			Uf[j*2-1][i*2] = 0.5625*U[j][i] + 0.1875*U[j][i+1] + 0.1875*U[j-1][i] + 0.0625*U[j-1][i+1];
			Uf[j*2][i*2-1] = 0.5625*U[j][i] + 0.1875*U[j][i-1] + 0.1875*U[j+1][i] + 0.0625*U[j+1][i-1];
			Uf[j*2-1][i*2-1] = 0.5625*U[j][i] + 0.1875*U[j][i-1] + 0.1875*U[j-1][i] + 0.0625*U[j-1][i-1];
		}
	}
	ghost(Uf);
	resizeMat(U, Uf.size(),Uf[0].size());
	AequalB(U,Uf);
	// return Uf;
}

// Point Gauss Seidel with over-relaxation
void PGS(vector<vector<double>> &U, vector<vector<double>> &S, double w)
{
	int jmax = U.size();
	int imax = U[0].size();
	double dx = xmax / (imax - 2);
	double dy = ymax / (jmax - 2);

	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			double dUji = pow(dx,2)/(2*(pow(dx,2)+pow(dy,2))) * (U[j + 1][i] + U[j - 1][i]) + pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * (U[j][i + 1] + U[j][i - 1]) - (pow(dx,2)*pow(dy,2))/(2*(pow(dx,2)+pow(dy,2)))*S[j][i] - U[j][i];
			U[j][i] += w * dUji;
		}
	}
	ghost(U);
}

vector<double> Vcycle(vector<vector<double>> &U, vector<vector<double>> &S, int nM, int nCI, vector<double> &L2)
{
	for (int n = 1; n < nM; n++)
	{
		PGS(U,S,w);
		PGS(U,S,w2);
		F2C(U);
	}
	for (int n = 0; n < nCI; n++)
	{
		PGS(U,S,w);
	}
	for (int n = 1; n < nM; n++)
	{
		// C2F_interp(U);
		C2F_inject(U);
		PGS(U,S,w);
		PGS(U,S,w2);
	}

	// W-cycle?
	// Second pass
	// double w2 = 0.8724;
	// for (int n = 1; n < nM; n++)
	// {
	// 	PGS(U,S,w2);
	// 	F2C(U);
	// }
	// for (int n = 1; n < nCI; n++)
	// {
	// 	PGS(U,S,w2);
	// }
	// for (int n = 1; n < nM; n++)
	// {
	// 	PGS(U,S,w2);
	// 	C2F_interp(U);
	// 	// C2F_inject(U);
	// }
	// PGS(U,S,w2);

	// L2.push_back (L2Norm(U));
	// return L2;
}

int main()
{
	vector<vector<double>> U(jmax, vector<double>(imax));
	vector<vector<double>> S(jmax, vector<double>(imax));
	init(U);
	source(S);

	int it = 0;
	vector<double> L2 (1, 1.0);
	int nM = 4;	// Number of mesh levels
	int nC = 2; // Number of coarsest mesh passes

	while (L2[it] > tol)
	{
	  it++;
		vector<vector<double>> Uprev(jmax, vector<double>(imax));
		AequalB(Uprev, U);

		Vcycle(U, S, nM, nC, L2);

		vector<vector<double>> delta = Msub(U, Uprev);
		L2.push_back (L2Norm(delta));
	}

	cout << "Iterations: " << it << endl;
	cout << "L2 Norm: " << L2[it] << endl;
	// string L2name = "L2_doublepass.dat";
	// vec1D2File(L2name, L2);

	return 0;
}





// ==========================
// Part 1 Code
// ==========================
// vector<vector<double>> Ucoarse(jmax/2+1, vector<double>(imax/2+1));
// init(Ucoarse);
// ghost(Ucoarse);
// vector<vector<double>> Uc = F2C(U);
// vector<vector<double>> Uf_inj = C2F_inject(Uc);
// vector<vector<double>> Uf_int = C2F_interp(Uc);
//
// vector<vector<double>> E_coarse = error(Uc, Ucoarse);
// vector<vector<double>> E_inject = error(Uf_inj, U);
// vector<vector<double>> E_interp = error(Uf_int, U);
//
// double L2_coarse = L2Norm(E_coarse);
// double L2_inject = L2Norm(E_inject);
// double L2_interp = L2Norm(E_interp);
//
// cout << "Transferred Coarse grid error: " << L2_coarse << endl;
// cout << "Transferred Fine Injection grid error: " << L2_inject << endl;
// cout << "Transferred Fine interpolation grid error: " << L2_interp << endl;

// ----------------------------
// Printing Grids to data files
// ----------------------------
// string OG_name = "data\\original.dat";
// vec2D2File(OG_name, U);
// string OGc_name = "data\\OG_coarse.dat";
// vec2D2File(OGc_name, Ucoarse);
// string Cname = "data\\coarse.dat";
// vec2D2File(Cname, Uc);
// string Ftname = "data\\fine_int.dat";
// vec2D2File(Ftname, Uf_int);
// string Fjname = "data\\fine_inj.dat";
// vec2D2File(Fjname, Uf_inj);
// string ECname = "data\\Ecoarse.dat";
// vec2D2File(ECname, E_coarse);
// string Ejname = "data\\Einject.dat";
// vec2D2File(Ejname, E_inject);
// string Etname = "data\\Einterp.dat";
// vec2D2File(Etname, E_interp);

// =================================
// Part 2 Code
// =================================
// while (L2 > tol)
// {
// 	for (int n = 0; n < 5; n++)
// 	{
// 		it++;
// 		PGS(U,S);
// 		F2C(U);
// 		for (int n = 0; n < 20; n++)
// 		{
// 			// it++;
// 			PGS(U,S);
// 		}
// 		// C2F_inject(U);
// 		C2F_interp(U);
// 		PGS(U,S);
// 		L2 = L2Norm(U);
// 	}
// }
