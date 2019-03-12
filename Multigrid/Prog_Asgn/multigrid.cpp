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

	U[0][0] = -U[1][0];
	U[jmax-1][0] = -U[jmax-1][1];
	U[0][imax-1] = -U[1][imax-1];
	U[jmax-1][imax-1] = -U[jmax-2][imax-1];

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

void source(vector<vector<double> > &S)
{
	int jmax = S.size();
	int imax = S[0].size();

	for (int j = 1; j < jmax - 1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 1; i < imax - 1; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			S[j][i] = 0;
			// S[j][i] = x*(1-x)*y*(1-y);
		}
	}
}

void zero(vector<vector<double> > &U)
{
	for (int j = 0; j < U.size(); j++)
	{
		for (int i = 0; i < U[0].size(); i++)
		{
			U[j][i] = 0;
		}
	}
}

// Fine-to-coarse transfer operator
// vector<vector<double>> F2C_restrict(vector<vector<double>> &U)
void F2C_restrict(vector<vector<double>> &U)
{
	vector<vector<double>> Uc(U.size()/2+1, vector<double>(U[0].size()/2+1));
	for (int j = 1; j < Uc.size()-1; j++)
	{
		for (int i = 1; i < Uc[j].size()-1; i++)
		{
			Uc[j][i] = 0.25 * (U[2*j][2*i] + U[2*j-1][2*i] + U[2*j][2*i-1] + U[2*j-1][2*i-1]);
		}
	}
	// ghost(Uc);
	resizeMat(U, Uc.size(),Uc[0].size());
	AequalB(U,Uc);
	// return Uc;
}

// // Coarse-to-fine transfer operator by injection
// vector<vector<double>> C2F_inject(vector<vector<double>> &U)
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
// vector<vector<double>> C2F_interp(vector<vector<double>> &U)
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
void PGSpoisson(vector<vector<double>> &U, vector<vector<double>> &S, double w)
{
	int jmax = U.size();
	int imax = U[0].size();
	double dx = xmax / (imax - 2);
	double dy = ymax / (jmax - 2);

	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			double dUji = pow(dx,2)/(2*(pow(dx,2)+pow(dy,2))) * (U[j + 1][i] + U[j - 1][i])
				+ pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * (U[j][i + 1] + U[j][i - 1])
				- (pow(dx,2)*pow(dy,2))/(2*(pow(dx,2)+pow(dy,2)))*S[j][i] - U[j][i];
			U[j][i] += w * dUji;
		}
	}
	ghost(U);
}

vector<vector<double> > residual(vector<vector<double> > &U, vector<vector<double> > &S)
{
	int jmax = U.size();
	int imax = U[0].size();
	double dx = xmax / (imax - 2);
	double dy = ymax / (jmax - 2);
	vector<vector<double> > R(U.size(), vector<double>(U[0].size()));	// Residual
	for (int j = 1; j < R.size()-1; j++)
	{
		for (int i = 1; i < R[0].size()-1; i++)
		{
			R[j][i] = S[j][i] - (U[j-1][i] - 2*U[j][i] + U[j+1][i])/pow(dy,2) - (U[j][i-1] - 2*U[j][i] + U[j][i+1])/pow(dx,2);
		}
	}
	return R;
}

vector<double> Vcycle(vector<vector<double> > &U, vector<vector<double> > &S, int NMeshes, int NCoarseIts, vector<double> &L2)
{
	PGSpoisson(U,S,w);
	vector<vector<vector<double> > > Utilda(NMeshes, vector<vector<double> >(jmax, vector<double>(imax)));
	vector<vector<vector<double> > > R(NMeshes, vector<vector<double> >(jmax, vector<double>(imax)));
	BequalA(U, Utilda[0]);
	R[0] = residual(Utilda[0],S);

	for (int n = 1; n < NMeshes; n++)
	{
		vector<vector<double> > Rint((jmax-2)/pow(2,n-1)+2, vector<double>((imax-2)/pow(2,n-1)+2));
		vector<vector<double> > Uint((jmax-2)/pow(2,n)+2, vector<double>((imax-2)/pow(2,n)+2));
		zero(Uint);
		AequalB(Rint, R[n-1]);
		F2C_restrict(Rint);

		PGSpoisson(Uint,Rint,w);
		PGSpoisson(Uint,Rint,w2);
		BequalA(Uint, Utilda[n]);
		Rint = residual(Uint, Rint);
		BequalA(Rint, R[n]);
	}

	vector<vector<double> > Rcoarse((jmax-2)/pow(2,NMeshes-2)+2, vector<double>((imax-2)/pow(2,NMeshes-2)+2));
	AequalB(Rcoarse, R[NMeshes-2]);
	F2C_restrict(Rcoarse);
	BequalA(Rcoarse,R[NMeshes-1]);
	vector<vector<double> > Ucoarse((jmax-2)/pow(2,NMeshes-1)+2, vector<double>((imax-2)/pow(2,NMeshes-1)+2));
	AequalB(Ucoarse, Utilda[NMeshes-1]);
	for (int n = 1; n < NCoarseIts; n++)
	{
		PGSpoisson(Ucoarse,Rcoarse,w);
		PGSpoisson(Ucoarse,Rcoarse,w2);
	}
	BequalA(Ucoarse, Utilda[NMeshes-1]);

	AequalB(R[0], S);

	for (int n = NMeshes-1; n > 0; n--)
	{
		vector<vector<double> > Uint((jmax-2)/pow(2,n)+2, vector<double>((imax-2)/pow(2,n)+2));
		AequalB(Uint, Utilda[n]);
		C2F_interp(Uint);
		MatAddAtoB(Utilda[n-1], Uint);
		vector<vector<double> > Rint((jmax-2)/pow(2,n-1)+2, vector<double>((imax-2)/pow(2,n-1)+2));
		AequalB(Rint, R[n-1]);
		PGSpoisson(Uint, Rint,w);
		PGSpoisson(Uint, Rint,w2);
		BequalA(Uint, Utilda[n-1]);
	}
	AequalB(U,Utilda[0]);

	L2.push_back (L2Norm(U));
	return L2;
}
//------------------------------------------------------------------------------
int main()
{
	vector<vector<double>> U(jmax, vector<double>(imax));
	vector<vector<double>> S(jmax, vector<double>(imax));
	init(U);
	source(S);

	int it = 0;
	vector<double> L2 (1, 1.0);
	int NMeshes = 4;	// Number of mesh levels
	int NCoarseIts = 2; // Number of passes on coarsest mesh

	while (L2[it] > tol)
	{
			it++;
			Vcycle(U,S,NMeshes,NCoarseIts,L2);
	}

	// while (L2[it] > tol)
	// {
	// 	it++;
	// 	PGSpoisson(U, S, w);
	// 	L2.push_back (L2Norm(U));
	// }

	// vector<vector<double>> Uprev(jmax, vector<double>(imax));
	// while (L2[it] > tol)
	// {
	//   it++;
	// 	AequalB(Uprev, U);
	// 	Vcycle(U, S, NMeshes, NCoarseIts, L2);
	//
	// 	vector<vector<double>> delta = Msub(U, Uprev);
	// 	L2.push_back (L2Norm(delta));
	// }


	cout << "Iterations: " << it << endl;
	cout << "L2 Norm: " << L2[it] << endl;
	string L2name = "L2_force_256.dat";
	vec1D2File(L2name, L2);

	return 0;
}


// ==========================
// Part 1 Code
// ==========================
// vector<vector<double>> Ucoarse(jmax/2+1, vector<double>(imax/2+1));
// init(Ucoarse);
// ghost(Ucoarse);
// vector<vector<double>> Uc = F2C_restrict(U);
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


//==============================================================================
// Code that I used to convince myself how to actually do this.
//==============================================================================
// while ((L2[it]) > tol)
// // while(it < 2)
// {
// 	it++;
// 	PGSpoisson(U,S,w);
// 	// Restriction
// 	//2nd level
// 	vector<vector<double> > Rh = residual(U,S);
// 	vector<vector<double> > R2h = F2C_restrict(Rh);
// 	vector<vector<double> > U2h((jmax-2)/2+2, vector<double>((imax-2)/2+2));
// 	// Start of third level
// 	PGSpoisson(U2h, R2h, w);
// 	R2h = residual(U2h, R2h);
// 	vector<vector<double> > R4h = F2C_restrict(R2h);
// 	vector<vector<double> > U4h((jmax-2)/4+2, vector<double>((imax-2)/4+2));
// 	// // Start of fourth level
// 	// PGSpoisson(U4h, R4h, w);
// 	// R4h = residual(U4h, R4h);
// 	// vector<vector<double> > R8h = F2C_restrict(R4h);
// 	// vector<vector<double> > U8h((jmax-2)/8+2, vector<double>((imax-2)/8+2));
// 	// // // Start of fifth level
// 	// PGSpoisson(U8h, R8h, w);
// 	// R8h = residual(U8h, R8h);
// 	// vector<vector<double> > R16h = F2C_restrict(R8h);
// 	// vector<vector<double> > U16h((jmax-2)/16+2, vector<double>((imax-2)/16+2));
//
// 	// Coarsest Mesh
// 	for (int n = 0; n < 2; n++)
// 	{
// 		PGSpoisson(U4h, R4h, w);
// 	}
//
// 	//Prolongation
// 	// vector<vector<double> > e8h = C2F_interp(U16h);
// 	// MatAddBtoA(U8h, e8h);
// 	// PGSpoisson(U8h, R8h, w);
// 	// vector<vector<double> > e4h = C2F_interp(U8h);
// 	// MatAddBtoA(U4h, e4h);
// 	// PGSpoisson(U4h, R4h, w);
// 	vector<vector<double> > e2h = C2F_inject(U4h);
// 	MatAddBtoA(U2h, e2h);
// 	PGSpoisson(U2h, R2h, w);
// 	vector<vector<double> > eh = C2F_inject(U2h);
// 	MatAddBtoA(U, eh);
// 	PGSpoisson(U, S, w);
//
// 	L2.push_back (L2Norm(U));
// 	// cout << "L2 Norm: " << L2[it] << endl;
// }


// ============================================================================
