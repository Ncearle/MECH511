//=================================================
// Source file for Energy Equation
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "thomas.h"
#include "vecmat_funcs.h"
// #include "xtensor/xarray.hpp"
// #include "xtensor/xio.hpp"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v)
{
	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		T[0][i] = -T[1][i];
		T[jmax-1][i] = 2. - T[jmax-2][i];
	}

	// Left and Right Ghost Cells
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		T[j][0] = y;
		T[j][imax-1] = T[j][imax-2];
	}
}

// Calculates the source term, given velocities
vector<vector<double> > source(vector<vector<double> > &u, vector<vector<double> > &v)
{
	vector<vector<double> > S(jmax, vector<double>(imax));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			S[j][i] = Ec / Re * (2 * pow((u[j][i + 1] - u[j][i - 1]) / (2 * dx), 2)
								+ 2 * pow((v[j + 1][i] - v[j - 1][i]) / (2 * dy), 2)
								+ pow((v[j][i + 1] - v[j][i - 1]) / (2 * dx) + (u[j + 1][i] - u[j - 1][i]) / (2 * dy), 2));
		}
	}
	return S;
}

// Initializes the domain and velocities
void init(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v)
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			T[j][i] = y;
			u[j][i] = 6.0*ubar*y*(1-y);
			v[j][i] = 0.;
		}
	}
	ghost(T, u, v);
}

// Calculates the gradient along the bottom wall
vector<double> grad(vector<vector<double> > &T)//, double y)
{
	vector<double> grad(imax-2);
	// int j = y * (jmax-2);
	for (int i = 1; i < imax-1; i++)
	{
		grad[i-1] = (T[1][i] - T[0][i])/dy;
	}
	return grad;
}

// Calculates the flux integral given the domain, velocities, and source term
vector<vector<double> > FI2C(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &S)
{
	vector<vector<double> > FI(jmax, vector<double>(imax));
	vector<vector<double> > con(jmax, vector<double>(imax));
	vector<vector<double> > dif(jmax, vector<double>(imax));

	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			// Convective term
			con[j][i] = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - v[j - 1][i] * T[j - 1][i]) / (2 * dy);

			// Diffusive term
			dif[j][i] = ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2)) / (Re * Pr);

			FI[j][i] = (con[j][i] + dif[j][i]) + S[j][i];
		}
	}
	return FI;
}

vector<double> forSub(vector<vector<double> > &L, vector<double> &b)
{
	int n = L.size();
	vector<double> y(n);

  for(int i = 0; i < n; i++)
	{
    y[i] = b[i];
    for(int j = 0; j < i; j++)
		{
        y[i] -= L[i][j] * y[j];
    }
    y[i] /= L[i][i];
   }
	return y;
}

vector<double> backSub(vector<vector<double> > &U, vector<double> &y)
{
	int n = U.size();
	vector<double> x(n);

	for (int i = n-1; i >=0; i--)
	{
		x[i] = y[i];
		for(int j = i+1; j < n; j++)
		{
			x[i] -= U[i][j]*x[j];
		}
		x[i] /= U[i][i];
	}
	return x;
}

//=================
// Explicit Schemes
//=================
// Explicit Euler Time advance
void EE(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &S)
{
	vector<vector<double> > T0(jmax, vector<double>(imax));
	int it = 0;
	double delta = 1.0;
	while (abs(delta) > tol)
	{
		it++;
		vector<vector<double> > FI = FI2C(T, u, v, S);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T0[j][i] = T[j][i];
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
		ghost(T, u, v);

		delta = maxChange(T0, T);
	}
	// printVec2D(T);
	//
	vector<vector<double> > ExT = exactTemp();
	vector<vector<double> > err = error(T, ExT);
	double L2 = L2Norm(err);

	cout << "Domain: " << imax-2 << " X " << jmax-2 << endl;
	cout << setprecision(6) << "L2 norm: " << L2 << endl;
	cout << "Iterations: " << it << endl;
	cout << "Timestep: " << dt << endl;
}

// Two Stage Runge Kutta time advance
void RK2(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &S)
{
	vector<vector<double> > T0(jmax, vector<double>(imax));
	int it = 0;
	double delta = 1.0;
	while (abs(delta) > tol)
	{
		it++;
		// Intermediate Step
		vector<vector<double> > FIint = FI2C(T, u, v, S);
		vector<vector<double> > Tint(jmax, vector<double>(imax));
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T0[j][i] = T[j][i];
				Tint[j][i] = T[j][i] + dt/2.0 * FIint[j][i];
			}
		}
		ghost(Tint, u, v);

		// Full Step
		vector<vector<double> > FI = FI2C(Tint, u, v, S);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
		ghost(T, u, v);

		delta = maxChange(T0, T);
	}
	// printVec2D(T);
	cout << setprecision(6) << delta << endl;
	vector<vector<double> > ExT = exactTemp();
	vector<vector<double> > err = error(T, ExT);
	double L2 = L2Norm(err);
	cout << setprecision(6) << "L2 norm: " << L2 << endl;
	cout << "Iterations: " << it << endl;
	cout << "Timestep: " << dt << endl;
	cout << "Tolerance: " << tol << endl;
}

//=================
// Implicit Schemes
//=================
vector<vector<double> > assembleLHS(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v)
{
	int dim = (imax-2)*(jmax-2);
	vector<vector<double> > LHS(dim, vector<double>(dim));
	vector<vector<double> > u_nogho = removeGhost(u);
	vector<vector<double> > v_nogho = removeGhost(v);
	vector<double> u_vec = mat2vec(u_nogho);
	vector<double> v_vec = mat2vec(v_nogho);

	double ax = dt / (Re*Pr*pow(dx, 2));
	double bx = dt / (2. * dx);
	double ay = dt / (Re*Pr*pow(dy, 2));
	double by = dt / (2. * dy);

	for (int i = 1; i < dim-1; i++)
	{
		LHS[i][i] = 1. + 2.*ax + 2.*ay;
		LHS[i][i+1] = u_vec[i+1]*bx - ax;
		LHS[i][i-1] = u_vec[i-1]*bx - ax;
		if(i % (imax-2) == 0 )
		{
			LHS[i-1][i] = 0.;
			LHS[i][i-1] = 0.;
		}
	}
	LHS[0][0] = 1. + 2.*ax + 2.*ay;
	LHS[dim-1][dim-1] = 1. + 2.*ax + 2.*ay;
	LHS[0][1] = u_vec[1]*bx - ax;
	LHS[dim-1][dim-2] = u_vec[dim-1]*bx - ax;

	for (int i = 0; i < LHS.size()-(imax-2); i++)
	{
		LHS[i][(imax-2)+i] = v_vec[imax-2+i]*by - ay;
		LHS[(imax-2)+i][i] = -v_vec[i]*by - ay;
	}

	return LHS;
}
// LU Factorisation
vector<double> LU(vector<vector<double> > &A, vector<double> &b)
{
	int n = A.size();
	vector<vector<double> > L = Id(n);
	vector<vector<double> > U = copyMat(A);

	for (int i = 0; i < n-1; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			L[j][i] = U[j][i]/U[i][i];
			for (int k = i; k < n; k++)
			{
				U[j][k] -= L[j][i]*U[i][k];
			}
		}
	}
	vector<double> y = forSub(L, b);
	vector<double> x = backSub(U, y);
	return x;
}
// Cholesky Factorisation
vector<double> Cholesky(vector<vector<double> > &A, vector<double> &b)
{
	int n = A.size();
	vector<vector<double> > L(n, vector<double>(n));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < (i+1); j++)
		{
      double s = 0;
      for (int k = 0; k < j; k++)
			{
          s += L[i][k] * L[j][k];
			}
			if(i == j)
				L[i][j] = sqrt(A[i][i] - s);
			else
				L[i][j] = (1.0 / L[j][j] * (A[i][j] - s));
		}
	}

	vector<double> y = forSub(L, b);
	vector<vector<double> > Lstar = transpose(L);
	vector<double> x = backSub(Lstar, y);
	return x;
}
// LU Factorisation with Partial Pivoting
vector<double> LUpivot(vector<vector<double> > &A, vector<double> &b)
{
	int n = A.size();
	vector<vector<double> > L = Id(n);
	vector<vector<double> > P = Id(n);
	vector<vector<double> > U = copyMat(A);

	for (int i = 0; i < n-1; i++)
	{
		int q = i;
		for (int k = i; k < n-1; k++)
		{
			if (abs(U[k][i] > U[i][i])) q = k;
		}
		swaprow(U, i, q);
		swaprow(L, i, q);
		swaprow(P, i, q);

		for (int j = i+1; j < n; j++)
		{
			L[j][i] = U[j][i]/U[i][i];
			for (int k = i; k < n; k++)
			{
				U[j][k] -= L[j][i]*U[i][k];
			}
		}
	}
	vector<double> y = forSub(L, b);
	vector<double> x = backSub(U, y);
	return x;
}
// Approximate Factorisation
vector<vector<double> > AppFac(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &S, vector<vector<double> >FI)
{
		// Matrix {[I] + dt*[Dx]}
		vector<vector<double> > DX(imax, vector<double>(3));
		vector<double> FIx(imax);
		vector<vector<double> > Ttilda(jmax);
		double ax = dt / (Re*Pr*pow(dx, 2));
		double bx = dt / (2 * dx);

		for (int j = 1; j < jmax - 1; j++)
		{
			for (int i = 1; i < imax - 1; i++)
			{
				DX[i][0] = -bx * u[j][i-1] - ax;
				DX[i][1] = 1 + 2 * ax;
				DX[i][2] = bx * u[j][i+1] - ax;
				FIx[i] = dt * FI[j][i];
			}
			DX[0][0] = 1;
			DX[0][1] = 1;
			DX[imax-1][1] = 1;
			DX[imax-1][2] = -1;

			FIx[0] = 0;
			FIx[imax-1] = 0;

			SolveThomas(DX, FIx, imax);
			Ttilda[j] = FIx;
		}

		// Matrix {[I] + dt*[Dy]}
		vector<vector<double> > DY(jmax, vector<double>(3));
		vector<double> Tty(jmax);		// Ttilda in vector form for each column
		vector<vector<double> > deltaT(imax);
		double ay = dt / (Re*Pr*pow(dy, 2));
		double by = dt / (2 * dy);

		for (int i = 1; i < imax -1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				DY[j][0] = -by * v[j-1][i] - ay;
				DY[j][1] = 1 + 2 * ay;
				DY[j][2] = by * v[j+1][i] - ay;
				Tty[j] = Ttilda[j][i];
			}
			DY[0][0] = 1;
			DY[0][1] = 1;
			DY[jmax-1][1] = 1;
			DY[jmax-1][2] = -1;

			Tty[0] = 0;
			Tty[jmax-1] = 0;

			SolveThomas(DY, Tty, jmax);
			deltaT[i] = Tty;
		}
	return deltaT;
}
// Least Squares solve
vector<double> leastSquares(vector<vector<double> > &A, vector<double> &b)
{
  int n = b.size();
  vector<vector<double> > Atran = transpose(A);
  vector<vector<double> > Asquare = MM(Atran, A);
  vector<double> bsquare = MVM(Atran, b);
  vector<double> x = LU(Asquare, bsquare);
  return x;
}
// GMRES
vector<double> GMRES(vector<vector<double> > &A, vector<double> &b)
{
  int n = b.size();
  vector<double> x0(n);	// Initial guess (zero vector)
	vector<double> r0 = residual(A, x0, b);
  vector<double> res = {MagV(r0)};
	vector<vector<double> > v;
  v.push_back (ScaV(1/res[0], r0));

  vector<double> Be1 = {res[0]};
  vector<vector<double> > H;

  int j = 0;
  cout << "GMRES: Iteration: " << j << "\t|\tResidual: " << res[0] << endl;
  while (res[j] > tol && j < n)
  {
    resizeMat(H, j+2, j+1);
    vector<double> w = MVM(A, v[j]);
    vector<double> wdiff(n);
    for (int i = 0; i <= j; i++)
    {
      H[i][j] = dot(w, v[i]);
      wdiff = ScaV(H[i][j], v[i]);
      w = Vsub(w, wdiff);
    }
    H[j+1][j] = MagV(w);
    Be1.push_back (0.);
    vector<double> y = leastSquares(H, Be1);
    v = transpose(v);
    vector<double> u = MVM(v, y);
    v = transpose(v);
    v.push_back (ScaV(1/H[j+1][j], w));

    zeroVec(x0);
    x0 = Vadd(x0, u);
    r0 = residual(A, x0, b);
    res.push_back (MagV(r0));
    j++;
    cout << "GMRES: Iteration: " << j << "\t|\tResidual: " << res[j] << endl;
  }
  return x0;
}
int main()
{
	int time=clock();
	vector<vector<double> > T(jmax, vector<double>(imax));
	vector<vector<double> > u(jmax, vector<double>(imax));
	vector<vector<double> > v(jmax, vector<double>(imax));
	init(T, u, v);
	vector<vector<double> > S = source(u,v);	// Source Term
	vector<vector<double> > T0 = copyMat(T);	// Initial
	vector<vector<double> > A = assembleLHS(T,u,v);
	vector<vector<double> > FI = FI2C(T, u, v, S);	// Only changes on every time loop
	vector<vector<double> > FI_nogho = removeGhost(FI);
	vector<double> b = mat2vec(FI_nogho);


	vector<vector<double> > x = AppFac(T, u, v, S, FI);
	// vector<vector<double> > sol = vec2mat(x);
	printVec2D(x);

	time = clock() - time;
	cout << "Domain size: " << (imax-2) << " X " << (jmax-2) << endl;
	cout << "time [sec]: " << time/double(CLOCKS_PER_SEC) << endl;

	string LHS = "LHS.dat";
	vec2D2File(LHS,A);

	string RHS = "RHS.dat";
	vec1D2File(RHS, b);




	// int it = 0;
	// double delta = 1.0;
	// cout << "Iteration: " << it << "\t|\tDelta: " << delta << endl;
	// while (abs(delta) > tol)
	// {
	// 	it++;
	// 	vector<vector<double> > FI = FI2C(T, u, v, S);	// Only changes on every time loop
	//
	// 	// vector<vector<double> >deltaT = AppFac(T, u, v, S, FI);
	// 	// // RK2(T, u, v);
	// 	// // EE(T, u, v);
	// 	// for (int j = 1; j < jmax-1; j++)
	// 	// {
	// 	// 	for (int i = 1; i < imax-1; i++)
	// 	// 	{
	// 	// 		T[j][i] += deltaT[i][j];
	// 	// 	}
	// 	// }
	//
	//
	// 	vector<vector<double> > FI_nogho = removeGhost(FI);
	// 	vector<double> b = mat2vec(FI_nogho);
	// 	b = ScaV(dt, b);
	// 	vector<double> x = LU(A, b);
	// 	// vector<double> x = Cholesky(A, b);
	// 	// vector<double> x = LUpivot(A, b);
	//
	// 	vector<vector<double> >deltaT = vec2mat(x);
	// 	for (int j = 0; j < jmax-2; j++)
	// 	{
	// 		for (int i = 0; i < imax-2; i++)
	// 		{
	// 			T[j+1][i+1] += deltaT[j][i];
	// 		}
	// 	}
	//
	// 	ghost(T, u, v);
	//
	// 	delta = maxChange(T0, T);
	// 	cout << "Iteration: " << it << "\t|\tDelta: " << delta << endl;
	// 	T0 = copyMat(T);
	// }
	//
	// vector<vector<double> > ExT = exactTemp();
	// cout << setprecision(6) << delta << endl;
	// vector<vector<double> > err = error(T, ExT);
	// double L2 = L2Norm(err);
	// cout << setprecision(6) << "L2 norm: " << L2 << endl;
	// cout << "Iterations: " << it << endl;
	// cout << "Timestep: " << dt << endl;
	//
	// time = clock() - time;
	// cout << "Domain size: " << (imax-2) << " X " << (jmax-2) << endl;
	// cout << "time [sec]: " << time/double(CLOCKS_PER_SEC) << endl;


	// string Tname = "T_dev_" + to_string(imax-2) + "x" + to_string(jmax-2) + ".dat";
	// cout << Tname << endl;
	// vec2D2File(Tname,T);
	//
	// int x = xmax;
	// int y = ymax;
	// string Gname = "Grad_x" + to_string(x) + "_w" + to_string(y) + ".dat";
	// cout << Gname << endl;
	// vector<double> G = grad(T);
	// vec1D2File(Gname, G);

	// printVec(G);

	return 0;
}
