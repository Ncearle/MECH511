//=================================================
// Test File
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "thomas.h"
#include "vecmat_funcs.h"

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

vector<vector<double> > preconditioner(vector<vector<double> > &A, string pre)
{
  int n = A.size();
  vector<vector<double> > Pinv(n, vector<double>(n));
  if (pre == "LR")
  {
    vector<vector<double> >  L = Id(n);
    vector<vector<double> >  R(n, vector<double>(n));
    A[0][0] = -1.;
    R[0][0] = A[0][0];
    for (int i = 1; i < n; i++)
    {
      L[i][i-1] = A[i][i-1]/R[i-1][i-1];
      R[i-1][i] = A[i-1][i];
      R[i][i] = A[i][i] - L[i][i-1]*R[i-1][i];
    }
    printVec2D(L);
    printVec2D(R);

    invertLowerTri(L);
    invertUpperTri(R);
    Pinv = MM(R, L);
    return Pinv;
  }
  else if (pre == "jacobi")
  {
    for (int i = 0; i < n; i++)
    {
      Pinv[i][i] = 1/A[i][i];
    }
  }
  else if (pre == "cholesky")
  {

  }
  else
  {
    cout << "Not a recognised preconditioner, solving with none" << endl;
    Pinv = Id(n);
  }

  return Pinv;
}

// GMRES w/ Preconditioning
vector<double> GMRESPre(vector<vector<double> > &A, vector<double> &b, string pre)
{
  int n = b.size();
  vector<double> x0(n);	// Initial guess (zero vector)
	vector<double> r0 = residual(A, x0, b);
  vector<double> res = {MagV(r0)};
	vector<vector<double> > v;
  v.push_back (ScaV(1/res[0], r0));
  vector<vector<double> > Pinv = preconditioner(A, "LR");
  vector<vector<double> > z;
  z.push_back (MVM(Pinv, v[0]));

  vector<double> Be1 = {res[0]};
  vector<vector<double> > H;

  int j = 0;
  cout << "GMRES: Iteration: " << j << "\t|\tResidual: " << res[0] << endl;
  while (res[j] > tol && j < 1)
  {
    resizeMat(H, j+2, j+1);
    vector<double> w = MVM(A, z[j]);
    vector<double> wdiff(n);
    for (int i = 0; i <= j; i++)
    {
      H[i][j] = dot(w, v[i]);
      wdiff = ScaV(H[i][j], v[i]);
      w = Vsub(w, wdiff);
    }
    H[j+1][j] = MagV(w);
    printVec2D(H);
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
  while (res[j] > tol && j < 2)
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
	printVec2D(v);
  return x0;
}

vector<vector<double> > assembleLHS(vector<vector<double> > &u, vector<vector<double> > &v)
{
	int N = (imax)*(jmax);
	vector<vector<double> > LHS(N, vector<double>(N));
	vector<double> u_vec = mat2vec(u);
	vector<double> v_vec = mat2vec(v);

	double ax = dt / (Re*Pr*pow(dx, 2));
	double bx = dt / (2. * dx);
	double ay = dt / (Re*Pr*pow(dy, 2));
	double by = dt / (2. * dy);

	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		LHS[i][i] = 1;
		LHS[i][i+imax] = 1;
		LHS[N-1-i][N-1-i] = 1;
		LHS[N-1-i][N-1-imax-i] = 1;
	}

	for (int j = 1; j < jmax-1; j++)
	{
		int n = j*imax;
		// Left and Right Ghost Cells
		LHS[n][n] = 1;
		LHS[n][n+1] = 1;
		LHS[n+imax-1][n+imax-1] = 1;
		LHS[n+imax-1][n+imax-2] = 1;

		for (int i = 1; i < imax-1; i++)
		{
			// Interior Cells
			LHS[n+i][n+i] = 1. + 2.*ax + 2.*ay;
			LHS[n+i][n+i+1] = u_vec[n+i+1]*bx - ax;
			LHS[n+i][n+i-1] = u_vec[n+i-1]*bx - ax;
			LHS[n+i][n+i+imax] = v_vec[n+i+imax]*by - ay;
			LHS[n+i][n+i-imax] = v_vec[n+i-imax]*by - ay;
		}
	}

	// Corner Cells = 0
	LHS[0][0] = 1.;
	LHS[N-1][N-1] = 1.;
	LHS[imax-1][imax-1] = 1.;
	LHS[N-imax][N-imax] = 1.;

	return LHS;
}

vector<double> assembleRHS(vector<vector<double> > &FI)
{
	int N = imax*jmax;
	vector<double> RHS(N);
	vector<double> FI_vec = mat2vec(FI);

	// Top and Bottom Boundaries
	for (int i = 1; i < imax-1; i++)
	{
		RHS[i] = 0;
		RHS[N-1-i] = 0;
	}

	for (int j = 1; j < jmax-1; j++)
	{
		int n = j*imax;
		// Left and Right Boundaries
		RHS[n] = 0;
		RHS[n+imax-1] = 0;

		for (int i = 1; i < imax-1; i++)
		{
			// Interior Cells
			RHS[n+i] = FI_vec[n+i];
		}
	}

	// Corner Cells = 0
	RHS[0] = 0.;
	RHS[imax-1] = 0.;
	RHS[N-imax] = 0.;
	RHS[N-1] = 0.;
	return RHS;
}

int main()
{

	vector<vector<double> > T(jmax, vector<double>(imax));
	vector<vector<double> > u(jmax, vector<double>(imax));
	vector<vector<double> > v(jmax, vector<double>(imax));
	init(T, u, v);
	vector<vector<double> > S = source(u, v);
	vector<vector<double> > FI = FI2C(T, u, v, S);

	vector<vector<double> > A = assembleLHS(u, v);
	vector<double> b = assembleRHS(FI);
	printVec2D(A);
	printVec(b);















  // vector<vector<double> > A = {{-2, 1, 0, 0, 0, 0, 0, 0, 0, 0},
  //                             {1, -2, 1, 0, 0, 0, 0, 0, 0, 0},
  //                             {0, 1, -2, 1, 0, 0, 0, 0, 0, 0},
  //                             {0, 0, 1, -2, 1, 0, 0, 0, 0, 0},
  //                             {0, 0, 0, 1, -2, 1, 0, 0, 0, 0},
  //                             {0, 0, 0, 0, 1, -2, 1, 0, 0, 0},
  //                             {0, 0, 0, 0, 0, 1, -2, 1, 0, 0},
  //                             {0, 0, 0, 0, 0, 0, 1, -2, 1, 0},
  //                             {0, 0, 0, 0, 0, 0, 0, 1, -2, 1},
  //                             {0, 0, 0, 0, 0, 0, 0, 0, 1, -2}};
	//
  // vector<double> b = {0, 0, 0, 0, 1, 5, 1, 0, 0, 0};
	//
	// vector<double> x = GMRES(A, b);

  // vector<vector<double> > P = preconditioner(A, "LR");

  // vector<double> x = GMRESPre(A, b, "LR");
  // cout << endl;
  // printVec(x);

  // vector<vector<double> > A = {{1, 3, 0, 5},
  //                             {3, 5, 0, 8},
  //                             {0, 0, 6, 12},
  //                             {5, 8, 12, 16}};
  // printVec2D(A);
  //
  // vector<vector<double> > L(A.size(), vector<double>(A[0].size()));
  // vector<vector<double> > U(A.size(), vector<double>(A[0].size()));
  // LUdecomp(A, L, U);
  // cout << "L:\n";
  // printVec2D(L);
  // invertLowerTri(L);
  // cout << "Linv:\n";
  // printVec2D(L);
  // cout << "U:\n";
  // printVec2D(U);
  // invertUpperTri(U);
  // cout << "Uinv:\n";
  // printVec2D(U);
  //
  // vector<vector<double> > Uinv = inverseLU(A);
  // printVec2D(Uinv);
  //
  // vector<vector<double> > B = MM(U, Uinv);
  // printVec2D(B);



  // vector<vector<double> > A = {{1}, {2}};
  // printVec2D(A);
  // resizeMat(A, 3, 2);
  // printVec2D(A);

  return 0;
}
