//=================================================
// Test File
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "thomas.h"
#include "vecmat_funcs.h"

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


int main()
{
  vector<vector<double> > A = {{-2, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                              {1, -2, 1, 0, 0, 0, 0, 0, 0, 0},
                              {0, 1, -2, 1, 0, 0, 0, 0, 0, 0},
                              {0, 0, 1, -2, 1, 0, 0, 0, 0, 0},
                              {0, 0, 0, 1, -2, 1, 0, 0, 0, 0},
                              {0, 0, 0, 0, 1, -2, 1, 0, 0, 0},
                              {0, 0, 0, 0, 0, 1, -2, 1, 0, 0},
                              {0, 0, 0, 0, 0, 0, 1, -2, 1, 0},
                              {0, 0, 0, 0, 0, 0, 0, 1, -2, 1},
                              {0, 0, 0, 0, 0, 0, 0, 0, 1, -2}};

  vector<double> b = {0, 0, 0, 0, 1, 5, 1, 0, 0, 0};

  // vector<vector<double> > P = preconditioner(A, "LR");

  vector<double> x = GMRESPre(A, b, "LR");
  cout << endl;
  printVec(x);

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
