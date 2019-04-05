//=================================================
// Test File
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "thomas.h"
#include "vecmat_funcs.h"

// GMRES
vector<double> GMRES(vector<vector<double> > &A, vector<double> &b)
{
  int n = b.size();
  vector<double> x0(n);	// Initial guess (zero vector)
	vector<double> r0 = residual(A, x0, b);
	double mag_r0 = MagV(r0);
	vector<vector<double> > v;
  v.push_back (ScaV(1/mag_r0, r0));

  vector<double> Be1 = {mag_r0};
  vector<vector<double> > H;

  int j = 0;
  while (mag_r0/n > tol && j < 2)
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
    v.push_back (ScaV(1/H[j+1][j], w));



    Be1.push_back (0.);
    j++;
  }
  printVec2D(v);
  printVec2D(H);
  return v[0];
}

int main()
{
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

  vector<vector<double> > A = {{1, 3, 0, 5},
                              {3, 5, 0, 8},
                              {0, 0, 6, 12},
                              {5, 8, 12, 16}};
  printVec2D(A);

  vector<vector<double> > L(A.size(), vector<double>(A[0].size()));
  vector<vector<double> > U(A.size(), vector<double>(A[0].size()));
  LUdecomp(A, L, U);
  cout << "L:\n";
  printVec2D(L);
  invertLowerTri(L);
  cout << "Linv:\n";
  printVec2D(L);
  cout << "U:\n";
  printVec2D(U);
  invertUpperTri(U);
  cout << "Uinv:\n";
  printVec2D(U);
  //
  // vector<vector<double> > Uinv = inverseLU(A);
  // printVec2D(Uinv);
  //
  // vector<vector<double> > B = MM(U, Uinv);
  // printVec2D(B);


  // vector<double> v1 = GMRES(A, b);
  // printVec(v1);

  // vector<vector<double> > A = {{1}, {2}};
  // printVec2D(A);
  // resizeMat(A, 3, 2);
  // printVec2D(A);

  return 0;
}
