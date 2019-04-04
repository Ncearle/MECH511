//=================================================
// Test File
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "thomas.h"
#include "vecmat_funcs.h"


void init(vector<vector<double> > &T, vector<vector<double> > &u, vector<vector<double> > &v)
{
	for (int j = 0; j < jmax-2; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax-2; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			T[j][i] = y;
			u[j][i] = 6.0*ubar*y*(1-y);
			v[j][i] = 0.;
		}
	}
	// ghost(T, u, v);
}

int main()
{
  vector<vector<double> > T((jmax-2), vector<double>(imax-2));
  vector<vector<double> > u((jmax-2), vector<double>(imax-2));
  vector<vector<double> > v((jmax-2), vector<double>(imax-2));
  init(T, u, v);

  printVec2D(T);
  vector<double> V = mat2vec(T);
  printVec(V);
  vector<vector<double> > M = vec2mat(V);
  printVec2D(M);


  return 0;
}
