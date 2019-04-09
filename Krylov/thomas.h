/* This definition must be changed if you want to use a mesh with more than 
   100 control volumes. */ 
#define NMAX 100

/* Solve a tri-diagonal system Ax = b.
 * 
 * Uses the Thomas algorithm, which is Gauss elimination and back
 * substitution specialized for a tri-diagonal matrix.
 *
 * Input:
 *  LHS: The matrix A.  The three columns are the three non-zero diagonals.
 *     (0: below main diag; 1: on main diag; 2: above main diag)
 *  RHS: A vector containing b.
 *  iSize: Number of equations.  For situations with BC's, those are included
 *     in the count.
 *
 * Output:
 *  LHS: Garbled.
 *  RHS: The solution x.
 */
void SolveThomas(vector<vector<double>> &LHS, vector<double> &RHS, int iSize)
{
  int i;
  /* This next line actually has no effect, but it -does- make clear that
     the values in those locations have no impact. */
  LHS[0][0] = LHS[iSize-1][2] = 0;
  /* Forward elimination */
  for (i = 0; i < iSize-1; i++) {
    LHS[i][2] /= LHS[i][1];
    RHS[i] /= LHS[i][1];
    LHS[i+1][1] -= LHS[i][2]*LHS[i+1][0];
    RHS[i+1] -= LHS[i+1][0]*RHS[i];
  }
  /* Last line of elimination */
  RHS[iSize-1] /= LHS[iSize-1][1];

  /* Back-substitution */
  for (i = iSize-2; i >= 0; i--) {
    RHS[i] -= RHS[i+1]*LHS[i][2];
  }
}

// /* A test program to confirm proper behavior. */
// #include <stdio.h>
// int thomas()
// {
//   double LHS[12][3], RHS[12];
//   int i;
  
//   for (i = 0; i <= 11; i++) {
//     LHS[i][0] = 1;
//     LHS[i][1] = 2+i;
//     LHS[i][2] = 3;
//     RHS[i]   = i;
//   }

//   SolveThomas(LHS, RHS, 12);
      
//   for (i = 0; i <= 11; i++) {
//     double result;
//     if (i == 0) result = RHS[0]*2 + RHS[1]*3;
//     else if (i == 11) result = RHS[10]*1 + RHS[11]*(2+11);
//     else result = RHS[i-1]*1 + RHS[i]*(2+i) + RHS[i+1]*3;
//     printf("%2d Soln: %10.6f \tRHS recomputed: %10.6f  \tError: %10.6G\n", i, RHS[i], result, result - i);
//   }

//   getchar();
//   return(0);
// }
