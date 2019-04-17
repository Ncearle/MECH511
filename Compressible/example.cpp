#include "constant.h"

int main(int argc, char* argv[])
{
    xarray<double> arr1
      {{1.0, 2.0, 3.0},
       {2.0, 5.0, 7.0},
       {2.0, 5.0, 7.0}};

    xarray<double> arr2
      {5.0,6.0,7.0};

    xarray<double> res = arr1 * arr2;

    cout << res << endl;

    return 0;
}
