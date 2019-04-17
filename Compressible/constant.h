//==================================================
// Header file containing constants and includes
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <ctime>
#include <string>
#include <chrono>
// #include "xtensor/xarray.hpp"
// #include "xtensor/xio.hpp"
// #include "xtensor/xview.hpp"

using namespace std;
using namespace std::chrono;
// using namespace xt;

// Domain Constants
extern constexpr int imax = 14;
extern constexpr int jmax = 6;

extern constexpr double xmax = 1.;
extern constexpr double ymax = 0.1;
extern constexpr double dx = xmax / (imax - 2);
extern constexpr double dy = ymax / (jmax - 2);

extern constexpr int ubar = 3;		// average velocity in x [m/s]
extern constexpr double dt = 0.05;
extern constexpr double tol = pow(10, -9); // Tolerance for max change, "nano"

// Constants
extern constexpr double pi = M_PI;	// Pi
extern constexpr int Re = 250.;		// Reynolds number
extern constexpr double B = 1.;
extern constexpr double w = 1.;

extern constexpr double R = 287.058; // Specific gas constant for air
extern constexpr double Cp = 1.006; // Specific heat constant for air
extern constexpr double Cv = 0.7171; // Specific heat constant for air
extern constexpr double gam = Cp/Cv;

extern constexpr int rhoL = 6;		// Initial left side density
extern constexpr int rhoR = 1;		// Initial right side density
extern constexpr int uL = 0;		// Initial left side velocity
extern constexpr int uR = 0;		// Initial right side velocity
extern constexpr int PL = 12;		// Initial left side pressure
extern constexpr int PR = 1;		// Initial right side pressure

extern constexpr int p = 6; // Precision for printing
