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

using namespace std;

// Domain Constants
extern constexpr int imax = 66;
extern constexpr int jmax = 66;

extern constexpr double xmin = 0.0;
extern constexpr double ymin = 0.0;
extern constexpr double xmax = 1.0;
extern constexpr double ymax = 1.0;
extern constexpr double dx = xmax / (imax - 2);
extern constexpr double dy = ymax / (jmax - 2);

// extern constexpr double dt = 0.05;

extern constexpr double tol = pow(10, -9); // Tolerance for max change, "nano"

// // Constants
extern constexpr double pi = M_PI;	// Pi
extern constexpr double w = 2.0/3; //0.5395;	// Overrelaxation constant
extern constexpr double w2 = 0;//0.8724; // Second Overrelaxation constant for double pass
extern constexpr int p = 4;		// Precision for printing to console
