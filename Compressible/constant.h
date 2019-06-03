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

using namespace std;
using namespace std::chrono;
// using namespace xt;

// Domain Constants
extern constexpr int imax = 204;
extern constexpr int jmax = 1;

extern constexpr double xmax = 1.0;
extern constexpr double ymax = 1.0;
extern constexpr double dx = xmax / (imax - 4);
extern constexpr double dy = ymax / (jmax);

extern constexpr double dt = 0.001;

// Constants
extern constexpr double pi = M_PI;	// Pi
extern constexpr double Rair = 1.0; // 287.058; // Specific gas constant for air
extern constexpr double Cp = 1.006; // Specific heat constant for air
extern constexpr double Cv = 0.7171; // Specific heat constant for air
extern constexpr double gam = 1.4; //Cp/Cv;

extern constexpr double rhoL = 6.0;		// Initial left side density
extern constexpr double rhoR = 1.0;		// Initial right side density
extern constexpr double uL = 0.0;		// Initial left side velocity
extern constexpr double uR = 0.0;		// Initial right side velocity
extern constexpr double PL = 12.0;		// Initial left side pressure
extern constexpr double PR = 1.0;		// Initial right side pressure

extern constexpr int p = 8; // Precision for printing
