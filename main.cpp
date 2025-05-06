#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <random>
#include <iostream>
#include <algorithm>  
#include <string.h>
#include <iomanip>
#include "adios2.h"
#include "mgard/compress_x.hpp"

// python example here: https://github.com/Ceyron/machine-learning-and-simulation/blob/main/english/simulation_scripts/lid_driven_cavity_python_simple.py 

// Enhanced parameter print ---------------------------------
void parameterPrint(int numPoints , float domainSize , float timeStep ,
	float kinematicViscoity , float density ,
	float horizontalVelocityTop , int NumPressurePoissonIterations ,
	float statbilitySafetyFactor) {
// Header decoration
	std::cout << "\n******************************* PARAMETERS *******************************\n";
	std::cout << "***************************************************************************\n";

	// Print the parameters with formatting
	std::cout << std::left << std::setw(40) << "Number of Points: " << std::setw(10) << numPoints << std::endl;
	std::cout << std::left << std::setw(40) << "Domain Size: " << std::setw(10) << domainSize << std::endl;
	std::cout << std::left << std::setw(40) << "Time Step: " << std::setw(10) << timeStep << std::endl;
	std::cout << std::left << std::setw(40) << "Kinematic Viscosity: " << std::setw(10) << kinematicViscoity << std::endl;
	std::cout << std::left << std::setw(40) << "Density: " << std::setw(10) << density << std::endl;
	std::cout << std::left << std::setw(40) << "Horizontal Velocity (Top): " << std::setw(10) << horizontalVelocityTop << std::endl;
	std::cout << std::left << std::setw(40) << "Pressure Poisson Iterations: " << std::setw(10) << NumPressurePoissonIterations << std::endl;
	std::cout << std::left << std::setw(40) << "Stability Safety Factor: " << std::setw(10) << statbilitySafetyFactor << std::endl;

	// Footer decoration
	std::cout << "***************************************************************************\n";
	std::cout << "***************************************************************************\n";
}
// ---------------------------------------------------------------------------

// index method
int idx(int i , int j , int N) {
	return i * N + j;
}


// difference methods -------------------------------------------
std::vector<double> centralDifferenceX(const std::vector<double>& f , int N , double dx) {
	std::vector<double> diff(N * N , 0.0);

	for (int i = 1; i < N - 1; ++i) {
		for (int j = 1; j < N - 1; ++j) {
			int center = idx(i , j , N);
			int left = idx(i , j - 1 , N);
			int right = idx(i , j + 1 , N);

			diff[center] = (f[right] - f[left]) / (2.0 * dx);
		}
	}

	return diff;
}

std::vector<double> centralDifferenceY(
	const std::vector<double>& f ,
	int rows ,
	int cols ,
	double elementLength
) {
	std::vector<double> diff(rows * cols , 0.0);

	for (int i = 1; i < rows - 1; ++i) {
		for (int j = 1; j < cols - 1; ++j) {
			int center = idx(i , j , cols);
			int up = idx(i + 1 , j , cols);
			int down = idx(i - 1 , j , cols);

			diff[center] = (f[up] - f[down]) / (2.0 * elementLength);
		}
	}

	return diff;
}
// --------------------------------------------------------------------------

// laplacian method
std::vector<double> laplacian(
	const std::vector<double>& f ,
	int rows ,
	int cols ,
	double elementLength
) {
	std::vector<double> laplacian(rows * cols , 0.0);

	for (int i = 1; i < rows - 1; ++i) {
		for (int j = 1; j < cols - 1; ++j) {
			int center = idx(i , j , cols);
			int up = idx(i + 1 , j , cols);
			int down = idx(i - 1 , j , cols);
			int left = idx(i , j - 1 , cols);
			int right = idx(i , j + 1 , cols);

			laplacian[center] = (f[up] + f[down] + f[left] + f[right] - 4.0 * f[center]) / (elementLength * elementLength);
		}
	}

	return laplacian;
}
// --------------------------------------------------------------------------

int main() {

	// parms setting ---------------------------------
	int numPoints = 41;
	float domainSize = 1.0;
	float timeStep = 0.001;
	float kinematicViscoity = 0.1;
	float density = 1.0;
	float horizontalVelocityTop = 1.0;
	int NumPressurePoissonIterations = 60;
	float statbilitySafetyFactor = 0.5;

	parameterPrint(numPoints , domainSize , timeStep ,
		kinematicViscoity , density ,
		horizontalVelocityTop , NumPressurePoissonIterations ,
		statbilitySafetyFactor);
	// end of parms  ---------------------------------	

	float elementLength = domainSize / (numPoints - 1);
	std::vector<double> x(numPoints);
	std::vector<double> y(numPoints);
	for (int i = 0; i < numPoints; ++i) {
		x[i] = i * elementLength;
	}
	y = x;
	std::vector<double> u_prev(numPoints * numPoints , 0.0);
	std::vector<double> v_prev(numPoints * numPoints , 0.0);
	std::vector<double> p_prev(numPoints * numPoints , 0.0);







	// moving foward in time
	std::vector<double> u_next = u_prev;
	std::vector<double> v_next = v_prev;
	std::vector<double> p_next = p_prev;


	return 0;

}
