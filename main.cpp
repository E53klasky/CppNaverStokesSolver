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

// probly have to test this
std::pair<std::vector<double> , std::vector<double>> centralDifference(std::vector<std::vector<double>>& fx ,
	std::vector<std::vector<double>>& fy , double elementLength) {
	std::vector<double> diffx(fx.size());
	std::vector<double> diffy()(fy.size());
	int numPoints = fx.size();

	for (size_t i = 0; i < numPoints; i++)
	{
		for (size_t j = 0; j < numPoints; j++)
		{
			int center = i * numPoints + j;
			int right = i * numPoints + (j + 1);
			int left = i * numPoints + (j - 1);
			diffx[center] = (f[right] - f[left]) / (2.0 * elementLength);
			diffyp[center] = (f[right] - f[left]) / (2.0 * elementLength);

		}

	}

	// this will be a central difference funcion 
	return std::make_pair(diffx , diffy);
}

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
