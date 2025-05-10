#include "mathFunc.h"
// TODO: at the 200 line convert the python to the 328 line + wright out the data 

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


int main() {

	// parms setting ---------------------------------
	int numPoints = 41;
	float domainSize = 1.0;
	float timeStep = 0.001;
	float kinematicViscoity = 0.1;
	float density = 1.0;
	float horizontalVelocityTop = 1.0;
	int NumPressurePoissonIterations = 500;
	float statbilitySafetyFactor = 0.5;


	parameterPrint(numPoints , domainSize , timeStep ,
		kinematicViscoity , density ,
		horizontalVelocityTop , NumPressurePoissonIterations ,
		statbilitySafetyFactor);
	// end of parms  ---------------------------------	

	double elementLength = domainSize / (numPoints - 1); // 1/40 = 0.025  
	std::vector<double> x(numPoints);
	std::vector<double> y(numPoints);
	for (int i = 0; i < numPoints; ++i) {
		x[i] = i * elementLength;
	}
	y = x;
	int size = numPoints * numPoints;
	std::vector<double> u_prev(size , 0.0);
	std::vector<double> v_prev(size , 0.0);
	std::vector<double> p_prev(size , 0.0);

	double maxPossibleTimeStep = 0.5 * elementLength * elementLength / kinematicViscoity;
	assert(maxPossibleTimeStep >= timeStep * maxPossibleTimeStep && "Time step is too large for stability! Killing the program.");
	// maybe add it for div * U ????????????? later todo 

	for (int iter = 0; iter < NumPressurePoissonIterations; ++iter) {
		// bug in here??? 
		std::vector<double> d_u_prev__d_x = centralDifferenceX(u_prev , size , elementLength);
		std::vector<double> d_u_prev__d_y = centralDifferenceY(u_prev , size , elementLength);
		std::vector<double> d_v_prev__d_x = centralDifferenceX(v_prev , size , elementLength);
		std::vector<double> d_v_prev__d_y = centralDifferenceY(v_prev , size , elementLength);
		std::vector<double> laplace__u_prev = laplacian(u_prev , size , elementLength);
		std::vector<double> laplace__v_prev = laplacian(v_prev , size , elementLength);

		// // Optional: simple progress print
		// if (iter % 100 == 0) {
		// 	std::cout << "Iteration " << iter << " / " << N_ITERATIONS << std::endl;
		// }
	}





			// moving foward in time
	std::vector<double> u_next = u_prev;
	std::vector<double> v_next = v_prev;
	std::vector<double> p_next = p_prev;


	return 0;

}
