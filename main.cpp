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
	assert(timeStep <= maxPossibleTimeStep && "Time step is too large for stability! Killing the program.");

	std::vector<double> duPrevdx(NumPressurePoissonIterations);
	std::vector<double> duPrevdy(NumPressurePoissonIterations);
	std::vector<double> dvPrevdx(NumPressurePoissonIterations);
	std::vector<double> dvPrevdy(NumPressurePoissonIterations);
	std::vector<double> uPrevLaplace(NumPressurePoissonIterations);
	std::vector<double> vPrevLaplace(NumPressurePoissonIterations);

	for (int iter = 0; iter < NumPressurePoissonIterations; ++iter) {
		duPrevdx = centralDifferenceX(u_prev , numPoints , elementLength);
		duPrevdy = centralDifferenceY(u_prev , numPoints , elementLength);
		dvPrevdx = centralDifferenceX(v_prev , numPoints , elementLength);
		dvPrevdy = centralDifferenceY(v_prev , numPoints , elementLength);
		uPrevLaplace = laplacian(u_prev , numPoints , elementLength);
		vPrevLaplace = laplacian(v_prev , numPoints , elementLength);
		if (iter % 100 == 0) {
			std::cout << "\n================== Iteration " << iter << " ==================\n";
			std::cout << std::setw(20) << std::left << "duPrevdx" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << duPrevdx[iter] << "\n";
			std::cout << std::setw(20) << std::left << "duPrevdy" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << duPrevdy[iter] << "\n";
			std::cout << std::setw(20) << std::left << "dvPrevdx" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << dvPrevdx[iter] << "\n";
			std::cout << std::setw(20) << std::left << "dvPrevdy" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << dvPrevdy[iter] << "\n";
			std::cout << std::setw(20) << std::left << "uPrevLaplace" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << uPrevLaplace[iter] << "\n";
			std::cout << std::setw(20) << std::left << "vPrevLaplace" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << vPrevLaplace[iter] << "\n";
			std::cout << "======================================================\n";
		}

			// Time stepping (placeholders — implement update logic)
		std::vector<double> u_next = u_prev;
		std::vector<double> v_next = v_prev;
		std::vector<double> p_next = p_prev;

	}




	return 0;

}
