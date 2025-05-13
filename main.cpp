#include "mathFunc.h"
#include <adios2.h>

// python example here: https://github.com/Ceyron/machine-learning-and-simulation/blob/main/english/simulation_scripts/lid_driven_cavity_python_simple.py 

// Enhanced parameter print ---------------------------------
void parameterPrint(int numPoints , float domainSize , float timeStep ,
	float kinematicViscosity , float density ,
	float horizontalVelocityTop , int NumPressurePoissonIterations ,
	float statbilitySafetyFactor) {
// Header decoration
	std::cout << "\n******************************* PARAMETERS *******************************\n";
	std::cout << "***************************************************************************\n";

	// Print the parameters with formatting
	std::cout << std::left << std::setw(40) << "Number of Points: " << std::setw(10) << numPoints << std::endl;
	std::cout << std::left << std::setw(40) << "Domain Size: " << std::setw(10) << domainSize << std::endl;
	std::cout << std::left << std::setw(40) << "Time Step: " << std::setw(10) << timeStep << std::endl;
	std::cout << std::left << std::setw(40) << "Kinematic Viscosity: " << std::setw(10) << kinematicViscosity << std::endl;
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
	float kinematicViscosity = 0.1;
	float density = 1.0;
	float horizontalVelocityTop = 1.0;
	int NumPressurePoissonIterations = 50;
	int numIterations = 500;
	float statbilitySafetyFactor = 0.5;


	parameterPrint(numPoints , domainSize , timeStep ,
		kinematicViscosity , density ,
		horizontalVelocityTop , NumPressurePoissonIterations ,
		statbilitySafetyFactor);
	// end of parms  ---------------------------------	

	double elementLength = domainSize / (numPoints - 1);
	std::vector<double> x(numPoints);
	std::vector<double> y(numPoints);

	// adios declaration ---------------------------------------------------------
	adios2::ADIOS adios;
	adios2::IO io = adios.DeclareIO("SimulationOutput");
	auto varU = io.DefineVariable<double>("u" , { numPoints , numPoints } , { 0,0 } , { numPoints , numPoints });
	auto varV = io.DefineVariable<double>("v" , { numPoints , numPoints } , { 0,0 } , { numPoints , numPoints });
	auto varP = io.DefineVariable<double>("p" , { numPoints , numPoints } , { 0,0 } , { numPoints , numPoints });

	std::string filename = "data.bp";
	adios2::Engine engine = io.Open(filename , adios2::Mode::Write);
	// end of adios declaration -------------------------------------------------

	for (int i = 0; i < numPoints; ++i) {
		x[i] = i * elementLength;
	}
	y = x;
	int size = numPoints * numPoints;
	std::vector<double> u_prev(size , 0.0);
	std::vector<double> v_prev(size , 0.0);
	std::vector<double> p_prev(size , 0.0);

	double maxPossibleTimeStep = 0;
	maxPossibleTimeStep = 0.5 * elementLength * elementLength / kinematicViscosity;
	if (timeStep > statbilitySafetyFactor * maxPossibleTimeStep) {
		std::cerr << "Time step is too large for stability! Killing the program." << std::endl;
		return 1;
	}

	for (int iter = 0; iter < numIterations; ++iter) {
		std::vector<double> duPrevdx = centralDifferenceX(u_prev , numPoints , elementLength);
		std::vector<double> duPrevdy = centralDifferenceY(u_prev , numPoints , elementLength);
		std::vector<double> dvPrevdx = centralDifferenceX(v_prev , numPoints , elementLength);
		std::vector<double> dvPrevdy = centralDifferenceY(v_prev , numPoints , elementLength);
		std::vector<double> uPrevLaplace = laplacian(u_prev , numPoints , elementLength);
		std::vector<double> vPrevLaplace = laplacian(v_prev , numPoints , elementLength);


		if (iter % 100 == 0) {
			std::cout << "\n================== Iteration " << iter << " ==================\n";
			int sampleIdx = idx(numPoints / 2 , numPoints / 2 , numPoints);
			std::cout << std::setw(20) << std::left << "duPrevdx sample" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << duPrevdx[sampleIdx] << "\n";
			std::cout << std::setw(20) << std::left << "duPrevdy sample" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << duPrevdy[sampleIdx] << "\n";
			std::cout << std::setw(20) << std::left << "dvPrevdx sample" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << dvPrevdx[sampleIdx] << "\n";
			std::cout << std::setw(20) << std::left << "dvPrevdy sample" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << dvPrevdy[sampleIdx] << "\n";
			std::cout << std::setw(20) << std::left << "uPrevLaplace sample" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << uPrevLaplace[sampleIdx] << "\n";
			std::cout << std::setw(20) << std::left << "vPrevLaplace sample" << std::setw(15) << std::right << std::fixed << std::setprecision(6) << vPrevLaplace[sampleIdx] << "\n";
			std::cout << "======================================================\n";
		}

		std::vector<double> u_tent(size , 0.0);
		std::vector<double> v_tent(size , 0.0);

		for (int i = 0; i < numPoints; ++i) {
			for (int j = 0; j < numPoints; ++j) {
				int idx_val = idx(i , j , numPoints);
				u_tent[idx_val] = u_prev[idx_val] + timeStep * (
					-(u_prev[idx_val] * duPrevdx[idx_val] +
						v_prev[idx_val] * duPrevdy[idx_val]) +
					kinematicViscosity * uPrevLaplace[idx_val]
					);

				v_tent[idx_val] = v_prev[idx_val] + timeStep * (
					-(u_prev[idx_val] * dvPrevdx[idx_val] +
						v_prev[idx_val] * dvPrevdy[idx_val]) +
					kinematicViscosity * vPrevLaplace[idx_val]
					);
			}
		}


		for (int j = 0; j < numPoints; ++j) {
			u_tent[idx(0 , j , numPoints)] = 0.0;
			v_tent[idx(0 , j , numPoints)] = 0.0;
		}

		for (int i = 0; i < numPoints; ++i) {
			u_tent[idx(i , 0 , numPoints)] = 0.0;
			v_tent[idx(i , 0 , numPoints)] = 0.0;
		}

		for (int i = 0; i < numPoints; ++i) {
			u_tent[idx(i , numPoints - 1 , numPoints)] = 0.0;
			v_tent[idx(i , numPoints - 1 , numPoints)] = 0.0;
		}

		for (int j = 0; j < numPoints; ++j) {
			u_tent[idx(numPoints - 1 , j , numPoints)] = horizontalVelocityTop;
			v_tent[idx(numPoints - 1 , j , numPoints)] = 0.0;
		}

		std::vector<double> duTentdx = centralDifferenceX(u_tent , numPoints , elementLength);
		std::vector<double> dvTentdy = centralDifferenceY(v_tent , numPoints , elementLength);

		std::vector<double> rhs(size , 0.0);

		for (int i = 1; i < numPoints - 1; ++i) {
			for (int j = 1; j < numPoints - 1; ++j) {
				int idx_val = idx(i , j , numPoints);
				rhs[idx_val] = density / timeStep * (
					duTentdx[idx_val] + dvTentdy[idx_val]
					);
			}
		}

		std::vector<double> p_next = p_prev;


		for (int poisson_iter = 0; poisson_iter < NumPressurePoissonIterations; ++poisson_iter) {
			std::vector<double> p_temp(size , 0.0);

			for (int i = 1; i < numPoints - 1; ++i) {
				for (int j = 1; j < numPoints - 1; ++j) {
					int center = idx(i , j , numPoints);
					int left = idx(i , j - 1 , numPoints);
					int right = idx(i , j + 1 , numPoints);
					int up = idx(i - 1 , j , numPoints);
					int down = idx(i + 1 , j , numPoints);

					p_temp[center] = 0.25 * (
						p_next[left] + p_next[right] +
						p_next[up] + p_next[down] -
						elementLength * elementLength * rhs[center]
						);
				}
			}


			for (int i = 0; i < numPoints; ++i) {
				p_temp[idx(i , 0 , numPoints)] = p_temp[idx(i , 1 , numPoints)];
			}

			for (int i = 0; i < numPoints; ++i) {
				p_temp[idx(i , numPoints - 1 , numPoints)] = p_temp[idx(i , numPoints - 2 , numPoints)];
			}

			for (int j = 0; j < numPoints; ++j) {
				p_temp[idx(0 , j , numPoints)] = p_temp[idx(1 , j , numPoints)];
			}

			for (int j = 0; j < numPoints; ++j) {
				p_temp[idx(numPoints - 1 , j , numPoints)] = 0.0;
			}

			p_next = p_temp;
		}

		std::vector<double> dpNextdx = centralDifferenceX(p_next , numPoints , elementLength);
		std::vector<double> dpNextdy = centralDifferenceY(p_next , numPoints , elementLength);

		std::vector<double> u_next(size , 0.0);
		std::vector<double> v_next(size , 0.0);

		// Velocity correction
		for (int i = 1; i < numPoints - 1; ++i) {
			for (int j = 1; j < numPoints - 1; ++j) {
				int idx_val = idx(i , j , numPoints);

				u_next[idx_val] = u_tent[idx_val] - timeStep / density * dpNextdx[idx_val];
				v_next[idx_val] = v_tent[idx_val] - timeStep / density * dpNextdy[idx_val];
			}
		}

		// Apply velocity boundary conditions
		for (int j = 0; j < numPoints; ++j) {
			u_next[idx(0 , j , numPoints)] = 0.0;
			v_next[idx(0 , j , numPoints)] = 0.0;
		}

		for (int i = 0; i < numPoints; ++i) {
			u_next[idx(i , 0 , numPoints)] = 0.0;
			v_next[idx(i , 0 , numPoints)] = 0.0;
		}

		for (int i = 0; i < numPoints; ++i) {
			u_next[idx(i , numPoints - 1 , numPoints)] = 0.0;
			v_next[idx(i , numPoints - 1 , numPoints)] = 0.0;
		}

		for (int j = 0; j < numPoints; ++j) {
			u_next[idx(numPoints - 1 , j , numPoints)] = horizontalVelocityTop;
			v_next[idx(numPoints - 1 , j , numPoints)] = 0.0;
		}


		u_prev = u_next;
		v_prev = v_next;
		p_prev = p_next;


		engine.BeginStep();
		engine.Put(varU , u_prev.data());
		engine.Put(varV , v_prev.data());
		engine.Put(varP , p_prev.data());
		engine.EndStep();

	}

	engine.Close();
	std::cout << "Simulation completed successfully" << std::endl;
	std::cout << "Data written to " << filename << " in the same place as the exe" << std::endl;
	return 0;
}