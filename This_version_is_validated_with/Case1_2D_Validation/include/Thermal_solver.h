#ifndef THERMAL_SOLVER_H
#define THERMAL_SOLVER_H

#include <sstream>  // string streams
#include <iostream> // for the use of 'cout'

#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"

//#if defined Flow_With_Thermal_Effect

class Thermal_solver
{
private:
	int az;
public:
	Thermal_solver();
	void initialize_t();
	void equilibrium_t(double, double, double, double); // compute the equilibrium populations of temperature from the fluid density and velocity
	void LBM_t(int); // perform LBM Thermal operations
	void BC_t(int); // perform boundary conditions for temperature operations
	void momenta_t(); // compute temperature from the populations
	bool check_residual();
	~Thermal_solver();
};

#endif

//#endif
