#ifndef THERMAL_SOLVER_H
#define THERMAL_SOLVER_H

#include <sstream>  // string streams
#include <iostream> // for the use of 'cout'
#include <vector>
using namespace std;

#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"

//#if defined Flow_With_Thermal_Effect

class Thermal_solver
{
private:
	int az;
    vector<int> cx_pop_T, cy_pop_T;
    vector<double> weight_T;//Thermal weights
    vector<int> bounce_back_pop_T;//bounce back direction of each direction (for thermal)
    vector<int> pop_to_top_T, pop_to_bottom_T, pop_to_left_T, pop_to_right_T;
    vector<int> pop_reverse_top_T, pop_reverse_bottom_T, pop_reverse_left_T, pop_reverse_right_T;

public:
	Thermal_solver();
	void initialize_t();
	void equilibrium_t(double, double, double, double); // compute the equilibrium populations of temperature from the fluid density and velocity
	void LBM_t(int); // perform LBM Thermal operations
	void BC_t(int); // perform boundary conditions for temperature operations
	void BC_t_OLD_VERSION(int);
	void momenta_t(); // compute temperature from the populations
	bool check_residual();
	~Thermal_solver();
};

#endif

//#endif
