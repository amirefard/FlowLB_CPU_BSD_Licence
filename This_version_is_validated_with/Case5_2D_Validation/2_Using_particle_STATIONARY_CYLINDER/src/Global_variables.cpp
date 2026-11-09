//#include "stdafx.h"
#include "Global_variables.h"

double ***pop_F = nullptr, ***pop_old_F = nullptr; // LBM populations (old and new)
double **density = nullptr; // fluid density
double **velocity_x = nullptr; // fluid velocity (x-component)
double **velocity_y = nullptr; // fluid velocity (y-component)
double **velocity_magnitude = nullptr; // fluid velocity magnitude
double **previous_velocity_magnitude = nullptr; // fluid velocity magnitude at previous time step

double **force_x = nullptr; // fluid force (x-component)
double **force_y = nullptr; // fluid force (y-component)
int **is_solid = nullptr;
int **solid_thermal = nullptr;
double **mdf_force_x = nullptr; // fluid force for MDF (x-component)
double **mdf_force_y = nullptr; // fluid force for MDF(y-component)
double **mdf_force_thermal = nullptr; // fluid force for MDF(thermal)
double **mdf_vel_x = nullptr; // vel_x for MDF(x-component)
double **mdf_vel_y = nullptr; // vel_y for MDF(y-component)
double **mdf_temperature = nullptr; // temperature for MDF
double **temperature = nullptr; // fluid temperature
double **previous_temperature = nullptr; // fluid temperature
double **force_thermal = nullptr; // force term for thermal lbm

double ***pop_T = nullptr, ***pop_old_T = nullptr; // LBM populations of temperature (old and new)
double pop_eq_F[9]; // equilibrium populations
double pop_eq_T[9]; // equilibrium populations of temperature

//solid_surface solid; // create solid surface tracer
//Inlet_outlets Inlets_Outlets;

int **boundary_cond = nullptr; //Primary boundary condition of the point: 0) Fluid 1) Constant Velocity 2) Constant Pressure 3) Zero Gradient 4) Stationary wall
double **boundary_vel_x = nullptr;//velocity x value at the boundary
double **boundary_vel_y = nullptr;//velocity y value at the boundary
double **boundary_pres = nullptr;//pressure value at the boundary

double **boundary_temp = nullptr;//temperature value at the boundary
double **boundary_temp_2 = nullptr;//temperature value at the boundary for edges
int **boundary_cond_temp = nullptr;//0) Fluid 1) Constant Temperature 2) Adiabatic
int **boundary_cond_temp_2 = nullptr;//0) Fluid 1) Constant Temperature 2) Adiabatic

int **edge_boundary = nullptr;
int **boundary_cond_2 = nullptr;//second boundary condition for wall
double **boundary_vel_x_2 = nullptr;//second velocity x value at the boundary (at edges only)
double **boundary_vel_y_2 = nullptr;//second velocity y value at the boundary (at edges only)
double **boundary_pres_2 = nullptr;//second pressure value at the boundary (at edges only)

double ***delta = nullptr;
int recovery_step=0;
