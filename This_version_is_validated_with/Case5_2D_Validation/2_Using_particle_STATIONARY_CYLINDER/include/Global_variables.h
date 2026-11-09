#ifndef ALBORZ_GLOBALVARIABLES_H_
#define ALBORZ_GLOBALVARIABLES_H_

#include <stdio.h>
#include "LB_Settings.h"
#include "LB_Parameters.h"

extern double ***pop_F, ***pop_old_F; // LBM populations (old and new)
extern double **density; // fluid density
extern double **velocity_x; // fluid velocity (x-component)
extern double **velocity_y; // fluid velocity (y-component)
extern double **velocity_magnitude; // fluid velocity magnitude
extern double **previous_velocity_magnitude; // fluid velocity magnitude at previous time step
extern double **force_x; // fluid force (x-component)
extern double **force_y; // fluid force (y-component)
extern int **is_solid;
extern int **solid_thermal;//Thermal Boundary condition for solid nodes. 0) Fluid, 1) Solid with constant temperature, 2) Solid adiabatic, 3) Solid conjugate
extern double **mdf_force_x; // fluid force for MDF (x-component)
extern double **mdf_force_y; // fluid force for MDF(y-component)
extern double **mdf_force_thermal; // fluid force for MDF(thermal)
extern double **mdf_vel_x; // vel_x for MDF(x-component)
extern double **mdf_vel_y; // vel_y for MDF(y-component)
extern double **mdf_temperature; // temperature for MDF
extern double **temperature; // fluid temperature
extern double **previous_temperature; // fluid temperature
extern double **force_thermal; // force term for thermal lbm

extern double ***pop_T, ***pop_old_T; // LBM populations of temperature (old and new)
extern double pop_eq_F[9]; // equilibrium populations
extern double pop_eq_T[9]; // equilibrium populations of temperature
//const double weight[9] = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. }; // lattice weights

//extern solid_surface solid; // create solid surface tracer
extern int **boundary_cond; // fluid density
//extern Inlet_outlets Inlets_Outlets;
extern double **boundary_vel_x;//velocity x value at the boundary
extern double **boundary_vel_y;//velocity y value at the boundary
extern double **boundary_pres;// pressure value at the boundary
extern int **edge_boundary;
extern int **boundary_cond_2;//second boundary condition for wall
extern double **boundary_vel_x_2;//second boundary condition for velocity x value at the boundary
extern double **boundary_vel_y_2;//second boundary condition for velocity y value at the boundary
extern double **boundary_pres_2;// second boundary condition for pressure value at the boundary

extern double **boundary_temp;//temperature value at the boundary
extern double **boundary_temp_2;//other temperature value at the boundary for edges
extern int **boundary_cond_temp;
extern int **boundary_cond_temp_2;

extern double ***delta;// delta distance for a solid boundary

extern int recovery_step;

//-------------------------------------------------------
//-------------------------------------------------------
#endif

