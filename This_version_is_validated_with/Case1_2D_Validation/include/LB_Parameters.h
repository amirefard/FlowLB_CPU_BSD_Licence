#ifndef LB_PARAMETERS_H_
#define LB_PARAMETERS_H_

#include "LB_Settings.h"

const int Nx = 1201;    // number of lattice nodes along the x-axis
const int Ny = 401;      // number of lattice nodes along the y-axis
const double tau = 0.6; // relaxation time
const int t_num = 2000; // number of time steps (running from 1 to t_num)
const int t_data = 100;  // disk write time step (data will be written to the disk every t_data step)
const int t_vtk = 100;  // VTK write time step
const int t_info = 100; // info time step (screen message will be printed every t_info step)
const int t_time = 100; // Report elapsed time on screen
const int t_recovery = 50000;
const double gravity = 0.000;// force density due to gravity (in positive x-direction). Just for fluid without particle.

#if defined Flow_With_Resolved_Particle
//// Particle Properties  ////
const int particle_num_nodes = 200; // number of surface nodes
const int num_particles = 1; // number of particles.
const double particle_radius = 25.0; // radius of a circular cylinder OR the radius along the X-axis for a "ELLIPSE"
const double particle_center_x = 400.0; // center position (x-component) of first particle
const double particle_center_y = 200.0; // center position (y-component) of first particle
const double particle_gravity_x = 0.000013611111111;//This force is applied to particles according to density ratio
const double particle_gravity_y = 0.000;  //This force is applied to particles according to density ratio
const double den_ratio = 1.25; //Particle to fluid density ratio. Note that in below simulations denisty of fluid is assumed to be 1.0
const double specific_heat_ratio = 1.0; //Particle to fluid specific heat ratio.
const double particle_radius_2 = 6.0; //Just used for "ELLIPSE": radius along the Y-axis
const double init_ang = 0.5; //Initial angle of ELLIPSE with respect to X-dir in clock-wise direction (Jahate Mosalasati)
const double offset_x = 51; //Distance of particles from each other in x direction.
const double offset_y = 0; //Distance of particles from each other in y direction.

////Particle collision parameters
const double cij = 0.002;
const double ep = 0.1;//stifness factor for particle-particle collision
const double Ep = 0.05;//stifness factor for particle-particle collision usually smaller than ep
const double epw = 0.1;//stifness factor for particle-wall collision
const double Epw = 0.05;//stifness factor for particle-wall collision usually smaller than epw
const double th = 2.0;// Thereshold
const double thw = 2.0;// Thereshold
#endif

#if defined Flow_With_Point_Particle
//// Particle Properties  ////
const int num_point_particles = 1; // number of particles.
const double point_particle_radius = 0.5; // radius of a circular cylinder OR the radius along the X-axis for a "ELLIPSE"
const double point_particle_center_x = 49; // center position (x-component) of first particle
const double point_particle_center_y = 70; // center position (y-component) of first particle
const double point_particle_gravity_x = 0.0001;//This force is applied to particles according to density ratio
const double point_particle_gravity_y = 0.000;  //This force is applied to particles according to density ratio
const double point_den_ratio = 5; //Particle to fluid density ratio. Note that in below simulations denisty of fluid is assumed to be 1.0
//const double specific_heat_ratio = 1.0; //Particle to fluid specific heat ratio.
//const double particle_radius_2 = 6.0; //Just used for "ELLIPSE": radius along the Y-axis
//const double init_ang = 0.0; //Initial angle of ELLIPSE with respect to X-dir in clock-wise direction (Jahate Mosalasati)
const double point_offset_x = 5; //Distance of particles from each other in x direction.
const double point_offset_y = 0; //Distance of particles from each other in y direction.
#endif

#if defined Flow_With_Thermal_Effect
const double prandtl = 0.7; //prandtl number
const double gbeta = 0.00007; //In NATURAL convection: Ra=Gr*Pr=gbeta*dT*H^3/(alpha_t*nu)==>gbeta=alpha_t*nu*Ra/(dT*H^3). If natural convection is not important, set gbeta=0.0. Negative gbeta=g works in positive direction
const double Tref = 0.0; //reference temperature in case of natural convection

#if defined Particle_Fixed_Temperature
const double particle_temperature = 1.0; //Particle temperature must be constant
#endif
#if defined Particle_Varying_Temperature
const double particle_temperature = 0.0; //Here temperature is the particle INITIAL temperature
const double particle_heat_source = 0.31416; //Heat source in the total volume of the particle (e.g. watts)
#endif
#endif

#if defined Flow_Without_Thermal_Effect
//// Do NOT change 4 below lines
const double prandtl = 1.0; //prandtl number
const double gbeta = 0.0;// Do NOT change
const double Tref = 0.0; //reference teperature in case of natural convection. Do NOT change
const double particle_temperature = 0.0; //Particle temperature must be constant. Do NOT change
#endif

#if defined Stop_by_Residual
const double residual_flow = 1e-8;
const double residual_thermal = 1e-8;
const int t_residual = 2000;
#endif

const int surface_num_nodes = 500;//number of surface nodes on a solid boundary

//-------------------------------------------------------
//-------------------------------------------------------
//--------------DO NOT MODIFY----------------------------
//-------------------------------------------------------
//-------------------------------------------------------
const double omega = 1. / tau; // relaxation frequency (inverse of relaxation time)
//#if defined Flow_With_Thermal_Effect
const double alpha_t = (tau - 0.5) / 3. / prandtl; // Thermal diffusion coefficient
const double tau_t = 3.* alpha_t + 0.5; // relaxation time of temperature
const double omega_t = 1. / tau_t; //For fluid nodes
#if defined Flow_Without_Particle
const double omega_ts = 1. / (3.* solid_to_fluid_conductivity* alpha_t + 0.5);//for solid nodes
#endif



#endif  //end of header
