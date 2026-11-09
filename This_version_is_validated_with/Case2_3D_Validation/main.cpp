// Welcome to ALBORZ 3D Ver 88.0
//Program written by: Amir Eshghinejadfard

// This program uses the D3Q19 lattice including forcing term.
// Rigid bottom and top walls are parallel to the x-axis (channel).
// particle(s) is positioned in the flow.
// All particles must have the same shape.
// This particle can be rigid and stationary/moving.
// The lattice velocities are defined according to the following scheme: D3Q19


// index:   0  1  2  3  4  5  6  7   8  9  10  11  12  13  14  15  16  17  18
// x:       0  1 -1  0  0  0  0  1  -1  1  -1  1   -1  1   -1   0   0   0   0
// y:       0  0  0  1 -1  0  0  1  1  -1  -1  0   0   0   0    1  -1   1  -1
// z:       0  0  0  0  0  1 -1  0  0   0  0   1   1   -1  -1   1   1  -1  -1

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************
//#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <vector>   // vector containers
#include <cmath>    // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream>  // file streams
#include <sstream>  // string streams
#include <cstdlib>  // standard library
#include <iomanip>  // For set precision. From ver 25
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <filesystem> // for std::filesystem
//#include <conio.h> // getch()
#ifdef __linux__
#include <sys/stat.h>//for mkdir
#else
#include "direct.h"//for mkdir
#endif
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
using namespace std; // permanently use the standard namespace

#define TRUE        1
#define FALSE       -1

/// *********************
/// SIMULATION PARAMETERS
/// *********************

// These are the relevant simulation parameters.
// They can be changed by the user.
// If a bottom or top wall shall move in negative x-direction, a negative velocity has to be specified.
// Moving walls and gravity can be switched on simultaneously.

#define Flow_With_Particle
//#define Flow_Without_Particle

//#define Flow_With_Thermal_Effect
#define Flow_Without_Thermal_Effect

// Condition of Read_Geometry_From_File is just for Flow_Without_Particle
// In this case you should:
// 1: convert the geometry to a matrix of 0 and 1 by the MATLAB program: image_to_matrix.m .This MATLAB rogram creates geometry.txt
// 2: put the geometry.txt in the same folder as where parameters.dat or fluid_data.dat are created.
// 3: Enter manually Nx and Ny values in ALBORZ program.
// Below condition can be used by different boundary conditions.
// If gravity has a non-zero value then Permeability is calculated and reported and saved to fluid_data.dat
// If you deactive below line you can still define some solid nodes in the geometry by putting is_solid=TRUE (See beginning of Flow_Without_Particle section)
// Important: Solid point have ZERO velocity in all directions

//#define Read_Geometry_From_File

// Force Calculation method.
// Guo: Z. Guo, C. Zheng, B. Shi, Discrete lattice effects on the forcing term in the lattice Boltzmann method, Phys. Rev. E 65 (046308) (2002) 16.
// Shan-Chen: Shan, X., Chen, H. (1993), Lattice Boltzmann model for simulating flows with multiple phase and components, Physical Review E, 47, 1815-1819
// Luo: L.-S. Luo, Lattice-Gas Automata and Lattice Boltzmann Equations for Two-Dimensional Hydrodynamics, Ph.D. thesis, Georgia Institute of Technology,1993.
// Guo model is a second order model but has less stability than Shan-Chen in low relaxation times. It is also somewhat slower.

#define Shan_Chen
//#define Guo
//#define Luo

// Select Single Relaxation Time (SRT) or Multi Relaxation Time (MRT)
// -MRT model is just adjusted for Guo force scheme.

#define SRT
//#define MRT

// Boundary types
// Exactly one of the following options has to be defined.
// X_Periodic
// - Boundary condition of X-Direction: Periodic inlet outlet
// X_velIn_ZeroGradOut
// - Boundary condition of X-Direction:  Constant velocity inlet, zero gradient Outlet.
// X_Wall
// - Boundary condition of X-Direction: Stationary Wall
// Z_Periodic
// - Boundary condition of Z-Direction: Periodic inlet outlet
// Z_Wall
// - Boundary condition of Z-Direction: Stationary or Moving Wall
// Y_Periodic
// - Boundary condition of Y-Direction: Periodic inlet outlet
// Y_Wall
// - Boundary condition of Y-Direction: Stationary or Moving Wall

//#define X_Periodic
#define X_Wall
//#define X_velIn_ZeroGradOut
//#define ZeroGradIn_ZeroGradOut


//#define Y_Periodic
#define Y_Wall

//#define Z_Periodic
#define Z_Wall

// Output visualisation file format:

//#define VTK_ASCII //Paraview file ASCII format
#define VTK_BINARY //Paraview file ASCII format
//#define CASE //Ensight file format

#if defined Flow_With_Thermal_Effect
// Temperature boundary condition for different walls
// Each wall can have only ONE of the boundary conditions
//#define Left_Fixed_Temp //left wall (X=0)
//#define Left_Adiabatic //left wall (X=0)
//#define Right_Fixed_Temp //right wall (X=Nx-1)
//#define Right_Adiabatic //right wall (X=Nx-1)
//#define Bottom_Fixed_Temp //bottom wall (Y=0)
//#define Bottom_Adiabatic //bottom wall (Y=0)
//#define Top_Fixed_Temp //Top wall (Y=Ny-1)
//#define Top_Adiabatic //Top wall (Y=Ny-1)
//#define Back_Fixed_Temp //back wall (Z=0)
//#define Back_Adiabatic //back wall (Z=0)
//#define Front_Fixed_Temp //Front wall (Z=Nz-1)
//#define Front_Adiabatic //Front wall (Z=Nz-1)
#endif

#if defined Flow_With_Thermal_Effect
#define Force_Thermal_1
//#define Force_Thermal_2

#define Equilibrium_Thermal_1  //The approach that equilibrium distribution function is defined: like paper Mohamad, A. A., M. El-Ganaoui, and R. Bennacer. "Lattice Boltzmann simulation of natural convection in an open ended cavity." International Journal of Thermal Sciences 48, no. 10 (2009): 1870-1875.
//#define Equilibrium_Thermal_2    //The approach that equilibrium distribution function is defined: Lai, Feng-Hsiang, and Yue-Tzu Yang. "Lattice Boltzmann simulation of natural convection heat transfer of Al 2 O 3/water nanofluids in a square enclosure." International Journal of Thermal Sciences 50, no. 10 (2011): 1930-1941.
#endif

#if defined Flow_With_Particle
/// Simulation types
// Exactly one of the following options has to be defined.
// STATIONARY_SPHERE
// - for simulation of, e.g., Karman vortex street
// - the sphere is kept in space, it does not move
// - In case of X_velIn_ZeroGradOut condition, Just Fixed particles (STATIONARY_SPHERE) are allowed
// STATIONARY_SPHEROID
// - the spheroid is kept in space, it does not move
// MOVING_SPHERE
// - for simulation of a moving sphere.
// - If you choose Direct_Forcing then the sphere CANNOT deform but moves in the fluid.
// - the sphere moves along with the flow.
// MOVING_SPHEROID
// - for simulation of a moving spheroid.

//#define STATIONARY_SPHERE
//#define STATIONARY_SPHEROID
//#define MOVING_SPHERE
#define MOVING_SPHEROID

// SDF: Single-Direct-Forcing.
// MDF: Multi-Direct-Forcing. Do NOT use it when there is NO particle inside the flow.
// use MDF in case of Flow_Without_Thermal_Effect
#define SDF
//#define MDF
const int mdfstep = 10;//Number of MDF steps. Adjust this line

//For description of added mass refer to page 88 and 89 of Kang PhD Thesis (2010) and Eq.18: Z.-G. Feng, E.E. Michaelides, Computers & Fluids 38 (2009) 370381
//#define Without_added_mass
#define With_added_mass


// -In Kang PhD thesis and Feng, Michaelides (paper 841, 2009) it is suggested to use corrected radius with diffuse interface scheme.

//#define Corrected_Radius

// Delta Function choice
// - You can choose Original delta function which was first developed by Krüger. //Removed at this version
// - You can choose 2 point function: Eq 67 Kang thesis (804)
// - You can choose 4 point function: Eq 68 Kang thesis (804)
// - You can choose Kernel function: Eq 17 Feng (806) Z.-G. Feng, E.E. Michaelides / Journal of Computational Physics 202 (2005) 2051
// - You can choose Modified_Two_Point_Delta: Yang, Xiaolei, et al. "A smoothing technique for discrete delta functions with application to immersed boundary method in moving boundary simulations." Journal of Computational Physics 228.20 (2009): 7821-7836.

//#define Two_Point_Delta
#define Four_Point_Delta
//#define Kernel_Delta
//#define Modified_Two_Point_Delta

// Collision model
// - select collision model Particle-Particle or Particle-Wall collision. (Just for sphere shape and just Direct forcing).
//#define Sphere_Collision_1_Feng   // Feng ZG, Michaelides 10E. Journal of Computational Physics,202 (2005), 20-51.
//#define Sphere_Collision_2_Feng_Verlet   // Use verlet list and avoid checking all particles. Feng ZG, Michaelides 10E. Journal of Computational Physics,202 (2005), 20-51.
#define Sphere_Collision_3_Feng_Lub_Verlet    //Lubrication and Feng repulsive models. Lubrication based on paper 952.

// Collision model for Spheroid
//// - select collision model for more than 1 moving Spheroidal particles or particle-wall collision.
//#define Spheroid_Collision_1_Feng
//#define Spheroid_Collision_2_Feng_Verlet
#define Spheroid_Collision_3_Feng_Lub_Verlet


#if defined Flow_With_Thermal_Effect
#define Particle_Fixed_Temperature
//#define Particle_Varying_Temperature
#endif

// - Random_Arrangement_ON means that location of 2nd until last particle is randomly defined by the program itself.
#define Random_Arrangement_OFF
//#define Random_Arrangement_ON

#endif


// Activate "Recovery_ON" if you want to continue the program from a recovery step (for example if the program has stopped due to electercity failure but you do NOT want to start from 1st iteration)
// In case of "Recovery_ON", initial values for fluid and particles are read from "data_fluid_recovery.dat" (See Initialize subroutine)
// Activate "Recovery_OFF" if you want to start the program from the 1st iteration.
// Only one of Recovery_OFF or Recovery_ON must be active.

#define Recovery_OFF
//#define Recovery_ON

#define ASCII_Recovery_File
//#define BINARY_Recovery_File

/// Fluid/lattice properties

const int Nx = 100; // number of lattice nodes along the x-axis (periodic)
const int Ny = 100; // number of lattice nodes along the y-axis (including two wall nodes)
const int Nz = 100; // number of lattice nodes along the z-axis (including two wall nodes)
const double tau = 0.8; // relaxation time
const int t_num = 500; // number of time steps (running from 1 to t_num)
const int t_data = 10; // disk write time step (data will be written to the disk every t_data step)
const int t_vtk = 100; // VTK write time step
const int t_info = 10; // info time step (screen message will be printed every t_info step)
const int t_time = 20; // Report elapsed time on screen
const int t_recovery = 200;
const double gravity = 0.000;// force density due to gravity (in positive x-direction). Just for flow without particle
const double wall_vel_bottom = 0.0; //velocity of the bottom wall(Y=0) (in positive x-direction)
const double wall_vel_top = 0.0; // velocity of the top wall(Y=Ny-1) (in positive x-direction)
const double wall_vel_back = 0.0;//velocity of the back wall(Z=0) (in positive x-direction)
const double wall_vel_front = 0.0;//velocity of the front wall(Z=Nz-1) (in positive x-direction)
const double inletvel = 0.0;//Inlet velocity. Just used in case of #define X_velIn_ZeroGradOut

#if defined Flow_With_Particle
							/// ******** PARTICLE PROPERTIES *********
const int num_particles = 2; // number of particles.
double particle_radius  = 5.0; // particle radius
const double particle_center_x = 65; // center position (x-component) of first particle
const double particle_center_y = 10.0; // center position (y-component) of first particle
const double particle_center_z = 10.0; // center position (z-component) of first particle
const double particle_gravity_x = 0.0004;//This force is applied to particles according to density ratio
const double particle_gravity_y = 0.00;//This force is applied to particles according to density ratio
const double particle_gravity_z = 0.00;//This force is applied to particles according to density ratio
const double den_ratio = 2; //Particle to fluid density ratio. Note that in below simulations density of fluid is assumed to be 1.0
const double offset_x = 16; //Distance of particles from each other in x direction (For Random_Arrangement_OFF).
const double offset_y = 0; //Distance of particles from each other in y direction (For Random_Arrangement_OFF).
const double offset_z = 0; //Distance of particles from each other in z direction (For Random_Arrangement_OFF).

double particle_radius_2 = 4.0; //Just used for "SPHEROID": radius along the Y-axis and Z-axis. The equation of the Spheroid is: (X^2/particle_radius^2) + (Y^2/particle_radius_2^2) + (Z^2/particle_radius_2^2) = 1. must NOT be exactly equal to particle_radius
#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
const int num_strip = 30; // Manual adjustment of total number of strips.
const double cir_factor = 2.5; //This parameter adjusts the number of points on each strip. Recommended: 0.6-2.0
#endif

							   //Particle collision parameters
const double r_buffer = 17;//Preferably = (r_cutoff+max(particle_radius,particle_radius_2))
const double r_cutoff = 12;//Preferably >= [2*max(particle_radius,particle_radius_2)+collision_threshold]
const double cij = 0.4;
const double ep = 2.0;//stifness factor for particle-particle collision
const double Ep = 0.1;//stifness factor for particle-particle collision usually smaller than ep
const double epw = 2.0;//stifness factor for particle-wall collision
const double Epw = 0.1;//stifness factor for particle-wall collision usually smaller than epw
const double th = 0.4;// Thereshold for particle-particle collision
const double thw = 1;// Thereshold for particle-wall collision
const double th_lub = 2.0;// Lubrication force thereshold (larger than th)
const double thw_lub = 2.0;// Lubrication force thereshold (larger than thw)
#endif

#if defined Flow_With_Thermal_Effect
const double specific_heat_ratio = 1.0; //Particle to fluid specific heat ratio.
const double prandtl = 0.7; //prandtl number
double gbeta = 0.00007335681; //In NATURAL convection: Ra=Gr*Pr=gbeta*dT*H^3/(alpha_t*nu)==>gbeta=alpha_t*nu*Ra/(dT*H^3). If natural convection is not important, set gbeta=0.0
double Tref = 0.0; //reference teperature in case of natural convection

#if defined Particle_Fixed_Temperature
double particle_temperature = 1.0; //Particle temperature must be constant
#endif
#if defined Particle_Varying_Temperature
double particle_temperature = 1.0; //Here temperature is the particle INITIAL temperature
double particle_heat_source = 0.0; //Heat source in the total volume of the particle (e.g. watts)
#endif

#if defined Left_Fixed_Temp //left wall (X=0)
double temp_left = 0.0;
#endif
#if defined Right_Fixed_Temp //right wall (X=Nx-1)
double temp_right = 0.0;
#endif
#if defined Bottom_Fixed_Temp //bottom wall (Y=0)
double temp_bot = 0.0;
#endif
#if defined Top_Fixed_Temp //Top wall (Y=Ny-1)
double temp_top = 0.0;
#endif
#if defined Back_Fixed_Temp //back wall (Z=0)
double temp_back = 0.0;
#endif
#if defined Front_Fixed_Temp //Front wall (Z=Nz-1)
double temp_front = 0.0;
#endif
#endif


//////////// Do NOT change  ///////////
#if defined Flow_Without_Thermal_Effect // Do NOT change
// Do NOT change 3 below lines
double gbeta = 0.0;
double Tref = 0.0; //reference teperature in case of natural convection
double particle_temperature = 0.0; //Particle temperature must be constant
#endif

#if defined Flow_With_Particle// Do not change
#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
int particle_num_nodes = int(M_PI / 3.*(3 * SQ(2 * particle_radius) + 1)); // Uhlmann (2005) recommendation for number of surface nodes for sphere
																		   //int particle_num_nodes = 2500; // Manual adjustment of number of surface nodes for sphere
#endif

#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
int particle_num_nodes = 6000; // Important. This parameter is just for initial adjustment. The actual particle_num_nodes for MOVING_SPHEROID will be specified in the main program.
#endif

std::vector<std::vector<int>> verlet_list(num_particles, std::vector<int>(5));//5 will change later
#endif
																			  /// *****************
																			  /// DECLARE VARIABLES
																			  /// *****************

																			  // The following code should not be modified when it is first used.

const double omega = 1. / tau; // relaxation frequency (inverse of relaxation time)

double ****pop, ****pop_old; // LBM populations (old and new)
double ***density; // fluid density
double ***velocity_x; // fluid velocity (x-component)
double ***velocity_y; // fluid velocity (y-component)
double ***velocity_z; // fluid velocity (z-component)
double ***force_x; // fluid force (x-component)
double ***force_y; // fluid force (y-component)
double ***force_z; // fluid force (z-component)
double ***mdf_force_x; // fluid force for MDF (x-component)
double ***mdf_force_y; // fluid force for MDF(y-component)
double ***mdf_force_z; // fluid force for MDF(z-component)
double ***mdf_vel_x; // vel_x for MDF(x-component)
double ***mdf_vel_y; // vel_y for MDF(y-component)
double ***mdf_vel_z; // vel_y for MDF(z-component)
int ***is_solid;
int* triangleset1;// Location of first triangle node. Just used for Particle_VTK to make a surface
int* triangleset2;// Location of first triangle node. Just used for Particle_VTK to make a surface
int* triangleset3;// Location of first triangle node. Just used for Particle_VTK to make a surface
int number_triangles;//number of triangles created for vtk file of particle
double pop_eq[19]; // equilibrium populations
const double weight[19] = { 1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36. }; // lattice weights

#if defined Flow_With_Thermal_Effect
const double alpha_t = (tau - 0.5) / 3. / prandtl; // Thermal diffusion coefficient
const double tau_t = 3.* alpha_t + 0.5; // relaxation time of temperature
const double omega_t = 1. / tau_t; // relaxation frequency (inverse of relaxation time for temperature)
#endif
double ***temperature; // fluid temperature
double ****pop_t, ****pop_old_t; // LBM populations of temperature (old and new)
double pop_eq_t[19]; // equilibrium populations of temperature
double ***force_thermal; // force term for thermal lbm
						 /// ******************
						 /// PARTICLE STRUCTURE
						 /// ******************

						 // The following code handles the object immersed in the flow.
						 // In the present implementation, only a single object can be put into the flow.

						 /// Structure for surface nodes
						 // Each node has a current x- and y-position and a reference x- and y-position.
#if defined Flow_With_Particle
struct node_struct {

	/// Constructor
	node_struct() {
		x = 0;
		y = 0;
		z = 0;
		vel_x = 0;
		vel_y = 0;
		vel_z = 0;
		force_x = 0;
		force_y = 0;
		force_z = 0;
		den = 0;
		temperature = particle_temperature;
		force_thermal = 0;
	}

	/// Elements
	double x; // current x-position
	double y; // current y-position
	double z; // current z-position
	double vel_x; // node velocity (x-component)
	double vel_y; // node velocity (y-component)
	double vel_z; // node velocity (z-component)
	double force_x; // node force (x-component)
	double force_y; // node force (y-component)
	double force_z; // node force (z-component)
	double den;//Density at boundary point
	double force_thermal; // node force due to thermal effect (it is not in Newton)
	double temperature; // node force due to thermal effect (it is not in Newton)
};

struct center_struct {

	/// Constructor
	center_struct() {
		x = 0;
		y = 0;
		z = 0;
		vel_x = 0;
		vel_y = 0;
		vel_z = 0;
		force_x = 0;
		force_y = 0;
		force_z = 0;
		den = 0;
		temperature = particle_temperature;
		force_thermal = 0;
		omgX = 0.0;
		omgY = 0.0;
		omgZ = 0.0;
		previous_vel_x = 0.0;
		previous_vel_y = 0.0;
		previous_vel_z = 0.0;
		tetaX = 0.0;
		tetaY = 0.0;
		tetaZ = 0.0;
		tempx[0] = 0.0;
		tempx[1] = 0.0;
		tempy[0] = 0.0;
		tempy[1] = 0.0;
		tempz[0] = 0.0;
		tempz[1] = 0.0;
		temp_omgX[0] = 0.0;
		temp_omgX[1] = 0.0;
		temp_omgY[0] = 0.0;
		temp_omgY[1] = 0.0;
		temp_omgZ[0] = 0.0;
		temp_omgZ[1] = 0.0;
		temp_part = 0.0;
		temp_temperat[0] = 0;
		temp_temperat[1] = 0;
		delta_teta_X = 0;
		delta_teta_Y = 0;
		delta_teta_Z = 0;
		Omg_prime_X = 0;
		Omg_prime_Y = 0;
		Omg_prime_Z = 0;
		Qu0 = 0;
		Qu1 = 0;
		Qu2 = 0;
		Qu3 = 0;
		displacement_verlet = 0;
		deter = 0;
		quater[0] = 0;
		quater[1] = 0;
		quater[2] = 0;
		quater[3] = 0;
		quater[4] = 0;
		quater[5] = 0;
		quater[6] = 0;
		quater[7] = 0;
		quater[8] = 0;
		quater[9] = 0;
		f3x = 0;
		f3y = 0;
		f3z = 0;
		f4x = 0;
		f4y = 0;
		f4z = 0;
		t1_x = 0.0;// Torque on particle center from fluid
		t1_y = 0.0;// Torque on particle center from fluid
		t1_z = 0.0;// Torque on particle center from fluid
		t2_x = 0.0;// Torque on particle center from fluid
		t2_y = 0.0;// Torque on particle center from fluid
		t2_z = 0.0;// Torque on particle center from fluid
		t3_x = 0.0;
		t3_y = 0.0;
		t3_z = 0.0;
	}

	/// Elements
	double x; // current x-position
	double y; // current y-position
	double z; // current z-position
	double vel_x; // center velocity (x-component)
	double vel_y; // center velocity (y-component)
	double vel_z; // center velocity (z-component)
	double force_x; // center force (x-component)
	double force_y; // center force (y-component)
	double force_z; // center force (z-component)
	double den;//Density at boundary point
	double omgX; // center rotational speed
	double omgY; // center rotational speed
	double omgZ; // center rotational speed
	double temperature; // node force due to thermal effect (it is not in Newton)
	double force_thermal; // center force due to thermal effect (it is not in Newton)
	double previous_vel_x; // Velocity x in previous time step
	double previous_vel_y; // Velocity y in previous time step
	double previous_vel_z; // Velocity z in previous time step
	double tetaX; // Angle x
	double tetaY; // Angle y
	double tetaZ; // Angle z
	double tempx[2];// Translational speed of particle at two previous time steps
	double tempy[2];// Translational speed of particle at two previous time steps
	double tempz[2];// Translational speed of particle at two previous time steps
	double temp_omgX[2]; // Rotational speed of particle at two previous time steps
	double temp_omgY[2]; // Rotational speed of particle at two previous time steps
	double temp_omgZ[2]; // Rotational speed of particle at two previous time steps
	double temp_part; // At first iteration it is Initial particle temperature and then take the value of particle temperature at pervious iteration
	double temp_temperat[2]; // Temperature of particle at two previous time steps
	double delta_teta_X, delta_teta_Y, delta_teta_Z; // Particle change in rotation angle
	double Omg_prime_X, Omg_prime_Y, Omg_prime_Z; //ELLIPSOID only: Rotational velocity in body-fixed coordinate
	double Qu0, Qu1, Qu2, Qu3; //Quaternions
	double displacement_verlet;//Displacement of particle sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2) between two Verlet list update
	double deter;//determinant of transformations matrix for spheroid
	double quater[10];
	double f3x;// Total Collision force on particle center due to particle-particle collision
	double f3y;
	double f3z;
	double f4x;// Total Collision force on particle center due to particle-wall collision
	double f4y;
	double f4z;
	double t1_x; // Torque on particle center from fluid
	double t1_y; // Torque on particle center from fluid
	double t1_z; // Torque on particle center from fluid
	double t2_x; // Torque on particle center from fluid
	double t2_y; // Torque on particle center from fluid
	double t2_z; // Torque on particle center from fluid
	double t3_x; // Torque on particle center due to collision with other particles or walls--> Just for Ellipsoid
	double t3_y; // Torque on particle center due to collision with other particles or walls--> Just for Ellipsoid
	double t3_z; // Torque on particle center due to collision with other particles or walls--> Just for Ellipsoid

};

struct particle_struct {

	/// Constructor
	particle_struct() {
		//num_nodes = particle_num_nodes;
		radius = particle_radius;
		radius_2 = particle_radius;
#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
		radius_2 = particle_radius_2;
#endif
		center = new center_struct[1];
		center[0].x = 0;    //Will be initialized and changed later. particle_center_x;
		center[0].y = 0;    //Will be initialized and changed later.particle_center_y;
		center[0].z = 0;    //Will be initialized and changed later.particle_center_z;
		node = new node_struct[particle_num_nodes];
		// The initial position of the object is set in main().
	}

	/// Elements
	//int num_nodes; // number of surface nodes
	double radius, radius_2; // object radius
	center_struct *center; // center node
	node_struct *node; // list of nodes
};
#endif
/// *****************
/// DECLARE FUNCTIONS
/// *****************

// The following functions are used in the simulation code.

void initialize(int &); // allocate memory and initialize variables
void LBM(int); // perform LBM operations
void momenta(); // compute fluid density and velocity from the populations
void equilibrium(double, double, double, double); // compute the equilibrium populations from the fluid density and velocity
void write_fluid_vtk(int); // write the fluid state to the disk as VTK file
void write_solid_vtk(int, int);// write SOLID points to a VTK file.
void write_fluid_data(int); // write fluid data to the disk (velocity/density)
void write_parametes();
void write_fluid_recovery(int); // write fluid data to the disk for case of recovery if program crashes
void read_fluid_recovery(int &); // read fluid initial values from a file if "Recovery_ON" is defined
double meq(double, double, double, double, int);// Equilibrium momenta
std::string numberToString(int number); //Convert number to string


void swap_endian(char *, size_t);
void swap_endian(char *, char *, size_t);
void write_bigendian(std::ofstream &, char *, size_t, size_t);
template <typename T>
void write_bigendian(std::ofstream &, T  *, size_t);

void equilibrium_t(double, double, double, double, double); // compute the equilibrium populations of temperature from the fluid density and velocity
void LBM_t(int); // perform LBM operations for temperature
void momenta_t(); // compute fluid temperature from the populations

#if defined Flow_With_Particle
void spread(particle_struct &); // spread node forces to fluid lattice
void interpolate(particle_struct &); // interpolate node velocities from fluid velocity
void write_particle_vtk(int, particle_struct, int); // write the particle state to the disk as VTK file
void write_particle_data(int, particle_struct, int, double, double, double); // write particle data to the disk (drag/lift, center position)
void rungekutta(double &, double &, double &, double, double, double, double, double, double, double &, double &, double &, double &); // Runge-Kutta at each time step for the spheroid rotation
void quaternion(double, double, double, double, int, double[]);//Calculation of quaternion matrix elements. Eq2.13 Article 876 Huang 2012
double determinant(int, double[]);
void comp_rot_speed(double, double &, double &, double &, int, double[], double, double, double);//update rotational speed. Just for SPHEROID
void collision_sphere(particle_struct[num_particles], int);//Collision forces for moving circular objects
void collision_spheroid(particle_struct[num_particles], int, double[], double[], double[]);//Collision forces for moving Spheroidal objects
void check_zero(double &, double);
void update_rot_spheroid(double, double &, double &, double &, double[], double, double, double, double, double, double);//update initial position of nodes if Euler_angles are not initially zero. Just for SPHEROID
void normal_vec_spheroid(double, double &, double &, double &, int, double[], double, double, double);//Calculates NORMAL vector at each surface point (in INERTIAL coordinate) of spheroid . Used to calculate collision forces
void update_verlet(particle_struct[num_particles]);
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MAIN
//-----------------------------------------------------------------------------------------------------------------------------
//
// This is the main function, containing the simulation initialization and the simulation loop.

int main(int argc, char** argv) {

	clock_t elapsed_time, elapsed_time1;
	elapsed_time = clock();
	/// ************
	/// PREPARATIONS
	/// ************
	int recovery_step; //If "Recovery_ON" is activated the program starts from step=recovery_step+1
	initialize(recovery_step); // allocate memory and initialize variables
	write_parametes();

#if defined Flow_With_Particle
	particle_struct particle[num_particles]; // create immersed object
	double ck, tetak, fik, fikold;//parameters for the uniform distribution of points on the sphere
	double particle_mass, particle_area, mom_inertia_x, mom_inertia_y, mom_inertia_z, parvol;
	// double area[particle_num_nodes];//area belonging to each surafce node.
	double* area = 0;//area belonging to each surafce node.
	area = new double[particle_num_nodes];
	double aa, cc, hh, zz;
	int margin;

	margin = 0000;//This parameter is helpful for STATIONARY_SPHERE and STATIONARY_SPHEROID. You can change this line if you want the particles to be
				  //... stationary until a certain iteration (margin) and the starts motion.

	double disp_x, disp_y, disp_z;
	double disp_max_1 = 0;// Max Displacement of particles
	double disp_max_2 = 0;// Second max of Displacement of particles

#if defined STATIONARY_SPHERE || defined STATIONARY_SPHEROID
	margin = t_num + 1;//Do not change this line
#endif
#endif

#if defined Recovery_OFF && defined Flow_With_Particle
					   //***************************************************


					   /// ****************************************
					   /// PARTICLE MASS AND MOMENT OF INERTIA
					   /// ****************************************
#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
	parvol = 4.0 / 3.0 * M_PI * particle_radius * particle_radius * particle_radius;//Particle volume
	particle_mass = (den_ratio*1.0) * parvol;//Particle mass
	particle_area = M_PI * (12.0 * SQ(particle_radius) + 1.0) / 3.0;	//Eq. 2.14 Lucci et. al. J. Fluid Mech. (2010), vol. 650, pp. 555.
	mom_inertia_x = 2.0 / 5.0 * particle_mass * SQ(particle_radius);//Moment of inertia
	mom_inertia_y = 2.0 / 5.0 * particle_mass * SQ(particle_radius);
	mom_inertia_z = 2.0 / 5.0 * particle_mass * SQ(particle_radius);
#endif

#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
	parvol = 4.0 / 3.0 * M_PI * particle_radius * particle_radius_2 * particle_radius_2;//Particle volume
	particle_mass = (den_ratio*1.0) * parvol;
	cc = particle_radius;
	aa = particle_radius_2;
	hh = sqrt(MAX((aa - cc), (cc - aa))*(aa + cc) / pow(cc, 4));
	zz = particle_radius * hh;
	if (cc > aa) {
		particle_area = 2 * M_PI * aa * 2 * (asin(zz) / 2. + (zz*sqrt(1 - SQ(zz))) / 2.) / hh;
	}
	else {
		particle_area = 2 * M_PI * aa * 2 * (asinh(zz) / 2. + (zz*sqrt(1 + SQ(zz))) / 2.) / hh;
	}
	//cout<<"particle area "<<particle_area<<endl;
	mom_inertia_x = 1.0 / 5.0 * particle_mass * (SQ(particle_radius_2) + SQ(particle_radius_2));
	mom_inertia_y = 1.0 / 5.0 * particle_mass * (SQ(particle_radius) + SQ(particle_radius_2));
	mom_inertia_z = 1.0 / 5.0 * particle_mass * (SQ(particle_radius) + SQ(particle_radius_2));
#endif

	//***************************************************
#endif  // End for Recovery_OFF && Flow_With_Particle

#if defined Corrected_Radius && defined Flow_With_Particle
    double corrected_radius, corrected_radius_2;
    corrected_radius   = pow((pow(particle_radius, 3)   + pow(particle_radius - 1.0, 3))   / 2.0, 1.0 / 3.0);
    corrected_radius_2 = pow((pow(particle_radius_2, 3) + pow(particle_radius_2 - 1.0, 3)) / 2.0, 1.0 / 3.0);

	particle_radius   = corrected_radius;//If corrected-radius is defined, from now on particle-radius represents corrected radius
#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
	particle_radius_2 = corrected_radius_2;//If corrected-radius is defined, from now on particle-radius_2 represents corrected radius
#endif
#endif


#if defined Flow_With_Particle
										   //***************************************************

										   /// ***********************************
										   /// CREATING OBJECTS (1st PARTICLE)
										   /// ***********************************
	if (num_particles > 0) {
		particle[0].center[0].x = particle_center_x;
		particle[0].center[0].y = particle_center_y;
		particle[0].center[0].z = particle_center_z;

#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
		// Distributing points on a sphere according to:
		//Eq. 2.10 , 2.11 Lucci et. al. J. Fluid Mech. (2010), vol. 650, pp. 555.
		for (int k = 0; k < particle_num_nodes; ++k) {
			ck = -1.0 + 2.0 * (double)k / (double)(particle_num_nodes - 1.0);
			tetak = acos(ck);
			if (k == 0 || k == particle_num_nodes - 1) {
				fik = 0.0;
				fikold = 0.0;
			}
			else
			{
				fik = fikold + 3.6 / sqrt(particle_num_nodes) / sqrt(1.0 - SQ(ck));
			}
			fikold = fik;
			particle[0].node[k].x = particle[0].center[0].x + particle_radius * cos(fik) * sin(tetak);
			particle[0].node[k].y = particle[0].center[0].y + particle_radius * sin(fik) * sin(tetak);
			particle[0].node[k].z = particle[0].center[0].z + particle_radius * cos(tetak);
			//Eq. 2.14 Lucci et. al. J. Fluid Mech. (2010), vol. 650, pp. 555.
			area[k] = M_PI * (12.0 * SQ(particle_radius) + 1.0) / 3.0 / (double)particle_num_nodes; // area belonging to a node
			if (particle[0].node[k].x > Nx - 0.001 || particle[0].node[k].x < 0.001) {
				printf("ERROR: Particle %d position is out of the domain \n", 1);
				return 0;
			}
			if (particle[0].node[k].y > Ny - 0.001 || particle[0].node[k].y < 0.001) {
				printf("ERROR: Particle %d position is out of the domain \n", 1);
				return 0;
			}
			if (particle[0].node[k].z > Nz - 0.001 || particle[0].node[k].z < 0.001) {
				printf("ERROR: Particle %d position is out of the domain \n", 1);
				return 0;
			}
		}
#endif

#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
		// Distributing points on a spheroid.
		int kkk, mmm;
		double perimeter;
		cc = particle_radius;
		aa = particle_radius_2;
		perimeter = M_PI * (3.0*(particle_radius + particle_radius_2) - sqrt((3.0*particle_radius + particle_radius_2)*(particle_radius + 3.0*particle_radius_2)));
		//num_strip = (int)(perimeter / 2.0 * particle_radius / particle_radius_2) + 2;//Total number of parallel strips (including top and bottom points)

		//double hei[num_strip];//height of each strip
		double* hei = 0;//height of each strip
		hei = new double[num_strip];

		int** nodes_in_hei;//It tells us at each strip height wich nodes are located (node number). Important: Just for VTK file
		nodes_in_hei = new int*[num_strip];
		for (int nps = 0; nps < num_strip; ++nps) {
			nodes_in_hei[nps] = new int[2000];//Certainly the number of nodes on one strip is less than 2000
		}

		int* po_per_str = 0;//Number of points per each strip
		po_per_str = new int[num_strip];

		//Important: In the following configuration, the particle_radius is initially along the X-axis and particle_radius_2 is initially along the Y and Z-axis. You can change it later by adjusting initial euler angles
		kkk = 0;
		particle[0].node[kkk].x = particle[0].center[0].x + particle_radius;
		particle[0].node[kkk].y = particle[0].center[0].y;
		particle[0].node[kkk].z = particle[0].center[0].z;
		kkk = kkk + 1;//kkk is the index for number of particles.
		hei[0] = particle[0].node[0].x;//hei is the elevation of each strip
		po_per_str[0] = 1;
		nodes_in_hei[0][0] = 0;

		for (int hf = 1; hf < num_strip - 1; hf++) {
			int ghk = 0;
			tetak = (double)(hf)* M_PI / (double)(num_strip - 1);

			////////// Method 1
			po_per_str[hf] = int(cir_factor * num_strip * sin(tetak));
			//////// End of Method 1

			////////// Method 2
			//double dist_to_cent = particle_radius_2*sqrt(2)*sin(tetak);
			//po_per_str[hf] = abs(int(2 * M_PI*dist_to_cent*cir_factor));
			//////// End of Method 2

			if (po_per_str[hf] < 2) {
				po_per_str[hf] = 2;
			}
			fik = 0.0;
			for (int hd = 1; hd <= po_per_str[hf]; hd++) {
				fik = fik + 2. * M_PI / po_per_str[hf];
				particle[0].node[kkk].y = particle[0].center[0].y + particle_radius_2 * cos(fik) * sin(tetak);
				particle[0].node[kkk].z = particle[0].center[0].z + particle_radius_2 * sin(fik) * sin(tetak);
				particle[0].node[kkk].x = particle[0].center[0].x + particle_radius * cos(tetak);
				hei[hf] = particle[0].node[kkk].x;
				nodes_in_hei[hf][ghk] = kkk;
				//cout << hf << " " << nodes_in_hei[hf][ghk] << " " << ghk << endl;
				ghk = ghk + 1;
				kkk = kkk + 1;
			}
		}
		particle[0].node[kkk].x = particle[0].center[0].x - particle_radius;//Last point
		particle[0].node[kkk].y = particle[0].center[0].y;
		particle[0].node[kkk].z = particle[0].center[0].z;
		po_per_str[num_strip - 1] = 1;
		nodes_in_hei[num_strip - 1][0] = kkk;//Index of point that is located on this strip (last strip)
		hei[num_strip - 1] = particle[0].node[kkk].x;//hei is the elevation of each strip
		particle_num_nodes = kkk + 1; //The actual particle_num_nodes for spheroid
									  //double area[particle_num_nodes];//If you activate this line you would have wrong results (Just remember not to define a parameter twice)

		mmm = 0;
		double sum, h1, h2, a1, a2;
		sum = 0.0;//Total area. It is to make a comparison with a Full Spheroid formula
		for (int hf = 0; hf < num_strip - 1; hf++) {
			/////////tetak = (double)(hf)* M_PI / (double)(num_strip - 1);
			for (int hd = 1; hd <= po_per_str[hf]; hd++) {
				if (hf != 0) {
					h1 = (hei[hf] - particle[0].center[0].x + hei[hf - 1] - particle[0].center[0].x) / 2. * hh;
				}
				else {
					h1 = cc *hh;//(hei[hf]-particle[0].center[0].x ) * hh;
				}
				h2 = (hei[hf] - particle[0].center[0].x + hei[hf + 1] - particle[0].center[0].x) / 2. * hh;
				if (cc > aa) {
					a1 = M_PI * aa * 2 * (asin(h1) / 2. + (h1*sqrt(1 - SQ(h1))) / 2.) / hh / po_per_str[hf];
					a2 = M_PI * aa * 2 * (asin(h2) / 2. + (h2*sqrt(1 - SQ(h2))) / 2.) / hh / po_per_str[hf];
				}
				else
				{
					a1 = M_PI * aa * 2 * (asinh(h1) / 2. + (h1*sqrt(1 + SQ(h1))) / 2.) / hh / po_per_str[hf];
					a2 = M_PI * aa * 2 * (asinh(h2) / 2. + (h2*sqrt(1 + SQ(h2))) / 2.) / hh / po_per_str[hf];
				}
				area[mmm] = abs(a1 - a2);
				//cout<<area[mmm]<<"  ";
				sum = sum + area[mmm];
				mmm = mmm + 1;
			}
		}
		area[mmm] = area[0];//The last point
		sum = sum + area[mmm];
		//cout<<"total area: "<<sum<<endl;

		/////////////////////////////// Below part is just for creating surface triangles for Particle_VTK
		for (int i = 0; i < po_per_str[1]; i++) {
			triangleset1[i] = 0;
			triangleset2[i] = nodes_in_hei[1][i];
			triangleset3[i] = nodes_in_hei[1][(i + 1) % po_per_str[1]];
		}
		// For Particle_vtk we want to create a surface by triangles.
		kkk = 1;
		int tindex = po_per_str[1];
		double dist_node;
		double min_dist;
		int min_node, ipp;
		double max_dist = -2;
		//cout << "SSSSSS " << endl;
		int ik;

		for (int hf = 1; hf < num_strip - 1; hf++) {
			for (int ij = 0; ij < po_per_str[hf]; ij++) {
				ik = hf - 1;
				min_dist = 2 * MAX(particle_radius, particle_radius_2);
				for (int ip = 0; ip < po_per_str[ik]; ip++) {
					dist_node = sqrt(SQ(particle[0].node[nodes_in_hei[hf][ij]].x - particle[0].node[nodes_in_hei[ik][ip]].x) + SQ(particle[0].node[nodes_in_hei[hf][ij]].y - particle[0].node[nodes_in_hei[ik][ip]].y) + SQ(particle[0].node[nodes_in_hei[hf][ij]].z - particle[0].node[nodes_in_hei[ik][ip]].z));
					if (ip == 0) {
						ipp = 0;
					}

					if (dist_node<min_dist) {
						min_node = nodes_in_hei[ik][ip];
						ipp = ip;
						min_dist = dist_node;
					}
				}

				triangleset2[tindex] = nodes_in_hei[ik][ipp];
				triangleset1[tindex] = kkk;
				triangleset3[tindex] = nodes_in_hei[ik][(ipp + 1) % po_per_str[ik]];
				//cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
				tindex = tindex + 1;
				triangleset2[tindex] = nodes_in_hei[ik][ipp];
				triangleset1[tindex] = kkk;
				triangleset3[tindex] = nodes_in_hei[ik][(ipp - 1 + po_per_str[ik]) % po_per_str[ik]];
				//cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
				tindex = tindex + 1;

				ik = hf + 1;
				min_dist = 2 * MAX(particle_radius, particle_radius_2);
				for (int ip = 0; ip < po_per_str[ik]; ip++) {
					dist_node = sqrt(SQ(particle[0].node[nodes_in_hei[hf][ij]].x - particle[0].node[nodes_in_hei[ik][ip]].x) + SQ(particle[0].node[nodes_in_hei[hf][ij]].y - particle[0].node[nodes_in_hei[ik][ip]].y) + SQ(particle[0].node[nodes_in_hei[hf][ij]].z - particle[0].node[nodes_in_hei[ik][ip]].z));
					if (ip == 0) {
						ipp = 0;
					}

					if (dist_node<min_dist) {
						min_node = nodes_in_hei[ik][ip];
						ipp = ip;
						min_dist = dist_node;
					}
				}

				triangleset2[tindex] = nodes_in_hei[ik][ipp];
				triangleset1[tindex] = kkk;
				triangleset3[tindex] = nodes_in_hei[ik][(ipp + 1) % po_per_str[ik]];
				//cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
				tindex = tindex + 1;
				triangleset2[tindex] = nodes_in_hei[ik][ipp];
				triangleset1[tindex] = kkk;
				triangleset3[tindex] = nodes_in_hei[ik][(ipp - 1 + po_per_str[ik]) % po_per_str[ik]];
				//cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
				tindex = tindex + 1;
				kkk = kkk + 1;
			}
		}
		number_triangles = tindex;

		for (int i = 0; i < po_per_str[num_strip - 2]; i++) {
			triangleset1[tindex + i] = particle_num_nodes - 1;
			triangleset2[tindex + i] = nodes_in_hei[num_strip - 2][i];
			triangleset3[tindex + i] = nodes_in_hei[num_strip - 2][(i + 1) % po_per_str[1]];
			number_triangles++;
		}
		//////////////////////////////////End of Particle_VTK part

		delete hei;
		delete po_per_str;
#endif //spheroid
	}

	/// ***********************************
	/// LOCATION OF 2nd UNTIL Nth PARTICLE
	/// ***********************************
#if defined Random_Arrangement_ON
	for (int i = 1; i < num_particles; ++i) {
		int part_position;
		srand(5323);
		double max_dist;
#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
		max_dist = particle[0].radius;
#endif
#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
		max_dist = MAX(particle[0].radius_2, particle[0].radius);
#endif
		do {
			part_position = TRUE;
			particle[i].center[0].x = rand() % Nx;
			particle[i].center[0].y = rand() % Ny;
			particle[i].center[0].z = rand() % Nz;
			printf("Particle i is %d and x, y, z are %5.3f , %5.3f, %5.3f \n", i + 1, particle[i].center[0].x, particle[i].center[0].y, particle[i].center[0].z);

			for (int j = 0; j < num_particles; ++j) {
				if (j != i) {
					double Dij = sqrt(SQ(particle[i].center[0].x - particle[j].center[0].x) + SQ(particle[i].center[0].y - particle[j].center[0].y) + SQ(particle[i].center[0].z - particle[j].center[0].z));
					//#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
					if (Dij < 2 * max_dist + 1) {
						part_position = FALSE;
						printf("Error, contact with particle %d \n", j);
						//	break;
					}
					//#endif
				}
			}

			if (particle[i].center[0].x > Nx - 1 - max_dist - 1 || particle[i].center[0].x < max_dist + 1) {
				part_position = FALSE;
				cout << "Error 1" << endl;
			}
			if (particle[i].center[0].y > Ny - 1 - max_dist - 1 || particle[i].center[0].y < max_dist + 1) {
				part_position = FALSE;
				cout << "Error 2" << endl;
			}
			if (particle[i].center[0].z > Nz - 1 - max_dist - 1 || particle[i].center[0].z < max_dist + 1) {
				part_position = FALSE;
				cout << "Error 3" << endl;
			}
		} while (part_position == FALSE);

		for (int n = 0; n < particle_num_nodes; ++n) {
			particle[i].node[n].x = particle[0].node[n].x - particle[0].center[0].x + particle[i].center[0].x;
			particle[i].node[n].y = particle[0].node[n].y - particle[0].center[0].y + particle[i].center[0].y;
			particle[i].node[n].z = particle[0].node[n].z - particle[0].center[0].z + particle[i].center[0].z;
		}
	}
#endif
#if defined Random_Arrangement_OFF
	for (int i = 1; i < num_particles; ++i) {
		particle[i].center[0].x = particle_center_x + i * offset_x;//Location of other particles is adjusted according to location of first one. You can change it.
		particle[i].center[0].y = particle_center_y + i * offset_y;
		particle[i].center[0].z = particle_center_z + i * offset_z;
		for (int n = 0; n < particle_num_nodes; ++n) {
			particle[i].node[n].x = particle[0].node[n].x - particle[0].center[0].x + particle[i].center[0].x;
			particle[i].node[n].y = particle[0].node[n].y - particle[0].center[0].y + particle[i].center[0].y;
			particle[i].node[n].z = particle[0].node[n].z - particle[0].center[0].z + particle[i].center[0].z;
			if (particle[i].node[n].x > Nx - 1 - 0.1 || particle[i].node[n].x < 0.1) {
				printf("ERROR: Particle %d position is out of the domain \n", i + 1);
				return 0;
			}
			if (particle[i].node[n].y > Ny - 1 - 0.1 || particle[i].node[n].y < 0.1) {
				printf("ERROR: Particle %d position is out of the domain \n", i + 1);
				return 0;
			}
			if (particle[i].node[n].z > Nz - 1 - 0.1 || particle[i].node[n].z < 0.1) {
				printf("ERROR: Particle %d position is out of the domain \n", i + 1);
				return 0;
			}
		}
	}
#endif

	//***************************************************
#endif  // End for Recovery_OFF && Flow_With_Particle

	/// ***********************************
	/// REPORT DERIVED PARAMETERS
	/// ***********************************
	/// Compute derived quantities
	const double D = Ny - 2; // inner channel diameter
	const double Dz = Nz - 2; // inner channel diameter
	const double nu = (tau - 0.5) / 3.0; // lattice viscosity
	const double umax = gravity / (2.0 * nu) * SQ(0.5 * D); // expected maximum velocity for Poiseuille flow without immersed object
	double Re;
	if (gravity != 0) {
		Re = D * 2.0 / 3.0*umax / nu; // Reynolds number for Poiseuille flow without immersed object.
	}

	/// Report derived parameters
	printf("Simulation Parameters \n");
	printf("===================== \n");
	printf("Nx                  = %d \n", Nx);
	printf("Ny                  = %d \n", Ny);
	printf("Nz                  = %d \n", Nz);
	printf("Viscosity           = %6.4f \n", nu);
	printf("Gravity             = %6.4f \n", gravity);


#if defined X_velIn_ZeroGradOut
	cout << "Inlet Velocity      = " << inletvel << endl;
#endif
#if defined	Flow_With_Particle
	printf("Number of Particles = %d \n", num_particles);
#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
	printf("Spheroid Particle  \n");
#endif
#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
	printf("Sphere Particle  \n");
#endif
#endif

	/// ********************************
	/// ********************************
	///         SIMULATION LOOP
	/// ********************************
	/// ********************************

#if defined Flow_With_Particle
	double fx_surf[num_particles], fy_surf[num_particles], fz_surf[num_particles]; // Total force on particle center
																				   //double tempx[num_particles][2], tempy[num_particles][2], tempz[num_particles][2]; // Translational speed of particle at two previous time steps
																				   //double temp_omgX[num_particles][2], temp_omgY[num_particles][2], temp_omgZ[num_particles][2]; // Rotational speed of particle at two previous time steps
																				   //double delta_teta_X[num_particles], delta_teta_Y[num_particles], delta_teta_Z[num_particles]; // Particle change in rotation angle
	double center_old_x, center_old_y, center_old_z; // Keeps particle center position before updating by Newton's equations
	int col_eval[num_particles];

	//double temp_part[num_particles]; // At first iteration it is Initial particle temperature and then take the value of particle temperature at pervious iteration
	double heat_surf[num_particles]; // Total heat on particle from fluid
									 //double temp_temperat[num_particles][2]; // Temperature of particle at two previous time steps

	double x_pos, y_pos, z_pos; //Temporary parameter

								//Next 11 lines parameters are just for MOVING_SPHEROID to model its rotation
								/////double Omg_prime_X[num_particles], Omg_prime_Y[num_particles], Omg_prime_Z[num_particles]; //Rotational velocity in body-fixed coordinate
								//double Qu0[num_particles], Qu1[num_particles], Qu2[num_particles], Qu3[num_particles]; //Quaternions
	double T_prime_X[num_particles], T_prime_Y[num_particles], T_prime_Z[num_particles];//Torques in body-fixed coordinate
	double euler_teta[num_particles], euler_fi[num_particles], euler_sai[num_particles]; //Euler angles
	double  sumq;
	double* x_pos_1; //Node x-location in BODY-fixed coordinate with respect to particle center. It does not change during time steps
	double* y_pos_1; //Node y-location in BODY-fixed coordinate with respect to particle center. It does not change during time steps
	double* z_pos_1; //Node z-location in BODY-fixed coordinate with respect to particle center. It does not change during time steps
	x_pos_1 = new double[particle_num_nodes];
	y_pos_1 = new double[particle_num_nodes];
	z_pos_1 = new double[particle_num_nodes];

	//***********************************

	/// ---------------------------------------------------------------------------------------------
	/// ----------------------------------- PARTICLE INITIALIZATION ---------------------------------
	/// ---------------------------------------------------------------------------------------------
	for (int np = 0; np < num_particles; ++np) {
		update_verlet(particle);
		center_old_x = 0.0;//position of particle center at previous iteration. used to update rotational position.
		center_old_y = 0.0;
		center_old_z = 0.0;
		particle[np].center[0].temp_omgX[0] = 0.0;
		particle[np].center[0].temp_omgX[1] = 0.0;
		particle[np].center[0].temp_omgY[0] = 0.0;
		particle[np].center[0].temp_omgY[1] = 0.0;
		particle[np].center[0].temp_omgZ[0] = 0.0;
		particle[np].center[0].temp_omgZ[1] = 0.0;
		particle[np].center[0].delta_teta_X = 0.0;//Change of Angle (between two iterarions)-Rotation aroud X axis
		particle[np].center[0].delta_teta_Y = 0.0;//Change of Angle (between two iterarions)-Rotation aroud Y axis
		particle[np].center[0].delta_teta_Z = 0.0;//Change of Angle (between two iterarions)-Rotation aroud Z axis
		particle[np].center[0].tempx[0] = 0.0;
		particle[np].center[0].tempx[1] = 0.0;
		particle[np].center[0].tempy[0] = 0.0;
		particle[np].center[0].tempy[1] = 0.0;
		particle[np].center[0].tempz[0] = 0.0;
		particle[np].center[0].tempz[1] = 0.0;

		particle[np].center[0].temp_part = particle_temperature; //Initial particle temperature

		particle[np].center[0].Omg_prime_X = 0.0;//Angular velocity in body-fixed coordinate. Just used for deformable-spheroid. Must be ZERO initially
		particle[np].center[0].Omg_prime_Y = 0.0;//Angular velocity in body-fixed coordinate. Just used for deformable-spheroid. Must be ZERO initially
		particle[np].center[0].Omg_prime_Z = 0.0;//Angular velocity in body-fixed coordinate. Just used for deformable-spheroid. Must be ZERO initially
		euler_teta[np] = 0.0; //Euler angles (radian). You can adjust them. Important: when all Euler angles are zero, particle_radius is along X-axis and particle_radius_2 is along Y and Z-axis.
		euler_fi[np] = 0.0; //Euler angles (radian). You can adjust them. Important: when all Euler angles are zero, particle_radius is along X-axis and particle_radius_2 is along Y and Z-axis.
		euler_sai[np] = 0.0; //Euler angles (radian). You can adjust them. Important: when all Euler angles are zero, particle_radius is along X-axis and particle_radius_2 is along Y and Z-axis.
							 //When euler_teta, euler_fi and euler_sai = 0, it means that evolution axis is along X-axis
							 //Evolution axis is the axis along particle_rasius and on the other hand it connects particle center to 1st Lagrangian point
		particle[np].center[0].Qu0 = cos(0.5*euler_teta[np]) * cos(0.5*(euler_fi[np] + euler_sai[np]));
		particle[np].center[0].Qu1 = sin(0.5*euler_teta[np]) * cos(0.5*(euler_fi[np] - euler_sai[np]));
		particle[np].center[0].Qu2 = sin(0.5*euler_teta[np]) * sin(0.5*(euler_fi[np] - euler_sai[np]));
		particle[np].center[0].Qu3 = cos(0.5*euler_teta[np]) * sin(0.5*(euler_fi[np] + euler_sai[np]));

#if defined MOVING_SPHEROID	|| defined STATIONARY_SPHEROID
		//// In the next 9 lines we update each Lagragian node position if the initial values of Euler angles are not zero. (It means particle rotation)
		quaternion(particle[np].center[0].Qu0, particle[np].center[0].Qu1, particle[np].center[0].Qu2, particle[np].center[0].Qu3, np, particle[np].center[0].quater);
		particle[np].center[0].deter = determinant(np, particle[np].center[0].quater);//determinant of Quaternion matrix
		for (int n = 0; n < particle_num_nodes; ++n) {
			x_pos_1[n] = particle[np].node[n].x - particle[np].center[0].x;
			y_pos_1[n] = particle[np].node[n].y - particle[np].center[0].y;
			z_pos_1[n] = particle[np].node[n].z - particle[np].center[0].z;
			update_rot_spheroid(particle[np].center[0].deter, particle[np].node[n].x, particle[np].node[n].y, particle[np].node[n].z, particle[np].center[0].quater, x_pos_1[n], y_pos_1[n], z_pos_1[n], particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Node location at INERTIAL coordinate is calculated based on the Quaternions
		}
#endif
	}
	/////////////////////////////////////////////


#if defined Recovery_ON
	if (num_particles > 0) {

		for (int npn = 0; npn < num_particles; ++npn) {

			//string npstr = static_cast<ostringstream*>(&(ostringstream() << npn + 1))->str();//convert number to string
            ostringstream oss;
            oss << npn + 1;
            string npstr = oss.str();//convert number to string
			string input_filename("Alborz_Results/data_particle_recovery_");
			input_filename += npstr;
			input_filename += ".dat";
			ifstream input_file;

			/// Open file
			input_file.open(input_filename.c_str());

			/// Read particle data
			input_file >> particle_mass;
			input_file >> parvol;
			input_file >> mom_inertia_x;
			input_file >> mom_inertia_y;
			input_file >> mom_inertia_z;
			input_file >> particle_num_nodes;
			input_file >> particle_area;

			/// Read data
			input_file >> particle[npn].center[0].x; // center position (x-component)
			input_file >> particle[npn].center[0].y; // center position (y-component)
			input_file >> particle[npn].center[0].z; // center position (y-component)
			input_file >> particle[npn].center[0].vel_x; // center velocity (x-component)
			input_file >> particle[npn].center[0].vel_y; // center velocity (y-component)
			input_file >> particle[npn].center[0].vel_z; // center velocity (y-component)
			input_file >> particle[npn].center[0].omgX; // Rotational speed
			input_file >> particle[npn].center[0].omgY; // Rotational speed
			input_file >> particle[npn].center[0].omgZ; // Rotational speed
			input_file >> particle[npn].center[0].tetaX; // Angle
			input_file >> particle[npn].center[0].tetaY; // Angle
			input_file >> particle[npn].center[0].tetaZ; // Angle
			input_file >> particle[npn].center[0].previous_vel_x;
			input_file >> particle[npn].center[0].previous_vel_y;
			input_file >> particle[npn].center[0].previous_vel_z;
			input_file >> particle[npn].center[0].tempx[0];
			input_file >> particle[npn].center[0].tempx[1];
			input_file >> particle[npn].center[0].tempy[0];
			input_file >> particle[npn].center[0].tempy[1];
			input_file >> particle[npn].center[0].tempz[0];
			input_file >> particle[npn].center[0].tempz[1];
			input_file >> particle[npn].center[0].temp_omgX[0];
			input_file >> particle[npn].center[0].temp_omgX[1];
			input_file >> particle[npn].center[0].temp_omgY[0];
			input_file >> particle[npn].center[0].temp_omgY[1];
			input_file >> particle[npn].center[0].temp_omgZ[0];
			input_file >> particle[npn].center[0].temp_omgZ[1];
			input_file >> particle[npn].center[0].Omg_prime_X;
			input_file >> particle[npn].center[0].Omg_prime_Y;
			input_file >> particle[npn].center[0].Omg_prime_Z;
			input_file >> particle[npn].center[0].Qu0;
			input_file >> particle[npn].center[0].Qu1;
			input_file >> particle[npn].center[0].Qu2;
			input_file >> particle[npn].center[0].Qu3;

			input_file >> particle[npn].center[0].temperature; // Angle with respect to Horizontal line
			input_file >> particle[npn].center[0].temp_temperat[0];
			input_file >> particle[npn].center[0].temp_temperat[1];
			input_file >> particle[npn].center[0].temp_part;

			for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
				input_file >> particle[npn].node[np1].x;
				input_file >> particle[npn].node[np1].y;
				input_file >> particle[npn].node[np1].z;
				//input_file >> particle[npn].node[np1].vel_x;
				//input_file >> particle[npn].node[np1].vel_y;
				//input_file >> particle[npn].node[np1].vel_z;
				//input_file >> particle[npn].node[np1].temperature;
#if defined MOVING_SPHEROID	|| defined STATIONARY_SPHEROID
				input_file >> x_pos_1[np1];
				input_file >> y_pos_1[np1];
				input_file >> z_pos_1[np1];
#endif
				input_file >> area[np1];
			}
			/// Close particle recovery data file
			input_file.close();
		}
	}
#endif   // Recovery_ON

	printf("Particle number of surface nodes = %d \n", particle_num_nodes);
	cout << endl;
	printf("STARTING SIMULATION ... \n");
	srand(1);

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

#if defined X_Wall
	start_X = 1;
	end_X = Nx - 2;
#endif

#if defined X_velIn_ZeroGradOut
	start_X = 1;
	end_X = Nx - 2;
#endif
#if defined ZeroGradIn_ZeroGradOut
	start_X = 1;
	end_X = Nx - 2;
#endif
	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

#if defined Y_Wall
	start_Y = 1;
	end_Y = Ny - 2;
#endif

	int start_Z, end_Z;
	start_Z = 0;
	end_Z = Nz - 1;

#if defined Z_Wall
	start_Z = 1;
	end_Z = Nz - 2;
#endif



	/// ---------------------------------------------------------------------------------------------
	/// -----------------------------------  MAIN TIME LOOP -----------------------------------------
	/// ---------------------------------------------------------------------------------------------
	for (int t = 1 + recovery_step; t <= t_num; ++t) { // START OF TIME STEPS (run over all times between 1 and t_num)

													   /// ---------------------------------------------------------------------------------------------
													   /// ----------------------------- LBM AND THERMAL LBM -------------------------------------------
													   /// ---------------------------------------------------------------------------------------------
		LBM(t);

#if defined Flow_With_Thermal_Effect
		LBM_t(t);
#endif
		/// ---------------------------------------------------------------------------------------------
		/// ----------------------------- END OF LBM AND THERMAL LBM ------------------------------------
		/// ---------------------------------------------------------------------------------------------

		/// ---------------------------------------------------------------------------------------------
		/// ----------------------------- EULERIAN FORCE INITIALIZATION ---------------------------------
		/// ---------------------------------------------------------------------------------------------
		/// Important: Forces at Eulerian points must be set to ZERO at each iteration before SPREAD. We
		/// do this at the end of MOMENTA subroutine.
		/// This is done so just to reduce computational cost & time. (We can set all forces Zero here
		/// but it is time consuming)
		/// ---------------------------------------------------------------------------------------------
		/// ---------------------------------------------------------------------------------------------

		/// ---------------------------------------------------------------------------------------------
		/// --------------------------- LAGRANGIAN FORCES CALCULATION AND SPREAD  -----------------------
		/// ---------------------------------------------------------------------------------------------

		/// ---------------------------------------------------------------------------------------------
		/// ----------------------------------------- SDF  ----------------------------------------------
		/// ---------------------------------------------------------------------------------------------
		if ((disp_max_1 + disp_max_2)>(r_buffer - r_cutoff))
		{
			update_verlet(particle);
			disp_max_1 = 0;
			disp_max_2 = 0;
			for (int ii = 0; ii < num_particles; ++ii) {
				particle[ii].center[0].displacement_verlet = 0;
			}
		}

#if defined SDF
		for (int np = 0; np < num_particles; ++np) {

			for (int nn = 0; nn < particle_num_nodes; ++nn) {
				particle[np].node[nn].force_x = 0;
				particle[np].node[nn].force_y = 0;
				particle[np].node[nn].force_z = 0;
				particle[np].node[nn].force_thermal = 0.0;
			}

			interpolate(particle[np]); // interpolate velocity from the Eulerian mesh to Lagrangian node

			heat_surf[np] = 0.0;

			fx_surf[np] = 0.0;
			fy_surf[np] = 0.0;
			fz_surf[np] = 0.0;
			particle[np].center[0].f3x = 0.0; //Total Force due to particle-particle collision
			particle[np].center[0].f3y = 0.0; //Total Force due to particle-particle collision
			particle[np].center[0].f3z = 0.0; //Total Force due to particle-particle collision
			particle[np].center[0].f4x = 0.0; //Total Force due to particle-wall collision
			particle[np].center[0].f4y = 0.0; //Total Force due to particle-wall collision
			particle[np].center[0].f4z = 0.0; //Total Force due to particle-wall collision

			particle[np].center[0].t1_x = 0.0;//Torque due to forces in X direction
			particle[np].center[0].t1_y = 0.0;//Torque due to forces in Y direction
			particle[np].center[0].t1_z = 0.0;//Torque due to forces in Z direction

			particle[np].center[0].t2_x = 0.0;//Torque due to forces in X direction
			particle[np].center[0].t2_y = 0.0;//Torque due to forces in Y direction
			particle[np].center[0].t2_z = 0.0;//Torque due to forces in Z direction

			particle[np].center[0].t3_x = 0.0;
			particle[np].center[0].t3_y = 0.0;
			particle[np].center[0].t3_z = 0.0;

			T_prime_X[np] = 0.0;
			T_prime_Y[np] = 0.0;
			T_prime_Z[np] = 0.0;

#if defined MOVING_SPHERE
			collision_sphere(particle, np);
#endif
#if defined MOVING_SPHEROID
			collision_spheroid(particle, np, x_pos_1, y_pos_1, z_pos_1);
#endif


			for (int n = 0; n < particle_num_nodes; ++n) {
				particle[np].node[n].force_x = (1.0) * (particle[np].center[0].previous_vel_x - particle[np].center[0].omgZ * (particle[np].node[n].y - particle[np].center[0].y) + particle[np].center[0].omgY * (particle[np].node[n].z - particle[np].center[0].z) - particle[np].node[n].vel_x) * area[n];//Surface force
				particle[np].node[n].force_y = (1.0) * (particle[np].center[0].previous_vel_y + particle[np].center[0].omgZ * (particle[np].node[n].x - particle[np].center[0].x) - particle[np].center[0].omgX * (particle[np].node[n].z - particle[np].center[0].z) - particle[np].node[n].vel_y) * area[n];//Surface force
				particle[np].node[n].force_z = (1.0) * (particle[np].center[0].previous_vel_z + particle[np].center[0].omgX * (particle[np].node[n].y - particle[np].center[0].y) - particle[np].center[0].omgY * (particle[np].node[n].x - particle[np].center[0].x) - particle[np].node[n].vel_z) * area[n];//Surface force

				particle[np].node[n].force_thermal = 1 * (particle[np].center[0].temp_part - particle[np].node[n].temperature) * area[n];//Surface heat force

				particle[np].center[0].t1_x += particle[np].node[n].force_z * (particle[np].node[n].y - particle[np].center[0].y);
				particle[np].center[0].t2_x += particle[np].node[n].force_y * (particle[np].node[n].z - particle[np].center[0].z);
				particle[np].center[0].t1_y += particle[np].node[n].force_z * (particle[np].node[n].x - particle[np].center[0].x);
				particle[np].center[0].t2_y += particle[np].node[n].force_x * (particle[np].node[n].z - particle[np].center[0].z);
				particle[np].center[0].t1_z += particle[np].node[n].force_x * (particle[np].node[n].y - particle[np].center[0].y);
				particle[np].center[0].t2_z += particle[np].node[n].force_y * (particle[np].node[n].x - particle[np].center[0].x);
				fx_surf[np] += particle[np].node[n].force_x;
				fy_surf[np] += particle[np].node[n].force_y;
				fz_surf[np] += particle[np].node[n].force_z;
				heat_surf[np] += particle[np].node[n].force_thermal;
			}
			spread(particle[np]); // spread forces from the Lagrangian to the Eulerian mesh
		}

#endif //End of:SDF
		/// ---------------------------------------------------------------------------------------------
		/// ----------------------------------------- END OF SDF  ---------------------------------------
		/// ---------------------------------------------------------------------------------------------m


		/// ---------------------------------------------------------------------------------------------
		/// --------------------------------- GUO MODEL  ------------------------------------------------
		/// ---------------------------------------------------------------------------------------------
#if defined Guo
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				for (int Z = start_Z; Z <= end_Z; ++Z) {
					velocity_x[X][Y][Z] += (force_x[X][Y][Z]) / 2.0;
					velocity_y[X][Y][Z] += (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) / 2.0;
					velocity_z[X][Y][Z] += (force_z[X][Y][Z]) / 2.0;
				}
			}
		}
#endif //End of: Guo
		/// ---------------------------------------------------------------------------------------------
		/// --------------------------------- END OF GUO MODEL  -----------------------------------------
		/// ---------------------------------------------------------------------------------------------

#if defined Force_Thermal_2
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				for (int Z = start_Z; Z <= end_Z; ++Z) {
					temperature[X][Y][Z] = temperature[X][Y][Z] + force_thermal[X][Y][Z] / 2.0;
				}
			}
		}
#endif

		/// ---------------------------------------------------------------------------------------------
		/// ------------------------- NEWTON'S EQUATIONS OF MOTION  -------------------------------------
		/// ---------------------------------------------------------------------------------------------
		for (int np = 0; np < num_particles; ++np) {
			check_zero(particle[np].center[0].t1_x, 1.0e-6);
			check_zero(particle[np].center[0].t1_y, 1.0e-6);
			check_zero(particle[np].center[0].t1_z, 1.0e-6);

			check_zero(particle[np].center[0].t2_x, 1.0e-6);
			check_zero(particle[np].center[0].t2_y, 1.0e-6);
			check_zero(particle[np].center[0].t2_z, 1.0e-6);


			////////// Updating Velocity and center position ///////
			if (t > margin) {
				if (t > margin + 1) {
#if defined With_added_mass
					particle[np].center[0].vel_x = particle[np].center[0].vel_x - fx_surf[np] / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3x / particle_mass + particle[np].center[0].f4x / particle_mass + 1.0 / den_ratio*(particle[np].center[0].vel_x - particle[np].center[0].tempx[0]);//Newton law
					particle[np].center[0].vel_y = particle[np].center[0].vel_y - fy_surf[np] / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3y / particle_mass + particle[np].center[0].f4y / particle_mass + 1.0 / den_ratio*(particle[np].center[0].vel_y - particle[np].center[0].tempy[0]);//Newton law
					particle[np].center[0].vel_z = particle[np].center[0].vel_z - fz_surf[np] / particle_mass + particle_gravity_z * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3z / particle_mass + particle[np].center[0].f4z / particle_mass + 1.0 / den_ratio*(particle[np].center[0].vel_z - particle[np].center[0].tempz[0]);//Newton law
#endif
#if defined Without_added_mass
					particle[np].center[0].vel_x = particle[np].center[0].vel_x - fx_surf[np] / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3x / particle_mass + particle[np].center[0].f4x / particle_mass;//Newton law
					particle[np].center[0].vel_y = particle[np].center[0].vel_y - fy_surf[np] / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3y / particle_mass + particle[np].center[0].f4y / particle_mass;//Newton law
					particle[np].center[0].vel_z = particle[np].center[0].vel_z - fz_surf[np] / particle_mass + particle_gravity_z * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3z / particle_mass + particle[np].center[0].f4z / particle_mass;//Newton law
#endif
				}
				else {
					particle[np].center[0].vel_x = particle[np].center[0].vel_x - fx_surf[np] / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3x / particle_mass + particle[np].center[0].f4x / particle_mass;//Newton law
					particle[np].center[0].vel_y = particle[np].center[0].vel_y - fy_surf[np] / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3y / particle_mass + particle[np].center[0].f4y / particle_mass;//Newton law
					particle[np].center[0].vel_z = particle[np].center[0].vel_z - fz_surf[np] / particle_mass + particle_gravity_z * (1.0 - 1.0 / den_ratio) + particle[np].center[0].f3z / particle_mass + particle[np].center[0].f4z / particle_mass;//Newton law
				}

				particle[np].center[0].tempx[0] = particle[np].center[0].tempx[1];
				particle[np].center[0].tempx[1] = particle[np].center[0].vel_x;
				particle[np].center[0].tempy[0] = particle[np].center[0].tempy[1];
				particle[np].center[0].tempy[1] = particle[np].center[0].vel_y;
				particle[np].center[0].tempz[0] = particle[np].center[0].tempz[1];
				particle[np].center[0].tempz[1] = particle[np].center[0].vel_z;

				center_old_x = particle[np].center[0].x;
				center_old_y = particle[np].center[0].y;
				center_old_z = particle[np].center[0].z;

				disp_x = 0.5 *(particle[np].center[0].vel_x + particle[np].center[0].tempx[0]);
				disp_y = 0.5 *(particle[np].center[0].vel_y + particle[np].center[0].tempy[0]);
				disp_z = 0.5 *(particle[np].center[0].vel_z + particle[np].center[0].tempz[0]);
				particle[np].center[0].x = particle[np].center[0].x + disp_x;
				particle[np].center[0].y = particle[np].center[0].y + disp_y;
				particle[np].center[0].z = particle[np].center[0].z + disp_z;

				particle[np].center[0].displacement_verlet += sqrt(SQ(disp_x) + SQ(disp_y) + SQ(disp_z));//Displacement between two Verlet list update steps
				if (particle[np].center[0].displacement_verlet>disp_max_1)
				{
					disp_max_2 = disp_max_1;
					disp_max_1 = particle[np].center[0].displacement_verlet;
				}
				else if ((particle[np].center[0].displacement_verlet>disp_max_2) && (particle[np].center[0].displacement_verlet<disp_max_1))
				{
					disp_max_2 = particle[np].center[0].displacement_verlet;
				}
			}
			////////////////////////////////////////////

			////////// Updating Angular Velocity ///////
#if defined MOVING_SPHERE  //Calculating rotational speed for sphere
			/////////// For X rotation
			if (t > margin) {
				if (t > margin + 1) {
#if defined With_added_mass
					particle[np].center[0].omgX = particle[np].center[0].omgX + 1.0 / mom_inertia_x * (particle[np].center[0].t2_x - particle[np].center[0].t1_x) + 1.0 / den_ratio*(particle[np].center[0].omgX - particle[np].center[0].temp_omgX[0]);//Eq. 100 Kang Thesis (804)
#endif
#if defined Without_added_mass
					particle[np].center[0].omgX = particle[np].center[0].omgX + 1.0 / mom_inertia_x * (particle[np].center[0].t2_x - particle[np].center[0].t1_x);//Eq. 100 Kang Thesis (804)
#endif
				}
				else {
					particle[np].center[0].omgX = particle[np].center[0].omgX + 1.0 / mom_inertia_x * (particle[np].center[0].t2_x - particle[np].center[0].t1_x);
				}

				check_zero(particle[np].center[0].omgX, 1.0e-12);

				particle[np].center[0].temp_omgX[0] = particle[np].center[0].temp_omgX[1];
				particle[np].center[0].temp_omgX[1] = particle[np].center[0].omgX;
				/////////// For Y rotation
				if (t > margin + 1) {
#if defined With_added_mass
					particle[np].center[0].omgY = particle[np].center[0].omgY + 1.0 / mom_inertia_y * (particle[np].center[0].t1_y - particle[np].center[0].t2_y) + 1.0 / den_ratio*(particle[np].center[0].omgY - particle[np].center[0].temp_omgY[0]);//Eq. 100 Kang Thesis (804)
#endif
#if defined Without_added_mass
					particle[np].center[0].omgY = particle[np].center[0].omgY + 1.0 / mom_inertia_y * (particle[np].center[0].t1_y - particle[np].center[0].t2_y);//Eq. 100 Kang Thesis (804)
#endif
				}
				else {
					particle[np].center[0].omgY = particle[np].center[0].omgY + 1.0 / mom_inertia_y * (particle[np].center[0].t1_y - particle[np].center[0].t2_y);
				}

				check_zero(particle[np].center[0].omgY, 1.0e-12);
				particle[np].center[0].temp_omgY[0] = particle[np].center[0].temp_omgY[1];
				particle[np].center[0].temp_omgY[1] = particle[np].center[0].omgY;
				/////////// For Z rotation
				if (t > margin + 1) {
#if defined With_added_mass
					particle[np].center[0].omgZ = particle[np].center[0].omgZ + 1.0 / mom_inertia_z * (particle[np].center[0].t1_z - particle[np].center[0].t2_z) + 1.0 / den_ratio*(particle[np].center[0].omgZ - particle[np].center[0].temp_omgZ[0]);//Eq. 100 Kang Thesis (804)
#endif
#if defined Without_added_mass
					particle[np].center[0].omgZ = particle[np].center[0].omgZ + 1.0 / mom_inertia_z * (particle[np].center[0].t1_z - particle[np].center[0].t2_z);//Eq. 100 Kang Thesis (804)
#endif
				}
				else {
					particle[np].center[0].omgZ = particle[np].center[0].omgZ + 1.0 / mom_inertia_z * (particle[np].center[0].t1_z - particle[np].center[0].t2_z);
				}

				check_zero(particle[np].center[0].omgZ, 1.0e-12);
				particle[np].center[0].temp_omgZ[0] = particle[np].center[0].temp_omgZ[1];
				particle[np].center[0].temp_omgZ[1] = particle[np].center[0].omgZ;
				//////////////////


				particle[np].center[0].delta_teta_X = 0.5 * (particle[np].center[0].omgX + particle[np].center[0].temp_omgX[0]);
				particle[np].center[0].delta_teta_Y = 0.5 * (particle[np].center[0].omgY + particle[np].center[0].temp_omgY[0]);
				particle[np].center[0].delta_teta_Z = 0.5 * (particle[np].center[0].omgZ + particle[np].center[0].temp_omgZ[0]);
				particle[np].center[0].tetaX = particle[np].center[0].tetaX + 0.5 * (particle[np].center[0].omgX + particle[np].center[0].temp_omgX[0]);
				particle[np].center[0].tetaY = particle[np].center[0].tetaY + 0.5 * (particle[np].center[0].omgY + particle[np].center[0].temp_omgY[0]);
				particle[np].center[0].tetaZ = particle[np].center[0].tetaZ + 0.5 * (particle[np].center[0].omgZ + particle[np].center[0].temp_omgZ[0]);

				// Applying the rotation matrix
				for (int n = 0; n < particle_num_nodes; ++n) {
					//Around Z-Axis rotation
					x_pos = particle[np].node[n].x - center_old_x;
					y_pos = particle[np].node[n].y - center_old_y;
					z_pos = particle[np].node[n].z - center_old_z;
					particle[np].node[n].x = particle[np].center[0].x + x_pos * cos(particle[np].center[0].delta_teta_Z) - y_pos * sin(particle[np].center[0].delta_teta_Z);
					particle[np].node[n].y = particle[np].center[0].y + x_pos * sin(particle[np].center[0].delta_teta_Z) + y_pos * cos(particle[np].center[0].delta_teta_Z);
					particle[np].node[n].z = particle[np].center[0].z + z_pos;
					//Around X-Axis rotation
					x_pos = particle[np].node[n].x - particle[np].center[0].x;
					y_pos = particle[np].node[n].y - particle[np].center[0].y;
					z_pos = particle[np].node[n].z - particle[np].center[0].z;
					particle[np].node[n].x = particle[np].center[0].x + x_pos;
					particle[np].node[n].y = particle[np].center[0].y + y_pos * cos(particle[np].center[0].delta_teta_X) - z_pos * sin(particle[np].center[0].delta_teta_X);
					particle[np].node[n].z = particle[np].center[0].z + y_pos * sin(particle[np].center[0].delta_teta_X) + z_pos * cos(particle[np].center[0].delta_teta_X);

					//Around Y-Axis rotation
					x_pos = particle[np].node[n].x - particle[np].center[0].x;
					y_pos = particle[np].node[n].y - particle[np].center[0].y;
					z_pos = particle[np].node[n].z - particle[np].center[0].z;
					particle[np].node[n].x = particle[np].center[0].x + x_pos * cos(particle[np].center[0].delta_teta_Y) + z_pos * sin(particle[np].center[0].delta_teta_Y);;
					particle[np].node[n].y = particle[np].center[0].y + y_pos;
					particle[np].node[n].z = particle[np].center[0].z - x_pos * sin(particle[np].center[0].delta_teta_Y) + z_pos * cos(particle[np].center[0].delta_teta_Y);

					/////////////////////////
					particle[np].node[n].vel_x = particle[np].center[0].vel_x;  // I think these 3 lines can be removed because node velocity is determined by INTERPOLATION
					particle[np].node[n].vel_y = particle[np].center[0].vel_y;
					particle[np].node[n].vel_z = particle[np].center[0].vel_z;
				}
			}

#endif  //End of rotation procedure for SPHERE


#if defined MOVING_SPHEROID   //Calculating rotation speed for spheroid. Using quaternion concept. See: H. Huang et al, J. Fluid Mech. (2012), vol. 692, pp. 369-394.
			if (t > margin) {
				quaternion(particle[np].center[0].Qu0, particle[np].center[0].Qu1, particle[np].center[0].Qu2, particle[np].center[0].Qu3, np, particle[np].center[0].quater);
				double torque_1 = (particle[np].center[0].t2_x - particle[np].center[0].t1_x - particle[np].center[0].t3_x);
				double torque_2 = (particle[np].center[0].t1_y - particle[np].center[0].t2_y - particle[np].center[0].t3_y);
				double torque_3 = (particle[np].center[0].t1_z - particle[np].center[0].t2_z - particle[np].center[0].t3_z);
				// Torque transformation from Inertial coordinate to body-fixed coordinate
				T_prime_X[np] = particle[np].center[0].quater[1] * torque_1 + particle[np].center[0].quater[2] * torque_2 + particle[np].center[0].quater[3] * torque_3;
				T_prime_Y[np] = particle[np].center[0].quater[4] * torque_1 + particle[np].center[0].quater[5] * torque_2 + particle[np].center[0].quater[6] * torque_3;
				T_prime_Z[np] = particle[np].center[0].quater[7] * torque_1 + particle[np].center[0].quater[8] * torque_2 + particle[np].center[0].quater[9] * torque_3;

				rungekutta(particle[np].center[0].Omg_prime_X, particle[np].center[0].Omg_prime_Y, particle[np].center[0].Omg_prime_Z, T_prime_X[np], T_prime_Y[np], T_prime_Z[np], mom_inertia_x, mom_inertia_y, mom_inertia_z, particle[np].center[0].Qu0, particle[np].center[0].Qu1, particle[np].center[0].Qu2, particle[np].center[0].Qu3);

				quaternion(particle[np].center[0].Qu0, particle[np].center[0].Qu1, particle[np].center[0].Qu2, particle[np].center[0].Qu3, np, particle[np].center[0].quater);

				sumq = sqrt(SQ(particle[np].center[0].Qu0) + SQ(particle[np].center[0].Qu1) + SQ(particle[np].center[0].Qu2) + SQ(particle[np].center[0].Qu3));
				particle[np].center[0].Qu0 = particle[np].center[0].Qu0 / sumq;
				particle[np].center[0].Qu1 = particle[np].center[0].Qu1 / sumq;
				particle[np].center[0].Qu2 = particle[np].center[0].Qu2 / sumq;
				particle[np].center[0].Qu3 = particle[np].center[0].Qu3 / sumq;
				particle[np].center[0].deter = determinant(np, particle[np].center[0].quater);//determinant of Quaternion matrix

				comp_rot_speed(particle[np].center[0].deter, particle[np].center[0].omgX, particle[np].center[0].omgY, particle[np].center[0].omgZ, np, particle[np].center[0].quater, particle[np].center[0].Omg_prime_X, particle[np].center[0].Omg_prime_Y, particle[np].center[0].Omg_prime_Z);// Angular velocity transformation from body-fixed coordinate to Inertial coordinate.

				if (t > margin + 1) {
#if defined With_added_mass
					particle[np].center[0].omgX = particle[np].center[0].omgX + 1.0 / den_ratio*(particle[np].center[0].omgX - particle[np].center[0].temp_omgX[0]);//Eq. 100 Kang Thesis (804)
					particle[np].center[0].omgY = particle[np].center[0].omgY + 1.0 / den_ratio*(particle[np].center[0].omgY - particle[np].center[0].temp_omgY[0]);//Eq. 100 Kang Thesis (804)
					particle[np].center[0].omgZ = particle[np].center[0].omgZ + 1.0 / den_ratio*(particle[np].center[0].omgZ - particle[np].center[0].temp_omgZ[0]);//Eq. 100 Kang Thesis (804)
					particle[np].center[0].Omg_prime_X = particle[np].center[0].quater[1] * (particle[np].center[0].omgX) + particle[np].center[0].quater[2] * (particle[np].center[0].omgY) + particle[np].center[0].quater[3] * (particle[np].center[0].omgZ);
					particle[np].center[0].Omg_prime_Y = particle[np].center[0].quater[4] * (particle[np].center[0].omgX) + particle[np].center[0].quater[5] * (particle[np].center[0].omgY) + particle[np].center[0].quater[6] * (particle[np].center[0].omgZ);
					particle[np].center[0].Omg_prime_Z = particle[np].center[0].quater[7] * (particle[np].center[0].omgX) + particle[np].center[0].quater[8] * (particle[np].center[0].omgY) + particle[np].center[0].quater[9] * (particle[np].center[0].omgZ);
#endif
				}

				particle[np].center[0].temp_omgX[0] = particle[np].center[0].temp_omgX[1];
				particle[np].center[0].temp_omgX[1] = particle[np].center[0].omgX;
				particle[np].center[0].temp_omgY[0] = particle[np].center[0].temp_omgY[1];
				particle[np].center[0].temp_omgY[1] = particle[np].center[0].omgY;
				particle[np].center[0].temp_omgZ[0] = particle[np].center[0].temp_omgZ[1];
				particle[np].center[0].temp_omgZ[1] = particle[np].center[0].omgZ;

				particle[np].center[0].delta_teta_X = 0.5 * (particle[np].center[0].omgX + particle[np].center[0].temp_omgX[0]);// For spheroid they are not important because node location is updated by quaternions
				particle[np].center[0].delta_teta_Y = 0.5 * (particle[np].center[0].omgY + particle[np].center[0].temp_omgY[0]);
				particle[np].center[0].delta_teta_Z = 0.5 * (particle[np].center[0].omgZ + particle[np].center[0].temp_omgZ[0]);
				particle[np].center[0].tetaX = particle[np].center[0].tetaX + 0.5 * (particle[np].center[0].omgX + particle[np].center[0].temp_omgX[0]);
				particle[np].center[0].tetaY = particle[np].center[0].tetaY + 0.5 * (particle[np].center[0].omgY + particle[np].center[0].temp_omgY[0]);
				particle[np].center[0].tetaZ = particle[np].center[0].tetaZ + 0.5 * (particle[np].center[0].omgZ + particle[np].center[0].temp_omgZ[0]);


				// Applying the rotation matrix
				for (int n = 0; n < particle_num_nodes; ++n) {
					//Updating nodes location
					update_rot_spheroid(particle[np].center[0].deter, particle[np].node[n].x, particle[np].node[n].y, particle[np].node[n].z, particle[np].center[0].quater, x_pos_1[n], y_pos_1[n], z_pos_1[n], particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Node location at INERTIAL coordinate is calculated based on the Quaternions


																																																																							   /////////////////////////
					particle[np].node[n].vel_x = particle[np].center[0].vel_x;  // I think these 3 lines can be removed because node velocity is determined by INTERPOLATION
					particle[np].node[n].vel_y = particle[np].center[0].vel_y;
					particle[np].node[n].vel_z = particle[np].center[0].vel_z;
				}
			}
#endif //End of rotation for spheroid
			////////////////////////////////////////////

			////////// Updating Temperature ///////
#if defined Particle_Varying_Temperature
			if (t > 1) {
				particle[np].center[0].temperature = particle[np].center[0].temperature - heat_surf[np] / (particle_mass*specific_heat_ratio) + 1.0 / (den_ratio*specific_heat_ratio)*(particle[np].center[0].temperature - particle[np].center[0].temp_temperat[0]) + particle_heat_source / (particle_mass*specific_heat_ratio);
			}
			else {
				particle[np].center[0].temperature = particle[np].center[0].temperature - heat_surf[np] / (particle_mass*specific_heat_ratio) + particle_heat_source / (particle_mass*specific_heat_ratio);
			}
#endif
			if (t == 1) {
				particle[np].center[0].temp_temperat[0] = particle_temperature;
				particle[np].center[0].temp_temperat[1] = particle_temperature;
			}
			particle[np].center[0].temp_temperat[0] = particle[np].center[0].temp_temperat[1];
			particle[np].center[0].temp_temperat[1] = particle[np].center[0].temperature;
			/////////////////////////////////////////

			if (t > margin) {
#if defined X_Periodic  //In periodic condition, if particle exits from one side it should enter from other side
				if (particle[np].center[0].x < 0) {
					disp_max_1 = 100000;// A large number
					disp_max_2 = 100000;
					particle[np].center[0].x += Nx;
					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						particle[np].node[np1].x += Nx;
					}
				}
				else if (particle[np].center[0].x >= Nx) {
					disp_max_1 = 100000;// A large number
					disp_max_2 = 100000;
					particle[np].center[0].x -= Nx;
					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						particle[np].node[np1].x -= Nx;
					}
				}
#endif

#if defined Z_Periodic //In periodic condition, if particle exits from one side it should enter from other side
				if (particle[np].center[0].z < 0) {
					disp_max_1 = 100000;// A large number
					disp_max_2 = 100000;
					particle[np].center[0].z += Nz;
					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						particle[np].node[np1].z += Nz;
					}
				}
				else if (particle[np].center[0].z > Nz) {
					disp_max_1 = 100000;// A large number
					disp_max_2 = 100000;
					particle[np].center[0].z -= Nz;
					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						particle[np].node[np1].z -= Nz;
					}
				}
#endif
#if defined Y_Periodic //In periodic condition, if particle exits from one side it should enter from other side
				if (particle[np].center[0].y < 0) {
					disp_max_1 = 100000;// A large number
					disp_max_2 = 100000;
					particle[np].center[0].y += Ny;
					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						particle[np].node[np1].y += Ny;
					}
				}
				else if (particle[np].center[0].y >= Ny) {
					disp_max_1 = 100000;// A large number
					disp_max_2 = 100000;
					particle[np].center[0].y -= Ny;
					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						particle[np].node[np1].y -= Ny;
					}
				}
#endif
			} //End of : if (t > margin)

			particle[np].center[0].previous_vel_x = particle[np].center[0].vel_x;
			particle[np].center[0].previous_vel_y = particle[np].center[0].vel_y;
			particle[np].center[0].previous_vel_z = particle[np].center[0].vel_z;

#if defined Particle_Varying_Temperature
			particle[np].center[0].temp_part = particle[np].center[0].temperature;
#endif

		} //End of: NEWTON'S EQUATIONS OF MOTION
		  /// ---------------------------------------------------------------------------------------------
		  /// ------------------------- END OF NEWTON'S EQUATIONS OF MOTION  ------------------------------
		  /// ---------------------------------------------------------------------------------------------

		  /// ****************************************
		  ///  WRITE FLUID AND PARTICLES DATA TO DISC
		  /// ****************************************
		if (t % t_data == 0) {
			for (int npn = 0; npn < num_particles; ++npn) {
				write_particle_data(t, particle[npn], npn, fx_surf[npn], fy_surf[npn], fz_surf[npn]);
			}
			write_fluid_data(t);//You can activate or deactivate this line.
		}


		/// ***************************
		///  WRITE RECOVERY DATA
		/// ***************************
		if (t % t_recovery == 0) {
			write_fluid_recovery(t); //Recommended: You can write fluid data into a "dat" file. In case of sudden program failure just activate "Recovery_ON" and run it so the program continues from a almost last iteration

			if (num_particles > 0) {

				for (int npn = 0; npn < num_particles; ++npn) {
				    std::string npstr = numberToString(npn + 1);
					//string npstr = static_cast<ostringstream*>(&(ostringstream() << npn + 1))->str();//convert number to string
					string output_filename("Alborz_Results/data_particle_recovery_");
					output_filename += npstr;
					output_filename += ".dat";
					ofstream output_file;

					/// Open file
					output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

																			   /// Write particle data
					output_file << std::fixed << std::setprecision(12) << particle_mass << " ";
					output_file << std::fixed << std::setprecision(12) << parvol << " ";
					output_file << std::fixed << std::setprecision(12) << mom_inertia_x << " ";
					output_file << std::fixed << std::setprecision(12) << mom_inertia_y << " ";
					output_file << std::fixed << std::setprecision(12) << mom_inertia_z << " ";
					output_file << std::fixed << std::setprecision(12) << particle_num_nodes << " ";
					output_file << std::fixed << std::setprecision(12) << particle_area << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].x << " "; // center position (x-component)
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].y << " "; // center position (y-component)
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].z << " "; // center position (y-component)
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].vel_x << " "; // center velocity (x-component)
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].vel_y << " "; // center velocity (y-component)
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].vel_z << " "; // center velocity (y-component)
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].omgX << " "; // Rotational speed
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].omgY << " "; // Rotational speed
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].omgZ << " "; // Rotational speed
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tetaX << " "; // Angle with respect to Horizontal line
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tetaY << " "; // Angle with respect to Horizontal line
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tetaZ << " "; // Angle with respect to Horizontal line
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].previous_vel_x << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].previous_vel_y << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].previous_vel_z << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tempx[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tempx[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tempy[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tempy[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tempz[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].tempz[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_omgX[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_omgX[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_omgY[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_omgY[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_omgZ[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_omgZ[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Omg_prime_X << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Omg_prime_Y << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Omg_prime_Z << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Qu0 << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Qu1 << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Qu2 << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].Qu3 << " ";

					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temperature << " "; // Temperature
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_temperat[0] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_temperat[1] << " ";
					output_file << std::fixed << std::setprecision(12) << particle[npn].center[0].temp_part << " ";

					for (int np1 = 0; np1 < particle_num_nodes; ++np1) {
						output_file << std::fixed << std::setprecision(12) << particle[npn].node[np1].x << " ";
						output_file << std::fixed << std::setprecision(12) << particle[npn].node[np1].y << " ";
						output_file << std::fixed << std::setprecision(12) << particle[npn].node[np1].z << " ";
						//output_file << std::fixed << std::setprecision(12) << particle[npn].node[np1].vel_x << " ";
						//output_file << std::fixed << std::setprecision(12) << particle[npn].node[np1].vel_y << " ";
						//output_file << std::fixed << std::setprecision(12) << particle[npn].node[np1].vel_z << " ";
						//output_file << std::fixed << std::setprecision(15) << particle[npn].node[np1].temperature << " ";
#if defined MOVING_SPHEROID	|| defined STATIONARY_SPHEROID
						output_file << std::fixed << std::setprecision(12) << x_pos_1[np1] << " ";
						output_file << std::fixed << std::setprecision(12) << y_pos_1[np1] << " ";
						output_file << std::fixed << std::setprecision(12) << z_pos_1[np1] << " ";
#endif
						output_file << std::fixed << std::setprecision(12) << area[np1] << " ";
					}
					/// Close particle recovery data file
					output_file.close();
				}

			}
		}


		/// ***************************
		///  WRITE VTK DATA
		/// ***************************
		if (t % t_vtk == 0) {
			write_fluid_vtk(t);
			for (int npn = 0; npn < num_particles; ++npn) {
				write_particle_vtk(t, particle[npn], npn);
			}
		}
		/// ***************************
		///  REPORT TIME ON SCREEN
		/// ***************************
		if (t % t_time == 0) {
			elapsed_time1 = clock() - elapsed_time;
			printf("The program took %f seconds from beginning .\n", ((float)elapsed_time1) / CLOCKS_PER_SEC);
		}
		/// ***************************
		///  REPORT ON COMMAND SCREEN
		/// ***************************
		/// Report end of time step

		if (t % t_info == 0) {
			cout << "Completed Time Step " << t << " in [1, " << t_num << "]" << endl;
			printf("===========================================================================\n");
		}

	}//End of time step loop

	 /// Report successful end of simulation

	printf("SIMULATION COMPLETED. THANKS \n");
	delete[] force_x;
	delete[] force_y;
	delete[] force_z;
	delete[] pop;
	delete[] density;
	delete[] velocity_x;
	delete[] velocity_y;
	delete[] velocity_z;
	delete[] pop_old;
	delete[] x_pos_1;
	delete[] y_pos_1;
	delete[] z_pos_1;
	delete[] area;
	delete[] mdf_force_x;
	delete[] mdf_force_y;
	delete[] mdf_force_z;
	delete[] mdf_vel_x;
	delete[] mdf_vel_y;
	delete[] mdf_vel_z;
	delete[] force_thermal;
	delete[] temperature;
	elapsed_time = clock() - elapsed_time;
	printf("THE PROGRAM TOOK %d CLICKS (%f seconds).\n", elapsed_time, ((double)elapsed_time) / CLOCKS_PER_SEC);
//	getch();
	std::cin.clear(); // reset any error flags
	std::cin.ignore(32767, '\n'); // ignore any characters in the input buffer until we find an enter character
	std::cin.get(); // get one more char from the user
#endif // END of Direct Forcing for Moving Spheres & Spheroids


#if defined Flow_Without_Particle
	cout << endl;
	printf("STARTING SIMULATION ... \n");

#if defined Read_Geometry_From_File
	// Read geometry from the file (geometry.txt). This file must be created by MATLAB program: image_to_matrix.m
	// Before using program specify Nx and Ny at the beginning of program according to geometry size.
	int sol;
	string output_filename("geometry");
	output_filename += ".txt";
	ifstream output_file;
	/// Open file
	output_file.open(output_filename.c_str());
	/// Read geometry data
	for (int X = 0; X<Nx; X++) {
		for (int Y = 0; Y<Ny; Y++) {
			for (int Z = 0; Z<Nz; Z++) {
				if (!output_file.eof()) {
					output_file >> sol;
					if (sol == 0) { // 0 in geometry.txt means SOLID
						is_solid[X][Y][Z] = TRUE;
					}
				}
			}
		}
	}

	/// Close geometry file
	output_file.close();

#endif

	//// In below lines specify which nodes you want to be solid. Hint: You can define the whole domain by periodic boundary condition and then define the walls by is_solid
	//// This is_solid approach can be just used with periodic or wall boundary condition.
	//// You can use below lines in addition to read_geometry_from_file or deactive below lines.
	//for (int X = 20; X <= 40; ++X) {//Solid nodes
	//	for (int Y = 20; Y <= 40; ++Y) {
	//		for (int Z = 20; Z <= 40; ++Z) {
	//			is_solid[X][Y][Z] = TRUE;
	//		}
	//	}
	//}

	// Counting total solid nodes.
	// Do not change below lines:
	int total_solid = 0;
	for (int X = 0; X < Nx; ++X) {//Solid nodes
		for (int Y = 0; Y < Ny; ++Y) {
			for (int Z = 0; Z < Nz; ++Z) {
				if (is_solid[X][Y][Z] == TRUE) {
					total_solid = total_solid + 1;
				}
			}
		}
	}

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

#if defined X_Wall
	start_X = 1;
	end_X = Nx - 2;
#endif

#if defined X_velIn_ZeroGradOut
	start_X = 1;
	end_X = Nx - 2;
#endif
#if defined ZeroGradIn_ZeroGradOut
	start_X = 1;
	end_X = Nx - 2;
#endif
	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

#if defined Y_Wall
	start_Y = 1;
	end_Y = Ny - 2;
#endif

	int start_Z, end_Z;
	start_Z = 0;
	end_Z = Nz - 1;

#if defined Z_Wall
	start_Z = 1;
	end_Z = Nz - 2;
#endif

	for (int t = 1 + recovery_step; t <= t_num; ++t) { // START OF TIME STEPS (run over all times between 1 and t_num)

													   /// ******************************
													   ///  LBM ALGORITHM
													   /// ******************************
		LBM(t);

		/// ******************************
		///  LBM ALGORITHM FOR TEMPERATURE
		/// ******************************
#if defined Flow_With_Thermal_Effect
		LBM_t(t);
#endif

		/// ****************
		///     GUO MODEL
		/// ****************
#if defined Guo
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				for (int Z = start_Z; Z <= end_Z; ++Z) {
					if (is_solid[X][Y][Z] == FALSE) {
						velocity_x[X][Y][Z] += gravity / 2.0;
						velocity_y[X][Y][Z] += gbeta*(temperature[X][Y][Z] - Tref) / 2.0; // Because gravity just acts in X direction
						velocity_z[X][Y][Z] += 0; // Because gravity just acts in X direction
					}
				}
			}
		}
#endif

#if defined Force_Thermal_2
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				for (int Z = start_Z; Z <= end_Z; ++Z) {
					temperature[X][Y][Z] = temperature[X][Y][Z] + force_thermal[X][Y][Z] / 2.0;
				}
			}
		}
#endif
		/// ****************************************
		///  WRITE FLUID AND PARTICLES DATA TO DISC
		/// ****************************************
		if (t % t_data == 0) {
			write_fluid_data(t);//You can write fluid data into a .dat file by activating this command
		}


		/// ***************************
		///  WRITE RECOVERY DATA
		/// ***************************
		if (t % t_recovery == 0) {
			write_fluid_recovery(t);
		}
		/// ***************************
		///  WRITE VTK DATA
		/// ***************************
		if (t % t_vtk == 0) {
			write_fluid_vtk(t);
		}
		if (t == t_vtk) {
			write_solid_vtk(t, total_solid);
		}
		/// ***************************
		///  REPORT TIME ON SCREEN
		/// ***************************
		if (t % t_time == 0) {
			elapsed_time1 = clock() - elapsed_time;
			printf("The program took %f seconds from beginning .\n", ((float)elapsed_time1) / CLOCKS_PER_SEC);
		}

		/// ***************************
		///  REPORT ON COMMAND SCREEN
		/// ***************************
		/// Report end of time step
		if (t % t_info == 0) {
			cout << "Completed Time Step " << t << " in [1, " << t_num << "]" << endl;
			printf("===========================================================================\n");
		}

	}//END OF TIME STEP LOOP

	 /// Report successful end of simulation
	printf("SIMULATION COMPLETED. THANKS \n");
	delete[] force_x;
	delete[] force_y;
	delete[] force_z;
	delete[] pop;
	delete[] density;
	delete[] velocity_x;
	delete[] velocity_y;
	delete[] velocity_z;
	delete[] pop_old;
	delete[] force_thermal;
	delete[] temperature;
	elapsed_time = clock() - elapsed_time;
	printf("THE PROGRAM TOOK %d CLICKS (%f seconds).\n", elapsed_time, ((double)elapsed_time) / CLOCKS_PER_SEC);
	getch();
	std::cin.clear(); // reset any error flags
	std::cin.ignore(32767, '\n'); // ignore any characters in the input buffer until we find an enter character
	std::cin.get(); // get one more char from the user
#endif // End of: Direct forcing and Flow without particle

	int num_r;
	//scanf_s("%i", &num_r);
	return 0;
} // end of main function

  //=============================================================================================================================
  //=============================================================================================================================
  //-----------------------------------------------------------------------------------------------------------------------------
  // Subroutine: INITIALIZE
  //-----------------------------------------------------------------------------------------------------------------------------
  //
  // The memory for lattice variables (populations, density, velocity, forces) is allocated.
  // The variables are initialized.

void initialize(int & rec_step) {

	/// Create folders, delete data file
	// Make sure that the VTK folders exist.
	// Old file data.dat is deleted, if existing.
	//int ignore; // ignore return value of system calls
	//ignore = system("mkdir Alborz_Results\\vtk_fluid"); // create folder if not existing
	//ignore = system("mkdir Alborz_Results\\vtk_particle"); // create folder if not existing
	//ignore = system("mkdir Alborz_Results\\Data_particle"); // create folder if not existing
	//ignore = system("mkdir -p vtk_particle"); // create folder if not existing


    std::string base_folder = "Alborz_Results";
#if defined Recovery_OFF
#ifdef __linux__
    if (std::filesystem::exists(base_folder)){
        // Command to delete folder recursively
        std::string command1 = "rm -r " + base_folder;

        // Execute the command
        int result1 = system(command1.c_str());
        if (result1 == 0) {
            // Command executed successfully
            std::cout << "Folder '" << base_folder << "' deleted successfully\n";
        } else {
            // Command failed
            std::cerr << "Failed to delete '" << base_folder << "' folder\n";
        }
    }
#else
    // Command to delete folder with /Q option
    std::string command1 = "rmdir /Q /S " + base_folder;
    int result1 = std::system(command1.c_str());

    // Check if the command executed successfully
    if (result1 == 0) {
        // Command executed successfully
        std::cout << "Folder '" << base_folder << "' deleted successfully\n";
    } else {
        // Command failed
        std::cerr << "Failed to delete folder '" << base_folder << "'\n";
    }
	/*std::string command1 = "del /Q ";
	std::string path1 = "Alborz_Results\\vtk_fluid\\*.*";
	system(command1.append(path1).c_str());

	std::string command2 = "del /Q ";
	std::string path2 = "Alborz_Results\\vtk_particle\\*.*";
	system(command2.append(path2).c_str());

	std::string command3 = "del /Q ";
	std::string path3 = "Alborz_Results\\Data_particle\\*.*";
	system(command3.append(path3).c_str());

	std::string command4 = "del /Q ";
	std::string path4 = "Alborz_Results\\*.*";
	system(command4.append(path4).c_str());*/
#endif
#endif

// Create folder
#ifdef __linux__
    if (std::filesystem::create_directory(base_folder)) {
        std::cout << "Folder '" << base_folder << "' created successfully\n";
        // Create subfolders
        std::string vtk_fluid_folder = base_folder + "/vtk_fluid";
        std::string vtk_particle_folder = base_folder + "/vtk_particle";
        std::string data_particle_folder = base_folder + "/Data_particle";
        std::filesystem::create_directory(vtk_fluid_folder);
        std::filesystem::create_directory(vtk_particle_folder);
        std::filesystem::create_directory(data_particle_folder);

    }
	//mkdir("FlowLB_Results", 0777); // Or what parameter you need here ...
	//mkdir("FlowLB_Results/vtk_fluid", 0777); // Or what parameter you need here ...
	//mkdir("FlowLB_Results/vtk_particle", 0777); // Or what parameter you need here ...
	//mkdir("FlowLB_Results/Data_particle", 0777); // Or what parameter you need here ...
#else //for Windows system
    if (!std::filesystem::exists(base_folder)){
        if (_mkdir(base_folder.c_str()) == 0 || errno == EEXIST) {
            std::cout << "Folder '" << base_folder << "' created successfully\n";
            std::string vtk_fluid_folder = base_folder + "\\vtk_fluid";
            std::string vtk_particle_folder = base_folder + "\\vtk_particle";
            std::string data_particle_folder = base_folder + "\\Data_particle";
            _mkdir(vtk_fluid_folder.c_str());
            _mkdir(vtk_particle_folder.c_str());
            _mkdir(data_particle_folder.c_str());
        }
       // _mkdir("Alborz_Results");
       // _mkdir("Alborz_Results/vtk_fluid");
       // _mkdir("Alborz_Results/vtk_particle");
       // _mkdir("Alborz_Results/Data_particle");
	}
#endif

	//#if defined Recovery_OFF
	//	// Below lines reomve existing data_particle and data_fluid and recovery files if "Recovery_OFF" is active
	//	stringstream myfilename;
	//#if defined Flow_With_Particle
	//	for (int npp = 0; npp < num_particles; ++npp) {
	//		string npstr = static_cast<ostringstream*>(&(ostringstream() << npp + 1))->str();//convert number to string
	//		string output_filename("data_particle_");
	//		output_filename += npstr;
	//		output_filename += ".dat";
	//		remove(output_filename.c_str());
	//	}
	//#endif
	//	remove("data_fluid.dat");
	//	remove("data_particle_recovery.dat");
	//	remove("data_fluid_recovery.dat");
	//#endif

	/// Allocate memory for the fluid density, velocity, and force

	density = new double**[Nx];
	velocity_x = new double**[Nx];
	velocity_y = new double**[Nx];
	velocity_z = new double**[Nx];
	force_x = new double**[Nx];
	force_y = new double**[Nx];
	force_z = new double**[Nx];
	is_solid = new int**[Nx];
	mdf_vel_x = new double**[Nx];
	mdf_vel_y = new double**[Nx];
	mdf_vel_z = new double**[Nx];
	mdf_force_x = new double**[Nx];
	mdf_force_y = new double**[Nx];
	mdf_force_z = new double**[Nx];

	temperature = new double**[Nx];
	force_thermal = new double**[Nx];

	triangleset1 = new int[24000];//Total number of triangles is particle_num_nodesx4. Therefore, it must be less than 24000
	triangleset2 = new int[24000];
	triangleset3 = new int[24000];


	for (int X = 0; X < Nx; ++X) {
		density[X] = new double*[Ny];
		velocity_x[X] = new double*[Ny];
		velocity_y[X] = new double*[Ny];
		velocity_z[X] = new double*[Ny];
		force_x[X] = new double*[Ny];
		force_y[X] = new double*[Ny];
		force_z[X] = new double*[Ny];
		is_solid[X] = new int*[Ny];
		mdf_vel_x[X] = new double*[Ny];
		mdf_vel_y[X] = new double*[Ny];
		mdf_vel_z[X] = new double*[Ny];
		mdf_force_x[X] = new double*[Ny];
		mdf_force_y[X] = new double*[Ny];
		mdf_force_z[X] = new double*[Ny];

		temperature[X] = new double*[Ny];
		force_thermal[X] = new double*[Ny];

		for (int Y = 0; Y < Ny; ++Y) {
			density[X][Y] = new double[Nz];
			velocity_x[X][Y] = new double[Nz];
			velocity_y[X][Y] = new double[Nz];
			velocity_z[X][Y] = new double[Nz];
			force_x[X][Y] = new double[Nz];
			force_y[X][Y] = new double[Nz];
			force_z[X][Y] = new double[Nz];
			is_solid[X][Y] = new int[Nz];
			mdf_vel_x[X][Y] = new double[Nz];
			mdf_vel_y[X][Y] = new double[Nz];
			mdf_vel_z[X][Y] = new double[Nz];
			mdf_force_x[X][Y] = new double[Nz];
			mdf_force_y[X][Y] = new double[Nz];
			mdf_force_z[X][Y] = new double[Nz];

			temperature[X][Y] = new double[Nz];
			force_thermal[X][Y] = new double[Nz];


		}
	}


	/// Initialize the fluid density and velocity
	// Start with unit density and zero velocity.
	rec_step = 0;
	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int Z = 0; Z < Nz; ++Z) {
				density[X][Y][Z] = 1;
				velocity_x[X][Y][Z] = 0.;
				velocity_y[X][Y][Z] = 0;
				velocity_z[X][Y][Z] = 0;
				force_x[X][Y][Z] = 0;
				force_y[X][Y][Z] = 0;
				force_z[X][Y][Z] = 0;
				is_solid[X][Y][Z] = FALSE;
				mdf_vel_x[X][Y][Z] = 0;
				mdf_vel_y[X][Y][Z] = 0;
				mdf_vel_z[X][Y][Z] = 0;
				mdf_force_x[X][Y][Z] = 0;
				mdf_force_y[X][Y][Z] = 0;
				mdf_force_z[X][Y][Z] = 0;

				temperature[X][Y][Z] = 0.;
				force_thermal[X][Y][Z] = 0;
			}
		}
	}


	/// Allocate memory for the populations

	pop = new double***[19];
	pop_old = new double***[19];

	pop_t = new double***[19]; // for temperature
	pop_old_t = new double***[19]; // for temperature

	for (int c_i = 0; c_i < 19; ++c_i) {
		pop[c_i] = new double**[Nx];
		pop_old[c_i] = new double**[Nx];
		pop_t[c_i] = new double**[Nx];
		pop_old_t[c_i] = new double**[Nx];
		for (int X = 0; X < Nx; ++X) {
			pop[c_i][X] = new double*[Ny];
			pop_old[c_i][X] = new double*[Ny];
			pop_t[c_i][X] = new double*[Ny];
			pop_old_t[c_i][X] = new double*[Ny];
			for (int Y = 0; Y < Ny; ++Y) {
				pop[c_i][X][Y] = new double[Nz];
				pop_old[c_i][X][Y] = new double[Nz];
				pop_t[c_i][X][Y] = new double[Nz];
				pop_old_t[c_i][X][Y] = new double[Nz];
				for (int Z = 0; Z < Nz; ++Z) {
					pop[c_i][X][Y][Z] = 0;
					pop_old[c_i][X][Y][Z] = 0;
					pop_t[c_i][X][Y][Z] = 0;
					pop_old_t[c_i][X][Y][Z] = 0;
				}
			}
		}
	}


	//#if defined Recovery_OFF
	/// Initialize the populations
	// Use the equilibrium populations corresponding to the initialized fluid density and velocity.

	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int Z = 0; Z < Nz; ++Z) {
				equilibrium(density[X][Y][Z], velocity_x[X][Y][Z], velocity_y[X][Y][Z], velocity_z[X][Y][Z]);
				equilibrium_t(temperature[X][Y][Z], density[X][Y][Z], velocity_x[X][Y][Z], velocity_y[X][Y][Z], velocity_z[X][Y][Z]);

				for (int c_i = 0; c_i < 19; ++c_i) {
					pop_old[c_i][X][Y][Z] = pop_eq[c_i];
					pop[c_i][X][Y][Z] = pop_eq[c_i];
					pop_old_t[c_i][X][Y][Z] = pop_eq_t[c_i];
					pop_t[c_i][X][Y][Z] = pop_eq_t[c_i];
				}
			}
		}
	}
	//#endif

#if defined Recovery_ON
	read_fluid_recovery(rec_step);
#endif

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: EQUILIBRIUM
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the equilibrium populations from the fluid density and velocity.
// It computes the equilibrium only at a specific lattice node: Function has to be called at each lattice node.
// The standard quadratic euilibrium is used.
// reminder: SQ(x) = x * x

void equilibrium(double den, double vel_x, double vel_y, double vel_z) {
	pop_eq[0] = weight[0] * den * (1 - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[1] = weight[1] * den * (1 + 3.0 * (vel_x)+4.5 * SQ(vel_x) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[2] = weight[2] * den * (1 + 3.0 * (-vel_x) + 4.5 * SQ(-vel_x) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[3] = weight[3] * den * (1 + 3.0 * (vel_y)+4.5 * SQ(vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[4] = weight[4] * den * (1 + 3.0 * (-vel_y) + 4.5 * SQ(-vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[5] = weight[5] * den * (1 + 3.0 * (vel_z)+4.5 * SQ(vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[6] = weight[6] * den * (1 + 3.0 * (-vel_z) + 4.5 * SQ(-vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[11] = weight[11] * den * (1 + 3.0 * (vel_x + vel_z) + 4.5 * SQ(vel_x + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[14] = weight[14] * den * (1 + 3.0 * (-vel_x - vel_z) + 4.5 * SQ(-vel_x - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[10] = weight[10] * den * (1 + 3.0 * (-vel_x - vel_y) + 4.5 * SQ(-vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[7] = weight[7] * den * (1 + 3.0 * (vel_x + vel_y) + 4.5 * SQ(vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[12] = weight[12] * den * (1 + 3.0 * (-vel_x + vel_z) + 4.5 * SQ(-vel_x + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[13] = weight[13] * den * (1 + 3.0 * (vel_x - vel_z) + 4.5 * SQ(vel_x - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[9] = weight[9] * den * (1 + 3.0 * (vel_x - vel_y) + 4.5 * SQ(vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[8] = weight[8] * den * (1 + 3.0 * (-vel_x + vel_y) + 4.5 * SQ(-vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[16] = weight[16] * den * (1 + 3.0 * (-vel_y + vel_z) + 4.5 * SQ(-vel_y + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[17] = weight[17] * den * (1 + 3.0 * (+vel_y - vel_z) + 4.5 * SQ(+vel_y - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[18] = weight[18] * den * (1 + 3.0 * (-vel_y - vel_z) + 4.5 * SQ(-vel_y - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq[15] = weight[15] * den * (1 + 3.0 * (+vel_y + vel_z) + 4.5 * SQ(+vel_y + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: EQUILIBRIUM-TEMPERATURE
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the equilibrium populations from the fluid temperature.
// It computes the equilibrium only at a specific lattice node: Function has to be called at each lattice node.
// The standard quadratic euilibrium is used.
void equilibrium_t(double tem, double den, double vel_x, double vel_y, double vel_z) {
#if defined Equilibrium_Thermal_1
	pop_eq_t[0] = weight[0] * tem * (1);
	pop_eq_t[1] = weight[1] * tem * (1 + 3 * (vel_x));
	pop_eq_t[2] = weight[2] * tem * (1 + 3.0 * (-vel_x));
	pop_eq_t[3] = weight[3] * tem * (1 + 3.0 * (vel_y));
	pop_eq_t[4] = weight[4] * tem * (1 + 3.0 * (-vel_y));
	pop_eq_t[5] = weight[5] * tem * (1 + 3.0 * (vel_z));
	pop_eq_t[6] = weight[6] * tem * (1 + 3.0 * (-vel_z));
	pop_eq_t[11] = weight[11] * tem * (1 + 3.0 * (vel_x + vel_z));
	pop_eq_t[14] = weight[14] * tem * (1 + 3.0 * (-vel_x - vel_z));
	pop_eq_t[10] = weight[10] * tem * (1 + 3.0 * (-vel_x - vel_y));
	pop_eq_t[7] = weight[7] * tem * (1 + 3.0 * (vel_x + vel_y));
	pop_eq_t[12] = weight[12] * tem * (1 + 3.0 * (-vel_x + vel_z));
	pop_eq_t[13] = weight[13] * tem * (1 + 3.0 * (vel_x - vel_z));
	pop_eq_t[9] = weight[9] * tem * (1 + 3.0 * (vel_x - vel_y));
	pop_eq_t[8] = weight[8] * tem * (1 + 3.0 * (-vel_x + vel_y));
	pop_eq_t[16] = weight[16] * tem * (1 + 3.0 * (-vel_y + vel_z));
	pop_eq_t[17] = weight[17] * tem * (1 + 3.0 * (+vel_y - vel_z));
	pop_eq_t[18] = weight[18] * tem * (1 + 3.0 * (-vel_y - vel_z));
	pop_eq_t[15] = weight[15] * tem * (1 + 3.0 * (+vel_y + vel_z));
#endif
#if defined Equilibrium_Thermal_2
	pop_eq_t[0] = weight[0] * tem* den * (1 - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[1] = weight[1] * tem* den * (1 + 3.0 * (vel_x)+4.5 * SQ(vel_x) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[2] = weight[2] * tem* den * (1 + 3.0 * (-vel_x) + 4.5 * SQ(-vel_x) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[3] = weight[3] * tem* den * (1 + 3.0 * (vel_y)+4.5 * SQ(vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[4] = weight[4] * tem* den * (1 + 3.0 * (-vel_y) + 4.5 * SQ(-vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[5] = weight[5] * tem* den * (1 + 3.0 * (vel_z)+4.5 * SQ(vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[6] = weight[6] * tem* den * (1 + 3.0 * (-vel_z) + 4.5 * SQ(-vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[11] = weight[11] * tem* den * (1 + 3.0 * (vel_x + vel_z) + 4.5 * SQ(vel_x + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[14] = weight[14] * tem* den * (1 + 3.0 * (-vel_x - vel_z) + 4.5 * SQ(-vel_x - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[10] = weight[10] * tem* den * (1 + 3.0 * (-vel_x - vel_y) + 4.5 * SQ(-vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[7] = weight[7] * tem* den * (1 + 3.0 * (vel_x + vel_y) + 4.5 * SQ(vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[12] = weight[12] * tem* den * (1 + 3.0 * (-vel_x + vel_z) + 4.5 * SQ(-vel_x + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[13] = weight[13] * tem* den * (1 + 3.0 * (vel_x - vel_z) + 4.5 * SQ(vel_x - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[9] = weight[9] * tem* den * (1 + 3.0 * (vel_x - vel_y) + 4.5 * SQ(vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[8] = weight[8] * tem* den * (1 + 3.0 * (-vel_x + vel_y) + 4.5 * SQ(-vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[16] = weight[16] * tem* den * (1 + 3.0 * (-vel_y + vel_z) + 4.5 * SQ(-vel_y + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[17] = weight[17] * tem* den * (1 + 3.0 * (+vel_y - vel_z) + 4.5 * SQ(+vel_y - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[18] = weight[18] * tem* den * (1 + 3.0 * (-vel_y - vel_z) + 4.5 * SQ(-vel_y - vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
	pop_eq_t[15] = weight[15] * tem* den * (1 + 3.0 * (+vel_y + vel_z) + 4.5 * SQ(+vel_y + vel_z) - 1.5 * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));

#endif
	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM
//-----------------------------------------------------------------------------------------------------------------------------
//
void LBM(int time) {

	/// Swap populations
	// The present code used old and new populations which are swapped at the beginning of each time step.
	// This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
	// This way, the old populations are not overwritten during propagation.
	// The resulting code is easier to write and to debug.
	// The memory requirement for the populations is twice as large.
	double ****swap_temp = pop_old;
	pop_old = pop;
	pop = swap_temp;
	double ssf = 1. / tau;
	double s[19] = { 0.0, 1.19, 1.4, 0.0, 1.2, 0.0, 1.2, 0.0, 1.2, ssf, 1.4, ssf, 1.4, ssf, ssf, ssf, 1.98, 1.98, 1.98 };//As sugested in page 103 of Mohamad book
																														 //double s[19] = {ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf};//If all coefficients are 1/tau then we have SRT
	double MEQ;
	double m[19];
	double deltam[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double mforce[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double fforce[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };



	/// Lattice Boltzmann equation
	// The lattice Boltzmann equation is solved in the following.
	// The algorithm includes
	// - computation of the lattice force
	// - combined collision and propagation (faster than first collision and then propagation)

	int start_X, end_X;//Start and End of LBM calculation for processor.
	start_X = 0;
	end_X = Nx - 1;

#if defined X_Wall
	start_X = 1;
	end_X = Nx - 2;
#endif
#if defined X_velIn_ZeroGradOut
	start_X = 1;
	end_X = Nx - 2;
#endif
#if defined ZeroGradIn_ZeroGradOut
	start_X = 1;
	end_X = Nx - 2;
#endif

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

#if defined Y_Wall
	start_Y = 1;
	end_Y = Ny - 2;
#endif

	int start_Z, end_Z;
	start_Z = 0;
	end_Z = Nz - 1;

#if defined Z_Wall
	start_Z = 1;
	end_Z = Nz - 2;
#endif


	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			for (int Z = start_Z; Z <= end_Z; ++Z) {
				if (is_solid[X][Y][Z] == FALSE) {

					/// Compute equilibrium
					// The equilibrium populations are computed.
					// Forces are coupled via Guo model or shan-chen model

#if defined Shan_Chen
					equilibrium(density[X][Y][Z], (velocity_x[X][Y][Z] + (force_x[X][Y][Z] + gravity) * tau), (velocity_y[X][Y][Z] + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * tau), (velocity_z[X][Y][Z] + (force_z[X][Y][Z]) * tau));
#endif
#if defined Guo || defined Luo
					equilibrium(density[X][Y][Z], velocity_x[X][Y][Z], velocity_y[X][Y][Z], velocity_z[X][Y][Z]);
#endif

					/// Compute new populations
					// This is the lattice Boltzmann equation (combined collision and propagation) including external forcing.
					// Periodicity of the lattice in x-direction is taken into account by the %-operator.

#if defined Shan_Chen
					pop[0][X][Y][Z] = pop_old[0][X][Y][Z] * (1.0 - omega) + pop_eq[0] * omega;
					pop[1][(X + 1) % Nx][Y][Z] = pop_old[1][X][Y][Z] * (1.0 - omega) + pop_eq[1] * omega;
					pop[2][(X - 1 + Nx) % Nx][Y][Z] = pop_old[2][X][Y][Z] * (1.0 - omega) + pop_eq[2] * omega;
					pop[3][X][(Y + 1) % Ny][Z] = pop_old[3][X][Y][Z] * (1.0 - omega) + pop_eq[3] * omega;
					pop[4][X][(Y - 1 + Ny) % Ny][Z] = pop_old[4][X][Y][Z] * (1.0 - omega) + pop_eq[4] * omega;
					pop[5][X][Y][(Z + 1) % Nz] = pop_old[5][X][Y][Z] * (1.0 - omega) + pop_eq[5] * omega;
					pop[6][X][Y][(Z - 1 + Nz) % Nz] = pop_old[6][X][Y][Z] * (1.0 - omega) + pop_eq[6] * omega;
					pop[11][(X + 1) % Nx][Y][(Z + 1) % Nz] = pop_old[11][X][Y][Z] * (1.0 - omega) + pop_eq[11] * omega;
					pop[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old[14][X][Y][Z] * (1.0 - omega) + pop_eq[14] * omega;
					pop[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old[10][X][Y][Z] * (1.0 - omega) + pop_eq[10] * omega;
					pop[7][(X + 1) % Nx][(Y + 1) % Ny][Z] = pop_old[7][X][Y][Z] * (1.0 - omega) + pop_eq[7] * omega;
					pop[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] = pop_old[12][X][Y][Z] * (1.0 - omega) + pop_eq[12] * omega;
					pop[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old[13][X][Y][Z] * (1.0 - omega) + pop_eq[13] * omega;
					pop[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old[9][X][Y][Z] * (1.0 - omega) + pop_eq[9] * omega;
					pop[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] = pop_old[8][X][Y][Z] * (1.0 - omega) + pop_eq[8] * omega;
					pop[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] = pop_old[16][X][Y][Z] * (1.0 - omega) + pop_eq[16] * omega;
					pop[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] = pop_old[17][X][Y][Z] * (1.0 - omega) + pop_eq[17] * omega;
					pop[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] = pop_old[18][X][Y][Z] * (1.0 - omega) + pop_eq[18] * omega;
					pop[15][X][(Y + 1) % Ny][(Z + 1) % Nz] = pop_old[15][X][Y][Z] * (1.0 - omega) + pop_eq[15] * omega;
#endif

#if defined Luo
					pop[0][X][Y][Z] = pop_old[0][X][Y][Z] * (1.0 - omega) + pop_eq[0] * omega;
					pop[1][(X + 1) % Nx][Y][Z] = pop_old[1][X][Y][Z] * (1.0 - omega) + pop_eq[1] * omega + 3 * weight[1] * (1 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[2][(X - 1 + Nx) % Nx][Y][Z] = pop_old[2][X][Y][Z] * (1.0 - omega) + pop_eq[2] * omega + 3 * weight[2] * (-1 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[3][X][(Y + 1) % Ny][Z] = pop_old[3][X][Y][Z] * (1.0 - omega) + pop_eq[3] * omega + 3 * weight[3] * (0 * (force_x[X][Y][Z] + gravity) + 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[4][X][(Y - 1 + Ny) % Ny][Z] = pop_old[4][X][Y][Z] * (1.0 - omega) + pop_eq[4] * omega + 3 * weight[4] * (0 * (force_x[X][Y][Z] + gravity) - 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[5][X][Y][(Z + 1) % Nz] = pop_old[5][X][Y][Z] * (1.0 - omega) + pop_eq[5] * omega + 3 * weight[5] * (0 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[6][X][Y][(Z - 1 + Nz) % Nz] = pop_old[6][X][Y][Z] * (1.0 - omega) + pop_eq[6] * omega + 3 * weight[6] * (0 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) - 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[11][(X + 1) % Nx][Y][(Z + 1) % Nz] = pop_old[11][X][Y][Z] * (1.0 - omega) + pop_eq[11] * omega + 3 * weight[11] * (1 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old[14][X][Y][Z] * (1.0 - omega) + pop_eq[14] * omega + 3 * weight[14] * (-1 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) - 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old[10][X][Y][Z] * (1.0 - omega) + pop_eq[10] * omega + 3 * weight[10] * (-1 * (force_x[X][Y][Z] + gravity) - 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[7][(X + 1) % Nx][(Y + 1) % Ny][Z] = pop_old[7][X][Y][Z] * (1.0 - omega) + pop_eq[7] * omega + 3 * weight[7] * (1 * (force_x[X][Y][Z] + gravity) + 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] = pop_old[12][X][Y][Z] * (1.0 - omega) + pop_eq[12] * omega + 3 * weight[12] * (-1 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old[13][X][Y][Z] * (1.0 - omega) + pop_eq[13] * omega + 3 * weight[13] * (1 * (force_x[X][Y][Z] + gravity) + 0 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) - 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old[9][X][Y][Z] * (1.0 - omega) + pop_eq[9] * omega + 3 * weight[9] * (1 * (force_x[X][Y][Z] + gravity) - 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] = pop_old[8][X][Y][Z] * (1.0 - omega) + pop_eq[8] * omega + 3 * weight[8] * (-1 * (force_x[X][Y][Z] + gravity) + 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 0 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] = pop_old[16][X][Y][Z] * (1.0 - omega) + pop_eq[16] * omega + 3 * weight[16] * (0 * (force_x[X][Y][Z] + gravity) - 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] = pop_old[17][X][Y][Z] * (1.0 - omega) + pop_eq[17] * omega + 3 * weight[17] * (0 * (force_x[X][Y][Z] + gravity) + 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) - 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] = pop_old[18][X][Y][Z] * (1.0 - omega) + pop_eq[18] * omega + 3 * weight[18] * (0 * (force_x[X][Y][Z] + gravity) - 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) - 1 * force_z[X][Y][Z])* density[X][Y][Z];
					pop[15][X][(Y + 1) % Ny][(Z + 1) % Nz] = pop_old[15][X][Y][Z] * (1.0 - omega) + pop_eq[15] * omega + 3 * weight[15] * (0 * (force_x[X][Y][Z] + gravity) + 1 * (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) + 1 * force_z[X][Y][Z])* density[X][Y][Z];
#endif

#if defined Guo
#if defined SRT
					pop[0][X][Y][Z] = pop_old[0][X][Y][Z] * (1.0 - omega) + pop_eq[0] * omega + (1.0 - 0.5*omega)* weight[0] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[1][(X + 1) % Nx][Y][Z] = pop_old[1][X][Y][Z] * (1.0 - omega) + pop_eq[1] * omega + (1.0 - 0.5*omega)* weight[1] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[2][(X - 1 + Nx) % Nx][Y][Z] = pop_old[2][X][Y][Z] * (1.0 - omega) + pop_eq[2] * omega + (1.0 - 0.5*omega)* weight[2] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[3][X][(Y + 1) % Ny][Z] = pop_old[3][X][Y][Z] * (1.0 - omega) + pop_eq[3] * omega + (1.0 - 0.5*omega)* weight[3] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[4][X][(Y - 1 + Ny) % Ny][Z] = pop_old[4][X][Y][Z] * (1.0 - omega) + pop_eq[4] * omega + (1.0 - 0.5*omega)* weight[4] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[5][X][Y][(Z + 1) % Nz] = pop_old[5][X][Y][Z] * (1.0 - omega) + pop_eq[5] * omega + (1.0 - 0.5*omega)* weight[5] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[6][X][Y][(Z - 1 + Nz) % Nz] = pop_old[6][X][Y][Z] * (1.0 - omega) + pop_eq[6] * omega + (1.0 - 0.5*omega)* weight[6] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (-1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[11][(X + 1) % Nx][Y][(Z + 1) % Nz] = pop_old[11][X][Y][Z] * (1.0 - omega) + pop_eq[11] * omega + (1.0 - 0.5*omega)* weight[11] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old[14][X][Y][Z] * (1.0 - omega) + pop_eq[14] * omega + (1.0 - 0.5*omega)* weight[14] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old[10][X][Y][Z] * (1.0 - omega) + pop_eq[10] * omega + (1.0 - 0.5*omega)* weight[10] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[7][(X + 1) % Nx][(Y + 1) % Ny][Z] = pop_old[7][X][Y][Z] * (1.0 - omega) + pop_eq[7] * omega + (1.0 - 0.5*omega)* weight[7] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] = pop_old[12][X][Y][Z] * (1.0 - omega) + pop_eq[12] * omega + (1.0 - 0.5*omega)* weight[12] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old[13][X][Y][Z] * (1.0 - omega) + pop_eq[13] * omega + (1.0 - 0.5*omega)* weight[13] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old[9][X][Y][Z] * (1.0 - omega) + pop_eq[9] * omega + (1.0 - 0.5*omega)* weight[9] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] = pop_old[8][X][Y][Z] * (1.0 - omega) + pop_eq[8] * omega + (1.0 - 0.5*omega)* weight[8] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] = pop_old[16][X][Y][Z] * (1.0 - omega) + pop_eq[16] * omega + (1.0 - 0.5*omega)* weight[16] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (-1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] = pop_old[17][X][Y][Z] * (1.0 - omega) + pop_eq[17] * omega + (1.0 - 0.5*omega)* weight[17] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] = pop_old[18][X][Y][Z] * (1.0 - omega) + pop_eq[18] * omega + (1.0 - 0.5*omega)* weight[18] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					pop[15][X][(Y + 1) % Ny][(Z + 1) % Nz] = pop_old[15][X][Y][Z] * (1.0 - omega) + pop_eq[15] * omega + (1.0 - 0.5*omega)* weight[15] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
#endif
#if defined MRT
					// Transformation from velocity space to moment space:
					m[0] = pop_old[0][X][Y][Z] + pop_old[1][X][Y][Z] + pop_old[2][X][Y][Z] + pop_old[3][X][Y][Z] + pop_old[4][X][Y][Z] + pop_old[5][X][Y][Z] + pop_old[6][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] + pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] + pop_old[13][X][Y][Z] + pop_old[14][X][Y][Z] + pop_old[15][X][Y][Z] + pop_old[16][X][Y][Z] + pop_old[17][X][Y][Z] + pop_old[18][X][Y][Z];
					m[1] = -30 * pop_old[0][X][Y][Z] - 11 * pop_old[1][X][Y][Z] - 11 * pop_old[2][X][Y][Z] - 11 * pop_old[3][X][Y][Z] - 11 * pop_old[4][X][Y][Z] - 11 * pop_old[5][X][Y][Z] - 11 * pop_old[6][X][Y][Z] + 8 * pop_old[7][X][Y][Z] + 8 * pop_old[8][X][Y][Z] + 8 * pop_old[9][X][Y][Z] + 8 * pop_old[10][X][Y][Z] + 8 * pop_old[11][X][Y][Z] + 8 * pop_old[12][X][Y][Z] + 8 * pop_old[13][X][Y][Z] + 8 * pop_old[14][X][Y][Z] + 8 * pop_old[15][X][Y][Z] + 8 * pop_old[16][X][Y][Z] + 8 * pop_old[17][X][Y][Z] + 8 * pop_old[18][X][Y][Z];
					m[2] = 12 * pop_old[0][X][Y][Z] - 4 * pop_old[1][X][Y][Z] - 4 * pop_old[2][X][Y][Z] - 4 * pop_old[3][X][Y][Z] - 4 * pop_old[4][X][Y][Z] - 4 * pop_old[5][X][Y][Z] - 4 * pop_old[6][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] + pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] + pop_old[13][X][Y][Z] + pop_old[14][X][Y][Z] + pop_old[15][X][Y][Z] + pop_old[16][X][Y][Z] + pop_old[17][X][Y][Z] + pop_old[18][X][Y][Z];
					m[3] = pop_old[1][X][Y][Z] - pop_old[2][X][Y][Z] + pop_old[7][X][Y][Z] - pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] - pop_old[10][X][Y][Z] + pop_old[11][X][Y][Z] - pop_old[12][X][Y][Z] + pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z];
					m[4] = -4 * pop_old[1][X][Y][Z] + 4 * pop_old[2][X][Y][Z] + pop_old[7][X][Y][Z] - pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] - pop_old[10][X][Y][Z] + pop_old[11][X][Y][Z] - pop_old[12][X][Y][Z] + pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z];
					m[5] = pop_old[3][X][Y][Z] - pop_old[4][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] - pop_old[9][X][Y][Z] - pop_old[10][X][Y][Z] + pop_old[15][X][Y][Z] - pop_old[16][X][Y][Z] + pop_old[17][X][Y][Z] - pop_old[18][X][Y][Z];
					m[6] = -4 * pop_old[3][X][Y][Z] + 4 * pop_old[4][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] - pop_old[9][X][Y][Z] - pop_old[10][X][Y][Z] + pop_old[15][X][Y][Z] - pop_old[16][X][Y][Z] + pop_old[17][X][Y][Z] - pop_old[18][X][Y][Z];
					m[7] = pop_old[5][X][Y][Z] - pop_old[6][X][Y][Z] + pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z] + pop_old[15][X][Y][Z] + pop_old[16][X][Y][Z] - pop_old[17][X][Y][Z] - pop_old[18][X][Y][Z];
					m[8] = -4 * pop_old[5][X][Y][Z] + 4 * pop_old[6][X][Y][Z] + pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z] + pop_old[15][X][Y][Z] + pop_old[16][X][Y][Z] - pop_old[17][X][Y][Z] - pop_old[18][X][Y][Z];
					m[9] = 2 * pop_old[1][X][Y][Z] + 2 * pop_old[2][X][Y][Z] - pop_old[3][X][Y][Z] - pop_old[4][X][Y][Z] - pop_old[5][X][Y][Z] - pop_old[6][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] + pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] + pop_old[13][X][Y][Z] + pop_old[14][X][Y][Z] - 2 * pop_old[15][X][Y][Z] - 2 * pop_old[16][X][Y][Z] - 2 * pop_old[17][X][Y][Z] - 2 * pop_old[18][X][Y][Z];
					m[10] = -4 * pop_old[1][X][Y][Z] - 4 * pop_old[2][X][Y][Z] + 2 * pop_old[3][X][Y][Z] + 2 * pop_old[4][X][Y][Z] + 2 * pop_old[5][X][Y][Z] + 2 * pop_old[6][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] + pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] + pop_old[13][X][Y][Z] + pop_old[14][X][Y][Z] - 2 * pop_old[15][X][Y][Z] - 2 * pop_old[16][X][Y][Z] - 2 * pop_old[17][X][Y][Z] - 2 * pop_old[18][X][Y][Z];
					m[11] = pop_old[3][X][Y][Z] + pop_old[4][X][Y][Z] - pop_old[5][X][Y][Z] - pop_old[6][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] - pop_old[11][X][Y][Z] - pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z];
					m[12] = -2 * pop_old[3][X][Y][Z] - 2 * pop_old[4][X][Y][Z] + 2 * pop_old[5][X][Y][Z] + 2 * pop_old[6][X][Y][Z] + pop_old[7][X][Y][Z] + pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] - pop_old[11][X][Y][Z] - pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z];
					m[13] = pop_old[7][X][Y][Z] - pop_old[8][X][Y][Z] - pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z];
					m[14] = pop_old[15][X][Y][Z] - pop_old[16][X][Y][Z] - pop_old[17][X][Y][Z] + pop_old[18][X][Y][Z];
					m[15] = pop_old[11][X][Y][Z] - pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] + pop_old[14][X][Y][Z];
					m[16] = pop_old[7][X][Y][Z] - pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] - pop_old[10][X][Y][Z] - pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] + pop_old[14][X][Y][Z];
					m[17] = -1 * pop_old[7][X][Y][Z] - pop_old[8][X][Y][Z] + pop_old[9][X][Y][Z] + pop_old[10][X][Y][Z] + pop_old[15][X][Y][Z] - pop_old[16][X][Y][Z] + pop_old[17][X][Y][Z] - pop_old[18][X][Y][Z];
					m[18] = pop_old[11][X][Y][Z] + pop_old[12][X][Y][Z] - pop_old[13][X][Y][Z] - pop_old[14][X][Y][Z] - pop_old[15][X][Y][Z] - pop_old[16][X][Y][Z] + pop_old[17][X][Y][Z] + pop_old[18][X][Y][Z];


					fforce[0] = weight[0] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[1] = weight[1] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[2] = weight[2] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[3] = weight[3] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[4] = weight[4] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[5] = weight[5] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[6] = weight[6] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (-1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[11] = weight[11] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[14] = weight[14] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[10] = weight[10] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[7] = weight[7] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[12] = weight[12] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[13] = weight[13] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (0 - velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[9] = weight[9] * ((force_x[X][Y][Z] + gravity)*(3 * (1 - velocity_x[X][Y][Z]) + 9.0 * 1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (1 * velocity_x[X][Y][Z] + -1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[8] = weight[8] * ((force_x[X][Y][Z] + gravity)*(3 * (-1 - velocity_x[X][Y][Z]) + 9.0 * -1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (-1 * velocity_x[X][Y][Z] + 1 * velocity_y[X][Y][Z])) + force_z[X][Y][Z] * (3 * (0 - velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[16] = weight[16] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (-1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[17] = weight[17] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[18] = weight[18] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (-1 - velocity_y[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (-1 - velocity_z[X][Y][Z]) + 9 * -1 * (-1 * velocity_y[X][Y][Z] + -1 * velocity_z[X][Y][Z])))*density[X][Y][Z];
					fforce[15] = weight[15] * ((force_x[X][Y][Z] + gravity)*(3 * (0 - velocity_x[X][Y][Z])) + (force_y[X][Y][Z] + gbeta*(temperature[X][Y][Z] - Tref)) * (3 * (1 - velocity_y[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])) + force_z[X][Y][Z] * (3 * (1 - velocity_z[X][Y][Z]) + 9 * 1 * (1 * velocity_y[X][Y][Z] + 1 * velocity_z[X][Y][Z])))*density[X][Y][Z];

					mforce[0] = fforce[0] + fforce[1] + fforce[2] + fforce[3] + fforce[4] + fforce[5] + fforce[6] + fforce[7] + fforce[8] + fforce[9] + fforce[10] + fforce[11] + fforce[12] + fforce[13] + fforce[14] + fforce[15] + fforce[16] + fforce[17] + fforce[18];
					mforce[1] = -30 * fforce[0] - 11 * fforce[1] - 11 * fforce[2] - 11 * fforce[3] - 11 * fforce[4] - 11 * fforce[5] - 11 * fforce[6] + 8 * fforce[7] + 8 * fforce[8] + 8 * fforce[9] + 8 * fforce[10] + 8 * fforce[11] + 8 * fforce[12] + 8 * fforce[13] + 8 * fforce[14] + 8 * fforce[15] + 8 * fforce[16] + 8 * fforce[17] + 8 * fforce[18];
					mforce[2] = 12 * fforce[0] - 4 * fforce[1] - 4 * fforce[2] - 4 * fforce[3] - 4 * fforce[4] - 4 * fforce[5] - 4 * fforce[6] + fforce[7] + fforce[8] + fforce[9] + fforce[10] + fforce[11] + fforce[12] + fforce[13] + fforce[14] + fforce[15] + fforce[16] + fforce[17] + fforce[18];
					mforce[3] = fforce[1] - fforce[2] + fforce[7] - fforce[8] + fforce[9] - fforce[10] + fforce[11] - fforce[12] + fforce[13] - fforce[14];
					mforce[4] = -4 * fforce[1] + 4 * fforce[2] + fforce[7] - fforce[8] + fforce[9] - fforce[10] + fforce[11] - fforce[12] + fforce[13] - fforce[14];
					mforce[5] = fforce[3] - fforce[4] + fforce[7] + fforce[8] - fforce[9] - fforce[10] + fforce[15] - fforce[16] + fforce[17] - fforce[18];
					mforce[6] = -4 * fforce[3] + 4 * fforce[4] + fforce[7] + fforce[8] - fforce[9] - fforce[10] + fforce[15] - fforce[16] + fforce[17] - fforce[18];
					mforce[7] = fforce[5] - fforce[6] + fforce[11] + fforce[12] - fforce[13] - fforce[14] + fforce[15] + fforce[16] - fforce[17] - fforce[18];
					mforce[8] = -4 * fforce[5] + 4 * fforce[6] + fforce[11] + fforce[12] - fforce[13] - fforce[14] + fforce[15] + fforce[16] - fforce[17] - fforce[18];
					mforce[9] = 2 * fforce[1] + 2 * fforce[2] - fforce[3] - fforce[4] - fforce[5] - fforce[6] + fforce[7] + fforce[8] + fforce[9] + fforce[10] + fforce[11] + fforce[12] + fforce[13] + fforce[14] - 2 * fforce[15] - 2 * fforce[16] - 2 * fforce[17] - 2 * fforce[18];
					mforce[10] = -4 * fforce[1] - 4 * fforce[2] + 2 * fforce[3] + 2 * fforce[4] + 2 * fforce[5] + 2 * fforce[6] + fforce[7] + fforce[8] + fforce[9] + fforce[10] + fforce[11] + fforce[12] + fforce[13] + fforce[14] - 2 * fforce[15] - 2 * fforce[16] - 2 * fforce[17] - 2 * fforce[18];
					mforce[11] = fforce[3] + fforce[4] - fforce[5] - fforce[6] + fforce[7] + fforce[8] + fforce[9] + fforce[10] - fforce[11] - fforce[12] - fforce[13] - fforce[14];
					mforce[12] = -2 * fforce[3] - 2 * fforce[4] + 2 * fforce[5] + 2 * fforce[6] + fforce[7] + fforce[8] + fforce[9] + fforce[10] - fforce[11] - fforce[12] - fforce[13] - fforce[14];
					mforce[13] = fforce[7] - fforce[8] - fforce[9] + fforce[10];
					mforce[14] = fforce[15] - fforce[16] - fforce[17] + fforce[18];
					mforce[15] = fforce[11] - fforce[12] - fforce[13] + fforce[14];
					mforce[16] = fforce[7] - fforce[8] + fforce[9] - fforce[10] - fforce[11] + fforce[12] - fforce[13] + fforce[14];
					mforce[17] = -1 * fforce[7] - fforce[8] + fforce[9] + fforce[10] + fforce[15] - fforce[16] + fforce[17] - fforce[18];
					mforce[18] = fforce[11] + fforce[12] - fforce[13] - fforce[14] - fforce[15] - fforce[16] + fforce[17] + fforce[18];

					for (int k = 0; k < 19; k++)
					{
						MEQ = meq(density[X][Y][Z], velocity_x[X][Y][Z], velocity_y[X][Y][Z], velocity_z[X][Y][Z], k);
						m[k] = m[k] - s[k] * (m[k] - MEQ) + (1 - 1. / 2.*s[k]) * mforce[k]; // relaxation
						mforce[k] = 0.0;
						fforce[k] = 0.0;
					}

					pop[0][X][Y][Z] = 0.0526315789473684 * m[0] - 0.0125313283208020	* m[1] + 0.0476190476190476 * m[2];
					pop[1][(X + 1) % Nx][Y][Z] = 0.0526315789473684	* m[0] - 0.00459482038429407* m[1] - 0.0158730158730159	* m[2] + 0.1 * m[3] - 0.1   * m[4] + 0.0555555555555556 * m[9] - 0.0555555555555556 * m[10];
					pop[2][(X - 1 + Nx) % Nx][Y][Z] = 0.0526315789473684	* m[0] - 0.00459482038429407* m[1] - 0.0158730158730159	* m[2] - 0.1 * m[3] + 0.1   * m[4] + 0.0555555555555556 * m[9] - 0.0555555555555556 * m[10];
					pop[3][X][(Y + 1) % Ny][Z] = 0.0526315789473684	* m[0] - 0.00459482038429407* m[1] - 0.0158730158730159	* m[2] + 0.1 * m[5] - 0.1   * m[6] - 0.0277777777777778 * m[9] + 0.0277777777777778 * m[10] + 0.0833333333333333 * m[11] - 0.0833333333333333 * m[12];
					pop[4][X][(Y - 1 + Ny) % Ny][Z] = 0.0526315789473684	* m[0] - 0.00459482038429407* m[1] - 0.0158730158730159 * m[2] - 0.1 * m[5] + 0.1   * m[6] - 0.0277777777777778 * m[9] + 0.0277777777777778 * m[10] + 0.0833333333333333 * m[11] - 0.0833333333333333 * m[12];
					pop[5][X][Y][(Z + 1) % Nz] = 0.0526315789473684	* m[0] - 0.00459482038429407* m[1] - 0.0158730158730159 * m[2] + 0.1 * m[7] - 0.1   * m[8] - 0.0277777777777778 * m[9] + 0.0277777777777778 * m[10] - 0.0833333333333333 * m[11] + 0.0833333333333333 * m[12];
					pop[6][X][Y][(Z - 1 + Nz) % Nz] = 0.0526315789473684	* m[0] - 0.00459482038429407* m[1] - 0.0158730158730159 * m[2] - 0.1 * m[7] + 0.1   * m[8] - 0.0277777777777778 * m[9] + 0.0277777777777778 * m[10] - 0.0833333333333333 * m[11] + 0.0833333333333333 * m[12];
					pop[7][(X + 1) % Nx][(Y + 1) % Ny][Z] = 0.0526315789473685	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] + 0.1 * m[3] + 0.025 * m[4] + 0.1 * m[5] + 0.025 * m[6] + 0.0277777777777778 * m[9] + 0.0138888888888889  * m[10] + 0.0833333333333333 * m[11] + 0.0416666666666667 * m[12] + 0.25 * m[13] + 0.125 * m[16] - 0.125 * m[17];
					pop[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] = 0.0526315789473685 * m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] - 0.1 * m[3] - 0.025 * m[4] + 0.1 * m[5] + 0.025 * m[6] + 0.0277777777777778	* m[9] + 0.0138888888888889  * m[10] + 0.0833333333333333 * m[11] + 0.0416666666666667 * m[12] - 0.25 * m[13] - 0.125 * m[16] - 0.125 * m[17];
					pop[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] = 0.0526315789473685	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] + 0.1 * m[3] + 0.025 * m[4] - 0.1 * m[5] - 0.025 * m[6] + 0.0277777777777778 * m[9] + 0.0138888888888889  * m[10] + 0.0833333333333333 * m[11] + 0.0416666666666667 * m[12] - 0.25 * m[13] + 0.125 * m[16] + 0.125 * m[17];
					pop[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] = 0.0526315789473685	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] - 0.1 * m[3] - 0.025 * m[4] - 0.1 * m[5] - 0.025 * m[6] + 0.0277777777777778 * m[9] + 0.0138888888888889  * m[10] + 0.0833333333333333 * m[11] + 0.0416666666666667 * m[12] + 0.25 * m[13] - 0.125 * m[16] + 0.125 * m[17];
					pop[11][(X + 1) % Nx][Y][(Z + 1) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] + 0.1 * m[3] + 0.025 * m[4] + 0.1 * m[7] + 0.025 * m[8] + 0.0277777777777778 * m[9] + 0.0138888888888889  * m[10] - 0.0833333333333333 * m[11] - 0.0416666666666667 * m[12] + 0.25 * m[15] - 0.125 * m[16] + 0.125 * m[18];
					pop[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] - 0.1 * m[3] - 0.025 * m[4] + 0.1 * m[7] + 0.025 * m[8] + 0.0277777777777778 * m[9] + 0.0138888888888889	 * m[10] - 0.0833333333333333 * m[11] - 0.0416666666666667 * m[12] - 0.25 * m[15] + 0.125 * m[16] + 0.125 * m[18];
					pop[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] = 0.0526315789473684 * m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] + 0.1 * m[3] + 0.025 * m[4] - 0.1 * m[7] - 0.025 * m[8] + 0.0277777777777778 * m[9] + 0.0138888888888889  * m[10] - 0.0833333333333333 * m[11] - 0.0416666666666667 * m[12] - 0.25 * m[15] - 0.125 * m[16] - 0.125 * m[18];
					pop[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] - 0.1 * m[3] - 0.025 * m[4] - 0.1 * m[7] - 0.025 * m[8] + 0.0277777777777778 * m[9] + 0.0138888888888889  * m[10] - 0.0833333333333333 * m[11] - 0.0416666666666667 * m[12] + 0.25 * m[15] + 0.125 * m[16] - 0.125 * m[18];
					pop[15][X][(Y + 1) % Ny][(Z + 1) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] + 0.1 * m[5] + 0.025 * m[6] + 0.1 * m[7] + 0.025 * m[8] - 0.0555555555555556 * m[9] - 0.0277777777777778  * m[10] + 0.25 * m[14] + 0.125 * m[17] - 0.125 * m[18];
					pop[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] - 0.1 * m[5] - 0.025 * m[6] + 0.1 * m[7] + 0.025 * m[8] - 0.0555555555555556 * m[9] - 0.0277777777777778  * m[10] - 0.25 * m[14] - 0.125 * m[17] - 0.125 * m[18];
					pop[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] + 0.1 * m[5] + 0.025 * m[6] - 0.1 * m[7] - 0.025 * m[8] - 0.0555555555555556 * m[9] - 0.0277777777777778	 * m[10] - 0.25 * m[14] + 0.125 * m[17] + 0.125 * m[18];
					pop[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] = 0.0526315789473684	* m[0] + 0.00334168755221387* m[1] + 0.00396825396825397* m[2] - 0.1 * m[5] - 0.025 * m[6] - 0.1 * m[7] - 0.025 * m[8] - 0.0555555555555556	* m[9] - 0.0277777777777778  * m[10] + 0.25 * m[14] - 0.125 * m[17] + 0.125 * m[18];
#endif //End of MRT
#endif //End of GUO
				} //if (is_solid[X][Y][Z] == FALSE)
			} //for (int X = start_X; X <= end_X; ++X)
		} // for (int Y = start_Y; Y <= end_Y; ++Y)
	} //  for (int Z = start_Z; Z <= end_Z; ++Z)


	  // Bounce back on SOLID nodes (Just for Flow_Without_Particle, because I have not defined collision between particles and solid objects.)
#if defined Flow_Without_Particle
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			for (int Z = start_Z; Z <= end_Z; ++Z) {
				if (is_solid[X][Y][Z] == FALSE) {

					if (is_solid[(X - 1 + Nx) % Nx][Y][Z] == TRUE) {
						pop[1][X][Y][Z] = pop[2][(X - 1 + Nx) % Nx][Y][Z];
					}
					if (is_solid[(X + 1) % Nx][Y][Z] == TRUE) {
						pop[2][X][Y][Z] = pop[1][(X + 1) % Nx][Y][Z];
					}
					if (is_solid[X][(Y - 1 + Ny) % Ny][Z] == TRUE) {
						pop[3][X][Y][Z] = pop[4][X][(Y - 1 + Ny) % Ny][Z];
					}
					if (is_solid[X][(Y + 1) % Ny][Z] == TRUE) {
						pop[4][X][Y][Z] = pop[3][X][(Y + 1) % Ny][Z];
					}
					if (is_solid[X][Y][(Z - 1 + Nz) % Nz] == TRUE) {
						pop[5][X][Y][Z] = pop[6][X][Y][(Z - 1 + Nz) % Nz];
					}
					if (is_solid[X][Y][(Z + 1) % Nz] == TRUE) {
						pop[6][X][Y][Z] = pop[5][X][Y][(Z + 1) % Nz];
					}
					if (is_solid[(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] == TRUE) {
						pop[7][X][Y][Z] = pop[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z];
					}
					if (is_solid[(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] == TRUE) {
						pop[8][X][Y][Z] = pop[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z];
					}
					if (is_solid[(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] == TRUE) {
						pop[9][X][Y][Z] = pop[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z];
					}
					if (is_solid[(X + 1) % Nx][(Y + 1) % Ny][Z] == TRUE) {
						pop[10][X][Y][Z] = pop[7][(X + 1) % Nx][(Y + 1) % Ny][Z];
					}
					if (is_solid[(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] == TRUE) {
						pop[11][X][Y][Z] = pop[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz];
					}
					if (is_solid[(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] == TRUE) {
						pop[12][X][Y][Z] = pop[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz];
					}
					if (is_solid[(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] == TRUE) {
						pop[13][X][Y][Z] = pop[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz];
					}
					if (is_solid[(X + 1) % Nx][Y][(Z + 1) % Nz] == TRUE) {
						pop[14][X][Y][Z] = pop[11][(X + 1) % Nx][Y][(Z + 1) % Nz];
					}
					if (is_solid[X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] == TRUE) {
						pop[15][X][Y][Z] = pop[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz];
					}
					if (is_solid[X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] == TRUE) {
						pop[16][X][Y][Z] = pop[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz];
					}
					if (is_solid[X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] == TRUE) {
						pop[17][X][Y][Z] = pop[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz];
					}
					if (is_solid[X][(Y + 1) % Ny][(Z + 1) % Nz] == TRUE) {
						pop[18][X][Y][Z] = pop[15][X][(Y + 1) % Ny][(Z + 1) % Nz];
					}
				}
			}
		}
	}
#endif


	/// Boundary-Conditions on walls
	// Due to the presence of the rigid walls at y = 0 and y = Ny - 1, the populations have to be bounced back.
	// Ladd's momentum correction term is included for moving walls (wall velocity parallel to x-axis).
	// Periodicity of the lattice in x-direction is taken into account via the %-operator if spesified.

	int Xtest;
	Xtest = 0; //This parameter is just for avoiding applying X-direc boundary condition multiple times. Therefore it just improves computational time and nothing else
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {

#if defined Z_Wall
			if (is_solid[X][Y][1] == FALSE) {
				/// back wall (z = 0)
				pop[5][X][Y][1] = pop[6][X][Y][0];
				pop[11][X][Y][1] = pop[14][(X - 1 + Nx) % Nx][Y][0] + 6 * weight[14] * density[X][Y][1] * wall_vel_back;
				pop[12][X][Y][1] = pop[13][(X + 1) % Nx][Y][0] - 6 * weight[13] * density[X][Y][1] * wall_vel_back;
				pop[16][X][Y][1] = pop[17][X][(Y + 1) % Ny][0];
				pop[15][X][Y][1] = pop[18][X][(Y - 1 + Ny) % Ny][0];
			}
			if (is_solid[X][Y][Nz - 2] == FALSE) {
				/// Front wall (z = Nz - 1)
				pop[6][X][Y][Nz - 2] = pop[5][X][Y][Nz - 1];
				pop[14][X][Y][Nz - 2] = pop[11][(X + 1) % Nx][Y][Nz - 1] - 6 * weight[11] * density[X][Y][Nz - 2] * wall_vel_front;
				pop[13][X][Y][Nz - 2] = pop[12][(X - 1 + Nx) % Nx][Y][Nz - 1] + 6 * weight[12] * density[X][Y][Nz - 2] * wall_vel_front;
				pop[17][X][Y][Nz - 2] = pop[16][X][(Y - 1 + Ny) % Ny][Nz - 1];
				pop[18][X][Y][Nz - 2] = pop[15][X][(Y + 1) % Ny][Nz - 1];
			}
#endif

#if defined X_velIn_ZeroGradOut || defined X_Wall || defined ZeroGradIn_ZeroGradOut
			if (Xtest == 0) {
				for (int Z = start_Z; Z <= end_Z; ++Z) {

#if defined X_velIn_ZeroGradOut

					if (is_solid[1][Y][Z] == FALSE) {
						/// left (x = 0)
						pop[1][1][Y][Z] = pop[2][0][Y][Z] + 6 * weight[2] * density[1][Y][Z] * inletvel;
						pop[11][1][Y][Z] = pop[14][0][Y][(Z - 1 + Nz) % Nz] + 6 * weight[14] * density[1][Y][Z] * inletvel;
						pop[7][1][Y][Z] = pop[10][0][(Y - 1 + Ny) % Ny][Z] + 6 * weight[10] * density[1][Y][Z] * inletvel;
						pop[13][1][Y][Z] = pop[12][0][Y][(Z + 1) % Nz] + 6 * weight[12] * density[1][Y][Z] * inletvel;
						pop[9][1][Y][Z] = pop[8][0][(Y + 1) % Ny][Z] + 6 * weight[8] * density[1][Y][Z] * inletvel;
					}

					if (is_solid[Nx - 2][Y][Z] == FALSE) {
						/// right (x = Nx - 1)
						pop[2][Nx - 2][Y][Z] = pop[2][Nx - 3][Y][Z];//zero gradient condition
						pop[14][Nx - 2][Y][Z] = pop[14][Nx - 3][Y][Z];//zero gradient condition
						pop[10][Nx - 2][Y][Z] = pop[10][Nx - 3][Y][Z];//zero gradient condition
						pop[12][Nx - 2][Y][Z] = pop[12][Nx - 3][Y][Z];//zero gradient condition
						pop[8][Nx - 2][Y][Z] = pop[8][Nx - 3][Y][Z];//zero gradient condition

																	//Second order
																	//pop[2][Nx - 2][Y][Z] = 2 * pop[2][Nx - 3][Y][Z] - pop[2][Nx - 4][Y][Z];//zero gradient condition
																	//pop[14][Nx - 2][Y][Z] = 2 * pop[14][Nx - 3][Y][Z] - pop[14][Nx - 4][Y][Z];//zero gradient condition
																	//pop[10][Nx - 2][Y][Z] = 2 * pop[10][Nx - 3][Y][Z] - pop[10][Nx - 4][Y][Z];//zero gradient condition
																	//pop[12][Nx - 2][Y][Z] = 2 * pop[12][Nx - 3][Y][Z] - pop[12][Nx - 4][Y][Z];//zero gradient condition
																	//pop[8][Nx - 2][Y][Z] = 2 * pop[8][Nx - 3][Y][Z] - pop[8][Nx - 4][Y][Z];//zero gradient condition
					}
#endif

#if defined ZeroGradIn_ZeroGradOut
					if (is_solid[1][Y][Z] == FALSE) {
						/// left (x = 0)
						pop[1][1][Y][Z] = pop[1][2][Y][Z];
						pop[11][1][Y][Z] = pop[11][2][Y][Z];
						pop[7][1][Y][Z] = pop[7][2][Y][Z];
						pop[13][1][Y][Z] = pop[13][2][Y][Z];
						pop[9][1][Y][Z] = pop[9][2][Y][Z];
						//pop[1][1][Y][Z] = 2 * pop[1][2][Y][Z] - pop[1][3][Y][Z];
						//pop[11][1][Y][Z] = 2 * pop[11][2][Y][Z] - pop[11][3][Y][Z];
						//pop[7][1][Y][Z] = 2 * pop[7][2][Y][Z] - pop[7][3][Y][Z];
						//pop[13][1][Y][Z] = 2 * pop[13][2][Y][Z] - pop[13][3][Y][Z];
						//pop[9][1][Y][Z] = 2 * pop[9][2][Y][Z] - pop[9][3][Y][Z];
					}
					if (is_solid[Nx - 2][Y][Z] == FALSE) {
						/// right (x = Nx - 1)
						pop[2][Nx - 2][Y][Z] = pop[2][Nx - 3][Y][Z];//zero gradient condition
						pop[14][Nx - 2][Y][Z] = pop[14][Nx - 3][Y][Z];//zero gradient condition
						pop[10][Nx - 2][Y][Z] = pop[10][Nx - 3][Y][Z];//zero gradient condition
						pop[12][Nx - 2][Y][Z] = pop[12][Nx - 3][Y][Z];//zero gradient condition
						pop[8][Nx - 2][Y][Z] = pop[8][Nx - 3][Y][Z];//zero gradient condition

																	//Second order
																	//pop[2][Nx - 2][Y][Z] = 2 * pop[2][Nx - 3][Y][Z] - pop[2][Nx - 4][Y][Z];//zero gradient condition
																	//pop[14][Nx - 2][Y][Z] = 2 * pop[14][Nx - 3][Y][Z] - pop[14][Nx - 4][Y][Z];//zero gradient condition
																	//pop[10][Nx - 2][Y][Z] = 2 * pop[10][Nx - 3][Y][Z] - pop[10][Nx - 4][Y][Z];//zero gradient condition
																	//pop[12][Nx - 2][Y][Z] = 2 * pop[12][Nx - 3][Y][Z] - pop[12][Nx - 4][Y][Z];//zero gradient condition
																	//pop[8][Nx - 2][Y][Z] = 2 * pop[8][Nx - 3][Y][Z] - pop[8][Nx - 4][Y][Z];//zero gradient condition
					}
#endif

#if defined X_Wall
					if (is_solid[1][Y][Z] == FALSE) {
						/// left (x = 0)
						pop[1][1][Y][Z] = pop[2][0][Y][Z];
						pop[11][1][Y][Z] = pop[14][0][Y][(Z - 1 + Nz) % Nz];
						pop[7][1][Y][Z] = pop[10][0][(Y - 1 + Ny) % Ny][Z];
						pop[13][1][Y][Z] = pop[12][0][Y][(Z + 1) % Nz];
						pop[9][1][Y][Z] = pop[8][0][(Y + 1) % Ny][Z];
					}
					if (is_solid[Nx - 2][Y][Z] == FALSE) {
						/// right (x = Nx - 1)
						pop[2][Nx - 2][Y][Z] = pop[1][Nx - 1][Y][Z];
						pop[14][Nx - 2][Y][Z] = pop[11][Nx - 1][Y][(Z + 1) % Nz];
						pop[10][Nx - 2][Y][Z] = pop[7][Nx - 1][(Y + 1) % Ny][Z];
						pop[12][Nx - 2][Y][Z] = pop[13][Nx - 1][Y][(Z - 1 + Nz) % Nz];
						pop[8][Nx - 2][Y][Z] = pop[9][Nx - 1][(Y - 1 + Ny) % Ny][Z];
					}
#endif
				}
			}
#endif
		}
		Xtest = 1;


#if defined Y_Wall
		for (int Z = start_Z; Z <= end_Z; ++Z) {
			if (is_solid[X][1][Z] == FALSE) {
				/// Bottom wall (y = 0)
				pop[3][X][1][Z] = pop[4][X][0][Z];
				pop[7][X][1][Z] = pop[10][(X - 1 + Nx) % Nx][0][Z] + 6 * weight[10] * density[X][1][Z] * wall_vel_bottom;
				pop[8][X][1][Z] = pop[9][(X + 1) % Nx][0][Z] - 6 * weight[9] * density[X][1][Z] * wall_vel_bottom;
				pop[17][X][1][Z] = pop[16][X][0][(Z + 1) % Nz];
				pop[15][X][1][Z] = pop[18][X][0][(Z - 1 + Nz) % Nz];
			}
			if (is_solid[X][Ny - 2][Z] == FALSE) {
				/// Top wall (y = Ny - 1)
				pop[4][X][Ny - 2][Z] = pop[3][X][Ny - 1][Z];
				pop[10][X][Ny - 2][Z] = pop[7][(X + 1) % Nx][Ny - 1][Z] - 6 * weight[7] * density[X][Ny - 2][Z] * wall_vel_top;
				pop[9][X][Ny - 2][Z] = pop[8][(X - 1 + Nx) % Nx][Ny - 1][Z] + 6 * weight[8] * density[X][Ny - 2][Z] * wall_vel_top;
				pop[16][X][Ny - 2][Z] = pop[17][X][Ny - 1][(Z - 1 + Nz) % Nz];
				pop[18][X][Ny - 2][Z] = pop[15][X][Ny - 1][(Z + 1) % Nz];
			}
		}
#endif
	}

	/// Compute fluid density and velocity
	// The fluid density and velocity are obtained from the populations.
	momenta();
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MOMENTA
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is not included here.
// It must be taken into account whenever the physical velocity is required.


void momenta() {


	int start_X, end_X;//Start and End of LBM calculation for processor.
	start_X = 0;
	end_X = Nx - 1;

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

	int start_Z, end_Z;
	start_Z = 0;
	end_Z = Nz - 1;

	//#if defined X_Wall
	//	start_X = 1;
	//	end_X = Nx - 2;
	//#endif
	//
	//#if defined X_velIn_ZeroGradOut
	//	start_X = 1;
	//	end_X = Nx - 2;
	//#endif
	//
	//#if defined ZeroGradIn_ZeroGradOut
	//	start_X = 1;
	//	end_X = Nx - 2;
	//#endif
	//
	//
	//#if defined Y_Wall
	//	start_Y = 1;
	//	end_Y = Ny - 2;
	//#endif
	//
	//
	//
	//#if defined Z_Wall
	//	start_Z = 1;
	//	end_Z = Nz - 2;
	//#endif


	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			for (int Z = start_Z; Z <= end_Z; ++Z) {
				//				if (is_solid[X][Y][Z] == FALSE){
				density[X][Y][Z] = pop[0][X][Y][Z] + pop[1][X][Y][Z] + pop[2][X][Y][Z] + pop[3][X][Y][Z] + pop[4][X][Y][Z] + pop[5][X][Y][Z] + pop[6][X][Y][Z] + pop[7][X][Y][Z] + pop[8][X][Y][Z] + pop[9][X][Y][Z] + pop[10][X][Y][Z] + pop[11][X][Y][Z] + pop[12][X][Y][Z] + pop[13][X][Y][Z] + pop[14][X][Y][Z] + pop[15][X][Y][Z] + pop[16][X][Y][Z] + pop[17][X][Y][Z] + pop[18][X][Y][Z];
				velocity_x[X][Y][Z] = (pop[1][X][Y][Z] - pop[2][X][Y][Z] + pop[7][X][Y][Z] - pop[8][X][Y][Z] + pop[9][X][Y][Z] - pop[10][X][Y][Z] + pop[11][X][Y][Z] - pop[12][X][Y][Z] + pop[13][X][Y][Z] - pop[14][X][Y][Z]) / density[X][Y][Z];
				velocity_y[X][Y][Z] = (pop[3][X][Y][Z] - pop[4][X][Y][Z] + pop[7][X][Y][Z] + pop[8][X][Y][Z] - pop[9][X][Y][Z] - pop[10][X][Y][Z] + pop[15][X][Y][Z] - pop[16][X][Y][Z] + pop[17][X][Y][Z] - pop[18][X][Y][Z]) / density[X][Y][Z];
				velocity_z[X][Y][Z] = (pop[5][X][Y][Z] - pop[6][X][Y][Z] + pop[11][X][Y][Z] + pop[12][X][Y][Z] - pop[13][X][Y][Z] - pop[14][X][Y][Z] + pop[15][X][Y][Z] + pop[16][X][Y][Z] - pop[17][X][Y][Z] - pop[18][X][Y][Z]) / density[X][Y][Z];
				force_x[X][Y][Z] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
				force_y[X][Y][Z] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
				force_z[X][Y][Z] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
										//				}
										//				else{
										//					velocity_x[X][Y][Z] = 0.0;
										//					velocity_y[X][Y][Z] = 0.0;
										//					velocity_z[X][Y][Z] = 0.0;

										//				}
			}
		}
	}

	//	for (int X = 0; X < Nx; ++X){
	//		for (int Y = 0; Y < Ny; ++Y){
	//			force_x[X][Y][0] = 0.0; //Back face
	//			force_x[X][Y][Nz - 1] = 0.0;//Front face
	//			force_y[X][Y][0] = 0.0;
	//			force_y[X][Y][Nz - 1] = 0.0;
	//			force_z[X][Y][0] = 0.0;
	//			force_z[X][Y][Nz - 1] = 0.0;
	//#if defined Z_Wall
	//			velocity_x[X][Y][0] = wall_vel_back;
	//			velocity_x[X][Y][Nz - 1] = wall_vel_front;
	//			velocity_y[X][Y][0] = 0.;
	//			velocity_y[X][Y][Nz - 1] = 0.;
	//			velocity_z[X][Y][0] = 0.;
	//			velocity_z[X][Y][Nz - 1] = 0.;
	//#endif
	//		}
	//	}
	//
	//	for (int Y = 0; Y < Ny; ++Y){
	//		for (int Z = 0; Z < Nz; ++Z){
	//			force_x[0][Y][Z] = 0.0;//Left wall
	//			force_x[Nx - 1][Y][Z] = 0.0;//Rgiht wall
	//			force_y[0][Y][Z] = 0.0;
	//			force_y[Nx - 1][Y][Z] = 0.0;
	//			force_z[0][Y][Z] = 0.0;
	//			force_z[Nx - 1][Y][Z] = 0.0;
	//#if defined X_Wall
	//			velocity_x[0][Y][Z] = 0.;
	//			velocity_x[Nx - 1][Y][Z] = 0.;
	//			velocity_y[0][Y][Z] = 0.;
	//			velocity_y[Nx - 1][Y][Z] = 0.;
	//			velocity_z[0][Y][Z] = 0.;
	//			velocity_z[Nx - 1][Y][Z] = 0.;
	//#endif
	//		}
	//	}
	//
	//	for (int X = 0; X < Nx; ++X){
	//		for (int Z = 0; Z < Nz; ++Z){
	//			force_x[X][0][Z] = 0.0;//Bottom wall
	//			force_x[X][Ny - 1][Z] = 0.0;//Top wall
	//			force_y[X][0][Z] = 0.0;
	//			force_y[X][Ny - 1][Z] = 0.0;
	//			force_z[X][0][Z] = 0.0;
	//			force_z[X][Ny - 1][Z] = 0.0;
	//#if defined Y_Wall
	//			velocity_x[X][0][Z] = wall_vel_bottom;
	//			velocity_x[X][Ny - 1][Z] = wall_vel_top;
	//			velocity_y[X][0][Z] = 0.;
	//			velocity_y[X][Ny - 1][Z] = 0.;
	//			velocity_z[X][0][Z] = 0.;
	//			velocity_z[X][Ny - 1][Z] = 0.;
	//#endif
	//		}
	//	}

	return;

}


//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM_T
//-----------------------------------------------------------------------------------------------------------------------------
//
#if defined Flow_With_Thermal_Effect
void LBM_t(int time) {

	/// Swap populations
	// The present code used old and new populations which are swapped at the beginning of each time step.
	// This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
	// This way, the old populations are not overwritten during propagation.
	// The resulting code is easier to write and to debug.
	// The memory requirement for the populations is twice as large.

	double ****swap_temp_t = pop_old_t;
	pop_old_t = pop_t;
	pop_t = swap_temp_t;

	/// Lattice Boltzmann equation
	// The lattice Boltzmann equation is solved in the following.
	// The algorithm includes
	// - computation of the lattice force
	// - combined collision and propagation (faster than first collision and then propagation)

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

	int start_Z, end_Z;
	start_Z = 0;
	end_Z = Nz - 1;

#if defined Left_Fixed_Temp || defined Left_Adiabatic
	start_X = 1;
#endif
#if defined Right_Fixed_Temp || defined Right_Adiabatic
	end_X = Nx - 2;
#endif
#if defined Bottom_Fixed_Temp || defined Bottom_Adiabatic
	start_Y = 1;
#endif
#if defined Top_Fixed_Temp || defined Top_Adiabatic
	end_Y = Ny - 2;
#endif
#if defined Back_Fixed_Temp || defined Back_Adiabatic
	start_Z = 1;
#endif
#if defined Front_Fixed_Temp || defined Front_Adiabatic
	end_Z = Nz - 2;
#endif

	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			for (int Z = start_Z; Z <= end_Z; ++Z) {

				/// Compute equilibrium
				// The equilibrium populations are computed.

				equilibrium_t(temperature[X][Y][Z], density[X][Y][Z], velocity_x[X][Y][Z], velocity_y[X][Y][Z], velocity_z[X][Y][Z]);

				/// Compute new populations
				// This is the lattice Boltzmann equation (combined collision and propagation).
				// Periodicity of the lattice in x-direction is taken into account by the %-operator.

#if defined Force_Thermal_1
				pop_t[0][X][Y][Z] = pop_old_t[0][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[0] * omega_t + weight[0] * force_thermal[X][Y][Z];
				pop_t[1][(X + 1) % Nx][Y][Z] = pop_old_t[1][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[1] * omega_t + weight[1] * force_thermal[X][Y][Z];
				pop_t[2][(X - 1 + Nx) % Nx][Y][Z] = pop_old_t[2][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[2] * omega_t + weight[2] * force_thermal[X][Y][Z];
				pop_t[3][X][(Y + 1) % Ny][Z] = pop_old_t[3][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[3] * omega_t + weight[3] * force_thermal[X][Y][Z];
				pop_t[4][X][(Y - 1 + Ny) % Ny][Z] = pop_old_t[4][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[4] * omega_t + weight[4] * force_thermal[X][Y][Z];
				pop_t[5][X][Y][(Z + 1) % Nz] = pop_old_t[5][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[5] * omega_t + weight[5] * force_thermal[X][Y][Z];
				pop_t[6][X][Y][(Z - 1 + Nz) % Nz] = pop_old_t[6][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[6] * omega_t + weight[6] * force_thermal[X][Y][Z];
				pop_t[11][(X + 1) % Nx][Y][(Z + 1) % Nz] = pop_old_t[11][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[11] * omega_t + weight[11] * force_thermal[X][Y][Z];
				pop_t[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old_t[14][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[14] * omega_t + weight[14] * force_thermal[X][Y][Z];
				pop_t[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old_t[10][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[10] * omega_t + weight[10] * force_thermal[X][Y][Z];
				pop_t[7][(X + 1) % Nx][(Y + 1) % Ny][Z] = pop_old_t[7][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[7] * omega_t + weight[7] * force_thermal[X][Y][Z];
				pop_t[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] = pop_old_t[12][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[12] * omega_t + weight[12] * force_thermal[X][Y][Z];
				pop_t[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old_t[13][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[13] * omega_t + weight[13] * force_thermal[X][Y][Z];
				pop_t[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old_t[9][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[9] * omega_t + weight[9] * force_thermal[X][Y][Z];
				pop_t[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] = pop_old_t[8][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[8] * omega_t + weight[8] * force_thermal[X][Y][Z];
				pop_t[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] = pop_old_t[16][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[16] * omega_t + weight[16] * force_thermal[X][Y][Z];
				pop_t[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] = pop_old_t[17][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[17] * omega_t + weight[17] * force_thermal[X][Y][Z];
				pop_t[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] = pop_old_t[18][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[18] * omega_t + weight[18] * force_thermal[X][Y][Z];
				pop_t[15][X][(Y + 1) % Ny][(Z + 1) % Nz] = pop_old_t[15][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[15] * omega_t + weight[15] * force_thermal[X][Y][Z];


#endif
#if defined Force_Thermal_2
				pop_t[0][X][Y][Z] = pop_old_t[0][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[0] * omega_t + (1 - 0.5*omega_t)*weight[0] * force_thermal[X][Y][Z];
				pop_t[1][(X + 1) % Nx][Y][Z] = pop_old_t[1][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[1] * omega_t + (1 - 0.5*omega_t)* weight[1] * force_thermal[X][Y][Z];
				pop_t[2][(X - 1 + Nx) % Nx][Y][Z] = pop_old_t[2][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[2] * omega_t + (1 - 0.5*omega_t)*weight[2] * force_thermal[X][Y][Z];
				pop_t[3][X][(Y + 1) % Ny][Z] = pop_old_t[3][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[3] * omega_t + (1 - 0.5*omega_t)*weight[3] * force_thermal[X][Y][Z];
				pop_t[4][X][(Y - 1 + Ny) % Ny][Z] = pop_old_t[4][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[4] * omega_t + (1 - 0.5*omega_t)*weight[4] * force_thermal[X][Y][Z];
				pop_t[5][X][Y][(Z + 1) % Nz] = pop_old_t[5][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[5] * omega_t + (1 - 0.5*omega_t)*weight[5] * force_thermal[X][Y][Z];
				pop_t[6][X][Y][(Z - 1 + Nz) % Nz] = pop_old_t[6][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[6] * omega_t + (1 - 0.5*omega_t)*weight[6] * force_thermal[X][Y][Z];
				pop_t[11][(X + 1) % Nx][Y][(Z + 1) % Nz] = pop_old_t[11][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[11] * omega_t + (1 - 0.5*omega_t)*weight[11] * force_thermal[X][Y][Z];
				pop_t[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old_t[14][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[14] * omega_t + (1 - 0.5*omega_t)* weight[14] * force_thermal[X][Y][Z];
				pop_t[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old_t[10][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[10] * omega_t + (1 - 0.5*omega_t)*weight[10] * force_thermal[X][Y][Z];
				pop_t[7][(X + 1) % Nx][(Y + 1) % Ny][Z] = pop_old_t[7][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[7] * omega_t + (1 - 0.5*omega_t)*weight[7] * force_thermal[X][Y][Z];
				pop_t[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz] = pop_old_t[12][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[12] * omega_t + (1 - 0.5*omega_t)*weight[12] * force_thermal[X][Y][Z];
				pop_t[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz] = pop_old_t[13][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[13] * omega_t + (1 - 0.5*omega_t)*weight[13] * force_thermal[X][Y][Z];
				pop_t[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z] = pop_old_t[9][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[9] * omega_t + (1 - 0.5*omega_t)*weight[9] * force_thermal[X][Y][Z];
				pop_t[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z] = pop_old_t[8][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[8] * omega_t + (1 - 0.5*omega_t)*weight[8] * force_thermal[X][Y][Z];
				pop_t[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz] = pop_old_t[16][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[16] * omega_t + (1 - 0.5*omega_t)*weight[16] * force_thermal[X][Y][Z];
				pop_t[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz] = pop_old_t[17][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[17] * omega_t + (1 - 0.5*omega_t)*weight[17] * force_thermal[X][Y][Z];
				pop_t[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz] = pop_old_t[18][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[18] * omega_t + (1 - 0.5*omega_t)*weight[18] * force_thermal[X][Y][Z];
				pop_t[15][X][(Y + 1) % Ny][(Z + 1) % Nz] = pop_old_t[15][X][Y][Z] * (1.0 - omega_t) + pop_eq_t[15] * omega_t + (1 - 0.5*omega_t)*weight[15] * force_thermal[X][Y][Z];
#endif

			}
		}
	}
	///////////////////////////
	/// Temperature boundary conditions on walls

	int Xtest;
	Xtest = 0; //This parameter is just for avoiding applying X-direc boundary condition multiple times. Therefore it just improves computational time and nothing else
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {

#if defined Back_Fixed_Temp
			pop_t[5][X][Y][1] = (weight[5] + weight[6])*temp_back - pop_t[6][X][Y][0];
			pop_t[11][X][Y][1] = (weight[11] + weight[14])*temp_back - pop_t[14][(X - 1 + Nx) % Nx][Y][0];
			pop_t[12][X][Y][1] = (weight[12] + weight[13])*temp_back - pop_t[13][(X + 1) % Nx][Y][0];
			pop_t[16][X][Y][1] = (weight[16] + weight[17])*temp_back - pop_t[17][X][(Y + 1) % Ny][0];
			pop_t[15][X][Y][1] = (weight[15] + weight[18])*temp_back - pop_t[18][X][(Y - 1 + Ny) % Ny][0];
#endif
#if defined Back_Adiabatic
			pop_t[5][X][Y][1] = pop_t[5][X][Y][2];
			pop_t[11][X][Y][1] = pop_t[11][X][Y][2];
			pop_t[12][X][Y][1] = pop_t[12][X][Y][2];
			pop_t[16][X][Y][1] = pop_t[16][X][Y][2];
			pop_t[15][X][Y][1] = pop_t[15][X][Y][2];

			//Second order
			//pop_t[5][X][Y][1] = 2 * pop_t[5][X][Y][2] - pop_t[5][X][Y][3];
			//pop_t[11][X][Y][1] = 2 * pop_t[11][X][Y][2] - pop_t[11][X][Y][3];
			//pop_t[12][X][Y][1] = 2 * pop_t[12][X][Y][2] - pop_t[12][X][Y][3];
			//pop_t[16][X][Y][1] = 2 * pop_t[16][X][Y][2] - pop_t[16][X][Y][3];
			//pop_t[15][X][Y][1] = 2 * pop_t[15][X][Y][2] - pop_t[15][X][Y][3];
#endif
#if defined Front_Fixed_Temp
			pop_t[6][X][Y][Nz - 2] = (weight[5] + weight[6])*temp_front - pop_t[5][X][Y][Nz - 1];
			pop_t[14][X][Y][Nz - 2] = (weight[11] + weight[14])*temp_front - pop_t[11][(X + 1) % Nx][Y][Nz - 1];
			pop_t[13][X][Y][Nz - 2] = (weight[12] + weight[13])*temp_front - pop_t[12][(X - 1 + Nx) % Nx][Y][Nz - 1];
			pop_t[17][X][Y][Nz - 2] = (weight[16] + weight[17])*temp_front - pop_t[16][X][(Y - 1 + Ny) % Ny][Nz - 1];
			pop_t[18][X][Y][Nz - 2] = (weight[15] + weight[18])*temp_front - pop_t[15][X][(Y + 1) % Ny][Nz - 1];
#endif
#if defined Front_Adiabatic
			pop_t[6][X][Y][Nz - 2] = pop_t[6][X][Y][Nz - 3];
			pop_t[14][X][Y][Nz - 2] = pop_t[14][X][Y][Nz - 3];
			pop_t[13][X][Y][Nz - 2] = pop_t[13][X][Y][Nz - 3];
			pop_t[17][X][Y][Nz - 2] = pop_t[17][X][Y][Nz - 3];
			pop_t[18][X][Y][Nz - 2] = pop_t[18][X][Y][Nz - 3];

			//Second order
			//pop_t[6][X][Y][Nz - 2] = 2 * pop_t[6][X][Y][Nz - 3] - pop_t[6][X][Y][Nz - 4];
			//pop_t[14][X][Y][Nz - 2] = 2 * pop_t[14][X][Y][Nz - 3] - pop_t[14][X][Y][Nz - 4];
			//pop_t[13][X][Y][Nz - 2] = 2 * pop_t[13][X][Y][Nz - 3] - pop_t[13][X][Y][Nz - 4];
			//pop_t[17][X][Y][Nz - 2] = 2 * pop_t[17][X][Y][Nz - 3] - pop_t[17][X][Y][Nz - 4];
			//pop_t[18][X][Y][Nz - 2] = 2 * pop_t[18][X][Y][Nz - 3] - pop_t[18][X][Y][Nz - 4];
#endif

#if defined Left_Fixed_Temp || defined Right_Fixed_Temp || defined Left_Adiabatic || defined Right_Adiabatic
			if (Xtest == 0) {
				for (int Z = start_Z; Z <= end_Z; ++Z) {

#if defined Left_Fixed_Temp
					/// left (x = 0)
					pop_t[1][1][Y][Z] = (weight[1] + weight[2])*temp_left - pop_t[2][0][Y][Z];
					pop_t[11][1][Y][Z] = (weight[11] + weight[14])*temp_left - pop_t[14][0][Y][(Z - 1 + Nz) % Nz];
					pop_t[7][1][Y][Z] = (weight[7] + weight[10])*temp_left - pop_t[10][0][(Y - 1 + Ny) % Ny][Z];
					pop_t[13][1][Y][Z] = (weight[13] + weight[12])*temp_left - pop_t[12][0][Y][(Z + 1) % Nz];
					pop_t[9][1][Y][Z] = (weight[9] + weight[8])*temp_left - pop_t[8][0][(Y + 1) % Ny][Z];
#endif
#if defined Left_Adiabatic
					pop_t[1][1][Y][Z] = pop_t[1][2][Y][Z];
					pop_t[11][1][Y][Z] = pop_t[11][2][Y][Z];
					pop_t[7][1][Y][Z] = pop_t[7][2][Y][Z];
					pop_t[13][1][Y][Z] = pop_t[13][2][Y][Z];
					pop_t[9][1][Y][Z] = pop_t[9][2][Y][Z];

					//Second order
					//pop_t[1][1][Y][Z] = 2 * pop_t[1][2][Y][Z] - pop_t[1][3][Y][Z];
					//pop_t[11][1][Y][Z] = 2 * pop_t[11][2][Y][Z] - pop_t[11][3][Y][Z];
					//pop_t[7][1][Y][Z] = 2 * pop_t[7][2][Y][Z] - pop_t[7][3][Y][Z];
					//pop_t[13][1][Y][Z] = 2 * pop_t[13][2][Y][Z] - pop_t[13][3][Y][Z];
					//pop_t[9][1][Y][Z] = 2 * pop_t[9][2][Y][Z] - pop_t[9][3][Y][Z];
#endif

#if defined Right_Fixed_Temp
					pop_t[2][Nx - 2][Y][Z] = (weight[1] + weight[2])*temp_right - pop_t[1][Nx - 1][Y][Z];
					pop_t[14][Nx - 2][Y][Z] = (weight[11] + weight[14])*temp_right - pop_t[11][Nx - 1][Y][(Z + 1) % Nz];
					pop_t[10][Nx - 2][Y][Z] = (weight[7] + weight[10])*temp_right - pop_t[7][Nx - 1][(Y + 1) % Ny][Z];
					pop_t[12][Nx - 2][Y][Z] = (weight[13] + weight[12])*temp_right - pop_t[13][Nx - 1][Y][(Z - 1 + Nz) % Nz];
					pop_t[8][Nx - 2][Y][Z] = (weight[9] + weight[8])*temp_right - pop_t[9][Nx - 1][(Y - 1 + Ny) % Ny][Z];
#endif
#if defined Right_Adiabatic
					//pop_t[2][Nx - 2][Y][Z] = pop_t[2][Nx - 3][Y][Z];
					//pop_t[14][Nx - 2][Y][Z] = pop_t[14][Nx - 3][Y][Z];
					//pop_t[10][Nx - 2][Y][Z] = pop_t[10][Nx - 3][Y][Z];
					//pop_t[12][Nx - 2][Y][Z] = pop_t[12][Nx - 3][Y][Z];
					//pop_t[8][Nx - 2][Y][Z] = pop_t[8][Nx - 3][Y][Z];

					//Second order
					pop_t[2][Nx - 2][Y][Z] = 2 * pop_t[2][Nx - 3][Y][Z] - pop_t[2][Nx - 4][Y][Z];
					pop_t[14][Nx - 2][Y][Z] = 2 * pop_t[14][Nx - 3][Y][Z] - pop_t[14][Nx - 4][Y][Z];
					pop_t[10][Nx - 2][Y][Z] = 2 * pop_t[10][Nx - 3][Y][Z] - pop_t[10][Nx - 4][Y][Z];
					pop_t[12][Nx - 2][Y][Z] = 2 * pop_t[12][Nx - 3][Y][Z] - pop_t[12][Nx - 4][Y][Z];
					pop_t[8][Nx - 2][Y][Z] = 2 * pop_t[8][Nx - 3][Y][Z] - pop_t[8][Nx - 4][Y][Z];
#endif
				}
			}
#endif
		}
		Xtest = 1;

		for (int Z = start_Z; Z <= end_Z; ++Z) {
#if defined Bottom_Fixed_Temp
			pop_t[3][X][1][Z] = (weight[3] + weight[4])*temp_bot - pop_t[4][X][0][Z];
			pop_t[7][X][1][Z] = (weight[7] + weight[10])*temp_bot - pop_t[10][(X - 1 + Nx) % Nx][0][Z];
			pop_t[8][X][1][Z] = (weight[8] + weight[9])*temp_bot - pop_t[9][(X + 1) % Nx][0][Z];
			pop_t[17][X][1][Z] = (weight[17] + weight[16])*temp_bot - pop_t[16][X][0][(Z + 1) % Nz];
			pop_t[15][X][1][Z] = (weight[15] + weight[18])*temp_bot - pop_t[18][X][0][(Z - 1 + Nz) % Nz];
#endif
#if defined Bottom_Adiabatic
			pop_t[3][X][1][Z] = pop_t[3][X][2][Z];
			pop_t[7][X][1][Z] = pop_t[7][X][2][Z];
			pop_t[8][X][1][Z] = pop_t[8][X][2][Z];
			pop_t[17][X][1][Z] = pop_t[17][X][2][Z];
			pop_t[15][X][1][Z] = pop_t[15][X][2][Z];

			//Second order
			//pop_t[3][X][1][Z] = 2 * pop_t[3][X][2][Z] - pop_t[3][X][3][Z];
			//pop_t[7][X][1][Z] = 2 * pop_t[7][X][2][Z] - pop_t[7][X][3][Z];
			//pop_t[8][X][1][Z] = 2 * pop_t[8][X][2][Z] - pop_t[8][X][3][Z];
			//pop_t[17][X][1][Z] = 2 * pop_t[17][X][2][Z] - pop_t[17][X][3][Z];
			//pop_t[15][X][1][Z] = 2 * pop_t[15][X][2][Z] - pop_t[15][X][3][Z];
#endif
#if defined Top_Fixed_Temp
			pop_t[4][X][Ny - 2][Z] = (weight[4] + weight[3])*temp_top - pop_t[3][X][Ny - 1][Z];
			pop_t[10][X][Ny - 2][Z] = (weight[10] + weight[7])*temp_top - pop_t[7][(X + 1) % Nx][Ny - 1][Z];
			pop_t[9][X][Ny - 2][Z] = (weight[9] + weight[8])*temp_top - pop_t[8][(X - 1 + Nx) % Nx][Ny - 1][Z];
			pop_t[16][X][Ny - 2][Z] = (weight[16] + weight[17])*temp_top - pop_t[17][X][Ny - 1][(Z - 1 + Nz) % Nz];
			pop_t[18][X][Ny - 2][Z] = (weight[18] + weight[15])*temp_top - pop_t[15][X][Ny - 1][(Z + 1) % Nz];
#endif
#if defined Top_Adiabatic
			pop_t[4][X][Ny - 2][Z] = pop_t[4][X][Ny - 3][Z];
			pop_t[10][X][Ny - 2][Z] = pop_t[10][X][Ny - 3][Z];
			pop_t[9][X][Ny - 2][Z] = pop_t[9][X][Ny - 3][Z];
			pop_t[16][X][Ny - 2][Z] = pop_t[16][X][Ny - 3][Z];
			pop_t[18][X][Ny - 2][Z] = pop_t[18][X][Ny - 3][Z];

			//Second order
			//pop_t[4][X][Ny - 2][Z] = 2 * pop_t[4][X][Ny - 3][Z] - pop_t[4][X][Ny - 4][Z];
			//pop_t[10][X][Ny - 2][Z] = 2 * pop_t[10][X][Ny - 3][Z] - pop_t[10][X][Ny - 4][Z];
			//pop_t[9][X][Ny - 2][Z] = 2 * pop_t[9][X][Ny - 3][Z] - pop_t[9][X][Ny - 4][Z];
			//pop_t[16][X][Ny - 2][Z] = 2 * pop_t[16][X][Ny - 3][Z] - pop_t[16][X][Ny - 4][Z];
			//pop_t[18][X][Ny - 2][Z] = 2 * pop_t[18][X][Ny - 3][Z] - pop_t[18][X][Ny - 4][Z];
#endif
		}
	}

	/// Compute fluid density and velocity
	// The fluid density and velocity are obtained from the populations.
	momenta_t();

	return;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MOMENTA_T
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the fluid temperature from the populations.
#if defined Flow_With_Thermal_Effect
void momenta_t() {

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

	int start_Z, end_Z;
	start_Z = 0;
	end_Z = Nz - 1;

	//#if defined Left_Fixed_Temp || defined Left_Adiabatic
	//	start_X = 1;
	//#endif
	//#if defined Right_Fixed_Temp || defined Right_Adiabatic
	//	end_X = Nx - 2;
	//#endif
	//#if defined Bottom_Fixed_Temp || defined Bottom_Adiabatic
	//	start_Y = 1;
	//#endif
	//#if defined Top_Fixed_Temp || defined Top_Adiabatic
	//	end_Y = Ny - 2;
	//#endif
	//#if defined Back_Fixed_Temp || defined Back_Adiabatic
	//	start_Z = 1;
	//#endif
	//#if defined Front_Fixed_Temp || defined Front_Adiabatic
	//	end_Z = Nz - 2;
	//#endif

	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			for (int Z = start_Z; Z <= end_Z; ++Z) {
#if defined Equilibrium_Thermal_1
				temperature[X][Y][Z] = pop_t[0][X][Y][Z] + pop_t[1][X][Y][Z] + pop_t[2][X][Y][Z] + pop_t[3][X][Y][Z] + pop_t[4][X][Y][Z] + pop_t[5][X][Y][Z] + pop_t[6][X][Y][Z] + pop_t[7][X][Y][Z] + pop_t[8][X][Y][Z] + pop_t[9][X][Y][Z] + pop_t[10][X][Y][Z] + pop_t[11][X][Y][Z] + pop_t[12][X][Y][Z] + pop_t[13][X][Y][Z] + pop_t[14][X][Y][Z] + pop_t[15][X][Y][Z] + pop_t[16][X][Y][Z] + pop_t[17][X][Y][Z] + pop_t[18][X][Y][Z];

#endif
#if defined Equilibrium_Thermal_2
				temperature[X][Y][Z] = (pop_t[0][X][Y][Z] + pop_t[1][X][Y][Z] + pop_t[2][X][Y][Z] + pop_t[3][X][Y][Z] + pop_t[4][X][Y][Z] + pop_t[5][X][Y][Z] + pop_t[6][X][Y][Z] + pop_t[7][X][Y][Z] + pop_t[8][X][Y][Z] + pop_t[9][X][Y][Z] + pop_t[10][X][Y][Z] + pop_t[11][X][Y][Z] + pop_t[12][X][Y][Z] + pop_t[13][X][Y][Z] + pop_t[14][X][Y][Z] + pop_t[15][X][Y][Z] + pop_t[16][X][Y][Z] + pop_t[17][X][Y][Z] + pop_t[18][X][Y][Z]) / density[X][Y][Z];
#endif
				force_thermal[X][Y][Z] = 0.0;
			}
		}
	}


	//	for (int X = 0; X < Nx; ++X){
	//		for (int Y = 0; Y < Ny; ++Y){
	//			force_thermal[X][Y][0] = 0.0; //Back face
	//			force_thermal[X][Y][Nz - 1] = 0.0;//Front face
	//#if defined Back_Fixed_Temp
	//			temperature[X][Y][0] = temp_back;
	//#endif
	//#if defined Front_Fixed_Temp
	//			temperature[X][Y][Nz - 1] = temp_front;
	//#endif
	//		}
	//	}
	//
	//	for (int Y = 0; Y < Ny; ++Y){
	//		for (int Z = 0; Z < Nz; ++Z){
	//			force_thermal[0][Y][Z] = 0.0;//Left wall
	//			force_thermal[Nx - 1][Y][Z] = 0.0;//Rgiht wall
	//
	//#if defined Left_Fixed_Temp
	//			temperature[0][Y][Z] = temp_left;
	//#endif
	//#if defined Right_Fixed_Temp
	//			temperature[Nx - 1][Y][Z] = temp_right;
	//#endif
	//		}
	//	}
	//
	//	for (int X = 0; X < Nx; ++X){
	//		for (int Z = 0; Z < Nz; ++Z){
	//			force_thermal[X][0][Z] = 0.0;//Bottom wall
	//			force_thermal[X][Ny - 1][Z] = 0.0;//Top wall
	//
	//#if defined Bottom_Fixed_Temp
	//			temperature[X][0][Z] = temp_bot;
	//#endif
	//#if defined Bottom_Fixed_Temp
	//			temperature[X][Ny - 1][Z] = temp_top;
	//#endif
	//		}
	//	}

	return;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MEQ
//-----------------------------------------------------------------------------------------------------------------------------
//calculation the equilibrium moment for MRT

double meq(double RHO, double U, double V, double W, int k)
{
	double x;
	switch (k)
	{
	case 0: {x = RHO; break; }
	case 1: {x = RHO*(-11 + 19 * (U*U + V*V + W*W)); break; }
	case 2: {x = RHO*(3 - 11. / 2.*(U*U + V*V + W*W)); break; }
	case 3: {x = RHO*U; break; }
	case 4: {x = RHO*(-2. / 3.)*U; break; }
	case 5: {x = RHO*V; break; }
	case 6: {x = RHO*(-2. / 3.)*V; break; }
	case 7: {x = RHO*W; break; }
	case 8: {x = RHO*(-2. / 3.)*W; break; }
	case 9: {x = RHO*(2 * U*U - V*V - W*W); break; }
	case 10: {x = -1. / 2.*RHO*(2.*U*U - V*V - W*W); break; }
	case 11: {x = RHO*(V*V - W*W); break; }
	case 12: {x = -1. / 2.*RHO*(V*V - W*W); break; }
	case 13: {x = RHO*(U*V); break; }
	case 14: {x = RHO*(V*W); break; }
	case 15: {x = RHO*(U*W); break; }
	case 16: {x = 0.0; break; }
	case 17: {x = 0.0; break; }
	case 18: {x = 0.0; break; }
	default: x = 0;
	}
	return x;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: UPDATE_VERLET
//-----------------------------------------------------------------------------------------------------------------------------
#if defined Flow_With_Particle
void update_verlet(particle_struct particle[num_particles]) {

	int ssu = 0;
	for (size_t i = 0; i < num_particles; i++) {
		verlet_list[i].resize(0);
		for (size_t j = 0; j < num_particles; j++) {

			if (j != i) {
				double Dij = sqrt(SQ(particle[i].center[0].x - particle[j].center[0].x) + SQ(particle[i].center[0].y - particle[j].center[0].y) + SQ(particle[i].center[0].z - particle[j].center[0].z));
				if (Dij<r_buffer)
				{
					verlet_list[i].push_back(j);
					ssu++;
				}
			}
		}
	}
	//cout << "total "<<ssu << endl;
	return;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_SPHERE
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for spherical objects
#if defined Flow_With_Particle
#if defined MOVING_SPHERE

void collision_sphere(particle_struct particle[num_particles], int np) {


#if defined Sphere_Collision_1_Feng
	//Collision Forces. JUST for Sphere shape
	double Dij;
	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dij = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y) + SQ(particle[np].center[0].z - particle[np2].center[0].z));

			if (Dij <= (2.0 * particle[np].radius + th) && Dij>(2.0 * particle[np].radius)) {
				particle[np].center[0].f3x += cij / ep * ((particle[np].center[0].x - particle[np2].center[0].x)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
				particle[np].center[0].f3y += cij / ep * ((particle[np].center[0].y - particle[np2].center[0].y)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
				particle[np].center[0].f3z += cij / ep * ((particle[np].center[0].z - particle[np2].center[0].z)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
			}
			if (Dij <= (2.0 * particle[np].radius)) {
				particle[np].center[0].f3x += (particle[np].center[0].x - particle[np2].center[0].x) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
				particle[np].center[0].f3y += (particle[np].center[0].y - particle[np2].center[0].y) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
				particle[np].center[0].f3z += (particle[np].center[0].z - particle[np2].center[0].z) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_Wall
	if (particle[np].center[0].y <= (Ny - 2) / 2) {//Particle-Y-Wall Collision section
		Rij = particle[np].center[0].y + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Ny - 1) - particle[np].center[0].y + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}
#endif
#if defined X_Wall
	if (particle[np].center[0].x <= (Nx - 2) / 2) {//Particle-Wall Collision section
		Rij = particle[np].center[0].x + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Nx - 1) - particle[np].center[0].x + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}
#endif
#if defined Z_Wall
	if (particle[np].center[0].z <= (Nz - 2) / 2) {//Particle-Wall Collision section
		Rij = particle[np].center[0].z + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Nz - 1) - particle[np].center[0].z + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
	}
#endif
	// END of Particle-Wall Collision
#endif //end of Sphere_Collision_1_Feng

#if defined Sphere_Collision_2_Feng_Verlet
	//Collision Forces. JUST for Sphere shape
	double Dij;
	for (int np2 = 0; np2 < verlet_list[np].size(); ++np2) {//Particle-Particle Collision
		int index_p = verlet_list[np][np2];
		if (index_p != np) {
			Dij = sqrt(SQ(particle[np].center[0].x - particle[index_p].center[0].x) + SQ(particle[np].center[0].y - particle[index_p].center[0].y) + SQ(particle[np].center[0].z - particle[index_p].center[0].z));

			if (Dij <= (2.0 * particle[np].radius + th) && Dij>(2.0 * particle[np].radius)) {
				particle[np].center[0].f3x += cij / ep * ((particle[np].center[0].x - particle[index_p].center[0].x)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
				particle[np].center[0].f3y += cij / ep * ((particle[np].center[0].y - particle[index_p].center[0].y)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
				particle[np].center[0].f3z += cij / ep * ((particle[np].center[0].z - particle[index_p].center[0].z)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
			}
			if (Dij <= (2.0 * particle[np].radius)) {
				particle[np].center[0].f3x += (particle[np].center[0].x - particle[index_p].center[0].x) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
				particle[np].center[0].f3y += (particle[np].center[0].y - particle[index_p].center[0].y) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
				particle[np].center[0].f3z += (particle[np].center[0].z - particle[index_p].center[0].z) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_Wall
	if (particle[np].center[0].y <= (Ny - 2) / 2) {//Particle-Y-Wall Collision section
		Rij = particle[np].center[0].y + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Ny - 1) - particle[np].center[0].y + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}
#endif
#if defined X_Wall
	if (particle[np].center[0].x <= (Nx - 2) / 2) {//Particle-Wall Collision section
		Rij = particle[np].center[0].x + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Nx - 1) - particle[np].center[0].x + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}
#endif
#if defined Z_Wall
	if (particle[np].center[0].z <= (Nz - 2) / 2) {//Particle-Wall Collision section
		Rij = particle[np].center[0].z + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Nz - 1) - particle[np].center[0].z + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - th) / th, 2) + cij / Epw / th*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
	}
#endif
	// END of Particle-Wall Collision
#endif //end of Sphere_Collision_2_Feng_Verlet
	//////////////////////////////////////////////////////

#if defined Sphere_Collision_3_Feng_Lub_Verlet
	//Collision Forces. JUST for Sphere shape
	double Dij;
	for (int np2 = 0; np2 < verlet_list[np].size(); ++np2) {//Particle-Particle Collision
		int index_p = verlet_list[np][np2];
		if (index_p != np) {
			Dij = sqrt(SQ(particle[np].center[0].x - particle[index_p].center[0].x) + SQ(particle[np].center[0].y - particle[index_p].center[0].y) + SQ(particle[np].center[0].z - particle[index_p].center[0].z));

			if (Dij <= (2.0 * particle[np].radius + th_lub) && Dij>(2.0 * particle[np].radius + th)) {
				double dot_pro;
				dot_pro = (particle[np].center[0].vel_x - particle[index_p].center[0].vel_x)*(particle[np].center[0].x - particle[index_p].center[0].x);
				particle[np].center[0].f3x += -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius*particle_radius / (particle_radius + particle_radius))  * (1. / (Dij - 2 * particle_radius) - 1. / th_lub) * dot_pro / Dij;
				dot_pro = (particle[np].center[0].vel_y - particle[index_p].center[0].vel_y)*(particle[np].center[0].y - particle[index_p].center[0].y);
				particle[np].center[0].f3y += -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius*particle_radius / (particle_radius + particle_radius))  * (1. / (Dij - 2 * particle_radius) - 1. / th_lub) * dot_pro / Dij;
				dot_pro = (particle[np].center[0].vel_z - particle[index_p].center[0].vel_z)*(particle[np].center[0].z - particle[index_p].center[0].z);
				particle[np].center[0].f3z += -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius*particle_radius / (particle_radius + particle_radius))  * (1. / (Dij - 2 * particle_radius) - 1. / th_lub) * dot_pro / Dij;

			}

			if (Dij <= (2.0 * particle[np].radius + th) && Dij>(2.0 * particle[np].radius)) {
				particle[np].center[0].f3x += cij / ep * ((particle[np].center[0].x - particle[index_p].center[0].x)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
				particle[np].center[0].f3y += cij / ep * ((particle[np].center[0].y - particle[index_p].center[0].y)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
				particle[np].center[0].f3z += cij / ep * ((particle[np].center[0].z - particle[index_p].center[0].z)) / Dij * pow((Dij - 2.0 * particle[np].radius - th) / th, 2);
			}
			if (Dij <= (2.0 * particle[np].radius)) {
				particle[np].center[0].f3x += (particle[np].center[0].x - particle[index_p].center[0].x) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
				particle[np].center[0].f3y += (particle[np].center[0].y - particle[index_p].center[0].y) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
				particle[np].center[0].f3z += (particle[np].center[0].z - particle[index_p].center[0].z) / Dij * (cij / ep * pow((Dij - 2.0 * particle[np].radius - th) / th, 2) + cij / Ep / th*(2.0*particle[np].radius - Dij));
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_Wall
	if (particle[np].center[0].y <= (Ny - 2) / 2) {//Particle-Y-Wall Collision section
		Rij = particle[np].center[0].y + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Ny - 1) - particle[np].center[0].y + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw_lub) && abs(Rij)>(2.0 * particle[np].radius + thw)) {
		double dot_pro = (particle[np].center[0].vel_y - 0)*Rij;
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius*particle_radius / (particle_radius + particle_radius))  * (1. / (abs(Rij) - 2 * particle_radius) - 1. / thw_lub) * dot_pro / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}


	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4z += 0.0;
	}
#endif

#if defined X_Wall
	if (particle[np].center[0].x <= (Nx - 2) / 2) {//Particle-Wall Collision section
		Rij = particle[np].center[0].x + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Nx - 1) - particle[np].center[0].x + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw_lub) && abs(Rij)>(2.0 * particle[np].radius + thw)) {
		double dot_pro = (particle[np].center[0].vel_x - 0)*Rij;
		particle[np].center[0].f4x += -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius*particle_radius / (particle_radius + particle_radius))  * (1. / (abs(Rij) - 2 * particle_radius) - 1. / thw_lub) * dot_pro / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += 0.0;
	}
#endif
#if defined Z_Wall
	if (particle[np].center[0].z <= (Nz - 2) / 2) {//Particle-Wall Collision section
		Rij = particle[np].center[0].z + particle[np].radius;
	}
	else
	{
		Rij = -1.0 * ((Nz - 1) - particle[np].center[0].z + particle[np].radius);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw_lub) && abs(Rij)>(2.0 * particle[np].radius + thw)) {
		double dot_pro = (particle[np].center[0].vel_z - 0)*Rij;
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius*particle_radius / (particle_radius + particle_radius))  * (1. / (abs(Rij) - 2 * particle_radius) - 1. / thw_lub) * dot_pro / abs(Rij);
	}

	if (abs(Rij) <= (2.0 * particle[np].radius + thw) && abs(Rij) > (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += cij / epw * pow((2.0 * particle[np].radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
	}
	if (abs(Rij) <= (2.0 * particle[np].radius)) {
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += 0.0;
		particle[np].center[0].f4z += (cij / epw *pow((abs(Rij) - 2.0 * particle[np].radius - th) / th, 2) + cij / Epw / th*(2.0*particle[np].radius - abs(Rij)))*Rij / abs(Rij);
	}
#endif
	// END of Particle-Wall Collision
#endif //end of Sphere_Collision_3_Feng_Lub_Verlet

	return;
}
#endif
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_SPHEROID
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for Spheroidal objects
#if defined Flow_With_Particle
#if defined MOVING_SPHEROID

void collision_spheroid(particle_struct particle[num_particles], int np, double xp[], double yp[], double zp[]) {

	//#if defined Spheroid_Collision_1
	//Collision Forces. for Spheroid shape
	double Dij, Dist, ffx, ffy, ffz, ffx4, ffy4, ffz4, norvec;
	double xnor_init, ynor_init, znor_init, xnor_slope, ynor_slope, znor_slope;
	double x_out, y_out, z_out;
	double next_point_x, next_point_y, next_point_z;
	double m_slope, n_slope, p_slope;
	double x_intersect1, y_intersect1, z_intersect1;
	double x_intersect2, y_intersect2, z_intersect2;
	double x_intersect_real, y_intersect_real, z_intersect_real;
	double tt_part1, tt_part2, tt_part3;
	double tt_intersect_1, tt_intersect_2;
	double x_0, y_0, z_0, slope_th, xp1, yp1, zp1;
	double a, b, radius_check;
	double x_0_inertial, y_0_inertial, z_0_inertial;
	double xnor_slope_old, ynor_slope_old, znor_slope_old;
	double Rij;
	double dist_to_wall;
	double M_curvature, N_curvature, R_curvature, fi_curvature;
	double x_intersect_real_inertial, y_intersect_real_inertial, z_intersect_real_inertial;
	double check_Y;
	double u_approach_x, u_approach_y, u_approach_z;
	double x_center_curve, y_center_curve, z_center_curve;
	double x_center_curve_inertial, y_center_curve_inertial, z_center_curve_inertial;
	double nor_slope_length;
	int mmm;

	slope_th = 0.001;//thereshold of slope. If less than this value, computation stops

	a = particle_radius;
	b = particle_radius_2;
	radius_check = MIN(SQ(a) / b, SQ(b) / a);
	std::vector<int> check_walls;

	//////// Particle-Wall collision  //////////////////
	//----------------------------------------------------//
#if defined X_Wall || defined Y_Wall || defined Z_Wall

	double thw_control;
#if defined Spheroid_Collision_1_Feng || defined Spheroid_Collision_2_Feng_Verlet
	thw_control = thw;
#endif
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
	thw_control = thw_lub;
#endif

#if defined Y_Wall
	if ((particle[np].center[0].y <= MAX(particle_radius, particle_radius_2) + thw_control)) {
		check_walls.push_back(1);//1 means lower Y-wall
	}
	if (particle[np].center[0].y >= Ny - 1 - (MAX(particle_radius, particle_radius_2) + thw_control)) {
		check_walls.push_back(2);//2 means upper Y-wall
	}
#endif
#if defined X_Wall
	if ((particle[np].center[0].x <= MAX(particle_radius, particle_radius_2) + thw_control)) {
		check_walls.push_back(3);//3 means lower X-wall
	}
	if (particle[np].center[0].x >= Nx - 1 - (MAX(particle_radius, particle_radius_2) + thw_control)) {
		check_walls.push_back(4);//4 means upper X-wall
	}
#endif
#if defined Z_Wall
	if ((particle[np].center[0].z <= MAX(particle_radius, particle_radius_2) + thw_control)) {
		check_walls.push_back(5);//5 means lower Z-wall
	}
	if (particle[np].center[0].z >= Nz - 1 - (MAX(particle_radius, particle_radius_2) + thw_control)) {
		check_walls.push_back(6);//6 means upper Z-wall
	}
#endif
	//////////////////////// Finding nearest point to wall and calculate the collision force  /////////////////////
	for (size_t i_ch = 0; i_ch < check_walls.size(); i_ch++) {
		dist_to_wall = 1000;

		switch (check_walls[i_ch])
		{
		case 1: {//Lower Y-wall
			xp1 = 0.0;//wall point in inertial coordinate minus particle center position
			yp1 = 0 - particle[np].center[0].y;
			zp1 = 0.0;

			xnor_slope = 0.0;
			ynor_slope = 1.0;
			znor_slope = 0.0;
			break;
		}
		case 2: {//Upper Y-wall
			xp1 = 0.0;//wall point in inertial coordinate minus particle center position
			yp1 = Ny - 1 - particle[np].center[0].y;
			zp1 = 0.0;

			xnor_slope = 0.0;
			ynor_slope = 1.0;
			znor_slope = 0.0;
			break;
		}
		case 3: {//Lower X-wall
			xp1 = 0 - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
			yp1 = 0.0;
			zp1 = 0.0;

			xnor_slope = 1.0;
			ynor_slope = 0.0;
			znor_slope = 0.0;
			break;
		}
		case 4: {//Upper X-wall
			xp1 = Nx - 1 - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
			yp1 = 0.0;
			zp1 = 0.0;

			xnor_slope = 1.0;
			ynor_slope = 0.0;
			znor_slope = 0.0;
			break;
		}
		case 5: {//Lower Z-wall
			xp1 = 0.0;//wall point in inertial coordinate minus particle center position
			yp1 = 0.0;
			zp1 = 0 - particle[np].center[0].z;

			xnor_slope = 0.0;
			ynor_slope = 0.0;
			znor_slope = 1.0;
			break;
		}
		case 6: {//Lower Z-wall
			xp1 = 0.0;//wall point in inertial coordinate minus particle center position
			yp1 = 0.0;
			zp1 = Nz - 1 - particle[np].center[0].z;

			xnor_slope = 0.0;
			ynor_slope = 0.0;
			znor_slope = 1.0;
			break;
		}
		default: std::cout << "Kein Antwort";
		}

		x_0 = 0;//means particle center position in body-fixed frame
		y_0 = 0;//means particle center position in body-fixed frame
		z_0 = 0;//means particle center position in body-fixed frame

		mmm = 0;

		do {
			x_out = particle[np].center[0].quater[1] * (xp1)+particle[np].center[0].quater[2] * (yp1)+particle[np].center[0].quater[3] * (zp1);//location of wall point in body-fixed frame
			y_out = particle[np].center[0].quater[4] * (xp1)+particle[np].center[0].quater[5] * (yp1)+particle[np].center[0].quater[6] * (zp1);
			z_out = particle[np].center[0].quater[7] * (xp1)+particle[np].center[0].quater[8] * (yp1)+particle[np].center[0].quater[9] * (zp1);
			m_slope = x_out - x_0;//slope of the line conecting the wall and the center of imaginary circle
			n_slope = y_out - y_0;
			p_slope = z_out - z_0;

			/// Intersection of ellipsoid and line, Obtained by MATLAB  ///
			tt_part1 = a*sqrt(SQ(a) * SQ(b) * SQ(n_slope) + SQ(a) *SQ(b)  * SQ(p_slope) - SQ(a) * SQ(n_slope) * SQ(z_0) + 2 * SQ(a) * n_slope*p_slope*y_0*z_0 - SQ(a * p_slope * y_0) + SQ(SQ(b))  * SQ(m_slope) - SQ(b* m_slope* y_0) - SQ(b * m_slope  * z_0) + 2 * SQ(b) * m_slope*n_slope*x_0*y_0 + 2 * SQ(b) * m_slope*p_slope*x_0*z_0 - SQ(b * n_slope  * x_0) - SQ(b * p_slope * x_0));
			tt_part2 = SQ(b) * m_slope*x_0 + SQ(a) * n_slope*y_0 + SQ(a) * p_slope*z_0;
			tt_part3 = (SQ(a  * n_slope) + SQ(a * p_slope) + SQ(b * m_slope));
			tt_intersect_1 = -(tt_part1 + tt_part2) / tt_part3;//Answers of equation
			tt_intersect_2 = -(-tt_part1 + tt_part2) / tt_part3;

			x_intersect1 = x_0 + m_slope*tt_intersect_1;//in body fixed frame
			y_intersect1 = y_0 + n_slope*tt_intersect_1;
			z_intersect1 = z_0 + p_slope*tt_intersect_1;

			x_intersect2 = x_0 + m_slope*tt_intersect_2;
			y_intersect2 = y_0 + n_slope*tt_intersect_2;
			z_intersect2 = z_0 + p_slope*tt_intersect_2;

			/// Find the nearst point
			if ((SQ(x_intersect1 - x_out) + SQ(y_intersect1 - y_out) + SQ(z_intersect1 - z_out)) < (SQ(x_intersect2 - x_out) + SQ(y_intersect2 - y_out) + SQ(z_intersect2 - z_out)))
			{
				x_intersect_real = x_intersect1;//in body fixed frame
				y_intersect_real = y_intersect1;//in body fixed frame
				z_intersect_real = z_intersect1;//in body fixed frame
			}
			else
			{
				x_intersect_real = x_intersect2;//in body fixed frame
				y_intersect_real = y_intersect2;//in body fixed frame
				z_intersect_real = z_intersect2;//in body fixed frame
			}
			xnor_slope_old = xnor_slope;
			ynor_slope_old = ynor_slope;
			znor_slope_old = znor_slope;

			xnor_slope = -2.*x_intersect_real / particle_radius;//slope of normal vector on spheroid in body_fixed frame. INWARD direction
			ynor_slope = -2.*y_intersect_real / particle_radius_2;
			znor_slope = -2.*z_intersect_real / particle_radius_2;

			nor_slope_length = (sqrt(SQ(xnor_slope) + SQ(ynor_slope) + SQ(znor_slope)));
			xnor_slope /= nor_slope_length;
			ynor_slope /= nor_slope_length;
			znor_slope /= nor_slope_length;

			x_0 = x_intersect_real + xnor_slope * radius_check;//center of an imaginary circle for next check. (body-fixed frame)
			y_0 = y_intersect_real + ynor_slope * radius_check;
			z_0 = z_intersect_real + znor_slope * radius_check;

			//finding the center of imaginary circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
			update_rot_spheroid(particle[np].center[0].deter, x_0_inertial, y_0_inertial, z_0_inertial, particle[np].center[0].quater, x_0, y_0, z_0, particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Find x_0 position in inertial frame

			switch (check_walls[i_ch])
			{
			case 1: {
				xp1 = x_0_inertial - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
				yp1 = 0 - particle[np].center[0].y;
				zp1 = z_0_inertial - particle[np].center[0].z;
				break;
			}
			case 2: {
				xp1 = x_0_inertial - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
				yp1 = Ny - 1 - particle[np].center[0].y;
				zp1 = z_0_inertial - particle[np].center[0].z;
				break;
			}
			case 3: {
				xp1 = 0 - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
				yp1 = y_0_inertial - particle[np].center[0].y;
				zp1 = z_0_inertial - particle[np].center[0].z;
				break;
			}
			case 4: {
				xp1 = Nx - 1 - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
				yp1 = y_0_inertial - particle[np].center[0].y;
				zp1 = z_0_inertial - particle[np].center[0].z;
				break;
			}
			case 5: {
				xp1 = x_0_inertial - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
				yp1 = y_0_inertial - particle[np].center[0].y;
				zp1 = 0 - particle[np].center[0].z;
				break;
			}
			case 6: {
				xp1 = x_0_inertial - particle[np].center[0].x;//wall point in inertial coordinate minus particle center position
				yp1 = y_0_inertial - particle[np].center[0].y;
				zp1 = Nz - 1 - particle[np].center[0].z;
				break;
			}
			default: cout << "Keine Antwort";
			}

			mmm++;
		} while (abs(xnor_slope - xnor_slope_old) >= slope_th || abs(ynor_slope - ynor_slope_old) >= slope_th || abs(znor_slope - znor_slope_old) >= slope_th);

		update_rot_spheroid(particle[np].center[0].deter, x_intersect_real_inertial, y_intersect_real_inertial, z_intersect_real_inertial, particle[np].center[0].quater, x_intersect_real, y_intersect_real, z_intersect_real, particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Find x_0 position in inertial frame

		switch (check_walls[i_ch])
		{
		case 1: {
			dist_to_wall = y_intersect_real_inertial - 0;
			break;
		}
		case 2: {
			dist_to_wall = Ny - 1 - y_intersect_real_inertial;
			break;
		}
		case 3: {
			dist_to_wall = x_intersect_real_inertial - 0;
			break;
		}
		case 4: {
			dist_to_wall = Nx - 1 - x_intersect_real_inertial;
			break;
		}
		case 5: {
			dist_to_wall = z_intersect_real_inertial - 0;
			break;
		}
		case 6: {
			dist_to_wall = Nz - 1 - z_intersect_real_inertial;
			break;
		}
		default: cout << "Keine Antwort";
		}

		if (dist_to_wall<thw_control) {
			//cout << time << endl;
			fi_curvature = asin(x_intersect_real / (sqrt(SQ(x_intersect_real) + SQ(y_intersect_real) + SQ(z_intersect_real))));
			//cout << fi_curvature * 180.0 / M_PI<<endl;
			R_curvature = SQ(a)*b / (SQ(a*sin(fi_curvature)) + SQ(b*cos(fi_curvature))); //Paper 1366

			x_center_curve = x_intersect_real + xnor_slope * R_curvature;//center of circle of curvature. (body-fixed frame)
			y_center_curve = y_intersect_real + ynor_slope * R_curvature;
			z_center_curve = z_intersect_real + znor_slope * R_curvature;

			//finding the center of curvature circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
			update_rot_spheroid(particle[np].center[0].deter, x_center_curve_inertial, y_center_curve_inertial, z_center_curve_inertial, particle[np].center[0].quater, x_center_curve, y_center_curve, z_center_curve, particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Find  position of circle of curvature in inertial frame

			switch (check_walls[i_ch])
			{
			case 1: {
				Rij = y_center_curve_inertial + R_curvature;
				break;
			}
			case 2: {
				Rij = -1.0 * ((Ny - 1) - y_center_curve_inertial + R_curvature);
				break;
			}
			case 3: {
				Rij = x_center_curve_inertial + R_curvature;

				break;
			}
			case 4: {
				//double ooo = -1.0 * ((Nx - 1) - x_intersect_real_inertial + 2 * R_curvature);
				Rij = -1.0 * ((Nx - 1) - x_center_curve_inertial + R_curvature);
				break;
			}
			case 5: {
				Rij = z_center_curve_inertial + R_curvature;

				break;
			}
			case 6: {
				Rij = -1.0 * ((Nz - 1) - z_center_curve_inertial + R_curvature);

				break;
			}
			default: std::cout << "Keine Antwort";
			}

			////////////// Force and torque calculation
			switch (check_walls[i_ch])
			{
			case 1: { // Lower Y-Wall
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (abs(Rij) <= (2.0 * R_curvature + thw_lub) && abs(Rij)>(2.0 * R_curvature + thw)) {
					u_approach_y = particle[np].center[0].vel_y - particle[np].center[0].omgZ * (x_center_curve_inertial - particle[np].center[0].x) + particle[np].center[0].omgX * (z_center_curve_inertial - particle[np].center[0].z);
					double dot_pro = (u_approach_y - 0)*Rij;
					ffx4 = 0.0;
					ffy4 = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature*R_curvature / (R_curvature + R_curvature))  * (1. / (abs(Rij) - 2 * R_curvature) - 1. / thw_lub) * dot_pro *Rij / SQ(Rij);
					ffz4 = 0.0;
				}
#endif
				if (abs(Rij) <= (2.0 * R_curvature + thw) && abs(Rij) > (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = cij / epw * pow((2.0 * R_curvature - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
					ffz4 = 0.0;
				}
				if (abs(Rij) <= (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = (cij / epw *pow((abs(Rij) - 2.0 * R_curvature - thw) / thw, 2) + cij / Epw / thw*(2.0*R_curvature - abs(Rij)))*Rij / abs(Rij);
					ffz4 = 0.0;
				}
				break;
			}
			case 2: {// Upper Y-Wall
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (abs(Rij) <= (2.0 * R_curvature + thw_lub) && abs(Rij)>(2.0 * R_curvature + thw)) {
					u_approach_y = particle[np].center[0].vel_y - particle[np].center[0].omgZ * (x_center_curve_inertial - particle[np].center[0].x) + particle[np].center[0].omgX * (z_center_curve_inertial - particle[np].center[0].z);
					double dot_pro = (u_approach_y - 0)*Rij;
					ffx4 = 0.0;
					ffy4 = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature*R_curvature / (R_curvature + R_curvature))  * (1. / (abs(Rij) - 2 * R_curvature) - 1. / thw_lub) * dot_pro *Rij / SQ(Rij);
					ffz4 = 0.0;
				}
#endif
				if (abs(Rij) <= (2.0 * R_curvature + thw) && abs(Rij) > (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = cij / epw * pow((2.0 * R_curvature - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
					ffz4 = 0.0;
				}
				if (abs(Rij) <= (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = (cij / epw *pow((abs(Rij) - 2.0 * R_curvature - thw) / thw, 2) + cij / Epw / thw*(2.0*R_curvature - abs(Rij)))*Rij / abs(Rij);
					ffz4 = 0.0;
				}
				break;
			}
			case 3: {// Lower X-Wall
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (abs(Rij) <= (2.0 * R_curvature + thw_lub) && abs(Rij)>(2.0 * R_curvature + thw)) {
					u_approach_x = particle[np].center[0].vel_x + particle[np].center[0].omgZ * (y_center_curve_inertial - particle[np].center[0].y) - particle[np].center[0].omgY * (z_center_curve_inertial - particle[np].center[0].z);
					double dot_pro = (u_approach_x - 0)*Rij;
					ffx4 = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature*R_curvature / (R_curvature + R_curvature))  * (1. / (abs(Rij) - 2 * R_curvature) - 1. / thw_lub) * dot_pro *Rij / SQ(Rij);
					ffy4 = 0.0;
					ffz4 = 0.0;
				}
#endif
				if (abs(Rij) <= (2.0 * R_curvature + thw) && abs(Rij) > (2.0 * R_curvature)) {
					ffx4 = cij / epw * pow((2.0 * R_curvature - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
					ffy4 = 0.0;
					ffz4 = 0.0;
				}
				if (abs(Rij) <= (2.0 * R_curvature)) {
					ffx4 = (cij / epw *pow((abs(Rij) - 2.0 * R_curvature - thw) / thw, 2) + cij / Epw / thw*(2.0*R_curvature - abs(Rij)))*Rij / abs(Rij);
					ffy4 = 0.0;
					ffz4 = 0.0;
				}
				break;
			}
			case 4: {// Upper X-Wall
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (abs(Rij) <= (2.0 * R_curvature + thw_lub) && abs(Rij)>(2.0 * R_curvature + thw)) {
					u_approach_x = particle[np].center[0].vel_x + particle[np].center[0].omgZ * (y_center_curve_inertial - particle[np].center[0].y) - particle[np].center[0].omgY * (z_center_curve_inertial - particle[np].center[0].z);
					//cout << "salam A " << x_intersect_real_inertial << " " << dist_to_wall << " " << " " << Rij << " " << R_curvature << " " << x_0_inertial<<endl;
					double dot_pro = (u_approach_x - 0)*Rij;
					ffx4 = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature*R_curvature / (R_curvature + R_curvature))  * (1. / (abs(Rij) - 2 * R_curvature) - 1. / thw_lub) * dot_pro *Rij / SQ(Rij);
					ffy4 = 0.0;
					ffz4 = 0.0;
				}
#endif
				if (abs(Rij) <= (2.0 * R_curvature + thw) && abs(Rij) > (2.0 * R_curvature)) {
					ffx4 = cij / epw * pow((2.0 * R_curvature - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
					ffy4 = 0.0;
					ffz4 = 0.0;
				}
				if (abs(Rij) <= (2.0 * R_curvature)) {
					ffx4 = (cij / epw *pow((abs(Rij) - 2.0 * R_curvature - thw) / thw, 2) + cij / Epw / thw*(2.0*R_curvature - abs(Rij)))*Rij / abs(Rij);
					ffy4 = 0.0;
					ffz4 = 0.0;
				}
				break;
			}
			case 5: {// Lower Z-Wall
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (abs(Rij) <= (2.0 * R_curvature + thw_lub) && abs(Rij)>(2.0 * R_curvature + thw)) {
					u_approach_z = particle[np].center[0].vel_z - particle[np].center[0].omgX * (y_center_curve_inertial - particle[np].center[0].y) + particle[np].center[0].omgY * (x_center_curve_inertial - particle[np].center[0].x);
					double dot_pro = (u_approach_z - 0)*Rij;
					ffx4 = 0.0;
					ffy4 = 0.0;
					ffz4 = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature*R_curvature / (R_curvature + R_curvature))  * (1. / (abs(Rij) - 2 * R_curvature) - 1. / thw_lub) * dot_pro *Rij / SQ(Rij);
				}
#endif
				if (abs(Rij) <= (2.0 * R_curvature + thw) && abs(Rij) > (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = 0.0;
					ffz4 = cij / epw * pow((2.0 * R_curvature - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
				}
				if (abs(Rij) <= (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = 0.0;
					ffz4 = (cij / epw *pow((abs(Rij) - 2.0 * R_curvature - thw) / thw, 2) + cij / Epw / thw*(2.0*R_curvature - abs(Rij)))*Rij / abs(Rij);
				}
				break;
			}
			case 6: {// Upper Z-Wall
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (abs(Rij) <= (2.0 * R_curvature + thw_lub) && abs(Rij)>(2.0 * R_curvature + thw)) {
					u_approach_z = particle[np].center[0].vel_z - particle[np].center[0].omgX * (y_center_curve_inertial - particle[np].center[0].y) + particle[np].center[0].omgY * (x_center_curve_inertial - particle[np].center[0].x);
					double dot_pro = (u_approach_z - 0)*Rij;
					ffx4 = 0.0;
					ffy4 = 0.0;
					ffz4 = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature*R_curvature / (R_curvature + R_curvature))  * (1. / (abs(Rij) - 2 * R_curvature) - 1. / thw_lub) * dot_pro *Rij / SQ(Rij);
				}
#endif
				if (abs(Rij) <= (2.0 * R_curvature + thw) && abs(Rij) > (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = 0.0;
					ffz4 = cij / epw * pow((2.0 * R_curvature - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
				}
				if (abs(Rij) <= (2.0 * R_curvature)) {
					ffx4 = 0.0;
					ffy4 = 0.0;
					ffz4 = (cij / epw *pow((abs(Rij) - 2.0 * R_curvature - thw) / thw, 2) + cij / Epw / thw*(2.0*R_curvature - abs(Rij)))*Rij / abs(Rij);
				}
				break;
			}
			default: std::cout << "Keine Antwort";
			}

			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].f4z += ffz4;

			particle[np].center[0].t3_x += ffy4 * (z_intersect_real_inertial - particle[np].center[0].z) - ffz4 * (y_intersect_real_inertial - particle[np].center[0].y);
			particle[np].center[0].t3_y += ffz4 * (x_intersect_real_inertial - particle[np].center[0].x) - ffx4 * (z_intersect_real_inertial - particle[np].center[0].z);
			particle[np].center[0].t3_z += ffx4 * (y_intersect_real_inertial - particle[np].center[0].y) - ffy4 * (x_intersect_real_inertial - particle[np].center[0].x);
			//cout << dist_to_wall << " " << particle[np].center[0].f4x << endl;

		}  // if (dist_to_wall<thw_control)
	}///for (size_t i_ch = 0; i_ch < check_walls.size(); i_ch++)
#endif  //Walls collision



	 //------------------------------------------------------------------
	 //------------------------------------------------------------------
	 //------------------------------------------------------------------
	 //-------------- Particle_Particle collision -----------------/////
	int max_par_index;
	int index_p;
	double par_par_dist;
	double xp1_on_p2, yp1_on_p2, zp1_on_p2;
	double x_0_of_p1, y_0_of_p1, z_0_of_p1;
	double x_out_on_p2, y_out_on_p2, z_out_on_p2;
	double m_slope_p1, n_slope_p1, p_slope_p1;
	double x_intersect1_of_p1, y_intersect1_of_p1, z_intersect1_of_p1;
	double x_intersect2_of_p1, y_intersect2_of_p1, z_intersect2_of_p1;
	double x_intersect_real_of_p1, y_intersect_real_of_p1, z_intersect_real_of_p1;
	double xnor_slope_old_of_p1, ynor_slope_old_of_p1, znor_slope_old_of_p1;
	double xnor_slope_of_p1, ynor_slope_of_p1, znor_slope_of_p1;
	double x_0_inertial_of_p1, y_0_inertial_of_p1, z_0_inertial_of_p1;
	double x_intersect_real_inertial_of_p1, y_intersect_real_inertial_of_p1, z_intersect_real_inertial_of_p1;
	double xp1_on_p1, yp1_on_p1, zp1_on_p1;
	double temp_slope_length;
	double xnor_slope_of_p2, ynor_slope_of_p2, znor_slope_of_p2;
	double x_0_of_p2, y_0_of_p2, z_0_of_p2;
	double x_out_on_p1, y_out_on_p1, z_out_on_p1;
	double m_slope_p2, n_slope_p2, p_slope_p2;
	double x_intersect1_of_p2, y_intersect1_of_p2, z_intersect1_of_p2;
	double x_intersect2_of_p2, y_intersect2_of_p2, z_intersect2_of_p2;
	double x_intersect_real_of_p2, y_intersect_real_of_p2, z_intersect_real_of_p2;
	double xnor_slope_old_of_p2, ynor_slope_old_of_p2, znor_slope_old_of_p2;
	double x_0_inertial_of_p2, y_0_inertial_of_p2, z_0_inertial_of_p2;
	double x_intersect_real_inertial_of_p2, y_intersect_real_inertial_of_p2, z_intersect_real_inertial_of_p2;
	double fi_curvature_p1, fi_curvature_p2;
	double R_curvature_p1, R_curvature_p2;
	double thp_control;
	double x_center_curve_p1, y_center_curve_p1, z_center_curve_p1;
	double x_center_curve_inertial_p1, y_center_curve_inertial_p1, z_center_curve_inertial_p1;
	double x_center_curve_inertial_p2, y_center_curve_inertial_p2, z_center_curve_inertial_p2;
	double x_center_curve_p2, y_center_curve_p2, z_center_curve_p2;
	double u_approach_x_p1, u_approach_y_p1, u_approach_z_p1;
	double u_approach_x_p2, u_approach_y_p2, u_approach_z_p2;

#if defined Spheroid_Collision_1_Feng || defined Spheroid_Collision_2_Feng_Verlet
	thp_control = th;
#endif
#if defined Spheroid_Collision_3_Feng_Lub_Verlet
	thp_control = th_lub;
#endif

	ffx = 0.0;
	ffy = 0.0;
	ffz = 0.0;

#if defined Spheroid_Collision_1_Feng
	max_par_index = num_particles;
#endif
#if defined Spheroid_Collision_2_Feng_Verlet || defined Spheroid_Collision_3_Feng_Lub_Verlet
	max_par_index = verlet_list[np].size();
#endif

	for (int np2 = 0; np2 < max_par_index; ++np2) {//Particle-Particle Collision
#if defined Spheroid_Collision_1_Feng
		index_p = np2;
#endif
#if defined Spheroid_Collision_2_Feng_Verlet || defined Spheroid_Collision_3_Feng_Lub_Verlet
		index_p = verlet_list[np][np2];
#endif
		if (index_p != np) {
			par_par_dist = 1000;

			xp1_on_p2 = particle[index_p].center[0].x - particle[np].center[0].x;//point on P2 in inertial coordinate minus particle P1 center position
			yp1_on_p2 = particle[index_p].center[0].y - particle[np].center[0].y;
			zp1_on_p2 = particle[index_p].center[0].z - particle[np].center[0].z;

			xnor_slope_of_p1 = particle[np].center[0].x - particle[index_p].center[0].x;//slope of line connecting nearest points
			ynor_slope_of_p1 = particle[np].center[0].y - particle[index_p].center[0].y;
			znor_slope_of_p1 = particle[np].center[0].z - particle[index_p].center[0].z;

			temp_slope_length = sqrt(SQ(particle[np].center[0].x - particle[index_p].center[0].x) + SQ(particle[np].center[0].y - particle[index_p].center[0].y) + SQ(particle[np].center[0].z - particle[index_p].center[0].z));
			xnor_slope_of_p1 /= temp_slope_length;
			ynor_slope_of_p1 /= temp_slope_length;
			znor_slope_of_p1 /= temp_slope_length;

			x_0_of_p1 = 0;//means particle center position in body-fixed frame
			y_0_of_p1 = 0;//means particle center position in body-fixed frame
			z_0_of_p1 = 0;//means particle center position in body-fixed frame

			x_0_of_p2 = 0;//means particle center position in body-fixed frame
			y_0_of_p2 = 0;//means particle center position in body-fixed frame
			z_0_of_p2 = 0;//means particle center position in body-fixed frame

			xnor_slope_of_p2 = particle[index_p].center[0].x - particle[np].center[0].x;//slope of line connecting nearest points
			ynor_slope_of_p2 = particle[index_p].center[0].y - particle[np].center[0].y;
			znor_slope_of_p2 = particle[index_p].center[0].z - particle[np].center[0].z;

			//initial set. Changes later
			temp_slope_length = sqrt(SQ(particle[index_p].center[0].x - particle[np].center[0].x) + SQ(particle[index_p].center[0].y - particle[np].center[0].y) + SQ(particle[index_p].center[0].z - particle[np].center[0].z));
			xnor_slope_of_p2 /= temp_slope_length;
			ynor_slope_of_p2 /= temp_slope_length;
			znor_slope_of_p2 /= temp_slope_length;

			mmm = 0;

			do {
				x_out_on_p2 = particle[np].center[0].quater[1] * (xp1_on_p2)+particle[np].center[0].quater[2] * (yp1_on_p2)+particle[np].center[0].quater[3] * (zp1_on_p2);//location of P2 point in body-fixed frame of P1
				y_out_on_p2 = particle[np].center[0].quater[4] * (xp1_on_p2)+particle[np].center[0].quater[5] * (yp1_on_p2)+particle[np].center[0].quater[6] * (zp1_on_p2);
				z_out_on_p2 = particle[np].center[0].quater[7] * (xp1_on_p2)+particle[np].center[0].quater[8] * (yp1_on_p2)+particle[np].center[0].quater[9] * (zp1_on_p2);
				m_slope_p1 = x_out_on_p2 - x_0_of_p1; //slope of the line conecting the point on P2 and the center of imaginary circle of P1
				n_slope_p1 = y_out_on_p2 - y_0_of_p1;
				p_slope_p1 = z_out_on_p2 - z_0_of_p1;

				/// Intersection of ellipsoid and line, Obtained by MATLAB  ///
				tt_part1 = a*sqrt(SQ(a) * SQ(b) * SQ(n_slope_p1) + SQ(a) *SQ(b)  * SQ(p_slope_p1) - SQ(a) * SQ(n_slope_p1) * SQ(z_0_of_p1) + 2 * SQ(a) * n_slope_p1*p_slope_p1*y_0_of_p1*z_0_of_p1 - SQ(a * p_slope_p1 * y_0_of_p1) + SQ(SQ(b))  * SQ(m_slope_p1) - SQ(b* m_slope_p1* y_0_of_p1) - SQ(b * m_slope_p1  * z_0_of_p1) + 2 * SQ(b) * m_slope_p1*n_slope_p1*x_0_of_p1*y_0_of_p1 + 2 * SQ(b) * m_slope_p1*p_slope_p1*x_0_of_p1*z_0_of_p1 - SQ(b * n_slope_p1  * x_0_of_p1) - SQ(b * p_slope_p1 * x_0_of_p1));
				tt_part2 = SQ(b) * m_slope_p1*x_0_of_p1 + SQ(a) * n_slope_p1*y_0_of_p1 + SQ(a) * p_slope_p1*z_0_of_p1;
				tt_part3 = (SQ(a  * n_slope_p1) + SQ(a * p_slope_p1) + SQ(b * m_slope_p1));
				tt_intersect_1 = -(tt_part1 + tt_part2) / tt_part3;//Answers of equation
				tt_intersect_2 = -(-tt_part1 + tt_part2) / tt_part3;

				x_intersect1_of_p1 = x_0_of_p1 + m_slope_p1*tt_intersect_1;//in body fixed frame of P1
				y_intersect1_of_p1 = y_0_of_p1 + n_slope_p1*tt_intersect_1;
				z_intersect1_of_p1 = z_0_of_p1 + p_slope_p1*tt_intersect_1;

				x_intersect2_of_p1 = x_0_of_p1 + m_slope_p1*tt_intersect_2;
				y_intersect2_of_p1 = y_0_of_p1 + n_slope_p1*tt_intersect_2;
				z_intersect2_of_p1 = z_0_of_p1 + p_slope_p1*tt_intersect_2;

				/// Find the nearst point
				if ((SQ(x_intersect1_of_p1 - x_out_on_p2) + SQ(y_intersect1_of_p1 - y_out_on_p2) + SQ(z_intersect1_of_p1 - z_out_on_p2)) < (SQ(x_intersect2_of_p1 - x_out_on_p2) + SQ(y_intersect2_of_p1 - y_out_on_p2) + SQ(z_intersect2_of_p1 - z_out_on_p2)))
				{
					x_intersect_real_of_p1 = x_intersect1_of_p1;//in body fixed frame
					y_intersect_real_of_p1 = y_intersect1_of_p1;//in body fixed frame
					z_intersect_real_of_p1 = z_intersect1_of_p1;//in body fixed frame
				}
				else
				{
					x_intersect_real_of_p1 = x_intersect2_of_p1;//in body fixed frame
					y_intersect_real_of_p1 = y_intersect2_of_p1;//in body fixed frame
					z_intersect_real_of_p1 = z_intersect2_of_p1;//in body fixed frame
				}
				xnor_slope_old_of_p1 = xnor_slope_of_p1;
				ynor_slope_old_of_p1 = ynor_slope_of_p1;
				znor_slope_old_of_p1 = znor_slope_of_p1;

				xnor_slope_of_p1 = -2.*x_intersect_real_of_p1 / particle_radius;//slope of normal vector on spheroid in body_fixed frame. INWARD direction
				ynor_slope_of_p1 = -2.*y_intersect_real_of_p1 / particle_radius_2;
				znor_slope_of_p1 = -2.*z_intersect_real_of_p1 / particle_radius_2;

				nor_slope_length = (sqrt(SQ(xnor_slope_of_p1) + SQ(ynor_slope_of_p1) + SQ(znor_slope_of_p1)));
				xnor_slope_of_p1 /= nor_slope_length;
				ynor_slope_of_p1 /= nor_slope_length;
				znor_slope_of_p1 /= nor_slope_length;

				x_0_of_p1 = x_intersect_real_of_p1 + xnor_slope_of_p1 * radius_check;//center of an imaginary circle for next check. (body-fixed frame)
				y_0_of_p1 = y_intersect_real_of_p1 + ynor_slope_of_p1 * radius_check;
				z_0_of_p1 = z_intersect_real_of_p1 + znor_slope_of_p1 * radius_check;

				//finding the center of imaginary circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
				update_rot_spheroid(particle[np].center[0].deter, x_0_inertial_of_p1, y_0_inertial_of_p1, z_0_inertial_of_p1, particle[np].center[0].quater, x_0_of_p1, y_0_of_p1, z_0_of_p1, particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Find x_0 position in inertial frame

																																																																			//finding the intersection point of P1 in INERTIAL frame (x_intersect_real_inertial_of_p1, y_intersect_real_inertial_of_p1, z_intersect_real_inertial_of_p1)
				update_rot_spheroid(particle[np].center[0].deter, x_intersect_real_inertial_of_p1, y_intersect_real_inertial_of_p1, z_intersect_real_inertial_of_p1, particle[np].center[0].quater, x_intersect_real_of_p1, y_intersect_real_of_p1, z_intersect_real_of_p1, particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Find x_0 position in inertial frame


																																																																																						  //////////----------------------------------------------------------------------//////////////////////
																																																																																						  //Now we start the process for P2
				xp1_on_p1 = x_intersect_real_inertial_of_p1 - particle[index_p].center[0].x;//point on P1 in inertial coordinate minus particle P2 center position
				yp1_on_p1 = y_intersect_real_inertial_of_p1 - particle[index_p].center[0].y;
				zp1_on_p1 = z_intersect_real_inertial_of_p1 - particle[index_p].center[0].z;

				x_out_on_p1 = particle[index_p].center[0].quater[1] * (xp1_on_p1)+particle[index_p].center[0].quater[2] * (yp1_on_p1)+particle[index_p].center[0].quater[3] * (zp1_on_p1);//location of P1 point in body-fixed frame of P2
				y_out_on_p1 = particle[index_p].center[0].quater[4] * (xp1_on_p1)+particle[index_p].center[0].quater[5] * (yp1_on_p1)+particle[index_p].center[0].quater[6] * (zp1_on_p1);
				z_out_on_p1 = particle[index_p].center[0].quater[7] * (xp1_on_p1)+particle[index_p].center[0].quater[8] * (yp1_on_p1)+particle[index_p].center[0].quater[9] * (zp1_on_p1);

				m_slope_p2 = x_out_on_p1 - x_0_of_p2; //slope of the line conecting the point on P1 and the center of imaginary circle of P2
				n_slope_p2 = y_out_on_p1 - y_0_of_p2;
				p_slope_p2 = z_out_on_p1 - z_0_of_p2;

				///// Intersection of ellipsoid and line, Obtained by MATLAB  ///
				tt_part1 = a*sqrt(SQ(a) * SQ(b) * SQ(n_slope_p2) + SQ(a) *SQ(b)  * SQ(p_slope_p2) - SQ(a) * SQ(n_slope_p2) * SQ(z_0_of_p2) + 2 * SQ(a) * n_slope_p2*p_slope_p2*y_0_of_p2*z_0_of_p2 - SQ(a * p_slope_p2 * y_0_of_p2) + SQ(SQ(b))  * SQ(m_slope_p2) - SQ(b* m_slope_p2* y_0_of_p2) - SQ(b * m_slope_p2  * z_0_of_p2) + 2 * SQ(b) * m_slope_p2*n_slope_p2*x_0_of_p2*y_0_of_p2 + 2 * SQ(b) * m_slope_p2*p_slope_p2*x_0_of_p2*z_0_of_p2 - SQ(b * n_slope_p2  * x_0_of_p2) - SQ(b * p_slope_p2 * x_0_of_p2));
				tt_part2 = SQ(b) * m_slope_p2*x_0_of_p2 + SQ(a) * n_slope_p2*y_0_of_p2 + SQ(a) * p_slope_p2*z_0_of_p2;
				tt_part3 = (SQ(a  * n_slope_p2) + SQ(a * p_slope_p2) + SQ(b * m_slope_p2));
				tt_intersect_1 = -(tt_part1 + tt_part2) / tt_part3;//Answers of equation
				tt_intersect_2 = -(-tt_part1 + tt_part2) / tt_part3;

				x_intersect1_of_p2 = x_0_of_p2 + m_slope_p2*tt_intersect_1;//in body fixed frame of P1
				y_intersect1_of_p2 = y_0_of_p2 + n_slope_p2*tt_intersect_1;
				z_intersect1_of_p2 = z_0_of_p2 + p_slope_p2*tt_intersect_1;

				x_intersect2_of_p2 = x_0_of_p2 + m_slope_p2*tt_intersect_2;
				y_intersect2_of_p2 = y_0_of_p2 + n_slope_p2*tt_intersect_2;
				z_intersect2_of_p2 = z_0_of_p2 + p_slope_p2*tt_intersect_2;

				///// Find the nearst point
				if ((SQ(x_intersect1_of_p2 - x_out_on_p1) + SQ(y_intersect1_of_p2 - y_out_on_p1) + SQ(z_intersect1_of_p2 - z_out_on_p1)) < (SQ(x_intersect2_of_p2 - x_out_on_p1) + SQ(y_intersect2_of_p2 - y_out_on_p1) + SQ(z_intersect2_of_p2 - z_out_on_p1)))
				{
					x_intersect_real_of_p2 = x_intersect1_of_p2;//in body fixed frame
					y_intersect_real_of_p2 = y_intersect1_of_p2;//in body fixed frame
					z_intersect_real_of_p2 = z_intersect1_of_p2;//in body fixed frame
				}
				else
				{
					x_intersect_real_of_p2 = x_intersect2_of_p2;//in body fixed frame
					y_intersect_real_of_p2 = y_intersect2_of_p2;//in body fixed frame
					z_intersect_real_of_p2 = z_intersect2_of_p2;//in body fixed frame
				}

				xnor_slope_old_of_p2 = xnor_slope_of_p2;
				ynor_slope_old_of_p2 = ynor_slope_of_p2;
				znor_slope_old_of_p2 = znor_slope_of_p2;

				xnor_slope_of_p2 = -2.*x_intersect_real_of_p2 / particle_radius;//slope of normal vector on spheroid in body_fixed frame. INWARD direction
				ynor_slope_of_p2 = -2.*y_intersect_real_of_p2 / particle_radius_2;
				znor_slope_of_p2 = -2.*z_intersect_real_of_p2 / particle_radius_2;

				nor_slope_length = (sqrt(SQ(xnor_slope_of_p2) + SQ(ynor_slope_of_p2) + SQ(znor_slope_of_p2)));
				xnor_slope_of_p2 /= nor_slope_length;
				ynor_slope_of_p2 /= nor_slope_length;
				znor_slope_of_p2 /= nor_slope_length;

				x_0_of_p2 = x_intersect_real_of_p2 + xnor_slope_of_p2 * radius_check;//center of an imaginary circle for next check. (body-fixed frame)
				y_0_of_p2 = y_intersect_real_of_p2 + ynor_slope_of_p2 * radius_check;
				z_0_of_p2 = z_intersect_real_of_p2 + znor_slope_of_p2 * radius_check;

				////finding the center of imaginary circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
				update_rot_spheroid(particle[index_p].center[0].deter, x_0_inertial_of_p2, y_0_inertial_of_p2, z_0_inertial_of_p2, particle[index_p].center[0].quater, x_0_of_p2, y_0_of_p2, z_0_of_p2, particle[index_p].center[0].x, particle[index_p].center[0].y, particle[index_p].center[0].z);// Find x_0 position in inertial frame

																																																																									 //finding the intersection point of P2 in INERTIAL frame (x_intersect_real_inertial_of_p2, y_intersect_real_inertial_of_p2, z_intersect_real_inertial_of_p2)
				update_rot_spheroid(particle[index_p].center[0].deter, x_intersect_real_inertial_of_p2, y_intersect_real_inertial_of_p2, z_intersect_real_inertial_of_p2, particle[index_p].center[0].quater, x_intersect_real_of_p2, y_intersect_real_of_p2, z_intersect_real_of_p2, particle[index_p].center[0].x, particle[index_p].center[0].y, particle[index_p].center[0].z);// Find x_0 position in inertial frame

				xp1_on_p2 = x_intersect_real_inertial_of_p2 - particle[np].center[0].x;
				yp1_on_p2 = y_intersect_real_inertial_of_p2 - particle[np].center[0].y;
				zp1_on_p2 = z_intersect_real_inertial_of_p2 - particle[np].center[0].z;

				mmm++;
			} while (abs(xnor_slope_of_p1 - xnor_slope_old_of_p1) >= slope_th || abs(ynor_slope_of_p1 - ynor_slope_old_of_p1) >= slope_th || abs(znor_slope_of_p1 - znor_slope_old_of_p1) >= slope_th || abs(xnor_slope_of_p2 - xnor_slope_old_of_p2) >= slope_th || abs(ynor_slope_of_p2 - ynor_slope_old_of_p2) >= slope_th || abs(znor_slope_of_p2 - znor_slope_old_of_p2) >= slope_th);

			par_par_dist = sqrt(SQ(x_intersect_real_inertial_of_p1 - x_intersect_real_inertial_of_p2) + SQ(y_intersect_real_inertial_of_p1 - y_intersect_real_inertial_of_p2) + SQ(z_intersect_real_inertial_of_p1 - z_intersect_real_inertial_of_p2));

			if (par_par_dist<thp_control) {

				fi_curvature_p1 = asin(x_intersect_real_of_p1 / (sqrt(SQ(x_intersect_real_of_p1) + SQ(y_intersect_real_of_p1) + SQ(z_intersect_real_of_p1))));
				fi_curvature_p2 = asin(x_intersect_real_of_p2 / (sqrt(SQ(x_intersect_real_of_p2) + SQ(y_intersect_real_of_p2) + SQ(z_intersect_real_of_p2))));

				R_curvature_p1 = SQ(a)*b / (SQ(a*sin(fi_curvature_p1)) + SQ(b*cos(fi_curvature_p1))); //Paper 1366
				R_curvature_p2 = SQ(a)*b / (SQ(a*sin(fi_curvature_p2)) + SQ(b*cos(fi_curvature_p2))); //Paper 1366

				x_center_curve_p1 = x_intersect_real_of_p1 + xnor_slope_of_p1 * R_curvature_p1;//center of circle of curvature. (body-fixed frame)
				y_center_curve_p1 = y_intersect_real_of_p1 + ynor_slope_of_p1 * R_curvature_p1;
				z_center_curve_p1 = z_intersect_real_of_p1 + znor_slope_of_p1 * R_curvature_p1;

				x_center_curve_p2 = x_intersect_real_of_p2 + xnor_slope_of_p2 * R_curvature_p2;//center of circle of curvature. (body-fixed frame)
				y_center_curve_p2 = y_intersect_real_of_p2 + ynor_slope_of_p2 * R_curvature_p2;
				z_center_curve_p2 = z_intersect_real_of_p2 + znor_slope_of_p2 * R_curvature_p2;

				//finding the center of curvature circle in INERTIAL frame
				update_rot_spheroid(particle[np].center[0].deter, x_center_curve_inertial_p1, y_center_curve_inertial_p1, z_center_curve_inertial_p1, particle[np].center[0].quater, x_center_curve_p1, y_center_curve_p1, z_center_curve_p1, particle[np].center[0].x, particle[np].center[0].y, particle[np].center[0].z);// Find  position of circle of curvature in inertial frame
				update_rot_spheroid(particle[index_p].center[0].deter, x_center_curve_inertial_p2, y_center_curve_inertial_p2, z_center_curve_inertial_p2, particle[index_p].center[0].quater, x_center_curve_p2, y_center_curve_p2, z_center_curve_p2, particle[index_p].center[0].x, particle[index_p].center[0].y, particle[index_p].center[0].z);// Find  position of circle of curvature in inertial frame

				Dij = sqrt(SQ(x_center_curve_inertial_p1 - x_center_curve_inertial_p2) + SQ(y_center_curve_inertial_p1 - y_center_curve_inertial_p2) + SQ(z_center_curve_inertial_p1 - z_center_curve_inertial_p2));

#if defined Spheroid_Collision_3_Feng_Lub_Verlet
				if (par_par_dist <= th_lub && par_par_dist>th) {
					u_approach_x_p1 = particle[np].center[0].vel_x + particle[np].center[0].omgZ * (y_center_curve_inertial_p1 - particle[np].center[0].y) - particle[np].center[0].omgY * (z_center_curve_inertial_p1 - particle[np].center[0].z);
					u_approach_y_p1 = particle[np].center[0].vel_y - particle[np].center[0].omgZ * (x_center_curve_inertial_p1 - particle[np].center[0].x) + particle[np].center[0].omgX * (z_center_curve_inertial_p1 - particle[np].center[0].z);
					u_approach_z_p1 = particle[np].center[0].vel_z - particle[np].center[0].omgX * (y_center_curve_inertial_p1 - particle[np].center[0].y) + particle[np].center[0].omgY * (x_center_curve_inertial_p1 - particle[np].center[0].x);

					u_approach_x_p2 = particle[index_p].center[0].vel_x + particle[np].center[0].omgZ * (y_center_curve_inertial_p2 - particle[index_p].center[0].y) - particle[np].center[0].omgY * (z_center_curve_inertial_p2 - particle[index_p].center[0].z);
					u_approach_y_p2 = particle[index_p].center[0].vel_y - particle[np].center[0].omgZ * (x_center_curve_inertial_p2 - particle[index_p].center[0].x) + particle[np].center[0].omgX * (z_center_curve_inertial_p2 - particle[index_p].center[0].z);
					u_approach_z_p2 = particle[index_p].center[0].vel_z - particle[np].center[0].omgX * (y_center_curve_inertial_p2 - particle[index_p].center[0].y) + particle[np].center[0].omgY * (x_center_curve_inertial_p2 - particle[index_p].center[0].x);

					double dot_pro = ((u_approach_x_p1 - u_approach_x_p2)*(x_center_curve_inertial_p1 - x_center_curve_inertial_p2) + (u_approach_y_p1 - u_approach_y_p2)*(y_center_curve_inertial_p1 - y_center_curve_inertial_p2) + (u_approach_z_p1 - u_approach_z_p2)*(z_center_curve_inertial_p1 - z_center_curve_inertial_p2));
					ffx = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature_p1*R_curvature_p2 / (R_curvature_p1 + R_curvature_p2))  * (1. / (par_par_dist)-1. / th_lub) * dot_pro * (x_center_curve_inertial_p1 - x_center_curve_inertial_p2) / SQ(Dij);
					ffy = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature_p1*R_curvature_p2 / (R_curvature_p1 + R_curvature_p2))  * (1. / (par_par_dist)-1. / th_lub) * dot_pro * (y_center_curve_inertial_p1 - y_center_curve_inertial_p2) / SQ(Dij);
					ffz = -6 * M_PI* (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature_p1*R_curvature_p2 / (R_curvature_p1 + R_curvature_p2))  * (1. / (par_par_dist)-1. / th_lub) * dot_pro * (z_center_curve_inertial_p1 - z_center_curve_inertial_p2) / SQ(Dij);

				}
#endif
				if (par_par_dist <= th && par_par_dist>0) {
					ffx = cij / ep * ((x_center_curve_inertial_p1 - x_center_curve_inertial_p2)) / Dij * pow((par_par_dist - th) / th, 2);
					ffy = cij / ep * ((y_center_curve_inertial_p1 - y_center_curve_inertial_p2)) / Dij * pow((par_par_dist - th) / th, 2);
					ffz = cij / ep * ((z_center_curve_inertial_p1 - z_center_curve_inertial_p2)) / Dij * pow((par_par_dist - th) / th, 2);
				}
				if (par_par_dist < 0) {
					ffx = (x_center_curve_inertial_p1 - x_center_curve_inertial_p2) / Dij * (cij / ep * pow((par_par_dist - th) / th, 2) + cij / Ep / th*(-1 * par_par_dist));
					ffy = (y_center_curve_inertial_p1 - y_center_curve_inertial_p2) / Dij * (cij / ep * pow((par_par_dist - th) / th, 2) + cij / Ep / th*(-1 * par_par_dist));
					ffz = (z_center_curve_inertial_p1 - z_center_curve_inertial_p2) / Dij * (cij / ep * pow((par_par_dist - th) / th, 2) + cij / Ep / th*(-1 * par_par_dist));
				}

				particle[np].center[0].f3x += ffx;
				particle[np].center[0].f3y += ffy;
				particle[np].center[0].f3z += ffz;

				particle[np].center[0].t3_x += ffy * (z_intersect_real_inertial_of_p1 - particle[np].center[0].z) - ffz * (y_intersect_real_inertial_of_p1 - particle[np].center[0].y);
				particle[np].center[0].t3_y += ffz * (x_intersect_real_inertial_of_p1 - particle[np].center[0].x) - ffx * (z_intersect_real_inertial_of_p1 - particle[np].center[0].z);
				particle[np].center[0].t3_z += ffx * (y_intersect_real_inertial_of_p1 - particle[np].center[0].y) - ffy * (x_intersect_real_inertial_of_p1 - particle[np].center[0].x);

			} //// if (par_par_dist<thp_control)
		}  ///// if (index_p != np)
	}// END of Particle-Particle Collision
	return;
}
#endif
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: RUNGEKUTTA
//-----------------------------------------------------------------------------------------------------------------------------
//this function solves 3+4 equations for "3:rotational speed in body-fixed coordinate" + "4:Quaternions:q0,q1,q2,q3"

void rungekutta(double & omega_prm_x, double & omega_prm_y, double & omega_prm_z,
	double Tor_prm_x, double Tor_prm_y, double Tor_prm_z,
	double mom_inr_x, double mom_inr_y, double mom_inr_z, double & Q0, double & Q1, double & Q2, double & Q3) {

	double nn, dh, f, k1, k2, k3, k4;
	double g, r1, r2, r3, r4;
	double h, s1, s2, s3, s4;
	double m, b1, b2, b3, b4;
	double n, c1, c2, c3, c4;
	double p, j1, j2, j3, j4;
	double q, t1, t2, t3, t4;
	double omg_prm_x_1, omg_prm_y_1, omg_prm_z_1;
	double Q0_1, Q1_1, Q2_1, Q3_1;
	double Tx = Tor_prm_x;
	double Ty = Tor_prm_y;
	double Tz = Tor_prm_z;
	double mom_x = mom_inr_x;
	double mom_y = mom_inr_y;
	double mom_z = mom_inr_z;

#define FF(omg_prm_y,omg_prm_z)  Tx/mom_x + omg_prm_y * omg_prm_z * (mom_y - mom_z)/mom_x
#define GG(omg_prm_x,omg_prm_z)  Ty/mom_y + omg_prm_x * omg_prm_z * (mom_z - mom_x)/mom_y
#define HH(omg_prm_x,omg_prm_y)  Tz/mom_z + omg_prm_x * omg_prm_y * (mom_x - mom_y)/mom_z
#define MM(Q1,Q2,Q3,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*(-Q1*omg_prm_x - Q2*omg_prm_y - Q3*omg_prm_z)
#define NN(Q0,Q2,Q3,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*( Q0*omg_prm_x - Q3*omg_prm_y + Q2*omg_prm_z)
#define PP(Q0,Q1,Q3,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*( Q3*omg_prm_x + Q0*omg_prm_y - Q1*omg_prm_z)
#define QQ(Q0,Q1,Q2,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*(-Q2*omg_prm_x + Q1*omg_prm_y + Q0*omg_prm_z)

	double omg_prm_x_0 = omega_prm_x;//Initial value
	double omg_prm_y_0 = omega_prm_y;//Initial value
	double omg_prm_z_0 = omega_prm_z;//Initial value
	double Q0_0 = Q0;//Initial value
	double Q1_0 = Q1;//Initial value
	double Q2_0 = Q2;//Initial value
	double Q3_0 = Q3;//Initial value

	dh = 1;//dh here must be time step which is in lattice calculations = 1

		   //printf("\nEnter the value of last point: ");
		   //scanf("%lf",&nn);
		   //for(double t=0; t<nn; t=t+dh)
		   //{

	f = FF(omg_prm_y_0, omg_prm_z_0);
	k1 = dh * f;
	g = GG(omg_prm_x_0, omg_prm_z_0);
	r1 = dh * g;
	h = HH(omg_prm_x_0, omg_prm_y_0);
	s1 = dh * h;
	m = MM(Q1_0, Q2_0, Q3_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
	b1 = dh * m;
	n = NN(Q0_0, Q2_0, Q3_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
	c1 = dh * n;
	p = PP(Q0_0, Q1_0, Q3_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
	j1 = dh * p;
	q = QQ(Q0_0, Q1_0, Q2_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
	t1 = dh * q;

	f = FF((omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
	k2 = dh * f;
	g = GG((omg_prm_x_0 + k1 / 2), (omg_prm_z_0 + s1 / 2));
	r2 = dh * g;
	h = HH((omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2));
	s2 = dh * h;
	m = MM((Q1_0 + c1 / 2), (Q2_0 + j1 / 2), (Q3_0 + t1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
	b2 = dh * m;
	n = NN((Q0_0 + b1 / 2), (Q2_0 + j1 / 2), (Q3_0 + t1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
	c2 = dh * n;
	p = PP((Q0_0 + b1 / 2), (Q1_0 + c1 / 2), (Q3_0 + t1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
	j2 = dh * p;
	q = QQ((Q0_0 + b1 / 2), (Q1_0 + c1 / 2), (Q2_0 + j1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
	t2 = dh * q;

	f = FF((omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
	k3 = dh * f;
	g = GG((omg_prm_x_0 + k2 / 2), (omg_prm_z_0 + s2 / 2));
	r3 = dh * g;
	h = HH((omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2));
	s3 = dh * h;
	m = MM((Q1_0 + c2 / 2), (Q2_0 + j2 / 2), (Q3_0 + t2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
	b3 = dh * m;
	n = NN((Q0_0 + b2 / 2), (Q2_0 + j2 / 2), (Q3_0 + t2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
	c3 = dh * n;
	p = PP((Q0_0 + b2 / 2), (Q1_0 + c2 / 2), (Q3_0 + t2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
	j3 = dh * p;
	q = QQ((Q0_0 + b2 / 2), (Q1_0 + c2 / 2), (Q2_0 + j2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
	t3 = dh * q;

	f = FF((omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
	k4 = dh * f;
	g = GG((omg_prm_x_0 + k3), (omg_prm_z_0 + s3));
	r4 = dh * g;
	h = HH((omg_prm_x_0 + k3), (omg_prm_y_0 + r3));
	s4 = dh * h;
	m = MM((Q1_0 + c3), (Q2_0 + j3), (Q3_0 + t3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
	b4 = dh * m;
	n = NN((Q0_0 + b3), (Q2_0 + j3), (Q3_0 + t3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
	c4 = dh * n;
	p = PP((Q0_0 + b3), (Q1_0 + c3), (Q3_0 + t3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
	j4 = dh * p;
	q = QQ((Q0_0 + b3), (Q1_0 + c3), (Q2_0 + j3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
	t4 = dh * q;

	omg_prm_x_1 = omg_prm_x_0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
	omg_prm_y_1 = omg_prm_y_0 + (r1 + 2 * r2 + 2 * r3 + r4) / 6.;
	omg_prm_z_1 = omg_prm_z_0 + (s1 + 2 * s2 + 2 * s3 + s4) / 6.;
	Q0_1 = Q0_0 + (b1 + 2 * b2 + 2 * b3 + b4) / 6.;
	Q1_1 = Q1_0 + (c1 + 2 * c2 + 2 * c3 + c4) / 6.;
	Q2_1 = Q2_0 + (j1 + 2 * j2 + 2 * j3 + j4) / 6.;
	Q3_1 = Q3_0 + (t1 + 2 * t2 + 2 * t3 + t4) / 6.;

	omg_prm_x_0 = omg_prm_x_1;
	omg_prm_y_0 = omg_prm_y_1;
	omg_prm_z_0 = omg_prm_z_1;
	Q0_0 = Q0_1;
	Q1_0 = Q1_1;
	Q2_0 = Q2_1;
	Q3_0 = Q3_1;
	// }
	omega_prm_x = omg_prm_x_1;
	omega_prm_y = omg_prm_y_1;
	omega_prm_z = omg_prm_z_1;
	Q0 = Q0_1;
	Q1 = Q1_1;
	Q2 = Q2_1;
	Q3 = Q3_1;
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: QUATERNION
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function calculates the elements of quaternion matrix.(Eq.2.13 Article 876: H.Huang et al. J. Fluid Mech. (2012), vol. 692, pp. 369-394.
void quaternion(double Q0, double Q1, double Q2, double Q3, int np, double quat[]) {


	quat[0] = 0.;//this value is useless
	quat[1] = SQ(Q0) + SQ(Q1) - SQ(Q2) - SQ(Q3);
	quat[2] = 2 * (Q1*Q2 + Q0*Q3);
	quat[3] = 2 * (Q1*Q3 - Q0*Q2);
	quat[4] = 2 * (Q1*Q2 - Q0*Q3);
	quat[5] = SQ(Q0) - SQ(Q1) + SQ(Q2) - SQ(Q3);
	quat[6] = 2 * (Q2*Q3 + Q0*Q1);
	quat[7] = 2 * (Q1*Q3 + Q0*Q2);
	quat[8] = 2 * (Q2*Q3 - Q0*Q1);
	quat[9] = SQ(Q0) - SQ(Q1) - SQ(Q2) + SQ(Q3);

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: DETERMINANT
//-----------------------------------------------------------------------------------------------------------------------------
//This function calculates the determinant of quaternion matrix (in order to obtain its inverse)

double determinant(int np, double quat[]) {
	double det;
	det = quat[1] * (quat[5] * quat[9] - quat[6] * quat[8])
		- quat[2] * (quat[4] * quat[9] - quat[6] * quat[7])
		+ quat[3] * (quat[4] * quat[8] - quat[5] * quat[7]);
	return det;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COMP_ROT_SPEED
//-----------------------------------------------------------------------------------------------------------------------------
// (This function transforms the rotational speed from body-fixed coordinate to Inertial coordinate
// This function is just for SPHEROID
void comp_rot_speed(double determi, double & omgX, double & omgY, double & omgZ, int n, double quat[], double omg_pr_x, double omg_pr_y, double omg_pr_z) {

	double RNI[10];
	double om_temp_x = omg_pr_x;
	double om_temp_y = omg_pr_y;
	double om_temp_z = omg_pr_z;

	RNI[0] = 0.;
	RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
	RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
	RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
	RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
	RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
	RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
	RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
	RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
	RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
	RNI[1] = RNI[1] / determi;
	RNI[2] = RNI[2] / determi;
	RNI[3] = RNI[3] / determi;
	RNI[4] = RNI[4] / determi;
	RNI[5] = RNI[5] / determi;
	RNI[6] = RNI[6] / determi;
	RNI[7] = RNI[7] / determi;
	RNI[8] = RNI[8] / determi;
	RNI[9] = RNI[9] / determi;

	omgX = RNI[1] * om_temp_x + RNI[2] * om_temp_y + RNI[3] * om_temp_z;
	omgY = RNI[4] * om_temp_x + RNI[5] * om_temp_y + RNI[6] * om_temp_z;
	omgZ = RNI[7] * om_temp_x + RNI[8] * om_temp_y + RNI[9] * om_temp_z;

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: UPDATE_ROT_SPHEROID
//-----------------------------------------------------------------------------------------------------------------------------
// If at t=0 Euler angles are not zero then the spheroid must be rotated. This function is part of caculating the Lagrangian node position in INERTIAL coordinate based on Euler angles and node location in BODY coordinate
// This function is just for SPHEROID
void update_rot_spheroid(double determi, double & X_inertial, double & Y_inertial, double & Z_inertial, double quat[], double X_body, double Y_body, double Z_body
	, double center_x, double center_y, double center_z) {
	// X_body is node X position in body coordinate minus particle center
	// X_inertial is node X position in INERTIAL coordinate

	double RNI[10];

	RNI[0] = 0.;
	RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
	RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
	RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
	RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
	RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
	RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
	RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
	RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
	RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
	RNI[1] = RNI[1] / determi;
	RNI[2] = RNI[2] / determi;
	RNI[3] = RNI[3] / determi;
	RNI[4] = RNI[4] / determi;
	RNI[5] = RNI[5] / determi;
	RNI[6] = RNI[6] / determi;
	RNI[7] = RNI[7] / determi;
	RNI[8] = RNI[8] / determi;
	RNI[9] = RNI[9] / determi;

	X_inertial = (RNI[1] * X_body + RNI[2] * Y_body + RNI[3] * Z_body) + center_x;
	Y_inertial = (RNI[4] * X_body + RNI[5] * Y_body + RNI[6] * Z_body) + center_y;
	Z_inertial = (RNI[7] * X_body + RNI[8] * Y_body + RNI[9] * Z_body) + center_z;

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: NORMAL_VEC_SPHEROID
//-----------------------------------------------------------------------------------------------------------------------------
// Calculates NORMAL vector at each node at INERTIAL coordinate based on the body coordinate valuse
void normal_vec_spheroid(double determi, double & X_inertial, double & Y_inertial, double & Z_inertial, int n, double quat[], double X_body, double Y_body, double Z_body) {
	// X_body is normal vector X index in body coordinate
	// X_inertial normal vector X index in INERTIAL coordinate

	double RNI[10];

	RNI[0] = 0.;
	RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
	RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
	RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
	RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
	RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
	RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
	RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
	RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
	RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
	RNI[1] = RNI[1] / determi;
	RNI[2] = RNI[2] / determi;
	RNI[3] = RNI[3] / determi;
	RNI[4] = RNI[4] / determi;
	RNI[5] = RNI[5] / determi;
	RNI[6] = RNI[6] / determi;
	RNI[7] = RNI[7] / determi;
	RNI[8] = RNI[8] / determi;
	RNI[9] = RNI[9] / determi;

	X_inertial = (RNI[1] * X_body + RNI[2] * Y_body + RNI[3] * Z_body);
	Y_inertial = (RNI[4] * X_body + RNI[5] * Y_body + RNI[6] * Z_body);
	Z_inertial = (RNI[7] * X_body + RNI[8] * Y_body + RNI[9] * Z_body);

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: CHECK_ZERO
//-----------------------------------------------------------------------------------------------------------------------------
// This function checks if a variable is less than epsilon it is set to zero
void check_zero(double & variab, double epsilon) {
	if (abs(variab)<epsilon) {
		variab = 0.0;
	}
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SREAD
//-----------------------------------------------------------------------------------------------------------------------------
//
// The node forces are spread to the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.
#if defined Flow_With_Particle
void spread(particle_struct & particle) {


	int x_int;
	int y_int;
	int z_int;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	int zmin;
	int zmax;
	double dist_x;
	double dist_y;
	double dist_z;
	double weight_x;
	double weight_y;
	double weight_z;

	/// Spread forces
	// Run over all object nodes.
	for (int n = 0; n < particle_num_nodes; ++n) {

		// Identify the lowest fluid lattice node in interpolation range.
		// 'Lowest' means: its x- and y-values are the smallest.
		// The other fluid nodes in range have coordinates
		// (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).

#if defined Kernel_Delta || defined Two_Point_Delta || defined Four_Point_Delta	|| defined Modified_Two_Point_Delta
		x_int = floor(particle.node[n].x);
		y_int = floor(particle.node[n].y);
		z_int = floor(particle.node[n].z);
		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
#if defined Kernel_Delta || defined Four_Point_Delta || defined Modified_Two_Point_Delta
		xmin = x_int - 2;
		xmax = x_int + 2;
		ymin = y_int - 2;
		ymax = y_int + 2;
		zmin = z_int - 2;
		zmax = z_int + 2;
#endif

#if defined Two_Point_Delta
		xmin = x_int - 1;
		xmax = x_int + 1;
		ymin = y_int - 1;
		ymax = y_int + 1;
		zmin = z_int - 1;
		zmax = z_int + 1;
#endif
#if defined X_Wall
		if (xmin <= 0) {
			xmin = 0;
		}
		if (xmax >= Nx - 1) {
			xmax = Nx - 1;
		}
#endif
#if defined Y_Wall
		if (ymin <= 0) {
			ymin = 0;
		}
		if (ymax >= Ny - 1) {
			ymax = Ny - 1;
		}
#endif
#if defined Z_Wall
		if (zmin <= 0) {
			zmin = 0;
		}
		if (zmax >= Nz - 1) {
			zmax = Nz - 1;
		}
#endif
#endif


		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
		for (int X = xmin; X <= xmax; ++X)
			for (int Y = ymin; Y <= ymax; ++Y)
				for (int Z = zmin; Z <= zmax; ++Z)
				{
					{
						{

#if defined Kernel_Delta || defined Two_Point_Delta || defined Four_Point_Delta || defined Modified_Two_Point_Delta
							// Compute distance between object node and fluid lattice node.
							dist_x = particle.node[n].x - X;
							dist_y = particle.node[n].y - Y;
							dist_z = particle.node[n].z - Z;
							// Compute interpolation weights for x- and y-direction based on the distance.

#if defined Kernel_Delta
							if (abs(dist_x) <= 2) {
								weight_x = 1. / 4.0*(1. + cos(M_PI*abs(dist_x) / 2));
							}
							else weight_x = 0.0;
							if (abs(dist_y) <= 2) {
								weight_y = 1. / 4.0*(1. + cos(M_PI*abs(dist_y) / 2));
							}
							else weight_y = 0.0;
							if (abs(dist_z) <= 2) {
								weight_z = 1. / 4.0*(1. + cos(M_PI*abs(dist_z) / 2));
							}
							else weight_z = 0.0;
#endif

#if defined Four_Point_Delta
							if (abs(dist_x) >= 0 && abs(dist_x) < 1) {
								weight_x = 1. / 8.0*(3. - 2.* abs(dist_x) + sqrt(1 + 4.*abs(dist_x) - 4.*SQ(dist_x)));
							}
							else if (abs(dist_x) >= 1 && abs(dist_x) < 2) {
								weight_x = 1. / 8.0*(5. - 2.* abs(dist_x) - sqrt(-7 + 12.*abs(dist_x) - 4.*SQ(dist_x)));
							}
							else weight_x = 0.0;

							if (abs(dist_y) >= 0 && abs(dist_y) < 1) {
								weight_y = 1. / 8.0*(3. - 2.* abs(dist_y) + sqrt(1 + 4.*abs(dist_y) - 4.*SQ(dist_y)));
							}
							else if (abs(dist_y) >= 1 && abs(dist_y) < 2) {
								weight_y = 1. / 8.0*(5. - 2.* abs(dist_y) - sqrt(-7 + 12.*abs(dist_y) - 4.*SQ(dist_y)));
							}
							else weight_y = 0.0;

							if (abs(dist_z) >= 0 && abs(dist_z) < 1) {
								weight_z = 1. / 8.0*(3. - 2.* abs(dist_z) + sqrt(1 + 4.*abs(dist_z) - 4.*SQ(dist_z)));
							}
							else if (abs(dist_z) >= 1 && abs(dist_z) < 2) {
								weight_z = 1. / 8.0*(5. - 2.* abs(dist_z) - sqrt(-7 + 12.*abs(dist_z) - 4.*SQ(dist_z)));
							}
							else weight_z = 0.0;
#endif

#if defined Two_Point_Delta
							if (abs(dist_x) <= 1) {
								weight_x = 1. - abs(dist_x);
							}
							else weight_x = 0.0;
							if (abs(dist_y) <= 1) {
								weight_y = 1. - abs(dist_y);
							}
							else weight_y = 0.0;
							if (abs(dist_z) <= 1) {
								weight_z = 1. - abs(dist_z);
							}
							else weight_z = 0.0;
#endif
#if defined Modified_Two_Point_Delta
							if (abs(dist_x) <= 0.5) {
								weight_x = 3. / 4. - (dist_x)*(dist_x);
							}
							else if (abs(dist_x) > 0.5 && abs(dist_x) <= 1.5) {
								weight_x = 9. / 8. - 3.*abs(dist_x) / 2. + SQ(dist_x) / 2.;
							}
							else weight_x = 0.0;

							if (abs(dist_y) <= 0.5) {
								weight_y = 3. / 4. - (dist_y)*(dist_y);
							}
							else if (abs(dist_y) > 0.5 && abs(dist_y) <= 1.5) {
								weight_y = 9. / 8. - 3.*abs(dist_y) / 2. + SQ(dist_y) / 2.;
							}
							else weight_y = 0.0;

							if (abs(dist_z) <= 0.5) {
								weight_z = 3. / 4. - (dist_z)*(dist_z);
							}
							else if (abs(dist_z) > 0.5 && abs(dist_z) <= 1.5) {
								weight_z = 9. / 8. - 3.*abs(dist_z) / 2. + SQ(dist_z) / 2.;
							}
							else weight_z = 0.0;

#endif
#endif

							force_x[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz] += (particle.node[n].force_x * weight_x * weight_y * weight_z);
							force_y[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz] += (particle.node[n].force_y * weight_x * weight_y * weight_z);
							force_z[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz] += (particle.node[n].force_z * weight_x * weight_y * weight_z);
#if defined Flow_With_Thermal_Effect
							force_thermal[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz] += (particle.node[n].force_thermal * weight_x * weight_y* weight_z);
#endif
						}
					}
				}
	}

	return;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: INTERPOLATE
//-----------------------------------------------------------------------------------------------------------------------------
//
// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) and other interpolation schemes are included in the present code.
#if defined Flow_With_Particle

void interpolate(particle_struct & particle) {

	int x_int;
	int y_int;
	int z_int;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	int zmin;
	int zmax;
	double dist_x;
	double dist_y;
	double dist_z;
	double weight_x;
	double weight_y;
	double weight_z;

	// Run over all object nodes.

	for (int n = 0; n < particle_num_nodes; ++n) {

		// Reset node velocity first since '+=' is used.
		particle.node[n].vel_x = 0;
		particle.node[n].vel_y = 0;
		particle.node[n].vel_z = 0;
		particle.node[n].den = 0;
		particle.node[n].temperature = 0;


		// Identify the lowest fluid lattice node in interpolation range (see spreading).

#if defined Kernel_Delta || defined Two_Point_Delta || defined Four_Point_Delta	|| defined Modified_Two_Point_Delta
		x_int = floor(particle.node[n].x);
		y_int = floor(particle.node[n].y);
		z_int = floor(particle.node[n].z);
		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
#if defined Kernel_Delta || defined Four_Point_Delta || defined Modified_Two_Point_Delta
		xmin = x_int - 2;
		xmax = x_int + 2;
		ymin = y_int - 2;
		ymax = y_int + 2;
		zmin = z_int - 2;
		zmax = z_int + 2;
#endif

#if defined Two_Point_Delta
		xmin = x_int - 1;
		xmax = x_int + 1;
		ymin = y_int - 1;
		ymax = y_int + 1;
		zmin = z_int - 1;
		zmax = z_int + 1;
#endif
#if defined X_Wall
		if (xmin <= 0) {
			xmin = 0;
		}
		if (xmax >= Nx - 1) {
			xmax = Nx - 1;
		}
#endif
#if defined Y_Wall

		if (ymin <= 0) {
			ymin = 0;
		}
		if (ymax >= Ny - 1) {
			ymax = Ny - 1;
		}
#endif
#if defined Z_Wall
		if (zmin <= 0) {
			zmin = 0;
		}
		if (zmax >= Nz - 1) {
			zmax = Nz - 1;
		}
#endif

#endif


		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
		for (int X = xmin; X <= xmax; ++X)
			for (int Y = ymin; Y <= ymax; ++Y)
				for (int Z = zmin; Z <= zmax; ++Z)
				{
					{
						{

#if defined Kernel_Delta || defined Two_Point_Delta || defined Four_Point_Delta	|| defined Modified_Two_Point_Delta
							// Compute distance between object node and fluid lattice node.
							dist_x = particle.node[n].x - X;
							dist_y = particle.node[n].y - Y;
							dist_z = particle.node[n].z - Z;
							// Compute interpolation weights for x- and y-direction based on the distance.

#if defined Kernel_Delta
							if (abs(dist_x) <= 2) {
								weight_x = 1. / 4.0*(1. + cos(M_PI*abs(dist_x) / 2));
							}
							else weight_x = 0.0;
							if (abs(dist_y) <= 2) {
								weight_y = 1. / 4.0*(1. + cos(M_PI*abs(dist_y) / 2));
							}
							else weight_y = 0.0;
							if (abs(dist_z) <= 2) {
								weight_z = 1. / 4.0*(1. + cos(M_PI*abs(dist_z) / 2));
							}
							else weight_z = 0.0;
#endif

#if defined Four_Point_Delta
							if (abs(dist_x) >= 0 && abs(dist_x) < 1) {
								weight_x = 1. / 8.0*(3. - 2.* abs(dist_x) + sqrt(1 + 4.*abs(dist_x) - 4.*SQ(dist_x)));
							}
							else if (abs(dist_x) >= 1 && abs(dist_x) < 2) {
								weight_x = 1. / 8.0*(5. - 2.* abs(dist_x) - sqrt(-7 + 12.*abs(dist_x) - 4.*SQ(dist_x)));
							}
							else weight_x = 0.0;

							if (abs(dist_y) >= 0 && abs(dist_y) < 1) {
								weight_y = 1. / 8.0*(3. - 2.* abs(dist_y) + sqrt(1 + 4.*abs(dist_y) - 4.*SQ(dist_y)));
							}
							else if (abs(dist_y) >= 1 && abs(dist_y) < 2) {
								weight_y = 1. / 8.0*(5. - 2.* abs(dist_y) - sqrt(-7 + 12.*abs(dist_y) - 4.*SQ(dist_y)));
							}
							else weight_y = 0.0;

							if (abs(dist_z) >= 0 && abs(dist_z) < 1) {
								weight_z = 1. / 8.0*(3. - 2.* abs(dist_z) + sqrt(1 + 4.*abs(dist_z) - 4.*SQ(dist_z)));
							}
							else if (abs(dist_z) >= 1 && abs(dist_z) < 2) {
								weight_z = 1. / 8.0*(5. - 2.* abs(dist_z) - sqrt(-7 + 12.*abs(dist_z) - 4.*SQ(dist_z)));
							}
							else weight_z = 0.0;
#endif

#if defined Two_Point_Delta
							if (abs(dist_x) <= 1) {
								weight_x = 1. - abs(dist_x);
							}
							else weight_x = 0.0;
							if (abs(dist_y) <= 1) {
								weight_y = 1. - abs(dist_y);
							}
							else weight_y = 0.0;
							if (abs(dist_z) <= 1) {
								weight_z = 1. - abs(dist_z);
							}
							else weight_z = 0.0;
#endif
#if defined Modified_Two_Point_Delta
							if (abs(dist_x) <= 0.5) {
								weight_x = 3. / 4. - (dist_x)*(dist_x);
							}
							else if (abs(dist_x) > 0.5 && abs(dist_x) <= 1.5) {
								weight_x = 9. / 8. - 3.*abs(dist_x) / 2. + SQ(dist_x) / 2.;
							}
							else weight_x = 0.0;

							if (abs(dist_y) <= 0.5) {
								weight_y = 3. / 4. - (dist_y)*(dist_y);
							}
							else if (abs(dist_y) > 0.5 && abs(dist_y) <= 1.5) {
								weight_y = 9. / 8. - 3.*abs(dist_y) / 2. + SQ(dist_y) / 2.;
							}
							else weight_y = 0.0;

							if (abs(dist_z) <= 0.5) {
								weight_z = 3. / 4. - (dist_z)*(dist_z);
							}
							else if (abs(dist_z) > 0.5 && abs(dist_z) <= 1.5) {
								weight_z = 9. / 8. - 3.*abs(dist_z) / 2. + SQ(dist_z) / 2.;
							}
							else weight_z = 0.0;

#endif
#endif

							particle.node[n].vel_x += ((velocity_x[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz]) * weight_x * weight_y * weight_z);
							particle.node[n].vel_y += ((velocity_y[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz]) * weight_x * weight_y * weight_z);
							particle.node[n].vel_z += ((velocity_z[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz]) * weight_x * weight_y * weight_z);
							particle.node[n].den += ((density[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz]) * weight_x * weight_y * weight_z);
#if defined Flow_With_Thermal_Effect
							particle.node[n].temperature += ((temperature[(X + Nx) % Nx][(Y + Ny) % Ny][(Z + Nz) % Nz]) * weight_x * weight_y* weight_z);
#endif
						}
					}
				}
	}
	return;
}
#endif

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SWAP_ENDIAN
//-----------------------------------------------------------------------------------------------------------------------------
void swap_endian(char * buffer, size_t size) {
	for (size_t b = 0; b < size / 2; b++) {
		char temp = buffer[b];
		buffer[b] = buffer[size - 1 - b];
		buffer[size - 1 - b] = temp;
	}
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SWAP_ENDIAN
//-----------------------------------------------------------------------------------------------------------------------------
void swap_endian(char * buffer, char * buffer2, size_t size) {
	for (size_t b = 0; b < size; b++) {
		buffer2[b] = buffer[size - 1 - b];
	}

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_BIGENDIAN
//-----------------------------------------------------------------------------------------------------------------------------

void write_bigendian(std::ofstream & file, char * buffer, size_t count, size_t size) {
	char * _buf = new char[size];
	for (size_t i = 0; i < count; i++) {
		swap_endian(buffer + (size*i), _buf, size);
		file.write(_buf, size);
	}
	delete[] _buf;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_BIGENDIAN
//-----------------------------------------------------------------------------------------------------------------------------
template <typename T>
void write_bigendian(std::ofstream & file, T * buffer, size_t count) {
	char * _buf = new char[sizeof(T)];
	for (size_t i = 0; i < count; i++) {
		swap_endian((char*)(buffer + i), _buf, sizeof(T));
		file.write(_buf, sizeof(T));
	}
	delete[] _buf;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// The fluid state is writen to a VTK file at each t_vtk step.
// The following data is written:
// - density difference (density - 1)
// - x-component of velocity
// - y-component of velocity
// The following code is designed in such a way that the file can be read by ParaView.

void write_fluid_vtk(int time) {


#if defined VTK_BINARY
	/// Create filename

	stringstream output_filename;
	output_filename << "Alborz_Results/vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "BINARY\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
	output_file << "X_COORDINATES " << Nx << " float\n";
	output_file.close();

	//{
	//	std::stringstream str_line1;
	//	str_line1 << "# vtk DataFile Version 3.0\n";
	//	std::string strstr_line1 = str_line1.str();
	//	output_file1.write(strstr_line1.c_str(), sizeof(char)*strstr_line1.size());

	//output_file1.write((char*)z1, sizeof(float)* 1);
	//}

	float XX[Nx];
	float YY[Ny];
	float ZZ[Nz];

	for (int X = 0; X < Nx; ++X) {
		XX[X] = X;
	}
	for (int Y = 0; Y < Ny; ++Y) {
		YY[Y] = Y;
	}
	for (int Z = 0; Z < Nz; ++Z) {
		ZZ[Z] = Z;
	}

	/// Open file for append
	std::stringstream str_line1;
	str_line1 << "Alborz_Results/vtk_fluid/fluid_t" << time << ".vtk";
	std::string strstr_line1 = str_line1.str();

	std::ofstream output_file1(strstr_line1, std::ios::app | std::ios::binary);

	//	std::ofstream output_file1("Alborz_Results/vtk_fluid/fluid_t200.vtk", std::ios::app | std::ios::binary);

	write_bigendian(output_file1, XX, Nx);

	output_file1 << "\nY_COORDINATES " << Ny << " float\n";

	write_bigendian(output_file1, YY, Ny);

	output_file1 << "\nZ_COORDINATES " << Nz << " float\n";

	write_bigendian(output_file1, ZZ, Nz);

	output_file1 << "\nPOINT_DATA " << Nx * Ny *Nz << "\n";

	/// Write density difference
	output_file1 << "SCALARS density_difference float 1\n";
	output_file1 << "LOOKUP_TABLE default\n";

	float mm[1];
	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				mm[0] = density[X][Y][Z] - 1;
				write_bigendian(output_file1, mm, 1);
			}
		}
	}

	output_file1 << "\nSCALARS temperature float 1\n";
	output_file1 << "LOOKUP_TABLE default\n";
	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				mm[0] = temperature[X][Y][Z];
				write_bigendian(output_file1, mm, 1);
			}
		}
	}

	float nn[2];
	///// Write velocity
	output_file1 << "VECTORS velocity_vector float\n";
	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				mm[0] = velocity_x[X][Y][Z];
				write_bigendian(output_file1, mm, 1);
				mm[0] = velocity_y[X][Y][Z];
				write_bigendian(output_file1, mm, 1);
				mm[0] = velocity_z[X][Y][Z];
				write_bigendian(output_file1, mm, 1);
			}
		}
	}

	output_file1.close();
#endif

#if defined VTK_ASCII

	/// Create filename

	stringstream output_filename;
	output_filename << "Alborz_Results/vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file

	output_file.open(output_filename.str().c_str());

	/// Write VTK header

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";

	//output_file << "DIMENSIONS " << Nx << " " << Ny - 2  << " " << Nz - 2 << "\n";
	output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";

	output_file << "X_COORDINATES " << Nx << " float\n";

	for (int X = 0; X < Nx; ++X) {
		output_file << X << " ";
	}

	output_file << "\n";
	output_file << "Y_COORDINATES " << Ny << " float\n";

	for (int Y = 0; Y < Ny; ++Y) {
		output_file << Y << " ";
	}

	output_file << "\n";
	output_file << "Z_COORDINATES " << Nz << " float\n";

	for (int Z = 0; Z < Nz; ++Z) {
		output_file << Z << " ";
	}


	//output_file << "\n";
	//output_file << "Z_COORDINATES " << 1 << " double\n";
	//output_file << 0 << "\n";
	output_file << "POINT_DATA " << Nx * Ny * Nz << "\n";

	/// Write density difference

	output_file << "SCALARS density_difference float 1\n";
	output_file << "LOOKUP_TABLE default\n";

	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				output_file << density[X][Y][Z] - 1 << "\n";
			}
		}
	}

	/// Write temperature

	output_file << "SCALARS temperature float 1\n";
	output_file << "LOOKUP_TABLE default\n";

	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				output_file << temperature[X][Y][Z] << "\n";
			}
		}
	}

	/// Write velocity

	output_file << "VECTORS velocity_vector float\n";

	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				//#if defined Shan_Chen
				//				output_file << velocity_x[X][Y][Z] + 0.5 * (force_x[X][Y][Z] + gravity) / density[X][Y][Z] << " " << velocity_y[X][Y][Z] + 0.5 * (force_y[X][Y][Z]) / density[X][Y][Z] << " " << velocity_z[X][Y][Z] + 0.5 * (force_z[X][Y][Z]) / density[X][Y][Z] << "\n";
				//#endif
				//#if (defined Guo || defined Luo)
				output_file << velocity_x[X][Y][Z] << " " << velocity_y[X][Y][Z] << " " << velocity_z[X][Y][Z] << "\n";
				//#endif
			}
		}
	}
	/// Close file

	output_file.close();
#endif
#if defined CASE
	unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);// number are shown as 2.0e02 instead of 2.0e002 (necessay for geo and data files)

																		  /////////// Writing CASE file (.case) ///////////////////
																		  /////////////////////////////////////////////////////////
																		  /// Create filename

	stringstream output_filename;
	output_filename << "Alborz_Results/vtk_fluid/fluid" << ".case";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# EnSight Gold Model\n";
	output_file << "# Created by ALBORZ\n\n";
	output_file << "FORMAT\n";
	output_file << "type:                   ensight gold\n\n";
	output_file << "GEOMETRY\n";
	output_file << "model:                  fluid.geo\n\n";
	output_file << "VARIABLE\n";
	output_file << "scalar per node:           1 density_difference fluid.density_difference******\n";
	output_file << "scalar per node:           1 temperature fluid.temperature******\n";
	output_file << "vector per node:           1 VELOC_VECTOR fluid.VELOC_VECTOR******\n\n";
	output_file << "TIME\n";
	output_file << "time set:               1\n";
	output_file << "number of steps:        " << t_num / t_vtk << "\n";
	output_file << "filename start number:  " << t_vtk << "\n";
	output_file << "filename increment:     " << t_vtk << "\n";
	output_file << "time values:\n";
	for (int itm = t_vtk; itm <= t_num; itm = itm + t_vtk) {
		output_file << itm << "\n";

	}
	/// Close file
	output_file.close();
	/////////// End of Writing CASE file (.case) ///////////////
	////////////////////////////////////////////////////////////

	if (time == t_vtk) {//Because .geo file do not change during iterations

						/////////// Writing GEOMETRY file (.geo) ///////////////////
						////////////////////////////////////////////////////////////
		stringstream output_filename1;
		output_filename1 << "Alborz_Results/vtk_fluid/fluid" << ".geo";//Geometry file
		ofstream output_file1;

		/// Open file
		output_file1.open(output_filename1.str().c_str());

		/// Write VTK header
		output_file1 << "EnSight Model Geometry File\n";
		output_file1 << "EnSight 10.1.1\n";
		output_file1 << "node id given\n";
		output_file1 << "element id given\n";
		output_file1 << "extents\n";
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << (double)Nx - 1. << "\n";
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << (double)Ny - 1. << "\n";
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << (double)Nz - 1. << "\n";
		output_file1 << "part\n";
		output_file1 << std::setw(10) << 1 << "\n";
		output_file1 << "FLUID\n";
		output_file1 << "coordinates\n";
		output_file1 << std::setw(10) << Nx*Ny*Nz << "\n";
		for (int itm = 1; itm <= Nx*Ny*Nz; itm++) {
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int Z = 0; Z < Nz; Z++) {
			for (int Y = 0; Y < Ny; Y++) {
				for (int X = 0; X < Nx; X++) {
					output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << (double)X << "\n";
				}
			}
		}
		for (int Z = 0; Z < Nz; Z++) {
			for (int Y = 0; Y < Ny; Y++) {
				for (int X = 0; X < Nx; X++) {
					output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << (double)Y << "\n";
				}
			}
		}

		for (int Z = 0; Z < Nz; Z++) {
			for (int Y = 0; Y < Ny; Y++) {
				for (int X = 0; X < Nx; X++) {
					output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << (double)Z << "\n";
				}
			}
		}
		output_file1 << "hexa8\n";
		output_file1 << std::setw(10) << (Nx - 1)*(Ny - 1)*(Nz - 1) << "\n";
		for (int itm = 1; itm <= (Nx - 1)*(Ny - 1)*(Nz - 1); itm++) {
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int K = 1; K < Nz; K++) {
			for (int J = 1; J < Ny; J++) {
				for (int I = 1; I < Nx; I++) {
					output_file1 << std::setw(10) << I + (J - 1)*Nx + (K - 1)*Nx*Ny << std::setw(10) << I + 1 + (J - 1)*Nx + (K - 1)*Nx*Ny << std::setw(10) << (I + 1) + (J*Nx) + (K - 1)*Nx*Ny << std::setw(10) << (I)+(J*Nx) + (K - 1)*Nx*Ny << std::setw(10) << I + (J - 1)*Nx + (K)*Nx*Ny << std::setw(10) << I + 1 + (J - 1)*Nx + (K)*Nx*Ny << std::setw(10) << I + 1 + (J)*Nx + (K)*Nx*Ny << std::setw(10) << I + (J)*Nx + (K)*Nx*Ny << "\n";
				}
			}
		}
		/// Close file
		output_file1.close();
		/////////// End of Writing GEOMETRY file (.geo) ///////////////////
		///////////////////////////////////////////////////////////////////

	}// End of if (time == t_vtk)


	 ///////////// Writing DENSITY_DIFFERENCE file ////////////////////////
	 ////////////////////////////////////////////////////////////
	stringstream output_filename6;
	output_filename6 << "Alborz_Results/vtk_fluid/fluid.density_difference" << setfill('0') << setw(6) << time;
	ofstream output_file6;

	/// Open file
	output_file6.open(output_filename6.str().c_str());

	output_file6 << "density_difference\n";
	output_file6 << "part\n";
	output_file6 << std::setw(10) << 1 << "\n";
	output_file6 << "coordinates\n";
	for (int Z = 0; Z < Nz; Z++) {
		for (int Y = 0; Y < Ny; Y++) {
			for (int X = 0; X < Nx; X++) {
				output_file6 << std::scientific << std::setw(12) << std::setprecision(5) << density[X][Y][Z] - 1. << "\n";
			}
		}
	}

	/// Close file
	output_file6.close();

	/////////// End of Writing DENSITY_DIFFERENCE file /////////////////
	////////////////////////////////////////////////////////////

	///////////// Writing TEMPERATURE file ////////////////////////
	/////////////////////////////////////////////////////////////////////
	stringstream output_filename2;
	output_filename2 << "Alborz_Results/vtk_fluid/fluid.temperature" << setfill('0') << setw(6) << time;
	ofstream output_file2;

	/// Open file
	output_file2.open(output_filename2.str().c_str());

	output_file2 << "temperature\n";
	output_file2 << "part\n";
	output_file2 << std::setw(10) << 1 << "\n";
	output_file2 << "coordinates\n";

	for (int Z = 0; Z < Nz; Z++) {
		for (int Y = 0; Y < Ny; Y++) {
			for (int X = 0; X < Nx; X++) {
				output_file2 << std::scientific << std::setw(12) << std::setprecision(5) << temperature[X][Y][Z] << "\n";
			}
		}
	}

	/// Close file
	output_file2.close();

	/////////// End of Writing TEMPERATURE file /////////////////
	////////////////////////////////////////////////////////////

	///////////// Writing VELOCITY_VECTOR file ////////////////////////
	/////////////////////////////////////////////////////////////////////
	stringstream output_filename5;
	output_filename5 << "Alborz_Results/vtk_fluid/fluid.VELOC_VECTOR" << setfill('0') << setw(6) << time;
	ofstream output_file5;

	/// Open file
	output_file5.open(output_filename5.str().c_str());

	output_file5 << "VELOC_VECTOR\n";
	output_file5 << "part\n";
	output_file5 << std::setw(10) << 1 << "\n";
	output_file5 << "coordinates\n";

	//#if (defined Guo || defined Luo)
	for (int Z = 0; Z < Nz; Z++) {
		for (int Y = 0; Y < Ny; Y++) {
			for (int X = 0; X < Nx; X++) {
				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_x[X][Y][Z] << "\n";
			}
		}
	}
	for (int Z = 0; Z < Nz; Z++) {
		for (int Y = 0; Y < Ny; Y++) {
			for (int X = 0; X < Nx; X++) {
				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_y[X][Y][Z] << "\n";
			}
		}
	}
	for (int Z = 0; Z < Nz; Z++) {
		for (int Y = 0; Y < Ny; Y++) {
			for (int X = 0; X < Nx; X++) {
				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_z[X][Y][Z] << "\n";
			}
		}
	}
	//#endif
	//#if defined Shan_Chen
	//	for (int Z = 0; Z < Nz; Z++){
	//		for (int Y = 0; Y < Ny; Y++){
	//			for (int X = 0; X < Nx; X++){
	//				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_x[X][Y][Z] + 0.5 * (force_x[X][Y][Z] + gravity) / density[X][Y][Z] << "\n";
	//			}
	//		}
	//	}
	//	for (int Z = 0; Z < Nz; Z++){
	//		for (int Y = 0; Y < Ny; Y++){
	//			for (int X = 0; X < Nx; X++){
	//				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_y[X][Y][Z] + 0.5 * (force_y[X][Y][Z]) / density[X][Y][Z] << "\n";
	//			}
	//		}
	//	}
	//	for (int Z = 0; Z < Nz; Z++){
	//		for (int Y = 0; Y < Ny; Y++){
	//			for (int X = 0; X < Nx; X++){
	//				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_z[X][Y][Z] + 0.5 * (force_z[X][Y][Z]) / density[X][Y][Z] << "\n";
	//			}
	//		}
	//	}
	//#endif

	/// Close file
	output_file5.close();

	/////////// End of Writing VELOCITY_VECTOR file /////////////
	////////////////////////////////////////////////////////////

#endif

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// The particle state (node positions) is writen to a VTK file at each t_vtk step.
// The following code is designed in such a way that the file can be read by ParaView.
#if defined Flow_With_Particle

void write_particle_vtk(int time, particle_struct particle, int npp) {
#if defined VTK_BINARY || defined VTK_ASCII

	/// Create filename

	stringstream output_filename;
	output_filename << "Alborz_Results/vtk_particle/particle_" << npp + 1 << "_t" << time << ".vtk";
	ofstream output_file;

	/// Open file

	output_file.open(output_filename.str().c_str());

	/// Write VTK header

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "particle_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET POLYDATA\n";

	/// Write node positions
	output_file << "POINTS " << particle_num_nodes << " float\n";

	for (int n = 0; n < particle_num_nodes; ++n) {
		output_file << particle.node[n].x << " " << particle.node[n].y << " " << particle.node[n].z << "\n";
	}

#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
	output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
	output_file << particle_num_nodes << "\n";
	for (int n = 0; n < particle_num_nodes; ++n) {
		output_file << n << "\n";
	}

	// Below part is useful for SPHERE to create surface using Glyph for canter in Paraview
	stringstream output_filename1;
	output_filename1 << "Alborz_Results/vtk_particle/particle_center_" << npp + 1 << "_t" << time << ".vtk";
	ofstream output_file1;

	/// Open file
	output_file1.open(output_filename1.str().c_str());

	/// Write VTK header
	output_file1 << "# vtk DataFile Version 3.0\n";
	output_file1 << "particle_state\n";
	output_file1 << "ASCII\n";
	output_file1 << "DATASET POLYDATA\n";

	output_file1 << "POINTS " << "2" << " float\n";
	output_file1 << particle.center[0].x << " " << particle.center[0].y << " " << particle.center[0].z << "\n";
	output_file1 << particle.center[0].x + 0.000001 << " " << particle.center[0].y << " " << particle.center[0].z << "\n";

	/// Close file
	output_file1.close();
#endif

#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
	// Create triangle surfaces
	output_file << "POLYGONS " << number_triangles << " " << 4 * number_triangles << "\n";
	for (int n = 0; n < number_triangles; ++n) {
		output_file << "3 " << triangleset1[n] << " " << triangleset2[n] << " " << triangleset3[n] << "\n";
	}
#endif
	/// Close file
	output_file.close();
#endif

#if defined CASE
	unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);// number are shown as 2.0e02 instead of 2.0e002 (necessay for geo and data files)


																		  /////////// Writing CASE file (.case) ///////////////////
																		  /////////////////////////////////////////////////////////
																		  /// Create filename

	stringstream output_filename;
	output_filename << "Alborz_Results/vtk_particle/particle_" << npp + 1 << ".case";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# EnSight Gold Model\n";
	output_file << "# Created by ALBORZ\n\n";
	output_file << "FORMAT\n";
	output_file << "type:                   ensight gold\n\n";
	output_file << "GEOMETRY\n";
	output_file << "model:                 1 particle_" << npp + 1 << "_******.geo\n\n";
	output_file << "TIME\n";
	output_file << "time set:               1\n";
	output_file << "number of steps:        " << t_num / t_vtk << "\n";
	output_file << "filename start number:  " << t_vtk << "\n";
	output_file << "filename increment:     " << t_vtk << "\n";
	output_file << "time values:\n";
	for (int itm = t_vtk; itm <= t_num; itm = itm + t_vtk) {
		output_file << itm << "\n";

	}
	/// Close file
	output_file.close();
	/////////// End of Writing CASE file (.case) ///////////////
	////////////////////////////////////////////////////////////


	/////////// Writing GEOMETRY file (.geo) ///////////////////
	////////////////////////////////////////////////////////////

	double min_node_x, min_node_y, min_node_z;
	double max_node_x, max_node_y, max_node_z;
	min_node_x = particle.node[0].x;
	min_node_y = particle.node[0].y;
	min_node_z = particle.node[0].z;
	max_node_x = particle.node[0].x;
	max_node_y = particle.node[0].y;
	max_node_z = particle.node[0].z;


	for (int n = 0; n < particle_num_nodes; ++n) {
		if (particle.node[n].x<min_node_x) {
			min_node_x = particle.node[n].x;
		}
		if (particle.node[n].y<min_node_y) {
			min_node_y = particle.node[n].y;
		}
		if (particle.node[n].z<min_node_z) {
			min_node_z = particle.node[n].z;
		}
		if (particle.node[n].x>max_node_x) {
			max_node_x = particle.node[n].x;
		}
		if (particle.node[n].y>max_node_y) {
			max_node_y = particle.node[n].y;
		}
		if (particle.node[n].z>max_node_z) {
			max_node_z = particle.node[n].z;
		}
	}
	stringstream output_filename1;
	output_filename1 << "Alborz_Results/vtk_particle/particle_" << npp + 1 << "_" << setfill('0') << setw(6) << time << ".geo";//Geometry file
	ofstream output_file1;

	/// Open file
	output_file1.open(output_filename1.str().c_str());

	/// Write VTK header
	output_file1 << "EnSight Model Geometry File\n";
	output_file1 << "EnSight 10.1.1\n";
	output_file1 << "node id given\n";
	output_file1 << "element id given\n";
	output_file1 << "extents\n";
	output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << min_node_x << std::scientific << std::setw(12) << std::setprecision(5) << max_node_x << "\n";
	output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << min_node_y << std::scientific << std::setw(12) << std::setprecision(5) << max_node_y << "\n";
	output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << min_node_z << std::scientific << std::setw(12) << std::setprecision(5) << max_node_z << "\n";
	output_file1 << "part\n";
	output_file1 << std::setw(10) << 1 << "\n";
	output_file1 << "Particle" << npp + 1 << "\n";
	output_file1 << "coordinates\n";
	output_file1 << std::setw(10) << particle_num_nodes << "\n";
	for (int itm = 1; itm <= particle_num_nodes; itm++) {
		output_file1 << std::setw(10) << itm << "\n";
	}

	for (int i = 0; i < particle_num_nodes; i++) {
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << particle.node[i].x << "\n";
	}
	for (int i = 0; i < particle_num_nodes; i++) {
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << particle.node[i].y << "\n";
	}
	for (int i = 0; i < particle_num_nodes; i++) {
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << particle.node[i].z << "\n";
	}


#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
	output_file1 << "tria3\n";
	output_file1 << std::setw(10) << number_triangles << "\n";
	for (int itm = 1; itm <= number_triangles; itm++) {
		output_file1 << std::setw(10) << itm << "\n";
	}

	// Create triangle surfaces
	for (int n = 0; n < number_triangles; ++n) {
		output_file1 << std::setw(10) << triangleset1[n] + 1 << std::setw(10) << triangleset2[n] + 1 << std::setw(10) << triangleset3[n] + 1 << "\n";
	}
#endif


	/// Close file
	output_file1.close();
	/////////// End of Writing GEOMETRY file (.geo) ///////////////////
	///////////////////////////////////////////////////////////////////

#endif//End of CASE file for Ensight
	return;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_SOLID_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// Location of solid nodes specified by is_solid=TRUE is written to a vtk file
void write_solid_vtk(int time, int tot_solid) {

#if defined VTK_BINARY || defined VTK_ASCII

	/// Create filename

	stringstream output_filename;
	output_filename << "Alborz_Results/vtk_particle/solid_nodes_" << 1 << "_t" << time << ".vtk";

	ofstream output_file;

	/// Open file

	output_file.open(output_filename.str().c_str());

	/// Write VTK header

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";

	output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";

	output_file << "X_COORDINATES " << Nx << " float\n";

	for (int X = 0; X < Nx; ++X) {
		output_file << X << " ";
	}

	output_file << "\n";
	output_file << "Y_COORDINATES " << Ny << " float\n";

	for (int Y = 0; Y < Ny; ++Y) {
		output_file << Y << " ";
	}

	output_file << "\n";
	output_file << "Z_COORDINATES " << Nz << " float\n";

	for (int Z = 0; Z < Nz; ++Z) {
		output_file << Z << " ";
	}

	output_file << "POINT_DATA " << Nx * Ny * Nz << "\n";

	/// Write solid_nodes

	output_file << "SCALARS solid_nodes float 1\n";
	output_file << "LOOKUP_TABLE default\n";

	for (int Z = 0; Z < Nz; ++Z) {
		for (int Y = 0; Y < Ny; ++Y) {
			for (int X = 0; X < Nx; ++X) {
				if (is_solid[X][Y][Z] == TRUE) {
					output_file << 0 << "\n";
				}
				else {
					output_file << 1 << "\n";
				}
			}
		}
	}

	/// Close file

	output_file.close();


	////////////////////////////////////
	////////////////////////////////////
	/// Create filename
	stringstream output_filename1;
	output_filename1 << "Alborz_Results/vtk_particle/solid_nodes_" << 2 << "_t" << time << ".vtk";
	ofstream output_file1;

	/// Open file
	output_file1.open(output_filename1.str().c_str());

	/// Write VTK header
	output_file1 << "# vtk DataFile Version 3.0\n";
	output_file1 << "particle_state\n";
	output_file1 << "ASCII\n";
	output_file1 << "DATASET POLYDATA\n";

	/// Write node positions
	output_file1 << "POINTS " << tot_solid << " float\n";
	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			for (int Z = 0; Z < Nz; ++Z) {
				if (is_solid[X][Y][Z] == TRUE) {
					output_file1 << X << " " << Y << " " << Z << " \n";
				}
			}
		}
	}

	/// Write vertices
	output_file1 << "VERTICES 1 " << tot_solid + 1 << "\n";
	output_file1 << tot_solid << " ";
	for (int n = 0; n < tot_solid; ++n) {
		output_file1 << n << " ";
	}

	/// Close file
	output_file1.close();
#endif
	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_DATA
//-----------------------------------------------------------------------------------------------------------------------------
//
// The following PARTICLE quantities are written to the disk at each t_data step:
// - drag and lift forces (x- and y-components of the force)
// - object center position (x- and y-components)
// - object center velocity (x- and y-components)
// The data file is readable by gnuplot
#if defined Flow_With_Particle

void write_particle_data(int time, particle_struct particle, int npp, double force_tot_x, double force_tot_y, double force_tot_z) {

	/// Create filename

    std::string npstr = numberToString(npp + 1);
	//string npstr = static_cast<ostringstream*>(&(ostringstream() << npp + 1))->str();//convert number to string
	string output_filename("Alborz_Results/Data_particle/data_particle_");
	output_filename += npstr;
	output_filename += ".dat";
	ofstream output_file;

	//string output_filename("data.dat");
	//ofstream output_file;


	/// Open file

	output_file.open(output_filename.c_str(), fstream::app);


	/// Compute quantities

	double vel_center_x = 0;
	double vel_center_y = 0;
	double vel_center_z = 0;


	/// Write data

	output_file << time << " "; // time step
	output_file << force_tot_x << " "; // drag force
	output_file << force_tot_y << " "; // lift force
	output_file << force_tot_z << " "; // lift force
	output_file << particle.center[0].x << " "; // center position (x-component)
	output_file << particle.center[0].y << " "; // center position (y-component)
	output_file << particle.center[0].z << " "; // center position (y-component)
	output_file << particle.center[0].vel_x << " "; // center velocity (x-component)
	output_file << particle.center[0].vel_y << " "; // center velocity (y-component)
	output_file << particle.center[0].vel_z << " "; // center velocity (y-component)
	output_file << particle.center[0].omgX << " "; // center velocity (x-component)
	output_file << particle.center[0].omgY << " "; // center velocity (y-component)
	output_file << particle.center[0].omgZ << " "; // center velocity (y-component)
#if defined Flow_Without_Thermal_Effect
	output_file << "\n";
#endif
#if defined Flow_With_Thermal_Effect
	output_file << particle.center[0].temperature << "\n"; // center velocity (y-component)
#endif
														   /// Close file
	output_file.close();

	return;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are written to the disk at each t_recovery step:
void write_fluid_recovery(int time) {

#if defined BINARY_Recovery_File

	std::stringstream str_line1;
	str_line1 << "Alborz_Results/data_fluid_recovery.dat";
	std::string strstr_line1 = str_line1.str();

	std::ofstream output_file1(strstr_line1, std::ios::out | std::ios::binary);

	float tt[1];
	tt[0] = time;
	//std::cout << tt[0];
	output_file1.write((char*)tt, sizeof(float) * 1);

	float mm[1];

	for (int X = 0; X < Nx; X++) {
		for (int Y = 0; Y < Ny; Y++) {
			for (int Z = 0; Z < Nz; Z++) {

				for (int c_i = 0; c_i < 19; ++c_i) {
					mm[0] = pop[c_i][X][Y][Z];
					output_file1.write((char*)mm, sizeof(float) * 1);

				}
				mm[0] = density[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = velocity_x[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = velocity_y[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = velocity_z[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = force_x[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = force_y[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = force_z[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);

#if defined Flow_With_Thermal_Effect
				for (int c_i = 0; c_i < 19; ++c_i) {
					mm[0] = pop_t[c_i][X][Y][Z];
					output_file1.write((char*)mm, sizeof(float) * 1);
				}
				mm[0] = temperature[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
				mm[0] = force_thermal[X][Y][Z];
				output_file1.write((char*)mm, sizeof(float) * 1);
#endif
			}
		}
	}

	output_file1.close();
#endif


#if defined ASCII_Recovery_File

	/// Create filename
	string output_filename("Alborz_Results/data_fluid_recovery");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

    /// Write data
	output_file << time << " "; // time step is written just as the first number in the recovery file.

	for (int X = 0; X<Nx; X++) {
		for (int Y = 0; Y<Ny; Y++) {
			for (int Z = 0; Z<Nz; Z++) {
				for (int c_i = 0; c_i < 19; ++c_i) {
					//	output_file << std::fixed << std::setprecision(15) << pop_old[c_i][X][Y][Z] << " ";
					output_file << std::fixed << std::setprecision(15) << pop[c_i][X][Y][Z] << " ";
				}
				output_file << std::fixed << std::setprecision(15) << density[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << velocity_x[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << velocity_y[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << velocity_z[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << force_x[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << force_y[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << force_z[X][Y][Z] << " ";
#if defined Flow_With_Thermal_Effect
				for (int c_i = 0; c_i < 19; ++c_i) {
					output_file << std::fixed << std::setprecision(15) << pop_t[c_i][X][Y][Z] << " ";
				}
				output_file << std::fixed << std::setprecision(15) << temperature[X][Y][Z] << " ";
				output_file << std::fixed << std::setprecision(15) << force_thermal[X][Y][Z] << "\n";
#endif
			}
		}
	}

	/// Close file
	output_file.close();
#endif

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: READ_FLUID_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are READ from the disk.
void read_fluid_recovery(int & Step) {

#if defined BINARY_Recovery_File
	std::stringstream str_line1;
	str_line1 << "Alborz_Results/data_fluid_recovery.dat";
	std::string strstr_line1 = str_line1.str();

	std::ifstream output_file1(strstr_line1, std::ios::in | std::ios::binary);

	float tt[1];
	output_file1.read((char*)tt, sizeof(float));

	Step = static_cast<int>(tt[0]);
	std::cout << "Recovery step is" << Step << endl;
	if (output_file1.is_open()) {
		for (int X = 0; X < Nx; X++) {
			for (int Y = 0; Y < Ny; Y++) {
				for (int Z = 0; Z < Nz; Z++) {

					for (int c_i = 0; c_i < 19; ++c_i) {
						output_file1.read((char*)tt, sizeof(float));
						pop[c_i][X][Y][Z] = tt[0];
					}
					output_file1.read((char*)tt, sizeof(float));
					density[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					velocity_x[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					velocity_y[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					velocity_z[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					force_x[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					force_y[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					force_z[X][Y][Z] = tt[0];
#if defined Flow_With_Thermal_Effect
					for (int c_i = 0; c_i < 19; ++c_i) {
						output_file1.read((char*)tt, sizeof(float));
						pop_t[c_i][X][Y][Z] = tt[0];
					}
					output_file1.read((char*)tt, sizeof(float));
					temperature[X][Y][Z] = tt[0];
					output_file1.read((char*)tt, sizeof(float));
					force_thermal[X][Y][Z] = tt[0];
#endif
				}
			}
		}
		output_file1.close();
	}

#endif


#if defined ASCII_Recovery_File
	/// Create filename
	string output_filename("Alborz_Results/data_fluid_recovery");
	output_filename += ".dat";
	ifstream output_file;//File is open for READING

						 /// Open file
	output_file.open(output_filename.c_str());

	/// Read data
	output_file >> Step;      // get first number from the file. The first number is: "saved recovery step"

	for (int X = 0; X<Nx; X++) {
		for (int Y = 0; Y<Ny; Y++) {
			for (int Z = 0; Z<Nz; Z++) {
				if (!output_file.eof()) {
					for (int c_i = 0; c_i < 19; ++c_i) {
						//	output_file >> pop_old[c_i][X][Y][Z];
						output_file >> pop[c_i][X][Y][Z];
					}
					output_file >> density[X][Y][Z];
					output_file >> velocity_x[X][Y][Z];
					output_file >> velocity_y[X][Y][Z];
					output_file >> velocity_z[X][Y][Z];
					output_file >> force_x[X][Y][Z];
					output_file >> force_y[X][Y][Z];
					output_file >> force_z[X][Y][Z];
#if defined Flow_With_Thermal_Effect
					for (int c_i = 0; c_i < 19; ++c_i) {
						output_file >> pop_t[c_i][X][Y][Z];
					}
					output_file >> temperature[X][Y][Z];
					output_file >> force_thermal[X][Y][Z];
#endif
				}
			}
		}
	}

	/// Close file
	output_file.close();
#endif
	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_DATA
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are written to the disk at each t_data step:

void write_fluid_data(int time) {

	/// Create filename

	string output_filename("Alborz_Results/data_fluid");
	output_filename += ".dat";
	ofstream output_file;


	/// Open file

	output_file.open(output_filename.c_str(), fstream::app);

	/// Write data

	output_file << time << " "; // time step

								//You can change the location at which you want to save the data
								//for (int X = 0; X<10; X++){
								//	output_file << velocity_x[X][40][20] << " ";
								//}
								//	double tgrad[Ny];
								//	double nusselt1, nusselt2 = 0;
								//	for (int Y = 0; Y<Ny; Y++){
								//		tgrad[Y] = (-3.0*temperature[1][Y][10] + 4.0*temperature[2][Y][10] - temperature[3][Y][10]) / 2.0;
								//		nusselt1 += tgrad[Y];
								//		tgrad[Y] = (temperature[1][Y][10] - temperature[2][Y][10]) / 1.0;
								//		nusselt2 += tgrad[Y];
								//	}
								//	output_file << nusselt1 << " ";
								//	output_file << nusselt2 << "\n";

								//Permeability calculation and storage
	double permeability;
	if (gravity>0) {
		double uav = 0.;
		for (int X = 0; X < Nx; ++X) {
			for (int Y = 0; Y < Ny; ++Y) {
				for (int Z = 0; Z < Nz; ++Z) {
					if (is_solid[X][Y][Z] == FALSE) {
						uav = uav + velocity_x[X][Y][Z]; //Because gravity is just applied in X direction
					}
				}
			}
		}
		uav = uav / (Nx * Ny*Nz);
		permeability = uav * (tau - 0.5) / 3.0 / gravity;
		printf("Flow permeability in Lattice units is %f .\n", permeability);
		output_file << permeability << " ";
	}
	//////////////////////////////////////////

	//output_file << density[10][40][20] << "\n";

	/// Close file

	output_file.close();

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARAMETES
//-----------------------------------------------------------------------------------------------------------------------------
//
// The Main simulation parameters will be written into a file:

void write_parametes() {

	/// Create filename

	string output_filename("Alborz_Results/parameters");
	output_filename += ".dat";
	ofstream output_file;


	/// Open file

	output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

															   /// Write data
#if defined	Flow_With_Particle
	if (num_particles > 0) {

#if defined STATIONARY_SPHERE
		output_file << "STATIONARY_SPHERE " << endl;
#endif
#if defined MOVING_SPHERE
		output_file << "MOVING_SPHERE " << endl;
#endif
#if defined MOVING_SPHEROID
		output_file << "MOVING_SPHEROID " << endl;
#endif
#if defined STATIONARY_SPHEROID
		output_file << "STATIONARY_SPHEROID " << endl;
#endif
	}
#endif
#if defined Shan_Chen
	output_file << "Shan Chen model " << endl;
	output_file << "Single Relaxation Time (SRT)" << endl;
#endif
#if defined Luo
	output_file << "Luo model";
	output_file << "- Single Relaxation Time (SRT)" << endl;
#endif
#if defined MDF
	output_file << " Multi Direct Forcing" << endl;
#endif
#if defined Guo
	output_file << "Guo model ";
#endif
#if defined Guo && defined MRT
	output_file << " - Multi Relaxation Time (MRT)" << endl;
#endif
#if defined Guo && defined SRT
	output_file << " - Single Relaxation Time (SRT)" << endl;
#endif
#if defined X_Periodic
	output_file << "Periodic in X direction " << endl;
#endif
#if defined X_velIn_ZeroGradOut
	output_file << "Constant Velocity Inlet & Zero Gradient Outlet in X direction " << endl;
	output_file << "Inlet velocity  : " << inletvel << endl;
#endif
#if defined ZeroGradIn_ZeroGradOut
	output_file << " Boundary condition:Zero Gradient Inlet & Zero Gradient Outlet in X direction " << endl;
#endif
#if defined X_Wall
	output_file << "Boundary consition in X-Dir: WALL " << endl;
#endif
#if defined Z_Wall
	output_file << "Boundary consition in Z-Dir: WALL " << endl;
#endif
#if defined Z_Periodic
	output_file << "Periodic in Z direction " << endl;
#endif
#if defined Y_Wall
	output_file << "Boundary consition in Y-Dir: WALL " << endl;
#endif
#if defined Y_Periodic
	output_file << "Periodic in Y direction " << endl;
#endif

	output_file << "Nx  : " << Nx << endl;
	output_file << "Ny  : " << Ny << endl;
	output_file << "Nz  : " << Nz << endl;
	output_file << "Number of iterations  : " << t_num << endl;
	output_file << "Wrtie to .dat file (iteration): " << t_data << endl;
	output_file << "Wrtie to .vtk file (iteration): " << t_vtk << endl;
	output_file << "Tau : " << tau << endl;
	output_file << "Lattice Viscosity : " << (2 * tau - 1) / 6 << endl;
	output_file << "Gravity in X direction : " << gravity << endl;
	output_file << "Bottom wall (Y=0)velocity : " << wall_vel_bottom << endl;
	output_file << "Top wall (Y=Ny-1)velocity    : " << wall_vel_top << endl;
	output_file << "Back wall (Z=0)velocity : " << wall_vel_back << endl;
	output_file << "Front wall (Y=Ny-1)velocity    : " << wall_vel_front << endl;
#if defined	Flow_With_Particle
	output_file << "Number of particles : " << num_particles << endl;
	if (num_particles > 0) {
		output_file << "Particle radius     : " << particle_radius << endl;
#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
		output_file << "Particle radius 2   : " << particle_radius_2 << endl;
#endif
		output_file << "1st particle X position : " << particle_center_x << endl;
		output_file << "1st particle Y position : " << particle_center_y << endl;
		output_file << "1st particle Z position : " << particle_center_z << endl;
#if defined Without_added_mass
		output_file << "Without added mass " << endl;
#endif
#if defined With_added_mass
		output_file << "With added mass " << endl;
#endif
		output_file << "Particle gravity in X dir  : " << particle_gravity_x << endl;
		output_file << "Particle gravity in Y dir  : " << particle_gravity_y << endl;
		output_file << "Particle gravity in Z dir  : " << particle_gravity_z << endl;
		output_file << "Particle to Fluid density ratio  : " << den_ratio << endl;
#if defined Two_Point_Delta
		output_file << "Two_Point_Delta   " << endl;
#endif
#if defined Modified_Two_Point_Delta
		output_file << "Modified_Two_Point_Delta   " << endl;
#endif
#if defined Kernel_Delta
		output_file << "Kernel_Delta   " << endl;
#endif
#if defined Four_Point_Delta
		output_file << "Four_Point_Delta   " << endl;
#endif

	}
#endif

	/// Close file

	output_file.close();

	return;
}


// Function to convert a number to string
std::string numberToString(int number) {
    std::ostringstream oss;
    oss << number;
    return oss.str();
}
