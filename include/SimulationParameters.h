#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H

#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <unordered_map>

#include "GlobalFunctions.h"

class SimulationSettings;//Forward declaration

class SimulationParameters  // It is a Singleton class: there can only be one instance of the class created throughout the lifetime of the program.
{
public:
    // Static member function to access the singleton instance
    static SimulationParameters& get_instance(const std::string& filename, const SimulationSettings& settingsLB);

    std::map<std::string, double> sortedParameters;

    void report_parameters_on_screen();

    // Accessor methods for the parameters (Get functions)
    INT32    get_domain_dimension_F() const {return domain_dimension_F;}
    INT32    get_stencil_directions_F() const {return stencil_directions_F;}
    INT32    get_domain_dimension_T() const {return domain_dimension_T;}
    INT32    get_stencil_directions_T() const {return stencil_directions_T;}
    std::string get_stencil_flow() const {return stencil_flow;}
    std::string get_stencil_thermal() const {return stencil_thermal;}
    INT32    get_Nx() const {return Nx;}
    INT32    get_Ny() const {return Ny;}
    INT32    get_Nz() const {return Nz;}
    double get_tau() const {return tau;}
    INT32    get_t_end() const {return t_end;}
    INT32    get_t_data() const {return t_data;}
    INT32    get_t_vtk() const {return t_vtk;}
    INT32    get_t_info() const {return t_info;}
    INT32    get_t_time() const {return t_time;}
    INT32    get_t_recovery() const {return t_recovery;}
    double get_gravity() const {return gravity;}
    double get_wall_vel_bottom() const {return wall_vel_bottom;}
    double get_wall_vel_top() const {return wall_vel_top;}
    double get_wall_vel_back() const {return wall_vel_back;}
    double get_wall_vel_front() const {return wall_vel_front;}
    double get_inletvel() const {return inletvel;}
    INT32    get_mdfstep() const {return mdfstep;}
    double get_omega() const {return omega;}
    double get_fluid_init_temp() const {return fluid_init_temp;}
    double get_gbeta_x() const {return gbeta_x;}
    double get_gbeta_y() const {return gbeta_y;}
    double get_gbeta_z() const {return gbeta_z;}
    double get_Tref() const {return Tref;}
    double get_prandtl() const {return prandtl;}
    double get_alpha_t() const {return alpha_t;}
    double get_tau_t() const {return tau_t;}
    double get_omega_t() const {return omega_t;}
    double get_omega_ts() const {return omega_ts;}
    double get_particle_specific_heat_ratio() const {return particle_specific_heat_ratio;}
    double get_temperature_left() const {return temperature_left;}
    double get_temperature_right() const {return temperature_right;}
    double get_temperature_bot() const {return temperature_bot;}
    double get_temperature_top() const {return temperature_top;}
    double get_temperature_back() const {return temperature_back;}
    double get_temperature_front() const {return temperature_front;}
    double get_heat_flux_left() const {return heat_flux_left;}
    double get_heat_flux_right() const {return heat_flux_right;}
    double get_heat_flux_bottom() const {return heat_flux_bottom;}
    double get_heat_flux_top() const {return heat_flux_top;}
    double get_heat_flux_back() const {return heat_flux_back;}
    double get_heat_flux_front() const {return heat_flux_front;}
    INT32  get_num_particles() const {return num_particles;}
    double get_particle_radius() const {return particle_radius;}
    double get_particle_center_x() const {return particle_center_x;}
    double get_particle_center_y() const {return particle_center_y;}
    double get_particle_center_z() const {return particle_center_z;}
    double get_particle_gravity_x() const {return particle_gravity_x;}
    double get_particle_gravity_y() const {return particle_gravity_y;}
    double get_particle_gravity_z() const {return particle_gravity_z;}
    double get_particle_den_ratio() const {return particle_den_ratio;}
    double get_particle_offset_x() const {return particle_offset_x;}
    double get_particle_offset_y() const {return particle_offset_y;}
    double get_particle_offset_z() const {return particle_offset_z;}
    double get_particle_intial_temperature() const {return particle_temperature;}
    double get_particle_heat_source() const {return particle_heat_source;}
    double get_particle_radius_2() const {return particle_radius_2;}
    INT32  get_particle_num_nodes() const {return particle_num_nodes;}
    INT32  get_particle_num_strip() const {return particle_num_strip;}
    double get_particle_cir_factor() const {return particle_cir_factor;}
    double get_ellipse_init_ang() const {return ellipse_init_ang;}
    double get_particle_r_buffer_verlet() const {return particle_r_buffer_verlet;}
    double get_particle_r_cutoff_verlet() const {return particle_r_cutoff_verlet;}
    double get_particle_cij() const {return particle_cij;}
    double get_particle_ep() const {return particle_ep;}
    double get_particle_Ep() const {return particle_Ep;}
    double get_particle_epw() const {return particle_epw;}
    double get_particle_Epw() const {return particle_Epw;}
    double get_particle_th() const {return particle_th;}
    double get_particle_thw() const {return particle_thw;}
    double get_particle_th_lub() const {return particle_th_lub;}
    double get_particle_thw_lub() const {return particle_thw_lub;}
    INT32  get_initial_num_buffer_layers_X() const {return initial_num_buffer_layers_X;}
    INT32  get_initial_num_buffer_layers_Y() const {return initial_num_buffer_layers_Y;}
    INT32  get_initial_num_buffer_layers_Z() const {return initial_num_buffer_layers_Z;}
    double get_solid_temperature() const {return solid_temperature;}
    double get_solid_to_fluid_conductivity() const {return solid_to_fluid_conductivity;}
    INT32  get_Stop_by_Residual() const {return Stop_by_Residual;}
    double get_residual_flow() const {return residual_flow;}
    double get_residual_thermal() const {return residual_thermal;}
    INT32  get_t_residual() const {return t_residual;}

private:
    /// Private constructor to prevent external instantiation
    SimulationParameters(const std::string& filename, const SimulationSettings&);

    //void set_simulation_settings(const SimulationSettings&);

    /// Destructor
    ~SimulationParameters();

    /// Prevent copy construction and assignment
    SimulationParameters(const SimulationParameters&) = delete;
    void operator=(const SimulationParameters&) = delete;

    void parse_file(const std::string &filename);
    void parse_line(const std::string &line);
    void parse_stencil(const std::string &stencil, std::pair<INT32, INT32> &result);
    template<typename T> void assign_param(const std::string& key, T& var);

    /// Private Parameters Declaration
    ///Settings
    // These const bool variables will be initialized as initializer list in constructor because they are const
    const bool Flow_With_Particle;
    const bool Flow_With_Thermal_Effect;
    const UINT32 Create_Solid;
    const bool Shan_Chen;
    const bool Guo;
    const bool Luo;
    const bool SRT;
    const bool MRT;
    const bool X_Periodic;
    const bool X_Wall;
    const bool X_velIn_ZeroGradOut;
    const bool X_velIn_PressureOut;
    const bool Y_Periodic;
    const bool Y_Wall;
    const bool Z_Periodic;
    const bool Z_Wall;
    const bool VTK_ASCII;
    const bool VTK_BINARY;
    const bool CASE;
    const bool Left_Fixed_Temp;
    const bool Left_Adiabatic;
    const bool Left_Fixed_Heat_Flux;
    const bool Right_Fixed_Temp;
    const bool Right_Adiabatic;
    const bool Right_Fixed_Heat_Flux;
    const bool Bottom_Fixed_Temp;
    const bool Bottom_Adiabatic;
    const bool Bottom_Fixed_Heat_Flux;
    const bool Top_Fixed_Temp;
    const bool Top_Adiabatic;
    const bool Top_Fixed_Heat_Flux;
    const bool Back_Fixed_Temp;
    const bool Back_Adiabatic;
    const bool Back_Fixed_Heat_Flux;
    const bool Front_Fixed_Temp;
    const bool Front_Adiabatic;
    const bool Front_Fixed_Heat_Flux;
    const bool Force_Thermal_1;
    const bool Force_Thermal_2;
    const bool Equilibrium_Thermal_1;
    const bool Equilibrium_Thermal_2;
    const bool STATIONARY_SPHERE;
    const bool STATIONARY_SPHEROID;
    const bool MOVING_SPHERE;
    const bool MOVING_SPHEROID;
    const bool SDF;
    const bool MDF;
    const bool Without_added_mass;
    const bool With_added_mass;
    const bool Corrected_Radius;
    const bool Two_Point_Delta;
    const bool Four_Point_Delta;
    const bool Kernel_Delta;
    const bool Modified_Two_Point_Delta;
    const bool Sphere_Collision_1_Feng;
    const bool Sphere_Collision_2_Feng_Verlet;
    const bool Sphere_Collision_3_Feng_Lub_Verlet;
    const bool Spheroid_Collision_1_Feng;
    const bool Spheroid_Collision_2_Feng_Verlet;
    const bool Spheroid_Collision_3_Feng_Lub_Verlet;
    const bool Particle_Fixed_Temperature;
    const bool Particle_Varying_Temperature;
    const bool Random_Arrangement_OFF;
    const bool Random_Arrangement_ON;
    const bool Recovery_OFF;
    const bool Recovery_ON;
    const bool ASCII_Recovery_File;
    const bool BINARY_Recovery_File;

    ///2D Settings
    const bool STATIONARY_CYLINDER;
    const bool STATIONARY_ELLIPSE;
    const bool MOVING_CYLINDER;
    const bool MOVING_ELLIPSE;
    const bool Circle_Collision_1;   // R. Glowinski et al, J. Comput. Phys., 169 (2001) 363-426.///Z.-G. Feng, E.E. Michaelides,  806(Journal of Computational Physics, 195, (2004), pp.602�628)
    const bool Circle_Collision_2;   // Wan & Turek, Int. J. Numer. Meth. Fluids, 51 (2006), 531-566./// Liang Wang et al., 814(Commun. Comput. Phys., 13 (2013), pp.1151-1172)
    const bool Circle_Collision_3;   //usually used. Feng ZG, Michaelides EE. Journal of Computational Physics,202 (2005), pp.20-51.
    const bool Ellipse_Collision_1;  // R. Glowinski et al, J. Comput. Phys., 169 (2001) 363�426.///Z.-G. Feng, E.E. Michaelides,  806(Journal of Computational Physics, 195, (2004), pp.602�628)
    const bool Ellipse_Collision_2;
    const bool Ellipse_Collision_3;
    const bool Ellipse_Collision_4;  //usually used.

    ///Variables
    std::unordered_map<std::string, double> parameters;
    std::unordered_map<std::string, double> parameters1;
    std::vector<std::string> insertion_order;
    std::pair<INT32, INT32> Flow_Stencil;
    std::pair<INT32, INT32> Thermal_Stencil;

    ///80 parameters = (70 read from parametersLB.txt) + (2 stencils - std::string) + (Calculated: omega, alpha_t, tau_t, omega_t) + (domain_dimension_F, stencil_directions_F, domain_dimension_T, stencil_directions_T)
    std::string stencil_flow;
    std::string stencil_thermal;
    INT32 Nx{}; // number of lattice nodes along the x-axis (periodic)
    INT32 Ny{}; // number of lattice nodes along the y-axis (including two wall nodes)
    INT32 Nz{}; // number of lattice nodes along the z-axis (including two wall nodes)
    double tau{}; // relaxation time
    INT32 t_end{}; // number of time steps (running from 1 to t_num)
    INT32 t_data{}; // disk write time step (data will be written to the disk every t_data step)
    INT32 t_vtk{}; // VTK write time step
    INT32 t_info{}; // info time step (screen message will be printed every t_info step)
    INT32 t_time{}; // Report elapsed time on screen
    INT32 t_recovery{}; // Time to store in recovery file
    double gravity{};// force density due to gravity (in positive x-direction). Just for flow without particle
    double wall_vel_bottom{}; //velocity of the bottom wall(Y=0) (in positive x-direction)
    double wall_vel_top{}; // velocity of the top wall(Y=Ny-1) (in positive x-direction)
    double wall_vel_back{};//velocity of the back wall(Z=0) (in positive x-direction)
    double wall_vel_front{};//velocity of the front wall(Z=Nz-1) (in positive x-direction)
    double inletvel{};//Inlet velocity. Just used in case of #define X_velIn_ZeroGradOut
    INT32  mdfstep{};//Number of MDF steps. Adjust this line
    double fluid_init_temp{};// Initial fluid temperature
    double gbeta_x{}; //In NATURAL convection: Ra=Gr*Pr=gbeta*dT*H^3/(alpha_t*nu)==>gbeta=alpha_t*nu*Ra/(dT*H^3). If natural convection is not important, set gbeta=0.0
    double gbeta_y{};
    double gbeta_z{};
    double Tref{}; //infinity temperature in case of natural convection
    double prandtl{}; //fluid Prandtl number
    double particle_specific_heat_ratio{}; //Particle to fluid specific heat ratio.
    double temperature_left{};
    double temperature_right{};
    double temperature_bot{};
    double temperature_top{};
    double temperature_back{};
    double temperature_front{};
    double heat_flux_left{};
    double heat_flux_right{};
    double heat_flux_bottom{};
    double heat_flux_top{};
    double heat_flux_back{};
    double heat_flux_front{};
    double solid_temperature{};
    double solid_to_fluid_conductivity{};
    INT32  Stop_by_Residual{};
    double residual_flow{};
    double residual_thermal{};
    INT32  t_residual{};
    /// PARTICLE PROPERTIES
    INT32 num_particles{}; // number of particles.
    double particle_radius{}; // particle radius
    double particle_center_x{}; // center position (x-component) of first particle
    double particle_center_y{}; // center position (y-component) of first particle
    double particle_center_z{}; // center position (z-component) of first particle
    double particle_gravity_x{};//This force is applied to particles according to density ratio
    double particle_gravity_y{};//This force is applied to particles according to density ratio
    double particle_gravity_z{};//This force is applied to particles according to density ratio
    double particle_den_ratio{}; //Particle to fluid density ratio. Note that in below simulations density of fluid is assumed to be 1.0
    double particle_offset_x{}; //Distance of particles from each other in x direction (For Random_Arrangement_OFF).
    double particle_offset_y{}; //Distance of particles from each other in y direction (For Random_Arrangement_OFF).
    double particle_offset_z{}; //Distance of particles from each other in z direction (For Random_Arrangement_OFF).
    double particle_temperature{}; //Particle temperature must be constant
    double particle_heat_source{}; //Heat source inside the particle
    INT32 particle_num_nodes{};// For 2D particle comes from the ParametersLB.txt (for 3D calculated)
    double particle_radius_2{}; //Just used for "SPHEROID": radius along the Y-axis and Z-axis. The equation of the Spheroid is: (X^2/particle_radius^2) + (Y^2/particle_radius_2^2) + (Z^2/particle_radius_2^2) = 1. must NOT be exactly equal to particle_radius
    INT32 particle_num_strip{}; // Manual adjustment of total number of strips.
    double particle_cir_factor{}; //This parameter adjusts the number of points on each strip. Recommended: 0.6-2.0
    double ellipse_init_ang{};//Initial angle of ELLIPSE with respect to X-dir in anti-clock-wise direction
    ///Particle collision parameters
    double particle_r_buffer_verlet{};//Preferably = (particle_r_cutoff + max(particle_radius,particle_radius_2))
    double particle_r_cutoff_verlet{};//Preferably >= [2 * max(particle_radius,particle_radius_2) + collision_threshold]
    double particle_cij{};
    double particle_ep{};//stiffness factor for particle-particle collision
    double particle_Ep{};//stiffness factor for particle-particle collision usually smaller than ep
    double particle_epw{};//stiffness factor for particle-wall collision
    double particle_Epw{};//stiffness factor for particle-wall collision usually smaller than epw
    double particle_th{};// Threshold for particle-particle collision
    double particle_thw{};// Threshold for particle-wall collision
    double particle_th_lub{};// Lubrication force threshold (larger than th)
    double particle_thw_lub{};// Lubrication force threshold (larger than thw)
    ///MPI parallel parameters
    INT32 initial_num_buffer_layers_X{};
    INT32 initial_num_buffer_layers_Y{};
    INT32 initial_num_buffer_layers_Z{};

    ///Calculated parameters
    double omega{}; // relaxation frequency (inverse of relaxation time)
    double alpha_t{}; // Thermal diffusion coefficient
    double tau_t{}; // relaxation time of temperature
    double omega_t{}; // relaxation frequency (inverse of relaxation time for temperature)
    double omega_ts{}; // relaxation frequency (inverse of relaxation time for temperature based on thermal diffusivity of solid nodes)
    INT32 domain_dimension_F{};//Domain dimensions
    INT32 stencil_directions_F{};// Lattice stencils for Flow (e.g. 19)
    INT32 domain_dimension_T{};//Domain dimensions
    INT32 stencil_directions_T{};// Lattice stencils for Flow (e.g. 19)
    double nu{};//kinematic viscosity
};

#endif // SIMULATIONPARAMETERS_H
