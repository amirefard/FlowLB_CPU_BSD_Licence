#ifndef FLOWSOLVER_H
#define FLOWSOLVER_H

#include <vector>   // vector containers
#include <sys/stat.h> //for mkdir
#include <iostream>
#if defined(_WIN32)
    #include <direct.h>   // Windows: _mkdir
#else
    #include <sys/stat.h> // Linux/macOS: mkdir
    #include <sys/types.h>
#endif

//#ifdef __linux__
//#include <sys/stat.h>//for mkdir
//#else
//#include "direct.h"//for mkdir
//#endif
#include <filesystem> // Only C++17, for std::filesystem
#include <fstream> // Include the header for file I/O operations
#include <iomanip> // for setprecision
#include <stdexcept>
using namespace std;//needed for cin, cout, vector

#include "GlobalFunctions.h"

// Forward declaration: tells the compiler “this class exists, but I don’t need to know its details yet. Needed when two classes use eachother.
class ParticleResolved;
class SimulationParameters;
class SimulationSettings;
class ParallelMPI;

//=============================================================================================================================
enum class FaceDirection
{
    PlusX,
    MinusX,
    PlusY,
    MinusY,
    PlusZ,
    MinusZ
};
//=============================================================================================================================
class FlowSolver
{
public:
    //FlowSolver();
    FlowSolver(const SimulationParameters&, const SimulationSettings&, const ParallelMPI&);// constructor
    FlowSolver(const FlowSolver& other);//copy constructor
    virtual ~FlowSolver();

    ///Methods
    //void set_simulation_settings(const SimulationSettings&);
    void set_simulation_parameters();
    void stensil_setup();
    void initialize_flow();// allocate memory and initialize variables
    double equilibrium_flow(double den, double vel_x, double vel_y, double vel_z, INT32 i);// compute the equilibrium populations from the fluid density and velocity
    void LBM_flow(INT32 time, const ParallelMPI&, const ParticleResolved* = nullptr);// perform LBM operations. Important: I used pointer for ParticleResolved because I want to set it nullptr for !Flow_With_Particle
    void exchange_flow_data_neighbors(const std::string&);
    void gather_fluid_macroscopic_data_to_master();
    void gather_fluid_macroscopic_data_to_master_optimized();
    void scatter_force_modified_indices(const ParallelMPI& paralLB);
    void scatter_fluid_force_from_master();
    void scatter_fluid_force_from_master_optimized(const ParallelMPI& paralLB);
    void bc_at_walls_flow(const ParallelMPI&);// BC on the walls
    void bc_at_solid_points_flow();// BC on solid points
    void momenta_flow();  // compute fluid density and velocity from the populations
    double meq_MRT(double, double, double, double, INT32);// Equilibrium momenta
    double equilibrium_thermal(double tem, double den, double vel_x, double vel_y, double vel_z, INT32 i);
    void LBM_thermal(INT32 time);// perform LBM operations for temperature
    void bc_at_walls_thermal(const ParallelMPI&);
    void bc_at_solid_points_thermal();// Thermal BC on solid points
    void momenta_thermal();// compute fluid temperature from the populations
    void write_fluid_vtk(INT32 time);// write the fluid state to the disk as VTK file
    void write_solid_vtk(INT32 time, INT32 tot_solid);// write SOLID points to a VTK file.
    void write_fluid_data(INT32); // write fluid data to the disk (velocity/density)
    void gather_fluid_recovery_data_to_master(); // Send pop_F, pop_T, etc to MASTER
    void compress_file(const std::string&);
    bool decompress_file(const std::string&, const std::string&);
    void write_fluid_recovery(INT32);// write fluid data to the disk for case of recovery if program crashes
    void read_fluid_recovery(INT32&);// read fluid initial values from a file if "Recovery_ON" is defined
    void read_fluid_initialize_data(INT32&, const INT32, const INT32, const INT32);// read fluid initial values from a file if "Read_Initialize_Data_From_File" is defined
    void scatter_fluid_recovery_data_from_master(INT32&);// share fluid recovery data from the master to other processors
    void Guo_part_two(const ParallelMPI&, const ParticleResolved* = nullptr);// Important: I used a pointer for ParticleResolved because I want to set it nullptr for !Flow_With_Particle
    void force_thermal_2_part_two();
    void set_Eulerian_forces_in_master_zero();
    void apply_sponge_layer_on_outlet_X(INT32);
    INT32 idx_3d_vector(INT32 x, INT32 y, INT32 z) const;       // index for 3d vectors like temperature, density whih are a function of X, Y, Z
    INT32 idx_4d_vector(INT32 q, INT32 id3) const;
    INT32 idx_4d_vector(INT32 q, INT32 x, INT32 y, INT32 z) const;// index for 4d vectors like pop, velocity. Function of X, Y, Z and 1 other parameter
    void decode_idx_3d(INT32 id3, INT32 &x, INT32 &y, INT32 &z) const;
    void read_solid_from_file(const std::string& filename, const ParallelMPI&);
    void create_solid_from_dimensions(const ParallelMPI& parallelLB);
    INT32 count_solid(const ParallelMPI&);
    void scatter_is_solid_points_from_master();
    bool check_residual();

    #if STORE_1D_VECTOR
    std::vector<double> density;
    std::vector<bool> is_solid;
    std::vector<INT16> solid_thermal;//Thermal BC for solid points
    std::vector<double> temperature;
    std::vector<double> previous_temperature;
    std::vector<double> force_thermal;// force term for thermal lbm
    std::vector<double> velocity;
    std::vector<double> previous_velocity;
    std::vector<double> force;
    std::vector<double> mdfVelocity;
    std::vector<double> mdfTemperature;
    std::vector<double> mdfForce;//fluid force for MDF
    std::vector<bool> cell_marked;
    std::vector<INT32> force_modified_indices;   // flat indices for 1D mode
    std::vector<std::vector<INT32>> proc_indices;
    #endif

    #if !STORE_1D_VECTOR
    vector<vector<vector<double>>> density;
    vector<vector<vector<bool>>> is_solid;
    vector<vector<vector<INT16>>> solid_thermal;
    vector<vector<vector<double>>> temperature;
    vector<vector<vector<double>>> previous_temperature;
    vector<vector<vector<double>>> force_thermal; // force term for thermal lbm
    vector<vector<vector<vector<double>>>> velocity;
    vector<vector<vector<vector<double>>>> previous_velocity;
    vector<vector<vector<vector<double>>>> force;
    vector<vector<vector<vector<double>>>> mdfVelocity;
    vector<vector<vector<vector<double>>>> mdfForce;//fluid force for MDF
    vector<vector<vector<double>>> mdfTemperature;
    vector<vector<vector<bool>>> cell_marked;
    std::vector<std::tuple<INT32,INT32,INT32>> force_modified_indices;   // flat indices for 1D mode
    std::vector<std::vector<std::tuple<INT32,INT32,INT32>>> proc_indices;
    #endif

    #if STORE_1D_VECTOR
    std::vector<double> pop_F;// Flow distribution function
    std::vector<double> pop_old_F;
    std::vector<double> pop_T;// Temperature distribution function
    std::vector<double> pop_old_T;
    #endif // STORE_1D_VECTOR

    #if !STORE_1D_VECTOR
    std::vector<vector<vector<vector<double>>>> pop_F; // Flow distribution function
    std::vector<vector<vector<vector<double>>>> pop_old_F;
    std::vector<vector<vector<vector<double>>>> pop_T;// Temperature distribution function
    std::vector<vector<vector<vector<double>>>> pop_old_T;
    #endif
//===========================================================================//
protected:

//===========================================================================//
private:
    double cs_2_Flow;// means c_s^2 for flow lattice
    double cs_2_Thermal;// means c_s^2 for thermal lattice
    vector<INT32> cx_pop_F, cy_pop_F, cz_pop_F;
    vector<INT32> cx_pop_T, cy_pop_T, cz_pop_T;
    vector<double> weight_F;//Flow weights
    vector<double> weight_T;//Thermal weights
       // const double weight[19] = { 1. / 3., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36. }; // lattice weights
    vector<INT32> bounce_back_pop_F;//bounce back direction of each direction (for flow)
    vector<INT32> bounce_back_pop_T;//bounce back direction of each direction (for thermal)
    vector<INT32> pop_to_top_F, pop_to_bottom_F, pop_to_left_F, pop_to_right_F, pop_to_front_F, pop_to_back_F;
    vector<INT32> pop_to_top_T, pop_to_bottom_T, pop_to_left_T, pop_to_right_T, pop_to_front_T, pop_to_back_T;
    vector<INT32> pop_reverse_top_F, pop_reverse_bottom_F, pop_reverse_left_F, pop_reverse_right_F, pop_reverse_front_F, pop_reverse_back_F;
    vector<INT32> pop_reverse_top_T, pop_reverse_bottom_T, pop_reverse_left_T, pop_reverse_right_T, pop_reverse_front_T, pop_reverse_back_T;

    double ssf_MRT{};
    std::vector<double> s_MRT;
    std::vector<double> fforce_MRT;
    std::vector<double> mforce_MRT;
    std::vector<double> m_MRT;
    std::vector<double> f_post_MRT;
    double M_D3Q19[19][19];//Transformation matrix for D3Q19
    double Minv_D3Q19[19][19];
	double M_D2Q9[9][9];//Transformation matrix for D2Q9
    double Minv_D2Q9[9][9];

    ///Settings
    const bool Flow_With_Particle;
    const bool Flow_With_Thermal_Effect;
    const bool X_Periodic;
    const bool Y_Periodic;
    const bool Z_Periodic;
    const bool Recovery_OFF;
    const bool Recovery_ON;
    const bool Shan_Chen;
    const bool Guo;
    const bool Luo;
    const bool SRT;
    const bool MRT;
    const bool X_Wall;
    const bool Y_Wall;
    const bool Z_Wall;
    const bool X_velIn_ZeroGradOut;
    const bool X_velIn_PressureOut;
    const bool X_PressureIn_PressureOut;
    const bool Parabolic_Velocity_Inlet;
    const bool Constant_Velocity_Inlet;
    const bool Equilibrium_Thermal_1;
    const bool Equilibrium_Thermal_2;
    const bool Force_Thermal_1;
    const bool Force_Thermal_2;
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
    const bool VTK_ASCII;
    const bool VTK_BINARY;
    const bool CASE;
    const bool ASCII_Recovery_File;
    const bool BINARY_Recovery_File;
    const INT16 Solid_Thermal_BC;


    ///Required parameters
    const std::string stencil_flow;
    const std::string stencil_thermal;
    const INT32 Nx;
    const INT32 Ny;
    const INT32 Nz;
    const double tau;
    const INT32 t_end;
    const INT32 t_data;
    const INT32 t_vtk;
    const INT32 t_info;
    const INT32 t_time;
    const INT32 t_recovery;
    const bool Stop_by_Residual;
    const double residual_flow;
    const double residual_thermal;
    const INT32 t_residual;
    const double gravity;
    const double wall_vel_bottom;
    const double wall_vel_top;
    const double wall_vel_back;
    const double wall_vel_front;
    const double inletvel;
    const INT32 mdfstep;
    const double fluid_init_temp;
    const double gbeta_x;
    const double gbeta_y;
    const double gbeta_z;
    const double Tref;
    const double prandtl;
    const double particle_specific_heat_ratio;
    const double temperature_left;
    const double temperature_right;
    const double temperature_bot;
    const double temperature_top;
    const double temperature_back;
    const double temperature_front;
    const double heat_flux_left;
    const double heat_flux_right;
    const double heat_flux_bottom;
    const double heat_flux_top;
    const double heat_flux_back;
    const double heat_flux_front;
    const INT32 domain_dimension_F;
    const INT32 stencil_directions_F;
    const INT32 domain_dimension_T;
    const INT32 stencil_directions_T;
    const double omega;
    const double alpha_t;
    const double tau_t;
    const double omega_t;
    const INT32 num_particles;
    const double solid_temperature;
    const double solid_to_fluid_conductivity;
    const double omega_ts;

    ///Parameters (params from ParallelMPI class)
    const INT32 this_processor_id;
    const INT32 num_processors_total;
    const std::vector<INT32> start_global_x;
    const std::vector<INT32> end_global_x;
    const std::vector<INT32> start_global_y;
    const std::vector<INT32> end_global_y;
    const std::vector<INT32> start_global_z;
    const std::vector<INT32> end_global_z;
    const std::vector<INT32> start_local_x;
    const std::vector<INT32> end_local_x;
    const std::vector<INT32> start_local_y;
    const std::vector<INT32> end_local_y;
    const std::vector<INT32> start_local_z;
    const std::vector<INT32> end_local_z;
    const std::vector<std::vector<INT32>> local_neighbors_left;
    const std::vector<std::vector<INT32>> local_neighbors_right;
    const std::vector<std::vector<INT32>> local_neighbors_bottom;
    const std::vector<std::vector<INT32>> local_neighbors_top;
    const std::vector<std::vector<INT32>> local_neighbors_front;
    const std::vector<std::vector<INT32>> local_neighbors_back;
    const std::vector<INT32> length_with_buffer_X;
    const std::vector<INT32> length_with_buffer_Y;
    const std::vector<INT32> length_with_buffer_Z;
    const std::vector<INT32>  num_buffer_layers_X;
    const std::vector<INT32>  num_buffer_layers_Y;
    const std::vector<INT32>  num_buffer_layers_Z;

    double numerator_residual_flow{0};
    double denominator_residual_flow{0};

    double numerator_residual_thermal{0};
    double denominator_residual_thermal{0};

    INT32 limitBound_left_X{}, limitBound_right_X{}, limitBound_bottom_Y{}, limitBound_top_Y{}, limitBound_back_Z{}, limitBound_front_Z{};
    INT32 total_solid{0};//total number of is_solid
    void get_do_LBM_thermal_omega_thermal(bool&, double&, INT32 X, INT32 Y, INT32 Z);
    void Apply_MovingWall_BC(const std::vector<INT32>& pop_to, const std::vector<INT32>& pop_reverse, INT32 X, INT32 Y, INT32 Z, double wall_ux, double wall_uy, double wall_uz);
    void Apply_Velocity_Inlet_Left_ZouHe_D2Q9(const ParallelMPI& paralLB, INT32 X, INT32 Y, INT32 Z, double ux, double uy);
    void Apply_ZeroGradient_Right_BC(INT32 X, INT32 Y, INT32 Z, const std::vector<INT32>& pop_reverse, const std::vector<INT32>& pop_to, INT32 normalX, INT32 normalY, INT32 normalZ);
    void Apply_NonEquilibrium_Extrapolation_Right_BC(INT32 X, INT32 Y, INT32 Z, const std::vector<INT32>& pop_reverse, const std::vector<INT32>& pop_to, INT32 normalX, INT32 normalY, INT32 normalZ);
    void Apply_Convective_Outlet_Right_BC(INT32 X, INT32 Y, INT32 Z, const std::vector<INT32>& pop_reverse, const std::vector<INT32>& pop_to, INT32 normalX, INT32 normalY, INT32 normalZ);
    void Apply_Pressure_Inlet_Left_Ginzburg(const std::vector<INT32>& pop_to, const std::vector<INT32>& pop_reverse, INT32 X, INT32 Y, INT32 Z, double rho_in);
    void Apply_Pressure_Outlet_Right_Ginzburg(const std::vector<INT32>& pop_to, const std::vector<INT32>& pop_reverse, INT32 X, INT32 Y, INT32 Z, INT32 normalX, INT32 normalY, INT32 normalZ, double rho_out);
    void Apply_Pressure_Outlet_Right_ZouHe_D2Q9(INT32 X, INT32 Y, INT32 Z, double rho_out);
    void Apply_Pressure_Inlet_Left_ZouHe_D2Q9(INT32 X, INT32 Y, INT32 Z, double rho_in);
    void Apply_Velocity_Inlet_Left_Ginzburg(const ParallelMPI& paralLB, const std::vector<INT32>& pop_to, const std::vector<INT32>& pop_reverse, INT32 X, INT32 Y, INT32 Z, double ux, double uy, double uz);
    void Process_Surface_Data(INT32 x_local_start, INT32 x_local_end, INT32 y_local_start, INT32 y_local_end, INT32 z_local_start,INT32 z_local_end, FaceDirection face_dir, std::vector<double>& theta_values, std::vector<double>& nusselt_values);

    void TransformPopulationToMomentSpace(const std::vector<double>& pop_old_F, std::vector<double>& m_MRT, const double* M, INT32 Q, INT32 X, INT32 Y, INT32 Z);
    void TransformForceToMomentSpace(const std::vector<double>& f, std::vector<double>& m, const double* M, int Q);
    void TransformToPopulationSpace(const std::vector<double>& m, std::vector<double>& f, const double* Minv, const INT32 Q);
    void InvertMatrix(const double* A, double* Ainv, int N);
    void StreamPostCollisionPopulations(const std::vector<double>& f_post, std::vector<double>&  pop_F, INT32 X, INT32 Y, INT32 Z);
    //double heat_flux_at_bc_over_lambda = 0.1; //q/k for the fixed flux BC, where q is the heat flux at the boundary and k is the conductivity of the fluid.
};

#endif // FLOWSOLVER_H
