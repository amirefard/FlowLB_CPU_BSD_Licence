#ifndef PARTICLERESOLVED_H
#define PARTICLERESOLVED_H

#include <iostream>
#include <vector>
#include <cmath>    // mathematical library (M_PI)
#include <sstream> // Include the <sstream> header for std::stringstream
#include <fstream> // Include the <fstream> header for std::ofstream
#include <mpi.h>
#include <algorithm> // For std::max,  std::min

#include "GlobalFunctions.h"

// Forward declaration: tells the compiler “this class exists, but I don’t need to know its details yet. Needed when two classes use eachother.
class FlowSolver;
class SimulationParameters;
class SimulationSettings;
class ParallelMPI;

//=============================================================================================================================
class ParticleResolved
{
public:
    //ParticleResolved();
    ParticleResolved(const SimulationParameters&, const SimulationSettings&, const ParallelMPI&);    //ctor
    ParticleResolved(const ParticleResolved& other);//copy constructor
    virtual ~ParticleResolved();

    ///Public Methods
    void initialize_particles();
    //void set_simulation_settings(const SimulationSettings&);
    void set_simulation_parameters(const SimulationParameters&);
    void direct_forcing(INT32, FlowSolver&, const ParallelMPI&);
    void newton_eqs_motion(INT32, FlowSolver&);
    void write_particle_vtk(INT32); // write the particle state to the disk as VTK file
    void write_particle_data(INT32); // write particle data to the disk (drag/lift, center position)
    void write_particle_recovery(INT32);
    void read_particle_recovery();
    void calculate_force_modified_indices(FlowSolver&, const ParallelMPI&);
    INT32  calculate_gbetaindex(INT32 X, INT32 Y, INT32 Z, const ParallelMPI&) const;
    void scatter_particles_center_data_from_MASTER();

private:

    INT32 num_mpi_processors{};//number of mpi processors
    INT32 mpi_processor_rank{};//rank of the mpi process
    std::vector<INT32> num_particles_per_mpi_processor;

    /// Settings
    const bool Flow_With_Thermal_Effect;
    const bool X_Periodic;
    const bool X_Wall;
    const bool Y_Periodic;
    const bool Y_Wall;
    const bool Z_Periodic;
    const bool Z_Wall;
    const bool VTK_ASCII;
    const bool VTK_BINARY;
    const bool CASE;
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

    ///2D Settings
    const bool STATIONARY_CYLINDER;
    const bool STATIONARY_ELLIPSE;
    const bool MOVING_CYLINDER;
    const bool MOVING_ELLIPSE;
    const bool Circle_Collision_1;
    const bool Circle_Collision_2;
    const bool Circle_Collision_3;
    const bool Ellipse_Collision_1;
    const bool Ellipse_Collision_2;
    const bool Ellipse_Collision_3;
    const bool Ellipse_Collision_4;

    ///Required parameters
    const INT32 t_end;
    const INT32 t_vtk;
    const INT32 Nx;
    const INT32 Ny;
    const INT32 Nz;
    const double tau;
    const INT32 mdfstep;
    const INT32 num_particles;
    double particle_radius; // Important: This parameter is not const becuase this class modifies it
    const double particle_center_x;
    const double particle_center_y;
    const double particle_center_z;
    const double particle_gravity_x;
    const double particle_gravity_y;
    const double particle_gravity_z;
    const double particle_den_ratio;
    const double particle_offset_x;
    const double particle_offset_y;
    const double particle_offset_z;
    const double particle_temperature;
    const double particle_heat_source;
    const double particle_specific_heat_ratio;
    double particle_radius_2;// Important: This parameter is not const becuase Particle class modifies it
    INT32 particle_num_nodes;// Important: This parameter is not const becuase Particle class modifies it
    const INT32 particle_num_strip;
    const double particle_cir_factor;
    const double ellipse_init_ang;
    double particle_r_buffer_verlet;// Important: This parameter is not const becuase Particle class modifies it
    double particle_r_cutoff_verlet;// Important: This parameter is not const becuase Particle class modifies it
    const double particle_cij;
    const double particle_ep;
    const double particle_Ep;
    const double particle_epw;
    const double particle_Epw;
    const double particle_th;
    const double particle_thw;
    const double particle_th_lub;
    const double particle_thw_lub;
    const INT32 num_processors_for_particles_only;
    const std::vector<INT32> length_with_buffer_X;
    const std::vector<INT32> length_with_buffer_Y;
    const std::vector<INT32> length_with_buffer_Z;
    const INT32 stencil_directions_F;
    const INT32 stencil_directions_T;

    ///Parameters (params from ParallelMPI class)
    const INT32 this_processor_id;
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
    const std::vector<INT32>  num_buffer_layers_X;
    const std::vector<INT32>  num_buffer_layers_Y;
    const std::vector<INT32>  num_buffer_layers_Z;

    ///Only created and used in ParticleResolved class:
    std::vector <INT32> triangleset1;//Total number of triangles is particle_num_nodesx4. Therefore, it must be less than 24000
    std::vector <INT32> triangleset2;
    std::vector <INT32> triangleset3;
    INT32 number_triangles{};//number of triangles created for vtk file of particle
    double corrected_particle_radius{}, corrected_particle_radius_2{};
    double uncorrected_particle_radius{}, uncorrected_particle_radius_2{};
    double all_particles_disp_max_1{};// Max Displacement of particles
    double all_particles_disp_max_2{};// Second max of Displacement of particles
    INT32 margin{};//This parameter is helpful for STATIONARY_SPHERE and STATIONARY_SPHEROID. You can change this line if you want the particles to be stationary until a certain iteration (margin) and the starts motion.


    std::vector<std::vector<INT32>> verlet_list;

    const double Kernel_Max_Dist = 2.0;
    const double Kernel_Resolution = 0.0001;
    const INT32 Kernel_Table_Size = static_cast<INT32>(Kernel_Max_Dist / Kernel_Resolution) + 1;
    std::vector<double> kernel_table;

    double total_displacement_x{0};   //Just used for keeping particle in center
    double partial_displacement_x{0}; //Just used for keeping particle in center

       struct node_struct
    {
        // Elements
        double x{};             // current x-position
        double y{};             // current y-position
        double z{};             // current z-position
        double vel_x{};         // node velocity (x-component)
        double vel_y{};         // node velocity (y-component)
        double vel_z{};         // node velocity (z-component)
        double force_x{};       // node force (x-component)
        double force_y{};       // node force (y-component)
        double force_z{};       // node force (z-component)
        double total_node_mdf_force_x{}; // total node force on each node after mdf (x-component)
        double total_node_mdf_force_y{}; // total node force on each node after mdf (y-component)
        double total_node_mdf_force_z{}; // total node force on each node after mdf (y-component)
        double den{};           // Density at boundary point
        double force_thermal{}; // node force due to thermal effect (it is not in Newton)
        double total_node_mdf_force_thermal; // total node force due to thermal effect at each node after all mdf steps (it is not in Newton)
        double temperature{};   // node force due to thermal effect (it is not in Newton)
        double x_pos_1{};       // Node x-location in BODY-fixed coordinate with respect to particle center. It does not change during time steps
        double y_pos_1{};       // Node y-location in BODY-fixed coordinate with respect to particle center. It does not change during time steps
        double z_pos_1{};       // Node z-location in BODY-fixed coordinate with respect to particle center. It does not change during time steps
        double area_node{};     // area belonging to each surface node.
    };
    struct center_struct
    {
        // Elements
        double x{};                         // current x-position
        double y{};                         // current y-position
        double z{};                         // current z-position
        double vel_x{};                     // center velocity (x-component)
        double vel_y{};                     // center velocity (y-component)
        double vel_z{};                     // center velocity (z-component)
        double force_x{};                   // center force (x-component)
        double force_y{};                   // center force (y-component)
        double force_z{};                   // center force (z-component)
        double den{};                       // Density at boundary point
        double omgX{};                      // center rotational speed
        double omgY{};                      // center rotational speed
        double omgZ{};                      // center rotational speed
        double temperature{};               // node force due to thermal effect (it is not in Newton)
        double force_thermal{};            // center force due to thermal effect (it is not in Newton)
        double previous_vel_x{};           // Velocity x in previous time step
        double previous_vel_y{};           // Velocity y in previous time step
        double previous_vel_z{};           // Velocity z in previous time step
        double tetaX{};                    // Angle x
        double tetaY{};                    // Angle y
        double tetaZ{};                    // Angle z
        double tempx[2] {};               // Translational speed of particle at two previous time steps
        double tempy[2] {};               // Translational speed of particle at two previous time steps
        double tempz[2] {};               // Translational speed of particle at two previous time steps
        double temp_omgX[2] {};          // Rotational speed of particle at two previous time steps
        double temp_omgY[2] {};          // Rotational speed of particle at two previous time steps
        double temp_omgZ[2] {};          // Rotational speed of particle at two previous time steps
        double temp_part{};              // At first iteration it is Initial particle temperature and then takes "previous" iteration temperature of particle
        double temp_temperat[2] {};    // Temperature of particle at two previous time steps
        double delta_teta_X{};        // Particle change in rotation angle
        double delta_teta_Y{};        // Particle change in rotation angle
        double delta_teta_Z{};        // Particle change in rotation angle
        double Omg_prime_X{};         // ELLIPSOID only: Rotational velocity in body-fixed coordinate
        double Omg_prime_Y{};         // ELLIPSOID only: Rotational velocity in body-fixed coordinate
        double Omg_prime_Z{};         // ELLIPSOID only: Rotational velocity in body-fixed coordinate
        double Qu0{};                // Quaternions
        double Qu1{};                // Quaternions
        double Qu2{};                // Quaternions
        double Qu3{};                // Quaternions
        double displacement_verlet{}; // Displacement of particle sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2) between two Verlet list update
        double deter{};              // determinant of transformations matrix for spheroid
        double quater[10] {};       // quater
        double f3x{};               // Total Collision force on particle center due to particle-particle collision
        double f3y{};               // Total Collision force on particle center due to particle-particle collision
        double f3z{};               // Total Collision force on particle center due to particle-particle collision
        double f4x{};               // Total Collision force on particle center due to particle-wall collision
        double f4y{};               // Total Collision force on particle center due to particle-wall collision
        double f4z{};               // Total Collision force on particle center due to particle-wall collision
        double t1_x{};             // Torque on particle center from fluid
        double t1_y{};             // Torque on particle center from fluid
        double t1_z{};             // Torque on particle center from fluid
        double t2_x{};             // Torque on particle center from fluid
        double t2_y{};             // Torque on particle center from fluid
        double t2_z{};             // Torque on particle center from fluid
        double t3_x{};             // Torque on particle center due to collision with other particles or walls--> Just for Ellipsoid
        double t3_y{};             // Torque on particle center due to collision with other particles or walls--> Just for Ellipsoid
        double t3_z{};             // Torque on particle center due to collision with other particles or walls--> Just for Ellipsoid
        double fx_surf{};          // Total force on particle center
        double fy_surf{};          // Total force on particle center
        double fz_surf{};          // Total force on particle center
        double center_old_x{};     // Keeps particle center position before updating by Newton's equations
        double center_old_y{};     // Keeps particle center position before updating by Newton's equations
        double center_old_z{};     // Keeps particle center position before updating by Newton's equations
        double heat_surf{};        // Total heat on particle from fluid
        double T_prime_X{};        // Torques in body-fixed coordinate for Spheroid
        double T_prime_Y{};        // Torques in body-fixed coordinate
        double T_prime_Z{};        // Torques in body-fixed coordinate
        double euler_teta{};       // Euler angles. (radian). You can adjust them. Important: when all Euler angles are zero, particle_radius is along X-axis and particle_radius_2 is along Y and Z-axis.
        double euler_fi{};         // Euler angles
        double euler_sai{};        // Euler angles
    };
    // Define particle_struct
    struct particle_struct
    {
        // Elements
        double radius{};                   // object radius
        double radius_2{};                 // object radius
        double particle_mass{};
        double particle_perimeter{};
        double mom_inertia_x{};
        double mom_inertia_y{};
        double mom_inertia_z{};
        double parvol{};
        std::vector<center_struct> center; // center node
        std::vector<node_struct> node;     // list of nodes
    };

    std::vector<particle_struct> particle_list;


    ///======== Private methods ===========///
    void DF_force_calc(INT32 npn);
    void total_force_torque_calc(INT32 npn);    // Calculates total force on the particle center (at the end of all mdfs)
    void interpolate(INT32, const FlowSolver&, const ParallelMPI&, INT32); // interpolate node velocities from fluid velocity
    void spread(INT32, FlowSolver&, const ParallelMPI&); // spread node forces to fluid lattice (Notice that FlowSolver is not const becuase it changes here.)
    void update_verlet();
    void check_zero(double&, double);
    void set_particle_forces_torques_zero(INT32);
    void collision_model(INT32 np);
    void collision_circle(INT32 np);
    void collision_ellipse(INT32 np);
    /// For spheres only:
    void collision_sphere(INT32);//Collision forces for moving circular objects
    /// 8 functions for spheroids only:
    void rungekutta(INT32); // Runge-Kutta at each time step for the spheroid rotation
    void quaternion(INT32); //Calculation of quaternion matrix elements. Eq2.13 Article 876 Huang 2012
    double determinant(INT32);
    void comp_rot_speed(INT32);//update rotational speed. Just for SPHEROID
    void update_rot_spheroid(INT32);//update initial position of nodes if Euler_angles are not initially zero. Just for SPHEROID
    void update_rot_spheroid_collision(double, double &, double &, double &, double[], double, double, double, double, double, double); // only used in collision
    void collision_spheroid(INT32);//Collision forces for moving Spheroidal objects
    void normal_vec_spheroid(double, double&, double&, double&, INT32, double[], double, double, double);//Calculates NORMAL vector at each surface point (in INERTIAL coordinate) of spheroid . Used to calculate collision forces

    INT32 idx_3d_vector(INT32 x, INT32 y, INT32 z);// index for 3d vectors like temperature, density whih are a function of X, Y, Z
    INT32 idx_4d_vector(INT32 q, INT32 x, INT32 y, INT32 z);// index for 4d vectors like pop, velocity. Function of X, Y, Z and 1 other parameter

    void init_kernel_table();//Construct the table for kernels (two-point, four-point, etc.)
    double kernel_weight(double dist);//Calculate the weight from the kernel_table

    //mpi methods:
    void distribute_particles_mpi(INT32);

    void keep_particle_at_center(INT32, FlowSolver&);

};

#endif // PARTICLERESOLVED_H
