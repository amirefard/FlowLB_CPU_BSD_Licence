#ifndef PARTICLE_SIM_H
#define PARTICLE_SIM_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>  // file stream
#include <sstream>  // string streams
#include <cstdlib>  // standard library
#include <iomanip>  // For set precision.
#include <iostream> // for the use of 'cout'
#include <vector>
#include <cstring>   // for memcpy


#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"

#if defined Flow_With_Resolved_Particle


class Particle_sim
{
private:
    double total_displacement_x{0};   //Just used for keeping particle in center
    double partial_displacement_x{0}; //Just used for keeping particle in center
    double center_old_x, center_old_y;
    double area[particle_num_nodes];
    double particle_mass, particle_area, mom_inertia;
    int margin;
    const double KERNEL_MAX_DIST = 2.0;
    const double KERNEL_RESOLUTION = 0.0001;
    const int KERNEL_TABLE_SIZE = static_cast<int>(KERNEL_MAX_DIST / KERNEL_RESOLUTION) + 1;
    std::vector<double> kernel_table;

    void init_kernel_table();//Construct the table for kernels (two-point, four-point, etc.)
    double kernel_weight(double dist);//Calculate the weight from the kernel_table
    void keep_particle_at_center(int);
public:
    struct node_struct {

        /// Constructor

        node_struct() {
            x = 0;
            y = 0;
            vel_x = 0;
            vel_y = 0;
            force_x = 0;
            force_y = 0;
            total_node_mdf_force_x = 0;
            total_node_mdf_force_y = 0;
            den = 0;
            temperature = particle_temperature;
            force_thermal = 0;
            total_node_mdf_force_thermal = 0;
        }

        /// Elements
        double x; // current x-position
        double y; // current y-position
        double vel_x; // node velocity (x-component)
        double vel_y; // node velocity (y-component)
        double force_x; // node force (x-component) at each mdf step
        double force_y; // node force (y-component) at each mdf step
        double total_node_mdf_force_x; // total node force on eah node after mdf (x-component)
        double total_node_mdf_force_y; // total node force on each node after mdf (y-component)
        double den;//Density at boundary point
        double force_thermal; // node force due to thermal effect (it is not in Newton)
        double total_node_mdf_force_thermal; // total node force due to thermal effect at each node after all mdf steps (it is not in Newton)
        double temperature; // node force due to thermal effect (it is not in Newton)
    };

    struct center_struct {

        /// Constructor
        center_struct() {
            x = 0;
            y = 0;
            vel_x = 0;
            vel_y = 0;
            force_x = 0;
            force_y = 0;
            den = 0;
            temperature = particle_temperature;
            force_thermal = 0;
            omg = 0.0;
            previous_vel_x = 0.0;
            previous_vel_y = 0.0;
            teta = 0.0;
            tempx[0] = 0.0;
            tempx[1] = 0.0;
            tempy[0] = 0.0;
            tempy[1] = 0.0;
            temp_omg[0] = 0.0;
            temp_omg[1] = 0.0;
            temp_temperat[0] = 0.0;
            temp_temperat[1] = 0.0;
            temp_part = particle_temperature;
            delta_teta = 0;
            fx_surf = 0;
            fy_surf = 0;
            t1_z = 0;
            t2_z = 0;
            heat_surf = 0;
            f3x = 0;
            f3y = 0;
            f4x = 0;
            f4y = 0;
        }

        /// Elements
        double x; // current x-position
        double y; // current y-position
        double vel_x; // center velocity (x-component)
        double vel_y; // center velocity (y-component)
        double force_x; // center force (x-component)
        double force_y; // center force (y-component)
        double den;//Density at boundary point
        double omg; // center rotational speed
        double temperature; // node force due to thermal effect (it is not in Newton)
        double force_thermal; // center force due to thermal effect (it is not in Newton)
        double previous_vel_x; // Velocity x in previous time step
        double previous_vel_y; // Velocity y in previous time step
        double teta; // Angle
        double tempx[2];// Translational speed of particle at two previous time steps
        double tempy[2];// Translational speed of particle at two previous time steps
        double temp_omg[2]; // Rotational speed of particle at two previous time steps
        double temp_temperat[2]; // Temperature of particle at two previous time steps
        double temp_part; // At first iteration it is Initial particle temperature and then take the value of particle temperature at pervious iteration
        double delta_teta;
        double heat_surf;
        double fx_surf, fy_surf;// Total force on particle center
        double t1_z, t2_z;// Torque on particle center from fluid
        double t3_z;//Torque due to COLLISION-forces. Just for ELLIPTICAL particles
        double f3x, f3y; // Total particle-particle Collision force on particle center
        double f4x, f4y; // Total particle-wall Collision force on particle center
    };


    /// Structure for object (either cylinder or ellipse)
    struct particle_struct {

        /// Constructor

        particle_struct() {
            num_nodes = particle_num_nodes;
            center = new center_struct[1];
            center[0].x = 0;    //Will be initialized and changed later. particle_center_x;
            center[0].y = 0;    //Will be initialized and changed later.particle_center_y;
            node = new node_struct[num_nodes];
            // The initial shape of the object is set in main().
        }

        /// Elements
        int num_nodes; // number of surface nodes
        //double radius; // object radius
        center_struct *center; // center node
        node_struct *node; // list of nodes
    };

    particle_struct particle[num_particles]; // create immersed object

    ///Methods
    Particle_sim();
    void Particle_initialize();
    void set_forces_torques_zero(particle_struct);//set all forces and torques of the particle = 0
    void interpolate(particle_struct, int);
    void DF_force_calc(particle_struct, int);
    void total_force_torque_calc(particle_struct, int);
    void spread(particle_struct);
    void collision_cir(particle_struct[num_particles], int);//Collision forces for moving circular objects
    void collision_ellipse(particle_struct[num_particles], int);//Collision forces for moving elliptical objects
    void update_location(particle_struct particle[num_particles], int t);
    void write_particle_vtk(int, particle_struct[num_particles]); // write the particle state to the disk as VTK file
    void write_particle_data(int, particle_struct[num_particles]); // write particle data to the disk (drag/lift, center position)
    void write_particle_recovery(int, particle_struct[num_particles]);
    void read_particle_recovery(int, particle_struct[num_particles]);

    ~Particle_sim();

};

#endif

#endif
