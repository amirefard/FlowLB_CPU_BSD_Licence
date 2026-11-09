

#if defined Flow_With_Point_Particle


#ifndef POINT_PARTICLE_H
#define POINT_PARTICLE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>  // string streams
#include <fstream>  // file stream
#include <cstdlib>  // standard library
#include <iostream> // for the use of 'cout'
#include <iomanip>  // For set precision. From ver 28
//#include <conio.h>

#include "Global_variables.h"
#include "LB_setting.h"
#include "LB_Macros.h"


class Point_Particle
{
private:
	double center_old_x, center_old_y;

	int margin;

public:


	struct point_particle_struct {

		/// Constructor
		point_particle_struct() {
			x = 0;
			y = 0;
			radius = 0;
			area = 0;
			front_area=0;
			mass = 0;
			mom_inertia = 0;
			den_ratio = 1;
			vel_x = 0;
			vel_y = 0;
			vel_fluid_x = 0;
			vel_fluid_y = 0;
			force_x = 0;
			force_y = 0;
			den = 0;
			//temperature = particle_temperature;
			//force_thermal = 0;
			//omg = 0.0;
			previous_vel_x = 0.0;
			previous_vel_y = 0.0;
			teta = 0.0;
			tempx[0] = 0.0;
			tempx[1] = 0.0;
			tempy[0] = 0.0;
			tempy[1] = 0.0;
			//temp_omg[0] = 0.0;
			//temp_omg[1] = 0.0;
			//temp_temperat[0] = 0.0;
			//temp_temperat[1] = 0.0;
			//temp_part = particle_temperature;
			//delta_teta = 0;
			//fx_surf = 0;
			//fy_surf = 0;
			//t1_z = 0;
			//t2_z = 0;
			//heat_surf = 0;
			//f3x = 0;
			//f3y = 0;
			//f4x = 0;
			//f4y = 0;
		}

		/// Elements
		double x; // current x-position
		double y; // current y-position
		double radius;
		double area;
		double front_area;
		double mass;
		double mom_inertia;
		double den_ratio;
		double vel_x; // center velocity (x-component)
		double vel_y; // center velocity (y-component)
		double vel_fluid_x;// Fluid velocity at the particle position
		double vel_fluid_y;
		double force_x; // center force (x-component)
		double force_y; // center force (y-component)
		double den;//Density at boundary point
		//double omg; // center rotational speed
		//double temperature; // node force due to thermal effect (it is not in Newton)
		//double force_thermal; // center force due to thermal effect (it is not in Newton)
		double previous_vel_x; // Velocity x in previous time step
		double previous_vel_y; // Velocity y in previous time step
		double teta; // Angle
		double tempx[2];// Translational speed of particle at two previous time steps
		double tempy[2];// Translational speed of particle at two previous time steps
		//double temp_omg[2]; // Rotational speed of particle at two previous time steps
		//double temp_temperat[2]; // Temperature of particle at two previous time steps
		//double temp_part; // At first iteration it is Initial particle temperature and then take the value of particle temperature at pervious iteration
		//double delta_teta;
		//double heat_surf;
		//double fx_surf, fy_surf;// Total force on particle center
		//double t1_z, t2_z;// Torque on particle center from fluid
		//double t3_z;//Torque due to COLLISION-forces. Just for ELLIPTICAL particles
		//double f3x, f3y; // Total particle-particle Collision force on particle center
		//double f4x, f4y; // Total particle-wall Collision force on particle center
	};
	point_particle_struct particle[num_point_particles]; // create immersed object


	public:
	Point_Particle();
	void Particle_initialize();
	void interpolate(point_particle_struct &);
	void force_calc(point_particle_struct & );
	void update_location(point_particle_struct particle[num_point_particles], int t);
	void write_particle_vtk(int, point_particle_struct particle[num_point_particles]);
	void write_particle_data(int, point_particle_struct particle[num_point_particles]); // write particle data to the disk (drag/lift, center position)



	~Point_Particle();
};

#endif
#endif
