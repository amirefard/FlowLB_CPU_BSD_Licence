//#include "stdafx.h"
#include "Point_Particle.h"


using namespace std;

#if defined Flow_With_Point_Particle

Point_Particle::Point_Particle()
{
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Particle_initialize()
//-----------------------------------------------------------------------------------------------------------------------------
void Point_Particle::Particle_initialize() {

	margin = 0;

#if defined Recovery_OFF

	/// ***********************************
	/// LOCATION OF PARTICLES
	/// ***********************************
	for (int i = 0; i < num_point_particles; ++i) {
		particle[i].radius = point_particle_radius;
		particle[i].den_ratio = point_den_ratio;

		particle[i].area = M_PI * SQ(particle[i].radius);
		particle[i].front_area = M_PI * 2 * particle[i].radius;
		particle[i].mass = (particle[i].den_ratio*1.0) * (M_PI*SQ(particle[i].radius));//Mass
		particle[i].mom_inertia = 0.5 * particle[i].mass * SQ(particle[i].radius);//Mass moment of inertia
	}

	particle[0].x = point_particle_center_x;
	particle[0].y = point_particle_center_y;

#if defined Random_Arrangement_OFF

	for (int i = 1; i < num_point_particles; ++i) {

		particle[i].x = point_particle_center_x + i * point_offset_x;//Location of other particles is adjusted according to location of first one. You can change it.
		particle[i].y = point_particle_center_y + i * point_offset_y;

		//cout << particle[0].x << endl;

		if (particle[i].x > Nx - 0.001 || particle[i].x < 0.001) {
			cout << "ERROR: Particle " << i + 1 << " position is out of the domain";
			return;
		}
		if (particle[i].y > Ny - 0.001 || particle[i].y < 0.001) {
			cout << "ERROR: Particle " << i + 1 << " position is out of the domain";
			return;
		}
}
#endif
#if defined Random_Arrangement_ON
	for (int i = 1; i < num_point_particles; ++i) {
		int part_position;
		srand(5323);
		double max_dist;
//#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
		max_dist = particle[0].radius;
//#endif
//#if defined MOVING_SPHEROID || defined STATIONARY_SPHEROID
//		max_dist = MAX(particle[0].radius_2, particle[0].radius);
//#endif
		do {
			part_position = TRUE;
			particle[i].x = rand() % Nx;
			particle[i].y = rand() % Ny;
			printf("Particle i is %d and x, y, z are %5.3f , %5.3f \n", i + 1, particle[i].x, particle[i].y);

			for (int j = 0; j < num_point_particles; ++j) {
				if (j != i) {
					double Dij = sqrt(SQ(particle[i].x - particle[j].x) + SQ(particle[i].y - particle[j].y));
					//#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
					if (Dij < 2 * max_dist + 1) {
						part_position = FALSE;
						printf("Error, contact with particle %d \n", j);
						//	break;
					}
					//#endif
				}
			}

			if (particle[i].x > Nx - 1 - max_dist - 1 || particle[i].x < max_dist + 1) {
				part_position = FALSE;
				cout << "Error 1" << endl;
			}
			if (particle[i].y > Ny - 1 - max_dist - 1 || particle[i].y < max_dist + 1) {
				part_position = FALSE;
				cout << "Error 2" << endl;
			}

		} while (part_position == FALSE);
	}
#endif
	//***************************************
#endif  // End for Recovery_OFF
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: INTERPOLATE
//-----------------------------------------------------------------------------------------------------------------------------
//
// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) and other interpolation schemes are included in the present code.
void Point_Particle::interpolate(point_particle_struct & particle) {

	int x_int;
	int y_int;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	double dist_x;
	double dist_y;
	double weight_x;
	double weight_y;

	// Run over all object nodes.
	//for (int n = 0; n < particle.num_nodes; ++n) {

		// Reset node velocity first since '+=' is used.
		particle.vel_fluid_x = 0;
		particle.vel_fluid_y = 0;
		//particle.node[n].den = 0;
		//particle.node[n].temperature = 0;

		x_int = floor(particle.x);
		y_int = floor(particle.y);

		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
#if defined Four_Point_Delta
		xmin = x_int - 2;
		xmax = x_int + 2;
		ymin = y_int - 2;
		ymax = y_int + 2;
#endif
#if defined Two_Point_delta
		xmin = x_int - 1;
		xmax = x_int + 1;
		ymin = y_int - 1;
		ymax = y_int + 1;
#endif


#if defined X_NOT_Periodic
		if (xmin <= 1) {
			xmin = 1;
		}
		if (xmax >= Nx - 2) {
			xmax = Nx - 2;
		}
#endif
#if defined Y_NOT_Periodic
		if (ymin <= 1) {
			ymin = 1;
		}
		if (ymax >= Ny - 2) {
			ymax = Ny - 2;
		}
#endif

		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
		for (int X = xmin; X <= xmax; ++X)
			for (int Y = ymin; Y <= ymax; ++Y)
			{
				{

					// Compute distance between object node and fluid lattice node.
					dist_x = particle.x - X;
					dist_y = particle.y - Y;

#if defined Two_Point_delta
					// Compute interpolation weights for x- and y-direction based on the distance.
					if (abs(dist_x) <= 1) {
						weight_x = 1. - abs(dist_x);
					}
					else weight_x = 0.0;
					if (abs(dist_y) <= 1) {
						weight_y = 1. - abs(dist_y);
					}
					else weight_y = 0.0;
#endif

#if defined Four_Point_Delta
					if (abs(dist_x) >= 0 && abs(dist_x)< 1) {
						weight_x = 1. / 8.0*(3. - 2.* abs(dist_x) + sqrt(1 + 4.*abs(dist_x) - 4.*SQ(dist_x)));
					}
					else if (abs(dist_x) >= 1 && abs(dist_x)< 2) {
						weight_x = 1. / 8.0*(5. - 2.* abs(dist_x) - sqrt(-7 + 12.*abs(dist_x) - 4.*SQ(dist_x)));
					}
					else weight_x = 0.0;

					if (abs(dist_y) >= 0 && abs(dist_y)< 1) {
						weight_y = 1. / 8.0*(3. - 2.* abs(dist_y) + sqrt(1 + 4.*abs(dist_y) - 4.*SQ(dist_y)));
					}
					else if (abs(dist_y) >= 1 && abs(dist_y)< 2) {
						weight_y = 1. / 8.0*(5. - 2.* abs(dist_y) - sqrt(-7 + 12.*abs(dist_y) - 4.*SQ(dist_y)));
					}
					else weight_y = 0.0;
#endif

					particle.vel_fluid_x += ((velocity_x[(X + Nx) % Nx][(Y + Ny) % Ny]) * weight_x * weight_y);
					particle.vel_fluid_y += ((velocity_y[(X + Nx) % Nx][(Y + Ny) % Ny]) * weight_x * weight_y);

				}
			}
		//cout << "hi " << particle.vel_fluid_y << endl;
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: force_calc
//-----------------------------------------------------------------------------------------------------------------------------
//void Particle_sim::force_calc(particle_struct particle, int npn){
void Point_Particle::force_calc(point_particle_struct & particle) {
	double C_d_x, Re_p_x;
	double C_d_y, Re_p_y;

	Re_p_x = 2 * particle.radius  * abs(particle.vel_fluid_x - particle.vel_x) / ((tau - 0.5) / 3.0);
	if (Re_p_x>0)
	{
		C_d_x = 10. / Re_p_x;  //Stokes regime
		//C_d_x = 24. / Re_p_x;  //Stokes regime for SPHERE
	}
	else
	{
		C_d_x = 0;
	}
	//cout << "Reynolds_p X " << Re_p_x << endl;
	Re_p_y = 2 * particle.radius  * abs(particle.vel_fluid_y - particle.vel_y) / ((tau - 0.5) / 3.0);
	if (Re_p_y>0)
	{
		C_d_y = 10. / Re_p_y;  //Stokes regime for CIRCLE
		//C_d_y = 24. / Re_p_y;  //Stokes regime for SPHERE

	}
	else
	{
		C_d_y = 0;
	}
	//cout << "rrr " << particle.vel_fluid_x << endl;
	particle.force_x = 0.5 * C_d_x * particle.front_area * (particle.vel_fluid_x - particle.vel_x) * abs(particle.vel_fluid_x - particle.vel_x);
	particle.force_y = 0.5 * C_d_y * particle.front_area * (particle.vel_fluid_y - particle.vel_y) * abs(particle.vel_fluid_y - particle.vel_y);
	//cout << "salam_aaa" << Re_p_y<<" "<<particle.force_y << endl;

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: update_location
//-----------------------------------------------------------------------------------------------------------------------------
void Point_Particle::update_location(point_particle_struct particle[num_point_particles], int t) {

	for (int npn = 0; npn < num_point_particles; ++npn) {

		////////// Updating Velocity and center position ///////
		///// Tracer
		//particle[npn].vel_x = particle[npn].vel_fluid_x;// particle[npn].vel_x + particle[npn].force_x / particle[npn].mass + point_particle_gravity_x * (1.0 - 1.0 / particle[npn].den_ratio);// +particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
		//particle[npn].vel_y = particle[npn].vel_fluid_y; //particle[npn].vel_y + particle[npn].force_y / particle[npn].mass + point_particle_gravity_y * (1.0 - 1.0 / particle[npn].den_ratio);// +particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law


	//	if (t > margin) {
		//	if (t > margin + 1) {
#if defined With_added_mass
		particle[npn].vel_x = particle[npn].vel_x + particle[npn].force_x / particle[npn].mass + point_particle_gravity_x * (1.0 - 1.0 / particle[npn].den_ratio) + 1.0 / particle[npn].den_ratio*(particle[npn].vel_x - particle[npn].tempx[0]);// +particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
		particle[npn].vel_y = particle[npn].vel_y + particle[npn].force_y / particle[npn].mass + point_particle_gravity_y * (1.0 - 1.0 / particle[npn].den_ratio) + 1.0 / particle[npn].den_ratio*(particle[npn].vel_y - particle[npn].tempy[0]);;// +particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
#endif
#if defined Without_added_mass
		/// Inertial
		particle[npn].vel_x += particle[npn].force_x / particle[npn].mass + point_particle_gravity_x * (1.0 - 1.0 / particle[npn].den_ratio);// +particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
		particle[npn].vel_y += particle[npn].force_y / particle[npn].mass + point_particle_gravity_y * (1.0 - 1.0 / particle[npn].den_ratio);// +particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
#endif

			//}
			//else {
			//	//particle[npn].center[0].vel_x = particle[npn].center[0].vel_x - particle[npn].center[0].fx_surf / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
			//	//particle[npn].center[0].vel_y = particle[npn].center[0].vel_y - particle[npn].center[0].fy_surf / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3y / particle_mass + particle[npn].center[0].f4y / particle_mass;//Newton law
			//}
	//	}
		if (t == 1) {
			particle[npn].tempx[0] = 0;
			particle[npn].tempx[1] = 0;
			particle[npn].tempy[0] = 0;
			particle[npn].tempy[1] = 0;
		}
		particle[npn].tempx[0] = particle[npn].tempx[1];
		particle[npn].tempx[1] = particle[npn].vel_x;
		particle[npn].tempy[0] = particle[npn].tempy[1];
		particle[npn].tempy[1] = particle[npn].vel_y;

		center_old_x = particle[npn].x;
		center_old_y = particle[npn].y;

		//particle[npn].x = particle[npn].x + 0.5*(particle[npn].vel_x + particle[npn].tempx[0]);
		//particle[npn].y = particle[npn].y + 0.5*(particle[npn].vel_y + particle[npn].tempy[0]);

		particle[npn].x = particle[npn].x + particle[npn].vel_x;
		particle[npn].y = particle[npn].y + particle[npn].vel_y;

////////////////////////////////////////////

//		////////// Updating Angular Velocity ///////
//		if (t > margin) {
//			if (t > margin + 1) {
//#if defined With_added_mass
//				particle[npn].center[0].omg = particle[npn].center[0].omg + 1.0 / (mom_inertia)*(particle[npn].center[0].t1_z - particle[npn].center[0].t2_z - particle[npn].center[0].t3_z) + 1.0 / den_ratio*(particle[npn].center[0].omg - particle[npn].center[0].temp_omg[0]);//Eq. 100 Kang Thesis (804)
//#endif
//#if defined Without_added_mass
//				particle[npn].center[0].omg = particle[npn].center[0].omg + 1.0 / (mom_inertia)*(particle[npn].center[0].t1_z - particle[npn].center[0].t2_z - particle[npn].center[0].t3_z);
//#endif
//			}
//			else {
//				particle[npn].center[0].omg = particle[npn].center[0].omg + 1.0 / (mom_inertia)*(particle[npn].center[0].t1_z - particle[npn].center[0].t2_z - particle[npn].center[0].t3_z);
//			}
//		}
//		if (t == 1) {
//			particle[npn].center[0].temp_omg[0] = 0.0;
//			particle[npn].center[0].temp_omg[1] = 0.0;
//		}
//		particle[npn].center[0].temp_omg[0] = particle[npn].center[0].temp_omg[1];
//		particle[npn].center[0].temp_omg[1] = particle[npn].center[0].omg;
//
//		particle[npn].center[0].teta = particle[npn].center[0].teta + 0.5 * (particle[npn].center[0].omg + particle[npn].center[0].temp_omg[0]);
//		particle[npn].center[0].delta_teta = 0.5 * (particle[npn].center[0].omg + particle[npn].center[0].temp_omg[0]);
//		////////////////////////////////////////////
//
//		////////// Updating Temperature ///////
//#if defined Particle_Varying_Temperature
//		if (t > 1) {
//			particle[npn].center[0].temperature = particle[npn].center[0].temperature - particle[npn].center[0].heat_surf / (particle_mass*specific_heat_ratio) + 1.0 / (den_ratio*specific_heat_ratio)*(particle[npn].center[0].temperature - particle[npn].center[0].temp_temperat[0]) + particle_heat_source / (particle_mass*specific_heat_ratio);
//		}
//		else {
//			particle[npn].center[0].temperature = particle[npn].center[0].temperature - particle[npn].center[0].heat_surf / (particle_mass*specific_heat_ratio) + particle_heat_source / (particle_mass*specific_heat_ratio);
//		}
//#endif
//		if (t == 1) {
//			particle[npn].center[0].temp_temperat[0] = particle_temperature;
//			particle[npn].center[0].temp_temperat[1] = particle_temperature;
//		}
//		particle[npn].center[0].temp_temperat[0] = particle[npn].center[0].temp_temperat[1];
//		particle[npn].center[0].temp_temperat[1] = particle[npn].center[0].temperature;
//		/////////////////////////////////////////
//
//		double x_pos, y_pos;
//		////////// Updating nodes position ///////
//		for (int n = 0; n < particle[npn].num_nodes; ++n) {
//			x_pos = particle[npn].node[n].x - center_old_x;
//			y_pos = particle[npn].node[n].y - center_old_y;
//			particle[npn].node[n].x = particle[npn].center[0].x + x_pos * cos(particle[npn].center[0].delta_teta) - y_pos * sin(particle[npn].center[0].delta_teta);
//			particle[npn].node[n].y = particle[npn].center[0].y + x_pos * sin(particle[npn].center[0].delta_teta) + y_pos * cos(particle[npn].center[0].delta_teta);
//
//			particle[npn].node[n].vel_x = particle[npn].center[0].vel_x;
//			particle[npn].node[n].vel_y = particle[npn].center[0].vel_y;
//			particle[npn].node[n].temperature = particle[npn].center[0].temperature;
//		}

#if defined X_Periodic
		if (particle[npn].x < 0) {
			particle[npn].x += Nx;

		}
		else if (particle[npn].x >= Nx) {
			particle[npn].x -= Nx;

		}
#endif
#if defined Y_Periodic
		if (particle[npn].y < 0) {
			particle[npn].y += Ny;

		}
		else if (particle[npn].y >= Ny) {
			particle[npn].y -= Ny;

		}
#endif
		particle[npn].previous_vel_x = particle[npn].vel_x;
		particle[npn].previous_vel_y = particle[npn].vel_y;

//#if defined Particle_Varying_Temperature
//		particle[npn].center[0].temp_part = particle[npn].center[0].temperature;
//#endif
	}
	return;
};

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// The particle state (node positions) is writen to a VTK file at each t_vtk step and can be read by ParaView

void Point_Particle::write_particle_vtk(int time, point_particle_struct particle[num_point_particles]) {

#if defined VTK_BINARY || defined VTK_ASCII

	for (int npn = 0; npn < num_point_particles; ++npn) {
		/// Create filename
		stringstream output_filename;
        output_filename << "LB_Results/vtk_particle/particle_" << npn + 1 << "_t" << time << ".vtk";
		ofstream output_file;

		/// Open file
		output_file.open(output_filename.str().c_str());

		/// Write VTK header
		output_file << "# vtk DataFile Version 3.0\n";
		output_file << "particle_state\n";
		output_file << "ASCII\n";
		output_file << "DATASET POLYDATA\n";

		/// Write node positions
		output_file << "POINTS " << 1 << " float\n";
			output_file << particle[npn].x << " " << particle[npn].y << " 0\n";
			//cout << time << " mmm" << particle[npn].x << endl;
			output_file << "VERTICES 1 " << 1 + 1 << "\n";
			output_file << 1 << "\n";
			for (int n = 0; n < 1; ++n) {
				output_file << n << "\n";
			}





		/// Write lines between neighboring nodes
		//output_file << "LINES " << particle_num_nodes << " " << 3 * particle_num_nodes << "\n";
		//for (int n = 0; n < particle_num_nodes; ++n) {
		//	output_file << "2 " << n << " " << (n + 1) % particle_num_nodes << "\n";
		//}

		/// Write vertices
		//output_file << "VERTICES 1 " << 1 + 1 << "\n";
		//output_file << 1 << " ";
	//	for (int n = 0; n < particle_num_nodes; ++n) {
		//	output_file << 0 << " ";
	//	}

		/// Close file
		output_file.close();
	}
#endif
#if defined CASE
	unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);// number are shown as 2.0e02 instead of 2.0e002 (necessay for geo and data files)


																		  /////////// Writing CASE file (.case) ///////////////////
																		  /////////////////////////////////////////////////////////
																		  /// Create filename
	for (int npn = 0; npn < num_particles; ++npn) {

		stringstream output_filename;
        output_filename << "LB_Results/vtk_particle/particle_" << npn + 1 << ".case";
		ofstream output_file;

		/// Open file
		output_file.open(output_filename.str().c_str());

		/// Write VTK header
		output_file << "# EnSight Gold Model\n";
        output_file << "# Created by LB\n\n";
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

		double min_node_x, min_node_y;
		double max_node_x, max_node_y;
		min_node_x = particle.node[0].x;
		min_node_y = particle.node[0].y;
		max_node_x = particle.node[0].x;
		max_node_y = particle.node[0].y;

		for (int n = 0; n < particle_num_nodes; ++n) {
			if (particle.node[n].x<min_node_x) {
				min_node_x = particle.node[n].x;
			}
			if (particle.node[n].y<min_node_y) {
				min_node_y = particle.node[n].y;
			}
			if (particle.node[n].x>max_node_x) {
				max_node_x = particle.node[n].x;
			}
			if (particle.node[n].y>max_node_y) {
				max_node_y = particle.node[n].y;
			}
		}
		stringstream output_filename1;
        output_filename1 << "LB_Results/vtk_particle/particle_" << npp + 1 << "_" << setfill('0') << setw(6) << time << ".geo";//Geometry file
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
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << "\n";
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
			output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << "\n";
		}
		output_file1 << "bar2\n";
		output_file1 << std::setw(10) << particle_num_nodes << "\n";
		for (int itm = 1; itm <= particle_num_nodes; itm++) {
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int I = 1; I < particle_num_nodes; I++) {
			output_file1 << std::setw(10) << I << std::setw(10) << I + 1 << "\n";
		}
		output_file1 << std::setw(10) << particle_num_nodes << std::setw(10) << 1 << "\n";

		/// Close file
		output_file1.close();
	}
	/////////// End of Writing GEOMETRY file (.geo) ///////////////////
	///////////////////////////////////////////////////////////////////


#endif//End of CASE file for Ensight

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_DATA
//-----------------------------------------------------------------------------------------------------------------------------
void Point_Particle::write_particle_data(int time, point_particle_struct particle[num_point_particles]) {
	for (int npn = 0; npn < num_point_particles; ++npn) {

		/// Create filename
		stringstream output_filename;
        output_filename << "LB_Results/Data_particle/data_particle_" << npn + 1 << ".dat";
		ofstream output_file;

		/// Open file
		output_file.open(output_filename.str().c_str(), ofstream::app);


		/// Write data
		output_file << time << " "; // time step
		output_file << particle[npn].force_x << " "; // drag force
		output_file << particle[npn].force_y << " "; // lift force
		output_file << particle[npn].x << " "; // center position (x-component)
		output_file << particle[npn].y << " "; // center position (y-component)
		output_file << particle[npn].vel_x << " "; // center velocity (x-component)
		output_file << particle[npn].vel_y << endl; // center velocity (y-component)

		/// Close file
		output_file.close();
	}
	return;
}


Point_Particle::~Point_Particle()
{
}

#endif
