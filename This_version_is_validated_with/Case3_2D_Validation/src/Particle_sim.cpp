#include "Particle_sim.h"
using namespace std;


//particle_struct particle[num_particles]; // create immersed object
//double area[particle_num_nodes];
//double particle_mass, particle_area, mom_inertia;
//double particle_area, mom_inertia;
//int margin;

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------

Particle_sim::Particle_sim()
{
    init_kernel_table();
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: init_kernel_table()
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::init_kernel_table()
{
    kernel_table.resize(KERNEL_TABLE_SIZE);
    for (int i = 0; i < KERNEL_TABLE_SIZE; ++i)
    {
        double r = i * KERNEL_RESOLUTION;
#if defined(Kernel_Delta)
        kernel_table[i] = (r <= 2.0) ? 0.25 * (1.0 + cos(M_PI * r / 2.0)) : 0.0;
#elif defined(Four_Point_Delta)
        if (r < 1.0)
            kernel_table[i] = 0.125 * (3.0 - 2.0 * r + sqrt(1 + 4 * r - 4 * r * r));
        else if (r >= 1 && r < 2.0)
            kernel_table[i] = 0.125 * (5.0 - 2.0 * r - sqrt(-7 + 12 * r - 4 * r * r));
        else
            kernel_table[i] = 0.0;
#elif defined(Two_Point_delta)
        kernel_table[i] = (r <= 1.0) ? (1.0 - r) : 0.0;
#endif
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: kernel_weight()
//-----------------------------------------------------------------------------------------------------------------------------
double Particle_sim::kernel_weight(double dist)
{
    double r = std::abs(dist);
    if (r > KERNEL_MAX_DIST)
        return 0.0;
    int idx = static_cast<int>(r / KERNEL_RESOLUTION);
    return kernel_table[idx];
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Particle_initialize()
//-----------------------------------------------------------------------------------------------------------------------------

void Particle_sim::Particle_initialize(){
	margin = 0000;//This parameter is helpful for STATIONARY_CYLINDER or STATIONARY_ELLIPSE. You can change this line if you want the particles to be
	//... stationary until a certain iteration (margin) and then starts motion.

#if defined STATIONARY_CYLINDER || defined STATIONARY_ELLIPSE
	margin = t_num + 1;//Do not change this line
#endif

#if defined Recovery_OFF

	/// ********************************
	/// CREATING OBJECTS (1st PARTICLE)
	/// ********************************
	//Defining particles initial positions. (Based on Shape)
	if (num_particles > 0) {
		particle[0].center[0].x = particle_center_x;
		particle[0].center[0].y = particle_center_y;
		for (int n = 0; n < particle[0].num_nodes; ++n) {
            #if defined STATIONARY_CYLINDER || defined MOVING_CYLINDER
			particle[0].node[n].x = particle[0].center[0].x + particle_radius * sin(2. * M_PI * (double)n / particle[0].num_nodes);
			particle[0].node[n].y = particle[0].center[0].y + particle_radius * cos(2. * M_PI * (double)n / particle[0].num_nodes);
			area[n] = 2.0 * M_PI * particle_radius / particle[0].num_nodes;//Perimeter divided by number of nodes for CIRCLE
            #endif
            #if defined STATIONARY_ELLIPSE || defined MOVING_ELLIPSE
			particle[0].node[n].x = particle[0].center[0].x + particle_radius   * cos(2. * M_PI * (double)n / particle[0].num_nodes);
			particle[0].node[n].y = particle[0].center[0].y + particle_radius_2 * sin(2. * M_PI * (double)n / particle[0].num_nodes);
            #endif
			if (particle[0].node[n].x > Nx - 0.001 || particle[0].node[n].x < 0.001){
				cout << "ERROR: Particle " << 1 << " position is out of the domain";
				return;
			}
			if (particle[0].node[n].y > Ny - 0.001 || particle[0].node[n].y < 0.001){
				cout << "ERROR: Particle " << 1 << " position is out of the domain";
				return;
			}
		}
	}

	/// ***********************************
	/// LOCATION OF 2nd UNTIL Nth PARTICLE
	/// ***********************************
	for (int i = 1; i < num_particles; ++i) {
		particle[i].center[0].x = particle_center_x + i * offset_x;//Location of other particles is adjusted according to location of first one. You can change it.
		particle[i].center[0].y = particle_center_y + i * offset_y;
		for (int n = 0; n < particle[i].num_nodes; ++n) {
			particle[i].node[n].x = particle[0].node[n].x - particle[0].center[0].x + particle[i].center[0].x;
			particle[i].node[n].y = particle[0].node[n].y - particle[0].center[0].y + particle[i].center[0].y;
			if (particle[i].node[n].x > Nx - 0.001 || particle[i].node[n].x < 0.001){
				cout << "ERROR: Particle " << i + 1 << " position is out of the domain";
				return;
			}
			if (particle[i].node[n].y > Ny - 0.001 || particle[i].node[n].y < 0.001){
				cout << "ERROR: Particle " << i + 1 << " position is out of the domain";
				return;
			}
		}
	}

	/// *****************************
	/// INITIAL ROTATION OF ELLIPSE
	/// *****************************
	// If the ellipse first radius (particle_radius) is not initially alligned along Y-direction it is adjusted here. (Like paper 885, Suzuki, Inamuro test case)
	// This angle is measured in counter clock wise direction
#if defined STATIONARY_ELLIPSE || defined MOVING_ELLIPSE
	double x_pos_1, y_pos_1;
	for (int i = 0; i < num_particles; ++i) {
		particle[i].center[0].teta = particle[i].center[0].teta + M_PI / 2. - init_ang;
		for (int n = 0; n < particle[i].num_nodes; ++n) {
			x_pos_1 = particle[i].node[n].x - particle[i].center[0].x;
			y_pos_1 = particle[i].node[n].y - particle[i].center[0].y;
			particle[i].node[n].x = particle[i].center[0].x + x_pos_1 * cos(init_ang) + y_pos_1 * sin(init_ang);
			particle[i].node[n].y = particle[i].center[0].y - x_pos_1 * sin(init_ang) + y_pos_1 * cos(init_ang);
		}
	}
#endif


	/// ********************************
	/// NODE AREA FOR ELLIPSE PARTICLE
	/// ********************************
	if (num_particles > 0){
		double dist1, dist2;
		double sumarea;
		sumarea = 0.0;
        #if defined MOVING_ELLIPSE || defined STATIONARY_ELLIPSE
		for (int n = 0; n < particle[0].num_nodes; ++n) {
			dist1 = sqrt(SQ(particle[0].node[(n + 1) % particle[0].num_nodes].x - particle[0].node[n].x) + SQ(particle[0].node[(n + 1) % particle[0].num_nodes].y - particle[0].node[n].y));
			dist2 = sqrt(SQ(particle[0].node[n].x - particle[0].node[(n - 1 + particle[0].num_nodes) % particle[0].num_nodes].x) + SQ(particle[0].node[n].y - particle[0].node[(n - 1 + particle[0].num_nodes) % particle[0].num_nodes].y));
			area[n] = (dist1 + dist2) / 2.;
			sumarea += area[n];
		}
		//cout<<"sum area"<<sumarea<<endl;
        #endif
	}


	/// ---------------------------------------------------------------------------------------------
	/// -------------------------- PARTICLE MASS AND MOMENT OF INERTIA ------------------------------
	/// ---------------------------------------------------------------------------------------------
#if defined STATIONARY_CYLINDER || defined MOVING_CYLINDER
	//particle_area = 2.0 * M_PI * particle_radius;//Perimeter
	particle_mass = (den_ratio*1.0) * (M_PI*SQ(particle_radius));//Mass
	mom_inertia = 0.5 * particle_mass * SQ(particle_radius);//Mass moment of inertia
#endif

#if defined STATIONARY_ELLIPSE || defined MOVING_ELLIPSE
	//particle_area = M_PI * (3.0*(particle_radius+particle_radius_2) - sqrt((3.0*particle_radius+particle_radius_2)*(particle_radius+3.0*particle_radius_2))) ;
	particle_mass = (den_ratio*1.0) * M_PI * particle_radius * particle_radius_2;
	mom_inertia = 0.25 * particle_mass * (SQ(particle_radius) + SQ(particle_radius_2));
#endif

	//***************************************
#endif  // End for Recovery_OFF && Flow_With_Resolved_Particle
	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SET_FORCES_TORQUES_ZERO
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::set_forces_torques_zero(particle_struct particle){//set all forces and torques of the particle = 0
    for (int nn = 0; nn < particle.num_nodes; ++nn) {
        particle.node[nn].force_x = 0;
        particle.node[nn].force_y = 0;
        particle.node[nn].force_thermal = 0.0;
        particle.node[nn].total_node_mdf_force_x = 0;
        particle.node[nn].total_node_mdf_force_y = 0;
        particle.node[nn].total_node_mdf_force_thermal = 0.0;
    }
    particle.center[0].fx_surf = 0.0;//Total Force x due to IBM
    particle.center[0].fy_surf = 0.0;//Total Force y due to IBM
    particle.center[0].heat_surf = 0.0;//Total heat due to TIBM
    particle.center[0].f3x = 0.0; //Total Force x due to particle-particle collision
    particle.center[0].f3y = 0.0; //Total Force y due to particle-particle collision
    particle.center[0].f4x = 0.0; //Forces x due to particle-wall collision
    particle.center[0].f4y = 0.0; //Forces y due to particle-wall collision
    particle.center[0].t1_z = 0.0;//Torque z due to forces in X direction
    particle.center[0].t2_z = 0.0;//Torque z due to forces in Y direction
    particle.center[0].t3_z = 0.0;//Torque z due to COLLISION-forces. Just for ELLIPTICAL particles
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: DF_FORCE_CALC
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::DF_force_calc(particle_struct particle, int npn){
	for (int n = 0; n < particle.num_nodes; ++n) {	//Compute Lagrangian Forces
        ///force at each mdf step
		particle.node[n].force_x       = (1) * (particle.center[0].previous_vel_x - particle.center[0].omg * (particle.node[n].y - particle.center[0].y) - particle.node[n].vel_x) * area[n];//Surface force
		particle.node[n].force_y       = (1) * (particle.center[0].previous_vel_y + particle.center[0].omg * (particle.node[n].x - particle.center[0].x) - particle.node[n].vel_y) * area[n];//Surface force
		particle.node[n].force_thermal = (1) * (particle.center[0].temp_part - particle.node[n].temperature) * area[n];//Surface force

		///total force on each node after each mdf step
		particle.node[n].total_node_mdf_force_x       += particle.node[n].force_x;
        particle.node[n].total_node_mdf_force_y       += particle.node[n].force_y;
		particle.node[n].total_node_mdf_force_thermal += particle.node[n].force_thermal;
	}

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: TOTAL_FORCE_TORQUE_CALC
//-----------------------------------------------------------------------------------------------------------------------------
// Calculates total force on the particle center (at the end of all mdfs)
void Particle_sim::total_force_torque_calc(particle_struct particle, int npn){
	for (int n = 0; n < particle.num_nodes; ++n) {	//Compute Lagrangian Forces
        particle.center[0].fx_surf += particle.node[n].total_node_mdf_force_x;
		particle.center[0].fy_surf += particle.node[n].total_node_mdf_force_y;

		particle.center[0].t1_z += particle.node[n].total_node_mdf_force_x * (particle.node[n].y - particle.center[0].y);
		particle.center[0].t2_z += particle.node[n].total_node_mdf_force_y * (particle.node[n].x - particle.center[0].x);

		particle.center[0].heat_surf += particle.node[n].total_node_mdf_force_thermal;
	}

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: update_location
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::update_location(particle_struct particle[num_particles], int t){

	for (int npn = 0; npn < num_particles; ++npn) {
        /*if (particle[npn].center[0].t1_z < 1.0e-6)
        {
            particle[npn].center[0].t1_z = 0.0;
        }
        if (particle[npn].center[0].t2_z < 1.0e-6)
        {
            particle[npn].center[0].t2_z = 0.0;
        }*/

		/// Updating Velocity and center position ///
		if (t > margin){
			if (t > margin + 1){
#if defined With_added_mass
				particle[npn].center[0].vel_x = particle[npn].center[0].vel_x - particle[npn].center[0].fx_surf / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass + 1.0 / den_ratio*(particle[npn].center[0].vel_x - particle[npn].center[0].tempx[0]);//Newton law
				particle[npn].center[0].vel_y = particle[npn].center[0].vel_y - particle[npn].center[0].fy_surf / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3y / particle_mass + particle[npn].center[0].f4y / particle_mass + 1.0 / den_ratio*(particle[npn].center[0].vel_y - particle[npn].center[0].tempy[0]);//Newton law
#endif
#if defined Without_added_mass
				particle[npn].center[0].vel_x = particle[npn].center[0].vel_x - particle[npn].center[0].fx_surf / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
				particle[npn].center[0].vel_y = particle[npn].center[0].vel_y - particle[npn].center[0].fy_surf / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3y / particle_mass + particle[npn].center[0].f4y / particle_mass;//Newton law
#endif
			}
			else{
				particle[npn].center[0].vel_x = particle[npn].center[0].vel_x - particle[npn].center[0].fx_surf / particle_mass + particle_gravity_x * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3x / particle_mass + particle[npn].center[0].f4x / particle_mass;//Newton law
				particle[npn].center[0].vel_y = particle[npn].center[0].vel_y - particle[npn].center[0].fy_surf / particle_mass + particle_gravity_y * (1.0 - 1.0 / den_ratio) + particle[npn].center[0].f3y / particle_mass + particle[npn].center[0].f4y / particle_mass;//Newton law
			}
		}
		if (t == 1){
			particle[npn].center[0].tempx[0] = 0;
			particle[npn].center[0].tempx[1] = 0;
			particle[npn].center[0].tempy[0] = 0;
			particle[npn].center[0].tempy[1] = 0;
		}
		particle[npn].center[0].tempx[0] = particle[npn].center[0].tempx[1];
		particle[npn].center[0].tempx[1] = particle[npn].center[0].vel_x;
		particle[npn].center[0].tempy[0] = particle[npn].center[0].tempy[1];
		particle[npn].center[0].tempy[1] = particle[npn].center[0].vel_y;

		center_old_x = particle[npn].center[0].x;
		center_old_y = particle[npn].center[0].y;
		particle[npn].center[0].x = particle[npn].center[0].x + 0.5*(particle[npn].center[0].vel_x + particle[npn].center[0].tempx[0]);
		particle[npn].center[0].y = particle[npn].center[0].y + 0.5*(particle[npn].center[0].vel_y + particle[npn].center[0].tempy[0]);
		////////////////////////////////////////////

		keep_particle_at_center(npn);///To keep particle at center. Does not work with Recovery_ON at the moment. Change write_particle_data() if needed.

		/// Updating Angular Velocity ///
		if (t > margin){
			if (t > margin + 1){
#if defined With_added_mass
				particle[npn].center[0].omg = particle[npn].center[0].omg + 1.0 / (mom_inertia)*(particle[npn].center[0].t1_z - particle[npn].center[0].t2_z - particle[npn].center[0].t3_z) + 1.0 / den_ratio*(particle[npn].center[0].omg - particle[npn].center[0].temp_omg[0]);//Eq. 100 Kang Thesis (804)
#endif
#if defined Without_added_mass
				particle[npn].center[0].omg = particle[npn].center[0].omg + 1.0 / (mom_inertia)*(particle[npn].center[0].t1_z - particle[npn].center[0].t2_z - particle[npn].center[0].t3_z);
#endif
			}
			else{
				particle[npn].center[0].omg = particle[npn].center[0].omg + 1.0 / (mom_inertia)*(particle[npn].center[0].t1_z - particle[npn].center[0].t2_z - particle[npn].center[0].t3_z);
			}
		}
		if (t == 1){
			particle[npn].center[0].temp_omg[0] = 0.0;
			particle[npn].center[0].temp_omg[1] = 0.0;
		}
		particle[npn].center[0].temp_omg[0] = particle[npn].center[0].temp_omg[1];
		particle[npn].center[0].temp_omg[1] = particle[npn].center[0].omg;

		particle[npn].center[0].teta = particle[npn].center[0].teta + 0.5 * (particle[npn].center[0].omg + particle[npn].center[0].temp_omg[0]);
		particle[npn].center[0].delta_teta = 0.5 * (particle[npn].center[0].omg + particle[npn].center[0].temp_omg[0]);
		////////////////////////////////////////////

		/// Updating Particle Temperature ///
#if defined Particle_Varying_Temperature
		if (t > 1){
#if defined With_added_heat
			particle[npn].center[0].temperature = particle[npn].center[0].temperature + Dt_LB * (-particle[npn].center[0].heat_surf / (particle_mass*specific_heat_ratio) + particle_heat_source / (particle_mass*specific_heat_ratio)) + 1.0 / (den_ratio*specific_heat_ratio) * (particle[npn].center[0].temperature - particle[npn].center[0].temp_temperat[0]);
#endif
#if defined Without_added_heat
			particle[npn].center[0].temperature = particle[npn].center[0].temperature + Dt_LB * (-particle[npn].center[0].heat_surf / (particle_mass*specific_heat_ratio) + particle_heat_source / (particle_mass*specific_heat_ratio));
#endif
		}
		else{
			particle[npn].center[0].temperature = particle[npn].center[0].temperature + Dt_LB * (-particle[npn].center[0].heat_surf / (particle_mass*specific_heat_ratio) + particle_heat_source / (particle_mass*specific_heat_ratio));
		}
#endif
		if (t == 1){
			particle[npn].center[0].temp_temperat[0] = particle_temperature;
			particle[npn].center[0].temp_temperat[1] = particle_temperature;
		}
		particle[npn].center[0].temp_temperat[0] = particle[npn].center[0].temp_temperat[1];
		particle[npn].center[0].temp_temperat[1] = particle[npn].center[0].temperature;
		/////////////////////////////////////////

		double x_pos, y_pos;
		////////// Updating nodes position ///////
		for (int n = 0; n < particle[npn].num_nodes; ++n) {
			x_pos = particle[npn].node[n].x - center_old_x;
			y_pos = particle[npn].node[n].y - center_old_y;
			particle[npn].node[n].x = particle[npn].center[0].x + x_pos * cos(particle[npn].center[0].delta_teta) - y_pos * sin(particle[npn].center[0].delta_teta);
			particle[npn].node[n].y = particle[npn].center[0].y + x_pos * sin(particle[npn].center[0].delta_teta) + y_pos * cos(particle[npn].center[0].delta_teta);

			particle[npn].node[n].vel_x = particle[npn].center[0].vel_x;
			particle[npn].node[n].vel_y = particle[npn].center[0].vel_y;
			particle[npn].node[n].temperature = particle[npn].center[0].temperature;
		}

#if defined X_Periodic
		if (particle[npn].center[0].x < 0) {
			particle[npn].center[0].x += Nx - 1;
			for (int np1 = 0; np1 < particle[npn].num_nodes; ++np1) {
				particle[npn].node[np1].x += Nx - 1;
			}
		}
		else if (particle[npn].center[0].x > Nx - 1) {
			particle[npn].center[0].x -= Nx - 1;
			for (int np1 = 0; np1 < particle[npn].num_nodes; ++np1) {
				particle[npn].node[np1].x -= Nx - 1;
			}
		}
#endif
#if defined Y_Periodic
		if (particle[npn].center[0].y < 0) {
			particle[npn].center[0].y += Ny - 1;
			for (int np1 = 0; np1 < particle[npn].num_nodes; ++np1) {
				particle[npn].node[np1].y += Ny - 1;
			}
		}
		else if (particle[npn].center[0].y > Ny - 1) {
			particle[npn].center[0].y -= Ny - 1;
			for (int np1 = 0; np1 < particle[np].num_nodes; ++np1) {
				particle[npn].node[np1].y -= Ny - 1;
			}
		}
#endif
		particle[npn].center[0].previous_vel_x = particle[npn].center[0].vel_x;
		particle[npn].center[0].previous_vel_y = particle[npn].center[0].vel_y;

#if defined Particle_Varying_Temperature
		particle[npn].center[0].temp_part = particle[npn].center[0].temperature;
#endif
	}
	return;
};
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: INTERPOLATE
//-----------------------------------------------------------------------------------------------------------------------------
//
// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) and other interpolation schemes are included in the present code.
void Particle_sim::interpolate(particle_struct particle, int mdf_iteration) {

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
	for (int n = 0; n < particle.num_nodes; ++n) {
		// Reset node velocity first since '+=' is used.
		particle.node[n].vel_x = 0;
		particle.node[n].vel_y = 0;
		particle.node[n].den = 0;
		particle.node[n].temperature = 0;

		// Identify the lowest fluid lattice node in interpolation range (see spreading).
		x_int = floor(particle.node[n].x);
		y_int = floor(particle.node[n].y);
		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
#if defined Kernel_Delta || defined Four_Point_Delta
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
		if (xmin <= 0){
			xmin = 0;
		}
		if (xmax >= Nx - 1){
			xmax = Nx - 1;
		}
#endif
#if defined Y_NOT_Periodic
		if (ymin <= 0){
			ymin = 0;
		}
		if (ymax >= Ny - 1){
			ymax = Ny - 1;
		}
#endif

		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
		for (int X = xmin; X <= xmax; ++X)
		for (int Y = ymin; Y <= ymax; ++Y)
		{
			{
				// Compute distance between object node and fluid lattice node.
				dist_x = particle.node[n].x - X;
				dist_y = particle.node[n].y - Y;

/*				// Compute interpolation weights for x- and y-direction based on the distance.
#if defined Kernel_Delta
				if (abs(dist_x) <= 2){
					weight_x = 1.0 / 4.0*(1. + cos(M_PI*abs(dist_x) / 2.0));
				}
				else weight_x = 0.0;
				if (abs(dist_y) <= 2){
					weight_y = 1.0 / 4.0*(1.0 + cos(M_PI*abs(dist_y) / 2.0));
				}
				else weight_y = 0.0;
#endif

#if defined Four_Point_Delta
				if (abs(dist_x) >= 0 && abs(dist_x)< 1){
					weight_x = 1. / 8.0*(3. - 2.* abs(dist_x) + sqrt(1 + 4.*abs(dist_x) - 4.*SQ(dist_x)));
				}
				else if (abs(dist_x) >= 1 && abs(dist_x)< 2){
					weight_x = 1. / 8.0*(5. - 2.* abs(dist_x) - sqrt(-7 + 12.*abs(dist_x) - 4.*SQ(dist_x)));
				}
				else weight_x = 0.0;

				if (abs(dist_y) >= 0 && abs(dist_y)< 1){
					weight_y = 1. / 8.0*(3. - 2.* abs(dist_y) + sqrt(1 + 4.*abs(dist_y) - 4.*SQ(dist_y)));
				}
				else if (abs(dist_y) >= 1 && abs(dist_y)< 2){
					weight_y = 1. / 8.0*(5. - 2.* abs(dist_y) - sqrt(-7 + 12.*abs(dist_y) - 4.*SQ(dist_y)));
				}
				else weight_y = 0.0;
#endif

#if defined Two_Point_delta
				if (abs(dist_x) <= 1){
					weight_x = 1. - abs(dist_x);
				}
				else weight_x = 0.0;
				if (abs(dist_y) <= 1){
					weight_y = 1. - abs(dist_y);
				}
				else weight_y = 0.0;
#endif*/
                weight_x = kernel_weight(dist_x);
                weight_y = kernel_weight(dist_y);

                const int Xnew = (X + Nx) % Nx;
                const int Ynew = (Y + Ny) % Ny;
                const double wxyz = weight_x * weight_y;

if (mdf_iteration == 0) {//First interpolation
				particle.node[n].vel_x += ((velocity_x[Xnew][Ynew]) * wxyz);
				particle.node[n].vel_y += ((velocity_y[Xnew][Ynew]) * wxyz);
				particle.node[n].den   += ((density[Xnew][Ynew])    * wxyz);
#if defined Flow_With_Thermal_Effect
				particle.node[n].temperature += ((temperature[Xnew][Ynew]) * wxyz);
#endif
}else{// Next mdf iterations
				particle.node[n].vel_x += ((mdf_vel_x[Xnew][Ynew]) * wxyz);
				particle.node[n].vel_y += ((mdf_vel_y[Xnew][Ynew]) * wxyz);
				particle.node[n].den   += ((density[Xnew][Ynew])   * wxyz);
#if defined Flow_With_Thermal_Effect
				particle.node[n].temperature += ((mdf_temperature[Xnew][Ynew]) * wxyz);
#endif
}
			}
		}
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
void Particle_sim::spread(particle_struct particle) {

	/// Reset forces (Now done in momenta and momenta_t)
	// This is necessary since '+=' is used afterwards.
	/*for(int X = 0; X < Nx; ++X) {
	for(int Y = 1; Y < Ny - 1; ++Y) {
        force_x[X][Y] = 0;
        force_y[X][Y] = 0;
        force_thermal[X][Y] = 0;
	}
	}*/

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
	for (int n = 0; n < particle.num_nodes; ++n) {

		// Identify the lowest fluid lattice node in interpolation range (see spreading).
		x_int = floor(particle.node[n].x);
		y_int = floor(particle.node[n].y);
		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
#if defined Kernel_Delta || defined Four_Point_Delta
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
		if (xmin <= 0){
			xmin = 0;
		}
		if (xmax >= Nx - 1){
			xmax = Nx - 1;
		}
#endif
#if defined Y_NOT_Periodic
		if (ymin <= 0){
			ymin = 0;
		}
		if (ymax >= Ny - 1){
			ymax = Ny - 1;
		}
#endif

		// Run over all neighboring fluid nodes.
		// In the case of the two-point interpolation, it is 2x2 fluid nodes.
		for (int X = xmin; X <= xmax; ++X)
		for (int Y = ymin; Y <= ymax; ++Y)
		{
			{
				// Compute distance between object node and fluid lattice node.
				dist_x = particle.node[n].x - X;
				dist_y = particle.node[n].y - Y;

		/*		// Compute interpolation weights for x- and y-direction based on the distance.
#if defined Kernel_Delta
				if (abs(dist_x) <= 2){
					weight_x = 1.0 / 4.0*(1. + cos(M_PI*abs(dist_x) / 2.0));
				}
				else weight_x = 0.0;
				if (abs(dist_y) <= 2){
					weight_y = 1.0 / 4.0*(1. + cos(M_PI*abs(dist_y) / 2.0));
				}
				else weight_y = 0.0;
#endif

#if defined Four_Point_Delta
				if (abs(dist_x) >= 0 && abs(dist_x)< 1){
					weight_x = 1. / 8.0*(3. - 2.* abs(dist_x) + sqrt(1 + 4.*abs(dist_x) - 4.*SQ(dist_x)));
				}
				else if (abs(dist_x) >= 1 && abs(dist_x)< 2){
					weight_x = 1. / 8.0*(5. - 2.* abs(dist_x) - sqrt(-7 + 12.*abs(dist_x) - 4.*SQ(dist_x)));
				}
				else weight_x = 0.0;

				if (abs(dist_y) >= 0 && abs(dist_y)< 1){
					weight_y = 1. / 8.0*(3. - 2.* abs(dist_y) + sqrt(1 + 4.*abs(dist_y) - 4.*SQ(dist_y)));
				}
				else if (abs(dist_y) >= 1 && abs(dist_y)< 2){
					weight_y = 1. / 8.0*(5. - 2.* abs(dist_y) - sqrt(-7 + 12.*abs(dist_y) - 4.*SQ(dist_y)));
				}
				else weight_y = 0.0;
#endif

#if defined Two_Point_delta
				if (abs(dist_x) <= 1){
					weight_x = 1. - abs(dist_x);
				}
				else weight_x = 0.0;
				if (abs(dist_y) <= 1){
					weight_y = 1. - abs(dist_y);
				}
				else weight_y = 0.0;
#endif*/

                weight_x = kernel_weight(dist_x);
                weight_y = kernel_weight(dist_y);

                const int Xnew = (X + Nx) % Nx;
                const int Ynew = (Y + Ny) % Ny;
                const double wxyz = weight_x * weight_y;

				force_x[Xnew][Ynew] += (particle.node[n].force_x * wxyz);
				force_y[Xnew][Ynew] += (particle.node[n].force_y * wxyz);
#if defined Flow_With_Thermal_Effect
				force_thermal[Xnew][Ynew] += (particle.node[n].force_thermal * wxyz);
#endif
			}
		}
	}

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: keep_particle_at_center
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::keep_particle_at_center(int npn){

    /// To keep particle at center (Does not work with Recovery_ON currently)
    const int displacement_before_shift = 5;// How many cell particle moves before shift of domain happens
    const int desire_cen_x  = int(particle_center_x);// To always keep the particle at center
    total_displacement_x   += 0.5 * (particle[npn].center[0].vel_x + particle[npn].center[0].tempx[0]);// Just for saving in write_particle_data
    partial_displacement_x += 0.5 * (particle[npn].center[0].vel_x + particle[npn].center[0].tempx[0]);

    std::vector<std::vector<double>> swap_vel_mov_x(Nx, std::vector<double>(Ny));// Vectors are defined on the heap
    std::vector<std::vector<double>> swap_vel_mov_y(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_den_mov(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_force_mov_x(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_force_mov_y(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_force_thermal(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_temperature(Nx, std::vector<double>(Ny));

    std::vector<std::vector<std::vector<double>>> swap_pop(
        stencil_directions_F, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny))
    );

    std::vector<std::vector<std::vector<double>>> swap_pop_t(
        stencil_directions_T, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny))
    );

    if (abs (partial_displacement_x) > displacement_before_shift){
        //double*** swap_pop           = pop;// Wrong: shallow copy
        //double**  swap_vel_mov_x     = velocity_x;// Wrong: shallow copy
        //double**  swap_vel_mov_y     = velocity_y;// Wrong: shallow copy
        //double**  swap_den_mov       = density;// Wrong: shallow copy
        //double**  swap_force_mov_x   = force_x;// Wrong: shallow copy
        //double**  swap_force_mov_y   = force_y;// Wrong: shallow copy
        //double**  swap_force_thermal = force_thermal;// Wrong: shallow copy
        //double**  swap_temperature   = temperature;// Wrong: shallow copy
        //double*** swap_pop_t         = pop_t;

        for (int i = 0; i < Nx; ++i) {
            std::copy(velocity_x[i], velocity_x[i] + Ny, swap_vel_mov_x[i].begin());
            std::copy(velocity_y[i], velocity_y[i] + Ny, swap_vel_mov_y[i].begin());
            std::copy(density[i], density[i] + Ny, swap_den_mov[i].begin());
            std::copy(force_x[i], force_x[i] + Ny, swap_force_mov_x[i].begin());
            std::copy(force_y[i], force_y[i] + Ny, swap_force_mov_y[i].begin());
            std::copy(force_thermal[i], force_thermal[i] + Ny, swap_force_thermal[i].begin());
            std::copy(temperature[i], temperature[i] + Ny, swap_temperature[i].begin());
        }

        for (int d = 0; d < stencil_directions_F; ++d) {
            for (int i = 0; i < Nx; ++i) {
                std::copy(pop_F[d][i], pop_F[d][i] + Ny, swap_pop[d][i].begin());
            }
        }

        for (int d = 0; d < stencil_directions_T; ++d) {
            for (int i = 0; i < Nx; ++i) {
                std::copy(pop_T[d][i], pop_T[d][i] + Ny, swap_pop_t[d][i].begin());
            }
        }

        int start_X = desire_cen_x - 600;
        int end_X   = desire_cen_x + 400;

        for (int X = start_X; X < end_X; ++X)
        {
            for (int Y = 0; Y < Ny; ++Y)
            {
                int copy_from_X = int(floor(particle[npn].center[0].x)) - desire_cen_x + X;

                velocity_x[X][Y]    = swap_vel_mov_x[copy_from_X][Y];
                velocity_y[X][Y]    = swap_vel_mov_y[copy_from_X][Y];
                force_x[X][Y]       = swap_force_mov_x[copy_from_X][Y];
                force_y[X][Y]       = swap_force_mov_y[copy_from_X][Y];
                density[X][Y]       = swap_den_mov[copy_from_X][Y];
                temperature[X][Y]   = swap_temperature[copy_from_X][Y];
                force_thermal[X][Y] = swap_force_thermal[copy_from_X][Y];
                for (int c_i = 0; c_i < stencil_directions_F; ++c_i) {
                    pop_F[c_i][X][Y]  = swap_pop[c_i][copy_from_X][Y];
                }
                for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
                    pop_T[c_i][X][Y] = swap_pop_t[c_i][copy_from_X][Y];
                }
            }
        }

        particle[npn].center[0].x -= floor(partial_displacement_x);//Now shift the particle position
        partial_displacement_x    -= floor(partial_displacement_x);//Now reset partial_displacement_x
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_CIR
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for circular objects
void Particle_sim::collision_cir(particle_struct particle[num_particles], int np) {

#if defined Circle_Collision_1
	//Collision Forces. JUST for Cylinder shape
	double Dij;
	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dij = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dij <= (2.0 * particle_radius + th)) {
				particle[np].center[0].f3x += cij / ep * ((particle[np].center[0].x - particle[np2].center[0].x)) / Dij * pow((Dij - 2.0 * particle_radius - th) / th, 2);
				particle[np].center[0].f3y += cij / ep * ((particle[np].center[0].y - particle[np2].center[0].y)) / Dij * pow((Dij - 2.0 * particle_radius - th) / th, 2);
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_NOT_Periodic
	if (particle[np].center[0].y <= (Ny - 2) / 2){//Particle-Y-Wall Collision section
		Rij = particle[np].center[0].y + particle_radius;
	}
	else
	{
		Rij = -1.0 * ((double)(Ny - 1) - particle[np].center[0].y + particle_radius);
	}

	if (abs(Rij) <= (2.0 * particle_radius + thw)){
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += cij / epw * pow((2 * particle_radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
	}
#endif

#if defined X_NOT_Periodic
	if (particle[np].center[0].x <= (Nx - 2) / 2){//Particle-X-Wall Collision section
		Rij = particle[np].center[0].x + particle_radius;
	}
	else
	{
		Rij = -1.0 * ((double)(Nx - 1) - particle[np].center[0].x + particle_radius);
	}

	if (abs(Rij) <= (2.0 * particle_radius + thw)){
		particle[np].center[0].f4x += cij / epw * pow((2.0 * particle_radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4y += 0;
	}
#endif	 // END of Particle-Wall Collision
#endif //end of Circle_Collision_1

#if defined Circle_Collision_2
	//Collision Forces. JUST for Cylinder shape
	//const double ep = 1.0;//stifness factor for particle-particle collision
	//const double epw = 0.5;//stifness factor for particle-wall collision
	//const double th = 1.0;// Thereshold
	double Dij;
	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dij = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dij <= (2.0*particle_radius)){
				particle[np].center[0].f3x += 1.0 / ep * ((particle[np].center[0].x - particle[np2].center[0].x)) * (2.0 * particle_radius - Dij);//Total particle-particle collision force
				particle[np].center[0].f3y += 1.0 / ep * ((particle[np].center[0].y - particle[np2].center[0].y)) * (2.0 * particle_radius - Dij);
			}
			if (Dij >(2 * particle_radius) && Dij <= (2.0*particle_radius + th)){
				particle[np].center[0].f3x += 1.0 / ep * ((particle[np].center[0].x - particle[np2].center[0].x)) * pow((2.0 * particle_radius - Dij + th), 2);
				particle[np].center[0].f3y += 1.0 / ep * ((particle[np].center[0].y - particle[np2].center[0].y)) * pow((2.0 * particle_radius - Dij + th), 2);
			}
		}
	} // END of Particle-Particle Collision

	double Rij;
#if defined Y_NOT_Periodic
	if (particle[np].center[0].y < (Ny - 2) / 2){//Particle-Wall Collision section
		Rij = particle[np].center[0].y + particle_radius;
	}
	else
	{
		Rij = -1.0 * ((Ny - 1) - particle[np].center[0].y + particle_radius);
	}
	if (abs(Rij) < 2.0 * particle_radius) {
		particle[np].center[0].f4x += 1.0 / epw * 0.;
		particle[np].center[0].f4y += 1.0 / epw * Rij * (2 * particle_radius - abs(Rij));
	}
	if (abs(Rij) >= 2 * particle_radius && abs(Rij) <= (2.0 * particle_radius + thw)){
		particle[np].center[0].f4x += 1.0 / epw * 0;
		particle[np].center[0].f4y += 1.0 / epw * Rij * pow((2 * particle_radius - abs(Rij) + thw), 2);
	}
#endif
#if defined X_NOT_Periodic
	if (particle[np].center[0].x < (Nx - 2) / 2){//Particle-Y-Wall Collision section
		Rij = particle[np].center[0].x + particle_radius;
	}
	else
	{
		Rij = -1 * ((Nx - 1) - particle[np].center[0].x + particle_radius);
	}
	if (abs(Rij) < 2 * particle_radius) {
		particle[np].center[0].f4x += 1 / epw * Rij * (2.0 * particle_radius - abs(Rij));
		particle[np].center[0].f4y += 1 / epw * 0;
	}
	if (abs(Rij) >= 2 * particle_radius && abs(Rij) <= (2 * particle_radius + thw)){
		particle[np].center[0].f4x += 1 / epw * Rij * pow((2.0 * particle_radius - abs(Rij) + thw), 2);
		particle[np].center[0].f4y += 1 / epw * 0.;
	}
#endif

	// END of Particle-Wall Collision
#endif //end of Circle_Collision_2

#if defined Circle_Collision_3
	//Collision Forces. JUST for Cylinder shape

	double Dij;
	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dij = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dij <= (2.0 * particle_radius + th) && Dij>(2.0 * particle_radius)) {
				particle[np].center[0].f3x += cij / ep * ((particle[np].center[0].x - particle[np2].center[0].x)) / Dij * pow((Dij - 2.0 * particle_radius - th) / th, 2);
				particle[np].center[0].f3y += cij / ep * ((particle[np].center[0].y - particle[np2].center[0].y)) / Dij * pow((Dij - 2.0 * particle_radius - th) / th, 2);
			}
			if (Dij <= (2.0 * particle_radius)) {
				particle[np].center[0].f3x += (particle[np].center[0].x - particle[np2].center[0].x) / Dij * (cij / ep *pow((Dij - 2.0 * particle_radius - th) / th, 2) + cij / Ep / th*(2 * particle_radius - Dij));
				particle[np].center[0].f3y += (particle[np].center[0].y - particle[np2].center[0].y) / Dij * (cij / ep *pow((Dij - 2.0 * particle_radius - th) / th, 2) + cij / Ep / th*(2 * particle_radius - Dij));
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_NOT_Periodic
	if (particle[np].center[0].y <= (Ny - 2) / 2){//Particle-Y-Wall Collision section
		Rij = particle[np].center[0].y + particle_radius;
	}
	else
	{
		Rij = -1.0 * ((Ny - 1) - particle[np].center[0].y + particle_radius);
	}

	if (abs(Rij) <= (2.0 * particle_radius + thw) && abs(Rij) > (2.0 * particle_radius)){
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += cij / epw * pow((2.0 * particle_radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
	}
	if (abs(Rij) <= (2 * particle_radius)){
		particle[np].center[0].f4x += 0.0;
		particle[np].center[0].f4y += (cij / epw *pow((abs(Rij) - 2.0 * particle_radius - thw) / thw, 2.0) + cij / Epw / thw*(2.0*particle_radius - abs(Rij)))*Rij / abs(Rij);
	}
#endif

#if defined X_NOT_Periodic
	if (particle[np].center[0].x <= (Nx - 2) / 2){//Particle-Wall Collision section
		Rij = particle[np].center[0].x + particle_radius;
	}
	else
	{
		Rij = -1.0 * ((Nx - 1) - particle[np].center[0].x + particle_radius);
	}

	if (abs(Rij) <= (2.0 * particle_radius + thw) && abs(Rij) > (2.0 * particle_radius)){
		particle[np].center[0].f4x += cij / epw * pow((2.0 * particle_radius - abs(Rij) + thw) / thw, 2)*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
	}
	if (abs(Rij) <= (2.0 * particle_radius)){
		particle[np].center[0].f4x += (cij / epw *pow((abs(Rij) - 2.0 * particle_radius - thw) / thw, 2) + cij / Epw / thw*(2.0*particle_radius - abs(Rij)))*Rij / abs(Rij);
		particle[np].center[0].f4y += 0.0;
	}
#endif
	// END of Particle-Wall Collision
#endif //end of Circle_Collision_3
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_ELLIPSE
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for Elliptical objects
void Particle_sim::collision_ellipse(particle_struct particle[num_particles], int np) {

#if defined Ellipse_Collision_1
	//Collision Forces. JUST for ELLIPSE shape
	double Dij, Dist, ffx, ffy, ffx4, ffy4, angle_1, angle_2, fas, xnor, ynor, norvec, xnor_init, ynor_init;
	ffx = 0.0;
	ffy = 0.0;

	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dist = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dist <= (2 * MAX(particle_radius, particle_radius_2) + th + 2)){
				for (int n = 0; n < particle[np].num_nodes; ++n) {
					angle_1 = atan((particle[np].node[0].y - particle[np].center[0].y) / (particle[np].node[0].x - particle[np].center[0].x));//It gives us the angle of particle np
					xnor = (2 * (particle[np].node[n].x - particle[np].center[0].x)*cos(angle_1) - 2 * (particle[np].node[n].y - particle[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
					ynor = (2 * (particle[np].node[n].x - particle[np].center[0].x)*sin(angle_1) + 2 * (particle[np].node[n].y - particle[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
					xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
					ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
					norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
					for (int m = 0; m < particle[np2].num_nodes; ++m) {
						Dij = sqrt(SQ(particle[np].node[n].x - particle[np2].node[m].x) + SQ(particle[np].node[n].y - particle[np2].node[m].y));
						angle_2 = atan((particle[np2].node[0].y - particle[np2].center[0].y) / (particle[np2].node[0].x - particle[np2].center[0].x));//It gives us the angle of particle np2
						fas = SQ((particle[np].node[n].x - particle[np2].center[0].x)*cos(angle_2) + (particle[np].node[n].y - particle[np2].center[0].y)*sin(angle_2)) / SQ(particle_radius) + SQ((particle[np].node[n].x - particle[np2].center[0].x)*sin(angle_2) - (particle[np].node[n].y - particle[np2].center[0].y)*cos(angle_2)) / SQ(particle_radius_2);

						if (Dij <= th  && fas > 1) { //It means that particles have NOT penetrated into each other
							// for circle// alwaya keep deactive// ffx = -1. * cij / ep * ((particle[np].node[n].x-particle[np].center[0].x))/sqrt(SQ(particle[np].node[n].x-particle[np].center[0].x)+SQ(particle[np].node[n].y-particle[np].center[0].y))  * pow((Dij - th) / th , 2);
							// for circle// alwaya keep deactive// ffy = -1. * cij / ep * ((particle[np].node[n].y-particle[np].center[0].y))/sqrt(SQ(particle[np].node[n].x-particle[np].center[0].x)+SQ(particle[np].node[n].y-particle[np].center[0].y))  * pow((Dij - th) / th , 2);
							ffx = -1. * cij / ep * xnor_init / norvec * pow((Dij - th) / th, 2);
							ffy = -1. * cij / ep * ynor_init / norvec * pow((Dij - th) / th, 2);
							particle[np].center[0].f3x += ffx;
							particle[np].center[0].f3y += ffy;
							particle[np].center[0].t3_z += ffx * (particle[np].node[n].y - particle[np].center[0].y) - ffy * (particle[np].node[n].x - particle[np].center[0].x);
						}
						/*if (Dij <= th && fas <= 1) {//It means that node n of particle np has penetrated to particle np2
						ffx = -1. * cij / Ep * xnor_init / norvec * pow((Dij - th) / th , 2);
						ffy = -1. * cij / Ep * ynor_init / norvec * pow((Dij - th) / th , 2);
						particle[np].center[0].f3x += ffx;
						particle[np].center[0].f3y += ffy;
						particle[np].center[0].t3_z += ffx * (particle[np].node[n].y-particle[np].center[0].y) - ffy * (particle[np].node[n].x-particle[np].center[0].x);
						}*/

					}
				}
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_NOT_Periodic
	if ((particle[np].center[0].y <= MAX(particle_radius, particle_radius_2) + thw + 2) || ((Ny - particle[np].center[0].y) <= MAX(particle_radius, particle_radius_2) + thw + 2)){
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			angle_1 = atan((particle[np].node[0].y - particle[np].center[0].y) / (particle[np].node[0].x - particle[np].center[0].x));//It gives us the angle of particle np
			xnor = (2 * (particle[np].node[n].x - particle[np].center[0].x)*cos(angle_1) - 2 * (particle[np].node[n].y - particle[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
			ynor = (2 * (particle[np].node[n].x - particle[np].center[0].x)*sin(angle_1) + 2 * (particle[np].node[n].y - particle[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
			xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
			ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
			norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
			if (particle[np].node[n].y <= (Ny - 1) / 2){//Particle-Y-Wall Collision section
				Rij = particle[np].node[n].y - 0;// 0 means y-location of lower wall
			}
			else
			{
				Rij = -1.0 * ((Ny - 1) - particle[np].node[n].y);// Ny-1 means y-location of upper wall
			}
			if (abs(Rij) <= (thw)){
				ffx4 = 0.0;
				ffy4 = -1. * cij / epw * ynor_init / norvec * pow((-abs(Rij) + thw) / thw, 2);
				particle[np].center[0].f4x += ffx4;
				particle[np].center[0].f4y += ffy4;
				particle[np].center[0].t3_z += ffx4 * (particle[np].node[n].y - particle[np].center[0].y) - ffy4 * (particle[np].node[n].x - particle[np].center[0].x);
			}
		}
	}
#endif
#if defined X_NOT_Periodic
	if ((particle[np].center[0].x <= MAX(particle_radius, particle_radius_2) + thw + 2) || ((Nx - particle[np].center[0].x) <= MAX(particle_radius, particle_radius_2) + thw + 2)){
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			angle_1 = atan((particle[np].node[0].y - particle[np].center[0].y) / (particle[np].node[0].x - particle[np].center[0].x));//It gives us the angle of particle np
			xnor = (2 * (particle[np].node[n].x - particle[np].center[0].x)*cos(angle_1) - 2 * (particle[np].node[n].y - particle[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
			ynor = (2 * (particle[np].node[n].x - particle[np].center[0].x)*sin(angle_1) + 2 * (particle[np].node[n].y - particle[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
			xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
			ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
			norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
			if (particle[np].node[n].x <= (Nx - 1) / 2.){//Particle-X-Wall Collision section
				Rij = particle[np].node[n].x - 0;// 0 means x-location of left wall
			}
			else
			{
				Rij = -1.0 * ((Nx - 1) - particle[np].node[n].x);// Nx-1 means x-location of right wall
			}

			if (abs(Rij) <= (thw)){
				// 	cout<<"angle "<<angle_1<<"xnor "<<xnor_init<<" norvec"<<norvec<<endl;
				ffx4 = -1. * cij / epw * xnor_init / norvec * pow((-abs(Rij) + thw) / thw, 2);
				ffy4 = 0;
				particle[np].center[0].f4x += ffx4;
				particle[np].center[0].f4y += ffy4;
				particle[np].center[0].t3_z += ffx4 * (particle[np].node[n].y - particle[np].center[0].y) - ffy4 * (particle[np].node[n].x - particle[np].center[0].x);
			}
		}
	}
#endif	 // END of Particle-Wall Collision
#endif //end of Ellipse_Collision_1

#if defined Ellipse_Collision_2
	//Collision Forces. JUST for ELLIPSE shape
	const double min_th = 0.0;// If distance is less than this value, Ep is used instead of ep
	const double min_thw = 0.0;// If distance is less than this value, Epw is used instead of epw
	double Dij, Dist, ffx, ffy, ffx4, ffy4, angle_1, angle_2, fas, xnor, ynor, norvec, xnor_init, ynor_init;
	ffx = 0.0;
	ffy = 0.0;
	int node_1;
	int node_2;
	double Dij_min;//Just to initially have a value greater than th

	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dist = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dist <= (2 * MAX(particle_radius, particle_radius_2) + th)){
				node_1 = 0;
				node_2 = 0;
				Dij_min = th + 1;//Just to initially have a value greater than th
				for (int n = 0; n < particle[np].num_nodes; ++n) {
					for (int m = 0; m < particle[np2].num_nodes; ++m) {
						Dij = sqrt(SQ(particle[np].node[n].x - particle[np2].node[m].x) + SQ(particle[np].node[n].y - particle[np2].node[m].y));
						if (Dij <= Dij_min) {
							Dij_min = Dij;
							node_1 = n;
							node_2 = m;
						}
					}
				}
				if (Dij_min >min_th && Dij_min <= th){
					angle_1 = atan((particle[np].node[0].y - particle[np].center[0].y) / (particle[np].node[0].x - particle[np].center[0].x));//It gives us the angle of particle np
					xnor = (2 * (particle[np].node[node_1].x - particle[np].center[0].x)*cos(angle_1) - 2 * (particle[np].node[node_1].y - particle[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
					ynor = (2 * (particle[np].node[node_1].x - particle[np].center[0].x)*sin(angle_1) + 2 * (particle[np].node[node_1].y - particle[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
					xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
					ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
					norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
					ffx = -1. * cij / ep * xnor_init / norvec * pow((Dij_min - th) / th, 2);
					ffy = -1. * cij / ep * ynor_init / norvec * pow((Dij_min - th) / th, 2);
					particle[np].center[0].f3x += ffx;
					particle[np].center[0].f3y += ffy;
					particle[np].center[0].t3_z += ffx * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy * (particle[np].node[node_1].x - particle[np].center[0].x);
				}
				if (Dij_min <= min_th){
					angle_1 = atan((particle[np].node[0].y - particle[np].center[0].y) / (particle[np].node[0].x - particle[np].center[0].x));//It gives us the angle of particle np
					xnor = (2 * (particle[np].node[node_1].x - particle[np].center[0].x)*cos(angle_1) - 2 * (particle[np].node[node_1].y - particle[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
					ynor = (2 * (particle[np].node[node_1].x - particle[np].center[0].x)*sin(angle_1) + 2 * (particle[np].node[node_1].y - particle[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
					xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
					ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
					norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
					ffx = -1. * cij / Ep * xnor_init / norvec * pow((Dij_min - th) / th, 2);
					ffy = -1. * cij / Ep * ynor_init / norvec * pow((Dij_min - th) / th, 2);
					particle[np].center[0].f3x += ffx;
					particle[np].center[0].f3y += ffy;
					particle[np].center[0].t3_z += ffx * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy * (particle[np].node[node_1].x - particle[np].center[0].x);
				}
			}
		}
	}// END of Particle-Particle Collision

	double Rij, Rij_min;;
#if defined Y_NOT_Periodic
	if ((particle[np].center[0].y <= MAX(particle_radius, particle_radius_2) + thw) || ((Ny - particle[np].center[0].y) <= MAX(particle_radius, particle_radius_2) + thw)){
		node_1 = 0;
		Rij_min = thw + 1;//Just to initially have a value greater than th
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			if (particle[np].node[n].y <= (Ny - 1) / 2.){//Particle-Y-Wall Collision section
				Rij = abs(particle[np].node[n].y - 0);// 0 means y-location of lower wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			else
			{
				Rij = -1.0 * abs((Ny - 1) - particle[np].node[n].y);// Ny-1 means y-location of upper wall	.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			if (abs(Rij) <= abs(Rij_min)){
				Rij_min = Rij;
				node_1 = n;
			}
		}
		if (abs(Rij_min) > min_thw && abs(Rij_min) <= thw){
			ffx4 = 0.0;
			ffy4 = 1. * cij / epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
		if (abs(Rij_min) <= min_thw){
			ffx4 = 0.0;
			ffy4 = 1. * cij / Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
	}
#endif
#if defined X_NOT_Periodic
	if ((particle[np].center[0].x <= MAX(particle_radius, particle_radius_2) + thw) || ((Nx - particle[np].center[0].x) <= MAX(particle_radius, particle_radius_2) + thw)){
		node_1 = 0;
		Rij_min = thw + 1;//Just to initially have a value greater than th
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			if (particle[np].node[n].x <= (Nx - 1) / 2.){//Particle-Y-Wall Collision section
				Rij = abs(particle[np].node[n].x - 0);// 0 means y-location of lower wall.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			else
			{
				Rij = -1.0 * abs((Nx - 1) - particle[np].node[n].x);// Ny-1 means y-location of upper wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			if (abs(Rij) <= abs(Rij_min)){
				Rij_min = Rij;
				node_1 = n;
			}
		}
		if (abs(Rij_min) > min_thw && abs(Rij_min) <= thw){
			ffx4 = 1. * cij / epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			ffy4 = 0.0;
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
		if (abs(Rij_min) <= min_thw){
			ffx4 = 1. * cij / Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			ffy4 = 0.0;
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
	}
#endif// END of Particle_X-Wall Collision

#endif
#if defined Ellipse_Collision_3
	//Collision Forces. JUST for ELLIPSE shape
	const double min_th = 0.0;// If distance is less than this value, Ep is used instead of ep
	const double min_thw = 0.0;// If distance is less than this value, Epw is used instead of epw
	double Dij, Dist, ffx, ffy, ffx4, ffy4, angle_2, fas;
	ffx = 0.0;
	ffy = 0.0;
	int node_1;
	int node_2;
	double Dij_min;//Just to initially have a value greater than th

	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dist = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dist <= (2 * MAX(particle_radius, particle_radius_2) + th)){
				node_1 = 0;
				node_2 = 0;
				Dij_min = th + 1;//Just to initially have a value greater than th
				for (int n = 0; n < particle[np].num_nodes; ++n) {
					for (int m = 0; m < particle[np2].num_nodes; ++m) {
						Dij = sqrt(SQ(particle[np].node[n].x - particle[np2].node[m].x) + SQ(particle[np].node[n].y - particle[np2].node[m].y));
						if (Dij <= Dij_min) {
							Dij_min = Dij;
							node_1 = n;
							node_2 = m;
						}
					}
				}

				if (Dij_min >min_th && Dij_min <= th){
					ffx = 1. * cij / ep * (particle[np].node[node_1].x - particle[np2].node[node_2].x) / Dij_min * pow((Dij_min - th) / th, 2);
					ffy = 1. * cij / ep * (particle[np].node[node_1].y - particle[np2].node[node_2].y) / Dij_min * pow((Dij_min - th) / th, 2);
					particle[np].center[0].f3x += ffx;
					particle[np].center[0].f3y += ffy;
					particle[np].center[0].t3_z += ffx * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy * (particle[np].node[node_1].x - particle[np].center[0].x);
				}
				if (Dij_min <= min_th){
					ffx = 1. * cij / Ep * (particle[np].node[node_1].x - particle[np2].node[node_2].x) / Dij_min * pow((Dij_min - th) / th, 2);
					ffy = 1. * cij / Ep * (particle[np].node[node_1].y - particle[np2].node[node_2].y) / Dij_min * pow((Dij_min - th) / th, 2);
					particle[np].center[0].f3x += ffx;
					particle[np].center[0].f3y += ffy;
					particle[np].center[0].t3_z += ffx * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy * (particle[np].node[node_1].x - particle[np].center[0].x);
				}
			}
		}
	}// END of Particle-Particle Collision

	double Rij, Rij_min;;
#if defined Y_NOT_Periodic
	if ((particle[np].center[0].y <= MAX(particle_radius, particle_radius_2) + thw) || ((Ny - particle[np].center[0].y) <= MAX(particle_radius, particle_radius_2) + thw)){
		node_1 = 0;
		Rij_min = thw + 1;//Just to initially have a value greater than th
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			if (particle[np].node[n].y <= (Ny - 1) / 2.){//Particle-Y-Wall Collision section
				Rij = abs(particle[np].node[n].y - 0);// 0 means y-location of lower wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			else
			{
				Rij = -1.0 * abs((Ny - 1) - particle[np].node[n].y);// Ny-1 means y-location of upper wall	.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			if (abs(Rij) <= abs(Rij_min)){
				Rij_min = Rij;
				node_1 = n;
			}
		}
		if (abs(Rij_min) > min_thw && abs(Rij_min) <= thw){
			ffx4 = 0.0;
			ffy4 = 1. * cij / epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
		if (abs(Rij_min) <= min_thw){
			ffx4 = 0.0;
			ffy4 = 1. * cij / Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
	}
#endif
#if defined X_NOT_Periodic
	if ((particle[np].center[0].x <= MAX(particle_radius, particle_radius_2) + thw) || ((Nx - particle[np].center[0].x) <= MAX(particle_radius, particle_radius_2) + thw)){
		node_1 = 0;
		Rij_min = thw + 1;//Just to initially have a value greater than th
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			if (particle[np].node[n].x <= (Nx - 1) / 2.){//Particle-Y-Wall Collision section
				Rij = abs(particle[np].node[n].x - 0);// 0 means y-location of lower wall.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			else
			{
				Rij = -1.0 * abs((Nx - 1) - particle[np].node[n].x);// Ny-1 means y-location of upper wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
			}
			if (abs(Rij) <= abs(Rij_min)){
				Rij_min = Rij;
				node_1 = n;
			}
		}
		if (abs(Rij_min) > min_thw && abs(Rij_min) <= thw){
			ffx4 = 1. * cij / epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			ffy4 = 0.0;
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
		if (abs(Rij_min) <= min_thw){
			ffx4 = 1. * cij / Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + thw) / thw, 2);
			ffy4 = 0.0;
			particle[np].center[0].f4x += ffx4;
			particle[np].center[0].f4y += ffy4;
			particle[np].center[0].t3_z += ffx4 * (particle[np].node[node_1].y - particle[np].center[0].y) - ffy4 * (particle[np].node[node_1].x - particle[np].center[0].x);
		}
	}
#endif// END of Particle_X-Wall Collision

#endif //End of Ellipse_Collision_3

#if defined Ellipse_Collision_4
	//Collision Forces. JUST for ELLIPSE shape
	double Dij, Dist, ffx, ffy, ffx4, ffy4;
	ffx = 0.0;
	ffy = 0.0;

	for (int np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
		if (np2 != np) {
			Dist = sqrt(SQ(particle[np].center[0].x - particle[np2].center[0].x) + SQ(particle[np].center[0].y - particle[np2].center[0].y));
			if (Dist <= (2 * MAX(particle_radius, particle_radius_2) + th + 2)){
				for (int n = 0; n < particle[np].num_nodes; ++n) {
					for (int m = 0; m < particle[np2].num_nodes; ++m) {
						Dij = sqrt(SQ(particle[np].node[n].x - particle[np2].node[m].x) + SQ(particle[np].node[n].y - particle[np2].node[m].y));
						if (Dij <= th) {
							ffx = 1. * cij / ep * (particle[np].node[n].x - particle[np2].node[m].x) / Dij * pow((Dij - th) / th, 2);
							ffy = 1. * cij / ep * (particle[np].node[n].y - particle[np2].node[m].y) / Dij * pow((Dij - th) / th, 2);
							particle[np].center[0].f3x += ffx;
							particle[np].center[0].f3y += ffy;
							particle[np].center[0].t3_z += ffx * (particle[np].node[n].y - particle[np].center[0].y) - ffy * (particle[np].node[n].x - particle[np].center[0].x);
						}


					}
				}
			}
		}
	}// END of Particle-Particle Collision

	double Rij;
#if defined Y_NOT_Periodic
	if ((particle[np].center[0].y <= MAX(particle_radius, particle_radius_2) + thw + 2) || ((Ny - particle[np].center[0].y) <= MAX(particle_radius, particle_radius_2) + thw + 2)){
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			if (particle[np].node[n].y <= (Ny - 1) / 2.){//Particle-Y-Wall Collision section
				Rij = particle[np].node[n].y - 0;// 0 means y-location of lower wall
			}
			else
			{
				Rij = -1.0 * ((Ny - 1) - particle[np].node[n].y);// Ny-1 means y-location of upper wall
			}
			if (abs(Rij) <= (thw)){
				ffx4 = 0.0;
				ffy4 = 1. * cij / epw * Rij / abs(Rij) * pow((-abs(Rij) + thw) / thw, 2);
				particle[np].center[0].f4x += ffx4;
				particle[np].center[0].f4y += ffy4;
				particle[np].center[0].t3_z += ffx4 * (particle[np].node[n].y - particle[np].center[0].y) - ffy4 * (particle[np].node[n].x - particle[np].center[0].x);
			}
		}
	}
#endif
#if defined X_NOT_Periodic
	if ((particle[np].center[0].x <= MAX(particle_radius, particle_radius_2) + thw + 2) || ((Nx - particle[np].center[0].x) <= MAX(particle_radius, particle_radius_2) + thw + 2)){
		for (int n = 0; n < particle[np].num_nodes; ++n) {
			if (particle[np].node[n].x <= (Nx - 1) / 2.){//Particle-X-Wall Collision section
				Rij = particle[np].node[n].x - 0;// 0 means x-location of left wall
			}
			else
			{
				Rij = -1.0 * ((Nx - 1) - particle[np].node[n].x);// Nx-1 means x-location of right wall
			}

			if (abs(Rij) <= (thw)){
				ffx4 = 1. * cij / epw * Rij / abs(Rij) * pow((-abs(Rij) + thw) / thw, 2);
				ffy4 = 0;
				particle[np].center[0].f4x += ffx4;
				particle[np].center[0].f4y += ffy4;
				particle[np].center[0].t3_z += ffx4 * (particle[np].node[n].y - particle[np].center[0].y) - ffy4 * (particle[np].node[n].x - particle[np].center[0].x);
			}
		}
	}
#endif	 // END of Particle-Wall Collision
#endif //end of Ellipse_Collision_4
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// The particle state (node positions) is writen to a VTK file at each t_vtk step and can be read by ParaView

void Particle_sim::write_particle_vtk(int time, particle_struct particle[num_particles]) {

#if defined VTK_ASCII || defined VTK_BINARY

	for (int npn = 0; npn < num_particles; ++npn) {

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
		output_file << "POINTS " << particle_num_nodes << " float\n";
		for (int n = 0; n < particle_num_nodes; ++n) {
			output_file << particle[npn].node[n].x << " " << particle[npn].node[n].y << " 0\n";
		}

		/// Write lines between neighboring nodes
		output_file << "LINES " << particle_num_nodes << " " << 3 * particle_num_nodes << "\n";
		for (int n = 0; n < particle_num_nodes; ++n) {
			output_file << "2 " << n << " " << (n + 1) % particle_num_nodes << "\n";
		}

		/// Write vertices
		output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
		output_file << particle_num_nodes << " ";
		for (int n = 0; n < particle_num_nodes; ++n) {
			output_file << n << " ";
		}

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
		for (int itm = t_vtk; itm <= t_num; itm = itm + t_vtk){
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
			if (particle.node[n].x<min_node_x)	{
				min_node_x = particle.node[n].x;
			}
			if (particle.node[n].y<min_node_y)	{
				min_node_y = particle.node[n].y;
			}
			if (particle.node[n].x>max_node_x)	{
				max_node_x = particle.node[n].x;
			}
			if (particle.node[n].y>max_node_y)	{
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
		for (int itm = 1; itm <= particle_num_nodes; itm++){
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int i = 0; i < particle_num_nodes; i++){
			output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << particle.node[i].x << "\n";
		}
		for (int i = 0; i < particle_num_nodes; i++){
			output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << particle.node[i].y << "\n";
		}
		for (int i = 0; i < particle_num_nodes; i++){
			output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << "\n";
		}
		output_file1 << "bar2\n";
		output_file1 << std::setw(10) << particle_num_nodes << "\n";
		for (int itm = 1; itm <= particle_num_nodes; itm++){
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int I = 1; I < particle_num_nodes; I++){
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
//
// The following PARTICLE quantities are written to the disk at each t_data step:
// - drag and lift forces (x- and y-components of the force)
// - object center position (x- and y-components)
// - object center velocity (x- and y-components)
// The data file is readable by gnuplot
void Particle_sim::write_particle_data(int time, particle_struct particle[num_particles]) {

	for (int npn = 0; npn < num_particles; ++npn) {

        /// Create filename
        ostringstream oss;
        oss << npn + 1;
        string npstr = oss.str();//convert number to string
        string output_filename("LB_Results/Data_particle/data_particle_");
        output_filename += npstr;
        output_filename += ".dat";
        //ofstream output_file;

            // Check if file already exists or is empty
        bool write_header = false;
        {
            std::ifstream infile(output_filename.c_str());
            if (!infile.good() || infile.peek() == std::ifstream::traits_type::eof()) {
                write_header = true; // file doesn’t exist or is empty
            }
        }

		/// Open file for append
        std::ofstream output_file(output_filename.c_str(), std::ios::app);

        if (write_header) {
            output_file << "#time, fx_surf, fy_surf, center.x, center.y, "
                           "vel_x, vel_y, omgZ, tetaZ, Q_1/4, Q_2/4, Q_3/4, Q_4/4, Q_tot\n";
        }

		/// Compute quantities

		double vel_center_x = 0;
		double vel_center_y = 0;

		/// Write data
		output_file << time << " "; // time step
		output_file << particle[npn].center[0].fx_surf << " "; // drag force
		output_file << particle[npn].center[0].fy_surf << " "; // lift force
        //output_file << particle[npn].center[0].x<< " "; // center position (x-component)
        output_file << particle_center_x + total_displacement_x<< " "; // center position (x-component): Activate for keep_particle_at_center instead of above line
		output_file << particle[npn].center[0].y << " "; // center position (y-component)
		output_file << particle[npn].center[0].vel_x << " "; // center velocity (x-component)
		output_file << particle[npn].center[0].vel_y << " "; // center velocity (y-component)
		output_file << particle[npn].center[0].omg << " "; // Rotational speed
		output_file << particle[npn].center[0].teta << " "; // Angle with respect to Horizontal line
        output_file << particle[npn].node[(int)(0 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
        output_file << particle[npn].node[(int)(1 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
        output_file << particle[npn].node[(int)(2 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
        output_file << particle[npn].node[(int)(3 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
        output_file << particle[npn].center[0].heat_surf << "\n"; //Total heat trasnferred
		/// Close file
		output_file.close();
	}
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::write_particle_recovery(int time, particle_struct particle[num_particles]) {

	for (int npn = 0; npn < num_particles; ++npn) {

	/// Create filename for particle data
    // create a stringstream and insert npn+1
    ostringstream oss;
    oss << (npn + 1);
    // get the string value
    string npstr = oss.str();
    string output_filename("LB_Results/data_particle_recovery_");
    output_filename += npstr;
    output_filename += ".dat";
    ofstream output_file(output_filename);
    /*string RRR = "HHHh";
    string npstr = (std::ostringstream() << npn + 1).str();//convert number to string
//    string npstr = static_cast<std::ostringstream*>(&(std::ostringstream() << npn + 1))->str();//convert number to string
    string output_filename("LB_Results/data_particle_recovery_");
	output_filename += npstr;
	output_filename += ".dat";
	ofstream output_file;*/

	/// Open file
	output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

	/// Write particle data
	output_file << std::fixed << std::setprecision(15) << particle_mass << " ";
	output_file << std::fixed << std::setprecision(15) << mom_inertia << " ";
	//for (int npn = 0; npn < num_particles; ++npn) {
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].x << " "; // center position (x-component)
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].y << " "; // center position (y-component)
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].vel_x << " "; // center velocity (x-component)
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].vel_y << " "; // center velocity (y-component)
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].omg << " "; // Rotational speed
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].teta << " "; // Angle with respect to Horizontal line
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].temperature << " "; // Temperature
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].previous_vel_x << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].previous_vel_y << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].tempx[0] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].tempx[1] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].tempy[0] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].tempy[1] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].temp_omg[0] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].temp_omg[1] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].temp_part << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].temp_temperat[0] << " ";
		output_file << std::fixed << std::setprecision(15) << particle[npn].center[0].temp_temperat[1] << " ";
		for (int np1 = 0; np1 < particle[npn].num_nodes; ++np1) {
			output_file << std::fixed << std::setprecision(15) << particle[npn].node[np1].x << " ";
			output_file << std::fixed << std::setprecision(15) << particle[npn].node[np1].y << " ";
			output_file << std::fixed << std::setprecision(15) << area[np1] << " ";
		}
	/// Close particle recovery data file
	output_file.close();
	}

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: READ_PARTICLE_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
void Particle_sim::read_particle_recovery(int time, particle_struct particle[num_particles]) {

	for (int npn = 0; npn < num_particles; ++npn) {

        // Convert integer to string
        ostringstream oss;
        oss << (npn + 1);
        string npstr = oss.str();
        // Build filename
        string output_filename("LB_Results/data_particle_recovery_");
        output_filename += npstr;
        output_filename += ".dat";
        // Open file for reading
        ifstream output_file(output_filename);
        /*string npstr = (ostringstream() << npn + 1).str();//convert number to string
        string output_filename("LB_Results/data_particle_recovery_");
		output_filename += npstr;
		output_filename += ".dat";
		ifstream output_file;*/

		/// Open file
		output_file.open(output_filename.c_str());
		/// Read particle data
		output_file >> particle_mass;
		output_file >> mom_inertia;
			/// Read data
			output_file >> particle[npn].center[0].x; // center position (x-component)
			output_file >> particle[npn].center[0].y; // center position (y-component)
			output_file >> particle[npn].center[0].vel_x; // center velocity (x-component)
			output_file >> particle[npn].center[0].vel_y; // center velocity (y-component)
			output_file >> particle[npn].center[0].omg; // Rotational speed
			output_file >> particle[npn].center[0].teta; // Angle with respect to Horizontal line
			output_file >> particle[npn].center[0].temperature; // Angle with respect to Horizontal line
			output_file >> particle[npn].center[0].previous_vel_x;
			output_file >> particle[npn].center[0].previous_vel_y;
			output_file >> particle[npn].center[0].tempx[0];
			output_file >> particle[npn].center[0].tempx[1];
			output_file >> particle[npn].center[0].tempy[0];
			output_file >> particle[npn].center[0].tempy[1];
			output_file >> particle[npn].center[0].temp_omg[0];
			output_file >> particle[npn].center[0].temp_omg[1];
			output_file >> particle[npn].center[0].temp_part;
			output_file >> particle[npn].center[0].temp_temperat[0];
			output_file >> particle[npn].center[0].temp_temperat[1];
			for (int np1 = 0; np1 < particle[npn].num_nodes; ++np1) {
				output_file >> particle[npn].node[np1].x;
				output_file >> particle[npn].node[np1].y;
				output_file >> area[np1];
			}
		/// Close particle recovery data file
		output_file.close();
	}
	return;
}
//=============================================================================================================================
//=============================================================================================================================
Particle_sim::~Particle_sim()
{
}

