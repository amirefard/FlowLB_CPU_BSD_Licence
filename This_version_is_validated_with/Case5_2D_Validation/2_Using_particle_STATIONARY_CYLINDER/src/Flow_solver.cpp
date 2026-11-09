#include "Flow_solver.h"



using namespace std;

Flow_solver::Flow_solver()
{
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: INITIALIZE
//-----------------------------------------------------------------------------------------------------------------------------
//
// The memory for lattice variables (populations, density, velocity, forces) is allocated.
// The variables are initialized.

void Flow_solver::initialize() {

    this->cx_pop_F  = vector<int>(stencil_directions_F, 0);
    this->cy_pop_F  = vector<int>(stencil_directions_F, 0);
    this->bounce_back_pop_F = vector<int>(stencil_directions_F, 0);

    if (stencil_directions_F == 9)
    {
        /// D2Q9 stencil:
        // index:   0  1  2  3  4  5  6  7  8
        // x:       0  1  0 -1  0  1 -1 -1  1
        // y:       0  0  1  0 -1  1  1 -1 -1
        // z:       0  0  0  0  0  0  0  0  0

        weight_F = {4. / 9.,
                    1. / 9., 1. / 9., 1. / 9., 1. / 9.,
                    1. / 36., 1. / 36., 1. / 36., 1. / 36.};

        cx_pop_F = {0, 1, 0, -1, 0, 1, -1, -1,  1};
        cy_pop_F = {0, 0, 1, 0, -1, 1,  1, -1, -1};
        bounce_back_pop_F = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    }

        for (int alpha = 0; alpha < stencil_directions_F; alpha++)
    {
        if (cx_pop_F[alpha] > 0 )
        {
            pop_to_right_F.push_back(alpha);
            pop_reverse_right_F.push_back(bounce_back_pop_F[alpha]);
        }
        if (cx_pop_F[alpha] < 0 )
        {
            pop_to_left_F.push_back(alpha);
            pop_reverse_left_F.push_back(bounce_back_pop_F[alpha]);
        }
        if (cy_pop_F[alpha] > 0)
        {
            pop_to_top_F.push_back(alpha);
            pop_reverse_top_F.push_back(bounce_back_pop_F[alpha]);
        //    cout << "AA " <<pop_to_top_F[pop_to_top_F.size() - 1] << " " << pop_reverse_top_F[pop_reverse_top_F.size() - 1] << endl;
        }
        if (cy_pop_F[alpha] < 0)
        {
            pop_to_bottom_F.push_back(alpha);
            pop_reverse_bottom_F.push_back(bounce_back_pop_F[alpha]);
        }
    }

#if defined __linux__ || defined __APPLE__
    mkdir("LB_Results", 0777); /* Or what parameter you need here ... */
    mkdir("LB_Results/vtk_fluid", 0777); /* Or what parameter you need here ... */
    mkdir("LB_Results/vtk_particle", 0777); /* Or what parameter you need here ... */
    mkdir("LB_Results/Data_particle", 0777); /* Or what parameter you need here ... */
#else
    _mkdir("LB_Results");
    _mkdir("LB_Results/vtk_fluid");
    _mkdir("LB_Results/vtk_particle");
    _mkdir("LB_Results/Data_particle");

#endif

#if defined Recovery_OFF
    system("rm -f LB_Results/vtk_fluid/*.*");
    system("rm -f LB_Results/vtk_particle/*.*");
    system("rm -f LB_Results/Data_particle/*.*");
    system("rm -f LB_Results/*.*");

	/*std::string command1 = "del /Q ";
    std::string path1 = "LB_Results\\vtk_fluid\\*.*";
	system(command1.append(path1).c_str());

	std::string command2 = "del /Q ";
    std::string path2 = "LB_Results\\vtk_particle\\*.*";
	system(command2.append(path2).c_str());

	std::string command3 = "del /Q ";
    std::string path3 = "LB_Results\\Data_particle\\*.*";
	system(command3.append(path3).c_str());

	std::string command4 = "del /Q ";
    std::string path4 = "LB_Results\\*.*";
	system(command4.append(path4).c_str());*/
#endif

	/// Allocate memory for the fluid density, velocity, and force
	density = new double*[Nx];
	velocity_x = new double*[Nx];
	velocity_y = new double*[Nx];
	force_x = new double*[Nx];
	force_y = new double*[Nx];
	is_solid = new int*[Nx];
	mdf_vel_x = new double*[Nx];
	mdf_vel_y = new double*[Nx];
    mdf_temperature = new double*[Nx];
	mdf_force_x = new double*[Nx];
	mdf_force_y = new double*[Nx];
    mdf_force_thermal = new double*[Nx];
	velocity_magnitude = new double*[Nx];
	previous_velocity_magnitude = new double*[Nx];

	boundary_cond = new int*[Nx];
	boundary_vel_x = new double*[Nx];
	boundary_vel_y = new double*[Nx];
	boundary_pres = new double*[Nx];
	boundary_cond_2 = new int*[Nx];
	edge_boundary = new int*[Nx];
	boundary_vel_x_2 = new double*[Nx];
	boundary_vel_y_2 = new double*[Nx];
	boundary_pres_2 = new double*[Nx];

	for (int X = 0; X < Nx; ++X) {
		density[X] = new double[Ny];
		velocity_x[X] = new double[Ny];
		velocity_y[X] = new double[Ny];
		force_x[X] = new double[Ny];
		force_y[X] = new double[Ny];
		is_solid[X] = new int[Ny];
		mdf_vel_x[X] = new double[Ny];
		mdf_vel_y[X] = new double[Ny];
        mdf_temperature[X] = new double[Ny];
		mdf_force_x[X] = new double[Ny];
		mdf_force_y[X] = new double[Ny];
        mdf_force_thermal[X] = new double[Ny];
		previous_velocity_magnitude[X] = new double[Ny];
		velocity_magnitude[X] = new double[Ny];

		boundary_cond[X] = new int[Ny];
		boundary_vel_x[X] = new double[Ny];
		boundary_vel_y[X] = new double[Ny];
		boundary_pres[X] = new double[Ny];
		boundary_cond_2[X] = new int[Ny];
		edge_boundary[X] = new int[Ny];
		boundary_vel_x_2[X] = new double[Ny];
		boundary_vel_y_2[X] = new double[Ny];
		boundary_pres_2[X] = new double[Ny];
	}


	/// Initialize the fluid density and velocity
	// Start with unit density and zero velocity.
//	rec_step = 0;
	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			density[X][Y] = 1.0;
			velocity_x[X][Y] = Fluid_init_vel_x;
			velocity_y[X][Y] = Fluid_init_vel_y;
			force_x[X][Y] = 0;
			force_y[X][Y] = 0;
			is_solid[X][Y] = FALSE;
			mdf_vel_x[X][Y] = 0.0;
			mdf_vel_y[X][Y] = 0.0;
            mdf_temperature[X][Y] = 0.0;
			mdf_force_x[X][Y] = 0;
			mdf_force_y[X][Y] = 0;
            mdf_force_thermal[X][Y] = 0.0;
			velocity_magnitude[X][Y] = 0.0;
			previous_velocity_magnitude[X][Y] = 0.0;

			boundary_cond[X][Y] = 0;//0 means fluid
			boundary_vel_x[X][Y] = 0;//0 means fluid
			boundary_vel_y[X][Y] = 0;//0 means fluid
			boundary_pres[X][Y] = 0;//0 means fluid
			boundary_cond_2[X][Y] = 0;
			edge_boundary[X][Y] = FALSE;
			boundary_vel_x_2[X][Y] = 0;//0 means fluid
			boundary_vel_y_2[X][Y] = 0;//0 means fluid
			boundary_pres_2[X][Y] = 0;//0 means fluid

			if ((X == 1) || (X == Nx - 2))
			{
#if defined X_NOT_Periodic
				boundary_cond[X][Y] = 7;//7 means solid on x boundary
				boundary_cond_2[X][Y] = 7;
#endif
			}

			if ((Y == 1) || (Y == Ny - 2))
			{
#if defined Y_NOT_Periodic
				boundary_cond[X][Y] = 8;//8 means solid on y boundary
				boundary_cond_2[X][Y] = 8;
#endif
			}

		}
	}

	/// Allocate memory for the populations

	pop_F = new double**[stencil_directions_F];
	pop_old_F = new double**[stencil_directions_F];
	delta = new double**[stencil_directions_F];

	for (int c_i = 0; c_i < stencil_directions_F; ++c_i) {
		pop_F[c_i] = new double*[Nx];
		pop_old_F[c_i] = new double*[Nx];
		delta[c_i] = new double*[Nx];
		for (int X = 0; X < Nx; ++X) {
			pop_F[c_i][X] = new double[Ny];
			pop_old_F[c_i][X] = new double[Ny];
			delta[c_i][X] = new double[Ny];
			for (int Y = 0; Y < Ny; ++Y) {
				pop_F[c_i][X][Y] = 0;
				pop_old_F[c_i][X][Y] = 0;
				delta[c_i][X][Y] = 0;
			}
		}
	}


	//#if defined Recovery_OFF
	/// Initialize the populations
	// Use the equilibrium populations corresponding to the initialized fluid density and velocity.
	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
			for (int c_i = 0; c_i < stencil_directions_F; ++c_i) {
				pop_old_F[c_i][X][Y] = pop_eq_F[c_i];
				pop_F[c_i][X][Y] = pop_eq_F[c_i];
			}
		}
	}
	//#endif

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

void Flow_solver::equilibrium(double den, double vel_x, double vel_y) {
    for (size_t c_i = 0; c_i < stencil_directions_F; ++c_i ){
        pop_eq_F[c_i] = weight_F[c_i] * den * (1 + 3.0 * (cx_pop_F[c_i] * vel_x + cy_pop_F[c_i] * vel_y)
                                             + 4.5 * SQ(cx_pop_F[c_i] * vel_x + cy_pop_F[c_i] * vel_y)
                                             - 1.5 * (SQ(vel_x) + SQ(vel_y)));
    }
	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: calculate_gbetaindex
//-----------------------------------------------------------------------------------------------------------------------------
int Flow_solver::calculate_gbetaindex(int X, int Y, const Particle_sim& Particle_Model){

    int gbindex = 1;//Will be 0 outside particle and 1 inside particle
    for (int npa = 0; npa < num_particles; npa++) {
#if defined STATIONARY_CYLINDER || defined MOVING_CYLINDER
        if (sqrt(SQ(X - Particle_Model.particle[npa].center[0].x) + SQ(Y - Particle_Model.particle[npa].center[0].y)) <= particle_radius) {
            gbindex = 0;
        }
#endif
#if defined STATIONARY_ELLIPSE || defined MOVING_ELLIPSE
        double alphaAngle = Particle_Model.particle[npa].center[0].teta - M_PI / 2.;
        if ((SQ((X - Particle_Model.particle[npa].center[0].x) *  cos(alphaAngle) + (Y - Particle_Model.particle[npa].center[0].y) * sin(alphaAngle)) / SQ(particle_radius) + SQ((X - Particle_Model.particle[npa].center[0].x) * sin(alphaAngle) - (Y - Particle_Model.particle[npa].center[0].y) * cos(alphaAngle)) / SQ(particle_radius_2)) <= 1) {
            gbindex = 0;
        }
        //cout << alphaAngle << endl;
        //if ((X + start_x[processor_id] - 2 == 110) && (Y == 1023)) { { cout << gbetaindex<<endl;  } }
#endif
    }
    return gbindex;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM
//-----------------------------------------------------------------------------------------------------------------------------
//
void Flow_solver::LBM(int time, const Particle_sim& Particle_Model) {

	/// Swap populations
	// The present code used old and new populations which are swapped at the beginning of each time step.
	// This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
	// This way, the old populations are not overwritten during propagation.
	// The resulting code is easier to write and to debug.
	// The memory requirement for the populations is twice as large.

	//double s[9];
	//s[7]=s[8]=1.0/tau; s[0]=s[3]=s[5]=0.0; s[4]=s[6]=8.*(2.-s[7])/(8.-s[7]); s[1]=1.6; s[2]=1.8; // relaxation rates for MRT Page 28 Guo BOOK
	double ssf = 1. / tau;
	double s[9] = { 1.0, 1.4, 1.4, 1.0, 1.2, 1.0, 1.2, ssf, ssf };//As sugested in page 103 of Mohamad book
	//double s[9] = {ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf, ssf};//If all coefficients are 1/tau then we have SRT
	double MEQ;
	double m[9];
	double mforce[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double fforce[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double ***swap_temp = pop_old_F;
	pop_old_F = pop_F;
	pop_F = swap_temp;
	double velv_x = 0.0;
	double velv_y = 0.0;

	/// Lattice Boltzmann equation
	// The lattice Boltzmann equation is solved in the following.
	// The algorithm includes
	// - computation of the lattice force
	// - combined collision and propagation (faster than first collision and then propagation)

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

#if defined X_NOT_Periodic
	start_X = 1;
	end_X = Nx - 2;
#endif

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

#if defined Y_NOT_Periodic
	start_Y = 1;
	end_Y = Ny - 2;
#endif

int gbetaindex = 1; //Inside particle =0, outside particle = 1
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE){

                #if defined Flow_With_Resolved_Particle && defined Flow_With_Thermal_Effect
                gbetaindex = calculate_gbetaindex(X, Y, Particle_Model);// To set gbetaindex=0 for fluid nodes inside the particle
                #endif

				/// Compute new populations
				// This is the lattice Boltzmann equation (combined collision and propagation) including external forcing.
				// Periodicity of the lattice in x-direction is taken into account by the %-operator.
///===================================================================================================================================== ///
///                                                 SHAN-CHEN MODEL                                                                      ///
///===================================================================================================================================== ///
#if defined Shan_Chen
                double forceX = force_x[X][Y] + gbetaindex * gbeta_x * (temperature[X][Y] - Tref) + gravity;
                double forceY = force_y[X][Y] + gbetaindex * gbeta_y * (temperature[X][Y] - Tref);
                double v0 = velocity_x[X][Y] + forceX * tau;
                double v1 = velocity_y[X][Y] + forceY * tau;
                double den = density[X][Y];
                equilibrium(den, v0, v1);
                for (int i = 0; i < stencil_directions_F; i++)
                {
                    int newX = (X + cx_pop_F[i] + Nx) % Nx;
                    int newY = (Y + cy_pop_F[i] + Ny) % Ny;
                    pop_F[i][newX][newY] = pop_old_F[i][X][Y] * (1.0 - omega) + pop_eq_F[i] * omega;
                }
#endif
///===================================================================================================================================== ///
///                                                    LUO MODEL                                                                         ///
///===================================================================================================================================== ///
#if defined Luo
                equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
                for (int i = 0; i < stencil_directions_F; i++)
                {
                    int newX = (X + cx_pop_F[i] + Nx) % Nx;
                    int newY = (Y + cy_pop_F[i] + Ny) % Ny;
                    double forceX = force_x[X][Y] + gbetaindex * gbeta_x * (temperature[X][Y] - Tref) + gravity;
                    double forceY = force_y[X][Y] + gbetaindex * gbeta_y * (temperature[X][Y] - Tref);
                    pop_F[i][newX][newY] = pop_old_F[i][X][Y] * (1.0 - omega) + pop_eq_F[i] * omega
                                       + 3 * weight_F[i] * (cx_pop_F[i] * forceX + cy_pop_F[i] * forceY);// * density[X][Y];
                }
#endif
///===================================================================================================================================== ///
///                                                    GUO MODEL                                                                         ///
///===================================================================================================================================== ///
#if defined Guo
#if defined SRT
                equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);

                for (int i = 0; i < stencil_directions_F; i++)
                {
                    int newX = (X + cx_pop_F[i] + Nx) % Nx;
                    int newY = (Y + cy_pop_F[i] + Ny) % Ny;
                    double forceX = force_x[X][Y] + gbetaindex * gbeta_x * (temperature[X][Y] - Tref) + gravity;
                    double forceY = force_y[X][Y] + gbetaindex * gbeta_y * (temperature[X][Y] - Tref);
                    pop_F[i][newX][newY] = pop_old_F[i][X][Y] * (1.0 - omega) + pop_eq_F[i] * omega
                                       + (1.0 - 0.5 * omega) * weight_F[i] * (forceX * (3 * (cx_pop_F[i] - velocity_x[X][Y]) + 9.0 * cx_pop_F[i] * (cx_pop_F[i] * velocity_x[X][Y] + cy_pop_F[i] * velocity_y[X][Y]))
                                               + forceY * (3 * (cy_pop_F[i] - velocity_y[X][Y]) + 9.0 * cy_pop_F[i] * (cx_pop_F[i] * velocity_x[X][Y] + cy_pop_F[i] * velocity_y[X][Y])));// * density[X][Y];
                }

#endif
#if defined MRT
                double forceX = force_x[X][Y] + gbetaindex * gbeta_x * (temperature[X][Y] - Tref) + gravity;
                double forceY = force_y[X][Y] + gbetaindex * gbeta_y * (temperature[X][Y] - Tref);
				// Transformation from velocity space to moment space:
				m[0] = pop_old_F[0][X][Y] + pop_old_F[1][X][Y] + pop_old_F[2][X][Y] + pop_old_F[3][X][Y] + pop_old_F[4][X][Y] + pop_old_F[5][X][Y] + pop_old_F[6][X][Y] + pop_old_F[7][X][Y] + pop_old_F[8][X][Y];
				m[1] = -4 * pop_old_F[0][X][Y] - pop_old_F[1][X][Y] - pop_old_F[2][X][Y] - pop_old_F[3][X][Y] - pop_old_F[4][X][Y] + 2 * (pop_old_F[5][X][Y] + pop_old_F[6][X][Y] + pop_old_F[7][X][Y] + pop_old_F[8][X][Y]);
				m[2] = 4 * pop_old_F[0][X][Y] - 2 * (pop_old_F[1][X][Y] + pop_old_F[2][X][Y] + pop_old_F[3][X][Y] + pop_old_F[4][X][Y]) + pop_old_F[5][X][Y] + pop_old_F[6][X][Y] + pop_old_F[7][X][Y] + pop_old_F[8][X][Y];
				m[3] = pop_old_F[1][X][Y] - pop_old_F[3][X][Y] + pop_old_F[5][X][Y] - pop_old_F[6][X][Y] - pop_old_F[7][X][Y] + pop_old_F[8][X][Y];
				m[4] = -2 * (pop_old_F[1][X][Y] - pop_old_F[3][X][Y]) + pop_old_F[5][X][Y] - pop_old_F[6][X][Y] - pop_old_F[7][X][Y] + pop_old_F[8][X][Y];
				m[5] = pop_old_F[2][X][Y] - pop_old_F[4][X][Y] + pop_old_F[5][X][Y] + pop_old_F[6][X][Y] - pop_old_F[7][X][Y] - pop_old_F[8][X][Y];
				m[6] = -2 * (pop_old_F[2][X][Y] - pop_old_F[4][X][Y]) + pop_old_F[5][X][Y] + pop_old_F[6][X][Y] - pop_old_F[7][X][Y] - pop_old_F[8][X][Y];
				m[7] = pop_old_F[1][X][Y] - pop_old_F[2][X][Y] + pop_old_F[3][X][Y] - pop_old_F[4][X][Y];
				m[8] = pop_old_F[5][X][Y] - pop_old_F[6][X][Y] + pop_old_F[7][X][Y] - pop_old_F[8][X][Y];

				fforce[0] = weight_F[0] * (forceX*(3 * (0 - velocity_x[X][Y])) + forceY * (3 * (0 - velocity_y[X][Y])));//*density[X][Y];
				fforce[1] = weight_F[1] * (forceX*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y])) + forceY * (3 * (0 - velocity_y[X][Y])));//*density[X][Y];
				fforce[2] = weight_F[2] * (forceX*(3 * (0 - velocity_x[X][Y])) + forceY * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (1 * velocity_y[X][Y])));//*density[X][Y];
				fforce[3] = weight_F[3] * (forceX*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y])) + forceY * (3 * (0 - velocity_y[X][Y])));//*density[X][Y];
				fforce[4] = weight_F[4] * (forceX*(3 * (0 - velocity_x[X][Y])) + forceY * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (-1 * velocity_y[X][Y])));//*density[X][Y];
				fforce[5] = weight_F[5] * (forceX*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))) + forceY * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))));//*density[X][Y];
				fforce[6] = weight_F[6] * (forceX*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))) + forceY * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (-1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))));//*density[X][Y];
				fforce[7] = weight_F[7] * (forceX*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))) + forceY * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))));//*density[X][Y];
				fforce[8] = weight_F[8] * (forceX*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))) + forceY * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))));//*density[X][Y];

				mforce[0] = fforce[0] + fforce[1] + fforce[2] + fforce[3] + fforce[4] + fforce[5] + fforce[6] + fforce[7] + fforce[8];
				mforce[1] = -4 * fforce[0] - fforce[1] - fforce[2] - fforce[3] - fforce[4] + 2 * (fforce[5] + fforce[6] + fforce[7] + fforce[8]);
				mforce[2] = 4 * fforce[0] - 2 * (fforce[1] + fforce[2] + fforce[3] + fforce[4]) + fforce[5] + fforce[6] + fforce[7] + fforce[8];
				mforce[3] = fforce[1] - fforce[3] + fforce[5] - fforce[6] - fforce[7] + fforce[8];
				mforce[4] = -2 * (fforce[1] - fforce[3]) + fforce[5] - fforce[6] - fforce[7] + fforce[8];
				mforce[5] = fforce[2] - fforce[4] + fforce[5] + fforce[6] - fforce[7] - fforce[8];
				mforce[6] = -2 * (fforce[2] - fforce[4]) + fforce[5] + fforce[6] - fforce[7] - fforce[8];
				mforce[7] = fforce[1] - fforce[2] + fforce[3] - fforce[4];
				mforce[8] = fforce[5] - fforce[6] + fforce[7] - fforce[8];

				for (int k = 0; k < 9; k++)
				{
					MEQ = meq(density[X][Y], velocity_x[X][Y], velocity_y[X][Y], k);
					m[k] = m[k] - s[k] * (m[k] - MEQ) + (1 - 1. / 2.*s[k]) * mforce[k]; // relaxation
					m[k] /= 36.; // Page 102 Mohamad Book
					mforce[k] = 0.0;
					fforce[k] = 0.0;
				}

				pop_F[0][X][Y] = 4 * m[0] - 4 * (m[1] - m[2]);
				pop_F[1][(X + 1) % Nx][Y] = 4 * m[0] - m[1] - 2 * m[2] + 6 * m[3] - 6 * m[4] + 9 * m[7];
				pop_F[2][X][(Y + 1) % Ny] = 4 * m[0] - m[1] - 2 * m[2] + 6 * m[5] - 6 * m[6] - 9 * m[7];
				pop_F[3][(X - 1 + Nx) % Nx][Y] = 4 * m[0] - m[1] - 2 * m[2] - 6 * m[3] + 6 * m[4] + 9 * m[7];
				pop_F[4][X][(Y - 1 + Ny) % Ny] = 4 * m[0] - m[1] - 2 * m[2] - 6 * m[5] + 6 * m[6] - 9 * m[7];
				pop_F[5][(X + 1) % Nx][(Y + 1) % Ny] = 4 * m[0] + 2 * m[1] + m[2] + 6 * m[3] + 3 * m[4] + 6 * m[5] + 3 * m[6] + 9 * m[8];
				pop_F[6][(X - 1 + Nx) % Nx][(Y + 1) % Ny] = 4 * m[0] + 2 * m[1] + m[2] - 6 * m[3] - 3 * m[4] + 6 * m[5] + 3 * m[6] - 9 * m[8];
				pop_F[7][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny] = 4 * m[0] + 2 * m[1] + m[2] - 6 * m[3] - 3 * m[4] - 6 * m[5] - 3 * m[6] + 9 * m[8];
				pop_F[8][(X + 1) % Nx][(Y - 1 + Ny) % Ny] = 4 * m[0] + 2 * m[1] + m[2] + 6 * m[3] + 3 * m[4] - 6 * m[5] - 3 * m[6] - 9 * m[8];
#endif
#endif
			}
		}
	}

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: BC
//-----------------------------------------------------------------------------------------------------------------------------
void Flow_solver::BC(int time) {

    int newX, newY;

    for (int ii = 0; ii < 4; ii++)
    {
        int start_BC_X{}, end_BC_X{}, start_BC_Y{}, end_BC_Y{};
        int tag;
        std::vector<int> pop_to, pop_reverse;
        std::string BC_type;

        pop_to.clear();
        pop_reverse.clear();
        BC_type = "";

        switch(ii)
        {
        case 0:
            /// Right BC
                pop_to      = pop_to_right_F;
                pop_reverse = pop_reverse_right_F;
                start_BC_X  = Nx - 2;
                end_BC_X    = Nx - 2;
                start_BC_Y  = 1;
                end_BC_Y    = Ny - 2;
                tag = ii;
                //cout << "Right, processor:" << processor_id << " "<<start_BC_X<< "-"<<end_BC_X << ", "<< start_BC_Y<< "-"<<end_BC_Y << ", "<< start_BC_Z<< "-"<<end_BC_Z << endl;
                #if defined X_Wall
                    BC_type = "wall";
                #endif
        break;
        case 1:
            /// Left BC
                pop_to      = pop_to_left_F;
                pop_reverse = pop_reverse_left_F;
                start_BC_X  = 1;
                end_BC_X    = 1;
                start_BC_Y  = 1;
                end_BC_Y    = Ny - 2;
                tag = ii;
                #if defined X_Wall
                    BC_type = "wall";
                #endif
        break;
        case 2:
            /// Top BC
                pop_to      = pop_to_top_F;
                pop_reverse = pop_reverse_top_F;
                start_BC_Y  = Ny - 2;
                end_BC_Y    = Ny - 2;
                start_BC_X  = 1;
                end_BC_X    = Nx - 2;
                tag = ii;
                #if defined Y_Wall
                    BC_type = "wall";
                #endif
            break;
        case 3:
            /// Bottom BC
                pop_to      = pop_to_bottom_F;
                pop_reverse = pop_reverse_bottom_F;
                start_BC_Y  = 1;
                end_BC_Y    = 1;
                start_BC_X  = 1;
                end_BC_X    = Nx - 2;
                tag = ii;
                #if defined Y_Wall
                    BC_type = "wall";
                #endif
            break;
        default:
            std::cout << "Invalid choice" << std::endl;
            break;
        }

        for (int alpha = 0; alpha < pop_to.size(); alpha++)
        {
            for (int X = start_BC_X; X <= end_BC_X; ++X)
            {
                for (int Y = start_BC_Y; Y <= end_BC_Y; ++Y)
                {
                    if (BC_type == "wall")
                    {
                        int pop_id = pop_to[alpha];
                        newX = X + cx_pop_F[pop_id];
                        newY = Y + cy_pop_F[pop_id];
                        pop_F[pop_reverse[alpha]][X][Y] = pop_F[pop_id][newX][newY];
                        //cout << "BB " << pop_to[alpha] << " " << pop_reverse[alpha] <<endl;
                    }
                }
            }
        }
    }


///##############################################################################################///

#if defined Flow_Without_Particle
		//// Bounce back on SOLID nodes (Just for Flow_Without_Particle, because I have not defined collision between particles and solid objects.)
		double min_dist, dist;
		int i_min;
		int x1, y1, x2, y2;
		int Xp, Yp, Xp2, Yp2;
#if defined Curved_boundary_treatment
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				if (is_solid[X][Y] == FALSE){
					for (int dir = 1; dir < stencil_directions_F; dir++)
					{
						if (is_solid[(X - cx_pop_F[dir] + Nx) % Nx][(Y - cy_pop_F[dir] + Ny) % Ny] == TRUE){
							Xp  = (X - cx_pop_F[dir] + Nx) % Nx;
							Yp  = (Y - cy_pop_F[dir] + Ny) % Ny;
							Xp2 = (X + cx_pop_F[dir] + Nx) % Nx;
							Yp2 = (Y + cy_pop_F[dir] + Ny) % Ny;

							if (time == recovery_step + 1)	{ // To avoid calculating delta multiple times

								x1 = X;
								y1 = Y;
								x2 = X - cx_pop_F[dir];
								y2 = Y - cy_pop_F[dir];
								min_dist = 10000;

								for (size_t i = 0; i < surface_num_nodes; i++)
								{
									if (solid.node_surf[i].x <= (MAX(x1, x2) + 1) && solid.node_surf[i].x >= (MIN(x1, x2) - 1) && solid.node_surf[i].y <= (MAX(y1, y2) + 1) && solid.node_surf[i].y >= (MIN(y1, y2) - 1))
									{
										dist = (sqrt(SQ(solid.node_surf[i].x - x1) + SQ(solid.node_surf[i].y - y1)) + sqrt(SQ(solid.node_surf[i].x - x2) + SQ(solid.node_surf[i].y - y2)));
										if (dist<min_dist)
										{
											min_dist = dist;
											i_min = i;
										}
									}
								}
								delta[oppose[dir]][X][Y] = abs(sqrt(SQ(x1 - solid.node_surf[i_min].x) + SQ(y1 - solid.node_surf[i_min].y)) / sqrt(SQ(x1 - x2) + SQ(y1 - y2)));
							} //End of: if (time == recovery_step+1)

#if defined Bouzidi
							if (delta[oppose[dir]][X][Y] >= 0.5)
							{
								pop_F[dir][X][Y] = 1.0 / (2.0*delta[oppose[dir]][X][Y])*pop_F[oppose[dir]][Xp][Yp] + (2.0*delta[oppose[dir]][X][Y] - 1) / (2.0*delta[oppose[dir]][X][Y])*pop_F[dir][Xp2][Yp2];
							}
							else
							{
								pop_F[dir][X][Y] = 2.0*delta[oppose[dir]][X][Y] * pop_F[oppose[dir]][Xp][Yp] + (1 - 2.0*delta[oppose[dir]][X][Y]) * pop_F[oppose[dir]][X][Y];
							}
#endif
#if defined Mei_Luo
							double temp_pop, u_bf_x, u_bf_y, zet;
							if (delta[oppose[dir]][X][Y] >= 0.5)
							{
								u_bf_x = (delta[oppose[dir]][X][Y] - 1) * velocity_x[X][Y] / delta[oppose[dir]][X][Y];
								u_bf_y = (delta[oppose[dir]][X][Y] - 1) * velocity_y[X][Y] / delta[oppose[dir]][X][Y];
								zet = (2 * delta[oppose[dir]][X][Y] - 1.0) / (tau);
							}
							else
							{
								u_bf_x = velocity_x[Xp2][Yp2];
								u_bf_y = velocity_y[Xp2][Yp2];
								zet = (2 * delta[oppose[dir]][X][Y] - 1.0) / (tau - 2);
							}
							temp_pop = weight_F[oppose[dir]] * density[X][Y] * (1 + 3 * (cx_pop_F[oppose[dir]] * u_bf_x + cy_pop_F[oppose[dir]] * u_bf_y) + 4.5*SQ(cx_pop_F[oppose[dir]] * velocity_x[X][Y] + cy_pop_F[oppose[dir]] * velocity_y[X][Y]) - 1.5*(SQ(velocity_x[X][Y]) + SQ(velocity_y[X][Y])));
							pop_F[dir][X][Y] = (1 - zet)*pop_F[oppose[dir]][Xp][Yp] + zet*temp_pop;
#endif
#if defined Filippova
							double temp_pop, u_bf_x, u_bf_y, zet;
							if (delta[oppose[dir]][X][Y] >= 0.5)
							{
								u_bf_x = (delta[oppose[dir]][X][Y] - 1) * velocity_x[X][Y] / delta[oppose[dir]][X][Y];
								u_bf_y = (delta[oppose[dir]][X][Y] - 1) * velocity_y[X][Y] / delta[oppose[dir]][X][Y];
								zet = (2 * delta[oppose[dir]][X][Y] - 1.0) / (tau);
							}
							else
							{
								u_bf_x = velocity_x[X][Y];
								u_bf_y = velocity_y[X][Y];
								zet = (2 * delta[oppose[dir]][X][Y] - 1.0) / (tau - 1);
							}
							temp_pop = weight_F[oppose[dir]] * density[X][Y] * (1 + 3 * (cx_pop_F[oppose[dir]] * u_bf_x + cy_pop_F[oppose[dir]] * u_bf_y) + 4.5*SQ(cx_pop_F[oppose[dir]] * velocity_x[X][Y] + cy_pop_F[oppose[dir]] * velocity_y[X][Y]) - 1.5*(SQ(velocity_x[X][Y]) + SQ(velocity_y[X][Y])));
							pop_F[dir][X][Y] = (1 - zet)*pop_F[oppose[dir]][Xp][Yp] + zet*temp_pop;
#endif

							//	pop_F[dir][X][Y] = pop_F[oppose[dir]][(X - cx_pop_F[dir] + Nx) % Nx][(Y - cy_pop_F[dir] + Ny) % Ny];
						}
					} //End of: for (int dir = 1; dir < stencil_directions_F; dir++)
				}
			}
		}
#endif

#if defined No_curved_boundary_treatment
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				if (is_solid[X][Y] == FALSE){
					for (int dir = 1; dir < stencil_directions_F; dir++)	{
						Xp = (X - cx_pop_F[dir] + Nx) % Nx;
						Yp = (Y - cy_pop_F[dir] + Ny) % Ny;
						if (is_solid[Xp][Yp] == TRUE){
							pop_F[dir][X][Y] = pop_F[oppose[dir]][Xp][Yp];
						}
					}
				}
			}
		}
#endif
#endif

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: BC_OLD_VERSION
//-----------------------------------------------------------------------------------------------------------------------------
void Flow_solver::BC_OLD_VERSION(int time) {
    /*

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

#if defined X_NOT_Periodic
	start_X = 1;
	end_X = Nx - 2;
#endif

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

#if defined Y_NOT_Periodic
	start_Y = 1;
	end_Y = Ny - 2;
#endif

	///////////////////////////
	/// Bounce-back on walls
	// Due to the presence of the rigid walls at y = 0 and y = Ny - 1, the populations have to be bounced back.
	// Ladd's momentum correction term is included for moving walls (wall velocity parallel to x-axis).
	// Periodicity of the lattice in x-direction is taken into account via the %-operator if spesified.

	int Xtest;
	Xtest = 0; //This parameter is just for avoiding applying X-direc boundary condition multiple times. Therefore it just improves computational time and nothing else
	double rhow, u_X,u_Y;
	int cond, bou;
	double b_vel_x, b_vel_y, b_pre;//velocity and pressure boundary conditions

#if defined X_NOT_Periodic
	for (int Y = start_Y; Y <= end_Y; ++Y) {

		/////////////// Left side
		if (is_solid[1][Y] == FALSE){

			cond = 1;
			if (edge_boundary[1][Y] == TRUE)
			{
				cond = 2;
			}

			for (size_t i = 1; i <= cond; i++)
			{
				bou = boundary_cond[1][Y];
				b_vel_x = boundary_vel_x[1][Y];
			//	b_vel_x = (4.0 * b_vel_x*((Y - 0.5) / (Ny - 2.))*(1 - (Y - 0.5) / (Ny - 2.)));
				b_vel_y = boundary_vel_y[1][Y];
				b_pre = boundary_pres[1][Y];

				if (i == 2)	{
					bou = boundary_cond_2[1][Y];
					b_vel_x = boundary_vel_x_2[1][Y];
			//		b_vel_x=(4.0 * b_vel_x*((Y-0.5) / (Ny - 2.))*(1 - (Y-0.5) / (Ny - 2.)));
					b_vel_y = boundary_vel_y_2[1][Y];
					b_pre = boundary_pres_2[1][Y];
				}

				switch (bou)//Left side
				{
				case 1:		{//Velocity boundary condition along X
#if defined Ladd_Velocity
								pop_F[1][1][Y] = pop_F[3][0][Y] + 6 * cx_pop_F[1] * weight_F[3] * density[1][Y] * b_vel_x + 6 * cy_pop_F[1] * weight_F[3] * density[1][Y] * b_vel_y;
								pop_F[5][1][Y] = pop_F[7][0][(Y - 1 + Ny) % Ny] + 6 * cx_pop_F[5] * weight_F[7] * density[1][Y] * b_vel_x + 6 * cy_pop_F[5] * weight_F[7] * density[1][Y] * b_vel_y;
								pop_F[8][1][Y] = pop_F[6][0][(Y + 1) % Ny] + 6 * cx_pop_F[8] * weight_F[6] * density[1][Y] * b_vel_x + 6 * cy_pop_F[8] * weight_F[6] * density[1][Y] * b_vel_y;

#endif
#if defined Zou_He_Velocity
#if defined Y_NOT_Periodic

								if (Y==1) 	{ //Left Bottom edge

									rhow = (density[2][1] + density[1][2]) / 2.0;

									pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*rhow * b_vel_x;
									pop_F[5][1][Y] = pop_F[7][1][Y] + 1. / 6.*rhow * b_vel_x + 1. / 6.*rhow * b_vel_y;
									pop_F[8][1][Y] = rhow / 2. - pop_F[3][1][Y] - pop_F[4][1][Y] - pop_F[7][1][Y] - 0.5*pop_F[0][1][Y] - 1. / 2.*rhow * b_vel_y - 1. / 3.*rhow * b_vel_x;
								}
								else if (Y == Ny - 2){//Left Top edge

									rhow = (density[2][Ny-2] + density[1][Ny-3]) / 2.0;

									pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*rhow * b_vel_x;
									pop_F[5][1][Y] = 1. / 2.*rhow - pop_F[2][1][Y] - pop_F[3][1][Y] - pop_F[6][1][Y] - 1. / 2.*pop_F[0][1][Y] - 1. / 3.*rhow * b_vel_x + 1. / 2.*rhow * b_vel_y;
									pop_F[8][1][Y] = pop_F[6][1][Y] + 1. / 6.*rhow * b_vel_x - 1. / 6.*rhow * b_vel_y;
								}
								else
								{
								rhow = (pop_F[0][1][Y] + pop_F[2][1][Y] + pop_F[4][1][Y] + 2 * (pop_F[3][1][Y] + pop_F[6][1][Y] + pop_F[7][1][Y])) / (1. - b_vel_x);
								pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*rhow * b_vel_x;
								pop_F[5][1][Y] = pop_F[7][1][Y] - 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*rhow * b_vel_x + 0.5*rhow*b_vel_y;
								pop_F[8][1][Y] = pop_F[6][1][Y] + 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*rhow * b_vel_x - 0.5*rhow*b_vel_y;
								}
#endif
#if defined Y_Periodic
								rhow = (pop_F[0][1][Y] + pop_F[2][1][Y] + pop_F[4][1][Y] + 2 * (pop_F[3][1][Y] + pop_F[6][1][Y] + pop_F[7][1][Y])) / (1. - b_vel_x);
								pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*rhow * b_vel_x;
								pop_F[5][1][Y] = pop_F[7][1][Y] - 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*rhow * b_vel_x + 0.5*rhow*b_vel_y;
								pop_F[8][1][Y] = pop_F[6][1][Y] + 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*rhow * b_vel_x - 0.5*rhow*b_vel_y;
#endif


#endif
								break;
				}

				case 2:{break; }//Velocity boundary condition along Y

				case 3: {//Pressure boundary condition along X
#if defined Ladd_Pressure
							pop_F[1][1][Y] = pop_F[3][0][Y] + 6 * cx_pop_F[1] * weight_F[3] * b_pre * velocity_x[1][Y];
							pop_F[5][1][Y] = pop_F[7][0][(Y - 1 + Ny) % Ny] + 6 * cx_pop_F[5] * weight_F[7] * b_pre * velocity_x[1][Y];
							pop_F[8][1][Y] = pop_F[6][0][(Y + 1) % Ny] + 6 * cx_pop_F[8] * weight_F[6] * b_pre * velocity_x[1][Y];
#endif
#if defined Zou_He_Pressure
#if defined Y_NOT_Periodic

							//Important: Y-velocity is assumed zero
							if (Y == 1) 	{ //Left bottom edge
								u_X = (velocity_x[2][1] + velocity_x[1][2]) / 2.0;

								pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*b_pre * u_X;
								pop_F[5][1][Y] = pop_F[7][1][Y] + 1. / 6.*b_pre * u_X;
								pop_F[8][1][Y] = b_pre / 2. - pop_F[3][1][Y] - pop_F[4][1][Y] - pop_F[7][1][Y] - 0.5*pop_F[0][1][Y] - 1. / 3.*b_pre * u_X;
							}
							else if (Y == Ny - 2){
								u_X = (velocity_x[2][Ny - 2] + velocity_x[1][Ny - 3]) / 2.0;

								pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*b_pre * u_X;
								pop_F[5][1][Y] = 1. / 2.*b_pre - pop_F[2][1][Y] - pop_F[3][1][Y] - pop_F[6][1][Y] - 1. / 2.*pop_F[0][1][Y] - 1. / 3.*b_pre * u_X;
								pop_F[8][1][Y] = pop_F[6][1][Y] + 1. / 6.*b_pre * u_X;
							}
							else{
								u_X = 1 - 1.0 / (b_pre)*(pop_F[0][1][Y] + pop_F[2][1][Y] + pop_F[4][1][Y] + 2 * (pop_F[3][1][Y] + pop_F[6][1][Y] + pop_F[7][1][Y]));
								pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*b_pre * (4.0 * u_X*(Y / (Ny - 2.))*(1 - Y / (Ny - 2.)));
								pop_F[5][1][Y] = pop_F[7][1][Y] - 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*b_pre * u_X;
								pop_F[8][1][Y] = pop_F[6][1][Y] + 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*b_pre * u_X;
							}
#endif
#if defined Y_Periodic
							u_X = 1 - 1.0 / (b_pre)*(pop_F[0][1][Y] + pop_F[2][1][Y] + pop_F[4][1][Y] + 2 * (pop_F[3][1][Y] + pop_F[6][1][Y] + pop_F[7][1][Y]));
							pop_F[1][1][Y] = pop_F[3][1][Y] + 2. / 3.*b_pre * (4.0 * u_X*(Y / (Ny - 2.))*(1 - Y / (Ny - 2.)));
							pop_F[5][1][Y] = pop_F[7][1][Y] - 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*b_pre * u_X;
							pop_F[8][1][Y] = pop_F[6][1][Y] + 0.5* (pop_F[2][1][Y] - pop_F[4][1][Y]) + 1. / 6.*b_pre * u_X;
#endif
#endif

							break;
				}
				case 4:{break; }//Pressure boundary condition along Y

				case 5: {//Zero gradient boundary condition in X direction
							//First order
							//pop_F[1][1][Y] = pop_F[1][2][Y];
							//pop_F[5][1][Y] = pop_F[5][2][Y];
							//pop_F[8][1][Y] = pop_F[8][2][Y];

							////Second order
							pop_F[1][1][Y] = 2 * pop_F[1][2][Y] - pop_F[1][3][Y];
							pop_F[5][1][Y] = 2 * pop_F[5][2][Y] - pop_F[5][3][Y];
							pop_F[8][1][Y] = 2 * pop_F[8][2][Y] - pop_F[8][3][Y];

							break;
				}

				case 6:{break; }

				case 7: {//Stationary wall boundary condition in X direction
#if defined Ladd_Stationary_Wall
							pop_F[1][1][Y] = pop_F[3][0][Y];
							pop_F[5][1][Y] = pop_F[7][0][(Y - 1 + Ny) % Ny];
							pop_F[8][1][Y] = pop_F[6][0][(Y + 1) % Ny];
#endif

#if defined Zou_He_Stationary_Wall
#if defined Y_NOT_Periodic

							if (Y == 1) 	{ //Left bottom edge
								rhow = (density[2][1] + density[1][2]) / 2.0;

								pop_F[1][1][Y] = pop_F[3][1][Y];
								pop_F[5][1][Y] = pop_F[7][1][Y];
								pop_F[8][1][Y] = rhow / 2. - pop_F[3][1][Y] - pop_F[4][1][Y] - pop_F[7][1][Y] - 0.5*pop_F[0][1][Y];
							}
							else if (Y == Ny - 2){//Left top edge
								rhow = (density[2][Ny - 2] + density[1][Ny - 3]) / 2.0;

								pop_F[1][1][Y] = pop_F[3][1][Y];
								pop_F[5][1][Y] = 1. / 2.*rhow - pop_F[2][1][Y] - pop_F[3][1][Y] - pop_F[6][1][Y] - 1. / 2.*pop_F[0][1][Y];
								pop_F[8][1][Y] = pop_F[6][1][Y];
							}
							else
							{

								pop_F[1][1][Y] = pop_F[3][0][Y];
								pop_F[5][1][Y] = pop_F[7][0][(Y - 1 + Ny) % Ny];
								pop_F[8][1][Y] = pop_F[6][0][(Y + 1) % Ny];
							}
#endif
#if defined Y_Periodic
							pop_F[1][1][Y] = pop_F[3][0][Y];
							pop_F[5][1][Y] = pop_F[7][0][(Y - 1 + Ny) % Ny];
							pop_F[8][1][Y] = pop_F[6][0][(Y + 1) % Ny];
#endif
#endif
							break;
				}

				case 8:{break; }

				default: cout << "uuu";
				}
			}
		}

		//////////////////////////////////////////////////////
		////Right wall
		if (is_solid[Nx - 2][Y] == FALSE){

			cond = 1;
			if (edge_boundary[Nx - 2][Y] == TRUE) {
				cond = 2;
			}

			for (size_t i = 1; i <= cond; i++)
			{
				bou = boundary_cond[Nx - 2][Y];
				b_vel_x = boundary_vel_x[Nx - 2][Y];
				b_vel_y = boundary_vel_y[Nx - 2][Y];
				b_pre = boundary_pres[Nx - 2][Y];

				if (i == 2)	{
					bou = boundary_cond_2[Nx - 2][Y];
					b_vel_x = boundary_vel_x_2[Nx - 2][Y];
					b_vel_y = boundary_vel_y_2[Nx - 2][Y];
					b_pre = boundary_pres_2[Nx - 2][Y];
				}

				///////////////////////////////
				switch (bou)//Right side
				{

				case 1:		{//Velocity boundary condition
#if defined Ladd_Velocity
								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 1][Y] + 6 * cx_pop_F[3] * weight_F[1] * density[Nx - 2][Y] * b_vel_x + 6 * cy_pop_F[3] * weight_F[1] * density[Nx - 2][Y] * b_vel_y;
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 1][(Y + 1) % Ny] + 6 * cx_pop_F[7] * weight_F[5] * density[Nx - 2][Y] * b_vel_x + 6 * cy_pop_F[7] * weight_F[5] * density[Nx - 2][Y] * b_vel_y;
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 1][(Y - 1 + Ny) % Ny] + 6 * cx_pop_F[6] * weight_F[8] * density[Nx - 2][Y] * b_vel_x + 6 * cy_pop_F[6] * weight_F[8] * density[Nx - 2][Y] * b_vel_y;
#endif
#if defined Zou_He_Velocity
#if defined Y_NOT_Periodic
								if (Y == 1) 	{ //Right Bottom edge
									rhow = (density[Nx - 3][1] + density[Nx - 2][2]) / 2.0;

									pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
									pop_F[7][Nx - 2][Y] = 0.5*rhow - pop_F[1][Nx - 2][Y] - pop_F[4][Nx - 2][Y] - pop_F[8][Nx - 2][Y] - 1. / 2.*pop_F[0][Nx - 2][Y] + 1. / 3.*rhow*b_vel_x - 1. / 2.*rhow*b_vel_y;
									pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y] - 1. / 6.*rhow*b_vel_x + 1. / 6.*rhow*b_vel_y;
								}
								else if(Y==Ny-2){//Right Top edge
									rhow = (density[Nx - 3][Ny-2] + density[Nx - 2][Ny-3]) / 2.0;

									pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
									pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y] - 1. / 6.*rhow*b_vel_x - 1. / 6.*rhow*b_vel_y;
									pop_F[6][Nx - 2][Y] = rhow / 2. - pop_F[1][Nx - 2][Y] - pop_F[2][Nx - 2][Y] - pop_F[5][Nx - 2][Y] - 1. / 2.*pop_F[0][Nx - 2][Y] + 1. / 2.*rhow*b_vel_y + 1. / 6.*rhow*b_vel_x;
								}
								else
								{
									rhow = (pop_F[0][Nx - 2][Y] + pop_F[2][Nx - 2][Y] + pop_F[4][Nx - 2][Y] + 2 * (pop_F[1][Nx - 2][Y] + pop_F[5][Nx - 2][Y] + pop_F[8][Nx - 2][Y])) / (1. + b_vel_x);
									pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
									pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y] + 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x - 0.5*rhow*b_vel_y;
									pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y] - 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x + 0.5*rhow*b_vel_y;
								}
#endif
#if defined Y_Periodic
								rhow = (pop_F[0][Nx - 2][Y] + pop_F[2][Nx - 2][Y] + pop_F[4][Nx - 2][Y] + 2 * (pop_F[1][Nx - 2][Y] + pop_F[5][Nx - 2][Y] + pop_F[8][Nx - 2][Y])) / (1. + b_vel_x);
								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y] + 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x - 0.5*rhow*b_vel_y;
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y] - 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x + 0.5*rhow*b_vel_y;
#endif
#endif
								break;
				}

				case 2:{break; }

				case 3: {//Pressure boundary condition
#if defined Ladd_Pressure
							pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 1][Y] + 6 * cx_pop_F[3] * weight_F[1] * b_pre * velocity_x[Nx - 2][Y];
							pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 1][(Y + 1) % Ny] + 6 * cx_pop_F[7] * weight_F[5] * b_pre * velocity_x[Nx - 2][Y];
							pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 1][(Y - 1 + Ny) % Ny] + 6 * cx_pop_F[6] * weight_F[8] * b_pre * velocity_x[Nx - 2][Y];
#endif
#if defined Zou_He_Pressure
#if defined Y_NOT_Periodic
							//Important: Y-velocity is assumed zero
							if (Y == 1) 	{ //Right Bottom edge
								u_X = (velocity_x[Nx - 3][1] + velocity_x[Nx - 2][2]) / 2.0;

								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
								pop_F[7][Nx - 2][Y] = 0.5*b_pre - pop_F[1][Nx - 2][Y] - pop_F[4][Nx - 2][Y] - pop_F[8][Nx - 2][Y] - 1. / 2.*pop_F[0][Nx - 2][Y] + 1. / 3.*b_pre*u_X;
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y] - 1. / 6.*b_pre*u_X;
							}
							else if (Y == Ny - 2){
								u_X = (velocity_x[Nx - 3][Ny - 2] + velocity_x[Nx - 2][Ny - 3]) / 2.0;

								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y] - 1. / 6.*b_pre*u_X;
								pop_F[6][Nx - 2][Y] = b_pre / 2. - pop_F[1][Nx - 2][Y] - pop_F[2][Nx - 2][Y] - pop_F[5][Nx - 2][Y] - 1. / 2.*pop_F[0][Nx - 2][Y] + 1. / 6.*b_pre*u_X;
							}
							else
							{
								u_X = (pop_F[0][Nx - 2][Y] + pop_F[2][Nx - 2][Y] + pop_F[4][Nx - 2][Y] + 2 * (pop_F[1][Nx - 2][Y] + pop_F[5][Nx - 2][Y] + pop_F[8][Nx - 2][Y])) / b_pre - 1.0;
								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y] + 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y] - 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
							}
#endif
#if defined Y_Periodic
							u_X = (pop_F[0][Nx - 2][Y] + pop_F[2][Nx - 2][Y] + pop_F[4][Nx - 2][Y] + 2 * (pop_F[1][Nx - 2][Y] + pop_F[5][Nx - 2][Y] + pop_F[8][Nx - 2][Y])) / b_pre - 1.0;
							pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
							pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y] + 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
							pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y] - 0.5* (pop_F[2][Nx - 2][Y] - pop_F[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
#endif
#endif
							break;
				}

				case 4:{break; }

				case 5: {//Zero gradient boundary condition

							//First order
							//pop_F[3][Nx - 2][Y] = pop_F[3][Nx - 3][Y];//zero gradient condition
							//pop_F[7][Nx - 2][Y] = pop_F[7][Nx - 3][Y];
							//pop_F[6][Nx - 2][Y] = pop_F[6][Nx - 3][Y];

							//Second order
							pop_F[3][Nx - 2][Y] = 2 * pop_F[3][Nx - 3][Y] - pop_F[3][Nx - 4][Y];//zero gradient condition
							pop_F[7][Nx - 2][Y] = 2 * pop_F[7][Nx - 3][Y] - pop_F[7][Nx - 4][Y];
							pop_F[6][Nx - 2][Y] = 2 * pop_F[6][Nx - 3][Y] - pop_F[6][Nx - 4][Y];
							break;
				}

				case 6:{break; }

				case 7: {//Stationary wall boundary condition
#if defined Ladd_Stationary_Wall

								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 1][Y];
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 1][(Y + 1) % Ny];
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 1][(Y - 1 + Ny) % Ny];
#endif
#if defined Zou_He_Stationary_Wall
#if defined Y_NOT_Periodic

							if (Y == 1) 	{ //Right Bottom edge
								rhow = (density[Nx - 3][1] + density[Nx - 2][2]) / 2.0;

								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y] ;
								pop_F[7][Nx - 2][Y] = 1. / 2. * rhow - pop_F[1][Nx - 2][Y] - pop_F[4][Nx - 2][Y] - pop_F[8][Nx - 2][Y] - 1. / 2. * pop_F[0][Nx - 2][Y];
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 2][Y];
							}
							else if (Y == Ny - 2){
								rhow = (density[Nx - 3][Ny - 2] + density[Nx - 2][Ny - 3]) / 2.0;

								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 2][Y];
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 2][Y];
								pop_F[6][Nx - 2][Y] = rhow / 2. - pop_F[1][Nx - 2][Y] - pop_F[2][Nx - 2][Y] - pop_F[5][Nx - 2][Y] - 1. / 2.*pop_F[0][Nx - 2][Y];
							}
							else
							{
								pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 1][Y];
								pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 1][(Y + 1) % Ny];
								pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 1][(Y - 1 + Ny) % Ny];
							}
#endif
#if defined Y_Periodic
							pop_F[3][Nx - 2][Y] = pop_F[1][Nx - 1][Y];
							pop_F[7][Nx - 2][Y] = pop_F[5][Nx - 1][(Y + 1) % Ny];
							pop_F[6][Nx - 2][Y] = pop_F[8][Nx - 1][(Y - 1 + Ny) % Ny];
#endif
#endif
							break;
				}

				case 8:{break; }

				default: cout << "fff";
				}
			}
		}
	}//End of: for (int X = start_X; X <= end_X; ++X)
#endif

////	/////////////////////////////////////////////
////	/////////////////////////////////////////////

	#if defined Y_NOT_Periodic
		for (int X = start_X; X <= end_X; ++X) {

			////Bottom boundary
			if (is_solid[X][1] == FALSE){

				cond = 1;
				if (edge_boundary[X][1] == TRUE) {
					cond = 2;
				}

				for (size_t i = 1; i <= cond; i++)
				{
					bou = boundary_cond[X][1];
					b_vel_x = boundary_vel_x[X][1];
					b_vel_y = boundary_vel_y[X][1];
					b_pre = boundary_pres[X][1];

					if (i == 2)	{
						bou = boundary_cond_2[X][1];
						b_vel_x = boundary_vel_x_2[X][1];
						b_vel_y = boundary_vel_y_2[X][1];
						b_pre = boundary_pres_2[X][1];
					}

					switch (bou)
					{

					case 1:{break; }

					case 2:
					{
#if defined Ladd_Velocity
							  pop_F[2][X][1] = pop_F[4][X][0] + 6 * cx_pop_F[2] * weight_F[4] * density[X][1] * b_vel_x + 6 * cy_pop_F[2] * weight_F[4] * density[X][1] * b_vel_y;
							  pop_F[5][X][1] = pop_F[7][(X - 1 + Nx) % Nx][0] + 6 * cx_pop_F[5] * weight_F[7] * density[X][1] * b_vel_x + 6 * cy_pop_F[5] * weight_F[7] * density[X][1] * b_vel_y;
							  pop_F[6][X][1] = pop_F[8][(X + 1) % Nx][0] + 6 * cx_pop_F[6] * weight_F[8] * density[X][1] * b_vel_x + 6 * cy_pop_F[6] * weight_F[8] * density[X][1] * b_vel_y;
#endif
#if defined Zou_He_Velocity
#if defined X_NOT_Periodic

							  if (X==1)	  {//Left Bottom edge
								  rhow = (density[2][1] + density[1][2]) / 2.0;

								  pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * rhow * b_vel_y;
								  pop_F[5][X][1] = pop_F[7][X][1] + 1. / 6.*rhow * b_vel_x + 1. / 6.*rhow * b_vel_y;
								  pop_F[6][X][1] = rhow / 2. - pop_F[3][X][1] - pop_F[4][X][1] - pop_F[7][X][1] - 1. / 2.*pop_F[0][X][1] - 1. / 2.*rhow * b_vel_x + 1. / 3.*rhow * b_vel_y;
							  }
							  else if (X==Nx-2)
							  {
								  rhow = (density[Nx - 2][2] + density[Nx - 3][1]) / 2.0;

								  pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * rhow * b_vel_y;
								  pop_F[5][X][1] = rhow / 2. - pop_F[1][X][1] - pop_F[4][X][1] - pop_F[8][X][1] - 1. / 2.*pop_F[0][X][1] + 1. / 2.*rhow * b_vel_x - 1. / 3.*rhow * b_vel_y;
								  pop_F[6][X][1] = pop_F[8][X][1] - 1. / 6.*rhow * b_vel_x + 1. / 6.*rhow * b_vel_y;
							  }
							  else
							  {
								  rhow = 1. / (1 - b_vel_y)*(pop_F[0][X][1] + pop_F[1][X][1] + pop_F[3][X][1] + 2 * (pop_F[4][X][1] + pop_F[7][X][1] + pop_F[8][X][1]));
								  pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * rhow * b_vel_y;
								  pop_F[5][X][1] = pop_F[7][X][1] - 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*rhow * b_vel_y + 1. / 2.*rhow * b_vel_x;
								  pop_F[6][X][1] = pop_F[8][X][1] + 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*rhow * b_vel_y - 1. / 2.*rhow * b_vel_x;
							  }
#endif
#if defined X_Periodic
							  rhow = 1. / (1 - b_vel_y)*(pop_F[0][X][1] + pop_F[1][X][1] + pop_F[3][X][1] + 2 * (pop_F[4][X][1] + pop_F[7][X][1] + pop_F[8][X][1]));
							  pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * rhow * b_vel_y;
							  pop_F[5][X][1] = pop_F[7][X][1] - 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*rhow * b_vel_y + 1. / 2.*rhow * b_vel_x;
							  pop_F[6][X][1] = pop_F[8][X][1] + 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*rhow * b_vel_y - 1. / 2.*rhow * b_vel_x;
#endif

#endif
							  break;
					}

					case 3:{break; }

					case 4: {
#if defined Ladd_Pressure
								pop_F[2][X][1] = pop_F[4][X][0] + 6 * cy_pop_F[2] * weight_F[4] * b_pre * velocity_y[X][1];
								pop_F[5][X][1] = pop_F[7][(X - 1 + Nx) % Nx][0] + 6 * cy_pop_F[5] * weight_F[7] * b_pre * velocity_y[X][1];
								pop_F[6][X][1] = pop_F[8][(X + 1) % Nx][0] + 6 * cy_pop_F[6] * weight_F[8] * b_pre * velocity_y[X][1];
#endif
#if defined Zou_He_Pressure
#if defined X_NOT_Periodic


								//Important: X-velocity is assumed zero
								if (X == 1)	{//Left bottom edge
									u_Y = (velocity_y[2][1] + velocity_y[1][2]) / 2.0;

									pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * b_pre * u_Y;
									pop_F[5][X][1] = pop_F[7][X][1] + 1. / 6.*b_pre * u_Y;
									pop_F[6][X][1] = b_pre / 2. - pop_F[3][X][1] - pop_F[4][X][1] - pop_F[7][X][1] - 1 / 2.*pop_F[0][X][1] + 1. / 3.*b_pre * u_Y;
								}
								else if (X == Nx - 2){//Right bottom edge
									u_Y = (velocity_y[Nx - 2][2] + velocity_y[Nx - 3][1]) / 2.0;

									pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * b_pre * u_Y;
									pop_F[5][X][1] = b_pre / 2. - pop_F[1][X][1] - pop_F[4][X][1] - pop_F[8][X][1] - 1. / 2.*pop_F[0][X][1] - 1. / 3.*b_pre * u_Y;
									pop_F[6][X][1] = pop_F[8][X][1] + 1. / 6.*b_pre * u_Y;
								}
								else
								{
								u_Y = 1.0 - 1. / b_pre*(pop_F[0][X][1] + pop_F[1][X][1] + pop_F[3][X][1] + 2 * (pop_F[4][X][1] + pop_F[7][X][1] + pop_F[8][X][1]));
								pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * b_pre * u_Y;
								pop_F[5][X][1] = pop_F[7][X][1] - 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*b_pre * u_Y;
								pop_F[6][X][1] = pop_F[8][X][1] + 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*b_pre * u_Y;
								}
#endif
#if defined X_Periodic
								u_Y = 1.0 - 1. / b_pre*(pop_F[0][X][1] + pop_F[1][X][1] + pop_F[3][X][1] + 2 * (pop_F[4][X][1] + pop_F[7][X][1] + pop_F[8][X][1]));
								pop_F[2][X][1] = pop_F[4][X][1] + 2. / 3.  * b_pre * u_Y;
								pop_F[5][X][1] = pop_F[7][X][1] - 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*b_pre * u_Y;
								pop_F[6][X][1] = pop_F[8][X][1] + 0.5*(pop_F[1][X][1] - pop_F[3][X][1]) + 1. / 6.*b_pre * u_Y;
#endif

#endif
								break;
					}
					case 5:{break; }

					case 6: {
								//First order
								//pop_F[2][X][1] = pop_F[2][X][2];
								//pop_F[5][X][1] = pop_F[5][X][2];
								//pop_F[6][X][1] = pop_F[6][X][2];

								//Second order
								pop_F[2][X][1] = 2 * pop_F[2][X][2] - pop_F[2][X][3];
								pop_F[5][X][1] = 2 * pop_F[5][X][2] - pop_F[5][X][3];
								pop_F[6][X][1] = 2 * pop_F[6][X][2] - pop_F[6][X][3];
								break;
					}
					case 7:{break; }

					case 8: {
#if defined Ladd_Stationary_Wall
								pop_F[2][X][1] = pop_F[4][X][0];
								pop_F[5][X][1] = pop_F[7][(X - 1 + Nx) % Nx][0];
								pop_F[6][X][1] = pop_F[8][(X + 1) % Nx][0];
#endif
#if defined Zou_He_Stationary_Wall
#if defined X_NOT_Periodic
								if (X == 1)	{ //Left Bottom edge
									rhow = (density[2][1]+density[1][2])/2.0;

									pop_F[2][X][1] = pop_F[4][X][1] ;
									pop_F[5][X][1] = pop_F[7][X][1];
									pop_F[6][X][1] = rhow / 2. - pop_F[3][X][1] - pop_F[4][X][1] - pop_F[7][X][1] - 1 / 2.*pop_F[0][X][1];
								}
								else if (X == Nx - 2)	{//Right bottom edge
									rhow = (density[Nx - 2][2] + density[Nx - 3][1]) / 2.0;

									pop_F[2][X][1] = pop_F[4][X][1];
									pop_F[5][X][1] = rhow / 2. - pop_F[1][X][1] - pop_F[4][X][1] - pop_F[8][X][1] - 1. / 2.*pop_F[0][X][1];
									pop_F[6][X][1] = pop_F[8][X][1];
								}
								else
								{
									//rhow = 1. / (1)*(pop_F[0][X][1] + pop_F[1][X][1] + pop_F[3][X][1] + 2 * (pop_F[4][X][1] + pop_F[7][X][1] + pop_F[8][X][1]));
									//pop_F[2][X][1] = pop_F[4][X][1] ;
									//pop_F[5][X][1] = pop_F[7][X][1] - 0.5*(pop_F[1][X][1] - pop_F[3][X][1]);
									//pop_F[6][X][1] = pop_F[8][X][1] + 0.5*(pop_F[1][X][1] - pop_F[3][X][1]);

									pop_F[2][X][1] = pop_F[4][X][0];
									pop_F[5][X][1] = pop_F[7][(X - 1 + Nx) % Nx][0];
									pop_F[6][X][1] = pop_F[8][(X + 1) % Nx][0];
								}
#endif
#if defined X_Periodic
								pop_F[2][X][1] = pop_F[4][X][0];
								pop_F[5][X][1] = pop_F[7][(X - 1 + Nx) % Nx][0];
								pop_F[6][X][1] = pop_F[8][(X + 1) % Nx][0];
#endif
#endif
								break;
					}
					default: cout << "nnn";
					}
				}
			}
			////////////////////////////////
			////////////////////////////////
			/// TOP surface
			////////////////////////////////

			if (is_solid[X][Ny - 2] == FALSE){
				cond = 1;
				if (edge_boundary[X][Ny - 2] == TRUE) {
					cond = 2;
				}

				for (size_t i = 1; i <= cond; i++)
				{
					bou = boundary_cond[X][Ny - 2];
					b_vel_x = boundary_vel_x[X][Ny - 2];
					b_vel_y = boundary_vel_y[X][Ny - 2];
					b_pre = boundary_pres[X][Ny - 2];
					if (i == 2)	{
						bou = boundary_cond_2[X][Ny - 2];
						b_vel_x = boundary_vel_x_2[X][Ny - 2];
						b_vel_y = boundary_vel_y_2[X][Ny - 2];
						b_pre = boundary_pres_2[X][Ny - 2];
					}

					switch (bou)
					{

					case 1:{break; }

					case 2:
					{
#if defined Ladd_Velocity
							  pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 1] + 6 * cx_pop_F[4] * weight_F[2] * density[X][Ny - 2] * b_vel_x + 6 * cy_pop_F[4] * weight_F[2] * density[X][Ny - 2] * b_vel_y;
							  pop_F[7][X][Ny - 2] = pop_F[5][(X + 1) % Nx][Ny - 1] + 6 * cx_pop_F[7] * weight_F[5] * density[X][Ny - 2] * b_vel_x + 6 * cy_pop_F[7] * weight_F[5] * density[X][Ny - 2] * b_vel_y;
							  pop_F[8][X][Ny - 2] = pop_F[6][(X - 1 + Nx) % Nx][Ny - 1] + 6 * cx_pop_F[8] * weight_F[6] * density[X][Ny - 2] * b_vel_x + 6 * cy_pop_F[8] * weight_F[6] * density[X][Ny - 2] * b_vel_y;
#endif
#if defined Zou_He_Velocity
#if defined X_NOT_Periodic
							  if (X==1)//Left Top corner
							  {
								  rhow = (density[2][Ny-2] + density[1][Ny-3]) / 2.0;

								  pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
								  pop_F[7][X][Ny - 2] = 1. / 2.*rhow - pop_F[2][X][Ny - 2] - pop_F[3][X][Ny - 2] - pop_F[6][X][Ny - 2] - 1. / 2.*pop_F[0][X][Ny - 2] - 1. / 2.*rhow*b_vel_x + 1. / 3.*rhow*b_vel_y;
								  pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2] + 1. / 6.*rhow*b_vel_x - 1. / 6.*rhow*b_vel_y;
							  }
							  else if (X==Nx-2)//Right Top corner
							  {
								  rhow = (density[Nx - 2][Ny - 3] + density[Nx - 3][Ny - 2]) / 2.0;

								  pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
								  pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2] - 1. / 6.*rhow*b_vel_x - 1. / 6.*rhow*b_vel_y;
								  pop_F[8][X][Ny - 2] = rhow / 2. - pop_F[1][X][Ny - 2] - pop_F[2][X][Ny - 2] - pop_F[5][X][Ny - 2] - 1. / 2.*pop_F[0][X][Ny - 2] + 1. / 2.*rhow*b_vel_x + 1. / 3.*rhow*b_vel_y;
							  }
							  else
							  {
							  rhow = 1. / (1.0 + b_vel_y)*(pop_F[0][X][Ny - 2] + pop_F[1][X][Ny - 2] + pop_F[3][X][Ny - 2] + 2 * (pop_F[2][X][Ny - 2] + pop_F[5][X][Ny - 2] + pop_F[6][X][Ny - 2]));
							  pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
							  pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2] + 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y - 1. / 2.*rhow*b_vel_x;
							  pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2] - 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y + 1. / 2.*rhow*b_vel_x;
							  }
#endif
#if defined X_Periodic
							  rhow = 1. / (1.0 + b_vel_y)*(pop_F[0][X][Ny - 2] + pop_F[1][X][Ny - 2] + pop_F[3][X][Ny - 2] + 2 * (pop_F[2][X][Ny - 2] + pop_F[5][X][Ny - 2] + pop_F[6][X][Ny - 2]));
							  pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
							  pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2] + 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y - 1. / 2.*rhow*b_vel_x;
							  pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2] - 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y + 1. / 2.*rhow*b_vel_x;
#endif

#endif
							  break;
					}

					case 3:{break; }

					case 4: {//Pressure boundary condition along Y
#if defined Ladd_Pressure

								pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 1] + 6 * cy_pop_F[4] * weight_F[2] * b_pre * velocity_y[X][Ny - 2];
								pop_F[7][X][Ny - 2] = pop_F[5][(X + 1) % Nx][Ny - 1] + 6 * cy_pop_F[7] * weight_F[5] * b_pre * velocity_y[X][Ny - 2];
								pop_F[8][X][Ny - 2] = pop_F[6][(X - 1 + Nx) % Nx][Ny - 1] + 6 * cy_pop_F[8] * weight_F[6] * b_pre * velocity_y[X][Ny - 2];
#endif
#if defined Zou_He_Pressure
#if defined X_NOT_Periodic

								//Important: X-velocity is assumed zero
								if (X == 1)
								{
									u_Y = (velocity_y[2][Ny - 2] + velocity_y[1][Ny - 3]) / 2.0;

									pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
									pop_F[7][X][Ny - 2] = 1. / 2.*b_pre - pop_F[2][X][Ny - 2] - pop_F[3][X][Ny - 2] - pop_F[6][X][Ny - 2] - 1. / 2.*pop_F[0][X][Ny - 2] + 1. / 3.*b_pre* u_Y;
									pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2] - 1. / 6.*b_pre*u_Y;
								}
								else if (X == Nx - 2)//Right Top corner
								{
									u_Y = (velocity_y[Nx - 2][Ny - 3] + velocity_y[Nx - 3][Ny - 2]) / 2.0;

									pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
									pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2] - 1. / 6.*b_pre*u_Y;
									pop_F[8][X][Ny - 2] = b_pre / 2. - pop_F[1][X][Ny - 2] - pop_F[2][X][Ny - 2] - pop_F[5][X][Ny - 2] - 1. / 2.*pop_F[0][X][Ny - 2] + 1. / 3.*b_pre*u_Y;
								}
								else
								{
								u_Y = (pop_F[0][X][Ny - 2] + pop_F[1][X][Ny - 2] + pop_F[3][X][Ny - 2] + 2 * (pop_F[2][X][Ny - 2] + pop_F[5][X][Ny - 2] + pop_F[6][X][Ny - 2])) / b_pre - 1.0;
								pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
								pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2] + 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
								pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2] - 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
								}
#endif
#if defined X_Periodic
								u_Y = (pop_F[0][X][Ny - 2] + pop_F[1][X][Ny - 2] + pop_F[3][X][Ny - 2] + 2 * (pop_F[2][X][Ny - 2] + pop_F[5][X][Ny - 2] + pop_F[6][X][Ny - 2])) / b_pre - 1.0;
								pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
								pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2] + 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
								pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2] - 0.5*(pop_F[1][X][Ny - 2] - pop_F[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
#endif

#endif
								break;
					}

					case 5:{break; }

					case 6: {
								//First order
								//pop_F[4][X][Ny - 2] = pop_F[4][X][Ny - 3];
								//pop_F[7][X][Ny - 2] = pop_F[7][X][Ny - 3];
								//pop_F[8][X][Ny - 2] = pop_F[8][X][Ny - 3];

								//Second order
								pop_F[4][X][Ny - 2] = 2 * pop_F[4][X][Ny - 3] - pop_F[4][X][Ny - 4];
								pop_F[7][X][Ny - 2] = 2 * pop_F[7][X][Ny - 3] - pop_F[7][X][Ny - 4];
								pop_F[8][X][Ny - 2] = 2 * pop_F[8][X][Ny - 3] - pop_F[8][X][Ny - 4];
								break;
					}

					case 7:{break; }

					case 8: {
#if defined Ladd_Stationary_Wall
								pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 1];
								pop_F[7][X][Ny - 2] = pop_F[5][(X + 1) % Nx][Ny - 1];
								pop_F[8][X][Ny - 2] = pop_F[6][(X - 1 + Nx) % Nx][Ny - 1];
#endif
#if defined Zou_He_Stationary_Wall
#if defined X_NOT_Periodic

								if (X == 1)	{//Left top corner
									rhow = (density[2][Ny - 2] + density[1][Ny - 3]) / 2.0;

									pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2];
									pop_F[7][X][Ny - 2] = 1. / 2.*rhow - pop_F[2][X][Ny - 2] - pop_F[3][X][Ny - 2] - pop_F[6][X][Ny - 2] - 1. / 2.*pop_F[0][X][Ny - 2];
									pop_F[8][X][Ny - 2] = pop_F[6][X][Ny - 2];
								}
								else if (X == Nx - 2)//Right Top corner
								{
									rhow = (density[Nx - 2][Ny - 3] + density[Nx - 3][Ny - 2]) / 2.0;

									pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 2] ;
									pop_F[7][X][Ny - 2] = pop_F[5][X][Ny - 2];
									pop_F[8][X][Ny - 2] = rhow / 2. - pop_F[1][X][Ny - 2] - pop_F[2][X][Ny - 2] - pop_F[5][X][Ny - 2] - 1. / 2.*pop_F[0][X][Ny - 2];
								}
								else
								{
									pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 1];
									pop_F[7][X][Ny - 2] = pop_F[5][(X + 1) % Nx][Ny - 1];
									pop_F[8][X][Ny - 2] = pop_F[6][(X - 1 + Nx) % Nx][Ny - 1];
								}
#endif
#if defined X_Periodic
								pop_F[4][X][Ny - 2] = pop_F[2][X][Ny - 1];
								pop_F[7][X][Ny - 2] = pop_F[5][(X + 1) % Nx][Ny - 1];
								pop_F[8][X][Ny - 2] = pop_F[6][(X - 1 + Nx) % Nx][Ny - 1];
#endif


#endif
								break;
					}
					default: cout << "mmm";
					}
				}
			}
		}//End of: for (int X = start_X; X <= end_X; ++X)
	#endif*/

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

void Flow_solver::momenta() {

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE)
            {
                previous_velocity_magnitude[X][Y] = velocity_magnitude[X][Y];
                density[X][Y] = 0;
                velocity_x[X][Y] = 0;
                velocity_y[X][Y] = 0;

                for (int c_i = 0; c_i < stencil_directions_F; ++c_i) {
                    density[X][Y]    += pop_F[c_i][X][Y];
                    velocity_x[X][Y] += cx_pop_F[c_i] * pop_F[c_i][X][Y];
                    velocity_y[X][Y] += cy_pop_F[c_i] * pop_F[c_i][X][Y];
                }
                #if defined Ladd_Pressure
                if (boundary_cond[X][Y] == 3 || boundary_cond[X][Y] == 4)
                {
                    density[X][Y] = boundary_pres[X][Y];
                }
                #endif
                velocity_x[X][Y] /= density[X][Y];
                velocity_y[X][Y] /= density[X][Y];
                velocity_magnitude[X][Y] = sqrt(SQ(velocity_x[X][Y]) + SQ(velocity_x[X][Y]));
                force_x[X][Y] = 0.0;//It is done because Eulerian forces must be zero at the beginning of each iteration
                force_y[X][Y] = 0.0;//It is done because Eulerian forces must be zero at the beginning of each iteration
			}else if (is_solid[X][Y] == TRUE) {
				density[X][Y] = 1.0;
                velocity_x[X][Y] = 0;
                velocity_y[X][Y] = 0;
                velocity_magnitude[X][Y] = 0;
			}
		}
	}

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MEQ
//-----------------------------------------------------------------------------------------------------------------------------
//calculation the equilibrium moment for MRT

double Flow_solver::meq(double RHO, double U, double V, int k)
{
	double x;
	switch (k)
	{
	case 0: {x = RHO; break; }
	case 1: {x = RHO*(-2 + 3 * (U*U + V*V)); break; }
	case 2: {x = RHO*(1 - 3 * (U*U + V*V)); break; }
	case 3: {x = RHO*U; break; }
	case 4: {x = -RHO*U; break; }
	case 5: {x = RHO*V; break; }
	case 6: {x = -RHO*V; break; }
	case 7: {x = RHO*(U*U - V*V); break; }
	case 8: {x = RHO*U*V; break; }
	default: x = 0;
	}
	return x;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: CHECK_RESIDUAL
//-----------------------------------------------------------------------------------------------------------------------------
#if defined Stop_by_Residual
bool Flow_solver::check_residual(){//Based on paper Applied Mathematical Modelling 39 (2015) 2436�2451
	bool st_cr = false;
	double vel_criteria = 0;
	double vel_criteria_2 = 0;

	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			if (is_solid[X][Y] == FALSE)			{
					vel_criteria += abs(velocity_magnitude[X][Y] - previous_velocity_magnitude[X][Y]);
					vel_criteria_2 +=  velocity_magnitude[X][Y];
			}
		}
	}

	if ((vel_criteria / vel_criteria_2)<residual_flow)	{
		st_cr = true;
		return st_cr;
	}
	return st_cr;

}
#endif
//=============================================================================================================================
//=============================================================================================================================
Flow_solver::~Flow_solver()
{
    // 3D arrays
    delete3D(pop_F, stencil_directions_F, Nx);
    delete3D(pop_old_F, stencil_directions_F, Nx);
    delete3D(pop_T, stencil_directions_T, Nx);
    delete3D(pop_old_T, stencil_directions_T, Nx);

    // 2D double arrays
    delete2D(density, Nx);
    delete2D(velocity_x, Nx);
    delete2D(velocity_y, Nx);
    delete2D(velocity_magnitude, Nx);
    delete2D(previous_velocity_magnitude, Nx);

    delete2D(force_x, Nx);
    delete2D(force_y, Nx);
    delete2D(mdf_force_x, Nx);
    delete2D(mdf_force_y, Nx);
    delete2D(mdf_force_thermal, Nx);
    delete2D(mdf_vel_x, Nx);
    delete2D(mdf_vel_y, Nx);
    delete2D(mdf_temperature, Nx);
    delete2D(temperature, Nx);
    delete2D(previous_temperature, Nx);
    delete2D(force_thermal, Nx);

    delete2D(boundary_vel_x, Nx);
    delete2D(boundary_vel_y, Nx);
    delete2D(boundary_pres, Nx);
    delete2D(boundary_temp, Nx);
    delete2D(boundary_temp_2, Nx);
    delete2D(boundary_vel_x_2, Nx);
    delete2D(boundary_vel_y_2, Nx);
    delete2D(boundary_pres_2, Nx);

    // 2D int arrays
    delete2D_int(is_solid, Nx);
    delete2D_int(solid_thermal, Nx);
    delete2D_int(boundary_cond, Nx);
    delete2D_int(boundary_cond_temp, Nx);
    delete2D_int(boundary_cond_temp_2, Nx);
    delete2D_int(edge_boundary, Nx);
    delete2D_int(boundary_cond_2, Nx);
}
