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
			velocity_x[X][Y] = 0.0;
			velocity_y[X][Y] = 0.0;
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

	pop = new double**[9];
	pop_old = new double**[9];
	delta = new double**[9];

	for (int c_i = 0; c_i < 9; ++c_i) {
		pop[c_i] = new double*[Nx];
		pop_old[c_i] = new double*[Nx];
		delta[c_i] = new double*[Nx];
		for (int X = 0; X < Nx; ++X) {
			pop[c_i][X] = new double[Ny];
			pop_old[c_i][X] = new double[Ny];
			delta[c_i][X] = new double[Ny];
			for (int Y = 0; Y < Ny; ++Y) {
				pop[c_i][X][Y] = 0;
				pop_old[c_i][X][Y] = 0;
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
			for (int c_i = 0; c_i < 9; ++c_i) {
				pop_old[c_i][X][Y] = pop_eq[c_i];
				pop[c_i][X][Y] = pop_eq[c_i];
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
	pop_eq[0] = weight[0] * den * (1 - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[1] = weight[1] * den * (1 + 3 * (vel_x)+4.5 * SQ(vel_x) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[2] = weight[2] * den * (1 + 3 * (vel_y)+4.5 * SQ(vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[3] = weight[3] * den * (1 + 3 * (-vel_x) + 4.5 * SQ(-vel_x) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[4] = weight[4] * den * (1 + 3 * (-vel_y) + 4.5 * SQ(-vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[5] = weight[5] * den * (1 + 3 * (vel_x + vel_y) + 4.5 * SQ(vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[6] = weight[6] * den * (1 + 3 * (-vel_x + vel_y) + 4.5 * SQ(-vel_x + vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[7] = weight[7] * den * (1 + 3 * (-vel_x - vel_y) + 4.5 * SQ(-vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	pop_eq[8] = weight[8] * den * (1 + 3 * (vel_x - vel_y) + 4.5 * SQ(vel_x - vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM
//-----------------------------------------------------------------------------------------------------------------------------
//
void Flow_solver::LBM(int time) {

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
	double ***swap_temp = pop_old;
	pop_old = pop;
	pop = swap_temp;
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


	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE){

				/// Compute new populations
				// This is the lattice Boltzmann equation (combined collision and propagation) including external forcing.
				// Periodicity of the lattice in x-direction is taken into account by the %-operator.
#if defined Shan_Chen
                double v0 = velocity_x[X][Y] + (force_x[X][Y] + gravity) * tau;
                double v1 = velocity_y[X][Y] + (force_y[X][Y] + gbeta * (temperature[X][Y] - Tref)) * tau;
                double den = density[X][Y];
                equilibrium(den, v0, v1);
                for (int i = 0; i < 9; i++)
                {
                    int newX = (X + c_x[i] + Nx) % Nx;
                    int newY = (Y + c_y[i] + Ny) % Ny;
                    pop[i][newX][newY] = pop_old[i][X][Y] * (1.0 - omega) + pop_eq[i] * omega;
                }

				/*pop[0][X][Y] = pop_old[0][X][Y] * (1 - omega) + pop_eq[0] * omega;
				pop[1][(X + 1) % Nx][Y] = pop_old[1][X][Y] * (1 - omega) + pop_eq[1] * omega;
				pop[2][X][(Y + 1) % Ny] = pop_old[2][X][Y] * (1 - omega) + pop_eq[2] * omega;
				pop[3][(X - 1 + Nx) % Nx][Y] = pop_old[3][X][Y] * (1 - omega) + pop_eq[3] * omega;
				pop[4][X][(Y - 1 + Ny) % Ny] = pop_old[4][X][Y] * (1 - omega) + pop_eq[4] * omega;
				pop[5][(X + 1) % Nx][(Y + 1) % Ny] = pop_old[5][X][Y] * (1 - omega) + pop_eq[5] * omega;
				pop[6][(X - 1 + Nx) % Nx][(Y + 1) % Ny] = pop_old[6][X][Y] * (1 - omega) + pop_eq[6] * omega;
				pop[7][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny] = pop_old[7][X][Y] * (1 - omega) + pop_eq[7] * omega;
				pop[8][(X + 1) % Nx][(Y - 1 + Ny) % Ny] = pop_old[8][X][Y] * (1 - omega) + pop_eq[8] * omega;*/
#endif
#if defined Luo
                equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
                for (int i = 0; i < 9; i++)
                {
                    int newX = (X + c_x[i] + Nx) % Nx;
                    int newY = (Y + c_y[i] + Ny) % Ny;
                    pop[i][newX][newY] = pop_old[i][X][Y] * (1.0 - omega) + pop_eq[i] * omega
                                       + 3 * weight[i] * (c_x[i] * (force_x[X][Y] + gravity)  + c_y[i] * (force_y[X][Y] + gbeta * (temperature[X][Y] - Tref)));// * density[X][Y];
                }
				/*pop[0][X][Y] = pop_old[0][X][Y] * (1 - omega) + pop_eq[0] * omega;
				pop[1][(X + 1) % Nx][Y] = pop_old[1][X][Y] * (1 - omega) + pop_eq[1] * omega + 3 * weight[1] * (1 * (force_x[X][Y] + gravity) + 0 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[2][X][(Y + 1) % Ny] = pop_old[2][X][Y] * (1 - omega) + pop_eq[2] * omega + 3 * weight[2] * (0 * (force_x[X][Y] + gravity) + 1 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[3][(X - 1 + Nx) % Nx][Y] = pop_old[3][X][Y] * (1 - omega) + pop_eq[3] * omega + 3 * weight[3] * (-1 * (force_x[X][Y] + gravity) + 0 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[4][X][(Y - 1 + Ny) % Ny] = pop_old[4][X][Y] * (1 - omega) + pop_eq[4] * omega + 3 * weight[4] * (0 * (force_x[X][Y] + gravity) - 1 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[5][(X + 1) % Nx][(Y + 1) % Ny] = pop_old[5][X][Y] * (1 - omega) + pop_eq[5] * omega + 3 * weight[5] * (1 * (force_x[X][Y] + gravity) + 1 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[6][(X - 1 + Nx) % Nx][(Y + 1) % Ny] = pop_old[6][X][Y] * (1 - omega) + pop_eq[6] * omega + 3 * weight[6] * (-1 * (force_x[X][Y] + gravity) + 1 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[7][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny] = pop_old[7][X][Y] * (1 - omega) + pop_eq[7] * omega + 3 * weight[7] * (-1 * (force_x[X][Y] + gravity) - 1 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];
				pop[8][(X + 1) % Nx][(Y - 1 + Ny) % Ny] = pop_old[8][X][Y] * (1 - omega) + pop_eq[8] * omega + 3 * weight[8] * (1 * (force_x[X][Y] + gravity) - 1 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)))* density[X][Y];*/
#endif
#if defined Guo
#if defined SRT
                equilibrium(density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);

                for (int i = 0; i < 9; i++)
                {
                    int newX = (X + c_x[i] + Nx) % Nx;
                    int newY = (Y + c_y[i] + Ny) % Ny;
                    pop[i][newX][newY] = pop_old[i][X][Y] * (1.0 - omega) + pop_eq[i] * omega
                                       + (1.0 - 0.5 * omega) * weight[i] * ((force_x[X][Y] + gravity) * (3 * (c_x[i] - velocity_x[X][Y]) + 9.0 * c_x[i] * (c_x[i] * velocity_x[X][Y] + c_y[i] * velocity_y[X][Y]))
                                               + (force_y[X][Y] + gbeta * (temperature[X][Y] - Tref)) * (3 * (c_y[i] - velocity_y[X][Y]) + 9.0 * c_y[i] * (c_x[i] * velocity_x[X][Y] + c_y[i] * velocity_y[X][Y])));// * density[X][Y];
                }

				/*pop[0][X][Y] = pop_old[0][X][Y] * (1 - omega) + pop_eq[0] * omega + (1 - 0.5*omega)*weight[0] * ((force_x[X][Y] + gravity)*(3 * (0 - velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (0 - velocity_y[X][Y])))*density[X][Y];
				pop[1][(X + 1) % Nx][Y] = pop_old[1][X][Y] * (1 - omega) + pop_eq[1] * omega + (1 - 0.5*omega)*weight[1] * ((force_x[X][Y] + gravity)*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (0 - velocity_y[X][Y])))*density[X][Y];
				pop[2][X][(Y + 1) % Ny] = pop_old[2][X][Y] * (1 - omega) + pop_eq[2] * omega + (1 - 0.5*omega)*weight[2] * ((force_x[X][Y] + gravity)*(3 * (0 - velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (1 * velocity_y[X][Y])))*density[X][Y];
				pop[3][(X - 1 + Nx) % Nx][Y] = pop_old[3][X][Y] * (1 - omega) + pop_eq[3] * omega + (1 - 0.5*omega)*weight[3] * ((force_x[X][Y] + gravity)*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (0 - velocity_y[X][Y])))*density[X][Y];
				pop[4][X][(Y - 1 + Ny) % Ny] = pop_old[4][X][Y] * (1 - omega) + pop_eq[4] * omega + (1 - 0.5*omega)*weight[4] * ((force_x[X][Y] + gravity)*(3 * (0 - velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (-1 * velocity_y[X][Y])))*density[X][Y];
				pop[5][(X + 1) % Nx][(Y + 1) % Ny] = pop_old[5][X][Y] * (1 - omega) + pop_eq[5] * omega + (1 - 0.5*omega)*weight[5] * ((force_x[X][Y] + gravity)*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))))*density[X][Y];
				pop[6][(X - 1 + Nx) % Nx][(Y + 1) % Ny] = pop_old[6][X][Y] * (1 - omega) + pop_eq[6] * omega + (1 - 0.5*omega)*weight[6] * ((force_x[X][Y] + gravity)*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (-1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))))*density[X][Y];
				pop[7][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny] = pop_old[7][X][Y] * (1 - omega) + pop_eq[7] * omega + (1 - 0.5*omega)*weight[7] * ((force_x[X][Y] + gravity)*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))))*density[X][Y];
				pop[8][(X + 1) % Nx][(Y - 1 + Ny) % Ny] = pop_old[8][X][Y] * (1 - omega) + pop_eq[8] * omega + (1 - 0.5*omega)*weight[8] * ((force_x[X][Y] + gravity)*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))))*density[X][Y];*/
#endif
#if defined MRT
				// Transformation from velocity space to moment space:
				m[0] = pop_old[0][X][Y] + pop_old[1][X][Y] + pop_old[2][X][Y] + pop_old[3][X][Y] + pop_old[4][X][Y] + pop_old[5][X][Y] + pop_old[6][X][Y] + pop_old[7][X][Y] + pop_old[8][X][Y];
				m[1] = -4 * pop_old[0][X][Y] - pop_old[1][X][Y] - pop_old[2][X][Y] - pop_old[3][X][Y] - pop_old[4][X][Y] + 2 * (pop_old[5][X][Y] + pop_old[6][X][Y] + pop_old[7][X][Y] + pop_old[8][X][Y]);
				m[2] = 4 * pop_old[0][X][Y] - 2 * (pop_old[1][X][Y] + pop_old[2][X][Y] + pop_old[3][X][Y] + pop_old[4][X][Y]) + pop_old[5][X][Y] + pop_old[6][X][Y] + pop_old[7][X][Y] + pop_old[8][X][Y];
				m[3] = pop_old[1][X][Y] - pop_old[3][X][Y] + pop_old[5][X][Y] - pop_old[6][X][Y] - pop_old[7][X][Y] + pop_old[8][X][Y];
				m[4] = -2 * (pop_old[1][X][Y] - pop_old[3][X][Y]) + pop_old[5][X][Y] - pop_old[6][X][Y] - pop_old[7][X][Y] + pop_old[8][X][Y];
				m[5] = pop_old[2][X][Y] - pop_old[4][X][Y] + pop_old[5][X][Y] + pop_old[6][X][Y] - pop_old[7][X][Y] - pop_old[8][X][Y];
				m[6] = -2 * (pop_old[2][X][Y] - pop_old[4][X][Y]) + pop_old[5][X][Y] + pop_old[6][X][Y] - pop_old[7][X][Y] - pop_old[8][X][Y];
				m[7] = pop_old[1][X][Y] - pop_old[2][X][Y] + pop_old[3][X][Y] - pop_old[4][X][Y];
				m[8] = pop_old[5][X][Y] - pop_old[6][X][Y] + pop_old[7][X][Y] - pop_old[8][X][Y];

				fforce[0] = weight[0] * ((force_x[X][Y] + gravity)*(3 * (0 - velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (0 - velocity_y[X][Y])));//*density[X][Y];
				fforce[1] = weight[1] * ((force_x[X][Y] + gravity)*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (0 - velocity_y[X][Y])));//*density[X][Y];
				fforce[2] = weight[2] * ((force_x[X][Y] + gravity)*(3 * (0 - velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (1 * velocity_y[X][Y])));//*density[X][Y];
				fforce[3] = weight[3] * ((force_x[X][Y] + gravity)*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (0 - velocity_y[X][Y])));//*density[X][Y];
				fforce[4] = weight[4] * ((force_x[X][Y] + gravity)*(3 * (0 - velocity_x[X][Y])) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (-1 * velocity_y[X][Y])));//*density[X][Y];
				fforce[5] = weight[5] * ((force_x[X][Y] + gravity)*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))));//*density[X][Y];
				fforce[6] = weight[6] * ((force_x[X][Y] + gravity)*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (1 - velocity_y[X][Y]) + 9 * 1 * (-1 * velocity_x[X][Y] + (1 * velocity_y[X][Y]))));//*density[X][Y];
				fforce[7] = weight[7] * ((force_x[X][Y] + gravity)*(3 * (-1 - velocity_x[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (-1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))));//*density[X][Y];
				fforce[8] = weight[8] * ((force_x[X][Y] + gravity)*(3 * (1 - velocity_x[X][Y]) + 9 * 1 * (1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))) + (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) * (3 * (-1 - velocity_y[X][Y]) + 9 * -1 * (1 * velocity_x[X][Y] + (-1 * velocity_y[X][Y]))));//*density[X][Y];

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

				pop[0][X][Y] = 4 * m[0] - 4 * (m[1] - m[2]);
				pop[1][(X + 1) % Nx][Y] = 4 * m[0] - m[1] - 2 * m[2] + 6 * m[3] - 6 * m[4] + 9 * m[7];
				pop[2][X][(Y + 1) % Ny] = 4 * m[0] - m[1] - 2 * m[2] + 6 * m[5] - 6 * m[6] - 9 * m[7];
				pop[3][(X - 1 + Nx) % Nx][Y] = 4 * m[0] - m[1] - 2 * m[2] - 6 * m[3] + 6 * m[4] + 9 * m[7];
				pop[4][X][(Y - 1 + Ny) % Ny] = 4 * m[0] - m[1] - 2 * m[2] - 6 * m[5] + 6 * m[6] - 9 * m[7];
				pop[5][(X + 1) % Nx][(Y + 1) % Ny] = 4 * m[0] + 2 * m[1] + m[2] + 6 * m[3] + 3 * m[4] + 6 * m[5] + 3 * m[6] + 9 * m[8];
				pop[6][(X - 1 + Nx) % Nx][(Y + 1) % Ny] = 4 * m[0] + 2 * m[1] + m[2] - 6 * m[3] - 3 * m[4] + 6 * m[5] + 3 * m[6] - 9 * m[8];
				pop[7][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny] = 4 * m[0] + 2 * m[1] + m[2] - 6 * m[3] - 3 * m[4] - 6 * m[5] - 3 * m[6] + 9 * m[8];
				pop[8][(X + 1) % Nx][(Y - 1 + Ny) % Ny] = 4 * m[0] + 2 * m[1] + m[2] + 6 * m[3] + 3 * m[4] - 6 * m[5] - 3 * m[6] - 9 * m[8];
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
								pop[1][1][Y] = pop[3][0][Y] + 6 * c_x[1] * weight[3] * density[1][Y] * b_vel_x + 6 * c_y[1] * weight[3] * density[1][Y] * b_vel_y;
								pop[5][1][Y] = pop[7][0][(Y - 1 + Ny) % Ny] + 6 * c_x[5] * weight[7] * density[1][Y] * b_vel_x + 6 * c_y[5] * weight[7] * density[1][Y] * b_vel_y;
								pop[8][1][Y] = pop[6][0][(Y + 1) % Ny] + 6 * c_x[8] * weight[6] * density[1][Y] * b_vel_x + 6 * c_y[8] * weight[6] * density[1][Y] * b_vel_y;

#endif
#if defined Zou_He_Velocity
#if defined Y_NOT_Periodic

								if (Y==1) 	{ //Left Bottom edge

									rhow = (density[2][1] + density[1][2]) / 2.0;

									pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*rhow * b_vel_x;
									pop[5][1][Y] = pop[7][1][Y] + 1. / 6.*rhow * b_vel_x + 1. / 6.*rhow * b_vel_y;
									pop[8][1][Y] = rhow / 2. - pop[3][1][Y] - pop[4][1][Y] - pop[7][1][Y] - 0.5*pop[0][1][Y] - 1. / 2.*rhow * b_vel_y - 1. / 3.*rhow * b_vel_x;
								}
								else if (Y == Ny - 2){//Left Top edge

									rhow = (density[2][Ny-2] + density[1][Ny-3]) / 2.0;

									pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*rhow * b_vel_x;
									pop[5][1][Y] = 1. / 2.*rhow - pop[2][1][Y] - pop[3][1][Y] - pop[6][1][Y] - 1. / 2.*pop[0][1][Y] - 1. / 3.*rhow * b_vel_x + 1. / 2.*rhow * b_vel_y;
									pop[8][1][Y] = pop[6][1][Y] + 1. / 6.*rhow * b_vel_x - 1. / 6.*rhow * b_vel_y;
								}
								else
								{
								rhow = (pop[0][1][Y] + pop[2][1][Y] + pop[4][1][Y] + 2 * (pop[3][1][Y] + pop[6][1][Y] + pop[7][1][Y])) / (1. - b_vel_x);
								pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*rhow * b_vel_x;
								pop[5][1][Y] = pop[7][1][Y] - 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*rhow * b_vel_x + 0.5*rhow*b_vel_y;
								pop[8][1][Y] = pop[6][1][Y] + 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*rhow * b_vel_x - 0.5*rhow*b_vel_y;
								}
#endif
#if defined Y_Periodic
								rhow = (pop[0][1][Y] + pop[2][1][Y] + pop[4][1][Y] + 2 * (pop[3][1][Y] + pop[6][1][Y] + pop[7][1][Y])) / (1. - b_vel_x);
								pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*rhow * b_vel_x;
								pop[5][1][Y] = pop[7][1][Y] - 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*rhow * b_vel_x + 0.5*rhow*b_vel_y;
								pop[8][1][Y] = pop[6][1][Y] + 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*rhow * b_vel_x - 0.5*rhow*b_vel_y;
#endif


#endif
								break;
				}

				case 2:{break; }//Velocity boundary condition along Y

				case 3: {//Pressure boundary condition along X
#if defined Ladd_Pressure
							pop[1][1][Y] = pop[3][0][Y] + 6 * c_x[1] * weight[3] * b_pre * velocity_x[1][Y];
							pop[5][1][Y] = pop[7][0][(Y - 1 + Ny) % Ny] + 6 * c_x[5] * weight[7] * b_pre * velocity_x[1][Y];
							pop[8][1][Y] = pop[6][0][(Y + 1) % Ny] + 6 * c_x[8] * weight[6] * b_pre * velocity_x[1][Y];
#endif
#if defined Zou_He_Pressure
#if defined Y_NOT_Periodic

							//Important: Y-velocity is assumed zero
							if (Y == 1) 	{ //Left bottom edge
								u_X = (velocity_x[2][1] + velocity_x[1][2]) / 2.0;

								pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*b_pre * u_X;
								pop[5][1][Y] = pop[7][1][Y] + 1. / 6.*b_pre * u_X;
								pop[8][1][Y] = b_pre / 2. - pop[3][1][Y] - pop[4][1][Y] - pop[7][1][Y] - 0.5*pop[0][1][Y] - 1. / 3.*b_pre * u_X;
							}
							else if (Y == Ny - 2){
								u_X = (velocity_x[2][Ny - 2] + velocity_x[1][Ny - 3]) / 2.0;

								pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*b_pre * u_X;
								pop[5][1][Y] = 1. / 2.*b_pre - pop[2][1][Y] - pop[3][1][Y] - pop[6][1][Y] - 1. / 2.*pop[0][1][Y] - 1. / 3.*b_pre * u_X;
								pop[8][1][Y] = pop[6][1][Y] + 1. / 6.*b_pre * u_X;
							}
							else{
								u_X = 1 - 1.0 / (b_pre)*(pop[0][1][Y] + pop[2][1][Y] + pop[4][1][Y] + 2 * (pop[3][1][Y] + pop[6][1][Y] + pop[7][1][Y]));
								pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*b_pre * (4.0 * u_X*(Y / (Ny - 2.))*(1 - Y / (Ny - 2.)));
								pop[5][1][Y] = pop[7][1][Y] - 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*b_pre * u_X;
								pop[8][1][Y] = pop[6][1][Y] + 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*b_pre * u_X;
							}
#endif
#if defined Y_Periodic
							u_X = 1 - 1.0 / (b_pre)*(pop[0][1][Y] + pop[2][1][Y] + pop[4][1][Y] + 2 * (pop[3][1][Y] + pop[6][1][Y] + pop[7][1][Y]));
							pop[1][1][Y] = pop[3][1][Y] + 2. / 3.*b_pre * (4.0 * u_X*(Y / (Ny - 2.))*(1 - Y / (Ny - 2.)));
							pop[5][1][Y] = pop[7][1][Y] - 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*b_pre * u_X;
							pop[8][1][Y] = pop[6][1][Y] + 0.5* (pop[2][1][Y] - pop[4][1][Y]) + 1. / 6.*b_pre * u_X;
#endif
#endif

							break;
				}
				case 4:{break; }//Pressure boundary condition along Y

				case 5: {//Zero gradient boundary condition in X direction
							//First order
							//pop[1][1][Y] = pop[1][2][Y];
							//pop[5][1][Y] = pop[5][2][Y];
							//pop[8][1][Y] = pop[8][2][Y];

							////Second order
							pop[1][1][Y] = 2 * pop[1][2][Y] - pop[1][3][Y];
							pop[5][1][Y] = 2 * pop[5][2][Y] - pop[5][3][Y];
							pop[8][1][Y] = 2 * pop[8][2][Y] - pop[8][3][Y];

							break;
				}

				case 6:{break; }

				case 7: {//Stationary wall boundary condition in X direction
#if defined Ladd_Stationary_Wall
							pop[1][1][Y] = pop[3][0][Y];
							pop[5][1][Y] = pop[7][0][(Y - 1 + Ny) % Ny];
							pop[8][1][Y] = pop[6][0][(Y + 1) % Ny];
#endif

#if defined Zou_He_Stationary_Wall
#if defined Y_NOT_Periodic

							if (Y == 1) 	{ //Left bottom edge
								rhow = (density[2][1] + density[1][2]) / 2.0;

								pop[1][1][Y] = pop[3][1][Y];
								pop[5][1][Y] = pop[7][1][Y];
								pop[8][1][Y] = rhow / 2. - pop[3][1][Y] - pop[4][1][Y] - pop[7][1][Y] - 0.5*pop[0][1][Y];
							}
							else if (Y == Ny - 2){//Left top edge
								rhow = (density[2][Ny - 2] + density[1][Ny - 3]) / 2.0;

								pop[1][1][Y] = pop[3][1][Y];
								pop[5][1][Y] = 1. / 2.*rhow - pop[2][1][Y] - pop[3][1][Y] - pop[6][1][Y] - 1. / 2.*pop[0][1][Y];
								pop[8][1][Y] = pop[6][1][Y];
							}
							else
							{

								pop[1][1][Y] = pop[3][0][Y];
								pop[5][1][Y] = pop[7][0][(Y - 1 + Ny) % Ny];
								pop[8][1][Y] = pop[6][0][(Y + 1) % Ny];
							}
#endif
#if defined Y_Periodic
							pop[1][1][Y] = pop[3][0][Y];
							pop[5][1][Y] = pop[7][0][(Y - 1 + Ny) % Ny];
							pop[8][1][Y] = pop[6][0][(Y + 1) % Ny];
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
								pop[3][Nx - 2][Y] = pop[1][Nx - 1][Y] + 6 * c_x[3] * weight[1] * density[Nx - 2][Y] * b_vel_x + 6 * c_y[3] * weight[1] * density[Nx - 2][Y] * b_vel_y;
								pop[7][Nx - 2][Y] = pop[5][Nx - 1][(Y + 1) % Ny] + 6 * c_x[7] * weight[5] * density[Nx - 2][Y] * b_vel_x + 6 * c_y[7] * weight[5] * density[Nx - 2][Y] * b_vel_y;
								pop[6][Nx - 2][Y] = pop[8][Nx - 1][(Y - 1 + Ny) % Ny] + 6 * c_x[6] * weight[8] * density[Nx - 2][Y] * b_vel_x + 6 * c_y[6] * weight[8] * density[Nx - 2][Y] * b_vel_y;
#endif
#if defined Zou_He_Velocity
#if defined Y_NOT_Periodic
								if (Y == 1) 	{ //Right Bottom edge
									rhow = (density[Nx - 3][1] + density[Nx - 2][2]) / 2.0;

									pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
									pop[7][Nx - 2][Y] = 0.5*rhow - pop[1][Nx - 2][Y] - pop[4][Nx - 2][Y] - pop[8][Nx - 2][Y] - 1. / 2.*pop[0][Nx - 2][Y] + 1. / 3.*rhow*b_vel_x - 1. / 2.*rhow*b_vel_y;
									pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y] - 1. / 6.*rhow*b_vel_x + 1. / 6.*rhow*b_vel_y;
								}
								else if(Y==Ny-2){//Right Top edge
									rhow = (density[Nx - 3][Ny-2] + density[Nx - 2][Ny-3]) / 2.0;

									pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
									pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y] - 1. / 6.*rhow*b_vel_x - 1. / 6.*rhow*b_vel_y;
									pop[6][Nx - 2][Y] = rhow / 2. - pop[1][Nx - 2][Y] - pop[2][Nx - 2][Y] - pop[5][Nx - 2][Y] - 1. / 2.*pop[0][Nx - 2][Y] + 1. / 2.*rhow*b_vel_y + 1. / 6.*rhow*b_vel_x;
								}
								else
								{
									rhow = (pop[0][Nx - 2][Y] + pop[2][Nx - 2][Y] + pop[4][Nx - 2][Y] + 2 * (pop[1][Nx - 2][Y] + pop[5][Nx - 2][Y] + pop[8][Nx - 2][Y])) / (1. + b_vel_x);
									pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
									pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y] + 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x - 0.5*rhow*b_vel_y;
									pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y] - 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x + 0.5*rhow*b_vel_y;
								}
#endif
#if defined Y_Periodic
								rhow = (pop[0][Nx - 2][Y] + pop[2][Nx - 2][Y] + pop[4][Nx - 2][Y] + 2 * (pop[1][Nx - 2][Y] + pop[5][Nx - 2][Y] + pop[8][Nx - 2][Y])) / (1. + b_vel_x);
								pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*rhow*b_vel_x;
								pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y] + 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x - 0.5*rhow*b_vel_y;
								pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y] - 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * rhow*b_vel_x + 0.5*rhow*b_vel_y;
#endif
#endif
								break;
				}

				case 2:{break; }

				case 3: {//Pressure boundary condition
#if defined Ladd_Pressure
							pop[3][Nx - 2][Y] = pop[1][Nx - 1][Y] + 6 * c_x[3] * weight[1] * b_pre * velocity_x[Nx - 2][Y];
							pop[7][Nx - 2][Y] = pop[5][Nx - 1][(Y + 1) % Ny] + 6 * c_x[7] * weight[5] * b_pre * velocity_x[Nx - 2][Y];
							pop[6][Nx - 2][Y] = pop[8][Nx - 1][(Y - 1 + Ny) % Ny] + 6 * c_x[6] * weight[8] * b_pre * velocity_x[Nx - 2][Y];
#endif
#if defined Zou_He_Pressure
#if defined Y_NOT_Periodic
							//Important: Y-velocity is assumed zero
							if (Y == 1) 	{ //Right Bottom edge
								u_X = (velocity_x[Nx - 3][1] + velocity_x[Nx - 2][2]) / 2.0;

								pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
								pop[7][Nx - 2][Y] = 0.5*b_pre - pop[1][Nx - 2][Y] - pop[4][Nx - 2][Y] - pop[8][Nx - 2][Y] - 1. / 2.*pop[0][Nx - 2][Y] + 1. / 3.*b_pre*u_X;
								pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y] - 1. / 6.*b_pre*u_X;
							}
							else if (Y == Ny - 2){
								u_X = (velocity_x[Nx - 3][Ny - 2] + velocity_x[Nx - 2][Ny - 3]) / 2.0;

								pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
								pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y] - 1. / 6.*b_pre*u_X;
								pop[6][Nx - 2][Y] = b_pre / 2. - pop[1][Nx - 2][Y] - pop[2][Nx - 2][Y] - pop[5][Nx - 2][Y] - 1. / 2.*pop[0][Nx - 2][Y] + 1. / 6.*b_pre*u_X;
							}
							else
							{
								u_X = (pop[0][Nx - 2][Y] + pop[2][Nx - 2][Y] + pop[4][Nx - 2][Y] + 2 * (pop[1][Nx - 2][Y] + pop[5][Nx - 2][Y] + pop[8][Nx - 2][Y])) / b_pre - 1.0;
								pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
								pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y] + 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
								pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y] - 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
							}
#endif
#if defined Y_Periodic
							u_X = (pop[0][Nx - 2][Y] + pop[2][Nx - 2][Y] + pop[4][Nx - 2][Y] + 2 * (pop[1][Nx - 2][Y] + pop[5][Nx - 2][Y] + pop[8][Nx - 2][Y])) / b_pre - 1.0;
							pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] - 2. / 3.*b_pre*u_X;
							pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y] + 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
							pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y] - 0.5* (pop[2][Nx - 2][Y] - pop[4][Nx - 2][Y]) - 1. / 6. * b_pre*u_X;
#endif
#endif
							break;
				}

				case 4:{break; }

				case 5: {//Zero gradient boundary condition

							//First order
							//pop[3][Nx - 2][Y] = pop[3][Nx - 3][Y];//zero gradient condition
							//pop[7][Nx - 2][Y] = pop[7][Nx - 3][Y];
							//pop[6][Nx - 2][Y] = pop[6][Nx - 3][Y];

							//Second order
							pop[3][Nx - 2][Y] = 2 * pop[3][Nx - 3][Y] - pop[3][Nx - 4][Y];//zero gradient condition
							pop[7][Nx - 2][Y] = 2 * pop[7][Nx - 3][Y] - pop[7][Nx - 4][Y];
							pop[6][Nx - 2][Y] = 2 * pop[6][Nx - 3][Y] - pop[6][Nx - 4][Y];
							break;
				}

				case 6:{break; }

				case 7: {//Stationary wall boundary condition
#if defined Ladd_Stationary_Wall

								pop[3][Nx - 2][Y] = pop[1][Nx - 1][Y];
								pop[7][Nx - 2][Y] = pop[5][Nx - 1][(Y + 1) % Ny];
								pop[6][Nx - 2][Y] = pop[8][Nx - 1][(Y - 1 + Ny) % Ny];
#endif
#if defined Zou_He_Stationary_Wall
#if defined Y_NOT_Periodic

							if (Y == 1) 	{ //Right Bottom edge
								rhow = (density[Nx - 3][1] + density[Nx - 2][2]) / 2.0;

								pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y] ;
								pop[7][Nx - 2][Y] = 1. / 2. * rhow - pop[1][Nx - 2][Y] - pop[4][Nx - 2][Y] - pop[8][Nx - 2][Y] - 1. / 2. * pop[0][Nx - 2][Y];
								pop[6][Nx - 2][Y] = pop[8][Nx - 2][Y];
							}
							else if (Y == Ny - 2){
								rhow = (density[Nx - 3][Ny - 2] + density[Nx - 2][Ny - 3]) / 2.0;

								pop[3][Nx - 2][Y] = pop[1][Nx - 2][Y];
								pop[7][Nx - 2][Y] = pop[5][Nx - 2][Y];
								pop[6][Nx - 2][Y] = rhow / 2. - pop[1][Nx - 2][Y] - pop[2][Nx - 2][Y] - pop[5][Nx - 2][Y] - 1. / 2.*pop[0][Nx - 2][Y];
							}
							else
							{
								pop[3][Nx - 2][Y] = pop[1][Nx - 1][Y];
								pop[7][Nx - 2][Y] = pop[5][Nx - 1][(Y + 1) % Ny];
								pop[6][Nx - 2][Y] = pop[8][Nx - 1][(Y - 1 + Ny) % Ny];
							}
#endif
#if defined Y_Periodic
							pop[3][Nx - 2][Y] = pop[1][Nx - 1][Y];
							pop[7][Nx - 2][Y] = pop[5][Nx - 1][(Y + 1) % Ny];
							pop[6][Nx - 2][Y] = pop[8][Nx - 1][(Y - 1 + Ny) % Ny];
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
							  pop[2][X][1] = pop[4][X][0] + 6 * c_x[2] * weight[4] * density[X][1] * b_vel_x + 6 * c_y[2] * weight[4] * density[X][1] * b_vel_y;
							  pop[5][X][1] = pop[7][(X - 1 + Nx) % Nx][0] + 6 * c_x[5] * weight[7] * density[X][1] * b_vel_x + 6 * c_y[5] * weight[7] * density[X][1] * b_vel_y;
							  pop[6][X][1] = pop[8][(X + 1) % Nx][0] + 6 * c_x[6] * weight[8] * density[X][1] * b_vel_x + 6 * c_y[6] * weight[8] * density[X][1] * b_vel_y;
#endif
#if defined Zou_He_Velocity
#if defined X_NOT_Periodic

							  if (X==1)	  {//Left Bottom edge
								  rhow = (density[2][1] + density[1][2]) / 2.0;

								  pop[2][X][1] = pop[4][X][1] + 2. / 3.  * rhow * b_vel_y;
								  pop[5][X][1] = pop[7][X][1] + 1. / 6.*rhow * b_vel_x + 1. / 6.*rhow * b_vel_y;
								  pop[6][X][1] = rhow / 2. - pop[3][X][1] - pop[4][X][1] - pop[7][X][1] - 1. / 2.*pop[0][X][1] - 1. / 2.*rhow * b_vel_x + 1. / 3.*rhow * b_vel_y;
							  }
							  else if (X==Nx-2)
							  {
								  rhow = (density[Nx - 2][2] + density[Nx - 3][1]) / 2.0;

								  pop[2][X][1] = pop[4][X][1] + 2. / 3.  * rhow * b_vel_y;
								  pop[5][X][1] = rhow / 2. - pop[1][X][1] - pop[4][X][1] - pop[8][X][1] - 1. / 2.*pop[0][X][1] + 1. / 2.*rhow * b_vel_x - 1. / 3.*rhow * b_vel_y;
								  pop[6][X][1] = pop[8][X][1] - 1. / 6.*rhow * b_vel_x + 1. / 6.*rhow * b_vel_y;
							  }
							  else
							  {
								  rhow = 1. / (1 - b_vel_y)*(pop[0][X][1] + pop[1][X][1] + pop[3][X][1] + 2 * (pop[4][X][1] + pop[7][X][1] + pop[8][X][1]));
								  pop[2][X][1] = pop[4][X][1] + 2. / 3.  * rhow * b_vel_y;
								  pop[5][X][1] = pop[7][X][1] - 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*rhow * b_vel_y + 1. / 2.*rhow * b_vel_x;
								  pop[6][X][1] = pop[8][X][1] + 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*rhow * b_vel_y - 1. / 2.*rhow * b_vel_x;
							  }
#endif
#if defined X_Periodic
							  rhow = 1. / (1 - b_vel_y)*(pop[0][X][1] + pop[1][X][1] + pop[3][X][1] + 2 * (pop[4][X][1] + pop[7][X][1] + pop[8][X][1]));
							  pop[2][X][1] = pop[4][X][1] + 2. / 3.  * rhow * b_vel_y;
							  pop[5][X][1] = pop[7][X][1] - 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*rhow * b_vel_y + 1. / 2.*rhow * b_vel_x;
							  pop[6][X][1] = pop[8][X][1] + 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*rhow * b_vel_y - 1. / 2.*rhow * b_vel_x;
#endif

#endif
							  break;
					}

					case 3:{break; }

					case 4: {
#if defined Ladd_Pressure
								pop[2][X][1] = pop[4][X][0] + 6 * c_y[2] * weight[4] * b_pre * velocity_y[X][1];
								pop[5][X][1] = pop[7][(X - 1 + Nx) % Nx][0] + 6 * c_y[5] * weight[7] * b_pre * velocity_y[X][1];
								pop[6][X][1] = pop[8][(X + 1) % Nx][0] + 6 * c_y[6] * weight[8] * b_pre * velocity_y[X][1];
#endif
#if defined Zou_He_Pressure
#if defined X_NOT_Periodic


								//Important: X-velocity is assumed zero
								if (X == 1)	{//Left bottom edge
									u_Y = (velocity_y[2][1] + velocity_y[1][2]) / 2.0;

									pop[2][X][1] = pop[4][X][1] + 2. / 3.  * b_pre * u_Y;
									pop[5][X][1] = pop[7][X][1] + 1. / 6.*b_pre * u_Y;
									pop[6][X][1] = b_pre / 2. - pop[3][X][1] - pop[4][X][1] - pop[7][X][1] - 1 / 2.*pop[0][X][1] + 1. / 3.*b_pre * u_Y;
								}
								else if (X == Nx - 2){//Right bottom edge
									u_Y = (velocity_y[Nx - 2][2] + velocity_y[Nx - 3][1]) / 2.0;

									pop[2][X][1] = pop[4][X][1] + 2. / 3.  * b_pre * u_Y;
									pop[5][X][1] = b_pre / 2. - pop[1][X][1] - pop[4][X][1] - pop[8][X][1] - 1. / 2.*pop[0][X][1] - 1. / 3.*b_pre * u_Y;
									pop[6][X][1] = pop[8][X][1] + 1. / 6.*b_pre * u_Y;
								}
								else
								{
								u_Y = 1.0 - 1. / b_pre*(pop[0][X][1] + pop[1][X][1] + pop[3][X][1] + 2 * (pop[4][X][1] + pop[7][X][1] + pop[8][X][1]));
								pop[2][X][1] = pop[4][X][1] + 2. / 3.  * b_pre * u_Y;
								pop[5][X][1] = pop[7][X][1] - 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*b_pre * u_Y;
								pop[6][X][1] = pop[8][X][1] + 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*b_pre * u_Y;
								}
#endif
#if defined X_Periodic
								u_Y = 1.0 - 1. / b_pre*(pop[0][X][1] + pop[1][X][1] + pop[3][X][1] + 2 * (pop[4][X][1] + pop[7][X][1] + pop[8][X][1]));
								pop[2][X][1] = pop[4][X][1] + 2. / 3.  * b_pre * u_Y;
								pop[5][X][1] = pop[7][X][1] - 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*b_pre * u_Y;
								pop[6][X][1] = pop[8][X][1] + 0.5*(pop[1][X][1] - pop[3][X][1]) + 1. / 6.*b_pre * u_Y;
#endif

#endif
								break;
					}
					case 5:{break; }

					case 6: {
								//First order
								//pop[2][X][1] = pop[2][X][2];
								//pop[5][X][1] = pop[5][X][2];
								//pop[6][X][1] = pop[6][X][2];

								//Second order
								pop[2][X][1] = 2 * pop[2][X][2] - pop[2][X][3];
								pop[5][X][1] = 2 * pop[5][X][2] - pop[5][X][3];
								pop[6][X][1] = 2 * pop[6][X][2] - pop[6][X][3];
								break;
					}
					case 7:{break; }

					case 8: {
#if defined Ladd_Stationary_Wall
								pop[2][X][1] = pop[4][X][0];
								pop[5][X][1] = pop[7][(X - 1 + Nx) % Nx][0];
								pop[6][X][1] = pop[8][(X + 1) % Nx][0];
#endif
#if defined Zou_He_Stationary_Wall
#if defined X_NOT_Periodic
								if (X == 1)	{ //Left Bottom edge
									rhow = (density[2][1]+density[1][2])/2.0;

									pop[2][X][1] = pop[4][X][1] ;
									pop[5][X][1] = pop[7][X][1];
									pop[6][X][1] = rhow / 2. - pop[3][X][1] - pop[4][X][1] - pop[7][X][1] - 1 / 2.*pop[0][X][1];
								}
								else if (X == Nx - 2)	{//Right bottom edge
									rhow = (density[Nx - 2][2] + density[Nx - 3][1]) / 2.0;

									pop[2][X][1] = pop[4][X][1];
									pop[5][X][1] = rhow / 2. - pop[1][X][1] - pop[4][X][1] - pop[8][X][1] - 1. / 2.*pop[0][X][1];
									pop[6][X][1] = pop[8][X][1];
								}
								else
								{
									//rhow = 1. / (1)*(pop[0][X][1] + pop[1][X][1] + pop[3][X][1] + 2 * (pop[4][X][1] + pop[7][X][1] + pop[8][X][1]));
									//pop[2][X][1] = pop[4][X][1] ;
									//pop[5][X][1] = pop[7][X][1] - 0.5*(pop[1][X][1] - pop[3][X][1]);
									//pop[6][X][1] = pop[8][X][1] + 0.5*(pop[1][X][1] - pop[3][X][1]);

									pop[2][X][1] = pop[4][X][0];
									pop[5][X][1] = pop[7][(X - 1 + Nx) % Nx][0];
									pop[6][X][1] = pop[8][(X + 1) % Nx][0];
								}
#endif
#if defined X_Periodic
								pop[2][X][1] = pop[4][X][0];
								pop[5][X][1] = pop[7][(X - 1 + Nx) % Nx][0];
								pop[6][X][1] = pop[8][(X + 1) % Nx][0];
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
							  pop[4][X][Ny - 2] = pop[2][X][Ny - 1] + 6 * c_x[4] * weight[2] * density[X][Ny - 2] * b_vel_x + 6 * c_y[4] * weight[2] * density[X][Ny - 2] * b_vel_y;
							  pop[7][X][Ny - 2] = pop[5][(X + 1) % Nx][Ny - 1] + 6 * c_x[7] * weight[5] * density[X][Ny - 2] * b_vel_x + 6 * c_y[7] * weight[5] * density[X][Ny - 2] * b_vel_y;
							  pop[8][X][Ny - 2] = pop[6][(X - 1 + Nx) % Nx][Ny - 1] + 6 * c_x[8] * weight[6] * density[X][Ny - 2] * b_vel_x + 6 * c_y[8] * weight[6] * density[X][Ny - 2] * b_vel_y;
#endif
#if defined Zou_He_Velocity
#if defined X_NOT_Periodic
							  if (X==1)//Left Top corner
							  {
								  rhow = (density[2][Ny-2] + density[1][Ny-3]) / 2.0;

								  pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
								  pop[7][X][Ny - 2] = 1. / 2.*rhow - pop[2][X][Ny - 2] - pop[3][X][Ny - 2] - pop[6][X][Ny - 2] - 1. / 2.*pop[0][X][Ny - 2] - 1. / 2.*rhow*b_vel_x + 1. / 3.*rhow*b_vel_y;
								  pop[8][X][Ny - 2] = pop[6][X][Ny - 2] + 1. / 6.*rhow*b_vel_x - 1. / 6.*rhow*b_vel_y;
							  }
							  else if (X==Nx-2)//Right Top corner
							  {
								  rhow = (density[Nx - 2][Ny - 3] + density[Nx - 3][Ny - 2]) / 2.0;

								  pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
								  pop[7][X][Ny - 2] = pop[5][X][Ny - 2] - 1. / 6.*rhow*b_vel_x - 1. / 6.*rhow*b_vel_y;
								  pop[8][X][Ny - 2] = rhow / 2. - pop[1][X][Ny - 2] - pop[2][X][Ny - 2] - pop[5][X][Ny - 2] - 1. / 2.*pop[0][X][Ny - 2] + 1. / 2.*rhow*b_vel_x + 1. / 3.*rhow*b_vel_y;
							  }
							  else
							  {
							  rhow = 1. / (1.0 + b_vel_y)*(pop[0][X][Ny - 2] + pop[1][X][Ny - 2] + pop[3][X][Ny - 2] + 2 * (pop[2][X][Ny - 2] + pop[5][X][Ny - 2] + pop[6][X][Ny - 2]));
							  pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
							  pop[7][X][Ny - 2] = pop[5][X][Ny - 2] + 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y - 1. / 2.*rhow*b_vel_x;
							  pop[8][X][Ny - 2] = pop[6][X][Ny - 2] - 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y + 1. / 2.*rhow*b_vel_x;
							  }
#endif
#if defined X_Periodic
							  rhow = 1. / (1.0 + b_vel_y)*(pop[0][X][Ny - 2] + pop[1][X][Ny - 2] + pop[3][X][Ny - 2] + 2 * (pop[2][X][Ny - 2] + pop[5][X][Ny - 2] + pop[6][X][Ny - 2]));
							  pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*rhow*b_vel_y;
							  pop[7][X][Ny - 2] = pop[5][X][Ny - 2] + 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y - 1. / 2.*rhow*b_vel_x;
							  pop[8][X][Ny - 2] = pop[6][X][Ny - 2] - 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*rhow*b_vel_y + 1. / 2.*rhow*b_vel_x;
#endif

#endif
							  break;
					}

					case 3:{break; }

					case 4: {//Pressure boundary condition along Y
#if defined Ladd_Pressure

								pop[4][X][Ny - 2] = pop[2][X][Ny - 1] + 6 * c_y[4] * weight[2] * b_pre * velocity_y[X][Ny - 2];
								pop[7][X][Ny - 2] = pop[5][(X + 1) % Nx][Ny - 1] + 6 * c_y[7] * weight[5] * b_pre * velocity_y[X][Ny - 2];
								pop[8][X][Ny - 2] = pop[6][(X - 1 + Nx) % Nx][Ny - 1] + 6 * c_y[8] * weight[6] * b_pre * velocity_y[X][Ny - 2];
#endif
#if defined Zou_He_Pressure
#if defined X_NOT_Periodic

								//Important: X-velocity is assumed zero
								if (X == 1)
								{
									u_Y = (velocity_y[2][Ny - 2] + velocity_y[1][Ny - 3]) / 2.0;

									pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
									pop[7][X][Ny - 2] = 1. / 2.*b_pre - pop[2][X][Ny - 2] - pop[3][X][Ny - 2] - pop[6][X][Ny - 2] - 1. / 2.*pop[0][X][Ny - 2] + 1. / 3.*b_pre* u_Y;
									pop[8][X][Ny - 2] = pop[6][X][Ny - 2] - 1. / 6.*b_pre*u_Y;
								}
								else if (X == Nx - 2)//Right Top corner
								{
									u_Y = (velocity_y[Nx - 2][Ny - 3] + velocity_y[Nx - 3][Ny - 2]) / 2.0;

									pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
									pop[7][X][Ny - 2] = pop[5][X][Ny - 2] - 1. / 6.*b_pre*u_Y;
									pop[8][X][Ny - 2] = b_pre / 2. - pop[1][X][Ny - 2] - pop[2][X][Ny - 2] - pop[5][X][Ny - 2] - 1. / 2.*pop[0][X][Ny - 2] + 1. / 3.*b_pre*u_Y;
								}
								else
								{
								u_Y = (pop[0][X][Ny - 2] + pop[1][X][Ny - 2] + pop[3][X][Ny - 2] + 2 * (pop[2][X][Ny - 2] + pop[5][X][Ny - 2] + pop[6][X][Ny - 2])) / b_pre - 1.0;
								pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
								pop[7][X][Ny - 2] = pop[5][X][Ny - 2] + 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
								pop[8][X][Ny - 2] = pop[6][X][Ny - 2] - 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
								}
#endif
#if defined X_Periodic
								u_Y = (pop[0][X][Ny - 2] + pop[1][X][Ny - 2] + pop[3][X][Ny - 2] + 2 * (pop[2][X][Ny - 2] + pop[5][X][Ny - 2] + pop[6][X][Ny - 2])) / b_pre - 1.0;
								pop[4][X][Ny - 2] = pop[2][X][Ny - 2] - 2. / 3.*b_pre*u_Y;
								pop[7][X][Ny - 2] = pop[5][X][Ny - 2] + 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
								pop[8][X][Ny - 2] = pop[6][X][Ny - 2] - 0.5*(pop[1][X][Ny - 2] - pop[3][X][Ny - 2]) - 1. / 6.*b_pre*u_Y;
#endif

#endif
								break;
					}

					case 5:{break; }

					case 6: {
								//First order
								//pop[4][X][Ny - 2] = pop[4][X][Ny - 3];
								//pop[7][X][Ny - 2] = pop[7][X][Ny - 3];
								//pop[8][X][Ny - 2] = pop[8][X][Ny - 3];

								//Second order
								pop[4][X][Ny - 2] = 2 * pop[4][X][Ny - 3] - pop[4][X][Ny - 4];
								pop[7][X][Ny - 2] = 2 * pop[7][X][Ny - 3] - pop[7][X][Ny - 4];
								pop[8][X][Ny - 2] = 2 * pop[8][X][Ny - 3] - pop[8][X][Ny - 4];
								break;
					}

					case 7:{break; }

					case 8: {
#if defined Ladd_Stationary_Wall
								pop[4][X][Ny - 2] = pop[2][X][Ny - 1];
								pop[7][X][Ny - 2] = pop[5][(X + 1) % Nx][Ny - 1];
								pop[8][X][Ny - 2] = pop[6][(X - 1 + Nx) % Nx][Ny - 1];
#endif
#if defined Zou_He_Stationary_Wall
#if defined X_NOT_Periodic

								if (X == 1)	{//Left top corner
									rhow = (density[2][Ny - 2] + density[1][Ny - 3]) / 2.0;

									pop[4][X][Ny - 2] = pop[2][X][Ny - 2];
									pop[7][X][Ny - 2] = 1. / 2.*rhow - pop[2][X][Ny - 2] - pop[3][X][Ny - 2] - pop[6][X][Ny - 2] - 1. / 2.*pop[0][X][Ny - 2];
									pop[8][X][Ny - 2] = pop[6][X][Ny - 2];
								}
								else if (X == Nx - 2)//Right Top corner
								{
									rhow = (density[Nx - 2][Ny - 3] + density[Nx - 3][Ny - 2]) / 2.0;

									pop[4][X][Ny - 2] = pop[2][X][Ny - 2] ;
									pop[7][X][Ny - 2] = pop[5][X][Ny - 2];
									pop[8][X][Ny - 2] = rhow / 2. - pop[1][X][Ny - 2] - pop[2][X][Ny - 2] - pop[5][X][Ny - 2] - 1. / 2.*pop[0][X][Ny - 2];
								}
								else
								{
									pop[4][X][Ny - 2] = pop[2][X][Ny - 1];
									pop[7][X][Ny - 2] = pop[5][(X + 1) % Nx][Ny - 1];
									pop[8][X][Ny - 2] = pop[6][(X - 1 + Nx) % Nx][Ny - 1];
								}
#endif
#if defined X_Periodic
								pop[4][X][Ny - 2] = pop[2][X][Ny - 1];
								pop[7][X][Ny - 2] = pop[5][(X + 1) % Nx][Ny - 1];
								pop[8][X][Ny - 2] = pop[6][(X - 1 + Nx) % Nx][Ny - 1];
#endif


#endif
								break;
					}
					default: cout << "mmm";
					}
				}
			}
		}//End of: for (int X = start_X; X <= end_X; ++X)
	#endif





		//// Bounce back on SOLID nodes (Just for Flow_Without_Particle, because I have not defined collision between particles and solid objects.)
		double min_dist, dist;
		int i_min;
		int x1, y1, x2, y2;
		int Xp, Yp, Xp2, Yp2;

#if defined Flow_Without_Particle
#if defined Curved_boundary_treatment
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				if (is_solid[X][Y] == FALSE){
					for (int dir = 1; dir < 9; dir++)
					{
						if (is_solid[(X - c_x[dir] + Nx) % Nx][(Y - c_y[dir] + Ny) % Ny] == TRUE){
							Xp = (X - c_x[dir] + Nx) % Nx;
							Yp = (Y - c_y[dir] + Ny) % Ny;
							Xp2 = (X + c_x[dir] + Nx) % Nx;
							Yp2 = (Y + c_y[dir] + Ny) % Ny;

							if (time == recovery_step + 1)	{ // To avoid calculating delta multiple times

								x1 = X;
								y1 = Y;
								x2 = X - c_x[dir];
								y2 = Y - c_y[dir];
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
								pop[dir][X][Y] = 1.0 / (2.0*delta[oppose[dir]][X][Y])*pop[oppose[dir]][Xp][Yp] + (2.0*delta[oppose[dir]][X][Y] - 1) / (2.0*delta[oppose[dir]][X][Y])*pop[dir][Xp2][Yp2];
							}
							else
							{
								pop[dir][X][Y] = 2.0*delta[oppose[dir]][X][Y] * pop[oppose[dir]][Xp][Yp] + (1 - 2.0*delta[oppose[dir]][X][Y]) * pop[oppose[dir]][X][Y];
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
							temp_pop = weight[oppose[dir]] * density[X][Y] * (1 + 3 * (c_x[oppose[dir]] * u_bf_x + c_y[oppose[dir]] * u_bf_y) + 4.5*SQ(c_x[oppose[dir]] * velocity_x[X][Y] + c_y[oppose[dir]] * velocity_y[X][Y]) - 1.5*(SQ(velocity_x[X][Y]) + SQ(velocity_y[X][Y])));
							pop[dir][X][Y] = (1 - zet)*pop[oppose[dir]][Xp][Yp] + zet*temp_pop;
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
							temp_pop = weight[oppose[dir]] * density[X][Y] * (1 + 3 * (c_x[oppose[dir]] * u_bf_x + c_y[oppose[dir]] * u_bf_y) + 4.5*SQ(c_x[oppose[dir]] * velocity_x[X][Y] + c_y[oppose[dir]] * velocity_y[X][Y]) - 1.5*(SQ(velocity_x[X][Y]) + SQ(velocity_y[X][Y])));
							pop[dir][X][Y] = (1 - zet)*pop[oppose[dir]][Xp][Yp] + zet*temp_pop;
#endif

							//	pop[dir][X][Y] = pop[oppose[dir]][(X - c_x[dir] + Nx) % Nx][(Y - c_y[dir] + Ny) % Ny];
						}
					} //End of: for (int dir = 1; dir < 9; dir++)
				}
			}
		}
#endif

#if defined No_curved_boundary_treatment
		for (int X = start_X; X <= end_X; ++X) {
			for (int Y = start_Y; Y <= end_Y; ++Y) {
				if (is_solid[X][Y] == FALSE){
					for (int dir = 1; dir < 9; dir++)	{
						Xp = (X - c_x[dir] + Nx) % Nx;
						Yp = (Y - c_y[dir] + Ny) % Ny;
						if (is_solid[Xp][Yp] == TRUE){
							pop[dir][X][Y] = pop[oppose[dir]][Xp][Yp];
						}
					}
				}
			}
		}
#endif
#endif
	/// Compute fluid density and velocity
	// The fluid density and velocity are obtained from the populations.
		momenta(); mkmk = 75;
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
			if (is_solid[X][Y] == FALSE)	{
			previous_velocity_magnitude[X][Y] = velocity_magnitude[X][Y];
			density[X][Y] = pop[0][X][Y] + pop[1][X][Y] + pop[2][X][Y] + pop[3][X][Y] + pop[4][X][Y] + pop[5][X][Y] + pop[6][X][Y] + pop[7][X][Y] + pop[8][X][Y];
#if defined Ladd_Pressure
			if (boundary_cond[X][Y] == 3 || boundary_cond[X][Y] == 4)
			{
				density[X][Y] = boundary_pres[X][Y];
			}
#endif

			velocity_x[X][Y] = (pop[1][X][Y] - pop[3][X][Y] + pop[5][X][Y] - pop[7][X][Y] + pop[8][X][Y] - pop[6][X][Y]) / density[X][Y];
			velocity_y[X][Y] = (pop[2][X][Y] - pop[4][X][Y] + pop[5][X][Y] - pop[7][X][Y] - pop[8][X][Y] + pop[6][X][Y]) / density[X][Y];
			velocity_magnitude[X][Y] = sqrt(SQ(velocity_x[X][Y]) + SQ(velocity_x[X][Y]));

			force_x[X][Y] = 0.0;//It is done because Eulerian forces must be zero at the beginning of each iteration
			force_y[X][Y] = 0.0;//It is done because Eulerian forces must be zero at the beginning of each iteration
			}
			if (is_solid[X][Y] == TRUE) {
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
	delete[] force_x;
	delete[] force_y;
	delete[] pop;
	delete[] density;
	delete[] velocity_x;
	delete[] velocity_y;
	delete[] pop_old;
	delete[] mdf_force_x;
	delete[] mdf_force_y;
    delete[] mdf_force_thermal;
	delete[] mdf_vel_x;
	delete[] mdf_vel_y;
    delete[] mdf_temperature;
	delete[] force_thermal;
	delete[] temperature;
	delete[] boundary_pres;
	delete[] boundary_vel_x;
	delete[] boundary_vel_y;
}
