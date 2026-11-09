//#include "stdafx.h"
#include "Thermal_solver.h"

using namespace std;


Thermal_solver::Thermal_solver()
{
}


////=============================================================================================================================
////=============================================================================================================================
////-----------------------------------------------------------------------------------------------------------------------------
//// Subroutine: INITIALIZE
////-----------------------------------------------------------------------------------------------------------------------------
////
//// The memory for lattice variables (populations, density, velocity, forces) is allocated.
//// The variables are initialized.
//
void Thermal_solver::initialize_t() {

    this->cx_pop_T  = vector<int>(stencil_directions_T, 0);
    this->cy_pop_T  = vector<int>(stencil_directions_T, 0);
    this->bounce_back_pop_T = vector<int>(stencil_directions_T, 0);

    if (stencil_directions_T == 9)
    {
        /// D2Q9 stencil:
        // index:   0  1  2  3  4  5  6  7  8
        // x:       0  1  0 -1  0  1 -1 -1  1
        // y:       0  0  1  0 -1  1  1 -1 -1
        // z:       0  0  0  0  0  0  0  0  0

        weight_T = {4. / 9.,
                    1. / 9., 1. / 9., 1. / 9., 1. / 9.,
                    1. / 36., 1. / 36., 1. / 36., 1. / 36.};

        cx_pop_T = {0, 1, 0, -1, 0, 1, -1, -1,  1};
        cy_pop_T = {0, 0, 1, 0, -1, 1,  1, -1, -1};
        bounce_back_pop_T = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    }
    if (stencil_directions_T == 5)
    {
        /// D2Q5 stencil:
        // index:   0  1  2  3  4
        // x:       0  1  0 -1  0
        // y:       0  0  1  0 -1
        weight_T = {1. / 3., 1. / 6., 1. / 6., 1. / 6., 1. / 6.};
        cx_pop_T = {0, 1, 0, -1, 0};
        cy_pop_T = {0, 0, 1, 0, -1};
        bounce_back_pop_T = {0, 3, 4, 1, 2};
    }

        for (int alpha = 0; alpha < stencil_directions_T; alpha++)
    {
        if (cx_pop_T[alpha] > 0 )
        {
            pop_to_right_T.push_back(alpha);
            pop_reverse_right_T.push_back(bounce_back_pop_T[alpha]);
        }
        if (cx_pop_T[alpha] < 0 )
        {
            pop_to_left_T.push_back(alpha);
            pop_reverse_left_T.push_back(bounce_back_pop_T[alpha]);
        }
        if (cy_pop_T[alpha] > 0)
        {
            pop_to_top_T.push_back(alpha);
            pop_reverse_top_T.push_back(bounce_back_pop_T[alpha]);
        }
        if (cy_pop_T[alpha] < 0)
        {
            pop_to_bottom_T.push_back(alpha);
            pop_reverse_bottom_T.push_back(bounce_back_pop_T[alpha]);
        }
    }


	temperature = new double*[Nx];
	force_thermal = new double*[Nx];
	previous_temperature = new double*[Nx];
	boundary_cond_temp = new int*[Nx];
	boundary_temp = new double*[Nx];
	boundary_cond_temp_2 = new int*[Nx];
	boundary_temp_2 = new double*[Nx];
	solid_thermal = new int*[Nx];

	for (int X = 0; X < Nx; ++X) {
		temperature[X] = new double[Ny];
		force_thermal[X] = new double[Ny];
		previous_temperature[X] = new double[Ny];
		boundary_temp[X] = new double[Ny];
		boundary_temp_2[X] = new double[Ny];
		boundary_cond_temp[X] = new int[Ny];
		boundary_cond_temp_2[X] = new int[Ny];
		solid_thermal[X] = new int[Ny];
	}


//	rec_step = 0;
	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			temperature[X][Y] = Fluid_init_temp;
			previous_temperature[X][Y] = Fluid_init_temp;
			force_thermal[X][Y] = 0;
			boundary_temp[X][Y] = 0;
			boundary_temp_2[X][Y] = 0;
			boundary_cond_temp[X][Y] = 0;
			boundary_cond_temp_2[X][Y] = 0;
			solid_thermal[X][Y] = 0;//All nodes are set as fluid initially

			if ((X == 1) || (X == Nx - 2))
			{
#if defined X_NOT_Periodic
				boundary_cond_temp[X][Y] = 5;//5 means zero temperature on x boundary
				boundary_cond_temp_2[X][Y] = 5;
#endif
			}

			if ((Y == 1) || (Y == Ny - 2))
			{
#if defined Y_NOT_Periodic
				boundary_cond_temp[X][Y] = 6;//6 means zero temperature on y boundary
				boundary_cond_temp_2[X][Y] = 6;
#endif
			}
		}
	}

	/// Allocate memory for the populations

	pop_T = new double**[stencil_directions_T]; // for temperature
	pop_old_T = new double**[stencil_directions_T]; // for temperature

	for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
		pop_T[c_i] = new double*[Nx];
		pop_old_T[c_i] = new double*[Nx];
		for (int X = 0; X < Nx; ++X) {
			pop_T[c_i][X] = new double[Ny];
			pop_old_T[c_i][X] = new double[Ny];
			for (int Y = 0; Y < Ny; ++Y) {
				pop_T[c_i][X][Y] = 0;
				pop_old_T[c_i][X][Y] = 0;
			}
		}
	}


	//#if defined Recovery_OFF
	/// Initialize the populations
	// Use the equilibrium populations corresponding to the initialized fluid density and velocity.
	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
			equilibrium_t(temperature[X][Y], density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);
			for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
				pop_old_T[c_i][X][Y] = pop_eq_T[c_i];
				pop_T[c_i][X][Y] = pop_eq_T[c_i];
			}
		}
	}
	//#endif

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

void Thermal_solver::equilibrium_t(double tem, double den, double vel_x, double vel_y) {
#if defined Equilibrium_Thermal_1
    for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
        pop_eq_T[c_i] = weight_T[c_i] * tem * (1.0 + 3.0 * (cx_pop_T[c_i] * vel_x + cy_pop_T[c_i] * vel_y));
    }
#endif
#if defined Equilibrium_Thermal_2
    for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
        pop_eq_T[c_i] = weight_T[c_i] * tem * den * (1.0 + 3.0 * (cx_pop_T[c_i] * vel_x + cy_pop_T[c_i] * vel_y) + 4.5 * SQ(cx_pop_T[c_i] * vel_x + cy_pop_T[c_i] * vel_y) - 1.5 * (SQ(vel_x) + SQ(vel_y)));
    }
#endif
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM_T
//-----------------------------------------------------------------------------------------------------------------------------
//
void Thermal_solver::LBM_t(int time) {

	/// Swap populations
	// The present code used old and new populations which are swapped at the beginning of each time step.
	// This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
	// This way, the old populations are not overwritten during propagation.
	// The resulting code is easier to write and to debug.
	// The memory requirement for the populations is twice as large.

	double ***swap_temp_t = pop_old_T;
	pop_old_T = pop_T;
	pop_T = swap_temp_t;

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
	double omega_thermal;
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (solid_thermal[X][Y] == 0 || solid_thermal[X][Y] == 3)		{//0 means fluid; 3 means conjugate heat transfer

				omega_thermal = omega_t;
#if defined Flow_Without_Particle
				if (solid_thermal[X][Y] == 3)	{//conjugate nodes
					omega_thermal = omega_ts;
				}
#endif
					/// Compute equilibrium
					// The equilibrium populations are computed.
					equilibrium_t(temperature[X][Y], density[X][Y], velocity_x[X][Y], velocity_y[X][Y]);

				/// Compute new populations
				// This is the lattice Boltzmann equation (combined collision and propagation).
				// Periodicity of the lattice in x-direction is taken into account by the %-operator.

#if defined Force_Thermal_1
                for (int i = 0; i < stencil_directions_T; i++)
                {
                    int newX = (X + cx_pop_T[i] + Nx) % Nx;
                    int newY = (Y + cy_pop_T[i] + Ny) % Ny;
                    pop_T[i][newX][newY] = pop_old_T[i][X][Y] * (1.0 - omega_thermal) + pop_eq_T[i] * omega_thermal + weight_T[i] * force_thermal[X][Y];;
                }
#endif
#if defined Force_Thermal_2
                for (int i = 0; i < stencil_directions_T; i++)
                {
                    int newX = (X + cx_pop_T[i] + Nx) % Nx;
                    int newY = (Y + cy_pop_T[i] + Ny) % Ny;
                    pop_T[i][newX][newY] = pop_old_T[i][X][Y] * (1.0 - omega_thermal) + pop_eq_T[i] * omega_thermal + (1.0 - 0.5 * omega_thermal) * weight_T[i] * force_thermal[X][Y];
                }
#endif
			}
		}
	}

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: BC_T
//-----------------------------------------------------------------------------------------------------------------------------

void Thermal_solver::BC_t(int time) {

    int surf_normal_X{}, surf_normal_Y{}, surf_normal_Z{};
    double temperature_at_bc{};

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
                pop_to      = pop_to_right_T;
                pop_reverse = pop_reverse_right_T;
                start_BC_X  = Nx - 2;
                end_BC_X    = Nx - 2;
                start_BC_Y  = 1;
                end_BC_Y    = Ny - 2;
                tag = ii;
                //cout << "Right, processor:" << processor_id << " "<<start_BC_X<< "-"<<end_BC_X << ", "<< start_BC_Y<< "-"<<end_BC_Y << ", "<< start_BC_Z<< "-"<<end_BC_Z << endl;
                #if defined Right_Fixed_Temp
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_right;
                #endif
                #if defined Right_Adiabatic
                    BC_type = "adiabatic";
                    surf_normal_X = -1;
                    surf_normal_Y = 0;
                #endif
        break;
        case 1:
            /// Left BC
                pop_to      = pop_to_left_T;
                pop_reverse = pop_reverse_left_T;
                start_BC_X  = 1;
                end_BC_X    = 1;
                start_BC_Y  = 1;
                end_BC_Y    = Ny - 2;
                tag = ii;
                #if defined Left_Fixed_Temp
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_left;
                #endif
                #if defined Left_Adiabatic
                    BC_type = "adiabatic";
                    surf_normal_X = 1;
                    surf_normal_Y = 0;
                #endif
        break;
        case 2:
            /// Top BC
                pop_to      = pop_to_top_T;
                pop_reverse = pop_reverse_top_T;
                start_BC_Y  = Ny - 2;
                end_BC_Y    = Ny - 2;
                start_BC_X  = 1;
                end_BC_X    = Nx - 2;
                tag = ii;
                #if defined Top_Fixed_Temp
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_top;
                #endif
                #if defined Top_Adiabatic
                    BC_type = "adiabatic";
                    surf_normal_X = 0;
                    surf_normal_Y = -1;
                #endif
            break;
        case 3:
            /// Bottom BC
                pop_to      = pop_to_bottom_T;
                pop_reverse = pop_reverse_bottom_T;
                start_BC_Y  = 1;
                end_BC_Y    = 1;
                start_BC_X  = 1;
                end_BC_X    = Nx - 2;
                tag = ii;
                #if defined Bottom_Fixed_Temp
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_bot;
                #endif
                #if defined  Bottom_Adiabatic
                    BC_type = "adiabatic";
                    surf_normal_X = 0;
                    surf_normal_Y = 1;
                #endif
            break;

        default:
            std::cout << "Invalid choice of BC Thermal" << std::endl;
            break;
        }

        for (int alpha = 0; alpha < pop_to.size(); alpha++)
        {
            for (int X = start_BC_X; X <= end_BC_X; ++X)
            {
                for (int Y = start_BC_Y; Y <= end_BC_Y; ++Y)
                {
                    if (BC_type == "fixed_temperature"){
                        int pop_id  = pop_to[alpha];
                        int pop_rev = pop_reverse[alpha];
                        pop_T[pop_rev][X][Y] = (weight_T[pop_rev] + weight_T[pop_id]) * temperature_at_bc - pop_T[pop_id][X + cx_pop_T[pop_id]][Y + cy_pop_T[pop_id]]; //Drichlet BC
                    }
                    if (BC_type == "adiabatic"){
                        int pop_id  = pop_to[alpha];
                        int pop_rev = pop_reverse[alpha];
                        pop_T[pop_rev][X][Y] = pop_T[pop_rev][X + surf_normal_X][Y + surf_normal_Y];//Neumann BC (first order)
                    }
                }
            }
        }
    }

///##############################################################################################///
#if defined Flow_Without_Particle
    	int Xp, Yp, Xp2,Yp2;
#if defined No_curved_boundary_treatment
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE){
				for (int dir = 1; dir < stencil_directions_T; dir++)	{
					Xp = (X - cx_pop_T[dir] + Nx) % Nx;
					Yp = (Y - cy_pop_T[dir] + Ny) % Ny;
					Xp2 = (X + cx_pop_T[dir] + Nx) % Nx;
					Yp2 = (Y + cy_pop_T[dir] + Ny) % Ny;
					if (is_solid[Xp][Yp] == TRUE){
						if (solid_thermal[Xp][Yp] == 1){//constant temp.
							pop_T[dir][X][Y] = (weight_T[dir] + weight_T[oppose[dir]])*solid_temperature - pop_T[oppose[dir]][Xp][Yp];
						}
						if (solid_thermal[Xp][Yp] == 2){//adiabatic
							pop_T[dir][X][Y] = pop_T[dir][Xp2][Yp2];
						}
					}
				}
			}
		}
	}
#endif
	double tb1, tb2, tbstar, ubstar_x, ubstar_y;
	double u_bf_x, u_bf_y, pop_eqff;
	double u_bff_x, u_bff_y, h_oppose_eq, h_oppose_noneq;

#if defined Curved_boundary_treatment
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE){
				for (int dir = 1; dir < stencil_directions_T; dir++)	{
					Xp = (X - cx_pop_T[dir] + Nx) % Nx;//wall
					Yp = (Y - cy_pop_T[dir] + Ny) % Ny;//wall
					Xp2 = (X + cx_pop_T[dir] + Nx) % Nx;//ff
					Yp2 = (Y + cy_pop_T[dir] + Ny) % Ny;//ff
					if (is_solid[Xp][Yp] == TRUE){

						if (solid_thermal[Xp][Yp] == 1){//constant temp.

							u_bf_x = (delta[oppose[dir]][X][Y] - 1) * velocity_x[X][Y] / delta[oppose[dir]][X][Y];
							u_bf_y = (delta[oppose[dir]][X][Y] - 1) * velocity_y[X][Y] / delta[oppose[dir]][X][Y];

							u_bff_x = ((delta[oppose[dir]][X][Y] - 1) * velocity_x[Xp2][Yp2]) / (1 + delta[oppose[dir]][X][Y]);
							u_bff_y = ((delta[oppose[dir]][X][Y] - 1) * velocity_y[Xp2][Yp2]) / (1 + delta[oppose[dir]][X][Y]);


							tb1 = (solid_temperature + (delta[oppose[dir]][X][Y] - 1)*temperature[X][Y]) / delta[oppose[dir]][X][Y];
							tb2 = (2 * solid_temperature + (delta[oppose[dir]][X][Y] - 1)*temperature[Xp2][Yp2]) / (1 + delta[oppose[dir]][X][Y]);
							if (delta[oppose[dir]][X][Y] >= 0.75)
							{
								tbstar = tb1;
								ubstar_x = u_bf_x;
								ubstar_y = u_bf_y;
								h_oppose_eq = weight_T[dir] * tbstar * (1 + 3 * (cx_pop_T[dir] * ubstar_x + cy_pop_T[dir] * ubstar_y));
								h_oppose_noneq = (pop_T[dir][X][Y] - h_oppose_eq) / (1 - omega_t);//pop_T[dir][X][Y] is the post collision one
							}
							else
							{
								tbstar = delta[oppose[dir]][X][Y]*tb1 + (1 - delta[oppose[dir]][X][Y])*tb2;
								ubstar_x = u_bf_x + (1 - delta[oppose[dir]][X][Y])*u_bff_x;
								ubstar_y = u_bf_y + (1 - delta[oppose[dir]][X][Y])*u_bff_y;
								h_oppose_eq = weight_T[dir] * tbstar * (1 + 3 * (cx_pop_T[dir] * ubstar_x + cy_pop_T[dir] * ubstar_y));
								pop_eqff = weight_T[dir] * temperature[Xp2][Yp2] * (1 + 3 * (cx_pop_T[dir] * velocity_x[Xp2][Yp2] + cy_pop_T[dir] * velocity_y[Xp2][Yp2]));
								h_oppose_noneq = delta[oppose[dir]][X][Y] * ((pop_T[dir][X][Y] - h_oppose_eq) / (1 - omega_t)) + (1 - delta[oppose[dir]][X][Y])*(pop_T[dir][Xp2][Yp2] - pop_eqff);
							}

							pop_T[dir][X][Y] = h_oppose_eq + h_oppose_noneq;
						}

						if (solid_thermal[Xp][Yp] == 2){//adiabatic
							pop_T[dir][X][Y] = pop_T[dir][Xp2][Yp2];
						}
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
// Subroutine: BC_T_OLD_VERSION
//-----------------------------------------------------------------------------------------------------------------------------

void Thermal_solver::BC_t_OLD_VERSION(int time) {

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

	int Xtest;
	Xtest = 0; //This parameter is just for avoiding applying X-direc boundary condition multiple times. Therefore it just improves computational time and nothing else

	int cond, bou;
	double b_temp;//temperature boundary conditions
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
				bou = boundary_cond_temp[1][Y];
				b_temp = boundary_temp[1][Y];


				if (i == 2)	{
					bou = boundary_cond_temp_2[1][Y];
					b_temp = boundary_temp_2[1][Y];
				}

				switch (bou)//Left side
				{
				case 1:		{//Constant Temperature condition along X
				    if (stencil_directions_T == 9)
                    {
                        pop_T[1][1][Y] = (weight_T[1] + weight_T[3])*b_temp - pop_T[3][0][Y];
                        pop_T[5][1][Y] = (weight_T[5] + weight_T[7])*b_temp - pop_T[7][0][(Y - 1 + Ny) % Ny];
                        pop_T[8][1][Y] = (weight_T[8] + weight_T[6])*b_temp - pop_T[6][0][(Y + 1) % Ny];
                    }else if (stencil_directions_T == 5){
                        pop_T[1][1][Y] = (weight_T[1] + weight_T[3])*b_temp - pop_T[3][0][Y];
                    }
                    break;
				}

				case 2:{break; }//Constant Temperature condition along Y

				case 3: {//Adiabatic condition along X
                    if (stencil_directions_T == 9)
                    {
                        pop_T[1][1][Y] = pop_T[1][2][Y];
                        pop_T[5][1][Y] = pop_T[5][2][Y];
                        pop_T[8][1][Y] = pop_T[8][2][Y];

                        //Second order
                        //pop_T[1][1][Y] = 2 * pop_T[1][2][Y] - pop_T[1][3][Y];
                        //pop_T[5][1][Y] = 2 * pop_T[5][2][Y] - pop_T[1][3][Y];
                        //pop_T[8][1][Y] = 2 * pop_T[8][2][Y] - pop_T[1][3][Y];
                    }else if (stencil_directions_T == 5){
                        pop_T[1][1][Y] = pop_T[1][2][Y];
                    }
                    break;
				}
				case 4:{break; }//Adiabatic condition along Y


				case 5: {//Zero Temperature condition along X
                    if (stencil_directions_T == 9)
                    {
                        pop_T[1][1][Y] = (weight_T[1] + weight_T[3])* 0.0 - pop_T[3][0][Y];
                        pop_T[5][1][Y] = (weight_T[5] + weight_T[7])* 0.0 - pop_T[7][0][(Y - 1 + Ny) % Ny];
                        pop_T[8][1][Y] = (weight_T[8] + weight_T[6])* 0.0 - pop_T[6][0][(Y + 1) % Ny];
                    }else if (stencil_directions_T == 5){
                        pop_T[1][1][Y] = (weight_T[1] + weight_T[3])* 0.0 - pop_T[3][0][Y];
                    }
                    break;
				}

				case 6:{break; }//Zero Temperature condition along Y

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
				bou = boundary_cond_temp[Nx - 2][Y];
				b_temp = boundary_temp[Nx - 2][Y];

				if (i == 2)	{
					bou = boundary_cond_temp_2[Nx - 2][Y];
					b_temp = boundary_temp_2[Nx - 2][Y];
				}

				switch (bou)//Right side
				{

				case 1:		{//Constant Temperature condition along X
                    if (stencil_directions_T == 9)
                    {
                        pop_T[3][Nx - 2][Y] = (weight_T[3] + weight_T[1]) * b_temp - pop_T[1][Nx - 1][Y];
                        pop_T[7][Nx - 2][Y] = (weight_T[7] + weight_T[5]) * b_temp - pop_T[5][Nx - 1][(Y + 1) % Ny];
                        pop_T[6][Nx - 2][Y] = (weight_T[6] + weight_T[8]) * b_temp - pop_T[8][Nx - 1][(Y - 1 + Ny) % Ny];
                    }else if (stencil_directions_T == 5){
                        pop_T[3][Nx - 2][Y] = (weight_T[3] + weight_T[1]) * b_temp - pop_T[1][Nx - 1][Y];
                    }
                    break;
				}

				case 2:{break; }//Constant Temperature condition along Y

				case 3: {//Adiabatic condition along X
                    if (stencil_directions_T == 9)
                    {
                        //pop_T[3][Nx - 2][Y] = pop_T[3][Nx - 3][Y];
                        //pop_T[7][Nx - 2][Y] = pop_T[7][Nx - 3][Y];
                        //pop_T[6][Nx - 2][Y] = pop_T[6][Nx - 3][Y];

                        //Second order
                        pop_T[3][Nx - 2][Y] = 2 * pop_T[3][Nx - 3][Y] - pop_T[3][Nx - 4][Y];
                        pop_T[7][Nx - 2][Y] = 2 * pop_T[7][Nx - 3][Y] - pop_T[7][Nx - 4][Y];
                        pop_T[6][Nx - 2][Y] = 2 * pop_T[6][Nx - 3][Y] - pop_T[6][Nx - 4][Y];
                    }else if (stencil_directions_T == 5){
                        pop_T[3][Nx - 2][Y] = 2 * pop_T[3][Nx - 3][Y] - pop_T[3][Nx - 4][Y];
                    }
                    break;
				}

				case 4:{break; }//Adiabatic condition along Y

				case 5: {//Zero Temperature condition along X
                    if (stencil_directions_T == 9)
                    {
                        pop_T[3][Nx - 2][Y] = (weight_T[3] + weight_T[1])*0 - pop_T[1][Nx - 1][Y];
                        pop_T[7][Nx - 2][Y] = (weight_T[7] + weight_T[5])*0 - pop_T[5][Nx - 1][(Y + 1) % Ny];
                        pop_T[6][Nx - 2][Y] = (weight_T[6] + weight_T[8])*0 - pop_T[8][Nx - 1][(Y - 1 + Ny) % Ny];
                    }else if (stencil_directions_T == 5){
                        pop_T[3][Nx - 2][Y] = (weight_T[3] + weight_T[1])*0 - pop_T[1][Nx - 1][Y];
                    }
                    break;
				}

				case 6:{break; }//Zero Temperature condition along Y

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
				bou = boundary_cond_temp[X][1];
				b_temp = boundary_temp[X][1];

				if (i == 2)	{
					bou = boundary_cond_temp_2[X][1];
					b_temp = boundary_temp_2[X][1];
				}

				switch (bou)
				{

				case 1:{break; }//Constant Temperature condition along X

				case 2:	{//Constant Temperature condition along Y
                    if (stencil_directions_T == 9)
                    {
                        pop_T[2][X][1] = (weight_T[2] + weight_T[4])*b_temp - pop_T[4][X][0];
                        pop_T[5][X][1] = (weight_T[5] + weight_T[7])*b_temp - pop_T[7][(X - 1 + Nx) % Nx][0];
                        pop_T[6][X][1] = (weight_T[6] + weight_T[8])*b_temp - pop_T[8][(X + 1) % Nx][0];
                    }else if (stencil_directions_T == 5){
                        pop_T[2][X][1] = (weight_T[2] + weight_T[4])*b_temp - pop_T[4][X][0];
                    }
                    break;
				}

				case 3:{break; }//Adiabatic condition along X

				case 4: {//Adiabatic condition along Y
                    if (stencil_directions_T == 9)
                    {
                        pop_T[2][X][1] = pop_T[2][X][2];
                        pop_T[5][X][1] = pop_T[5][X][2];
                        pop_T[6][X][1] = pop_T[6][X][2];

                        //Second order
                        //pop_T[2][X][1] = 2 * pop_T[2][X][2] - pop_T[2][X][3];
                        //pop_T[5][X][1] = 2 * pop_T[5][X][2] - pop_T[5][X][3];
                        //pop_T[6][X][1] = 2 * pop_T[6][X][2] - pop_T[6][X][3];
                    }else if (stencil_directions_T == 5){
                        pop_T[2][X][1] = pop_T[2][X][2];
                    }
                    break;
				}

				case 5:{break; }//Zero Temperature condition along X

				case 6: {//Zero Temperature condition along Y
                    if (stencil_directions_T == 9)
                    {
                        pop_T[2][X][1] = (weight_T[2] + weight_T[4])*0 - pop_T[4][X][0];
                        pop_T[5][X][1] = (weight_T[5] + weight_T[7])*0 - pop_T[7][(X - 1 + Nx) % Nx][0];
                        pop_T[6][X][1] = (weight_T[6] + weight_T[8])*0 - pop_T[8][(X + 1) % Nx][0];
                    }else if (stencil_directions_T == 5){
                        pop_T[2][X][1] = (weight_T[2] + weight_T[4])*0 - pop_T[4][X][0];
                    }
                    break;
				}
				default: cout << "nnn";
				}
			}
		}

		/// TOP surface
		if (is_solid[X][Ny - 2] == FALSE){
			cond = 1;
			if (edge_boundary[X][Ny - 2] == TRUE) {
				cond = 2;
			}

			for (size_t i = 1; i <= cond; i++)
			{
				bou = boundary_cond_temp[X][Ny - 2];
				b_temp = boundary_temp[X][Ny - 2];

				if (i == 2)	{
					bou = boundary_cond_temp_2[X][Ny - 2];
					b_temp = boundary_temp_2[X][Ny - 2];
				}

				switch (bou)
				{

				case 1:{break; }

				case 2:	{
                    if (stencil_directions_T == 9)
                    {
                        pop_T[4][X][Ny - 2] = (weight_T[4] + weight_T[2])*b_temp - pop_T[2][X][Ny - 1];
                        pop_T[7][X][Ny - 2] = (weight_T[7] + weight_T[5])*b_temp - pop_T[5][(X + 1) % Nx][Ny - 1];
                        pop_T[8][X][Ny - 2] = (weight_T[8] + weight_T[6])*b_temp - pop_T[6][(X - 1 + Nx) % Nx][Ny - 1];
                    }else if (stencil_directions_T == 5){
                        pop_T[4][X][Ny - 2] = (weight_T[4] + weight_T[2])*b_temp - pop_T[2][X][Ny - 1];
                    }
                    break;
				}

				case 3:{break; }//Adiabatic condition along X

				case 4: {//Adiabatic condition along Y
                    if (stencil_directions_T == 9)
                    {
                        pop_T[4][X][Ny - 2] = pop_T[4][X][Ny - 3];
                        pop_T[7][X][Ny - 2] = pop_T[7][X][Ny - 3];
                        pop_T[8][X][Ny - 2] = pop_T[8][X][Ny - 3];

                        //Second order
                        //pop_T[4][X][Ny - 2] = 2 * pop_T[4][X][Ny - 3] - pop_T[4][X][Ny - 4];
                        //pop_T[7][X][Ny - 2] = 2 * pop_T[7][X][Ny - 3] - pop_T[7][X][Ny - 4];
                        //pop_T[8][X][Ny - 2] = 2 * pop_T[8][X][Ny - 3] - pop_T[8][X][Ny - 4];
                    }else if (stencil_directions_T == 5){
                        pop_T[4][X][Ny - 2] = pop_T[4][X][Ny - 3];
                    }
                    break;
				}

				case 5:{break; }//Zero Temperature condition along X

				case 6: {//Zero Temperature condition along Y
                    if (stencil_directions_T == 9)
                    {
                        pop_T[4][X][Ny - 2] = (weight_T[4] + weight_T[2])*0 - pop_T[2][X][Ny - 1];
                        pop_T[7][X][Ny - 2] = (weight_T[7] + weight_T[5])*0 - pop_T[5][(X + 1) % Nx][Ny - 1];
                        pop_T[8][X][Ny - 2] = (weight_T[8] + weight_T[6])*0 - pop_T[6][(X - 1 + Nx) % Nx][Ny - 1];
                    }else if (stencil_directions_T == 5){
                        pop_T[4][X][Ny - 2] = (weight_T[4] + weight_T[2])*0 - pop_T[2][X][Ny - 1];
                    }
                    break;
				}
				default: cout << "mmm";
				}
			}
		}
	}//End of: for (int X = start_X; X <= end_X; ++X)
#endif

	int Xp, Yp, Xp2,Yp2;

#if defined Flow_Without_Particle
#if defined No_curved_boundary_treatment
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE){
				for (int dir = 1; dir < stencil_directions_T; dir++)	{
					Xp = (X - cx_pop_T[dir] + Nx) % Nx;
					Yp = (Y - cy_pop_T[dir] + Ny) % Ny;
					Xp2 = (X + cx_pop_T[dir] + Nx) % Nx;
					Yp2 = (Y + cy_pop_T[dir] + Ny) % Ny;
					if (is_solid[Xp][Yp] == TRUE){
						if (solid_thermal[Xp][Yp] == 1){//constant temp.
							pop_T[dir][X][Y] = (weight_T[dir] + weight_T[oppose[dir]])*solid_temperature - pop_T[oppose[dir]][Xp][Yp];
						}
						if (solid_thermal[Xp][Yp] == 2){//adiabatic
							pop_T[dir][X][Y] = pop_T[dir][Xp2][Yp2];
						}
					}
				}
			}
		}
	}
#endif
	double tb1, tb2, tbstar, ubstar_x, ubstar_y;
	double u_bf_x, u_bf_y, pop_eqff;
	double u_bff_x, u_bff_y, h_oppose_eq, h_oppose_noneq;

#if defined Curved_boundary_treatment
	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			if (is_solid[X][Y] == FALSE){
				for (int dir = 1; dir < stencil_directions_T; dir++)	{
					Xp = (X - cx_pop_T[dir] + Nx) % Nx;//wall
					Yp = (Y - cy_pop_T[dir] + Ny) % Ny;//wall
					Xp2 = (X + cx_pop_T[dir] + Nx) % Nx;//ff
					Yp2 = (Y + cy_pop_T[dir] + Ny) % Ny;//ff
					if (is_solid[Xp][Yp] == TRUE){

						if (solid_thermal[Xp][Yp] == 1){//constant temp.

							u_bf_x = (delta[oppose[dir]][X][Y] - 1) * velocity_x[X][Y] / delta[oppose[dir]][X][Y];
							u_bf_y = (delta[oppose[dir]][X][Y] - 1) * velocity_y[X][Y] / delta[oppose[dir]][X][Y];

							u_bff_x = ((delta[oppose[dir]][X][Y] - 1) * velocity_x[Xp2][Yp2]) / (1 + delta[oppose[dir]][X][Y]);
							u_bff_y = ((delta[oppose[dir]][X][Y] - 1) * velocity_y[Xp2][Yp2]) / (1 + delta[oppose[dir]][X][Y]);


							tb1 = (solid_temperature + (delta[oppose[dir]][X][Y] - 1)*temperature[X][Y]) / delta[oppose[dir]][X][Y];
							tb2 = (2 * solid_temperature + (delta[oppose[dir]][X][Y] - 1)*temperature[Xp2][Yp2]) / (1 + delta[oppose[dir]][X][Y]);
							if (delta[oppose[dir]][X][Y] >= 0.75)
							{
								tbstar = tb1;
								ubstar_x = u_bf_x;
								ubstar_y = u_bf_y;
								h_oppose_eq = weight_T[dir] * tbstar * (1 + 3 * (cx_pop_T[dir] * ubstar_x + cy_pop_T[dir] * ubstar_y));
								h_oppose_noneq = (pop_T[dir][X][Y] - h_oppose_eq) / (1 - omega_t);//pop_T[dir][X][Y] is the post collision one
							}
							else
							{
								tbstar = delta[oppose[dir]][X][Y]*tb1 + (1 - delta[oppose[dir]][X][Y])*tb2;
								ubstar_x = u_bf_x + (1 - delta[oppose[dir]][X][Y])*u_bff_x;
								ubstar_y = u_bf_y + (1 - delta[oppose[dir]][X][Y])*u_bff_y;
								h_oppose_eq = weight_T[dir] * tbstar * (1 + 3 * (cx_pop_T[dir] * ubstar_x + cy_pop_T[dir] * ubstar_y));
								pop_eqff = weight_T[dir] * temperature[Xp2][Yp2] * (1 + 3 * (cx_pop_T[dir] * velocity_x[Xp2][Yp2] + cy_pop_T[dir] * velocity_y[Xp2][Yp2]));
								h_oppose_noneq = delta[oppose[dir]][X][Y] * ((pop_T[dir][X][Y] - h_oppose_eq) / (1 - omega_t)) + (1 - delta[oppose[dir]][X][Y])*(pop_T[dir][Xp2][Yp2] - pop_eqff);
							}

							pop_T[dir][X][Y] = h_oppose_eq + h_oppose_noneq;
						}

						if (solid_thermal[Xp][Yp] == 2){//adiabatic
							pop_T[dir][X][Y] = pop_T[dir][Xp2][Yp2];
						}
					}
				}
			}
		}
	}
#endif
#endif



	/// Compute fluid temperature
	// The fluid temperature is obtained from the populations.
	momenta_t();

	return;*/
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MOMENTA_T
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the fluid temperature from the populations.
void Thermal_solver::momenta_t() {

	int start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
	start_X = 0;
	end_X = Nx - 1;

	int start_Y, end_Y;
	start_Y = 0;
	end_Y = Ny - 1;

	for (int X = start_X; X <= end_X; ++X) {
		for (int Y = start_Y; Y <= end_Y; ++Y) {
			previous_temperature[X][Y] = temperature[X][Y];
            temperature[X][Y] = 0.0;
#if defined Equilibrium_Thermal_1
            for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
                temperature[X][Y] += pop_T[c_i][X][Y];
            }
#endif
#if defined Equilibrium_Thermal_2
            for (int c_i = 0; c_i < stencil_directions_T; ++c_i) {
                temperature[X][Y] += pop_T[c_i][X][Y];
            }
            temperature[X][Y] /= density[X][Y];
#endif
			force_thermal[X][Y] = 0.0;
		}
	}

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: CHECK_RESIDUAL
//-----------------------------------------------------------------------------------------------------------------------------
#if defined Stop_by_Residual
bool Thermal_solver::check_residual(){
	bool st_cr = false;
	double temp_criteria = 0;
	double temp_criteria_1 = 0;

	for (int X = 0; X < Nx; ++X) {
		for (int Y = 0; Y < Ny; ++Y) {
		//	if (is_solid[X][Y] == FALSE)			{

					temp_criteria += abs(temperature[X][Y] - previous_temperature[X][Y]);
					temp_criteria_1 += abs(temperature[X][Y]);

		//	}
		}
	}
	if ((temp_criteria / temp_criteria_1)<residual_thermal)	{
		st_cr = true;
		return st_cr;
	}
	return st_cr;
}
#endif
//=============================================================================================================================
//=============================================================================================================================
Thermal_solver::~Thermal_solver()
{
}
