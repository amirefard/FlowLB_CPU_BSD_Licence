//#include "stdafx.h"
#include "Inlet_outlets.h"

using namespace std;

Inlet_outlets::Inlet_outlets()
{

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SET_BOUNDARIES
//-----------------------------------------------------------------------------------------------------------------------------
//
// Set the location of inlet and outlets

void Inlet_outlets::set_boundaries(){

#if defined X_NOT_Periodic || defined Y_NOT_Periodic

    string output_filename("LB_Results/parameters");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::app);

	int b_type;
	int low_y, high_y, low_x, high_x;
	double vel_x, vel_y, pres;
	int X, Y;
	int bbtype[4];
	bbtype[0] = Left_in_out;
	bbtype[1] = Right_in_out;
	bbtype[2] = Bot_in_out;
	bbtype[3] = Top_in_out;
	int start_i = 0;
	int end_i = 3;

	int Xmin=1;
	int Xmax=Nx-2;
	int Ymin = 1;
	int Ymax = Ny - 2;

#if defined X_Periodic
	start_i = 2;
	Xmin = 0;
	Xmax = Nx - 1;
#endif
#if defined Y_Periodic
	end_i = 1;
	Ymin = 0;
	Ymax = Ny - 1;
#endif
	std::cout << endl;
	string st;

	for (size_t kk = start_i; kk <= end_i; kk++)	{
		switch (kk)
		{
		case 0: {st = "X=1"; break; }
		case 1: {st = "X=Nx-2"; break; }
		case 2: {st = "Y=1"; break; }
		case 3: {st = "Y=Ny-2"; break; }
		default: st = "X=0";
		}
		std::cout << "Number of Inlet/Outlet at " << st << " is: " << bbtype[kk] << endl;

		for (size_t i = 0; i < bbtype[kk]; i++)
		{
			std::cout << "Please specify the boundary type of ";
			switch (((i + 1) % 10))
			{
			case 1: {std::cout << "1st inlet/outlet at " << st << endl; break; }
			case 2: {std::cout << "2nd inlet/outlet at " << st << endl; break; }
			case 3: {std::cout << "3rd inlet/outlet at " << st << endl; break; }
			default: std::cout << i + 1 << "th inlet/outlet at " << st << endl;
			}
			std::cout << "1) Velocity 2) Pressure 3) Zero Gradient: ";
			std::cin >> b_type;

			while (b_type > 3 || b_type < 1)
			{
				std::cout << "1) Velocity 2) Pressure 3) Zero Gradient: ";
				std::cin >> b_type;
			}

			if (kk < 2)
			{
				std::cout << "Please specify lower Y location of the boundary (" << Ymin << "-" << Ymax << "): ";
				std::cin >> low_y;
				std::cout << "Please specify higher Y location of the boundary (" << Ymin << "-" << Ymax << "): ";
				std::cin >> high_y;
			}
			else
			{
				std::cout << "Please specify lower X location of the boundary (" << Xmin << "-" << Xmax << "): ";
				std::cin >> low_x;
				std::cout << "Please specify higher X location of the boundary (" << Xmin << "-" << Xmax << "): ";
				std::cin >> high_x;
			}

			if (b_type == 1)
			{
				std::cout << "Enter velocity_X value: ";
				std::cin >> vel_x;
				std::cout << "Enter velocity_Y value: ";
				std::cin >> vel_y;
				output_file << "Boundary No. " << i + 1 << " at " << st << ": Velocity" << endl;
				output_file << "Vel_x: " << vel_x << ", Vel_y: " << vel_y << endl;
			}
			if (b_type == 2)
			{
				std::cout << "Enter Pressure value: ";
				std::cin >> pres;
				output_file << "Boundary No. " << i + 1 << " at " << st << ": Pressure" << endl;
				output_file << "Pressure: " << pres << endl;
			}


			if (b_type == 3)
			{
				output_file << "Boundary No. " << i + 1 << " at " << st << ": Zero Gradient" << endl;
			}

			if (kk < 2)
			{
				output_file << "Lower Y value: " << low_y << ", Higher Y value: " << high_y << endl;
			}
			else
			{
				output_file << "Lower X value: " << low_x << ", Higher X value: " << high_x << endl << endl;
			}

			switch (kk)
			{

			case 0: {X = 1; break; }
			case 1: {X = Nx - 2; break; }//Nx-2 because there is a buffer layer
			case 2: {Y = 1; break; }
			case 3: {Y = Ny - 2; break; }
			default: { cout << "uizu"; }
			}

			if ((kk == 0) || (kk == 1))
			{
				for (int io = low_y; io <= high_y; io++)
				{

					if (b_type == 1)
					{
						boundary_cond[X][io] = 1;
						boundary_cond_2[X][io] = 1;
						boundary_vel_x[X][io] = vel_x;
						boundary_vel_y[X][io] = vel_y;
					}
					if (b_type == 2)
					{
						boundary_cond[X][io] = 3;
						boundary_cond_2[X][io] = 3;
						boundary_pres[X][io] = pres;
					}
					if (b_type == 3)
					{
						boundary_cond[X][io] = 5;//5 means zero gradient in X direction
						boundary_cond_2[X][io] = 5;
					//	boundary_pres[X][io] = pres;
					}
				}
			}

			if ((kk == 2) || (kk == 3))
			{
				for (int oi = low_x; oi <= high_x; oi++)
				{

					if (b_type == 1)
					{
						boundary_cond[oi][Y] = 2;
						boundary_cond_2[oi][Y] = 2;
						boundary_vel_x[oi][Y] = vel_x;
						boundary_vel_y[oi][Y] = vel_y;
					}
					if (b_type == 2)
					{
						boundary_cond[oi][Y] = 4;
						boundary_cond_2[oi][Y] = 4;
						boundary_pres[oi][Y] = pres;
					}
					if (b_type == 3)
					{
						boundary_cond[oi][Y] = 6;//6 means zero gradient in Y direction
						boundary_cond_2[oi][Y] = 6;
					}
				}
			}

			std::cout << endl;
		}
	}

#if defined X_NOT_Periodic && defined Y_NOT_Periodic
		int Y1 = 1;
		int Y2 = Ny - 2;

		for (size_t ik = 0; ik <= 1; ik++)
		{
			int XX = Nx - 2;
			if (ik == 1)
			{
				XX = 1;
			}

			if (boundary_cond[XX][Y1] != boundary_cond[XX][Y1+1])
			{
				edge_boundary[XX][Y1] = TRUE;
				boundary_cond_2[XX][Y1] = boundary_cond[XX][Y1 + 1];

				if (boundary_cond[XX][Y1 + 1] == 1)
				{
					boundary_vel_x_2[XX][Y1] = boundary_vel_x[XX][Y1 + 1];
					boundary_vel_y_2[XX][Y1] = boundary_vel_y[XX][Y1 + 1];
				}
				if (boundary_cond[XX][Y1 + 1] == 3)
				{
					boundary_pres_2[XX][Y1] = boundary_pres[XX][Y1 + 1];
				}
			}

			if (boundary_cond[XX][Y2] != boundary_cond[XX][Y2-1])
			{
				edge_boundary[XX][Y2] = TRUE;
				boundary_cond_2[XX][Y2] = boundary_cond[XX][Y2-1];

				if (boundary_cond[XX][Y2 - 1] == 1)
				{
					boundary_vel_x_2[XX][Y2] = boundary_vel_x[XX][Y2 - 1];
					boundary_vel_y_2[XX][Y2] = boundary_vel_y[XX][Y2 - 1];
				}
				if (boundary_cond[XX][Y2 - 1] == 3)
				{
					boundary_pres_2[XX][Y2] = boundary_pres[XX][Y2 - 1];
				}
			}
		}
#endif

#if defined Y_NOT_Periodic && defined X_NOT_Periodic
		for (size_t ik = 0; ik <= 1; ik++)
		{
			int YY = Ny - 2;
			if (ik == 1)
			{
				YY = 1;
			}

			if (boundary_cond[1][YY] != boundary_cond[2][YY])
			{
				edge_boundary[1][YY] = TRUE;
				boundary_cond_2[1][YY] = boundary_cond[2][YY];

				if (boundary_cond[2][YY] == 2)
				{
					boundary_vel_x_2[1][YY] = boundary_vel_x[2][YY];
					boundary_vel_y_2[1][YY] = boundary_vel_y[2][YY];
				}

				if (boundary_cond[2][YY] == 4)
				{
					boundary_pres_2[1][YY] = boundary_pres[2][YY];
				}
			}

			if (boundary_cond[Nx - 2][YY] != boundary_cond[Nx - 3][YY])	{
				edge_boundary[Nx - 2][YY] = TRUE;
				boundary_cond_2[Nx - 2][YY] = boundary_cond[Nx - 3][YY];

				if (boundary_cond[Nx - 3][YY] == 2)
				{
					boundary_vel_x_2[Nx - 2][YY] = boundary_vel_x[Nx - 3][YY];
					boundary_vel_y_2[Nx - 2][YY] = boundary_vel_y[Nx - 3][YY];
				}
				if (boundary_cond[Nx - 3][YY] == 4)
				{
					boundary_pres_2[Nx - 2][YY] = boundary_pres[Nx - 3][YY];
				}
			}

		}
#endif

	std::cout << endl;

	output_file.close();

#endif
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SET_BOUNDARIES_THERMAL
//-----------------------------------------------------------------------------------------------------------------------------
//
// Set the location of temperature boundary condition

void Inlet_outlets::set_boundaries_thermal(){
#if defined Flow_With_Thermal_Effect

#if defined X_NOT_Periodic || defined Y_NOT_Periodic

    string output_filename("LB_Results/parameters");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::app);

	int b_type;
	int low_y, high_y, low_x, high_x;
	double temp_bound;
	int X, Y;
	int bbtype[4];
	bbtype[0] = Left_non_zero_temperature;
	bbtype[1] = Right_non_zero_temperature;
	bbtype[2] = Bot_non_zero_temperature;
	bbtype[3] = Top_non_zero_temperature;

	int start_i = 0;
	int end_i = 3;

	int Xmin = 1;
	int Xmax = Nx - 2;
	int Ymin = 1;
	int Ymax = Ny - 2;

#if defined X_Periodic
	start_i = 2;
	Xmin = 0;
	Xmax = Nx - 1;
#endif
#if defined Y_Periodic
	end_i = 1;
	Ymin = 0;
	Ymax = Ny - 1;
#endif
	std::cout << endl;
	string st;

	for (size_t kk = start_i; kk <= end_i; kk++)	{
		switch (kk)
		{
		case 0: {st = "X=1"; break; }
		case 1: {st = "X=Nx-2"; break; }
		case 2: {st = "Y=1"; break; }
		case 3: {st = "Y=Ny-2"; break; }
		default: st = "X=0";
		}
		std::cout << "Number of Non-zero temperature boundaries at " << st << " is: " << bbtype[kk] << endl;

		for (size_t i = 0; i < bbtype[kk]; i++)
		{
			std::cout << "Please specify the boundary type of ";
			switch (((i + 1) % 10))
			{
			case 1: {std::cout << "1st non-zero temperature boundary at " << st << endl; break; }
			case 2: {std::cout << "2nd non-zero temperature boundary at " << st << endl; break; }
			case 3: {std::cout << "3rd non-zero temperature boundary at " << st << endl; break; }
			default: std::cout << i + 1 << "th non-zero temperature boundary at " << st << endl;
			}
			std::cout << "1) Constant Temperature 2) Adiabatic: ";
			std::cin >> b_type;

			while (b_type > 3 || b_type < 1)
			{
				std::cout << "1) Constant Temperature 2) Adiabatic: ";
				std::cin >> b_type;
			}

			if (kk < 2)
			{
				std::cout << "Please specify lower Y location of the boundary (" << Ymin << "-" << Ymax << "): ";
				std::cin >> low_y;
				std::cout << "Please specify higher Y location of the boundary (" << Ymin << "-" << Ymax << "): ";
				std::cin >> high_y;
			}
			else
			{
				std::cout << "Please specify lower X location of the boundary (" << Xmin << "-" << Xmax << "): ";
				std::cin >> low_x;
				std::cout << "Please specify higher X location of the boundary (" << Xmin << "-" << Xmax << "): ";
				std::cin >> high_x;
			}

			if (b_type == 1)
			{
				std::cout << "Enter temperature value: ";
				std::cin >> temp_bound;
				output_file << "Temperature Boundary No. " << i + 1 << " at " << st << ": Constant Temperature" << endl;
				output_file << "Temperature: " << temp_bound << endl;
			}


			if (b_type == 2)
			{
				output_file << "Boundary No. " << i + 1 << " at " << st << ": Adiabatic" << endl;
			}

			if (kk < 2)
			{
				output_file << "Lower Y value: " << low_y << ", Higher Y value: " << high_y << endl;
			}
			else
			{
				output_file << "Lower X value: " << low_x << ", Higher X value: " << high_x << endl << endl;
			}

			switch (kk)
			{

			case 0: {X = 1; break; }
			case 1: {X = Nx - 2; break; }//Nx-2 because there is a buffer layer
			case 2: {Y = 1; break; }
			case 3: {Y = Ny - 2; break; }
			default: { cout << "uizu"; }
			}

			if ((kk == 0) || (kk == 1))
			{
				for (int io = low_y; io <= high_y; io++)
				{

					if (b_type == 1)
					{
						boundary_cond_temp[X][io] = 1;
						boundary_cond_temp_2[X][io] = 1;
						boundary_temp[X][io] = temp_bound;
					}
					if (b_type == 2)
					{
						boundary_cond_temp[X][io] = 3;
						boundary_cond_temp_2[X][io] = 3;
					}
				}
			}

			if ((kk == 2) || (kk == 3))
			{
				for (int oi = low_x; oi <= high_x; oi++)
				{

					if (b_type == 1)
					{
						boundary_cond_temp[oi][Y] = 2;
						boundary_cond_temp_2[oi][Y] = 2;
						boundary_temp[oi][Y] = temp_bound;
					}
					if (b_type == 2)
					{
						boundary_cond_temp[oi][Y] = 4;
						boundary_cond_temp_2[oi][Y] = 4;
					}

				}
			}
			std::cout << endl;
		}
	}

#if defined X_NOT_Periodic && defined Y_NOT_Periodic
	int Y1 = 1;
	int Y2 = Ny - 2;

	for (size_t ik = 0; ik <= 1; ik++)
	{
		int XX = Nx - 2;
		if (ik == 1)
		{
			XX = 1;
		}

		if (boundary_cond_temp[XX][Y1] != boundary_cond_temp[XX][Y1 + 1])
		{
			//edge_boundary[XX][Y1] = TRUE;
			boundary_cond_temp_2[XX][Y1] = boundary_cond_temp[XX][Y1 + 1];

			if (boundary_cond_temp[XX][Y1 + 1] == 1)
			{
				boundary_temp_2[XX][Y1] = boundary_temp[XX][Y1 + 1];
			}
		}

		if (boundary_cond_temp[XX][Y2] != boundary_cond_temp[XX][Y2 - 1])
		{
			//edge_boundary[XX][Y2] = TRUE;
			boundary_cond_temp_2[XX][Y2] = boundary_cond_temp[XX][Y2 - 1];

			if (boundary_cond_temp[XX][Y2 - 1] == 1)
			{
				boundary_temp_2[XX][Y2] = boundary_temp[XX][Y2 - 1];
			}
		}
	}
#endif

#if defined Y_NOT_Periodic && defined X_NOT_Periodic
	for (size_t ik = 0; ik <= 1; ik++)
	{
		int YY = Ny - 2;
		if (ik == 1)
		{
			YY = 1;
		}

		if (boundary_cond_temp[1][YY] != boundary_cond_temp[2][YY])
		{
		//	edge_boundary[1][YY] = TRUE;
			boundary_cond_temp_2[1][YY] = boundary_cond_temp[2][YY];

			if (boundary_cond_temp[2][YY] == 2)
			{
				boundary_temp_2[1][YY] = boundary_temp[2][YY];
			}
		}

		if (boundary_cond_temp[Nx - 2][YY] != boundary_cond_temp[Nx - 3][YY])	{
		//	edge_boundary[Nx - 2][YY] = TRUE;
			boundary_cond_temp_2[Nx - 2][YY] = boundary_cond_temp[Nx - 3][YY];

			if (boundary_cond_temp[Nx - 3][YY] == 2)
			{
				boundary_temp_2[Nx - 2][YY] = boundary_temp[Nx - 3][YY];
			}
		}

	}
#endif

	std::cout << endl;

	output_file.close();

#endif
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
Inlet_outlets::~Inlet_outlets()
{
}
