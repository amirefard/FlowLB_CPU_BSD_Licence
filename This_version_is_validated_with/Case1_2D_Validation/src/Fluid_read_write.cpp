//#include "stdafx.h"
#include "Fluid_read_write.h"


using namespace std; // permanently use the standard namespace

Fluid_read_write::Fluid_read_write()
{
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SWAP_ENDIAN
//-----------------------------------------------------------------------------------------------------------------------------
void Fluid_read_write::swap_endian(char * buffer, size_t size) {
	for (size_t b = 0; b < size / 2; b++) {
		char temp = buffer[b];
		buffer[b] = buffer[size - 1 - b];
		buffer[size - 1 - b] = temp;
	}
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SWAP_ENDIAN
//-----------------------------------------------------------------------------------------------------------------------------
void Fluid_read_write::swap_endian(char * buffer, char * buffer2, size_t size) {
	for (size_t b = 0; b < size; b++) {
		buffer2[b] = buffer[size - 1 - b];
	}

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_BIGENDIAN
//-----------------------------------------------------------------------------------------------------------------------------

void Fluid_read_write::write_bigendian(std::ofstream & file, char * buffer, size_t count, size_t size) {
	char * _buf = new char[size];
	for (size_t i = 0; i < count; i++) {
		swap_endian(buffer + (size*i), _buf, size);
		file.write(_buf, size);
	}
	delete[] _buf;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_BIGENDIAN
//-----------------------------------------------------------------------------------------------------------------------------
template <typename T>
void Fluid_read_write::write_bigendian(std::ofstream & file, T * buffer, size_t count) {
	char * _buf = new char[sizeof(T)];
	for (size_t i = 0; i < count; i++) {
		swap_endian((char*)(buffer + i), _buf, sizeof(T));
		file.write(_buf, sizeof(T));
	}
	delete[] _buf;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// The fluid state is writen to a VTK file at each t_vtk step.
// The following data is written:
// - density difference (density - 1)
// - x-component of velocity
// - y-component of velocity
// The following code is designed in such a way that the file can be read by ParaView.
void Fluid_read_write::write_fluid_vtk(int time) {

#if defined VTK_BINARY
	/// Create filename

	stringstream output_filename;
    output_filename << "LB_Results/vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "BINARY\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << Nx << " " << Ny << " 1" << "\n";
	output_file << "X_COORDINATES " << Nx << " float\n";
	output_file.close();

	//{
	//	std::stringstream str_line1;
	//	str_line1 << "# vtk DataFile Version 3.0\n";
	//	std::string strstr_line1 = str_line1.str();
	//	output_file1.write(strstr_line1.c_str(), sizeof(char)*strstr_line1.size());

	//output_file1.write((char*)z1, sizeof(float)* 1);
	//}

	float XX[Nx];
	float YY[Ny];
	float ZZ[1];

	for (int X = 0; X < Nx; ++X) {
		XX[X] = X;
	}
	for (int Y = 0; Y < Ny; ++Y) {
		YY[Y] = Y;
	}
	for (int Z = 0; Z < 1; ++Z) {
		ZZ[Z] = Z;
	}

	/// Open file for append
	std::stringstream str_line1;
    str_line1 << "LB_Results/vtk_fluid/fluid_t" << time << ".vtk";
	std::string strstr_line1 = str_line1.str();

	std::ofstream output_file1(strstr_line1, std::ios::app | std::ios::binary);

//	std::ofstream output_file1("LB_Results/vtk_fluid/fluid_t200.vtk", std::ios::app | std::ios::binary);

	write_bigendian(output_file1, XX, Nx);

	output_file1 << "\nY_COORDINATES " << Ny << " float\n";

	write_bigendian(output_file1, YY, Ny);

	output_file1 << "\nZ_COORDINATES " << 1 << " float\n";

	write_bigendian(output_file1, ZZ, 1);

	output_file1 << "\nPOINT_DATA " << Nx * Ny << "\n";

	/// Write density difference
	output_file1 << "SCALARS density_difference float 1\n";
	output_file1 << "LOOKUP_TABLE default\n";

	float mm[1];
	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			mm[0] = density[X][Y]-1;
			write_bigendian(output_file1, mm, 1);
		}
	}

	output_file1 << "\nSCALARS temperature float 1\n";
	output_file1 << "LOOKUP_TABLE default\n";
	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			mm[0] = temperature[X][Y];
			write_bigendian(output_file1, mm, 1);
		}
	}

	float nn[2];
	///// Write velocity
	output_file1 << "VECTORS velocity_vector float\n";
	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			mm[0] = velocity_x[X][Y];
			write_bigendian(output_file1, mm, 1);
			mm[0] = velocity_y[X][Y];
			write_bigendian(output_file1, mm, 1);
			mm[0] = 0;
			write_bigendian(output_file1, mm, 1);
		}
	}

	output_file1.close();
#endif


#if defined VTK_ASCII
	/// Create filename

	stringstream output_filename;
    output_filename << "LB_Results/vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << Nx << " " << Ny << " 1" << "\n";
	output_file << "X_COORDINATES " << Nx << " float\n";

	for (int X = 0; X < Nx; ++X) {
		output_file << X << " ";
	}

	output_file << "\n";
	output_file << "Y_COORDINATES " << Ny << " float\n";
	for (int Y = 0; Y < Ny; ++Y) {
		output_file << Y << " ";
	}

	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << Nx * Ny << "\n";

	///// Write density difference
	output_file << "SCALARS density_difference float 1\n";
	output_file << "LOOKUP_TABLE default\n";

	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			output_file << density[X][Y] - 1 << "\n";
		}
	}

	/// Write temperature
	output_file << "SCALARS temperature float 1\n";
	output_file << "LOOKUP_TABLE default\n";

	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			output_file << temperature[X][Y] << "\n";
		}
	}

	/// Write velocity
	output_file << "VECTORS velocity_vector float\n";

	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			//#if defined Shan_Chen
			//			output_file << velocity_x[X][Y] + 0.5 * (force_x[X][Y] + gravity) / density[X][Y] << " " << velocity_y[X][Y] + 0.5 * (force_y[X][Y] + gbeta*(temperature[X][Y] - Tref)) / density[X][Y] << " 0\n";
			//#endif
			//#if (defined Guo || defined Luo)
			output_file << velocity_x[X][Y] << " " << velocity_y[X][Y] << " 0\n";
			//#endif
		}
	}

	/// Close file
	output_file.close();
#endif

#if defined CASE
	unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);// number are shown as 2.0e02 instead of 2.0e002 (necessay for geo and data files)


	/////////// Writing CASE file (.case) ///////////////////
	/////////////////////////////////////////////////////////
	/// Create filename

	stringstream output_filename;
    output_filename << "LB_Results/vtk_fluid/fluid" << ".case";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# EnSight Gold Model\n";
    output_file << "# Created by LB\n\n";
	output_file << "FORMAT\n";
	output_file << "type:                   ensight gold\n\n";
	output_file << "GEOMETRY\n";
	output_file << "model:                  fluid.geo\n\n";
	output_file << "VARIABLE\n";
	output_file << "scalar per node:           1 density_difference fluid.density_difference******\n";
	output_file << "scalar per node:           1 temperature fluid.temperature******\n";
	output_file << "vector per node:           1 VELOC_VECTOR fluid.VELOC_VECTOR******\n\n";
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

	if (time == t_vtk){//Because .geo file do not change during iterations

		/////////// Writing GEOMETRY file (.geo) ///////////////////
		////////////////////////////////////////////////////////////
		stringstream output_filename1;
        output_filename1 << "LB_Results/vtk_fluid/fluid" << ".geo";//Geometry file
		ofstream output_file1;

		/// Open file
		output_file1.open(output_filename1.str().c_str());

		/// Write VTK header
		output_file1 << "EnSight Model Geometry File\n";
		output_file1 << "EnSight 10.1.1\n";
		output_file1 << "node id given\n";
		output_file1 << "element id given\n";
		output_file1 << "extents\n";
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << (double)Nx - 1. << "\n";
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << (double)Ny - 1. << "\n";
		output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0. << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << "\n";
		output_file1 << "part\n";
		output_file1 << std::setw(10) << 1 << "\n";
		output_file1 << "FLUID\n";
		output_file1 << "coordinates\n";
		output_file1 << std::setw(10) << Nx*Ny << "\n";
		for (int itm = 1; itm <= Nx*Ny; itm++){
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int Y = 0; Y < Ny; Y++){
			for (int X = 0; X < Nx; X++){
				output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << (double)X << "\n";
			}
		}
		for (int Y = 0; Y < Ny; Y++){
			for (int X = 0; X < Nx; X++){
				output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << (double)Y << "\n";
			}
		}
		for (int Y = 0; Y < Ny; Y++){
			for (int X = 0; X < Nx; X++){
				output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << "\n";
			}
		}
		output_file1 << "quad4\n";
		output_file1 << std::setw(10) << (Nx - 1)*(Ny - 1) << "\n";
		for (int itm = 1; itm <= (Nx - 1)*(Ny - 1); itm++){
			output_file1 << std::setw(10) << itm << "\n";
		}

		for (int J = 1; J < Ny; J++){
			for (int I = 1; I < Nx; I++){
				output_file1 << std::setw(10) << I + (J - 1)*Nx << std::setw(10) << I + 1 + (J - 1)*Nx << std::setw(10) << (I + 1) + (J*Nx) << std::setw(10) << (I)+(J*Nx) << "\n";

			}
		}
		/// Close file
		output_file1.close();
		/////////// End of Writing GEOMETRY file (.geo) ///////////////////
		///////////////////////////////////////////////////////////////////

	}// End of if (time == t_vtk)


	///////////// Writing DENSITY_DIFFERENCE file ////////////////////////
	/////////////////////////////////////////////////////////////////////
	stringstream output_filename4;
    output_filename4 << "LB_Results/vtk_fluid/fluid.density_difference" << setfill('0') << setw(6) << time;
	ofstream output_file4;

	/// Open file
	output_file4.open(output_filename4.str().c_str());

	output_file4 << "density_difference\n";
	output_file4 << "part\n";
	output_file4 << std::setw(10) << 1 << "\n";
	output_file4 << "coordinates\n";

	for (int Y = 0; Y < Ny; Y++){
		for (int X = 0; X < Nx; X++){
			output_file4 << std::scientific << std::setw(12) << std::setprecision(5) << density[X][Y] - 1. << "\n";
		}
	}

	/// Close file
	output_file4.close();

	/////////// End of Writing DENSITY_DIFFERENCE file /////////////////
	////////////////////////////////////////////////////////////

	///////////// Writing TEMPERATURE file ////////////////////////
	/////////////////////////////////////////////////////////////////////
	stringstream output_filename6;
    output_filename6 << "LB_Results/vtk_fluid/fluid.temperature" << setfill('0') << setw(6) << time;
	ofstream output_file6;

	/// Open file
	output_file6.open(output_filename6.str().c_str());

	output_file6 << "temperature\n";
	output_file6 << "part\n";
	output_file6 << std::setw(10) << 1 << "\n";
	output_file6 << "coordinates\n";

	for (int Y = 0; Y < Ny; Y++){
		for (int X = 0; X < Nx; X++){
			output_file6 << std::scientific << std::setw(12) << std::setprecision(5) << temperature[X][Y] << "\n";
		}
	}

	/// Close file
	output_file6.close();

	/////////// End of Writing TEMPERATURE file /////////////////
	////////////////////////////////////////////////////////////

	///////////// Writing VELOCITY_VECTOR file ////////////////////////
	/////////////////////////////////////////////////////////////////////
	stringstream output_filename5;
    output_filename5 << "LB_Results/vtk_fluid/fluid.VELOC_VECTOR" << setfill('0') << setw(6) << time;
	ofstream output_file5;

	/// Open file
	output_file5.open(output_filename5.str().c_str());

	output_file5 << "VELOC_VECTOR\n";
	output_file5 << "part\n";
	output_file5 << std::setw(10) << 1 << "\n";
	output_file5 << "coordinates\n";

	for (int Y = 0; Y < Ny; Y++){
		for (int X = 0; X < Nx; X++){
			output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_x[X][Y] << "\n";
		}
	}
	for (int Y = 0; Y < Ny; Y++){
		for (int X = 0; X < Nx; X++){
			output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_y[X][Y] << "\n";
		}
	}

	for (int Y = 0; Y < Ny; Y++){
		for (int X = 0; X < Nx; X++){
			output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << "\n";
		}
	}

	/// Close file
	output_file5.close();

	/////////// End of Writing VELOCITY_VECTOR file /////////////
	////////////////////////////////////////////////////////////


#endif

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_SOLID_VTK
//-----------------------------------------------------------------------------------------------------------------------------
// Location of solid nodes specified by is_solid=TRUE
void Fluid_read_write::write_solid_vtk(int time, int tot_solid) {

	/// Create filename
	stringstream output_filename;
    output_filename << "LB_Results/vtk_particle/solid_nodes_" << 1 << "_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());


	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "solid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << Nx << " " << Ny << " 1" << "\n";
	output_file << "X_COORDINATES " << Nx << " float\n";

	for (int X = 0; X < Nx; ++X) {
		output_file << X << " ";
	}

	output_file << "\n";
	output_file << "Y_COORDINATES " << Ny << " float\n";
	for (int Y = 0; Y < Ny; ++Y) {
		output_file << Y << " ";
	}

	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << Nx * Ny << "\n";

	///// Write density difference
	output_file << "SCALARS solid int 1\n";
	output_file << "LOOKUP_TABLE default\n";

	for (int Y = 0; Y < Ny; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			output_file << is_solid[X][Y] << "\n";
		}
	}


	/// Close file
	output_file.close();


	////////////////// Second file

	/// Create filename
	stringstream output_filename2;
    output_filename2 << "LB_Results/vtk_particle/solid_nodes_" << 2 << "_t" << time << ".vtk";
	ofstream output_file2;

	/// Open file
	output_file2.open(output_filename2.str().c_str());

	/// Write VTK header
	output_file2 << "# vtk DataFile Version 3.0\n";
	output_file2 << "particle_state\n";
	output_file2 << "ASCII\n";
	output_file2 << "DATASET POLYDATA\n";

	 //Write node positions
	output_file2 << "POINTS " << tot_solid << " float\n";
	for (int Y = 1; Y < Ny - 1; ++Y) {
		for (int X = 0; X < Nx; ++X) {
			if (is_solid[X][Y] == TRUE)	{
				output_file2 << X << " " << Y << " 0\n";
			}
		}
	}

	///// Write vertices
	output_file2 << "VERTICES 1 " << tot_solid + 1 << "\n";
	output_file2 << tot_solid << " ";
	for (int n = 0; n < tot_solid; ++n) {
		output_file2 << n << " ";
	}

	/// Close file
	output_file2.close();


	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_DATA
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are written to the disk at each t_data step:
void Fluid_read_write::write_fluid_data(int time) {

	/// Create filename
    string output_filename("LB_Results/data_fluid");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::app);

	/// Write data
	output_file << time << " "; // time step

//	double l2_error = -10000;
//	double error=0;
//	for (size_t X = 1; X < Nx - 1; X++)
//	{
//		for (size_t Y = 1; Y < Ny - 1; Y++)
//		{
//			error = 0;
//			error = SQ(velocity_x[X][Y] - (4.0 * 0.01*((Y - 0.5) / (Ny - 2.))*(1 - (Y - 0.5) / (Ny - 2.)))) + SQ(velocity_y[X][Y] - 0);
//			if (error>l2_error)
//			{
//				l2_error = error;
//			}

//		}
//	}
//	l2_error = sqrt(l2_error) / ((Nx - 1)*(Ny - 1));
//	output_file << l2_error << endl;

	//You can change the location at which you want to save the data

	//	double tgrad;
	//	for (int X = 1; X<Nx - 1; X++) {
	//		tgrad = (-3.0*temperature[X][1] + 4.0*temperature[X][2] - temperature[X][3]) / 2.0;
	//		output_file << tgrad << " ";
	//	}
	//	output_file << "\n";

	//	double tbulk;
	//	for (int X = 1; X<Nx - 1; X++) {
	//		double uu = 0.0;
	//		double uutt = 0.0;
	//		for (int Y = 1; Y < Ny - 1; Y++) {
	//			uu = uu + velocity_x[X][Y];
	//			uutt = uutt + velocity_x[X][Y] * temperature[X][Y];
	//		}

	//		tbulk = uutt / uu;
	//		output_file << tbulk << " ";
	//	}

	output_file << "\n";


	//Permeability calculation and storage
	double permeability;
	if (gravity>0) {
		double uav = 0.;
		for (int X = 0; X < Nx; ++X) {
			for (int Y = 0; Y < Ny; ++Y) {
				if (is_solid[X][Y] == FALSE)	{
					uav = uav + velocity_x[X][Y];
				}
			}
		}
		uav = uav / (Nx * Ny);
		permeability = uav * (tau - 0.5) / 3 / gravity;
		printf("Flow permeability in Lattice units is %f .\n", permeability);
		output_file << permeability << " ";
	}
	//////////////////////////////////////////

	//output_file << velocity_x[10][10] << "\n";

	/// Close file
	output_file.close();

	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are written to the disk at each t_recovery step:
void Fluid_read_write::write_fluid_recovery(int time) {


#if defined BINARY_Recovery_File

	std::stringstream str_line1;
    str_line1 << "LB_Results/data_fluid_recovery.dat";
	std::string strstr_line1 = str_line1.str();

	std::ofstream output_file1(strstr_line1, std::ios::out | std::ios::binary);

	float tt[1];
	tt[0] = time;
	//std::cout << tt[0];
	output_file1.write((char*)tt, sizeof(float)*1);

	float mm[1];

	for (int X = 0; X<Nx; X++){
		for (int Y = 0; Y<Ny; Y++){
			for (int c_i = 0; c_i < 9; ++c_i) {
				mm[0] = pop[c_i][X][Y];
				output_file1.write((char*)mm, sizeof(float)* 1);

			}
			mm[0] = density[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);
			mm[0] = velocity_x[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);
			mm[0] = velocity_y[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);
			mm[0] = force_x[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);
			mm[0] = force_y[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);

#if defined Flow_With_Thermal_Effect
			for (int c_i = 0; c_i < 9; ++c_i) {
				mm[0] = pop_t[c_i][X][Y];
				output_file1.write((char*)mm, sizeof(float)* 1);
			}
			mm[0] = temperature[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);
			mm[0] = force_thermal[X][Y];
			output_file1.write((char*)mm, sizeof(float)* 1);
#endif
		}
	}

	output_file1.close();
#endif


#if defined ASCII_Recovery_File
	/// Create filename
    string output_filename("LB_Results/data_fluid_recovery");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

	/// Write data
	output_file << time << " "; // time step is written just as the first number in the recovery file.

	for (int X = 0; X<Nx; X++){
		for (int Y = 0; Y<Ny; Y++){
			for (int c_i = 0; c_i < 9; ++c_i) {
				//	output_file << std::fixed << std::setprecision(15) << pop_old[c_i][X][Y] << " ";
				output_file << std::fixed << std::setprecision(15) << pop[c_i][X][Y] << " ";
			}
			output_file << std::fixed << std::setprecision(15) << density[X][Y] << " ";
			output_file << std::fixed << std::setprecision(15) << velocity_x[X][Y] << " ";
			output_file << std::fixed << std::setprecision(15) << velocity_y[X][Y] << " ";
			output_file << std::fixed << std::setprecision(15) << force_x[X][Y] << " ";
			output_file << std::fixed << std::setprecision(15) << force_y[X][Y] << " ";
#if defined Flow_With_Thermal_Effect
			for (int c_i = 0; c_i < 9; ++c_i) {
				output_file << std::fixed << std::setprecision(15) << pop_t[c_i][X][Y] << " ";
			}
			output_file << std::fixed << std::setprecision(15) << temperature[X][Y] << " ";
			output_file << std::fixed << std::setprecision(15) << force_thermal[X][Y] << "\n";
#endif
		}
	}

	/// Close file
	output_file.close();
#endif
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: READ_FLUID_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are READ from the disk.
void Fluid_read_write::read_fluid_recovery(int & Step) {

#if defined BINARY_Recovery_File
	std::stringstream str_line1;
    str_line1 << "LB_Results/data_fluid_recovery.dat";
	std::string strstr_line1 = str_line1.str();

	std::ifstream output_file1(strstr_line1, std::ios::in | std::ios::binary);

	float tt[1];
	output_file1.read((char*)tt,sizeof(float));

	Step=static_cast<int>(tt[0]);
	//std::cout << Step << endl;
	if (output_file1.is_open())	{
		for (int X = 0; X < Nx; X++) {
			for (int Y = 0; Y < Ny; Y++) {
				for (int c_i = 0; c_i < 9; ++c_i) {
					output_file1.read((char*)tt, sizeof(float));
					pop[c_i][X][Y] = tt[0];
				}
				output_file1.read((char*)tt, sizeof(float));
				density[X][Y] = tt[0];
				output_file1.read((char*)tt, sizeof(float));
				velocity_x[X][Y] = tt[0];
				output_file1.read((char*)tt, sizeof(float));
				velocity_y[X][Y] = tt[0];
				output_file1.read((char*)tt, sizeof(float));
				force_x[X][Y] = tt[0];
				output_file1.read((char*)tt, sizeof(float));
				force_y[X][Y] = tt[0];
#if defined Flow_With_Thermal_Effect
				for (int c_i = 0; c_i < 9; ++c_i) {
					output_file1.read((char*)tt, sizeof(float));
					pop_t[c_i][X][Y] = tt[0];
				}
				output_file1.read((char*)tt, sizeof(float));
				temperature[X][Y] = tt[0];
				output_file1.read((char*)tt, sizeof(float));
				force_thermal[X][Y] = tt[0];
#endif
			}
		}

		output_file1.close();
	}

#endif

#if defined ASCII_Recovery_File
	/// Create filename
    string output_filename("LB_Results/data_fluid_recovery");
	output_filename += ".dat";
	ifstream output_file;//File is open for READING

	/// Open file
	output_file.open(output_filename.c_str());


	output_file >> Step;        // get first number from the file. The first number is: "saved recovery step"

	//while ( ! output_file.eof() ) {

	for (int X = 0; X<Nx; X++) {
		for (int Y = 0; Y<Ny; Y++) {
			if (!output_file.eof()) {
				for (int c_i = 0; c_i < 9; ++c_i) {
					output_file >> pop[c_i][X][Y];
				}
				output_file >> density[X][Y];
				output_file >> velocity_x[X][Y];
				output_file >> velocity_y[X][Y];
				output_file >> force_x[X][Y];
				output_file >> force_y[X][Y];
#if defined Flow_With_Thermal_Effect
				for (int c_i = 0; c_i < 9; ++c_i) {
					output_file >> pop_t[c_i][X][Y];
				}
				output_file >> temperature[X][Y];
				output_file >> force_thermal[X][Y];
#endif
			}
		}
	}

	/// Close file
	output_file.close();
#endif
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARAMETES
//-----------------------------------------------------------------------------------------------------------------------------
//
// The Main simulation parameters will be written into a file named data_param.
void Fluid_read_write::write_parametes() {

	/// Create filename
    string output_filename("LB_Results/parameters");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

	/// Write data
#if defined Flow_With_Resolved_Particle
#if defined STATIONARY_CYLINDER
	output_file << "STATIONARY_CYLINDER " << endl;
#endif
#if defined MOVING_CYLINDER
	output_file << "MOVING_CYLINDER " << endl;
#endif
#if defined STATIONARY_ELLIPSE
	output_file << "STATIONARY_ELLIPSE " << endl;
#endif
#if defined MOVING_ELLIPSE
	output_file << "MOVING_ELLIPSE  " << endl;
#endif
#endif
#if defined Shan_Chen
	output_file << "Shan Chen model ";
	output_file << "- Single Relaxation Time (SRT)" << endl;
#endif
#if defined Guo
	output_file << "Guo model";
#endif
#if defined Luo
	output_file << "Luo model";
	output_file << "- Single Relaxation Time (SRT)" << endl;
#endif
#if defined Guo && defined MRT
	output_file << " - Multi Relaxation Time (MRT)" << endl;
#endif
#if defined Guo && defined SRT
	output_file << " - Single Relaxation Time (SRT)" << endl;
#endif
#if defined MDF
	output_file << " Multi Direct Forcing" << endl;
#endif
#if defined X_Periodic
	output_file << "Boundary condition: Periodic in X direction " << endl;
#endif
#if defined X_NOT_Periodic
	output_file << "Boundary condition: NOT Periodic in X direction " << endl;
#endif

#if defined Y_Periodic
	output_file << "Periodic in Y direction " << endl;
#endif
#if defined Y_NOT_Periodic
	output_file << "Boundary condition: NOT Periodic in Y direction " << endl;
#endif


	output_file << "Nx  : " << Nx << endl;
	output_file << "Ny  : " << Ny << endl;
	output_file << "Number of iterations  : " << t_num << endl;
	output_file << "Wrtie to .dat file (iteration): " << t_data << endl;
	output_file << "Wrtie to .vtk file (iteration): " << t_vtk << endl;
	output_file << "Wrtie to recovery file (iteration): " << t_recovery << endl;
	output_file << "Tau : " << tau << endl;
	output_file << "Lattice Viscosity : " << (2 * tau - 1) / 6 << endl;
	output_file << "Gravity : " << gravity << endl;
//#if defined Y_Wall
//	output_file << "Bottom wall (Y=0) velocity : " << wall_vel_bottom << endl;
//	output_file << "Top wall (Y=Ny-1) velocity : " << wall_vel_top << endl;
//#endif

#if defined Flow_With_Resolved_Particle
	output_file << "Number of particles : " << num_particles << endl;

	if (num_particles > 0) {
		output_file << "Particle radius     : " << particle_radius << endl;
#if defined STATIONARY_ELLIPSE || defined MOVING_ELLIPSE
		output_file << "Particle radius 2   : " << particle_radius_2 << endl;
#endif
		output_file << "First particle X position : " << particle_center_x << endl;
		output_file << "First particle Y position : " << particle_center_y << endl;
#if defined Without_added_mass
		output_file << "Without added mass " << endl;
#endif
#if defined With_added_mass
		output_file << "With added mass " << endl;
#endif
		output_file << "Particle number of nodes  : " << particle_num_nodes << endl;
		output_file << "Particle gravity in X dir  : " << particle_gravity_x << endl;
		output_file << "Particle gravity in Y dir  : " << particle_gravity_y << endl;
		output_file << "Particle to Fluid density ratio  : " << den_ratio << endl;
#if defined Two_Point_delta
		output_file << "Two_Point_delta   " << endl;
#endif
#if defined Original_Delta
		output_file << "Original_Delta   " << endl;
#endif
#if defined Four_Point_Delta
		output_file << "Four_Point_Delta   " << endl;
#endif
#if defined Kernel_Delta
		output_file << "Kernel_Delta   " << endl;
#endif
#if  defined MOVING_CYLINDER
#if defined Circle_Collision_1
		output_file << "Collision model for CIRCLULAR particle: Collision 1 " << endl;
#endif
#if defined Circle_Collision_2
		output_file << "Collision model for CIRCLULAR particle: Collision 2 " << endl;
#endif
#if defined Circle_Collision_3
		output_file << "Collision model for CIRCLULAR particle: Collision 3 " << endl;
#endif
#endif
#if defined MOVING_ELLIPSE
#if defined Ellipse_Collision_1
		output_file << "Collision model for ELLIPTICAL particle: Collision 1 " << endl;
#endif
#endif
	}
#endif
	output_file << endl;
	/// Close file
	output_file.close();

	return;
}
//=============================================================================================================================
//=============================================================================================================================

Fluid_read_write::~Fluid_read_write()
{
}
