// LB_Solver 2D Ver 104.0 Without-MPI
// This is 2D immersed boundary lattice Boltzmann method code.
// It uses the D2Q9 lattice including forcing term.
// Bottom and top walls are parallel to the x-axis (channel).
// Bottom and top walls can be stationary or moving.
// particles are positioned in the flow.
// All particles must have the same shape.
// This particle can be rigid and stationary/moving.
// The lattice velocities are defined according to the following scheme:
// index:   0  1  2  3  4  5  6  7  8
// ----------------------------------
// x:       0 +1  0 -1  0 +1 -1 -1 +1
// y:       0  0 +1  0 -1 +1 +1 -1 -1
//
// 6 2 5  ^y
//  \|/   |   x
// 3-0-1   --->
//  /|\
// 7 4 8

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************
//#include "stdafx.h"  //This line should be active just when you use Visual-Studio
#define _USE_MATH_DEFINES
#include <vector>   // vector containers
#include <cmath>    // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream>  // file stream
#include <sstream>  // string streams
#include <cstdlib>  // standard library
#include <iomanip>  // For set precision. From ver 28
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

using namespace std; // permanently use the standard namespace

#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"
#include "Flow_solver.h"
#include "Thermal_solver.h"
#include "Fluid_read_write.h"
                     //#include "Solid_node.h"
                     //#include "solid_surface.h"
//#include "Inlet_outlets.h"

Flow_solver Flow_Field;
Thermal_solver Temperature_Field;
Fluid_read_write Flow_read_write;//To read and write fluid data in the file
//Inlet_outlets Inlets_Outlets;

#if defined Flow_With_Resolved_Particle
//#include "node_struct.h"
//#include "center_struct.h"
//#include "particle_struct.h"
#include "Particle_sim.h"
Particle_sim Particle_Model;
#endif

#if defined Flow_With_Point_Particle
#include "Point_Particle.h"
Point_Particle Point_Particle_Model;
#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MAIN
//-----------------------------------------------------------------------------------------------------------------------------
//
// This is the main function, containing the simulation initialization and the simulation loop.

int main(int argc, char** argv) {

    clock_t elapsed_time, elapsed_time1;
    elapsed_time = clock();

    /// ---------------------------------------------------------------------------------------------
    /// ------------------------------ REPORT PARAMETERS ON SCREEN ----------------------------------
    /// ---------------------------------------------------------------------------------------------
    const double nu = (tau - 0.5) / 3; // lattice viscosity

    cout << "Simulation Parameters" << endl;
    cout << "=====================" << endl;
    cout << "Nx                  = " << Nx << endl;
    cout << "Ny                  = " << Ny << endl;
    cout << "nu (viscosity)      = " << nu << endl;
    cout << "Gravity             = " << gravity << endl;

#if defined Flow_With_Resolved_Particle
    cout << "Number of Particles = " << num_particles << endl;
    if (num_particles > 0) {
        cout << "Particle Diameter   = " << particle_radius * 2 << endl;
#if defined MOVING_ELLIPSE || defined STATIONARY_ELLIPSE
        cout << "Particle Diameter 2 = " << particle_radius_2 * 2 << endl;
#endif
    }
#endif
    cout << endl;
    /// ---------------------------------------------------------------------------------------------
    /// --------------------------- END OF REPORT PARAMETERS ON SCREEN ------------------------------
    /// ---------------------------------------------------------------------------------------------

    /// **************************
    /// PREPARATIONS
    /// **************************
    //int recovery_step = 0; //If "Recovery_ON" is activated the program starts from step=recovery_step+1
    Flow_Field.initialize(); // allocate memory and initialize variables
    Temperature_Field.initialize_t(); // allocate memory and initialize variables
    Flow_read_write.write_parametes();
    //Inlets_Outlets.set_boundaries();
    //Inlets_Outlets.set_boundaries_thermal();

#if defined Recovery_ON
    Flow_read_write.read_fluid_recovery(recovery_step);
#endif

    /// ********************************
    /// ********************************
    ///         SIMULATION LOOP
    /// ********************************
    /// ********************************

    cout << "STARTING SIMULATION ..." << endl;
    srand(1);

#if defined Flow_With_Resolved_Particle
    Particle_Model.Particle_initialize();

    /// ---------------------------------------------------------------------------------------------
    /// ---------------------- READ PARTICLE RECOVERY DATA FROM FILE  -------------------------------
    /// ---------------------------------------------------------------------------------------------
#if defined Recovery_ON
    Particle_Model.read_particle_recovery(t, particle);
#endif
    /// ---------------------------------------------------------------------------------------------
    /// --------------------- END OF READ PARTICLE RECOVERY DATA FROM FILE  -------------------------
    /// ---------------------------------------------------------------------------------------------

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

    for (int t = 1 + recovery_step; t <= t_num; ++t) { // START OF TIME STEPS (run over all times between 1 and t_num)

        /// ---------------------------------------------------------------------------------------------
        /// ----------------------------- LBM & THERMAL LBM  --------------------------------------------
        /// ---------------------------------------------------------------------------------------------
        Flow_Field.LBM(t, Particle_Model);
        Flow_Field.BC(t);
        Flow_Field.momenta();


#if defined Flow_With_Thermal_Effect
        Temperature_Field.LBM_t(t);
        Temperature_Field.BC_t(t);
        Temperature_Field.momenta_t();
#endif
        /// ---------------------------------------------------------------------------------------------
        /// ----------------------------- END OF LBM & THERMAL LBM  -------------------------------------
        /// ---------------------------------------------------------------------------------------------


        /// ---------------------------------------------------------------------------------------------
        /// ----------------------------- EULERIAN FORCE INITIALIZATION ---------------------------------
        /// ---------------------------------------------------------------------------------------------
        /// Important: Forces at Eulerian points must be set to ZERO at each iteration before SPREAD. We
        /// do this at the end of MOMENTA subroutine.
        /// This is done so just to reduce computational cost & time. (We can set all forces Zero here
        /// but it is time consuming)
        /// ---------------------------------------------------------------------------------------------
        /// ---------------------------------------------------------------------------------------------


        /// ---------------------------------------------------------------------------------------------
        /// ---------------------------------- DIRECT FORCING  ------------------------------------------
        /// ---------------------------------------------------------------------------------------------
        for (int npr2 = 0; npr2 < num_particles; ++npr2) {
            Particle_Model.set_forces_torques_zero(Particle_Model.particle[npr2]);
            Particle_Model.interpolate(Particle_Model.particle[npr2], 0); // interpolate velocity from the Eulerian mesh to Lagrangian node

    #if defined MOVING_CYLINDER
            Particle_Model.collision_cir(Particle_Model.particle, npr2);
    #endif
    #if defined MOVING_ELLIPSE
            Particle_Model.collision_ellipse(Particle_Model.particle, npr2);
    #endif

    #if defined SDF
            Particle_Model.DF_force_calc(Particle_Model.particle[npr2], npr2);//Force and heat calculation at each node
            Particle_Model.spread(Particle_Model.particle[npr2]); // Spread forces from the Lagrangian to the Eulerian mesh
    #endif // defined SDF

    #if defined MDF
            int start_x_force_zero = MAX(0, Particle_Model.particle[npr2].center[0].x      - MAX(particle_radius, particle_radius_2) - 4);//instead of start_X
            int end_x_force_zero   = MIN(Nx - 1, Particle_Model.particle[npr2].center[0].x + MAX(particle_radius, particle_radius_2) + 4);//instead of end_X
            int start_y_force_zero = MAX(0, Particle_Model.particle[npr2].center[0].y      - MAX(particle_radius, particle_radius_2) - 4);//instead of start_Y
            int end_y_force_zero   = MIN(Ny - 1, Particle_Model.particle[npr2].center[0].y + MAX(particle_radius, particle_radius_2) + 4);//instead of end_Y

            for (int df_iter = 1; df_iter <= mdfstep; ++df_iter){//mdf iterations. mdfstep = 1 means SDF

                Particle_Model.DF_force_calc(Particle_Model.particle[npr2], npr2);//Force and heat calculation at each node

                Particle_Model.spread(Particle_Model.particle[npr2]); // Spread forces from the Lagrangian to the Eulerian mesh

                for (int X = start_x_force_zero; X <= end_x_force_zero; ++X) {
                    for (int Y = start_y_force_zero; Y <= end_y_force_zero; ++Y) {
                        mdf_vel_x[X][Y] = velocity_x[X][Y] + force_x[X][Y] * Dt_LB;// (See Luo et al. PRE 066709 (2007) http://dx.doi.org/10.1103/PhysRevE.76.066709)
                        mdf_vel_y[X][Y] = velocity_y[X][Y] + force_y[X][Y] * Dt_LB;
                        #if defined Flow_With_Thermal_Effect
                            mdf_temperature[X][Y] = temperature[X][Y] + force_thermal[X][Y] * Dt_LB;
                        #endif

                    }
                }
                Particle_Model.interpolate(Particle_Model.particle[npr2], df_iter); // interpolate velocity from the Eulerian mesh to Lagrangian node
            }
    #endif // defined MDF

            Particle_Model.total_force_torque_calc(Particle_Model.particle[npr2], npr2);//Total force and heat calculation exerted on the "center"
        }
        /// ---------------------------------------------------------------------------------------------
        /// ------------------------------- END OF DIRECT FORCING   -------------------------------------
        /// ---------------------------------------------------------------------------------------------

        /// ---------------------------------------------------------------------------------------------
        /// --------------------------------- GUO MODEL  ------------------------------------------------
        /// ---------------------------------------------------------------------------------------------
#if defined Guo
        for (int X = start_X; X <= end_X; ++X) {
            for (int Y = start_Y; Y <= end_Y; ++Y) {
                int gbetaindex = 1;
                #if defined Flow_With_Resolved_Particle
                gbetaindex = Flow_Field.calculate_gbetaindex(X, Y, Particle_Model);
                #endif
                double forceX = force_x[X][Y] + gbetaindex * gbeta_x * (temperature[X][Y] - Tref) + gravity;
                double forceY = force_y[X][Y] + gbetaindex * gbeta_y * (temperature[X][Y] - Tref);
                velocity_x[X][Y] += forceX / 2.0;
                velocity_y[X][Y] += forceY / 2.0;
            }
        }
#endif
        /// ---------------------------------------------------------------------------------------------
        /// --------------------------------- END OF GUO MODEL  -----------------------------------------
        /// ---------------------------------------------------------------------------------------------

#if defined Force_Thermal_2
        for (int X = start_X; X <= end_X; ++X) {
            for (int Y = start_Y; Y <= end_Y; ++Y) {
                temperature[X][Y] = temperature[X][Y] + force_thermal[X][Y] / 2.0;
            }
        }
#endif

        Particle_Model.update_location(Particle_Model.particle, t); //Update particle location and temperature

        if (t % t_data == 0) {
            Particle_Model.write_particle_data(t, Particle_Model.particle);//This file is used for post processing
            Flow_read_write.write_fluid_data(t);//You can write fluid data into a .dat file by activating this command
        }

        if (t % t_recovery == 0) {
            Flow_read_write.write_fluid_recovery(t); //Recommended: You can write fluid data into a "dat" file. In case of sudden program failure just activate "Recovery_ON" and run it so the program continues from a almost last iteration
            Particle_Model.write_particle_recovery(t, Particle_Model.particle);
        }


        /// ***************************
        ///  WRITE VTK DATA
        /// ***************************
        if (t % t_vtk == 0) {
            Flow_read_write.write_fluid_vtk(t);
            Particle_Model.write_particle_vtk(t, Particle_Model.particle);
        }

        /// ***************************
        ///  REPORT TIME ON SCREEN
        /// ***************************
        if (t % t_time == 0) {
            elapsed_time1 = clock() - elapsed_time;
            printf("The program took %f seconds from beginning .\n", ((float)elapsed_time1) / CLOCKS_PER_SEC);

        }

        /// ***************************
        ///  REPORT ON COMMAND SCREEN
        /// ***************************
        /// Report end of time step
        if (t % t_info == 0) {
            cout << "Completed Time Step " << t << " in [1, " << t_num << "]" << endl;
            printf("===========================================================================\n");
        }

        /// ***************************
        ///  REPORT RESIDUAL ON COMMAND SCREEN
        /// ***************************
        /// Stop program by residual
#if defined Stop_by_Residual
        if (t % t_residual == 0) {
            bool crit_flow, crit_thermal;
            crit_flow = Flow_Field.check_residual();

#if defined Flow_Without_Thermal_Effect
            if (crit_flow == true)
            {
                cout << "The program is finished by RESIDUAL condition. " << endl;
                printf("===========================================================================\n");
                int num_r;
                scanf("%i", &num_r);
                exit(0);
            }
        }
#endif
#if defined Flow_With_Thermal_Effect
        crit_thermal = Temperature_Field.check_residual();

        if (crit_flow == true && crit_thermal == true)
        {
            cout << "The program is finished by RESIDUAL condition. " << endl;
            printf("===========================================================================\n");
            int num_r;
            scanf("%i", &num_r);
            exit(0);
        }
    }
#endif
#endif

}//END OF TIME STEP LOOP

 /// Report successful end of simulation
cout << "SIMULATION COMPLETED. THANKS" << endl;
elapsed_time = clock() - elapsed_time;
printf("THE PROGRAM TOOK %d CLICKS (%f seconds).\n", elapsed_time, ((float)elapsed_time) / CLOCKS_PER_SEC);
std::cin.clear(); // reset any error flags
std::cin.ignore(32767, '\n'); // ignore any characters in the input buffer until we find an enter character
std::cin.get(); // get one more char from the user#endif // End of: Direct forcing and Flow without particle
#endif // END of Direct Forcing && define Flow_With_Resolved_Particle

//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================
//=============================================================================================================================

#if defined Flow_Without_Particle

#if defined Read_Geometry_From_File
                // Read geometry from the file (geometry.txt). This file must be created by MATLAB program: image_to_matrix.m
                // Before using program specify Nx and Ny at the beginning of program according to geometry size.
int sol;
string output_filename("geometry");
output_filename += ".txt";
ifstream output_file;
/// Open file
output_file.open(output_filename.c_str());
/// Read geometry data
for (int X = 0; X<Nx; X++) {
    for (int Y = 0; Y<Ny; Y++) {
        if (!output_file.eof()) {
            output_file >> sol;
            if (sol == 0) {
                is_solid[X][Y] = TRUE;
            }
        }
    }
}
/// Close geometry file
output_file.close();
#endif

//// In below lines specify which nodes you want to be solid. Hint: You can define the whole domain by periodic boundary condition and then define the walls by is_solid
//// This is_solid approach can be just used with periodic or wall boundary condition.
//// You can use below lines in addition to read_geometry_from_file or deactive below lines.
double cir_center_x = 20;
double cir_center_y = 20;
double cir_radius = 6;

//for (int X = 0; X < Nx; ++X) {//Define Solid nodes
//	for (int Y = 0; Y < Ny; ++Y) {
//		if ((SQ(X - cir_center_x) + SQ(Y - cir_center_y)) <= SQ(cir_radius))	{
//			is_solid[X][Y] = TRUE;
//		}
//	}
//}


for (int n = 0; n < surface_num_nodes; ++n) {// Location of nodes on the surface of solid (Used for curved boundary treatment)
    Flow_Field.solid.node_surf[n].x = cir_center_x + cir_radius * sin(2. * M_PI * (double)n / Flow_Field.solid.num_surf_nodes);
    Flow_Field.solid.node_surf[n].y = cir_center_y + cir_radius * cos(2. * M_PI * (double)n / Flow_Field.solid.num_surf_nodes);
}


// Counting total solid nodes.
// Do NOT change below lines:
int total_solid = 0;
for (int X = 0; X < Nx; ++X) {//Solid nodes. Do NOT change
    for (int Y = 0; Y < Ny; ++Y) {
        if (is_solid[X][Y] == TRUE)
        {
            total_solid = total_solid + 1;
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

#if defined Solid_Constant_Temperature//Do not change
            solid_thermal[X][Y] = 1;//Solid nodes are considered as solid with constant temperature
#endif
#if defined Solid_Adiabatic//Do not change
            solid_thermal[X][Y] = 2;//Solid nodes are considered as solid with adiabatic condition
#endif
#if defined Solid_Conjugate//Do not change
            solid_thermal[X][Y] = 3;//Solid nodes are considered as solid with conjuagte heat transfer
#endif

        }
    }
}

int start_X, end_X;
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

for (int t = 1 + recovery_step; t <= t_num; ++t) { // START OF TIME STEPS (run over all times between 1 and t_num)

                                                   /// ******************************
                                                   ///  LBM ALGORITHM
                                                   /// ******************************
    Flow_Field.LBM(t);
    Flow_Field.BC(t);
    /// ******************************
    ///  LBM ALGORITHM FOR TEMPERATURE
    /// ******************************
#if defined Flow_With_Thermal_Effect
    Temperature_Field.LBM_t(t);
    Temperature_Field.BC_t(t);
#endif

    //****************
    //    GUO MODEL
    //****************
#if defined Guo
    for (int X = start_X; X <= end_X; ++X) {
        for (int Y = start_Y; Y <= end_Y; ++Y) {
            if (is_solid[X][Y] == FALSE) {
                velocity_x[X][Y] += (gravity) / 2.0;
                velocity_y[X][Y] += gbeta*(temperature[X][Y] - Tref) / 2.0;// Natural convection only acts in X-direction
            }
        }
    }
#endif

#if defined Force_Thermal_2
    for (int X = start_X; X <= end_X; ++X) {
        for (int Y = start_Y; Y <= end_Y; ++Y) {
            temperature[X][Y] = temperature[X][Y] + force_thermal[X][Y] / 2.0;
        }
    }
#endif


    /// ****************************************
    ///  WRITE FLUID AND PARTICLES DATA TO DISC
    /// ****************************************
    if (t % t_data == 0) {
        Flow_read_write.write_fluid_data(t);//You can write fluid data into a .dat file by activating this command
    }


    /// ***************************
    ///  WRITE RECOVERY DATA
    /// ***************************
    if (t % t_recovery == 0) {
        Flow_read_write.write_fluid_recovery(t);
    }
    /// ***************************
    ///  WRITE VTK DATA
    /// ***************************
    if (t % t_vtk == 0) {
        Flow_read_write.write_fluid_vtk(t);
    }
    if (t == t_vtk) {
        Flow_read_write.write_solid_vtk(t, total_solid);
    }

    /// ***************************
    ///  REPORT TIME ON SCREEN
    /// ***************************
    if (t % t_time == 0) {
        elapsed_time1 = clock() - elapsed_time;
        printf("The program took %f seconds from beginning .\n", ((float)elapsed_time1) / CLOCKS_PER_SEC);

    }

    /// ***************************
    ///  REPORT ON COMMAND SCREEN
    /// ***************************
    /// Report end of time step
    if (t % t_info == 0) {
        cout << "Completed Time Step " << t << " in [1, " << t_num << "]" << endl;
        printf("===========================================================================\n");
    }

    /// ***************************
    ///  REPORT RESIDUAL ON COMMAND SCREEN
    /// ***************************
    /// Stop program by residual
#if defined Stop_by_Residual
    if (t % t_residual == 0) {
        bool crit_flow, crit_thermal;
        crit_flow = Flow_Field.check_residual();

#if defined Flow_Without_Thermal_Effect
        if (crit_flow == true)
        {
            cout << "The program is finished by RESIDUAL condition. " << endl;
            printf("===========================================================================\n");
            int num_r;
            scanf("%i", &num_r);
            exit(0);
        }
    }
#endif
#if defined Flow_With_Thermal_Effect
    crit_thermal = Temperature_Field.check_residual();

    if (crit_flow == true && crit_thermal == true)
    {
        cout << "The program is finished by RESIDUAL condition. " << endl;
        printf("===========================================================================\n");
        int num_r;
        scanf("%i", &num_r);
        exit(0);
    }
}
#endif
#endif

    }//END OF TIME STEP LOOP

     /// Report successful end of simulation
    cout << "SIMULATION COMPLETED. THANKS" << endl;
    elapsed_time = clock() - elapsed_time;
    printf("THE PROGRAM TOOK %d CLICKS (%f seconds).\n", elapsed_time, ((float)elapsed_time) / CLOCKS_PER_SEC);

    std::cin.clear(); // reset any error flags
    std::cin.ignore(32767, '\n'); // ignore any characters in the input buffer until we find an enter character
    std::cin.get(); // get one more char from the user#endif // End of: Direct forcing and Flow without particle
#endif


#if defined Flow_With_Point_Particle

#if defined Read_Geometry_From_File
                    // Read geometry from the file (geometry.txt). This file must be created by MATLAB program: image_to_matrix.m
                    // Before using program specify Nx and Ny at the beginning of program according to geometry size.
    int sol;
    string output_filename("geometry");
    output_filename += ".txt";
    ifstream output_file;
    /// Open file
    output_file.open(output_filename.c_str());
    /// Read geometry data
    for (int X = 0; X<Nx; X++) {
        for (int Y = 0; Y<Ny; Y++) {
            if (!output_file.eof()) {
                output_file >> sol;
                if (sol == 0) {
                    is_solid[X][Y] = TRUE;
                }
            }
        }
    }
    /// Close geometry file
    output_file.close();
#endif

    //// In below lines specify which nodes you want to be solid. Hint: You can define the whole domain by periodic boundary condition and then define the walls by is_solid
    //// This is_solid approach can be just used with periodic or wall boundary condition.
    //// You can use below lines in addition to read_geometry_from_file or deactive below lines.
    double cir_center_x = 20;
    double cir_center_y = 20;
    double cir_radius = 6;

    //for (int X = 0; X < Nx; ++X) {//Define Solid nodes
    //	for (int Y = 0; Y < Ny; ++Y) {
    //		if ((SQ(X - cir_center_x) + SQ(Y - cir_center_y)) <= SQ(cir_radius))	{
    //			is_solid[X][Y] = TRUE;
    //		}
    //	}
    //}


    for (int n = 0; n < surface_num_nodes; ++n) {// Location of nodes on the surface of solid (Used for curved boundary treatment)
        Flow_Field.solid.node_surf[n].x = cir_center_x + cir_radius * sin(2. * M_PI * (double)n / Flow_Field.solid.num_surf_nodes);
        Flow_Field.solid.node_surf[n].y = cir_center_y + cir_radius * cos(2. * M_PI * (double)n / Flow_Field.solid.num_surf_nodes);
    }


    // Counting total solid nodes.
    // Do NOT change below lines:
    int total_solid = 0;
    for (int X = 0; X < Nx; ++X) {//Solid nodes. Do NOT change
        for (int Y = 0; Y < Ny; ++Y) {
            if (is_solid[X][Y] == TRUE)
            {
                total_solid = total_solid + 1;
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

#if defined Solid_Constant_Temperature//Do not change
                solid_thermal[X][Y] = 1;//Solid nodes are considered as solid with constant temperature
#endif
#if defined Solid_Adiabatic//Do not change
                solid_thermal[X][Y] = 2;//Solid nodes are considered as solid with adiabatic condition
#endif
#if defined Solid_Conjugate//Do not change
                solid_thermal[X][Y] = 3;//Solid nodes are considered as solid with conjuagte heat transfer
#endif

            }
        }
    }

    Point_Particle_Model.Particle_initialize();

    /// ---------------------------------------------------------------------------------------------
    /// ---------------------- READ PARTICLE RECOVERY DATA FROM FILE  -------------------------------
    /// ---------------------------------------------------------------------------------------------
//#if defined Recovery_ON
//	Particle_Model.read_particle_recovery(t, particle);
//#endif

    int start_X, end_X;
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
    Point_Particle_Model.write_particle_vtk(0, Point_Particle_Model.particle);

    for (int t = 1 + recovery_step; t <= t_num; ++t) { // START OF TIME STEPS (run over all times between 1 and t_num)
         /// ******************************
         ///  LBM ALGORITHM
         /// ******************************
        Flow_Field.LBM(t);
        Flow_Field.BC(t);
        /// ******************************
        ///  LBM ALGORITHM FOR TEMPERATURE
        /// ******************************
//#if defined Flow_With_Thermal_Effect
//		Temperature_Field.LBM_t(t);
//		Temperature_Field.BC_t(t);
//#endif


        for (int npr2 = 0; npr2 < num_point_particles; ++npr2) {

                Point_Particle_Model.particle[npr2].force_x = 0;
                Point_Particle_Model.particle[npr2].force_y = 0;
                //Point_Particle_Model.particle[npr2].node[nn].force_thermal = 0.0;

                Point_Particle_Model.interpolate(Point_Particle_Model.particle[npr2]); // interpolate velocity from the Eulerian mesh to Lagrangian node

                //Point_Particle_Model.particle[npr2].fx_surf = 0.0;
                //Point_Particle_Model.particle[npr2].fy_surf = 0.0;
                //Point_Particle_Model.particle[npr2].center[0].heat_surf = 0.0;
                //Point_Particle_Model.particle[npr2].center[0].f3x = 0.0; //Total Force due to particle-particle collision
                //Point_Particle_Model.particle[npr2].center[0].f3y = 0.0; //Total Force due to particle-particle collision
                //Point_Particle_Model.particle[npr2].center[0].f4x = 0.0; //Forces due to particle-wall collision
                //Point_Particle_Model.particle[npr2].center[0].f4y = 0.0; //Forces due to particle-wall collision
                //Point_Particle_Model.particle[npr2].center[0].t1_z = 0.0;//Torque due to forces in X direction
                //Point_Particle_Model.particle[npr2].center[0].t2_z = 0.0;//Torque due to forces in Y direction
                //Point_Particle_Model.particle[npr2].center[0].t3_z = 0.0;//Torque due to COLLISION-forces. Just for ELLIPTICAL particles

//#if defined MOVING_CYLINDER
//			Particle_Model.collision_cir(Particle_Model.particle, npr2);
//#endif
//#if defined MOVING_ELLIPSE
//			Particle_Model.collision_ellipse(Particle_Model.particle, npr2);
//#endif
                double mm;
                Point_Particle_Model.force_calc(Point_Particle_Model.particle[npr2]);//Force and heat calculation at each node
                //cout << "salam" << Point_Particle_Model.particle[0].force_x << endl;

            //Particle_Model.spread(Particle_Model.particle[npr2]); // Spread forces from the Lagrangian to the Eulerian mesh
        }

        //****************
        //    GUO MODEL
        //****************
#if defined Guo
        for (int X = start_X; X <= end_X; ++X) {
            for (int Y = start_Y; Y <= end_Y; ++Y) {
                if (is_solid[X][Y] == FALSE) {
                    velocity_x[X][Y] += (gravity) / 2.0;
                    velocity_y[X][Y] += gbeta*(temperature[X][Y] - Tref) / 2.0;// Natural convection only acts in X-direction
                }
            }
        }
#endif

//#if defined Force_Thermal_2
//		for (int X = start_X; X <= end_X; ++X) {
//			for (int Y = start_Y; Y <= end_Y; ++Y) {
//				temperature[X][Y] = temperature[X][Y] + force_thermal[X][Y] / 2.0;
//			}
//		}
//#endif
        //for (int npr2 = 0; npr2 < num_point_particles; ++npr2) {

        //	Point_Particle_Model.particle[npr2].x = Point_Particle_Model.particle[npr2].x + Point_Particle_Model.particle[npr2].vel_fluid_x;
        //	Point_Particle_Model.particle[npr2].y = Point_Particle_Model.particle[npr2].y + Point_Particle_Model.particle[npr2].vel_fluid_y;
        //}
        Point_Particle_Model.update_location(Point_Particle_Model.particle, t); //Update particle location and temperature


        //cout << t << " ssss " << Point_Particle_Model.particle[0].x << endl;
        /// ****************************************
        ///  WRITE FLUID AND PARTICLES DATA TO DISC
        /// ****************************************
        if (t % t_data == 0) {
            Flow_read_write.write_fluid_data(t);//You can write fluid data into a .dat file by activating this command
            Point_Particle_Model.write_particle_data(t, Point_Particle_Model.particle);
        }


        /// ***************************
        ///  WRITE RECOVERY DATA
        /// ***************************
        if (t % t_recovery == 0) {
            Flow_read_write.write_fluid_recovery(t);
        }
        /// ***************************
        ///  WRITE VTK DATA
        /// ***************************
        if (t % t_vtk == 0) {
            Flow_read_write.write_fluid_vtk(t);
            Point_Particle_Model.write_particle_vtk(t, Point_Particle_Model.particle);

        }
        if (t == t_vtk) {
            Flow_read_write.write_solid_vtk(t, total_solid);
        }

        /// ***************************
        ///  REPORT TIME ON SCREEN
        /// ***************************
        if (t % t_time == 0) {
            elapsed_time1 = clock() - elapsed_time;
            printf("The program took %f seconds from beginning .\n", ((float)elapsed_time1) / CLOCKS_PER_SEC);

        }

        /// ***************************
        ///  REPORT ON COMMAND SCREEN
        /// ***************************
        /// Report end of time step
        if (t % t_info == 0) {
            cout << "Completed Time Step " << t << " in [1, " << t_num << "]" << endl;
            printf("===========================================================================\n");
        }

        /// ***************************
        ///  REPORT RESIDUAL ON COMMAND SCREEN
        /// ***************************
        /// Stop program by residual
#if defined Stop_by_Residual
        if (t % t_residual == 0) {
            bool crit_flow, crit_thermal;
            crit_flow = Flow_Field.check_residual();

#if defined Flow_Without_Thermal_Effect
            if (crit_flow == true)
            {
                cout << "The program is finished by RESIDUAL condition. " << endl;
                printf("===========================================================================\n");
                int num_r;
                scanf("%i", &num_r);
                exit(0);
            }
        }
#endif
#if defined Flow_With_Thermal_Effect
        crit_thermal = Temperature_Field.check_residual();

        if (crit_flow == true && crit_thermal == true)
        {
            cout << "The program is finished by RESIDUAL condition. " << endl;
            printf("===========================================================================\n");
            int num_r;
            scanf("%i", &num_r);
            exit(0);
        }
    }
#endif
#endif

    }//END OF TIME STEP LOOP

     /// Report successful end of simulation
    cout << "SIMULATION COMPLETED. THANKS" << endl;
    elapsed_time = clock() - elapsed_time;
    printf("THE PROGRAM TOOK %d CLICKS (%f seconds).\n", elapsed_time, ((float)elapsed_time) / CLOCKS_PER_SEC);

    std::cin.clear(); // reset any error flags
    std::cin.ignore(32767, '\n'); // ignore any characters in the input buffer until we find an enter character
    std::cin.get(); // get one more char from the user#endif // End of: Direct forcing and Flow without particle
#endif
    int num_r;
    scanf("%i", &num_r);
    return 0;
} // end of MAIN function



