/// FlowLB_2D_3D.cpp : This file contains the 'main' function.
/// Program originally developed by: Amir E. Fard

/// *********************
/// PREPROCESSOR COMMANDS
/// *********************
#define _USE_MATH_DEFINES
#include <vector>   // vector containers
#include <cmath>    // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream>  // file streams
#include <sstream>  // string streams
#include <cstdlib>  // standard library
#include <iomanip>  // For set precision. From ver 25
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <mpi.h>
using namespace std; // permanently use the standard namespace

//=============================================================================================================================
//=============================================================================================================================
// CLASSES
#include "GlobalFunctions.h"
#include "SimulationSettings.h"
#include "SimulationParameters.h"
#include "ParticleResolved.h"
#include "FlowSolver.h"
#include "ParallelMPI.h"
#include "RunFlowWithParticle.h"
#include "RunFlowWithoutParticle.h"

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: MAIN
//-----------------------------------------------------------------------------------------------------------------------------
//
// This is the main function, containing the simulation initialization and the simulation loop.

int main(int argc, char** argv)
{
    setbuf(stdout, NULL); // Disable buffering for stdout to ensure immediate output in printing the output
///===================================================================================================================================== ///
///                                        DO NOT CHANGE (UNLESS NEEDED)                                                                 ///
///===================================================================================================================================== ///
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;

    std::string address_of_settings_file   = "config/settings.txt";
    std::string address_of_parameters_file = "config/params.txt";

    /// Initialize the MPI environment (parallelLB)
    ParallelMPI parallelLB(argc, argv);//Class of parallel structure
    INT32 this_processor_id = parallelLB.get_processor_id();
    INT32 num_processors_total = parallelLB.get_num_processors_total();

    /// Initialize the Simulation Settings (settingsLB)
    SimulationSettings& settingsLB = SimulationSettings::get_instance(address_of_settings_file);// Get the singleton instance of SimulationSettings with a specific filename
    bool Flow_With_Particle    = settingsLB.get_bool_setting("Flow_With_Particle");

    /// Initialize the Simulation Parameters (paramsLB)
    SimulationParameters& paramsLB = SimulationParameters::get_instance(address_of_parameters_file, settingsLB);// Get the singleton instance of SimulationParameters with a specific filename
    if (this_processor_id == MASTER)
    {
        printf("\nWelcome to 'FlowLB' for Multiphase Flow Simulations. \n\n");
        printf("============= Parameters Report =========== \n");
        paramsLB.report_parameters_on_screen();
    }

    /// Setup the MPI domain decomposition
    parallelLB.prepareMpiSetup(paramsLB, settingsLB);
    //std::vector<bool> is_processor_for_particle_only = parallelLB.get_is_processor_for_particle_only();// Should be after prepareMpiSetup. Before here it is false for all processors

///===================================================================================================================================== ///
///                                                Flow_With_Particle                                                                    ///
///===================================================================================================================================== ///
    if (Flow_With_Particle)
    {
      if (this_processor_id == MASTER) cout << "[Main] Running Flow With Particle Case\n" << endl;
        run_flow_with_particle(paramsLB, settingsLB, parallelLB);
    }// End of: if (Flow_With_Particle)
///===================================================================================================================================== ///
///                                              !Flow_With_Particle                                                                   ///
///===================================================================================================================================== ///
    else if (!Flow_With_Particle)
    {
        if (this_processor_id == MASTER) cout << "[Main] Running Flow Without Particle Case\n" << endl;
        run_flow_without_particle(paramsLB, settingsLB, parallelLB);
    }// End of: if (!Flow_With_Particle)
///===================================================================================================================================== ///

    return 0;
} // end of main function

