#include "RunFlowWithoutParticle.h"
#include "GlobalFunctions.h"
#include "SimulationSettings.h"
#include "SimulationParameters.h"
#include "ParallelMPI.h"
#include "FlowSolver.h"
//=============================================================================================================================
void run_flow_without_particle(const SimulationParameters& paramsLB, const SimulationSettings& settingsLB, const ParallelMPI& parallelLB)
{
    bool Flow_With_Thermal_Effect = settingsLB.get_bool_setting("Flow_With_Thermal_Effect");
    bool Guo                      = settingsLB.get_bool_setting("Guo");
    bool Force_Thermal_2          = settingsLB.get_bool_setting("Force_Thermal_2");
    bool Recovery_ON              = settingsLB.get_bool_setting("Recovery_ON");
    UINT32 Create_Solid           = settingsLB.get_int_setting("Create_Solid");
    bool Read_Initialize_Data_From_File = settingsLB.get_bool_setting("Read_Initialize_Data_From_File");
    INT16 Solid_Thermal_BC    = settingsLB.get_int_setting("Solid_Thermal_BC");

    INT32 this_processor_id    = parallelLB.get_processor_id();
    INT32 num_processors_total = parallelLB.get_num_processors_total();// Should be after prepareMpiSetup. Before here it is false for all processors

    INT32 t_end      = paramsLB.get_t_end();
    INT32 t_data     = paramsLB.get_t_data();
    INT32 t_vtk      = paramsLB.get_t_vtk();
    INT32 t_info     = paramsLB.get_t_info();
    INT32 t_time     = paramsLB.get_t_time();
    INT32 t_recovery = paramsLB.get_t_recovery();
    INT32 Stop_by_Residual  = paramsLB.get_Stop_by_Residual();
    INT32 t_residual        = paramsLB.get_t_residual();

    std::vector<bool> is_processor_for_particle_only = parallelLB.get_is_processor_for_particle_only();// Should be after prepareMpiSetup. Before here it is false for all processors

    clock_t elapsed_time_0, elapsed_time_1, elapsed_time_2, elapsed_time_3;
    elapsed_time_0 = clock();
    INT32 recovery_step{}; //If "Recovery_ON" is activated the program starts from step=recovery_step+1

     /// Initialize the Flow
    FlowSolver flowSolverLB(paramsLB, settingsLB, parallelLB);//object of flowSolverLB class
    //Not needed here. flowSolverLB.initialize_flow(); // initialize variables for flow solver

    INT32 total_is_solid_points{0};
    /// Define solid nodes (is_solid) by MASTER processor
    if (this_processor_id == MASTER){
        if (Create_Solid == 0)
        {
            /// Do nothing
        } else if (Create_Solid == 1)///Create is_solid here
        {
            /// Important: Go to create_solid_from_dimensions and set is_solid locations
            flowSolverLB.create_solid_from_dimensions(parallelLB);
        } else if (Create_Solid == 2)
        {
            flowSolverLB.read_solid_from_file("geometry.txt", parallelLB);
        }

        total_is_solid_points = flowSolverLB.count_solid(parallelLB);// counts total number of is_solid
    }

    flowSolverLB.scatter_is_solid_points_from_master();// send is_solid from MASTER to others

    /// If Recovery_ON is active
    if (Recovery_ON)
    {
        if (this_processor_id == MASTER){
            flowSolverLB.read_fluid_recovery(recovery_step);// recovery_step is read from data_fluid_recovery
        }
        flowSolverLB.scatter_fluid_recovery_data_from_master(recovery_step);//send recovery data from MASTER to other processors
    }

    /// If Read_Initialize_Data_From_File is active
    if (Read_Initialize_Data_From_File){
        if (this_processor_id == MASTER){
            flowSolverLB.read_fluid_initialize_data(recovery_step, 2, 2, 2);
            recovery_step = 0;
        }
        flowSolverLB.scatter_fluid_recovery_data_from_master(recovery_step);//send recovery data from MASTER to other processors
    }

    /// Initialize populations
    flowSolverLB.initialize_flow(); // create folders and initialize populations (pop_F, pop_T) for flow solver

    if (this_processor_id == MASTER){printf("\n=========================== STARTING SIMULATION ===========================\n");}//srand(1);

    elapsed_time_1 = clock();
///===================================================================================================================================== ///
///                                                    MAIN TIME LOOP                                                                    ///
///===================================================================================================================================== ///
    for (INT32 t = 1 + recovery_step; t <= t_end; ++t)   // START OF TIME STEPS (run over all times between 1 and t_end)
    {
        /// ----------------------------- LBM AND THERMAL LBM -------------------------------------------
        if (this_processor_id != MASTER || num_processors_total == 1)
        {
            flowSolverLB.LBM_flow(t, parallelLB, nullptr);
            flowSolverLB.exchange_flow_data_neighbors("flow");
            flowSolverLB.bc_at_walls_flow(parallelLB);
            flowSolverLB.bc_at_solid_points_flow();
            flowSolverLB.momenta_flow();
        }

        if (Flow_With_Thermal_Effect)
        {
            if (this_processor_id != MASTER || num_processors_total == 1)
            {
                flowSolverLB.LBM_thermal(t);
                flowSolverLB.exchange_flow_data_neighbors("thermal");
                flowSolverLB.bc_at_walls_thermal(parallelLB);
                flowSolverLB.bc_at_solid_points_thermal();
                flowSolverLB.momenta_thermal();
            }
        }

        if (this_processor_id == MASTER || num_processors_total == 1)
        {
            flowSolverLB.set_Eulerian_forces_in_master_zero();
        }

        /// --------------------------------- GUO MODEL  ------------------------------------------------
        if (this_processor_id != MASTER || num_processors_total == 1)
        {
            if (Guo)
            {
                flowSolverLB.Guo_part_two(parallelLB, nullptr);
            }//End of: Guo
        }
        /// --------------------------------- FORCE_THERMAL_2  -----------------------------------------
        if (this_processor_id != MASTER || num_processors_total == 1)
        {
            if (Flow_With_Thermal_Effect)
            {
                if (Force_Thermal_2)
                {
                    flowSolverLB.force_thermal_2_part_two();
                }
            }
        }
        // It is better to gather_flow_data_to_master after Guo_part_two or force_thermal_2_part_two
        // (just to update velocity or temperature data for VTK and Fluid Data. No effect on code results.)

        if (t % t_data == 0 || t % t_vtk == 0 || t == 1){
            flowSolverLB.gather_fluid_macroscopic_data_to_master();// send density, temperature and velocity and previous_velocity and previous_temperature from slaves to MASTER
        }
        /// -------------------------------- WRITE FLUID DATA TO DISC  ------------------------------------
        if (t % t_data == 0)
        {
            if (this_processor_id == MASTER){
                flowSolverLB.write_fluid_data(t);
            }
        }

        /// ------------------------- WRITE RECOVERY DATA  -------------------------------------------------
        if (t % t_recovery == 0)
        {
            flowSolverLB.gather_fluid_recovery_data_to_master();// Send fluid recovery data to MASTER
            if (this_processor_id == MASTER){
                flowSolverLB.write_fluid_recovery(t); //Recommended: You can write fluid data into a "dat" file. In case of sudden program failure just activate "Recovery_ON" and run it so the program continues from a almost last iteration
            }
        }

        /// ------------------------- WRITE VTK DATA  ------------------------------------------------------
        if (t % t_vtk == 0 || t == 1)
        {
            if (this_processor_id == MASTER){
                flowSolverLB.write_fluid_vtk(t);
                flowSolverLB.write_solid_vtk(t, total_is_solid_points);
            }
        }

        /// ------------------------- REPORT TIME ON SCREEN  -----------------------------------------------
        if (t % t_time == 0)
        {
            if (this_processor_id == MASTER){
                elapsed_time_2 = clock();
                printf("The program took %f seconds from start of iterations.\n", ((float)elapsed_time_2 - elapsed_time_1) / CLOCKS_PER_SEC);
            }
        }

        if (t % t_info == 0)
        {
            if (this_processor_id == MASTER){
                printf("Completed Time Step %d in [1, %d]\n", t, t_end);
                printf("===========================================================================\n");
            }
        }

        ///-------------------------------- RESIDUAL CHECK -----------------------------------------------
        if (Stop_by_Residual && (t_residual > 0) && (t % t_residual == 0)) {

            bool crit_res = 0;
            crit_res = flowSolverLB.check_residual();

            if (this_processor_id == MASTER) {
                if (crit_res == 1) {
                    std::cout << "The program is finished by RESIDUAL condition.\n";
                    std::cout.flush();
                }
            }
            MPI_Bcast(&crit_res, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
            //all ranks exit consistently:
            if (crit_res) {
                MPI_Finalize();
                std::exit(0);
            }
        }// End of residual check

    }//End of time step loop
///===================================================================================================================================== ///
///                                                    END OF MAIN TIME LOOP                                                             ///
///===================================================================================================================================== ///

    /// Report successful end of simulation
    if (this_processor_id == MASTER)
    {
        printf("SIMULATION COMPLETED. THANKS \n");
        elapsed_time_3 = clock();
        printf("The program took %f seconds from beginning.\n", ((double)elapsed_time_3 - elapsed_time_0) / CLOCKS_PER_SEC);
        printf("The program took %f seconds from start of iterations.\n", ((double)elapsed_time_3 - elapsed_time_1) / CLOCKS_PER_SEC);
    }
}//End of: run_flow_without_particle
