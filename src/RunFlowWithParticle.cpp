#include "RunFlowWithParticle.h"
#include "GlobalFunctions.h"
#include "SimulationSettings.h"
#include "SimulationParameters.h"
#include "ParallelMPI.h"
#include "FlowSolver.h"
#include "ParticleResolved.h"
//=============================================================================================================================
void run_flow_with_particle(const SimulationParameters& paramsLB, const SimulationSettings& settingsLB, const ParallelMPI& parallelLB)
{
    bool Flow_With_Thermal_Effect = settingsLB.get_bool_setting("Flow_With_Thermal_Effect");
    bool Guo                      = settingsLB.get_bool_setting("Guo");
    bool Force_Thermal_2          = settingsLB.get_bool_setting("Force_Thermal_2");
    bool Recovery_ON              = settingsLB.get_bool_setting("Recovery_ON");
    bool Read_Initialize_Data_From_File = settingsLB.get_bool_setting("Read_Initialize_Data_From_File");

    INT32 this_processor_id    = parallelLB.get_processor_id();
    INT32 num_processors_total = parallelLB.get_num_processors_total();

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
    //Not needed here. flowSolverLB.initialize_flow(); // allocate memory and initialize variables for flow solver

    /// If Recovery_ON is active
    if (Recovery_ON)
    {
        if (this_processor_id == MASTER){
            flowSolverLB.read_fluid_recovery(recovery_step);
        }
        flowSolverLB.scatter_fluid_recovery_data_from_master(recovery_step);//send recovery data from MASTER to other processors
    }

    /// If Read_Initialize_Data_From_File is active
    if (Read_Initialize_Data_From_File){
        if (this_processor_id == MASTER){
            flowSolverLB.read_fluid_initialize_data(recovery_step, 5, 5, 5);
            recovery_step = 0;
        }
        flowSolverLB.scatter_fluid_recovery_data_from_master(recovery_step);//send recovery data from MASTER to other processors
    }

    /// Initialize populations
    flowSolverLB.initialize_flow(); // create folders and initialize populations (pop_F, pop_T) for flow solver

    /// Initialize Particles
    ParticleResolved particleResolvedLB(paramsLB, settingsLB, parallelLB);//object of ParticleResolved class
    particleResolvedLB.initialize_particles();
    if (Recovery_ON)
    {
        particleResolvedLB.read_particle_recovery();
    }   // Recovery_ON


    ///---------------------------------------------------------------------------///
    if (this_processor_id == MASTER){printf("\n=========================== STARTING SIMULATION ===========================\n");}//srand(1);

    elapsed_time_1 = clock();
///===================================================================================================================================== ///
///                                                    MAIN TIME LOOP                                                                    ///
///===================================================================================================================================== ///
    for (INT32 t = 1 + recovery_step; t <= t_end; ++t)   // START OF TIME STEPS (run over all times between 1 and t_end)
    {

        /// ----------------------------- LBM AND THERMAL LBM -------------------------------------------
        if (!is_processor_for_particle_only[this_processor_id])
        {
            flowSolverLB.LBM_flow(t, parallelLB, &particleResolvedLB);
            flowSolverLB.exchange_flow_data_neighbors("flow");
            flowSolverLB.bc_at_walls_flow(parallelLB);
            flowSolverLB.momenta_flow();
        }

        if (Flow_With_Thermal_Effect)
        {
            if (!is_processor_for_particle_only[this_processor_id])
            {
                flowSolverLB.LBM_thermal(t);
                flowSolverLB.exchange_flow_data_neighbors("thermal");
                flowSolverLB.bc_at_walls_thermal(parallelLB);
                flowSolverLB.momenta_thermal();
            }
        }

        if (is_processor_for_particle_only[this_processor_id] || num_processors_total == 1)
        {
            flowSolverLB.set_Eulerian_forces_in_master_zero();
        }

        ///The next 5 lines are only needed if you use scatter_...._optimized(parallelLB); (not needed for scatter_fluid_force_from_master())
        if (is_processor_for_particle_only[this_processor_id] || num_processors_total == 1)//equal to: if (this_processor_id == MASTER) for num_processors_for_particles_only=1 and num_processors_total>1
        {
            particleResolvedLB.calculate_force_modified_indices(flowSolverLB, parallelLB);// To calculate the indices of Eulerian nodes that are close to particle
        }
        flowSolverLB.scatter_force_modified_indices(parallelLB);// To send the indices of modified flow points (e.g. because of interpolation, spread)

        /// Gather velocity and temperature data to Master (needed for interpolation();)
        //flowSolverLB.gather_flow_data_to_master();
        flowSolverLB.gather_fluid_macroscopic_data_to_master_optimized();

        /// ----------------------------- EULERIAN FORCE INITIALIZATION ---------------------------------
        /// Important: Forces at Eulerian points must be set to ZERO at each iteration before SPREAD. We
        /// do this at the end of MOMENTA subroutine.
        /// This is done so just to reduce computational cost & time. (We can set all forces Zero here
        /// but it is time consuming)

        /// --------------------------- LAGRANGIAN FORCES CALCULATION AND SPREAD  -----------------------
        if (is_processor_for_particle_only[this_processor_id] || num_processors_total == 1)//equal to: if (this_processor_id == MASTER) for num_processors_for_particles_only=1 and num_processors_total>1
        {
            particleResolvedLB.direct_forcing(t, flowSolverLB, parallelLB);
        }

        /// ------------ SCATTER FLOW FORCES FROM PARTICLE PROCESSORS TO FLOW PROCESSORS  ---------------
        flowSolverLB.scatter_fluid_force_from_master_optimized(parallelLB);// Only force of Eulerian points around the particle are send to flow computing processors (Huge computational time saving)
        //flowSolverLB.scatter_fluid_force_from_master();// Sends all Eulerian points (even if no force computation has happened) to flow computing processors

        /// --------------------------------- GUO MODEL  ------------------------------------------------
        if (!is_processor_for_particle_only[this_processor_id])
        {
            if (Guo)
            {
                flowSolverLB.Guo_part_two(parallelLB, &particleResolvedLB);
            }//End of: Guo
        }
        /// --------------------------------- FORCE_THERMAL_2  -----------------------------------------
        if (!is_processor_for_particle_only[this_processor_id])
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

        /// ------------------------- NEWTON EQ.: UPDATE PARTICLE LOCATION AND ORIENTATION  ----------------
        if (is_processor_for_particle_only[this_processor_id] || num_processors_total == 1)
        {
            particleResolvedLB.newton_eqs_motion(t, flowSolverLB);
        }

        /// ----------- SCATTER PARTICLE DATA FROM MASTER (NEEDED FOR calculate_gbetaindex)  ---------------
        particleResolvedLB.scatter_particles_center_data_from_MASTER();
        /*if (this_processor_id == 5 && t == 2000){
            std::cout << "AAAA " << particleResolvedLB.particle_list[0].center[0].x << endl;
        }*/
        if (t % t_data == 0 || t % t_vtk == 0 || t == 1){
            flowSolverLB.gather_fluid_macroscopic_data_to_master();// send density, temperature and velocity from slaves to MASTER
        }
        /// ------------------------- WRITE FLUID AND PARTICLES DATA TO DISC  ------------------------------
        if (t % t_data == 0)
        {
            if (this_processor_id == MASTER){
                flowSolverLB.write_fluid_data(t);
                particleResolvedLB.write_particle_data(t);
            }
        }

        /// ------------------------- WRITE RECOVERY DATA  -------------------------------------------------
        if (t % t_recovery == 0)
        {
            flowSolverLB.gather_fluid_recovery_data_to_master();// Send fluid recovery data to MASTER
            if (this_processor_id == MASTER){
                flowSolverLB.write_fluid_recovery(t); //Recommended: You can write fluid data into a "dat" file. In case of sudden program failure just activate "Recovery_ON" and run it so the program continues from a almost last iteration
                particleResolvedLB.write_particle_recovery(t);
            }
        }

        /// ------------------------- WRITE VTK DATA  ------------------------------------------------------
        if (t % t_vtk == 0 || t == 1)
        {
            if (this_processor_id == MASTER){
                flowSolverLB.write_fluid_vtk(t);
                particleResolvedLB.write_particle_vtk(t);
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

}//End of: run_flow_with_particle
