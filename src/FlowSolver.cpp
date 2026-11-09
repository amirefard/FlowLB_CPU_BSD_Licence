#include "FlowSolver.h"
#include "ParticleResolved.h"
#include "SimulationParameters.h"
#include "SimulationSettings.h"
#include "ParallelMPI.h"
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Constructor
//-----------------------------------------------------------------------------------------------------------------------------
FlowSolver::FlowSolver(const SimulationParameters& paramsLB, const SimulationSettings& settingsLB, const ParallelMPI& paralLB): // Constructor with parameter
    //Settings from here
    Flow_With_Particle(settingsLB.get_bool_setting("Flow_With_Particle")),
    Flow_With_Thermal_Effect(settingsLB.get_bool_setting("Flow_With_Thermal_Effect")),
    X_Periodic(settingsLB.get_bool_setting("X_Periodic")),
    Y_Periodic(settingsLB.get_bool_setting("Y_Periodic")),
    Z_Periodic(settingsLB.get_bool_setting("Z_Periodic")),
    X_Wall(settingsLB.get_bool_setting("X_Wall")),
    Y_Wall(settingsLB.get_bool_setting("Y_Wall")),
    Z_Wall(settingsLB.get_bool_setting("Z_Wall")),
    X_velIn_ZeroGradOut(settingsLB.get_bool_setting("X_velIn_ZeroGradOut")),
    X_velIn_PressureOut(settingsLB.get_bool_setting("X_velIn_PressureOut")),
    X_PressureIn_PressureOut(settingsLB.get_bool_setting("X_PressureIn_PressureOut")),
    Parabolic_Velocity_Inlet(settingsLB.get_bool_setting("Parabolic_Velocity_Inlet")),
    Constant_Velocity_Inlet(settingsLB.get_bool_setting("Constant_Velocity_Inlet")),
    Recovery_OFF(settingsLB.get_bool_setting("Recovery_OFF")),
    Recovery_ON(settingsLB.get_bool_setting("Recovery_ON")),
    Shan_Chen(settingsLB.get_bool_setting("Shan_Chen")),
    Guo(settingsLB.get_bool_setting("Guo")),
    Luo(settingsLB.get_bool_setting("Luo")),
    SRT(settingsLB.get_bool_setting("SRT")),
    MRT(settingsLB.get_bool_setting("MRT")),
    Equilibrium_Thermal_1(settingsLB.get_bool_setting("Equilibrium_Thermal_1")),
    Equilibrium_Thermal_2(settingsLB.get_bool_setting("Equilibrium_Thermal_2")),
    Force_Thermal_1(settingsLB.get_bool_setting("Force_Thermal_1")),
    Force_Thermal_2(settingsLB.get_bool_setting("Force_Thermal_2")),
    Left_Fixed_Temp(settingsLB.get_bool_setting("Left_Fixed_Temp")),
    Left_Adiabatic(settingsLB.get_bool_setting("Left_Adiabatic")),
    Left_Fixed_Heat_Flux(settingsLB.get_bool_setting("Left_Fixed_Heat_Flux")),
    Right_Fixed_Temp(settingsLB.get_bool_setting("Right_Fixed_Temp")),
    Right_Adiabatic(settingsLB.get_bool_setting("Right_Adiabatic")),
    Right_Fixed_Heat_Flux(settingsLB.get_bool_setting("Right_Fixed_Heat_Flux")),
    Bottom_Fixed_Temp(settingsLB.get_bool_setting("Bottom_Fixed_Temp")),
    Bottom_Adiabatic(settingsLB.get_bool_setting("Bottom_Adiabatic")),
    Bottom_Fixed_Heat_Flux(settingsLB.get_bool_setting("Bottom_Fixed_Heat_Flux")),
    Top_Fixed_Temp(settingsLB.get_bool_setting("Top_Fixed_Temp")),
    Top_Adiabatic(settingsLB.get_bool_setting("Top_Adiabatic")),
    Top_Fixed_Heat_Flux(settingsLB.get_bool_setting("Top_Fixed_Heat_Flux")),
    Back_Fixed_Temp(settingsLB.get_bool_setting("Back_Fixed_Temp")),
    Back_Adiabatic(settingsLB.get_bool_setting("Back_Adiabatic")),
    Back_Fixed_Heat_Flux(settingsLB.get_bool_setting("Back_Fixed_Heat_Flux")),
    Front_Fixed_Temp(settingsLB.get_bool_setting("Front_Fixed_Temp")),
    Front_Adiabatic(settingsLB.get_bool_setting("Front_Adiabatic")),
    Front_Fixed_Heat_Flux(settingsLB.get_bool_setting("Front_Fixed_Heat_Flux")),
    VTK_ASCII(settingsLB.get_bool_setting("VTK_ASCII")),
    VTK_BINARY(settingsLB.get_bool_setting("VTK_BINARY")),
    CASE(settingsLB.get_bool_setting("CASE")),
    ASCII_Recovery_File(settingsLB.get_bool_setting("ASCII_Recovery_File")),
    BINARY_Recovery_File(settingsLB.get_bool_setting("BINARY_Recovery_File")),
    Solid_Thermal_BC(settingsLB.get_int_setting("Solid_Thermal_BC")),
    ///Parameters from SimulationParameters
    stencil_flow(paramsLB.get_stencil_flow()),
    stencil_thermal(paramsLB.get_stencil_thermal()),
    Nx(paramsLB.get_Nx()),
    Ny(paramsLB.get_Ny()),
    Nz(paramsLB.get_Nz()),
    tau(paramsLB.get_tau()),
    t_end(paramsLB.get_t_end()),
    t_data(paramsLB.get_t_data()),
    t_vtk(paramsLB.get_t_vtk()),
    t_info(paramsLB.get_t_info()),
    t_time(paramsLB.get_t_time()),
    t_recovery(paramsLB.get_t_recovery()),
    Stop_by_Residual(paramsLB.get_Stop_by_Residual()),
    residual_flow(paramsLB.get_residual_flow()),
    residual_thermal(paramsLB.get_residual_thermal()),
    t_residual(paramsLB.get_t_residual()),
    gravity(paramsLB.get_gravity()),
    wall_vel_bottom(paramsLB.get_wall_vel_bottom()),
    wall_vel_top(paramsLB.get_wall_vel_top()),
    wall_vel_back(paramsLB.get_wall_vel_back()),
    wall_vel_front(paramsLB.get_wall_vel_front()),
    inletvel(paramsLB.get_inletvel()),
    mdfstep(paramsLB.get_mdfstep()),
    fluid_init_temp(paramsLB.get_fluid_init_temp()),
    gbeta_x(paramsLB.get_gbeta_x()),
    gbeta_y(paramsLB.get_gbeta_y()),
    gbeta_z(paramsLB.get_gbeta_z()),
    Tref(paramsLB.get_Tref()),
    prandtl(paramsLB.get_prandtl()),
    particle_specific_heat_ratio(paramsLB.get_particle_specific_heat_ratio()),
    temperature_left(paramsLB.get_temperature_left()),
    temperature_right(paramsLB.get_temperature_right()),
    temperature_bot(paramsLB.get_temperature_bot()),
    temperature_top(paramsLB.get_temperature_top()),
    temperature_back(paramsLB.get_temperature_back()),
    temperature_front(paramsLB.get_temperature_front()),
    heat_flux_left(paramsLB.get_heat_flux_left()),
    heat_flux_right(paramsLB.get_heat_flux_right()),
    heat_flux_bottom(paramsLB.get_heat_flux_bottom()),
    heat_flux_top(paramsLB.get_heat_flux_top()),
    heat_flux_back(paramsLB.get_heat_flux_back()),
    heat_flux_front(paramsLB.get_heat_flux_front()),    
    domain_dimension_F(paramsLB.get_domain_dimension_F()),
    stencil_directions_F(paramsLB.get_stencil_directions_F()),
    domain_dimension_T(paramsLB.get_domain_dimension_T()),
    stencil_directions_T(paramsLB.get_stencil_directions_T()),
    omega(paramsLB.get_omega()),
    alpha_t(paramsLB.get_alpha_t()),
    tau_t(paramsLB.get_tau_t()),
    omega_t(paramsLB.get_omega_t()),
    num_particles(paramsLB.get_num_particles()),
    solid_temperature(paramsLB.get_solid_temperature()),
    solid_to_fluid_conductivity(paramsLB.get_solid_to_fluid_conductivity()),
    omega_ts(paramsLB.get_omega_ts()),

    ///Parameters from ParallelMPI
    this_processor_id(paralLB.get_processor_id()),
    num_processors_total(paralLB.get_num_processors_total()),
    start_global_x(paralLB.get_start_global_x()),
    end_global_x(paralLB.get_end_global_x()),
    start_global_y(paralLB.get_start_global_y()),
    end_global_y(paralLB.get_end_global_y()),
    start_global_z(paralLB.get_start_global_z()),
    end_global_z(paralLB.get_end_global_z()),
    start_local_x(paralLB.get_start_local_x()),
    end_local_x(paralLB.get_end_local_x()),
    start_local_y(paralLB.get_start_local_y()),
    end_local_y(paralLB.get_end_local_y()),
    start_local_z(paralLB.get_start_local_z()),
    end_local_z(paralLB.get_end_local_z()),
    local_neighbors_left(paralLB.get_local_neighbors_left()),
    local_neighbors_right(paralLB.get_local_neighbors_right()),
    local_neighbors_bottom(paralLB.get_local_neighbors_bottom()),
    local_neighbors_top(paralLB.get_local_neighbors_top()),
    local_neighbors_front(paralLB.get_local_neighbors_front()),
    local_neighbors_back(paralLB.get_local_neighbors_back()),
    length_with_buffer_X(paralLB.get_length_with_buffer_X()),
    length_with_buffer_Y(paralLB.get_length_with_buffer_Y()),
    length_with_buffer_Z(paralLB.get_length_with_buffer_Z()),
    num_buffer_layers_X(paralLB.get_num_buffer_layers_X()),
    num_buffer_layers_Y(paralLB.get_num_buffer_layers_Y()),
    num_buffer_layers_Z(paralLB.get_num_buffer_layers_Z())
{
    //set_simulation_settings(settingsLB);

    set_simulation_parameters();

    stensil_setup();

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: set_simulation_settings
//-----------------------------------------------------------------------------------------------------------------------------
/*void FlowSolver::set_simulation_settings(const SimulationSettings& settingsLB)
{
    Flow_With_Particle = settingsLB.get_bool_setting("Flow_With_Particle");
    Flow_With_Thermal_Effect = settingsLB.get_bool_setting("Flow_With_Thermal_Effect");
}*/
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: set_simulation_parameters
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::set_simulation_parameters()
{
    /*domain_dimension_F = params.get_domain_dimension_F();
    stencil_num_directions_F = params.get_stencil_directions_F();
    domain_dimension_T = params.get_domain_dimension_T();
    stencil_num_directions_T = params.get_stencil_directions_T();
    stencil_flow = params.get_stencil_flow();
    stencil_thermal = params.get_stencil_thermal();
    Nx = params.get_Nx();
    Ny = params.get_Ny();
    Nz = params.get_Nz();
    etc .....        */

    limitBound_left_X   = 0; //I recommend 0 (0: bc and lbm is from 0-Nx-1, 1: bc and lbm is from 1-Nx-2 for non-periodic case)
    limitBound_right_X  = 0; //I recommend 0
    limitBound_bottom_Y = 0; //I recommend 0
    limitBound_top_Y    = 0; //I recommend 0
    limitBound_back_Z   = 0; //I recommend 0
    limitBound_front_Z  = 0; //I recommend 0
     ///if you want the same result as Alborz in case of non periodic BC activate the next 6 lines. Because in ALBORZ bc and lbm are from 1:Nx-2, 1:Ny-2, 1:Nz-2 for non-periodic case
    /*if (!X_Periodic && start_global_x[this_processor_id] == 0)     {limitBound_left_X  = 1;}
    if (!X_Periodic && end_global_x[this_processor_id]   == Nx - 1){limitBound_right_X = 1;}
    if (!Y_Periodic && start_global_y[this_processor_id] == 0)     {limitBound_bottom_Y = 1;}
    if (!Y_Periodic && end_global_y[this_processor_id]   == Ny - 1){limitBound_top_Y = 1;}
    if (!Z_Periodic && start_global_z[this_processor_id] == 0)     {limitBound_back_Z = 1;}
    if (!Z_Periodic && end_global_z[this_processor_id]   == Nz - 1){limitBound_front_Z = 1;}*/
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: idx_3d_vector
//-----------------------------------------------------------------------------------------------------------------------------
INT32 FlowSolver::idx_3d_vector(INT32 x, INT32 y, INT32 z) const
{
    return (x * this->length_with_buffer_Y[this_processor_id] + y) * this->length_with_buffer_Z[this_processor_id] + z;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: idx_4d_vector(INT32 q, INT32 x, INT32 y, INT32 z)
//-----------------------------------------------------------------------------------------------------------------------------
INT32 FlowSolver::idx_4d_vector(INT32 q, INT32 x, INT32 y, INT32 z) const
{
    return ((q * this->length_with_buffer_X[this_processor_id] + x) * this->length_with_buffer_Y[this_processor_id] + y) * this->length_with_buffer_Z[this_processor_id] + z;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: decode_idx_3d(INT32 id3, INT32 &x, INT32 &y, INT32 &z)
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::decode_idx_3d(INT32 id3, INT32 &x, INT32 &y, INT32 &z) const
{
    const INT32 Ny_here = this->length_with_buffer_Y[this_processor_id];
    const INT32 Nz_here = this->length_with_buffer_Z[this_processor_id];

    z = id3 % Nz_here;        // remainder gives z
    id3 /= Nz_here;
    y = id3 % Ny_here;        // remainder gives y
    x = id3 / Ny_here;        // remaining quotient gives x
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: idx_4d_vector(INT32 q, INT32 id3)
//-----------------------------------------------------------------------------------------------------------------------------
/// Gets q, id3 and returns id4
INT32 FlowSolver::idx_4d_vector(INT32 q, INT32 id3) const
{
    INT32 x, y, z;
    decode_idx_3d(id3, x, y, z);
    return idx_4d_vector(q, x, y, z);
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: stensil_setup
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::stensil_setup()
{
    /*INT32 length_with_buffer_X = end_local_x[this_processor_id] - start_local_x[this_processor_id] + 1;
    INT32 length_with_buffer_Y = end_local_y[this_processor_id] - start_local_y[this_processor_id] + 1;
    INT32 length_with_buffer_Z = end_local_z[this_processor_id] - start_local_z[this_processor_id] + 1;*/

    force_modified_indices.clear();
    proc_indices.clear();

    #if STORE_1D_VECTOR
    density.resize      (1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 1.0);
    temperature.resize  (1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], fluid_init_temp);
    force_thermal.resize(1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    is_solid.resize     (1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], false);
    solid_thermal.resize(1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0);
    mdfTemperature.resize(1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);

    velocity.resize   (3 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    force.resize      (3 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    mdfVelocity.resize(3 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    mdfForce.resize   (3 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    cell_marked.resize(1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], false);
    previous_velocity.resize   (3 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    previous_temperature.resize(1 * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], fluid_init_temp);
    #endif

    #if !STORE_1D_VECTOR
    density       = vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 1.0)));
    temperature   = vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], fluid_init_temp)));
    force_thermal = vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0)));
    is_solid      = vector<vector<vector<bool>>>(length_with_buffer_X[this_processor_id], vector<vector<bool>>(length_with_buffer_Y[this_processor_id], vector<bool>(length_with_buffer_Z[this_processor_id], false)));
    solid_thermal = vector<vector<vector<INT16>>>(length_with_buffer_X[this_processor_id], vector<vector<INT16>>(length_with_buffer_Y[this_processor_id], vector<INT16>(length_with_buffer_Z[this_processor_id], 0)));
    velocity      = vector<vector<vector<vector<double>>>>(3, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    force         = vector<vector<vector<vector<double>>>>(3, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    mdfVelocity   = vector<vector<vector<vector<double>>>>(3, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    mdfForce      = vector<vector<vector<vector<double>>>>(3, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    cell_marked   = vector<vector<vector<bool>>>(length_with_buffer_X[this_processor_id], vector<vector<bool>>(length_with_buffer_Y[this_processor_id], vector<bool>(length_with_buffer_Z[this_processor_id], false)));
    mdfTemperature= vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0)));
    previous_velocity    = vector<vector<vector<vector<double>>>>(3, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    previous_temperature = vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], fluid_init_temp)));
    #endif

    /// Flow populations
    #if !STORE_1D_VECTOR
    this->pop_F     = vector<vector<vector<vector<double>>>>(stencil_directions_F, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    this->pop_old_F = vector<vector<vector<vector<double>>>>(stencil_directions_F, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    #endif
    #if STORE_1D_VECTOR
    pop_F.resize    (stencil_directions_F * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    pop_old_F.resize(stencil_directions_F * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    #endif
    this->cx_pop_F  = vector<INT32>(stencil_directions_F, 0);
    this->cy_pop_F  = vector<INT32>(stencil_directions_F, 0);
    this->cz_pop_F  = vector<INT32>(stencil_directions_F, 0);
    this->bounce_back_pop_F = vector<INT32>(stencil_directions_F, 0);

    /// Thermal populations
    #if !STORE_1D_VECTOR
    this->pop_T     = vector<vector<vector<vector<double>>>>(stencil_directions_T, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    this->pop_old_T = vector<vector<vector<vector<double>>>>(stencil_directions_T, vector<vector<vector<double>>>(length_with_buffer_X[this_processor_id], vector<vector<double>>(length_with_buffer_Y[this_processor_id], vector<double>(length_with_buffer_Z[this_processor_id], 0.0))));
    #endif
    #if STORE_1D_VECTOR
    pop_T.resize    (stencil_directions_T * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    pop_old_T.resize(stencil_directions_T * length_with_buffer_X[this_processor_id] * length_with_buffer_Y[this_processor_id] * length_with_buffer_Z[this_processor_id], 0.0);
    #endif
    this->cx_pop_T  = vector<INT32>(stencil_directions_T, 0);
    this->cy_pop_T  = vector<INT32>(stencil_directions_T, 0);
    this->cz_pop_T  = vector<INT32>(stencil_directions_T, 0);
    this->bounce_back_pop_T = vector<INT32>(stencil_directions_T, 0);

    // ########################### FLOW STENCIL SETUP ########################### //
    // ================ D2Q9 Flow stencil setup ================ //
    if (stencil_flow == "D2Q9")
    {
        cs_2_Flow = 1.0/3.0;

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
        cz_pop_F = {0, 0, 0, 0,  0, 0,  0,  0,  0};
        bounce_back_pop_F = {0, 3, 4, 1, 2, 7, 8, 5, 6};

        ///MRT variables
        ssf_MRT = 1. / tau;
	    s_MRT  = {1.0, 1.4, 1.4, 1.0, 1.2, 1.0, 1.2, ssf_MRT, ssf_MRT};//As sugested in page 103 of Mohamad book

        double M[9][9] =
        {
            { 1,  1,  1,  1,  1,  1,  1,  1,  1},

            {-4, -1, -1, -1, -1,  2,  2,  2,  2},

            { 4, -2, -2, -2, -2,  1,  1,  1,  1},

            { 0,  1,  0, -1,  0,  1, -1, -1,  1},

            { 0, -2,  0,  2,  0,  1, -1, -1,  1},

            { 0,  0,  1,  0, -1,  1,  1, -1, -1},

            { 0,  0, -2,  0,  2,  1,  1, -1, -1},

            { 0,  1, -1,  1, -1,  0,  0,  0,  0},

            { 0,  0,  0,  0,  0,  1, -1,  1, -1}
        };

        std::copy(&M[0][0], &M[0][0] + 9 * 9, &M_D2Q9[0][0]);

        double MInv[9][9] =
        {
            { 1.0/9.0,  -1.0/9.0,   1.0/9.0,    0.0,        0.0,        0.0,        0.0,        0.0,        0.0 },

            { 1.0/9.0,  -1.0/36.0, -1.0/18.0,   1.0/6.0,   -1.0/6.0,    0.0,        0.0,        1.0/4.0,    0.0 },

            { 1.0/9.0,  -1.0/36.0, -1.0/18.0,   0.0,        0.0,        1.0/6.0,   -1.0/6.0,   -1.0/4.0,    0.0 },

            { 1.0/9.0,  -1.0/36.0, -1.0/18.0,  -1.0/6.0,    1.0/6.0,    0.0,        0.0,        1.0/4.0,    0.0 },

            { 1.0/9.0,  -1.0/36.0, -1.0/18.0,   0.0,        0.0,       -1.0/6.0,    1.0/6.0,   -1.0/4.0,    0.0 },

            { 1.0/9.0,   1.0/18.0,  1.0/36.0,   1.0/6.0,    1.0/12.0,   1.0/6.0,    1.0/12.0,   0.0,         1.0/4.0 },

            { 1.0/9.0,   1.0/18.0,  1.0/36.0,  -1.0/6.0,   -1.0/12.0,   1.0/6.0,    1.0/12.0,   0.0,        -1.0/4.0 },

            { 1.0/9.0,   1.0/18.0,  1.0/36.0,  -1.0/6.0,   -1.0/12.0,  -1.0/6.0,   -1.0/12.0,   0.0,         1.0/4.0 },

            { 1.0/9.0,   1.0/18.0,  1.0/36.0,   1.0/6.0,    1.0/12.0,  -1.0/6.0,   -1.0/12.0,   0.0,        -1.0/4.0 }
        };

        std::copy(&MInv[0][0], &MInv[0][0] + 9 * 9, &this->Minv_D2Q9[0][0]);

        ///InvertMatrix(&M_D2Q9[0][0], &Minv_D2Q9[0][0], 9); // If you want to calculate the inverse of M by code instead of hardcoding it, you can use this line.
    }  
    
    // ================ D3Q19 Flow stencil setup ================ //
    if (stencil_flow == "D3Q19")
    {
        cs_2_Flow = 1.0/3.0;

        /// D3Q19 stencil:
        // index:   0  1  2  3  4  5  6  7   8  9  10  11  12  13  14  15  16  17  18
        // x:       0  1 -1  0  0  0  0  1  -1  1  -1  1   -1  1   -1   0   0   0   0
        // y:       0  0  0  1 -1  0  0  1  1  -1  -1  0   0   0   0    1  -1   1  -1
        // z:       0  0  0  0  0  1 -1  0  0   0  0   1   1   -1  -1   1   1  -1  -1
        weight_F =
        {
            1.0 / 3.0,
            1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
        };
        cx_pop_F = {0, 1, -1,  0, 0,  0,  0,  1,  -1,  1,  -1,  1,  -1,  1,   -1,  0,  0,   0,  0};
        cy_pop_F = {0, 0, 0,   1, -1, 0,  0,  1,  1,   -1, -1,  0,  0,   0,   0,   1,  -1,  1,  -1};
        cz_pop_F = {0, 0, 0,   0, 0,  1,  -1, 0,  0,   0,  0,   1,  1,   -1,  -1,  1,  1,   -1, -1};
        bounce_back_pop_F = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};

        ///MRT variables
        ssf_MRT = 1. / tau;
        s_MRT = { 0.0, 1.19, 1.4, 0.0, 1.2, 0.0, 1.2, 0.0, 1.2, ssf_MRT, 1.4, ssf_MRT, 1.4, ssf_MRT, ssf_MRT, ssf_MRT, 1.98, 1.98, 1.98 };//As sugested in page 103 of Mohamad book
        //s_MRT = {ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT};//If all coefficients are 1/tau then we have SRT
        double M[19][19] =
        {
            {  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},

            {-30,-11,-11,-11,-11,-11,-11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8},

            { 12,-4,-4,-4,-4,-4,-4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},

            { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0},

            { 0,-4, 4, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0},

            { 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1},

            { 0, 0, 0,-4, 4, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1},

            { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1},

            { 0, 0, 0, 0, 0,-4, 4, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1},

            { 0, 2, 2,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2},

            { 0,-4,-4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2},

            { 0, 0, 0, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},

            { 0, 0, 0,-2,-2, 2, 2, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1},

            { 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0},

            { 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1, 0, 0, 0, 0, 1,-1, 1,-1},

            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1,-1,-1, 1, 1}
        };

        std::copy(&M[0][0], &M[0][0] + 19 * 19, &M_D3Q19[0][0]);

        double MInv[19][19] =
        {
            {1.0/19.0, -5.0/399.0,  1.0/21.0,   0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0,          0},

            {1.0/19.0, -11.0/2394.0,-1.0/63.0,  1.0/10.0, -1.0/10.0,   0,          0,          0,          0,          1.0/18.0, -1.0/18.0,  0,          0,          0,          0,          0,          0,          0,          0},

            {1.0/19.0, -11.0/2394.0,-1.0/63.0, -1.0/10.0,  1.0/10.0,   0,          0,          0,          0,          1.0/18.0, -1.0/18.0,  0,          0,          0,          0,          0,          0,          0,          0},

            {1.0/19.0, -11.0/2394.0,-1.0/63.0,  0,          0,          1.0/10.0, -1.0/10.0,   0,          0,         -1.0/36.0,  1.0/36.0,  1.0/12.0, -1.0/12.0,  0,          0,          0,          0,          0,          0},

            {1.0/19.0, -11.0/2394.0,-1.0/63.0,  0,          0,         -1.0/10.0,  1.0/10.0,   0,          0,         -1.0/36.0,  1.0/36.0,  1.0/12.0, -1.0/12.0,  0,          0,          0,          0,          0,          0},

            {1.0/19.0, -11.0/2394.0,-1.0/63.0,  0,          0,          0,          0,          1.0/10.0, -1.0/10.0, -1.0/36.0,  1.0/36.0, -1.0/12.0,  1.0/12.0,  0,          0,          0,          0,          0,          0},

            {1.0/19.0, -11.0/2394.0,-1.0/63.0,  0,          0,          0,          0,         -1.0/10.0,  1.0/10.0, -1.0/36.0,  1.0/36.0, -1.0/12.0,  1.0/12.0,  0,          0,          0,          0,          0,          0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  1.0/10.0,  1.0/40.0,   1.0/10.0,  1.0/40.0,   0,          0,          1.0/36.0,  1.0/72.0,  1.0/12.0,  1.0/24.0,  1.0/4.0,   0,          0,          1.0/8.0,  -1.0/8.0,   0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0,   1.0/10.0,  1.0/40.0,   0,          0,          1.0/36.0,  1.0/72.0,  1.0/12.0,  1.0/24.0, -1.0/4.0,   0,          0,         -1.0/8.0,  -1.0/8.0,   0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  1.0/10.0,  1.0/40.0,  -1.0/10.0, -1.0/40.0,   0,          0,          1.0/36.0,  1.0/72.0,  1.0/12.0,  1.0/24.0, -1.0/4.0,   0,          0,          1.0/8.0,   1.0/8.0,   0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0,  -1.0/10.0, -1.0/40.0,   0,          0,          1.0/36.0,  1.0/72.0,  1.0/12.0,  1.0/24.0,  1.0/4.0,   0,          0,         -1.0/8.0,   1.0/8.0,   0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  1.0/10.0,  1.0/40.0,   0,          0,          1.0/10.0,  1.0/40.0,   1.0/36.0,  1.0/72.0, -1.0/12.0, -1.0/24.0,  0,          1.0/4.0,   0,         -1.0/8.0,   0,          1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0,   0,          0,          1.0/10.0,  1.0/40.0,   1.0/36.0,  1.0/72.0, -1.0/12.0, -1.0/24.0,  0,         -1.0/4.0,   0,          1.0/8.0,   0,          1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  1.0/10.0,  1.0/40.0,   0,          0,         -1.0/10.0, -1.0/40.0,   1.0/36.0,  1.0/72.0, -1.0/12.0, -1.0/24.0,  0,         -1.0/4.0,   0,         -1.0/8.0,   0,         -1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0,   0,          0,         -1.0/10.0, -1.0/40.0,   1.0/36.0,  1.0/72.0, -1.0/12.0, -1.0/24.0,  0,          1.0/4.0,   0,          1.0/8.0,   0,         -1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  0,          0,          1.0/10.0,  1.0/40.0,   1.0/10.0,  1.0/40.0,  -1.0/18.0, -1.0/36.0,  0,          0,          0,          0,          1.0/4.0,   0,          1.0/8.0,  -1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  0,          0,         -1.0/10.0, -1.0/40.0,   1.0/10.0,  1.0/40.0,  -1.0/18.0, -1.0/36.0,  0,          0,          0,          0,         -1.0/4.0,   0,         -1.0/8.0,  -1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  0,          0,          1.0/10.0,  1.0/40.0,  -1.0/10.0, -1.0/40.0,  -1.0/18.0, -1.0/36.0,  0,          0,          0,          0,         -1.0/4.0,   0,          1.0/8.0,   1.0/8.0},

            {1.0/19.0,  4.0/1197.0, 1.0/252.0,  0,          0,         -1.0/10.0, -1.0/40.0,  -1.0/10.0, -1.0/40.0,  -1.0/18.0, -1.0/36.0,  0,          0,          0,          0,          1.0/4.0,   0,         -1.0/8.0,   1.0/8.0}
        };

        std::copy(&MInv[0][0], &MInv[0][0] + 19 * 19, &this->Minv_D3Q19[0][0]);
        ///InvertMatrix(&M_D3Q19[0][0], &Minv_D3Q19[0][0], 19);// If you want to calculate the inverse of M by code instead of hardcoding it, you can use this line.
    }
    
    // ================ D3Q27 Flow stencil setup ================ //
    if (stencil_flow == "D3Q27")
    {
        cs_2_Flow = 1.0/3.0;

        /// D3Q27 stencil
        /// index:  0   1   2   3   4   5   6   7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26
        /// x:      0   1  -1   0   0   0   0   1   -1    1   -1    1   -1    1   -1    0    0    0    0    1   -1    1   -1    1   -1    1   -1
        /// y:      0   0   0   1  -1   0   0   1    1   -1   -1    0    0    0    0    1   -1    1   -1    1    1   -1   -1    1    1   -1   -1
        /// z:      0   0   0   0   0   1  -1   0    0    0    0    1    1   -1   -1    1    1   -1   -1    1    1    1    1   -1   -1   -1   -1

        weight_F =
        {
            8.0 / 27.0,

            2.0 / 27.0, 2.0 / 27.0,
            2.0 / 27.0, 2.0 / 27.0,
            2.0 / 27.0, 2.0 / 27.0,

            1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0,

            1.0 / 216.0, 1.0 / 216.0,
            1.0 / 216.0, 1.0 / 216.0,
            1.0 / 216.0, 1.0 / 216.0,
            1.0 / 216.0, 1.0 / 216.0
        };

        cx_pop_F =
        {
             0,
             1,-1, 0, 0, 0, 0,
             1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0,
             1,-1, 1,-1, 1,-1, 1,-1
        };

        cy_pop_F =
        {
             0,
             0, 0, 1,-1, 0, 0,
             1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1,
             1, 1,-1,-1, 1, 1,-1,-1
        };

        cz_pop_F =
        {
             0,
             0, 0, 0, 0, 1,-1,
             0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1,
             1, 1, 1, 1,-1,-1,-1,-1
        };

        bounce_back_pop_F =
        {
             0,
             2, 1,
             4, 3,
             6, 5,

            10, 9,
             8, 7,

            14,13,
            12,11,

            18,17,
            16,15,

            26,25,
            24,23,
            22,21,
            20,19
        };

        /// MRT variables
        ssf_MRT = 1.0 / tau;

        // For now use SRT-equivalent MRT
        s_MRT =
        {
            ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT,
            ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT,
            ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT,
            ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT,
            ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT, ssf_MRT,
            ssf_MRT, ssf_MRT
        };

        // NOTE:
        // D3Q27 MRT matrix M and inverse Minv are substantially larger.
        // In practice many solvers:
        // 1) use central-moment/cumulant formulations
        // 2) use SRT for D3Q27
        // 3) generate M automatically

        // Recommended temporary approach:
        // use D3Q27 with BGK/SRT collision until full MRT basis is added.
    }

    // ############################# THERMAL STENCIL SETUP ############################# //
    // ================ D2Q5 Thermal stencil setup ================ //
    if (stencil_thermal == "D2Q5")
    {
        cs_2_Thermal = 1.0/3.0;

        weight_T = {1. / 3., 1. / 6., 1. / 6., 1. / 6., 1. / 6.};
        cx_pop_T = {0, 1, 0, -1, 0};
        cy_pop_T = {0, 0, 1, 0, -1};
        cz_pop_T = {0, 0, 0, 0,  0};
        bounce_back_pop_T = {0, 3, 4, 1, 2};
    }
    // ================ D2Q9 Thermal stencil setup ================ //
    if (stencil_thermal == "D2Q9")
    {
        cs_2_Thermal = 1.0/3.0;

        weight_T = {4. / 9.,
                    1. / 9., 1. / 9., 1. / 9., 1. / 9.,
                    1. / 36., 1. / 36., 1. / 36., 1. / 36.};
        cx_pop_T = {0, 1, 0, -1, 0, 1, -1, -1,  1};
        cy_pop_T = {0, 0, 1, 0, -1, 1,  1, -1, -1};
        cz_pop_T = {0, 0, 0, 0,  0, 0,  0,  0,  0};
        bounce_back_pop_T = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    }
    // ================ D3Q7 Thermal stencil setup ================ //
    if (stencil_thermal == "D3Q7")
    {
        cs_2_Thermal = 1.0/4.0;

        weight_T =
        {
            1.0/4.0,
            1.0/8.0, 1.0/8.0,
            1.0/8.0, 1.0/8.0,
            1.0/8.0, 1.0/8.0
        };

        cx_pop_T =
        {
            0,
            1,-1,
            0, 0,
            0, 0
        };

        cy_pop_T =
        {
            0,
            0, 0,
            1,-1,
            0, 0
        };

        cz_pop_T =
        {
            0,
            0, 0,
            0, 0,
            1,-1
        };

        bounce_back_pop_T =
        {
            0,
            2,1,
            4,3,
            6,5
        };
    }
    // ================ D3Q19 Thermal stencil setup ================ //
    if (stencil_thermal == "D3Q19")
    {
        weight_T =
        {
            1.0 / 3.0,
            1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
        };
        cx_pop_T = {0, 1, -1,  0, 0,  0,  0,  1,  -1,  1,  -1,  1,  -1,  1,   -1,  0,  0,   0,  0};
        cy_pop_T = {0, 0, 0,   1, -1, 0,  0,  1,  1,   -1, -1,  0,  0,   0,   0,   1,  -1,  1,  -1};
        cz_pop_T = {0, 0, 0,   0, 0,  1,  -1, 0,  0,   0,  0,   1,  1,   -1,  -1,  1,  1,   -1, -1};
        bounce_back_pop_T = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
    }
    // ============================================================= //


    for (INT32 alpha = 0; alpha < stencil_directions_F; alpha++)
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
        if (cz_pop_F[alpha] > 0 )
        {
            pop_to_front_F.push_back(alpha);
            pop_reverse_front_F.push_back(bounce_back_pop_F[alpha]);
        }
        if (cz_pop_F[alpha] < 0 )
        {
            pop_to_back_F.push_back(alpha);
            pop_reverse_back_F.push_back(bounce_back_pop_F[alpha]);
        }
    }

    for (INT32 alpha = 0; alpha < stencil_directions_T; alpha++)
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
        if (cz_pop_T[alpha] > 0 )
        {
            pop_to_front_T.push_back(alpha);
            pop_reverse_front_T.push_back(bounce_back_pop_T[alpha]);
        }
        if (cz_pop_T[alpha] < 0 )
        {
            pop_to_back_T.push_back(alpha);
            pop_reverse_back_T.push_back(bounce_back_pop_T[alpha]);
        }
    }

    this->fforce_MRT = vector<double>(stencil_directions_F, 0.0);
    this->mforce_MRT = vector<double>(stencil_directions_F, 0.0);
    this->m_MRT      = vector<double>(stencil_directions_F, 0.0);

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: equilibrium_flow
//-----------------------------------------------------------------------------------------------------------------------------
double FlowSolver::equilibrium_flow(double den, double vel_x, double vel_y, double vel_z, INT32 i)
{
    double pop_eq_f{};
#if !STORE_1D_VECTOR
    pop_eq_f = weight_F[i] * den * (1 + 1.0 / cs_2_Flow * (cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z)
                                + 1.0 / (2.0 * cs_2_Flow * cs_2_Flow) * SQ(cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z)
                                - 1.0 / (2.0 * cs_2_Flow) * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
#endif
#if STORE_1D_VECTOR
    pop_eq_f = weight_F[i] * den * (1 + 1.0 / cs_2_Flow * (cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z)
                           + 1.0 / (2.0 * cs_2_Flow * cs_2_Flow) * SQ(cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z)
                           - 1.0 / (2.0 * cs_2_Flow) * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
#endif
    return pop_eq_f;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: meq_MRT
//-----------------------------------------------------------------------------------------------------------------------------
//calculation the equilibrium moment for MRT
double FlowSolver::meq_MRT(double RHO, double U, double V, double W, INT32 k)
{
    double x;
    switch (k)
    {
    case 0:
    {
        x = RHO;
        break;
    }
    case 1:
    {
        x = RHO * (-11 + 19 * (U * U + V * V + W * W));
        break;
    }
    case 2:
    {
        x = RHO * (3 - 11. / 2. * (U * U + V * V + W * W));
        break;
    }
    case 3:
    {
        x = RHO * U;
        break;
    }
    case 4:
    {
        x = RHO * (-2. / 3.) * U;
        break;
    }
    case 5:
    {
        x = RHO * V;
        break;
    }
    case 6:
    {
        x = RHO * (-2. / 3.) * V;
        break;
    }
    case 7:
    {
        x = RHO * W;
        break;
    }
    case 8:
    {
        x = RHO * (-2. / 3.) * W;
        break;
    }
    case 9:
    {
        x = RHO * (2 * U * U - V * V - W * W);
        break;
    }
    case 10:
    {
        x = -1. / 2. * RHO * (2. * U * U - V * V - W * W);
        break;
    }
    case 11:
    {
        x = RHO * (V * V - W * W);
        break;
    }
    case 12:
    {
        x = -1. / 2. * RHO * (V * V - W * W);
        break;
    }
    case 13:
    {
        x = RHO * (U * V);
        break;
    }
    case 14:
    {
        x = RHO * (V * W);
        break;
    }
    case 15:
    {
        x = RHO * (U * W);
        break;
    }
    case 16:
    {
        x = 0.0;
        break;
    }
    case 17:
    {
        x = 0.0;
        break;
    }
    case 18:
    {
        x = 0.0;
        break;
    }
    default:
        x = 0;
    }
    return x;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: initialize_flow
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::initialize_flow()
{
    /*INT32 length_with_buffer_X = paral.get_end_local_x()[paral.get_processor_id()] - paral.get_start_local_x()[paral.get_processor_id()] + 1;
    INT32 length_with_buffer_Y = paral.get_end_local_y()[paral.get_processor_id()] - paral.get_start_local_y()[paral.get_processor_id()] + 1;
    INT32 length_with_buffer_Z = paral.get_end_local_z()[paral.get_processor_id()] - paral.get_start_local_z()[paral.get_processor_id()] + 1;*/

    /// Create folders, delete data file
    // Make sure that the VTK folders exist.
    // Old file data.dat is deleted, if existing.
    //INT32 ignore; // ignore return value of system calls
    //ignore = system("mkdir FlowLB_Results\\vtk_fluid"); // create folder if not existing
    //ignore = system("mkdir FlowLB_Results\\vtk_particle"); // create folder if not existing
    //ignore = system("mkdir FlowLB_Results\\Data_particle"); // create folder if not existing
    //ignore = system("mkdir -p vtk_particle"); // create folder if not existing

    std::string base_folder = "FlowLB_Results";
    if (this_processor_id == MASTER)
    {
        if (Recovery_OFF)
        {
            #ifdef _WIN32 //for Windows system
            // Command to delete folder with /Q option
            std::string command1 = "rmdir /Q /S " + base_folder;
            INT32 result1 = std::system(command1.c_str());

            // Check if the command executed successfully.
            if (result1 == 0)
            {
                // Command executed successfully
                std::cout << "\nFolder '" << base_folder << "' deleted successfully.\n";
            }
            else
            {
                // Command failed
                std::cerr << "\nFailed to delete folder '" << base_folder << "'\n";
            }
            /*std::string command1 = "del /Q ";
            std::string path1 = "FlowLB_Results\\vtk_fluid\\*.*";
            system(command1.append(path1).c_str());*/
            #else   //#ifdef __linux__
            if (std::filesystem::exists(base_folder))
            {
                // Command to delete folder recursively
                std::string command1 = "rm -r " + base_folder;

                // Execute the command
                INT32 result1 = system(command1.c_str());
                if (result1 == 0)
                {
                    // Command executed successfully
                    std::cout << "\nFolder '" << base_folder << "' deleted successfully.\n";
                }
                else
                {
                    // Command failed
                    std::cerr << "\nFailed to delete '" << base_folder << "' folder\n";
                }
            }
            #endif
        }

        /// Create folder
        #if _WIN32 //for Windows system
        if (_mkdir(base_folder.c_str()) == 0 || errno == EEXIST)
        {
            std::cout << "Folder '" << base_folder << "' created successfully.\n";
            std::string vtk_fluid_folder = base_folder + "\\vtk_fluid";
            std::string vtk_particle_folder = base_folder + "\\vtk_particle";
            std::string data_particle_folder = base_folder + "\\Data_particle";
            _mkdir(vtk_fluid_folder.c_str());
            _mkdir(vtk_particle_folder.c_str());
            _mkdir(data_particle_folder.c_str());
        }
        #else  // #ifdef __linux__
         if (std::filesystem::create_directory(base_folder))
        {
            std::cout << "Folder '" << base_folder << "' created successfully.\n";
            // Create subfolders
            std::string vtk_fluid_folder = base_folder + "/vtk_fluid";
            std::string vtk_particle_folder = base_folder + "/vtk_particle";
            std::string data_particle_folder = base_folder + "/Data_particle";
            std::filesystem::create_directory(vtk_fluid_folder);
            std::filesystem::create_directory(vtk_particle_folder);
            std::filesystem::create_directory(data_particle_folder);
        }
        #endif
    }// End of MASTER process

    const INT32 start_x = start_local_x[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id];

    // Use the equilibrium populations corresponding to the initialized fluid density and velocity.
    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                {
                    double popEq_F{};
                    #if !STORE_1D_VECTOR
                    popEq_F   = equilibrium_flow(density[X][Y][Z], velocity[0][X][Y][Z], velocity[1][X][Y][Z], velocity[2][X][Y][Z], c_i);
                    pop_old_F[c_i][X][Y][Z] = popEq_F;
                    pop_F[c_i][X][Y][Z]     = popEq_F;
                    #endif
                    #if STORE_1D_VECTOR
                    INT32 index_pop      = idx_4d_vector(c_i, X, Y, Z);
                    popEq_F   = equilibrium_flow(density[idx_3d_vector(X, Y, Z)], velocity[idx_4d_vector(0, X, Y, Z)], velocity[idx_4d_vector(1, X, Y, Z)], velocity[idx_4d_vector(2, X, Y, Z)], c_i);
                    pop_old_F[index_pop] = popEq_F;
                    pop_F[index_pop]     = popEq_F;
                    #endif
                }
                for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                {
                    double popEq_T{};
                    #if !STORE_1D_VECTOR
                    popEq_T   = equilibrium_thermal(temperature[X][Y][Z], density[X][Y][Z], velocity[0][X][Y][Z], velocity[1][X][Y][Z], velocity[2][X][Y][Z], c_i);
                    pop_old_T[c_i][X][Y][Z] = popEq_T;
                    pop_T[c_i][X][Y][Z]     = popEq_T;
                    #endif
                    #if STORE_1D_VECTOR
                    INT32 index_pop      = idx_4d_vector(c_i, X, Y, Z);
                    popEq_T   = equilibrium_thermal(temperature[idx_3d_vector(X, Y, Z)],density[idx_3d_vector(X, Y, Z)], velocity[idx_4d_vector(0, X, Y, Z)], velocity[idx_4d_vector(1, X, Y, Z)], velocity[idx_4d_vector(2, X, Y, Z)], c_i);
                    pop_old_T[index_pop] = popEq_T;
                    pop_T[index_pop]     = popEq_T;
                    #endif
                }
            }
        }
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM_flow
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::LBM_flow(INT32 time, const ParallelMPI& paralLB, const ParticleResolved* Particle_Model)/// Important: I used a pointer for ParticleResolved because I want to set it nullptr for !Flow_With_Particle
{
    /*INT32 num_buffer_layers_X = paral.get_num_buffer_layers_X();
    INT32 num_buffer_layers_Y = paral.get_num_buffer_layers_Y();
    INT32 num_buffer_layers_Z = paral.get_num_buffer_layers_Z();*/

    /// Swap populations
    // The present code used old and new populations which are swapped at the beginning of each time step.
    // This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
    // This way, the old populations are not overwritten during propagation.
    // The resulting code is easier to write and to debug.
    // The memory requirement for the populations is twice as large.

    /*vector<vector<vector<vector<double>>>> swap_temp;
    swap_temp = pop_old_F;
    pop_old_F = pop_F;
    pop_F = swap_temp;*/
    std::swap(pop_F, pop_old_F); //much (almost 2.5 times) faster than the 3 above lines.


  /*  INT32 start_X, end_X;//Start and End of LBM calculation for processor.
    start_X = 0;    end_X = Nx - 1;
    if (X_Wall)    {  start_X = 1; end_X = Nx - 2; }
    if (X_velIn_ZeroGradOut)    {start_X = 1; end_X = Nx - 2; }
    if (ZeroGradIn_ZeroGradOut) {start_X = 1; end_X = Nx - 2; }
    INT32 start_Y, end_Y;
    start_Y = 0;  end_Y = Ny - 1;
    if (Y_Wall)    {start_Y = 1; end_Y = Ny - 2;}
    INT32 start_Z, end_Z;
    start_Z = 0;   end_Z = Nz - 1;
    if (Z_Wall)    {start_Z = 1; end_Z = Nz - 2;
    }*/

    /// Lattice Boltzmann equation
    // The lattice Boltzmann equation is solved in the following.
    // The algorithm includes
    // - computation of the lattice force
    // - combined collision and propagation (faster than first collision and then propagation)

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;

    double v0, v1, v2, rho;// This together with private () increases the openmp speed
    INT32 newX, newY, newZ;
    INT32 id3, id4_0, id4_1, id4_2;
    INT32 gbetaindex; //Inside particle =0, outside particle = 1

    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                double popEq_F{};
                id3 = idx_3d_vector(X, Y, Z);
                #if !STORE_1D_VECTOR
                if (!is_solid[X][Y][Z])
                #endif
                #if STORE_1D_VECTOR
                if (!is_solid[id3])
                #endif
                {
                gbetaindex = 1;
                if (Particle_Model && Flow_With_Particle && Flow_With_Thermal_Effect){
                    gbetaindex = Particle_Model->calculate_gbetaindex(X, Y, Z, paralLB);// To set gbetaindex=0 for fluid nodes inside the particle
                }

                    /// Compute new populations
                    // This is the lattice Boltzmann equation (combined collision and propagation) including external forcing.
                    // Periodicity of the lattice in x-direction is taken into account by the %-operator.
///===================================================================================================================================== ///
///                                                 SHAN-CHEN MODEL                                                                      ///
///===================================================================================================================================== ///
                    if (Shan_Chen)
                    {
                        for (INT32 i = 0; i < stencil_directions_F; i++)
                        {
                            #if !STORE_1D_VECTOR
                            double forceX = force[0][X][Y][Z] + gbetaindex * gbeta_x * (temperature[X][Y][Z] - Tref) + gravity;
                            double forceY = force[1][X][Y][Z] + gbetaindex * gbeta_y * (temperature[X][Y][Z] - Tref);
                            double forceZ = force[2][X][Y][Z] + gbetaindex * gbeta_z * (temperature[X][Y][Z] - Tref);
                            v0 = velocity[0][X][Y][Z] + forceX * tau;
                            v1 = velocity[1][X][Y][Z] + forceY * tau;
                            v2 = velocity[2][X][Y][Z] + forceZ * tau;
                            rho = density[X][Y][Z];

                            popEq_F = equilibrium_flow(rho, v0, v1, v2, i);

                            newX = X + cx_pop_F[i];//(X + cx_pop_F[i] + Nx) % Nx;
                            newY = Y + cy_pop_F[i];//(Y + cy_pop_F[i] + Ny) % Ny;
                            newZ = Z + cz_pop_F[i];//(Z + cz_pop_F[i] + Nz) % Nz;
                            pop_F[i][newX][newY][newZ] = pop_old_F[i][X][Y][Z] * (1.0 - omega) + popEq_F * omega;
                            #endif

                            #if STORE_1D_VECTOR
                            id4_0 = idx_4d_vector(0, X, Y, Z);
                            id4_1 = idx_4d_vector(1, X, Y, Z);
                            id4_2 = idx_4d_vector(2, X, Y, Z);
                            double forceX = force[id4_0] + gbetaindex * gbeta_x * (temperature[id3] - Tref) + gravity;
                            double forceY = force[id4_1] + gbetaindex * gbeta_y * (temperature[id3] - Tref);
                            double forceZ = force[id4_2] + gbetaindex * gbeta_z * (temperature[id3] - Tref);
                            v0    = velocity[id4_0] + forceX * tau;
                            v1    = velocity[id4_1] + forceY * tau;
                            v2    = velocity[id4_2] + forceZ * tau;
                            rho   = density[id3];

                            INT32 id_old = idx_4d_vector(i, X, Y, Z);
                            popEq_F = equilibrium_flow(rho, v0, v1, v2, i);

                            newX = X + cx_pop_F[i];//(X + cx_pop_F[i] + Nx) % Nx;
                            newY = Y + cy_pop_F[i];//(Y + cy_pop_F[i] + Ny) % Ny;
                            newZ = Z + cz_pop_F[i];//(Z + cz_pop_F[i] + Nz) % Nz;
                            INT32 id_new  = idx_4d_vector(i, newX, newY, newZ);
                            pop_F[id_new] = pop_old_F[id_old] * (1.0 - omega) + popEq_F * omega;
                            #endif
                        }
                    }
///===================================================================================================================================== ///
///                                                    LUO MODEL                                                                         ///
///===================================================================================================================================== ///
                    if (Luo)
                    {
                        for (INT32 i = 0; i < stencil_directions_F; i++)
                        {
                            #if !STORE_1D_VECTOR
                            popEq_F = equilibrium_flow(density[X][Y][Z], velocity[0][X][Y][Z], velocity[1][X][Y][Z], velocity[2][X][Y][Z], i);
                            newX = X + cx_pop_F[i];//(X + cx_pop_F[i] + Nx) % Nx;
                            newY = Y + cy_pop_F[i];//(Y + cy_pop_F[i] + Ny) % Ny;
                            newZ = Z + cz_pop_F[i];//(Z + cz_pop_F[i] + Nz) % Nz;
                            double forceX = force[0][X][Y][Z] + gbetaindex * gbeta_x * (temperature[X][Y][Z] - Tref) + gravity;
                            double forceY = force[1][X][Y][Z] + gbetaindex * gbeta_y * (temperature[X][Y][Z] - Tref);
                            double forceZ = force[2][X][Y][Z] + gbetaindex * gbeta_z * (temperature[X][Y][Z] - Tref);
                            pop_F[i][newX][newY][newZ] = pop_old_F[i][X][Y][Z] * (1.0 - omega) + popEq_F * omega
                                                         + 3 * weight_F[i] * (cx_pop_F[i] * forceX  + cy_pop_F[i] * forceY + cz_pop_F[i] * forceZ);
                            #endif
                            #if STORE_1D_VECTOR
                            id4_0 = idx_4d_vector(0, X, Y, Z);
                            id4_1 = idx_4d_vector(1, X, Y, Z);
                            id4_2 = idx_4d_vector(2, X, Y, Z);
                            double forceX = force[id4_0] + gbetaindex * gbeta_x * (temperature[id3] - Tref) + gravity;
                            double forceY = force[id4_1] + gbetaindex * gbeta_y * (temperature[id3] - Tref);
                            double forceZ = force[id4_2] + gbetaindex * gbeta_z * (temperature[id3] - Tref);
                            v0 = velocity[id4_0];
                            v1 = velocity[id4_1];
                            v2 = velocity[id4_2];

                            INT32 id_old = idx_4d_vector(i, X, Y, Z);
                            popEq_F = equilibrium_flow(density[id3], v0, v1, v2, i);

                            newX = X + cx_pop_F[i];//(X + cx_pop_F[i] + Nx) % Nx;
                            newY = Y + cy_pop_F[i];//(Y + cy_pop_F[i] + Ny) % Ny;
                            newZ = Z + cz_pop_F[i];//(Z + cz_pop_F[i] + Nz) % Nz;
                            INT32 id_new = idx_4d_vector(i, newX, newY, newZ);
                            pop_F[id_new] = pop_old_F[id_old] * (1.0 - omega) + popEq_F * omega
                                                         + 3 * weight_F[i] * (cx_pop_F[i] * forceX  + cy_pop_F[i] * forceY + cz_pop_F[i] * forceZ);
                            #endif
                        }
                    }
///===================================================================================================================================== ///
///                                                    GUO MODEL                                                                         ///
///===================================================================================================================================== ///
                    if (Guo)
                    {
                        if (SRT)
                        {
                            for (INT32 i = 0; i < stencil_directions_F; i++)
                            {
                                #if !STORE_1D_VECTOR
                                popEq_F = equilibrium_flow(density[X][Y][Z], velocity[0][X][Y][Z], velocity[1][X][Y][Z], velocity[2][X][Y][Z], i);
                                newX = X + cx_pop_F[i];//(X + cx_pop_F[i] + Nx) % Nx;
                                newY = Y + cy_pop_F[i];//(Y + cy_pop_F[i] + Ny) % Ny;
                                newZ = Z + cz_pop_F[i];//(Z + cz_pop_F[i] + Nz) % Nz;
                                double forceX = force[0][X][Y][Z] + gbetaindex * gbeta_x * (temperature[X][Y][Z] - Tref) + gravity;
                                double forceY = force[1][X][Y][Z] + gbetaindex * gbeta_y * (temperature[X][Y][Z] - Tref);
                                double forceZ = force[2][X][Y][Z] + gbetaindex * gbeta_z * (temperature[X][Y][Z] - Tref);
                                pop_F[i][newX][newY][newZ] = pop_old_F[i][X][Y][Z] * (1.0 - omega) + popEq_F * omega
                                                             + (1.0 - 0.5 * omega) * weight_F[i] * (forceX * (3 * (cx_pop_F[i] - velocity[0][X][Y][Z]) + 9.0 * cx_pop_F[i] * (cx_pop_F[i] * velocity[0][X][Y][Z] + cy_pop_F[i] * velocity[1][X][Y][Z] + cz_pop_F[i] * velocity[2][X][Y][Z]))
                                                                                                  + forceY * (3 * (cy_pop_F[i] - velocity[1][X][Y][Z]) + 9.0 * cy_pop_F[i] * (cx_pop_F[i] * velocity[0][X][Y][Z] + cy_pop_F[i] * velocity[1][X][Y][Z] + cz_pop_F[i] * velocity[2][X][Y][Z]))
                                                                                                  + forceZ * (3 * (cz_pop_F[i] - velocity[2][X][Y][Z]) + 9.0 * cz_pop_F[i] * (cx_pop_F[i] * velocity[0][X][Y][Z] + cy_pop_F[i] * velocity[1][X][Y][Z] + cz_pop_F[i] * velocity[2][X][Y][Z]))   );
                                #endif
                                #if STORE_1D_VECTOR
                                id4_0 = idx_4d_vector(0, X, Y, Z);
                                id4_1 = idx_4d_vector(1, X, Y, Z);
                                id4_2 = idx_4d_vector(2, X, Y, Z);
                                double forceX = force[id4_0] + gbetaindex * gbeta_x * (temperature[id3] - Tref) + gravity;
                                double forceY = force[id4_1] + gbetaindex * gbeta_y * (temperature[id3] - Tref);
                                double forceZ = force[id4_2] + gbetaindex * gbeta_z * (temperature[id3] - Tref);
                                v0 = velocity[id4_0];
                                v1 = velocity[id4_1];
                                v2 = velocity[id4_2];

                                INT32 id_old = idx_4d_vector(i, X, Y, Z);
                                popEq_F = equilibrium_flow(density[id3], v0, v1, v2, i);

                                newX = X + cx_pop_F[i];//(X + cx_pop_F[i] + Nx) % Nx;
                                newY = Y + cy_pop_F[i];//(Y + cy_pop_F[i] + Ny) % Ny;
                                newZ = Z + cz_pop_F[i];//(Z + cz_pop_F[i] + Nz) % Nz;
                                INT32 id_new = idx_4d_vector(i, newX, newY, newZ);
                                pop_F[id_new] = pop_old_F[id_old] * (1.0 - omega) + popEq_F * omega
                                                             + (1.0 - 0.5 * omega) * weight_F[i] * (forceX * (3 * (cx_pop_F[i] - velocity[id4_0]) + 9.0 * cx_pop_F[i] * (cx_pop_F[i] * velocity[id4_0] + cy_pop_F[i] * velocity[id4_1] + cz_pop_F[i] * velocity[id4_2]))
                                                                                                  + forceY * (3 * (cy_pop_F[i] - velocity[id4_1]) + 9.0 * cy_pop_F[i] * (cx_pop_F[i] * velocity[id4_0] + cy_pop_F[i] * velocity[id4_1] + cz_pop_F[i] * velocity[id4_2]))
                                                                                                  + forceZ * (3 * (cz_pop_F[i] - velocity[id4_2]) + 9.0 * cz_pop_F[i] * (cx_pop_F[i] * velocity[id4_0] + cy_pop_F[i] * velocity[id4_1] + cz_pop_F[i] * velocity[id4_2])) );
                                #endif
                            }
                        }
                        /// ================ MRT ================ ///
                        if (MRT)
                        {
                            #if !STORE_1D_VECTOR
                            std::cerr << "MRT with 4D vector is not implemented yet. Please set STORE_1D_VECTOR to 1 in the header file and recompile the code." << std::endl;
                            #endif

                            #if STORE_1D_VECTOR
                            // f -> m
                            if (stencil_flow == "D3Q19"){
                                TransformPopulationToMomentSpace(pop_old_F, m_MRT, &M_D3Q19[0][0], 19, X, Y, Z);
                            }else if (stencil_flow == "D2Q9"){
                                TransformPopulationToMomentSpace(pop_old_F, m_MRT, &M_D2Q9[0][0], 9, X, Y, Z);
                            }
                            double forceX = force[id4_0] + gbetaindex * gbeta_x * (temperature[id3] - Tref) + gravity;
                            double forceY = force[id4_1] + gbetaindex * gbeta_y * (temperature[id3] - Tref);
                            double forceZ = force[id4_2] + gbetaindex * gbeta_z * (temperature[id3] - Tref);
                            
                            for (INT32 i = 0; i < stencil_directions_F; i++)
                            {
                                double vel_x = velocity[id4_0];
                                double vel_y = velocity[id4_1];
                                double vel_z = velocity[id4_2];
                                fforce_MRT[i] = weight_F[i] * (forceX * (3 * (cx_pop_F[i] - vel_x) + 9.0 * cx_pop_F[i] * (cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z))
                                                             + forceY * (3 * (cy_pop_F[i] - vel_y) + 9.0 * cy_pop_F[i] * (cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z))
                                                             + forceZ * (3 * (cz_pop_F[i] - vel_z) + 9.0 * cz_pop_F[i] * (cx_pop_F[i] * vel_x + cy_pop_F[i] * vel_y + cz_pop_F[i] * vel_z)));                                
                            }

                            // mforce_MRT = M * fforce_MRT
                            if (stencil_flow == "D3Q19"){
                                TransformForceToMomentSpace(fforce_MRT, mforce_MRT, &M_D3Q19[0][0], 19);
                            }else if (stencil_flow == "D2Q9"){
                                TransformForceToMomentSpace(fforce_MRT, mforce_MRT, &M_D2Q9[0][0], 9);
                            }
                            // relaxation of m_MRT
                            for (INT32 k = 0; k < stencil_directions_F; k++)
                            {
                                double MEQ = meq_MRT(density[id3], velocity[id4_0], velocity[id4_1], velocity[id4_2], k);
                                m_MRT[k] = m_MRT[k] - s_MRT[k] * (m_MRT[k] - MEQ) + (1 - 1. / 2. * s_MRT[k]) * mforce_MRT[k]; // relaxation
                                mforce_MRT[k] = 0.0;
                                fforce_MRT[k] = 0.0;
                            }

                            // m -> f*
                            if (stencil_flow == "D3Q19"){
                                TransformToPopulationSpace(m_MRT, f_post_MRT, &Minv_D3Q19[0][0], stencil_directions_F);
                            }else if (stencil_flow == "D2Q9"){
                                TransformToPopulationSpace(m_MRT, f_post_MRT, &Minv_D2Q9[0][0], stencil_directions_F);
                            }

                            // stream f* -> f
                            StreamPostCollisionPopulations(f_post_MRT, pop_F, X, Y, Z);
                            #endif // STORE_1D_VECTOR
                        } //End of MRT
                    } //End of GUO
                } //if (is_solid[X][Y][Z] == FALSE)
            } //for (INT32 Z = start_Z; Z <= end_Z; ++Z)
        } // for (INT32 Y = start_Y; Y <= end_Y; ++Y)
    } // for (INT32 X = start_X; X <= end_X; ++X)

  //  return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: TransformPopulationToMomentSpace
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::TransformPopulationToMomentSpace(
    const std::vector<double>& pop_old_F,
    std::vector<double>& m_MRT,
    const double* M,
    INT32 Q,
    INT32 X,
    INT32 Y,
    INT32 Z)
{
    for (INT32 k = 0; k < Q; k++)
    {
        m_MRT[k] = 0.0;

        for (INT32 i = 0; i < Q; i++)
        {
            m_MRT[k] +=
                M[k * Q + i]
                * pop_old_F[idx_4d_vector(i, X, Y, Z)];
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: TransformForceToMomentSpace
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::TransformForceToMomentSpace(
    const std::vector<double>& f,
    std::vector<double>& m,
    const double* M,
    int Q)
{
    for (int k = 0; k < Q; ++k)
    {
        m[k] = 0.0;

        for (int i = 0; i < Q; ++i)
        {
            m[k] += M[k * Q + i] * f[i];
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: TransformToPopulationSpace
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::TransformToPopulationSpace(
    const std::vector<double>& m,
    std::vector<double>& f,
    const double* Minv,
    const INT32 Q)
{
    for (INT32 i = 0; i < Q; i++)
    {
        f[i] = 0.0;

        for (INT32 k = 0; k < Q; k++)
        {
            f[i] += Minv[i * Q + k] * m[k];
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: InvertMatrix
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::InvertMatrix(
    const double* A,
    double* Ainv,
    int N)
{
    std::vector<double> aug(N * 2 * N);

    // Build augmented matrix [A | I]
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            aug[i * (2 * N) + j] = A[i * N + j];
            aug[i * (2 * N) + N + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan elimination with partial pivoting
    for (int i = 0; i < N; i++)
    {
        // Find pivot row
        int pivot_row = i;
        double max_val = std::abs(aug[i * (2 * N) + i]);

        for (int k = i + 1; k < N; k++)
        {
            double val = std::abs(aug[k * (2 * N) + i]);

            if (val > max_val)
            {
                max_val = val;
                pivot_row = k;
            }
        }

        // Check singularity
        if (max_val < 1e-14)
        {
            throw std::runtime_error("Matrix is singular");
        }

        // Swap rows if needed
        if (pivot_row != i)
        {
            for (int j = 0; j < 2 * N; j++)
            {
                std::swap(
                    aug[i * (2 * N) + j],
                    aug[pivot_row * (2 * N) + j]);
            }
        }

        // Normalize pivot row
        double pivot = aug[i * (2 * N) + i];

        for (int j = 0; j < 2 * N; j++)
        {
            aug[i * (2 * N) + j] /= pivot;
        }

        // Eliminate other rows
        for (int k = 0; k < N; k++)
        {
            if (k == i) continue;

            double factor = aug[k * (2 * N) + i];

            for (int j = 0; j < 2 * N; j++)
            {
                aug[k * (2 * N) + j]
                    -= factor * aug[i * (2 * N) + j];
            }
        }
    }

    // Extract inverse
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Ainv[i * N + j]
                = aug[i * (2 * N) + N + j];
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: StreamPostCollisionPopulations
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::StreamPostCollisionPopulations(
    const std::vector<double>& f_post,
    std::vector<double>& pop_F,
    INT32 X,
    INT32 Y,
    INT32 Z)
{
    for (INT32 i = 0; i < stencil_directions_F; i++)
    {
        pop_F[idx_4d_vector(i, X + cx_pop_F[i], Y + cy_pop_F[i], Z + cz_pop_F[i])] = f_post[i];
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: exchange_flow_data_neighbors
//-----------------------------------------------------------------------------------------------------------------------------
//void FlowSolver::exchange_flow_data_neighbors(const ParallelMPI& paral, const std::string& flow_or_thermal)
void FlowSolver::exchange_flow_data_neighbors(const std::string& flow_or_thermal)
{
    //INT32 processor_id = paral.get_processor_id();

    /*INT32 num_buffer_layers_X = paral.get_num_buffer_layers_X();
    INT32 num_buffer_layers_Y = paral.get_num_buffer_layers_Y();
    INT32 num_buffer_layers_Z = paral.get_num_buffer_layers_Z();*/

    /*std::vector<std::vector<INT32>> local_neighbors_left   = paral.get_local_neighbors_left();
    std::vector<std::vector<INT32>> local_neighbors_right  = paral.get_local_neighbors_right();
    std::vector<std::vector<INT32>> local_neighbors_bottom = paral.get_local_neighbors_bottom();
    std::vector<std::vector<INT32>> local_neighbors_top    = paral.get_local_neighbors_top();
    std::vector<std::vector<INT32>> local_neighbors_front  = paral.get_local_neighbors_front();
    std::vector<std::vector<INT32>> local_neighbors_back   = paral.get_local_neighbors_back();*/

    /*INT32 start_local_x = paral.get_start_local_x()[this_processor_id];
    INT32 end_local_x   = paral.get_end_local_x()[this_processor_id];
    INT32 start_local_y = paral.get_start_local_y()[this_processor_id];
    INT32 end_local_y   = paral.get_end_local_y()[this_processor_id];
    INT32 start_local_z = paral.get_start_local_z()[this_processor_id];
    INT32 end_local_z   = paral.get_end_local_z()[this_processor_id];*/


    MPI_Status status;
    MPI_Request request_right, request_left;

    for (INT32 ii = 0; ii < 6; ii++)
    {
        INT32 sending_to_neighbor, receiving_from_neighbor;
        INT32 start_send_X, end_send_X, start_send_Y, end_send_Y, start_send_Z, end_send_Z;
        INT32 start_receive_X, end_receive_X, start_receive_Y, end_receive_Y, start_receive_Z, end_receive_Z;
        INT32 tag;
        MPI_Request request_send, request_left;
        std::vector<INT32> pop_to;
        std::vector<double> send_buffer, recv_buffer;

        pop_to.clear();

        switch(ii)
        {
        case 0:
            {
            /// Send to right
            sending_to_neighbor     = local_neighbors_right[this_processor_id][0]; // Assuming at most one neighbor in the right direction
            receiving_from_neighbor = local_neighbors_left[this_processor_id][0];
            if (flow_or_thermal == "flow"){pop_to = pop_to_right_F;}else if(flow_or_thermal == "thermal"){pop_to = pop_to_right_T;}
            start_send_X    = end_local_x[this_processor_id] - num_buffer_layers_X[this_processor_id] + 1;
            end_send_X      = end_local_x[this_processor_id] - num_buffer_layers_X[this_processor_id] + 1;
            start_send_Y    = 0;
            end_send_Y      = end_local_y[this_processor_id];
            start_send_Z    = 0;
            end_send_Z      = end_local_z[this_processor_id];
            start_receive_X = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
            end_receive_X   = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
            start_receive_Y = 0;
            end_receive_Y   = end_local_y[this_processor_id];
            start_receive_Z = 0;
            end_receive_Z   = end_local_z[this_processor_id];
            tag = ii;
            }
            break;
        case 1:
            {
            /// Send to left
            sending_to_neighbor     = local_neighbors_left[this_processor_id][0];
            receiving_from_neighbor = local_neighbors_right[this_processor_id][0];
            if (flow_or_thermal == "flow"){pop_to = pop_to_left_F;}else if(flow_or_thermal == "thermal"){pop_to = pop_to_left_T;}
            start_send_X    = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] - 1;
            end_send_X      = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] - 1;
            start_send_Y    = 0;
            end_send_Y      = end_local_y[this_processor_id];
            start_send_Z    = 0;
            end_send_Z      = end_local_z[this_processor_id];
            start_receive_X = end_local_x[this_processor_id] - num_buffer_layers_X[this_processor_id];
            end_receive_X   = end_local_x[this_processor_id] - num_buffer_layers_X[this_processor_id];
            start_receive_Y = 0;
            end_receive_Y   = end_local_y[this_processor_id];
            start_receive_Z = 0;
            end_receive_Z   = end_local_z[this_processor_id];
            tag = ii;
            }
            break;
        case 2:
            {
            /// Send to top
            sending_to_neighbor     = local_neighbors_top[this_processor_id][0]; // Assuming at most one neighbor in the right direction
            receiving_from_neighbor = local_neighbors_bottom[this_processor_id][0];
            if (flow_or_thermal == "flow"){pop_to = pop_to_top_F;}else if(flow_or_thermal == "thermal"){pop_to = pop_to_top_T;}
            start_send_Y    = end_local_y[this_processor_id] - num_buffer_layers_Y[this_processor_id] + 1;
            end_send_Y      = end_local_y[this_processor_id] - num_buffer_layers_Y[this_processor_id] + 1;
            start_send_X    = 0;
            end_send_X      = end_local_x[this_processor_id];
            start_send_Z    = 0;
            end_send_Z      = end_local_z[this_processor_id];
            start_receive_Y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
            end_receive_Y   = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
            start_receive_X = 0;
            end_receive_X   = end_local_x[this_processor_id];
            start_receive_Z = 0;
            end_receive_Z   = end_local_z[this_processor_id];
            tag = ii;
            }
            break;
        case 3:
            {
            /// Send to bottom
            sending_to_neighbor     = local_neighbors_bottom[this_processor_id][0];
            receiving_from_neighbor = local_neighbors_top[this_processor_id][0];
            if (flow_or_thermal == "flow"){pop_to = pop_to_bottom_F;}else if(flow_or_thermal == "thermal"){pop_to = pop_to_bottom_T;}
            start_send_Y    = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] - 1;
            end_send_Y      = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] - 1;
            start_send_X    = 0;
            end_send_X      = end_local_x[this_processor_id];
            start_send_Z    = 0;
            end_send_Z      = end_local_z[this_processor_id];
            start_receive_Y = end_local_y[this_processor_id] - num_buffer_layers_Y[this_processor_id];
            end_receive_Y   = end_local_y[this_processor_id] - num_buffer_layers_Y[this_processor_id];
            start_receive_X = 0;
            end_receive_X   = end_local_x[this_processor_id];
            start_receive_Z = 0;
            end_receive_Z   = end_local_z[this_processor_id];
            tag = ii;
            }
            break;
        case 4:
            {
            /// Send to front
            sending_to_neighbor     = local_neighbors_front[this_processor_id][0]; // Assuming at most one neighbor in the right direction
            receiving_from_neighbor = local_neighbors_back[this_processor_id][0];
            if (flow_or_thermal == "flow"){pop_to = pop_to_front_F;}else if(flow_or_thermal == "thermal"){pop_to = pop_to_front_T;}
            start_send_Z    = end_local_z[this_processor_id] - num_buffer_layers_Z[this_processor_id] + 1;
            end_send_Z      = end_local_z[this_processor_id] - num_buffer_layers_Z[this_processor_id] + 1;
            start_send_Y    = 0;
            end_send_Y      = end_local_y[this_processor_id];
            start_send_X    = 0;
            end_send_X      = end_local_x[this_processor_id];
            start_receive_Z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
            end_receive_Z   = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
            start_receive_Y = 0;
            end_receive_Y   = end_local_y[this_processor_id];
            start_receive_X = 0;
            end_receive_X   = end_local_x[this_processor_id];
            tag = ii;
            }
            break;
        case 5:
            {
            /// Send to back
            sending_to_neighbor     = local_neighbors_back[this_processor_id][0];
            receiving_from_neighbor = local_neighbors_front[this_processor_id][0];
            if (flow_or_thermal == "flow"){pop_to = pop_to_back_F;}else if(flow_or_thermal == "thermal"){pop_to = pop_to_back_T;}
            start_send_Z    = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] - 1;
            end_send_Z      = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] - 1;
            start_send_Y    = 0;
            end_send_Y      = end_local_y[this_processor_id];
            start_send_X    = 0;
            end_send_X      = end_local_x[this_processor_id];
            start_receive_Z = end_local_z[this_processor_id] - num_buffer_layers_Z[this_processor_id];
            end_receive_Z   = end_local_z[this_processor_id] - num_buffer_layers_Z[this_processor_id];
            start_receive_Y = 0;
            end_receive_Y   = end_local_y[this_processor_id];
            start_receive_X = 0;
            end_receive_X   = end_local_x[this_processor_id];
            tag = ii;
            }
            break;
        default:
            std::cout << "Invalid choice" << std::endl;
            break;
        }

        send_buffer.clear();
        recv_buffer.clear();

        for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
        {
            for (INT32 X = start_send_X; X <= end_send_X; ++X)
            {
                for (INT32 Y = start_send_Y; Y <= end_send_Y; ++Y)
                {
                    for (INT32 Z = start_send_Z; Z <= end_send_Z; ++Z)
                    {
                        if (flow_or_thermal == "flow")
                        {
                            #if !STORE_1D_VECTOR
                            send_buffer.push_back(pop_F[pop_to[alpha]][X][Y][Z]);
                            #endif
                            #if STORE_1D_VECTOR
                            INT32 id4 = idx_4d_vector(pop_to[alpha], X, Y, Z);
                            send_buffer.push_back(pop_F[id4]);
                            #endif
                        }
                        else if (flow_or_thermal == "thermal")
                        {
                            #if !STORE_1D_VECTOR
                            send_buffer.push_back(pop_T[pop_to[alpha]][X][Y][Z]);
                            #endif
                            #if STORE_1D_VECTOR
                            INT32 id4 = idx_4d_vector(pop_to[alpha], X, Y, Z);
                            send_buffer.push_back(pop_T[id4]);
                            #endif
                        }
                    }
                }
            }
        }

        //std::cout << "Processor " << processor_id << " sending to " << sending_to_neighbor << " with buffer size: " << send_buffer.size() << std::endl;

        recv_buffer.resize(send_buffer.size());

        // Use MPI_Sendrecv to send and receive in a single call
        MPI_Sendrecv(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, sending_to_neighbor, tag,
                     recv_buffer.data(), recv_buffer.size(), MPI_DOUBLE, receiving_from_neighbor, tag,
                     MPI_COMM_WORLD, &status);


        /* MPI_Isend(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, sending_to_neighbor, tag, MPI_COMM_WORLD, &request_right);//Send to the right

         MPI_Irecv(recv_buffer.data(), recv_buffer.size(), MPI_DOUBLE, receiving_from_neighbor, tag, MPI_COMM_WORLD, &request_left);
             // Wait for receive operation to complete
         MPI_Wait(&request_left, MPI_STATUS_IGNORE);*/

        // Integrate received data into pop_F (or pop_T)
        INT32 index = 0;
        for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
        {
            for (INT32 X = start_receive_X; X <= end_receive_X; ++X)
            {
                for (INT32 Y = start_receive_Y; Y <= end_receive_Y; ++Y)
                {
                    for (INT32 Z = start_receive_Z; Z <= end_receive_Z; ++Z)
                    {
                        if (flow_or_thermal == "flow")
                        {
                            #if !STORE_1D_VECTOR
                            pop_F[pop_to[alpha]][X][Y][Z] = recv_buffer[index++];
                            #endif
                            #if STORE_1D_VECTOR
                            INT32 id4 = idx_4d_vector(pop_to[alpha], X, Y, Z);
                            pop_F[id4] = recv_buffer[index++];
                            #endif // STORE_1D_VECTOR

                        }
                        else if (flow_or_thermal == "thermal")
                        {
                            #if !STORE_1D_VECTOR
                            pop_T[pop_to[alpha]][X][Y][Z] = recv_buffer[index++];
                            #endif
                            #if STORE_1D_VECTOR
                            INT32 id4 = idx_4d_vector(pop_to[alpha], X, Y, Z);
                            pop_T[id4] = recv_buffer[index++];
                            #endif // STORE_1D_VECTOR
                        }
                    }
                }
            }
        }
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: bc_at_walls_flow
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::bc_at_walls_flow(const ParallelMPI& paralLB)
{

    INT32 newX, newY, newZ;
    for (INT32 ii = 0; ii < 6; ii++)
    {
        INT32 start_BC_X{}, end_BC_X{}, start_BC_Y{}, end_BC_Y{}, start_BC_Z{}, end_BC_Z{};
        INT32 tag;
        std::vector<INT32> pop_to, pop_reverse;
        std::string BC_type;
        double velocity_along_x{0};
        INT32 surf_normal_X{}, surf_normal_Y{},surf_normal_Z{};

        pop_to.clear();
        pop_reverse.clear();
        BC_type = "";

        switch(ii)
        {
        case 0:
            /// Right BC
            if (end_global_x[this_processor_id] == this->Nx - 1)
            {
                pop_to      = pop_to_right_F;
                pop_reverse = pop_reverse_right_F;
                start_BC_X  = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                //cout << "Right, processor:" << processor_id << " "<<start_BC_X<< "-"<<end_BC_X << ", "<< start_BC_Y<< "-"<<end_BC_Y << ", "<< start_BC_Z<< "-"<<end_BC_Z << endl;
                if (X_Wall)
                {
                    BC_type = "wall";
                    velocity_along_x = 0.0;
                }
                if (X_velIn_ZeroGradOut)
                {
                    BC_type = "zero_grad";
                    surf_normal_X = -1;
                    surf_normal_Y = 0;
                    surf_normal_Z = 0;
                }
                if (X_velIn_PressureOut)
                {
                    BC_type = "press_out";
                }
                if (X_PressureIn_PressureOut)
                {
                    BC_type = "press_out";
                }
            }
        break;
        case 1:
            /// Left BC
            if (start_global_x[this_processor_id] == 0)
            {
                pop_to      = pop_to_left_F;
                pop_reverse = pop_reverse_left_F;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
                end_BC_X    = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                if (X_Wall)
                {
                    BC_type = "wall";
                    velocity_along_x = 0.0;
                }
                if (X_velIn_ZeroGradOut){
                    BC_type = "velocity";
                    velocity_along_x = inletvel;
                }
                if (X_velIn_PressureOut){
                    BC_type = "velocity";
                    velocity_along_x = inletvel;
                }
                if (X_PressureIn_PressureOut){
                    BC_type = "press_in";
                }
            }
        break;
        case 2:
            /// Top BC
            if (end_global_y[this_processor_id] == this->Ny - 1)
            {
                pop_to      = pop_to_top_F;
                pop_reverse = pop_reverse_top_F;
                start_BC_Y  = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                if (Y_Wall)
                {
                    BC_type = "wall";
                    velocity_along_x = wall_vel_top;
                }
            }
            break;
        case 3:
            /// Bottom BC
            if (start_global_y[this_processor_id] == 0)
            {
                pop_to      = pop_to_bottom_F;
                pop_reverse = pop_reverse_bottom_F;
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
                end_BC_Y    = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                if (Y_Wall)
                {
                    BC_type = "wall";
                    velocity_along_x = wall_vel_bottom;
                }
            }
            break;
        case 4:
            /// Front BC
            if (end_global_z[this_processor_id] == this->Nz - 1)
            {
                pop_to      = pop_to_front_F;
                pop_reverse = pop_reverse_front_F;
                start_BC_Z  = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                tag = ii;
                if (Z_Wall)
                {
                    BC_type = "wall";
                    velocity_along_x = wall_vel_front;
                }
            }
            tag = ii;
            break;
        case 5:
            /// Back BC
            if (start_global_z[this_processor_id] == 0)
            {
                pop_to      = pop_to_back_F;
                pop_reverse = pop_reverse_back_F;
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
                end_BC_Z    = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                tag = ii;
                if (Z_Wall)
                {
                    BC_type = "wall";
                    velocity_along_x = wall_vel_back;
                }
            }
            tag = ii;
            break;
        default:
            std::cout << "Invalid choice" << std::endl;
            break;
        }

        for (INT32 X = start_BC_X; X <= end_BC_X; ++X)
        {
            for (INT32 Y = start_BC_Y; Y <= end_BC_Y; ++Y)
            {
                for (INT32 Z = start_BC_Z; Z <= end_BC_Z; ++Z)
                {
                    //=========================== WALL ===========================//
                    if (BC_type == "wall")
                    {
                        Apply_MovingWall_BC(pop_to, pop_reverse, X, Y, Z, velocity_along_x, 0.0, 0.0);
                    }
                    //========================= VELOCITY =========================//
                    if (BC_type == "velocity") // Zou-He method for velocity inlet
                    {
                        if (stencil_flow == "D2Q9")
                        {
                            //Apply_Velocity_Inlet_Left_ZouHe_D2Q9(paralLB, X, Y, Z, velocity_along_x, 0.0);
                            Apply_Velocity_Inlet_Left_Ginzburg(paralLB, pop_to, pop_reverse, X, Y, Z, velocity_along_x, 0.0, 0.0);
                            //Apply_MovingWall_BC(pop_to, pop_reverse, X, Y, Z, velocity_along_x, 0.0, 0.0);// Not good for velocity inlet, better for moving wall
                        }
                    }
                    //====================== ZERO-GRADIENT =======================//
                    if (BC_type == "zero_grad")
                    {
                        //Apply_ZeroGradient_Right_BC(X, Y, Z, pop_reverse, pop_to, surf_normal_X, surf_normal_Y, surf_normal_Z);
                        Apply_NonEquilibrium_Extrapolation_Right_BC(X, Y, Z, pop_reverse, pop_to, surf_normal_X, surf_normal_Y, surf_normal_Z);// Better than first order zero gradient.
                        //Apply_Convective_Outlet_Right_BC(X, Y, Z, pop_reverse, pop_to, surf_normal_X, surf_normal_Y, surf_normal_Z);
                    }
                    //======================== PRESSURE-IN =======================//
                    if (BC_type == "press_in")
                    {
                        double rho_in = 1.0002; // Example density for pressure inlet
                        Apply_Pressure_Inlet_Left_Ginzburg(pop_to, pop_reverse, X, Y, Z, rho_in);
                        //Apply_Pressure_Inlet_Left_ZouHe_D2Q9(X, Y, Z, rho_in);
                    }
                    //======================= PRESSURE-OUT =======================//
                    if (BC_type == "press_out")
                    {
                        double rho_out = 1.0;
                        Apply_Pressure_Outlet_Right_Ginzburg(pop_to, pop_reverse, X, Y, Z, surf_normal_X, surf_normal_Y, surf_normal_Z, rho_out);// Ginzburg, Verhaeghe, Communications in Computational Physics (2008) -See Eq (5.53) of Krueger LBM book
                        //Apply_Pressure_Outlet_Right_ZouHe_D2Q9(X, Y, Z, rho_out);
                    }
                    //============================================================//
                }//Z
            }//Y
        }//X
    }// End of: for (INT32 ii = 0; ii < 6; ii++)

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_MovingWall_BC
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_MovingWall_BC(const std::vector<INT32>& pop_to, const std::vector<INT32>& pop_reverse, INT32 X, INT32 Y, INT32 Z, double wall_ux, double wall_uy, double wall_uz)
{
    #if STORE_1D_VECTOR
    INT32 id3 = idx_3d_vector(X, Y, Z);

    double rho = density[id3];

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor fluid node
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        INT32 id_boundary = idx_4d_vector(pop_rev, X, Y, Z);

        INT32 id_fluid = idx_4d_vector(pop, newX, newY, newZ);

        // c_i . u_wall
        double cu =
            cx_pop_F[pop] * wall_ux
            +
            cy_pop_F[pop] * wall_uy
            +
            cz_pop_F[pop] * wall_uz;

        // Moving wall bounce-back
        pop_F[id_boundary] = pop_F[id_fluid] -  (2.0 / cs_2_Flow) * weight_F[pop] * rho * cu;
    }
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
        double rho = density[X][Y][Z];

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor fluid node
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        // c_i . u_wall
        double cu =
            cx_pop_F[pop] * wall_ux
            +
            cy_pop_F[pop] * wall_uy
            +
            cz_pop_F[pop] * wall_uz;

        // Moving wall bounce-back
        pop_F[pop_rev][X][Y][Z] =
            pop_F[pop][newX][newY][newZ]
            -
            (2.0 / cs_2_Flow)
            *
            weight_F[pop]
            *
            rho
            *
            cu;
    }
    #endif // !STORE_1D_VECTOR

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Velocity_Inlet_Left_ZouHe_D2Q9
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_Velocity_Inlet_Left_ZouHe_D2Q9(
    const ParallelMPI& paralLB,
    INT32 X,
    INT32 Y,
    INT32 Z,
    double ux_mean,
    double uy = 0.0)
{
    double ux;
    if (Constant_Velocity_Inlet)
    {
        ux = ux_mean;
    }else if (Parabolic_Velocity_Inlet)
    {
        double H = Ny - 1;
        double y_global = paralLB.get_Yglobal_from_Ylocal(Y, this_processor_id);
        ux = 6.0 * ux_mean * y_global * (H - y_global) / (H * H); //full developed parabolic profile for 2D velocity inlet, can be modified to uniform or other profiles as needed
    }

    #if STORE_1D_VECTOR
    // --- Known populations ---
    double f0 = pop_F[idx_4d_vector(0, X, Y, Z)];
    double f2 = pop_F[idx_4d_vector(2, X, Y, Z)];
    double f3 = pop_F[idx_4d_vector(3, X, Y, Z)];
    double f4 = pop_F[idx_4d_vector(4, X, Y, Z)];
    double f6 = pop_F[idx_4d_vector(6, X, Y, Z)];
    double f7 = pop_F[idx_4d_vector(7, X, Y, Z)];

    // --- Density from Zou-He constraint ---
    double rho = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - ux);

    // --- Unknown populations ---
    // f1
    pop_F[idx_4d_vector(1, X, Y, Z)] = f3 + (2.0 / 3.0) * rho * ux;

    // f5
    pop_F[idx_4d_vector(5, X, Y, Z)] = f7 - 0.5 * (f2 - f4) + (1.0 / 6.0) * rho * ux + 0.5 * rho * uy;

    // f8
    pop_F[idx_4d_vector(8, X, Y, Z)] = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * rho * ux - 0.5 * rho * uy;
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
        // --- Known populations ---
    double f0 = pop_F[0][X][Y][Z];
    double f2 = pop_F[2][X][Y][Z];
    double f3 = pop_F[3][X][Y][Z];
    double f4 = pop_F[4][X][Y][Z];
    double f6 = pop_F[6][X][Y][Z];
    double f7 = pop_F[7][X][Y][Z];

    // --- Density from Zou-He constraint ---
    double rho =
        (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7))
        /
        (1.0 - ux);

    // --- Unknown populations ---

    // f1
    pop_F[1][X][Y][Z] =
        f3
        +
        (2.0 / 3.0) * rho * ux;

    // f5
    pop_F[5][X][Y][Z] =
        f7
        -
        0.5 * (f2 - f4)
        +
        (1.0 / 6.0) * rho * ux
        +
        0.5 * rho * uy;

    // f8
    pop_F[8][X][Y][Z] =
        f6
        +
        0.5 * (f2 - f4)
        +
        (1.0 / 6.0) * rho * ux
        -
        0.5 * rho * uy;
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Velocity_Inlet_Left_Ginzburg
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_Velocity_Inlet_Left_Ginzburg(
    const ParallelMPI& paralLB,
    const std::vector<INT32>& pop_to,
    const std::vector<INT32>& pop_reverse,
    INT32 X,
    INT32 Y,
    INT32 Z,
    double ux_mean,
    double uy,
    double uz)
{
    double ux;
    if (Constant_Velocity_Inlet)
    {
        ux = ux_mean;
    }else if (Parabolic_Velocity_Inlet)
    {
        double H = Ny - 1;
        double y_global = paralLB.get_Yglobal_from_Ylocal(Y, this_processor_id);
        ux = 6.0 * ux_mean * y_global * (H - y_global) / (H * H); //full developed parabolic profile for 2D velocity inlet, can be modified to uniform or other profiles as needed
    }

    #if STORE_1D_VECTOR
    // Density extrapolated from adjacent fluid node
    double rho = density[idx_3d_vector(X + 1, Y, Z)];

    // u^2
    double u2 =
        ux * ux
        +
        uy * uy
        +
        uz * uz;

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor fluid node
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        // Population indices
        INT32 id_boundary =
            idx_4d_vector(pop_rev, X, Y, Z);

        INT32 id_fluid =
            idx_4d_vector(pop, newX, newY, newZ);

        // e_i . u
        double eu =
            cx_pop_F[pop] * ux
            +
            cy_pop_F[pop] * uy
            +
            cz_pop_F[pop] * uz;

        // 2*f_eq
        double feq =
            2.0
            *
            weight_F[pop]
            *
            rho
            *
            (
                1.0
                +
                eu * eu / (2.0 * cs_2_Flow * cs_2_Flow)
                -
                u2 / (2.0 * cs_2_Flow)
            );

        // Ginzburg-Verhaeghe velocity BC
        pop_F[id_boundary] =
            -pop_F[id_fluid]
            +
            feq;
    }

    #endif // STORE_1D_VECTOR

    //==================================================================//

    #if !STORE_1D_VECTOR

    // Density extrapolated from adjacent fluid node
    double rho =
        density[X + 1][Y][Z];

    // u^2
    double u2 =
        ux * ux
        +
        uy * uy
        +
        uz * uz;

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor fluid node
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        // e_i . u
        double eu =
            cx_pop_F[pop] * ux
            +
            cy_pop_F[pop] * uy
            +
            cz_pop_F[pop] * uz;

        // 2*f_eq
        double feq =
            2.0
            *
            weight_F[pop]
            *
            rho
            *
            (
                1.0
                +
                eu * eu / (2.0 * cs_2_Flow * cs_2_Flow)
                -
                u2 / (2.0 * cs_2_Flow)
            );

        // Ginzburg-Verhaeghe velocity BC
        pop_F[pop_rev][X][Y][Z] =
            -pop_F[pop][newX][newY][newZ]
            +
            feq;
    }

    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_ZeroGradient_Right_BC
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_ZeroGradient_Right_BC(INT32 X, INT32 Y, INT32 Z, const std::vector<INT32>& pop_reverse, const std::vector<INT32>& pop_to, INT32 normalX, INT32 normalY, INT32 normalZ)
{
    #if STORE_1D_VECTOR
    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop_rev = pop_reverse[alpha];

        INT32 id_boundary = idx_4d_vector(pop_rev, X, Y, Z);

        INT32 id_fluid = idx_4d_vector(pop_rev, X + normalX, Y + normalY, Z + normalZ);

        // First-order extrapolation
        pop_F[id_boundary] = pop_F[id_fluid];

        // Optional second-order:
        /*
        INT32 id_fluid2 =
            idx_4d_vector(
                pop_rev,
                X + 2 * normalX,
                Y + 2 * normalY,
                Z + 2 * normalZ);

        pop_F[id_boundary] =
            2.0 * pop_F[id_fluid]
            - pop_F[id_fluid2];
        */
    }
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop_rev = pop_reverse[alpha];

        // First-order extrapolation
        pop_F[pop_rev][X][Y][Z] =
            pop_F[pop_rev]
                 [X + normalX]
                 [Y + normalY]
                 [Z + normalZ];

        // Optional second-order:
        /*
        pop_F[pop_rev][X][Y][Z] =
            2.0 *
            pop_F[pop_rev]
                 [X + normalX]
                 [Y + normalY]
                 [Z + normalZ]
            -
            pop_F[pop_rev]
                 [X + 2 * normalX]
                 [Y + 2 * normalY]
                 [Z + 2 * normalZ];
        */
    }
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_NonEquilibrium_Extrapolation_Right_BC
//-----------------------------------------------------------------------------------------------------------------------------
//See Z. L. Guo, C. G. Zheng, B. C. Shi, "Non-equilibrium extrapolation method for velocity and pressure boundary conditions in the lattice Boltzmann method"
///Journal: Chinese Physics, Volume 11, Number 4, pp. 366–374 (2002)
void FlowSolver::Apply_NonEquilibrium_Extrapolation_Right_BC(INT32 X, INT32 Y, INT32 Z, const std::vector<INT32>& pop_reverse, const std::vector<INT32>& pop_to, INT32 normalX, INT32 normalY, INT32 normalZ)
{
    // Fluid node adjacent to boundary
    INT32 Xf = X + normalX;
    INT32 Yf = Y + normalY;
    INT32 Zf = Z + normalZ;

    #if STORE_1D_VECTOR
    // Boundary macroscopic quantities
    double rho_b = 1;//density[idx_3d_vector(X, Y, Z)];

    double ux_b = velocity[idx_4d_vector(0, Xf, Yf, Zf)];

    double uy_b = velocity[idx_4d_vector(1, Xf, Yf, Zf)];

    double uz_b = velocity[idx_4d_vector(2, Xf, Yf, Zf)];

    // Interior fluid macroscopic quantities
    double rho_f = density[idx_3d_vector(Xf, Yf, Zf)];

    double ux_f = velocity[idx_4d_vector(0, Xf, Yf, Zf)];

    double uy_f = velocity[idx_4d_vector(1, Xf, Yf, Zf)];

    double uz_f = velocity[idx_4d_vector(2, Xf, Yf, Zf)];

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop_rev = pop_reverse[alpha];

        INT32 id_boundary = idx_4d_vector(pop_rev, X, Y, Z);
        INT32 id_fluid    = idx_4d_vector(pop_rev, Xf, Yf, Zf);

        // Equilibrium at boundary node
        double feq_b = equilibrium_flow(rho_b, ux_b, uy_b, uz_b, pop_rev);

        // Equilibrium at adjacent fluid node
        double feq_f = equilibrium_flow(rho_f, ux_f, uy_f, uz_f, pop_rev);

        // NNEQ extrapolation
        pop_F[id_boundary] = feq_b + (pop_F[id_fluid] - feq_f);
    }
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
    // Boundary macroscopic quantities
    double rho_b = 1;//density[X][Y][Z];

    double ux_b = velocity[0][Xf][Yf][Zf];

    double uy_b = velocity[1][Xf][Yf][Zf];

    double uz_b = velocity[2][Xf][Yf][Zf];

    // Interior fluid macroscopic quantities
    double rho_f = density[Xf][Yf][Zf];

    double ux_f = velocity[0][Xf][Yf][Zf];

    double uy_f = velocity[1][Xf][Yf][Zf];

    double uz_f = velocity[2][Xf][Yf][Zf];

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop_rev = pop_reverse[alpha];

        // Equilibrium at boundary node
        double feq_b =
            equilibrium_flow(
                rho_b,
                ux_b,
                uy_b,
                uz_b,
                pop_rev);

        // Equilibrium at adjacent fluid node
        double feq_f =
            equilibrium_flow(
                rho_f,
                ux_f,
                uy_f,
                uz_f,
                pop_rev);

        // NNEQ extrapolation
        pop_F[pop_rev][X][Y][Z] =
            feq_b
            +
            (
                pop_F[pop_rev][Xf][Yf][Zf]
                -
                feq_f
            );
    }    
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Convective_Outlet_Right_BC
//-----------------------------------------------------------------------------------------------------------------------------
//Z. Yang, T. Tang, et al. "Lattice Boltzmann outflow treatments: Convective condition and others"
//Journal: Computers & Mathematics with Applications, Volume 65, Issue 10, Pages 1604–1616 (2013)
void FlowSolver::Apply_Convective_Outlet_Right_BC(INT32 X, INT32 Y, INT32 Z, 
    const std::vector<INT32>& pop_reverse, const std::vector<INT32>& pop_to, 
    INT32 normalX, INT32 normalY, INT32 normalZ)
{
    // Fluid node adjacent to boundary (inside domain)
    INT32 Xf = X + normalX;
    INT32 Yf = Y + normalY;
    INT32 Zf = Z + normalZ;

    // Convective velocity (use previous time step values!)
    double u_conv = 0.0;

#if STORE_1D_VECTOR
    u_conv = velocity[idx_4d_vector(0, Xf, Yf, Zf)];           // ux at interior node
    // Optional: more stable average
    // u_conv = (velocity[idx_4d_vector(0, Xf, Yf, Zf)] + velocity[idx_4d_vector(0, Xf-1, Yf, Zf)]) * 0.5;

#else
    u_conv = velocity[0][Xf][Yf][Zf];   // ux at interior node (previous step)
    // Optional average:
    // u_conv = (velocity[0][Xf][Yf][Zf] + velocity[0][Xf-1][Yf][Zf]) * 0.5;
#endif

    // Clamp convective velocity for stability (optional but recommended)
    //if (u_conv < 0.0) u_conv = 0.0;
    //if (u_conv > 0.1) u_conv = 0.1;   // safety for your low Mach case

    // ==============================================================
#if STORE_1D_VECTOR
    for (INT32 alpha = 0; alpha < stencil_directions_F; alpha++)   // Apply to ALL populations (important!)
    {
        INT32 id_boundary = idx_4d_vector(alpha, X, Y, Z);
        INT32 id_fluid    = idx_4d_vector(alpha, Xf, Yf, Zf);

        // Simple first-order convective BC
        pop_F[id_boundary] = pop_F[id_boundary] - u_conv * (pop_F[id_boundary] - pop_F[id_fluid]);
    }
#endif // STORE_1D_VECTOR
    // ==============================================================
#if !STORE_1D_VECTOR
    for (INT32 alpha = 0; alpha < stencil_directions_F; alpha++)   // Apply to ALL populations
    {
        // Simple first-order convective BC
        pop_F[alpha][X][Y][Z] = pop_F[alpha][X][Y][Z] - u_conv * (pop_F[alpha][X][Y][Z] - pop_F[alpha][Xf][Yf][Zf]);
    }
#endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Pressure_Inlet_Left_Ginzburg
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_Pressure_Inlet_Left_Ginzburg(const std::vector<INT32>& pop_to, const std::vector<INT32>& pop_reverse, INT32 X, INT32 Y, INT32 Z, double rho_in)
{
    #if STORE_1D_VECTOR
    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor fluid node
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        // Population indices
        INT32 id_boundary =
            idx_4d_vector(pop_rev, X, Y, Z);

        INT32 id_fluid =
            idx_4d_vector(pop, newX, newY, newZ);

        // Velocity at boundary node
        // Krüger Eq. (5.53)
        double ux = velocity[idx_4d_vector(0, X, Y, Z)];
        double uy = 0.0;
        double uz = 0.0;

        // e_i . u
        double eu =
            cx_pop_F[pop] * ux
            +
            cy_pop_F[pop] * uy
            +
            cz_pop_F[pop] * uz;

        // u^2
        double u2 =
            ux * ux
            +
            uy * uy
            +
            uz * uz;

        // Equilibrium contribution
        double feq =
            2.0
            *
            weight_F[pop]
            *
            rho_in
            *
            (
                1.0
                +
                eu * eu / (2.0 * cs_2_Flow * cs_2_Flow)
                -
                u2 / (2.0 * cs_2_Flow)
            );

        // Ginzburg-Verhaeghe pressure BC
        pop_F[id_boundary] =
            -pop_F[id_fluid]
            +
            feq;
    }
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
    // Velocity at boundary node
    // Krüger Eq. (5.53)
    double ux = velocity[0][X][Y][Z];
    double uy = 0.0;
    double uz = 0.0;

    // u^2
    double u2 =
        ux * ux
        +
        uy * uy
        +
        uz * uz;

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor fluid node
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        // e_i . u
        double eu =
            cx_pop_F[pop] * ux
            +
            cy_pop_F[pop] * uy
            +
            cz_pop_F[pop] * uz;

        // Equilibrium contribution
        double feq =
            2.0
            *
            weight_F[pop]
            *
            rho_in
            *
            (
                1.0
                +
                eu * eu / (2.0 * cs_2_Flow * cs_2_Flow)
                -
                u2 / (2.0 * cs_2_Flow)
            );

        // Ginzburg-Verhaeghe pressure BC
        pop_F[pop_rev][X][Y][Z] =
            -pop_F[pop][newX][newY][newZ]
            +
            feq;
    }
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Pressure_Inlet_Left_ZouHe_D2Q9
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_Pressure_Inlet_Left_ZouHe_D2Q9(
    INT32 X,
    INT32 Y,
    INT32 Z,
    double rho_in)
{
    #if STORE_1D_VECTOR
    // --- Known populations ---
    double f0 = pop_F[idx_4d_vector(0, X, Y, Z)];
    double f2 = pop_F[idx_4d_vector(2, X, Y, Z)];
    double f3 = pop_F[idx_4d_vector(3, X, Y, Z)];
    double f4 = pop_F[idx_4d_vector(4, X, Y, Z)];
    double f6 = pop_F[idx_4d_vector(6, X, Y, Z)];
    double f7 = pop_F[idx_4d_vector(7, X, Y, Z)];

    // --- Velocity from Zou-He constraint ---
    double ux =
        1.0
        -
        (
            f0 + f2 + f4
            +
            2.0 * (f3 + f6 + f7)
        ) / rho_in;

    double uy = 0.0;

    // --- Unknown populations ---

    // f1
    pop_F[idx_4d_vector(1, X, Y, Z)] =
        f3
        +
        (2.0 / 3.0) * rho_in * ux;

    // f5
    pop_F[idx_4d_vector(5, X, Y, Z)] =
        f7
        -
        0.5 * (f2 - f4)
        +
        (1.0 / 6.0) * rho_in * ux
        +
        0.5 * rho_in * uy;

    // f8
    pop_F[idx_4d_vector(8, X, Y, Z)] =
        f6
        +
        0.5 * (f2 - f4)
        +
        (1.0 / 6.0) * rho_in * ux
        -
        0.5 * rho_in * uy;
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
    // --- Known populations ---
    double f0 = pop_F[0][X][Y][Z];
    double f2 = pop_F[2][X][Y][Z];
    double f3 = pop_F[3][X][Y][Z];
    double f4 = pop_F[4][X][Y][Z];
    double f6 = pop_F[6][X][Y][Z];
    double f7 = pop_F[7][X][Y][Z];

    // --- Velocity from Zou-He constraint ---
    double ux =
        1.0
        -
        (
            f0 + f2 + f4
            +
            2.0 * (f3 + f6 + f7)
        ) / rho_in;

    double uy = 0.0;

    // --- Unknown populations ---

    // f1
    pop_F[1][X][Y][Z] =
        f3
        +
        (2.0 / 3.0) * rho_in * ux;

    // f5
    pop_F[5][X][Y][Z] =
        f7
        -
        0.5 * (f2 - f4)
        +
        (1.0 / 6.0) * rho_in * ux
        +
        0.5 * rho_in * uy;

    // f8
    pop_F[8][X][Y][Z] =
        f6
        +
        0.5 * (f2 - f4)
        +
        (1.0 / 6.0) * rho_in * ux
        -
        0.5 * rho_in * uy;
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Pressure_Outlet_Right_Ginzburg
// Ginzburg, Verhaeghe, Communications in Computational Physics (2008) -See Eq (5.53) of Krueger LBM book
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_Pressure_Outlet_Right_Ginzburg(
    const std::vector<INT32>& pop_to,
    const std::vector<INT32>& pop_reverse,
    INT32 X,
    INT32 Y,
    INT32 Z,
    INT32 normalX,
    INT32 normalY,
    INT32 normalZ,
    double rho_out)
{
    // Adjacent interior fluid node
    INT32 Xf = X + normalX;
    INT32 Yf = Y + normalY;
    INT32 Zf = Z + normalZ;

    #if STORE_1D_VECTOR
    // Velocity extrapolated from fluid node
    double ux = velocity[idx_4d_vector(0, Xf, Yf, Zf)];
    double uy = 0.0;
    double uz = 0.0;

    double u2 = ux * ux + uy * uy + uz * uz;

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor node in streaming direction
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        INT32 id_boundary =
            idx_4d_vector(pop_rev, X, Y, Z);

        INT32 id_fluid =
            idx_4d_vector(pop, newX, newY, newZ);

        // e_i . u
        double eu =
            cx_pop_F[pop] * ux
            +
            cy_pop_F[pop] * uy
            +
            cz_pop_F[pop] * uz;

        // 2*f_eq
        double feq = 2.0 * weight_F[pop] * rho_out * (1.0 + eu * eu / (2.0 * cs_2_Flow * cs_2_Flow) - u2 / (2.0 * cs_2_Flow));

        // Ginzburg-Verhaeghe BC
        pop_F[id_boundary] = -pop_F[id_fluid] + feq;
    }
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
    // Velocity extrapolated from fluid node
    double ux = velocity[0][Xf][Yf][Zf];
    double uy = 0.0;
    double uz = 0.0;

    double u2 =
        ux * ux
        +
        uy * uy
        +
        uz * uz;

    for (INT32 alpha = 0; alpha < pop_to.size(); alpha++)
    {
        INT32 pop     = pop_to[alpha];
        INT32 pop_rev = pop_reverse[alpha];

        // Neighbor node in streaming direction
        INT32 newX = X + cx_pop_F[pop];
        INT32 newY = Y + cy_pop_F[pop];
        INT32 newZ = Z + cz_pop_F[pop];

        // e_i . u
        double eu =
            cx_pop_F[pop] * ux
            +
            cy_pop_F[pop] * uy
            +
            cz_pop_F[pop] * uz;

        // 2*f_eq
        double feq =
            2.0
            *
            weight_F[pop]
            *
            rho_out
            *
            (
                1.0
                +
                eu * eu / (2.0 * cs_2_Flow * cs_2_Flow)
                -
                u2 / (2.0 * cs_2_Flow)
            );

        // Ginzburg-Verhaeghe BC
        pop_F[pop_rev][X][Y][Z] =
            -pop_F[pop][newX][newY][newZ]
            +
            feq;
    }
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Apply_Pressure_Outlet_Right_ZouHe_D2Q9
// Zou, He, Journal of Computational Physics (1997) -See Eq (4.10) of Krueger LBM book
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Apply_Pressure_Outlet_Right_ZouHe_D2Q9(
    INT32 X,
    INT32 Y,
    INT32 Z,
    double rho_out)
{
    #if STORE_1D_VECTOR
    // --- Known populations ---
    double f0 = pop_F[idx_4d_vector(0, X, Y, Z)];
    double f1 = pop_F[idx_4d_vector(1, X, Y, Z)];
    double f2 = pop_F[idx_4d_vector(2, X, Y, Z)];
    double f4 = pop_F[idx_4d_vector(4, X, Y, Z)];
    double f5 = pop_F[idx_4d_vector(5, X, Y, Z)];
    double f8 = pop_F[idx_4d_vector(8, X, Y, Z)];

    // --- Velocity from Zou-He constraint ---
    double ux =
        -1.0
        +
        (
            f0 + f2 + f4
            +
            2.0 * (f1 + f5 + f8)
        ) / rho_out;

    double uy = 0.0;

    // --- Unknown populations ---

    // f3
    pop_F[idx_4d_vector(3, X, Y, Z)] =
        f1
        -
        (2.0 / 3.0) * rho_out * ux;

    // f6
    pop_F[idx_4d_vector(6, X, Y, Z)] =
        f8
        +
        0.5 * (f4 - f2)
        -
        (1.0 / 6.0) * rho_out * ux
        -
        0.5 * rho_out * uy;

    // f7
    pop_F[idx_4d_vector(7, X, Y, Z)] =
        f5
        +
        0.5 * (f2 - f4)
        -
        (1.0 / 6.0) * rho_out * ux
        +
        0.5 * rho_out * uy;
    #endif // STORE_1D_VECTOR
    //==================================================================//
    #if !STORE_1D_VECTOR
    // --- Known populations ---
    double f0 = pop_F[0][X][Y][Z];
    double f1 = pop_F[1][X][Y][Z];
    double f2 = pop_F[2][X][Y][Z];
    double f4 = pop_F[4][X][Y][Z];
    double f5 = pop_F[5][X][Y][Z];
    double f8 = pop_F[8][X][Y][Z];

    // --- Velocity from Zou-He constraint ---
    double ux =
        -1.0
        +
        (
            f0 + f2 + f4
            +
            2.0 * (f1 + f5 + f8)
        ) / rho_out;

    double uy = 0.0;

    // --- Unknown populations ---

    // f3
    pop_F[3][X][Y][Z] =
        f1
        -
        (2.0 / 3.0) * rho_out * ux;

    // f6
    pop_F[6][X][Y][Z] =
        f8
        +
        0.5 * (f4 - f2)
        -
        (1.0 / 6.0) * rho_out * ux
        -
        0.5 * rho_out * uy;

    // f7
    pop_F[7][X][Y][Z] =
        f5
        +
        0.5 * (f2 - f4)
        -
        (1.0 / 6.0) * rho_out * ux
        +
        0.5 * rho_out * uy;
    #endif // !STORE_1D_VECTOR
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: bc_at_solid_points_flow
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::bc_at_solid_points_flow()
{
        ///Now for is_solid BC
    if (!Flow_With_Particle){
        const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
        const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
        const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
        const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
        const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
        const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;
        for (INT32 X = start_x; X <= end_x; ++X) {
            for (INT32 Y = start_y; Y <= end_y; ++Y) {
                for (INT32 Z = start_z; Z <= end_z; ++Z) {
                    #if STORE_1D_VECTOR
                    INT32 id3 = idx_3d_vector(X, Y, Z);
                    if (is_solid[id3] == FALSE){
                        for (INT32 dir = 1; dir < stencil_directions_F; dir++)// dire=0 is not needed
                        {
                            const INT32 Xp = (X - cx_pop_F[dir]);
                            const INT32 Yp = (Y - cy_pop_F[dir]);
                            const INT32 Zp = (Z - cz_pop_F[dir]);
                            INT32 id3_p = idx_3d_vector(Xp, Yp, Zp);
                            INT32 id4   = idx_4d_vector(dir, X, Y, Z);
                            INT32 id4_p = idx_4d_vector(bounce_back_pop_F[dir], Xp, Yp, Zp);
                            if (is_solid[id3_p] == TRUE)
                            {
                                pop_F[id4] = pop_F[id4_p];
                            }
                        }
                    }
                    #endif
                    #if !STORE_1D_VECTOR
                    if (is_solid[X][Y][Z] == FALSE){
                        for (INT32 dir = 1; dir < stencil_directions_F; dir++)// dire=0 is not needed
                        {
                            const INT32 Xp = (X - cx_pop_F[dir]);
                            const INT32 Yp = (Y - cy_pop_F[dir]);
                            const INT32 Zp = (Z - cz_pop_F[dir]);
                            if (is_solid[Xp][Yp][Zp] == TRUE)
                            {
                                pop_F[dir][X][Y][Z] = pop_F[bounce_back_pop_F[dir]][Xp][Yp][Zp];
                            }
                        }
                    }
                    #endif
                }
            }
        }
    }


   /* const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;

    INT32 newX, newY, newZ;

    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                INT32 id3      = idx_3d_vector(X, Y, Z);
                #if !STORE_1D_VECTOR
                if (is_solid[X][Y][Z])
                    #endif
                #if STORE_1D_VECTOR
                if (is_solid[id3])
                #endif // STORE_1D_VECTOR
                {
                    for (INT32 alpha = 0; alpha < stencil_directions_F; alpha++)
                    {
                        newX = X - cx_pop_F[alpha];
                        newY = Y - cy_pop_F[alpha];
                        newZ = Z - cz_pop_F[alpha];
                        INT32 id3_new = idx_3d_vector(newX, newY, newZ);
                        #if !STORE_1D_VECTOR
                        if (!is_solid[newX][newY][newZ])
                        #endif
                        #if STORE_1D_VECTOR
                        if (!is_solid[id3_new])
                        #endif
                        {
                            #if !STORE_1D_VECTOR
                            pop_F[bounce_back_pop_F[alpha]][X][Y][Z] = pop_F[alpha][X][Y][Z];//Full-way Bounce Back
                            #endif
                            #if STORE_1D_VECTOR
                            INT32 id4      = idx_4d_vector(alpha                   , X, Y, Z);
                            INT32 id4_back = idx_4d_vector(bounce_back_pop_F[alpha], X, Y, Z);
                            pop_F[id4_back] = pop_F[id4];//Full-way Bounce Back
                            #endif
                        }
                    }
                }
            }
        }
    }*/

    /*
    INT32 start_X, end_X;//Start and End of LBM calculation for processor.
    start_X = 0;
    end_X = Nx - 1;

    if (X_Wall)
    {
        start_X = 1;
        end_X = Nx - 2;
    }
    if (X_velIn_ZeroGradOut)
    {
        start_X = 1;
        end_X = Nx - 2;
    }
    if (ZeroGradIn_ZeroGradOut)
    {
        start_X = 1;
        end_X = Nx - 2;
    }

    INT32 start_Y, end_Y;
    start_Y = 0;
    end_Y = Ny - 1;

    if (Y_Wall)
    {
        start_Y = 1;
        end_Y = Ny - 2;
    }

    INT32 start_Z, end_Z;
    start_Z = 0;
    end_Z = Nz - 1;

    if (Z_Wall)
    {
        start_Z = 1;
        end_Z = Nz - 2;
    }

    // Bounce back on SOLID nodes (Just for !Flow_With_Particle, because I have not defined collision between particles and solid objects.)
    if (!Flow_With_Particle)
    {
        for (INT32 X = start_X; X <= end_X; ++X)
        {
            for (INT32 Y = start_Y; Y <= end_Y; ++Y)
            {
                for (INT32 Z = start_Z; Z <= end_Z; ++Z)
                {
                    if (!is_solid[X][Y][Z])
                    {

                        if (is_solid[(X - 1 + Nx) % Nx][Y][Z])
                        {
                            pop_F[1][X][Y][Z] = pop_F[2][(X - 1 + Nx) % Nx][Y][Z];
                        }
                        if (is_solid[(X + 1) % Nx][Y][Z])
                        {
                            pop_F[2][X][Y][Z] = pop_F[1][(X + 1) % Nx][Y][Z];
                        }
                        if (is_solid[X][(Y - 1 + Ny) % Ny][Z])
                        {
                            pop_F[3][X][Y][Z] = pop_F[4][X][(Y - 1 + Ny) % Ny][Z];
                        }
                        if (is_solid[X][(Y + 1) % Ny][Z])
                        {
                            pop_F[4][X][Y][Z] = pop_F[3][X][(Y + 1) % Ny][Z];
                        }
                        if (is_solid[X][Y][(Z - 1 + Nz) % Nz])
                        {
                            pop_F[5][X][Y][Z] = pop_F[6][X][Y][(Z - 1 + Nz) % Nz];
                        }
                        if (is_solid[X][Y][(Z + 1) % Nz])
                        {
                            pop_F[6][X][Y][Z] = pop_F[5][X][Y][(Z + 1) % Nz];
                        }
                        if (is_solid[(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z])
                        {
                            pop_F[7][X][Y][Z] = pop_F[10][(X - 1 + Nx) % Nx][(Y - 1 + Ny) % Ny][Z];
                        }
                        if (is_solid[(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z])
                        {
                            pop_F[8][X][Y][Z] = pop_F[9][(X + 1) % Nx][(Y - 1 + Ny) % Ny][Z];
                        }
                        if (is_solid[(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z])
                        {
                            pop_F[9][X][Y][Z] = pop_F[8][(X - 1 + Nx) % Nx][(Y + 1) % Ny][Z];
                        }
                        if (is_solid[(X + 1) % Nx][(Y + 1) % Ny][Z])
                        {
                            pop_F[10][X][Y][Z] = pop_F[7][(X + 1) % Nx][(Y + 1) % Ny][Z];
                        }
                        if (is_solid[(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz])
                        {
                            pop_F[11][X][Y][Z] = pop_F[14][(X - 1 + Nx) % Nx][Y][(Z - 1 + Nz) % Nz];
                        }
                        if (is_solid[(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz])
                        {
                            pop_F[12][X][Y][Z] = pop_F[13][(X + 1) % Nx][Y][(Z - 1 + Nz) % Nz];
                        }
                        if (is_solid[(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz])
                        {
                            pop_F[13][X][Y][Z] = pop_F[12][(X - 1 + Nx) % Nx][Y][(Z + 1) % Nz];
                        }
                        if (is_solid[(X + 1) % Nx][Y][(Z + 1) % Nz])
                        {
                            pop_F[14][X][Y][Z] = pop_F[11][(X + 1) % Nx][Y][(Z + 1) % Nz];
                        }
                        if (is_solid[X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz])
                        {
                            pop_F[15][X][Y][Z] = pop_F[18][X][(Y - 1 + Ny) % Ny][(Z - 1 + Nz) % Nz];
                        }
                        if (is_solid[X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz])
                        {
                            pop_F[16][X][Y][Z] = pop_F[17][X][(Y + 1) % Ny][(Z - 1 + Nz) % Nz];
                        }
                        if (is_solid[X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz])
                        {
                            pop_F[17][X][Y][Z] = pop_F[16][X][(Y - 1 + Ny) % Ny][(Z + 1) % Nz];
                        }
                        if (is_solid[X][(Y + 1) % Ny][(Z + 1) % Nz])
                        {
                            pop_F[18][X][Y][Z] = pop_F[15][X][(Y + 1) % Ny][(Z + 1) % Nz];
                        }
                    }
                }
            }
        }
    }*/

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: momenta_flow
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is not included here.
// It must be taken into account whenever the physical velocity is required.
void FlowSolver::momenta_flow()
{
    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    //apply_sponge_layer_on_outlet_X(30); //It is done before calculating momenta, because it is necessary to update the populations at the outlet before calculating momenta. Otherwise, the velocity and density at the outlet will not be updated according to the new populations after applying sponge layer.

    numerator_residual_flow   = 0;
    denominator_residual_flow = 0;

    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                #if !STORE_1D_VECTOR
                density[X][Y][Z] = 0;
                previous_velocity[0][X][Y][Z] = velocity[0][X][Y][Z];
                previous_velocity[1][X][Y][Z] = velocity[1][X][Y][Z];
                previous_velocity[2][X][Y][Z] = velocity[2][X][Y][Z];
                velocity[0][X][Y][Z] = 0;
                velocity[1][X][Y][Z] = 0;
                velocity[2][X][Y][Z] = 0;
                force[0][X][Y][Z] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
                force[1][X][Y][Z] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
                force[2][X][Y][Z] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
                if (!is_solid[X][Y][Z]){
                    for (INT32 i = 0; i < stencil_directions_F; i++)
                    {
                        double f = pop_F[i][X][Y][Z];
                        density[X][Y][Z]     += f;
                        velocity[0][X][Y][Z] += cx_pop_F[i] * f;
                        velocity[1][X][Y][Z] += cy_pop_F[i] * f;
                        velocity[2][X][Y][Z] += cz_pop_F[i] * f;
                    }
                    velocity[0][X][Y][Z] /= density[X][Y][Z];
                    velocity[1][X][Y][Z] /= density[X][Y][Z];
                    velocity[2][X][Y][Z] /= density[X][Y][Z];
                }else{
                    density[X][Y][Z]     = 1.0;
                    velocity[0][X][Y][Z] = 0.0;
                    velocity[1][X][Y][Z] = 0.0;
                    velocity[2][X][Y][Z] = 0.0;
                }
                numerator_residual_flow += pow(velocity[0][X][Y][Z] - previous_velocity[0][X][Y][Z], 2);
                numerator_residual_flow += pow(velocity[1][X][Y][Z] - previous_velocity[1][X][Y][Z], 2);
                numerator_residual_flow += pow(velocity[2][X][Y][Z] - previous_velocity[2][X][Y][Z], 2);
                denominator_residual_flow += pow(previous_velocity[0][X][Y][Z], 2);
                denominator_residual_flow += pow(previous_velocity[1][X][Y][Z], 2);
                denominator_residual_flow += pow(previous_velocity[2][X][Y][Z], 2);
                #endif
                #if STORE_1D_VECTOR
                INT32 id3 = idx_3d_vector(X, Y, Z);
                INT32 id4_0 = idx_4d_vector(0, X, Y, Z);
                INT32 id4_1 = idx_4d_vector(1, X, Y, Z);
                INT32 id4_2 = idx_4d_vector(2, X, Y, Z);
                density[id3] = 0;
                previous_velocity[id4_0] = velocity[id4_0];
                previous_velocity[id4_1] = velocity[id4_1];
                previous_velocity[id4_2] = velocity[id4_2];
                velocity[id4_0] = 0;
                velocity[id4_1] = 0;
                velocity[id4_2] = 0;
                force[id4_0] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
                force[id4_1] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
                force[id4_2] = 0.0; //It is done because Eulerian forces must be zero at the beginning of each iteration
                if (!is_solid[id3]){
                    for (INT32 i = 0; i < stencil_directions_F; i++)
                    {
                        INT32 id4_pop = idx_4d_vector(i, X, Y, Z);
                        double f = pop_F[id4_pop];
                        density[id3]     += f;
                        velocity[id4_0] += cx_pop_F[i] * f;
                        velocity[id4_1] += cy_pop_F[i] * f;
                        velocity[id4_2] += cz_pop_F[i] * f;
                    }
                    velocity[id4_0] /= density[id3];
                    velocity[id4_1] /= density[id3];
                    velocity[id4_2] /= density[id3];
                }else{
                    density[id3]    = 1.0;
                    velocity[id4_0] = 0.0;
                    velocity[id4_1] = 0.0;
                    velocity[id4_2] = 0.0;
                }
                numerator_residual_flow += pow(velocity[id4_0] - previous_velocity[id4_0], 2);
                numerator_residual_flow += pow(velocity[id4_1] - previous_velocity[id4_1], 2);
                numerator_residual_flow += pow(velocity[id4_2] - previous_velocity[id4_2], 2);
                denominator_residual_flow += pow(previous_velocity[id4_0], 2);
                denominator_residual_flow += pow(previous_velocity[id4_1], 2);
                denominator_residual_flow += pow(previous_velocity[id4_2], 2);
                #endif // STORE_1D_VECTOR
            }
        }
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: apply_sponge_layer_on_outlet_X
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::apply_sponge_layer_on_outlet_X(INT32 sponge_width){
/// This function applies a sponge layer at the outlet (x = Nx - 1) to prevent reflections. It relaxes the populations towards equilibrium based on local velocity and density.
    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    if (end_global_x[this_processor_id] == Nx - 1)
    {
        //cout << this_processor_id <<endl;
        int x_start_sponge = end_local_x[this_processor_id] - num_buffer_layers_X[this_processor_id] - sponge_width;
        for (INT32 X = x_start_sponge; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    for (INT32 i = 0; i < stencil_directions_F; i++)
                    {
                        double popEq_F{};
                        #if STORE_1D_VECTOR
                        int id3 = idx_3d_vector(X, Y, Z);
                        int id4_0 = idx_4d_vector(0, X, Y, Z);
                        int id4_1 = idx_4d_vector(1, X, Y, Z);
                        int id4_2 = idx_4d_vector(2, X, Y, Z);
                        double v0    = velocity[id4_0];
                        double v1    = velocity[id4_1];
                        double v2    = velocity[id4_2];
                        double rho   = density[id3];

                        INT32 id_pop = idx_4d_vector(i, X, Y, Z);
                        popEq_F = equilibrium_flow(rho, v0, v1, v2, i);

                        double alpha = (X - x_start_sponge) / (double)sponge_width;
                        pop_F[id_pop] = (1.0 - alpha) * pop_F[id_pop] + alpha * popEq_F;
                        #endif
                    }
                }
            }
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: scatter_force_modified_indices
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::scatter_force_modified_indices(const ParallelMPI& paralLB){
    /// Send indices of modified flow points (due to interpolation, spread) which are
    /// ... calculated in ParticleResolved::calculate_force_modified_indices from MASTER to other processors who owns that points

    // For global coordinates
    std::vector<INT32> send_coord_counts(num_processors_total, 0);
    std::vector<INT32> coord_displs(num_processors_total, 0);
    std::vector<INT32> send_coords_flat; // master only

    // Per-processor cell lists
    //proc_indices.clear();// To empty proc_indices (It has been already cleared in set_Eulerian_forces_in_master_zero(). Otherwise this line must be active.
    proc_indices.resize(num_processors_total);

    // MASTER PROCESS — packs and scatters
    if (this_processor_id == MASTER)
    {
        /// STEP 1: Create per-processor index lists
        //cout << force_modified_indices.size() <<endl;
        #if STORE_1D_VECTOR
        for (INT32 id3 : force_modified_indices)
        {
            INT32 X_local, Y_local, Z_local;
            decode_idx_3d(id3, X_local, Y_local, Z_local);

            INT32 X_global = paralLB.get_Xglobal_from_Xlocal(X_local, this_processor_id); //X_local - num_buffer_layers_X[this_processor_id];
            INT32 Y_global = paralLB.get_Yglobal_from_Ylocal(Y_local, this_processor_id); //Y_local - num_buffer_layers_Y[this_processor_id];
            INT32 Z_global = paralLB.get_Zglobal_from_Zlocal(Z_local, this_processor_id); //Z_local - num_buffer_layers_Z[this_processor_id];

            INT32 proc = paralLB.get_owner_processor_for_cell(X_global, Y_global, Z_global);
            if (proc >= 0)
                proc_indices[proc].push_back(id3);
        }
        #endif
        #if !STORE_1D_VECTOR
        for (auto [X_local, Y_local, Z_local] : force_modified_indices)
        {
            INT32 X_global = paralLB.get_Xglobal_from_Xlocal(X_local, this_processor_id); //X_local - num_buffer_layers_X[this_processor_id];
            INT32 Y_global = paralLB.get_Yglobal_from_Ylocal(Y_local, this_processor_id); //Y_local - num_buffer_layers_Y[this_processor_id];
            INT32 Z_global = paralLB.get_Zglobal_from_Zlocal(Z_local, this_processor_id); //Z_local - num_buffer_layers_Z[this_processor_id];

            INT32 proc = paralLB.get_owner_processor_for_cell(X_global, Y_global, Z_global);
            if (proc >= 0)
                proc_indices[proc].emplace_back(X_local, Y_local, Z_local);
        }
        #endif

        /// STEP 2: Flatten coordinates
        INT32 coord_offset = 0;
        for (INT32 p = 0; p < num_processors_total; ++p)
        {
           // cout << " proc_indices[p].size() for processor: " << p << " " <<proc_indices[p].size() <<endl;
            send_coord_counts[p] = static_cast<INT32>(proc_indices[p].size()) * 3; // X,Y,Z
            coord_displs[p] = coord_offset;
            coord_offset += send_coord_counts[p];
        }

        send_coords_flat.reserve(coord_offset);
        for (INT32 p = 0; p < num_processors_total; ++p)
        {
            #if STORE_1D_VECTOR
            for (INT32 id3 : proc_indices[p]) {
                INT32 X_local, Y_local, Z_local;
                decode_idx_3d(id3, X_local, Y_local, Z_local);

                INT32 X_global = paralLB.get_Xglobal_from_Xlocal(X_local, this_processor_id); //X_local - num_buffer_layers_X[this_processor_id];
                INT32 Y_global = paralLB.get_Yglobal_from_Ylocal(Y_local, this_processor_id); //Y_local - num_buffer_layers_Y[this_processor_id];
                INT32 Z_global = paralLB.get_Zglobal_from_Zlocal(Z_local, this_processor_id); //Z_local - num_buffer_layers_Z[this_processor_id];

                send_coords_flat.push_back(X_global);
                send_coords_flat.push_back(Y_global);
                send_coords_flat.push_back(Z_global);
            }
            #endif
            #if !STORE_1D_VECTOR
            for (auto [X_local, Y_local, Z_local] : proc_indices[p]) {
                INT32 X_global = paralLB.get_Xglobal_from_Xlocal(X_local, this_processor_id); //X_local - num_buffer_layers_X[this_processor_id];
                INT32 Y_global = paralLB.get_Yglobal_from_Ylocal(Y_local, this_processor_id); //Y_local - num_buffer_layers_Y[this_processor_id];
                INT32 Z_global = paralLB.get_Zglobal_from_Zlocal(Z_local, this_processor_id); //Z_local - num_buffer_layers_Z[this_processor_id];

                send_coords_flat.push_back(X_global);
                send_coords_flat.push_back(Y_global);
                send_coords_flat.push_back(Z_global);
            }
            #endif
        }
    }

    /// ---------------- STEP 3: BROADCAST COUNT INFO ----------------
    MPI_Bcast(send_coord_counts.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(coord_displs.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// ---------------- RECEIVE BUFFERS ----------------
    INT32 num_coords = send_coord_counts[this_processor_id];
    std::vector<INT32> recv_coords(num_coords);

    /// ---------------- STEP 4: CATTER COORDINATES ----------------
    MPI_Scatterv(this_processor_id == MASTER ? send_coords_flat.data() : nullptr,
                 send_coord_counts.data(), coord_displs.data(), MPI_INT,
                 recv_coords.data(), num_coords, MPI_INT,
                 MASTER, MPI_COMM_WORLD);

   /// ---------------- STEP 5: CONVERT GLOBAL TO LOCAL & UNPACK ----------------
    INT32 num_cells = num_coords / 3;
    force_modified_indices.clear();
    force_modified_indices.reserve(num_cells);//Allocates memory to hold at least num_cells elements without actually adding elements. The size() remains whatever it was before.

    INT32 counter = 0;
    for (INT32 i = 0; i < num_cells; ++i)
    {
        INT32 X_global = recv_coords[counter++];
        INT32 Y_global = recv_coords[counter++];
        INT32 Z_global = recv_coords[counter++];

        INT32 x_local = paralLB.get_Xlocal_from_Xglobal(X_global, this_processor_id); //X_global - start_global_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
        INT32 y_local = paralLB.get_Ylocal_from_Yglobal(Y_global, this_processor_id); //Y_global - start_global_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
        INT32 z_local = paralLB.get_Zlocal_from_Zglobal(Z_global, this_processor_id); //Z_global - start_global_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];

        #if STORE_1D_VECTOR
        force_modified_indices.push_back(idx_3d_vector(x_local, y_local, z_local));
        #endif
        #if !STORE_1D_VECTOR
        force_modified_indices.emplace_back(x_local, y_local, z_local);
        #endif
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: gather_flow_data_to_master_optimized
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::gather_fluid_macroscopic_data_to_master_optimized(){

    // 1) Prepare local send_buffer (every rank builds it, including MASTER)
    INT32 send_count = static_cast<INT32>(force_modified_indices.size()) * 5; // 3 doubles per cell (vx,vy,vz) + density + temperature
    std::vector<double> send_buffer;
    send_buffer.reserve(send_count);
    send_buffer.resize(send_count);

    if (this_processor_id != MASTER) {
        INT32 sidx = 0;
        #if STORE_1D_VECTOR
        for (size_t id3 : force_modified_indices) {
            size_t id4_0 = idx_4d_vector(0, id3);
            size_t id4_1 = idx_4d_vector(1, id3);
            size_t id4_2 = idx_4d_vector(2, id3);

            send_buffer[sidx++] = density[id3];
            send_buffer[sidx++] = velocity[id4_0];
            send_buffer[sidx++] = velocity[id4_1];
            send_buffer[sidx++] = velocity[id4_2];
            send_buffer[sidx++] = temperature[id3];
        }
        #else
        for (auto & tpl : force_modified_indices) {
            INT32 x = std::get<0>(tpl);
            INT32 y = std::get<1>(tpl);
            INT32 z = std::get<2>(tpl);

            send_buffer[sidx++] = density[x][y][z];
            send_buffer[sidx++] = velocity[0][x][y][z];
            send_buffer[sidx++] = velocity[1][x][y][z];
            send_buffer[sidx++] = velocity[2][x][y][z];
            send_buffer[sidx++] = temperature[x][y][z];
        }
        #endif
    }

    // 2) Gather counts at MASTER
    std::vector<INT32> recv_counts(num_processors_total, 0);
    MPI_Gather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MASTER, MPI_COMM_WORLD);

   // 3) MASTER computes displacements & total receive size
    std::vector<INT32> displs(num_processors_total, 0);
    std::vector<double> recv_buffer;
    if (this_processor_id == MASTER) {
        INT32 total_recv = 0;
        for (INT32 p = 0; p < num_processors_total; ++p) {
            displs[p] = total_recv;
            total_recv += recv_counts[p];
        }
        recv_buffer.resize(total_recv);
      //  cout << "total_recv: " << total_recv <<endl;
    }

    // 4) Gatherv all partial velocity buffers to MASTER
    MPI_Gatherv(send_buffer.data(), send_count, MPI_DOUBLE,
                this_processor_id == MASTER ? recv_buffer.data() : nullptr,
                recv_counts.data(), displs.data(), MPI_DOUBLE,
                MASTER, MPI_COMM_WORLD);

    // 5) MASTER unpacks velocities into global arrays using proc_indices ordering.
    if (this_processor_id == MASTER) {
        INT32 offset = 0;
        for (INT32 p = 1; p < num_processors_total; ++p) {
            const INT32 entries = static_cast<INT32>(proc_indices[p].size()); // number of cells that were assigned to proc p
          //  cout << "VVVVVVVVVVV " << p << " " << proc_indices[p].size() <<endl;
            for (INT32 j = 0; j < entries; ++j) {
                // recv order matches the order proc p had in its force_modified_indices
                double den= recv_buffer[offset++];
                double vx = recv_buffer[offset++];
                double vy = recv_buffer[offset++];
                double vz = recv_buffer[offset++];
                double temper = recv_buffer[offset++];

                #if STORE_1D_VECTOR
                // proc_indices[p] holds flat id3 (master local indices including buffers)
                size_t id3_master_local = proc_indices[p][j];
                INT32 X_local, Y_local, Z_local;
                decode_idx_3d(id3_master_local, X_local, Y_local, Z_local);

                // store into global (MASTER) arrays using global indices
                // When using STORE_1D_VECTOR your global storage for velocity is flattened with idx_4d_vector(q, X_global, Y_global, Z_global)
                size_t id4_vx = idx_4d_vector(0, X_local, Y_local, Z_local); // The same as idx_4d_vector(0, id3_master_local);
                size_t id4_vy = idx_4d_vector(1, X_local, Y_local, Z_local); // The same as idx_4d_vector(1, id3_master_local);
                size_t id4_vz = idx_4d_vector(2, X_local, Y_local, Z_local); // The same as idx_4d_vector(2, id3_master_local);

                density[id3_master_local] = den;
                velocity[id4_vx] = vx;
                velocity[id4_vy] = vy;
                velocity[id4_vz] = vz;
                temperature[id3_master_local] = temper;
                #endif
                #if !STORE_1D_VECTOR
                // proc_indices[p] holds tuples (X_local, Y_local, Z_local) as stored on MASTER earlier
                auto [X_local, Y_local, Z_local] = proc_indices[p][j];

                // MASTER's velocity arrays are indexed by global coordinates
                density[X_local][Y_local][Z_local]     = den;
                velocity[0][X_local][Y_local][Z_local] = vx;
                velocity[1][X_local][Y_local][Z_local] = vy;
                velocity[2][X_local][Y_local][Z_local] = vz;
                temperature[X_local][Y_local][Z_local] = temper;
                #endif
            } // j (cells for processor p)
        } // p (processors)
    } // if MASTER

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: gather_flow_data_to_master
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::gather_fluid_macroscopic_data_to_master(){

    /// This function gathers velocity data in MASTER to be used for interpolation in ParticleResolved

    /// Approach 1: push_back is removed in this function to make the process 10% faster compared to Approach 2 (below)
    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    // -----------------------------
    // 1) Preallocate and fill send_buffer (non-master only)
    // -----------------------------
    std::vector<double> send_buffer;
    INT32 send_count = 0;

    if (this_processor_id != MASTER) {
        const INT32 nx = end_x - start_x + 1;
        const INT32 ny = end_y - start_y + 1;
        const INT32 nz = end_z - start_z + 1;
        send_count = nx * ny * nz * 5; // density + 3 velocity comps + temperature

        send_buffer.resize(send_count);

        INT32 idx = 0;
        for (INT32 X = start_x; X <= end_x; ++X) {
            for (INT32 Y = start_y; Y <= end_y; ++Y) {
                for (INT32 Z = start_z; Z <= end_z; ++Z) {
                    #if !STORE_1D_VECTOR
                    send_buffer[idx++] = density[X][Y][Z];
                    send_buffer[idx++] = velocity[0][X][Y][Z];
                    send_buffer[idx++] = velocity[1][X][Y][Z];
                    send_buffer[idx++] = velocity[2][X][Y][Z];
                    send_buffer[idx++] = temperature[X][Y][Z];
                    #endif
                    #if STORE_1D_VECTOR
                    const size_t id3   = idx_3d_vector(X, Y, Z);
                    const size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    const size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    const size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    send_buffer[idx++] = density[id3];
                    send_buffer[idx++] = velocity[id4_0];
                    send_buffer[idx++] = velocity[id4_1];
                    send_buffer[idx++] = velocity[id4_2];
                    send_buffer[idx++] = temperature[id3];
                    #endif
                }
            }
        }
    }

    // -----------------------------
    // 2) Gather send_count from all processors
    // -----------------------------
    std::vector<INT32> recv_counts(num_processors_total, 0);
    MPI_Gather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    // -----------------------------
    // 3) On MASTER: compute displacements and total receive count
    // -----------------------------
    std::vector<INT32> displs(num_processors_total, 0);
    std::vector<double> recv_buffer;

    if (this_processor_id == MASTER) {
        INT32 total_recv_count = 0;
        for (INT32 i = 0; i < num_processors_total; ++i) {
            displs[i] = total_recv_count;
            total_recv_count += recv_counts[i];
        }
        recv_buffer.resize(total_recv_count);
    }

    // -----------------------------
    // 4) Gather all data to MASTER
    // -----------------------------
    MPI_Gatherv(send_buffer.data(), send_count, MPI_DOUBLE,
                recv_buffer.data(), recv_counts.data(), displs.data(),
                MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    // -----------------------------
    // 5) MASTER unpacks received data
    // -----------------------------
    if (this_processor_id == MASTER) {
        INT32 index = 0;
        for (INT32 proc = 1; proc < num_processors_total; ++proc) {  // start from 1: master didn’t send
            const INT32 x0 = start_global_x[proc] + num_buffer_layers_X[this_processor_id];
            const INT32 x1 = end_global_x[proc]   + num_buffer_layers_X[this_processor_id];
            const INT32 y0 = start_global_y[proc] + num_buffer_layers_Y[this_processor_id];
            const INT32 y1 = end_global_y[proc]   + num_buffer_layers_Y[this_processor_id];
            const INT32 z0 = start_global_z[proc] + num_buffer_layers_Z[this_processor_id];
            const INT32 z1 = end_global_z[proc]   + num_buffer_layers_Z[this_processor_id];

            for (INT32 X = x0; X <= x1; ++X) {
                for (INT32 Y = y0; Y <= y1; ++Y) {
                    for (INT32 Z = z0; Z <= z1; ++Z) {
                        #if !STORE_1D_VECTOR
                        density[X][Y][Z]     = recv_buffer[index++];
                        velocity[0][X][Y][Z] = recv_buffer[index++];
                        velocity[1][X][Y][Z] = recv_buffer[index++];
                        velocity[2][X][Y][Z] = recv_buffer[index++];
                        temperature[X][Y][Z] = recv_buffer[index++];
                        #endif
                        #if STORE_1D_VECTOR
                        const size_t id3   = idx_3d_vector(X, Y, Z);
                        const size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        const size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        const size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        density[id3]     = recv_buffer[index++];
                        velocity[id4_0]  = recv_buffer[index++];
                        velocity[id4_1]  = recv_buffer[index++];
                        velocity[id4_2]  = recv_buffer[index++];
                        temperature[id3] = recv_buffer[index++];
                        #endif
                    }
                }
            }
        }
    }

    /// Approach 2: push_back is used in this function to make send_buffer (functional but much slower than Approach 1)
 /*   const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X;
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X;
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y;
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y;
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z;
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z;

    /// Flatten the velocity data of other processors into a 1D buffer
    std::vector<double> send_buffer;

    if (this_processor_id != MASTER)//Important: MASTER should not send its data (that are zero) to itself
    {
        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    #if !STORE_1D_VECTOR
                    send_buffer.push_back(density[X][Y][Z]);
                    send_buffer.push_back(velocity[0][X][Y][Z]);
                    send_buffer.push_back(velocity[1][X][Y][Z]);
                    send_buffer.push_back(velocity[2][X][Y][Z]);
                    send_buffer.push_back(temperature[X][Y][Z]);
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    send_buffer.push_back(density[id3]);
                    send_buffer.push_back(velocity[id4_0]);
                    send_buffer.push_back(velocity[id4_1]);
                    send_buffer.push_back(velocity[id4_2]);
                    send_buffer.push_back(temperature[id3]);
                    #endif // STORE_1D_VECTOR

                }
            }
        }
    }

    /// Synchronize before proceeding
    //MPI_Barrier(MPI_COMM_WORLD);

    INT32 send_count = send_buffer.size();
    std::vector<INT32> recv_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);

    /// Gather the sizes of the buffers
    MPI_Gather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Calculate displacements on the master processor
    if (this_processor_id == MASTER) {
        displs[0] = 0;
        for (INT32 i = 1; i < num_processors_total; ++i) {
            displs[i] = displs[i-1] + recv_counts[i-1];
        }
    }

    /// Gather all the data to the master processor
    std::vector<double> recv_buffer;
    if (this_processor_id == MASTER) {
        INT32 total_recv_count = 0;
        for (INT32 count : recv_counts) {
            total_recv_count += count;
        }
        recv_buffer.resize(total_recv_count);
    }

    MPI_Gatherv(send_buffer.data(), send_count, MPI_DOUBLE,
            recv_buffer.data(), recv_counts.data(), displs.data(),
            MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (this_processor_id == 0) {
        INT32 index = 0;
        for (INT32 proc = 1; proc < num_processors_total; ++proc) {
            for (INT32 X = start_global_x[proc] + num_buffer_layers_X; X <= end_global_x[proc] + num_buffer_layers_X; ++X) {
                for (INT32 Y = start_global_y[proc] + num_buffer_layers_Y; Y <= end_global_y[proc] + num_buffer_layers_Y; ++Y) {
                    for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z; Z <= end_global_z[proc] + num_buffer_layers_Z; ++Z) {
                    //for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z; Z <= end_global_z[proc] + num_buffer_layers_Z; ++Z) {
                        #if !STORE_1D_VECTOR
                        density[X][Y][Z]     = recv_buffer[index++];
                        velocity[0][X][Y][Z] = recv_buffer[index++];
                        velocity[1][X][Y][Z] = recv_buffer[index++];
                        velocity[2][X][Y][Z] = recv_buffer[index++];
                        temperature[X][Y][Z] = recv_buffer[index++];
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        density[id3]     = recv_buffer[index++];
                        velocity[id4_0]  = recv_buffer[index++];
                        velocity[id4_1]  = recv_buffer[index++];
                        velocity[id4_2]  = recv_buffer[index++];
                        temperature[id3] = recv_buffer[index++];
                        #endif // STORE_1D_VECTOR
                    }
                }
            }
        }

        /// Print the velocity at a specific point
      //  INT32 x = 10, y = 20, z = 30;
       // std::cout << "Velocity at (" << x << ", " << y << ", " << z << "): " << "vx = " << velocity[0][x][y][z] << ", " << "vy = " << velocity[1][x][y][z] << ", " << "vz = " << velocity[2][x][y][z] << std::endl;
    }*/
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: scatter_fluid_force_from_master_optimized
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::scatter_fluid_force_from_master_optimized(const ParallelMPI& paralLB){
/// Data are only sent to processors whose force is modified in calculate_force_modified_indices and scatter_force_modified_indices function.

    // Each processor will receive 4 doubles per cell: Fx, Fy, Fz, Fthermal
    std::vector<INT32> send_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);

    INT32 total_send_count = 0;

    // MASTER PROCESS — packs and scatters
    if (this_processor_id == MASTER)
    {
        /// STEP 1: Prepare send_counts for forces (4 values per cell)
        total_send_count = 0;
        for (INT32 p = 0; p < num_processors_total; ++p)
        {
            INT32 num_cells = static_cast<INT32>(proc_indices[p].size());
            send_counts[p] = num_cells * 4; // Fx,Fy,Fz,Fthermal
            displs[p] = total_send_count;
            total_send_count += send_counts[p];
        }
    }

    /// ---------------- STEP 2: BROADCAST COUNT INFO ----------------
    MPI_Bcast(send_counts.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);

    INT32 recv_count = send_counts[this_processor_id];
    std::vector<double> recv_buffer(recv_count);

    std::vector<double> send_buffer;
    if (this_processor_id == MASTER) send_buffer.resize(total_send_count);


   /// ---------------- STEP 3: PACK FORCE DATA ON MASTER ----------------
    if (this_processor_id == MASTER)
    {
        INT32 offset = 0;
        for (INT32 p = 0; p < num_processors_total; ++p)
        {
            #if STORE_1D_VECTOR
            for (size_t id3 : proc_indices[p]) {
                size_t id4_0 = idx_4d_vector(0, id3);
                size_t id4_1 = idx_4d_vector(1, id3);
                size_t id4_2 = idx_4d_vector(2, id3);

                send_buffer[offset++] = force[id4_0];
                send_buffer[offset++] = force[id4_1];
                send_buffer[offset++] = force[id4_2];
                send_buffer[offset++] = force_thermal[id3];
            }
            #endif
            #if !STORE_1D_VECTOR
            for (auto [x, y, z] : proc_indices[p]) {
                send_buffer[offset++] = force[0][x][y][z];
                send_buffer[offset++] = force[1][x][y][z];
                send_buffer[offset++] = force[2][x][y][z];
                send_buffer[offset++] = force_thermal[x][y][z];
            }
            #endif
        }
    }

    /// ---------------- STEP 4: SCATTER FORCES ----------------//
    MPI_Scatterv(this_processor_id == MASTER ? send_buffer.data() : nullptr,
                 send_counts.data(), displs.data(), MPI_DOUBLE,
                 recv_buffer.data(), recv_count, MPI_DOUBLE,
                 MASTER, MPI_COMM_WORLD);

    /// ---------------- STEP 5: UPDATE LOCAL FORCES ----------------//
    INT32 counter = 0;
    #if STORE_1D_VECTOR
    for (size_t id3 : force_modified_indices) {
        size_t id4_0 = idx_4d_vector(0, id3);
        size_t id4_1 = idx_4d_vector(1, id3);
        size_t id4_2 = idx_4d_vector(2, id3);

        force[id4_0]       = recv_buffer[counter++];
        force[id4_1]       = recv_buffer[counter++];
        force[id4_2]       = recv_buffer[counter++];
        force_thermal[id3] = recv_buffer[counter++];
    }
    #endif
    #if !STORE_1D_VECTOR
    for (auto [x, y, z] : force_modified_indices) {
        force[0][x][y][z]       = recv_buffer[counter++];
        force[1][x][y][z]       = recv_buffer[counter++];
        force[2][x][y][z]       = recv_buffer[counter++];
        force_thermal[x][y][z]  = recv_buffer[counter++];
    }
    #endif
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: scatter_fluid_force_from_master
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::scatter_fluid_force_from_master(){

    /// Approach 1: Does not use push_back and is therefore faster than Approach2 (below and commented)
    // ---------- 1) calculate send_counts and displs (MASTER) ----------
    std::vector<INT32> send_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);
    INT32 total_size = 0;

    if (this_processor_id == MASTER)
    {
        for (INT32 proc = 0; proc < num_processors_total; ++proc) {
            INT32 x_size = end_global_x[proc] - start_global_x[proc] + 1;
            INT32 y_size = end_global_y[proc] - start_global_y[proc] + 1;
            INT32 z_size = end_global_z[proc] - start_global_z[proc] + 1;
            send_counts[proc] = x_size * y_size * z_size * 4; // 4 values per cell
            displs[proc] = total_size;
            total_size += send_counts[proc];
        }
    }

    // Broadcast send_counts and displs
    MPI_Bcast(send_counts.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);


    // ---------- 2) MASTER: allocate send_buffer once ----------
    std::vector<double> send_buffer;
    if (this_processor_id == MASTER) {
        send_buffer.resize(total_size); // preallocate exact size
    }

    // ---------- 3) MASTER: fill send_buffer by writing into known offsets ----------
    if (this_processor_id == MASTER) {
        // We'll fill each proc's chunk at send_buffer[displs[proc] + local_index]
        // This avoids push_back and allocations.
        for (INT32 proc = 0; proc < num_processors_total; ++proc) {
            const INT32 base_disp = displs[proc];
            INT32 local_index = 0;

            const INT32 x0 = start_global_x[proc] + num_buffer_layers_X[this_processor_id];
            const INT32 x1 = end_global_x[proc]   + num_buffer_layers_X[this_processor_id];
            const INT32 y0 = start_global_y[proc] + num_buffer_layers_Y[this_processor_id];
            const INT32 y1 = end_global_y[proc]   + num_buffer_layers_Y[this_processor_id];
            const INT32 z0 = start_global_z[proc] + num_buffer_layers_Z[this_processor_id];
            const INT32 z1 = end_global_z[proc]   + num_buffer_layers_Z[this_processor_id];

            // Iterate and write directly to send_buffer[base_disp + local_index]
            for (INT32 X = x0; X <= x1; ++X) {
                for (INT32 Y = y0; Y <= y1; ++Y) {
                    // If possible, move some index arithmetic out: e.g. compute idx for X,Y once per Z loop
                    for (INT32 Z = z0; Z <= z1; ++Z) {
                        #if !STORE_1D_VECTOR
                        send_buffer[base_disp + local_index++] = force[0][X][Y][Z];// Much faster than push_back
                        send_buffer[base_disp + local_index++] = force[1][X][Y][Z];
                        send_buffer[base_disp + local_index++] = force[2][X][Y][Z];
                        send_buffer[base_disp + local_index++] = force_thermal[X][Y][Z];
                        #endif
                        #if STORE_1D_VECTOR
                        // minimize calls to idx functions by computing base ids once if possible
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);// id4_0 + 1 ;
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);//id4_0 + 2;
                        send_buffer[base_disp + local_index++] = force[id4_0];
                        send_buffer[base_disp + local_index++] = force[id4_1];
                        send_buffer[base_disp + local_index++] = force[id4_2];
                        send_buffer[base_disp + local_index++] = force_thermal[id3];
                        #endif
                    }
                }
            }
            // Sanity check (optional)
            // if (local_index != send_counts[proc]) { std::cerr << "pack size mismatch proc " << proc << "\n"; }
        }
       // std::cout << "Master: send_buffer size = " << send_buffer.size() << std::endl;
    }

    // ---------- 4) Scatterv ----------
    INT32 recv_count = send_counts[this_processor_id];
    std::vector<double> recv_buffer(recv_count);

    MPI_Scatterv(this_processor_id == MASTER ? send_buffer.data() : nullptr,
                 send_counts.data(), displs.data(), MPI_DOUBLE,
                 recv_buffer.data(), recv_count, MPI_DOUBLE,
                 MASTER, MPI_COMM_WORLD);

    // ---------- 5) unpack on receivers (same as you had, but use recv_count if needed) ----------
    if (this_processor_id != MASTER)
    {
        const INT32 x_size = end_global_x[this_processor_id] - start_global_x[this_processor_id] + 1;
        const INT32 y_size = end_global_y[this_processor_id] - start_global_y[this_processor_id] + 1;
        const INT32 z_size = end_global_z[this_processor_id] - start_global_z[this_processor_id] + 1;

        INT32 index = 0;
        for (INT32 X = num_buffer_layers_X[this_processor_id]; X < x_size + num_buffer_layers_X[this_processor_id]; ++X) {
            for (INT32 Y = num_buffer_layers_Y[this_processor_id]; Y < y_size + num_buffer_layers_Y[this_processor_id]; ++Y) {
                for (INT32 Z = num_buffer_layers_Z[this_processor_id]; Z < z_size + num_buffer_layers_Z[this_processor_id]; ++Z) {
                    #if !STORE_1D_VECTOR
                    force[0][X][Y][Z] = recv_buffer[index++];
                    force[1][X][Y][Z] = recv_buffer[index++];
                    force[2][X][Y][Z] = recv_buffer[index++];
                    force_thermal[X][Y][Z] = recv_buffer[index++];
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    force[id4_0] = recv_buffer[index++];
                    force[id4_1] = recv_buffer[index++];
                    force[id4_2] = recv_buffer[index++];
                    force_thermal[id3] = recv_buffer[index++];
                    #endif
                }
            }
        }
    }

    /// Approach 2: Similar to Approach 1 but slower because of using push_back
   /* std::vector<double> send_buffer;

    /// Flatten the force data into a 1D buffer (send_buffer)
     if (this_processor_id == MASTER) {
        INT32 index = 0;
        for (INT32 proc = 0; proc < num_processors_total; ++proc) {
            for (INT32 X = start_global_x[proc] + num_buffer_layers_X; X <= end_global_x[proc] + num_buffer_layers_X; ++X) {
                for (INT32 Y = start_global_y[proc] + num_buffer_layers_Y; Y <= end_global_y[proc] + num_buffer_layers_Y; ++Y) {
                    for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z; Z <= end_global_z[proc] + num_buffer_layers_Z; ++Z) {
                   // for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z; Z <= end_global_z[proc] + num_buffer_layers_Z; ++Z) {
                        #if !STORE_1D_VECTOR
                        send_buffer.push_back(force[0][X][Y][Z]);
                        send_buffer.push_back(force[1][X][Y][Z]);
                        send_buffer.push_back(force[2][X][Y][Z]);
                        send_buffer.push_back(force_thermal[X][Y][Z]);
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        send_buffer.push_back(force[id4_0]);
                        send_buffer.push_back(force[id4_1]);
                        send_buffer.push_back(force[id4_2]);
                        send_buffer.push_back(force_thermal[id3]);
                        #endif // STORE_1D_VECTOR
                    }
                }
            }
        }
        //std::cout << "Master: send_buffer size = " << send_buffer.size() << std::endl;
    }

    /// Calculate send counts and displacements for each processor
    std::vector<INT32> send_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);
    if (this_processor_id == MASTER)//only MASTER calculates send_counts and displs to save time
    {
        INT32 index = 0;
        for (INT32 proc = 0; proc < num_processors_total; ++proc)
        {
            INT32 x_size = end_global_x[proc] - start_global_x[proc] + 1;
            INT32 y_size = end_global_y[proc] - start_global_y[proc] + 1;
            INT32 z_size = end_global_z[proc] - start_global_z[proc] + 1;
            //std::cout << "Proc " << proc << ": " << x_size << " " << y_size << " " << z_size << std::endl;
            send_counts[proc] = x_size * y_size * z_size * 4;// 3 components for force + 1 for force_thermal
            displs[proc] = index;
            index += send_counts[proc];
        }
    }

    /// Broadcast send_counts and displs from the master to all processors
    MPI_Bcast(send_counts.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Scatter the data from the master processor to all processors
    INT32 recv_count = send_counts[this_processor_id];
    std::vector<double> recv_buffer(recv_count);
    //std::cout << "Processor " << this_processor_id << ": recv_count = " << recv_count << std::endl;

    MPI_Scatterv(this_processor_id == MASTER ? send_buffer.data() : nullptr, send_counts.data(), displs.data(), MPI_DOUBLE,
              recv_buffer.data(), send_counts[this_processor_id], MPI_DOUBLE,
              MASTER, MPI_COMM_WORLD);

    /// Reshape the received data back into 3D arrays on each processor
    if (this_processor_id != MASTER) {
        INT32 x_size = end_global_x[this_processor_id] - start_global_x[this_processor_id] + 1;
        INT32 y_size = end_global_y[this_processor_id] - start_global_y[this_processor_id] + 1;
        INT32 z_size = end_global_z[this_processor_id] - start_global_z[this_processor_id] + 1;

        INT32 index = 0;
        for (INT32 X = num_buffer_layers_X; X < x_size + num_buffer_layers_X; ++X) {
            for (INT32 Y = num_buffer_layers_Y; Y < y_size + num_buffer_layers_Y; ++Y) {
                for (INT32 Z = num_buffer_layers_Z; Z < z_size + num_buffer_layers_Z; ++Z) {
                    #if !STORE_1D_VECTOR
                    force[0][X][Y][Z] = recv_buffer[index++];
                    force[1][X][Y][Z] = recv_buffer[index++];
                    force[2][X][Y][Z] = recv_buffer[index++];
                    force_thermal[X][Y][Z] = recv_buffer[index++];
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    force[id4_0] = recv_buffer[index++];
                    force[id4_1] = recv_buffer[index++];
                    force[id4_2] = recv_buffer[index++];
                    force_thermal[id3] = recv_buffer[index++];
                    #endif // STORE_1D_VECTOR
                }
            }
        }
    }*/

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: set_forces_zero
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::set_Eulerian_forces_in_master_zero()
{
    force_modified_indices.clear();
    proc_indices.clear();// To empty proc_indices

    #if STORE_1D_VECTOR
    std::fill(force.begin(), force.end(), 0.0);//fastest solution for 1D vector
    std::fill(force_thermal.begin(), force_thermal.end(), 0.0);
    std::fill(cell_marked.begin(), cell_marked.end(), false);
    #endif

    #if !STORE_1D_VECTOR
    for (INT32 X = start_local_x[this_processor_id]; X <= end_local_x[this_processor_id]; ++X)
    {
        for (INT32 Y = start_local_y[this_processor_id]; Y <= end_local_y[this_processor_id]; ++Y)
        {
            for (INT32 Z = start_local_z[this_processor_id]; Z <= end_local_z[this_processor_id]; ++Z)
            {
                force[0][X][Y][Z] = 0;
                force[1][X][Y][Z] = 0;
                force[2][X][Y][Z] = 0;
                force_thermal[X][Y][Z] = 0;
            }
        }
    }
    #endif
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Guo_part_two
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Guo_part_two(const ParallelMPI& paralLB, const ParticleResolved* Particle_Model){

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                INT32 gbetaindex = 1;
                if (Particle_Model && Flow_With_Particle && Flow_With_Thermal_Effect){
                    gbetaindex = Particle_Model->calculate_gbetaindex(X, Y, Z, paralLB);// To set gbetaindex=0 for fluid nodes inside the particle
                }
                #if !STORE_1D_VECTOR
                double forceX = force[0][X][Y][Z] + gbetaindex * gbeta_x * (temperature[X][Y][Z] - Tref) + gravity;
                double forceY = force[1][X][Y][Z] + gbetaindex * gbeta_y * (temperature[X][Y][Z] - Tref);
                double forceZ = force[2][X][Y][Z] + gbetaindex * gbeta_z * (temperature[X][Y][Z] - Tref);
                velocity[0][X][Y][Z] += forceX / 2.0;
                velocity[1][X][Y][Z] += forceY / 2.0;
                velocity[2][X][Y][Z] += forceZ / 2.0;
                #endif
                #if STORE_1D_VECTOR
                size_t id3   = idx_3d_vector(X, Y, Z);
                size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                double forceX = force[id4_0] + gbetaindex * gbeta_x * (temperature[id3] - Tref) + gravity;
                double forceY = force[id4_1] + gbetaindex * gbeta_y * (temperature[id3] - Tref);
                double forceZ = force[id4_2] + gbetaindex * gbeta_z * (temperature[id3] - Tref);
                velocity[id4_0] += forceX / 2.0;
                velocity[id4_1] += forceY / 2.0;
                velocity[id4_2] += forceZ / 2.0;
                #endif // STORE_1D_VECTOR
            }
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: equilibrium_thermal
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the equilibrium populations from the fluid temperature.
// It computes the equilibrium only at a specific lattice node: Function has to be called at each lattice node.
// The standard quadratic euilibrium is used.
double FlowSolver::equilibrium_thermal(double tem, double den, double vel_x, double vel_y, double vel_z, INT32 i)
{
    double pop_eq_T{};
    /// Equilibrium_Thermal_1
    if (Equilibrium_Thermal_1)
    {
        #if !STORE_1D_VECTOR
        pop_eq_T = weight_T[i] * tem * (1 + 1.0 / cs_2_Thermal * (cx_pop_T[i] * vel_x + cy_pop_T[i] * vel_y + cz_pop_T[i] * vel_z));
        #endif // !STORE_1D_VECTOR
        #if STORE_1D_VECTOR
        pop_eq_T = weight_T[i] * tem * (1 + 1.0 / cs_2_Thermal * (cx_pop_T[i] * vel_x + cy_pop_T[i] * vel_y + cz_pop_T[i] * vel_z));
        #endif
    }
    /// Equilibrium_Thermal_2: Recommended according to page 309 of Krueger's book
    if (Equilibrium_Thermal_2)
    {
        #if !STORE_1D_VECTOR
         pop_eq_T = weight_T[i] * tem * (1 + 1.0 / cs_2_Thermal * (cx_pop_T[i] * vel_x + cy_pop_T[i] * vel_y + cz_pop_T[i] * vel_z) + 1.0 / (2.0 * cs_2_Thermal * cs_2_Thermal) * SQ(cx_pop_T[i] * vel_x + cy_pop_T[i] * vel_y + cz_pop_T[i] * vel_z)
                              - 1.0 / (2.0 * cs_2_Thermal) * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
        #endif
        #if STORE_1D_VECTOR
         pop_eq_T= weight_T[i] * tem * (1 + 1.0 / cs_2_Thermal * (cx_pop_T[i] * vel_x + cy_pop_T[i] * vel_y + cz_pop_T[i] * vel_z) 
                                          + 1.0 / (2.0 * cs_2_Thermal * cs_2_Thermal) * SQ(cx_pop_T[i] * vel_x + cy_pop_T[i] * vel_y + cz_pop_T[i] * vel_z)
                                         - 1.0 / (2.0 * cs_2_Thermal) * (SQ(vel_x) + SQ(vel_y) + SQ(vel_z)));
        #endif
    }
    return pop_eq_T;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: LBM_thermal
//-----------------------------------------------------------------------------------------------------------------------------
//
//#if defined Flow_With_Thermal_Effect
void FlowSolver::LBM_thermal(INT32 time)
{

    /// Swap populations
    // The present code used old and new populations which are swapped at the beginning of each time step.
    // This is sometimes called 'double-buffered' or 'ping-pong' algorithm.
    // This way, the old populations are not overwritten during propagation.
    // The resulting code is easier to write and to debug.
    // The memory requirement for the populations is twice as large.

    /*vector<vector<vector<vector<double>>>> swap_temp_t;
    swap_temp_t = pop_old_T;
    pop_old_T = pop_T;
    pop_T = swap_temp_t;*/
    std::swap(pop_T, pop_old_T); //much (almost 2.5 times) faster than the 3 above lines.


    /// Lattice Boltzmann equation
    // The lattice Boltzmann equation is solved in the following.
    // The algorithm includes
    // - computation of the lattice force
    // - combined collision and propagation (faster than first collision and then propagation)

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    //============================================================//
    double source_thermal{0};
    /*double T_wall = 1.0;
    double tgrad{0};
    size_t id4_velx;
    double local_Nusselt{};
    double added_q = 0;
    for (int X = start_x; X <= end_x; X++){
        size_t id1 = idx_3d_vector(X, 1, 1);
        size_t id2 = idx_3d_vector(X, 2, 1);

        tgrad = (-3.0 * T_wall + 4.0 * temperature[id1] - temperature[id2]) / 2.0; //second-order accurate
        //tgrad   = (T_wall - temperature[id1]) / 0.5;// 0.5 is due two half way bounce back (First order accurate)). Interesting: This does not work properly here but for nusselt number.
        //tgrad = (temperature[idx_3d_vector(X, 1, 1)] - temperature[idx_3d_vector(X, 2, 1)]);
        added_q += -alpha_t * tgrad;
    }
    source_thermal =  2 * added_q / ((Nx - 1) * (Ny - 1));// 2 because of upper and lower wall, (Nx-1)*(Ny-1) is the number of fluid nodes in x-y plane
    */
    //============================================================//

    INT32 newX{}, newY{}, newZ{};
    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                bool do_LBM_thermal;
                double omega_thermal;
                double popEq_T{};
                get_do_LBM_thermal_omega_thermal(do_LBM_thermal, omega_thermal, X, Y, Z);

                /// Compute new populations
                // This is the lattice Boltzmann equation (combined collision and propagation).
                // Periodicity of the lattice in x-direction is taken into account by the %-operator.
///===================================================================================================================================== ///
///                                                 Force_Thermal_1                                                                      ///
///===================================================================================================================================== ///
                if (Force_Thermal_1 && do_LBM_thermal)
                {
                    for (INT32 i = 0; i < stencil_directions_T; i++)
                    {
                        #if !STORE_1D_VECTOR
                        newX = X + cx_pop_T[i];//(X + cx_pop_T[i] + Nx) % Nx
                        newY = Y + cy_pop_T[i];//(Y + cy_pop_T[i] + Ny) % Ny;
                        newZ = Z + cz_pop_T[i];//(Z + cz_pop_T[i] + Nz) % Nz;
                        popEq_T        = equilibrium_thermal(temperature[X][Y][Z], density[X][Y][Z], velocity[0][X][Y][Z], velocity[1][X][Y][Z], velocity[2][X][Y][Z], i);
                        pop_T[i][newX][newY][newZ] =  pop_old_T[i][X][Y][Z] * (1.0 - omega_thermal) + popEq_T * omega_thermal + weight_T[i] * force_thermal[X][Y][Z];
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        size_t id4_pop_i     = idx_4d_vector(i, X, Y, Z);
                        popEq_T   = equilibrium_thermal(temperature[id3], density[id3], velocity[id4_0], velocity[id4_1], velocity[id4_2], i);
                        newX = X + cx_pop_T[i];//(X + cx_pop_T[i] + Nx) % Nx
                        newY = Y + cy_pop_T[i];//(Y + cy_pop_T[i] + Ny) % Ny;
                        newZ = Z + cz_pop_T[i];//(Z + cz_pop_T[i] + Nz) % Nz;
                        size_t id4_pop_new_i = idx_4d_vector(i, newX, newY, newZ);
                        pop_T[id4_pop_new_i] = pop_old_T[id4_pop_i] * (1.0 - omega_thermal) + popEq_T * omega_thermal + weight_T[i] * force_thermal[id3] + weight_T[i] * source_thermal;
                        #endif
                    }
                }
///===================================================================================================================================== ///
///                                                 Force_Thermal_2                                                                      ///
///===================================================================================================================================== ///
                if (Force_Thermal_2 && do_LBM_thermal)
                {
                    for (INT32 i = 0; i < stencil_directions_T; i++)
                    {
                        #if !STORE_1D_VECTOR
                        newX = X + cx_pop_T[i];//(X + cx_pop_T[i] + Nx) % Nx
                        newY = Y + cy_pop_T[i];//(Y + cy_pop_T[i] + Ny) % Ny;
                        newZ = Z + cz_pop_T[i];//(Z + cz_pop_T[i] + Nz) % Nz;
                        popEq_T = equilibrium_thermal(temperature[X][Y][Z], density[X][Y][Z], velocity[0][X][Y][Z], velocity[1][X][Y][Z], velocity[2][X][Y][Z], i);
                        pop_T[i][newX][newY][newZ] = pop_old_T[i][X][Y][Z] * (1.0 - omega_thermal) + popEq_T * omega_thermal + (1 - 0.5 * omega_thermal) * weight_T[i] * force_thermal[X][Y][Z];
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        size_t id4_pop_i = idx_4d_vector(i, X, Y, Z);
                        popEq_T = equilibrium_thermal(temperature[id3], density[id3], velocity[id4_0], velocity[id4_1], velocity[id4_2], i);
                        newX = X + cx_pop_T[i];//(X + cx_pop_T[i] + Nx) % Nx
                        newY = Y + cy_pop_T[i];//(Y + cy_pop_T[i] + Ny) % Ny;
                        newZ = Z + cz_pop_T[i];//(Z + cz_pop_T[i] + Nz) % Nz;
                        size_t id4_pop_new_i = idx_4d_vector(i, newX, newY, newZ);
                        pop_T[id4_pop_new_i] = pop_old_T[id4_pop_i] * (1.0 - omega_thermal) + popEq_T * omega_thermal + (1 - 0.5 * omega_thermal) * weight_T[i] * force_thermal[id3];
                        #endif
                    }//End of for (INT32 i = 0; i < stencil_directions_T; i++)
                }// End of Force_Thermal_2
            }
        }
    }
}
//#endif
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_do_LBM_thermal_omega_thermal
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::get_do_LBM_thermal_omega_thermal(bool& do_LBM_th, double& omega_th, INT32 X, INT32 Y, INT32 Z){
    do_LBM_th = true;
    omega_th  = omega_t;//default value (uses fluid conductivity)

    if (!Flow_With_Particle){
        #if !STORE_1D_VECTOR
        if (solid_thermal[X][Y][Z] == 1 || solid_thermal[X][Y][Z] == 2){//1 means constant temperature; 2 means adiabatic
            do_LBM_th = false;//no thermal inside solid nodes is needed
        }
        if (solid_thermal[X][Y][Z] == 0){//0 means fluid
            do_LBM_th = true;
            omega_th = omega_t;
        }
        if (solid_thermal[X][Y][Z] == 3){//3 means conjugate heat transfer
            do_LBM_th = true;
            omega_th = omega_ts;
        }
        #endif
        #if STORE_1D_VECTOR
        size_t id3   = idx_3d_vector(X, Y, Z);
        if (solid_thermal[id3] == 1 || solid_thermal[id3] == 2){//1 means constant temperature; 2 means adiabatic
            do_LBM_th = false;//no thermal inside solid nodes is needed
        }
        if (solid_thermal[id3] == 0){//0 means fluid
            do_LBM_th = true;
            omega_th = omega_t;
        }
        if (solid_thermal[id3] == 3){//3 means conjugate heat transfer
            do_LBM_th = true;
            omega_th = omega_ts;
        }
        #endif
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: bc_at_walls_thermal
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::bc_at_walls_thermal(const ParallelMPI& paralLB)
{

    INT32 surf_normal_X{}, surf_normal_Y{}, surf_normal_Z{};
    double temperature_at_bc{};
    double heat_flux_at_bc{};

    for (INT32 ii = 0; ii < 6; ii++)
    {
        INT32 start_BC_X{}, end_BC_X{}, start_BC_Y{}, end_BC_Y{}, start_BC_Z{}, end_BC_Z{};
        INT32 tag;
        std::vector<INT32> pop_to, pop_reverse;
        std::string BC_type;

        pop_to.clear();
        pop_reverse.clear();
        BC_type = "";

        switch(ii)
        {
        case 0:
            /// Right BC
            if (end_global_x[this_processor_id] == this->Nx - 1)
            {
                pop_to      = pop_to_right_T;
                pop_reverse = pop_reverse_right_T;
                start_BC_X  = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                //cout << "Right, processor:" << processor_id << " "<<start_BC_X<< "-"<<end_BC_X << ", "<< start_BC_Y<< "-"<<end_BC_Y << ", "<< start_BC_Z<< "-"<<end_BC_Z << endl;
                if (Right_Fixed_Temp)
                {
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_right;
                }
                if (Right_Adiabatic)
                {
                    BC_type = "adiabatic";
                    surf_normal_X = -1;
                    surf_normal_Y = 0;
                    surf_normal_Z = 0;
                }
                if (Right_Fixed_Heat_Flux)
                {
                    BC_type = "fixed_heat_flux";
                    heat_flux_at_bc = heat_flux_right;
                    surf_normal_X = -1;
                    surf_normal_Y = 0;
                    surf_normal_Z = 0;
                }
            }
        break;
        case 1:
            /// Left BC
            if (start_global_x[this_processor_id] == 0)
            {
                pop_to      = pop_to_left_T;
                pop_reverse = pop_reverse_left_T;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
                end_BC_X    = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                if (Left_Fixed_Temp)
                {
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_left;
                }
                if (Left_Adiabatic)
                {
                    BC_type = "adiabatic";
                    surf_normal_X = 1;
                    surf_normal_Y = 0;
                    surf_normal_Z = 0;
                }
                if (Left_Fixed_Heat_Flux)
                {
                    BC_type = "fixed_heat_flux";
                    heat_flux_at_bc = heat_flux_left;
                    surf_normal_X = 1;
                    surf_normal_Y = 0;
                    surf_normal_Z = 0;
                }
            }
        break;
        case 2:
            /// Top BC
            if (end_global_y[this_processor_id] == this->Ny - 1)
            {
                pop_to      = pop_to_top_T;
                pop_reverse = pop_reverse_top_T;
                start_BC_Y  = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                if (Top_Fixed_Temp)
                {
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_top;
                }
                if (Top_Adiabatic)
                {
                    BC_type = "adiabatic";
                    surf_normal_X = 0;
                    surf_normal_Y = -1;
                    surf_normal_Z = 0;
                }
                if (Top_Fixed_Heat_Flux)
                {
                    BC_type = "fixed_heat_flux";
                    heat_flux_at_bc = heat_flux_top;
                    surf_normal_X = 0;
                    surf_normal_Y = -1;
                    surf_normal_Z = 0;
                }
            }
            break;
        case 3:
            /// Bottom BC
            if (start_global_y[this_processor_id] == 0)
            {
                pop_to      = pop_to_bottom_T;
                pop_reverse = pop_reverse_bottom_T;
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
                end_BC_Y    = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
                tag = ii;
                if (Bottom_Fixed_Temp)
                {
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_bot;
                }
                if (Bottom_Adiabatic)
                {
                    BC_type = "adiabatic";
                    surf_normal_X = 0;
                    surf_normal_Y = 1;
                    surf_normal_Z = 0;
                }
                if (Bottom_Fixed_Heat_Flux)
                {
                    BC_type = "fixed_heat_flux";
                    heat_flux_at_bc = heat_flux_bottom;
                    surf_normal_X = 0;
                    surf_normal_Y = 1;
                    surf_normal_Z = 0;
                }
            }
            break;
        case 4:
            /// Front BC
            if (end_global_z[this_processor_id] == this->Nz - 1)
            {
                pop_to      = pop_to_front_T;
                pop_reverse = pop_reverse_front_T;
                start_BC_Z  = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;
                end_BC_Z    = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                tag = ii;
                if (Front_Fixed_Temp)
                {
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_front;
                }
                if (Front_Adiabatic)
                {
                    BC_type = "adiabatic";
                    surf_normal_X = 0;
                    surf_normal_Y = 0;
                    surf_normal_Z = -1;
                }
                if (Front_Fixed_Heat_Flux)
                {
                    BC_type = "fixed_heat_flux";
                    heat_flux_at_bc = heat_flux_front;
                    surf_normal_X = 0;
                    surf_normal_Y = 0;
                    surf_normal_Z = -1;
                }
            }
            tag = ii;
            break;
        case 5:
            /// Back BC
            if (start_global_z[this_processor_id] == 0)
            {
                pop_to      = pop_to_back_T;
                pop_reverse = pop_reverse_back_T;
                start_BC_Z  = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
                end_BC_Z    = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
                start_BC_X  = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
                end_BC_X    = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
                start_BC_Y  = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
                end_BC_Y    = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
                tag = ii;
                if (Back_Fixed_Temp)
                {
                    BC_type = "fixed_temperature";
                    temperature_at_bc = temperature_back;
                }
                if (Back_Adiabatic)
                {
                    BC_type = "adiabatic";
                    surf_normal_X = 0;
                    surf_normal_Y = 0;
                    surf_normal_Z = 1;
                }
                if (Back_Fixed_Heat_Flux)
                {
                    BC_type = "fixed_heat_flux";
                    heat_flux_at_bc = heat_flux_back;
                    surf_normal_X = 0;
                    surf_normal_Y = 0;
                    surf_normal_Z = 1;
                }
            }
            tag = ii;
            break;
        default:
            std::cout << "Invalid choice" << std::endl;
            break;
        }
        //===============================================================//
        for (INT32 X = start_BC_X; X <= end_BC_X; ++X)
        {
            for (INT32 Y = start_BC_Y; Y <= end_BC_Y; ++Y)
            {
                for (INT32 Z = start_BC_Z; Z <= end_BC_Z; ++Z)
                {
                    //===============================================================//
                    if (BC_type == "fixed_temperature"){
                        for (size_t alpha = 0; alpha < pop_to.size(); alpha++)
                        {
                            INT32 pop     = pop_to[alpha];
                            INT32 pop_rev = pop_reverse[alpha];
                            #if !STORE_1D_VECTOR
                            //cout << "BB " << pop_to[alpha] << " " << pop_reverse[alpha] <<endl;
                            pop_T[pop_rev][X][Y][Z] = (weight_T[pop_rev] + weight_T[pop]) * temperature_at_bc - pop_T[pop][X + cx_pop_T[pop]][Y + cy_pop_T[pop]][Z + cz_pop_T[pop]];
                            #endif
                            #if STORE_1D_VECTOR
                            size_t id4_pop = idx_4d_vector(pop, X + cx_pop_T[pop], Y + cy_pop_T[pop], Z + cz_pop_T[pop]);
                            size_t id4_rev = idx_4d_vector(pop_rev, X, Y, Z);
                            pop_T[id4_rev] = (weight_T[pop_rev] + weight_T[pop]) * temperature_at_bc - pop_T[id4_pop];
                            #endif
                        }
                    }//End of fixed temperature
                    //===============================================================//
                    if (BC_type == "adiabatic"){
                        for (size_t alpha = 0; alpha < pop_to.size(); alpha++)
                        {
                            INT32 pop     = pop_to[alpha];
                            INT32 pop_rev = pop_reverse[alpha];
                            #if !STORE_1D_VECTOR
                            pop_T[pop_rev][X][Y][Z] = pop_T[pop_rev][X + surf_normal_X][Y + surf_normal_Y][Z + surf_normal_Z];//first order
                            //pop_T[pop_rev][X][Y][Z] = 2 * pop_T[pop_rev][X + 2 * surf_normal_X][Y + 2 * surf_normal_Y][Z + 2* surf_normal_Z] - pop_T[pop_rev][X + surf_normal_X][Y + surf_normal_Y][Z + surf_normal_Z];// second order
                            #endif
                            #if STORE_1D_VECTOR
                            size_t id4_pop = idx_4d_vector(pop_rev, X + surf_normal_X, Y + surf_normal_Y, Z + surf_normal_Z);
                            size_t id4_rev = idx_4d_vector(pop_rev, X, Y, Z);
                            pop_T[id4_rev] = pop_T[id4_pop];//first order
                            #endif
                        }
                    }// End of adiabatic
                    //===============================================================//
                    if (BC_type == "fixed_heat_flux"){
                        for (size_t alpha = 0; alpha < pop_to.size(); alpha++)
                        {
                            INT32 pop     = pop_to[alpha];
                            INT32 pop_rev = pop_reverse[alpha];
                            #if !STORE_1D_VECTOR
                            //double T_wall_at_bc = temperature[X + surf_normal_X][Y + surf_normal_Y][Z + surf_normal_Z] + (heat_flux_at_bc / alpha_t); //temperature at wall can be computed from the temperature at the first node and the heat flux at the wall (using Fourier's law)
                            double T_wall_at_bc = (2.0 * heat_flux_at_bc / alpha_t - temperature[X + 2 *surf_normal_X][Y + 2 * surf_normal_Y][Z + 2 * surf_normal_Z] + 4.0 * temperature[X + surf_normal_X][Y + surf_normal_Y][Z + surf_normal_Z]) / 3.0; //temperature at wall can be computed from the temperature at the first node and the heat flux at the wall (using Fourier's law)
                            pop_T[pop_rev][X][Y][Z] = (weight_T[pop_rev] + weight_T[pop]) * T_wall_at_bc - pop_T[pop][X + cx_pop_T[pop]][Y + cy_pop_T[pop]][Z + cz_pop_T[pop]];
                            #endif                            
                            #if STORE_1D_VECTOR
                            size_t id4_pop = idx_4d_vector(pop, X + cx_pop_T[pop], Y + cy_pop_T[pop], Z + cz_pop_T[pop]);
                            size_t id4_rev = idx_4d_vector(pop_rev, X, Y, Z);
                            size_t id3_1   = idx_3d_vector(X + surf_normal_X, Y + surf_normal_Y, Z + surf_normal_Z);
                            size_t id3_2   = idx_3d_vector(X + 2 * surf_normal_X, Y + 2 * surf_normal_Y, Z + 2 * surf_normal_Z);
                            //double T_wall_at_bc = temperature[id3_1] + (heat_flux_at_bc / alpha_t); //First order. temperature at wall can be computed from the temperature at the first node and the heat flux at the wall (using Fourier's law)
                            double T_wall_at_bc = (2.0 * heat_flux_at_bc / alpha_t - temperature[id3_2] + 4.0 * temperature[id3_1]) / 3.0; //Second order.
                            pop_T[id4_rev] = (weight_T[pop_rev] + weight_T[pop]) * T_wall_at_bc - pop_T[id4_pop];
                            #endif
                        }
                    }// End of fixed heat flux
                    //===============================================================//
                }//Z
            }//Y
        }//X
    }//End of for (INT32 ii = 0; ii < 6; ii++)

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: bc_at_solid_points_thermal
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::bc_at_solid_points_thermal()
{
    if (!Flow_With_Particle){
        const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id] + limitBound_left_X;
        const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id] - limitBound_right_X;
        const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id] + limitBound_bottom_Y;
        const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id] - limitBound_top_Y;
        const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id] + limitBound_back_Z;
        const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id] - limitBound_front_Z;

        INT32 Xp, Yp, Zp, Xp2, Yp2, Zp2;
        //#if defined No_curved_boundary_treatment
        for (INT32 X = start_x; X <= end_x; ++X) {
            for (INT32 Y = start_y; Y <= end_y; ++Y) {
                for (INT32 Z = start_z; Z <= end_z; ++Z) {
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    if (is_solid[id3] == FALSE){
                        for (INT32 dir = 1; dir < stencil_directions_T; dir++)	{
                            Xp = (X - cx_pop_T[dir]);
                            Yp = (Y - cy_pop_T[dir]);
                            Zp = (Z - cz_pop_T[dir]);
                            Xp2 = (X + cx_pop_T[dir]);
                            Yp2 = (Y + cy_pop_T[dir]);
                            Zp2 = (Z + cz_pop_T[dir]);
                            size_t id3_p  = idx_3d_vector(Xp, Yp, Zp);
                            size_t id4    = idx_4d_vector(dir, X, Y, Z);
                            size_t id4_p  = idx_4d_vector(bounce_back_pop_T[dir], Xp, Yp, Zp);
                            size_t id4_p2 = idx_4d_vector(dir, Xp2, Yp2, Zp2);
                            if (is_solid[id3_p] == TRUE){
                                if (solid_thermal[id3_p] == 1){//constant temp.
                                    pop_T[id4] = (weight_T[dir] + weight_T[bounce_back_pop_T[dir]]) * solid_temperature - pop_T[id4_p];
                                }
                                if (solid_thermal[id3_p] == 2){//adiabatic
                                    pop_T[id4] = pop_T[id4_p2];
                                }
                            }
                        }
                    }
                    #endif

                    #if !STORE_1D_VECTOR
                    if (is_solid[X][Y][Z] == FALSE){
                        for (INT32 dir = 1; dir < stencil_directions_T; dir++)	{
                            Xp = (X - cx_pop_T[dir]);
                            Yp = (Y - cy_pop_T[dir]);
                            Zp = (Z - cz_pop_T[dir]);
                            Xp2 = (X + cx_pop_T[dir]);
                            Yp2 = (Y + cy_pop_T[dir]);
                            Zp2 = (Z + cz_pop_T[dir]);
                            if (is_solid[Xp][Yp][Zp] == TRUE){
                                if (solid_thermal[Xp][Yp][Zp] == 1){//constant temp.
                                    pop_T[dir][X][Y][Z] = (weight_T[dir] + weight_T[bounce_back_pop_T[dir]])*solid_temperature - pop_T[bounce_back_pop_T[dir]][Xp][Yp][Zp];
                                }
                                if (solid_thermal[Xp][Yp][Zp] == 2){//adiabatic
                                    pop_T[dir][X][Y][Z] = pop_T[dir][Xp2][Yp2][Zp2];
                                }
                            }
                        }
                    }
                    #endif
                }
            }
        }
        //#endif /// End of: No_curved_boundary_treatment
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: momenta_thermal
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function computes the fluid temperature from the populations.
void FlowSolver::momenta_thermal()
{

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    numerator_residual_thermal   = 0;
    denominator_residual_thermal = 0;

    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                double sumT{0};
                #if !STORE_1D_VECTOR
                previous_temperature[X][Y][Z] = temperature[X][Y][Z];
                for (INT32 P = 0; P < stencil_directions_T; ++P)
                {
                    sumT += pop_T[P][X][Y][Z];
                }
                if (Equilibrium_Thermal_1)
                {
                    temperature[X][Y][Z] = sumT;

                }
                if (Equilibrium_Thermal_2)
                {
                    temperature[X][Y][Z] = sumT;
                }
                force_thermal[X][Y][Z] = 0.0;
                numerator_residual_thermal   += pow(temperature[X][Y][Z] - previous_temperature[X][Y][Z], 2);
                denominator_residual_thermal += pow(previous_temperature[X][Y][Z], 2);
                #endif
                #if STORE_1D_VECTOR
                size_t id3 = idx_3d_vector(X, Y, Z);
                previous_temperature[id3] = temperature[id3];
                for (INT32 P = 0; P < stencil_directions_T; ++P)
                {
                    size_t id4_pop = idx_4d_vector(P, X, Y, Z);
                    sumT += pop_T[id4_pop];
                }
                if (Equilibrium_Thermal_1)
                {
                    temperature[id3] = sumT;

                }
                if (Equilibrium_Thermal_2)
                {
                    temperature[id3] = sumT;
                }
                force_thermal[id3] = 0.0;
                numerator_residual_thermal   += pow(temperature[id3] - previous_temperature[id3], 2);
                denominator_residual_thermal += pow(previous_temperature[id3], 2);
                #endif // STORE_1D_VECTOR
            }
        }
    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: force_thermal_2_part_two
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::force_thermal_2_part_two(){

    /*INT32 num_buffer_layers_X = paral.get_num_buffer_layers_X();
    INT32 num_buffer_layers_Y = paral.get_num_buffer_layers_Y();
    INT32 num_buffer_layers_Z = paral.get_num_buffer_layers_Z();*/

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

       /* INT32 start_X, end_X;//Start and End of LBM calculation for processor. (May differ from processor X range, based on the boundary condition)
    start_X = 0;
    end_X = Nx - 1;

    if (X_Wall)
    {
        start_X = 1;
        end_X = Nx - 2;
    }

    if (X_velIn_ZeroGradOut)
    {
        start_X = 1;
        end_X = Nx - 2;
    }
    if (ZeroGradIn_ZeroGradOut)
    {
        start_X = 1;
        end_X = Nx - 2;
    }
    INT32 start_Y, end_Y;
    start_Y = 0;
    end_Y = Ny - 1;

    if (Y_Wall)
    {
        start_Y = 1;
        end_Y = Ny - 2;
    }

    INT32 start_Z, end_Z;
    start_Z = 0;
    end_Z = Nz - 1;

    if (Z_Wall)
    {
        start_Z = 1;
        end_Z = Nz - 2;
    }*/

    for (INT32 X = start_x; X <= end_x; ++X)
    {
        for (INT32 Y = start_y; Y <= end_y; ++Y)
        {
            for (INT32 Z = start_z; Z <= end_z; ++Z)
            {
                #if !STORE_1D_VECTOR
                temperature[X][Y][Z] = temperature[X][Y][Z] + force_thermal[X][Y][Z] / 2.0;
                #endif
                #if STORE_1D_VECTOR
                size_t id3 = idx_3d_vector(X, Y, Z);
                temperature[id3] = temperature[id3] + force_thermal[id3] / 2.0;
                #endif // STORE_1D_VECTOR
            }
        }
    }
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

void FlowSolver::write_fluid_vtk(INT32 time)
{

    /*INT32 num_buffer_layers_X = paral.get_num_buffer_layers_X();
    INT32 num_buffer_layers_Y = paral.get_num_buffer_layers_Y();
    INT32 num_buffer_layers_Z = paral.get_num_buffer_layers_Z();*/

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    if (VTK_BINARY)
    {
        /// Create filename

        stringstream output_filename;
        output_filename << "FlowLB_Results/vtk_fluid/fluid_t" << time << ".vtk";
        ofstream output_file;

        /// Open file
        output_file.open(output_filename.str().c_str());

        /// Write VTK header
        output_file << "# vtk DataFile Version 3.0\n";
        output_file << "fluid_state\n";
        output_file << "BINARY\n";
        output_file << "DATASET RECTILINEAR_GRID\n";
        output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
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
        float ZZ[Nz];

        for (INT32 X = 0; X < Nx; ++X)
        {
            XX[X] = X;
        }
        for (INT32 Y = 0; Y < Ny; ++Y)
        {
            YY[Y] = Y;
        }
        for (INT32 Z = 0; Z < Nz; ++Z)
        {
            ZZ[Z] = Z;
        }

        /// Open file for append
        std::stringstream str_line1;
        str_line1 << "FlowLB_Results/vtk_fluid/fluid_t" << time << ".vtk";
        std::string strstr_line1 = str_line1.str();

        std::ofstream output_file1(strstr_line1, std::ios::app | std::ios::binary);

        //	std::ofstream output_file1("FlowLB_Results/vtk_fluid/fluid_t200.vtk", std::ios::app | std::ios::binary);

        write_bigendian(output_file1, XX, Nx);

        output_file1 << "\nY_COORDINATES " << Ny << " float\n";

        write_bigendian(output_file1, YY, Ny);

        output_file1 << "\nZ_COORDINATES " << Nz << " float\n";

        write_bigendian(output_file1, ZZ, Nz);

        output_file1 << "\nPOINT_DATA " << Nx * Ny * Nz << "\n";

        /// Write density difference
        output_file1 << "SCALARS density_difference float 1\n";
        output_file1 << "LOOKUP_TABLE default\n";

        float mm[1];
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    mm[0] = density[X][Y][Z] - 1;
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    mm[0] = density[id3] - 1;
                    #endif // STORE_1D_VECTOR
                    write_bigendian(output_file1, mm, 1);
                }
            }
        }

        output_file1 << "\nSCALARS temperature float 1\n";
        output_file1 << "LOOKUP_TABLE default\n";
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    mm[0] = temperature[X][Y][Z];
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    mm[0] = temperature[id3];
                    #endif // STORE_1D_VECTOR
                    write_bigendian(output_file1, mm, 1);
                }
            }
        }

        //float nn[2];
        ///// Write velocity
        output_file1 << "VECTORS velocity_vector float\n";
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    mm[0] = velocity[0][X][Y][Z];
                    write_bigendian(output_file1, mm, 1);
                    mm[0] = velocity[1][X][Y][Z];
                    write_bigendian(output_file1, mm, 1);
                    mm[0] = velocity[2][X][Y][Z];
                    write_bigendian(output_file1, mm, 1);
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    mm[0] = velocity[id4_0];
                    write_bigendian(output_file1, mm, 1);
                    mm[0] = velocity[id4_1];
                    write_bigendian(output_file1, mm, 1);
                    mm[0] = velocity[id4_2];
                    write_bigendian(output_file1, mm, 1);
                    #endif // STORE_1D_VECTOR
                }
            }
        }

        output_file1.close();
    }

    if (VTK_ASCII)
    {

        /// Create filename

        stringstream output_filename;
        output_filename << "FlowLB_Results/vtk_fluid/fluid_t" << time << ".vtk";
        ofstream output_file;

        /// Open file

        output_file.open(output_filename.str().c_str());

        /// Write VTK header

        output_file << "# vtk DataFile Version 3.0\n";
        output_file << "fluid_state\n";
        output_file << "ASCII\n";
        output_file << "DATASET RECTILINEAR_GRID\n";

        //output_file << "DIMENSIONS " << Nx << " " << Ny - 2  << " " << Nz - 2 << "\n";
        output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";

        output_file << "X_COORDINATES " << Nx << " float\n";

        for (INT32 X = 0; X < Nx; ++X)
        {
            output_file << X << " ";
        }

        output_file << "\n";
        output_file << "Y_COORDINATES " << Ny << " float\n";

        for (INT32 Y = 0; Y < Ny; ++Y)
        {
            output_file << Y << " ";
        }

        output_file << "\n";
        output_file << "Z_COORDINATES " << Nz << " float\n";

        for (INT32 Z = 0; Z < Nz; ++Z)
        {
            output_file << Z << " ";
        }


        //output_file << "\n";
        //output_file << "Z_COORDINATES " << 1 << " double\n";
        //output_file << 0 << "\n";
        output_file << "POINT_DATA " << Nx * Ny * Nz << "\n";

        /// Write density difference

        output_file << "SCALARS density_difference float 1\n";
        output_file << "LOOKUP_TABLE default\n";

        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file << density[X][Y][Z] - 1 << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    output_file << density[id3] - 1 << "\n";
                    #endif
                }
            }
        }

        /// Write temperature

        output_file << "SCALARS temperature float 1\n";
        output_file << "LOOKUP_TABLE default\n";

        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file << temperature[X][Y][Z] << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    output_file << temperature[id3] << "\n";
                    #endif // STORE_1D_VECTOR
                }
            }
        }

        /// Write velocity

        output_file << "VECTORS velocity_vector float\n";

        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    //#if defined Shan_Chen
                    //				output_file << velocity_x[X][Y][Z] + 0.5 * (force_x[X][Y][Z] + gravity) / density[X][Y][Z] << " " << velocity_y[X][Y][Z] + 0.5 * (force_y[X][Y][Z]) / density[X][Y][Z] << " " << velocity_z[X][Y][Z] + 0.5 * (force_z[X][Y][Z]) / density[X][Y][Z] << "\n";
                    //#endif
                    //#if (defined Guo || defined Luo)
                    #if !STORE_1D_VECTOR
                    output_file << velocity[0][X][Y][Z] << " " << velocity[1][X][Y][Z] << " " << velocity[2][X][Y][Z] << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    output_file << velocity[id4_0] << " " << velocity[id4_1] << " " << velocity[id4_2] << "\n";
                    #endif // STORE_1D_VECTOR

                    //#endif
                }
            }
        }
        /// Close file

        output_file.close();
    }

    if (CASE)
    {
        //unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);// number are shown as 2.0e02 instead of 2.0e002 (necessay for geo and data files). This command only works in Windows

        /////////// Writing CASE file (.case) ///////////////////
        /////////////////////////////////////////////////////////
        /// Create filename

        stringstream output_filename;
        output_filename << "FlowLb_Results/vtk_fluid/fluid" << ".case";
        ofstream output_file;

        /// Open file
        output_file.open(output_filename.str().c_str());

        /// Write VTK header
        output_file << "# EnSight Gold Model\n";
        output_file << "# Created by FlowLB\n\n";
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
        output_file << "number of steps:        " << t_end / t_vtk << "\n";
        output_file << "filename start number:  " << t_vtk << "\n";
        output_file << "filename increment:     " << t_vtk << "\n";
        output_file << "time values:\n";
        for (INT32 itm = t_vtk; itm <= t_end; itm = itm + t_vtk)
        {
            output_file << itm << "\n";

        }
        /// Close file
        output_file.close();
        /////////// End of Writing CASE file (.case) ///////////////
        ////////////////////////////////////////////////////////////

        if (time == t_vtk)  //Because .geo file do not change during iterations
        {

            /////////// Writing GEOMETRY file (.geo) ///////////////////
            ////////////////////////////////////////////////////////////
            stringstream output_filename1;
            output_filename1 << "FlowLB_Results/vtk_fluid/fluid" << ".geo";//Geometry file
            ofstream output_file1;

            /// Open file
            output_file1.open(output_filename1.str().c_str());

            /// Write VTK header
            output_file1 << "EnSight Model Geometry File\n";
            output_file1 << "EnSight 10.1.1\n";
            output_file1 << "node id given\n";
            output_file1 << "element id given\n";
            output_file1 << "extents\n";
            output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << (double)Nx - 1. << "\n";
            output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << (double)Ny - 1. << "\n";
            output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << 0.0 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << (double)Nz - 1. << "\n";
            output_file1 << "part\n";
            output_file1 << std::setw(10) << 1 << "\n";
            output_file1 << "FLUID\n";
            output_file1 << "coordinates\n";
            output_file1 << std::setw(10) << Nx * Ny * Nz << "\n";
            for (INT32 itm = 1; itm <= Nx * Ny * Nz; itm++)
            {
                output_file1 << std::setw(10) << itm << "\n";
            }

            for (INT32 Z = 0; Z < Nz; Z++)
            {
                for (INT32 Y = 0; Y < Ny; Y++)
                {
                    for (INT32 X = 0; X < Nx; X++)
                    {
                        output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << (double)X << "\n";
                    }
                }
            }
            for (INT32 Z = 0; Z < Nz; Z++)
            {
                for (INT32 Y = 0; Y < Ny; Y++)
                {
                    for (INT32 X = 0; X < Nx; X++)
                    {
                        output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << (double)Y << "\n";
                    }
                }
            }

            for (INT32 Z = 0; Z < Nz; Z++)
            {
                for (INT32 Y = 0; Y < Ny; Y++)
                {
                    for (INT32 X = 0; X < Nx; X++)
                    {
                        output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << (double)Z << "\n";
                    }
                }
            }
            output_file1 << "hexa8\n";
            output_file1 << std::setw(10) << (Nx - 1) * (Ny - 1) * (Nz - 1) << "\n";
            for (INT32 itm = 1; itm <= (Nx - 1) * (Ny - 1) * (Nz - 1); itm++)
            {
                output_file1 << std::setw(10) << itm << "\n";
            }

            for (INT32 K = 1; K < Nz; K++)
            {
                for (INT32 J = 1; J < Ny; J++)
                {
                    for (INT32 I = 1; I < Nx; I++)
                    {
                        output_file1 << std::setw(10) << I + (J - 1) * Nx + (K - 1) * Nx * Ny << std::setw(10) << I + 1 + (J - 1) * Nx + (K - 1) * Nx * Ny << std::setw(10) << (I + 1) + (J * Nx) + (K - 1) * Nx * Ny << std::setw(10) << (I)+(J * Nx) + (K - 1) * Nx * Ny << std::setw(10) << I + (J - 1) * Nx + (K)*Nx * Ny << std::setw(10) << I + 1 + (J - 1) * Nx + (K)*Nx * Ny << std::setw(10) << I + 1 + (J)*Nx + (K)*Nx * Ny << std::setw(10) << I + (J)*Nx + (K)*Nx * Ny << "\n";
                    }
                }
            }
            /// Close file
            output_file1.close();
            /////////// End of Writing GEOMETRY file (.geo) ///////////////////
            ///////////////////////////////////////////////////////////////////

        }// End of if (time == t_vtk)


        ///////////// Writing DENSITY_DIFFERENCE file ////////////////////////
        ////////////////////////////////////////////////////////////
        stringstream output_filename6;
        output_filename6 << "FlowLB_Results/vtk_fluid/fluid.density_difference" << setfill('0') << setw(6) << time;
        ofstream output_file6;

        /// Open file
        output_file6.open(output_filename6.str().c_str());

        output_file6 << "density_difference\n";
        output_file6 << "part\n";
        output_file6 << std::setw(10) << 1 << "\n";
        output_file6 << "coordinates\n";
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file6 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << density[X][Y][Z] - 1. << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    output_file6 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << density[id3] - 1. << "\n";
                    #endif // STORE_1D_VECTOR
                }
            }
        }

        /// Close file
        output_file6.close();

        /////////// End of Writing DENSITY_DIFFERENCE file /////////////////
        ////////////////////////////////////////////////////////////

        ///////////// Writing TEMPERATURE file ////////////////////////
        /////////////////////////////////////////////////////////////////////
        stringstream output_filename2;
        output_filename2 << "FlowLB_Results/vtk_fluid/fluid.temperature" << setfill('0') << setw(6) << time;
        ofstream output_file2;

        /// Open file
        output_file2.open(output_filename2.str().c_str());

        output_file2 << "temperature\n";
        output_file2 << "part\n";
        output_file2 << std::setw(10) << 1 << "\n";
        output_file2 << "coordinates\n";

        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file2 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << temperature[X][Y][Z] << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    output_file2 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << temperature[id3] << "\n";
                    #endif // STORE_1D_VECTOR
                }
            }
        }

        /// Close file
        output_file2.close();

        /////////// End of Writing TEMPERATURE file /////////////////
        ////////////////////////////////////////////////////////////

        ///////////// Writing VELOCITY_VECTOR file ////////////////////////
        /////////////////////////////////////////////////////////////////////
        stringstream output_filename5;
        output_filename5 << "FlowLB_Results/vtk_fluid/fluid.VELOC_VECTOR" << setfill('0') << setw(6) << time;
        ofstream output_file5;

        /// Open file
        output_file5.open(output_filename5.str().c_str());

        output_file5 << "VELOC_VECTOR\n";
        output_file5 << "part\n";
        output_file5 << std::setw(10) << 1 << "\n";
        output_file5 << "coordinates\n";

        //#if (defined Guo || defined Luo)
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << velocity[0][X][Y][Z] << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << velocity[id4_0] << "\n";
                    #endif
                }
            }
        }
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << velocity[1][X][Y][Z] << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << velocity[id4_1] << "\n";
                    #endif
                }
            }
        }
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << velocity[2][X][Y][Z] << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id4_2 = idx_4d_vector(1, X, Y, Z);
                    output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << velocity[id4_2] << "\n";
                    #endif
                }
            }
        }
        //#endif
        //#if defined Shan_Chen
        //	for (INT32 Z = 0; Z < Nz; Z++){
        //		for (INT32 Y = 0; Y < Ny; Y++){
        //			for (INT32 X = 0; X < Nx; X++){
        //				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_x[X][Y][Z] + 0.5 * (force_x[X][Y][Z] + gravity) / density[X][Y][Z] << "\n";
        //			}
        //		}
        //	}
        //	for (INT32 Z = 0; Z < Nz; Z++){
        //		for (INT32 Y = 0; Y < Ny; Y++){
        //			for (INT32 X = 0; X < Nx; X++){
        //				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_y[X][Y][Z] + 0.5 * (force_y[X][Y][Z]) / density[X][Y][Z] << "\n";
        //			}
        //		}
        //	}
        //	for (INT32 Z = 0; Z < Nz; Z++){
        //		for (INT32 Y = 0; Y < Ny; Y++){
        //			for (INT32 X = 0; X < Nx; X++){
        //				output_file5 << std::scientific << std::setw(12) << std::setprecision(5) << velocity_z[X][Y][Z] + 0.5 * (force_z[X][Y][Z]) / density[X][Y][Z] << "\n";
        //			}
        //		}
        //	}
        //#endif

        /// Close file
        output_file5.close();

        /////////// End of Writing VELOCITY_VECTOR file /////////////
        ////////////////////////////////////////////////////////////

    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_SOLID_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// Location of solid nodes specified by is_solid=TRUE is written to a vtk file
void FlowSolver::write_solid_vtk(INT32 time, INT32 tot_solid)
{
    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    if (VTK_BINARY ||  VTK_ASCII)
    {
        /*

        /// Create filename

        stringstream output_filename;
        output_filename << "FlowLB_Results/vtk_particle/solid_nodes_" << 1 << "_t" << time << ".vtk";

        ofstream output_file;

        /// Open file

        output_file.open(output_filename.str().c_str());

        /// Write VTK header

        output_file << "# vtk DataFile Version 3.0\n";
        output_file << "fluid_state\n";
        output_file << "ASCII\n";
        output_file << "DATASET RECTILINEAR_GRID\n";

        output_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";

        output_file << "X_COORDINATES " << Nx << " float\n";

        for (INT32 X = 0; X < Nx; ++X)
        {
            output_file << X << " ";
        }

        output_file << "\n";
        output_file << "Y_COORDINATES " << Ny << " float\n";

        for (INT32 Y = 0; Y < Ny; ++Y)
        {
            output_file << Y << " ";
        }

        output_file << "\n";
        output_file << "Z_COORDINATES " << Nz << " float\n";

        for (INT32 Z = 0; Z < Nz; ++Z)
        {
            output_file << Z << " ";
        }

        output_file << "POINT_DATA " << Nx * Ny * Nz << "\n";

        /// Write solid_nodes

        output_file << "SCALARS solid_nodes float 1\n";
        output_file << "LOOKUP_TABLE default\n";


        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    if (is_solid[X][Y][Z])
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    if (is_solid[id3])
                    #endif
                    {
                        output_file << 0 << "\n";
                    }
                    else
                    {
                        output_file << 1 << "\n";
                    }
                }
            }
        }

        /// Close file

        output_file.close();*/


        ////////////////////////////////////
        ////////////////////////////////////
        /// Create filename
        stringstream output_filename1;
        output_filename1 << "FlowLB_Results/vtk_particle/solid_nodes" << "_t" << time << ".vtk";
        ofstream output_file1;

        /// Open file
        output_file1.open(output_filename1.str().c_str());

        /// Write VTK header
        output_file1 << "# vtk DataFile Version 3.0\n";
        output_file1 << "particle_state\n";
        output_file1 << "ASCII\n";
        output_file1 << "DATASET POLYDATA\n";

        /// Write node positions
        output_file1 << "POINTS " << tot_solid << " float\n";
        for (INT32 Z = start_z; Z <= end_z; ++Z)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 X = start_x; X <= end_x; ++X)
                {
                    #if !STORE_1D_VECTOR
                    if (is_solid[X][Y][Z])
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    if (is_solid[id3])
                    #endif
                    {
                        output_file1 << X - num_buffer_layers_X[this_processor_id] << " " << Y - num_buffer_layers_Y[this_processor_id] << " " << Z - num_buffer_layers_Z[this_processor_id] << " \n";
                    }
                }
            }
        }

        /// Write vertices
        output_file1 << "VERTICES 1 " << tot_solid + 1 << "\n";
        output_file1 << tot_solid << " ";
        for (INT32 n = 0; n < tot_solid; ++n)
        {
            output_file1 << n << " ";
        }

        /// Close file
        output_file1.close();
    }
    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Process_Surface_Data: Calculate dimensionless temperature and Nusselt number at the surface
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::Process_Surface_Data(
    INT32 x_local_start,
    INT32 x_local_end,
    INT32 y_local_start,
    INT32 y_local_end,
    INT32 z_local_start,
    INT32 z_local_end,
    FaceDirection face_dir,
    std::vector<double>& theta_values,
    std::vector<double>& nusselt_values)
{
    theta_values.clear();
    nusselt_values.clear();

    //----------------------------------------
    // Constants
    //----------------------------------------
    // Important: Adjust these constants according to the specific case and the definition of Nusselt number in your context. 
    // You may also need to modify the calculation of the "Nusselt number" and "non-dimensional temperature" (70 lines below) based on the specific definition and requirements of your case.
    double T_wall = 1.0;
    double T_bulk = 0.0;
    double T_ref;
    double D_h = (Ny - 1);
    double reference_heat_flux = 0.1;//W/m^2 


    //----------------------------------------
    // Normal direction
    //----------------------------------------

    INT32 nx = 0;
    INT32 ny = 0;
    INT32 nz = 0;

    switch(face_dir)
    {
        case FaceDirection::PlusX:  nx =  1; break;
        case FaceDirection::MinusX: nx = -1; break;

        case FaceDirection::PlusY:  ny =  1; break;
        case FaceDirection::MinusY: ny = -1; break;

        case FaceDirection::PlusZ:  nz =  1; break;
        case FaceDirection::MinusZ: nz = -1; break;
    }

    //----------------------------------------
    // Loop over surface
    //----------------------------------------
    short int x_step = (x_local_end >= x_local_start) ? +1 : -1;
    short int y_step = (y_local_end >= y_local_start) ? +1 : -1;
    short int z_step = (z_local_end >= z_local_start) ? +1 : -1;
#if STORE_1D_VECTOR
        for (INT32 Z = z_local_start; (z_step > 0) ? (Z <= z_local_end) : (Z >= z_local_end); Z += z_step)
        {
            for (INT32 Y = y_local_start; (y_step > 0) ? (Y <= y_local_end) : (Y >= y_local_end); Y += y_step)
            {
                for (INT32 X = x_local_start; (x_step > 0) ? (X <= x_local_end) : (X >= x_local_end); X += x_step)
                {
                T_ref = T_wall;// T_wall or temperature[idx_3d_vector(X, Y, Z)]; You can choose T_ref as the wall temperature or the local temperature at the surface node, depending on your case and definition of Nusselt number.
                //----------------------------------------
                // Using FD gradient
                //----------------------------------------
                INT32 id3_1 = idx_3d_vector(X, Y, Z);

                INT32 id3_2 = idx_3d_vector(X + nx, Y + ny, Z + nz);

                INT32 id3_3 = idx_3d_vector(X + 2 * nx, Y + 2 * ny, Z + 2 * nz);

                double tgrad_fd = (-3.0 * temperature[id3_1] + 4.0 * temperature[id3_2] - temperature[id3_3]) / 2.0;

                double local_nusselt_fd = tgrad_fd / (T_ref - T_bulk) * D_h;

                //----------------------------------------
                // Using NEQ gradient
                //----------------------------------------

                INT32 id4_0 = idx_4d_vector(0, X, Y, Z);
                INT32 id4_1 = idx_4d_vector(1, X, Y, Z);
                INT32 id4_2 = idx_4d_vector(2, X, Y, Z);

                double tgrad_neq = 0.0;

                for (INT32 P = 0; P < stencil_directions_T; ++P)
                {
                    INT32 id4_pop = idx_4d_vector(P, X, Y, Z);

                    double popEq_T =  equilibrium_thermal(temperature[id3_1], density[id3_1], velocity[id4_0], velocity[id4_1], velocity[id4_2], P);

                    tgrad_neq += (nx * cx_pop_T[P] + ny * cy_pop_T[P] + nz * cz_pop_T[P]) * (pop_T[id4_pop] - popEq_T);
                }

                tgrad_neq *= 1.0 / (cs_2_Thermal * tau_t);// See Eq 21 of Tao et al (2020)

                double local_nusselt_neq = tgrad_neq / (T_ref - T_bulk) * D_h;

                //----------------------------------------
                // Dimensionless temperature
                //----------------------------------------
                double theta = (temperature[id3_1] - T_bulk) / D_h / (reference_heat_flux / alpha_t);

                //----------------------------------------
                // Store
                //----------------------------------------

                theta_values.push_back(theta);

                nusselt_values.push_back(local_nusselt_neq);

                //----------------------------------------
                // Debug
                //----------------------------------------

               // std::cout << "X = " << X  << " , Y = " << Y << " , Z = " << Z << " , Nu = " << local_nusselt_neq << std::endl;
            }
        }
    }
    #endif
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_DATA
//-----------------------------------------------------------------------------------------------------------------------------
//
// The FLOW quantities (velocity, density) are written to the disk at each t_data step:

void FlowSolver::write_fluid_data(INT32 time)
{
    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    /// Create filename
	string output_filename("FlowLB_Results/data_fluid");
	output_filename += ".dat";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.c_str(), fstream::app);

	/// Write data
	output_file << time << " "; // time step
    
    //=====================================================================================//
    //                                   TEST CASE 6                                       //
    //=====================================================================================//
    /// Nusselt Number along the vertical wall (for natural convection in a cavity)
	/*double tgrad;
	double nussle_left1 = 0;
	double nussle_left2 = 0;
	for (INT32 Y = 1; Y < Ny - 1; Y++) {
        int x_left_local = 1;
		nussle_left1 += (1 - temperature[idx_3d_vector(x_left_local, Y, 1)]) / 0.5;//0.5 is due to half way
		nussle_left2 += (-3.0*temperature[idx_3d_vector(x_left_local, Y, 1)] + 4.0*temperature[idx_3d_vector(x_left_local + 1, Y, 1)] - temperature[idx_3d_vector(x_left_local + 2, Y, 1)]) / 2.0;
	}
	output_file << nussle_left1 << " " << nussle_left2 << " ";

	double nussle_right1 = 0;
	double nussle_right2 = 0;

	for (INT32 Y = 1; Y < Ny - 1; Y++) {
        int x_right_local = Nx;
		nussle_right1 += -1 * (0 - temperature[idx_3d_vector(x_right_local, Y, 1)]) / 0.5;//0.5 is due to half way
		nussle_right2 += -1 * (-3.0*temperature[idx_3d_vector(x_right_local, Y, 1)] + 4.0*temperature[idx_3d_vector(x_right_local - 1, Y, 1)] - temperature[idx_3d_vector(x_right_local - 2, Y, 1)]) / 2.0;
	}
	output_file << nussle_right1 << " " << nussle_right2 << " ";
    output_file << "\n";   */
    
    // ======== Or use Process_Surface_Data function =========
    /*std::vector<double> theta_left;
    std::vector<double> nu_left;
    INT32 x_local_start = 0 + 1;
    INT32 x_local_end   = 0 + 1;
    INT32 y_local_start = 0 + 1;
    INT32 y_local_end   = Ny - 1 + 1;
    INT32 z_local_start = 1;
    INT32 z_local_end   = 1;
    double avg_nu_left = 0;
    Process_Surface_Data(x_local_start, x_local_end, y_local_start, y_local_end, z_local_start, z_local_end, FaceDirection::PlusX, theta_left, nu_left);
    for (size_t i = 0; i < nu_left.size(); ++i)
    {
        avg_nu_left += nu_left[i];
        //output_file << nu_left[i] << ", ";
        //output_file << theta_left[i] << ", ";
    }
    avg_nu_left /= nu_left.size();
    output_file << avg_nu_left << " ";
    output_file << "\n";*/

    //=====================================================================================//
    //                                   TEST CASE 10                                      //
    //=====================================================================================//
    /// For channel with heat transfer where the Nu=7.54. (Test Case 10)
    /*#if STORE_1D_VECTOR
    double T_wall = 1.0;
    double D_h = 2 * (Ny - 2);//Hydraulic diameter
    double tgrad[Nx];
    size_t id3, id4_velx;
    double local_Nusselt{};
    for (int X = start_x; X <= end_x; X++){
        double T_bulk  = 0;
        double sum_u_T = 0;
        double sum_u   = 0;
        for (int Y = start_y; Y <= end_y; Y++){
            id3      = idx_3d_vector(X, Y, 1);
            id4_velx = idx_4d_vector(0, X, Y, 1);
            sum_u   += velocity[id4_velx];
            sum_u_T += velocity[id4_velx] * temperature[id3];
        }
        T_bulk   = sum_u_T / (sum_u);
        id3      = idx_3d_vector(X, 1, 1);
        tgrad[X] = (T_wall - temperature[id3]) / 0.5;// 0.5 is due two half way bounce back
        //cout << "X = " << X << " , T_bulk = " << T_bulk << " , tgrad[X]= " << tgrad[X] << endl;
        local_Nusselt = tgrad[X] / (T_wall - T_bulk) * D_h;
        output_file << local_Nusselt << " ";
    }
    #endif
	output_file << "\n";
	output_file << "\n";*/

    // ======== Or use of Process_Surface_Data function =========
    /*double T_wall = 1.0;
    size_t id3, id4_velx;
    std::vector<double> theta_bottom;
    std::vector<double> nu_bottom;
    INT32 x_local_start = 0 + 1;
    INT32 x_local_end   = Nx - 1 + 1;
    INT32 y_local_start = 0 + 1;
    INT32 y_local_end   = 0 + 1;
    INT32 z_local_start = 1;
    INT32 z_local_end   = 1;
    Process_Surface_Data(x_local_start, x_local_end, y_local_start, y_local_end, z_local_start, z_local_end, FaceDirection::PlusY, theta_bottom, nu_bottom);
    for (size_t i = 0; i < nu_bottom.size(); ++i)
    {
        int X = i + 1;
        double T_bulk  = 0;
        double sum_u_T = 0;
        double sum_u   = 0;
        for (int Y = start_y; Y <= end_y; Y++){
            id3      = idx_3d_vector(X, Y, 1);
            id4_velx = idx_4d_vector(0, X, Y, 1);
            sum_u   += velocity[id4_velx];
            sum_u_T += velocity[id4_velx] * temperature[id3];
        }
        T_bulk   = sum_u_T / (sum_u);
        double local_Nusselt = nu_bottom[i] / (T_wall - T_bulk) * 2; //2 because D_h = 2*(Ny-2) in this case.
        output_file << local_Nusselt << ", ";
    }
    output_file << "\n";
	output_file << "\n";*/

    //=====================================================================================//
    //                                   TEST CASE 11                                      //
    //=====================================================================================//
    /// For channel with constant heat fluxwhere the Nu=8.235. (Test Case 11)

    // ======== use of Process_Surface_Data function =========
    /*size_t id3, id4_velx;
    std::vector<double> theta_bottom;
    std::vector<double> nu_bottom;
    INT32 x_local_start = 0 + 1;
    INT32 x_local_end   = Nx - 1 + 1;
    INT32 y_local_start = 0 + 1;
    INT32 y_local_end   = 0 + 1;
    INT32 z_local_start = 1;
    INT32 z_local_end   = 1;
    Process_Surface_Data(x_local_start, x_local_end, y_local_start, y_local_end, z_local_start, z_local_end, FaceDirection::PlusY, theta_bottom, nu_bottom);
    for (size_t i = 0; i < nu_bottom.size(); ++i)
    {
        int X = i + 1;
        double T_bulk  = 0;
        double sum_u_T = 0;
        double sum_u   = 0;
        for (int Y = start_y; Y <= end_y; Y++){
            id3      = idx_3d_vector(X, Y, 1);
            id4_velx = idx_4d_vector(0, X, Y, 1);
            sum_u   += velocity[id4_velx];
            sum_u_T += velocity[id4_velx] * temperature[id3];
        }
        T_bulk   = sum_u_T / (sum_u);
        double T_wall = temperature[idx_3d_vector(X, 1, 1)];
        double local_Nusselt = nu_bottom[i] / (T_wall - T_bulk) * 2; //2 because D_h = 2*(Ny-2) in this case.
        output_file << local_Nusselt << ", ";
    }
    output_file << "\n";
	output_file << "\n";*/

    //=====================================================================================//
    //                                   PERMEABILITY                                      //
    //=====================================================================================//
    ///Permeability calculation and storage
    /*double permeability;
    if (gravity > 0)
    {
        double uav = 0.;
        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    #if STORE_1D_VECTOR
                    INT32 id3   = idx_3d_vector(X, Y, Z);
                    INT32 id4_0 = idx_4d_vector(0, X, Y, Z);
                    if (!is_solid[id3])
                    {
                        uav = uav + velocity[id4_0]; //Because gravity is just applied in X direction
                    }
                    #endif
                }
            }
        }
        uav = uav / (Nx * Ny * Nz);
        permeability = uav * (tau - 0.5) / 3.0 / gravity;
        printf("Flow permeability in Lattice units is %f .\n", permeability);
        output_file << permeability << " ";
    }*/
    //======================================================================================//

	output_file << "\n";

    /// Close file
    output_file.close();

    return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: gather_fluid_recovery_data_to_master
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::gather_fluid_recovery_data_to_master(){

    /// Flatten the velocity data of other processors into a 1D buffer
    std::vector<double> send_buffer;

    if (this_processor_id != MASTER)//Important: MASTER should not send its data (that are zero) to itself
    {
        const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
        const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
        const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
        const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
        const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
        const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];
        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    #if !STORE_1D_VECTOR
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            send_buffer.push_back(pop_F[c_i][X][Y][Z]);
                        }
                        send_buffer.push_back(density[X][Y][Z]);
                        send_buffer.push_back(velocity[0][X][Y][Z]);
                        send_buffer.push_back(velocity[1][X][Y][Z]);
                        send_buffer.push_back(velocity[2][X][Y][Z]);
                        send_buffer.push_back(force[0][X][Y][Z]);
                        send_buffer.push_back(force[1][X][Y][Z]);
                        send_buffer.push_back(force[2][X][Y][Z]);
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            send_buffer.push_back(pop_T[c_i][X][Y][Z]);
                        }
                        send_buffer.push_back(temperature[X][Y][Z]);
                        send_buffer.push_back(force_thermal[X][Y][Z]);
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        send_buffer.push_back(pop_F[id4_pop]);
                    }
                    send_buffer.push_back(density[id3]);
                    send_buffer.push_back(velocity[id4_0]);
                    send_buffer.push_back(velocity[id4_1]);
                    send_buffer.push_back(velocity[id4_2]);
                    send_buffer.push_back(force[id4_0]);
                    send_buffer.push_back(force[id4_1]);
                    send_buffer.push_back(force[id4_2]);
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        send_buffer.push_back(pop_T[id4_pop]);
                    }
                    send_buffer.push_back(temperature[id3]);
                    send_buffer.push_back(force_thermal[id3]);
                    #endif
                }
            }
        }
    }

    /// Synchronize before proceeding
    //MPI_Barrier(MPI_COMM_WORLD);

    INT32 send_count = send_buffer.size();
    std::vector<INT32> recv_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);

    /// Gather the sizes of the buffers
    MPI_Gather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Calculate displacements on the master processor
    if (this_processor_id == MASTER) {
        displs[0] = 0;
        for (INT32 i = 1; i < num_processors_total; ++i) {
            displs[i] = displs[i-1] + recv_counts[i-1];
        }
    }

    /// Gather all the data to the master processor
    std::vector<double> recv_buffer;
    if (this_processor_id == MASTER) {
        INT32 total_recv_count = 0;
        for (INT32 count : recv_counts) {
            total_recv_count += count;
        }
        recv_buffer.resize(total_recv_count);
    }

    MPI_Gatherv(send_buffer.data(), send_count, MPI_DOUBLE,
            recv_buffer.data(), recv_counts.data(), displs.data(),
            MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (this_processor_id == MASTER) {
        INT32 index = 0;
        for (INT32 proc = 1; proc < num_processors_total; ++proc) {
            for (INT32 X = start_global_x[proc] + num_buffer_layers_X[this_processor_id]; X <= end_global_x[proc] + num_buffer_layers_X[this_processor_id]; ++X) {
                for (INT32 Y = start_global_y[proc] + num_buffer_layers_Y[this_processor_id]; Y <= end_global_y[proc] + num_buffer_layers_Y[this_processor_id]; ++Y) {
                    for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z[this_processor_id]; Z <= end_global_z[proc] + num_buffer_layers_Z[this_processor_id]; ++Z) {
                        #if !STORE_1D_VECTOR
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            pop_F[c_i][X][Y][Z] = recv_buffer[index++];
                        }
                        density[X][Y][Z]     = recv_buffer[index++];
                        velocity[0][X][Y][Z] = recv_buffer[index++];
                        velocity[1][X][Y][Z] = recv_buffer[index++];
                        velocity[2][X][Y][Z] = recv_buffer[index++];
                        force[0][X][Y][Z]    = recv_buffer[index++];
                        force[1][X][Y][Z]    = recv_buffer[index++];
                        force[2][X][Y][Z]    = recv_buffer[index++];
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            pop_T[c_i][X][Y][Z] = recv_buffer[index++];
                        }
                        temperature[X][Y][Z]   = recv_buffer[index++];
                        force_thermal[X][Y][Z] = recv_buffer[index++];
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            pop_F[id4_pop] = recv_buffer[index++];
                        }
                        density[id3]    = recv_buffer[index++];
                        velocity[id4_0] = recv_buffer[index++];
                        velocity[id4_1] = recv_buffer[index++];
                        velocity[id4_2] = recv_buffer[index++];
                        force[id4_0]    = recv_buffer[index++];
                        force[id4_1]    = recv_buffer[index++];
                        force[id4_2]    = recv_buffer[index++];
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            pop_T[id4_pop] = recv_buffer[index++];
                        }
                        temperature[id3]   = recv_buffer[index++];
                        force_thermal[id3] = recv_buffer[index++];
                        #endif
                    }
                }
            }
        }
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: compress_file
//-----------------------------------------------------------------------------------------------------------------------------
// Function to compress a file
void FlowSolver::compress_file(const std::string& filename) {
    std::string command = "gzip " + filename;
    INT32 result = std::system(command.c_str());
    if (result != 0) {
        std::cerr << "Failed to compress file: " << filename << std::endl;
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: decompress_file
//-----------------------------------------------------------------------------------------------------------------------------
bool FlowSolver::decompress_file(const std::string& compressed_filename, const std::string& decompressed_filename) {
    std::string command = "gzip -dc " + compressed_filename + " > " + decompressed_filename;
    INT32 result = std::system(command.c_str());
    return result == 0;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_FLUID_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
// The FLOW quantities (velocity, density) are written to the disk at each t_recovery step:
void FlowSolver::write_fluid_recovery(INT32 time)
{

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    if (BINARY_Recovery_File)
    {

        std::stringstream str_line1;
        str_line1 << "FlowLB_Results/data_fluid_recovery.dat";
        std::string strstr_line1 = str_line1.str();

        std::ofstream output_file1(strstr_line1, std::ios::out | std::ios::binary);

        float tt[1];
        tt[0] = time;
        //std::cout << tt[0];
        output_file1.write((char*)tt, sizeof(float) * 1);

        float mm[1];

        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    #if !STORE_1D_VECTOR

                    mm[0] = density[X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = velocity[0][X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = velocity[1][X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = velocity[2][X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = temperature[X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);

                    /*mm[0] = force[0][X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = force[1][X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = force[2][X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = force_thermal[X][Y][Z];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        mm[0] = pop_F[c_i][X][Y][Z];
                        output_file1.write((char*)mm, sizeof(float) * 1);
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        mm[0] = pop_T[c_i][X][Y][Z];
                        output_file1.write((char*)mm, sizeof(float) * 1);
                    }*/
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);

                    mm[0] = density[id3];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = velocity[id4_0];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = velocity[id4_1];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = velocity[id4_2];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = temperature[id3];
                    output_file1.write((char*)mm, sizeof(float) * 1);

                    /*mm[0] = force[id4_0];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = force[id4_1];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = force[id4_2];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    mm[0] = force_thermal[id3];
                    output_file1.write((char*)mm, sizeof(float) * 1);
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        mm[0] = pop_F[id4_pop];
                        output_file1.write((char*)mm, sizeof(float) * 1);
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        mm[0] = pop_T[id4_pop];
                        output_file1.write((char*)mm, sizeof(float) * 1);
                    }*/
                    #endif
                }
            }
        }

        output_file1.close();

        std::string compressed_filename = "FlowLB_Results/data_fluid_recovery.dat.gz";
        remove(compressed_filename.c_str());//remove previous compressed file
        compress_file(strstr_line1);
    }


    if (ASCII_Recovery_File)
    {

        /// Create filename
        string output_filename("FlowLB_Results/data_fluid_recovery");
        output_filename += ".dat";
        ofstream output_file;

        /// Open file
        output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

        /// Write data
        output_file << time << " "; // time step is written just as the first number in the recovery file.

        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    #if !STORE_1D_VECTOR
                    output_file << std::fixed << std::setprecision(15) << density[X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << velocity[0][X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << velocity[1][X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << velocity[2][X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << temperature[X][Y][Z] << " ";

                    /*output_file << std::fixed << std::setprecision(15) << force[0][X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << force[1][X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << force[2][X][Y][Z] << " ";
                    output_file << std::fixed << std::setprecision(15) << force_thermal[X][Y][Z] << " ";
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        output_file << std::fixed << std::setprecision(15) << pop_F[c_i][X][Y][Z] << " ";
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        output_file << std::fixed << std::setprecision(15) << pop_T[c_i][X][Y][Z] << " ";
                    }*/

                    output_file << "\n";
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    output_file << std::fixed << std::setprecision(15) << density[id3] << " ";
                    output_file << std::fixed << std::setprecision(15) << velocity[id4_0] << " ";
                    output_file << std::fixed << std::setprecision(15) << velocity[id4_1] << " ";
                    output_file << std::fixed << std::setprecision(15) << velocity[id4_2] << " ";
                    output_file << std::fixed << std::setprecision(15) << temperature[id3] << " ";

                    /*output_file << std::fixed << std::setprecision(15) << force[id4_0] << " ";
                    output_file << std::fixed << std::setprecision(15) << force[id4_1] << " ";
                    output_file << std::fixed << std::setprecision(15) << force[id4_2] << " ";
                    output_file << std::fixed << std::setprecision(15) << force_thermal[id3] << " ";
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        output_file << std::fixed << std::setprecision(15) << pop_F[id4_pop] << " ";
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        output_file << std::fixed << std::setprecision(15) << pop_T[id4_pop] << " ";
                    }*/

                    output_file << "\n";
                    #endif // STORE_1D_VECTOR
                }
            }
        }

        /// Close file
        output_file.close();
    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: READ_FLUID_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
// The FLOW quantities (density, 3 velocities, temperature) are READ from the disk to initialize the flow.
void FlowSolver::read_fluid_recovery(INT32& Step)
{
    std::cout << "\nReading data from Recovery File. Please be patient ... "<< std::endl;

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    if (BINARY_Recovery_File)
    {
        std::string compressed_filename = "FlowLB_Results/data_fluid_recovery.dat.gz";
        std::string decompressed_filename = "FlowLB_Results/data_fluid_recovery.dat";

        // Decompress the file
        if (!decompress_file(compressed_filename, decompressed_filename)) {
            std::cerr << "Failed to decompress file or there is no .gz file" << std::endl;
        }

        std::ifstream input_file1(decompressed_filename, std::ios::in | std::ios::binary);

        float tt[1];
        if (input_file1.is_open())
        {
            input_file1.read((char*)tt, sizeof(float));
        }else{
            std::cerr << "Failed to open data_fluid_recovery.dat" << std::endl;
        }

        Step = static_cast<INT32>(tt[0]);
        std::cout << "Recovery step is: " << Step << endl;

        if (input_file1.is_open())
        {
            for (INT32 X = start_x; X <= end_x; ++X)
            {
                for (INT32 Y = start_y; Y <= end_y; ++Y)
                {
                    for (INT32 Z = start_z; Z <= end_z; ++Z)
                    {
                        #if !STORE_1D_VECTOR
                        input_file1.read((char*)tt, sizeof(float));
                        density[X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        velocity[0][X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        velocity[1][X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        velocity[2][X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        temperature[X][Y][Z] = tt[0];

                        /*input_file1.read((char*)tt, sizeof(float));
                        force[0][X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        force[1][X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        force[2][X][Y][Z] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        force_thermal[X][Y][Z] = tt[0];
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            input_file1.read((char*)tt, sizeof(float));
                            pop_F[c_i][X][Y][Z] = tt[0];
                        }
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            input_file1.read((char*)tt, sizeof(float));
                            pop_T[c_i][X][Y][Z] = tt[0];
                        }*/
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        input_file1.read((char*)tt, sizeof(float));
                        density[id3] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        velocity[id4_0] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        velocity[id4_1] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        velocity[id4_2] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        temperature[id3] = tt[0];

                        /*input_file1.read((char*)tt, sizeof(float));
                        force[id4_0] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        force[id4_1] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        force[id4_2] = tt[0];
                        input_file1.read((char*)tt, sizeof(float));
                        force_thermal[id3] = tt[0];
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            input_file1.read((char*)tt, sizeof(float));
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            pop_F[id4_pop] = tt[0];
                        }
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            input_file1.read((char*)tt, sizeof(float));
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            pop_T[id4_pop] = tt[0];
                        }*/
                        #endif
                    }
                }
            }
            input_file1.close();

            remove(decompressed_filename.c_str());
        }

    }


    if (ASCII_Recovery_File)
    {
        /// Create filename
        string input_filename("FlowLB_Results/data_fluid_recovery");
        input_filename += ".dat";
        ifstream input_file;//File is open for READING

        /// Open file
        input_file.open(input_filename.c_str());

        /// Read data
        input_file >> Step;      // get first number from the file. The first number is: "saved recovery step"
        std::cout << "Recovery step is: " << Step << endl;

        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    if (!input_file.eof())
                    {
                        #if !STORE_1D_VECTOR
                        input_file >> density[X][Y][Z];
                        input_file >> velocity[0][X][Y][Z];
                        input_file >> velocity[1][X][Y][Z];
                        input_file >> velocity[2][X][Y][Z];
                        input_file >> temperature[X][Y][Z];

                        /*input_file >> force[0][X][Y][Z];
                        input_file >> force[1][X][Y][Z];
                        input_file >> force[2][X][Y][Z];
                        input_file >> force_thermal[X][Y][Z];
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            //	input_file >> pop_old[c_i][X][Y][Z];
                            input_file >> pop_F[c_i][X][Y][Z];
                        }
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            input_file >> pop_T[c_i][X][Y][Z];
                        }*/
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        input_file >> density[id3];
                        input_file >> velocity[id4_0];
                        input_file >> velocity[id4_1];
                        input_file >> velocity[id4_2];
                        input_file >> temperature[id3];

                        /*input_file >> force[id4_0];
                        input_file >> force[id4_1];
                        input_file >> force[id4_2];
                        input_file >> force_thermal[id3];
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            input_file >> pop_F[id4_pop];
                        }
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            input_file >> pop_T[id4_pop];
                        }*/
                        #endif
                    }
                }
            }
        }

        /// Close file
        input_file.close();
    }
    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: read_fluid_initialize_data
//-----------------------------------------------------------------------------------------------------------------------------
// The FLOW quantities (density, 3 velocities, temperature) are READ from the disk to initialize the flow.
void FlowSolver::read_fluid_initialize_data(INT32& Step,
                                            const INT32 NxA,
                                            const INT32 NyA,
                                            const INT32 NzA)
{
    std::cout << "\nReading initialize data and interpolating ...\n";

    const INT32 start_x = start_local_x[this_processor_id] + num_buffer_layers_X[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id]   - num_buffer_layers_X[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id] + num_buffer_layers_Y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id]   - num_buffer_layers_Y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id] + num_buffer_layers_Z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id]   - num_buffer_layers_Z[this_processor_id];

    INT32 NxB = end_x - start_x + 1;
    INT32 NyB = end_y - start_y + 1;
    INT32 NzB = end_z - start_z + 1;

    size_t sizeA = (size_t)NxA * NyA * NzA;

    std::vector<float> rhoA(sizeA);
    std::vector<float> uA(sizeA);
    std::vector<float> vA(sizeA);
    std::vector<float> wA(sizeA);
    std::vector<float> TA(sizeA);

    auto idA = [&](INT32 x, INT32 y, INT32 z)
    {
        return (size_t)x + NxA * ((size_t)y + NyA * (size_t)z);
    };

    /// open file
    std::string input_filename("data_fluid_initialize.dat");
    std::ifstream input_file(input_filename.c_str());

    if(!input_file.is_open())
    {
        std::cerr << "Failed to open initialize file\n";
        exit(1);
    }

    /// read time
    float time_value;
    input_file >> time_value;

    Step = static_cast<INT32>(time_value);

    std::cout << "Reading time from initialize file = " << time_value << std::endl;

    /// read file grid data
    for(INT32 z=0; z<NzA; ++z)
    {
        for(INT32 y=0; y<NyA; ++y)
        {
            for(INT32 x=0; x<NxA; ++x)
            {
            size_t id = idA(x,y,z);

            input_file >> rhoA[id];
            input_file >> uA[id];
            input_file >> vA[id];
            input_file >> wA[id];
            input_file >> TA[id];
            //cout << x << " " << y << " " << z << " " << rhoA[id] << " " << uA[id] << " " << vA[id] << " " << wA[id] << " " << TA[id] <<endl;
            }
        }
    }

    input_file.close();

    /// scaling factors
    double rx = double(NxA-1) / double(NxB-1);
    double ry = double(NyA-1) / double(NyB-1);
    double rz = double(NzA-1) / double(NzB-1);

    /// trilinear interpolation
    // Trilinear interpolation function
    auto interp = [&](const std::vector<float>& f,
                      INT32 x0, INT32 x1,
                      INT32 y0, INT32 y1,
                      INT32 z0, INT32 z1,
                      double dx, double dy, double dz)
    {
        float c000 = f[idA(x0,y0,z0)];
        float c100 = f[idA(x1,y0,z0)];
        float c010 = f[idA(x0,y1,z0)];
        float c110 = f[idA(x1,y1,z0)];
        float c001 = f[idA(x0,y0,z1)];
        float c101 = f[idA(x1,y0,z1)];
        float c011 = f[idA(x0,y1,z1)];
        float c111 = f[idA(x1,y1,z1)];

        float c00 = c000*(1-dx) + c100*dx;
        float c10 = c010*(1-dx) + c110*dx;
        float c01 = c001*(1-dx) + c101*dx;
        float c11 = c011*(1-dx) + c111*dx;

        float c0 = c00*(1-dy) + c10*dy;
        float c1 = c01*(1-dy) + c11*dy;

        return c0*(1-dz) + c1*dz;
    };

    // Interpolate onto simulation grid
    for(INT32 X=start_x; X<=end_x; ++X)
    {
        for(INT32 Y=start_y; Y<=end_y; ++Y)
        {
            for(INT32 Z=start_z; Z<=end_z; ++Z)
            {
            // Map simulation indices to source grid coordinates
            double xA = (NxB>1) ? (double)(X - start_x) * (NxA-1) / (NxB-1) : 0.0;
            double yA = (NyB>1) ? (double)(Y - start_y) * (NyA-1) / (NyB-1) : 0.0;
            double zA = (NzB>1) ? (double)(Z - start_z) * (NzA-1) / (NzB-1) : 0.0;

            // Neighboring indices in source grid
            INT32 x0 = std::max(0, std::min((INT32)floor(xA), NxA-1));
            INT32 y0 = std::max(0, std::min((INT32)floor(yA), NyA-1));
            INT32 z0 = std::max(0, std::min((INT32)floor(zA), NzA-1));

            INT32 x1 = std::min(x0+1, NxA-1);
            INT32 y1 = std::min(y0+1, NyA-1);
            INT32 z1 = std::min(z0+1, NzA-1);

            double dx = xA - x0;
            double dy = yA - y0;
            double dz = zA - z0;

            // Interpolated values
            float rho = interp(rhoA,x0,x1,y0,y1,z0,z1,dx,dy,dz);
            float u   = interp(uA  ,x0,x1,y0,y1,z0,z1,dx,dy,dz);
            float v   = interp(vA  ,x0,x1,y0,y1,z0,z1,dx,dy,dz);
            float w   = interp(wA  ,x0,x1,y0,y1,z0,z1,dx,dy,dz);
            float T   = interp(TA  ,x0,x1,y0,y1,z0,z1,dx,dy,dz);

            #if !STORE_1D_VECTOR
            density[X][Y][Z]     = rho;
            velocity[0][X][Y][Z] = u;
            velocity[1][X][Y][Z] = v;
            velocity[2][X][Y][Z] = w;
            temperature[X][Y][Z] = T;
            #endif
            #if STORE_1D_VECTOR
            size_t id3  = idx_3d_vector(X, Y, Z);
            size_t id40 = idx_4d_vector(0, X, Y, Z);
            size_t id41 = idx_4d_vector(1, X, Y, Z);
            size_t id42 = idx_4d_vector(2, X, Y, Z);

            density[id3]     = rho;
            velocity[id40]   = u;
            velocity[id41]   = v;
            velocity[id42]   = w;
            temperature[id3] = T;
            //cout << X << " " << Y << " " << Z << ": " << rho << " " << velocity[id40] << " " << velocity[id41] << " " << velocity[id42] << " " << T << endl;
            #endif
            }
        }
    }

    std::cout << "Interpolation completed.\n";
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: scatter_fluid_recovery_data_from_master
// Share fluid recovery data from MASTER to other processors
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::scatter_fluid_recovery_data_from_master(INT32& Step)
{
    //INT32 processor_id = paral.get_processor_id();
    if (this_processor_id == MASTER) {std::cout << "Sharing recovery data from MASTER to all processors ..."<< std::endl;}

    std::vector<double> send_buffer;

    MPI_Bcast(&Step, 1, MPI_INT, MASTER, MPI_COMM_WORLD);// Share recovery time step with other processors

    /// Flatten the fluid data into a 1D buffer (send_buffer)
     if (this_processor_id == MASTER) {
        INT32 index = 0;
        for (INT32 proc = 0; proc < num_processors_total; ++proc) {
            for (INT32 X = start_global_x[proc] + num_buffer_layers_X[this_processor_id]; X <= end_global_x[proc] + num_buffer_layers_X[this_processor_id]; ++X) {
                for (INT32 Y = start_global_y[proc] + num_buffer_layers_Y[this_processor_id]; Y <= end_global_y[proc] + num_buffer_layers_Y[this_processor_id]; ++Y) {
                    for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z[this_processor_id]; Z <= end_global_z[proc] + num_buffer_layers_Z[this_processor_id]; ++Z) {
                        #if !STORE_1D_VECTOR
                        send_buffer.push_back(density[X][Y][Z]);
                        send_buffer.push_back(velocity[0][X][Y][Z]);
                        send_buffer.push_back(velocity[1][X][Y][Z]);
                        send_buffer.push_back(velocity[2][X][Y][Z]);
                        send_buffer.push_back(temperature[X][Y][Z]);

                        send_buffer.push_back(force[0][X][Y][Z]);
                        send_buffer.push_back(force[1][X][Y][Z]);
                        send_buffer.push_back(force[2][X][Y][Z]);
                        send_buffer.push_back(force_thermal[X][Y][Z]);
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            send_buffer.push_back(pop_F[c_i][X][Y][Z]);
                        }
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            send_buffer.push_back(pop_T[c_i][X][Y][Z]);
                        }
                        #endif
                        #if STORE_1D_VECTOR

                        size_t id3   = idx_3d_vector(X, Y, Z);
                        size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                        size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                        size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                        send_buffer.push_back(density[id3]);
                        send_buffer.push_back(velocity[id4_0]);
                        send_buffer.push_back(velocity[id4_1]);
                        send_buffer.push_back(velocity[id4_2]);
                        send_buffer.push_back(temperature[id3]);

                        send_buffer.push_back(force[id4_0]);
                        send_buffer.push_back(force[id4_1]);
                        send_buffer.push_back(force[id4_2]);
                        send_buffer.push_back(force_thermal[id3]);
                        for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                        {
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            send_buffer.push_back(pop_F[id4_pop]);
                        }
                        for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                        {
                            size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                            send_buffer.push_back(pop_T[id4_pop]);
                        }
                        #endif
                    }
                }
            }
        }
        //std::cout << "Master: send_buffer size = " << send_buffer.size() << std::endl;
    }

    /// Calculate send counts and displacements for each processor
    std::vector<INT32> send_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);
    if (this_processor_id == MASTER)//only MASTER calculates send_counts and displs to save time
    {
        INT32 index = 0;
        for (INT32 proc = 0; proc < num_processors_total; ++proc)
        {
            INT32 x_size = end_global_x[proc] - start_global_x[proc] + 1;
            INT32 y_size = end_global_y[proc] - start_global_y[proc] + 1;
            INT32 z_size = end_global_z[proc] - start_global_z[proc] + 1;
            send_counts[proc] = x_size * y_size * z_size * (stencil_directions_F + 7 + stencil_directions_T + 2);
            //send_counts[proc] = x_size * y_size * z_size * (4 + 1);//density, 3 velocities, temperature
            displs[proc] = index;
            index += send_counts[proc];
        }
    }

    /// Broadcast send_counts and displs from the master to all processors
    MPI_Bcast(send_counts.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Scatter the data from the master processor to all processors
    INT32 recv_count = send_counts[this_processor_id];
    std::vector<double> recv_buffer(recv_count);
    //std::cout << "Processor " << this_processor_id << ": recv_count = " << recv_count << std::endl;

    MPI_Scatterv(this_processor_id == MASTER ? send_buffer.data() : nullptr, send_counts.data(), displs.data(), MPI_DOUBLE,
              recv_buffer.data(), send_counts[this_processor_id], MPI_DOUBLE,
              MASTER, MPI_COMM_WORLD);

    /// Reshape the received data back into 3D arrays on each processor
    if (this_processor_id != MASTER) {
        INT32 x_size = end_global_x[this_processor_id] - start_global_x[this_processor_id] + 1;
        INT32 y_size = end_global_y[this_processor_id] - start_global_y[this_processor_id] + 1;
        INT32 z_size = end_global_z[this_processor_id] - start_global_z[this_processor_id] + 1;

        INT32 index = 0;
        for (INT32 X = num_buffer_layers_X[this_processor_id]; X < x_size + num_buffer_layers_X[this_processor_id]; ++X) {
            for (INT32 Y = num_buffer_layers_Y[this_processor_id]; Y < y_size + num_buffer_layers_Y[this_processor_id]; ++Y) {
                for (INT32 Z = num_buffer_layers_Z[this_processor_id]; Z < z_size + num_buffer_layers_Z[this_processor_id]; ++Z) {
                    #if !STORE_1D_VECTOR

                    density[X][Y][Z]     = recv_buffer[index++];
                    velocity[0][X][Y][Z] = recv_buffer[index++];
                    velocity[1][X][Y][Z] = recv_buffer[index++];
                    velocity[2][X][Y][Z] = recv_buffer[index++];
                    temperature[X][Y][Z] = recv_buffer[index++];

                    force[0][X][Y][Z]    = recv_buffer[index++];
                    force[1][X][Y][Z]    = recv_buffer[index++];
                    force[2][X][Y][Z]    = recv_buffer[index++];
                    force_thermal[X][Y][Z] = recv_buffer[index++];
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        pop_F[c_i][X][Y][Z] = recv_buffer[index++];
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        pop_T[c_i][X][Y][Z] = recv_buffer[index++];
                    }
                    #endif
                    #if STORE_1D_VECTOR

                    size_t id3   = idx_3d_vector(X, Y, Z);
                    size_t id4_0 = idx_4d_vector(0, X, Y, Z);
                    size_t id4_1 = idx_4d_vector(1, X, Y, Z);
                    size_t id4_2 = idx_4d_vector(2, X, Y, Z);
                    density[id3]     = recv_buffer[index++];
                    velocity[id4_0] = recv_buffer[index++];
                    velocity[id4_1] = recv_buffer[index++];
                    velocity[id4_2] = recv_buffer[index++];
                    temperature[id3]   = recv_buffer[index++];

                    force[id4_0]    = recv_buffer[index++];
                    force[id4_1]    = recv_buffer[index++];
                    force[id4_2]    = recv_buffer[index++];
                    force_thermal[id3] = recv_buffer[index++];
                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        pop_F[id4_pop] = recv_buffer[index++];
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i)
                    {
                        size_t id4_pop = idx_4d_vector(c_i, X, Y, Z);
                        pop_T[id4_pop] = recv_buffer[index++];
                    }
                    #endif
                }
            }
        }
    }


    /*const INT32 start_x = start_local_x[this_processor_id];
    const INT32 end_x   = end_local_x[this_processor_id];
    const INT32 start_y = start_local_y[this_processor_id];
    const INT32 end_y   = end_local_y[this_processor_id];
    const INT32 start_z = start_local_z[this_processor_id];
    const INT32 end_z   = end_local_z[this_processor_id];

    // Use the equilibrium populations corresponding to the initialized fluid density and velocity.
    if (this_processor_id == 1){
        for (INT32 X = start_x; X <= end_x; ++X)
        {
            for (INT32 Y = start_y; Y <= end_y; ++Y)
            {
                for (INT32 Z = start_z; Z <= end_z; ++Z)
                {
                    if (velocity[idx_4d_vector(0, X, Y, Z)] > 0.00001)
                    cout << this_processor_id << " " << X << " " << Y << " " << Z << " " << velocity[idx_4d_vector(0, X, Y, Z)] << " " << velocity[idx_4d_vector(1, X, Y, Z)] << endl;
                }
            }
        }
    }*/


}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: read_solid_from_file
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::read_solid_from_file(const std::string& filename, const ParallelMPI& parallelLB){
    // Read geometry from the file (geometry.txt). This file must be created by MATLAB program: image_to_matrix.m
    // Before using program specify Nx and Ny at the beginning of program according to geometry size.
    INT32 sol;
    string output_filename(filename);
    ifstream output_file;
    /// Open file
    output_file.open(output_filename.c_str());
    /// Read geometry data
    for (INT32 X = 0; X < Nx; X++) {
        for (INT32 Y = 0; Y < Ny; Y++) {
            for (INT32 Z = 0; Z < Nz; Z++) {
                if (!output_file.eof()) {
                    output_file >> sol;
                    if (sol == 0) { // 0 in geometry.txt means SOLID
                        const INT32 local_X = parallelLB.get_Xlocal_from_Xglobal(X, this_processor_id);
                        const INT32 local_Y = parallelLB.get_Ylocal_from_Yglobal(Y, this_processor_id);
                        const INT32 local_Z = parallelLB.get_Zlocal_from_Zglobal(Z, this_processor_id);

                        #if !STORE_1D_VECTOR
                            is_solid[local_X][local_Y][local_Z] = TRUE;
                        #endif
                        #if STORE_1D_VECTOR
                            size_t id3 = idx_3d_vector(local_X, local_Y, local_Z);
                            is_solid[id3] = TRUE;
                        #endif // STORE_1D_VECTOR

                    }
                }
            }
        }
    }

    /// Close geometry file
    output_file.close();
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: create_solid_from_dimensions
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::create_solid_from_dimensions(const ParallelMPI& parallelLB){

    double cir_center_x = 100 + 1;
    double cir_center_y = 120 + 1;
    double cir_radius   = 40;

    for (INT32 X = 0 + 1; X < Nx + 1; ++X) {//Define Solid nodes
    	for (INT32 Y = 0 + 1; Y < Ny + 1; ++Y) {
            for (INT32 Z = 1; Z < 2; ++Z) {
                if ((SQ(X - cir_center_x) + SQ(Y - cir_center_y)) <= SQ(cir_radius))	{
                    #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    is_solid[id3] = TRUE;
                    #endif
                    #if !STORE_1D_VECTOR
                    is_solid[X][Y][Z] = TRUE;
                    #endif
                }
            }
    	}
    }


    /*const INT32 start_solid_global_x = 500;//User defined
    const INT32 end_solid_global_x   = 570;//User defined
    const INT32 start_solid_global_y = 80;//User defined
    const INT32 end_solid_global_y   = 120;//User defined
    const INT32 start_solid_global_z = 0;//User defined
    const INT32 end_solid_global_z   = 0;//User defined

    const INT32 start_solid_local_x  = parallelLB.get_Xlocal_from_Xglobal(start_solid_global_x, this_processor_id);
    const INT32 end_solid_local_x    = parallelLB.get_Xlocal_from_Xglobal(end_solid_global_x,   this_processor_id);
    const INT32 start_solid_local_y  = parallelLB.get_Ylocal_from_Yglobal(start_solid_global_y, this_processor_id);
    const INT32 end_solid_local_y    = parallelLB.get_Ylocal_from_Yglobal(end_solid_global_y,   this_processor_id);
    const INT32 start_solid_local_z  = parallelLB.get_Zlocal_from_Zglobal(start_solid_global_z, this_processor_id);
    const INT32 end_solid_local_z    = parallelLB.get_Zlocal_from_Zglobal(end_solid_global_z,   this_processor_id);

    if (start_solid_global_x > Nx || end_solid_global_x > Nx || start_solid_global_y > Ny || end_solid_global_y > Ny || start_solid_global_z > Nz || end_solid_global_z > Nz) {
        throw std::runtime_error("\n\n------------- Error: Solid dimensions in create_solid_from_dimensions() is larger than Nx,y,z or smaller than 0! ----------------- \n");
    }

    for (INT32 X = start_solid_local_x; X <= end_solid_local_x; ++X) {//Solid nodes
        for (INT32 Y = start_solid_local_y; Y <= end_solid_local_y; ++Y) {
            for (INT32 Z = start_solid_local_z; Z <= end_solid_local_z; ++Z) {
                #if !STORE_1D_VECTOR
                    is_solid[X][Y][Z] = TRUE;
                #endif
                #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    is_solid[id3] = TRUE;
                #endif
            }
        }
    }*/


}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: count_solid
//-----------------------------------------------------------------------------------------------------------------------------
INT32 FlowSolver::count_solid(const ParallelMPI& parallelLB){
    const INT32 start_solid_global_x = 0;
    const INT32 end_solid_global_x   = this->Nx;
    const INT32 start_solid_global_y = 0;
    const INT32 end_solid_global_y   = this->Ny;
    const INT32 start_solid_global_z = 0;
    const INT32 end_solid_global_z   = this->Nz;

    const INT32 start_solid_local_x  = parallelLB.get_Xlocal_from_Xglobal(start_solid_global_x, this_processor_id);
    const INT32 end_solid_local_x    = parallelLB.get_Xlocal_from_Xglobal(end_solid_global_x,   this_processor_id);
    const INT32 start_solid_local_y  = parallelLB.get_Ylocal_from_Yglobal(start_solid_global_y, this_processor_id);
    const INT32 end_solid_local_y    = parallelLB.get_Ylocal_from_Yglobal(end_solid_global_y,   this_processor_id);
    const INT32 start_solid_local_z  = parallelLB.get_Zlocal_from_Zglobal(start_solid_global_z, this_processor_id);
    const INT32 end_solid_local_z    = parallelLB.get_Zlocal_from_Zglobal(end_solid_global_z,   this_processor_id);

    INT32 total_solid{0};
    for (INT32 X = start_solid_local_x; X <= end_solid_local_x; ++X) {//Solid nodes
        for (INT32 Y = start_solid_local_y; Y <= end_solid_local_y; ++Y) {
            for (INT32 Z = start_solid_local_z; Z <= end_solid_local_z; ++Z) {
                #if !STORE_1D_VECTOR
                    if (is_solid[X][Y][Z] == TRUE) {
                        total_solid = total_solid + 1;
                        if (Solid_Thermal_BC == 1){//Do not change
                            solid_thermal[X][Y][Z] = 1;//Solid_Constant_Temperature.
                        }
                        if (Solid_Thermal_BC == 2){//Do not change
                            solid_thermal[X][Y][Z] = 2;//Solid_Adiabatic
                        }
                        if (Solid_Thermal_BC == 3){//Do not change
                            solid_thermal[X][Y][Z] = 3;//Solid_Conjugate
                        }
                    }
                #endif
                #if STORE_1D_VECTOR
                    size_t id3 = idx_3d_vector(X, Y, Z);
                    if (is_solid[id3] == TRUE) {
                        total_solid = total_solid + 1;
                        if (Solid_Thermal_BC == 1){//Do not change
                            solid_thermal[id3] = 1;//Solid_Constant_Temperature.
                        }
                        if (Solid_Thermal_BC == 2){//Do not change
                            solid_thermal[id3] = 2;//Solid_Adiabatic
                        }
                        if (Solid_Thermal_BC == 3){//Do not change
                            solid_thermal[id3] = 3;//Solid_Conjugate
                        }
                    }
                #endif // STORE_1D_VECTOR
            }
        }
    }
    return total_solid;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: scatter_is_solid_points_from_master
// Share is_solid data from MASTER to other processors
//-----------------------------------------------------------------------------------------------------------------------------
void FlowSolver::scatter_is_solid_points_from_master()
{
    std::vector<double> send_buffer;

  //  MPI_Bcast(&Step, 1, MPI_INT, MASTER, MPI_COMM_WORLD);// Share recovery time step with other processors

    /// Flatten the is_solid data into a 1D buffer (send_buffer)
     if (this_processor_id == MASTER) {
        std::cout << "\nScattering is_solid and solid_thermal from MASTER to SLAVES." << std::endl;
        INT32 index = 0;
        for (INT32 proc = 0; proc < num_processors_total; ++proc) {
            for (INT32 X = start_global_x[proc] + num_buffer_layers_X[this_processor_id]; X <= end_global_x[proc] + num_buffer_layers_X[this_processor_id]; ++X) {
                for (INT32 Y = start_global_y[proc] + num_buffer_layers_Y[this_processor_id]; Y <= end_global_y[proc] + num_buffer_layers_Y[this_processor_id]; ++Y) {
                    for (INT32 Z = start_global_z[proc] + num_buffer_layers_Z[this_processor_id]; Z <= end_global_z[proc] + num_buffer_layers_Z[this_processor_id]; ++Z) {
                        #if !STORE_1D_VECTOR
                        send_buffer.push_back(is_solid[X][Y][Z]);
                        send_buffer.push_back(solid_thermal[X][Y][Z]);
                        #endif
                        #if STORE_1D_VECTOR
                        size_t id3 = idx_3d_vector(X, Y, Z);
                        send_buffer.push_back(is_solid[id3]);
                        send_buffer.push_back(solid_thermal[id3]);
                        #endif
                    }
                }
            }
        }
        //std::cout << "Master: send_buffer size = " << send_buffer.size() << std::endl;
    }

    /// Calculate send counts and displacements for each processor
    std::vector<INT32> send_counts(num_processors_total, 0);
    std::vector<INT32> displs(num_processors_total, 0);
    if (this_processor_id == MASTER)//only MASTER calculates send_counts and displs to save time
    {
        INT32 index = 0;
        for (INT32 proc = 0; proc < num_processors_total; ++proc)
        {
            INT32 x_size = end_global_x[proc] - start_global_x[proc] + 1;
            INT32 y_size = end_global_y[proc] - start_global_y[proc] + 1;
            INT32 z_size = end_global_z[proc] - start_global_z[proc] + 1;
            //std::cout << "Proc " << proc << ": " << x_size << " " << y_size << " " << z_size << std::endl;
            send_counts[proc] = x_size * y_size * z_size * 2;//2 means is_solid and solid_thermal
            displs[proc] = index;
            index += send_counts[proc];
        }
    }

    /// Broadcast send_counts and displs from the master to all processors
    MPI_Bcast(send_counts.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(displs.data(), num_processors_total, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Scatter the data from the master processor to all processors
    INT32 recv_count = send_counts[this_processor_id];
    std::vector<double> recv_buffer(recv_count);

    MPI_Scatterv(this_processor_id == MASTER ? send_buffer.data() : nullptr, send_counts.data(), displs.data(), MPI_DOUBLE,
              recv_buffer.data(), send_counts[this_processor_id], MPI_DOUBLE,
              MASTER, MPI_COMM_WORLD);

    /// Reshape the received data back into 3D arrays on each processor
    if (this_processor_id != MASTER) {
        INT32 x_size = end_global_x[this_processor_id] - start_global_x[this_processor_id] + 1;
        INT32 y_size = end_global_y[this_processor_id] - start_global_y[this_processor_id] + 1;
        INT32 z_size = end_global_z[this_processor_id] - start_global_z[this_processor_id] + 1;

        INT32 index = 0;
        for (INT32 X = num_buffer_layers_X[this_processor_id]; X < x_size + num_buffer_layers_X[this_processor_id]; ++X) {
            for (INT32 Y = num_buffer_layers_Y[this_processor_id]; Y < y_size + num_buffer_layers_Y[this_processor_id]; ++Y) {
                for (INT32 Z = num_buffer_layers_Z[this_processor_id]; Z < z_size + num_buffer_layers_Z[this_processor_id]; ++Z) {
                    #if !STORE_1D_VECTOR
                    is_solid[X][Y][Z] = recv_buffer[index++];
                    solid_thermal[X][Y][Z] = recv_buffer[index++];
                    #endif
                    #if STORE_1D_VECTOR
                    size_t id3    = idx_3d_vector(X, Y, Z);
                    is_solid[id3] = recv_buffer[index++];
                    solid_thermal[id3] = recv_buffer[index++];
                    #endif
                }
            }
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: check_residual
//-----------------------------------------------------------------------------------------------------------------------------
bool FlowSolver::check_residual()
{
    bool converged = false;

    double global_numerator_residual_flow{0};
    double global_denominator_residual_flow{0};
    double global_numerator_residual_thermal{0};
    double global_denominator_residual_thermal{0};
    double global_residual_calculated_flow{0};
    double global_residual_calculated_thermal{0};

    MPI_Reduce(&numerator_residual_flow,      &global_numerator_residual_flow,      1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&denominator_residual_flow,    &global_denominator_residual_flow,    1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&numerator_residual_thermal,   &global_numerator_residual_thermal,   1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&denominator_residual_thermal, &global_denominator_residual_thermal, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    global_residual_calculated_flow    = sqrt(global_numerator_residual_flow    / (global_denominator_residual_flow    + 1e-12));
    global_residual_calculated_thermal = sqrt(global_numerator_residual_thermal / (global_denominator_residual_thermal + 1e-12));

    if (this_processor_id == MASTER) {
        if (!Flow_With_Thermal_Effect) {
            if (global_residual_calculated_flow < residual_flow) {
                converged = true;
            }
        }
        if (Flow_With_Thermal_Effect) {
            if (global_residual_calculated_flow < residual_flow && global_residual_calculated_thermal < residual_thermal) {
                converged = true;
            }
        }
    }
    return converged;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Destructor
//-----------------------------------------------------------------------------------------------------------------------------
FlowSolver::~FlowSolver()
{
    //dtor
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Copy Constructor
//-----------------------------------------------------------------------------------------------------------------------------
FlowSolver::FlowSolver(const FlowSolver& other) = default; // = default is needed because I initialized all settings and parameters as initializer list in the main constructor and they all were "const"

