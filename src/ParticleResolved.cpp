#include "ParticleResolved.h"
#include "SimulationParameters.h"
#include "SimulationSettings.h"
#include "FlowSolver.h"
#include "ParallelMPI.h"
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Constructor
//-----------------------------------------------------------------------------------------------------------------------------

ParticleResolved::ParticleResolved(const SimulationParameters& paramsLB, const SimulationSettings& settingsLB, const ParallelMPI& paralLB):
    ///Settings from here
    Flow_With_Thermal_Effect(settingsLB.get_bool_setting("Flow_With_Thermal_Effect")),
    X_Periodic(settingsLB.get_bool_setting("X_Periodic")),
    X_Wall(settingsLB.get_bool_setting("X_Wall")),
    Y_Periodic(settingsLB.get_bool_setting("Y_Periodic")),
    Y_Wall(settingsLB.get_bool_setting("Y_Wall")),
    Z_Periodic(settingsLB.get_bool_setting("Z_Periodic")),
    Z_Wall(settingsLB.get_bool_setting("Z_Wall")),
    VTK_ASCII(settingsLB.get_bool_setting("VTK_ASCII")),
    VTK_BINARY(settingsLB.get_bool_setting("VTK_BINARY")),
    CASE(settingsLB.get_bool_setting("CASE")),
    STATIONARY_SPHERE(settingsLB.get_bool_setting("STATIONARY_SPHERE")),
    STATIONARY_SPHEROID(settingsLB.get_bool_setting("STATIONARY_SPHEROID")),
    MOVING_SPHERE(settingsLB.get_bool_setting("MOVING_SPHERE")),
    MOVING_SPHEROID(settingsLB.get_bool_setting("MOVING_SPHEROID")),
    SDF(settingsLB.get_bool_setting("SDF")),
    MDF(settingsLB.get_bool_setting("MDF")),
    Without_added_mass(settingsLB.get_bool_setting("Without_added_mass")),
    With_added_mass(settingsLB.get_bool_setting("With_added_mass")),
    Corrected_Radius(settingsLB.get_bool_setting("Corrected_Radius")),
    Two_Point_Delta(settingsLB.get_bool_setting("Two_Point_Delta")),
    Four_Point_Delta(settingsLB.get_bool_setting("Four_Point_Delta")),
    Kernel_Delta(settingsLB.get_bool_setting("Kernel_Delta")),
    Modified_Two_Point_Delta(settingsLB.get_bool_setting("Modified_Two_Point_Delta")),
    Sphere_Collision_1_Feng(settingsLB.get_bool_setting("Sphere_Collision_1_Feng")),
    Sphere_Collision_2_Feng_Verlet(settingsLB.get_bool_setting("Sphere_Collision_2_Feng_Verlet")),
    Sphere_Collision_3_Feng_Lub_Verlet(settingsLB.get_bool_setting("Sphere_Collision_3_Feng_Lub_Verlet")),
    Spheroid_Collision_1_Feng(settingsLB.get_bool_setting("Spheroid_Collision_1_Feng")),
    Spheroid_Collision_2_Feng_Verlet(settingsLB.get_bool_setting("Spheroid_Collision_2_Feng_Verlet")),
    Spheroid_Collision_3_Feng_Lub_Verlet(settingsLB.get_bool_setting("Spheroid_Collision_3_Feng_Lub_Verlet")),
    Particle_Fixed_Temperature(settingsLB.get_bool_setting("Particle_Fixed_Temperature")),
    Particle_Varying_Temperature(settingsLB.get_bool_setting("Particle_Varying_Temperature")),
    Random_Arrangement_OFF(settingsLB.get_bool_setting("Random_Arrangement_OFF")),
    Random_Arrangement_ON(settingsLB.get_bool_setting("Random_Arrangement_ON")),
    ///2D Settings
    STATIONARY_CYLINDER(settingsLB.get_bool_setting("STATIONARY_CYLINDER")),
    STATIONARY_ELLIPSE(settingsLB.get_bool_setting("STATIONARY_ELLIPSE")),
    MOVING_CYLINDER(settingsLB.get_bool_setting("MOVING_CYLINDER")),
    MOVING_ELLIPSE(settingsLB.get_bool_setting("MOVING_ELLIPSE")),
    Circle_Collision_1(settingsLB.get_bool_setting("Circle_Collision_1")),
    Circle_Collision_2(settingsLB.get_bool_setting("Circle_Collision_2")),
    Circle_Collision_3(settingsLB.get_bool_setting("Circle_Collision_3")),
    Ellipse_Collision_1(settingsLB.get_bool_setting("Ellipse_Collision_1")),
    Ellipse_Collision_2(settingsLB.get_bool_setting("Ellipse_Collision_2")),
    Ellipse_Collision_3(settingsLB.get_bool_setting("Ellipse_Collision_3")),
    Ellipse_Collision_4(settingsLB.get_bool_setting("Ellipse_Collision_4")),

    ///Parameters from here
    t_end(paramsLB.get_t_end()),
    t_vtk(paramsLB.get_t_vtk()),
    Nx(paramsLB.get_Nx()),
    Ny(paramsLB.get_Ny()),
    Nz(paramsLB.get_Nz()),
    tau(paramsLB.get_tau()),
    mdfstep(paramsLB.get_mdfstep()),
    num_particles(paramsLB.get_num_particles()),
    particle_radius(paramsLB.get_particle_radius()),
    particle_center_x(paramsLB.get_particle_center_x()),
    particle_center_y(paramsLB.get_particle_center_y()),
    particle_center_z(paramsLB.get_particle_center_z()),
    particle_gravity_x(paramsLB.get_particle_gravity_x()),
    particle_gravity_y(paramsLB.get_particle_gravity_y()),
    particle_gravity_z(paramsLB.get_particle_gravity_z()),
    particle_den_ratio(paramsLB.get_particle_den_ratio()),
    particle_offset_x(paramsLB.get_particle_offset_x()),
    particle_offset_y(paramsLB.get_particle_offset_y()),
    particle_offset_z(paramsLB.get_particle_offset_z()),
    particle_temperature(paramsLB.get_particle_intial_temperature()),
    particle_heat_source(paramsLB.get_particle_heat_source()),
    particle_specific_heat_ratio(paramsLB.get_particle_specific_heat_ratio()),
    particle_radius_2(paramsLB.get_particle_radius_2()),
    particle_num_nodes(paramsLB.get_particle_num_nodes()),
    particle_num_strip(paramsLB.get_particle_num_strip()),
    particle_cir_factor(paramsLB.get_particle_cir_factor()),
    ellipse_init_ang(paramsLB.get_ellipse_init_ang()),
    particle_r_buffer_verlet(paramsLB.get_particle_r_buffer_verlet()),
    particle_r_cutoff_verlet(paramsLB.get_particle_r_cutoff_verlet()),
    particle_cij(paramsLB.get_particle_cij()),
    particle_ep(paramsLB.get_particle_ep()),
    particle_Ep(paramsLB.get_particle_Ep()),
    particle_epw(paramsLB.get_particle_epw()),
    particle_Epw(paramsLB.get_particle_Epw()),
    particle_th(paramsLB.get_particle_th()),
    particle_thw(paramsLB.get_particle_thw()),
    particle_th_lub(paramsLB.get_particle_th_lub()),
    particle_thw_lub(paramsLB.get_particle_thw_lub()),
    stencil_directions_F(paramsLB.get_stencil_directions_F()),
    stencil_directions_T(paramsLB.get_stencil_directions_T()),
    ///Parameters from ParallelMPI
    num_processors_for_particles_only(paralLB.get_num_processors_for_particles_only()),
    this_processor_id(paralLB.get_processor_id()),
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
    length_with_buffer_X(paralLB.get_length_with_buffer_X()),
    length_with_buffer_Y(paralLB.get_length_with_buffer_Y()),
    length_with_buffer_Z(paralLB.get_length_with_buffer_Z()),
    num_buffer_layers_X(paralLB.get_num_buffer_layers_X()),
    num_buffer_layers_Y(paralLB.get_num_buffer_layers_Y()),
    num_buffer_layers_Z(paralLB.get_num_buffer_layers_Z())
{
    MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_processors);// Get the number of mpi processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_processor_rank);  // Get the rank of the mpi process

    //set_simulation_settings(settings);
    set_simulation_parameters(paramsLB);
    init_kernel_table();// Construct  kernel tables
    if (mpi_processor_rank == 0){
        distribute_particles_mpi(num_mpi_processors);/// Not used in the current version
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: idx_3d_vector
//-----------------------------------------------------------------------------------------------------------------------------
INT32 ParticleResolved::idx_3d_vector(INT32 x, INT32 y, INT32 z)
{
    return (x * this->length_with_buffer_Y[this_processor_id] + y) * this->length_with_buffer_Z[this_processor_id] + z;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: idx_4d_vector
//-----------------------------------------------------------------------------------------------------------------------------
INT32 ParticleResolved::idx_4d_vector(INT32 q, INT32 x, INT32 y, INT32 z)
{
    return ((q * this->length_with_buffer_X[this_processor_id] + x) * this->length_with_buffer_Y[this_processor_id] + y) * this->length_with_buffer_Z[this_processor_id] + z;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: set_simulation_settings
//-----------------------------------------------------------------------------------------------------------------------------
/*void ParticleResolved::set_simulation_settings(const SimulationSettings& settingsLB)
{
    Flow_With_Particle = settingsLB.get_bool_setting("Flow_With_Particle");
    Flow_With_Thermal_Effect = settingsLB.get_bool_setting("Flow_With_Thermal_Effect");
    etc.....

}*/
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: init_kernel_table()
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::init_kernel_table()
{
    kernel_table.resize(Kernel_Table_Size);
    for (INT32 i = 0; i < Kernel_Table_Size; ++i)
    {
        double r = i * Kernel_Resolution;
        if (Kernel_Delta){
            kernel_table[i] = (r <= 2.0) ? 0.25 * (1.0 + cos(M_PI * r / 2.0)) : 0.0;
        }
        else if(Four_Point_Delta)
        {
            if (r < 1.0)
                kernel_table[i] = 0.125 * (3.0 - 2.0 * r + sqrt(1 + 4 * r - 4 * r * r));
            else if (r >= 1 && r < 2.0)
                kernel_table[i] = 0.125 * (5.0 - 2.0 * r - sqrt(-7 + 12 * r - 4 * r * r));
            else
                kernel_table[i] = 0.0;
        }
        else if (Two_Point_Delta)
        {
            kernel_table[i] = (r <= 1.0) ? (1.0 - r) : 0.0;
        }
        else if (Modified_Two_Point_Delta)
        {
            if (r <= 0.5)
            {
                kernel_table[i] = 3. / 4. - SQ(r);
            }
            else if (r > 0.5 && r <= 1.5)
            {
                kernel_table[i] = 9. / 8. - 3. * r / 2. + SQ(r) / 2.;
            }
            else kernel_table[i] = 0.0;
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: kernel_weight()
//-----------------------------------------------------------------------------------------------------------------------------
double ParticleResolved::kernel_weight(double dist)
{
    double r = std::abs(dist);
    if (r > Kernel_Max_Dist)
        return 0.0;
    INT32 idx = static_cast<INT32>(r / Kernel_Resolution);
    return kernel_table[idx];
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: set_simulation_parameters
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::set_simulation_parameters(const SimulationParameters& params)
{
    /*domain_dimension_F = params.get_domain_dimension_F();
    stencil_directions_F = params.get_stencil_directions_F();
    domain_dimension_T = params.get_domain_dimension_T();
    stencil_directions_T = params.get_stencil_directions_T();
    Nx = params.get_Nx();
    Ny = params.get_Ny();
    Nz = params.get_Nz();
    etc....
*/
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: distribute_particles_mpi
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::distribute_particles_mpi(INT32 num_proc) {
    // Calculate the base number of particles per processor
    INT32 base_particles_per_processor = num_particles / num_proc;

    // Calculate the number of processors that will get extra particles
    INT32 num_processors_with_extra = num_particles % num_proc;

    // Distribute particles among processors
    num_particles_per_mpi_processor.resize(num_proc);
    for (INT32 i = 0; i < num_proc; ++i) {
        num_particles_per_mpi_processor[i] = base_particles_per_processor;
        if (i < num_processors_with_extra) {
            num_particles_per_mpi_processor[i]++; // Add one extra particle to each processor until the remainder is exhausted
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: initialize_particles
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::initialize_particles()
{
    if (MOVING_SPHERE || STATIONARY_SPHERE || STATIONARY_CYLINDER || MOVING_CYLINDER)
    {
        if (particle_r_cutoff_verlet < (2 * particle_radius + std::max(particle_th_lub, particle_th)))
        {
            particle_r_cutoff_verlet = 2 * particle_radius + std::max(particle_th_lub, particle_th) + 1;
        }
        if (particle_r_buffer_verlet < (particle_r_cutoff_verlet + particle_radius))
        {
            particle_r_buffer_verlet = (particle_r_cutoff_verlet + particle_radius) + 1;
        }
    }
    if (MOVING_SPHEROID || STATIONARY_SPHEROID || STATIONARY_ELLIPSE || MOVING_ELLIPSE)
    {
        if (particle_r_cutoff_verlet < (2 * std::max(particle_radius, particle_radius_2) + std::max(particle_th_lub, particle_th)))
        {
            particle_r_cutoff_verlet = 2 * std::max(particle_radius, particle_radius_2) + std::max(particle_th_lub, particle_th) + 1;
        }
        if (particle_r_buffer_verlet < (particle_r_cutoff_verlet + std::max(particle_radius, particle_radius_2)))
        {
            particle_r_buffer_verlet = (particle_r_cutoff_verlet + std::max(particle_radius, particle_radius_2)) + 1;
        }
    }

    triangleset1.resize(24000, 0);
    triangleset2.resize(24000, 0);
    triangleset3.resize(24000, 0);

    // Resize verlet_list to have num_particles elements, each with size 0
    verlet_list.resize(num_particles);

    if (STATIONARY_SPHERE || STATIONARY_SPHEROID || STATIONARY_CYLINDER || STATIONARY_ELLIPSE)
    {
        margin = t_end + 1;//Do not change this line
    }

    uncorrected_particle_radius   = particle_radius;//Do not change this line
    uncorrected_particle_radius_2 = particle_radius_2;//Do not change this line
    if (Corrected_Radius){
        corrected_particle_radius   = pow((pow(particle_radius,   3) + pow(particle_radius   - 1.0, 3)) / 2.0, 1.0 / 3.0);
        corrected_particle_radius_2 = pow((pow(particle_radius_2, 3) + pow(particle_radius_2 - 1.0, 3)) / 2.0, 1.0 / 3.0);
        particle_radius    = corrected_particle_radius;
        particle_radius_2  = corrected_particle_radius_2;
    }

    //Calculate particle_num_nodes
    if (STATIONARY_CYLINDER || STATIONARY_ELLIPSE || MOVING_CYLINDER || MOVING_ELLIPSE)
    {
        // Already read from ParametersLB.txt in the Constructor
    }
    if (MOVING_SPHERE || STATIONARY_SPHERE)
    {
        particle_num_nodes = INT32(M_PI / 3. * (3 * SQ(2 * uncorrected_particle_radius) + 1)); // Uhlmann (2005) recommendation for number of surface nodes for sphere
    }

    if (MOVING_SPHEROID || STATIONARY_SPHEROID)
    {
        particle_num_nodes = 6000; // Important. This parameter is just for initial adjustment.
        //The actual particle_num_nodes for SPHEROID will be specified later and used to "resize" the node elements
    }

    // Resize the vector
    INT32 num_verlet_entries = 0; // will change later
    verlet_list.resize(num_particles, std::vector<INT32>(num_verlet_entries));

    // Resize particle vector
    particle_list.resize(num_particles);

    // Initialize each particle
    for (INT32 i = 0; i < num_particles; ++i)
    {
        particle_list[i].node.resize(particle_num_nodes);//May resize later for Spheroid

        if (STATIONARY_CYLINDER || MOVING_CYLINDER)
        {
            particle_list[i].radius             = particle_radius;// Initialize radius
            particle_list[i].radius_2           = particle_radius;// Initialize radius_2
            particle_list[i].parvol             = M_PI * SQ(uncorrected_particle_radius);//Particle volume based on uncorrected_radius
            particle_list[i].particle_perimeter = 2.0 * M_PI * particle_radius;//Perimeter
            particle_list[i].particle_mass      = (particle_den_ratio * 1.0) * particle_list[i].parvol;//Mass
            particle_list[i].mom_inertia_x      = 0;//Moment of inertia
            particle_list[i].mom_inertia_y      = 0;
            particle_list[i].mom_inertia_z      = 0.5 * particle_list[i].particle_mass * SQ(uncorrected_particle_radius);//Mass moment of inertia
        }
        if (STATIONARY_ELLIPSE || MOVING_ELLIPSE)
        {
            particle_list[i].radius             = particle_radius;// Initialize radius
            particle_list[i].radius_2           = particle_radius_2;// Initialize radius_2
            particle_list[i].particle_perimeter = M_PI * (3 * (particle_radius + particle_radius_2) - sqrt((3 * particle_radius + particle_radius_2) * (particle_radius +  3 * particle_radius_2))) ;//Perimeter
            particle_list[i].parvol             = M_PI * uncorrected_particle_radius * uncorrected_particle_radius_2;//Particle volume based on uncorrected_radius
            particle_list[i].particle_mass      = (particle_den_ratio * 1.0) * particle_list[i].parvol;
            particle_list[i].mom_inertia_x      = 0;//Moment of inertia
            particle_list[i].mom_inertia_y      = 0;
            particle_list[i].mom_inertia_z      = 0.25 * particle_list[i].particle_mass * (SQ(uncorrected_particle_radius) + SQ(uncorrected_particle_radius_2));//Mass moment of inertia
        }

        if (MOVING_SPHERE || STATIONARY_SPHERE)
        {
            particle_list[i].radius             = particle_radius;// Initialize radius
            particle_list[i].radius_2           = particle_radius;// Initialize radius_2
            particle_list[i].parvol             = 4.0 / 3.0 * M_PI * pow(uncorrected_particle_radius , 3);//Particle volume based on uncorrected_radius
            particle_list[i].particle_mass      = (particle_den_ratio * 1.0) * particle_list[i].parvol;//Particle mass
            particle_list[i].particle_perimeter = M_PI * (12.0 * SQ(uncorrected_particle_radius) + 1.0) / 3.0;	//Eq. 2.14 Lucci et. al. J. Fluid Mech. (2010), vol. 650, pp. 555.
            particle_list[i].mom_inertia_x      = 2.0 / 5.0 * particle_list[i].particle_mass * SQ(uncorrected_particle_radius);//Moment of inertia
            particle_list[i].mom_inertia_y      = 2.0 / 5.0 * particle_list[i].particle_mass * SQ(uncorrected_particle_radius);
            particle_list[i].mom_inertia_z      = 2.0 / 5.0 * particle_list[i].particle_mass * SQ(uncorrected_particle_radius);
        }
        if (MOVING_SPHEROID || STATIONARY_SPHEROID)
        {
            particle_list[i].radius        = particle_radius;// Initialize radius
            particle_list[i].radius_2      = particle_radius_2;// Initialize radius_2
            particle_list[i].parvol        = 4.0 / 3.0 * M_PI * uncorrected_particle_radius * SQ(uncorrected_particle_radius_2);//Particle volume based on uncorrected_radius and uncorrected_radius_2
            particle_list[i].particle_mass = (particle_den_ratio * 1.0) * particle_list[i].parvol;
            double aa{}, cc{}, hh{}, zz{};
            cc = uncorrected_particle_radius;
            aa = uncorrected_particle_radius_2;
            hh = sqrt(std::max((aa - cc), (cc - aa)) * (aa + cc) / pow(cc, 4));
            zz = uncorrected_particle_radius * hh;
            if (cc > aa)
            {
                particle_list[i].particle_perimeter = 2 * M_PI * aa * 2 * (asin(zz) / 2. + (zz * sqrt(1 - SQ(zz))) / 2.) / hh;
            }
            else
            {
                particle_list[i].particle_perimeter = 2 * M_PI * aa * 2 * (asinh(zz) / 2. + (zz * sqrt(1 + SQ(zz))) / 2.) / hh;
            }
            particle_list[i].mom_inertia_x = 1.0 / 5.0 * particle_list[i].particle_mass * (SQ(uncorrected_particle_radius_2) + SQ(uncorrected_particle_radius_2));
            particle_list[i].mom_inertia_y = 1.0 / 5.0 * particle_list[i].particle_mass * (SQ(uncorrected_particle_radius)   + SQ(uncorrected_particle_radius_2));
            particle_list[i].mom_inertia_z = 1.0 / 5.0 * particle_list[i].particle_mass * (SQ(uncorrected_particle_radius)   + SQ(uncorrected_particle_radius_2));
        }

        // Resize center vector
        particle_list[i].center.resize(1);
        particle_list[i].center[0].temperature = particle_temperature;
        particle_list[i].center[0].temp_part = particle_temperature; //Initial particle temperature
        particle_list[i].center[0].Qu0 = cos(0.5 * particle_list[i].center[0].euler_teta) * cos(0.5 * (particle_list[i].center[0].euler_fi + particle_list[i].center[0].euler_sai));
        particle_list[i].center[0].Qu1 = sin(0.5 * particle_list[i].center[0].euler_teta) * cos(0.5 * (particle_list[i].center[0].euler_fi - particle_list[i].center[0].euler_sai));
        particle_list[i].center[0].Qu2 = sin(0.5 * particle_list[i].center[0].euler_teta) * sin(0.5 * (particle_list[i].center[0].euler_fi - particle_list[i].center[0].euler_sai));
        particle_list[i].center[0].Qu3 = cos(0.5 * particle_list[i].center[0].euler_teta) * sin(0.5 * (particle_list[i].center[0].euler_fi + particle_list[i].center[0].euler_sai));
    }

    /// ***********************************
    /// POSITION OF 1st PARTICLE
    /// ***********************************
    if (num_particles > 0)
    {
        particle_list[0].center[0].x = particle_center_x;
        particle_list[0].center[0].y = particle_center_y;
        particle_list[0].center[0].z = particle_center_z;

        if (STATIONARY_CYLINDER || MOVING_CYLINDER){
            particle_list[0].center[0].z = 0;
            for (INT32 n = 0; n < particle_num_nodes; ++n) {
                particle_list[0].node[n].x = particle_list[0].center[0].x + particle_radius * sin(2. * M_PI * (double)n / particle_num_nodes);
                particle_list[0].node[n].y = particle_list[0].center[0].y + particle_radius * cos(2. * M_PI * (double)n / particle_num_nodes);
                particle_list[0].node[n].z = 0.0;
                particle_list[0].node[n].area_node = 2.0 * M_PI * particle_radius / particle_num_nodes;//Perimeter divided by number of nodes for CIRCLE
                if (particle_list[0].node[n].x > Nx - 0.001 || particle_list[0].node[n].x < 0.001){
                    cout << "ERROR: Particle " << 1 << " position is out of the domain";
                    return;
                }
                if (particle_list[0].node[n].y > Ny - 0.001 || particle_list[0].node[n].y < 0.001){
                    cout << "ERROR: Particle " << 1 << " position is out of the domain";
                    return;
                }
           	}
		}

        if (STATIONARY_ELLIPSE || MOVING_ELLIPSE) {
            particle_list[0].center[0].z = 0;
            //Node points distribution
            for (INT32 n = 0; n < particle_num_nodes; ++n) {
                particle_list[0].node[n].x = particle_list[0].center[0].x + particle_radius   * cos(2. * M_PI * (double)n / particle_num_nodes);
                particle_list[0].node[n].y = particle_list[0].center[0].y + particle_radius_2 * sin(2. * M_PI * (double)n / particle_num_nodes);
                particle_list[0].node[n].z = 0.0;
                if (particle_list[0].node[n].x > Nx - 0.001 || particle_list[0].node[n].x < 0.001){
                    cout << "ERROR: Particle " << 1 << " position is out of the domain";
                    return;
                }
                if (particle_list[0].node[n].y > Ny - 0.001 || particle_list[0].node[n].y < 0.001){
                    cout << "ERROR: Particle " << 1 << " position is out of the domain";
                    return;
                }
            }
            //Ellipse node area calculation
            double dist1, dist2;
            double sum_perimeter;//Total perimeter of ellipse
            sum_perimeter = 0.0;
            for (INT32 n = 0; n < particle_num_nodes; ++n) {
                dist1 = sqrt(SQ(particle_list[0].node[(n + 1) % particle_num_nodes].x - particle_list[0].node[n].x) + SQ(particle_list[0].node[(n + 1) % particle_num_nodes].y - particle_list[0].node[n].y));
                dist2 = sqrt(SQ(particle_list[0].node[n].x - particle_list[0].node[(n - 1 + particle_num_nodes) % particle_num_nodes].x) + SQ(particle_list[0].node[n].y - particle_list[0].node[(n - 1 + particle_num_nodes) % particle_num_nodes].y));
                particle_list[0].node[n].area_node = (dist1 + dist2) / 2.;
                sum_perimeter += particle_list[0].node[n].area_node;
            }
            //cout << sum_perimeter <<endl;
            ///Initial orientation of ellipse
            double x_pos_1, y_pos_1;
            particle_list[0].center[0].tetaZ = M_PI / 2. - ellipse_init_ang;
            for (INT32 n = 0; n < particle_num_nodes; ++n) {
                x_pos_1 = particle_list[0].node[n].x - particle_list[0].center[0].x;
                y_pos_1 = particle_list[0].node[n].y - particle_list[0].center[0].y;
                particle_list[0].node[n].x = particle_list[0].center[0].x + x_pos_1 * cos(ellipse_init_ang) + y_pos_1 * sin(ellipse_init_ang);
                particle_list[0].node[n].y = particle_list[0].center[0].y - x_pos_1 * sin(ellipse_init_ang) + y_pos_1 * cos(ellipse_init_ang);
                particle_list[0].node[n].z = 0.0;
            }
        }

        // Distributing points on Sphere surface, Now only for particle index = 0 (first particle).
        double ck, tetak, fik, fikold;//parameters for the uniform distribution of points on the sphere
        if (MOVING_SPHERE ||  STATIONARY_SPHERE)
        {
            // Distributing points on a sphere according to:
            //Eq. 2.10 , 2.11 Lucci et. al. J. Fluid Mech. (2010), vol. 650, pp. 555.
            for (INT32 k = 0; k < particle_num_nodes; ++k)
            {
                ck = -1.0 + 2.0 * (double)k / (double)(particle_num_nodes - 1.0);
                tetak = acos(ck);
                if (k == 0 || k == particle_num_nodes - 1)
                {
                    fik = 0.0;
                    fikold = 0.0;
                }
                else
                {
                    fik = fikold + 3.6 / sqrt(particle_num_nodes) / sqrt(1.0 - SQ(ck));
                }
                fikold = fik;
                particle_list[0].node[k].x = particle_list[0].center[0].x + particle_radius * cos(fik) * sin(tetak);
                particle_list[0].node[k].y = particle_list[0].center[0].y + particle_radius * sin(fik) * sin(tetak);
                particle_list[0].node[k].z = particle_list[0].center[0].z + particle_radius * cos(tetak);
                //Eq. 2.14 Lucci et. al. J. Fluid Mech. (2010), vol. 650, pp. 555.
                particle_list[0].node[k].area_node = M_PI * (12.0 * SQ(particle_radius) + 1.0) / 3.0 / (double)particle_num_nodes; // area belonging to a node
                if (particle_list[0].node[k].x > Nx - 0.001 || particle_list[0].node[k].x < 0.001)
                {
                    printf("ERROR: Particle %d position is out of the domain \n", 1);
                    break;
                }
                if (particle_list[0].node[k].y > Ny - 0.001 || particle_list[0].node[k].y < 0.001)
                {
                    printf("ERROR: Particle %d position is out of the domain \n", 1);
                    break;
                }
                if (particle_list[0].node[k].z > Nz - 0.001 || particle_list[0].node[k].z < 0.001)
                {
                    printf("ERROR: Particle %d position is out of the domain \n", 1);
                    break;
                }
            }
        }

        if (MOVING_SPHEROID || STATIONARY_SPHEROID)
        {
            // Distributing points on a spheroid.
            double aa{}, cc{}, hh{}, zz{};
            cc = particle_radius;
            aa = particle_radius_2;
            hh = sqrt(std::max((aa - cc), (cc - aa)) * (aa + cc) / pow(cc, 4));
            zz = particle_radius * hh;
            INT32 kkk{}, mmm{};
            double perimeter{};
            perimeter = M_PI * (3.0 * (particle_radius + particle_radius_2) - sqrt((3.0 * particle_radius + particle_radius_2) * (particle_radius + 3.0 * particle_radius_2)));
            //num_strip = (INT32)(perimeter / 2.0 * particle_radius / particle_radius_2) + 2;//Total number of parallel strips (including top and bottom points)

            //double hei[num_strip];//height of each strip
            vector <double> hei(particle_num_strip);//height of each strip

            vector<vector<INT32>> nodes_in_hei(particle_num_strip);//It tells us at each strip height wich nodes are located (node number). Important: Just for VTK file
            // Resize each inner vector to have 2000 elements
            for (INT32 i = 0; i < particle_num_strip; ++i)
            {
                nodes_in_hei[i].resize(2000);
            }
            /*INT32** nodes_in_hei;//It tells us at each strip height wich nodes are located (node number). Important: Just for VTK file
            nodes_in_hei = new INT32*[num_strip];
            for (INT32 nps = 0; nps < num_strip; ++nps) {
                nodes_in_hei[nps] = new INT32[2000];//Certainly the number of nodes on one strip is less than 2000
            }*/

            vector <INT32> po_per_str(particle_num_strip);//Number of points per each strip

            //Important: In the following configuration, the particle_radius is initially along the X-axis and particle_radius_2 is initially along the Y and Z-axis. You can change it later by adjusting initial euler angles
            kkk = 0;
            particle_list[0].node[kkk].x = particle_list[0].center[0].x + particle_radius;
            particle_list[0].node[kkk].y = particle_list[0].center[0].y;
            particle_list[0].node[kkk].z = particle_list[0].center[0].z;
            kkk = kkk + 1;//kkk is the index for number of particles.
            hei[0] = particle_list[0].node[0].x;//hei is the elevation of each strip
            po_per_str[0] = 1;
            nodes_in_hei[0][0] = 0;


            for (INT32 hf = 1; hf < particle_num_strip - 1; hf++)
            {
                INT32 ghk = 0;
                tetak = (double)(hf)*M_PI / (double)(particle_num_strip - 1);

                ////////// Method 1
                po_per_str[hf] = INT32(particle_cir_factor * particle_num_strip * sin(tetak));
                //////// End of Method 1

                ////////// Method 2
                //double dist_to_cent = particle_radius_2*sqrt(2)*sin(tetak);
                //po_per_str[hf] = abs(INT32(2 * M_PI*dist_to_cent*cir_factor));
                //////// End of Method 2

                if (po_per_str[hf] < 2)
                {
                    po_per_str[hf] = 2;
                }
                fik = 0.0;
                for (INT32 hd = 1; hd <= po_per_str[hf]; hd++)
                {
                    fik = fik + 2. * M_PI / po_per_str[hf];
                    particle_list[0].node[kkk].y = particle_list[0].center[0].y + particle_radius_2 * cos(fik) * sin(tetak);
                    particle_list[0].node[kkk].z = particle_list[0].center[0].z + particle_radius_2 * sin(fik) * sin(tetak);
                    particle_list[0].node[kkk].x = particle_list[0].center[0].x + particle_radius * cos(tetak);
                    hei[hf] = particle_list[0].node[kkk].x;
                    //nodes_in_hei[hf].push_back(kkk);// means nodes_in_hei[hf][ghk] = kkk;
                    nodes_in_hei[hf][ghk] = kkk;

                    //cout << hf << " " << nodes_in_hei[hf][ghk] << " " << ghk << endl;
                    ghk = ghk + 1;
                    kkk = kkk + 1;
                }
            }
            particle_list[0].node[kkk].x = particle_list[0].center[0].x - particle_radius;//Last point
            particle_list[0].node[kkk].y = particle_list[0].center[0].y;
            particle_list[0].node[kkk].z = particle_list[0].center[0].z;
            po_per_str[particle_num_strip - 1] = 1;
            nodes_in_hei[particle_num_strip - 1][0] = kkk;//Index of point that is located on this strip (last strip)
            hei[particle_num_strip - 1] = particle_list[0].node[kkk].x;//hei is the elevation of each strip
            particle_num_nodes = kkk + 1; //The actual particle_num_nodes for spheroid

            mmm = 0;
            double sum, h1, h2, a1, a2;
            sum = 0.0;//Total area. It is to make a comparison with a Full Spheroid formula
            for (INT32 hf = 0; hf < particle_num_strip - 1; hf++)
            {
                /////////tetak = (double)(hf)* M_PI / (double)(num_strip - 1);
                for (INT32 hd = 1; hd <= po_per_str[hf]; hd++)
                {
                    if (hf != 0)
                    {
                        h1 = (hei[hf] - particle_list[0].center[0].x + hei[hf - 1] - particle_list[0].center[0].x) / 2. * hh;
                    }
                    else
                    {
                        h1 = cc * hh;//(hei[hf]-particle_list[0].center[0].x ) * hh;
                    }
                    h2 = (hei[hf] - particle_list[0].center[0].x + hei[hf + 1] - particle_list[0].center[0].x) / 2. * hh;
                    if (cc > aa)
                    {
                        a1 = M_PI * aa * 2 * (asin(h1) / 2. + (h1 * sqrt(1 - SQ(h1))) / 2.) / hh / po_per_str[hf];
                        a2 = M_PI * aa * 2 * (asin(h2) / 2. + (h2 * sqrt(1 - SQ(h2))) / 2.) / hh / po_per_str[hf];
                    }
                    else
                    {
                        a1 = M_PI * aa * 2 * (asinh(h1) / 2. + (h1 * sqrt(1 + SQ(h1))) / 2.) / hh / po_per_str[hf];
                        a2 = M_PI * aa * 2 * (asinh(h2) / 2. + (h2 * sqrt(1 + SQ(h2))) / 2.) / hh / po_per_str[hf];
                    }
                    particle_list[0].node[mmm].area_node = abs(a1 - a2);
                    sum = sum + particle_list[0].node[mmm].area_node;
                    mmm = mmm + 1;
                }
            }
            particle_list[0].node[mmm].area_node = particle_list[0].node[0].area_node;//The last point
            sum = sum + particle_list[0].node[mmm].area_node;
            //cout<<"total area: "<<sum<<endl;

            /////////////////////////////// Below part is just for creating surface triangles for Particle_VTK
            for (INT32 i = 0; i < po_per_str[1]; i++)
            {
                triangleset1[i] = 0;
                triangleset2[i] = nodes_in_hei[1][i];
                triangleset3[i] = nodes_in_hei[1][(i + 1) % po_per_str[1]];
            }
            // For Particle_vtk we want to create a surface by triangles.
            kkk = 1;
            INT32 tindex = po_per_str[1];
            double dist_node{};
            double min_dist{};
            INT32 min_node{}, ipp{};
            //double max_dist = -2;
            INT32 ik;

            for (INT32 hf = 1; hf < particle_num_strip - 1; hf++)
            {
                for (INT32 ij = 0; ij < po_per_str[hf]; ij++)
                {
                    ik = hf - 1;
                    min_dist = 2 * std::max(particle_radius, particle_radius_2);
                    for (INT32 ip = 0; ip < po_per_str[ik]; ip++)
                    {
                        dist_node = sqrt(SQ(particle_list[0].node[nodes_in_hei[hf][ij]].x - particle_list[0].node[nodes_in_hei[ik][ip]].x) + SQ(particle_list[0].node[nodes_in_hei[hf][ij]].y - particle_list[0].node[nodes_in_hei[ik][ip]].y) + SQ(particle_list[0].node[nodes_in_hei[hf][ij]].z - particle_list[0].node[nodes_in_hei[ik][ip]].z));
                        if (ip == 0)
                        {
                            ipp = 0;
                        }

                        if (dist_node < min_dist)
                        {
                            min_node = nodes_in_hei[ik][ip];
                            ipp = ip;
                            min_dist = dist_node;
                        }
                    }

                    triangleset2[tindex] = nodes_in_hei[ik][ipp];
                    triangleset1[tindex] = kkk;
                    triangleset3[tindex] = nodes_in_hei[ik][(ipp + 1) % po_per_str[ik]];
                    //cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
                    tindex = tindex + 1;
                    triangleset2[tindex] = nodes_in_hei[ik][ipp];
                    triangleset1[tindex] = kkk;
                    triangleset3[tindex] = nodes_in_hei[ik][(ipp - 1 + po_per_str[ik]) % po_per_str[ik]];
                    //cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
                    tindex = tindex + 1;

                    ik = hf + 1;
                    min_dist = 2 * std::max(particle_radius, particle_radius_2);
                    for (INT32 ip = 0; ip < po_per_str[ik]; ip++)
                    {
                        dist_node = sqrt(SQ(particle_list[0].node[nodes_in_hei[hf][ij]].x - particle_list[0].node[nodes_in_hei[ik][ip]].x) + SQ(particle_list[0].node[nodes_in_hei[hf][ij]].y - particle_list[0].node[nodes_in_hei[ik][ip]].y) + SQ(particle_list[0].node[nodes_in_hei[hf][ij]].z - particle_list[0].node[nodes_in_hei[ik][ip]].z));
                        if (ip == 0)
                        {
                            ipp = 0;
                        }

                        if (dist_node < min_dist)
                        {
                            min_node = nodes_in_hei[ik][ip];
                            ipp = ip;
                            min_dist = dist_node;
                        }
                    }

                    triangleset2[tindex] = nodes_in_hei[ik][ipp];
                    triangleset1[tindex] = kkk;
                    triangleset3[tindex] = nodes_in_hei[ik][(ipp + 1) % po_per_str[ik]];
                    //cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
                    tindex = tindex + 1;
                    triangleset2[tindex] = nodes_in_hei[ik][ipp];
                    triangleset1[tindex] = kkk;
                    triangleset3[tindex] = nodes_in_hei[ik][(ipp - 1 + po_per_str[ik]) % po_per_str[ik]];
                    //cout << tindex << " " << triangleset1[tindex] << " " << triangleset2[tindex] << " " << triangleset3[tindex] << endl;
                    tindex = tindex + 1;
                    kkk = kkk + 1;
                }
            }
            number_triangles = tindex;

            for (INT32 i = 0; i < po_per_str[particle_num_strip - 2]; i++)
            {
                triangleset1[tindex + i] = particle_num_nodes - 1;
                triangleset2[tindex + i] = nodes_in_hei[particle_num_strip - 2][i];
                triangleset3[tindex + i] = nodes_in_hei[particle_num_strip - 2][(i + 1) % po_per_str[1]];
                number_triangles++;
            }
            triangleset1.resize(number_triangles);///AAAA new
            triangleset2.resize(number_triangles);
            triangleset3.resize(number_triangles);
            //////////////////////////////////End of Particle_VTK part
        } //spheroid
    }
    /// ***********************************
    /// ***********************************
    // Resize node vector and initialize each node in the particle
    for (INT32 i = 0; i < num_particles; ++i)
    {
        particle_list[i].node.resize(particle_num_nodes);// Mainly to remove the extra unnecessary nodes of Spheroids
        for (INT32 j = 0; j < particle_num_nodes; ++j)
        {
            // Initialize node temperature
            particle_list[i].node[j].temperature = particle_temperature; // Initialize temperature at nodes. Later will come from "Interpolate"
        }
    }

    /// ***********************************
    /// LOCATION OF 2nd UNTIL Nth PARTICLE
    /// ***********************************
    if (Random_Arrangement_ON)
    {
        for (INT32 i = 1; i < num_particles; ++i)
        {
            bool part_position{};
            srand(5323);
            double max_dist{};
            if (MOVING_SPHERE || STATIONARY_SPHERE || STATIONARY_CYLINDER || MOVING_CYLINDER)
            {
                max_dist = particle_list[0].radius;
            }
            if (MOVING_SPHEROID || STATIONARY_SPHEROID || STATIONARY_ELLIPSE || MOVING_ELLIPSE)
            {
                max_dist = std::max(particle_list[0].radius_2, particle_list[0].radius);
            }
            do
            {
                part_position = true;
                particle_list[i].center[0].x = rand() % Nx;
                particle_list[i].center[0].y = rand() % Ny;
                particle_list[i].center[0].z = rand() % Nz;
                if (STATIONARY_CYLINDER || MOVING_CYLINDER || STATIONARY_ELLIPSE || MOVING_ELLIPSE) // Only for 2D particles
                {
                    particle_list[i].center[0].z = 0.0;
                }
                printf("Particle i is %d and x, y, z are %5.3f , %5.3f, %5.3f \n", i + 1, particle_list[i].center[0].x, particle_list[i].center[0].y, particle_list[i].center[0].z);

                for (INT32 j = 0; j < num_particles; ++j)
                {
                    if (j != i)
                    {
                        double Dij = sqrt(SQ(particle_list[i].center[0].x - particle_list[j].center[0].x) + SQ(particle_list[i].center[0].y - particle_list[j].center[0].y) + SQ(particle_list[i].center[0].z - particle_list[j].center[0].z));
                        //#if defined MOVING_SPHERE || defined STATIONARY_SPHERE
                        if (Dij < 2 * max_dist + 1)
                        {
                            part_position = false;
                            printf("Error, contact with particle %d \n", j);
                            //	break;
                        }
                        //#endif
                    }
                }

                if (particle_list[i].center[0].x > Nx - 1 - max_dist - 1 || particle_list[i].center[0].x < max_dist + 1)
                {
                    part_position = false;
                    cout << "Error 1" << endl;
                }
                if (particle_list[i].center[0].y > Ny - 1 - max_dist - 1 || particle_list[i].center[0].y < max_dist + 1)
                {
                    part_position = false;
                    cout << "Error 2" << endl;
                }
                if (particle_list[i].center[0].z > Nz - 1 - max_dist - 1 || particle_list[i].center[0].z < max_dist + 1)
                {
                    part_position = false;
                    cout << "Error 3" << endl;
                }
            }
            while (part_position == false);

            for (INT32 n = 0; n < particle_num_nodes; ++n)
            {
                particle_list[i].node[n].area_node = particle_list[0].node[n].area_node;
                particle_list[i].node[n].x = particle_list[0].node[n].x - particle_list[0].center[0].x + particle_list[i].center[0].x;
                particle_list[i].node[n].y = particle_list[0].node[n].y - particle_list[0].center[0].y + particle_list[i].center[0].y;
                particle_list[i].node[n].z = particle_list[0].node[n].z - particle_list[0].center[0].z + particle_list[i].center[0].z;
            }
        }
    }// End of Random_Arrangement_ON
    if (Random_Arrangement_OFF)
    {
        for (INT32 i = 1; i < num_particles; ++i)
        {
            particle_list[i].center[0].x = particle_center_x + i * particle_offset_x;//Location of other particles is adjusted according to location of first one. You can change it.
            particle_list[i].center[0].y = particle_center_y + i * particle_offset_y;
            particle_list[i].center[0].z = particle_center_z + i * particle_offset_z;
            for (INT32 n = 0; n < particle_num_nodes; ++n)
            {
                particle_list[i].node[n].area_node = particle_list[0].node[n].area_node;
                particle_list[i].node[n].x = particle_list[0].node[n].x - particle_list[0].center[0].x + particle_list[i].center[0].x;
                particle_list[i].node[n].y = particle_list[0].node[n].y - particle_list[0].center[0].y + particle_list[i].center[0].y;
                particle_list[i].node[n].z = particle_list[0].node[n].z - particle_list[0].center[0].z + particle_list[i].center[0].z;
                if (particle_list[i].node[n].x > Nx - 1 - 0.1 || particle_list[i].node[n].x < 0.1)
                {
                    printf("ERROR: Particle %d position is out of the domain \n", i + 1);
                    break;
                }
                if (particle_list[i].node[n].y > Ny - 1 - 0.1 || particle_list[i].node[n].y < 0.1)
                {
                    printf("ERROR: Particle %d position is out of the domain \n", i + 1);
                    break;
                }
                if (particle_list[i].node[n].z > Nz - 1 - 0.1 || particle_list[i].node[n].z < 0.1)
                {
                    printf("ERROR: Particle %d position is out of the domain \n", i + 1);
                    break;
                }
            }
        }
    }// End of Random_Arrangement_OFF
    /// ***********************************
    /// ***********************************

    for (INT32 np = 0; np < num_particles; ++np)
    {
        update_verlet();
        if (MOVING_SPHEROID	|| STATIONARY_SPHEROID)
        {
            //// In the next 9 lines we update each Lagrangian node position if the initial values of Euler angles are not zero. (It means particle rotation)
            quaternion(np);
            particle_list[np].center[0].deter = determinant(np);//determinant of Quaternion matrix
            for (INT32 n = 0; n < particle_num_nodes; ++n)
            {
                particle_list[np].node[n].x_pos_1 = particle_list[np].node[n].x - particle_list[np].center[0].x;
                particle_list[np].node[n].y_pos_1 = particle_list[np].node[n].y - particle_list[np].center[0].y;
                particle_list[np].node[n].z_pos_1 = particle_list[np].node[n].z - particle_list[np].center[0].z;
            }
            update_rot_spheroid(np);// Node location at INERTIAL coordinate is calculated based on the Quaternions
        }
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SET_PARTICLE_FORCES_TORQUES_ZERO
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::set_particle_forces_torques_zero(INT32 np){//set all forces and torques of the particle = 0

    for (INT32 nn = 0; nn < particle_num_nodes; ++nn) {
        particle_list[np].node[nn].force_x = 0;
        particle_list[np].node[nn].force_y = 0;
        particle_list[np].node[nn].force_z = 0;
        particle_list[np].node[nn].force_thermal = 0.0;
        particle_list[np].node[nn].total_node_mdf_force_x = 0;
        particle_list[np].node[nn].total_node_mdf_force_y = 0;
        particle_list[np].node[nn].total_node_mdf_force_z = 0;
        particle_list[np].node[nn].total_node_mdf_force_thermal = 0.0;
    }

    particle_list[np].center[0].heat_surf = 0.0;

    particle_list[np].center[0].fx_surf = 0.0;
    particle_list[np].center[0].fy_surf = 0.0;
    particle_list[np].center[0].fz_surf = 0.0;
    particle_list[np].center[0].f3x = 0.0; //Total Force due to particle-particle collision
    particle_list[np].center[0].f3y = 0.0; //Total Force due to particle-particle collision
    particle_list[np].center[0].f3z = 0.0; //Total Force due to particle-particle collision
    particle_list[np].center[0].f4x = 0.0; //Total Force due to particle-wall collision
    particle_list[np].center[0].f4y = 0.0; //Total Force due to particle-wall collision
    particle_list[np].center[0].f4z = 0.0; //Total Force due to particle-wall collision

    particle_list[np].center[0].t1_x = 0.0;//Torque due to forces in X direction
    particle_list[np].center[0].t1_y = 0.0;//Torque due to forces in Y direction
    particle_list[np].center[0].t1_z = 0.0;//Torque due to forces in Z direction

    particle_list[np].center[0].t2_x = 0.0;//Torque due to forces in X direction
    particle_list[np].center[0].t2_y = 0.0;//Torque due to forces in Y direction
    particle_list[np].center[0].t2_z = 0.0;//Torque due to forces in Z direction

    particle_list[np].center[0].t3_x = 0.0;//Torque due to collision in X direction
    particle_list[np].center[0].t3_y = 0.0;//Torque due to collision in Y direction
    particle_list[np].center[0].t3_z = 0.0;//Torque due to collision in Z direction

    particle_list[np].center[0].T_prime_X = 0.0;
    particle_list[np].center[0].T_prime_Y = 0.0;
    particle_list[np].center[0].T_prime_Z = 0.0;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: DIRECT_FORCING
//-----------------------------------------------------------------------------------------------------------------------------
INT32 ParticleResolved::calculate_gbetaindex(INT32 X, INT32 Y, INT32 Z, const ParallelMPI& paralLB) const
{
    ///Checks if the point (X, Y, Z) is inside the ellipse/circle or not
    INT32 X_global = paralLB.get_Xglobal_from_Xlocal(X, this_processor_id);// X + start_global_x[this_processor_id] - num_buffer_layers_X[this_processor_id];
    INT32 Y_global = paralLB.get_Yglobal_from_Ylocal(Y, this_processor_id);//Y + start_global_y[this_processor_id] - num_buffer_layers_Y[this_processor_id];
    INT32 Z_global = paralLB.get_Zglobal_from_Zlocal(Z, this_processor_id);//Z + start_global_z[this_processor_id] - num_buffer_layers_Z[this_processor_id];

    INT32 gbindex = 1;//Will be 0 outside particle and 1 inside particle
    for (INT32 npa = 0; npa < num_particles; npa++)
    {
        if (STATIONARY_CYLINDER || MOVING_CYLINDER){
            if (sqrt(SQ(X_global - particle_list[npa].center[0].x) + SQ(Y_global- particle_list[npa].center[0].y)) <= particle_radius) {
                gbindex = 0;
            }
        }
        if (STATIONARY_ELLIPSE || MOVING_ELLIPSE){
            double alphaAngle = particle_list[npa].center[0].tetaZ - M_PI / 2.;
            if ((SQ((X_global - particle_list[npa].center[0].x) *  cos(alphaAngle) + (Y_global - particle_list[npa].center[0].y) * sin(alphaAngle)) / SQ(particle_radius) + SQ((X_global - particle_list[npa].center[0].x) * sin(alphaAngle) - (Y_global - particle_list[npa].center[0].y) * cos(alphaAngle)) / SQ(particle_radius_2)) <= 1) {
                gbindex = 0;
            }
            //cout << alphaAngle << endl;
            //if ((X + start_x[processor_id] - 2 == 110) && (Y == 1023)) { { cout << gbetaindex<<endl;  } }
        }
    }
    return gbindex;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: DIRECT_FORCING
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::direct_forcing(INT32 time, FlowSolver& flowSolverLB, const ParallelMPI& paralLB)//Can modify FlowSolver properties (like Eulerian forcees)
{

    if ((all_particles_disp_max_1 + all_particles_disp_max_2) > (particle_r_buffer_verlet - particle_r_cutoff_verlet))
    {
        update_verlet();
        all_particles_disp_max_1 = 0;
        all_particles_disp_max_2 = 0;
        for (INT32 ii = 0; ii < num_particles; ++ii)
        {
            particle_list[ii].center[0].displacement_verlet = 0;
        }
    }

    if (SDF)
    {
        for (INT32 np = 0; np < num_particles; ++np)
        {
            set_particle_forces_torques_zero(np);

            interpolate(np, flowSolverLB, paralLB ,0); // interpolate velocity from the Eulerian mesh to Lagrangian node

            collision_model(np);

            DF_force_calc(np);

            spread(np, flowSolverLB, paralLB); // spread forces from the Lagrangian to the Eulerian mesh

            total_force_torque_calc(np);//Total force and heat calculation exerted on the "center"
        }

    } //End of:SDF

    if (MDF)
    {
        for (INT32 np = 0; np < num_particles; ++np)
        {
            set_particle_forces_torques_zero(np);

            interpolate(np, flowSolverLB, paralLB, 0); // interpolate velocity from the Eulerian mesh to Lagrangian node

            collision_model(np);

            for (INT32 df_iter = 1; df_iter <= mdfstep; ++df_iter){//mdf iterations. mdfstep = 1 means SDF

                DF_force_calc(np);

                #if STORE_1D_VECTOR
                std::fill(flowSolverLB.mdfVelocity.begin(), flowSolverLB.mdfVelocity.end(), 0.0);
                #endif

                spread(np, flowSolverLB, paralLB); // spread forces from the Lagrangian to the Eulerian mesh

                INT32 start_x = MAX(0,      particle_list[np].center[0].x - MAX(particle_radius, particle_radius_2) - 4) + num_buffer_layers_X[this_processor_id];//instead of start_X
                INT32 end_x   = MIN(Nx - 1, particle_list[np].center[0].x + MAX(particle_radius, particle_radius_2) + 4) + num_buffer_layers_X[this_processor_id];//instead of end_X
                INT32 start_y = MAX(0,      particle_list[np].center[0].y - MAX(particle_radius, particle_radius_2) - 4) + num_buffer_layers_Y[this_processor_id];//instead of start_Y
                INT32 end_y   = MIN(Ny - 1, particle_list[np].center[0].y + MAX(particle_radius, particle_radius_2) + 4) + num_buffer_layers_Y[this_processor_id];//instead of end_Y
                INT32 start_z = MAX(0,      particle_list[np].center[0].z - MAX(particle_radius, particle_radius_2) - 4) + num_buffer_layers_Z[this_processor_id];//instead of start_Z
                INT32 end_z   = MIN(Nz - 1, particle_list[np].center[0].z + MAX(particle_radius, particle_radius_2) + 4) + num_buffer_layers_Z[this_processor_id];//instead of end_Z

                for (INT32 X = start_x; X <= end_x; ++X) {
                    for (INT32 Y = start_y; Y <= end_y; ++Y) {
                        for (INT32 Z = start_z; Z <= end_z; ++Z) {
                            #if !STORE_1D_VECTOR
                            flowSolverLB.mdfVelocity[0][X][Y][Z] = flowSolverLB.velocity[0][X][Y][Z] + flowSolverLB.force[0][X][Y][Z] * Dt_LB;// (See Luo et al. PRE 066709 (2007) http://dx.doi.org/10.1103/PhysRevE.76.066709)
                            flowSolverLB.mdfVelocity[1][X][Y][Z] = flowSolverLB.velocity[1][X][Y][Z] + flowSolverLB.force[1][X][Y][Z] * Dt_LB;
                            flowSolverLB.mdfVelocity[2][X][Y][Z] = flowSolverLB.velocity[2][X][Y][Z] + flowSolverLB.force[2][X][Y][Z] * Dt_LB;
                            if (Flow_With_Thermal_Effect)
                            {
                                flowSolverLB.mdfTemperature[X][Y][Z] = flowSolverLB.temperature[X][Y][Z]+ flowSolverLB.force_thermal[X][Y][Z] * Dt_LB;
                            }
                            #endif
                            #if STORE_1D_VECTOR
                            INT32 id3 = idx_3d_vector(X, Y, Z);
                            INT32 id4_0 = idx_4d_vector(0, X, Y, Z);
                            INT32 id4_1 = idx_4d_vector(1, X, Y, Z);
                            INT32 id4_2 = idx_4d_vector(2, X, Y, Z);
                            flowSolverLB.mdfVelocity[id4_0] = flowSolverLB.velocity[id4_0] + flowSolverLB.force[id4_0] * Dt_LB;// (See Luo et al. PRE 066709 (2007) http://dx.doi.org/10.1103/PhysRevE.76.066709)
                            flowSolverLB.mdfVelocity[id4_1] = flowSolverLB.velocity[id4_1] + flowSolverLB.force[id4_1] * Dt_LB;
                            flowSolverLB.mdfVelocity[id4_2] = flowSolverLB.velocity[id4_2] + flowSolverLB.force[id4_2] * Dt_LB;
                            if (Flow_With_Thermal_Effect)
                            {
                                flowSolverLB.mdfTemperature[id3]= flowSolverLB.temperature[id3]+ flowSolverLB.force_thermal[id3] * Dt_LB;
                            }
                            #endif
                        }
                    }
                }
                interpolate(np, flowSolverLB, paralLB, df_iter); // interpolate velocity from the Eulerian mesh to Lagrangian node
            }
            total_force_torque_calc(np);//Total force and heat calculation exerted on the "center"
        }
    } //End of:MDF
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: DF_force_calc
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::DF_force_calc(INT32 npn){
    for (INT32 n = 0; n < particle_num_nodes; ++n)// Compute IBM Force
    {
        particle_list[npn].node[n].force_x       = (1.0) * (particle_list[npn].center[0].previous_vel_x - particle_list[npn].center[0].omgZ * (particle_list[npn].node[n].y - particle_list[npn].center[0].y) + particle_list[npn].center[0].omgY * (particle_list[npn].node[n].z - particle_list[npn].center[0].z) - particle_list[npn].node[n].vel_x) * particle_list[npn].node[n].area_node;//Surface force
        particle_list[npn].node[n].force_y       = (1.0) * (particle_list[npn].center[0].previous_vel_y + particle_list[npn].center[0].omgZ * (particle_list[npn].node[n].x - particle_list[npn].center[0].x) - particle_list[npn].center[0].omgX * (particle_list[npn].node[n].z - particle_list[npn].center[0].z) - particle_list[npn].node[n].vel_y) * particle_list[npn].node[n].area_node;//Surface force
        particle_list[npn].node[n].force_z       = (1.0) * (particle_list[npn].center[0].previous_vel_z + particle_list[npn].center[0].omgX * (particle_list[npn].node[n].y - particle_list[npn].center[0].y) - particle_list[npn].center[0].omgY * (particle_list[npn].node[n].x - particle_list[npn].center[0].x) - particle_list[npn].node[n].vel_z) * particle_list[npn].node[n].area_node;//Surface force
        particle_list[npn].node[n].force_thermal = (1.0) * (particle_list[npn].center[0].temp_part - particle_list[npn].node[n].temperature) * particle_list[npn].node[n].area_node;//Surface heat force

        ///total force on each node after each mdf step
		particle_list[npn].node[n].total_node_mdf_force_x       += particle_list[npn].node[n].force_x;
        particle_list[npn].node[n].total_node_mdf_force_y       += particle_list[npn].node[n].force_y;
        particle_list[npn].node[n].total_node_mdf_force_z       += particle_list[npn].node[n].force_z;
		particle_list[npn].node[n].total_node_mdf_force_thermal += particle_list[npn].node[n].force_thermal;
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: TOTAL_FORCE_TORQUE_CALC
//-----------------------------------------------------------------------------------------------------------------------------
// Calculates total force on the particle center (at the end of all mdfs)
void ParticleResolved::total_force_torque_calc(INT32 npn){

	for (INT32 n = 0; n < particle_num_nodes; ++n) {

        particle_list[npn].center[0].fx_surf   += particle_list[npn].node[n].total_node_mdf_force_x;
		particle_list[npn].center[0].fy_surf   += particle_list[npn].node[n].total_node_mdf_force_y;
		particle_list[npn].center[0].fz_surf   += particle_list[npn].node[n].total_node_mdf_force_z;
		particle_list[npn].center[0].heat_surf += particle_list[npn].node[n].total_node_mdf_force_thermal;

        particle_list[npn].center[0].t1_x += particle_list[npn].node[n].total_node_mdf_force_z * (particle_list[npn].node[n].y - particle_list[npn].center[0].y);
        particle_list[npn].center[0].t2_x += particle_list[npn].node[n].total_node_mdf_force_y * (particle_list[npn].node[n].z - particle_list[npn].center[0].z);
        particle_list[npn].center[0].t1_y += particle_list[npn].node[n].total_node_mdf_force_z * (particle_list[npn].node[n].x - particle_list[npn].center[0].x);
        particle_list[npn].center[0].t2_y += particle_list[npn].node[n].total_node_mdf_force_x * (particle_list[npn].node[n].z - particle_list[npn].center[0].z);
        particle_list[npn].center[0].t1_z += particle_list[npn].node[n].total_node_mdf_force_x * (particle_list[npn].node[n].y - particle_list[npn].center[0].y);
        particle_list[npn].center[0].t2_z += particle_list[npn].node[n].total_node_mdf_force_y * (particle_list[npn].node[n].x - particle_list[npn].center[0].x);

	}

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: newton_eqs_motion
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::newton_eqs_motion(INT32 time, FlowSolver& flowSolverLB)
{

    for (INT32 np = 0; np < num_particles; ++np)
    {
        check_zero(particle_list[np].center[0].t1_x, 1.0e-6);
        check_zero(particle_list[np].center[0].t1_y, 1.0e-6);
        check_zero(particle_list[np].center[0].t1_z, 1.0e-6);

        check_zero(particle_list[np].center[0].t2_x, 1.0e-6);
        check_zero(particle_list[np].center[0].t2_y, 1.0e-6);
        check_zero(particle_list[np].center[0].t2_z, 1.0e-6);

        /// Update Translational Velocity and center position
        if (time > margin)
        {
            if (time > margin + 1)
            {
                if (With_added_mass)
                {
                    particle_list[np].center[0].vel_x = particle_list[np].center[0].vel_x - particle_list[np].center[0].fx_surf / particle_list[np].particle_mass + particle_gravity_x * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3x / particle_list[np].particle_mass + particle_list[np].center[0].f4x / particle_list[np].particle_mass + 1.0 / particle_den_ratio * (particle_list[np].center[0].vel_x - particle_list[np].center[0].tempx[0]);//Newton law
                    particle_list[np].center[0].vel_y = particle_list[np].center[0].vel_y - particle_list[np].center[0].fy_surf / particle_list[np].particle_mass + particle_gravity_y * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3y / particle_list[np].particle_mass + particle_list[np].center[0].f4y / particle_list[np].particle_mass + 1.0 / particle_den_ratio * (particle_list[np].center[0].vel_y - particle_list[np].center[0].tempy[0]);//Newton law
                    particle_list[np].center[0].vel_z = particle_list[np].center[0].vel_z - particle_list[np].center[0].fz_surf / particle_list[np].particle_mass + particle_gravity_z * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3z / particle_list[np].particle_mass + particle_list[np].center[0].f4z / particle_list[np].particle_mass + 1.0 / particle_den_ratio * (particle_list[np].center[0].vel_z - particle_list[np].center[0].tempz[0]);//Newton law
                }
                if (Without_added_mass)
                {
                    particle_list[np].center[0].vel_x = particle_list[np].center[0].vel_x - particle_list[np].center[0].fx_surf / particle_list[np].particle_mass + particle_gravity_x * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3x / particle_list[np].particle_mass + particle_list[np].center[0].f4x / particle_list[np].particle_mass;//Newton law
                    particle_list[np].center[0].vel_y = particle_list[np].center[0].vel_y - particle_list[np].center[0].fy_surf / particle_list[np].particle_mass + particle_gravity_y * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3y / particle_list[np].particle_mass + particle_list[np].center[0].f4y / particle_list[np].particle_mass;//Newton law
                    particle_list[np].center[0].vel_z = particle_list[np].center[0].vel_z - particle_list[np].center[0].fz_surf / particle_list[np].particle_mass + particle_gravity_z * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3z / particle_list[np].particle_mass + particle_list[np].center[0].f4z / particle_list[np].particle_mass;//Newton law
                }
            }
            else
            {
                particle_list[np].center[0].vel_x = particle_list[np].center[0].vel_x - particle_list[np].center[0].fx_surf / particle_list[np].particle_mass + particle_gravity_x * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3x / particle_list[np].particle_mass + particle_list[np].center[0].f4x / particle_list[np].particle_mass;//Newton law
                particle_list[np].center[0].vel_y = particle_list[np].center[0].vel_y - particle_list[np].center[0].fy_surf / particle_list[np].particle_mass + particle_gravity_y * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3y / particle_list[np].particle_mass + particle_list[np].center[0].f4y / particle_list[np].particle_mass;//Newton law
                particle_list[np].center[0].vel_z = particle_list[np].center[0].vel_z - particle_list[np].center[0].fz_surf / particle_list[np].particle_mass + particle_gravity_z * (1.0 - 1.0 / particle_den_ratio) + particle_list[np].center[0].f3z / particle_list[np].particle_mass + particle_list[np].center[0].f4z / particle_list[np].particle_mass;//Newton law
            }

            particle_list[np].center[0].tempx[0] = particle_list[np].center[0].tempx[1];
            particle_list[np].center[0].tempx[1] = particle_list[np].center[0].vel_x;
            particle_list[np].center[0].tempy[0] = particle_list[np].center[0].tempy[1];
            particle_list[np].center[0].tempy[1] = particle_list[np].center[0].vel_y;
            particle_list[np].center[0].tempz[0] = particle_list[np].center[0].tempz[1];
            particle_list[np].center[0].tempz[1] = particle_list[np].center[0].vel_z;

            particle_list[np].center[0].center_old_x = particle_list[np].center[0].x;
            particle_list[np].center[0].center_old_y = particle_list[np].center[0].y;
            particle_list[np].center[0].center_old_z = particle_list[np].center[0].z;

            double disp_x{}, disp_y{}, disp_z{};
            disp_x = 0.5 * (particle_list[np].center[0].vel_x + particle_list[np].center[0].tempx[0]);
            disp_y = 0.5 * (particle_list[np].center[0].vel_y + particle_list[np].center[0].tempy[0]);
            disp_z = 0.5 * (particle_list[np].center[0].vel_z + particle_list[np].center[0].tempz[0]);
            particle_list[np].center[0].x = particle_list[np].center[0].x + disp_x;
            particle_list[np].center[0].y = particle_list[np].center[0].y + disp_y;
            particle_list[np].center[0].z = particle_list[np].center[0].z + disp_z;

            particle_list[np].center[0].displacement_verlet += sqrt(SQ(disp_x) + SQ(disp_y) + SQ(disp_z));//Displacement between two Verlet list update steps
            if (particle_list[np].center[0].displacement_verlet > all_particles_disp_max_1)
            {
                all_particles_disp_max_2 = all_particles_disp_max_1;
                all_particles_disp_max_1 = particle_list[np].center[0].displacement_verlet;
            }
            else if ((particle_list[np].center[0].displacement_verlet > all_particles_disp_max_2) && (particle_list[np].center[0].displacement_verlet < all_particles_disp_max_1))
            {
                all_particles_disp_max_2 = particle_list[np].center[0].displacement_verlet;
            }
        }

		//keep_particle_at_center(np, flowSolverLB);///To keep particle at center (only along the X direction). Does not work with Recovery_ON at the moment. Change write_particle_data() if needed.

        ////////// Updating Angular Velocity ///////
        ///Rotation of ellipse and cylinder
        if (MOVING_ELLIPSE || MOVING_CYLINDER)   //Calculating rotational speed for sphere
        {
            if (time > margin){
                if (time > margin + 1){
                    if (With_added_mass)
                    {
                        particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / (particle_list[np].mom_inertia_z)*(particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z - particle_list[np].center[0].t3_z) + 1.0 / particle_den_ratio*(particle_list[np].center[0].omgZ - particle_list[np].center[0].temp_omgZ[0]);//Eq. 100 Kang Thesis (804)
                    }
                    if (Without_added_mass)
                    {
                        particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / (particle_list[np].mom_inertia_z)*(particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z - particle_list[np].center[0].t3_z);
                    }
                }
                else{
                    particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / (particle_list[np].mom_inertia_z)*(particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z - particle_list[np].center[0].t3_z);
                }
            }
            if (time == 1){
                particle_list[np].center[0].temp_omgZ[0] = 0.0;
                particle_list[np].center[0].temp_omgZ[1] = 0.0;
            }
            particle_list[np].center[0].temp_omgZ[0] = particle_list[np].center[0].temp_omgZ[1];
            particle_list[np].center[0].temp_omgZ[1] = particle_list[np].center[0].omgZ;

            particle_list[np].center[0].tetaZ = particle_list[np].center[0].tetaZ + 0.5 * (particle_list[np].center[0].omgZ + particle_list[np].center[0].temp_omgZ[0]);
            particle_list[np].center[0].delta_teta_Z = 0.5 * (particle_list[np].center[0].omgZ + particle_list[np].center[0].temp_omgZ[0]);

            // Applying the rotation matrix
            for (INT32 n = 0; n < particle_num_nodes; ++n)
            {
                //Around Z-Axis rotation
                double x_pos = particle_list[np].node[n].x - particle_list[np].center[0].center_old_x;
                double y_pos = particle_list[np].node[n].y - particle_list[np].center[0].center_old_y;
                double z_pos = particle_list[np].node[n].z - particle_list[np].center[0].center_old_z;
                particle_list[np].node[n].x = particle_list[np].center[0].x + x_pos * cos(particle_list[np].center[0].delta_teta_Z) - y_pos * sin(particle_list[np].center[0].delta_teta_Z);
                particle_list[np].node[n].y = particle_list[np].center[0].y + x_pos * sin(particle_list[np].center[0].delta_teta_Z) + y_pos * cos(particle_list[np].center[0].delta_teta_Z);
                particle_list[np].node[n].z = 0;//Only for 2D particles. particle_list[np].center[0].z + z_pos;
            }
        } // End of ellipse and circle rotation

        ///Rotation of sphere
        ////////// Updating Angular Velocity ///////
        if (MOVING_SPHERE)   //Calculating rotational speed for sphere
        {
            /////////// For X rotation
            if (time > margin)
            {
                if (time > margin + 1)
                {
                    if (With_added_mass)
                    {
                        particle_list[np].center[0].omgX = particle_list[np].center[0].omgX + 1.0 / particle_list[np].mom_inertia_x * (particle_list[np].center[0].t2_x - particle_list[np].center[0].t1_x) + 1.0 / particle_den_ratio * (particle_list[np].center[0].omgX - particle_list[np].center[0].temp_omgX[0]);//Eq. 100 Kang Thesis (804)
                    }
                    if (Without_added_mass)
                    {
                        particle_list[np].center[0].omgX = particle_list[np].center[0].omgX + 1.0 / particle_list[np].mom_inertia_x * (particle_list[np].center[0].t2_x - particle_list[np].center[0].t1_x);//Eq. 100 Kang Thesis (804)
                    }
                }
                else
                {
                    particle_list[np].center[0].omgX = particle_list[np].center[0].omgX + 1.0 / particle_list[np].mom_inertia_x * (particle_list[np].center[0].t2_x - particle_list[np].center[0].t1_x);
                }

                check_zero(particle_list[np].center[0].omgX, 1.0e-12);

                particle_list[np].center[0].temp_omgX[0] = particle_list[np].center[0].temp_omgX[1];
                particle_list[np].center[0].temp_omgX[1] = particle_list[np].center[0].omgX;
                /////////// For Y rotation
                if (time > margin + 1)
                {
                    if (With_added_mass)
                    {
                        particle_list[np].center[0].omgY = particle_list[np].center[0].omgY + 1.0 / particle_list[np].mom_inertia_y * (particle_list[np].center[0].t1_y - particle_list[np].center[0].t2_y) + 1.0 / particle_den_ratio * (particle_list[np].center[0].omgY - particle_list[np].center[0].temp_omgY[0]);//Eq. 100 Kang Thesis (804)
                    }
                    if (Without_added_mass)
                    {
                        particle_list[np].center[0].omgY = particle_list[np].center[0].omgY + 1.0 / particle_list[np].mom_inertia_y * (particle_list[np].center[0].t1_y - particle_list[np].center[0].t2_y);//Eq. 100 Kang Thesis (804)
                    }
                }
                else
                {
                    particle_list[np].center[0].omgY = particle_list[np].center[0].omgY + 1.0 / particle_list[np].mom_inertia_y * (particle_list[np].center[0].t1_y - particle_list[np].center[0].t2_y);
                }

                check_zero(particle_list[np].center[0].omgY, 1.0e-12);
                particle_list[np].center[0].temp_omgY[0] = particle_list[np].center[0].temp_omgY[1];
                particle_list[np].center[0].temp_omgY[1] = particle_list[np].center[0].omgY;
                /////////// For Z rotation
                if (time > margin + 1)
                {
                    if (With_added_mass)
                    {
                        particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / particle_list[np].mom_inertia_z * (particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z) + 1.0 / particle_den_ratio * (particle_list[np].center[0].omgZ - particle_list[np].center[0].temp_omgZ[0]);//Eq. 100 Kang Thesis (804)
                    }
                    if (Without_added_mass)
                    {
                        particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / particle_list[np].mom_inertia_z * (particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z);//Eq. 100 Kang Thesis (804)
                    }
                }
                else
                {
                    particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / particle_list[np].mom_inertia_z * (particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z);
                }

                check_zero(particle_list[np].center[0].omgZ, 1.0e-12);
                particle_list[np].center[0].temp_omgZ[0] = particle_list[np].center[0].temp_omgZ[1];
                particle_list[np].center[0].temp_omgZ[1] = particle_list[np].center[0].omgZ;
                //////////////////

                particle_list[np].center[0].delta_teta_X = 0.5 * (particle_list[np].center[0].omgX + particle_list[np].center[0].temp_omgX[0]);
                particle_list[np].center[0].delta_teta_Y = 0.5 * (particle_list[np].center[0].omgY + particle_list[np].center[0].temp_omgY[0]);
                particle_list[np].center[0].delta_teta_Z = 0.5 * (particle_list[np].center[0].omgZ + particle_list[np].center[0].temp_omgZ[0]);
                particle_list[np].center[0].tetaX = particle_list[np].center[0].tetaX + 0.5 * (particle_list[np].center[0].omgX + particle_list[np].center[0].temp_omgX[0]);
                particle_list[np].center[0].tetaY = particle_list[np].center[0].tetaY + 0.5 * (particle_list[np].center[0].omgY + particle_list[np].center[0].temp_omgY[0]);
                particle_list[np].center[0].tetaZ = particle_list[np].center[0].tetaZ + 0.5 * (particle_list[np].center[0].omgZ + particle_list[np].center[0].temp_omgZ[0]);

                // Applying the rotation matrix
                for (INT32 n = 0; n < particle_num_nodes; ++n)
                {
                    //Around Z-Axis rotation
                    double x_pos = particle_list[np].node[n].x - particle_list[np].center[0].center_old_x;
                    double y_pos = particle_list[np].node[n].y - particle_list[np].center[0].center_old_y;
                    double z_pos = particle_list[np].node[n].z - particle_list[np].center[0].center_old_z;
                    particle_list[np].node[n].x = particle_list[np].center[0].x + x_pos * cos(particle_list[np].center[0].delta_teta_Z) - y_pos * sin(particle_list[np].center[0].delta_teta_Z);
                    particle_list[np].node[n].y = particle_list[np].center[0].y + x_pos * sin(particle_list[np].center[0].delta_teta_Z) + y_pos * cos(particle_list[np].center[0].delta_teta_Z);
                    particle_list[np].node[n].z = particle_list[np].center[0].z + z_pos;
                    //Around X-Axis rotation
                    x_pos = particle_list[np].node[n].x - particle_list[np].center[0].x;
                    y_pos = particle_list[np].node[n].y - particle_list[np].center[0].y;
                    z_pos = particle_list[np].node[n].z - particle_list[np].center[0].z;
                    particle_list[np].node[n].x = particle_list[np].center[0].x + x_pos;
                    particle_list[np].node[n].y = particle_list[np].center[0].y + y_pos * cos(particle_list[np].center[0].delta_teta_X) - z_pos * sin(particle_list[np].center[0].delta_teta_X);
                    particle_list[np].node[n].z = particle_list[np].center[0].z + y_pos * sin(particle_list[np].center[0].delta_teta_X) + z_pos * cos(particle_list[np].center[0].delta_teta_X);

                    //Around Y-Axis rotation
                    x_pos = particle_list[np].node[n].x - particle_list[np].center[0].x;
                    y_pos = particle_list[np].node[n].y - particle_list[np].center[0].y;
                    z_pos = particle_list[np].node[n].z - particle_list[np].center[0].z;
                    particle_list[np].node[n].x = particle_list[np].center[0].x + x_pos * cos(particle_list[np].center[0].delta_teta_Y) + z_pos * sin(particle_list[np].center[0].delta_teta_Y);;
                    particle_list[np].node[n].y = particle_list[np].center[0].y + y_pos;
                    particle_list[np].node[n].z = particle_list[np].center[0].z - x_pos * sin(particle_list[np].center[0].delta_teta_Y) + z_pos * cos(particle_list[np].center[0].delta_teta_Y);

                    /////////////////////////
                    /*particle_list[np].node[n].vel_x = particle_list[np].center[0].vel_x;  // I think these 3 lines can be removed because node velocity is determined by INTERPOLATION
                    particle_list[np].node[n].vel_y = particle_list[np].center[0].vel_y;
                    particle_list[np].node[n].vel_z = particle_list[np].center[0].vel_z;*/
                }
            }

        }  //End of rotation procedure for SPHERE

        ///Rotation of spheroid
        if (MOVING_SPHEROID)    //Calculating rotation speed for spheroid. Using quaternion concept. See: H. Huang et al, J. Fluid Mech. (2012), vol. 692, pp. 369-394.
        {
            if (time > margin)
            {
                quaternion(np);
                double torque_1 = (particle_list[np].center[0].t2_x - particle_list[np].center[0].t1_x - particle_list[np].center[0].t3_x);
                double torque_2 = (particle_list[np].center[0].t1_y - particle_list[np].center[0].t2_y - particle_list[np].center[0].t3_y);
                double torque_3 = (particle_list[np].center[0].t1_z - particle_list[np].center[0].t2_z - particle_list[np].center[0].t3_z);
                // Torque transformation from Inertial coordinate to body-fixed coordinate
                particle_list[np].center[0].T_prime_X = particle_list[np].center[0].quater[1] * torque_1 + particle_list[np].center[0].quater[2] * torque_2 + particle_list[np].center[0].quater[3] * torque_3;
                particle_list[np].center[0].T_prime_Y = particle_list[np].center[0].quater[4] * torque_1 + particle_list[np].center[0].quater[5] * torque_2 + particle_list[np].center[0].quater[6] * torque_3;
                particle_list[np].center[0].T_prime_Z = particle_list[np].center[0].quater[7] * torque_1 + particle_list[np].center[0].quater[8] * torque_2 + particle_list[np].center[0].quater[9] * torque_3;

                rungekutta(np);

                quaternion(np);

                double sumq = sqrt(SQ(particle_list[np].center[0].Qu0) + SQ(particle_list[np].center[0].Qu1) + SQ(particle_list[np].center[0].Qu2) + SQ(particle_list[np].center[0].Qu3));
                particle_list[np].center[0].Qu0 = particle_list[np].center[0].Qu0 / sumq;
                particle_list[np].center[0].Qu1 = particle_list[np].center[0].Qu1 / sumq;
                particle_list[np].center[0].Qu2 = particle_list[np].center[0].Qu2 / sumq;
                particle_list[np].center[0].Qu3 = particle_list[np].center[0].Qu3 / sumq;
                particle_list[np].center[0].deter = determinant(np);//determinant of Quaternion matrix

                comp_rot_speed(np);// Angular velocity transformation from body-fixed coordinate to Inertial coordinate.

                if (time > margin + 1)
                {
                    if (With_added_mass)
                    {
                        particle_list[np].center[0].omgX = particle_list[np].center[0].omgX + 1.0 / particle_den_ratio * (particle_list[np].center[0].omgX - particle_list[np].center[0].temp_omgX[0]);//Eq. 100 Kang Thesis (804)
                        particle_list[np].center[0].omgY = particle_list[np].center[0].omgY + 1.0 / particle_den_ratio * (particle_list[np].center[0].omgY - particle_list[np].center[0].temp_omgY[0]);//Eq. 100 Kang Thesis (804)
                        particle_list[np].center[0].omgZ = particle_list[np].center[0].omgZ + 1.0 / particle_den_ratio * (particle_list[np].center[0].omgZ - particle_list[np].center[0].temp_omgZ[0]);//Eq. 100 Kang Thesis (804)
                        particle_list[np].center[0].Omg_prime_X = particle_list[np].center[0].quater[1] * (particle_list[np].center[0].omgX) + particle_list[np].center[0].quater[2] * (particle_list[np].center[0].omgY) + particle_list[np].center[0].quater[3] * (particle_list[np].center[0].omgZ);
                        particle_list[np].center[0].Omg_prime_Y = particle_list[np].center[0].quater[4] * (particle_list[np].center[0].omgX) + particle_list[np].center[0].quater[5] * (particle_list[np].center[0].omgY) + particle_list[np].center[0].quater[6] * (particle_list[np].center[0].omgZ);
                        particle_list[np].center[0].Omg_prime_Z = particle_list[np].center[0].quater[7] * (particle_list[np].center[0].omgX) + particle_list[np].center[0].quater[8] * (particle_list[np].center[0].omgY) + particle_list[np].center[0].quater[9] * (particle_list[np].center[0].omgZ);
                    }
                }

                particle_list[np].center[0].temp_omgX[0] = particle_list[np].center[0].temp_omgX[1];
                particle_list[np].center[0].temp_omgX[1] = particle_list[np].center[0].omgX;
                particle_list[np].center[0].temp_omgY[0] = particle_list[np].center[0].temp_omgY[1];
                particle_list[np].center[0].temp_omgY[1] = particle_list[np].center[0].omgY;
                particle_list[np].center[0].temp_omgZ[0] = particle_list[np].center[0].temp_omgZ[1];
                particle_list[np].center[0].temp_omgZ[1] = particle_list[np].center[0].omgZ;

                particle_list[np].center[0].delta_teta_X = 0.5 * (particle_list[np].center[0].omgX + particle_list[np].center[0].temp_omgX[0]);// For spheroid they are not important because node location is updated by quaternions
                particle_list[np].center[0].delta_teta_Y = 0.5 * (particle_list[np].center[0].omgY + particle_list[np].center[0].temp_omgY[0]);
                particle_list[np].center[0].delta_teta_Z = 0.5 * (particle_list[np].center[0].omgZ + particle_list[np].center[0].temp_omgZ[0]);
                particle_list[np].center[0].tetaX = particle_list[np].center[0].tetaX + 0.5 * (particle_list[np].center[0].omgX + particle_list[np].center[0].temp_omgX[0]);
                particle_list[np].center[0].tetaY = particle_list[np].center[0].tetaY + 0.5 * (particle_list[np].center[0].omgY + particle_list[np].center[0].temp_omgY[0]);
                particle_list[np].center[0].tetaZ = particle_list[np].center[0].tetaZ + 0.5 * (particle_list[np].center[0].omgZ + particle_list[np].center[0].temp_omgZ[0]);


                // Applying the rotation matrix
                //Updating nodes location
                update_rot_spheroid(np);// Node location at INERTIAL coordinate is calculated based on the Quaternions

                /*for (INT32 n = 0; n < particle_num_nodes; ++n) {
                	particle_list[np].node[n].vel_x = particle_list[np].center[0].vel_x;  // I think these 3 lines can be removed because node velocity is determined by INTERPOLATION
                	particle_list[np].node[n].vel_y = particle_list[np].center[0].vel_y;
                	particle_list[np].node[n].vel_z = particle_list[np].center[0].vel_z;
                }*/
            }
        } //End of rotation for spheroid
        ////////////////////////////////////////////

        /// Updating Temperature
        if (Particle_Fixed_Temperature)
        {
            particle_list[np].center[0].temperature += 0; //No change in temperature
        }
        if (Particle_Varying_Temperature)
        {
            if (time > 1)
            {
                particle_list[np].center[0].temperature = particle_list[np].center[0].temperature - particle_list[np].center[0].heat_surf / (particle_list[np].particle_mass * particle_specific_heat_ratio) + 1.0 / (particle_den_ratio * particle_specific_heat_ratio) * (particle_list[np].center[0].temperature - particle_list[np].center[0].temp_temperat[0]) + particle_heat_source / (particle_list[np].particle_mass * particle_specific_heat_ratio);
            }
            else
            {
                particle_list[np].center[0].temperature = particle_list[np].center[0].temperature - particle_list[np].center[0].heat_surf / (particle_list[np].particle_mass * particle_specific_heat_ratio) + particle_heat_source / (particle_list[np].particle_mass * particle_specific_heat_ratio);
            }
        }
        if (time == 1)
        {
            particle_list[np].center[0].temp_temperat[0] = particle_temperature;
            particle_list[np].center[0].temp_temperat[1] = particle_temperature;
        }
        particle_list[np].center[0].temp_temperat[0] = particle_list[np].center[0].temp_temperat[1];
        particle_list[np].center[0].temp_temperat[1] = particle_list[np].center[0].temperature;
        /////////////////////////////////////////

        if (time > margin)
        {
            if (X_Periodic)   //In periodic condition, if particle exits from one side it should enter from other side
            {
                if (particle_list[np].center[0].x < 0)
                {
                    all_particles_disp_max_1 = 100000;// A large number
                    all_particles_disp_max_2 = 100000;
                    particle_list[np].center[0].x += Nx - 1;
                    for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
                    {
                        particle_list[np].node[np1].x += Nx - 1;
                    }
                }
                else if (particle_list[np].center[0].x > Nx - 1)
                {
                    all_particles_disp_max_1 = 100000;// A large number
                    all_particles_disp_max_2 = 100000;
                    particle_list[np].center[0].x -= Nx - 1;
                    for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
                    {
                        particle_list[np].node[np1].x -= Nx - 1;
                    }
                }
            }

            if (Z_Periodic)  //In periodic condition, if particle exits from one side it should enter from other side
            {
                if (particle_list[np].center[0].z < 0)
                {
                    all_particles_disp_max_1 = 100000;// A large number
                    all_particles_disp_max_2 = 100000;
                    particle_list[np].center[0].z += Nz - 1;
                    for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
                    {
                        particle_list[np].node[np1].z += Nz - 1;
                    }
                }
                else if (particle_list[np].center[0].z > Nz - 1)
                {
                    all_particles_disp_max_1 = 100000;// A large number
                    all_particles_disp_max_2 = 100000;
                    particle_list[np].center[0].z -= Nz - 1;
                    for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
                    {
                        particle_list[np].node[np1].z -= Nz - 1;
                    }
                }
            }
            if (Y_Periodic)  //In periodic condition, if particle exits from one side it should enter from other side
            {
                if (particle_list[np].center[0].y < 0)
                {
                    all_particles_disp_max_1 = 100000;// A large number
                    all_particles_disp_max_2 = 100000;
                    particle_list[np].center[0].y += Ny - 1;
                    for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
                    {
                        particle_list[np].node[np1].y += Ny - 1;
                    }
                }
                else if (particle_list[np].center[0].y > Ny - 1)
                {
                    all_particles_disp_max_1 = 100000;// A large number
                    all_particles_disp_max_2 = 100000;
                    particle_list[np].center[0].y -= Ny - 1;
                    for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
                    {
                        particle_list[np].node[np1].y -= Ny - 1;
                    }
                }
            }
        } //End of : if (t > margin)

        particle_list[np].center[0].previous_vel_x = particle_list[np].center[0].vel_x;
        particle_list[np].center[0].previous_vel_y = particle_list[np].center[0].vel_y;
        particle_list[np].center[0].previous_vel_z = particle_list[np].center[0].vel_z;

        if (Particle_Varying_Temperature)
        {
            particle_list[np].center[0].temp_part = particle_list[np].center[0].temperature;
        }

    } //End of: NEWTON'S EQUATIONS OF MOTION
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: INTERPOLATE
//-----------------------------------------------------------------------------------------------------------------------------
//
// The node velocities are interpolated from the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) and other interpolation schemes are included in the present code.

void ParticleResolved::interpolate(INT32 np, const FlowSolver& flowSolverLB, const ParallelMPI& paralLB, INT32 mdf_iteration)
{
    INT32 x_int{};
    INT32 y_int{};
    INT32 z_int{};
    INT32 xmin{};
    INT32 xmax{};
    INT32 ymin{};
    INT32 ymax{};
    INT32 zmin{};
    INT32 zmax{};
    double dist_x{};
    double dist_y{};
    double dist_z{};
    double weight_x{};
    double weight_y{};
    double weight_z{};

    // Run over all object nodes.

    for (INT32 n = 0; n < particle_num_nodes; ++n)
    {

        // Reset node velocity first since '+=' is used.
        particle_list[np].node[n].vel_x = 0;
        particle_list[np].node[n].vel_y = 0;
        particle_list[np].node[n].vel_z = 0;
        particle_list[np].node[n].den = 0;
        particle_list[np].node[n].temperature = 0;

        // Identify the lowest fluid lattice node in interpolation range (see spreading).

            x_int = floor(particle_list[np].node[n].x);
            y_int = floor(particle_list[np].node[n].y);
            z_int = floor(particle_list[np].node[n].z);
            // Run over all neighboring fluid nodes.
            // In the case of the two-point interpolation, it is 2x2 fluid nodes.
            if (Kernel_Delta || Four_Point_Delta || Modified_Two_Point_Delta)
            {
                xmin = x_int - 2;
                xmax = x_int + 2;
                ymin = y_int - 2;
                ymax = y_int + 2;
                zmin = z_int - 2;
                zmax = z_int + 2;
            }

            if (Two_Point_Delta)
            {
                xmin = x_int - 1;
                xmax = x_int + 1;
                ymin = y_int - 1;
                ymax = y_int + 1;
                zmin = z_int - 1;
                zmax = z_int + 1;
            }
            if (!X_Periodic)
            {
                if (xmin <= 0)
                {
                    xmin = 0;
                }
                if (xmax >= Nx - 1)
                {
                    xmax = Nx - 1;
                }
            }
            if (!Y_Periodic)
            {

                if (ymin <= 0)
                {
                    ymin = 0;
                }
                if (ymax >= Ny - 1)
                {
                    ymax = Ny - 1;
                }
            }
            if (!Z_Periodic)
            {
                if (zmin <= 0)
                {
                    zmin = 0;
                }
                if (zmax >= Nz - 1)
                {
                    zmax = Nz - 1;
                }
            }



        // Run over all neighboring fluid nodes.
        // In the case of the two-point interpolation, it is 2x2 fluid nodes.
        for (INT32 X = xmin; X <= xmax; ++X)
        {
            for (INT32 Y = ymin; Y <= ymax ; ++Y)
            {
                for (INT32 Z = zmin; Z <= zmax; ++Z)
                {
                    // Compute distance between object node and fluid lattice node.
                    dist_x = particle_list[np].node[n].x - X;
                    dist_y = particle_list[np].node[n].y - Y;
                    dist_z = particle_list[np].node[n].z - Z;
                    // Compute interpolation weights for x- and y-direction based on the distance.

                    #if UsePreCalculatedKernels
                    weight_x = kernel_weight(dist_x);
                    weight_y = kernel_weight(dist_y);
                    weight_z = kernel_weight(dist_z);
                    #else //Calculate Kernels here
                    if (Kernel_Delta)
                    {
                        if (abs(dist_x) <= 2)
                        {
                            weight_x = 0.25 * (1. + cos(M_PI * abs(dist_x) / 2));
                        }
                        else weight_x = 0.0;
                        if (abs(dist_y) <= 2)
                        {
                            weight_y = 0.25 * (1. + cos(M_PI * abs(dist_y) / 2));
                        }
                        else weight_y = 0.0;
                        if (abs(dist_z) <= 2)
                        {
                            weight_z = 0.25 * (1. + cos(M_PI * abs(dist_z) / 2));
                        }
                        else weight_z = 0.0;
                    }

                    if (Four_Point_Delta)
                    {
                        if (abs(dist_x) >= 0 && abs(dist_x) < 1)
                        {
                            weight_x = 0.125 * (3. - 2. * abs(dist_x) + sqrt(1 + 4. * abs(dist_x) - 4. * SQ(dist_x)));
                        }
                        else if (abs(dist_x) >= 1 && abs(dist_x) < 2)
                        {
                            weight_x = 0.125 * (5. - 2. * abs(dist_x) - sqrt(-7 + 12. * abs(dist_x) - 4. * SQ(dist_x)));
                        }
                        else weight_x = 0.0;

                        if (abs(dist_y) >= 0 && abs(dist_y) < 1)
                        {
                            weight_y = 0.125 * (3. - 2. * abs(dist_y) + sqrt(1 + 4. * abs(dist_y) - 4. * SQ(dist_y)));
                        }
                        else if (abs(dist_y) >= 1 && abs(dist_y) < 2)
                        {
                            weight_y = 0.125 * (5. - 2. * abs(dist_y) - sqrt(-7 + 12. * abs(dist_y) - 4. * SQ(dist_y)));
                        }
                        else weight_y = 0.0;

                        if (abs(dist_z) >= 0 && abs(dist_z) < 1)
                        {
                            weight_z = 0.125 * (3. - 2. * abs(dist_z) + sqrt(1 + 4. * abs(dist_z) - 4. * SQ(dist_z)));
                        }
                        else if (abs(dist_z) >= 1 && abs(dist_z) < 2)
                        {
                            weight_z = 0.125 * (5. - 2. * abs(dist_z) - sqrt(-7 + 12. * abs(dist_z) - 4. * SQ(dist_z)));
                        }
                        else weight_z = 0.0;
                    }

                    if (Two_Point_Delta)
                    {
                        if (abs(dist_x) <= 1)
                        {
                            weight_x = 1. - abs(dist_x);
                        }
                        else weight_x = 0.0;
                        if (abs(dist_y) <= 1)
                        {
                            weight_y = 1. - abs(dist_y);
                        }
                        else weight_y = 0.0;
                        if (abs(dist_z) <= 1)
                        {
                            weight_z = 1. - abs(dist_z);
                        }
                        else weight_z = 0.0;
                    }
                    if (Modified_Two_Point_Delta)
                    {
                        if (abs(dist_x) <= 0.5)
                        {
                            weight_x = 3. / 4. - (dist_x) * (dist_x);
                        }
                        else if (abs(dist_x) > 0.5 && abs(dist_x) <= 1.5)
                        {
                            weight_x = 9. / 8. - 3. * abs(dist_x) / 2. + SQ(dist_x) / 2.;
                        }
                        else weight_x = 0.0;

                        if (abs(dist_y) <= 0.5)
                        {
                            weight_y = 3. / 4. - (dist_y) * (dist_y);
                        }
                        else if (abs(dist_y) > 0.5 && abs(dist_y) <= 1.5)
                        {
                            weight_y = 9. / 8. - 3. * abs(dist_y) / 2. + SQ(dist_y) / 2.;
                        }
                        else weight_y = 0.0;

                        if (abs(dist_z) <= 0.5)
                        {
                            weight_z = 3. / 4. - (dist_z) * (dist_z);
                        }
                        else if (abs(dist_z) > 0.5 && abs(dist_z) <= 1.5)
                        {
                            weight_z = 9. / 8. - 3. * abs(dist_z) / 2. + SQ(dist_z) / 2.;
                        }
                        else weight_z = 0.0;
                    }
                    #endif

                    INT32 Xnew = paralLB.get_Xlocal_from_Xglobal((X + Nx) % Nx, this_processor_id);//(X + Nx) % Nx + num_buffer_layers_X[this_processor_id];
                    INT32 Ynew = paralLB.get_Ylocal_from_Yglobal((Y + Ny) % Ny, this_processor_id);//(Y + Ny) % Ny + num_buffer_layers_Y[this_processor_id];
                    INT32 Znew = paralLB.get_Zlocal_from_Zglobal((Z + Nz) % Nz, this_processor_id);//(Z + Nz) % Nz + num_buffer_layers_Z[this_processor_id];

                    const double wxyz = weight_x * weight_y * weight_z;

                    if (mdf_iteration == 0) {//First interpolation
                        #if !STORE_1D_VECTOR
                        particle_list[np].node[n].vel_x += ((flowSolverLB.velocity[0][Xnew][Ynew][Znew]) * wxyz);
                        particle_list[np].node[n].vel_y += ((flowSolverLB.velocity[1][Xnew][Ynew][Znew]) * wxyz);
                        particle_list[np].node[n].vel_z += ((flowSolverLB.velocity[2][Xnew][Ynew][Znew]) * wxyz);
                        particle_list[np].node[n].den   += ((flowSolverLB.density[Xnew][Ynew][Znew])     * wxyz);
                        if (Flow_With_Thermal_Effect)
                        {
                            particle_list[np].node[n].temperature += ((flowSolverLB.temperature[Xnew][Ynew][Znew]) * wxyz);
                        }
                        #endif
                        #if STORE_1D_VECTOR
                        INT32 id3, id4_0, id4_1, id4_2;
                        id3   = idx_3d_vector(Xnew, Ynew, Znew);
                        id4_0 = idx_4d_vector(0, Xnew, Ynew, Znew);
                        id4_1 = idx_4d_vector(1, Xnew, Ynew, Znew);
                        id4_2 = idx_4d_vector(2, Xnew, Ynew, Znew);
                        particle_list[np].node[n].vel_x += ((flowSolverLB.velocity[id4_0]) * wxyz);
                        particle_list[np].node[n].vel_y += ((flowSolverLB.velocity[id4_1]) * wxyz);
                        particle_list[np].node[n].vel_z += ((flowSolverLB.velocity[id4_2]) * wxyz);
                        particle_list[np].node[n].den += ((flowSolverLB.density[id3]) * wxyz);
                        if (Flow_With_Thermal_Effect)
                        {
                            particle_list[np].node[n].temperature += ((flowSolverLB.temperature[id3]) * wxyz);
                        }
                        #endif
                    }else{// Next mdf iterations
                        #if !STORE_1D_VECTOR
                        particle_list[np].node[n].vel_x += ((flowSolverLB.mdfVelocity[0][Xnew][Ynew][Znew]) * wxyz);
                        particle_list[np].node[n].vel_y += ((flowSolverLB.mdfVelocity[1][Xnew][Ynew][Znew]) * wxyz);
                        particle_list[np].node[n].vel_z += ((flowSolverLB.mdfVelocity[2][Xnew][Ynew][Znew]) * wxyz);
                        particle_list[np].node[n].den   += ((flowSolverLB.density[Xnew][Ynew][Znew])        * wxyz);
                        if (Flow_With_Thermal_Effect)
                        {
                            particle_list[np].node[n].temperature += ((flowSolverLB.mdfTemperature[Xnew][Ynew][Znew]) * wxyz);
                        }
                        #endif
                        #if STORE_1D_VECTOR
                        INT32 id3, id4_0, id4_1, id4_2;
                        id3   = idx_3d_vector(Xnew, Ynew, Znew);
                        id4_0 = idx_4d_vector(0, Xnew, Ynew, Znew);
                        id4_1 = idx_4d_vector(1, Xnew, Ynew, Znew);
                        id4_2 = idx_4d_vector(2, Xnew, Ynew, Znew);
                        particle_list[np].node[n].vel_x += ((flowSolverLB.mdfVelocity[id4_0]) * wxyz);
                        particle_list[np].node[n].vel_y += ((flowSolverLB.mdfVelocity[id4_1]) * wxyz);
                        particle_list[np].node[n].vel_z += ((flowSolverLB.mdfVelocity[id4_2]) * wxyz);
                        particle_list[np].node[n].den   += ((flowSolverLB.density[id3])       * wxyz);
                        if (Flow_With_Thermal_Effect)
                        {
                            particle_list[np].node[n].temperature += ((flowSolverLB.mdfTemperature[id3]) * wxyz);
                        }
                        #endif
                    }
                }
            }
        }
    }
    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: SREAD
//-----------------------------------------------------------------------------------------------------------------------------
//
// The node forces are spread to the fluid nodes via IBM.
// The two-point interpolation stencil (bi-linear interpolation) is used in the present code.
// It may be replaced by a higher-order interpolation.
void ParticleResolved::spread(INT32 np, FlowSolver& flowSolverLB, const ParallelMPI& paralLB)
{
    INT32 x_int{};
    INT32 y_int{};
    INT32 z_int{};
    INT32 xmin{};
    INT32 xmax{};
    INT32 ymin{};
    INT32 ymax{};
    INT32 zmin{};
    INT32 zmax{};
    double dist_x{};
    double dist_y{};
    double dist_z{};
    double weight_x{};
    double weight_y{};
    double weight_z{};

    /// Spread forces
    // Run over all object nodes.
    for (INT32 n = 0; n < particle_num_nodes; ++n)
    {

        // Identify the lowest fluid lattice node in interpolation range.
        // 'Lowest' means: its x- and y-values are the smallest.
        // The other fluid nodes in range have coordinates
        // (x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).

            x_int = floor(particle_list[np].node[n].x);
            y_int = floor(particle_list[np].node[n].y);
            z_int = floor(particle_list[np].node[n].z);
            // Run over all neighboring fluid nodes.
            // In the case of the two-point interpolation, it is 2x2 fluid nodes.
            if (Kernel_Delta || Four_Point_Delta || Modified_Two_Point_Delta)
            {
                xmin = x_int - 2;
                xmax = x_int + 2;
                ymin = y_int - 2;
                ymax = y_int + 2;
                zmin = z_int - 2;
                zmax = z_int + 2;
            }

            if (Two_Point_Delta)
            {
                xmin = x_int - 1;
                xmax = x_int + 1;
                ymin = y_int - 1;
                ymax = y_int + 1;
                zmin = z_int - 1;
                zmax = z_int + 1;
            }
            if (!X_Periodic)
            {
                if (xmin <= 0)
                {
                    xmin = 0;
                }
                if (xmax >= Nx - 1)
                {
                    xmax = Nx - 1;
                }
            }
            if (!Y_Periodic)
            {
                if (ymin <= 0)
                {
                    ymin = 0;
                }
                if (ymax >= Ny - 1)
                {
                    ymax = Ny - 1;
                }
            }
            if (!Z_Periodic)
            {
                if (zmin <= 0)
                {
                    zmin = 0;
                }
                if (zmax >= Nz - 1)
                {
                    zmax = Nz - 1;
                }
            }


        // Run over all neighboring fluid nodes.
        // In the case of the two-point interpolation, it is 2x2 fluid nodes.
        for (INT32 X = xmin; X <= xmax; ++X)
        {
            for (INT32 Y = ymin; Y <= ymax; ++Y)
            {
                for (INT32 Z = zmin; Z <= zmax; ++Z)
                {
                    // Compute distance between object node and fluid lattice node.
                    dist_x = particle_list[np].node[n].x - X;
                    dist_y = particle_list[np].node[n].y - Y;
                    dist_z = particle_list[np].node[n].z - Z;
                    // Compute interpolation weights for x- and y-direction based on the distance.
                    #if UsePreCalculatedKernels
                    weight_x = kernel_weight(dist_x);
                    weight_y = kernel_weight(dist_y);
                    weight_z = kernel_weight(dist_z);
                    #else // UsePreCalculatedKernels
                    if (Kernel_Delta)
                    {
                        if (abs(dist_x) <= 2)
                        {
                            weight_x = 0.25 * (1. + cos(M_PI * abs(dist_x) / 2));
                        }
                        else weight_x = 0.0;
                        if (abs(dist_y) <= 2)
                        {
                            weight_y = 0.25 * (1. + cos(M_PI * abs(dist_y) / 2));
                        }
                        else weight_y = 0.0;
                        if (abs(dist_z) <= 2)
                        {
                            weight_z = 0.25 * (1. + cos(M_PI * abs(dist_z) / 2));
                        }
                        else weight_z = 0.0;
                    }

                    if (Four_Point_Delta)
                    {
                        if (abs(dist_x) >= 0 && abs(dist_x) < 1)
                        {
                            weight_x = 0.125 * (3. - 2. * abs(dist_x) + sqrt(1 + 4. * abs(dist_x) - 4. * SQ(dist_x)));
                        }
                        else if (abs(dist_x) >= 1 && abs(dist_x) < 2)
                        {
                            weight_x = 0.125 * (5. - 2. * abs(dist_x) - sqrt(-7 + 12. * abs(dist_x) - 4. * SQ(dist_x)));
                        }
                        else weight_x = 0.0;

                        if (abs(dist_y) >= 0 && abs(dist_y) < 1)
                        {
                            weight_y = 0.125 * (3. - 2. * abs(dist_y) + sqrt(1 + 4. * abs(dist_y) - 4. * SQ(dist_y)));
                        }
                        else if (abs(dist_y) >= 1 && abs(dist_y) < 2)
                        {
                            weight_y = 0.125 * (5. - 2. * abs(dist_y) - sqrt(-7 + 12. * abs(dist_y) - 4. * SQ(dist_y)));
                        }
                        else weight_y = 0.0;

                        if (abs(dist_z) >= 0 && abs(dist_z) < 1)
                        {
                            weight_z = 0.125 * (3. - 2. * abs(dist_z) + sqrt(1 + 4. * abs(dist_z) - 4. * SQ(dist_z)));
                        }
                        else if (abs(dist_z) >= 1 && abs(dist_z) < 2)
                        {
                            weight_z = 0.125 * (5. - 2. * abs(dist_z) - sqrt(-7 + 12. * abs(dist_z) - 4. * SQ(dist_z)));
                        }
                        else weight_z = 0.0;
                    }

                    if (Two_Point_Delta)
                    {
                        if (abs(dist_x) <= 1)
                        {
                            weight_x = 1. - abs(dist_x);
                        }
                        else weight_x = 0.0;
                        if (abs(dist_y) <= 1)
                        {
                            weight_y = 1. - abs(dist_y);
                        }
                        else weight_y = 0.0;
                        if (abs(dist_z) <= 1)
                        {
                            weight_z = 1. - abs(dist_z);
                        }
                        else weight_z = 0.0;
                    }
                    if (Modified_Two_Point_Delta)
                    {
                        if (abs(dist_x) <= 0.5)
                        {
                            weight_x = 3. / 4. - (dist_x) * (dist_x);
                        }
                        else if (abs(dist_x) > 0.5 && abs(dist_x) <= 1.5)
                        {
                            weight_x = 9. / 8. - 3. * abs(dist_x) / 2. + SQ(dist_x) / 2.;
                        }
                        else weight_x = 0.0;

                        if (abs(dist_y) <= 0.5)
                        {
                            weight_y = 3. / 4. - (dist_y) * (dist_y);
                        }
                        else if (abs(dist_y) > 0.5 && abs(dist_y) <= 1.5)
                        {
                            weight_y = 9. / 8. - 3. * abs(dist_y) / 2. + SQ(dist_y) / 2.;
                        }
                        else weight_y = 0.0;

                        if (abs(dist_z) <= 0.5)
                        {
                            weight_z = 3. / 4. - (dist_z) * (dist_z);
                        }
                        else if (abs(dist_z) > 0.5 && abs(dist_z) <= 1.5)
                        {
                            weight_z = 9. / 8. - 3. * abs(dist_z) / 2. + SQ(dist_z) / 2.;
                        }
                        else weight_z = 0.0;
                    }
                    #endif

                    INT32 Xnew = paralLB.get_Xlocal_from_Xglobal((X + Nx) % Nx, this_processor_id);//(X + Nx) % Nx + num_buffer_layers_X[this_processor_id];
                    INT32 Ynew = paralLB.get_Ylocal_from_Yglobal((Y + Ny) % Ny, this_processor_id);//(Y + Ny) % Ny + num_buffer_layers_Y[this_processor_id];
                    INT32 Znew = paralLB.get_Zlocal_from_Zglobal((Z + Nz) % Nz, this_processor_id);//(Z + Nz) % Nz + num_buffer_layers_Z[this_processor_id];

                    const double wxyz = weight_x * weight_y * weight_z;

                    #if !STORE_1D_VECTOR
                    flowSolverLB.force[0][Xnew][Ynew][Znew] += (particle_list[np].node[n].force_x * wxyz);
                    flowSolverLB.force[1][Xnew][Ynew][Znew] += (particle_list[np].node[n].force_y * wxyz);
                    flowSolverLB.force[2][Xnew][Ynew][Znew] += (particle_list[np].node[n].force_z * wxyz);
                    if (Flow_With_Thermal_Effect)
                    {
                        flowSolverLB.force_thermal[Xnew][Ynew][Znew] += (particle_list[np].node[n].force_thermal * wxyz);
                    }
                    #endif
                    #if STORE_1D_VECTOR
                    INT32 id3, id4_0, id4_1, id4_2;
                    id3   = idx_3d_vector(Xnew, Ynew, Znew);
                    id4_0 = idx_4d_vector(0, Xnew, Ynew, Znew);
                    id4_1 = idx_4d_vector(1, Xnew, Ynew, Znew);
                    id4_2 = idx_4d_vector(2, Xnew, Ynew, Znew);
                    flowSolverLB.force[id4_0] += (particle_list[np].node[n].force_x * wxyz);
                    flowSolverLB.force[id4_1] += (particle_list[np].node[n].force_y * wxyz);
                    flowSolverLB.force[id4_2] += (particle_list[np].node[n].force_z * wxyz);
                    if (Flow_With_Thermal_Effect)
                    {
                        flowSolverLB.force_thermal[id3] += (particle_list[np].node[n].force_thermal * wxyz);
                    }
                    #endif // STORE_1D_VECTOR

                }
            }
        }
    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: keep_particle_at_center
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::keep_particle_at_center(INT32 npn, FlowSolver& flowSolverLB){

    /// To keep particle at center along the X direction (Does not work with Recovery_ON currently)
    const INT32 displacement_before_shift = 5;// How many cell particle moves before shift of domain happens
    const INT32 desire_cen_x  = INT32(particle_center_x);// To always keep the particle at center
    total_displacement_x   += 0.5 * (particle_list[npn].center[0].vel_x + particle_list[npn].center[0].tempx[0]);// Just for saving in write_particle_data
    partial_displacement_x += 0.5 * (particle_list[npn].center[0].vel_x + particle_list[npn].center[0].tempx[0]);

    #if !STORE_1D_VECTOR

    std::vector<std::vector<std::vector<double>>> swap_den_mov      (length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER])));
    std::vector<std::vector<std::vector<double>>> swap_force_thermal(length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER])));
    std::vector<std::vector<std::vector<double>>> swap_temperature  (length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER])));

    std::vector<std::vector<std::vector<std::vector<double>>>> swap_vel_mov  (3, std::vector<std::vector<std::vector<double>>>(length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER]))));
    std::vector<std::vector<std::vector<std::vector<double>>>> swap_force_mov(3, std::vector<std::vector<std::vector<double>>>(length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER]))));

    std::vector<std::vector<std::vector<std::vector<double>>>> swap_pop  (stencil_directions_F, std::vector<std::vector<std::vector<double>>>(length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER]))));
    std::vector<std::vector<std::vector<std::vector<double>>>> swap_pop_t(stencil_directions_T, std::vector<std::vector<std::vector<double>>>(length_with_buffer_X[MASTER], std::vector<std::vector<double>>(length_with_buffer_Y[MASTER], std::vector<double>(length_with_buffer_Z[MASTER]))));

    if (abs (partial_displacement_x) > displacement_before_shift){

        /// Make a copy of data
        swap_den_mov       = flowSolverLB.density;
        swap_force_thermal = flowSolverLB.force_thermal;
        swap_temperature   = flowSolverLB.temperature;
        swap_vel_mov       = flowSolverLB.velocity;
        swap_force_mov     = flowSolverLB.force;
        swap_pop           = flowSolverLB.pop_F;
        swap_pop_t         = flowSolverLB.pop_T;

        INT32 start_X = desire_cen_x - 500;
        INT32 end_X   = desire_cen_x + 400;

        for (INT32 X = start_X; X < end_X; ++X)
        {
            for (INT32 Y = 0; Y < length_with_buffer_Y[MASTER]; ++Y)
            {
                for (INT32 Z = 0; Z < length_with_buffer_Z[MASTER]; ++Z)
                {
                    INT32 copy_from_X = INT32(floor(particle_list[npn].center[0].x)) - desire_cen_x + X;

                    flowSolverLB.velocity[0][X][Y][Z]   = swap_vel_mov[0][copy_from_X][Y][Z];
                    flowSolverLB.velocity[1][X][Y][Z]   = swap_vel_mov[1][copy_from_X][Y][Z];
                    flowSolverLB.velocity[2][X][Y][Z]   = swap_vel_mov[2][copy_from_X][Y][Z];
                    flowSolverLB.force[0][X][Y][Z]      = swap_force_mov[0][copy_from_X][Y][Z];
                    flowSolverLB.force[1][X][Y][Z]      = swap_force_mov[1][copy_from_X][Y][Z];
                    flowSolverLB.force[2][X][Y][Z]      = swap_force_mov[2][copy_from_X][Y][Z];
                    flowSolverLB.density[X][Y][Z]       = swap_den_mov[copy_from_X][Y][Z];
                    flowSolverLB.temperature[X][Y][Z]   = swap_temperature[copy_from_X][Y][Z];
                    flowSolverLB.force_thermal[X][Y][Z] = swap_force_thermal[copy_from_X][Y][Z];

                    for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i) {
                        flowSolverLB.pop_F[c_i][X][Y][Z]  = swap_pop[c_i][copy_from_X][Y][Z];
                    }
                    for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i) {
                        flowSolverLB.pop_T[c_i][X][Y][Z] = swap_pop_t[c_i][copy_from_X][Y][Z];
                    }
                }
            }
        }

        particle_list[npn].center[0].x -= floor(partial_displacement_x);//Now shift the particle position
        partial_displacement_x         -= floor(partial_displacement_x);//Now reset partial_displacement_x

    } //End of: if (abs (partial_displacement_x) > displacement_before_shift)

    #endif  // End of: !STORE_1D_VECTOR




   /* std::vector<std::vector<double>> swap_vel_mov_x(Nx, std::vector<double  >(Ny));// Vectors are defined on the heap
    std::vector<std::vector<double>> swap_vel_mov_y(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_den_mov(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_force_mov_x(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_force_mov_y(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_force_thermal(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> swap_temperature(Nx, std::vector<double>(Ny));

    std::vector<std::vector<std::vector<double>>> swap_pop(
        stencil_directions_F, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny))
    );

    std::vector<std::vector<std::vector<double>>> swap_pop_t(
        stencil_directions_T, std::vector<std::vector<double>>(Nx, std::vector<double>(Ny))
    );

    if (abs (partial_displacement_x) > displacement_before_shift){

        for (INT32 i = 0; i < Nx; ++i) {
            std::copy(velocity_x[i], velocity_x[i] + Ny, swap_vel_mov_x[i].begin());
            std::copy(velocity_y[i], velocity_y[i] + Ny, swap_vel_mov_y[i].begin());
            std::copy(density[i], density[i] + Ny, swap_den_mov[i].begin());
            std::copy(force_x[i], force_x[i] + Ny, swap_force_mov_x[i].begin());
            std::copy(force_y[i], force_y[i] + Ny, swap_force_mov_y[i].begin());
            std::copy(force_thermal[i], force_thermal[i] + Ny, swap_force_thermal[i].begin());
            std::copy(temperature[i], temperature[i] + Ny, swap_temperature[i].begin());
        }

        for (INT32 d = 0; d < stencil_directions_F; ++d) {
            for (INT32 i = 0; i < Nx; ++i) {
                std::copy(pop_F[d][i], pop_F[d][i] + Ny, swap_pop[d][i].begin());
            }
        }

        for (INT32 d = 0; d < stencil_directions_T; ++d) {
            for (INT32 i = 0; i < Nx; ++i) {
                std::copy(pop_T[d][i], pop_T[d][i] + Ny, swap_pop_t[d][i].begin());
            }
        }

        INT32 start_X = desire_cen_x - 500;
        INT32 end_X   = desire_cen_x + 400;

        for (INT32 X = start_X; X < end_X; ++X)
        {
            for (INT32 Y = 0; Y < Ny; ++Y)
            {
                INT32 copy_from_X = INT32(floor(particle[npn].center[0].x)) - desire_cen_x + X;

                velocity_x[X][Y]    = swap_vel_mov_x[copy_from_X][Y];
                velocity_y[X][Y]    = swap_vel_mov_y[copy_from_X][Y];
                force_x[X][Y]       = swap_force_mov_x[copy_from_X][Y];
                force_y[X][Y]       = swap_force_mov_y[copy_from_X][Y];
                density[X][Y]       = swap_den_mov[copy_from_X][Y];
                temperature[X][Y]   = swap_temperature[copy_from_X][Y];
                force_thermal[X][Y] = swap_force_thermal[copy_from_X][Y];
                for (INT32 c_i = 0; c_i < stencil_directions_F; ++c_i) {
                    pop_F[c_i][X][Y]  = swap_pop[c_i][copy_from_X][Y];
                }
                for (INT32 c_i = 0; c_i < stencil_directions_T; ++c_i) {
                    pop_T[c_i][X][Y] = swap_pop_t[c_i][copy_from_X][Y];
                }
            }
        }

        particle[npn].center[0].x -= floor(partial_displacement_x);//Now shift the particle position
        partial_displacement_x    -= floor(partial_displacement_x);//Now reset partial_displacement_x
    }*/
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: calculate_force_modified_indices
//-----------------------------------------------------------------------------------------------------------------------------
//
/// To calculate force_modified_indices
void ParticleResolved::calculate_force_modified_indices(FlowSolver& flowSolverLB, const ParallelMPI& paralLB)
{
    INT32 x_int{};
    INT32 y_int{};
    INT32 z_int{};
    INT32 xmin{};
    INT32 xmax{};
    INT32 ymin{};
    INT32 ymax{};
    INT32 zmin{};
    INT32 zmax{};

    // Run over all object nodes.
    for (INT32 np = 0; np < num_particles; ++np)
    {
        for (INT32 n = 0; n < particle_num_nodes; ++n)
        {

            if (Kernel_Delta || Two_Point_Delta || Four_Point_Delta	|| Modified_Two_Point_Delta)
            {
                x_int = floor(particle_list[np].node[n].x);
                y_int = floor(particle_list[np].node[n].y);
                z_int = floor(particle_list[np].node[n].z);
                // Run over all neighboring fluid nodes.
                // In the case of the two-point interpolation, it is 2x2 fluid nodes.
                if (Kernel_Delta || Four_Point_Delta || Modified_Two_Point_Delta)
                {
                    xmin = x_int - 2;
                    xmax = x_int + 2;
                    ymin = y_int - 2;
                    ymax = y_int + 2;
                    zmin = z_int - 2;
                    zmax = z_int + 2;
                }

                if (Two_Point_Delta)
                {
                    xmin = x_int - 1;
                    xmax = x_int + 1;
                    ymin = y_int - 1;
                    ymax = y_int + 1;
                    zmin = z_int - 1;
                    zmax = z_int + 1;
                }
                if (!X_Periodic)
                {
                    if (xmin <= 0)
                    {
                        xmin = 0;
                    }
                    if (xmax >= Nx - 1)
                    {
                        xmax = Nx - 1;
                    }
                }
                if (!Y_Periodic)
                {
                    if (ymin <= 0)
                    {
                        ymin = 0;
                    }
                    if (ymax >= Ny - 1)
                    {
                        ymax = Ny - 1;
                    }
                }
                if (!Z_Periodic)
                {
                    if (zmin <= 0)
                    {
                        zmin = 0;
                    }
                    if (zmax >= Nz - 1)
                    {
                        zmax = Nz - 1;
                    }
                }
            }


            // Run over all neighboring fluid nodes.
            // In the case of the two-point interpolation, it is 2x2 fluid nodes.
            for (INT32 X = xmin; X <= xmax; ++X)
            {
                for (INT32 Y = ymin; Y <= ymax; ++Y)
                {
                    for (INT32 Z = zmin; Z <= zmax; ++Z)
                    {
                        INT32 Xnew = paralLB.get_Xlocal_from_Xglobal((X + Nx) % Nx, this_processor_id);//(X + Nx) % Nx + num_buffer_layers_X[this_processor_id];
                        INT32 Ynew = paralLB.get_Ylocal_from_Yglobal((Y + Ny) % Ny, this_processor_id);//(Y + Ny) % Ny + num_buffer_layers_Y[this_processor_id];
                        INT32 Znew = paralLB.get_Zlocal_from_Zglobal((Z + Nz) % Nz, this_processor_id);//(Z + Nz) % Nz + num_buffer_layers_Z[this_processor_id];

                        #if !STORE_1D_VECTOR
                        flowSolverLB.force_modified_indices.emplace_back(Xnew, Ynew, Znew);
                        #endif
                        #if STORE_1D_VECTOR
                        INT32 id3;
                        id3   = idx_3d_vector(Xnew, Ynew, Znew);
                        flowSolverLB.force_modified_indices.push_back(id3);
                        #endif // STORE_1D_VECTOR

                    }
                }
            }

        }// for (INT32 n = 0; n < particle_num_nodes; ++n)
    }//for (INT32 np = 0; np < num_particles; ++np)

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: scatter_particles_center_data_from_MASTER
//-----------------------------------------------------------------------------------------------------------------------------


// Assume your structs are already defined:
// struct center_struct { ... };
// struct particle_struct { std::vector<center_struct> center; std::vector<node_struct> node; ... };
// std::vector<particle_struct> particle_list;

//--------------------------------------------------
// Update particle centers from MASTER to all workers
//--------------------------------------------------
void ParticleResolved::scatter_particles_center_data_from_MASTER()//, INT32 MASTER, MPI_Comm comm)
{
    // Make a contiguous buffer of center_struct
    std::vector<center_struct> centers(particle_list.size());

    if (this_processor_id == MASTER) {
        // Fill buffer with current centers
        for (INT32 i = 0; i < particle_list.size(); ++i) {
            centers[i] = particle_list[i].center[0]; // assuming center vector size = 1
        }
    }

    // Broadcast the centers buffer as raw bytes
    MPI_Bcast(centers.data(),
              static_cast<INT32>(centers.size() * sizeof(center_struct)),
              MPI_BYTE,
              MASTER,
              MPI_COMM_WORLD);

    // On all non-master ranks, update local particle_list
    if (this_processor_id != MASTER) {
        for (INT32 i = 0; i < particle_list.size(); ++i) {
            particle_list[i].center[0] = centers[i];
        }
    }
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_VTK
//-----------------------------------------------------------------------------------------------------------------------------
//
// The particle state (node positions) is writen to a VTK file at each t_vtk step.
// The following code is designed in such a way that the file can be read by ParaView.

void ParticleResolved::write_particle_vtk(INT32 time)
{

    for (INT32 npp = 0; npp < num_particles; ++npp)
    {

        if (VTK_BINARY || VTK_ASCII)
        {

            /// Create filename

            stringstream output_filename;
            output_filename << "FlowLB_Results/vtk_particle/particle_" << npp + 1 << "_t" << time << ".vtk";
            ofstream output_file;

            /// Open file

            output_file.open(output_filename.str().c_str());

            /// Write VTK header

            output_file << "# vtk DataFile Version 3.0\n";
            output_file << "particle_state\n";
            output_file << "ASCII\n";
            output_file << "DATASET POLYDATA\n";

            /// Write node positions
            output_file << "POINTS " << particle_num_nodes << " float\n";

            for (INT32 n = 0; n < particle_num_nodes; ++n)
            {
                output_file << particle_list[npp].node[n].x << " " << particle_list[npp].node[n].y << " " << particle_list[npp].node[n].z << "\n";
            }

            if (MOVING_SPHERE || STATIONARY_SPHERE)
            {
                output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
                output_file << particle_num_nodes << "\n";
                for (INT32 n = 0; n < particle_num_nodes; ++n)
                {
                    output_file << n << "\n";
                }

                // Below part is useful for SPHERE to create surface using Glyph for canter in Paraview
                stringstream output_filename1;
                output_filename1 << "FlowLB_Results/vtk_particle/particle_center_" << npp + 1 << "_t" << time << ".vtk";
                ofstream output_file1;

                /// Open file
                output_file1.open(output_filename1.str().c_str());

                /// Write VTK header
                output_file1 << "# vtk DataFile Version 3.0\n";
                output_file1 << "particle_state\n";
                output_file1 << "ASCII\n";
                output_file1 << "DATASET POLYDATA\n";

                output_file1 << "POINTS " << "2" << " float\n";
                output_file1 << particle_list[npp].center[0].x << " " << particle_list[npp].center[0].y << " " << particle_list[npp].center[0].z << "\n";
                output_file1 << particle_list[npp].center[0].x + 0.000001 << " " << particle_list[npp].center[0].y << " " << particle_list[npp].center[0].z << "\n";

                /// Close file
                output_file1.close();
            }

            if (MOVING_SPHEROID || STATIONARY_SPHEROID)
            {
                // Create triangle surfaces
                output_file << "POLYGONS " << number_triangles << " " << 4 * number_triangles << "\n";
                for (INT32 n = 0; n < number_triangles; ++n)
                {
                    output_file << "3 " << triangleset1[n] << " " << triangleset2[n] << " " << triangleset3[n] << "\n";
                }
            }

            if (STATIONARY_CYLINDER || STATIONARY_ELLIPSE || MOVING_CYLINDER || MOVING_ELLIPSE)
            {
                /// Write lines between neighboring nodes
                output_file << "LINES " << particle_num_nodes << " " << 3 * particle_num_nodes << "\n";
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    output_file << "2 " << n << " " << (n + 1) % particle_num_nodes << "\n";
                }

                /// Write vertices
                output_file << "VERTICES 1 " << particle_num_nodes + 1 << "\n";
                output_file << particle_num_nodes << " ";
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    output_file << n << " ";
                }
            }
            /// Close file
            output_file.close();
        }

        if (CASE)
        {
            //unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);// number are shown as 2.0e02 instead of 2.0e002 (necessay for geo and data files)


            /////////// Writing CASE file (.case) ///////////////////
            /////////////////////////////////////////////////////////
            /// Create filename

            stringstream output_filename;
            output_filename << "FlowLB_Results/vtk_particle/particle_" << npp + 1 << ".case";
            ofstream output_file;

            /// Open file
            output_file.open(output_filename.str().c_str());

            /// Write VTK header
            output_file << "# EnSight Gold Model\n";
            output_file << "# Created by FlowLB\n\n";
            output_file << "FORMAT\n";
            output_file << "type:                   ensight gold\n\n";
            output_file << "GEOMETRY\n";
            output_file << "model:                 1 particle_" << npp + 1 << "_******.geo\n\n";
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


            /////////// Writing GEOMETRY file (.geo) ///////////////////
            ////////////////////////////////////////////////////////////

            double min_node_x, min_node_y, min_node_z;
            double max_node_x, max_node_y, max_node_z;
            min_node_x = particle_list[npp].node[0].x;
            min_node_y = particle_list[npp].node[0].y;
            min_node_z = particle_list[npp].node[0].z;
            max_node_x = particle_list[npp].node[0].x;
            max_node_y = particle_list[npp].node[0].y;
            max_node_z = particle_list[npp].node[0].z;


            for (INT32 n = 0; n < particle_num_nodes; ++n)
            {
                if (particle_list[npp].node[n].x < min_node_x)
                {
                    min_node_x = particle_list[npp].node[n].x;
                }
                if (particle_list[npp].node[n].y < min_node_y)
                {
                    min_node_y = particle_list[npp].node[n].y;
                }
                if (particle_list[npp].node[n].z < min_node_z)
                {
                    min_node_z = particle_list[npp].node[n].z;
                }
                if (particle_list[npp].node[n].x > max_node_x)
                {
                    max_node_x = particle_list[npp].node[n].x;
                }
                if (particle_list[npp].node[n].y > max_node_y)
                {
                    max_node_y = particle_list[npp].node[n].y;
                }
                if (particle_list[npp].node[n].z > max_node_z)
                {
                    max_node_z = particle_list[npp].node[n].z;
                }
            }
            stringstream output_filename1;
            output_filename1 << "FlowLB_Results/vtk_particle/particle_" << npp + 1 << "_" << setfill('0') << setw(6) << time << ".geo";//Geometry file
            ofstream output_file1;

            /// Open file
            output_file1.open(output_filename1.str().c_str());

            /// Write VTK header
            output_file1 << "EnSight Model Geometry File\n";
            output_file1 << "EnSight 10.1.1\n";
            output_file1 << "node id given\n";
            output_file1 << "element id given\n";
            output_file1 << "extents\n";
            output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << min_node_x << std::scientific << std::setw(12) << std::fixed << std::setprecision(5) << max_node_x << "\n";
            output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << min_node_y << std::scientific << std::setw(12) << std::fixed << std::setprecision(5) << max_node_y << "\n";
            output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << min_node_z << std::scientific << std::setw(12) << std::fixed << std::setprecision(5) << max_node_z << "\n";
            output_file1 << "part\n";
            output_file1 << std::setw(10) << 1 << "\n";
            output_file1 << "Particle" << npp + 1 << "\n";
            output_file1 << "coordinates\n";
            output_file1 << std::setw(10) << particle_num_nodes << "\n";
            for (INT32 itm = 1; itm <= particle_num_nodes; itm++)
            {
                output_file1 << std::setw(10) << itm << "\n";
            }

            for (INT32 i = 0; i < particle_num_nodes; i++)
            {
                output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << particle_list[npp].node[i].x << "\n";
            }
            for (INT32 i = 0; i < particle_num_nodes; i++)
            {
                output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << particle_list[npp].node[i].y << "\n";
            }
            for (INT32 i = 0; i < particle_num_nodes; i++)
            {
                output_file1 << std::scientific << std::setw(12) << std::setprecision(5) << std::fixed << particle_list[npp].node[i].z << "\n";
            }


            if (MOVING_SPHEROID || STATIONARY_SPHEROID)
            {
                output_file1 << "tria3\n";
                output_file1 << std::setw(10) << number_triangles << "\n";
                for (INT32 itm = 1; itm <= number_triangles; itm++)
                {
                    output_file1 << std::setw(10) << itm << "\n";
                }

                // Create triangle surfaces
                for (INT32 n = 0; n < number_triangles; ++n)
                {
                    output_file1 << std::setw(10) << triangleset1[n] + 1 << std::setw(10) << triangleset2[n] + 1 << std::setw(10) << triangleset3[n] + 1 << "\n";
                }
            }

            if (STATIONARY_CYLINDER || STATIONARY_ELLIPSE || MOVING_CYLINDER || MOVING_ELLIPSE)
            {
                output_file1 << "bar2\n";
                output_file1 << std::setw(10) << particle_num_nodes << "\n";
                for (INT32 itm = 1; itm <= particle_num_nodes; itm++){
                    output_file1 << std::setw(10) << itm << "\n";
                }

                for (INT32 I = 1; I < particle_num_nodes; I++){
                    output_file1 << std::setw(10) << I << std::setw(10) << I + 1 << "\n";
                }
                output_file1 << std::setw(10) << particle_num_nodes << std::setw(10) << 1 << "\n";
            }

            /// Close file
            output_file1.close();
            /////////// End of Writing GEOMETRY file (.geo) ///////////////////
            ///////////////////////////////////////////////////////////////////

        }//End of CASE file for Ensight
    }
    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_DATA
//-----------------------------------------------------------------------------------------------------------------------------
//
// The following PARTICLE quantities are written to the disk at each t_data step:
// - drag and lift forces (x- and y-components of the force)
// - object center position (x- and y-components)
// - object center velocity (x- and y-components)
// The data file is readable by gnuplot

void ParticleResolved::write_particle_data(INT32 time)
{

    for (INT32 npp = 0; npp < num_particles; ++npp)
    {

        /// Create filename
        ostringstream oss;
        oss << npp + 1;
        string npstr = oss.str();//convert number to string
        string output_filename("FlowLB_Results/Data_particle/data_particle_");
        output_filename += npstr;
        output_filename += ".dat";
        //ofstream output_file;

            // Check if file already exists or is empty
        bool write_header = false;
        {
            std::ifstream infile(output_filename.c_str());
            if (!infile.good() || infile.peek() == std::ifstream::traits_type::eof()) {
                write_header = true; // file doesn’t exist or is empty
            }
        }

        /// Open file for append
        std::ofstream output_file(output_filename.c_str(), std::ios::app);

        if (write_header) {
            if (!Flow_With_Thermal_Effect)
            {
                output_file << "# time, fx_surf, fy_surf, fz_surf, center.x, center.y, center.z, "
               "vel_x, vel_y, vel_z, omgX, omgY, omgZ, tetaX, tetaY, tetaZ \n";
            }
            if (Flow_With_Thermal_Effect)
            {
                output_file << "# time, fx_surf, fy_surf, fz_surf, center.x, center.y, center.z, "
               "vel_x, vel_y, vel_z, omgX, omgY, omgZ, tetaX, tetaY, tetaZ, Q_1/4, Q_2/4, Q_3/4, Q_4/4, Q_tot, Temperature\n";
            }
        }

        /// Write data
        output_file << time << " "; // time step
        output_file << particle_list[npp].center[0].fx_surf << " "; // drag force
        output_file << particle_list[npp].center[0].fy_surf << " "; // lift force
        output_file << particle_list[npp].center[0].fz_surf << " "; // lift force
        output_file << particle_list[npp].center[0].x << " "; // center position (x-component)
        //output_file << particle_center_x + total_displacement_x<< " "; // center position (x-component): Activate for keep_particle_at_center instead of above line
        output_file << particle_list[npp].center[0].y << " "; // center position (y-component)
        output_file << particle_list[npp].center[0].z << " "; // center position (y-component)
        output_file << particle_list[npp].center[0].vel_x << " "; // center velocity (x-component)
        output_file << particle_list[npp].center[0].vel_y << " "; // center velocity (y-component)
        output_file << particle_list[npp].center[0].vel_z << " "; // center velocity (y-component)
        output_file << particle_list[npp].center[0].omgX << " "; // center velocity (x-component)
        output_file << particle_list[npp].center[0].omgY << " "; // center velocity (y-component)
        output_file << particle_list[npp].center[0].omgZ << " "; // center velocity (y-component)
        output_file << particle_list[npp].center[0].tetaX << " "; // Angle with respect to Horizontal line
        output_file << particle_list[npp].center[0].tetaY << " "; // Angle with respect to Horizontal line
        output_file << particle_list[npp].center[0].tetaZ << " "; // Angle with respect to Horizontal line

        if (!Flow_With_Thermal_Effect)
        {
            output_file << "\n";
        }
        if (Flow_With_Thermal_Effect)
        {
            output_file << particle_list[npp].node[(int)(0 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
            output_file << particle_list[npp].node[(int)(1 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
            output_file << particle_list[npp].node[(int)(2 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
            output_file << particle_list[npp].node[(int)(3 * particle_num_nodes / 4)].total_node_mdf_force_thermal << " "; //Heat transfer at a specific node
            output_file << particle_list[npp].center[0].heat_surf << " "; //Total heat trasnferred
            output_file << particle_list[npp].center[0].temperature << "\n"; // center velocity (y-component)
        }
        /// Close file
        output_file.close();

    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: WRITE_PARTICLE_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::write_particle_recovery(INT32 time)
{

    if (num_particles > 0)
    {
        for (INT32 npn = 0; npn < num_particles; ++npn)
        {
            ostringstream oss;
            oss << npn + 1;
            string npstr = oss.str();//convert number to string
            string output_filename("FlowLB_Results/data_particle_recovery_");
            output_filename += npstr;
            output_filename += ".dat";
            ofstream output_file;

            /// Open file
            output_file.open(output_filename.c_str(), fstream::trunc); // using "trunc" means previous file (if exists) is deleted

            /// Write particle data
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].particle_mass << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].parvol << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].mom_inertia_x << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].mom_inertia_y << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].mom_inertia_z << " ";
            output_file << std::fixed << std::setprecision(12) << particle_num_nodes << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].particle_perimeter << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].x << " "; // center position (x-component)
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].y << " "; // center position (y-component)
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].z << " "; // center position (y-component)
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].vel_x << " "; // center velocity (x-component)
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].vel_y << " "; // center velocity (y-component)
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].vel_z << " "; // center velocity (y-component)
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].omgX << " "; // Rotational speed
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].omgY << " "; // Rotational speed
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].omgZ << " "; // Rotational speed
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tetaX << " "; // Angle with respect to Horizontal line
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tetaY << " "; // Angle with respect to Horizontal line
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tetaZ << " "; // Angle with respect to Horizontal line
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].previous_vel_x << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].previous_vel_y << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].previous_vel_z << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tempx[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tempx[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tempy[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tempy[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tempz[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].tempz[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_omgX[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_omgX[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_omgY[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_omgY[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_omgZ[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_omgZ[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Omg_prime_X << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Omg_prime_Y << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Omg_prime_Z << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Qu0 << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Qu1 << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Qu2 << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].Qu3 << " ";

            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temperature << " "; // Temperature
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_temperat[0] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_temperat[1] << " ";
            output_file << std::fixed << std::setprecision(12) << particle_list[npn].center[0].temp_part << " ";

            for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
            {
                output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].x << " ";
                output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].y << " ";
                output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].z << " ";
                //output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].vel_x << " ";
                //output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].vel_y << " ";
                //output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].vel_z << " ";
                //output_file << std::fixed << std::setprecision(15) << particle_list[npn].node[np1].temperature << " ";
                if (MOVING_SPHEROID	|| STATIONARY_SPHEROID)
                {
                    output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].x_pos_1 << " ";
                    output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].y_pos_1 << " ";
                    output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].z_pos_1 << " ";
                }
                output_file << std::fixed << std::setprecision(12) << particle_list[npn].node[np1].area_node<< " ";
            }
            /// Close particle recovery data file
            output_file.close();
        }
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: READ_PARTICLE_RECOVERY
//-----------------------------------------------------------------------------------------------------------------------------
//
void ParticleResolved::read_particle_recovery()
{

    if (num_particles > 0)
    {
        for (INT32 npn = 0; npn < num_particles; ++npn)
        {
            ostringstream oss;
            oss << npn + 1;
            string npstr = oss.str();//convert number to string
            string input_filename("FlowLB_Results/data_particle_recovery_");
            input_filename += npstr;
            input_filename += ".dat";
            ifstream input_file;

            /// Open file
            input_file.open(input_filename.c_str());

            /// Read particle data
            input_file >> particle_list[npn].particle_mass;
            input_file >> particle_list[npn].parvol;
            input_file >> particle_list[npn].mom_inertia_x;
            input_file >> particle_list[npn].mom_inertia_y;
            input_file >> particle_list[npn].mom_inertia_z;
            input_file >> particle_num_nodes;
            input_file >> particle_list[npn].particle_perimeter;

            /// Read data
            input_file >> particle_list[npn].center[0].x; // center position (x-component)
            input_file >> particle_list[npn].center[0].y; // center position (y-component)
            input_file >> particle_list[npn].center[0].z; // center position (y-component)
            input_file >> particle_list[npn].center[0].vel_x; // center velocity (x-component)
            input_file >> particle_list[npn].center[0].vel_y; // center velocity (y-component)
            input_file >> particle_list[npn].center[0].vel_z; // center velocity (y-component)
            input_file >> particle_list[npn].center[0].omgX; // Rotational speed
            input_file >> particle_list[npn].center[0].omgY; // Rotational speed
            input_file >> particle_list[npn].center[0].omgZ; // Rotational speed
            input_file >> particle_list[npn].center[0].tetaX; // Angle
            input_file >> particle_list[npn].center[0].tetaY; // Angle
            input_file >> particle_list[npn].center[0].tetaZ; // Angle
            input_file >> particle_list[npn].center[0].previous_vel_x;
            input_file >> particle_list[npn].center[0].previous_vel_y;
            input_file >> particle_list[npn].center[0].previous_vel_z;
            input_file >> particle_list[npn].center[0].tempx[0];
            input_file >> particle_list[npn].center[0].tempx[1];
            input_file >> particle_list[npn].center[0].tempy[0];
            input_file >> particle_list[npn].center[0].tempy[1];
            input_file >> particle_list[npn].center[0].tempz[0];
            input_file >> particle_list[npn].center[0].tempz[1];
            input_file >> particle_list[npn].center[0].temp_omgX[0];
            input_file >> particle_list[npn].center[0].temp_omgX[1];
            input_file >> particle_list[npn].center[0].temp_omgY[0];
            input_file >> particle_list[npn].center[0].temp_omgY[1];
            input_file >> particle_list[npn].center[0].temp_omgZ[0];
            input_file >> particle_list[npn].center[0].temp_omgZ[1];
            input_file >> particle_list[npn].center[0].Omg_prime_X;
            input_file >> particle_list[npn].center[0].Omg_prime_Y;
            input_file >> particle_list[npn].center[0].Omg_prime_Z;
            input_file >> particle_list[npn].center[0].Qu0;
            input_file >> particle_list[npn].center[0].Qu1;
            input_file >> particle_list[npn].center[0].Qu2;
            input_file >> particle_list[npn].center[0].Qu3;

            input_file >> particle_list[npn].center[0].temperature; // Angle with respect to Horizontal line
            input_file >> particle_list[npn].center[0].temp_temperat[0];
            input_file >> particle_list[npn].center[0].temp_temperat[1];
            input_file >> particle_list[npn].center[0].temp_part;

            for (INT32 np1 = 0; np1 < particle_num_nodes; ++np1)
            {
                input_file >> particle_list[npn].node[np1].x;
                input_file >> particle_list[npn].node[np1].y;
                input_file >> particle_list[npn].node[np1].z;
                //input_file >> particle_list[npn].node[np1].vel_x;
                //input_file >> particle_list[npn].node[np1].vel_y;
                //input_file >> particle_list[npn].node[np1].vel_z;
                //input_file >> particle_list[npn].node[np1].temperature;
                if (MOVING_SPHEROID	|| STATIONARY_SPHEROID)
                {
                    input_file >> particle_list[npn].node[np1].x_pos_1;
                    input_file >> particle_list[npn].node[np1].y_pos_1;
                    input_file >> particle_list[npn].node[np1].z_pos_1;
                }
                input_file >> particle_list[npn].node[np1].area_node;
            }
            /// Close particle recovery data file
            input_file.close();
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: RUNGEKUTTA
//-----------------------------------------------------------------------------------------------------------------------------
//this function solves 3+4 equations for "3:rotational speed in body-fixed coordinate" + "4:Quaternions:q0,q1,q2,q3"

void ParticleResolved::rungekutta(INT32 np)
{
//double& omega_prm_x, double& omega_prm_y, double& omega_prm_z,	double Tor_prm_x, double Tor_prm_y, double Tor_prm_z,
    //double mom_inr_x, double mom_inr_y, double mom_inr_z, double& Q0, double& Q1, double& Q2, double& Q3

    double dh, f, k1, k2, k3, k4;
    double g, r1, r2, r3, r4;
    double h, s1, s2, s3, s4;
    double m, b1, b2, b3, b4;
    double n, c1, c2, c3, c4;
    double p, j1, j2, j3, j4;
    double q, t1, t2, t3, t4;
    double omg_prm_x_1, omg_prm_y_1, omg_prm_z_1;
    double Q0_1, Q1_1, Q2_1, Q3_1;
    double Tx = particle_list[np].center[0].T_prime_X;
    double Ty = particle_list[np].center[0].T_prime_Y;
    double Tz = particle_list[np].center[0].T_prime_Z;
    double mom_x = particle_list[np].mom_inertia_x;
    double mom_y = particle_list[np].mom_inertia_y;
    double mom_z = particle_list[np].mom_inertia_z;

#define FF(omg_prm_y,omg_prm_z)  Tx/mom_x + omg_prm_y * omg_prm_z * (mom_y - mom_z)/mom_x
#define GG(omg_prm_x,omg_prm_z)  Ty/mom_y + omg_prm_x * omg_prm_z * (mom_z - mom_x)/mom_y
#define HH(omg_prm_x,omg_prm_y)  Tz/mom_z + omg_prm_x * omg_prm_y * (mom_x - mom_y)/mom_z
#define MM(Q1,Q2,Q3,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*(-Q1*omg_prm_x - Q2*omg_prm_y - Q3*omg_prm_z)
#define NN(Q0,Q2,Q3,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*( Q0*omg_prm_x - Q3*omg_prm_y + Q2*omg_prm_z)
#define PP(Q0,Q1,Q3,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*( Q3*omg_prm_x + Q0*omg_prm_y - Q1*omg_prm_z)
#define QQ(Q0,Q1,Q2,omg_prm_x,omg_prm_y,omg_prm_z)  0.5*(-Q2*omg_prm_x + Q1*omg_prm_y + Q0*omg_prm_z)

    double omg_prm_x_0 = particle_list[np].center[0].Omg_prime_X;//Initial value
    double omg_prm_y_0 = particle_list[np].center[0].Omg_prime_Y;//Initial value
    double omg_prm_z_0 = particle_list[np].center[0].Omg_prime_Z;//Initial value
    double Q0_0 = particle_list[np].center[0].Qu0;//Initial value
    double Q1_0 = particle_list[np].center[0].Qu1;//Initial value
    double Q2_0 = particle_list[np].center[0].Qu2;//Initial value
    double Q3_0 = particle_list[np].center[0].Qu3;//Initial value

    dh = 1;//dh here must be time step which is in lattice calculations = 1

    //double nn;
    //printf("\nEnter the value of last point: ");
    //scanf("%lf",&nn);
    //for(double t=0; t<nn; t=t+dh)
    //{

    f = FF(omg_prm_y_0, omg_prm_z_0);
    k1 = dh * f;
    g = GG(omg_prm_x_0, omg_prm_z_0);
    r1 = dh * g;
    h = HH(omg_prm_x_0, omg_prm_y_0);
    s1 = dh * h;
    m = MM(Q1_0, Q2_0, Q3_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
    b1 = dh * m;
    n = NN(Q0_0, Q2_0, Q3_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
    c1 = dh * n;
    p = PP(Q0_0, Q1_0, Q3_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
    j1 = dh * p;
    q = QQ(Q0_0, Q1_0, Q2_0, omg_prm_x_0, omg_prm_y_0, omg_prm_z_0);
    t1 = dh * q;

    f = FF((omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
    k2 = dh * f;
    g = GG((omg_prm_x_0 + k1 / 2), (omg_prm_z_0 + s1 / 2));
    r2 = dh * g;
    h = HH((omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2));
    s2 = dh * h;
    m = MM((Q1_0 + c1 / 2), (Q2_0 + j1 / 2), (Q3_0 + t1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
    b2 = dh * m;
    n = NN((Q0_0 + b1 / 2), (Q2_0 + j1 / 2), (Q3_0 + t1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
    c2 = dh * n;
    p = PP((Q0_0 + b1 / 2), (Q1_0 + c1 / 2), (Q3_0 + t1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
    j2 = dh * p;
    q = QQ((Q0_0 + b1 / 2), (Q1_0 + c1 / 2), (Q2_0 + j1 / 2), (omg_prm_x_0 + k1 / 2), (omg_prm_y_0 + r1 / 2), (omg_prm_z_0 + s1 / 2));
    t2 = dh * q;

    f = FF((omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
    k3 = dh * f;
    g = GG((omg_prm_x_0 + k2 / 2), (omg_prm_z_0 + s2 / 2));
    r3 = dh * g;
    h = HH((omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2));
    s3 = dh * h;
    m = MM((Q1_0 + c2 / 2), (Q2_0 + j2 / 2), (Q3_0 + t2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
    b3 = dh * m;
    n = NN((Q0_0 + b2 / 2), (Q2_0 + j2 / 2), (Q3_0 + t2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
    c3 = dh * n;
    p = PP((Q0_0 + b2 / 2), (Q1_0 + c2 / 2), (Q3_0 + t2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
    j3 = dh * p;
    q = QQ((Q0_0 + b2 / 2), (Q1_0 + c2 / 2), (Q2_0 + j2 / 2), (omg_prm_x_0 + k2 / 2), (omg_prm_y_0 + r2 / 2), (omg_prm_z_0 + s2 / 2));
    t3 = dh * q;

    f = FF((omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
    k4 = dh * f;
    g = GG((omg_prm_x_0 + k3), (omg_prm_z_0 + s3));
    r4 = dh * g;
    h = HH((omg_prm_x_0 + k3), (omg_prm_y_0 + r3));
    s4 = dh * h;
    m = MM((Q1_0 + c3), (Q2_0 + j3), (Q3_0 + t3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
    b4 = dh * m;
    n = NN((Q0_0 + b3), (Q2_0 + j3), (Q3_0 + t3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
    c4 = dh * n;
    p = PP((Q0_0 + b3), (Q1_0 + c3), (Q3_0 + t3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
    j4 = dh * p;
    q = QQ((Q0_0 + b3), (Q1_0 + c3), (Q2_0 + j3), (omg_prm_x_0 + k3), (omg_prm_y_0 + r3), (omg_prm_z_0 + s3));
    t4 = dh * q;

    omg_prm_x_1 = omg_prm_x_0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
    omg_prm_y_1 = omg_prm_y_0 + (r1 + 2 * r2 + 2 * r3 + r4) / 6.;
    omg_prm_z_1 = omg_prm_z_0 + (s1 + 2 * s2 + 2 * s3 + s4) / 6.;
    Q0_1 = Q0_0 + (b1 + 2 * b2 + 2 * b3 + b4) / 6.;
    Q1_1 = Q1_0 + (c1 + 2 * c2 + 2 * c3 + c4) / 6.;
    Q2_1 = Q2_0 + (j1 + 2 * j2 + 2 * j3 + j4) / 6.;
    Q3_1 = Q3_0 + (t1 + 2 * t2 + 2 * t3 + t4) / 6.;

    omg_prm_x_0 = omg_prm_x_1;
    omg_prm_y_0 = omg_prm_y_1;
    omg_prm_z_0 = omg_prm_z_1;
    Q0_0 = Q0_1;
    Q1_0 = Q1_1;
    Q2_0 = Q2_1;
    Q3_0 = Q3_1;
    // }
    particle_list[np].center[0].Omg_prime_X = omg_prm_x_1;
    particle_list[np].center[0].Omg_prime_Y = omg_prm_y_1;
    particle_list[np].center[0].Omg_prime_Z = omg_prm_z_1;
    particle_list[np].center[0].Qu0 = Q0_1;
    particle_list[np].center[0].Qu1 = Q1_1;
    particle_list[np].center[0].Qu2 = Q2_1;
    particle_list[np].center[0].Qu3 = Q3_1;
    return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: QUATERNION
//-----------------------------------------------------------------------------------------------------------------------------
//
// This function calculates the elements of quaternion matrix.(Eq.2.13 Article 876: H.Huang et al. J. Fluid Mech. (2012), vol. 692, pp. 369-394.
void ParticleResolved::quaternion(INT32 np)
{
    double Q0 = particle_list[np].center[0].Qu0;
    double Q1 = particle_list[np].center[0].Qu1;
    double Q2 = particle_list[np].center[0].Qu2;
    double Q3 = particle_list[np].center[0].Qu3;

    double quat[10];
    quat[0] = 0.;//this value is useless
    quat[1] = SQ(Q0) + SQ(Q1) - SQ(Q2) - SQ(Q3);
    quat[2] = 2 * (Q1 * Q2 + Q0 * Q3);
    quat[3] = 2 * (Q1 * Q3 - Q0 * Q2);
    quat[4] = 2 * (Q1 * Q2 - Q0 * Q3);
    quat[5] = SQ(Q0) - SQ(Q1) + SQ(Q2) - SQ(Q3);
    quat[6] = 2 * (Q2 * Q3 + Q0 * Q1);
    quat[7] = 2 * (Q1 * Q3 + Q0 * Q2);
    quat[8] = 2 * (Q2 * Q3 - Q0 * Q1);
    quat[9] = SQ(Q0) - SQ(Q1) - SQ(Q2) + SQ(Q3);

    for (INT32 i = 0; i < 10; i++)
    {
        particle_list[np].center[0].quater[i] = quat[i];
    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: DETERMINANT
//-----------------------------------------------------------------------------------------------------------------------------
//This function calculates the determinant of quaternion matrix (in order to obtain its inverse)

double ParticleResolved::determinant(INT32 np)
{
    double quat[10];
    for (INT32 i = 0; i < 10; i++)
    {
        quat[i] = particle_list[np].center[0].quater[i];
    }

    double det;
    det = quat[1] * (quat[5] * quat[9] - quat[6] * quat[8])
          - quat[2] * (quat[4] * quat[9] - quat[6] * quat[7])
          + quat[3] * (quat[4] * quat[8] - quat[5] * quat[7]);
    return det;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COMP_ROT_SPEED
//-----------------------------------------------------------------------------------------------------------------------------
// (This function transforms the rotational speed from body-fixed coordinate to Inertial coordinate
// This function is just for SPHEROID
void ParticleResolved::comp_rot_speed(INT32 np)
{
//void ParticleResolved::comp_rot_speed(double determi, double& omgX, double& omgY, double& omgZ, INT32 n, double quat[], double omg_pr_x, double omg_pr_y, double omg_pr_z) {

    double RNI[10];
    double om_temp_x = particle_list[np].center[0].Omg_prime_X;
    double om_temp_y = particle_list[np].center[0].Omg_prime_Y;
    double om_temp_z = particle_list[np].center[0].Omg_prime_Z;
    double determi   = particle_list[np].center[0].deter;
    double quat[10];
    for (INT32 i = 0; i < 10; i++)
    {
        quat[i] = particle_list[np].center[0].quater[i];
    }

    RNI[0] = 0.;
    RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
    RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
    RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
    RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
    RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
    RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
    RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
    RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
    RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
    RNI[1] = RNI[1] / determi;
    RNI[2] = RNI[2] / determi;
    RNI[3] = RNI[3] / determi;
    RNI[4] = RNI[4] / determi;
    RNI[5] = RNI[5] / determi;
    RNI[6] = RNI[6] / determi;
    RNI[7] = RNI[7] / determi;
    RNI[8] = RNI[8] / determi;
    RNI[9] = RNI[9] / determi;

    particle_list[np].center[0].omgX = RNI[1] * om_temp_x + RNI[2] * om_temp_y + RNI[3] * om_temp_z;
    particle_list[np].center[0].omgY = RNI[4] * om_temp_x + RNI[5] * om_temp_y + RNI[6] * om_temp_z;
    particle_list[np].center[0].omgZ = RNI[7] * om_temp_x + RNI[8] * om_temp_y + RNI[9] * om_temp_z;

    return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: UPDATE_ROT_SPHEROID
//-----------------------------------------------------------------------------------------------------------------------------
// If at t=0 Euler angles are not zero then the spheroid must be rotated. This function is part of caculating the Lagrangian node position in INERTIAL coordinate based on Euler angles and node location in BODY coordinate
// This function is just for SPHEROID
void ParticleResolved::update_rot_spheroid(INT32 np)
{
//void update_rot_spheroid(double determi, double& X_inertial, double& Y_inertial, double& Z_inertial, double quat[], double X_body, double Y_body, double Z_body, double center_x, double center_y, double center_z) {
    // X_body is node X position in body coordinate minus particle center
    // X_inertial is node X position in INERTIAL coordinate

    double quat[10];
    for (INT32 i = 0; i < 10; i++)
    {
        quat[i] = particle_list[np].center[0].quater[i];
    }
    double determi = particle_list[np].center[0].deter;
    double X_body{};
    double Y_body{};
    double Z_body{};

    double RNI[10];

    RNI[0] = 0.;
    RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
    RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
    RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
    RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
    RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
    RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
    RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
    RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
    RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
    RNI[1] = RNI[1] / determi;
    RNI[2] = RNI[2] / determi;
    RNI[3] = RNI[3] / determi;
    RNI[4] = RNI[4] / determi;
    RNI[5] = RNI[5] / determi;
    RNI[6] = RNI[6] / determi;
    RNI[7] = RNI[7] / determi;
    RNI[8] = RNI[8] / determi;
    RNI[9] = RNI[9] / determi;

    for (INT32 n = 0; n < particle_num_nodes; ++n)
    {
        X_body = particle_list[np].node[n].x_pos_1;
        Y_body = particle_list[np].node[n].y_pos_1;
        Z_body = particle_list[np].node[n].z_pos_1;
        particle_list[np].node[n].x = (RNI[1] * X_body + RNI[2] * Y_body + RNI[3] * Z_body) + particle_list[np].center[0].x;
        particle_list[np].node[n].y = (RNI[4] * X_body + RNI[5] * Y_body + RNI[6] * Z_body) + particle_list[np].center[0].y;
        particle_list[np].node[n].z = (RNI[7] * X_body + RNI[8] * Y_body + RNI[9] * Z_body) + particle_list[np].center[0].z;
    }

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: UPDATE_ROT_SPHEROID_COLLISION
//-----------------------------------------------------------------------------------------------------------------------------
// If at t=0 Euler angles are not zero then the spheroid must be rotated. This function is part of caculating the Lagrangian node position in INERTIAL coordinate based on Euler angles and node location in BODY coordinate
// This function is just for SPHEROID
void ParticleResolved::update_rot_spheroid_collision(double determi, double & X_inertial, double & Y_inertial, double & Z_inertial, double quat[], double X_body, double Y_body, double Z_body, double center_x, double center_y, double center_z)
{
    // X_body is node X position in body coordinate minus particle center
    // X_inertial is node X position in INERTIAL coordinate

    double RNI[10];

    RNI[0] = 0.;
    RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
    RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
    RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
    RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
    RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
    RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
    RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
    RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
    RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
    RNI[1] = RNI[1] / determi;
    RNI[2] = RNI[2] / determi;
    RNI[3] = RNI[3] / determi;
    RNI[4] = RNI[4] / determi;
    RNI[5] = RNI[5] / determi;
    RNI[6] = RNI[6] / determi;
    RNI[7] = RNI[7] / determi;
    RNI[8] = RNI[8] / determi;
    RNI[9] = RNI[9] / determi;

    X_inertial = (RNI[1] * X_body + RNI[2] * Y_body + RNI[3] * Z_body) + center_x;
    Y_inertial = (RNI[4] * X_body + RNI[5] * Y_body + RNI[6] * Z_body) + center_y;
    Z_inertial = (RNI[7] * X_body + RNI[8] * Y_body + RNI[9] * Z_body) + center_z;

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: NORMAL_VEC_SPHEROID
//-----------------------------------------------------------------------------------------------------------------------------
// Calculates NORMAL vector at each node at INERTIAL coordinate based on the body coordinate valuse
void ParticleResolved::normal_vec_spheroid(double determi, double& X_inertial, double& Y_inertial, double& Z_inertial, INT32 n, double quat[], double X_body, double Y_body, double Z_body)
{
    // X_body is normal vector X index in body coordinate
    // X_inertial normal vector X index in INERTIAL coordinate

    double RNI[10];

    RNI[0] = 0.;
    RNI[1] = quat[5] * quat[9] - quat[6] * quat[8];
    RNI[2] = -quat[2] * quat[9] + quat[3] * quat[8];
    RNI[3] = quat[2] * quat[6] - quat[3] * quat[5];
    RNI[4] = -quat[4] * quat[9] + quat[6] * quat[7];
    RNI[5] = quat[1] * quat[9] - quat[3] * quat[7];
    RNI[6] = -quat[1] * quat[6] + quat[3] * quat[4];
    RNI[7] = quat[4] * quat[8] - quat[5] * quat[7];
    RNI[8] = -quat[1] * quat[8] + quat[2] * quat[7];
    RNI[9] = quat[1] * quat[5] - quat[2] * quat[4];
    RNI[1] = RNI[1] / determi;
    RNI[2] = RNI[2] / determi;
    RNI[3] = RNI[3] / determi;
    RNI[4] = RNI[4] / determi;
    RNI[5] = RNI[5] / determi;
    RNI[6] = RNI[6] / determi;
    RNI[7] = RNI[7] / determi;
    RNI[8] = RNI[8] / determi;
    RNI[9] = RNI[9] / determi;

    X_inertial = (RNI[1] * X_body + RNI[2] * Y_body + RNI[3] * Z_body);
    Y_inertial = (RNI[4] * X_body + RNI[5] * Y_body + RNI[6] * Z_body);
    Z_inertial = (RNI[7] * X_body + RNI[8] * Y_body + RNI[9] * Z_body);

    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: CHECK_ZERO
//-----------------------------------------------------------------------------------------------------------------------------
// This function checks if a variable is less than epsilon it is set to zero
void ParticleResolved::check_zero(double& variab, double epsilon)
{
    if (abs(variab) < epsilon)
    {
        variab = 0.0;
    }
    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISION_MODEL
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::collision_model(INT32 np){
    if (MOVING_SPHERE)
    {
        collision_sphere(np);
    }
    if (MOVING_SPHEROID)
    {
        collision_spheroid(np);
    }
    if (MOVING_CYLINDER)
    {
        collision_circle(np);
    }
    if (MOVING_ELLIPSE)
    {
        collision_ellipse(np);
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_CIR
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for circular objects
void ParticleResolved::collision_circle(INT32 np) {
    if (Circle_Collision_1)
    {
        ///Particle-particle Collision Forces. JUST for Cylinder shape
        double Dij;
        for (INT32 np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
            if (np2 != np) {
                Dij = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dij <= (2.0 * particle_radius + particle_th)) {
                    particle_list[np].center[0].f3x += particle_cij / particle_ep * ((particle_list[np].center[0].x - particle_list[np2].center[0].x)) / Dij * pow((Dij - 2.0 * particle_radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3y += particle_cij / particle_ep * ((particle_list[np].center[0].y - particle_list[np2].center[0].y)) / Dij * pow((Dij - 2.0 * particle_radius - particle_th) / particle_th, 2);
                }
            }
        }// END of Particle-Particle Collision

        /// Particle-Wall collision
        double Rij;
        if (Y_Wall)
        {
            if (particle_list[np].center[0].y <= (Ny - 2) / 2){//Particle-Y-Wall Collision section
                Rij = particle_list[np].center[0].y + particle_radius;
            }
            else
            {
                Rij = -1.0 * ((double)(Ny - 1) - particle_list[np].center[0].y + particle_radius);
            }
            if (abs(Rij) <= (2.0 * particle_radius + particle_thw)){
                particle_list[np].center[0].f4x += 0.0;
                particle_list[np].center[0].f4y += particle_cij / particle_epw * pow((2 * particle_radius - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
            }
        }// END of Particle-Y-Wall Collision

        if (X_Wall)
        {
            if (particle_list[np].center[0].x <= (Nx - 2) / 2){//Particle-X-Wall Collision section
                Rij = particle_list[np].center[0].x + particle_radius;
            }
            else
            {
                Rij = -1.0 * ((double)(Nx - 1) - particle_list[np].center[0].x + particle_radius);
            }
            if (abs(Rij) <= (2.0 * particle_radius + particle_thw)){
                particle_list[np].center[0].f4x += particle_cij / particle_epw * pow((2.0 * particle_radius - abs(Rij) + particle_thw) / particle_thw, 2)*Rij / abs(Rij);
                particle_list[np].center[0].f4y += 0;
            }
        }	 // END of Particle-X-Wall Collision
    }// End of Circle_Collision_1
///********************************************************************************//
    if (Circle_Collision_2)
    {
        ///Particle-particle Collision Forces. JUST for Cylinder shape
        double Dij;
        for (INT32 np2 = 0; np2 < num_particles; ++np2)
        {
            if (np2 != np) {
                Dij = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dij <= (2.0*particle_radius)){
                    particle_list[np].center[0].f3x += 1.0 / particle_ep * ((particle_list[np].center[0].x - particle_list[np2].center[0].x)) * (2.0 * particle_radius - Dij);//Total particle-particle collision force
                    particle_list[np].center[0].f3y += 1.0 / particle_ep * ((particle_list[np].center[0].y - particle_list[np2].center[0].y)) * (2.0 * particle_radius - Dij);
                }
                if (Dij >(2 * particle_radius) && Dij <= (2.0*particle_radius + particle_th)){
                    particle_list[np].center[0].f3x += 1.0 / particle_ep * ((particle_list[np].center[0].x - particle_list[np2].center[0].x)) * pow((2.0 * particle_radius - Dij + particle_th), 2);
                    particle_list[np].center[0].f3y += 1.0 / particle_ep * ((particle_list[np].center[0].y - particle_list[np2].center[0].y)) * pow((2.0 * particle_radius - Dij + particle_th), 2);
                }
            }
        } // END of Particle-Particle Collision

        /// Particle-Wall collision
        double Rij;
        if (Y_Wall)
        {
            if (particle_list[np].center[0].y < (Ny - 2) / 2){//Particle-Wall Collision section
                Rij = particle_list[np].center[0].y + particle_radius;
            }
            else
            {
                Rij = -1.0 * ((Ny - 1) - particle_list[np].center[0].y + particle_radius);
            }
            if (abs(Rij) < 2.0 * particle_radius) {
                particle_list[np].center[0].f4x += 1.0 / particle_epw * 0.;
                particle_list[np].center[0].f4y += 1.0 / particle_epw * Rij * (2 * particle_radius - abs(Rij));
            }
            if (abs(Rij) >= 2 * particle_radius && abs(Rij) <= (2.0 * particle_radius + particle_thw)){
                particle_list[np].center[0].f4x += 1.0 / particle_epw * 0;
                particle_list[np].center[0].f4y += 1.0 / particle_epw * Rij * pow((2 * particle_radius - abs(Rij) + particle_thw), 2);
            }
        }// END of Particle-Y-Wall Collision

        if (X_Wall)
        {
            if (particle_list[np].center[0].x < (Nx - 2) / 2){//Particle-Y-Wall Collision section
                Rij = particle_list[np].center[0].x + particle_radius;
            }
            else
            {
                Rij = -1 * ((Nx - 1) - particle_list[np].center[0].x + particle_radius);
            }
            if (abs(Rij) < 2 * particle_radius) {
                particle_list[np].center[0].f4x += 1 / particle_epw * Rij * (2.0 * particle_radius - abs(Rij));
                particle_list[np].center[0].f4y += 1 / particle_epw * 0;
            }
            if (abs(Rij) >= 2 * particle_radius && abs(Rij) <= (2 * particle_radius + particle_thw)){
                particle_list[np].center[0].f4x += 1 / particle_epw * Rij * pow((2.0 * particle_radius - abs(Rij) + particle_thw), 2);
                particle_list[np].center[0].f4y += 1 / particle_epw * 0.;
            }
        }// END of Particle-X-Wall Collision
    } //end of Circle_Collision_2
///********************************************************************************//
    if (Circle_Collision_3)
    {
        ///Particle-particle Collision Forces. JUST for Cylinder shape
        double Dij;
        for (INT32 np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
            if (np2 != np) {
                Dij = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dij <= (2.0 * particle_radius + particle_th) && Dij>(2.0 * particle_radius)) {
                    particle_list[np].center[0].f3x += particle_cij / particle_ep * ((particle_list[np].center[0].x - particle_list[np2].center[0].x)) / Dij * pow((Dij - 2.0 * particle_radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3y += particle_cij / particle_ep * ((particle_list[np].center[0].y - particle_list[np2].center[0].y)) / Dij * pow((Dij - 2.0 * particle_radius - particle_th) / particle_th, 2);
                }
                if (Dij <= (2.0 * particle_radius)) {
                    particle_list[np].center[0].f3x += (particle_list[np].center[0].x - particle_list[np2].center[0].x) / Dij * (particle_cij / particle_ep *pow((Dij - 2.0 * particle_radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2 * particle_radius - Dij));
                    particle_list[np].center[0].f3y += (particle_list[np].center[0].y - particle_list[np2].center[0].y) / Dij * (particle_cij / particle_ep *pow((Dij - 2.0 * particle_radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2 * particle_radius - Dij));
                }
            }
        }// END of Particle-Particle Collision

        /// Particle-Wall collision
        double Rij;
        if (Y_Wall)
        {
            if (particle_list[np].center[0].y <= (Ny - 2) / 2){//Particle-Y-Wall Collision section
                Rij = particle_list[np].center[0].y + particle_radius;
            }
            else
            {
                Rij = -1.0 * ((Ny - 1) - particle_list[np].center[0].y + particle_radius);
            }

            if (abs(Rij) <= (2.0 * particle_radius + particle_thw) && abs(Rij) > (2.0 * particle_radius)){
                particle_list[np].center[0].f4x += 0.0;
                particle_list[np].center[0].f4y += particle_cij / particle_epw * pow((2.0 * particle_radius - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
            }
            if (abs(Rij) <= (2 * particle_radius)){
                particle_list[np].center[0].f4x += 0.0;
                particle_list[np].center[0].f4y += (particle_cij / particle_epw *pow((abs(Rij) - 2.0 * particle_radius - particle_thw) / particle_thw, 2.0) + particle_cij / particle_Epw / particle_thw * (2.0*particle_radius - abs(Rij)))*Rij / abs(Rij);
            }
        }// END of Particle-Y-Wall Collision

        if (X_Wall)
        {
            if (particle_list[np].center[0].x <= (Nx - 2) / 2){//Particle-Wall Collision section
                Rij = particle_list[np].center[0].x + particle_radius;
            }
            else
            {
                Rij = -1.0 * ((Nx - 1) - particle_list[np].center[0].x + particle_radius);
            }

            if (abs(Rij) <= (2.0 * particle_radius + particle_thw) && abs(Rij) > (2.0 * particle_radius)){
                particle_list[np].center[0].f4x += particle_cij / particle_epw * pow((2.0 * particle_radius - abs(Rij) + particle_thw) / particle_thw, 2)*Rij / abs(Rij);
                particle_list[np].center[0].f4y += 0.0;
            }
            if (abs(Rij) <= (2.0 * particle_radius)){
                particle_list[np].center[0].f4x += (particle_cij / particle_epw *pow((abs(Rij) - 2.0 * particle_radius - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0*particle_radius - abs(Rij)))*Rij / abs(Rij);
                particle_list[np].center[0].f4y += 0.0;
            }
        } // END of Particle-X-Wall Collision
    }//end of Circle_Collision_3

	return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_ELLIPSE
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for Elliptical objects
void ParticleResolved::collision_ellipse(INT32 np) {

    if (Ellipse_Collision_1)
    {
        //Collision Forces. JUST for ELLIPSE shape
        double Dij, Dist, ffx, ffy, ffx4, ffy4, angle_1, angle_2, fas, xnor, ynor, norvec, xnor_init, ynor_init;
        ffx = 0.0;
        ffy = 0.0;
        ///Particle-particle Collision Forces.
        for (INT32 np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
            if (np2 != np) {
                Dist = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dist <= (2 * std::max(particle_radius, particle_radius_2) + particle_th + 2)){
                    for (INT32 n = 0; n < particle_num_nodes; ++n) {
                        angle_1 = atan((particle_list[np].node[0].y - particle_list[np].center[0].y) / (particle_list[np].node[0].x - particle_list[np].center[0].x));//It gives us the angle of particle np
                        xnor = (2 * (particle_list[np].node[n].x - particle_list[np].center[0].x)*cos(angle_1) - 2 * (particle_list[np].node[n].y - particle_list[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
                        ynor = (2 * (particle_list[np].node[n].x - particle_list[np].center[0].x)*sin(angle_1) + 2 * (particle_list[np].node[n].y - particle_list[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
                        xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
                        ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
                        norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
                        for (INT32 m = 0; m < particle_num_nodes; ++m) {
                            Dij = sqrt(SQ(particle_list[np].node[n].x - particle_list[np2].node[m].x) + SQ(particle_list[np].node[n].y - particle_list[np2].node[m].y));
                            angle_2 = atan((particle_list[np2].node[0].y - particle_list[np2].center[0].y) / (particle_list[np2].node[0].x - particle_list[np2].center[0].x));//It gives us the angle of particle np2
                            fas = SQ((particle_list[np].node[n].x - particle_list[np2].center[0].x)*cos(angle_2) + (particle_list[np].node[n].y - particle_list[np2].center[0].y)*sin(angle_2)) / SQ(particle_radius) + SQ((particle_list[np].node[n].x - particle_list[np2].center[0].x)*sin(angle_2) - (particle_list[np].node[n].y - particle_list[np2].center[0].y)*cos(angle_2)) / SQ(particle_radius_2);

                            if (Dij <= particle_th  && fas > 1) { //It means that particles have NOT penetrated into each other
                                // for circle// alwaya keep deactive// ffx = -1. * cij / ep * ((particle_list[np].node[n].x-particle_list[np].center[0].x))/sqrt(SQ(particle_list[np].node[n].x-particle_list[np].center[0].x)+SQ(particle_list[np].node[n].y-particle_list[np].center[0].y))  * pow((Dij - th) / th , 2);
                                // for circle// alwaya keep deactive// ffy = -1. * cij / ep * ((particle_list[np].node[n].y-particle_list[np].center[0].y))/sqrt(SQ(particle_list[np].node[n].x-particle_list[np].center[0].x)+SQ(particle_list[np].node[n].y-particle_list[np].center[0].y))  * pow((Dij - th) / th , 2);
                                ffx = -1. * particle_cij / particle_ep * xnor_init / norvec * pow((Dij - particle_th) / particle_th, 2);
                                ffy = -1. * particle_cij / particle_ep * ynor_init / norvec * pow((Dij - particle_th) / particle_th, 2);
                                particle_list[np].center[0].f3x += ffx;
                                particle_list[np].center[0].f3y += ffy;
                                particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[n].y - particle_list[np].center[0].y) - ffy * (particle_list[np].node[n].x - particle_list[np].center[0].x);
                            }
                            /*if (Dij <= th && fas <= 1) {//It means that node n of particle np has penetrated to particle np2
                            ffx = -1. * cij / Ep * xnor_init / norvec * pow((Dij - th) / th , 2);
                            ffy = -1. * cij / Ep * ynor_init / norvec * pow((Dij - th) / th , 2);
                            particle_list[np].center[0].f3x += ffx;
                            particle_list[np].center[0].f3y += ffy;
                            particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[n].y-particle_list[np].center[0].y) - ffy * (particle_list[np].node[n].x-particle_list[np].center[0].x);
                            }*/

                        }
                    }
                }
            }
        }// END of Particle-Particle Collision

        double Rij;
        if (Y_Wall)
        {
            if ((particle_list[np].center[0].y <= std::max(particle_radius, particle_radius_2) + particle_thw + 2) || ((Ny - particle_list[np].center[0].y) <= std::max(particle_radius, particle_radius_2) + particle_thw + 2)){
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    angle_1 = atan((particle_list[np].node[0].y - particle_list[np].center[0].y) / (particle_list[np].node[0].x - particle_list[np].center[0].x));//It gives us the angle of particle np
                    xnor = (2 * (particle_list[np].node[n].x - particle_list[np].center[0].x)*cos(angle_1) - 2 * (particle_list[np].node[n].y - particle_list[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
                    ynor = (2 * (particle_list[np].node[n].x - particle_list[np].center[0].x)*sin(angle_1) + 2 * (particle_list[np].node[n].y - particle_list[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
                    xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
                    ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
                    norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
                    if (particle_list[np].node[n].y <= (Ny - 1) / 2){//Particle-Y-Wall Collision section
                        Rij = particle_list[np].node[n].y - 0;// 0 means y-location of lower wall
                    }
                    else
                    {
                        Rij = -1.0 * ((Ny - 1) - particle_list[np].node[n].y);// Ny-1 means y-location of upper wall
                    }
                    if (abs(Rij) <= (particle_thw)){
                        ffx4 = 0.0;
                        ffy4 = -1. * particle_cij / particle_epw * ynor_init / norvec * pow((-abs(Rij) + particle_thw) / particle_thw, 2);
                        particle_list[np].center[0].f4x += ffx4;
                        particle_list[np].center[0].f4y += ffy4;
                        particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[n].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[n].x - particle_list[np].center[0].x);
                    }
                }
            }
        }

        if (X_Wall)
        {
            if ((particle_list[np].center[0].x <= std::max(particle_radius, particle_radius_2) + particle_thw + 2) || ((Nx - particle_list[np].center[0].x) <= std::max(particle_radius, particle_radius_2) + particle_thw + 2)){
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    angle_1 = atan((particle_list[np].node[0].y - particle_list[np].center[0].y) / (particle_list[np].node[0].x - particle_list[np].center[0].x));//It gives us the angle of particle np
                    xnor = (2 * (particle_list[np].node[n].x - particle_list[np].center[0].x)*cos(angle_1) - 2 * (particle_list[np].node[n].y - particle_list[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
                    ynor = (2 * (particle_list[np].node[n].x - particle_list[np].center[0].x)*sin(angle_1) + 2 * (particle_list[np].node[n].y - particle_list[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
                    xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
                    ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
                    norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
                    if (particle_list[np].node[n].x <= (Nx - 1) / 2.){//Particle-X-Wall Collision section
                        Rij = particle_list[np].node[n].x - 0;// 0 means x-location of left wall
                    }
                    else
                    {
                        Rij = -1.0 * ((Nx - 1) - particle_list[np].node[n].x);// Nx-1 means x-location of right wall
                    }

                    if (abs(Rij) <= (particle_thw)){
                        // 	cout<<"angle "<<angle_1<<"xnor "<<xnor_init<<" norvec"<<norvec<<endl;
                        ffx4 = -1. * particle_cij / particle_epw * xnor_init / norvec * pow((-abs(Rij) + particle_thw) / particle_thw, 2);
                        ffy4 = 0;
                        particle_list[np].center[0].f4x += ffx4;
                        particle_list[np].center[0].f4y += ffy4;
                        particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[n].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[n].x - particle_list[np].center[0].x);
                    }
                }
            }
        }	 // END of Particle-Wall Collision
    } //end of Ellipse_Collision_1
///********************************************************************************//
    if (Ellipse_Collision_2)
    {
        //Collision Forces. JUST for ELLIPSE shape
        const double min_th = 0.0;// If distance is less than this value, Ep is used instead of ep
        const double min_thw = 0.0;// If distance is less than this value, Epw is used instead of epw
        double Dij, Dist, ffx, ffy, ffx4, ffy4, angle_1, angle_2, fas, xnor, ynor, norvec, xnor_init, ynor_init;
        ffx = 0.0;
        ffy = 0.0;
        INT32 node_1;
        INT32 node_2;
        double Dij_min;//Just to initially have a value greater than th

        for (INT32 np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
            if (np2 != np) {
                Dist = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dist <= (2 * std::max(particle_radius, particle_radius_2) + particle_th)){
                    node_1 = 0;
                    node_2 = 0;
                    Dij_min = particle_th + 1;//Just to initially have a value greater than th
                    for (INT32 n = 0; n < particle_num_nodes; ++n) {
                        for (INT32 m = 0; m < particle_num_nodes; ++m) {
                            Dij = sqrt(SQ(particle_list[np].node[n].x - particle_list[np2].node[m].x) + SQ(particle_list[np].node[n].y - particle_list[np2].node[m].y));
                            if (Dij <= Dij_min) {
                                Dij_min = Dij;
                                node_1 = n;
                                node_2 = m;
                            }
                        }
                    }
                    if (Dij_min >min_th && Dij_min <= particle_th){
                        angle_1 = atan((particle_list[np].node[0].y - particle_list[np].center[0].y) / (particle_list[np].node[0].x - particle_list[np].center[0].x));//It gives us the angle of particle np
                        xnor = (2 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x)*cos(angle_1) - 2 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
                        ynor = (2 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x)*sin(angle_1) + 2 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
                        xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
                        ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
                        norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
                        ffx = -1. * particle_cij / particle_ep * xnor_init / norvec * pow((Dij_min - particle_th) / particle_th, 2);
                        ffy = -1. * particle_cij / particle_ep * ynor_init / norvec * pow((Dij_min - particle_th) / particle_th, 2);
                        particle_list[np].center[0].f3x += ffx;
                        particle_list[np].center[0].f3y += ffy;
                        particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                    }
                    if (Dij_min <= min_th){
                        angle_1 = atan((particle_list[np].node[0].y - particle_list[np].center[0].y) / (particle_list[np].node[0].x - particle_list[np].center[0].x));//It gives us the angle of particle np
                        xnor = (2 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x)*cos(angle_1) - 2 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y)*sin(angle_1)) / SQ(particle_radius);
                        ynor = (2 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x)*sin(angle_1) + 2 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y)*cos(angle_1)) / SQ(particle_radius_2);
                        xnor_init = xnor*cos(angle_1) + ynor*sin(angle_1);
                        ynor_init = -1 * xnor*sin(angle_1) + ynor*cos(angle_1);
                        norvec = sqrt(SQ(xnor_init) + SQ(ynor_init));//  It is LENGTH of normal vector at each surface point of ellipse
                        ffx = -1. * particle_cij / particle_Ep * xnor_init / norvec * pow((Dij_min - particle_th) / particle_th, 2);
                        ffy = -1. * particle_cij / particle_Ep * ynor_init / norvec * pow((Dij_min - particle_th) / particle_th, 2);
                        particle_list[np].center[0].f3x += ffx;
                        particle_list[np].center[0].f3y += ffy;
                        particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                    }
                }
            }
        }// END of Particle-Particle Collision

        double Rij, Rij_min;;
        if (Y_Wall)
        {
            if ((particle_list[np].center[0].y <= std::max(particle_radius, particle_radius_2) + particle_thw) || ((Ny - particle_list[np].center[0].y) <= std::max(particle_radius, particle_radius_2) + particle_thw))
            {
                node_1 = 0;
                Rij_min = particle_thw + 1;//Just to initially have a value greater than th
                for (INT32 n = 0; n < particle_num_nodes; ++n)
                {
                    if (particle_list[np].node[n].y <= (Ny - 1) / 2.)
                    {//Particle-Y-Wall Collision section
                        Rij = abs(particle_list[np].node[n].y - 0);// 0 means y-location of lower wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    else
                    {
                        Rij = -1.0 * abs((Ny - 1) - particle_list[np].node[n].y);// Ny-1 means y-location of upper wall	.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    if (abs(Rij) <= abs(Rij_min)){
                        Rij_min = Rij;
                        node_1 = n;
                    }
                }
                if (abs(Rij_min) > min_thw && abs(Rij_min) <= particle_thw)
                {
                    ffx4 = 0.0;
                    ffy4 = 1. * particle_cij / particle_epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
                if (abs(Rij_min) <= min_thw)
                {
                    ffx4 = 0.0;
                    ffy4 = 1. * particle_cij / particle_Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
            }
        }// END of Particle_Y-Wall Collision

        if (X_Wall)
        {
            if ((particle_list[np].center[0].x <= std::max(particle_radius, particle_radius_2) + particle_thw) || ((Nx - particle_list[np].center[0].x) <= std::max(particle_radius, particle_radius_2) + particle_thw))
            {
                node_1 = 0;
                Rij_min = particle_thw + 1;//Just to initially have a value greater than th
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    if (particle_list[np].node[n].x <= (Nx - 1) / 2.){//Particle-Y-Wall Collision section
                        Rij = abs(particle_list[np].node[n].x - 0);// 0 means y-location of lower wall.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    else
                    {
                        Rij = -1.0 * abs((Nx - 1) - particle_list[np].node[n].x);// Ny-1 means y-location of upper wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    if (abs(Rij) <= abs(Rij_min)){
                        Rij_min = Rij;
                        node_1 = n;
                    }
                }
                if (abs(Rij_min) > min_thw && abs(Rij_min) <= particle_thw){
                    ffx4 = 1. * particle_cij / particle_epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    ffy4 = 0.0;
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
                if (abs(Rij_min) <= min_thw){
                    ffx4 = 1. * particle_cij / particle_Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    ffy4 = 0.0;
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
            }
        }// END of Particle_X-Wall Collision
    }//end of Ellipse_Collision_2
///********************************************************************************//
    if (Ellipse_Collision_3)
    {
        //Collision Forces. JUST for ELLIPSE shape
        const double min_th = 0.0;// If distance is less than this value, Ep is used instead of ep
        const double min_thw = 0.0;// If distance is less than this value, Epw is used instead of epw
        double Dij, Dist, ffx, ffy, ffx4, ffy4, angle_2, fas;
        ffx = 0.0;
        ffy = 0.0;
        INT32 node_1;
        INT32 node_2;
        double Dij_min;//Just to initially have a value greater than th

        for (INT32 np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
            if (np2 != np) {
                Dist = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dist <= (2 * std::max(particle_radius, particle_radius_2) + particle_th)){
                    node_1 = 0;
                    node_2 = 0;
                    Dij_min = particle_th + 1;//Just to initially have a value greater than th
                    for (INT32 n = 0; n < particle_num_nodes; ++n) {
                        for (INT32 m = 0; m < particle_num_nodes; ++m) {
                            Dij = sqrt(SQ(particle_list[np].node[n].x - particle_list[np2].node[m].x) + SQ(particle_list[np].node[n].y - particle_list[np2].node[m].y));
                            if (Dij <= Dij_min) {
                                Dij_min = Dij;
                                node_1 = n;
                                node_2 = m;
                            }
                        }
                    }

                    if (Dij_min >min_th && Dij_min <= particle_th){
                        ffx = 1. * particle_cij / particle_ep * (particle_list[np].node[node_1].x - particle_list[np2].node[node_2].x) / Dij_min * pow((Dij_min - particle_th) / particle_th, 2);
                        ffy = 1. * particle_cij / particle_ep * (particle_list[np].node[node_1].y - particle_list[np2].node[node_2].y) / Dij_min * pow((Dij_min - particle_th) / particle_th, 2);
                        particle_list[np].center[0].f3x += ffx;
                        particle_list[np].center[0].f3y += ffy;
                        particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                    }
                    if (Dij_min <= min_th){
                        ffx = 1. * particle_cij / particle_Ep * (particle_list[np].node[node_1].x - particle_list[np2].node[node_2].x) / Dij_min * pow((Dij_min - particle_th) / particle_th, 2);
                        ffy = 1. * particle_cij / particle_Ep * (particle_list[np].node[node_1].y - particle_list[np2].node[node_2].y) / Dij_min * pow((Dij_min - particle_th) / particle_th, 2);
                        particle_list[np].center[0].f3x += ffx;
                        particle_list[np].center[0].f3y += ffy;
                        particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                    }
                }
            }
        }// END of Particle-Particle Collision

        double Rij, Rij_min;;
        if (Y_Wall)
        {
            if ((particle_list[np].center[0].y <= std::max(particle_radius, particle_radius_2) + particle_thw) || ((Ny - particle_list[np].center[0].y) <= std::max(particle_radius, particle_radius_2) + particle_thw)){
                node_1 = 0;
                Rij_min = particle_thw + 1;//Just to initially have a value greater than th
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    if (particle_list[np].node[n].y <= (Ny - 1) / 2.){//Particle-Y-Wall Collision section
                        Rij = abs(particle_list[np].node[n].y - 0);// 0 means y-location of lower wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    else
                    {
                        Rij = -1.0 * abs((Ny - 1) - particle_list[np].node[n].y);// Ny-1 means y-location of upper wall	.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    if (abs(Rij) <= abs(Rij_min)){
                        Rij_min = Rij;
                        node_1 = n;
                    }
                }
                if (abs(Rij_min) > min_thw && abs(Rij_min) <= particle_thw){
                    ffx4 = 0.0;
                    ffy4 = 1. * particle_cij / particle_epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
                if (abs(Rij_min) <= min_thw){
                    ffx4 = 0.0;
                    ffy4 = 1. * particle_cij / particle_Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
            }
        }// END of Particle_Y-Wall Collision
        if (X_Wall)
        {
            if ((particle_list[np].center[0].x <= std::max(particle_radius, particle_radius_2) + particle_thw) || ((Nx - particle_list[np].center[0].x) <= std::max(particle_radius, particle_radius_2) + particle_thw)){
                node_1 = 0;
                Rij_min = particle_thw + 1;//Just to initially have a value greater than th
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    if (particle_list[np].node[n].x <= (Nx - 1) / 2.){//Particle-Y-Wall Collision section
                        Rij = abs(particle_list[np].node[n].x - 0);// 0 means y-location of lower wall.Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    else
                    {
                        Rij = -1.0 * abs((Nx - 1) - particle_list[np].node[n].x);// Ny-1 means y-location of upper wall. Important: why abs is used here? If particle penetrates to the wall, abs causes true force direction
                    }
                    if (abs(Rij) <= abs(Rij_min)){
                        Rij_min = Rij;
                        node_1 = n;
                    }
                }
                if (abs(Rij_min) > min_thw && abs(Rij_min) <= particle_thw){
                    ffx4 = 1. * particle_cij / particle_epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    ffy4 = 0.0;
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
                if (abs(Rij_min) <= min_thw){
                    ffx4 = 1. * particle_cij / particle_Epw * Rij_min / abs(Rij_min) * pow((-abs(Rij_min) + particle_thw) / particle_thw, 2);
                    ffy4 = 0.0;
                    particle_list[np].center[0].f4x += ffx4;
                    particle_list[np].center[0].f4y += ffy4;
                    particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[node_1].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[node_1].x - particle_list[np].center[0].x);
                }
            }
        }// END of Particle_X-Wall Collision
    } //End of Ellipse_Collision_3
///********************************************************************************//
    if (Ellipse_Collision_4)
    {
        //Collision Forces. JUST for ELLIPSE shape
        double Dij, Dist, ffx, ffy, ffx4, ffy4;
        ffx = 0.0;
        ffy = 0.0;

        for (INT32 np2 = 0; np2 < num_particles; ++np2) {//Particle-Particle Collision
            if (np2 != np) {
                Dist = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y));
                if (Dist <= (2 * std::max(particle_radius, particle_radius_2) + particle_th + 2)){
                    for (INT32 n = 0; n < particle_num_nodes; ++n) {
                        for (INT32 m = 0; m < particle_num_nodes; ++m) {
                            Dij = sqrt(SQ(particle_list[np].node[n].x - particle_list[np2].node[m].x) + SQ(particle_list[np].node[n].y - particle_list[np2].node[m].y));
                            if (Dij <= particle_th) {
                                ffx = 1. * particle_cij / particle_ep * (particle_list[np].node[n].x - particle_list[np2].node[m].x) / Dij * pow((Dij - particle_th) / particle_th, 2);
                                ffy = 1. * particle_cij / particle_ep * (particle_list[np].node[n].y - particle_list[np2].node[m].y) / Dij * pow((Dij - particle_th) / particle_th, 2);
                                particle_list[np].center[0].f3x += ffx;
                                particle_list[np].center[0].f3y += ffy;
                                particle_list[np].center[0].t3_z += ffx * (particle_list[np].node[n].y - particle_list[np].center[0].y) - ffy * (particle_list[np].node[n].x - particle_list[np].center[0].x);
                            }


                        }
                    }
                }
            }
        }// END of Particle-Particle Collision

        double Rij;
        if (Y_Wall)
        {
            if ((particle_list[np].center[0].y <= std::max(particle_radius, particle_radius_2) + particle_thw + 2) || ((Ny - particle_list[np].center[0].y) <= std::max(particle_radius, particle_radius_2) + particle_thw + 2)){
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    if (particle_list[np].node[n].y <= (Ny - 1) / 2.){//Particle-Y-Wall Collision section
                        Rij = particle_list[np].node[n].y - 0;// 0 means y-location of lower wall
                    }
                    else
                    {
                        Rij = -1.0 * ((Ny - 1) - particle_list[np].node[n].y);// Ny-1 means y-location of upper wall
                    }
                    if (abs(Rij) <= (particle_thw)){
                        ffx4 = 0.0;
                        ffy4 = 1. * particle_cij / particle_epw * Rij / abs(Rij) * pow((-abs(Rij) + particle_thw) / particle_thw, 2);
                        particle_list[np].center[0].f4x += ffx4;
                        particle_list[np].center[0].f4y += ffy4;
                        particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[n].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[n].x - particle_list[np].center[0].x);
                    }
                }
            }
        }
        if (X_Wall)
        {
            if ((particle_list[np].center[0].x <= std::max(particle_radius, particle_radius_2) + particle_thw + 2) || ((Nx - particle_list[np].center[0].x) <= std::max(particle_radius, particle_radius_2) + particle_thw + 2))
            {
                for (INT32 n = 0; n < particle_num_nodes; ++n) {
                    if (particle_list[np].node[n].x <= (Nx - 1) / 2.){//Particle-X-Wall Collision section
                        Rij = particle_list[np].node[n].x - 0;// 0 means x-location of left wall
                    }
                    else
                    {
                        Rij = -1.0 * ((Nx - 1) - particle_list[np].node[n].x);// Nx-1 means x-location of right wall
                    }

                    if (abs(Rij) <= (particle_thw)){
                        ffx4 = 1. * particle_cij / particle_epw * Rij / abs(Rij) * pow((-abs(Rij) + particle_thw) / particle_thw, 2);
                        ffy4 = 0;
                        particle_list[np].center[0].f4x += ffx4;
                        particle_list[np].center[0].f4y += ffy4;
                        particle_list[np].center[0].t3_z += ffx4 * (particle_list[np].node[n].y - particle_list[np].center[0].y) - ffy4 * (particle_list[np].node[n].x - particle_list[np].center[0].x);
                    }
                }
            }
        }	 // END of Particle-Wall Collision
    } //end of Ellipse_Collision_4
	return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_SPHERE
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for spherical objects

void ParticleResolved::collision_sphere(INT32 np)
{

    if (Sphere_Collision_1_Feng) // Verlet list is not used. Possibility of collision with ALL particles is checked.
    {
        /// Particle-particle collision forces. (Sphere)
        double Dij;
        for (INT32 np2 = 0; np2 < num_particles; ++np2)  //Particle-Particle Collision
        {
            if (np2 != np)
            {
                Dij = sqrt(SQ(particle_list[np].center[0].x - particle_list[np2].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[np2].center[0].y) + SQ(particle_list[np].center[0].z - particle_list[np2].center[0].z));

                if (Dij <= (2.0 * particle_list[np].radius + particle_th) && Dij > (2.0 * particle_list[np].radius))
                {
                    particle_list[np].center[0].f3x += particle_cij / particle_ep * ((particle_list[np].center[0].x - particle_list[np2].center[0].x)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3y += particle_cij / particle_ep * ((particle_list[np].center[0].y - particle_list[np2].center[0].y)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3z += particle_cij / particle_ep * ((particle_list[np].center[0].z - particle_list[np2].center[0].z)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                }
                if (Dij <= (2.0 * particle_list[np].radius))
                {
                    particle_list[np].center[0].f3x += (particle_list[np].center[0].x - particle_list[np2].center[0].x) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                    particle_list[np].center[0].f3y += (particle_list[np].center[0].y - particle_list[np2].center[0].y) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                    particle_list[np].center[0].f3z += (particle_list[np].center[0].z - particle_list[np2].center[0].z) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                }
            }
        }// END of Particle-Particle Collision

        /// Particle-Wall Collision
        for (INT32 dir = 1; dir <= 3; dir++)
        {
            double Rij{};
            double center_particle{};
            double Ndomain{};
            double vel_particle{};
            double dot_pro{};
            double collision_wall_force{};
            switch(dir)
            {
                case 1:
                {
                    center_particle = particle_list[np].center[0].x;
                    Ndomain = Nx;
                    vel_particle = particle_list[np].center[0].vel_x;
                    break;
                }
                case 2:
                {
                    center_particle = particle_list[np].center[0].y;
                    Ndomain = Ny;
                    vel_particle = particle_list[np].center[0].vel_y;
                    break;
                }
                case 3:
                {
                    center_particle = particle_list[np].center[0].z;
                    Ndomain = Nz;
                    vel_particle = particle_list[np].center[0].vel_z;
                    break;
                }
                default:
                {
                    std::cout << "Invalid dir" << std::endl;
                    break;
                }
            }

            if (center_particle <= (Ndomain - 2) / 2)  //Calculate Rij
            {
                Rij = center_particle + particle_list[np].radius;
            }
            else
            {
                Rij = -1.0 * ((Ndomain - 1) - center_particle + particle_list[np].radius);
            }

            //Calculate wall collision force
            if (abs(Rij) <= (2.0 * particle_list[np].radius + particle_thw) && abs(Rij) > (2.0 * particle_list[np].radius))
            {
                collision_wall_force = particle_cij / particle_epw * pow((2.0 * particle_list[np].radius - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
            }
            if (abs(Rij) <= (2.0 * particle_list[np].radius))
            {
                collision_wall_force = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * particle_list[np].radius - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * particle_list[np].radius - abs(Rij))) * Rij / abs(Rij);
            }

            switch(dir)
            {
                case 1:
                {
                    if (X_Wall){
                        particle_list[np].center[0].f4x += collision_wall_force;
                    }
                    break;
                }
                case 2:
                {
                    if (Y_Wall){
                        particle_list[np].center[0].f4y += collision_wall_force;
                    }
                    break;
                }
                case 3:
                {
                    if (Z_Wall){
                        particle_list[np].center[0].f4z += collision_wall_force;
                    }
                    break;
                }
                default:
                {
                    std::cout << "Invalid dir" << std::endl;
                    break;
                }
            }
        }// END of Particle-Wall Collision
    } //end of Sphere_Collision_1_Feng

    ///############################################################################///
    if (Sphere_Collision_2_Feng_Verlet) // Verlet list is used. Possibility of collision is only checked with those in the Verlet list.
    {
        /// Particle-particle collision forces. (Sphere)
        double Dij;
        for (INT32 np2 = 0; np2 < verlet_list[np].size(); ++np2)  //Particle-Particle Collision
        {
            INT32 index_p = verlet_list[np][np2];
            if (index_p != np)
            {
                Dij = sqrt(SQ(particle_list[np].center[0].x - particle_list[index_p].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[index_p].center[0].y) + SQ(particle_list[np].center[0].z - particle_list[index_p].center[0].z));

                if (Dij <= (2.0 * particle_list[np].radius + particle_th) && Dij > (2.0 * particle_list[np].radius))
                {
                    particle_list[np].center[0].f3x += particle_cij / particle_ep * ((particle_list[np].center[0].x - particle_list[index_p].center[0].x)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3y += particle_cij / particle_ep * ((particle_list[np].center[0].y - particle_list[index_p].center[0].y)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3z += particle_cij / particle_ep * ((particle_list[np].center[0].z - particle_list[index_p].center[0].z)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                }
                if (Dij <= (2.0 * particle_list[np].radius))
                {
                    particle_list[np].center[0].f3x += (particle_list[np].center[0].x - particle_list[index_p].center[0].x) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                    particle_list[np].center[0].f3y += (particle_list[np].center[0].y - particle_list[index_p].center[0].y) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                    particle_list[np].center[0].f3z += (particle_list[np].center[0].z - particle_list[index_p].center[0].z) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                }
            }
        }// END of Particle-Particle Collision

        /// Particle-Wall Collision
        for (INT32 dir = 1; dir <= 3; dir++)
        {
            double Rij{};
            double center_particle{};
            double Ndomain{};
            double vel_particle{};
            double dot_pro{};
            double collision_wall_force{};
            switch(dir)
            {
                case 1:
                {
                    center_particle = particle_list[np].center[0].x;
                    Ndomain = Nx;
                    vel_particle = particle_list[np].center[0].vel_x;
                    break;
                }
                case 2:
                {
                    center_particle = particle_list[np].center[0].y;
                    Ndomain = Ny;
                    vel_particle = particle_list[np].center[0].vel_y;
                    break;
                }
                case 3:
                {
                    center_particle = particle_list[np].center[0].z;
                    Ndomain = Nz;
                    vel_particle = particle_list[np].center[0].vel_z;
                    break;
                }
                default:
                {
                    std::cout << "Invalid dir" << std::endl;
                    break;
                }
            }

            if (center_particle <= (Ndomain - 2) / 2)  //Calculate Rij
            {
                Rij = center_particle + particle_list[np].radius;
            }
            else
            {
                Rij = -1.0 * ((Ndomain - 1) - center_particle + particle_list[np].radius);
            }

            //Calculate wall collision force
            if (abs(Rij) <= (2.0 * particle_list[np].radius + particle_thw) && abs(Rij) > (2.0 * particle_list[np].radius))
            {
                collision_wall_force = particle_cij / particle_epw * pow((2.0 * particle_list[np].radius - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
            }
            if (abs(Rij) <= (2.0 * particle_list[np].radius))
            {
                collision_wall_force = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * particle_list[np].radius - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * particle_list[np].radius - abs(Rij))) * Rij / abs(Rij);
            }

            switch(dir)
            {
                case 1:
                {
                    if (X_Wall){
                        particle_list[np].center[0].f4x += collision_wall_force;
                    }
                    break;
                }
                case 2:
                {
                    if (Y_Wall){
                        particle_list[np].center[0].f4y += collision_wall_force;
                    }
                    break;
                }
                case 3:
                {
                    if (Z_Wall){
                        particle_list[np].center[0].f4z += collision_wall_force;
                    }
                    break;
                }
                default:
                {
                    std::cout << "Invalid dir" << std::endl;
                    break;
                }
            }
        }// END of Particle-Wall Collision
    } //end of Sphere_Collision_2_Feng_Verlet
    ///############################################################################///

    if (Sphere_Collision_3_Feng_Lub_Verlet)
    {
        /// Particle-particle collision forces. (Sphere)
        double Dij{};
        for (INT32 np2 = 0; np2 < verlet_list[np].size(); ++np2)  //Particle-Particle Collision
        {
            INT32 index_p = verlet_list[np][np2];
            if (index_p != np)
            {
                Dij = sqrt(SQ(particle_list[np].center[0].x - particle_list[index_p].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[index_p].center[0].y) + SQ(particle_list[np].center[0].z - particle_list[index_p].center[0].z));

                if (Dij <= (2.0 * particle_list[np].radius + particle_th_lub) && Dij > (2.0 * particle_list[np].radius + particle_th))
                {
                    double dot_pro;
                    dot_pro = (particle_list[np].center[0].vel_x - particle_list[index_p].center[0].vel_x) * (particle_list[np].center[0].x - particle_list[index_p].center[0].x);
                    particle_list[np].center[0].f3x += -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius * particle_radius / (particle_radius + particle_radius)) * (1. / (Dij - 2 * particle_radius) - 1. / particle_th_lub) * dot_pro / Dij;
                    dot_pro = (particle_list[np].center[0].vel_y - particle_list[index_p].center[0].vel_y) * (particle_list[np].center[0].y - particle_list[index_p].center[0].y);
                    particle_list[np].center[0].f3y += -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius * particle_radius / (particle_radius + particle_radius)) * (1. / (Dij - 2 * particle_radius) - 1. / particle_th_lub) * dot_pro / Dij;
                    dot_pro = (particle_list[np].center[0].vel_z - particle_list[index_p].center[0].vel_z) * (particle_list[np].center[0].z - particle_list[index_p].center[0].z);
                    particle_list[np].center[0].f3z += -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius * particle_radius / (particle_radius + particle_radius)) * (1. / (Dij - 2 * particle_radius) - 1. / particle_th_lub) * dot_pro / Dij;

                }

                if (Dij <= (2.0 * particle_list[np].radius + particle_th) && Dij > (2.0 * particle_list[np].radius))
                {
                    particle_list[np].center[0].f3x += particle_cij / particle_ep * ((particle_list[np].center[0].x - particle_list[index_p].center[0].x)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3y += particle_cij / particle_ep * ((particle_list[np].center[0].y - particle_list[index_p].center[0].y)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                    particle_list[np].center[0].f3z += particle_cij / particle_ep * ((particle_list[np].center[0].z - particle_list[index_p].center[0].z)) / Dij * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2);
                }
                if (Dij <= (2.0 * particle_list[np].radius))
                {
                    particle_list[np].center[0].f3x += (particle_list[np].center[0].x - particle_list[index_p].center[0].x) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                    particle_list[np].center[0].f3y += (particle_list[np].center[0].y - particle_list[index_p].center[0].y) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                    particle_list[np].center[0].f3z += (particle_list[np].center[0].z - particle_list[index_p].center[0].z) / Dij * (particle_cij / particle_ep * pow((Dij - 2.0 * particle_list[np].radius - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (2.0 * particle_list[np].radius - Dij));
                }
            }
        }// END of Particle-Particle Collision

        /// Particle-Wall Collision
        for (INT32 dir = 1; dir <= 3; dir++)
        {
            double Rij{};
            double center_particle{};
            double Ndomain{};
            double vel_particle{};
            double dot_pro{};
            double collision_wall_force{};
            switch(dir)
            {
                case 1:
                {
                    center_particle = particle_list[np].center[0].x;
                    Ndomain = Nx;
                    vel_particle = particle_list[np].center[0].vel_x;
                    break;
                }
                case 2:
                {
                    center_particle = particle_list[np].center[0].y;
                    Ndomain = Ny;
                    vel_particle = particle_list[np].center[0].vel_y;
                    break;
                }
                case 3:
                {
                    center_particle = particle_list[np].center[0].z;
                    Ndomain = Nz;
                    vel_particle = particle_list[np].center[0].vel_z;
                    break;
                }
                default:
                {
                    std::cout << "Invalid dir" << std::endl;
                    break;
                }
            }

            if (center_particle <= (Ndomain - 2) / 2)  //Calculate Rij
            {
                Rij = center_particle + particle_list[np].radius;
            }
            else
            {
                Rij = -1.0 * ((Ndomain - 1) - center_particle + particle_list[np].radius);
            }

            //Calculate wall collision force
            if (abs(Rij) <= (2.0 * particle_list[np].radius + particle_thw_lub) && abs(Rij) > (2.0 * particle_list[np].radius + particle_thw))
            {
                dot_pro = (vel_particle - 0) * Rij;
                collision_wall_force = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(particle_radius * particle_radius / (particle_radius + particle_radius)) * (1. / (abs(Rij) - 2 * particle_radius) - 1. / particle_thw_lub) * dot_pro / abs(Rij);
            }
            if (abs(Rij) <= (2.0 * particle_list[np].radius + particle_thw) && abs(Rij) > (2.0 * particle_list[np].radius))
            {
                collision_wall_force = particle_cij / particle_epw * pow((2.0 * particle_list[np].radius - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
            }
            if (abs(Rij) <= (2.0 * particle_list[np].radius))
            {
                collision_wall_force = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * particle_list[np].radius - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * particle_list[np].radius - abs(Rij))) * Rij / abs(Rij);
            }

            switch(dir)
            {
                case 1:
                {
                    if (X_Wall){
                        particle_list[np].center[0].f4x += collision_wall_force;
                    }
                    break;
                }
                case 2:
                {
                    if (Y_Wall){
                        particle_list[np].center[0].f4y += collision_wall_force;
                    }
                    break;
                }
                case 3:
                {
                    if (Z_Wall){
                        particle_list[np].center[0].f4z += collision_wall_force;
                    }
                    break;
                }
                default:
                {
                    std::cout << "Invalid dir" << std::endl;
                    break;
                }
            }
        }// END of Particle-Wall Collision
    } //end of Sphere_Collision_3_Feng_Lub_Verlet

    return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: COLLISON_SPHEROID
//-----------------------------------------------------------------------------------------------------------------------------
// calculation the collision forces for Spheroidal objects
void ParticleResolved::collision_spheroid(INT32 np)
{

    //#if defined Spheroid_Collision_1
    //Collision Forces. for Spheroid shape
    double Dij{}, Dist{}, ffx{}, ffy{}, ffz{}, ffx4{}, ffy4{}, ffz4{}, norvec{};
    double xnor_init{}, ynor_init{}, znor_init{}, xnor_slope{}, ynor_slope{}, znor_slope{};
    double x_out{}, y_out{}, z_out{};
    double next_point_x{}, next_point_y{}, next_point_z{};
    double m_slope{}, n_slope{}, p_slope{};
    double x_intersect1{}, y_intersect1{}, z_intersect1{};
    double x_intersect2{}, y_intersect2{}, z_intersect2{};
    double x_intersect_real{}, y_intersect_real{}, z_intersect_real{};
    double tt_part1{}, tt_part2{}, tt_part3{};
    double tt_intersect_1{}, tt_intersect_2{};
    double x_0{}, y_0{}, z_0{}, slope_th{}, xp1{}, yp1{}, zp1{};
    double a{}, b{}, radius_check{};
    double x_0_inertial{}, y_0_inertial{}, z_0_inertial{};
    double xnor_slope_old{}, ynor_slope_old{}, znor_slope_old{};
    double Rij{};
    double dist_to_wall{};
    double M_curvature{}, N_curvature{}, R_curvature{}, fi_curvature{};
    double x_intersect_real_inertial{}, y_intersect_real_inertial{}, z_intersect_real_inertial{};
    double check_Y{};
    double u_approach_x{}, u_approach_y{}, u_approach_z{};
    double x_center_curve{}, y_center_curve{}, z_center_curve{};
    double x_center_curve_inertial{}, y_center_curve_inertial{}, z_center_curve_inertial{};
    double nor_slope_length{};
    INT32 mmm{};

    slope_th = 0.001;//thereshold of slope. If less than this value, computation stops

    a = particle_radius;
    b = particle_radius_2;
    radius_check = std::min(SQ(a) / b, SQ(b) / a);
    std::vector<INT32> check_walls;

    //////// Particle-Wall collision  //////////////////
    //----------------------------------------------------//
    if (X_Wall || Y_Wall || Z_Wall)
    {

        double thw_control{};
        if (Spheroid_Collision_1_Feng || Spheroid_Collision_2_Feng_Verlet)
        {
            thw_control = particle_thw;
        }
        if (Spheroid_Collision_3_Feng_Lub_Verlet)
        {
            thw_control = particle_thw_lub;
        }

        if (Y_Wall)
        {
            if ((particle_list[np].center[0].y <= std::max(particle_radius, particle_radius_2) + thw_control))
            {
                check_walls.push_back(1);//1 means lower Y-wall
            }
            if (particle_list[np].center[0].y >= Ny - 1 - (std::max(particle_radius, particle_radius_2) + thw_control))
            {
                check_walls.push_back(2);//2 means upper Y-wall
            }
        }
        if (X_Wall)
        {
            if ((particle_list[np].center[0].x <= std::max(particle_radius, particle_radius_2) + thw_control))
            {
                check_walls.push_back(3);//3 means lower X-wall
            }
            if (particle_list[np].center[0].x >= Nx - 1 - (std::max(particle_radius, particle_radius_2) + thw_control))
            {
                check_walls.push_back(4);//4 means upper X-wall
            }
        }
        if (Z_Wall)
        {
            if ((particle_list[np].center[0].z <= std::max(particle_radius, particle_radius_2) + thw_control))
            {
                check_walls.push_back(5);//5 means lower Z-wall
            }
            if (particle_list[np].center[0].z >= Nz - 1 - (std::max(particle_radius, particle_radius_2) + thw_control))
            {
                check_walls.push_back(6);//6 means upper Z-wall
            }
        }
        //////////////////////// Finding nearest point to wall and calculate the collision force  /////////////////////
        for (INT32 i_ch = 0; i_ch < check_walls.size(); i_ch++)
        {
            dist_to_wall = 1000;

            switch (check_walls[i_ch])
            {
            case 1:  //Lower Y-wall
            {
                xp1 = 0.0;//wall point in inertial coordinate minus particle center position
                yp1 = 0 - particle_list[np].center[0].y;
                zp1 = 0.0;

                xnor_slope = 0.0;
                ynor_slope = 1.0;
                znor_slope = 0.0;
                break;
            }
            case 2:  //Upper Y-wall
            {
                xp1 = 0.0;//wall point in inertial coordinate minus particle center position
                yp1 = Ny - 1 - particle_list[np].center[0].y;
                zp1 = 0.0;

                xnor_slope = 0.0;
                ynor_slope = 1.0;
                znor_slope = 0.0;
                break;
            }
            case 3:  //Lower X-wall
            {
                xp1 = 0 - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                yp1 = 0.0;
                zp1 = 0.0;

                xnor_slope = 1.0;
                ynor_slope = 0.0;
                znor_slope = 0.0;
                break;
            }
            case 4:  //Upper X-wall
            {
                xp1 = Nx - 1 - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                yp1 = 0.0;
                zp1 = 0.0;

                xnor_slope = 1.0;
                ynor_slope = 0.0;
                znor_slope = 0.0;
                break;
            }
            case 5:  //Lower Z-wall
            {
                xp1 = 0.0;//wall point in inertial coordinate minus particle center position
                yp1 = 0.0;
                zp1 = 0 - particle_list[np].center[0].z;

                xnor_slope = 0.0;
                ynor_slope = 0.0;
                znor_slope = 1.0;
                break;
            }
            case 6:  //Lower Z-wall
            {
                xp1 = 0.0;//wall point in inertial coordinate minus particle center position
                yp1 = 0.0;
                zp1 = Nz - 1 - particle_list[np].center[0].z;

                xnor_slope = 0.0;
                ynor_slope = 0.0;
                znor_slope = 1.0;
                break;
            }
            default:
                std::cout << "Kein Antwort";
            }

            x_0 = 0;//means particle center position in body-fixed frame
            y_0 = 0;//means particle center position in body-fixed frame
            z_0 = 0;//means particle center position in body-fixed frame

            mmm = 0;

            do
            {
                x_out = particle_list[np].center[0].quater[1] * (xp1)+particle_list[np].center[0].quater[2] * (yp1)+particle_list[np].center[0].quater[3] * (zp1);//location of wall point in body-fixed frame
                y_out = particle_list[np].center[0].quater[4] * (xp1)+particle_list[np].center[0].quater[5] * (yp1)+particle_list[np].center[0].quater[6] * (zp1);
                z_out = particle_list[np].center[0].quater[7] * (xp1)+particle_list[np].center[0].quater[8] * (yp1)+particle_list[np].center[0].quater[9] * (zp1);
                m_slope = x_out - x_0;//slope of the line conecting the wall and the center of imaginary circle
                n_slope = y_out - y_0;
                p_slope = z_out - z_0;

                /// Intersection of ellipsoid and line, Obtained by MATLAB  ///
                tt_part1 = a * sqrt(SQ(a) * SQ(b) * SQ(n_slope) + SQ(a) * SQ(b) * SQ(p_slope) - SQ(a) * SQ(n_slope) * SQ(z_0) + 2 * SQ(a) * n_slope * p_slope * y_0 * z_0 - SQ(a * p_slope * y_0) + SQ(SQ(b)) * SQ(m_slope) - SQ(b * m_slope * y_0) - SQ(b * m_slope * z_0) + 2 * SQ(b) * m_slope * n_slope * x_0 * y_0 + 2 * SQ(b) * m_slope * p_slope * x_0 * z_0 - SQ(b * n_slope * x_0) - SQ(b * p_slope * x_0));
                tt_part2 = SQ(b) * m_slope * x_0 + SQ(a) * n_slope * y_0 + SQ(a) * p_slope * z_0;
                tt_part3 = (SQ(a * n_slope) + SQ(a * p_slope) + SQ(b * m_slope));
                tt_intersect_1 = -(tt_part1 + tt_part2) / tt_part3;//Answers of equation
                tt_intersect_2 = -(-tt_part1 + tt_part2) / tt_part3;

                x_intersect1 = x_0 + m_slope * tt_intersect_1;//in body fixed frame
                y_intersect1 = y_0 + n_slope * tt_intersect_1;
                z_intersect1 = z_0 + p_slope * tt_intersect_1;

                x_intersect2 = x_0 + m_slope * tt_intersect_2;
                y_intersect2 = y_0 + n_slope * tt_intersect_2;
                z_intersect2 = z_0 + p_slope * tt_intersect_2;

                /// Find the nearst point
                if ((SQ(x_intersect1 - x_out) + SQ(y_intersect1 - y_out) + SQ(z_intersect1 - z_out)) < (SQ(x_intersect2 - x_out) + SQ(y_intersect2 - y_out) + SQ(z_intersect2 - z_out)))
                {
                    x_intersect_real = x_intersect1;//in body fixed frame
                    y_intersect_real = y_intersect1;//in body fixed frame
                    z_intersect_real = z_intersect1;//in body fixed frame
                }
                else
                {
                    x_intersect_real = x_intersect2;//in body fixed frame
                    y_intersect_real = y_intersect2;//in body fixed frame
                    z_intersect_real = z_intersect2;//in body fixed frame
                }
                xnor_slope_old = xnor_slope;
                ynor_slope_old = ynor_slope;
                znor_slope_old = znor_slope;

                xnor_slope = -2. * x_intersect_real / particle_radius;//slope of normal vector on spheroid in body_fixed frame. INWARD direction
                ynor_slope = -2. * y_intersect_real / particle_radius_2;
                znor_slope = -2. * z_intersect_real / particle_radius_2;

                nor_slope_length = (sqrt(SQ(xnor_slope) + SQ(ynor_slope) + SQ(znor_slope)));
                xnor_slope /= nor_slope_length;
                ynor_slope /= nor_slope_length;
                znor_slope /= nor_slope_length;

                x_0 = x_intersect_real + xnor_slope * radius_check;//center of an imaginary circle for next check. (body-fixed frame)
                y_0 = y_intersect_real + ynor_slope * radius_check;
                z_0 = z_intersect_real + znor_slope * radius_check;

                //finding the center of imaginary circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
                update_rot_spheroid_collision(particle_list[np].center[0].deter, x_0_inertial, y_0_inertial, z_0_inertial, particle_list[np].center[0].quater, x_0, y_0, z_0, particle_list[np].center[0].x, particle_list[np].center[0].y, particle_list[np].center[0].z);// Find x_0 position in inertial frame

                switch (check_walls[i_ch])
                {
                case 1:
                {
                    xp1 = x_0_inertial - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                    yp1 = 0 - particle_list[np].center[0].y;
                    zp1 = z_0_inertial - particle_list[np].center[0].z;
                    break;
                }
                case 2:
                {
                    xp1 = x_0_inertial - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                    yp1 = Ny - 1 - particle_list[np].center[0].y;
                    zp1 = z_0_inertial - particle_list[np].center[0].z;
                    break;
                }
                case 3:
                {
                    xp1 = 0 - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                    yp1 = y_0_inertial - particle_list[np].center[0].y;
                    zp1 = z_0_inertial - particle_list[np].center[0].z;
                    break;
                }
                case 4:
                {
                    xp1 = Nx - 1 - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                    yp1 = y_0_inertial - particle_list[np].center[0].y;
                    zp1 = z_0_inertial - particle_list[np].center[0].z;
                    break;
                }
                case 5:
                {
                    xp1 = x_0_inertial - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                    yp1 = y_0_inertial - particle_list[np].center[0].y;
                    zp1 = 0 - particle_list[np].center[0].z;
                    break;
                }
                case 6:
                {
                    xp1 = x_0_inertial - particle_list[np].center[0].x;//wall point in inertial coordinate minus particle center position
                    yp1 = y_0_inertial - particle_list[np].center[0].y;
                    zp1 = Nz - 1 - particle_list[np].center[0].z;
                    break;
                }
                default:
                    cout << "Keine Antwort";
                }

                mmm++;
            }
            while (abs(xnor_slope - xnor_slope_old) >= slope_th || abs(ynor_slope - ynor_slope_old) >= slope_th || abs(znor_slope - znor_slope_old) >= slope_th);

            update_rot_spheroid_collision(particle_list[np].center[0].deter, x_intersect_real_inertial, y_intersect_real_inertial, z_intersect_real_inertial, particle_list[np].center[0].quater, x_intersect_real, y_intersect_real, z_intersect_real, particle_list[np].center[0].x, particle_list[np].center[0].y, particle_list[np].center[0].z);// Find x_0 position in inertial frame

            switch (check_walls[i_ch])
            {
            case 1:
            {
                dist_to_wall = y_intersect_real_inertial - 0;
                break;
            }
            case 2:
            {
                dist_to_wall = Ny - 1 - y_intersect_real_inertial;
                break;
            }
            case 3:
            {
                dist_to_wall = x_intersect_real_inertial - 0;
                break;
            }
            case 4:
            {
                dist_to_wall = Nx - 1 - x_intersect_real_inertial;
                break;
            }
            case 5:
            {
                dist_to_wall = z_intersect_real_inertial - 0;
                break;
            }
            case 6:
            {
                dist_to_wall = Nz - 1 - z_intersect_real_inertial;
                break;
            }
            default:
                cout << "Keine Antwort";
            }

            if (dist_to_wall < thw_control)
            {
                //cout << time << endl;
                fi_curvature = asin(x_intersect_real / (sqrt(SQ(x_intersect_real) + SQ(y_intersect_real) + SQ(z_intersect_real))));
                //cout << fi_curvature * 180.0 / M_PI<<endl;
                R_curvature = SQ(a) * b / (SQ(a * sin(fi_curvature)) + SQ(b * cos(fi_curvature))); //Paper 1366

                x_center_curve = x_intersect_real + xnor_slope * R_curvature;//center of circle of curvature. (body-fixed frame)
                y_center_curve = y_intersect_real + ynor_slope * R_curvature;
                z_center_curve = z_intersect_real + znor_slope * R_curvature;

                //finding the center of curvature circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
                update_rot_spheroid_collision(particle_list[np].center[0].deter, x_center_curve_inertial, y_center_curve_inertial, z_center_curve_inertial, particle_list[np].center[0].quater, x_center_curve, y_center_curve, z_center_curve, particle_list[np].center[0].x, particle_list[np].center[0].y, particle_list[np].center[0].z);// Find  position of circle of curvature in inertial frame

                switch (check_walls[i_ch])
                {
                case 1:
                {
                    Rij = y_center_curve_inertial + R_curvature;
                    break;
                }
                case 2:
                {
                    Rij = -1.0 * ((Ny - 1) - y_center_curve_inertial + R_curvature);
                    break;
                }
                case 3:
                {
                    Rij = x_center_curve_inertial + R_curvature;

                    break;
                }
                case 4:
                {
                    //double ooo = -1.0 * ((Nx - 1) - x_intersect_real_inertial + 2 * R_curvature);
                    Rij = -1.0 * ((Nx - 1) - x_center_curve_inertial + R_curvature);
                    break;
                }
                case 5:
                {
                    Rij = z_center_curve_inertial + R_curvature;

                    break;
                }
                case 6:
                {
                    Rij = -1.0 * ((Nz - 1) - z_center_curve_inertial + R_curvature);

                    break;
                }
                default:
                    std::cout << "Keine Antwort";
                }

                ////////////// Force and torque calculation
                switch (check_walls[i_ch])
                {
                case 1:   // Lower Y-Wall
                {
                    if (Spheroid_Collision_3_Feng_Lub_Verlet)
                    {
                        if (abs(Rij) <= (2.0 * R_curvature + particle_thw_lub) && abs(Rij) > (2.0 * R_curvature + particle_thw))
                        {
                            u_approach_y = particle_list[np].center[0].vel_y - particle_list[np].center[0].omgZ * (x_center_curve_inertial - particle_list[np].center[0].x) + particle_list[np].center[0].omgX * (z_center_curve_inertial - particle_list[np].center[0].z);
                            double dot_pro = (u_approach_y - 0) * Rij;
                            ffx4 = 0.0;
                            ffy4 = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature * R_curvature / (R_curvature + R_curvature)) * (1. / (abs(Rij) - 2 * R_curvature) - 1. / particle_thw_lub) * dot_pro * Rij / SQ(Rij);
                            ffz4 = 0.0;
                        }
                    }
                    if (abs(Rij) <= (2.0 * R_curvature + particle_thw) && abs(Rij) > (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = particle_cij / particle_epw * pow((2.0 * R_curvature - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
                        ffz4 = 0.0;
                    }
                    if (abs(Rij) <= (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * R_curvature - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * R_curvature - abs(Rij))) * Rij / abs(Rij);
                        ffz4 = 0.0;
                    }
                    break;
                }
                case 2:  // Upper Y-Wall
                {
                    if (Spheroid_Collision_3_Feng_Lub_Verlet)
                    {
                        if (abs(Rij) <= (2.0 * R_curvature + particle_thw_lub) && abs(Rij) > (2.0 * R_curvature + particle_thw))
                        {
                            u_approach_y = particle_list[np].center[0].vel_y - particle_list[np].center[0].omgZ * (x_center_curve_inertial - particle_list[np].center[0].x) + particle_list[np].center[0].omgX * (z_center_curve_inertial - particle_list[np].center[0].z);
                            double dot_pro = (u_approach_y - 0) * Rij;
                            ffx4 = 0.0;
                            ffy4 = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature * R_curvature / (R_curvature + R_curvature)) * (1. / (abs(Rij) - 2 * R_curvature) - 1. / particle_thw_lub) * dot_pro * Rij / SQ(Rij);
                            ffz4 = 0.0;
                        }
                    }
                    if (abs(Rij) <= (2.0 * R_curvature + particle_thw) && abs(Rij) > (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = particle_cij / particle_epw * pow((2.0 * R_curvature - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
                        ffz4 = 0.0;
                    }
                    if (abs(Rij) <= (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * R_curvature - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * R_curvature - abs(Rij))) * Rij / abs(Rij);
                        ffz4 = 0.0;
                    }
                    break;
                }
                case 3:  // Lower X-Wall
                {
                    if (Spheroid_Collision_3_Feng_Lub_Verlet)
                    {
                        if (abs(Rij) <= (2.0 * R_curvature + particle_thw_lub) && abs(Rij) > (2.0 * R_curvature + particle_thw))
                        {
                            u_approach_x = particle_list[np].center[0].vel_x + particle_list[np].center[0].omgZ * (y_center_curve_inertial - particle_list[np].center[0].y) - particle_list[np].center[0].omgY * (z_center_curve_inertial - particle_list[np].center[0].z);
                            double dot_pro = (u_approach_x - 0) * Rij;
                            ffx4 = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature * R_curvature / (R_curvature + R_curvature)) * (1. / (abs(Rij) - 2 * R_curvature) - 1. / particle_thw_lub) * dot_pro * Rij / SQ(Rij);
                            ffy4 = 0.0;
                            ffz4 = 0.0;
                        }
                    }
                    if (abs(Rij) <= (2.0 * R_curvature + particle_thw) && abs(Rij) > (2.0 * R_curvature))
                    {
                        ffx4 = particle_cij / particle_epw * pow((2.0 * R_curvature - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
                        ffy4 = 0.0;
                        ffz4 = 0.0;
                    }
                    if (abs(Rij) <= (2.0 * R_curvature))
                    {
                        ffx4 = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * R_curvature - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * R_curvature - abs(Rij))) * Rij / abs(Rij);
                        ffy4 = 0.0;
                        ffz4 = 0.0;
                    }
                    break;
                }
                case 4:  // Upper X-Wall
                {
                    if (Spheroid_Collision_3_Feng_Lub_Verlet)
                    {
                        if (abs(Rij) <= (2.0 * R_curvature + particle_thw_lub) && abs(Rij) > (2.0 * R_curvature + particle_thw))
                        {
                            u_approach_x = particle_list[np].center[0].vel_x + particle_list[np].center[0].omgZ * (y_center_curve_inertial - particle_list[np].center[0].y) - particle_list[np].center[0].omgY * (z_center_curve_inertial - particle_list[np].center[0].z);
                            //cout << "salam A " << x_intersect_real_inertial << " " << dist_to_wall << " " << " " << Rij << " " << R_curvature << " " << x_0_inertial<<endl;
                            double dot_pro = (u_approach_x - 0) * Rij;
                            ffx4 = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature * R_curvature / (R_curvature + R_curvature)) * (1. / (abs(Rij) - 2 * R_curvature) - 1. / particle_thw_lub) * dot_pro * Rij / SQ(Rij);
                            ffy4 = 0.0;
                            ffz4 = 0.0;
                        }
                    }
                    if (abs(Rij) <= (2.0 * R_curvature + particle_thw) && abs(Rij) > (2.0 * R_curvature))
                    {
                        ffx4 = particle_cij / particle_epw * pow((2.0 * R_curvature - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
                        ffy4 = 0.0;
                        ffz4 = 0.0;
                    }
                    if (abs(Rij) <= (2.0 * R_curvature))
                    {
                        ffx4 = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * R_curvature - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * R_curvature - abs(Rij))) * Rij / abs(Rij);
                        ffy4 = 0.0;
                        ffz4 = 0.0;
                    }
                    break;
                }
                case 5:  // Lower Z-Wall
                {
                    if (Spheroid_Collision_3_Feng_Lub_Verlet)
                    {
                        if (abs(Rij) <= (2.0 * R_curvature + particle_thw_lub) && abs(Rij) > (2.0 * R_curvature + particle_thw))
                        {
                            u_approach_z = particle_list[np].center[0].vel_z - particle_list[np].center[0].omgX * (y_center_curve_inertial - particle_list[np].center[0].y) + particle_list[np].center[0].omgY * (x_center_curve_inertial - particle_list[np].center[0].x);
                            double dot_pro = (u_approach_z - 0) * Rij;
                            ffx4 = 0.0;
                            ffy4 = 0.0;
                            ffz4 = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature * R_curvature / (R_curvature + R_curvature)) * (1. / (abs(Rij) - 2 * R_curvature) - 1. / particle_thw_lub) * dot_pro * Rij / SQ(Rij);
                        }
                    }
                    if (abs(Rij) <= (2.0 * R_curvature + particle_thw) && abs(Rij) > (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = 0.0;
                        ffz4 = particle_cij / particle_epw * pow((2.0 * R_curvature - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
                    }
                    if (abs(Rij) <= (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = 0.0;
                        ffz4 = (particle_cij / particle_epw * pow((abs(Rij) - 2.0 * R_curvature - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * R_curvature - abs(Rij))) * Rij / abs(Rij);
                    }
                    break;
                }
                case 6:  // Upper Z-Wall
                {
                    if (Spheroid_Collision_3_Feng_Lub_Verlet)
                    {
                        if (abs(Rij) <= (2.0 * R_curvature + particle_thw_lub) && abs(Rij) > (2.0 * R_curvature + particle_thw))
                        {
                            u_approach_z = particle_list[np].center[0].vel_z - particle_list[np].center[0].omgX * (y_center_curve_inertial - particle_list[np].center[0].y) + particle_list[np].center[0].omgY * (x_center_curve_inertial - particle_list[np].center[0].x);
                            double dot_pro = (u_approach_z - 0) * Rij;
                            ffx4 = 0.0;
                            ffy4 = 0.0;
                            ffz4 = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature * R_curvature / (R_curvature + R_curvature)) * (1. / (abs(Rij) - 2 * R_curvature) - 1. / particle_thw_lub) * dot_pro * Rij / SQ(Rij);
                        }
                    }
                    if (abs(Rij) <= (2.0 * R_curvature + particle_thw) && abs(Rij) > (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = 0.0;
                        ffz4 = particle_cij / particle_epw * pow((2.0 * R_curvature - abs(Rij) + particle_thw) / particle_thw, 2) * Rij / abs(Rij);
                    }
                    if (abs(Rij) <= (2.0 * R_curvature))
                    {
                        ffx4 = 0.0;
                        ffy4 = 0.0;
                        ffz4 = (particle_cij /particle_epw * pow((abs(Rij) - 2.0 * R_curvature - particle_thw) / particle_thw, 2) + particle_cij / particle_Epw / particle_thw * (2.0 * R_curvature - abs(Rij))) * Rij / abs(Rij);
                    }
                    break;
                }
                default:
                    std::cout << "Keine Antwort";
                }

                particle_list[np].center[0].f4x += ffx4;
                particle_list[np].center[0].f4y += ffy4;
                particle_list[np].center[0].f4z += ffz4;

                particle_list[np].center[0].t3_x += ffy4 * (z_intersect_real_inertial - particle_list[np].center[0].z) - ffz4 * (y_intersect_real_inertial - particle_list[np].center[0].y);
                particle_list[np].center[0].t3_y += ffz4 * (x_intersect_real_inertial - particle_list[np].center[0].x) - ffx4 * (z_intersect_real_inertial - particle_list[np].center[0].z);
                particle_list[np].center[0].t3_z += ffx4 * (y_intersect_real_inertial - particle_list[np].center[0].y) - ffy4 * (x_intersect_real_inertial - particle_list[np].center[0].x);
                //cout << dist_to_wall << " " << particle_list[np].center[0].f4x << endl;

            }  // if (dist_to_wall<thw_control)
        }///for (INT32 i_ch = 0; i_ch < check_walls.size(); i_ch++)
    }  //Walls collision



    //------------------------------------------------------------------
    //------------------------------------------------------------------
    //------------------------------------------------------------------
    //-------------- Particle_Particle collision -----------------/////
    INT32 max_par_index{};
    INT32 index_p{};
    double par_par_dist{};
    double xp1_on_p2{}, yp1_on_p2{}, zp1_on_p2{};
    double x_0_of_p1{}, y_0_of_p1{}, z_0_of_p1{};
    double x_out_on_p2{}, y_out_on_p2{}, z_out_on_p2{};
    double m_slope_p1{}, n_slope_p1{}, p_slope_p1{};
    double x_intersect1_of_p1{}, y_intersect1_of_p1{}, z_intersect1_of_p1{};
    double x_intersect2_of_p1{}, y_intersect2_of_p1{}, z_intersect2_of_p1{};
    double x_intersect_real_of_p1{}, y_intersect_real_of_p1{}, z_intersect_real_of_p1{};
    double xnor_slope_old_of_p1{}, ynor_slope_old_of_p1{}, znor_slope_old_of_p1{};
    double xnor_slope_of_p1{}, ynor_slope_of_p1{}, znor_slope_of_p1{};
    double x_0_inertial_of_p1{}, y_0_inertial_of_p1{}, z_0_inertial_of_p1{};
    double x_intersect_real_inertial_of_p1{}, y_intersect_real_inertial_of_p1{}, z_intersect_real_inertial_of_p1{};
    double xp1_on_p1{}, yp1_on_p1{}, zp1_on_p1{};
    double temp_slope_length{};
    double xnor_slope_of_p2{}, ynor_slope_of_p2{}, znor_slope_of_p2{};
    double x_0_of_p2{}, y_0_of_p2{}, z_0_of_p2{};
    double x_out_on_p1{}, y_out_on_p1{}, z_out_on_p1{};
    double m_slope_p2{}, n_slope_p2{}, p_slope_p2{};
    double x_intersect1_of_p2{}, y_intersect1_of_p2{}, z_intersect1_of_p2{};
    double x_intersect2_of_p2{}, y_intersect2_of_p2{}, z_intersect2_of_p2{};
    double x_intersect_real_of_p2{}, y_intersect_real_of_p2{}, z_intersect_real_of_p2{};
    double xnor_slope_old_of_p2{}, ynor_slope_old_of_p2{}, znor_slope_old_of_p2{};
    double x_0_inertial_of_p2{}, y_0_inertial_of_p2{}, z_0_inertial_of_p2{};
    double x_intersect_real_inertial_of_p2{}, y_intersect_real_inertial_of_p2{}, z_intersect_real_inertial_of_p2{};
    double fi_curvature_p1{}, fi_curvature_p2{};
    double R_curvature_p1{}, R_curvature_p2{};
    double thp_control{};
    double x_center_curve_p1{}, y_center_curve_p1{}, z_center_curve_p1{};
    double x_center_curve_inertial_p1{}, y_center_curve_inertial_p1{}, z_center_curve_inertial_p1{};
    double x_center_curve_inertial_p2{}, y_center_curve_inertial_p2{}, z_center_curve_inertial_p2{};
    double x_center_curve_p2{}, y_center_curve_p2{}, z_center_curve_p2{};
    double u_approach_x_p1{}, u_approach_y_p1{}, u_approach_z_p1{};
    double u_approach_x_p2{}, u_approach_y_p2{}, u_approach_z_p2{};

    if (Spheroid_Collision_1_Feng || Spheroid_Collision_2_Feng_Verlet)
    {
        thp_control = particle_th;
    }
    if (Spheroid_Collision_3_Feng_Lub_Verlet)
    {
        thp_control = particle_th_lub;
    }

    ffx = 0.0;
    ffy = 0.0;
    ffz = 0.0;

    if (Spheroid_Collision_1_Feng)
    {
        max_par_index = num_particles;
    }
    if (Spheroid_Collision_2_Feng_Verlet || Spheroid_Collision_3_Feng_Lub_Verlet)
    {
        max_par_index = verlet_list[np].size();
    }

    for (INT32 np2 = 0; np2 < max_par_index; ++np2)  //Particle-Particle Collision
    {
        if (Spheroid_Collision_1_Feng)
        {
            index_p = np2;
        }
        if (Spheroid_Collision_2_Feng_Verlet || Spheroid_Collision_3_Feng_Lub_Verlet)
        {
            index_p = verlet_list[np][np2];
        }
        if (index_p != np)
        {
            par_par_dist = 1000;

            xp1_on_p2 = particle_list[index_p].center[0].x - particle_list[np].center[0].x;//point on P2 in inertial coordinate minus particle P1 center position
            yp1_on_p2 = particle_list[index_p].center[0].y - particle_list[np].center[0].y;
            zp1_on_p2 = particle_list[index_p].center[0].z - particle_list[np].center[0].z;

            xnor_slope_of_p1 = particle_list[np].center[0].x - particle_list[index_p].center[0].x;//slope of line connecting nearest points
            ynor_slope_of_p1 = particle_list[np].center[0].y - particle_list[index_p].center[0].y;
            znor_slope_of_p1 = particle_list[np].center[0].z - particle_list[index_p].center[0].z;

            temp_slope_length = sqrt(SQ(particle_list[np].center[0].x - particle_list[index_p].center[0].x) + SQ(particle_list[np].center[0].y - particle_list[index_p].center[0].y) + SQ(particle_list[np].center[0].z - particle_list[index_p].center[0].z));
            xnor_slope_of_p1 /= temp_slope_length;
            ynor_slope_of_p1 /= temp_slope_length;
            znor_slope_of_p1 /= temp_slope_length;

            x_0_of_p1 = 0;//means particle center position in body-fixed frame
            y_0_of_p1 = 0;//means particle center position in body-fixed frame
            z_0_of_p1 = 0;//means particle center position in body-fixed frame

            x_0_of_p2 = 0;//means particle center position in body-fixed frame
            y_0_of_p2 = 0;//means particle center position in body-fixed frame
            z_0_of_p2 = 0;//means particle center position in body-fixed frame

            xnor_slope_of_p2 = particle_list[index_p].center[0].x - particle_list[np].center[0].x;//slope of line connecting nearest points
            ynor_slope_of_p2 = particle_list[index_p].center[0].y - particle_list[np].center[0].y;
            znor_slope_of_p2 = particle_list[index_p].center[0].z - particle_list[np].center[0].z;

            //initial set. Changes later
            temp_slope_length = sqrt(SQ(particle_list[index_p].center[0].x - particle_list[np].center[0].x) + SQ(particle_list[index_p].center[0].y - particle_list[np].center[0].y) + SQ(particle_list[index_p].center[0].z - particle_list[np].center[0].z));
            xnor_slope_of_p2 /= temp_slope_length;
            ynor_slope_of_p2 /= temp_slope_length;
            znor_slope_of_p2 /= temp_slope_length;

            mmm = 0;

            do
            {
                x_out_on_p2 = particle_list[np].center[0].quater[1] * (xp1_on_p2)+particle_list[np].center[0].quater[2] * (yp1_on_p2)+particle_list[np].center[0].quater[3] * (zp1_on_p2);//location of P2 point in body-fixed frame of P1
                y_out_on_p2 = particle_list[np].center[0].quater[4] * (xp1_on_p2)+particle_list[np].center[0].quater[5] * (yp1_on_p2)+particle_list[np].center[0].quater[6] * (zp1_on_p2);
                z_out_on_p2 = particle_list[np].center[0].quater[7] * (xp1_on_p2)+particle_list[np].center[0].quater[8] * (yp1_on_p2)+particle_list[np].center[0].quater[9] * (zp1_on_p2);
                m_slope_p1 = x_out_on_p2 - x_0_of_p1; //slope of the line conecting the point on P2 and the center of imaginary circle of P1
                n_slope_p1 = y_out_on_p2 - y_0_of_p1;
                p_slope_p1 = z_out_on_p2 - z_0_of_p1;

                /// Intersection of ellipsoid and line, Obtained by MATLAB  ///
                tt_part1 = a * sqrt(SQ(a) * SQ(b) * SQ(n_slope_p1) + SQ(a) * SQ(b) * SQ(p_slope_p1) - SQ(a) * SQ(n_slope_p1) * SQ(z_0_of_p1) + 2 * SQ(a) * n_slope_p1 * p_slope_p1 * y_0_of_p1 * z_0_of_p1 - SQ(a * p_slope_p1 * y_0_of_p1) + SQ(SQ(b)) * SQ(m_slope_p1) - SQ(b * m_slope_p1 * y_0_of_p1) - SQ(b * m_slope_p1 * z_0_of_p1) + 2 * SQ(b) * m_slope_p1 * n_slope_p1 * x_0_of_p1 * y_0_of_p1 + 2 * SQ(b) * m_slope_p1 * p_slope_p1 * x_0_of_p1 * z_0_of_p1 - SQ(b * n_slope_p1 * x_0_of_p1) - SQ(b * p_slope_p1 * x_0_of_p1));
                tt_part2 = SQ(b) * m_slope_p1 * x_0_of_p1 + SQ(a) * n_slope_p1 * y_0_of_p1 + SQ(a) * p_slope_p1 * z_0_of_p1;
                tt_part3 = (SQ(a * n_slope_p1) + SQ(a * p_slope_p1) + SQ(b * m_slope_p1));
                tt_intersect_1 = -(tt_part1 + tt_part2) / tt_part3;//Answers of equation
                tt_intersect_2 = -(-tt_part1 + tt_part2) / tt_part3;

                x_intersect1_of_p1 = x_0_of_p1 + m_slope_p1 * tt_intersect_1;//in body fixed frame of P1
                y_intersect1_of_p1 = y_0_of_p1 + n_slope_p1 * tt_intersect_1;
                z_intersect1_of_p1 = z_0_of_p1 + p_slope_p1 * tt_intersect_1;

                x_intersect2_of_p1 = x_0_of_p1 + m_slope_p1 * tt_intersect_2;
                y_intersect2_of_p1 = y_0_of_p1 + n_slope_p1 * tt_intersect_2;
                z_intersect2_of_p1 = z_0_of_p1 + p_slope_p1 * tt_intersect_2;

                /// Find the nearst point
                if ((SQ(x_intersect1_of_p1 - x_out_on_p2) + SQ(y_intersect1_of_p1 - y_out_on_p2) + SQ(z_intersect1_of_p1 - z_out_on_p2)) < (SQ(x_intersect2_of_p1 - x_out_on_p2) + SQ(y_intersect2_of_p1 - y_out_on_p2) + SQ(z_intersect2_of_p1 - z_out_on_p2)))
                {
                    x_intersect_real_of_p1 = x_intersect1_of_p1;//in body fixed frame
                    y_intersect_real_of_p1 = y_intersect1_of_p1;//in body fixed frame
                    z_intersect_real_of_p1 = z_intersect1_of_p1;//in body fixed frame
                }
                else
                {
                    x_intersect_real_of_p1 = x_intersect2_of_p1;//in body fixed frame
                    y_intersect_real_of_p1 = y_intersect2_of_p1;//in body fixed frame
                    z_intersect_real_of_p1 = z_intersect2_of_p1;//in body fixed frame
                }
                xnor_slope_old_of_p1 = xnor_slope_of_p1;
                ynor_slope_old_of_p1 = ynor_slope_of_p1;
                znor_slope_old_of_p1 = znor_slope_of_p1;

                xnor_slope_of_p1 = -2. * x_intersect_real_of_p1 / particle_radius;//slope of normal vector on spheroid in body_fixed frame. INWARD direction
                ynor_slope_of_p1 = -2. * y_intersect_real_of_p1 / particle_radius_2;
                znor_slope_of_p1 = -2. * z_intersect_real_of_p1 / particle_radius_2;

                nor_slope_length = (sqrt(SQ(xnor_slope_of_p1) + SQ(ynor_slope_of_p1) + SQ(znor_slope_of_p1)));
                xnor_slope_of_p1 /= nor_slope_length;
                ynor_slope_of_p1 /= nor_slope_length;
                znor_slope_of_p1 /= nor_slope_length;

                x_0_of_p1 = x_intersect_real_of_p1 + xnor_slope_of_p1 * radius_check;//center of an imaginary circle for next check. (body-fixed frame)
                y_0_of_p1 = y_intersect_real_of_p1 + ynor_slope_of_p1 * radius_check;
                z_0_of_p1 = z_intersect_real_of_p1 + znor_slope_of_p1 * radius_check;

                //finding the center of imaginary circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
                update_rot_spheroid_collision(particle_list[np].center[0].deter, x_0_inertial_of_p1, y_0_inertial_of_p1, z_0_inertial_of_p1, particle_list[np].center[0].quater, x_0_of_p1, y_0_of_p1, z_0_of_p1, particle_list[np].center[0].x, particle_list[np].center[0].y, particle_list[np].center[0].z);// Find x_0 position in inertial frame

                //finding the intersection point of P1 in INERTIAL frame (x_intersect_real_inertial_of_p1, y_intersect_real_inertial_of_p1, z_intersect_real_inertial_of_p1)
                update_rot_spheroid_collision(particle_list[np].center[0].deter, x_intersect_real_inertial_of_p1, y_intersect_real_inertial_of_p1, z_intersect_real_inertial_of_p1, particle_list[np].center[0].quater, x_intersect_real_of_p1, y_intersect_real_of_p1, z_intersect_real_of_p1, particle_list[np].center[0].x, particle_list[np].center[0].y, particle_list[np].center[0].z);// Find x_0 position in inertial frame


                //////////----------------------------------------------------------------------//////////////////////
                //Now we start the process for P2
                xp1_on_p1 = x_intersect_real_inertial_of_p1 - particle_list[index_p].center[0].x;//point on P1 in inertial coordinate minus particle P2 center position
                yp1_on_p1 = y_intersect_real_inertial_of_p1 - particle_list[index_p].center[0].y;
                zp1_on_p1 = z_intersect_real_inertial_of_p1 - particle_list[index_p].center[0].z;

                x_out_on_p1 = particle_list[index_p].center[0].quater[1] * (xp1_on_p1)+particle_list[index_p].center[0].quater[2] * (yp1_on_p1)+particle_list[index_p].center[0].quater[3] * (zp1_on_p1);//location of P1 point in body-fixed frame of P2
                y_out_on_p1 = particle_list[index_p].center[0].quater[4] * (xp1_on_p1)+particle_list[index_p].center[0].quater[5] * (yp1_on_p1)+particle_list[index_p].center[0].quater[6] * (zp1_on_p1);
                z_out_on_p1 = particle_list[index_p].center[0].quater[7] * (xp1_on_p1)+particle_list[index_p].center[0].quater[8] * (yp1_on_p1)+particle_list[index_p].center[0].quater[9] * (zp1_on_p1);

                m_slope_p2 = x_out_on_p1 - x_0_of_p2; //slope of the line conecting the point on P1 and the center of imaginary circle of P2
                n_slope_p2 = y_out_on_p1 - y_0_of_p2;
                p_slope_p2 = z_out_on_p1 - z_0_of_p2;

                ///// Intersection of ellipsoid and line, Obtained by MATLAB  ///
                tt_part1 = a * sqrt(SQ(a) * SQ(b) * SQ(n_slope_p2) + SQ(a) * SQ(b) * SQ(p_slope_p2) - SQ(a) * SQ(n_slope_p2) * SQ(z_0_of_p2) + 2 * SQ(a) * n_slope_p2 * p_slope_p2 * y_0_of_p2 * z_0_of_p2 - SQ(a * p_slope_p2 * y_0_of_p2) + SQ(SQ(b)) * SQ(m_slope_p2) - SQ(b * m_slope_p2 * y_0_of_p2) - SQ(b * m_slope_p2 * z_0_of_p2) + 2 * SQ(b) * m_slope_p2 * n_slope_p2 * x_0_of_p2 * y_0_of_p2 + 2 * SQ(b) * m_slope_p2 * p_slope_p2 * x_0_of_p2 * z_0_of_p2 - SQ(b * n_slope_p2 * x_0_of_p2) - SQ(b * p_slope_p2 * x_0_of_p2));
                tt_part2 = SQ(b) * m_slope_p2 * x_0_of_p2 + SQ(a) * n_slope_p2 * y_0_of_p2 + SQ(a) * p_slope_p2 * z_0_of_p2;
                tt_part3 = (SQ(a * n_slope_p2) + SQ(a * p_slope_p2) + SQ(b * m_slope_p2));
                tt_intersect_1 = -(tt_part1 + tt_part2) / tt_part3;//Answers of equation
                tt_intersect_2 = -(-tt_part1 + tt_part2) / tt_part3;

                x_intersect1_of_p2 = x_0_of_p2 + m_slope_p2 * tt_intersect_1;//in body fixed frame of P1
                y_intersect1_of_p2 = y_0_of_p2 + n_slope_p2 * tt_intersect_1;
                z_intersect1_of_p2 = z_0_of_p2 + p_slope_p2 * tt_intersect_1;

                x_intersect2_of_p2 = x_0_of_p2 + m_slope_p2 * tt_intersect_2;
                y_intersect2_of_p2 = y_0_of_p2 + n_slope_p2 * tt_intersect_2;
                z_intersect2_of_p2 = z_0_of_p2 + p_slope_p2 * tt_intersect_2;

                ///// Find the nearst point
                if ((SQ(x_intersect1_of_p2 - x_out_on_p1) + SQ(y_intersect1_of_p2 - y_out_on_p1) + SQ(z_intersect1_of_p2 - z_out_on_p1)) < (SQ(x_intersect2_of_p2 - x_out_on_p1) + SQ(y_intersect2_of_p2 - y_out_on_p1) + SQ(z_intersect2_of_p2 - z_out_on_p1)))
                {
                    x_intersect_real_of_p2 = x_intersect1_of_p2;//in body fixed frame
                    y_intersect_real_of_p2 = y_intersect1_of_p2;//in body fixed frame
                    z_intersect_real_of_p2 = z_intersect1_of_p2;//in body fixed frame
                }
                else
                {
                    x_intersect_real_of_p2 = x_intersect2_of_p2;//in body fixed frame
                    y_intersect_real_of_p2 = y_intersect2_of_p2;//in body fixed frame
                    z_intersect_real_of_p2 = z_intersect2_of_p2;//in body fixed frame
                }

                xnor_slope_old_of_p2 = xnor_slope_of_p2;
                ynor_slope_old_of_p2 = ynor_slope_of_p2;
                znor_slope_old_of_p2 = znor_slope_of_p2;

                xnor_slope_of_p2 = -2. * x_intersect_real_of_p2 / particle_radius;//slope of normal vector on spheroid in body_fixed frame. INWARD direction
                ynor_slope_of_p2 = -2. * y_intersect_real_of_p2 / particle_radius_2;
                znor_slope_of_p2 = -2. * z_intersect_real_of_p2 / particle_radius_2;

                nor_slope_length = (sqrt(SQ(xnor_slope_of_p2) + SQ(ynor_slope_of_p2) + SQ(znor_slope_of_p2)));
                xnor_slope_of_p2 /= nor_slope_length;
                ynor_slope_of_p2 /= nor_slope_length;
                znor_slope_of_p2 /= nor_slope_length;

                x_0_of_p2 = x_intersect_real_of_p2 + xnor_slope_of_p2 * radius_check;//center of an imaginary circle for next check. (body-fixed frame)
                y_0_of_p2 = y_intersect_real_of_p2 + ynor_slope_of_p2 * radius_check;
                z_0_of_p2 = z_intersect_real_of_p2 + znor_slope_of_p2 * radius_check;

                ////finding the center of imaginary circle in INERTIAL frame (x_0_inertial, y_0_inertial, z_0_inertial)
                update_rot_spheroid_collision(particle_list[index_p].center[0].deter, x_0_inertial_of_p2, y_0_inertial_of_p2, z_0_inertial_of_p2, particle_list[index_p].center[0].quater, x_0_of_p2, y_0_of_p2, z_0_of_p2, particle_list[index_p].center[0].x, particle_list[index_p].center[0].y, particle_list[index_p].center[0].z);// Find x_0 position in inertial frame

                //finding the intersection point of P2 in INERTIAL frame (x_intersect_real_inertial_of_p2, y_intersect_real_inertial_of_p2, z_intersect_real_inertial_of_p2)
                update_rot_spheroid_collision(particle_list[index_p].center[0].deter, x_intersect_real_inertial_of_p2, y_intersect_real_inertial_of_p2, z_intersect_real_inertial_of_p2, particle_list[index_p].center[0].quater, x_intersect_real_of_p2, y_intersect_real_of_p2, z_intersect_real_of_p2, particle_list[index_p].center[0].x, particle_list[index_p].center[0].y, particle_list[index_p].center[0].z);// Find x_0 position in inertial frame

                xp1_on_p2 = x_intersect_real_inertial_of_p2 - particle_list[np].center[0].x;
                yp1_on_p2 = y_intersect_real_inertial_of_p2 - particle_list[np].center[0].y;
                zp1_on_p2 = z_intersect_real_inertial_of_p2 - particle_list[np].center[0].z;

                mmm++;
            }
            while (abs(xnor_slope_of_p1 - xnor_slope_old_of_p1) >= slope_th || abs(ynor_slope_of_p1 - ynor_slope_old_of_p1) >= slope_th || abs(znor_slope_of_p1 - znor_slope_old_of_p1) >= slope_th || abs(xnor_slope_of_p2 - xnor_slope_old_of_p2) >= slope_th || abs(ynor_slope_of_p2 - ynor_slope_old_of_p2) >= slope_th || abs(znor_slope_of_p2 - znor_slope_old_of_p2) >= slope_th);

            par_par_dist = sqrt(SQ(x_intersect_real_inertial_of_p1 - x_intersect_real_inertial_of_p2) + SQ(y_intersect_real_inertial_of_p1 - y_intersect_real_inertial_of_p2) + SQ(z_intersect_real_inertial_of_p1 - z_intersect_real_inertial_of_p2));

            if (par_par_dist < thp_control)
            {

                fi_curvature_p1 = asin(x_intersect_real_of_p1 / (sqrt(SQ(x_intersect_real_of_p1) + SQ(y_intersect_real_of_p1) + SQ(z_intersect_real_of_p1))));
                fi_curvature_p2 = asin(x_intersect_real_of_p2 / (sqrt(SQ(x_intersect_real_of_p2) + SQ(y_intersect_real_of_p2) + SQ(z_intersect_real_of_p2))));

                R_curvature_p1 = SQ(a) * b / (SQ(a * sin(fi_curvature_p1)) + SQ(b * cos(fi_curvature_p1))); //Paper 1366
                R_curvature_p2 = SQ(a) * b / (SQ(a * sin(fi_curvature_p2)) + SQ(b * cos(fi_curvature_p2))); //Paper 1366

                x_center_curve_p1 = x_intersect_real_of_p1 + xnor_slope_of_p1 * R_curvature_p1;//center of circle of curvature. (body-fixed frame)
                y_center_curve_p1 = y_intersect_real_of_p1 + ynor_slope_of_p1 * R_curvature_p1;
                z_center_curve_p1 = z_intersect_real_of_p1 + znor_slope_of_p1 * R_curvature_p1;

                x_center_curve_p2 = x_intersect_real_of_p2 + xnor_slope_of_p2 * R_curvature_p2;//center of circle of curvature. (body-fixed frame)
                y_center_curve_p2 = y_intersect_real_of_p2 + ynor_slope_of_p2 * R_curvature_p2;
                z_center_curve_p2 = z_intersect_real_of_p2 + znor_slope_of_p2 * R_curvature_p2;

                //finding the center of curvature circle in INERTIAL frame
                update_rot_spheroid_collision(particle_list[np].center[0].deter, x_center_curve_inertial_p1, y_center_curve_inertial_p1, z_center_curve_inertial_p1, particle_list[np].center[0].quater, x_center_curve_p1, y_center_curve_p1, z_center_curve_p1, particle_list[np].center[0].x, particle_list[np].center[0].y, particle_list[np].center[0].z);// Find  position of circle of curvature in inertial frame
                update_rot_spheroid_collision(particle_list[index_p].center[0].deter, x_center_curve_inertial_p2, y_center_curve_inertial_p2, z_center_curve_inertial_p2, particle_list[index_p].center[0].quater, x_center_curve_p2, y_center_curve_p2, z_center_curve_p2, particle_list[index_p].center[0].x, particle_list[index_p].center[0].y, particle_list[index_p].center[0].z);// Find  position of circle of curvature in inertial frame

                Dij = sqrt(SQ(x_center_curve_inertial_p1 - x_center_curve_inertial_p2) + SQ(y_center_curve_inertial_p1 - y_center_curve_inertial_p2) + SQ(z_center_curve_inertial_p1 - z_center_curve_inertial_p2));

                if (Spheroid_Collision_3_Feng_Lub_Verlet)
                {
                    if (par_par_dist <= particle_th_lub && par_par_dist > particle_th)
                    {
                        u_approach_x_p1 = particle_list[np].center[0].vel_x + particle_list[np].center[0].omgZ * (y_center_curve_inertial_p1 - particle_list[np].center[0].y) - particle_list[np].center[0].omgY * (z_center_curve_inertial_p1 - particle_list[np].center[0].z);
                        u_approach_y_p1 = particle_list[np].center[0].vel_y - particle_list[np].center[0].omgZ * (x_center_curve_inertial_p1 - particle_list[np].center[0].x) + particle_list[np].center[0].omgX * (z_center_curve_inertial_p1 - particle_list[np].center[0].z);
                        u_approach_z_p1 = particle_list[np].center[0].vel_z - particle_list[np].center[0].omgX * (y_center_curve_inertial_p1 - particle_list[np].center[0].y) + particle_list[np].center[0].omgY * (x_center_curve_inertial_p1 - particle_list[np].center[0].x);

                        u_approach_x_p2 = particle_list[index_p].center[0].vel_x + particle_list[np].center[0].omgZ * (y_center_curve_inertial_p2 - particle_list[index_p].center[0].y) - particle_list[np].center[0].omgY * (z_center_curve_inertial_p2 - particle_list[index_p].center[0].z);
                        u_approach_y_p2 = particle_list[index_p].center[0].vel_y - particle_list[np].center[0].omgZ * (x_center_curve_inertial_p2 - particle_list[index_p].center[0].x) + particle_list[np].center[0].omgX * (z_center_curve_inertial_p2 - particle_list[index_p].center[0].z);
                        u_approach_z_p2 = particle_list[index_p].center[0].vel_z - particle_list[np].center[0].omgX * (y_center_curve_inertial_p2 - particle_list[index_p].center[0].y) + particle_list[np].center[0].omgY * (x_center_curve_inertial_p2 - particle_list[index_p].center[0].x);

                        double dot_pro = ((u_approach_x_p1 - u_approach_x_p2) * (x_center_curve_inertial_p1 - x_center_curve_inertial_p2) + (u_approach_y_p1 - u_approach_y_p2) * (y_center_curve_inertial_p1 - y_center_curve_inertial_p2) + (u_approach_z_p1 - u_approach_z_p2) * (z_center_curve_inertial_p1 - z_center_curve_inertial_p2));
                        ffx = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature_p1 * R_curvature_p2 / (R_curvature_p1 + R_curvature_p2)) * (1. / (par_par_dist)-1. / particle_th_lub) * dot_pro * (x_center_curve_inertial_p1 - x_center_curve_inertial_p2) / SQ(Dij);
                        ffy = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature_p1 * R_curvature_p2 / (R_curvature_p1 + R_curvature_p2)) * (1. / (par_par_dist)-1. / particle_th_lub) * dot_pro * (y_center_curve_inertial_p1 - y_center_curve_inertial_p2) / SQ(Dij);
                        ffz = -6 * M_PI * (1.0) * ((tau - 0.5) / 3.0) * SQ(R_curvature_p1 * R_curvature_p2 / (R_curvature_p1 + R_curvature_p2)) * (1. / (par_par_dist)-1. / particle_th_lub) * dot_pro * (z_center_curve_inertial_p1 - z_center_curve_inertial_p2) / SQ(Dij);

                    }
                }
                if (par_par_dist <= particle_th && par_par_dist > 0)
                {
                    ffx = particle_cij / particle_ep * ((x_center_curve_inertial_p1 - x_center_curve_inertial_p2)) / Dij * pow((par_par_dist - particle_th) / particle_th, 2);
                    ffy = particle_cij / particle_ep * ((y_center_curve_inertial_p1 - y_center_curve_inertial_p2)) / Dij * pow((par_par_dist - particle_th) / particle_th, 2);
                    ffz = particle_cij / particle_ep * ((z_center_curve_inertial_p1 - z_center_curve_inertial_p2)) / Dij * pow((par_par_dist - particle_th) / particle_th, 2);
                }
                if (par_par_dist < 0)
                {
                    ffx = (x_center_curve_inertial_p1 - x_center_curve_inertial_p2) / Dij * (particle_cij / particle_ep * pow((par_par_dist - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (-1 * par_par_dist));
                    ffy = (y_center_curve_inertial_p1 - y_center_curve_inertial_p2) / Dij * (particle_cij / particle_ep * pow((par_par_dist - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (-1 * par_par_dist));
                    ffz = (z_center_curve_inertial_p1 - z_center_curve_inertial_p2) / Dij * (particle_cij / particle_ep * pow((par_par_dist - particle_th) / particle_th, 2) + particle_cij / particle_Ep / particle_th * (-1 * par_par_dist));
                }

                particle_list[np].center[0].f3x += ffx;
                particle_list[np].center[0].f3y += ffy;
                particle_list[np].center[0].f3z += ffz;

                particle_list[np].center[0].t3_x += ffy * (z_intersect_real_inertial_of_p1 - particle_list[np].center[0].z) - ffz * (y_intersect_real_inertial_of_p1 - particle_list[np].center[0].y);
                particle_list[np].center[0].t3_y += ffz * (x_intersect_real_inertial_of_p1 - particle_list[np].center[0].x) - ffx * (z_intersect_real_inertial_of_p1 - particle_list[np].center[0].z);
                particle_list[np].center[0].t3_z += ffx * (y_intersect_real_inertial_of_p1 - particle_list[np].center[0].y) - ffy * (x_intersect_real_inertial_of_p1 - particle_list[np].center[0].x);

            } //// if (par_par_dist<thp_control)
        }  ///// if (index_p != np)
    }// END of Particle-Particle Collision
    return;
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: UPDATE_VERLET
//-----------------------------------------------------------------------------------------------------------------------------
void ParticleResolved::update_verlet()
{

    INT32 ssu = 0;
    for (INT32 i = 0; i < num_particles; i++)
    {
        verlet_list[i].resize(0);
        for (INT32 j = 0; j < num_particles; j++)
        {
            if (j != i)
            {
                double Dij = sqrt(SQ(particle_list[i].center[0].x - particle_list[j].center[0].x) + SQ(particle_list[i].center[0].y - particle_list[j].center[0].y) + SQ(particle_list[i].center[0].z - particle_list[j].center[0].z));
                if (Dij < particle_r_buffer_verlet)
                {
                    verlet_list[i].push_back(j);
                    ssu++;
                }
            }
        }
    }
    //cout << "total "<<ssu << endl;
    return;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Copy Constructor
//-----------------------------------------------------------------------------------------------------------------------------
ParticleResolved::ParticleResolved(const ParticleResolved& other) = default; // = default is needed because I initialized all settings and parameters as initializer list in the main constructor and they all were "const"

//=============================================================================================================================
//=============================================================================================================================
ParticleResolved::~ParticleResolved()
{
    //dtor
}
