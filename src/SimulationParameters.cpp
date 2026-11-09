#include "SimulationParameters.h"
#include "SimulationSettings.h"
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Constructor
//-----------------------------------------------------------------------------------------------------------------------------
SimulationParameters::SimulationParameters(const std::string& filename, const SimulationSettings& settingsLB):
    Flow_With_Particle(settingsLB.get_bool_setting("Flow_With_Particle")),
    Flow_With_Thermal_Effect(settingsLB.get_bool_setting("Flow_With_Thermal_Effect")),
    Create_Solid(settingsLB.get_int_setting("Create_Solid")),
    Shan_Chen(settingsLB.get_bool_setting("Shan_Chen")),
    Guo(settingsLB.get_bool_setting("Guo")),
    Luo(settingsLB.get_bool_setting("Luo")),
    SRT(settingsLB.get_bool_setting("SRT")),
    MRT(settingsLB.get_bool_setting("MRT")),
    X_Periodic(settingsLB.get_bool_setting("X_Periodic")),
    X_Wall(settingsLB.get_bool_setting("X_Wall")),
    X_velIn_ZeroGradOut(settingsLB.get_bool_setting("X_velIn_ZeroGradOut")),
    X_velIn_PressureOut(settingsLB.get_bool_setting("X_velIn_PressureOut")),
    Y_Periodic(settingsLB.get_bool_setting("Y_Periodic")),
    Y_Wall(settingsLB.get_bool_setting("Y_Wall")),
    Z_Periodic(settingsLB.get_bool_setting("Z_Periodic")),
    Z_Wall(settingsLB.get_bool_setting("Z_Wall")),
    VTK_ASCII(settingsLB.get_bool_setting("VTK_ASCII")),
    VTK_BINARY(settingsLB.get_bool_setting("VTK_BINARY")),
    CASE(settingsLB.get_bool_setting("CASE")),
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
    Force_Thermal_1(settingsLB.get_bool_setting("Force_Thermal_1")),
    Force_Thermal_2(settingsLB.get_bool_setting("Force_Thermal_2")),
    Equilibrium_Thermal_1(settingsLB.get_bool_setting("Equilibrium_Thermal_1")),
    Equilibrium_Thermal_2(settingsLB.get_bool_setting("Equilibrium_Thermal_2")),
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
    Recovery_OFF(settingsLB.get_bool_setting("Recovery_OFF")),
    Recovery_ON(settingsLB.get_bool_setting("Recovery_ON")),
    ASCII_Recovery_File(settingsLB.get_bool_setting("ASCII_Recovery_File")),
    BINARY_Recovery_File(settingsLB.get_bool_setting("BINARY_Recovery_File")),
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
    Ellipse_Collision_4(settingsLB.get_bool_setting("Ellipse_Collision_4"))
{

    parse_file(filename);

    // Assign parsed values to member variables, use default values if not found
    domain_dimension_F   = Flow_Stencil.first;
    stencil_directions_F = Flow_Stencil.second;
    domain_dimension_T   = Thermal_Stencil.first;
    stencil_directions_T = Thermal_Stencil.second;

    assign_param("Nx", Nx);
    assign_param("Ny", Ny);
    assign_param("Nz", Nz);
    assign_param("tau", tau);

    assign_param("t_end", t_end);
    assign_param("t_data", t_data);
    assign_param("t_vtk", t_vtk);
    assign_param("t_info", t_info);
    assign_param("t_time", t_time);
    assign_param("t_recovery", t_recovery);

    assign_param("Stop_by_Residual", Stop_by_Residual);
    assign_param("residual_flow", residual_flow);
    assign_param("residual_thermal", residual_thermal);
    assign_param("t_residual", t_residual);

    assign_param("gravity", gravity);
    assign_param("wall_vel_bottom", wall_vel_bottom);
    assign_param("wall_vel_top", wall_vel_top);
    assign_param("wall_vel_back", wall_vel_back);
    assign_param("wall_vel_front", wall_vel_front);
    assign_param("inletvel", inletvel);

    assign_param("mdfstep", mdfstep);

    assign_param("fluid_init_temp", fluid_init_temp);
    assign_param("gbeta_x", gbeta_x);
    assign_param("gbeta_y", gbeta_y);
    assign_param("gbeta_z", gbeta_z);

    assign_param("Tref", Tref);
    assign_param("prandtl", prandtl);
    assign_param("particle_specific_heat_ratio", particle_specific_heat_ratio);
    assign_param("temperature_left", temperature_left);
    assign_param("temperature_right", temperature_right);
    assign_param("temperature_bot", temperature_bot);
    assign_param("temperature_top", temperature_top);
    assign_param("temperature_front", temperature_front);
    assign_param("temperature_back", temperature_back);
    assign_param("heat_flux_left", heat_flux_left);
    assign_param("heat_flux_right", heat_flux_right);
    assign_param("heat_flux_bottom", heat_flux_bottom);
    assign_param("heat_flux_top", heat_flux_top);
    assign_param("heat_flux_front", heat_flux_front);
    assign_param("heat_flux_back", heat_flux_back);

    assign_param("num_particles", num_particles);
    assign_param("particle_radius", particle_radius);
    assign_param("particle_center_x", particle_center_x);
    assign_param("particle_center_y", particle_center_y);
    assign_param("particle_center_z", particle_center_z);
    assign_param("particle_gravity_x", particle_gravity_x);
    assign_param("particle_gravity_y", particle_gravity_y);
    assign_param("particle_gravity_z", particle_gravity_z);
    assign_param("particle_den_ratio", particle_den_ratio);
    assign_param("particle_offset_x", particle_offset_x);
    assign_param("particle_offset_y", particle_offset_y);
    assign_param("particle_offset_z", particle_offset_z);
    assign_param("particle_temperature", particle_temperature);
    assign_param("particle_heat_source", particle_heat_source);

    assign_param("particle_radius_2", particle_radius_2);
    assign_param("particle_num_nodes", particle_num_nodes);
    assign_param("particle_num_strip", particle_num_strip);
    assign_param("particle_cir_factor", particle_cir_factor);
    assign_param("ellipse_init_ang", ellipse_init_ang);

    assign_param("particle_r_buffer_verlet", particle_r_buffer_verlet);
    assign_param("particle_r_cutoff_verlet", particle_r_cutoff_verlet);
    assign_param("particle_cij", particle_cij);
    assign_param("particle_ep", particle_ep);
    assign_param("particle_Ep", particle_Ep);
    assign_param("particle_epw", particle_epw);
    assign_param("particle_Epw", particle_Epw);
    assign_param("particle_th", particle_th);
    assign_param("particle_thw", particle_thw);
    assign_param("particle_th_lub", particle_th_lub);
    assign_param("particle_thw_lub", particle_thw_lub);

    assign_param("num_buffer_layers_X", initial_num_buffer_layers_X);
    assign_param("num_buffer_layers_Y", initial_num_buffer_layers_Y);
    assign_param("num_buffer_layers_Z", initial_num_buffer_layers_Z);

    assign_param("solid_temperature", solid_temperature);
    assign_param("solid_to_fluid_conductivity", solid_to_fluid_conductivity);

    //assign_param("num_processors_for_particles_only", num_processors_for_particles_only);

    ///Handling errors
    if (domain_dimension_F == 2 && Nz != 1) {
        throw std::runtime_error("Error =======> Nz must be 1 for a 2D domain. <=======");
    }
    if (domain_dimension_F != domain_dimension_T) {
        throw std::runtime_error("Error =======> domain_dimension_F must be equal to domain_dimension_T. <=======");
    }
    if (domain_dimension_F == 2 && !Z_Periodic) {
        throw std::runtime_error("Error =======> Z direction must be periodic BC for a 2D domain. <=======");
    }
    if (domain_dimension_F == 2 && abs(particle_center_z - 0.0) > 1e-6) {
        throw std::runtime_error("Error =======> particle_center_z must be 0 for a 2D domain. <=======");
    }

   /* Nx = static_cast<INT32>(parameters["Nx"]);
    Ny = static_cast<INT32>(parameters["Ny"]);
    Nz = static_cast<INT32>(parameters["Nz"]);
    tau = parameters["tau"];
    t_end = static_cast<INT32>(parameters["t_end"]);
    t_data = static_cast<INT32>(parameters["t_data"]);
    t_vtk = static_cast<INT32>(parameters["t_vtk"]);
    t_info = static_cast<INT32>(parameters["t_info"]);
    t_time = static_cast<INT32>(parameters["t_time"]);
    t_recovery = static_cast<INT32>(parameters["t_recovery"]);
    gravity = parameters["gravity"];
    wall_vel_bottom = parameters["wall_vel_bottom"];
    wall_vel_top = parameters["wall_vel_top"];
    wall_vel_back = parameters["wall_vel_back"];
    wall_vel_front = parameters["wall_vel_front"];
    inletvel = parameters["inletvel"];

    mdfstep = static_cast<INT32>(parameters["mdfstep"]);
    etc .....      */

    // Initialize the sortedParameters map with the contents of the parameters map
    sortedParameters.insert(parameters.begin(), parameters.end());

    omega = 1.0 / tau; // relaxation frequency (inverse of relaxation time)
    if (Flow_With_Thermal_Effect)
    {
        alpha_t = (tau - 0.5) / 3. / prandtl; // Thermal diffusion coefficient
        tau_t = 3. * alpha_t + 0.5; // relaxation time of temperature
        omega_t = 1.0 / tau_t; // relaxation frequency (inverse of relaxation time for temperature)
    }

    if (!Flow_With_Thermal_Effect)  // Do NOT change this part
    {
        /// Do NOT change below lines
        fluid_init_temp = 0.0;
        gbeta_x = 0.0;
        gbeta_y = 0.0;
        gbeta_z = 0.0;
        Tref = 0.0; //reference temperature in case of natural convection
        particle_temperature = 0.0; //Particle temperature must be constant
    }

    if (!Flow_With_Particle){
        omega_ts = 1. / (3.* solid_to_fluid_conductivity * alpha_t + 0.5);//for solid nodes
    }

    /// Compute derived quantities
    double Dia = Ny - 2; // inner channel diameter
    nu = (tau - 0.5) / 3.0; // lattice viscosity
    double umax = gravity / (2.0 * nu) * SQ(0.5 * Dia); // expected maximum velocity for Poiseuille flow without immersed object
    double Re{};
    if (gravity != 0)
    {
        Re = Dia * 2.0 / 3.0 * umax / nu; // Reynolds number for Poiseuille flow without immersed object.
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: assign_param
//-----------------------------------------------------------------------------------------------------------------------------
template<typename T> void SimulationParameters::assign_param(const std::string& key, T& var) {
    auto it = parameters.find(key);
    if (it != parameters.end()) {
        var = static_cast<T>(it->second);
    } else {
        std::cerr << "Error: parameter '" << key << "' not found in parametersLB file!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: load_parameters
//-----------------------------------------------------------------------------------------------------------------------------
SimulationParameters& SimulationParameters::get_instance(const std::string& filename, const SimulationSettings& settingsLB)
{
    static SimulationParameters instance(filename, settingsLB);
    return instance;
}


//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: parse_file
//-----------------------------------------------------------------------------------------------------------------------------
void SimulationParameters::parse_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open/find file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Remove CR from Windows line endings
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        if (line.empty() || line[0] == '#') continue; // Skip comments and empty lines
        parse_line(line);
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: parse_line
//-----------------------------------------------------------------------------------------------------------------------------
void SimulationParameters::parse_line(const std::string &line) {
    std::regex rgx(R"((\S+)\s*=\s*(\S+))");
    std::smatch match;
    std::string::const_iterator searchStart(line.cbegin());

    while (std::regex_search(searchStart, line.cend(), match, rgx)) {
        std::string key = match[1];
        std::string value = match[2];

        if (key == "Flow_Stencil") {
            parse_stencil(value, Flow_Stencil);
            stencil_flow = value;
         //   std::cout<< "Stencil of Flow is: "<< stencil_flow <<std::endl;
        } else if (key == "Thermal_Stencil") {
            parse_stencil(value, Thermal_Stencil);
            stencil_thermal = value;
         //   std::cout<< "Stencil of Thermal is: "<< stencil_thermal <<std::endl;
        } else {
            if (parameters.find(key) == parameters.end()) {
                insertion_order.push_back(key);
            }
            parameters[key] = std::stod(value);
        }

        searchStart = match.suffix().first;
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: parse_stencil
//-----------------------------------------------------------------------------------------------------------------------------
void SimulationParameters::parse_stencil(const std::string &stencil, std::pair<INT32, INT32> &result) {
    std::regex stencil_rgx(R"(D(\d+)Q(\d+))");
    std::smatch stencil_match;
    if (std::regex_match(stencil, stencil_match, stencil_rgx)) {
        result.first = std::stoi(stencil_match[1]);
        result.second = std::stoi(stencil_match[2]);
    } else {
        std::cerr << "Invalid stencil format: " << stencil << std::endl;
    }
}

//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: report_parameters_on_screen
//-----------------------------------------------------------------------------------------------------------------------------
void SimulationParameters::report_parameters_on_screen() {

    /*for (const auto& pair : paramsLB.sortedParameters)    // Iterate over the sorted map and print the key-value pairs
    {
        std::cout << "Key: " << pair.first << " = " << pair.second << std::endl;
    }*/

    /// Print parameters
    std::cout << "Flow_Stencil = " << "D" << domain_dimension_F << "Q" << stencil_directions_F << std::endl;
    std::cout << "Thermal_Stencil = " << "D" << domain_dimension_T << "Q" << stencil_directions_T << std::endl;
    INT32 count = 1;
    for (const auto& key : insertion_order) {
        std::cout << count << "- " << key << " = " << parameters[key] << '\n';
        count++;
    }
}
//=============================================================================================================================
//=============================================================================================================================
SimulationParameters::~SimulationParameters()
{
    //dtor
}
