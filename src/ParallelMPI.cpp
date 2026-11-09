#include "ParallelMPI.h"
#include "SimulationSettings.h"
#include "SimulationParameters.h"
//=============================================================================================================================
//=============================================================================================================================
ParallelMPI::ParallelMPI(int argc, char** argv)
{
    //ctor
    MPI_Init(&argc, &argv);

    // Get the number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors_total);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &this->processor_id);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    INT32 name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    //printf("Hello world from processor %s, rank %d out of %d processors\n",  processor_name, processor_id, this->num_processors_total);
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: prepareMpiSetup
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::prepareMpiSetup(const SimulationParameters& params, const SimulationSettings& settings)
{
    set_simulation_parameters(params);

    if (this->num_processors_total == 1)
    {
        this->num_processors_for_particles_only = 0;//Do NOT change. Because one processor should handle both the flow and particles
        if (processor_id == MASTER)
        {
            printf("\nNOTE: num_processors_for_particles_only changed to 0 because num_processors_total = 1.\n");
        }
    }
    is_processor_for_particle_only.resize(num_processors_total, false);//initialize
    if (this->processor_id < this->num_processors_for_particles_only)
    {
        is_processor_for_particle_only[this->processor_id] = true;
    }

    start_global_x.resize(this->num_processors_total);
    end_global_x.resize(this->num_processors_total);
    start_global_y.resize(this->num_processors_total);
    end_global_y.resize(this->num_processors_total);
    start_global_z.resize(this->num_processors_total);
    end_global_z.resize(this->num_processors_total);

    start_local_x.resize(this->num_processors_total);
    end_local_x.resize(this->num_processors_total);
    start_local_y.resize(this->num_processors_total);
    end_local_y.resize(this->num_processors_total);
    start_local_z.resize(this->num_processors_total);
    end_local_z.resize(this->num_processors_total);

    num_buffer_layers_X.resize(this->num_processors_total);
    num_buffer_layers_Y.resize(this->num_processors_total);
    num_buffer_layers_Z.resize(this->num_processors_total);

    length_with_buffer_X.resize(this->num_processors_total);
    length_with_buffer_Y.resize(this->num_processors_total);
    length_with_buffer_Z.resize(this->num_processors_total);

    local_neighbors_left.resize(this->num_processors_total);
    local_neighbors_right.resize(this->num_processors_total);
    local_neighbors_bottom.resize(this->num_processors_total);
    local_neighbors_top.resize(this->num_processors_total);
    local_neighbors_front.resize(this->num_processors_total);
    local_neighbors_back.resize(this->num_processors_total);

    if (this->processor_id == MASTER)
    {
        findOptimalProcessorDistribution(this->num_processors_total - this->num_processors_for_particles_only);// I send num_processors_total - 1 because Master is excluded. MASTER will only handle particles unless num_processors_total is 1.
        setupMpiDomain();//find start_x, y, z and end_x, y, z. Also find neighbors
    }

    /// Broadcast setup details to all processors
    MPI_Bcast(&num_processors_X, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&num_processors_Y, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&num_processors_Z, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Broadcast number of buffer layers
    MPI_Bcast(num_buffer_layers_X.data(), num_buffer_layers_X.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(num_buffer_layers_Y.data(), num_buffer_layers_Y.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(num_buffer_layers_Z.data(), num_buffer_layers_Z.size(), MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Broadcast start_x, end_x, start_y, end_y, start_z, end_z, and local_neighbors_left
    MPI_Bcast(start_global_x.data(), start_global_x.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(end_global_x.data(), end_global_x.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(start_global_y.data(), start_global_y.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(end_global_y.data(), end_global_y.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(start_global_z.data(), start_global_z.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(end_global_z.data(), end_global_z.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(start_local_x.data(), start_local_x.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(end_local_x.data(), end_local_x.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(start_local_y.data(), start_local_y.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(end_local_y.data(), end_local_y.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(start_local_z.data(), start_local_z.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(end_local_z.data(), end_local_z.size(), MPI_INT, MASTER, MPI_COMM_WORLD);

    MPI_Bcast(length_with_buffer_X.data(), length_with_buffer_X.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(length_with_buffer_Y.data(), length_with_buffer_Y.size(), MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(length_with_buffer_Z.data(), length_with_buffer_Z.size(), MPI_INT, MASTER, MPI_COMM_WORLD);

    /// Broadcast local_neighbors_left, local_neighbors_right, local_neighbors_bottom, local_neighbors_top, local_neighbors_front, local_neighbors_back
    broadcast2DVector(local_neighbors_left);
    broadcast2DVector(local_neighbors_right);
    broadcast2DVector(local_neighbors_bottom);
    broadcast2DVector(local_neighbors_top);
    broadcast2DVector(local_neighbors_front);
    broadcast2DVector(local_neighbors_back);

    if (this->processor_id == MASTER)
    {
        cout << "\n================= Printing MPI data by Processor 0 ================="<<endl;
        printAllNeighbors();
    }

}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: broadcast2DVector
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::broadcast2DVector(std::vector<std::vector<INT32>>& vec) {
    INT32 num_elements = vec.size();
    MPI_Bcast(&num_elements, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    for (INT32 i = 0; i < num_elements; ++i) {
        INT32 subvector_size = vec[i].size();
        MPI_Bcast(&subvector_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        vec[i].resize(subvector_size);
        if (subvector_size > 0) {
            MPI_Bcast(vec[i].data(), subvector_size, MPI_INT, MASTER, MPI_COMM_WORLD);
        }
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: findOptimalProcessorDistribution
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::findOptimalProcessorDistribution(INT32 total_processors) {
    if (total_processors > 0) {
        INT32 best_px = total_processors, best_py = 1, best_pz = 1;
        double best_score = std::numeric_limits<double>::max(); // Initialize to maximum value

        for (INT32 px = total_processors; px >= 1; --px) { // Start with larger px
            if (total_processors % px != 0) continue;
            INT32 remaining_processors = total_processors / px;
            for (INT32 py = remaining_processors; py >= 1; --py) { // Start with larger py
                if (remaining_processors % py != 0) continue;
                INT32 pz = remaining_processors / py;

                // Calculate the dimensions of the subdomains
                double subdomain_x = static_cast<double>(this->Nx) / px;
                double subdomain_y = static_cast<double>(this->Ny) / py;
                double subdomain_z = static_cast<double>(this->Nz) / pz;

                // Score based on how close the subdomains are to being cubic (or balanced)
                double score = computeScore(subdomain_x, subdomain_y, subdomain_z);

                // Prefer configurations with smaller score (closer to cubic)
                if (score < best_score) {
                    best_score = score;
                    best_px = px;
                    best_py = py;
                    best_pz = pz;
                }
            }
        }

        // Store the values directly in class private variables
        num_processors_X = best_px;
        num_processors_Y = best_py;
        num_processors_Z = best_pz;
    } else {//MASTER will handle both particles and fluid
        num_processors_X = 1;
        num_processors_Y = 1;
        num_processors_Z = 1;
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: computeScore
//-----------------------------------------------------------------------------------------------------------------------------
double ParallelMPI::computeScore(double subdomain_x, double subdomain_y, double subdomain_z) {
    double score;
    //Method 1: Calculate the deviation from being cubic
    score = std::max({subdomain_x / subdomain_y, subdomain_y / subdomain_x}) *
            std::max({subdomain_x / subdomain_z, subdomain_z / subdomain_x}) *
            std::max({subdomain_y / subdomain_z, subdomain_z / subdomain_y});

    return score;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: set_simulation_parameters
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::set_simulation_parameters(const SimulationParameters& params)
{
    this->Nx = params.get_Nx();
    this->Ny = params.get_Ny();
    this->Nz = params.get_Nz();

    this->initial_num_buffer_layers_X = params.get_initial_num_buffer_layers_X();
    this->initial_num_buffer_layers_Y = params.get_initial_num_buffer_layers_Y();
    this->initial_num_buffer_layers_Z = params.get_initial_num_buffer_layers_Z();

    /// At this version only one processor (MASTER) handles the particles
    this->num_processors_for_particles_only = 1;// params.get_num_processors_for_particles_only();// This parameter may change later in prepareMpiSetup
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: setupMpiDomain
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::setupMpiDomain()
{
    // Distribute domain among processors
    INT32 x_size = this->Nx / num_processors_X;
    INT32 y_size = this->Ny / num_processors_Y;
    INT32 z_size = this->Nz / num_processors_Z;

    INT32 x_remainder = this->Nx % num_processors_X;
    INT32 y_remainder = this->Ny % num_processors_Y;
    INT32 z_remainder = this->Nz % num_processors_Z;

    INT32 idx{};

    /// Find the start and end of all processors that are ONLY for particles.
    for (INT32 i = 0; i < this->num_processors_for_particles_only; i++){// If num_processors_for_particles_only = 0 then this loop is not executed
        num_buffer_layers_X[i] = initial_num_buffer_layers_X;
        num_buffer_layers_Y[i] = initial_num_buffer_layers_Y;
        num_buffer_layers_Z[i] = initial_num_buffer_layers_Z;

        // The Master (id = 0) has the whole domain length
        start_global_x[i] = 0;
        end_global_x[i]   = this->Nx - 1;
        start_global_y[i] = 0;
        end_global_y[i]   = this->Ny - 1;
        start_global_z[i] = 0;
        end_global_z[i]   = this->Nz - 1;
        start_local_x[i]  = 0;//Always 0
        end_local_x[i]    = end_global_x[i] - start_global_x[i] + 2 * num_buffer_layers_X[i];
        start_local_y[i]  = 0;//Always 0
        end_local_y[i]    = end_global_y[i] - start_global_y[i] + 2 * num_buffer_layers_Y[i];
        start_local_z[i]  = 0;//Always 0
        end_local_z[i]    = end_global_z[i] - start_global_z[i] + 2 * num_buffer_layers_Z[i];
    }

    /// Find the start and end of all processors that are not for particles.
    idx = this->num_processors_for_particles_only;// Start from idx=num_processors_for_particles (thus excluding the master processor (processor 0) that deals only with particles)
    for (INT32 px = 0; px < num_processors_X; ++px) {
        for (INT32 py = 0; py < num_processors_Y; ++py) {
            for (INT32 pz = 0; pz < num_processors_Z; ++pz) {

                num_buffer_layers_X[idx] = initial_num_buffer_layers_X;// Non-zero number of buffer layers for flow processors
                num_buffer_layers_Y[idx] = initial_num_buffer_layers_Y;// Non-zero number of buffer layers for flow processors
                num_buffer_layers_Z[idx] = initial_num_buffer_layers_Z;// Non-zero number of buffer layers for flow processors

                INT32 xs = px * x_size + std::min(px, x_remainder);
                INT32 xe = xs + x_size + (px < x_remainder ? 1 : 0) - 1;

                INT32 ys = py * y_size + std::min(py, y_remainder);
                INT32 ye = ys + y_size + (py < y_remainder ? 1 : 0) - 1;

                INT32 zs = pz * z_size + std::min(pz, z_remainder);
                INT32 ze = zs + z_size + (pz < z_remainder ? 1 : 0) - 1;

                start_global_x[idx] = xs;
                end_global_x[idx]   = xe;
                start_global_y[idx] = ys;
                end_global_y[idx]   = ye;
                start_global_z[idx] = zs;
                end_global_z[idx]   = ze;

                start_local_x[idx]  = 0;//Always 0
                end_local_x[idx]    = end_global_x[idx] - start_global_x[idx] + 2 * num_buffer_layers_X[idx];
                start_local_y[idx]  = 0;//Always 0
                end_local_y[idx]    = end_global_y[idx] - start_global_y[idx] + 2 * num_buffer_layers_Y[idx];
                start_local_z[idx]  = 0;//Always 0
                end_local_z[idx]    = end_global_z[idx] - start_global_z[idx] + 2 * num_buffer_layers_Z[idx];

                // Calculate neighboring processor IDs with periodic BC
               /* INT32 left_neighbor = (px == 0) ? idx + num_processors_X - 1 : idx - 1;
                INT32 right_neighbor = (px == num_processors_X - 1) ? idx - num_processors_X + 1 : idx + 1;
                INT32 bottom_neighbor = (py == 0) ? idx + num_processors_X * (num_processors_Y - 1) : idx - num_processors_X;
                INT32 top_neighbor = (py == num_processors_Y - 1) ? idx - num_processors_X * (num_processors_Y - 1) : idx + num_processors_X;
                INT32 front_neighbor = (pz == 0) ? idx + num_processors_X * num_processors_Y * (num_processors_Z - 1) : idx - num_processors_X * num_processors_Y;
                INT32 back_neighbor = (pz == num_processors_Z - 1) ? idx - num_processors_X * num_processors_Y * (num_processors_Z - 1) : idx + num_processors_X * num_processors_Y;*/

                INT32 left_neighbor  = ((px - 1 + num_processors_X) % num_processors_X) * num_processors_Y * num_processors_Z + py * num_processors_Z + pz + num_processors_for_particles_only;
                INT32 right_neighbor = ((px + 1) % num_processors_X) * num_processors_Y * num_processors_Z + py * num_processors_Z + pz + num_processors_for_particles_only;
                INT32 bottom_neighbor = px * num_processors_Y * num_processors_Z + ((py - 1 + num_processors_Y) % num_processors_Y) * num_processors_Z + pz + num_processors_for_particles_only;
                INT32 top_neighbor = px * num_processors_Y * num_processors_Z + ((py + 1) % num_processors_Y) * num_processors_Z + pz + num_processors_for_particles_only;
                INT32 front_neighbor = px * num_processors_Y * num_processors_Z + py * num_processors_Z + ((pz - 1 + num_processors_Z) % num_processors_Z) + num_processors_for_particles_only;
                INT32 back_neighbor = px * num_processors_Y * num_processors_Z + py * num_processors_Z + ((pz + 1) % num_processors_Z) + num_processors_for_particles_only;

             //   if (left_neighbor != idx)
                    {local_neighbors_left[idx].push_back(left_neighbor);}
             //   if (right_neighbor != idx)
                {local_neighbors_right[idx].push_back(right_neighbor);}
             //   if (bottom_neighbor != idx)
                 {local_neighbors_bottom[idx].push_back(bottom_neighbor);}
             //   if (top_neighbor != idx)
                {local_neighbors_top[idx].push_back(top_neighbor);}
              //  if (front_neighbor != idx)
                 {local_neighbors_front[idx].push_back(front_neighbor);}
             //   if (back_neighbor != idx)
                  {local_neighbors_back[idx].push_back(back_neighbor);}

                idx++;//go to the next processor
            }
        }
    }

    for (INT32 id = 0; id < this->num_processors_total; ++id) {
        length_with_buffer_X[id] = end_global_x[id] - start_global_x[id] + 1 + 2 * num_buffer_layers_X[id];
        length_with_buffer_Y[id] = end_global_y[id] - start_global_y[id] + 1 + 2 * num_buffer_layers_Y[id];
        length_with_buffer_Z[id] = end_global_z[id] - start_global_z[id] + 1 + 2 * num_buffer_layers_Z[id];
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: printNeighborsForProcessor
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::printNeighborsForProcessor(INT32 id) {
    std::stringstream ss;
    if (id == 0){//To avoid printing this line multiple times
        ss << "num_processors_X: " << num_processors_X <<", num_processors_Y: " << num_processors_Y << ", num_processors_Z: "<< num_processors_Z << std::endl;
        ss << "===================================================================="<<endl;
    }
    ss << "Processor: " << id << "\n";
    ss << "  is_processor_for_particle_only: " << (is_processor_for_particle_only[id] ? "true" : "false") << "\n";
    ss << "  start_global_x: " << start_global_x[id] << ", end_global_x: " << end_global_x[id] << "\n";
    ss << "  start_global_y: " << start_global_y[id] << ", end_global_y: " << end_global_y[id] << "\n";
    ss << "  start_global_z: " << start_global_z[id] << ", end_global_z: " << end_global_z[id] << "\n";
    ss << "  start_local_x:  " << start_local_x[id]  << ", end_local_x:  " << end_local_x[id] << "\n";
    ss << "  start_local_y:  " << start_local_y[id]  << ", end_local_y:  " << end_local_y[id] << "\n";
    ss << "  start_local_z:  " << start_local_z[id]  << ", end_local_z:  " << end_local_z[id] << "\n";


    ss << "Neighbors:\n";
    /// Left neighbors
    ss << "  Left: ";
    if (local_neighbors_left[id].size() > 0) {
        for (INT32 i = 0; i < local_neighbors_left[id].size(); ++i) {
            ss << local_neighbors_left[id][i] << " ";
        }
    }
    ss << "\n";

    /// Right neighbors
    ss << "  Right: ";
    if (local_neighbors_right[id].size() > 0) {
        for (INT32 i = 0; i < local_neighbors_right[id].size(); ++i) {
            ss << local_neighbors_right[id][i] << " ";
        }
    }
    ss << "\n";

    /// Bottom neighbors
    ss << "  Bottom: ";
    if (local_neighbors_bottom[id].size() > 0) {
        for (INT32 i = 0; i < local_neighbors_bottom[id].size(); ++i) {
            ss << local_neighbors_bottom[id][i] << " ";
        }
    }
    ss << "\n";

    /// Top neighbors
    ss << "  Top: ";
    if (local_neighbors_top[id].size() > 0) {
        for (INT32 i = 0; i < local_neighbors_top[id].size(); ++i) {
            ss << local_neighbors_top[id][i] << " ";
        }
    }
    ss << "\n";

    /// Front neighbors
    ss << "  Front: ";
    if (local_neighbors_front[id].size() > 0) {
        for (INT32 i = 0; i < local_neighbors_front[id].size(); ++i) {
            ss << local_neighbors_front[id][i] << " ";
        }
    }
    ss << "\n";

    /// Back neighbors
    ss << "  Back: ";
    if (local_neighbors_back[id].size() > 0) {
        for (INT32 i = 0; i < local_neighbors_back[id].size(); ++i) {
            ss << local_neighbors_back[id][i] << " ";
        }
    }
    ss << "\n============\n";

    std::cout << ss.str();
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: printAllNeighbors
//-----------------------------------------------------------------------------------------------------------------------------
void ParallelMPI::printAllNeighbors() {
    for (INT32 id = 0; id < this->num_processors_total; ++id) { // Start from 1 to exclude the master processor (processor 0)
        printNeighborsForProcessor(id);
    }
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_owner_processor_for_cell
//-----------------------------------------------------------------------------------------------------------------------------
INT32 ParallelMPI::get_owner_processor_for_cell(INT32 X, INT32 Y, INT32 Z) const {// Important: X, Y, Z are Global (not local) positions. So can change between 0_Nx-1, 0_Ny-1, 0_Nz-1. Only checks within fluid nodes.
    for (INT32 p = this->num_processors_for_particles_only; p < num_processors_total; ++p) {
        if (X >= start_global_x[p] && X <= end_global_x[p] &&
            Y >= start_global_y[p] && Y <= end_global_y[p] &&
            Z >= start_global_z[p] && Z <= end_global_z[p]) {
            return p;
        }
    }
    std::cout << "Could not find owner processor." <<std:: endl;
    return -1; // Should not happen
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_Xglobal_from_Xlocal
//-----------------------------------------------------------------------------------------------------------------------------
///Convert Xlocal to Xglobal
INT32 ParallelMPI::get_Xglobal_from_Xlocal(INT32 Xlocal, INT32 processor_id) const
{
    INT32 Xglobal{};
    Xglobal = Xlocal + start_global_x[processor_id] - num_buffer_layers_X[processor_id];
    return Xglobal;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_Yglobal_from_Ylocal
//-----------------------------------------------------------------------------------------------------------------------------
///Convert Ylocal to Yglobal
INT32 ParallelMPI::get_Yglobal_from_Ylocal(INT32 Ylocal, INT32 processor_id) const
{
    INT32 Yglobal{};
    Yglobal = Ylocal + start_global_y[processor_id] - num_buffer_layers_Y[processor_id];
    return Yglobal;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_Zglobal_from_Zlocal
//-----------------------------------------------------------------------------------------------------------------------------
///Convert Zlocal to Zglobal
INT32 ParallelMPI::get_Zglobal_from_Zlocal(INT32 Zlocal, INT32 processor_id) const
{
    INT32 Zglobal{};
    Zglobal = Zlocal + start_global_z[processor_id] - num_buffer_layers_Z[processor_id];
    return Zglobal;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_Xlocal_from_Xglobal
//-----------------------------------------------------------------------------------------------------------------------------
///Convert Xglobal to Xlocal
INT32 ParallelMPI::get_Xlocal_from_Xglobal(INT32 Xglobal, INT32 processor_id) const
{
    INT32 Xlocal{};
    Xlocal = Xglobal - start_global_x[processor_id] + num_buffer_layers_X[processor_id];
    return Xlocal;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_Ylocal_from_Yglobal
//-----------------------------------------------------------------------------------------------------------------------------
///Convert Xglobal to Xlocal
INT32 ParallelMPI::get_Ylocal_from_Yglobal(INT32 Yglobal, INT32 processor_id) const
{
    INT32 Ylocal{};
    Ylocal = Yglobal - start_global_y[processor_id] + num_buffer_layers_Y[processor_id];
    return Ylocal;
}
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: get_Zlocal_from_Zglobal
//-----------------------------------------------------------------------------------------------------------------------------
///Convert Zglobal to Zlocal
INT32 ParallelMPI::get_Zlocal_from_Zglobal(INT32 Zglobal, INT32 processor_id) const
{
    INT32 Zlocal{};
    Zlocal = Zglobal - start_global_z[processor_id] + num_buffer_layers_Z[processor_id];
    return Zlocal;
}
//=============================================================================================================================
//=============================================================================================================================
ParallelMPI::~ParallelMPI()
{
    //dtor
    // Finalize MPI
    MPI_Finalize();
}
