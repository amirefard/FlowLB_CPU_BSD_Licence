#ifndef PARALLELMPI_H
#define PARALLELMPI_H

#include <mpi.h>
#include <iostream>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;

#include "GlobalFunctions.h"

//Forward declaration
class SimulationParameters;
class SimulationSettings;

//=============================================================================================================================
class ParallelMPI
{
public:
    ParallelMPI(int argc, char** argv);
    virtual ~ParallelMPI();

    void prepareMpiSetup(const SimulationParameters&, const SimulationSettings&);
    void exchangeData()
    {
        // Implement communication routines here
        // For example, exchanging boundary data between neighboring processors
    }
    INT32 get_processor_id() const {return processor_id;}
    INT32 get_num_processors_for_particles_only() const {return num_processors_for_particles_only;}
    INT32 get_num_processors_total() const {return num_processors_total;}
    std::vector<INT32> get_num_buffer_layers_X() const {return num_buffer_layers_X;}
    std::vector<INT32> get_num_buffer_layers_Y() const {return num_buffer_layers_Y;}
    std::vector<INT32> get_num_buffer_layers_Z() const {return num_buffer_layers_Z;}
    //INT32 get_num_processors_for_particles_only() const {return num_processors_for_particles_only;}
    std::vector<bool> get_is_processor_for_particle_only() const {return is_processor_for_particle_only;}

    std::vector<INT32> get_start_global_x() const {return start_global_x;}
    std::vector<INT32> get_end_global_x() const {return end_global_x;}
    std::vector<INT32> get_start_global_y() const {return start_global_y;}
    std::vector<INT32> get_end_global_y() const {return end_global_y;}
    std::vector<INT32> get_start_global_z() const {return start_global_z;}
    std::vector<INT32> get_end_global_z() const {return end_global_z;}
    std::vector<INT32> get_start_local_x() const {return start_local_x;}
    std::vector<INT32> get_end_local_x() const {return end_local_x;}
    std::vector<INT32> get_start_local_y() const {return start_local_y;}
    std::vector<INT32> get_end_local_y() const {return end_local_y;}
    std::vector<INT32> get_start_local_z() const {return start_local_z;}
    std::vector<INT32> get_end_local_z() const {return end_local_z;}
    std::vector<INT32> get_length_with_buffer_X() const {return length_with_buffer_X;}
    std::vector<INT32> get_length_with_buffer_Y() const {return length_with_buffer_Y;}
    std::vector<INT32> get_length_with_buffer_Z() const {return length_with_buffer_Z;}
    INT32 get_Xglobal_from_Xlocal(INT32 Xlocal, INT32 id) const;
    INT32 get_Yglobal_from_Ylocal(INT32 Ylocal, INT32 id) const;
    INT32 get_Zglobal_from_Zlocal(INT32 Zlocal, INT32 id) const;
    INT32 get_Xlocal_from_Xglobal(INT32 Xglobal, INT32 id) const;
    INT32 get_Ylocal_from_Yglobal(INT32 Yglobal, INT32 id) const;
    INT32 get_Zlocal_from_Zglobal(INT32 Zglobal, INT32 id) const;


    INT32 get_owner_processor_for_cell(INT32 X, INT32 Y, INT32 Z) const;

    std::vector<std::vector<INT32>> get_local_neighbors_left() const {return local_neighbors_left;}
    std::vector<std::vector<INT32>> get_local_neighbors_right() const {return local_neighbors_right;}
    std::vector<std::vector<INT32>> get_local_neighbors_bottom() const {return local_neighbors_bottom;}
    std::vector<std::vector<INT32>> get_local_neighbors_top() const {return local_neighbors_top;}
    std::vector<std::vector<INT32>> get_local_neighbors_front() const {return local_neighbors_front;}
    std::vector<std::vector<INT32>> get_local_neighbors_back() const {return local_neighbors_back;}

protected:

private:
    ///Parameters
    INT32 Nx{};
    INT32 Ny{};
    INT32 Nz{};
    INT32 initial_num_buffer_layers_X{};// Different from num_buffer_layers_X[]
    INT32 initial_num_buffer_layers_Y{};
    INT32 initial_num_buffer_layers_Z{};
    INT32 num_processors_for_particles_only{};//To only handle the particles. This may change in prepareMpiSetup

    ///For this class only
    INT32 num_processors_total{};
    INT32 num_processors_X{}, num_processors_Y{}, num_processors_Z{};
    INT32 processor_id{};
    std::vector<bool> is_processor_for_particle_only;
    std::vector<INT32> num_buffer_layers_X, num_buffer_layers_Y, num_buffer_layers_Z;
    std::vector<INT32> length_with_buffer_X, length_with_buffer_Y, length_with_buffer_Z;// Total length in each direction including buffer layers
    std::vector<INT32> start_global_x, end_global_x, start_global_y, end_global_y, start_global_z, end_global_z;
    std::vector<INT32> start_local_x, end_local_x, start_local_y, end_local_y, start_local_z, end_local_z;
    ///list of neighbors
    std::vector<std::vector<INT32>> local_neighbors_left;
    std::vector<std::vector<INT32>> local_neighbors_right;
    std::vector<std::vector<INT32>> local_neighbors_bottom;
    std::vector<std::vector<INT32>> local_neighbors_top;
    std::vector<std::vector<INT32>> local_neighbors_front;
    std::vector<std::vector<INT32>> local_neighbors_back;

    void set_simulation_parameters(const SimulationParameters&);
    void findOptimalProcessorDistribution(INT32);
    void setupMpiDomain();
    double computeScore(double, double, double);
    void printNeighborsForProcessor(INT32);
    void printAllNeighbors();
    void broadcast2DVector(std::vector<std::vector<INT32>>& vec);

};

#endif // PARALLELMPI_H
