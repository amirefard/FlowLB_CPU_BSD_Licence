#ifndef FLOW_SOLVER_H
#define FLOW_SOLVER_H

#include <cmath>
#include <sstream>  // string streams
#include <iostream> // for the use of 'cout'
#include <vector>
#include <string>
using namespace std;
#if defined __linux__ || defined __APPLE__
#include <sys/stat.h>//for mkdir
#else
#include "direct.h"//for mkdir
#endif


#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"
#include "Particle_sim.h"


class Flow_solver
{
private:
    int az;
    vector<int> cx_pop_F, cy_pop_F;
    vector<double> weight_F;//Thermal weights
    vector<int> bounce_back_pop_F;//bounce back direction of each direction (for thermal)
    vector<int> pop_to_top_F, pop_to_bottom_F, pop_to_left_F, pop_to_right_F;
    vector<int> pop_reverse_top_F, pop_reverse_bottom_F, pop_reverse_left_F, pop_reverse_right_F;

public:

    struct solid_node {
        /// Constructor
        solid_node() {
            x = 0;
            y = 0;
        }
        /// Elements
        double x; // current x-position
        double y; // current y-position
    };

    struct solid_surface {
        /// Constructor
        solid_surface() {
            num_surf_nodes = surface_num_nodes;
            node_surf = new solid_node[num_surf_nodes];
        }

        /// Elements
        int num_surf_nodes; // number of surface nodes on solid structure
        solid_node *node_surf; // center node
    };

    solid_surface solid; // create solid surface tracer

    Flow_solver();
    void initialize();
    void equilibrium(double, double, double); // compute the equilibrium populations from the fluid density and velocity
    #if defined Flow_With_Resolved_Particle
        void LBM(int, const Particle_sim&); // perform LBM operations
    #endif
    #if defined Flow_Without_Particle
        void LBM(int);
    #endif
    void BC(int); // perform boundary conditions operations
    void BC_OLD_VERSION(int);
    void momenta(); // compute fluid density and velocity from the populations
    double meq(double, double, double, int);
    bool check_residual();
    #if defined Flow_With_Resolved_Particle
        int calculate_gbetaindex(int X, int Y, const Particle_sim&);
    #endif
    ~Flow_solver();
};

#endif
