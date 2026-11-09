#ifndef FLOW_SOLVER_H
#define FLOW_SOLVER_H

//#include "stdafx.h"
#include <cmath>
#include <sstream>  // string streams
#include <iostream> // for the use of 'cout'
#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"
//#include "Particle_sim.h"

#if defined __linux__ || defined __APPLE__
#include <sys/stat.h>//for mkdir
#else
#include "direct.h"//for mkdir
#endif

class Flow_solver
{
private:
    int az;

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

    int mkmk = 5;

public:
    Flow_solver();
    void initialize();
    void equilibrium(double, double, double); // compute the equilibrium populations from the fluid density and velocity
    void LBM(int); // perform LBM operations
    void BC(int); // perform boundary conditions operations
    void momenta(); // compute fluid density and velocity from the populations
    double meq(double, double, double, int);
    bool check_residual();
    ~Flow_solver();
};

#endif
