#ifndef RUNFLOWWITHOUTPARTICLE_H
#define RUNFLOWWITHOUTPARTICLE_H

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
// CLASSES
class SimulationParameters;
class SimulationSettings;
class ParallelMPI;

//=============================================================================================================================
void run_flow_without_particle(const SimulationParameters&, const SimulationSettings&, const ParallelMPI&);

#endif // RUNFLOWWITHOUTPARTICLE_H
