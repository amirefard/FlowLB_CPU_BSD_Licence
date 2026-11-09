#ifndef INLET_OUTLETS_H
#define INLET_OUTLETS_H
#include "LB_Settings.h"
#include "Global_variables.h"
#include "LB_Macros.h"

#include <iostream> // for the use of 'cout'
#include <fstream>  // file stream
#include <sstream>  // string streams
#include <string>

class Inlet_outlets
{
public:
	Inlet_outlets();
	void set_boundaries();
	void set_boundaries_thermal();
	~Inlet_outlets();
};

#endif
