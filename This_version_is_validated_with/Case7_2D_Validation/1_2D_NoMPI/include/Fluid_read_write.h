#ifndef FLUID_READ_WRITE_H
#define FLUID_READ_WRITE_H

#include <stdint.h>
//#include <conio.h>
#include <iostream> // for the use of 'cout'
#include <cmath>
#include <fstream>  // file stream
#include <sstream>  // string streams
#include <cstdlib>  // standard library
#include <iomanip>  // For set precision.
#include <iostream>
#include <string>
#include <stdint.h>//to use int32_t

#include "LB_Macros.h"
#include "LB_Settings.h"
#include "Global_variables.h"

class Fluid_read_write
{
public:

    Fluid_read_write();
    void swap_endian(char *, size_t);
    void swap_endian(char *, char *, size_t);
    void write_bigendian(std::ofstream & , char * , size_t , size_t );
    template <typename T>
    void write_bigendian(std::ofstream & , T  *, size_t );
    void write_fluid_vtk(int);
    void write_solid_vtk(int, int);
    void write_fluid_data(int);
    void write_fluid_recovery(int);
    void read_fluid_recovery(int &);
    void write_parametes();
    ~Fluid_read_write();
};

#endif
