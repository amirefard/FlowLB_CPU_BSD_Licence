// Here are general functions that all classes have access to.
#ifndef GLOBALFUNCTIONS_H
#define GLOBALFUNCTIONS_H

#include <cstddef> // Include the header for size_t
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <cstdint>// for int32_t, etc
constexpr double M_PI = 3.14159265358979323846;

using INT32 = std::int32_t;//signed 32-bit
using INT16 = std::int16_t;//signed 16-bit

using UINT32 = std::uint32_t;  // unsigned 32-bit
using UINT16 = std::uint16_t;  // unsigned 16-bit

///***********************************************************************
/// You can set the following two macros
///***********************************************************************
#define STORE_1D_VECTOR 1   //Set 0 or 1. To store pops and velocity etc as 1D vector (much faster)

#define UsePreCalculatedKernels 1    // To use pre-calculated kernel functions (2-point, 4-point, etc) instead of calculating on the fly.

///***********************************************************************
/// Do not change the following Macros
///***********************************************************************
#define MASTER 0   //Do not change. processor_id of first task

#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define TRUE  1
#define FALSE 0
///***********************************************************************
/// Global functions
///***********************************************************************
void swap_endian(char* buffer, INT32 size);
void swap_endian(char* buffer3, char* buffer2, INT32 size);
void write_bigendian(std::ofstream& file, char* buffer, INT32 count, INT32 size);

template <typename T>
void write_bigendian(std::ofstream& file, T* buffer, INT32 count);
std::string numberToString(INT32 number);
std::vector<std::string> readConfigFile(const std::string& filename);

const double Dt_LB = 1.0;//Do not change
const double Dh_LB = 1.0;//Do not change

#endif // GLOBALFUNCTIONS_H
