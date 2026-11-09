#ifndef SIMULATIONSETTINGS_H
#define SIMULATIONSETTINGS_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <algorithm> // Add this header for the remove_if function

#include "GlobalFunctions.h"

class SimulationSettings
{
public:
    // Destructor
    ~SimulationSettings();

    // Static member function to access the singleton instance
    static SimulationSettings& get_instance(const std::string& filename);
    bool get_bool_setting(const std::string& name) const;
    INT32  get_int_setting (const std::string& name) const;

    std::map<std::string, INT32> read_settings(const std::string& filename);

private:
    // Private constructor to prevent instantiation
    SimulationSettings(const std::string& filename);
    //bool string_to_bool(const std::string& s);
    INT32 string_to_int(const std::string& s);

    std::map<std::string, INT32> settings1;
};

#endif // SIMULATIONSETTINGS_H
