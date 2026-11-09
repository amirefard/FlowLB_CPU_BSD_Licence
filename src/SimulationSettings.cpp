#include "SimulationSettings.h"
//=============================================================================================================================
//=============================================================================================================================
//-----------------------------------------------------------------------------------------------------------------------------
// Subroutine: Constructor
//-----------------------------------------------------------------------------------------------------------------------------
SimulationSettings::SimulationSettings(const std::string& filename)
{
   settings1 = read_settings(filename);
   //std::cout <<"Salam " << settings1["!Flow_With_Thermal_Effect"] <<std::endl;
}
//=============================================================================================================================
//=============================================================================================================================
SimulationSettings::~SimulationSettings()
{
    /*if (file.is_open())
    {
        file.close();
    }*/
}
//=============================================================================================================================
//=============================================================================================================================
SimulationSettings& SimulationSettings::get_instance(const std::string& filename)
{
    static SimulationSettings instance(filename); // Use the provided filename
    return instance;
}
//=============================================================================================================================
//=============================================================================================================================
/*std::string SimulationSettings::get_setting(const std::string& name) const
{
    auto it = settings1.find(name);
    if (it != settings1.end())
    {
        return it->second;
    }
    else
    {
        std::cerr << "Error: Setting " << name << " not found" << std::endl;
        return ""; // Or throw an exception
    }
}*/
//=============================================================================================================================
//=============================================================================================================================
/*INT32 SimulationSettings::get_int_setting(const std::string& name) const
{
    std::string value = get_setting(name);
    try
    {
        return std::stoi(value);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Failed to convert setting " << name << " to INT32" << std::endl;
        return 0; // Or throw an exception
    }
}*/
//=============================================================================================================================
//=============================================================================================================================
bool SimulationSettings::get_bool_setting(const std::string& name) const
{
    auto it = settings1.find(name);
    if (it != settings1.end())
    {
        return it->second;
    }
    else
    {
        std::cerr << "Error: Setting " << name << " not found in settingsLB file" << std::endl;
        std::exit(EXIT_FAILURE);
        return false; // Or throw an exception
    }

    /*std::string value = get_setting(name);

    if (value == "1")
    {
        return true;
    }
    else if (value == "0")
    {
        return false;
    }
    else
    {
        throw std::invalid_argument("Invalid string for boolean conversion");
    }*/
}
//=============================================================================================================================
//=============================================================================================================================
INT32 SimulationSettings::get_int_setting(const std::string& name) const
{
    auto it = settings1.find(name);
    if (it != settings1.end())
    {
        return it->second;
    }
    else
    {
        std::cerr << "Error: Setting " << name << " not found in settingsLB file" << std::endl;
        std::exit(EXIT_FAILURE);
        return false; // Or throw an exception
    }
}
//=============================================================================================================================
//=============================================================================================================================
// bool SimulationSettings::string_to_bool(const std::string& s) {
//     if (s == "1") return true;
//     if (s == "0") return false;
//     throw std::invalid_argument("Invalid string for bool conversion: " + s);
// }
//=============================================================================================================================
//=============================================================================================================================
INT32 SimulationSettings::string_to_int(const std::string& s) {
    try {
        return std::stoi(s);
    } catch (const std::exception& e) {
        throw std::invalid_argument("Invalid string for INT32 conversion: " + s);
    }
}
//=============================================================================================================================
//=============================================================================================================================
std::map<std::string, INT32> SimulationSettings::read_settings(const std::string& filename)
{
    std::map<std::string, INT32> settings;

    std::ifstream file(filename);
    if (!file.is_open()) {
        //throw std::runtime_error("********************* Unable to find/open file: "+ filename + " *********************\n");
        std::cerr << "Error: Unable to open/find file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line))
    {
        // Remove CR from Windows line endings
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string key, value;

        // Read key
        if (std::getline(ss, key, '='))
        {
            // Remove whitespace from key
            key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());

            // Read value
            if (std::getline(ss, value))
            {
                // Remove whitespace from value
                value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());

                // Insert key-value pair into settings map
                settings[key] = string_to_int(value);
            }
            else
            {
                std::cerr << "Error: Missing value for setting " << key << std::endl;
            //    return false;
            }
        }
        else
        {
            std::cerr << "Error: Invalid line: " << line << std::endl;
        //    return false;
        }
    }
    if (file.is_open())
    {
        file.close();
    }
    return settings;
}


