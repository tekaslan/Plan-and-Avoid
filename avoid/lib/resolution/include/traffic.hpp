#pragma once

#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <cassert>
#include <utility>
#include <nlopt.hpp>
#include <chrono>
#include <numeric>

#include "structs.hpp"

class Traffic {
public:
    void loadDynamicTraffic(const std::filesystem::path& filedir);
    void turnRadius(const size_t& id,
                    const size_t& t,
                    double& R,
                    const DynamicTraffic& airtraffic) const;
    AircraftType type(const size_t& id) const;
    std::string typeString(const size_t& id) const;

    DynamicTraffic data;
};