#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "geo.h"
#include "dubins.h"
#include "aclm.h"

// Separation parameters
inline constexpr double RNP1           = 1.0;   // NM
inline constexpr double RNP03          = 0.3;   // NM
inline constexpr double ADSB_ERR_P     = 0.05;  // NM
inline constexpr double IAF_ALT        = 2500;  // FT
inline constexpr double MSA            = 2500;  // FT
inline constexpr double BARO_ERR_A5k   = 210;   // FT
inline constexpr double BARO_ERR_B5k   = 160;   // FT
inline constexpr double BARO_ALT_LIM   = 5000;  // FT
inline constexpr double HORZ_CLEARANCE = 3;     // NM
inline constexpr double VERT_CLEARANCE = 500;   // FT

struct Hold {
    const char* name;
    Pos pos;
};

inline constexpr std::array<Hold, 8> Holds = {{
    {"ALWYZ", {.lat = 38.9835, .lon = -77.65467, .alt = 10000.0, .hdg = 119}},
    {"MOEJO", {.lat = 38.5760, .lon = -77.11450, .alt = 10000.0, .hdg = 28}},
    {"HALLL", {.lat = 38.5190, .lon = -77.00950, .alt = 10000.0, .hdg = 53}},
    {"SKILS", {.lat = 39.5033, .lon = -76.63117, .alt = 10000.0, .hdg = 195}},
    {"TRISH", {.lat = 39.5035, .lon = -76.37383, .alt = 8000.0,  .hdg = 225}},
    {"BILIT", {.lat = 38.7558, .lon = -76.06750, .alt = 11000.0, .hdg = 277}},
    {"IRONS", {.lat = 38.5312, .lon = -77.10617, .alt = 4500.0,  .hdg = 18}},
    {"OOURR", {.lat = 38.7328, .lon = -77.8655,  .alt = 13000.0, .hdg = 79}}
}};

struct FSSample {
    Pos pos{};
    double gs = 0.0;      // Groundspeed [kts]
    int64_t ts = 0;       // Timestamp [unix]
    bool valid = false;   // Validity flag
};

enum class AircraftType : uint8_t {
    GA         = 1,
    Turbojet   = 2,
    Helicopter = 3
};

enum class PathType : int {
    Search     = 0,
    Dubins     = 1
};

struct FlightSeries {
    std::vector<FSSample> samples;
    AircraftType type;
};

struct DynamicTraffic {
    std::string datafile;
    char label[10] = {};
    uint32_t nflights = 0;
    uint32_t nsteps   = 0;
    uint32_t median   = 0;
    uint32_t p90      = 0;
    uint32_t peak     = 0;
    std::vector<FlightSeries> flights;
};

struct ConflictRow {
    int t = 0;      // Timestamp / index

    double Olat = 0.0;
    double Olon = 0.0;
    double Oalt = 0.0;
    double Ohdg = 0.0;

    double Ilat = 0.0;
    double Ilon = 0.0;
    double Ialt = 0.0;
    double Ihdg = 0.0;

    int Iid = 0;
    double JCPA = 0.0;
};

struct RiskMetrics {
    double integral = 0.0;
    double worst    = 0.0;
};

struct ResolutionLogRow {
    size_t case_id{};
    size_t intruder_id{};
    std::string intruder_type{};
    size_t conflict_time{};
    std::string advisory_type{};
    size_t exit_flag{};
    size_t opt_exit_flag{};
    double onsetTime{};
    double conflictTime{};
    RiskMetrics plan_traffic_risk{};
    RiskMetrics adv_traffic_risk{};
    RiskMetrics adv_ownship_risk{};

    std::optional<double> V_star{};
    std::optional<double> t0_star{};
    std::optional<double> dT_star{};

    std::optional<double> dH_star{};
    std::optional<double> T_star{};
    std::optional<double> gamma0_star{};
    std::optional<double> gamma1_star{};

    std::optional<double> theta_star{};
    std::optional<double> tf_star{};
    std::optional<double> length_star{};

    std::optional<double> gamma_GA_star{};
    std::optional<size_t> hold_idx_star{};

    std::optional<double> g1_star{};
    std::optional<double> g2_star{};
    std::optional<double> f_star{};
    std::optional<bool> feasible{};
    std::optional<double> runtime_ms{};
};