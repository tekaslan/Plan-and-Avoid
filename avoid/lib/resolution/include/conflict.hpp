#pragma once

#include <filesystem>
#include <fstream>
#include <tuple>
#include <vector>

#include "geo.h"
#include "aclm.h"
#include "dubins.h"
#include "structs.hpp"
#include "traffic.hpp"

class Conflict {
public:

    // Wind conditions
    double windSpeed = 0.0;
    double windDirection = 0.0;

    // Emergency landing solution runtime
    double pathRuntime = 0.0; // [s]
    void getPathRuntime(const std::filesystem::path& filedir);

    // Datalink time delay
    double tdatalink   = 0.5; // [s]

    // Advisory execution onset
    size_t onsetTime   = 0; // [s]

    // Ego airspeed
    double egoSpeedNom  = 90.0; // [kts]
    double egoSpeedFast = 95.0; // [kts]
    double egoSpeedSlow = 85.0; // [kts]

    std::vector<FSSample> emergencyPath;
    std::vector<FSSample> emergencyPath_fast;
    std::vector<FSSample> emergencyPath_slow;

    enum class Optimizer {
        local  = 0,
        global = 1
    };
    
    struct OptimSettings {
        Optimizer method = Optimizer::local;
        unsigned long seed;
        double xtol = 1e-4;
        double ftol = 1e-4;
        double nonlconTol = 1e-6;
        size_t maxFunEval = 200;
        double tlimit = 3.0;
    };

    template <typename ParamsT>
    struct AdvisoryState {
        int exitFlag = 0;
        ParamsT optimal{};
        double fval = 0.0;
        double nonlcon = 0.0;
        double nonlconDomino = 0.0;
        bool feasible = false;
        double runtime = 0.0;
    };

    std::vector<size_t> resolve(size_t case_id,
                                const std::filesystem::path& dir,
                                std::ofstream& caselog,
                                std::ofstream& resolog) const;

public:
    GeoOpt geoopt{.model = WGS84};
    DubinsOpt dubinsopt{.verbose = 0};

    double tlimit = 3.0;
    Optimizer opt = Optimizer::local;
    int population = 4;

    void loadEmergencyPath(const std::filesystem::path& filedir);
    FSSample pathInstance(size_t idx) const { return emergencyPath.at(idx); }
    std::vector<FSSample> path() const { return emergencyPath; }
    size_t pathSize() const { return emergencyPath.size(); }

    int loadConflictLog(const std::filesystem::path& filedir);
    ConflictRow conflictInstance(size_t idx) const;
    size_t conflictCount() const { return log.size(); }
    std::vector<std::tuple<size_t, size_t, ConflictRow>> unique() const;
    std::vector<FSSample> retimeEmergencyPathBySpeed(double egoSpeedKts,
                                                  double timestep = 1.0);

    void computeSeparation(const Pos& ownship,
                           const Pos& intruder,
                           double& dh,
                           double& dv) const;

    double CPACost(const double& deltah, const double& deltav) const;

    RiskMetrics trafficRisk(const ConflictRow& instance,
                            const std::vector<FSSample>& samples,
                            const DynamicTraffic& airtraffic) const;
    RiskMetrics ownshipRisk(const std::vector<FSSample>& ownship,
                            const DynamicTraffic& airtraffic) const;
    double ownshipMargin(const std::vector<FSSample>& ownship,
                        const DynamicTraffic& airtraffic) const;
    double ownshipMarginSingle(const std::vector<FSSample>& ownship,
                                const std::vector<FSSample>& intruder) const;

    double ownshipMarginRobustTimeWindow(const std::vector<FSSample>& ownship,
                                        const std::vector<FSSample>& intruder,
                                        double VminKts,
                                        double VnomKts,
                                        double VmaxKts,
                                        double timestep) const;           

    airspaceHeatmap prohibited{
        .LAT_SIZE = 250,
        .LON_SIZE = 250,
        .ALT_SIZE = 2,
        .grid_on_cells = 0,
        .gridfile = "lib/airtraffic/data/nofly_grid.bin",
        .datafile = "lib/airtraffic/data/nofly.bin"
    };

private:
    std::vector<ConflictRow> log;
};

// Free helpers
std::vector<size_t> loadLogCSV(const std::filesystem::path& p);
void copyPositions(Pos& ownship, Pos& intruder, const ConflictRow& c);
void printCasesHeader(std::ofstream& log);
void printCaseRow(std::ofstream& log,
                  size_t case_id,
                  const std::string& case_status,
                  size_t n_conflicts,
                  size_t n_resolved,
                  size_t n_failed,
                  double runtime_ms,
                  RiskMetrics risk);
void printResolutionsHeader(std::ofstream& log);
void printResolutionRow(std::ofstream& log, const ResolutionLogRow& row);
void printOptional(std::ofstream& log, const std::optional<double>& x);
void printOptionalBool(std::ofstream& log, const std::optional<bool>& x);