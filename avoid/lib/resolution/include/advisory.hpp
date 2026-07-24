#pragma once

#include <filesystem>
#include <optional>
#include <utility>
#include <vector>

#include <nlopt.hpp>
#include "conflict.hpp"

class Halt : public Conflict {
public:
    bool isTakeoff(const ConflictRow& conflict,
                   const DynamicTraffic& airtraffic) const;

    std::vector<FSSample> run(const ConflictRow& instance,
                              const DynamicTraffic& airtraffic,
                              const std::filesystem::path& caseDir,
                              ResolutionLogRow& resolutionLogRow);

    double runtime = 0.0;
};

class HeliSpeed : public Conflict {
public:
    struct Params {
        double t0 = 0.0;
        double deltaV = 0.0;
        double dVMax = 0.0;
        double maxGS = 0.0;
        double thold = 0.0;
        double tholdMax = 0.0;
        double acc = 0.0;
        double decel = 0.0;
        RiskMetrics trafficInteractionRisk{};
        RiskMetrics ownshipInteractionRisk{};
    };

    struct OptimizerInput {
        const HeliSpeed* self = nullptr;
        const Conflict* conflict = nullptr;
        const ConflictRow* instance = nullptr;
        DynamicTraffic* airtraffic = nullptr;
        Params base{};
        Conflict::OptimSettings settings{};
    };

    Conflict::AdvisoryState<Params> state{};
    
    std::vector<FSSample> speedCommand(const ConflictRow& instance,
                                       Params& params,
                                       DynamicTraffic& airtraffic,
                                       std::optional<std::filesystem::path> caseDir = std::nullopt) const;

    std::vector<FSSample> timeShiftTrajectory(const ConflictRow& instance,
                                              Params& params,
                                              DynamicTraffic& airtraffic,
                                              std::optional<std::filesystem::path> caseDir = std::nullopt) const;

    static double nonlcon(const std::vector<double>& x,
                          std::vector<double>& grad,
                          void* data);

    static double nonlconDomino(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data);

    static double cost(const std::vector<double>& x,
                       std::vector<double>& grad,
                       void* data);

    void optimize(const Params& params,
                  const ConflictRow& instance,
                  const Conflict& conflict,
                  DynamicTraffic& airtraffic,
                  const Conflict::OptimSettings& settings);

    std::vector<FSSample> run(const ConflictRow& instance,
                              const Conflict& conflict,
                              DynamicTraffic& airtraffic,
                              const std::filesystem::path& dir,
                              ResolutionLogRow& resolutionLogRow);
};

class Speed : public Conflict {
public:
    struct Params {
        double deltaV = 0.0;
        double dVMax = 0.0;
        std::pair<double, double> vaBounds = {0.0, 0.0};
        double avgGS = 0.0;
        double t0 = 0.0;
        double thold = 0.0;
        double tholdMax = 0.0;
        double acc = 0.0;
        double decel = 0.0;
        size_t validt0 = 0;
        size_t validtf = 0;
        RiskMetrics trafficInteractionRisk{};
        RiskMetrics ownshipInteractionRisk{};
    };

    struct OptimizerInput {
        const Speed* self = nullptr;
        const Conflict* conflict = nullptr;
        const ConflictRow* instance = nullptr;
        const DynamicTraffic* airtraffic = nullptr;
        Params base{};
        Conflict::OptimSettings settings{};
    };

    Conflict::AdvisoryState<Params> state{};

    std::pair<double, double> groundToAirspeed(const double groundspeed,
                                                const double course,
                                                const double windspeed,
                                                const double windDirection) const;

    double airToGroundspeed(double airspeed,
                            double course,
                            double windspeed,
                            double windDirection) const;

    void cumulativeDistance(double timestep,
                            const std::vector<FSSample>& samples,
                            std::vector<double>& cumulDistance) const;

    std::vector<FSSample> speedCommand(const ConflictRow& instance,
                                       Params& params,
                                       const DynamicTraffic& airtraffic,
                                       std::optional<std::filesystem::path> caseDir = std::nullopt) const;

    std::vector<FSSample> timeShiftTrajectory(const ConflictRow& instance,
                                              Params& params,
                                              const DynamicTraffic& airtraffic,
                                              std::optional<std::filesystem::path> caseDir = std::nullopt) const;

    static double nonlcon(const std::vector<double>& x,
                          std::vector<double>& grad,
                          void* data);

    static double nonlconDomino(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data);

    static double cost(const std::vector<double>& x,
                       std::vector<double>& grad,
                       void* data);

    void optimize(const Params& params,
                  const ConflictRow& instance,
                  const Conflict& conflict,
                  const DynamicTraffic& airtraffic,
                  const Conflict::OptimSettings& settings);

    std::vector<FSSample> run(const ConflictRow& instance,
                              const Conflict& conflict,
                              const Traffic& traffic,
                              const std::filesystem::path& dir,
                              ResolutionLogRow& resolutionLogRow);
};

class Altitude : public Conflict {
public:
    struct Params {
        double deltaH = 0.0;
        double deltaHMax = 0.0;
        double t0 = 0.0;
        double thold = 0.0;
        double tholdMax = 0.0;
        double gamma0 = 0.0;
        double gamma1 = 0.0;
        size_t idx0 = 0;
        size_t idxf = 0;
        std::pair<double, double> gamma0bounds{};
        std::pair<double, double> gamma1bounds{};
        std::pair<double, double> t0bounds{};
        RiskMetrics trafficInteractionRisk{};
        RiskMetrics ownshipInteractionRisk{};
    };

    struct OptimizerInput {
        const Altitude* self = nullptr;
        const ConflictRow* instance = nullptr;
        const Conflict* conflict = nullptr;
        const DynamicTraffic* airtraffic = nullptr;
        Params base{};
        Conflict::OptimSettings settings{};
    };

    Conflict::AdvisoryState<Params> state{};

    std::vector<FSSample> run(const ConflictRow& instance,
                              const Conflict& conflict,
                              const DynamicTraffic& airtraffic,
                              const std::filesystem::path& dir,
                              ResolutionLogRow& resolutionLogRow);

    double verticalRate(const std::vector<FSSample>& samples,
                        size_t start,
                        size_t deltaT) const;

    int isClimb(const std::vector<FSSample>& samples,
                size_t start,
                size_t deltaT) const;

    std::vector<FSSample> shiftVerticalProfile(const ConflictRow& instance,
                                               const Params& params,
                                               const DynamicTraffic& airtraffic,
                                               std::optional<std::filesystem::path> caseDir = std::nullopt) const;

    static double nonlcon(const std::vector<double>& x,
                          std::vector<double>& grad,
                          void* data);

    static double nonlconDomino(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data);

    static double cost(const std::vector<double>& x,
                       std::vector<double>& grad,
                       void* data);

    void optimize(const Params& params,
                  const ConflictRow& instance,
                  const Conflict& conflict,
                  const DynamicTraffic& airtraffic,
                  const Conflict::OptimSettings& settings);
};

class Extend : public Conflict {
public:
    struct Params {
        double aux = 0;
        double theta = 60.0;
        double t0 = 0.0;
        double tf = 0.0;
        double tmax = 0.0;
        double length = 0.0;
        RiskMetrics trafficInteractionRisk{};
        RiskMetrics ownshipInteractionRisk{};
    };

    struct OptimizerInput {
        Extend* self = nullptr;
        const ConflictRow* instance = nullptr;
        const Conflict* conflict = nullptr;
        const DynamicTraffic* airtraffic = nullptr;
        Params base{};
        Conflict::OptimSettings settings{};
    };

    DubinsPath path[2]{};
    Conflict::AdvisoryState<Params> state{};

    double mapTheta(double z);

    Pos interpPos(const std::pair<Pos, double>& x1,
                  const std::pair<Pos, double>& x2,
                  double t) const;

    void buildSTurnPath(DubinsPath& dubins,
                        const Params& params,
                        const Conflict& conflict);

    void generateSTurnPath(const ConflictRow& instance,
                        const Params& params,
                        const Conflict& conflict,
                        const DynamicTraffic& airtraffic);

    std::vector<FSSample> DubinsCoordinates(const Conflict& conflict,
                                            const double groundspeed);

    std::vector<FSSample> divert(const ConflictRow& instance,
                                 const Params& params,
                                 const Conflict& conflict,
                                 const DynamicTraffic& traffic,
                                 std::optional<std::filesystem::path> caseDir = std::nullopt);

    static double nonlcon(const std::vector<double>& x,
                          std::vector<double>& grad,
                          void* data);

    static double nonlconDomino(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data);

    static double cost(const std::vector<double>& x,
                       std::vector<double>& grad,
                       void* data);

    void optimize(const Params& params,
                  const ConflictRow& instance,
                  const Conflict& conflict,
                  const DynamicTraffic& airtraffic,
                  const Conflict::OptimSettings& settings);

    std::vector<FSSample> run(const ConflictRow& instance,
                              const Conflict& conflict,
                              const DynamicTraffic& airtraffic,
                              const std::filesystem::path& dir,
                              ResolutionLogRow& resolutionLogRow);
};

class Divert : public Conflict {
public:
    struct Params {
        double gamma = 0.0;
        double t0 = 0.0;
        std::pair<double, double> gammaBounds{};
        std::pair<double, double> t0Bounds{};
        size_t hold_idx = 0;
        RiskMetrics trafficInteractionRisk{};
        RiskMetrics ownshipInteractionRisk{};
    };

    struct OptimizerInput {
        Divert* self = nullptr;
        const ConflictRow* instance = nullptr;
        const Conflict* conflict = nullptr;
        const DynamicTraffic* airtraffic = nullptr;
        Params base{};
        Conflict::OptimSettings settings{};
    };

    DubinsPath path{};
    Conflict::AdvisoryState<Params> state{};

    void computePath(const ConflictRow& instance,
                     const Params& params,
                     const Conflict& conflict,
                     const DynamicTraffic& airtraffic);

    std::vector<FSSample> pathToHold(const ConflictRow& instance,
                                     const Params& params,
                                     const Conflict& conflict,
                                     const DynamicTraffic& airtraffic,
                                     std::optional<std::filesystem::path> caseDir = std::nullopt);

    static double nonlcon(const std::vector<double>& x,
                          std::vector<double>& grad,
                          void* data);

    static double nonlconDomino(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data);

    static double cost(const std::vector<double>& x,
                       std::vector<double>& grad,
                       void* data);

    void optimize(const Params& params,
                  const ConflictRow& instance,
                  const Conflict& conflict,
                  const DynamicTraffic& airtraffic,
                  const Conflict::OptimSettings& settings);

    std::vector<FSSample> run(const ConflictRow& instance,
                              const Conflict& conflict,
                              const DynamicTraffic& airtraffic,
                              const std::filesystem::path& dir,
                              ResolutionLogRow& resolutionLogRow);
};

bool invalidPos(const Pos& p);