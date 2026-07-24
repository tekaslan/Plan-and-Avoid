#include "advisory.hpp"

static FSSample NaNFSSample() {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    FSSample p{};
    p.pos.lat = NaN; p.pos.lon = NaN; p.pos.alt = NaN; p.pos.hdg = NaN;
    return p;
}

bool invalidPos(const Pos& p) {
    return !std::isfinite(p.lat) || !std::isfinite(p.lon) || !std::isfinite(p.alt);
}

/*
    Identifies whether flight data has take-off phase
*/
bool Halt::isTakeoff(const ConflictRow& instance,
                    const DynamicTraffic& airtraffic) const{

    // Initial active time index
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    auto it = std::find_if(samples.begin(), samples.end(), [](const FSSample& s){return s.valid;});
    if (it == samples.end()) return false;
    size_t init = std::distance(samples.begin(), it);

    // Check altitude change
    size_t dt = 30; // 30 seconds window
    if (init + dt >= samples.size()) return false;
    double dh = samples[init + dt].pos.alt - samples[init].pos.alt; // [ft]
    double vrate = dh / dt * 60; // [ft/min]

    return ((init > 0) && (vrate > 200));
}


std::vector<FSSample> Halt::run(const ConflictRow& instance,
                                const DynamicTraffic& airtraffic,
                                const std::filesystem::path& casedir,
                                ResolutionLogRow& ResolutionLogRow){
    
    auto start = std::chrono::high_resolution_clock::now();

    if (!isTakeoff(instance, airtraffic)) return {};

    // Modify air traffic
    const auto& plan = airtraffic.flights[instance.Iid].samples;
    std::vector<FSSample> newTraj(airtraffic.nsteps, NaNFSSample());
    for (size_t i = 0; i < airtraffic.nsteps; ++i) {
        newTraj[i].ts = i;
    }

    // Runtime
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    runtime = elapsed.count() * 1000;

    // Log resolution
    ResolutionLogRow.advisory_type  = "Halt T/O";
    ResolutionLogRow.feasible       = 1;
    ResolutionLogRow.runtime_ms     = runtime;
    ResolutionLogRow.exit_flag      = 0;

    // Log NAN coordinates to halt takeoff
    std::filesystem::path dir = casedir / (std::to_string(instance.Iid) + "_trajectory.csv");
    std::ofstream log(dir);
    for (size_t i = 0; i < airtraffic.nsteps; ++i) {
        log << plan[i].pos.lat << ","
            <<  plan[i].pos.lon << ","
            <<  plan[i].pos.alt << ","
            << plan[i].pos.hdg << ","
            << newTraj[i].pos.lat << ","
            << newTraj[i].pos.lon << ","
            << newTraj[i].pos.alt << ","
            << newTraj[i].pos.hdg << "\n";
    }
    log.close();

    return newTraj;
}

/*
    HELICOPTER SPEED REGULATION
*/

std::vector<FSSample> HeliSpeed::speedCommand(
    const ConflictRow& instance,
    Params& params,
    DynamicTraffic& airtraffic,
    std::optional<std::filesystem::path> casedir) const
{
    // Get the planned intruder path
    if ((size_t)instance.Iid >= airtraffic.flights.size()) {
        throw std::invalid_argument("Invalid intruder ID.");
    }

    const auto& samples = airtraffic.flights[(size_t)instance.Iid].samples;

    if (samples.empty()) {
        throw std::invalid_argument("Intruder trajectory is empty.");
    }

    std::vector<FSSample> TrajShifted(samples.begin(), samples.end());

    const size_t N = samples.size();

    // Pre-allocate speed shift array
    std::vector<double> GSShift(N, 0.0);

    // Compute advisory start index
    size_t t0 = static_cast<size_t>(std::max(0.0, params.t0));
    t0 = std::min(t0, N - 1);

    size_t shiftStart = std::max(t0, onsetTime);
    shiftStart = std::min(shiftStart, N - 1);

    const size_t conflictTime = std::min(static_cast<size_t>(instance.t), N - 1);

    if (shiftStart > conflictTime) {
        throw std::runtime_error("shiftStart > time of conflict.");
    }

    // Optimizer-requested ground-speed offset [kts].
    // This may not be physically achievable because gs is saturated.
    const double deltaVCmd = params.deltaV;

    // Ramp offset from 0 to deltaVCmd
    double dv_cmd = 0.0;
    size_t holdStart = N;

    for (size_t i = shiftStart; i < N; ++i) {

        if (!std::isfinite(samples[i].gs)) {
            continue;
        }

        const double diff = deltaVCmd - dv_cmd;

        if (diff > 0.0) {
            dv_cmd += std::min(params.acc, diff);
        } else {
            dv_cmd += std::max(params.decel, diff);
        }

        // Saturated speed command.
        // VTOL may hover, but ground speed cannot be negative.
        const double gsNom = samples[i].gs;
        const double gsAdv = std::clamp(gsNom + dv_cmd, 0.0, params.maxGS);

        TrajShifted[i].gs = gsAdv;
        GSShift[i] = gsAdv - gsNom;

        if (std::abs(dv_cmd - deltaVCmd) < 1e-9) {
            holdStart = i + 1;
            break;
        }
    }

    if (holdStart == N) {
        return {};
    }

    // Hold offset
    const size_t holdEnd = std::min(N, holdStart + static_cast<size_t>(params.thold));

    for (size_t i = holdStart; i < holdEnd; ++i) {

        if (!std::isfinite(samples[i].gs)) {
            continue;
        }

        const double gsNom = samples[i].gs;
        const double gsAdv = std::clamp(gsNom + deltaVCmd, 0.0, params.maxGS);

        TrajShifted[i].gs = gsAdv;
        GSShift[i] = gsAdv - gsNom;
    }

    // Recover offset from deltaVCmd to 0
    dv_cmd = deltaVCmd;
    size_t recoveryEnd = N;

    for (size_t i = holdEnd; i < N; ++i) {

        if (!std::isfinite(samples[i].gs)) {
            continue;
        }

        const double diff = 0.0 - dv_cmd;

        if (diff > 0.0) {
            dv_cmd += std::min(params.acc, diff);
        } else {
            dv_cmd += std::max(params.decel, diff);
        }

        const double gsNom = samples[i].gs;
        const double gsAdv = std::clamp(gsNom + dv_cmd, 0.0, params.maxGS);

        TrajShifted[i].gs = gsAdv;
        GSShift[i] = gsAdv - gsNom;

        if (std::abs(dv_cmd) < 1e-9) {
            recoveryEnd = i + 1;
            break;
        }
    }

    if (recoveryEnd == N) {
        return {};
    }

    // Compute physically achieved deltaV after saturation
    double actualDeltaV = 0.0;

    // Prefer the mean achieved offset during the hold segment because this is
    // the steady part of the advisory.
    double sumShift = 0.0;
    size_t nShift = 0;

    for (size_t i = holdStart; i < holdEnd; ++i) {
        if (std::isfinite(GSShift[i])) {
            sumShift += GSShift[i];
            ++nShift;
        }
    }

    if (nShift > 0) {
        actualDeltaV = sumShift / static_cast<double>(nShift);
    } else {
        // Fallback: use largest achieved signed shift over the modified region.
        double bestAbs = 0.0;

        for (size_t i = shiftStart; i < recoveryEnd; ++i) {
            if (!std::isfinite(GSShift[i])) {
                continue;
            }

            if (std::abs(GSShift[i]) > bestAbs) {
                bestAbs = std::abs(GSShift[i]);
                actualDeltaV = GSShift[i];
            }
        }
    }

    // Replace requested deltaV with the physically achieved one.
    // This makes logged V_star and parameter statistics reflect the actual
    // speed change instead of the optimizer request.
    params.deltaV = actualDeltaV;

    // Write to file
    if (casedir) {
        std::filesystem::path dir =
            casedir.value() / ("groundspeed" + std::to_string(instance.Iid) + ".csv");

        std::ofstream file(dir);

        file << "gs_nom,gs_adv,dv_actual,dv_requested,dv_reported\n";

        for (size_t i = 0; i < N; ++i) {
            file << samples[i].gs << ","
                 << TrajShifted[i].gs << ","
                 << GSShift[i] << ","
                 << deltaVCmd << ","
                 << params.deltaV << std::endl;
        }

        file.close();
    }

    return TrajShifted;
}


std::vector<FSSample> HeliSpeed::timeShiftTrajectory(
    const ConflictRow& instance,
    Params& params,
    DynamicTraffic& airtraffic,
    std::optional<std::filesystem::path> casedir) const
{
    const double timestep = 1.0;

    // Get the planned trajectory
    if ((size_t)instance.Iid >= airtraffic.flights.size()) {
        throw std::invalid_argument("Invalid intruder ID.");
    }

    const auto& samples = airtraffic.flights[(size_t)instance.Iid].samples;

    if (samples.empty()) {
        return {};
    }

    const size_t N = samples.size();

    // Planned cumulative distance
    Speed speedReg;
    std::vector<double> plannedCumulDistance(N, 0.0);
    speedReg.cumulativeDistance(timestep, samples, plannedCumulDistance);

    const double plannedEnd = plannedCumulDistance.back();

    // Check monotonicity
    for (size_t k = 1; k < N; ++k) {
        assert(plannedCumulDistance[k] >= plannedCumulDistance[k - 1]);
    }

    // Change the speed profile
    auto TrajShifted = speedCommand(instance, params, airtraffic, casedir);

    if (TrajShifted.empty()) {
        return {};
    }

    // Shifted cumulative distance
    std::vector<double> shiftedCumulDistance(N, 0.0);
    speedReg.cumulativeDistance(timestep, TrajShifted, shiftedCumulDistance);

    // Check monotonicity
    for (size_t k = 1; k < N; ++k) {
        assert(shiftedCumulDistance[k] >= shiftedCumulDistance[k - 1]);
    }

    // New trajectory initialized as planned trajectory
    auto newTraj = samples;

    // Shift start time index
    size_t t0 = static_cast<size_t>(std::max(0.0, params.t0));
    t0 = std::min(t0, N - 1);

    size_t shiftStart = std::max(t0, onsetTime);
    shiftStart = std::min(shiftStart, N - 1);

    // Interpolate trajectory samples after shift
    for (size_t i = shiftStart; i < N; ++i) {

        const double d = shiftedCumulDistance[i];

        if (d >= plannedEnd) {

            // Reached/passed end of planned trajectory
            auto nansample = NaNFSSample();

            for (size_t k = i; k < N; ++k) {
                newTraj[k] = nansample;
                newTraj[k].ts = k;
            }

            break;
        }

        auto it = std::lower_bound(
            plannedCumulDistance.begin() + 1,
            plannedCumulDistance.end(),
            d
        );

        if (it == plannedCumulDistance.end()) {
            break;
        }

        size_t j = std::max<size_t>(
            1,
            static_cast<size_t>(std::distance(plannedCumulDistance.begin(), it))
        );

        double r = 0.0;

        if (plannedCumulDistance[j] > plannedCumulDistance[j - 1]) {
            r = (d - plannedCumulDistance[j - 1]) /
                (plannedCumulDistance[j] - plannedCumulDistance[j - 1]);
        }

        const Pos before = samples[j - 1].pos;
        const Pos after  = samples[j].pos;

        Pos shifted;
        shifted.lat = before.lat + r * (after.lat - before.lat);
        shifted.lon = before.lon + r * (after.lon - before.lon);
        shifted.alt = before.alt + r * (after.alt - before.alt);
        shifted.hdg = after.hdg;

        newTraj[i].pos = shifted;
        newTraj[i].gs  = TrajShifted[i].gs;
        newTraj[i].ts  = i;
    }

    // Log coordinates
    if (casedir) {
        std::filesystem::path dir =
            casedir.value() / (std::to_string(instance.Iid) + "_trajectory.csv");

        std::ofstream plan(dir);

        for (size_t i = 0; i < N; ++i) {
            plan << samples[i].pos.lat << ","
                 << samples[i].pos.lon << ","
                 << samples[i].pos.alt << ","
                 << samples[i].pos.hdg << ","
                 << samples[i].gs << ","
                 << newTraj[i].pos.lat << ","
                 << newTraj[i].pos.lon << ","
                 << newTraj[i].pos.alt << ","
                 << newTraj[i].pos.hdg << ","
                 << newTraj[i].gs << "\n";
        }

        plan.close();
    }

    return newTraj;
}


double HeliSpeed::nonlcon(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data)
{
    // Intentionally unused argument
    (void)grad;

    auto* d = static_cast<OptimizerInput*>(data);

    auto& params     = d->base;
    auto& airtraffic = *d->airtraffic;
    auto& conflict   = *d->conflict;
    auto& instance   = *d->instance;

    // Decision vector:
    // x[0] = ground-speed offset [kts]
    // x[1] = speed-change initial time [s]
    // x[2] = hold time [s]
    params.deltaV = x[0];
    params.t0     = x[1];
    params.thold  = x[2];

    const auto newTraj = d->self->timeShiftTrajectory(instance, params, airtraffic);

    if (newTraj.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    return conflict.ownshipMarginSingle(conflict.path(), newTraj);
}


double HeliSpeed::nonlconDomino(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data)
{
    // Intentionally unused argument
    (void)grad;

    auto* d = static_cast<OptimizerInput*>(data);

    auto& params     = d->base;
    auto& conflict   = d->conflict;
    auto& instance   = *d->instance;
    auto& airtraffic = *d->airtraffic;

    // Decision vector:
    // x[0] = ground-speed offset [kts]
    // x[1] = speed-change initial time [s]
    // x[2] = hold time [s]
    params.deltaV = x[0];
    params.t0     = x[1];
    params.thold  = x[2];

    const auto shifted = d->self->timeShiftTrajectory(instance, params, airtraffic);

    if (shifted.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    auto metric = conflict->trafficRisk(instance, shifted, airtraffic);

    return std::max(metric.integral - params.trafficInteractionRisk.integral,
                    metric.worst    - params.trafficInteractionRisk.worst);
}


double HeliSpeed::cost(const std::vector<double>& x,
                                 std::vector<double>& grad,
                                 void* data) {
    // Intentionally unused argument
    (void)grad;

    auto* d = static_cast<OptimizerInput*>(data);

    auto& settings = d->settings;
    auto& params   = d->base;
    auto& instance = *d->instance;

    const double speedScale = std::max(params.maxGS, 1.0);
    const double timeScale  = std::max(static_cast<double>(instance.t), 1.0);
    const double holdScale  = std::max(params.tholdMax, 1.0);

    // x[0] = deltaV [kts]
    // x[1] = t0 [s]
    // x[2] = thold [s]
    const double cV = x[0] / speedScale;
    const double c0 = x[1] / timeScale;
    const double cT = x[2] / holdScale;

    double cost = 0.0;
    cost += cV * cV;
    cost += c0 * c0;
    cost += cT * cT;

    // Penalty formulation for local optimizer
    if (settings.method == Conflict::Optimizer::local) {
        const double g1 = nonlcon(x, grad, data);
        const double g2 = nonlconDomino(x, grad, data);

        cost += 1e6 * std::max(0.0, g1);
        cost += 1e6 * std::max(0.0, g2);
    }

    return cost;
}


void HeliSpeed::optimize(
    const Params& params,
    const ConflictRow& instance,
    const Conflict& conflict,
    DynamicTraffic& airtraffic,
    const Conflict::OptimSettings& settings)
{
    nlopt::srand(settings.seed);

    // Decision vector x = [deltaV, t0, thold]
    constexpr int n = 3;

    // Create optimizer object
    nlopt::opt opt(nlopt::GN_ISRES, n);
    nlopt::opt local_opt(nlopt::LN_COBYLA, n);

    if (settings.method == Conflict::Optimizer::local) {
        opt = nlopt::opt(nlopt::GN_MLSL, n);
        local_opt = nlopt::opt(nlopt::LN_COBYLA, n);

        local_opt.set_xtol_rel(settings.xtol);
        local_opt.set_ftol_rel(settings.ftol);

        opt.set_local_optimizer(local_opt);
        opt.set_population(conflict.population);
    }

    // Bounds
    // deltaV is allowed to be negative enough to produce hover.
    // Actual commanded ground speed is clamped in speedCommand:
    //      0 <= gs_adv <= params.maxGS
    const double dVMax = params.dVMax;

    opt.set_lower_bounds({
        -dVMax,
        static_cast<double>(onsetTime),
        0.0
    });

    opt.set_upper_bounds({
        dVMax,
        static_cast<double>(instance.t),
        params.tholdMax
    });

    // Input argument structure.
    // Optimizer callbacks mutate optimIn.base, not the original params.
    OptimizerInput optimIn{this, &conflict, &instance, &airtraffic, params, settings};

    // Objective
    opt.set_min_objective(cost, &optimIn);

    // Constraints: g(x) <= 0
    if (settings.method == Conflict::Optimizer::global) {
        opt.add_inequality_constraint(nonlcon,       &optimIn, settings.nonlconTol);
        opt.add_inequality_constraint(nonlconDomino, &optimIn, settings.nonlconTol);
    }

    // Optimizer settings
    std::vector<double> xtol_abs = {1.0, 3.0, 3.0};   // 1 kt, 3 s, 3 s

    opt.set_xtol_abs(xtol_abs);
    opt.set_xtol_rel(settings.xtol);
    opt.set_ftol_rel(settings.ftol);
    opt.set_maxeval(settings.maxFunEval);

    if (std::isfinite(settings.tlimit)) {
        opt.set_maxtime(settings.tlimit);
    }

    // Initial guess
    std::vector<double> x = {
        params.deltaV,
        params.t0,
        params.thold
    };

    // Clamp initial guess into bounds
    x[0] = std::clamp(x[0], -dVMax, dVMax);
    x[1] = std::clamp(x[1], static_cast<double>(onsetTime), static_cast<double>(instance.t));
    x[2] = std::clamp(x[2], 0.0, params.tholdMax);

    double minf = std::numeric_limits<double>::infinity();

    try {
        nlopt::result r = opt.optimize(x, minf);

        state.exitFlag = static_cast<int>(r);

        // Store optimized parameters
        state.optimal = optimIn.base;
        state.optimal.deltaV = x[0];
        state.optimal.t0     = x[1];
        state.optimal.thold  = x[2];

        state.fval = minf;

        std::vector<double> grad;
        state.nonlcon       = nonlcon(x, grad, &optimIn);
        state.nonlconDomino = nonlconDomino(x, grad, &optimIn);

        state.feasible =
            (state.nonlcon <= settings.nonlconTol) &&
            (state.nonlconDomino <= settings.nonlconTol);

    } catch (const std::exception& e) {
        std::cerr << "[opt] exception: " << e.what() << "\n" << std::flush;

        state.optimal = params;
        state.fval = std::numeric_limits<double>::infinity();
        state.feasible = false;
        state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();

    } catch (...) {
        std::cerr << "[opt] unknown exception\n" << std::flush;

        state.optimal = params;
        state.fval = std::numeric_limits<double>::infinity();
        state.feasible = false;
        state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();
    }
}


std::vector<FSSample> HeliSpeed::run(
    const ConflictRow& instance,
    const Conflict& conflict,
    DynamicTraffic& airtraffic,
    const std::filesystem::path& dir,
    ResolutionLogRow& ResolutionLogRow)
{
    // Find intruder's first active time step
    if ((size_t)instance.Iid >= airtraffic.flights.size()) {
        return {};
    }

    const auto& samples = airtraffic.flights[(size_t)instance.Iid].samples;

    auto it = std::find_if(
        samples.begin(),
        samples.end(),
        [](const FSSample& s) { return s.valid; }
    );

    if (it == samples.end()) {
        return {};
    }

    size_t idx = static_cast<size_t>(std::distance(samples.begin(), it));

    // Optimizer settings
    Conflict::OptimSettings settings;
    settings.method = conflict.opt;
    settings.xtol = 1e-2;
    settings.ftol = 1e-3;
    settings.nonlconTol = 1e-3;
    settings.maxFunEval = 30000;
    settings.tlimit = conflict.tlimit;
    settings.seed = 1995;

    // Parameters
    Params params;

    params.tholdMax = static_cast<double>(conflict.pathSize()); // [s]

    // VTOL speed limits
    params.maxGS = 150.0;   // [kts]
    params.dVMax = 200.0;   // [kts], allows slowdown to hover and speed-up

    // Offset rate limits
    params.acc   = 1.0;     // [kts/s]
    params.decel = -1.0;    // [kts/s]

    // Plan traffic risk
    params.trafficInteractionRisk = conflict.trafficRisk(instance, samples, airtraffic);

    // Initial guess
    params.deltaV = 5.0;                               // [kts]
    params.t0     = static_cast<double>(std::max(idx, onsetTime));
    params.thold  = 30.0;                              // [s]

    // Optimize
    auto start = std::chrono::high_resolution_clock::now();

    optimize(params, instance, conflict, airtraffic, settings);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    state.runtime = elapsed.count() * 1000.0;

    if (!state.feasible) {
        return {};
    }

    const auto newTraj = timeShiftTrajectory(
        instance,
        state.optimal,
        airtraffic,
        dir
    );

    if (newTraj.empty()) {
        return {};
    }

    // Log resolution
    ResolutionLogRow.advisory_type = "HeliSpeed";
    ResolutionLogRow.opt_exit_flag = state.exitFlag;
    ResolutionLogRow.onsetTime     = onsetTime;
    ResolutionLogRow.conflictTime  = instance.t;

    // Positive: speed up, negative: slow down/hover
    ResolutionLogRow.V_star  = state.optimal.deltaV;
    ResolutionLogRow.t0_star = state.optimal.t0;
    ResolutionLogRow.dT_star = state.optimal.thold;

    ResolutionLogRow.feasible   = state.feasible;
    ResolutionLogRow.g1_star    = state.nonlcon;
    ResolutionLogRow.g2_star    = state.nonlconDomino;
    ResolutionLogRow.f_star     = state.fval;
    ResolutionLogRow.runtime_ms = state.runtime;
    ResolutionLogRow.exit_flag  = 1;

    ResolutionLogRow.plan_traffic_risk = params.trafficInteractionRisk;
    ResolutionLogRow.adv_traffic_risk  = conflict.trafficRisk(instance, newTraj, airtraffic);

    return newTraj;
}

/*
    SPEED REGULATION
*/
std::pair<double, double> Speed::groundToAirspeed(double groundspeed,
                                         double course,
                                         double windspeed,
                                         double windDirection) const {
    
    std::pair<double, double> out{};

    // Convert to radians
    double chi = course * DEG_2_RAD;

    // Wind direction is "from" convert to "to"
    double wind_to = (windDirection + 180.0) * DEG_2_RAD;

    // Ground velocity vector
    double vg_x = groundspeed * std::sin(chi);  // East component
    double vg_y = groundspeed * std::cos(chi);  // North component

    // Wind velocity vector
    double vw_x = windspeed * std::sin(wind_to); // East component
    double vw_y = windspeed * std::cos(wind_to); // North component

    // Airspeed magnitude
    double va_x = vg_x - vw_x;  // East component
    double va_y = vg_y - vw_y;  // North component
    out.first = std::sqrt(va_x * va_x + va_y * va_y);

    // Heading
    double psi = std::atan2(va_x, va_y) * RAD_2_DEG;  // atan2(E,N)
    out.second = wrapTo360(psi);

    return out;
}

double Speed::airToGroundspeed(double airspeed,
                                         double course,
                                         double windspeed,
                                         double windDirection) const {
    // Convert to radians
    double chi = course * DEG_2_RAD;

    // Wind direction is "from" convert to "to"
    double wind_to = (windDirection + 180.0) * DEG_2_RAD;

    // Ground velocity vector
    double va_x = airspeed * std::sin(chi); // East component
    double va_y = airspeed * std::cos(chi); // North compoennt

    // Wind velocity vector
    double vw_x = windspeed * std::sin(wind_to); // East component
    double vw_y = windspeed * std::cos(wind_to); // North component

    // Airspeed vector
    double vg_x = va_x + vw_x;
    double vg_y = va_y + vw_y;

    return std::sqrt(vg_x * vg_x + vg_y * vg_y);
}

std::vector<FSSample> Speed::speedCommand(const ConflictRow& instance,
                                                    Params& params,
                                                    const DynamicTraffic& airtraffic,
                                                    std::optional<std::filesystem::path> casedir) const{
     
    // Get the planned intruder path
    if ((size_t) instance.Iid >= airtraffic.flights.size()) throw std::invalid_argument("Invalid intruder ID."); 
    const auto& samples = airtraffic.flights[(size_t) instance.Iid].samples;
    if (samples.empty()) throw std::invalid_argument("Intruder tajectory is empty.");
    std::vector<FSSample> TrajShifted(samples.begin(), samples.end());

    // Get the planned ground speed trajectory
    const size_t N = samples.size();

    // Speed change start
    double dv_cmd = 0.0;
    const double deltaV = params.deltaV;

    // Ramp to deltaV
    const size_t iEnd = params.validtf + 1;   // exclusive
    size_t t0 = std::clamp(static_cast<size_t>(params.t0), params.validt0, params.validtf);
    size_t tcatch = iEnd;
    for (size_t i = t0; i < iEnd; ++i){

        // Difference
        double diff = deltaV - dv_cmd;

        // Acceleration/deceleration
        if (diff > 0.0) {
            dv_cmd += std::min(params.acc, diff);
        } else {
            dv_cmd += std::max(params.decel, diff);
        }

        // New ground speed
        TrajShifted[i].gs = samples[i].gs + dv_cmd;
    
        if (std::abs(dv_cmd - deltaV) < 1e-9) {
            tcatch = i + 1;   // next sample starts hold
            break;
        }
    }

    // Speed change hold
    size_t holdEnd = std::min(iEnd, tcatch + static_cast<size_t>(params.thold));
    for (size_t i = tcatch; i < holdEnd; ++i) {
        TrajShifted[i].gs = samples[i].gs + deltaV;
    }

    // Speed change recover
    dv_cmd = deltaV;
    size_t recoveryEnd = iEnd;
    for (size_t i = holdEnd; i < iEnd; ++i) {

        // Difference
        double diff = 0.0 - dv_cmd;

        // Acceleration/deceleration
        if (diff > 0.0) {
            dv_cmd += std::min(params.acc, diff);
        } else {
            dv_cmd += std::max(params.decel, diff);
        }

        // Add to the plan
        TrajShifted[i].gs = samples[i].gs + dv_cmd;

        if (std::abs(dv_cmd) < 1e-9) {
            recoveryEnd = i;
            break;
        }
    }
    if (recoveryEnd == iEnd) return {};

    // Compute the average plan speed within the advisory window
    double sumGS = 0;
    for (size_t i = t0; i < recoveryEnd; i++){
        sumGS += samples[i].gs;
    }
    params.avgGS = sumGS / (recoveryEnd - t0);

    // Airspeed feasibility enforcement
    for (size_t i = t0; i < recoveryEnd; i++) {
        double airspeed = groundToAirspeed(TrajShifted[i].gs,
                                        TrajShifted[i].pos.hdg,
                                        windSpeed,
                                        windDirection).first;

        if (airspeed < params.vaBounds.first || airspeed > params.vaBounds.second) {
            double corrected = std::clamp(airspeed, params.vaBounds.first, params.vaBounds.second);
            TrajShifted[i].gs = corrected;
        }
    }

    // Write to file
    if (casedir){
        std::filesystem::path dir =  casedir.value() / ("groundspeed" + std::to_string(instance.Iid) + ".csv");
        std::ofstream file(dir);
        for (size_t i = 0; i < N; ++i){
            file << samples[i].gs << "," << TrajShifted[i].gs << std::endl ;
        }
        file.close();
    }

    return TrajShifted;
}

void Speed::cumulativeDistance(const double timestep,
                                        const std::vector<FSSample>& samples,
                                        std::vector<double>& cumulDistance) const {
    
    // Get ground speed trajectory array size
    size_t N = samples.size();
    
    // NaN to Zero conversion
    auto v = [](double x){ return std::isnan(x) ? 0.0 : x; };

    // Find cumulative distance at each time step using trapezoidal integration
    for (size_t i = 0; i < N - 1; i++){

        const double v0 = std::max(0.0, v(samples[i].gs));     // Lower bound at zero to avoid negative reported ADS-B speed
        const double v1 = std::max(0.0, v(samples[i+1].gs));   // Lower bound at zero to avoid negative reported ADS-B speed
        const double ds = 0.5 * (v0 + v1) * KTS_2_FTS * FT_2_NM * timestep;

        if (std::isnan(cumulDistance[i])) cumulDistance[i] = 0.0;
        cumulDistance[i+1] = cumulDistance[i] + ds;
    }
}

std::vector<FSSample> Speed::timeShiftTrajectory(const ConflictRow& instance,
                                                            Params& params,
                                                            const DynamicTraffic& airtraffic,
                                                            std::optional<std::filesystem::path> casedir) const {
    double timestep = 1;

    // Get the planned cumulative distance array
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    const size_t N = samples.size();
   
    std::vector<double> plannedCumulDistance(N, 0.0);
    cumulativeDistance(timestep, samples, plannedCumulDistance);
    const double plannedEnd = plannedCumulDistance.back();

    // Check monotonicity
    for (size_t k = 1; k < N; ++k) {
        assert(plannedCumulDistance[k] >= plannedCumulDistance[k-1]);
    }

    // Change the speed profile in-place
    auto TrajShifted = speedCommand(instance,
                                    params,
                                    airtraffic,
                                    casedir);
    if (TrajShifted.empty()) return {};

    std::vector<double> shiftedCumulDistance(N, 0.0);
    cumulativeDistance(timestep, TrajShifted, shiftedCumulDistance);

    // Check monotonicity
    for (size_t k = 1; k < N; ++k) {
        assert(shiftedCumulDistance[k] >= shiftedCumulDistance[k-1]);
    }

    // Change the trajectory
    auto newTraj = samples;
    
    // Shift start time index
    size_t t0 = static_cast<size_t>(params.t0);
    t0 = std::min(t0, N);
    
    // Interpolate trajectory samples after shift
    for (size_t i = t0; i < samples.size(); i++){

        // Get the shifted cumulative distance
        double d = shiftedCumulDistance[i];

        if (d >= plannedEnd){

            // Reached/passed end of planned trajectory
            auto nansample = NaNFSSample();
            nansample.ts = i;
            for (size_t k = i; k < N; ++k) {
                newTraj[k] = nansample;
                newTraj[k].ts = k;
            }
            break;
        }

        // Search for index
        auto it = std::lower_bound(plannedCumulDistance.begin() + 1, plannedCumulDistance.end(), d);
        if (it == plannedCumulDistance.end()) break;
        size_t j = static_cast<size_t>(std::distance(plannedCumulDistance.begin(), it));

        // Ratio
        double r = (d - plannedCumulDistance[j-1])/(plannedCumulDistance[j] - plannedCumulDistance[j-1]);

        // Get the planned trajectory samples
        Pos before = samples[j-1].pos;
        Pos after  = samples[j].pos;
        Pos shifted;
        shifted.lat = before.lat + r * (after.lat - before.lat);
        shifted.lon = before.lon + r * (after.lon - before.lon);
        shifted.alt = before.alt + r * (after.alt - before.alt);
        shifted.hdg = wrapTo360(before.hdg + r * (after.hdg - before.hdg));

        newTraj[i].pos = shifted;
        newTraj[i].ts = i;
        newTraj[i].gs = TrajShifted[i].gs;
    }

    // Log coordinates
	if (casedir){
        std::filesystem::path dir =  casedir.value() / (std::to_string(instance.Iid) + "_trajectory.csv");
        std::ofstream plan(dir);
        for (size_t i = 0; i < samples.size(); ++i) {
            plan << samples[i].pos.lat << ","
                << samples[i].pos.lon << ","
                << samples[i].pos.alt << ","
                << samples[i].pos.hdg << ","
                << samples[i].gs << ","
                << newTraj[i].pos.lat << ","
                << newTraj[i].pos.lon << ","
                << newTraj[i].pos.alt << ","
                << newTraj[i].pos.hdg << ","
                << newTraj[i].gs << "\n";
        }
		plan.close();
	}

    return newTraj;
}

double Speed::nonlcon(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data) {

    // Intentionally unused
    (void) grad;

    // Get optimization data
    auto* d = static_cast<OptimizerInput*>(data);

    auto& params     = d->base;
    auto& airtraffic = *d->airtraffic;
    auto& conflict   = *d->conflict;
    auto& instance   = *d->instance;

    // Update optimization parameters
    params.deltaV = x[0];  // target ground speed change [kts]
    params.t0     = x[1];  // advisory onset time [s]
    params.thold  = x[2];  // hold duration [s]

    // Generate advised intruder trajectory
    const auto newTraj = d->self->timeShiftTrajectory(instance,
                                                      params,
                                                      airtraffic);

    if (newTraj.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    // Ego-speed conflict constraint
    const double g = conflict.ownshipMarginRobustTimeWindow(conflict.emergencyPath,
                                                            newTraj,
                                                            conflict.egoSpeedSlow,
                                                            conflict.egoSpeedNom,
                                                            conflict.egoSpeedFast,
                                                            1);

    if (!std::isfinite(g)) {
        return std::numeric_limits<double>::infinity();
    }

    return g;
}

double Speed::nonlconDomino(const std::vector<double>& x,
                                    std::vector<double>& grad,
                                    void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& conflict   = d->conflict;
    auto& instance   = *d->instance;
    auto& airtraffic = *d->airtraffic;
	params.deltaV    = x[0]; // First optimization parameter, target ground speed change [kts]
	params.t0 	     = x[1]; // Second optimization parameter, initial speed change timestamp [s]
    params.thold 	 = x[2]; // Second optimization parameter, speed change duration [s]
    
	// Shift the trajectory in-place
    const auto shifted = d->self->timeShiftTrajectory(instance, params, airtraffic);

    // Get the new risk value                                                        
    auto metric = conflict->trafficRisk(instance, shifted, airtraffic);
    return std::max(metric.integral - params.trafficInteractionRisk.integral,
                    metric.worst    - params.trafficInteractionRisk.worst);

}

double Speed::cost(const std::vector<double>& x,
                            std::vector<double>& grad,
                            void* data) {
	// Intentioanlly unused arguments
	(void) grad;

    // Get base parameters and pointers
    auto* d         = static_cast<OptimizerInput*>(data);
	auto& params    = d->base;
    auto& settings  = d->settings;
    auto& instance  = *d->instance;

	// Quadratic summation
    auto cost = (x[0]/params.dVMax) * (x[0]/params.dVMax);    // Ground speed deviation penalty
    cost += (x[1])/((double) instance.t) * (x[1])/((double) instance.t); // Late reaction time penalty
    cost += (x[2]/params.tholdMax)*(x[2]/params.tholdMax);  // Speed hold time penalty
    if (settings.method == Conflict::Optimizer::local){
        auto g1 = nonlcon(x, grad, data);
        auto g2 = nonlconDomino(x, grad, data);
        cost += 1e6*std::max(0.0, g1) + 1e6*std::max(0.0, g2);
    }
    return cost;
}

void Speed::optimize(const Params& params,
                                const ConflictRow& instance,
                                const Conflict& conflict,
                                const DynamicTraffic& airtraffic,
                                const Conflict::OptimSettings& settings){

    nlopt::srand(settings.seed);

	// Decision vector x = [deltaV, t0, thold]
    constexpr int n = 3;
	
	// Create optimizer object
    nlopt::opt opt(nlopt::GN_ISRES, n);   // default global optimization
    nlopt::opt local_opt(nlopt::LN_COBYLA, n);
    if (settings.method == Conflict::Optimizer::local){
        opt = nlopt::opt(nlopt::GN_MLSL, n);
        local_opt = nlopt::opt(nlopt::LN_COBYLA, n);
        local_opt.set_xtol_rel(settings.xtol);
        local_opt.set_ftol_rel(settings.ftol);
        opt.set_local_optimizer(local_opt);
        opt.set_population(conflict.population);
    }
	
	// Bounds
    opt.set_lower_bounds({-params.dVMax, (double) onsetTime, 0.0});
    opt.set_upper_bounds({params.dVMax, (double) instance.t, params.tholdMax});
	
	// Input argument structure
	OptimizerInput optimIn{this, &conflict, &instance, &airtraffic, params, settings};

	// Set the objective function
	opt.set_min_objective(cost, &optimIn);

	// Set the inequality constraint: g(x) <= 0
    if (settings.method == Conflict::Optimizer::global){
        opt.add_inequality_constraint(nonlcon, &optimIn, settings.nonlconTol);
        opt.add_inequality_constraint(nonlconDomino, &optimIn, settings.nonlconTol);
    }
	
	// Optimizer setting
    std::vector<double> xtol_abs = {1.0, 3.0, 3.0};   // 1 kt, 3 s
    opt.set_xtol_abs(xtol_abs);
	opt.set_xtol_rel(settings.xtol);
    opt.set_ftol_rel(settings.ftol);
    opt.set_maxeval(settings.maxFunEval);
    if (std::isfinite(settings.tlimit)) opt.set_maxtime(settings.tlimit);

	// Initial guess
    std::vector<double> x = {
        params.deltaV,
        params.t0,
        params.thold
    };

    // Clamp into bounds
    x[0] = std::clamp(x[0], -params.dVMax, params.dVMax);
    x[1] = std::clamp(x[1], (double) onsetTime, (double) instance.t);
    x[2] = std::clamp(x[2], 0.0, params.tholdMax);

	// Initialize the cost value
    double minf = std::numeric_limits<double>::infinity();

	try {
		nlopt::result r = opt.optimize(x, minf);
		state.exitFlag = (int)r;

		// Store
		state.optimal.deltaV 	= x[0];
		state.optimal.t0 	    = x[1];
        state.optimal.thold 	= x[2];
		state.optimal.acc		= params.acc;
		state.optimal.decel		= params.decel;
		state.optimal.tholdMax  = params.tholdMax;
        state.optimal.validt0   = params.validt0;
        state.optimal.validtf   = params.validtf;
        state.optimal.vaBounds  = params.vaBounds;

		state.fval = minf;
		std::vector<double> grad;
        state.nonlcon = nonlcon(x, grad, &optimIn);
        state.nonlconDomino = nonlconDomino(x, grad, &optimIn);
        state.feasible = ((state.nonlcon <= settings.nonlconTol) && state.nonlconDomino <= settings.nonlconTol);

	} catch (const std::exception& e) {
		std::cerr << "[opt] exception: " << e.what() << "\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();

	} catch (...) {
		std::cerr << "[opt] unknown exception\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();
	}
}

std::vector<FSSample> Speed::run(const ConflictRow& instance,
                                        const Conflict& conflict,
                                        const Traffic& traffic,
                                        const std::filesystem::path& dir,
                                        ResolutionLogRow& ResolutionLogRow){

    // Optimizer settings
    Conflict::OptimSettings settings;
    settings.method = conflict.opt;
    settings.xtol = 1e-2;
    settings.ftol = 1e-3;
    settings.nonlconTol = 1e-3;
    settings.maxFunEval = 30000;
    settings.tlimit = conflict.tlimit;
    settings.seed = 1995;

    // Find intruder's first active time step
    Params params;
    const auto& samples = traffic.data.flights[instance.Iid].samples;
    auto it = std::find_if(samples.begin(), samples.end(), [](const FSSample& s){return s.valid;});
    if (it == samples.end()) return {};
    params.validt0 = std::distance(samples.begin(), it);
    auto itr = std::find_if(samples.rbegin(), samples.rend(), [](const FSSample& s){return s.valid;});
    if (itr == samples.rend()) return {};
    params.validtf = std::distance(samples.begin(), itr.base()) - 1;

    // Constant parameters
    params.deltaV      = -30.0; // [kts]
    params.t0          = static_cast<double>(std::max(std::max(onsetTime, params.validt0), (size_t) instance.t));
    params.thold       = 30.0;
    params.acc         = 1.0;  // [kts/s]
    params.decel       = -1.0; // [kts/s]

    // Plan risks
    params.trafficInteractionRisk = conflict.trafficRisk(instance, samples, traffic.data);
    params.ownshipInteractionRisk = conflict.ownshipRisk(conflict.path(), traffic.data);

    // Airspeed limits
    auto type = traffic.type(instance.Iid);
    if (type == AircraftType::GA){    
        params.vaBounds = {70.0, 150.0};
    } else if (type == AircraftType::Turbojet){
        params.vaBounds = {150.0, 350.0};
    }

    // Bounds
    params.dVMax       = 30.0; // [kts]
    params.tholdMax    = static_cast<double>(conflict.pathSize());
    
    // Optimize
    auto start = std::chrono::high_resolution_clock::now();
    optimize(params, instance, conflict, traffic.data, settings);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    state.runtime = elapsed.count() * 1000;
    if (!state.feasible) return {};

    const auto newTraj = timeShiftTrajectory(instance,
                                            state.optimal,
                                            traffic.data,
                                            dir);

    // Log resolution
    ResolutionLogRow.advisory_type  = "Speed";
    ResolutionLogRow.opt_exit_flag  = state.exitFlag;
    ResolutionLogRow.onsetTime      = onsetTime;
    ResolutionLogRow.conflictTime   = instance.t;
    ResolutionLogRow.V_star         = state.optimal.deltaV;
    ResolutionLogRow.t0_star        = state.optimal.t0;
    ResolutionLogRow.dT_star        = state.optimal.thold;
    ResolutionLogRow.feasible       = state.feasible;
    ResolutionLogRow.g1_star        = state.nonlcon;
    ResolutionLogRow.g2_star        = state.nonlconDomino;
    ResolutionLogRow.f_star         = state.fval;
    ResolutionLogRow.runtime_ms     = state.runtime;
    ResolutionLogRow.exit_flag      = 1;
    ResolutionLogRow.plan_traffic_risk = params.trafficInteractionRisk;
    ResolutionLogRow.adv_traffic_risk  = conflict.trafficRisk(instance, newTraj, traffic.data);

    return newTraj;
}

/*
    ALTITUDE REGULATION
*/

double Altitude::verticalRate(const std::vector<FSSample>& samples,
                                        const size_t start,
                                        const size_t deltaT) const {

    assert(deltaT > 0 && "Vertical rate check end time <= start time.");

    // Altitudes [ft]
    auto alt_start = samples[start].pos.alt;
    auto alt_end = samples[start + deltaT].pos.alt;
    if (!std::isfinite(alt_start) || !std::isfinite(alt_end)) return NAN;
    
    // Vertical rate [ft/min]
    auto verticalRate = (alt_end - alt_start) / deltaT * 60;

    return verticalRate;
}


int Altitude::isClimb(const std::vector<FSSample>& samples,
                                const size_t start,
                                const size_t deltaT) const {

    // Vertical rate [ft/min]
    auto vRate = verticalRate(samples, start, deltaT);
    if (!std::isfinite(vRate)) return NAN;

    if (vRate > 50) return 1;
    else if ((vRate >= -50) && (vRate < 50)) return 0;
    else return -1;
}

std::vector<FSSample> Altitude::shiftVerticalProfile(const ConflictRow& instance,
                                                                const Params& params,
                                                                const DynamicTraffic& airtraffic,
                                                                std::optional<std::filesystem::path> casedir) const{
    
    // Conflict time and intruder ID
    const auto ID = instance.Iid;

    // Intruder trajectory
    const auto& samples = airtraffic.flights[ID].samples;
    const auto type     = airtraffic.flights[ID].type;
    const size_t N      = samples.size();
    auto shifted        = samples;

    // Profile shift times
    const size_t twindow = 10;
    size_t shiftStart = std::max(onsetTime, static_cast<size_t>(params.t0));
    assert(shiftStart > 0 && "Altitude shift start time must be positive.");
    if (shiftStart + twindow >= N) return {};

    // Planned and commanded vertical rates [ft/min]
    const double vRatePlan = verticalRate(samples, shiftStart, twindow);
    const double vRateCmd0 = (samples[shiftStart].gs * KTS_2_FTS) * std::sin(params.gamma0 * DEG_2_RAD) * 60.0;   // [ft/min]

    // Relative rate during deviation-build phase
    const double relRate0 = vRateCmd0 - vRatePlan;
    if (!std::isfinite(relRate0) || std::fabs(relRate0) < 1e-9) return {};
    if (params.deltaH * relRate0 <= 0.0) return {};   // must move toward desired deviation

    // Build deviation
    // Grow the altitude offset until |dh| reaches |deltaH|
    size_t holdStart = N;
    double dh = 0.0;

    for (size_t i = shiftStart; i < N; ++i) {
        shifted[i].pos.alt = shifted[i - 1].pos.alt + vRateCmd0 / 60.0;

        dh = shifted[i].pos.alt - samples[i].pos.alt;
        if (!std::isfinite(dh)) return {};

        // Reached target deviation?
        if (std::fabs(dh) >= std::fabs(params.deltaH)) {
            holdStart = i + 1;   // next sample starts hold
            break;
        }
    }
    if (holdStart >= N) return {};

    // Hold
    size_t holdEnd = holdStart + static_cast<size_t>(params.thold);
    if (holdEnd >= N) return {};

    for (size_t i = holdStart; i < holdEnd; ++i) {
        shifted[i].pos.alt = samples[i].pos.alt + params.deltaH;
    }

    // Recovery
    const double vRateCmd1 = (samples[shiftStart].gs * KTS_2_FTS) * std::sin(params.gamma1 * DEG_2_RAD) * 60.0;   // [ft/min]
    const double relRate1 = vRateCmd1 - vRatePlan;
    if (!std::isfinite(relRate1) || std::fabs(relRate1) < 1e-9) return {};

    double dh0 = shifted[holdEnd - 1].pos.alt - samples[holdEnd - 1].pos.alt;
    if (!std::isfinite(dh0)) return {};
    if (dh0 * relRate1 >= 0.0) return {};   // Recovery must reduce deviation

    const double htol = 10.0;
    size_t recoveryEnd = N;
    dh = dh0;
    for (size_t i = holdEnd; i < N; ++i) {
        shifted[i].pos.alt = shifted[i - 1].pos.alt + vRateCmd1 / 60.0;

        dh = shifted[i].pos.alt - samples[i].pos.alt;
        if (!std::isfinite(dh)) return {};

        if (std::fabs(dh) <= htol) {
            recoveryEnd = i + 1;
            break;
        }
    }
    if (recoveryEnd == N) return {};
    if (!std::isfinite(dh) || std::fabs(dh) > htol) return {};

    // Altitude sanity check
    double hmin_h = 300, hmin = 2500; // Min. descent altitude: 300 ft if helicopter, 2500 ft otherwise
    double h0 = (type == AircraftType::Helicopter) ? hmin_h : hmin; 
    for (size_t i = shiftStart; i < recoveryEnd; i++){
        if (shifted[i].pos.alt <= h0) return {};
    }

    // Log
    if (casedir){
        std::filesystem::path dir =  casedir.value() / (std::to_string(instance.Iid) + "_trajectory.csv");
        std::ofstream file(dir);
        for (size_t i = 0; i < N; ++i) {
            file << samples[i].pos.lat << ","
                << samples[i].pos.lon << ","
                << samples[i].pos.alt << ","
                << samples[i].pos.hdg << ","
                << samples[i].gs << ","
                << shifted[i].pos.lat << ","
                << shifted[i].pos.lon << ","
                << shifted[i].pos.alt << ","
                << shifted[i].pos.hdg << ","
                << shifted[i].gs << "\n";
        }
		file.close();
    }

    return shifted;
}

double Altitude::nonlcon(const std::vector<double>& x,
                                    std::vector<double>& grad,
                                    void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& airtraffic = *d->airtraffic;
    auto& conflict   = *d->conflict;
    auto& instance   = *d->instance;
	params.deltaH = x[0];	// First optimization parameter, altitude deviation [ft]
	params.thold  = x[1]; // Second optimization parameter, deviation duration [s]
    params.gamma0 = x[2]; // Third optimization parameter, departure angle [deg]
    params.gamma1 = x[3]; // Fourth optimization parameter, recovery angle [deg]
    params.t0     = x[4]; // Fifth optimization parameter, initial execution timestamp [s]
	
    // Change the altitude profile
    const auto newTraj = d->self->shiftVerticalProfile(instance, params, airtraffic);

    // Ego-speed conflict constraint
    const double g = conflict.ownshipMarginRobustTimeWindow(conflict.emergencyPath,
                                                            newTraj,
                                                            conflict.egoSpeedSlow,
                                                            conflict.egoSpeedNom,
                                                            conflict.egoSpeedFast,
                                                            1);

    if (!std::isfinite(g)) {
        return std::numeric_limits<double>::infinity();
    }
    return g;
}

double Altitude::nonlconDomino(const std::vector<double>& x,
                                        std::vector<double>& grad,
                                        void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto params     = d->base;
    auto conflict   = *d->conflict;
    auto airtraffic = *d->airtraffic;
    auto instance   = *d->instance;
	params.deltaH = x[0]; // First optimization parameter, altitude deviation [ft]
	params.thold  = x[1]; // Second optimization parameter, deviation duration [s]
    params.gamma0 = x[2]; // Third optimization parameter, departure angle [deg]
    params.gamma1 = x[3]; // Fourth optimization parameter, recovery angle [deg]
    params.t0     = x[4]; // Fifth optimization parameter, initial execution timestamp [s]

    // Create a new altitude profile
    const auto shifted = d->self->shiftVerticalProfile(instance, params, airtraffic);
    if (shifted.empty()) return std::numeric_limits<double>::infinity();

    // Get the new risk
	auto metric = conflict.trafficRisk(instance, shifted, airtraffic);
    
    // Return the planned and diversion risk difference
    return std::max(metric.integral - params.trafficInteractionRisk.integral,
                    metric.worst    - params.trafficInteractionRisk.worst);
}

double Altitude::cost(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& settings   = d->settings;
    
	// Quadratic summation
    double J1 = std::pow(x[0]/params.deltaHMax,2);
    double J2 = std::pow(x[1]/params.tholdMax,2);
    double J3 = std::pow((x[2])/(params.gamma0bounds.second - params.gamma0bounds.first),2);
    double J4 = std::pow((x[3])/(params.gamma1bounds.second - params.gamma1bounds.first),2);
    double J5 = std::pow((x[4])/(params.t0bounds.second),2);
    double cost = J1 + J2 + J3 + J4 + J5;

    if (settings.method == Conflict::Optimizer::local){
        double g1 = nonlcon(x, grad, data);
        double g2 = nonlconDomino(x, grad, data);
        cost +=  1e6*std::max(0.0, g1) + 1e6*std::max(0.0, g2);
    }

    return cost;
}

void Altitude::optimize(const Params& params,
                                const ConflictRow& instance,
                                const Conflict& conflict,
                                const DynamicTraffic& airtraffic,
                                const Conflict::OptimSettings& settings){
    nlopt::srand(settings.seed);

	// Decision vector x = [deltaH, thold, gamma0, gamma1, t0]
    constexpr int n = 5;
	
	// Create optimizer object
    nlopt::opt opt(nlopt::GN_ISRES, n);   // default global optimization
    nlopt::opt local_opt(nlopt::LN_COBYLA, n);
    if (settings.method == Conflict::Optimizer::local){
        opt = nlopt::opt(nlopt::GN_MLSL, n);
        local_opt = nlopt::opt(nlopt::LN_COBYLA, n);
        local_opt.set_xtol_rel(settings.xtol);
        local_opt.set_ftol_rel(settings.ftol);
        opt.set_local_optimizer(local_opt);
        opt.set_population(conflict.population);
    }

	// Bounds
    opt.set_lower_bounds({-params.deltaHMax, 0.0, params.gamma0bounds.first, params.gamma1bounds.first, params.t0bounds.first});
    opt.set_upper_bounds({params.deltaHMax, params.tholdMax, params.gamma0bounds.second, params.gamma1bounds.second, params.t0bounds.second});
	
	// Input argument structure
	OptimizerInput optimIn{this, &instance, &conflict, &airtraffic, params, settings};

	// Set the objective function
	opt.set_min_objective(cost, &optimIn);

    // Set the inequality constraint: g(x) <= 0
    if (settings.method == Conflict::Optimizer::global){
        opt.add_inequality_constraint(nonlcon, &optimIn, settings.nonlconTol);
        opt.add_inequality_constraint(nonlconDomino, &optimIn, settings.nonlconTol);
    }
	
	// Optimizer setting
    std::vector<double> xtol_abs = {50.0, 3.0, 0.1, 0.1, 3.0};   // 1 kt, 3 s
    opt.set_xtol_abs(xtol_abs);
	opt.set_xtol_rel(settings.xtol);
    opt.set_ftol_rel(settings.ftol);
    opt.set_maxeval(settings.maxFunEval);
    if (std::isfinite(settings.tlimit)) opt.set_maxtime(settings.tlimit);

	// Initial guess
    std::vector<double> x = {
        params.deltaH,
        params.thold,
        params.gamma0,
        params.gamma1,
        params.t0
    };

    // Clamp into bounds
    x[0] = std::clamp(x[0], -params.deltaHMax, params.deltaHMax);
    x[1] = std::clamp(x[1], 0.0, params.tholdMax);
    x[2] = std::clamp(x[2], params.gamma0bounds.first, params.gamma0bounds.second);
    x[3] = std::clamp(x[3], params.gamma1bounds.first, params.gamma1bounds.second);
    x[4] = std::clamp(x[4], params.t0bounds.first,     params.t0bounds.second);

	// Initialize the cost value
    double minf = std::numeric_limits<double>::infinity();

	try {
		nlopt::result r = opt.optimize(x, minf);
		state.exitFlag = (int)r;

		// Store
        state.optimal           = params;
		state.optimal.deltaH 	= x[0];
		state.optimal.thold 	= x[1];
        state.optimal.gamma0    = x[2];
        state.optimal.gamma1    = x[3];
		state.optimal.t0        = x[4];

		state.fval = minf;
		std::vector<double> grad;
        state.nonlcon = nonlcon(x, grad, &optimIn);
        state.nonlconDomino = nonlconDomino(x, grad, &optimIn);
        state.feasible = (state.nonlcon <= settings.nonlconTol) && (state.nonlconDomino <= settings.nonlconTol);

	} catch (const std::exception& e) {
		std::cerr << "[opt] exception: " << e.what() << "\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();

	} catch (...) {
		std::cerr << "[opt] unknown exception\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();
	}
}

std::vector<FSSample> Altitude::run(const ConflictRow& instance,
                                            const Conflict& conflict,
                                            const DynamicTraffic& airtraffic,
                                            const std::filesystem::path& dir,
                                            ResolutionLogRow& ResolutionLogRow){

    // Planned trajectory
    const auto& samples = airtraffic.flights[instance.Iid].samples;

    auto it = std::find_if(samples.begin(), samples.end(), [](const FSSample& s){return s.valid;});
    auto rit = std::find_if(samples.rbegin(), samples.rend(), [](const FSSample& s){return s.valid;});

    // Optimizer settings
    Conflict::OptimSettings settings;
    settings.method = conflict.opt;
    settings.xtol = 1e-2;
    settings.ftol = 1e-3;
    settings.nonlconTol = 1e-3;
    settings.maxFunEval = 30000;
    settings.tlimit = conflict.tlimit;
    settings.seed = 1995;

    // Intruder type
    auto type = airtraffic.flights[instance.Iid].type;

    // Constant parameters and initial guess
    Params params;
    params.deltaHMax = 1000; // [ft]
    params.tholdMax = static_cast<double>(conflict.pathSize()); // [s]
    params.t0       = static_cast<double>(onsetTime);    // s
    params.idx0     = std::distance(samples.begin(), it);
    params.idxf     = std::distance(samples.begin(), rit.base()) - 1;
    params.t0bounds = {static_cast<double>(std::max(onsetTime, params.idx0)), static_cast<double>(instance.t)};
    
    if (type == AircraftType::Helicopter){
        params.gamma0bounds = {-15.0, 15.0}; // [deg]
        params.gamma1bounds = {-15.0, 15.0}; // [deg]
        params.gamma0 = 5.0;
        params.gamma1 = -5.0;
        params.deltaH = 500.0;
        params.thold  = 60.0;
    } else {
        params.gamma0bounds = {-5.0, 3.0}; // [deg]
        params.gamma1bounds = {-5.0, 3.0}; // [deg]
        params.gamma0 = 3.0;
        params.gamma1 = -3.0;
        params.deltaH = 500.0;
        params.thold  = 60.0;
    }

    // Plan risk metrics
    params.trafficInteractionRisk = conflict.trafficRisk(instance, samples, airtraffic);  
    params.ownshipInteractionRisk = conflict.ownshipRisk(conflict.path(), airtraffic);  

    // Optimize
    auto start = std::chrono::high_resolution_clock::now();
    optimize(params, instance, conflict, airtraffic, settings);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    state.runtime = elapsed.count() * 1000;
    if (!state.feasible) return {};

    // Out
    const auto newTraj = shiftVerticalProfile(instance, state.optimal, airtraffic, dir);

    // Log resolution
    ResolutionLogRow.advisory_type  = "Altitude";
    ResolutionLogRow.opt_exit_flag  = state.exitFlag;
    ResolutionLogRow.onsetTime      = onsetTime;
    ResolutionLogRow.conflictTime   = instance.t;
    ResolutionLogRow.dH_star        = state.optimal.deltaH;
    ResolutionLogRow.T_star         = state.optimal.thold;
    ResolutionLogRow.t0_star        = state.optimal.t0;
    ResolutionLogRow.gamma0_star    = state.optimal.gamma0;
    ResolutionLogRow.gamma1_star    = state.optimal.gamma1;
    ResolutionLogRow.feasible       = state.feasible;
    ResolutionLogRow.g1_star        = state.nonlcon;
    ResolutionLogRow.g2_star        = state.nonlconDomino;
    ResolutionLogRow.f_star         = state.fval;
    ResolutionLogRow.runtime_ms     = state.runtime;
    ResolutionLogRow.exit_flag      = 2;
    ResolutionLogRow.plan_traffic_risk = params.trafficInteractionRisk;
    ResolutionLogRow.adv_traffic_risk  = conflict.trafficRisk(instance, newTraj, airtraffic);

    return newTraj;
}


/*
    LATERAL DIVERSION
*/
Pos Extend::interpPos(const std::pair<Pos,double>& x1,
                                const std::pair<Pos,double>& x2,
                                double t) const {
    const Pos& p1 = x1.first;
    const Pos& p2 = x2.first;
    double t1 = x1.second;
    double t2 = x2.second;

    if (t2 == t1) return p1;

    double alpha = (t - t1) / (t2 - t1);

    Pos out;

    out.lat = p1.lat + alpha * (p2.lat - p1.lat);
    out.lon = p1.lon + alpha * (p2.lon - p1.lon);
    out.alt = p1.alt + alpha * (p2.alt - p1.alt);

    // Heading wrap-safe interpolation
    double dh = wrapTo360(p2.hdg - p1.hdg);
    if (dh > 180.0) dh -= 360.0;

    out.hdg = wrapTo360(p1.hdg + alpha * dh);

    return out;
}

void Extend::buildSTurnPath(DubinsPath& dubins,
                                      const Params& params,
                                      const Conflict& conflict){


    shortestDubins(&dubins, &conflict.dubinsopt);

    // Set the initial intermediate waypoint along the straight segment
    struct Pos interWaypoint0, interWaypoint;
    interWaypoint0.lat = 0.5*(dubins.traj[1].wpt.pos.lat + dubins.traj[2].wpt.pos.lat);
    interWaypoint0.lon = 0.5*(dubins.traj[1].wpt.pos.lon + dubins.traj[2].wpt.pos.lon);
    interWaypoint0.alt = 0.5*(dubins.traj[1].wpt.pos.alt + dubins.traj[2].wpt.pos.alt);
    interWaypoint0.hdg = dubins.traj[1].wpt.hdg;

    // Get the initial straight segment length
    double straightLength, course;
    geo_dist(&dubins.traj[1].wpt.pos, &dubins.traj[2].wpt.pos, &straightLength, &course, &conflict.geoopt);

    // Set theta
    double theta = std::fabs(params.theta); // [deg]
    
    // Get an intermediate waypoint
    int extendTo = (params.theta > 0) ? 1 : -1;
    auto geoopt = conflict.geoopt;
    getIntermediateWaypoint(&interWaypoint0, &interWaypoint, straightLength, theta, extendTo, &geoopt);

    // Create two Dubins path structures to store paths
    for (size_t i = 0; i < 2; i++) {
        Traj_InitArray(path[i].traj, 4);   // Allocates memory for trajectory structure
        path[i].size = 2;                  // Sets path structure size to 2, indicating the structure holds two Dubins paths (S-Turn)
    }

    // Copy the original initial and final waypoints to the Sturn path structures
    Traj_CopyAll(dubins.traj, path[0].traj);
    path[0].traj[3].wpt.pos = interWaypoint;
    path[0].traj[3].wpt.pos.alt = dubins.traj[0].wpt.pos.alt;
    path[0].traj[3].wpt.hdg = interWaypoint.hdg;

    Traj_CopyAll(dubins.traj, path[1].traj);
    path[1].traj[0].wpt.pos = interWaypoint;
    path[1].traj[0].wpt.pos.alt = dubins.traj[0].wpt.pos.alt;
    path[1].traj[0].wpt.hdg = interWaypoint.hdg;

    // leg 2: intermediate -> final
    path[1].traj[0].wpt.pos = interWaypoint;
    path[1].traj[0].wpt.hdg = interWaypoint.hdg;
    path[1].traj[3] = dubins.traj[3];

    // Build the path
    shortestDubins(&path[0], &conflict.dubinsopt);
    shortestDubins(&path[1], &conflict.dubinsopt);
}

double Extend::mapTheta(double z)
{
    z = std::clamp(z, -1.0, 1.0);

    const double thetaL_min = -180.0;
    const double thetaL_max = -120.0;
    const double thetaR_min =  120.0;
    const double thetaR_max =  180.0;

    if (z < 0.0) {
        // z = 0  -> -180 deg
        // z = -1 -> -120 deg
        return thetaL_min + std::abs(z) * (thetaL_max - thetaL_min);
    } else {
        // z = 0 -> 180 deg
        // z = 1 -> 120 deg
        return thetaR_max - z * (thetaR_max - thetaR_min);
    }
}


void Extend::generateSTurnPath(const ConflictRow& instance,
                                        const Params& params,
                                        const Conflict& conflict,
                                        const DynamicTraffic& airtraffic){
    
    // Intruder planned trajectory
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    const size_t N = samples.size();
    if (N == 0) throw std::runtime_error("Empty intruder trajectory.");

    // Initial waypoint of the diversion
    size_t tstart = static_cast<size_t>(std::lround(params.t0));
    size_t tfinal = static_cast<size_t>(std::lround(params.tf));

    // Turn radius
    Traffic traffic;
    double R;
    traffic.turnRadius(instance.Iid, tstart, R, airtraffic);

    // Initialize the path structure
    DubinsPath dubins;
    Traj_InitArray(dubins.traj, 4);
    dubins.traj[0].wpt.pos = samples[tstart].pos;
    dubins.traj[0].wpt.hdg = samples[tstart].pos.hdg;
    
    dubins.traj[3].wpt.pos = samples[tfinal].pos;
    dubins.traj[3].wpt.hdg = samples[tfinal].pos.hdg;
    
    dubins.traj[0].wpt.rad = R * FT_2_NM;
    dubins.traj[2].wpt.rad = R * FT_2_NM;

    // Initialize zero flight path angle
    for (size_t i = 0; i < 4; i++){
        dubins.traj[i].wpt.gam = 0;
    }

    // Compute S-Turn path
    buildSTurnPath(dubins, params, conflict);
}

std::vector<FSSample> Extend::DubinsCoordinates(const Conflict& conflict,
                                                         const double groundspeed){

    // Get S-Turn coordinates
    int sampleSize1, sampleSize2;
    Pos *sturn1 = getDubinsCoordinatesWithFixedTimeStep(&path[0], nullptr, &sampleSize1, &(conflict.geoopt), groundspeed); // Initial waypoint to the intermediate waypoint
    if (!sturn1 || sampleSize1 <= 0) return {};
    Pos *sturn2 = getDubinsCoordinatesWithFixedTimeStep(&path[1], &sturn1[sampleSize1-1], &sampleSize2, &(conflict.geoopt), groundspeed); // Intermediate waypoint to the final waypoint
    if (!sturn2 || sampleSize2 <= 0) {
        free(sturn1);
        return {};
    }

    // Concatenate samples
    size_t size1 = static_cast<size_t>(sampleSize1);
    size_t size2 = static_cast<size_t>(sampleSize2);
    size_t numSamples = size1 + size2 - 1;
    std::vector<FSSample> coordinates(numSamples, FSSample{});

    for (size_t i = 0; i < numSamples; i++) {
        if (i < size1) {
            coordinates[i].pos = sturn1[i];
        } else {
            coordinates[i].pos = sturn2[i - sampleSize1 + 1];
        }
        coordinates[i].gs = groundspeed;
    }
    free(sturn1); sturn1 = nullptr;
    free(sturn2); sturn2 = nullptr;

    return coordinates;
}

std::vector<FSSample> Extend::divert(const ConflictRow& instance,
                                                const Params& params,
                                                const Conflict& conflict,
                                                const DynamicTraffic& airtraffic,
                                                std::optional<std::filesystem::path> casedir) {

    // Compute diversion path
    generateSTurnPath(instance,
                    params,
                    conflict,
                    airtraffic);

    // Planned trajectory
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    if (samples.empty()) return {};

    size_t tstart = static_cast<size_t>(std::lround(params.t0));
    size_t tfinal = static_cast<size_t>(std::lround(params.tf));

    tstart = std::clamp(tstart, static_cast<size_t>(0), samples.size() - 1);
    tfinal = std::clamp(tfinal, static_cast<size_t>(0), samples.size() - 1);
                
    // Get S-Turn coordinates
    const auto gs = airtraffic.flights[instance.Iid].samples[tstart].gs;
    std::vector<FSSample> coordinates = DubinsCoordinates(conflict, gs);
    if (coordinates.empty()) return {};

    // Sanity check, each Dubins coordinate must be valid.
    auto it = std::find_if(coordinates.begin(), coordinates.end(),
                        [](const FSSample& s){return !std::isfinite(s.pos.lat);});
    if (it != coordinates.end()) return {};
    
    // Enforce prohibited zone avoidance
    it = std::find_if(coordinates.begin(), coordinates.end(),
                            [&conflict](const FSSample& s){
                                double tmp = getAirTrafficDensity(&s.pos, &conflict.prohibited);
                                return tmp > 0.5;});
    if (it != coordinates.end()) return {};

    // Store the new trajectory
    std::vector<FSSample> newTraj;
    newTraj.reserve(tstart + coordinates.size() + (samples.size() - (tfinal + 1)));

    auto itStart = samples.begin() + static_cast<std::ptrdiff_t>(tstart);
    auto itFinal = samples.begin() + static_cast<std::ptrdiff_t>(tfinal + 1);

    newTraj.insert(newTraj.end(), samples.begin(), itStart);
    newTraj.insert(newTraj.end(), coordinates.begin(), coordinates.end());
    newTraj.insert(newTraj.end(), itFinal, samples.end());

    // Project planned altitude onto diversion segment
    const double alt0 = samples[tstart].pos.alt;
    const double alt1 = samples[tfinal].pos.alt;

    for (size_t k = 0; k < coordinates.size(); ++k) {
        double alpha = (coordinates.size() == 1) ? 0.0 : static_cast<double>(k) / static_cast<double>(coordinates.size() - 1);
        newTraj[tstart + k].pos.alt = alt0 + alpha * (alt1 - alt0);
        newTraj[tstart + k].gs = gs;
    }

    // Log
    if (casedir){
        std::filesystem::path dir =  casedir.value() / (std::to_string(instance.Iid) + "_trajectory.csv");
        std::ofstream log(dir);
        auto N = samples.size();

        // Entire trajectory log: planned and diversion
        for (size_t i = 0; i < newTraj.size(); i++){

            if (i < N)
                log << samples[i].pos.lat << ", "
                    << samples[i].pos.lon << ", "
                    << samples[i].pos.alt << ", "
                    << samples[i].pos.hdg << ", "
                    << samples[i].gs << ", "
                    << newTraj[i].pos.lat << ", "
                    << newTraj[i].pos.lon << ", "
                    << newTraj[i].pos.alt << ", "
                    << newTraj[i].pos.hdg << ", "
                    << newTraj[i].gs << "\n";
            else 
                log << NAN << ", "
                    << NAN << ", "
                    << NAN << ", "
                    << NAN << ", "
                    << NAN << ", "
                    << newTraj[i].pos.lat << ", "
                    << newTraj[i].pos.lon << ", "
                    << newTraj[i].pos.alt << ", "
                    << newTraj[i].pos.hdg << ", "
                    << newTraj[i].gs << "\n";
        }
        log.close();

        // Diversion segment log: only the S-Turn segment
        dir =  casedir.value() / (std::to_string(instance.Iid) + "_diversion.csv");
        std::ofstream log_div(dir);
        for (size_t i = 0; i < 2; i++){
            for (size_t j = 0; j < 4; j++){
                log_div << path[i].traj[j].wpt.pos.lat << ", "
                        << path[i].traj[j].wpt.pos.lon << ", "
                        << path[i].traj[j].wpt.pos.alt << ", "
                        << path[i].traj[j].wpt.hdg << ", "
                        << path[i].traj[j].wpt.rad << ", "
                        << path[i].traj[j].wpt.gam << ", "
                        << path[i].type << "\n";
            }
        }
        log_div.close();
    }
    return newTraj;
}

double Extend::nonlcon(const std::vector<double>& x,
                                std::vector<double>& grad,
                                void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& conflict   = *d->conflict;
    auto& airtraffic = *d->airtraffic;
    auto& instance   = *d->instance;
    auto* self      = d->self;
    params.aux      = x[0];
	params.theta    = self->mapTheta(x[0]);
    params.t0       = x[1]; // Second optimization parameter, Dubins initial waypoint timestamp [s]
    params.tf       = x[2]; // Third optimization parameter, Dubins final waypoint timestap [s]

    size_t tstart = static_cast<size_t>(std::lround(params.t0));
    size_t tfinal = static_cast<size_t>(std::lround(params.tf));
    
    // Get the initial conflict instance
    const auto newTraj = self->divert(instance,
                                    params,
                                    conflict,
                                    airtraffic);

    if (newTraj.empty()) return std::numeric_limits<double>::infinity();

    // Enforce Minimum Safe Altitude to Dubins waypoints
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    double galt = std::max(MSA - samples[tstart].pos.alt,
                           MSA - samples[tfinal].pos.alt);

    // Enforce flight path angle
    double hdist = self->path[0].hdist + self->path[1].hdist;
    double dh    = samples[tstart].pos.alt - samples[tfinal].pos.alt;
    double gamma = std::atan(dh/(hdist*NM_2_FT)) * RAD_2_DEG;
    double ggam  = std::abs(gamma) - 3.0;

    // Robust ego-speed conflict constraint
    const double g = conflict.ownshipMarginRobustTimeWindow(conflict.emergencyPath,
                                                            newTraj,
                                                            conflict.egoSpeedSlow,
                                                            conflict.egoSpeedNom,
                                                            conflict.egoSpeedFast,
                                                            1);

    if (!std::isfinite(g)) {
        return std::numeric_limits<double>::infinity();
    }

	// NLopt enforces g(x) <= 0
    return std::max(g, std::max(galt, ggam));
}

double Extend::nonlconDomino(const std::vector<double>& x,
                                        std::vector<double>& grad,
                                        void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& conflict   = *d->conflict;
    auto& airtraffic = *d->airtraffic;
    auto& instance   = *d->instance;
    auto* self      = d->self;
    params.aux      = x[0];
	params.theta    = self->mapTheta(x[0]);	// First optimization parameter, S-Turn theta angle [deg]
    params.t0       = x[1]; // Second optimization param eter, Dubins initial waypoint timestap [s]
    params.tf       = x[2]; // Third optimization parameter, Dubins final waypoint timestap [s]

    // Divert
    const auto newTraj = self->divert(instance,
                                    params,
                                    conflict,
                                    airtraffic);
    if (newTraj.empty()) return std::numeric_limits<double>::infinity();

    // Get the new risk
	auto metric = conflict.trafficRisk(instance, newTraj, airtraffic);
    
    // Return the planned and diversion risk difference
    return std::max(metric.integral - params.trafficInteractionRisk.integral,
                    metric.worst    - params.trafficInteractionRisk.worst);
}

double Extend::cost(const std::vector<double>& x,
                            std::vector<double>& grad,
                            void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization settings
    auto& settings = d->settings;
    auto& conflict = d->conflict;

	// Quadratic summation
    double T = static_cast<double>(conflict->pathSize());
    double cost  = x[0]*x[0];
    cost += std::max(0.0, std::fabs(x[2] - x[1])/T);

    if (settings.method == Conflict::Optimizer::local){
        auto g1 = nonlcon(x, grad, data);
        auto g2 = nonlconDomino(x, grad, data);
        cost += 1e6*std::max(0.0, g1) + 1e6*std::max(0.0, g2);
    }

    return cost;
}

void Extend::optimize(const Params& params,
                                const ConflictRow& instance,
                                const Conflict& conflict,
                                const DynamicTraffic& airtraffic,
                                const Conflict::OptimSettings& settings){
    nlopt::srand(settings.seed);

	// Decision vector x = [aux_theta, t0, tf]
    constexpr int n = 3;
	
	// Create optimizer object
    nlopt::opt opt(nlopt::GN_ISRES, n);   // default global optimization
    nlopt::opt local_opt(nlopt::LN_COBYLA, n);
    if (settings.method == Conflict::Optimizer::local){
        opt = nlopt::opt(nlopt::GN_MLSL, n);
        local_opt = nlopt::opt(nlopt::LN_COBYLA, n);
        local_opt.set_xtol_rel(settings.xtol);
        local_opt.set_ftol_rel(settings.ftol);
        opt.set_local_optimizer(local_opt);
        opt.set_population(conflict.population);
    }

	// Bounds
    const size_t N = airtraffic.flights[instance.Iid].samples.size();
    opt.set_lower_bounds({0.0, (double) onsetTime, 0.0});
    opt.set_upper_bounds({1.0,  (double) instance.t, (double) (N-1)});

	// Input argument structure
	OptimizerInput optimIn{this, &instance, &conflict, &airtraffic, params, settings};

	// Set the objective function
	opt.set_min_objective(cost, &optimIn);

	// Set the inequality constraint: g(x) <= 0
    if (settings.method == Conflict::Optimizer::global){
        opt.add_inequality_constraint(nonlcon, &optimIn, settings.nonlconTol);
        opt.add_inequality_constraint(nonlconDomino, &optimIn, settings.nonlconTol);
    }
	
	// Optimizer setting
    std::vector<double> xtol_abs = {1e-2, 3.0, 3.0};
    opt.set_xtol_abs(xtol_abs);
	opt.set_xtol_rel(settings.xtol);
    opt.set_ftol_rel(settings.ftol);
    opt.set_maxeval(settings.maxFunEval);
    if (std::isfinite(settings.tlimit)) opt.set_maxtime(settings.tlimit);

	// Initial guess
    std::vector<double> x = {params.theta, params.t0, params.tf};

    // Clamp into bounds
    x[0] = std::clamp(x[0], -1.0, 1.0);
    x[1] = std::clamp(x[1], (double) onsetTime, (double) instance.t);
    x[2] = std::clamp(x[2], 0.0, (double) (N-1));

	// Initialize the cost value
    double minf = std::numeric_limits<double>::infinity();

	try {
		nlopt::result r = opt.optimize(x, minf);
		state.exitFlag = (int)r;

		// Store
        state.optimal       = params;
		state.optimal.aux 	= x[0];
        state.optimal.theta = mapTheta(x[0]);
        state.optimal.t0 	= x[1];
		state.optimal.tf    = x[2];
        state.optimal.length = path[0].hdist + path[1].hdist;

		state.fval = minf;
		std::vector<double> grad;
        state.nonlcon = nonlcon(x, grad, &optimIn);
        state.nonlconDomino = nonlconDomino(x, grad, &optimIn);
        state.feasible = (state.nonlcon <= settings.nonlconTol) && ((state.nonlconDomino <= settings.nonlconTol));

	} catch (const std::exception& e) {
		std::cerr << "[opt] exception: " << e.what() << "\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();

	} catch (...) {
		std::cerr << "[opt] unknown exception\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();
	}
}

std::vector<FSSample> Extend::run(const ConflictRow& instance,
                                            const Conflict& conflict,
                                            const DynamicTraffic& airtraffic,
                                            const std::filesystem::path& dir,
                                            ResolutionLogRow& ResolutionLogRow){

    // Find intruders active airborne time indices
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    const size_t N = samples.size();

    // Define initial parameter set
    Params params;
    params.aux         = 0.3;
    params.t0          = 0.0;// static_cast<double>(instance.t);
    params.tf          = static_cast<double>(std::min(instance.t - 120, (int) N-1));
    params.tmax        = static_cast<double>(N);

    // Plan risk metrics
    params.trafficInteractionRisk = conflict.trafficRisk(instance, samples, airtraffic);
    params.ownshipInteractionRisk = conflict.ownshipRisk(conflict.path(), airtraffic);

    // Optimizer settings
    Conflict::OptimSettings settings;
    settings.method = conflict.opt;
    settings.xtol = 1e-2;
    settings.ftol = 1e-3;
    settings.nonlconTol = 1e-3;
    settings.maxFunEval = 30000;
    settings.tlimit = conflict.tlimit;
    settings.seed = 1995;

    // Optimize
    auto start = std::chrono::high_resolution_clock::now();
    optimize(params, instance, conflict, airtraffic, settings);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    state.runtime = elapsed.count() * 1000;
    if (!state.feasible) return {};
    
    // Out
    const auto newTraj  = divert(instance, state.optimal, conflict, airtraffic, dir);

    // Log resolution
    ResolutionLogRow.advisory_type  = "Diversion";
    ResolutionLogRow.opt_exit_flag  = state.exitFlag;
    ResolutionLogRow.onsetTime      = onsetTime;
    ResolutionLogRow.conflictTime   = instance.t;
    ResolutionLogRow.theta_star     = state.optimal.theta;
    ResolutionLogRow.t0_star        = state.optimal.t0;
    ResolutionLogRow.tf_star        = state.optimal.tf;
    ResolutionLogRow.length_star    = state.optimal.length;
    ResolutionLogRow.feasible       = state.feasible;
    ResolutionLogRow.g1_star        = state.nonlcon;
    ResolutionLogRow.g2_star        = state.nonlconDomino;
    ResolutionLogRow.f_star         = state.fval;
    ResolutionLogRow.runtime_ms     = state.runtime;
    ResolutionLogRow.exit_flag      = 2;
    ResolutionLogRow.plan_traffic_risk = params.trafficInteractionRisk;
    ResolutionLogRow.adv_traffic_risk  = conflict.trafficRisk(instance, newTraj, airtraffic);

    return newTraj;
}

/*
    Diversion
*/
void Divert::computePath(const ConflictRow& instance,
                        const Params& params,
                        const Conflict& conflict,
                        const DynamicTraffic& airtraffic){

    // Intruder trajectory
    const auto& samples = airtraffic.flights[instance.Iid].samples;

    // Diversion time
    size_t tstart = static_cast<size_t>(params.t0);

    // Initial waypoint
    auto init = samples[tstart].pos;

    // Final waypoint
    Pos hold = Holds[params.hold_idx].pos;

    // Turn radius
    Traffic traffic;
    double R;
    traffic.turnRadius(instance.Iid, tstart, R, airtraffic);
    R *= FT_2_NM;

    // Dubins path
    Traj_InitArray(path.traj, 4);
    path.traj[0].wpt.pos = init;
    path.traj[0].wpt.hdg = init.hdg;
    path.traj[3].wpt.pos = hold;
    path.traj[3].wpt.hdg = hold.hdg;

    path.traj[0].wpt.rad = R;
    path.traj[2].wpt.rad = R;

    for (size_t i = 0; i < 4; i++){
        path.traj[i].wpt.gam = params.gamma;
    }

    // Determine the path type
    double bearing, d;
    const FSSample ownship = conflict.pathInstance(instance.t);
    geo_dist(&init, &ownship.pos, &d, &bearing, &conflict.geoopt);
    double rel = bearing - init.hdg;
    while (rel > 180.0) rel -= 360.0;
    while (rel <= -180.0) rel += 360.0;
    int type = 0;
    if (rel > 0){
        type = 3;   // LSL
    } else {
        type = 0;   // RSR
    }

    // shortestDubins(&path, &conflict.dubinsopt);
    getDubinsWithType(&path, type, &conflict.dubinsopt);
    path.hdist = 0;
    for (size_t i = 0; i < 3; i++){
        path.hdist += path.traj[i].hdist;
    }
}

std::vector<FSSample> Divert::pathToHold(const ConflictRow& instance,
                                        const Params& params,
                                        const Conflict& conflict,
                                        const DynamicTraffic& airtraffic,
                                        std::optional<std::filesystem::path> casedir){

    // Intruder's planned trajectory
    const auto& samples = airtraffic.flights[instance.Iid].samples;

    // Compute Dubins path to hold
    computePath(instance,
                params,
                conflict,
                airtraffic);
                
    // Get Dubins coordinates
    size_t tstart = static_cast<size_t>(params.t0);
    const double groundspeed = airtraffic.flights[instance.Iid].samples[tstart].gs;
    int sampleSize;
    Pos *coordinates = getDubinsCoordinatesWithFixedTimeStep(&path, nullptr, &sampleSize, &(conflict.geoopt), groundspeed);
    if (sampleSize <= 0 || !std::isfinite(sampleSize)){
        throw std::runtime_error("Divert: Coordinate sample size must be positive.");
    }
    
    // Store Dubins coordinates in a trajectory array
    size_t N = static_cast<size_t>(sampleSize);
    std::vector<FSSample> dubins(N, FSSample{});
    for (size_t i = 0; i < N; i++){
        dubins[i].pos = coordinates[i];
        dubins[i].gs = groundspeed;
        dubins[i].ts = tstart + i;
    }
    free(coordinates); coordinates = nullptr;

    // Sanity check, each Dubins coordinate must be valid.
    auto it = std::find_if(dubins.begin(), dubins.end(),
                        [](const FSSample& s){return !std::isfinite(s.pos.lat);});
    if (it != dubins.end()) return {};

    // Enforce prohibited zone avoidance
    it = std::find_if(dubins.begin(), dubins.end(),
                            [&conflict](const FSSample& s){
                            double tmp = getAirTrafficDensity(&s.pos, &conflict.prohibited);
                            return tmp > 0.5;});
    if (it != dubins.end()) return {};

    // Flatten altitude profile at hold altitude
    it = std::find_if(dubins.begin(), dubins.end(),
                        [params](const FSSample& s){return s.pos.alt >= Holds[params.hold_idx].pos.alt;});
    if (it != dubins.end()){
        size_t idx = std::distance(dubins.begin(), it);
        for (size_t i = idx; i < N; i++){
            dubins[i].pos.alt = Holds[params.hold_idx].pos.alt;
        }
    }

    // Store the new trajectory
    std::vector<FSSample> newTraj;
    newTraj.reserve(tstart + N - 1);
    auto itStart = samples.begin() + static_cast<std::ptrdiff_t>(tstart);
    newTraj.insert(newTraj.end(), samples.begin(), itStart);
    newTraj.insert(newTraj.end(), dubins.begin(), dubins.end());

    // Log
    if (casedir){
        std::filesystem::path dir =  casedir.value() / (std::to_string(instance.Iid) + "_trajectory.csv");
        std::ofstream log(dir);
        auto N = samples.size();
        for (size_t i = 0; i < newTraj.size(); i++){

            if (i < N)
                log << samples[i].pos.lat << ", "
                    << samples[i].pos.lon << ", "
                    << samples[i].pos.alt << ", "
                    << samples[i].pos.hdg << ", "
                    << samples[i].gs << ", "
                    << newTraj[i].pos.lat << ", "
                    << newTraj[i].pos.lon << ", "
                    << newTraj[i].pos.alt << ", "
                    << newTraj[i].pos.hdg << ", "
                    << newTraj[i].gs << "\n";

            else 
                log << NAN << ", "
                    << NAN << ", "
                    << NAN << ", "
                    << NAN << ", "
                    << NAN << ", "
                    << newTraj[i].pos.lat << ", "
                    << newTraj[i].pos.lon << ", "
                    << newTraj[i].pos.alt << ", "
                    << newTraj[i].pos.hdg << ", "
                    << newTraj[i].gs << "\n";
        }
        log.close();
    }
    return newTraj;
}

double Divert::nonlcon(const std::vector<double>& x,
                        std::vector<double>& grad,
                        void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& conflict   = *d->conflict;
    auto& airtraffic = *d->airtraffic;
    auto& instance   = *d->instance;
    auto* self       = d->self;
	params.gamma    = x[0];	 // First optimization parameter, flight path angle [deg]
    params.t0       = x[1];  // Second optimization parameter, initial advisory execution timestamp [s]
    
    // Get the initial conflict instance
    const auto newTraj = self->pathToHold(instance, params, conflict, airtraffic);
    if (newTraj.empty()) return std::numeric_limits<double>::infinity();

    // Ego-speed conflict constraint
    const double g = conflict.ownshipMarginRobustTimeWindow(conflict.emergencyPath,
                                                            newTraj,
                                                            conflict.egoSpeedSlow,
                                                            conflict.egoSpeedNom,
                                                            conflict.egoSpeedFast,
                                                            1);

    if (!std::isfinite(g)) {
        return std::numeric_limits<double>::infinity();
    }
    return g;
}

double Divert::nonlconDomino(const std::vector<double>& x,
                            std::vector<double>& grad,
                            void* data){
    
	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto params     = d->base;
    auto conflict   = *d->conflict;
    auto airtraffic = *d->airtraffic;
    auto instance   = *d->instance;
    auto* self      = d->self;
	params.gamma    = x[0];	 // First optimization parameter, flight path angle [deg]
    params.t0       = x[1];  // Second optimization parameter, initial advisory execution timestamp [s]

    // Divert
    const auto newTraj = self->pathToHold(instance,
                                    params,
                                    conflict,
                                    airtraffic);

    if (newTraj.empty()) return std::numeric_limits<double>::infinity();

    // Get the new risk
	auto metric = conflict.trafficRisk(instance, newTraj, airtraffic);
    
    // Return the planned and diversion risk difference
    return std::max(metric.integral - params.trafficInteractionRisk.integral,
                    metric.worst    - params.trafficInteractionRisk.worst);
}

double Divert::cost(const std::vector<double>& x,
                    std::vector<double>& grad,
                    void* data){

	// Intentioanlly unused arguments
	(void) grad;
   
	// Get base parameters and pointers
    auto* d = static_cast<OptimizerInput*>(data);
    
    // Get optimization parameters
    auto& params     = d->base;
    auto& settings   = d->settings;

    // Cost
    double cost1 = std::max(0.0, std::fabs(x[0]) - 3.0) / 3.0;
    double cost2 = params.t0/params.t0Bounds.second;
    if (settings.method == Conflict::Optimizer::local){
        auto g1 = nonlcon(x, grad, data);
        auto g2 = nonlconDomino(x, grad, data);
        return cost1 + cost2 + 1e6*std::max(0.0,g1) + 1e6*std::max(0.0,g2);
    } else {
        return cost1*cost1 + cost2*cost2;
    }
}

void Divert::optimize(const Params& params,
                        const ConflictRow& instance,
                        const Conflict& conflict,
                        const DynamicTraffic& airtraffic,
                        const Conflict::OptimSettings& settings){
    nlopt::srand(settings.seed);

	// Optimization vector x = [gamma, t0]
    constexpr int n = 2;
	
	// Create optimizer object
    nlopt::opt opt(nlopt::GN_ISRES, n);   // default global optimization
    nlopt::opt local_opt(nlopt::LN_COBYLA, n);
    if (settings.method == Conflict::Optimizer::local){
        opt = nlopt::opt(nlopt::GN_MLSL, n);
        local_opt = nlopt::opt(nlopt::LN_COBYLA, n);
        local_opt.set_xtol_rel(settings.xtol);
        local_opt.set_ftol_rel(settings.ftol);
        opt.set_local_optimizer(local_opt);
        opt.set_population(conflict.population);
    }

	// Bounds
    opt.set_lower_bounds({params.gammaBounds.first, params.t0Bounds.first});
    opt.set_upper_bounds({params.gammaBounds.second, params.t0Bounds.second});
	
	// Input argument structure
	OptimizerInput optimIn{this, &instance, &conflict, &airtraffic, params, settings};

	// Set the objective function
	opt.set_min_objective(cost, &optimIn);

	// Set the inequality constraint: g(x) <= 0
    if (settings.method == Conflict::Optimizer::global){
        opt.add_inequality_constraint(nonlcon, &optimIn, settings.nonlconTol);
        opt.add_inequality_constraint(nonlconDomino, &optimIn, settings.nonlconTol);
    }
	
	// Optimizer setting
    std::vector<double> xtol_abs = {0.1, 3.0};   // 1 kt, 3 s
    opt.set_xtol_abs(xtol_abs);
	opt.set_xtol_rel(settings.xtol);
    opt.set_ftol_rel(settings.ftol);
    opt.set_maxeval(settings.maxFunEval);
    if (std::isfinite(settings.tlimit)) opt.set_maxtime(settings.tlimit);

	// Initial guess
    std::vector<double> x = {params.gamma, params.t0};

    // Clamp into bounds
    x[0] = std::clamp(x[0], params.gammaBounds.first, params.gammaBounds.second);
    x[1] = std::clamp(x[1], params.t0Bounds.first, params.t0Bounds.second);

	// Initialize the cost value
    double minf = std::numeric_limits<double>::infinity();

	try {
		nlopt::result r = opt.optimize(x, minf);
		state.exitFlag = (int)r;

		// Store
        state.optimal           = params;
		state.optimal.gamma 	= x[0];
        state.optimal.t0 	    = x[1];
        state.optimal.hold_idx  = params.hold_idx;

		state.fval = minf;
		std::vector<double> grad;
        state.nonlcon = nonlcon(x, grad, &optimIn);
        state.nonlconDomino = nonlconDomino(x, grad, &optimIn);
        state.feasible = (state.nonlcon <= settings.nonlconTol) && ((state.nonlconDomino <= settings.nonlconTol));

	} catch (const std::exception& e) {
		std::cerr << "[opt] exception: " << e.what() << "\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();

	} catch (...) {
		std::cerr << "[opt] unknown exception\n" << std::flush;

		state.optimal = params;
		state.fval = std::numeric_limits<double>::infinity();
		state.feasible = false;
		state.nonlcon = std::numeric_limits<double>::infinity();
        state.nonlconDomino = std::numeric_limits<double>::infinity();
	}
}

std::vector<FSSample> Divert::run(const ConflictRow& instance,
                                    const Conflict& conflict,
                                    const DynamicTraffic& airtraffic,
                                    const std::filesystem::path& dir,
                                    ResolutionLogRow& ResolutionLogRow){
    
    // Find intruders active airborne time indices
    const auto& samples = airtraffic.flights[instance.Iid].samples;
    auto it = std::find_if(samples.begin(), samples.end(), [](const FSSample& s){return s.valid;}); // Forward iterator pointing the first valid sample
    if (it == samples.end()) return {};
    size_t idx0 = std::distance(samples.begin(), it);

    // Define initial parameter set
    Params params;
    params.gamma       = 3.0;
    params.t0          = std::max(onsetTime, idx0);
    params.gammaBounds = {0.0,5.0};
    params.t0Bounds    = {params.t0, (double) instance.t};

    // Plan risk metrics
    params.trafficInteractionRisk = conflict.trafficRisk(instance, samples, airtraffic);
    params.ownshipInteractionRisk = conflict.ownshipRisk(conflict.path(), airtraffic);

    // Optimizer settings
    Conflict::OptimSettings settings;
    settings.method = conflict.opt;
    settings.xtol = 1e-2;
    settings.ftol = 1e-3;
    settings.nonlconTol = 1e-3;
    settings.maxFunEval = 30000;
    settings.tlimit = conflict.tlimit;
    settings.seed = 1995;

    // Find distances to holding points
    // from the first conflict instance position
    size_t nHolds = Holds.size();
    std::vector<double> d(nHolds,0.0);
    const Pos init = airtraffic.flights[instance.Iid].samples[instance.t].pos;
    double bearing;
    for (size_t i = 0; i < nHolds; i++){
        geo_dist(&init, &(Holds[i].pos), &(d[i]), &bearing, &conflict.geoopt);
    }
    std::vector<size_t> idx(nHolds);        // Index array
    std::iota(idx.begin(), idx.end(), 0);   // Fill the array 0, 1, ..., nHolds
    std::sort(idx.begin(), idx.end(), [&d](size_t i, size_t j) { return d[i] < d[j]; }); // argsort

    // Optimize
    for (size_t i = 0; i < nHolds; i++){
        params.hold_idx = idx[i];
        auto start = std::chrono::high_resolution_clock::now();
        optimize(params, instance, conflict, airtraffic, settings);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        state.runtime = elapsed.count() * 1000;
        if (state.feasible) break;
    }
    if (!state.feasible) return {};

    // Out
    const auto newTraj = pathToHold(instance, state.optimal, conflict, airtraffic, dir);
    
    // Log resolution
    ResolutionLogRow.advisory_type  = "Divert";
    ResolutionLogRow.opt_exit_flag  = state.exitFlag;
    ResolutionLogRow.onsetTime      = onsetTime;
    ResolutionLogRow.conflictTime   = instance.t;
    ResolutionLogRow.gamma_GA_star  = state.optimal.gamma;
    ResolutionLogRow.t0_star        = state.optimal.t0;
    ResolutionLogRow.hold_idx_star  = state.optimal.hold_idx;
    ResolutionLogRow.length_star    = path.hdist;
    ResolutionLogRow.feasible       = state.feasible;
    ResolutionLogRow.g1_star        = state.nonlcon;
    ResolutionLogRow.g2_star        = state.nonlconDomino;
    ResolutionLogRow.f_star         = state.fval;
    ResolutionLogRow.runtime_ms     = state.runtime;
    ResolutionLogRow.exit_flag      = 4;
    ResolutionLogRow.plan_traffic_risk = params.trafficInteractionRisk;
    ResolutionLogRow.adv_traffic_risk  = conflict.trafficRisk(instance, newTraj, airtraffic);
    
    return newTraj;
}