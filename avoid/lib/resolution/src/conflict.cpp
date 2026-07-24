
#include "conflict.hpp"
#include "advisory.hpp"

/*
    = = = = = = = = = = = = = =
    CONFLICT CLASS FUNCTIONS
    = = = = = = = = = = = = = =
*/
std::vector<size_t> Conflict::resolve(const size_t case_id,
                                      const std::filesystem::path& dir,
                                      std::ofstream& caselog,
                                      std::ofstream& resolog) const{
    
    // Load dynamic traffic
    Traffic traffic;
    std::filesystem::path trafficpath =  std::filesystem::current_path() / "lib" / "airtraffic" / "data" / "VFR_HI.bin";
    traffic.loadDynamicTraffic(trafficpath);

    // Air traffic data
    auto& airtraffic = traffic.data;

    // Get the first conflict instances of unique intruders
    std::vector<std::tuple<size_t,size_t,ConflictRow>> conflicts = unique();

    // Exit flags
    // 0: Halt takeoff
    // 1: Speed regulation
    // 2: Altitude regulation
    // 3: Lateral diversion
    // 4: Go around and hold
    // 5: No resolution is found
    // 6: No conflict
    // 7: Initialized in conflict
    std::vector<size_t> exitflag(conflicts.size(), 5);

    // Resolve
    double runtime = 0.0;
    RiskMetrics postResolutionRisk{};
    for (size_t i = 0; i < conflicts.size(); ++i) {

        // Define advisory subclasses
        Halt halt;
        HeliSpeed heliSpeed;
        Speed speed;
        Altitude altitude;
        Extend extend;
        Divert divert;

        // Define initial advisory onset timestamp
        if (tlimit < 10.0){
            heliSpeed.onsetTime             = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink + this->tlimit);
            speed.onsetTime       = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink + this->tlimit);
            altitude.onsetTime    = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink + this->tlimit);
            extend.onsetTime                = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink + this->tlimit);
            divert.onsetTime              = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink + this->tlimit);
        } else {
            heliSpeed.onsetTime             = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink);
            speed.onsetTime       = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink);
            altitude.onsetTime    = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink);
            extend.onsetTime                = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink);
            divert.onsetTime              = static_cast<size_t>(this->pathRuntime + 2*this->tdatalink);
        }

        // Set wind conditions for airspeed check
        speed.windSpeed     = windSpeed;
        speed.windDirection = windDirection;

        // Get a conflict instance
        const auto& [ID, t, instance] = conflicts[i];

        // Console out
        printf("\t%s:%2lu%s %d ---> ", "ID", ID, ", type:", traffic.type(ID));

        // Initialize resolution log structure
        struct ResolutionLogRow ResolutionLogRow{};
        ResolutionLogRow.case_id        = case_id;
        ResolutionLogRow.intruder_id    = ID;
        ResolutionLogRow.intruder_type  = traffic.typeString(ID);

        // Detect initialization in conflict
        const auto& samples = airtraffic.flights[instance.Iid].samples;
        auto it = std::find_if(samples.begin(), samples.end(), [](const FSSample& s){return s.valid;});
        bool initInConflict = false;
        if (it != samples.end()) {
            size_t idx0 = std::distance(samples.begin(), it);
            double alt0 = it->pos.alt;

            if (instance.t < 10) initInConflict = true;
            
            if (!initInConflict){
                size_t nConflict = conflictCount();
                for (size_t j = 0; j < nConflict; ++j) {
                    auto row = conflictInstance(j);
                    if (row.Iid == instance.Iid && (row.t - static_cast<int>(idx0) < 10 && alt0 > 200)) {
                        initInConflict = true;
                    }
                    break;
                }
            }
        }

        // Check if random case is initialized in a conflict
        if (initInConflict){
            for (size_t j = 0; j < conflicts.size(); j++){
                exitflag[j] = 7;
            }
            ResolutionLogRow.exit_flag = exitflag[i];
            ResolutionLogRow.advisory_type = "InitConflict";
            printCaseRow(caselog, case_id, "InitConflict", conflicts.size(), 0, conflicts.size(), runtime, postResolutionRisk);
            std::cout << "Initialized in conflict.\n";
            return exitflag;
        }

        // Halt take-off
        std::vector<FSSample> newTraj{};
        newTraj = halt.run(instance, airtraffic, dir, ResolutionLogRow);
        if (!newTraj.empty()){
            runtime += halt.runtime;
            airtraffic.flights[ID].samples = std::move(newTraj);
            exitflag[i] = 0;
            printResolutionRow(resolog, ResolutionLogRow);
            printf("%s %zu\n", "Exit flag: ", exitflag[i]);
            continue;
        }

        // If the intruder is a fixed-wing or a helicopter
        // Try speed regulation
        if (traffic.type(ID) == AircraftType::Helicopter){
            newTraj = heliSpeed.run(instance, *this, airtraffic, dir, ResolutionLogRow);
            runtime += heliSpeed.state.runtime;
        } else {
            runtime = 0.0;
            newTraj = speed.run(instance, *this, traffic, dir, ResolutionLogRow);
        }
        if (!newTraj.empty()){
            runtime += speed.state.runtime;
            airtraffic.flights[ID].samples = std::move(newTraj);
            exitflag[i] = 1;
        } else {

            // Try altitude regulation
            newTraj = altitude.run(instance, *this, airtraffic, dir, ResolutionLogRow);
            if (!newTraj.empty()){
                runtime += altitude.state.runtime;
                airtraffic.flights[ID].samples = std::move(newTraj);
                exitflag[i] = 2;
            } else {

                // Try lateral diversion
                newTraj = extend.run(instance, *this, airtraffic, dir, ResolutionLogRow);
                if (!newTraj.empty()){
                    runtime += extend.state.runtime;
                    airtraffic.flights[ID].samples = std::move(newTraj);
                    exitflag[i] = 3;
                } else {
                        // Last resort divert to hold
                        newTraj = divert.run(instance, *this, airtraffic, dir, ResolutionLogRow);
                        if (!newTraj.empty()){
                            runtime += divert.state.runtime;
                            airtraffic.flights[ID].samples = std::move(newTraj);
                            exitflag[i] = 4;
                        } else {
                            exitflag[i] = 5;
                        }
                    }
                }
            }
            printf("%s %zu\n", "Exit flag: ", exitflag[i]);
            printResolutionRow(resolog, ResolutionLogRow);
        }

    // Log case results
    size_t n_resolved = 0;
    size_t n_failed   = 0;
    for (auto f : exitflag){
        if (f == 5) n_failed++;
        else if (f != 5) n_resolved++;
    }
    postResolutionRisk = ownshipRisk(path(), airtraffic);
    printCaseRow(caselog, case_id, "Processed", conflicts.size(), n_resolved, n_failed, runtime, postResolutionRisk);

    return exitflag;
}

/*
    Reads emergency landing path generation runtime
*/
void Conflict::getPathRuntime(const std::filesystem::path& filedir){

    // Open path generation runtime log file
    std::ifstream file(filedir);
    if (!file.is_open()){
        std::cerr << "Error opening path runtime log file." << std::endl;
    }

    std::string line;
    std::string lastLine;

    // Skip header
    std::getline(file, line);

    // Read and keep only the last non-empty line
    while (std::getline(file, line)){
        if (!line.empty()){
            lastLine = line;
        }
    }
    if (lastLine.empty()) return;

    // Parse last line
    std::stringstream ss(lastLine);
    std::string cell;

    for (size_t i = 0; i < 4; i++){
        std::getline(ss, cell, ',');
    }

    pathRuntime = std::stod(cell) / 1000.0;
}

/*
    Loads emergency landing path
*/
void Conflict::loadEmergencyPath(const std::filesystem::path& filedir){

    // Clear previously loaded path
    emergencyPath.clear();

    // Open emergency landing solution csv
    std::ifstream file(filedir);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening emergency landing solution file.");
    }

    // Create a string to read lines
    std::string line;

    // Skip header
    std::getline(file, line);

    // Read each line
    double t = 0;
    while (std::getline(file, line)) {
        
        // Convert line to stream
        std::stringstream ss(line);

        // Create a string to store
        std::string cell;

        // Create a temporary path structure
        FSSample s{};

        // Push values into the structure
        std::getline(ss, cell, ',');    // Latitude
        s.pos.lat = std::stod(cell);

        std::getline(ss, cell, ',');    // Longitude
        s.pos.lon = std::stod(cell);

        std::getline(ss, cell, ',');    // Altitude
        s.pos.alt = std::stod(cell);

        std::getline(ss, cell, ',');    // Heading
        s.pos.hdg = std::stod(cell);

        std::getline(ss, cell, ',');    // Timestamp
        s.ts = t;
        s.pos.t = t;

        s.valid = 1;

        emergencyPath.push_back(s);

        // Update timestamp
        t++;
    }
}

/*
    Loads conflict log file
*/
int Conflict::loadConflictLog(const std::filesystem::path& filedir) {

    // Reset previous contents
    log.clear();

    // Open conflict log file
    std::ifstream file(filedir);
    if (!file.is_open()) {
        std::cerr << "Error opening conflict log file." << std::endl;
    }

    // Create a string
    std::string line;

    // Skip header
    std::getline(file, line);

    // For each line
    while (std::getline(file, line)) {

        // Convert line to a stream
        std::stringstream ss(line);

        // Create a string to store
        std::string cell;

        // Create a conflict structure
        ConflictRow r{};

        // Push values into log
        std::getline(ss, cell, ',');
        r.t = std::stoi(cell);

        // Column 1
        std::getline(ss, cell, ',');
        r.Olat = std::stod(cell);

        // Column 2
        std::getline(ss, cell, ',');
        r.Olon = std::stod(cell);

        // Column 3
        std::getline(ss, cell, ',');
        r.Oalt = std::stod(cell);

        // Column 4
        std::getline(ss, cell, ',');
        r.Ohdg = std::stod(cell);

        // Column 5
        std::getline(ss, cell, ',');
        r.Ilat = std::stod(cell);

        // Column 6
        std::getline(ss, cell, ',');
        r.Ilon = std::stod(cell);

        // Column 7
        std::getline(ss, cell, ',');
        r.Ialt = std::stod(cell);

        // Column 8
        std::getline(ss, cell, ',');
        r.Ihdg = std::stod(cell);

        // Column 9
        std::getline(ss, cell, ',');
        r.Iid = std::stoi(cell);

        // Column 10
        std::getline(ss, cell, ',');
        r.JCPA = std::stod(cell);

        log.push_back(r);
    }

    file.close();

    if (log.size() == 0) return 0;
    else return 1;
}

ConflictRow Conflict::conflictInstance(size_t idx) const{
    return log.at(idx);
}

std::vector<std::tuple<size_t,size_t,ConflictRow>> Conflict::unique() const {

    const auto nConflict = conflictCount();

    std::set<size_t> seen;
    std::vector<std::tuple<size_t, size_t, ConflictRow>> unique;

    for (size_t i = 0; i < nConflict; i++) {

        const auto& instance = conflictInstance(i);
        const size_t id = instance.Iid;
        const size_t t  = instance.t;

        // Insert if not seen before
        if (seen.insert(id).second){
            unique.emplace_back(id, t, instance);
        }
    }

    return unique;
}

std::vector<FSSample> Conflict::retimeEmergencyPathBySpeed(double egoSpeedKts,
                                                           double timestep) {
    const size_t Nin = emergencyPath.size();

    if (Nin < 2) return emergencyPath;
    if (egoSpeedKts <= 0.0) return {};
    if (timestep <= 0.0) return {};

    // Build cumulative geometric distance along emergency path
    std::vector<double> sPath(Nin, 0.0);

    for (size_t k = 1; k < Nin; ++k) {

        const auto& p0 = emergencyPath[k-1].pos;
        const auto& p1 = emergencyPath[k].pos;

        double ds = 0.0;
        double course = 0.0;

        geo_dist(&p0, &p1, &ds, &course, &geoopt);

        if (!std::isfinite(ds) || ds < 0.0) {
            ds = 0.0;
        }

        sPath[k] = sPath[k-1] + ds;
    }

    const double sEnd = sPath.back();

    if (sEnd <= 1e-12) return emergencyPath;

    for (size_t k = 1; k < Nin; ++k) {
        assert(sPath[k] >= sPath[k-1]);
    }

    

    // Compute new duration from assumed ego speed
    const double egoSpeedNMps = egoSpeedKts / 3600.0;  // [NM/s]
    const double tEnd = sEnd / egoSpeedNMps;           // [s]

    const size_t Nout = static_cast<size_t>(std::ceil(tEnd / timestep)) + 1;

    std::vector<FSSample> egoTimed;
    egoTimed.reserve(Nout);

    // Resample fixed geometry by distance traveled
    for (size_t i = 0; i < Nout; ++i) {

        double t = static_cast<double>(i) * timestep;
        double s = egoSpeedNMps * t;

        FSSample sample = emergencyPath.front();

        // Make sure the final sample lands exactly at the final path point.
        if (s >= sEnd || i == Nout - 1) {

            sample = emergencyPath.back();
            sample.ts = tEnd;
            sample.gs = egoSpeedKts;

            egoTimed.push_back(sample);
            break;
        }

        auto it = std::lower_bound(sPath.begin() + 1, sPath.end(), s);

        if (it == sPath.end()) {

            sample = emergencyPath.back();
            sample.ts = tEnd;
            sample.gs = egoSpeedKts;

            egoTimed.push_back(sample);
            break;
        }

        const size_t j = static_cast<size_t>(
            std::distance(sPath.begin(), it)
        );

        const double dsSeg = sPath[j] - sPath[j-1];

        double ratio = 0.0;

        if (dsSeg > 1e-12) {
            ratio = (s - sPath[j-1]) / dsSeg;
        }

        ratio = std::clamp(ratio, 0.0, 1.0);

        const Pos& before = emergencyPath[j-1].pos;
        const Pos& after  = emergencyPath[j].pos;

        Pos p;
        p.lat = before.lat + ratio * (after.lat - before.lat);
        p.lon = before.lon + ratio * (after.lon - before.lon);
        p.alt = before.alt + ratio * (after.alt - before.alt);

        sample.pos = p;
        sample.ts  = t;
        sample.gs  = egoSpeedKts;

        egoTimed.push_back(sample);
    }

    return egoTimed;
}

void Conflict::computeSeparation(const Pos& ownship, const Pos& intruder, double& dh, double& dv) const {

    // Get horizontal separation
    double course;
    geo_dist(&ownship, &intruder, &dh, &course, &geoopt);

    if (intruder.alt > IAF_ALT)
        dh -= ADSB_ERR_P + RNP1;
    else
        dh -= ADSB_ERR_P + RNP03;

    // Get vertical separation
    if (intruder.alt > BARO_ALT_LIM)
        dv = fabs(ownship.alt - intruder.alt) - 2*BARO_ERR_A5k;
    else
        dv = fabs(ownship.alt - intruder.alt) - 2*BARO_ERR_B5k;
}

double Conflict::CPACost(const double& deltah, const double& deltav) const{

    if ((deltah >= HORZ_CLEARANCE) || (deltav >= VERT_CLEARANCE)) return 0;

    // Compute the horizontal margin
    double mh = std::max(0.0, HORZ_CLEARANCE - deltah);

    // Compute the vertical margin
    double mv = std::max(0.0, VERT_CLEARANCE - deltav);

    // Horizontal penalty
    double Jh;
    double kh = 1;
    Jh = 2 / (1 + std::exp(kh * mh));

    // Vertical penalty
    double Jv;
    double kv = 0.005;
    Jv = 2 / (1 + std::exp(kv * mv));

    return Jh*Jv;
}

RiskMetrics Conflict::trafficRisk(const ConflictRow& instance,
                                const std::vector<FSSample>& samples,
                                const DynamicTraffic& airtraffic) const{


    // Number of time steps
    const size_t N = samples.size();

    // Pre-allocate out
    RiskMetrics out{};

    // Integrate CPA along the action
    constexpr double dt = 1.0;  // Time step
    double prev = 0.0;  // Previous time step risk
    for (size_t i = 0; i < N; i++){

        // Query position
        const Pos& q = samples[i].pos;
        if (invalidPos(q)) continue;

        // Compute the absolute maximum CPA cost considering nearby traffic
        double Jmax = 0.0;
        for (uint32_t aircraftID = 0; aircraftID < airtraffic.nflights; aircraftID++){

            // If itself, pass
            if (aircraftID == (uint32_t) instance.Iid) continue;
            
            const auto& intrSamples = airtraffic.flights[aircraftID].samples;
            if (i >= intrSamples.size()) continue;

            const Pos& intruder = intrSamples[i].pos;
            if (invalidPos(intruder)) continue;

            // Compute CPA for air traffic with ID: aircraftID
            double deltah, deltav;
            computeSeparation(q, intruder, deltah, deltav);
            double J = CPACost(deltah, deltav);

            // Check if it's higher, if so assign it to the maximum cost
            Jmax = std::max(Jmax, J);
        }

        // Store the cost
        const double curr = Jmax;
        out.worst = std::max(out.worst, curr);

        // Trapezoidal integration
        if (i > 0) {
            out.integral += 0.5 * (prev + curr) * dt;
        }

        // Assign the current risk as the previous
        prev = curr;
    }

    return out;
}

double Conflict::ownshipMarginSingle(const std::vector<FSSample>& ownship,
                                     const std::vector<FSSample>& intruder) const
{
    double worst = -std::numeric_limits<double>::infinity();
    bool any_valid = false;

    const size_t N = std::min(ownship.size(), intruder.size());
    for (size_t i = 0; i < N; ++i) {
        const Pos& q = ownship[i].pos;
        const Pos& p = intruder[i].pos;
        if (invalidPos(q) || invalidPos(p)) continue;

        any_valid = true;

        double dh, dv;
        computeSeparation(q, p, dh, dv);

        const double m = std::min(HORZ_CLEARANCE - dh,
                                  VERT_CLEARANCE - dv);
        worst = std::max(worst, m);
    }

    return any_valid ? worst : -1.0;
}

double Conflict::ownshipMarginRobustTimeWindow(const std::vector<FSSample>& ownship,
                                               const std::vector<FSSample>& intruder,
                                               double VminKts,
                                               double VnomKts,
                                               double VmaxKts,
                                               double timestep) const
{
    double worst = -std::numeric_limits<double>::infinity();
    bool any_valid = false;

    if (ownship.empty() || intruder.empty()) {
        return -1.0;
    }

    if (VminKts <= 0.0 || VnomKts <= 0.0 || VmaxKts <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    if (!(VminKts <= VnomKts && VnomKts <= VmaxKts)) {
        return std::numeric_limits<double>::infinity();
    }

    if (timestep <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    const size_t Nint = intruder.size();
    const size_t Nego = ownship.size();

    for (size_t i = 0; i < Nint; ++i) {

        const Pos& p = intruder[i].pos;

        if (invalidPos(p)) {
            continue;
        }

        const double t = static_cast<double>(i) * timestep;

        // Ego speed uncertainty mapped into nominal ego-time window.
        const double tEarliest = t * VnomKts / VmaxKts;
        const double tLatest   = t * VnomKts / VminKts;

        size_t k0 = static_cast<size_t>(std::floor(tEarliest / timestep));
        size_t k1 = static_cast<size_t>(std::ceil (tLatest   / timestep));

        if (k0 >= Nego) {
            continue;
        }

        k1 = std::min(k1, Nego - 1);

        for (size_t k = k0; k <= k1; ++k) {

            const Pos& q = ownship[k].pos;

            if (invalidPos(q)) {
                continue;
            }

            any_valid = true;

            double dh = 0.0;
            double dv = 0.0;

            computeSeparation(q, p, dh, dv);

            const double m = std::min(HORZ_CLEARANCE - dh,
                                      VERT_CLEARANCE - dv);

            worst = std::max(worst, m);
        }
    }

    return any_valid ? worst : -1.0;
}

double Conflict::ownshipMargin(const std::vector<FSSample>& ownship,
                               const DynamicTraffic& airtraffic) const {

    double worst = -std::numeric_limits<double>::infinity();
    bool any_valid = false;

    const size_t N = ownship.size();
    for (size_t i = 0; i < N; ++i) {

        const Pos& q = ownship[i].pos;
        if (invalidPos(q)) continue;

        for (uint32_t aircraftID = 0; aircraftID < airtraffic.nflights; ++aircraftID) {
            const auto& intrSamples = airtraffic.flights[aircraftID].samples;
            if (i >= intrSamples.size()) continue;

            const Pos& intruder = intrSamples[i].pos;
            if (invalidPos(intruder)) continue;

            any_valid = true;

            double dh, dv;
            computeSeparation(q, intruder, dh, dv);

            const double m = std::min(HORZ_CLEARANCE - dh,
                                      VERT_CLEARANCE - dv);

            worst = std::max(worst, m);
        }
    }

    return any_valid ? worst : -1.0;
}

RiskMetrics Conflict::ownshipRisk(const std::vector<FSSample>& ownship,
                                  const DynamicTraffic& airtraffic) const{

    // Number of time steps
    const size_t N = ownship.size();

    // Pre-allocate out
    RiskMetrics out{};

    // Integrate CPA along the action
    constexpr double dt = 1.0;  // Time step
    double prev = 0.0;  // Previous time step risk
    for (size_t i = 0; i < N; i++){

        // Query position
        const Pos& q = ownship[i].pos;
        if (invalidPos(q)) continue;

        // Compute the absolute maximum CPA cost considering nearby traffic
        double Jmax = 0.0;
        for (uint32_t aircraftID = 0; aircraftID < airtraffic.nflights; aircraftID++){
            
            const auto& intrSamples = airtraffic.flights[aircraftID].samples;
            if (i >= intrSamples.size()) continue;

            const Pos& intruder = intrSamples[i].pos;
            if (invalidPos(intruder)) continue;

            // Compute CPA for air traffic with ID: aircraftID
            double deltah, deltav;
            computeSeparation(q, intruder, deltah, deltav);
            double J = CPACost(deltah, deltav);

            // Check if it's higher, if so assign it to the maximum cost
            Jmax = std::max(Jmax, J);
        }

        // Store the cost
        const double curr = Jmax;
        out.worst = std::max(out.worst, curr);

        // Trapezoidal integration
        if (i > 0) {
            out.integral += 0.5 * (prev + curr) * dt;
        }

        // Assign the current risk as the previous
        prev = curr;
    }

    return out;
}

void printCasesHeader(std::ofstream& log){
    log << "case_id,"
        << "case_status,"
        << "n_conflicts,"
        << "n_resolved,"
        << "n_failed,"
        << "total_runtime_ms,"
        << "adv_ownship_int_cpa,"
        << "adv_ownship_worst_cpa\n";
}

void printCaseRow(std::ofstream& log,
                  size_t case_id,
                  const std::string& case_status,
                  size_t n_conflicts,
                  size_t n_resolved,
                  size_t n_failed,
                  double runtime_ms,
                  RiskMetrics risk){
    log << case_id << ","
        << case_status << ","
        << n_conflicts << ","
        << n_resolved << ","
        << n_failed << ","
        << runtime_ms << ","
        << risk.integral << ","
        << risk.worst << std::endl;
}

void printResolutionsHeader(std::ofstream& log) {

    // Case identifiers
    log << "case_id,"
        << "intruder_id,"
        << "intruder_type,"
        << "advisory_type,"
        << "advisory_exit_flag,"
        << "optimization_exit_flag,"
        << "onsetTime_s,"
        << "conflictTime_s,"
        << "plan_traffic_int_cpa,"
        << "adv_traffic_int_cpa,"
        << "plan_traffic_worst_cpa,"
        << "adv_traffic_worst_cpa,"

        // Speed regulation optimization parameters
        << "V_kts,"
        << "t0_s,"
        << "T_s,"

        // Altitude regulation optimization parameters
        << "H_ft,"
        << "t0_s,"
        << "T_s,"
        << "gamma0_deg,"
        << "gamma1_deg,"

        // Lateral diversion optimization parameters
        << "theta_deg,"
        << "t0_s,"
        << "tf_s,"
        << "length_nm,"

        // Hold optimization parameters
        << "gamma_deg,"
        << "length_nm_hold,"
        << "t0_s,"
        << "hold_idx,"

        // Common advisory optimization variables
        << "g1_star,"
        << "g2_star,"
        << "f_star,"
        << "feasible,"
        << "runtime_ms" << std::endl;
}
void printResolutionRow(std::ofstream& log, const ResolutionLogRow& row){

    // Case identifiers
    log << row.case_id << ","
        << row.intruder_id << ","
        << row.intruder_type << ","
        << row.advisory_type << ","
        << row.exit_flag << ","
        << row.opt_exit_flag << ","
        << row.onsetTime << ","
        << row.conflictTime << ","

        // Traffic risk metrics
        << row.plan_traffic_risk.integral << ","
        << row.adv_traffic_risk.integral << ","
        << row.plan_traffic_risk.worst << ","
        << row.adv_traffic_risk.worst << ",";

    // Speed regulation optimization parameters
    if (row.advisory_type == "HeliSpeed" || row.advisory_type == "Speed"){
        printOptional(log, row.V_star);       log << ",";
        printOptional(log, row.t0_star);       log << ",";
        printOptional(log, row.dT_star);       log << ",";
    } else {
        log << "," << "," << ",";
    }

    // Altitude regulation optimization parameters
    if (row.advisory_type == "Altitude"){
        printOptional(log, row.dH_star);       log << ",";
        printOptional(log, row.t0_star);       log << ",";
        printOptional(log, row.T_star);        log << ",";
        printOptional(log, row.gamma0_star);   log << ",";
        printOptional(log, row.gamma1_star);   log << ",";
    } else {
        log << "," << "," << "," << "," << ",";
    }

    // Lateral diversion optimization parameters
    if (row.advisory_type == "Diversion"){
        printOptional(log, row.theta_star);   log << ",";
        printOptional(log, row.t0_star);      log << ",";
        printOptional(log, row.tf_star);      log << ",";
        printOptional(log, row.length_star);  log << ",";
    } else {
        log << "," << "," << "," << ",";
    }

    // Hold optimization parameters
    if (row.advisory_type == "GoAround"){
        printOptional(log, row.gamma_GA_star);  log << ",";
        printOptional(log, row.length_star);  log << ",";
        printOptional(log, row.t0_star);log << ",";
        printOptional(log, row.hold_idx_star);  log << ",";
    } else {
        log << "," << "," << "," << ",";
    }

    // Common advisory optimization variables
    printOptional(log, row.g1_star);       log << ",";
    printOptional(log, row.g2_star);       log << ",";
    printOptional(log, row.f_star);        log << ",";
    printOptionalBool(log, row.feasible);  log << ",";
    printOptional(log, row.runtime_ms);

    log << "\n";
}


void printOptional(std::ofstream& log, const std::optional<double>& x){
    if (x) log << *x;
}
void printOptionalBool(std::ofstream& log, const std::optional<bool>& x){
    if (x) log << static_cast<int>(*x);
}

/*
    Reads emergency landing planner log file
*/
std::vector<size_t> loadLogCSV(const std::filesystem::path& p){

    // Read log
    std::ifstream logfile(p);
    if (!logfile) throw std::runtime_error("Log file cannot be opened.");

    // Create a string
    std::string line;

    // Skip header
    std::getline(logfile, line);

    // Create exit flag array
    std::vector<size_t> exitflag;

    // Read solution flags
    while (std::getline(logfile, line)){

        // Convert line to a stream
        std::stringstream ss(line);

        // Create a string to store
        std::string cell;

        // Skip the first column
        std::getline(ss, cell, ',');

        // Push values into log
        std::getline(ss, cell, ',');
        exitflag.push_back(std::stoull(cell));

    }
    return exitflag;
}