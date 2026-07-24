#include "traffic.hpp"

template<typename T>
void read_exact(std::ifstream& file, T& value, const char* errMsg)
{
    file.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!file) {
        throw std::runtime_error(errMsg);
    }
}

/*
    Loads dynamic traffic data
*/
void Traffic::loadDynamicTraffic(const std::filesystem::path& filedir){

    // Open binary data file
    std::ifstream file(filedir, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening dynamic traffic binary file.");
    }

    // Read magic
    char magic[8] = {};
    file.read(magic, 8);
    if (std::memcmp(magic, "ATD1BIN", 7) != 0) {
        throw std::runtime_error("Not an ATD1BIN DynamicTraffic file.");
    }

    // Read data label
    file.read(data.label, 10);
    if (!file) {
        throw std::runtime_error("Dynamic traffic label read failed.");
    }

    // Read number of flights and number of time steps within the data file
    read_exact(file, data.nflights, "Dynamic traffic nflights read failed.");
    read_exact(file, data.nsteps, "Dynamic traffic nsteps read failed.");

    // Read stats
    read_exact(file, data.median, "Dynamic traffic median read failed.");
    read_exact(file, data.p90,    "Dynamic traffic p90 read failed.");
    read_exact(file, data.peak,   "Dynamic traffic peak read failed.");

    // Resize flights + sample vectors
    data.flights.resize(data.nflights);
    for (auto& fs : data.flights) {
        fs.samples.resize(data.nsteps);
    }

    // Read aircraft types
    uint32_t meta_bytes = 0;
    read_exact(file, meta_bytes, "meta_bytes read failed.");
    if (meta_bytes != data.nflights) {
        throw std::runtime_error("Unexpected metadata size.");
    }

    std::vector<uint8_t> type(data.nflights, 0);
    if (meta_bytes > 0) {
        if (meta_bytes < data.nflights) {
            throw std::runtime_error("meta_bytes < nflights (corrupt file).");
        }

        file.read(reinterpret_cast<char*>(type.data()), static_cast<std::streamsize>(type.size() * sizeof(uint8_t)));
        if (!file) throw std::runtime_error("Aircraft type[] read failed.");
    }

    // Reusable buffers (allocated once, reused per flight)
    std::vector<float>   lat(data.nsteps), lon(data.nsteps), alt(data.nsteps), track(data.nsteps), gs(data.nsteps);
    std::vector<int64_t> ts(data.nsteps);
    
    // Read data
    for (uint32_t i = 0; i < data.nflights; i++){
        file.read(reinterpret_cast<char*>(lat.data()), lat.size() * sizeof(float));
        file.read(reinterpret_cast<char*>(lon.data()), lon.size() * sizeof(float));
        file.read(reinterpret_cast<char*>(alt.data()), alt.size() * sizeof(float));
        file.read(reinterpret_cast<char*>(track.data()), track.size() * sizeof(float));
        file.read(reinterpret_cast<char*>(gs.data()), gs.size() * sizeof(float));
        file.read(reinterpret_cast<char*>(ts.data()), ts.size() * sizeof(int64_t));

        if (!file) throw std::runtime_error("DynamicTraffic flight data array read failed.");

        // Reference to sample
        auto& samples = data.flights[i].samples;

        // Write aircraft types
        for (size_t i = 0; i < data.nflights; i++){
            data.flights[i].type = static_cast<AircraftType>(type[i]);
        }

        // Write samples
        for (uint32_t k = 0; k < data.nsteps; ++k) {
            samples[k].pos.lat = lat[k];
            samples[k].pos.lon = lon[k];
            samples[k].pos.alt = alt[k];
            samples[k].pos.hdg = track[k];
            samples[k].gs      = gs[k];
            samples[k].ts      = ts[k];
            samples[k].valid   = !std::isnan(samples[k].pos.lat);
        }
    }
    file.close();
}

/*
    Prints intruder type
*/
AircraftType Traffic::type(const size_t& ID) const {
    return data.flights[ID].type;
}

std::string Traffic::typeString(const size_t& ID) const {

    auto t = type(ID);
    if (t == AircraftType::GA) return "GA";
    else if (t == AircraftType::Turbojet) return "TJ";
    else if (t == AircraftType::Helicopter) return "H";
    else return "XX";
}

/*
    Intruder turn radius estimation
*/
void Traffic::turnRadius(const size_t& ID,
                        const size_t& t,
                        double& R,
                        const DynamicTraffic& airtraffic) const{
    

    // Aircraft type
    auto type = airtraffic.flights[ID].type;
    
    // Get intruder trajectory samples
    auto samples = airtraffic.flights[ID].samples;
    
    // Groundspeed
    auto groundspeed = samples[t].gs * KTS_2_FTS; // [ft/s]

    // Turn radius
    const double g  = 9.81; // Gravitational acceleration [m/s^2]
    const double mu = 30;   // Bank angle [deg]
    if (type == AircraftType::Helicopter){ // If helicopter
        R = groundspeed*groundspeed / ((g * M_2_FT) * std::tan(mu*DEG_2_RAD)); // Minimum turn radius at 30 deg of bank angle
    } else {
        R = groundspeed / (3 * DEG_2_RAD); // Standard turn rate of 3 deg/s
    }
}