# Plan and Avoid: Aircraft Trajectory Planning and Coordination

This repository provides a C/C++ implementation of a **Conflict-Aware Gradient-guided Search for Aircraft Path Planning** in constrained airspace, alongside a reactive **Conflict Resolution Advisory** framework. The system enables path planning and coordination with conflict avoidance and minimizes disruptions to multi-agent activity during priority aerial operations.

Developed at the **Autonomous Aerospace Systems Laboratory (A2Sys)**, Kevin T. Crofton Department of Aerospace and Ocean Engineering, Virginia Tech.

---

## 📂 Repository Structure

The codebase is modularized into two primary subsystems:

- **`plan/`**: The core strategic path planner. Contains the Gradient-Guided 4D Search algorithms, global airspace risk modeling, landing site ranking, and initial trajectory generation.
- **`avoid/`**: The tactical resolution module. Handles localized conflict detection, dynamic air traffic avoidance, and real-time advisory generation (via the `resolution` library).

---

## 🔑 Key Features

- **Gradient-Guided 4D Search**  
  Enhances real-time emergency landing path planning with gradient-based expansion in 4D state space using tree search.
- **Air Traffic Deconfliction & Tactical Avoidance**  
  Incorporates traffic-aware risk to avoid congested regions during emergency descent, with local conflict resolution.
- **3D Airspace Risk Modeling**  
  Considers dense arrival/departure corridors, urban air corridors, and restricted zones.
- **ADS-B Driven Risk Estimation**  
  Uses historical Automatic Dependent Surveillance–Broadcast (ADS-B) data for probabilistic traffic density risk.
- **Proximity-Based 3D Heatmaps**  
  Constructs spatial risk zones around critical urban air corridors and no-fly areas via computational geometry.
- **Cumulative Time-Exposure Metric**  
  Quantifies conflict risk based on path duration through air traffic/multi-agent activity.
- **Look-Ahead Heuristic in Congested Airspace**  
  Improves search efficiency by forecasting future expansion penalties in traffic-heavy zones.
- **Landing Site Selection for Operational Continuity**  
  Ranks candidate sites based on minimal disruption to active traffic corridors.

---

## 📊 Datasets Provided

The repository includes curated datasets for Cessna 182 gliding aircraft performance, as well as urban airspace and population density in the Washington, D.C. area. These datasets support simulation and benchmarking of emergency landing scenarios in congested environments.

### 1️⃣ Gliding Aircraft Performance (`data/aircraft/`)
- **Cessna 182 Gliding Parameter Set and Flight Path Angle Look-Up Table**  
  1) Straight glide descent angles as a function of wind conditions.
  2) Turn glide descent angles with formal forward invariance guarantee as a function of wind conditions.

### 2️⃣ Ground Risk Data (`data/census/`)
- **Washington D.C. Area Census Dataset**  
  Provided as shapefiles (e.g., `dcarea_all_census.shp`), includes polygonal population zones used to compute ground risk exposure along descent paths using spatial R-tree queries.

### 3️⃣ Airport and Runway Data (`lib/aclm/data/`)
- **Washington D.C. Airport and Runway Geometry**  
  Binary file (`landing_sites.bin`) containing runway endpoints, headings, and elevations for major airports in the region. Used to initialize candidate landing sites.

### 4️⃣ Dynamic Air Traffic Data (`lib/airtraffic/data/`)
- **Washington D.C. Airspace Air Traffic Data**  
  Binary files containing spatiotemporal dynamic traffic activity used for planning and conflict resolution. Data file names start with `VFR` or `IFR` depending on the sampled flight conditions, and are distinguished with `LOW`, `MID`, and `HI` to indicate the traffic density (e.g., `VFR_HI.bin` or `IFR_LOW.bin`).

### 5️⃣ Helicopter Route Data (`lib/airtraffic/data/`)
- **Urban Helicopter Route Geometry and Heatmaps**  
  Binary files containing low-altitude helicopter route shapes and associated airspace density heatmaps. Modeled using proximity-based risk for mixed-use corridors.

### 6️⃣ No-Fly Zone Data (`lib/airtraffic/data/`)
- **Restricted and Prohibited Area Geometry and Heatmaps**  
  Binary files representing permanent and temporary no-fly zones in the Washington D.C. airspace. Includes 3D volumetric heatmaps for risk-aware planning.

### 7️⃣ Commercial Airport Corridor Heatmaps (`data/matlab/`)
- **Ronald Reagan Washington National Airport (DCA)**  
  Heatmaps capturing historical departure and arrival flows based on ADS-B data (e.g., `dca_corridors`). Used to evaluate air traffic deconfliction and trajectory feasibility.

> **Note:** All datasets are distributed across the `plan/data/`, `plan/lib/*/data/`, and `avoid/lib/*/data/` directories. Heatmap binary files are preprocessed and ready to be utilized by the planner modules.

---

## 🔧 Build Instructions

> **📖 Detailed Configuration Guide:** Please refer to the **`manual.pdf`** included in this repository for comprehensive details on how to shape the `Makefile` and properly configure the software before running.

### 1. Dependencies

This project uses standard C/C++ libraries and a few geospatial dependencies. On macOS, install them via Homebrew:

`brew install proj spatialindex shapelib`

### 2. Build and Run

The provided test cases simulate joint airspace and ground risk-aware emergency landing planning for an engine-out Cessna 182 over Washington, D.C. Once you have reviewed the `manual.pdf` and configured your Makefiles, build the respective modules:

**For the Planner:**
`cd plan`
`make && ./plan`

**For the Avoidance Module:**
`cd avoid`
`make && ./avoid`

### 3. Post-process
Results and trajectory outputs are written to the respective `out/` folder.

### 4. Custom Case Configuration
To run a different case, input your emergency state coordinates and heading into the configuration file (`plan/aclm.cfg`). 
* Make sure it is a reachable case defined within the modeled airspace environment (e.g., Washington D.C. by default).
* Goal and touchdown states may be left blank as they are dynamically set by the landing site selection module.
* *Note:* The default Cessna 182 dataset supports a maximum of 8 m/s (~ 15.6 knots) of wind speed.

### 5. Software Execution Exit Flags
- `0`: The search solver has found a solution to the best landing site.  
- `1`: The search solver has found a solution to an alternate landing site.  
- `2`: The search solver did not converge within the maximum allowed state expansions. A fallback solution from the Dubins solver has been returned. Increasing `MAX_ITER` in `aclm.cfg` may help find a search-based solution.
- `4`: The search found a solution, but the Dubins benchmark was unavailable or invalid.
- `-1`: An unreachable case has been detected.  
- `-2`: The emergency was initialized near or inside a prohibited area.  
- `-3`: The search open list (priority queue) is empty. A search-based solution was not found. Check the results folder for a Dubins fallback solution if the case is reachable.

---

## 📚 Methodological Foundations

- **Plan-and-Avoid Framework**  
  Tekaslan, H. E. and Atkins, E. M. and Neogi, N., “Plan–and–Avoid: Real-Time Aircraft Trajectory Coordination in a Multi-Agent Environment,” arXiv:XXXXXXXX. Available: http://arxiv.org/abs/XXXXXXXXX.

- **Gradient-guided Search Path Planner**  
  Tekaslan, H. E. and Atkins, E. M., “Gradient-Guided Search for Aircraft Contingency Landing Planning,” *Drones*, Vol. 9, No. 9, 2025. https://doi.org/10.3390/drones9090642.
  
- **Static Airspace Risk Modeling**  
  Tekaslan, H. E. and Atkins, E. M. (2026). Airspace-aware Contingency Landing Planning. arXiv:2602.07074. Available: http://arxiv.org/abs/2602.07074
  
- **Glide Descent Angle Derivation**  
  Tekaslan, H. E. and Atkins, E. M., “Vehicle-to-Vehicle Approach to Assured Aircraft Emergency Road Landings,” *AIAA Journal of Guidance, Control, and Dynamics*, Vol. 48, No. 8, pp. 1800–1817, 2025. https://doi.org/10.2514/1.G008803.
  Tekaslan, H. E. and Atkins, E. M., “Airspeed Forward-Invariance for Unpowered Fixed-Wing Aircraft,” arXiv:2604.22860. Available: http://arxiv.org/abs/2604.22860

---

## 📝 Citation

If you use this software in part or in full for your research or projects, please cite:

Tekaslan, H. E. and Atkins, E. M. (2026). Airspace-aware Contingency Landing Planning. arXiv:2602.07074. Available: http://arxiv.org/abs/2602.07074

`@misc{paa,`
`      title={Plan–and–Avoid: Real-Time Aircraft Trajectory Coordination in a Multi-Agent Environment},` 
`      author={H. Emre Tekaslan and Ella M. Atkins and Natasha Neogi},`
`      year={2026},`
`      eprint={XXXXXX},`
`      archivePrefix={arXiv},`
`      primaryClass={cs.RO},`
`      url={https://arxiv.org/abs/XXXXXX},` 
`}`
