/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   #
#    Gradient-guided Search for Assured Contingency Landing Management              #
#    Copyright (C) 2026  Huseyin Emre Tekaslan                                      #
#                                                                                   #
#    This program is free software: you can redistribute it and/or modify           #
#    it under the terms of the GNU General Public License as published by           #
#    the Free Software Foundation, either version 3 of the License, or              #
#    (at your option) any later version.                                            #
#                                                                                   #
#    This program is distributed in the hope that it will be useful,                #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  #
#    GNU General Public License for more details.                                   #
#                                                                                   #
#    You should have received a copy of the GNU General Public License              #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.         #
#                                                                                   #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   %
#    Assured Contingency Landing Management Main Functions                          %
#    Multi-Agent Conflict-Aware                                                     %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : H. Emre Tekaslan (tekaslan@vt.edu)                                   %
#    Date    : July 2026                                                         %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "aclm.h"

int runEmergencyPlanning(SearchProblem *problem, char * casefolder, char * configfile) {

    printf("*******************************************************************************************\n");
    printf("* Multi-Agent Conflict-Aware                                                              *\n");
    printf("* Aircraft Contingency Landing Planner                                                    *\n");
    printf("* Using Gradient-guided 4D Discrete Search and 3D Dubins Solver                           *\n");
    printf("* Copyright (C) 2026 Huseyin Emre Tekaslan, (tekaslan@vt.edu) @ A2Sys Lab, Virginia Tech  *\n");
    printf("* This program comes with ABSOLUTELY NO WARRANTY.                                         *\n");
    printf("* This is free software, and you are welcome to redistribute it under certain conditions. *\n");
    printf("*******************************************************************************************\n");

    // Initialize exit flag
    problem->exitFlag = 9;

    // CFG File directory
    char cfgdir[256];  // Adjust size as needed
    sprintf(cfgdir, "aclm.cfg");

    // Get the solver type
    readSolverType(problem, cfgdir);

    // Load problem parameters
    setUpProblem(problem, cfgdir);

    // Initial emergency state is the state where the search begins
    problem->InitialSearch = problem->Initial;

    // Load airspace heatmaps
    // This is a GGS-ACLMv1 feature.
    loadHeatmaps(problem);
    loadDynamicTrafficParams(problem, cfgdir);
    loadDynamicTraffic(problem);
    
    // Select a landing site, if none of the landing sites are infeasible, return failure.
    int feasible = selectLandingSite(problem);
    if (feasible){
        problem->exitFlag = -1;  // None of the landing sites is reachable.
        return EXIT_FAILURE;
    }

    // Log time
    FILE * file;
    char dir[300];
    sprintf(dir,"%s/runtime.csv", casefolder);
    file = fopen(dir, "w");
    if (!file){
        perror("Runtime log csv cannot be opened.\n");
        return EXIT_FAILURE;
    }
    fprintf(file, "Rank Index, Site ID, Exit flag, Runtime [ms]\n");
    
    // Iterate path search over ranked landing sites
    problem->siteIndex = -1;
    while (((problem->exitFlag == 2 || problem->exitFlag == -3) || problem->siteIndex == -1) && (problem->siteIndex < problem->numRanked - 1)) {

        // Update landing site ordered index
        problem->siteIndex++;

        // Initialize landing site-related distance parameters
        initDistanceParams(problem, problem->siteIndex);

        // Write the touchdown and goal states to the config file
        editCFG(configfile, casefolder, &problem->InitialSearch, &problem->Goal, &problem->Touchdown);

        // Recompute the goal state normal vector
        if (problem->siteIndex == 0) {
            problem->goalNormal = normalVector(&problem->Goal, &problem->GeoOpt);
        } else {
            free(problem->goalNormal); problem->goalNormal = NULL;
            problem->goalNormal = normalVector(&problem->Goal, &problem->GeoOpt);

            destroyPriorityQueue(problem->pq); problem->pq = NULL;
            destroyClosedList(problem->explored); problem->explored = NULL;

            if (problem->altitudeModification) {
                free(problem->altitudeModification);
                problem->altitudeModification = NULL;
            }
        }

        // Initialize the search
        initializeSearch(problem);
        if (problem->exitFlag == -2) return EXIT_FAILURE;    // Initialized in a prohibited airspace

        // Run the gradient-guided search emergency landing path planning
        runSearch(problem);

        // Update timer
        problem->searchRuntime[problem->siteIndex] = problem->totalSearchRuntime;

        // Log runtime
        fprintf(file, "%d, %d, %d, %.6f\n", problem->siteIndex, problem->rankedSites[problem->siteIndex]->ident, problem->exitFlag, problem->searchRuntime[problem->siteIndex]);
    }
    fclose(file);

    // Update path planning trial count
    if (problem->exitFlag == 0 && problem->siteIndex > 0) problem->exitFlag =  1;

    // Update total time
    problem->totalSearchRuntime = 0;
    for (int i = 0; i < problem->siteIndex; i++){
        problem->totalSearchRuntime += problem->searchRuntime[i];
    }

    return EXIT_SUCCESS;
    
}

/*
    Sets up the emergency path planning problem
*/
void setUpProblem(SearchProblem *problem, const char* cfgdir) {

    // Avoid dangling pointers
    problem->ac = NULL;
    problem->census = NULL;
    problem->merged = NULL;
    problem->traffic = NULL;
    problem->heli = NULL;
    problem->prohibited = NULL;
    problem->dynTraffic = NULL;
    problem->pq = NULL;
    problem->explored = NULL;

    problem->HoldingPoint = NULL;
    problem->holdInbound = NULL;
    problem->holdOutbound = NULL;

    problem->path2Hold = NULL;
    problem->finalTurn = NULL;

    problem->DubinsOpt = NULL;
    problem->path = NULL;
    problem->pathActions = NULL;

    problem->altitudeModification = NULL;
    problem->nearest = NULL;
    problem->bestSite = NULL;

    // Load initial, goal, and touchdown states
    loadInitialState(problem, cfgdir);
    loadGoalState(problem, cfgdir);
    loadTouchdownState(problem, cfgdir);

    // Load aircraft parameters
    problem->ac = (Aircraft *) malloc(sizeof(Aircraft));
    loadAircraftParams(problem, cfgdir);

    // Load Earth model
    loadEarthModelParams(problem, cfgdir);

    // Load population risk computation parameters 
    char *shapefile = loadRiskComputationParams(problem, cfgdir);

    // Load search parameters
    loadSearchParams(problem, cfgdir);

    // Load flight path angle datasets
    char ac_name[10] = "c182";
    readGammaStraight(problem->ac, ac_name);
    readGammaTurn(problem->ac);

    // Load census data
    problem->census = initCensus(shapefile);   // Population density query dataset
    free(shapefile); shapefile = NULL;
    
    // Get min. segment length
    problem->lmin = problem->ac->turnRadius * (tan(problem->deltaCourse[4]/2 * DEG_2_RAD) + tan(problem->deltaCourse[4]/2 * DEG_2_RAD));

    // Discretization parameters
    problem->dAltitude =  floor((360.0 / problem->deltaCourse[4]) * problem->lmin * tan(problem->ac->gammaOptTurn*DEG_2_RAD) / 10.0) * 10.0;;   // [ft]

    // Get the outer solution identification radius
    problem->identRadius = 3*problem->ac->turnRadius;

    // Goal state normal vector
    problem->goalNormal = normalVector(&problem->Goal, &problem->GeoOpt);

    // Allocate memory for the final turn Dubins path
    problem->finalTurn = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));
    Traj_InitArray(problem->finalTurn->traj,4);

    problem->DubinsOpt = (struct DubinsOpt *) malloc(sizeof(struct DubinsOpt));
    problem->DubinsOpt->verbose = 0;
    problem->DubinsOpt->trajopt.verbose = 0;
    problem->DubinsOpt->trajopt.geoopt.verbose = 0;
    problem->DubinsOpt->trajopt.geoopt.model = WGS84;

    // Timers
    problem->totalSearchRuntime = 0;    // Reset timer
    problem->totalRiskRuntime = 0;      // Reset timer

}

/*
    Loads airspace heatmaps into structures
    Grid sizes are hard-coded in the function.
*/
void loadHeatmaps(SearchProblem *problem) {

    // Air traffic heatmap
    problem->traffic = (airspaceHeatmap * ) calloc(1,sizeof(airspaceHeatmap));
    problem->traffic->LAT_SIZE = NLAT;
    problem->traffic->LON_SIZE = NLON;
    problem->traffic->ALT_SIZE = NALT;
    problem->traffic->grid_on_cells = 1;
    problem->traffic->gridfile = "lib/airtraffic/data/traffic_density_grid.bin";
    problem->traffic->datafile = "lib/airtraffic/data/traffic_density.bin";
    readHeatMap(problem->traffic);

    // Helicopter routes heatmap
    problem->heli = (airspaceHeatmap * ) calloc(1,sizeof(airspaceHeatmap));
    if (!problem->heli) {
        perror("problem->heli Invalid memory allocation.\n");
        return;
    }
    problem->heli->LAT_SIZE = 250;
    problem->heli->LON_SIZE = 250;
    problem->heli->ALT_SIZE = 43;
    problem->heli->grid_on_cells = 0;
    problem->heli->gridfile = "lib/airtraffic/data/heligrid.bin";
    problem->heli->datafile = "lib/airtraffic/data/helicost.bin";
    readHeatMap(problem->heli);

    // Prohibited zone heatmap
    problem->prohibited = (airspaceHeatmap * ) calloc(1,sizeof(airspaceHeatmap));
    problem->prohibited->LAT_SIZE = 250;
    problem->prohibited->LON_SIZE = 250;
    problem->prohibited->ALT_SIZE = 2;
    problem->prohibited->grid_on_cells = 0;
    problem->prohibited->gridfile = "lib/airtraffic/data/nofly_grid.bin";
    problem->prohibited->datafile = "lib/airtraffic/data/nofly.bin";
    readHeatMap(problem->prohibited);

}

void loadInitialState(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];

    if (!cfg) {
        perror("Cannot open the cfg file..\n");
        exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "INIT_LAT", 8) == 0) {
            sscanf(line, "INIT_LAT = %lf", &problem->Initial.lat);
        } else if (strncmp(line, "INIT_LON", 8) == 0) {
            sscanf(line, "INIT_LON = %lf", &problem->Initial.lon);
        } else if (strncmp(line, "INIT_ALT", 8) == 0) {
            sscanf(line, "INIT_ALT = %lf", &problem->Initial.alt);
        } else if (strncmp(line, "INIT_HDG", 8) == 0) {
            sscanf(line, "INIT_HDG = %lf", &problem->Initial.hdg);
        }
    }
}

void loadGoalState(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];

    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "GOAL_LAT", 8) == 0) {
            sscanf(line, "GOAL_LAT = %lf", &problem->Goal.lat);
        } else if (strncmp(line, "GOAL_LON", 8) == 0) {
            sscanf(line, "GOAL_LON = %lf", &problem->Goal.lon);
        } else if (strncmp(line, "GOAL_ALT", 8) == 0) {
            sscanf(line, "GOAL_ALT = %lf", &problem->Goal.alt);
        } else if (strncmp(line, "GOAL_HDG", 8) == 0) {
            sscanf(line, "GOAL_HDG = %lf", &problem->Goal.hdg);
        }
    }
    fclose(cfg);
}

void loadTouchdownState(SearchProblem *problem, const char* cfgdir) {
    FILE *cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}

    char line[300];

    while (fgets(line, sizeof(line), cfg)) {
        // Remove newline character if present
        line[strcspn(line, "\n")] = '\0';

        if (strncmp(line, "TD_LAT", 6) == 0) {
            if (sscanf(line, "TD_LAT = %lf", &problem->Touchdown.lat) != 1) {
                fprintf(stderr, "Error reading TD_LAT\n");
            }
        } else if (strncmp(line, "TD_LON", 6) == 0) {
            if (sscanf(line, "TD_LON = %lf", &problem->Touchdown.lon) != 1) {
                fprintf(stderr, "Error reading TD_LON\n");
            }
        } else if (strncmp(line, "TD_ALT", 6) == 0) {
            if (sscanf(line, "TD_ALT = %lf", &problem->Touchdown.alt) != 1) {
                fprintf(stderr, "Error reading TD_ALT\n");
            }
        } else if (strncmp(line, "TD_HDG", 6) == 0) {
            if (sscanf(line, "TD_HDG = %lf", &problem->Touchdown.hdg) != 1) {
                fprintf(stderr, "Error reading TD_HDG\n");
            }
        }
    }

    fclose(cfg);
}
void loadAircraftParams(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];

    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "TURN_RAD", 8) == 0) {
            sscanf(line, "TURN_RAD = %lf", &problem->ac->turnRadius);
        } else if (strncmp(line, "GAMMA_BG_TURN", 13) == 0) {
            sscanf(line, "GAMMA_BG_TURN = %lf", &problem->ac->gammaBGturn);
        } else if (strncmp(line, "GAMMA_OPT_TURN", 14) == 0) {
            sscanf(line, "GAMMA_OPT_TURN = %lf", &problem->ac->gammaOptTurn);
        } else if (strncmp(line, "AIRSPEED", 8) == 0) {
            sscanf(line, "AIRSPEED = %lf", &problem->ac->airspeed);
        }
    }

    // Convert airspeed to ft/s
    problem->ac->airspeed *= KT_2_MS*M_2_FT;

    fclose(cfg);
}

void loadEarthModelParams(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];

    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "EARTH_MODEL", 11) == 0) {
            sscanf(line, "EARTH_MODEL = %u", &problem->GeoOpt.model);
        } else if (strncmp(line, "WIND_SPEED", 10) == 0) {
            sscanf(line, "WIND_SPEED = %lf", &problem->windSpeed);
        } else if (strncmp(line, "WIND_DIR", 8) == 0) {
            sscanf(line, "WIND_DIR = %lf", &problem->windDirection);
        }
    }
    problem->GeoOpt.verbose = 0;

    fclose(cfg);
}

char *loadRiskComputationParams(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];
    char *shapefile = (char *) malloc(300*sizeof(char));

    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "ASRISK_W_HMAX", 13) == 0) {
            sscanf(line, "ASRISK_W_HMAX = %lf", &problem->asrisk_w_hmax);
        } else if (strncmp(line, "ASRISK_W_HMIN", 13) == 0) {
            sscanf(line, "ASRISK_W_HMIN = %lf", &problem->asrisk_w_hmin);
        } else if (strncmp(line, "AIRSPACE_INTEG_DT", 17) == 0) {
            sscanf(line, "AIRSPACE_INTEG_DT = %lf", &problem->airspaceRiskTimeStep);
        } else if (strncmp(line, "GROUND_INTEG_DT", 15) == 0) {
            sscanf(line, "GROUND_INTEG_DT = %lf", &problem->groundRiskTimeStep);
        } else if (strncmp(line, "GROUND_RISK_WEIGHT", 18) == 0) {
            sscanf(line, "GROUND_RISK_WEIGHT = %lf", &problem->w_gp);
        } else if (strncmp(line, "CROSSOVER_H", 11) == 0) {
            sscanf(line, "CROSSOVER_H = %lf", &problem->crossoverAlt);
        } else if (strncmp(line, "MIN_POP", 7) == 0) {
            sscanf(line, "MIN_POP = %lf", &problem->minPopulation);
        } else if (strncmp(line, "MAX_POP", 7) == 0) {
            sscanf(line, "MAX_POP = %lf", &problem->maxPopulation);
        } else if (strncmp(line, "CENSUS_SHPFILE_DIR", 18) == 0) {
            sscanf(line, "CENSUS_SHPFILE_DIR = %s", shapefile);
        }
    }

    // printf("shapefile %s\n", *shapefile);

    // Airspace risk weight
    problem->w_ga = 1 - problem->w_gp;

    // Spatial step for integration sampling
    problem->airspace_dx = problem->airspaceRiskTimeStep*problem->ac->airspeed;
    problem->ground_dx = problem->groundRiskTimeStep*problem->ac->airspeed;

    fclose(cfg);

    return shapefile;
}

void loadSearchParams(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];

    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "MAX_ITER", 8) == 0) {
            sscanf(line, "MAX_ITER = %zu", &problem->maxIter);
        } else if (strncmp(line, "DISCRETIZE_HDG", 14) == 0) {
            sscanf(line, "DISCRETIZE_HDG = %lf", &problem->dPsi);
        } else if (strncmp(line, "DISCRETIZE_TIME", 15) == 0) {
            sscanf(line, "DISCRETIZE_TIME = %d", &problem->dTime);
        } else if (strncmp(line, "ALLOWED_HDG_CHG", 15) == 0) {
            sscanf(line, "ALLOWED_HDG_CHG = %lf, %lf, %lf, %lf, %lf", &problem->deltaCourse[0],
                                                                    &problem->deltaCourse[1],
                                                                    &problem->deltaCourse[2],
                                                                    &problem->deltaCourse[3],
                                                                    &problem->deltaCourse[4]);
        } else if (strncmp(line, "ADAPTIVE_SEGM_LEN_H", 19) == 0) {
            sscanf(line, "ADAPTIVE_SEGM_LEN_H = %lf", &problem->adaptiveLenghtAltitude);
        }
    }
    fclose(cfg);
}

void loadDynamicTrafficParams(SearchProblem *problem, const char* cfgdir) {

    // Open the config file
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}
    char line[300];

    // If dynamic air traffic is enabled:
    if (problem->solver == 3){

        // Allocate memory for dynamic traffic structure
        problem->dynTraffic = calloc(1, sizeof(DynamicTraffic));
        problem->dynTraffic->datafile =  (char *) malloc(300*sizeof(char));
        if (!problem->dynTraffic || !problem->dynTraffic->datafile) {
            fprintf(stderr, "Dynamic traffic memory allocation failed.\n");
            exit(EXIT_FAILURE);
        }

        while (fgets(line, sizeof(line), cfg)) {
            if (strncmp(line, "DYNAMIC_TRAFFIC_DATAFILE", 24) == 0) {
                sscanf(line, "DYNAMIC_TRAFFIC_DATAFILE = %s", problem->dynTraffic->datafile);
            } else if (strncmp(line, "RNP_TER", 7) == 0) {
                sscanf(line, "RNP_TER = %lf", &problem->RNP_TER);
            } else if (strncmp(line, "RNP_FIN", 7) == 0) {
                sscanf(line, "RNP_FIN = %lf", &problem->RNP_FIN);
            } else if (strncmp(line, "ERR_ADSB_P", 10) == 0) {
                sscanf(line, "ERR_ADSB_P = %lf", &problem->ERR_ADSB_P);
            } else if (strncmp(line, "ERR_ADSB_V", 10) == 0) {
                sscanf(line, "ERR_ADSB_V = %lf", &problem->ERR_ADSB_V);    
            } else if (strncmp(line, "DELAY_ADSB", 10) == 0) {
                sscanf(line, "DELAY_ADSB = %lf", &problem->DELAY_ADSB);
            } else if (strncmp(line, "ERR_BARO_A5000", 14) == 0) {
                sscanf(line, "ERR_BARO_A5000 = %lf", &problem->ERR_BARO_A5000);
            } else if (strncmp(line, "ERR_BARO_B5000", 14) == 0) {
                sscanf(line, "ERR_BARO_B5000 = %lf", &problem->ERR_BARO_B5000);
            } else if (strncmp(line, "HORZ_CLEARANCE", 14) == 0) {
                sscanf(line, "HORZ_CLEARANCE = %lf", &problem->HORZ_CLEARANCE);
            } else if (strncmp(line, "VERT_CLEARANCE", 14) == 0) {
                sscanf(line, "VERT_CLEARANCE = %lf", &problem->VERT_CLEARANCE);
            } else if (strncmp(line, "IAF_ALT", 7) == 0) {
                sscanf(line, "IAF_ALT = %lf", &problem->IAF_Alt);
            }
        }
    }

    fclose(cfg);
}

/* 
    Returns flight path angle array for the straight segment
    Output = [gamma_bg, gamma_vfe, gamma_opt]
*/
double *getOptimalGamma(double *course, SearchProblem *problem) {

    /* Angular difference of wind and straight path segment */
    double angDiff;
    angDiff = wrapTo360(problem->windDirection - *course);

    // Temporary wind speed and direction arrays
    double windDir[9]   = {0, 45, 90, 135, 180, 225, 270, 315, 360};    // [deg]
    double windSpeed[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};                  // [m/s]

    // Get indices
    int windDirIdx      = 0;
    int windSpeedIdx    = 0;

    /* Get the speed index */
    for (; (windSpeedIdx < 9) && (windSpeed[windSpeedIdx] < problem->windSpeed*KT_2_MS) && (windSpeed[windSpeedIdx+1] < problem->windSpeed*KT_2_MS) ; windSpeedIdx++);
    /* Get the direction index */
    for (; (windDirIdx < 9) && (windDir[windDirIdx] < angDiff) && (windDir[windDirIdx+1] < angDiff); windDirIdx++);

    // Interpolate w.r.t to direction
    double tmpLookUp[9];
    for (int i = 0; i < 9; i++) {
        tmpLookUp[i] = (angDiff - windDir[windDirIdx])*(problem->ac->gamma_bg[i][windDirIdx+1] - problem->ac->gamma_bg[i][windDirIdx])/(windDir[windDirIdx+1] - windDir[windDirIdx]) + problem->ac->gamma_bg[i][windDirIdx];
    }

    // Interpolate w.r.t to speed
    double gamma_bg = (tmpLookUp[windSpeedIdx+1] - tmpLookUp[windSpeedIdx])*
                        (problem->windSpeed*KT_2_MS - windSpeed[windSpeedIdx]) + tmpLookUp[windSpeedIdx];

    
    // Interpolate w.r.t to direction
    for (int i = 0; i < 9; i++) {
        tmpLookUp[i] = (angDiff - windDir[windDirIdx])*(problem->ac->gamma_vfe[i][windDirIdx+1] - problem->ac->gamma_vfe[i][windDirIdx])/(windDir[windDirIdx+1] - windDir[windDirIdx]) + problem->ac->gamma_vfe[i][windDirIdx];
    }

    // Interpolate w.r.t to speed
    double gamma_vfe = (tmpLookUp[windSpeedIdx+1] - tmpLookUp[windSpeedIdx])*
                        (problem->windSpeed*KT_2_MS - windSpeed[windSpeedIdx]) + tmpLookUp[windSpeedIdx];

    // Interpolate w.r.t to direction
    for (int i = 0; i < 9; i++) {
        tmpLookUp[i] = (angDiff - windDir[windDirIdx])*(problem->ac->gamma_opt[i][windDirIdx+1] - problem->ac->gamma_opt[i][windDirIdx])/(windDir[windDirIdx+1] - windDir[windDirIdx]) + problem->ac->gamma_opt[i][windDirIdx];
    }

    // Interpolate w.r.t to speed
    double gamma_opt = (tmpLookUp[windSpeedIdx+1] - tmpLookUp[windSpeedIdx])*
                        (problem->windSpeed*KT_2_MS - windSpeed[windSpeedIdx]) + tmpLookUp[windSpeedIdx];

    // Flight path angles
    double *gamma = (double *) malloc(3*sizeof(double));
    if (gamma == NULL) {
        fprintf(stderr, "Memory allocation for gamma failed.\n");
        return NULL;
    }

    gamma[0] = -gamma_bg;
    gamma[1] = -gamma_vfe;
    gamma[2] = -gamma_opt;

    return gamma;
}

/*
    Returns the optimal flight path angle for turns
    to maintain the airspeed
*/
double getOptimalGammaTurn(Node *node, SearchProblem *problem) {
    // Dataset indices
    size_t i = 0, j = 0, k = 0;

    // Relative wind angle β = windDir - heading, wrapped to [-180,180)
    double angDiff = wrapTo180(problem->windDirection - node->state.hdg);

    // Wind direction bins: -180:15:165 (24 bins)
    double windDir[NDIR];   // NDIR = 24
    windDir[0] = -180.0;
    for (int n = 1; n < NDIR; n++) {
        windDir[n] = windDir[n-1] + 15.0; // with 15 degrees intervals
    }

    // Find lower index j0
    for (j = 0; j < NDIR-1; j++) {
        if (angDiff >= windDir[j] && angDiff < windDir[j+1]) break;
    }
    int j0 = (int)j;
    int j1 = (j0+1 < NDIR) ? j0+1 : 0; // wrap across seam

    // Fraction ty between windDir[j0] and windDir[j1]
    double beta0 = windDir[j0];
    double beta1 = (j1 == 0 ? 180.0 : windDir[j1]);
    double ty = (angDiff - beta0) / (beta1 - beta0);
    if (ty < 0) ty = 0; if (ty > 1) ty = 1;

    // Wind speed bins: 0..8 m/s (9 bins)
    double windSpeed[NW];   // NW = 9
    for (int n = 0; n < NW; n++) {
        windSpeed[n] = (double) n;
    }

    double Wq = problem->windSpeed * KT_2_MS;
    for (i = 0; i < NW - 1; i++) {
        if (Wq >= windSpeed[i] && Wq < windSpeed[i+1]) break;
    }
    int i0 = (int) i;
    int i1 = (i0+1 < NW) ? i0 + 1 : NW - 1;

    double tx = (Wq - windSpeed[i0]) / (windSpeed[i1] - windSpeed[i0]);
    if (tx < 0) tx = 0; if (tx > 1) tx = 1;

    // Heading change index for deltaCourse \in {-90,-75,...,90}
    double headingChange = node->action->deltaCourse;

    // Round to nearest multiple of 15
    int step = (int) lround(headingChange / 15.0);

    // Range check
    if (step < -6 || step > 6) {
        return NAN;  // invalid
    }

    // Map [-6..6] -> [0..12]
    if (step > 0) k = step + 5;
    else k = step + 6;

    // Bilinear interpolation from 4 surrounding bins
    double g00 = problem->ac->gamma_turn[i0][j0][k];
    double g01 = problem->ac->gamma_turn[i0][j1][k];
    double g10 = problem->ac->gamma_turn[i1][j0][k];
    double g11 = problem->ac->gamma_turn[i1][j1][k];
    

    double g0 = g00*(1.0 - ty) + g01*ty;
    double g1 = g10*(1.0 - ty) + g11*ty;
    double gamma_turn = g0*(1.0 - tx) + g1*tx;

    return fabs(gamma_turn);
}

/*
    Returns the optimal flight path angle for turns
    to maintain the airspeed
*/
double getViableGamma(Node *node, SearchProblem *problem) {
    // Dataset indices
    size_t i = 0, j = 0, k = 0;

    // Relative wind angle β = windDir - heading, wrapped to [-180,180)
    double angDiff = wrapTo180(problem->windDirection - node->state.hdg);

    // Wind direction bins: -180:15:165 (24 bins)
    double windDir[NDIR];   // NDIR = 24
    windDir[0] = -180.0;
    for (int n = 1; n < NDIR; n++) {
        windDir[n] = windDir[n-1] + 15.0; // with 15 degrees intervals
    }

    // Find lower index j0
    for (j = 0; j < NDIR-1; j++) {
        if (angDiff >= windDir[j] && angDiff < windDir[j+1]) break;
    }
    int j0 = (int)j;
    int j1 = (j0+1 < NDIR) ? j0+1 : 0; // wrap across seam

    // fraction ty between windDir[j0] and windDir[j1]
    double beta0 = windDir[j0];
    double beta1 = (j1 == 0 ? 180.0 : windDir[j1]);
    double ty = (angDiff - beta0) / (beta1 - beta0);
    if (ty < 0) ty = 0; if (ty > 1) ty = 1;

    // Wind speed bins: 0..8 m/s (9 bins)
    double windSpeed[NW];   // NW = 9
    for (int n = 0; n < NW; n++) {
        windSpeed[n] = (double) n;
    }

    double Wq = problem->windSpeed * KT_2_MS;
    for (i = 0; i < NW - 1; i++) {
        if (Wq >= windSpeed[i] && Wq < windSpeed[i+1]) break;
    }
    int i0 = (int) i;
    int i1 = (i0+1 < NW) ? i0 + 1 : NW - 1;

    double tx = (Wq - windSpeed[i0]) / (windSpeed[i1] - windSpeed[i0]);
    if (tx < 0) tx = 0; if (tx > 1) tx = 1;

    // Heading change index for deltaCourse \in {-90,-75,...,90}
    double headingChange = node->action->deltaCourse;

    // Round to nearest multiple of 15
    int step = (int) lround(headingChange / 15.0);

    // Range check
    if (step < -6 || step > 6) {
        return NAN;  // invalid
    }

    // Map [-6..6] → [0..12]
    if (step > 0) k = step + 5;
    else k = step + 6;

    // Bilinear interpolation from 4 surrounding bins
    double g00 = problem->ac->gamma_turn[i0][j0][k];
    double g01 = problem->ac->gamma_turn[i0][j1][k];
    double g10 = problem->ac->gamma_turn[i1][j0][k];
    double g11 = problem->ac->gamma_turn[i1][j1][k];

    double g0 = g00*(1.0 - ty) + g01*ty;
    double g1 = g10*(1.0 - ty) + g11*ty;
    double gamma = g0*(1.0 - tx) + g1*tx;

    return fabs(gamma);
}

/*
    Reads flight path angle [deg] datasets
    for straight gliding segments
*/
int readGammaStraight(struct Aircraft *ac, char *ac_name) {
    // Straight best-glide 
    FILE * file;
    char dir[300];
    sprintf(dir,"data/aircraft/gamma_straight.dat");
    file = fopen(dir,"r");

    if (!file) {
        printf("Can't open gamma_straight.dat...\n");
        return EXIT_FAILURE;
    }

    char line[300];
    fgets(line,300,file);
    for (int i = 0; i < 9; fgets(line,300,file), i++) {
        sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&ac->gamma_bg[i][0], &ac->gamma_bg[i][1],
                                                          &ac->gamma_bg[i][2], &ac->gamma_bg[i][3],
                                                          &ac->gamma_bg[i][4], &ac->gamma_bg[i][5],
                                                          &ac->gamma_bg[i][6], &ac->gamma_bg[i][7],
                                                          &ac->gamma_bg[i][8]);
        
    }
    fclose(file);

    // Straight steepest glide corresponding to the maximum flap extended airspeed
    sprintf(dir,"data/aircraft/gamma_straight_100kts.dat");
    file = fopen(dir,"r");

    if (!file) {
        printf("Can't open gamma_straight.dat...\n");
        return EXIT_FAILURE;
    }

    fgets(line,300,file);
    for (int i = 0; i < 9; fgets(line,300,file), i++) {
        sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&ac->gamma_vfe[i][0], &ac->gamma_vfe[i][1],
                                                          &ac->gamma_vfe[i][2], &ac->gamma_vfe[i][3],
                                                          &ac->gamma_vfe[i][4], &ac->gamma_vfe[i][5],
                                                          &ac->gamma_vfe[i][6], &ac->gamma_vfe[i][7],
                                                          &ac->gamma_vfe[i][8]);
        
    }
    fclose(file);

    // Straight optimal glide corresponding to the reference airspeed
    sprintf(dir,"data/aircraft/gammaOpt.dat");
    file = fopen(dir,"r");

    if (!file) {
        printf("Can't open gammaOpt.dat...\n");
        return EXIT_FAILURE;
    }

    fgets(line,300,file);
    for (int i = 0; i < 9; fgets(line,300,file), i++) {
        sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&ac->gamma_opt[i][0], &ac->gamma_opt[i][1],
                                                          &ac->gamma_opt[i][2], &ac->gamma_opt[i][3],
                                                          &ac->gamma_opt[i][4], &ac->gamma_opt[i][5],
                                                          &ac->gamma_opt[i][6], &ac->gamma_opt[i][7],
                                                          &ac->gamma_opt[i][8]);
        
    }
    fclose(file);

	return EXIT_SUCCESS;
}

static int read_exact(void *dst, size_t sz, size_t n, FILE *f) {
    return fread(dst, sz, n, f) == n ? 0 : -1;
}

/*
    Reads flight path angle [deg] datasets
    for turning segments
*/
static int readGammaTurn(Aircraft *ac) {
    if (!ac) return -100;

    // Flight path angle data file
    const char *path = "data/aircraft/viableGammaTurn.bin";

    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "fopen failed for '%s': %s\n", path, strerror(errno));
        return 1;
    }
    if (!f) return -1;

    unsigned char magic[4];
    if (read_exact(magic,1,4,f) || memcmp(magic,"GAMM",4)!=0) { fclose(f); return -2; }

    uint32_t ver, nW, nDir, nCase;
    if (read_exact(&ver,sizeof(uint32_t),1,f) || ver!=1)         { fclose(f); return -3; }
    if (read_exact(&nW, sizeof(uint32_t),1,f))                   { fclose(f); return -4; }
    if (read_exact(&nDir,sizeof(uint32_t),1,f))                  { fclose(f); return -5; }
    if (read_exact(&nCase,sizeof(uint32_t),1,f))                 { fclose(f); return -6; }

    if (nW!=NW || nDir!=NDIR || nCase!=NCASE) { fclose(f); return -7; }

    float *W   = (float*)malloc(nW   * sizeof(float));
    float *Dir = (float*)malloc(nDir * sizeof(float));
    float *dCh = (float*)malloc(nCase* sizeof(float));
    if (!W || !Dir || !dCh) { fclose(f); free(W); free(Dir); free(dCh); return -8; }

    if (read_exact(W,  sizeof(float), nW,   f) ||
        read_exact(Dir,sizeof(float), nDir, f) ||
        read_exact(dCh,sizeof(float), nCase,f)) {
        fclose(f); free(W); free(Dir); free(dCh); return -9;
    }

    // Read data: case-major, then W-major, Dir-fastest (float32)
    // For each case c, for each wind i, read NDIR floats for all directions
    for (uint32_t c=0; c<nCase; ++c) {
        for (uint32_t i=0; i<nW; ++i) {
            float row[NDIR];
            if (read_exact(row, sizeof(float), nDir, f)) { fclose(f); free(W); free(Dir); free(dCh); return -13; }
            for (uint32_t j=0; j<nDir; ++j) {
                ac->gamma_turn[i][j][c] = (double)row[j];  // promote to double
            }
        }
    }

    fclose(f);
    free(W); free(Dir); free(dCh);
    return 0;
}

/*
    Reads the user-defined solver type from the config file
*/
void readSolverType(SearchProblem *problem, const char* cfgdir) {
    FILE * cfg;
    cfg = fopen(cfgdir, "r");
    if (!cfg) { fprintf(stderr, "Config file cannot be opened.\n"); exit(EXIT_FAILURE);}

    char line[300];
    while (fgets(line, sizeof(line), cfg)) {
        if (strncmp(line, "SOLVER", 6) == 0) {
            sscanf(line, "SOLVER = %d", (int*) &problem->solver);
        }
    }
    fclose(cfg);
}

void editCFG(char *origCfg, char *caseFolder, struct Pos *initial, struct Pos *goal, struct Pos *touchdown)
{
    char newCfgPath[300];

    // Construct the new file path inside the case folder
    snprintf(newCfgPath, sizeof(newCfgPath), "%s/aclm.cfg", caseFolder);

    // Open the original CFG file for reading
    FILE *orig = fopen(origCfg, "r");
    if (orig == NULL) {
        perror("Error opening original configuration file");
        return;
    }

    // Open the new CFG file for writing
    FILE *newFile = fopen(newCfgPath, "w");
    if (newFile == NULL) {
        perror("Error creating new configuration file");
        fclose(orig);
        return;
    }

    char line[300];

    // Copy and modify content
    while (fgets(line, sizeof(line), orig)) {
        if (strncmp(line, "INIT_LAT", 8) == 0) {
            fprintf(newFile, "INIT_LAT = %.8f\n", initial->lat);
        } else if (strncmp(line, "INIT_LON", 8) == 0) {
            fprintf(newFile, "INIT_LON = %.8f\n", initial->lon);
        } else if (strncmp(line, "INIT_ALT", 8) == 0) {
            fprintf(newFile, "INIT_ALT = %.2f\n", initial->alt);
        } else if (strncmp(line, "INIT_HDG", 8) == 0) {
            fprintf(newFile, "INIT_HDG = %.1f\n", initial->hdg);
        
        } else if (strncmp(line, "GOAL_LAT", 8) == 0) {
            fprintf(newFile, "GOAL_LAT = %.8f\n", goal->lat);
        } else if (strncmp(line, "GOAL_LON", 8) == 0) {
            fprintf(newFile, "GOAL_LON = %.8f\n", goal->lon);
        } else if (strncmp(line, "GOAL_ALT", 8) == 0) {
            fprintf(newFile, "GOAL_ALT = %.2f\n", goal->alt);
        } else if (strncmp(line, "GOAL_HDG", 8) == 0) {
            fprintf(newFile, "GOAL_HDG = %.1f\n", goal->hdg);
        
        } else if (strncmp(line, "TD_LAT", 6) == 0) {
            fprintf(newFile, "TD_LAT = %.8f\n", touchdown->lat);
        } else if (strncmp(line, "TD_LON", 6) == 0) {
            fprintf(newFile, "TD_LON = %.8f\n", touchdown->lon);
        } else if (strncmp(line, "TD_ALT", 6) == 0) {
            fprintf(newFile, "TD_ALT = %.2f\n", touchdown->alt);
        } else if (strncmp(line, "TD_HDG", 6) == 0) {
            fprintf(newFile, "TD_HDG = %.1f\n", touchdown->hdg);
        }

        else {
            fputs(line, newFile); // Copy other lines unchanged
        }
    }

    // Close files
    fclose(orig);
    fclose(newFile);
}

void copyCFGfile(const char *sourcePath, const char *destinationPath) {
    FILE *sourceFile, *destinationFile;
    char buffer[1024];  // Buffer to store data
    size_t bytesRead;

    // Open the source file
    sourceFile = fopen(sourcePath, "rb");
    if (!sourceFile) {
        perror("Error opening source file");
        return;
    }

    // Open the destination file
    destinationFile = fopen(destinationPath, "wb");
    if (!destinationFile) {
        perror("Error opening destination file");
        fclose(sourceFile);
        return;
    }

    // Copy the file contents
    while ((bytesRead = fread(buffer, 1, sizeof(buffer), sourceFile)) > 0) {
        fwrite(buffer, 1, bytesRead, destinationFile);
    }

    // Close files
    fclose(sourceFile);
    fclose(destinationFile);
}

/*
    Reads and stores binary landing site dataset
*/
landingSite* readLandingSites(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("Failed to open file");
        return NULL;
    }

    landingSite *sites = NULL;

    double buffer[9];
    while (fread(buffer, sizeof(double), 9, fp) == 9) {
        landingSite *node = malloc(sizeof(landingSite));
        if (!node) {
            perror("Memory allocation failed");
            fclose(fp);
            return sites;
        }

        node->ident   = buffer[0];
        node->pos.lat = buffer[1];
        node->pos.lon = buffer[2];
        node->pos.alt = buffer[3];
        node->pos.hdg = buffer[4];
        node->length  = buffer[5];
        node->width   = buffer[6];
        node->commercial = buffer[7];
        node->military = buffer[8];
        
        // Insert at the front
        node->next = sites;
        sites = node;
    }

    fclose(fp);
    return sites;
}

/*
    Utility-based landing site selection
*/
int selectLandingSite(SearchProblem *problem) {

    // Optimal flight path angle
    double course = 0;
    double *gammaArray = getOptimalGamma(&course, problem);
    problem->gammaBG = gammaArray[0];
    problem->gammaVFE = gammaArray[1];
    problem->gammaOpt = gammaArray[2];
    free(gammaArray); gammaArray = NULL;

    // Allocate memory for nearest landing sites and read them
    problem->nearest = malloc(sizeof(struct landingSite));
    char filename[300];
    strcpy(filename, "lib/aclm/data/landing_sites.bin");
    problem->nearest = readLandingSites(filename);
    if (!problem->nearest) return EXIT_FAILURE;

    // Temporary landing site structure for iteration
    landingSite *options;
    options = problem->nearest;
    
    // Get absolute maximums and minimums of the landing site set
    double maxRunwayLength = 0, maxRunwayWidth = 0, maxGamma = 0, maxHeadwind = 0;
    struct Pos approachFix;
    while (options) {

        // Get 1 NM approach fix
        double d = 1, course = wrapTo360(options->pos.hdg - 180);
        geo_npos(&options->pos, &approachFix, &d, &course, &problem->GeoOpt);
        approachFix.alt = options->pos.alt + (d*NM_2_FT)*tan(problem->gammaOpt*DEG_2_RAD);
        
        // Great-circle distance to the landing site
        geo_dist(&problem->Initial, &approachFix, &d, &course, &problem->GeoOpt);

        // Altitude difference
        double dh = problem->Initial.alt - approachFix.alt;

        // Estimate straight gliding angle
        options->gamma = atan(dh/(d*NM_2_FT))*RAD_2_DEG;

        // Relative wind speeds
        options->headwind = problem->windSpeed*cos(wrapTo360(options->pos.hdg - problem->windDirection)*DEG_2_RAD);
        options->crosswind = problem->windSpeed*sin(wrapTo360(options->pos.hdg - problem->windDirection)*DEG_2_RAD);

        if (options->gamma > maxGamma) maxGamma = options->gamma;
        if (options->length > maxRunwayLength) maxRunwayLength = options->length;
        if (options->width > maxRunwayWidth) maxRunwayWidth = options->width;
        if (options->headwind > maxHeadwind) maxHeadwind = options->headwind;

        options = options->next;
    }
    
    // Count reachable sites
    options = problem->nearest;
    int nReachable = 0;
    while (options) {
        if (options->gamma >= problem->gammaBG) nReachable++;
        options = options->next;
    }
    problem->numRanked = nReachable;

    // If none of the flight path angles to the landing sites are feasible, do not continue.
    if (!nReachable) return EXIT_FAILURE;

    // Allocate array of pointers
    problem->rankedSites = malloc(nReachable * sizeof(landingSite*));

    // Allocate array for runtime (calloc for 0 initialization)
    problem->searchRuntime = calloc(nReachable, sizeof(double));

    // Populate and compute utility
    options = problem->nearest;
    int idx = 0;
    double maxU = 0;
    while (options) {
        if (options->gamma < problem->gammaBG) {
            options->U = -1;  // Mark unreachable
        } else {
            double Ugamma   = options->gamma / maxGamma;
            double Ulength  = options->length / maxRunwayLength;
            double Uwidth   = options->width / maxRunwayWidth;
            double Ucom     = 1 - options->commercial;
            double Umil     = 1 - options->military;
            double Uwind    = options->headwind/maxHeadwind;

            // Total utility
            options->U = (0.5*Ugamma + 0.05*Ulength + 0.05*Uwidth + 0.15*Ucom + 0.15*Umil + 0.1*Uwind); // Version 2

            problem->rankedSites[idx++] = options;
        }
        options = options->next;
    }
    
    // Sort landing sites based on the utility
    qsort(problem->rankedSites, problem->numRanked, sizeof(landingSite*), cmp);

    return EXIT_SUCCESS;
}

// Sort by utility descending
int cmp(const void *a, const void *b) {
    double ua = (*(landingSite**)a)->U;
    double ub = (*(landingSite**)b)->U;
    return (ua < ub) - (ua > ub);  // descending
}

/*
    Initializes distance parameters related to a chosen landing site
*/
void initDistanceParams(SearchProblem * problem, int siteIndex){

    // Assign touchdown state
    problem->Touchdown = problem->rankedSites[siteIndex]->pos;

    // Get bearing from initial state to the landing site
    double d, course;
    geo_dist(&problem->Initial, &problem->rankedSites[siteIndex]->pos, &d, &course, &problem->GeoOpt);

    // Optimal flight path angle to the bestlanding site
    double *gammaArray = getOptimalGamma(&course, problem);
    problem->gammaBG = gammaArray[0];
    problem->gammaVFE = gammaArray[1];
    problem->gammaOpt = gammaArray[2];
    free(gammaArray); gammaArray = NULL;

    // Compute 1 NM approach fix of the best landing site
    d = 1;
    course = wrapTo360(problem->rankedSites[siteIndex]->pos.hdg - 180);
    geo_npos(&problem->rankedSites[siteIndex]->pos, &problem->Goal, &d, &course, &problem->GeoOpt);
    problem->Goal.alt = problem->rankedSites[siteIndex]->pos.alt + (d*NM_2_FT) * tan(problem->gammaBG*DEG_2_RAD);
    problem->Goal.hdg = problem->rankedSites[siteIndex]->pos.hdg;

    // If the default optimal gamma is infeasible, set it to the gamma to the best landing site
    if (problem->rankedSites[siteIndex]->gamma < problem->gammaOpt) {
        problem->gammaOpt = problem->rankedSites[siteIndex]->gamma;
        problem->limitedR = 1;
        printf("- Limited reachability detected.");
    }

    // Altitude to lose
    problem->alt_to_lose = problem->Initial.alt -  problem->Goal.alt;

    // Best-glide range
    problem->bestGlideRange = problem->alt_to_lose / tan(problem->gammaBG*DEG_2_RAD);
}


double groundRisk(SearchProblem *problem, struct DubinsPath *dubins, char *filename, int type)
{   

    // Pre-allocate sampling coordinate structure
    int numSamples;
    struct Pos * coordinates = NULL;
    
    // Get samples based on path type (0: Search-based, 1: Dubins-based)
    if (!type) {
        searchIntegrationCoordinates(problem, problem->ground_dx, &coordinates, &numSamples);
    } else {
        if (dubins->size < 2){
            coordinates = getDubinsCoordinates(dubins, problem->ground_dx, &numSamples, &problem->GeoOpt, NULL);
        } else {

            // Get S-Turn coordinates
            int sampleSize1, sampleSize2;
            struct Pos *sturn1 = getDubinsCoordinates(&dubins[0], problem->ground_dx, &sampleSize1, &problem->GeoOpt, NULL); // Initial waypoint to the intermediate waypoint
            struct Pos *sturn2 = getDubinsCoordinates(&dubins[1], problem->ground_dx, &sampleSize2, &problem->GeoOpt, NULL); // Intermediate waypoint to the final waypoint

            // Concatenate samples
            numSamples = sampleSize1 + sampleSize2 - 1;
            coordinates = (struct Pos *) malloc(numSamples*sizeof(struct Pos));
            for (int i = 0; i < numSamples; i++) {
                if (i < sampleSize1) {
                    coordinates[i] = sturn1[i];
                } else {
                    coordinates[i] = sturn2[i - sampleSize1 + 1];
                }
            }
            free(sturn1); sturn1 = NULL;
            free(sturn2); sturn2 = NULL;
        }
    }

    // Initialize population density samples
    double * c = (double *) malloc(numSamples * sizeof(double));
    if (!c) {
        perror("malloc failed for c");
        free(coordinates);
        return NAN;
    }

    // Write coordinates to a file
    FILE * file;
    if (filename){
        char dir[200];
        sprintf(dir, "%s/census_data_%d.csv", filename, type);
        file = fopen(dir, "w");
    }
    
    // Calculate population cost
    double density;
    double I = 0;
    for (int i = 0; i < numSamples; ++i) {  
        
        // Sample population density below the crossover altitude
        if (coordinates[i].alt < problem->crossoverAlt){
            density = getPopulationDensityRTree(problem->census, coordinates[i].lat, coordinates[i].lon, problem->maxPopulation);
            c[i] = (density - problem->minPopulation)/(problem->maxPopulation - problem->minPopulation);
        } else {
            c[i] = 0;
        }

        // Assert the sampled density is within the closed interval [0,1]
        if (c[i] < 0 || c[i] > 1) {
            printf("Population sample: %.4f\n", c[i]);
            free(coordinates); coordinates = NULL;
            free(c); c = NULL;
            exit(EXIT_FAILURE);
        } else if (isnan(c[i])) {
            c[i] = 0;
        }

        // Trapezoidal integration
        if (i > 0) {
            I += 0.5 * (c[i - 1] + c[i]) * problem->groundRiskTimeStep;
        }

        if (filename) fprintf(file, "%.6f, %.6f, %.6f, %.6f, %6f\n", coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, c[i], I);
    }

    if (filename) fclose(file);

    // Flight time
    if (!type) problem->flightTime_gp = ((numSamples - 1) * problem->ground_dx) / problem->ac->airspeed;

    // Free memory
    free(coordinates); coordinates = NULL;
    free(c); c = NULL;

    // Overflown population risk
    return I;
}

struct DubinsPath *minRiskDubins(struct DubinsPath * dubins, SearchProblem * problem, double * totalRuntime, double * riskRuntime)
{
    // Initialize Dubins structure
    Traj_InitArray(dubins->traj, 4);
    dubins->traj[0].wpt.pos = problem->InitialSearch;
    dubins->traj[0].wpt.hdg = problem->InitialSearch.hdg;
    dubins->traj[3].wpt.pos = problem->Goal;
    dubins->traj[3].wpt.hdg = problem->Goal.hdg;

    dubins->traj[0].wpt.gam = -problem->ac->gammaOptTurn;
    dubins->traj[1].wpt.gam = -problem->gammaOpt;
    dubins->traj[2].wpt.gam = -problem->ac->gammaOptTurn;

    dubins->traj[0].wpt.rad = problem->ac->turnRadius*FT_2_NM;
    dubins->traj[1].wpt.rad = 0;
    dubins->traj[2].wpt.rad = problem->ac->turnRadius*FT_2_NM;
        
    // Compute the initial min. risk Dubins path
    dubins->feasible = true;
    struct DubinsPath *bestDubins = getMinRiskDubinsPath(dubins, problem, totalRuntime, riskRuntime);

    return bestDubins;
}

double _airspaceRisk(SearchProblem *problem, struct DubinsPath *dubins, char *filename, int type) {   

    // Pre-allocate sampling coordinate structure
    int numSamples;
    struct Pos * coordinates = NULL;
    
    // Get samples based on path type (0: Search-based, 1: Dubins-based)
    if (!type) {
        searchIntegrationCoordinates(problem, problem->airspace_dx, &coordinates, &numSamples);
    } else {
        if (dubins->size < 2){
            coordinates = getDubinsCoordinates(dubins, problem->airspace_dx, &numSamples, &problem->GeoOpt, NULL);
        } else {

            // Get S-Turn coordinates
            int sampleSize1, sampleSize2;
            struct Pos *sturn1 = getDubinsCoordinates(&dubins[0], problem->airspace_dx, &sampleSize1, &problem->GeoOpt, NULL); // Initial waypoint to the intermediate waypoint
            struct Pos *sturn2 = getDubinsCoordinates(&dubins[1], problem->airspace_dx, &sampleSize2, &problem->GeoOpt, NULL); // Intermediate waypoint to the final waypoint

            // Concatenate samples
            numSamples = sampleSize1 + sampleSize2 - 1;
            coordinates = (struct Pos *) malloc(numSamples*sizeof(struct Pos));
            for (int i = 0; i < numSamples; i++) {
                if (i < sampleSize1) {
                    coordinates[i] = sturn1[i];
                } else {
                    coordinates[i] = sturn2[i - sampleSize1 + 1];
                }
            }
            free(sturn1); sturn1 = NULL;
            free(sturn2); sturn2 = NULL;
        }
    }


    // Initialize airspace risk sampling array
    double * c = (double *) malloc(numSamples * sizeof(double));
    if (!c) {
        perror("malloc failed for c");
        free(coordinates);
        return NAN;
    }

    // Write integration details to a file
    FILE * file;
    if (filename){
        char dir[200];
        sprintf(dir, "%s/airspace_samples.csv", filename);
        file = fopen(dir, "w");
    }
    
    // Calculate airspace risk
    double density;
    double I = 0;
    for (int i = 0; i < numSamples; ++i) {   

        // Sample airspace density
        double traffic_cost = getAirTrafficDensity(&coordinates[i], problem->traffic);
        double heli_cost = getAirTrafficDensity(&(coordinates[i]), problem->heli);
        double prohibitedArea_cost = getAirTrafficDensity(&(coordinates[i]), problem->prohibited);

        if (prohibitedArea_cost > 0.5) {
            I = INFINITY;
            break;
        } else {
            c[i] = 0.5*traffic_cost + 0.25*heli_cost + 0.25*prohibitedArea_cost;
        }

        // Assert the sampled density is within the closed interval [0,1]
        if (c[i] < 0) {
            printf("airspaceRisk type %d - Airspace density sample: %.4f\n", type, c[i]);
            free(coordinates); coordinates = NULL;
            free(c); c = NULL;
            exit(EXIT_FAILURE);
        } else if (isnan(c[i])) {
            c[i] = 0;
        }
    
        // Trapezoidal integration
        if (i > 0) {
            I += (0.5 * (c[i - 1] + c[i]) * problem->airspaceRiskTimeStep);
        }
        
        if (filename) fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, traffic_cost, heli_cost, prohibitedArea_cost, c[i], I);
    }
        
    if (filename) fclose(file);

    // Free memory
    free(coordinates); coordinates = NULL;
    free(c); c = NULL;

    // Time averaging
    if (!type) problem->flightTime_ga = problem->airspaceRiskTimeStep*numSamples;

    return I;
}
 
/*
    Trapezoidal integration for airspace risk computation
    Uses fixed integration time step of 1 second
    Discretizes path, samples risk, and integrates
    Integration details are written in "airspaceRisk0.csv" for search-based path
    Integration details are written in "airspaceRisk1.csv" for Dubins-based path
    Conflict predictions are written in "conflicts0.csv" for search-based path
    Conflict predictions are written in "conflicts1.csv" for Dubins-based path
*/
double airspaceRisk(SearchProblem *problem, struct DubinsPath *dubins, char *filename, int type){   

    // Pre-allocate sampling coordinate structure
    int numSamplesSearch, numSamplesDubins, numSamples;
    struct Pos * searchcoord = NULL; // Stores search-based path coordinates
    struct Pos * coordinates = NULL; // Stores the coordinates of the entire path (Search + final turn Dubins or minimum Risk Dubins including S-turn)
    double * groundspeedSearch = NULL;
    
    // Sampling frequency [s]
    double dt_sample = 1;
    if (!type) {
        searchIntegrationCoordinatesWithFixedTimeStep(problem, dt_sample, &searchcoord, &groundspeedSearch, &numSamplesSearch); // Search-based path samples
        struct Pos *finalTurn = getDubinsCoordinatesWithFixedTimeStep(problem->finalTurn, NULL, &numSamplesDubins, &problem->GeoOpt, NULL, problem); // Final turn Dubins path samples

        // Concatenate samples
        numSamples = numSamplesSearch + numSamplesDubins - 1;
        coordinates = (struct Pos *) malloc(numSamples*sizeof(struct Pos));
        for (int i = 0; i < numSamples; i++) {
            if (i < numSamplesSearch) { 
                coordinates[i] = searchcoord[i];
            } else {
                coordinates[i] = finalTurn[i - numSamplesSearch + 1];
            }
        }
        free(finalTurn); finalTurn = NULL;

    } else {
        if (dubins->size < 2){
            coordinates = getDubinsCoordinatesWithFixedTimeStep(dubins, NULL, &numSamples, &problem->GeoOpt, NULL, problem);
        } else {

            // Get S-Turn coordinates
            int sampleSize1, sampleSize2;
            struct Pos *sturn1 = getDubinsCoordinatesWithFixedTimeStep(&dubins[0], NULL, &sampleSize1, &problem->GeoOpt, NULL, problem); // Initial waypoint to the intermediate waypoint
            struct Pos *sturn2 = getDubinsCoordinatesWithFixedTimeStep(&dubins[1], &sturn1[sampleSize1-1], &sampleSize2, &problem->GeoOpt, NULL, problem); // Intermediate waypoint to the final waypoint

            // Concatenate samples
            numSamples = sampleSize1 + sampleSize2 - 1;
            coordinates = (struct Pos *) malloc(numSamples*sizeof(struct Pos));
            for (int i = 0; i < numSamples; i++) {
                if (i < sampleSize1) {
                    coordinates[i] = sturn1[i];
                } else {
                    coordinates[i] = sturn2[i - sampleSize1 + 1];
                }
            }
            free(sturn1); sturn1 = NULL;
            free(sturn2); sturn2 = NULL;
        }
    }

    if (numSamples > problem->dynTraffic->nsteps){
        free(coordinates);
        return NAN;
    } 

    // Initialize airspace risk sampling array
    double * c = (double *) malloc(numSamples * sizeof(double));
    if (!c) {
        perror("malloc failed for c");
        free(coordinates);
        return NAN;
    }

    // Write integration details to a file
    FILE * file = NULL;
    if (filename){
        char dir[200];
        snprintf(dir, sizeof(dir), "%s/airspaceRisk%d.csv", filename, type);
        file = fopen(dir, "w");
        if (!file) {
            perror("fopen airspaceRisk.csv failed.\n");
        }
        fprintf(file, "Lat [deg], Lon [deg], Alt [ft], Groundspeed [kts], JCPA, Integral JCPA\n");
    }

    FILE * file_conflict = NULL;
    if (filename){
        char dir[200];
        snprintf(dir, sizeof(dir), "%s/conflicts%d.csv", filename, type);
        file_conflict = fopen(dir, "w");
        if (!file_conflict){
            perror("fopen conflicts.csv failed.\n");
        }
        fprintf(file_conflict, "Time [s], OLat [deg], OLon [deg], OAlt [ft], OCourse [deg], ILat [deg], ILon [deg], IAlt [ft], ICourse [deg], ID, JCPA\n");
    }

    // Initialize time index
    int timeIDX = 0;

    // Integrate CPA along the action
    double I = 0;
    for (int i = 0; i < numSamples; i++){

        // Compute the absolute maximum CPA cost considering nearby traffic
        double Jmax = 0;
        for (size_t aircraftID = 0; aircraftID < (size_t) problem->dynTraffic->nflights; aircraftID++){

            // Compute CPA for air traffic with ID: aircraftID
            double J = JCPA(&coordinates[i], aircraftID, timeIDX, problem);

            // Check if it's higher, if so assign it to the maximum cost
            if (J > Jmax) Jmax = J;

            // If J > 0 log conflict
            if (filename && (J > 0)){
                FSSample s = problem->dynTraffic->flights[aircraftID].sample[timeIDX];
                fprintf(file_conflict, "%d, %.6f, %.6f, %.2f, %.2f, %.6f, %.6f, %.2f, %.2f, %zu, %.6f\n",\
                                        timeIDX, coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, coordinates[i].hdg,\
                                        s.pos.lat, s.pos.lon, s.pos.alt, s.pos.hdg, aircraftID, J);
            }
        }

        // Store the cost
        c[i] = Jmax;

        // Update time index
        timeIDX++;

        // Trapezoidal integration
        if (i > 0) {
            I += 0.5 * (c[i - 1] + c[i])*dt_sample;
        }

        // Log airspace risk samples and costs
        if (file){
            double gs = (groundspeedSearch ? groundspeedSearch[i] : NAN);
            fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, gs, c[i], I);
        }
    }
        
    if (filename){
        fclose(file);
        fclose(file_conflict);
    }

    // Free memory
    free(coordinates); coordinates = NULL;
    free(c); c = NULL;

    // Time averaging
    if (!type) problem->flightTime_ga = problem->airspaceRiskTimeStep*numSamples;

    return I;
}


/*
    DYNAMIC MEMORY HANDLING FUNCTIONS
*/
void freePath(struct Path *path) {
    if (!path) return;

    // Free each action in the nodes array
    if (path->nodes) {
        for (size_t i = 0; i < path->depth; i++) {
            if (path->nodes[i].action != NULL) {
                free(path->nodes[i].action);
                path->nodes[i].action = NULL;
            }
        }

        // Free the node array itself
        free(path->nodes);
        path->nodes = NULL;
    }

    // Finally, free the path structure
    free(path);
}

void freeHeatmap(struct airspaceHeatmap * hm) {

    if (!hm) return;

    if (hm->heatmap){
        free(hm->heatmap);
        hm->heatmap = NULL;
    }

    if (hm->lat_grid){
        free(hm->lat_grid);
        hm->lat_grid = NULL;
    }

    if (hm->lon_grid){
        free(hm->lon_grid);
        hm->lon_grid = NULL;
    }

    if (hm->alt_grid){
        free(hm->alt_grid);
        hm->alt_grid = NULL;
    }

    free(hm); //hm = NULL;
}

void freeDynamicTraffic(struct DynamicTraffic *t){

    if (!t) return;

    // Free datafile
    if (t->datafile) {
        free(t->datafile);
        t->datafile = NULL;
    }

    // Free flights
    if (t->flights) {

        // Free samples
        for (uint32_t i = 0; i < t->nflights; i++) {
            if (t->flights[i].sample) {
                free(t->flights[i].sample);
                t->flights[i].sample = NULL;
            }
        }
        free(t->flights);
        t->flights = NULL;
    }

    free(t);
}


void freeRankedSites(struct landingSite **rankedSites, int count) {
    if (!rankedSites) return;

    free(rankedSites); rankedSites = NULL;
}

void freeProblem(SearchProblem *problem) {

    if (!problem) return;

    // printPointerAddress(problem);

    // Free simple pointers
    if (problem->ac) {
        free(problem->ac); problem->ac = NULL;
    }
   
    if (problem->goalNormal) {
        free(problem->goalNormal);
        problem->goalNormal = NULL;
    }
    if (problem->noFlyZoneVerticesLon) {
        free(problem->noFlyZoneVerticesLon);
        problem->noFlyZoneVerticesLon = NULL;
    }
    if (problem->noFlyZoneVerticesLat) {
        free(problem->noFlyZoneVerticesLat);
        problem->noFlyZoneVerticesLat = NULL;
    }
    if (problem->altitudeModification) {
        free(problem->altitudeModification);
        problem->altitudeModification = NULL;
    }

    // Census pointers
    if (problem->census) {
        freeCensus(problem->census);
        problem->census = NULL;
    }

    if (problem->merged) {
        freeCensus(problem->merged);
        problem->merged = NULL;
    }

    // Lists
    destroyPriorityQueue(problem->pq);
    problem->pq = NULL;
    
    destroyClosedList(problem->explored);
    problem->explored = NULL;

    // Free DubinsPath structs if dynamically allocated
    if (problem->path2Hold) {
        free(problem->path2Hold);
        problem->path2Hold = NULL;
    }

    if (problem->finalTurn) {
        free(problem->finalTurn);
        problem->finalTurn = NULL;
    }

    // Free Dubins options
    if (problem->DubinsOpt) {
        free(problem->DubinsOpt);
        problem->DubinsOpt = NULL;
    }

    // Free path node structure
    if (problem->path) {
        if (problem->path->nodes) {
            for (size_t i = 0; i < problem->path->depth; ++i) {
                if (problem->path->nodes[i].action) {
                    if (problem->path->nodes[i].action != NULL) {
                        free(problem->path->nodes[i].action);
                        problem->path->nodes[i].action = NULL;
                    }
                }
            }
            free(problem->path->nodes);
            problem->path->nodes = NULL;
        }
    
        free(problem->path);
        problem->path = NULL;
    }

    // Free Action array
    if (problem->pathActions) {
        free(problem->pathActions);
        problem->pathActions = NULL;
    }

    // Free landing site list
    while (problem->nearest) {
        struct landingSite *temp = problem->nearest;
        problem->nearest = problem->nearest->next;
        free(temp); temp = NULL;
    }
    if (problem->numRanked){
        free(problem->rankedSites); problem->rankedSites = NULL;
    }
    
    // Free heatmaps
    if (problem->traffic) {
        freeHeatmap(problem->traffic);
        problem->traffic = NULL;
    }
    
    if (problem->heli) {
        freeHeatmap(problem->heli);
        problem->heli = NULL;
    }
    
    if (problem->prohibited) {
        freeHeatmap(problem->prohibited);
        problem->prohibited = NULL;
    }

    // Free dynamic traffic
    if (problem->dynTraffic){
        freeDynamicTraffic(problem->dynTraffic);
        problem->dynTraffic = NULL;
    }
    
    // printPointerAddress(problem);

    // Free the SearchProblem structure
    free(problem);

    // problem = NULL;
}

void printPointerAddress(SearchProblem *problem){
    printf("%-50s %p\n", "Address problem->ac:", problem->ac);
    printf("%-50s %p\n", "Address problem->goalNormal:", problem->goalNormal);
    printf("%-50s %p\n", "Address problem->noFlyZoneVerticesLon:", problem->noFlyZoneVerticesLon);
    printf("%-50s %p\n", "Address problem->noFlyZoneVerticesLat:", problem->noFlyZoneVerticesLat);
    printf("%-50s %p\n", "Address problem->altitudeModification:", problem->altitudeModification);
    printf("%-50s %p\n", "Address problem->census:", problem->census);
    printf("%-50s %p\n", "Address problem->merged:", problem->merged);
    printf("%-50s %p\n", "Address problem->pq:", problem->pq);
    printf("%-50s %p\n", "Address problem->explored:", problem->explored);
    printf("%-50s %p\n", "Address problem->path2Hold:", problem->path2Hold);
    printf("%-50s %p\n", "Address problem->finalTurn:", problem->finalTurn);
    printf("%-50s %p\n", "Address problem->DubinsOpt:", problem->DubinsOpt);
    printf("%-50s %p\n", "Address problem->path:", problem->path);
    printf("%-50s %p\n", "Address problem->pathActions:", problem->pathActions);
    printf("%-50s %p\n", "Address problem->nearest:", problem->nearest);
    printf("%-50s %p\n", "Address problem->rankedSites:", problem->rankedSites);
    printf("%-50s %p\n", "Address problem->traffic:", problem->traffic);
    printf("%-50s %p\n", "Address problem->heli:", problem->heli);
    printf("%-50s %p\n", "Address problem->prohibited:", problem->prohibited);
    printf("%-50s %p\n", "Address problem->dynTraffic:", problem->dynTraffic);
}
