/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   #
#    Gradient-guided Search for Assured Contingency Landing Management              #
#    Copyright (C) 2025  Huseyin Emre Tekaslan                                      #
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

#pragma once

#ifdef __cplusplus
extern "C"{
#endif

#include "search.h"
#include "geo_structs.h"
#include "airtraffic.h"
#include "math_utils.h"
#include "node.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

// struct GeoOpt;
typedef struct PriorityQueue PriorityQueue;
typedef struct closedList closedList;
typedef struct Node Node;

// Grid sizes of the optimal gamma data for straight flight
#define NDIR_STR 9

// Grid sizes of the optimal turn gamma data
#define NW   9                 // wind speeds 0:1:8 m/s
#define NDIR 24                // relWindDeg -180:15:165
#define NCASE 12               // deltaChi {-90:15:90} ~~~excluding {0}

// Grid sizes of static air traffic data
#define NLAT 31              
#define NLON 31
#define NALT 11

typedef struct Aircraft {
    double turnRadius;      // ft
    double airspeed;        // ft/s
    double gammaOptTurn;    // [deg]
    double gammaBGturn;     // [deg]
    double gammaOpt;        // [deg]
    double gamma_bg[NW][NDIR_STR];  // [deg]
    double gamma_vfe[NW][NDIR_STR]; // [deg]
    double gamma_opt[NW][NDIR_STR];  // [deg]
    double gamma_turn[NW][NDIR][NCASE]; // [deg] -> [windSpeed][relWind][deltaChiCase]
} Aircraft;

// Solver type structure
enum SolverType
{
    SEARCH,
    DUBINS_SEARCH,
    SEARCH_AIRSPACE,
    SEARCH_AIRSPACE_DYN
};

// Landing site structure
typedef struct landingSite
{
    int ident;      // Ident. integer corresponding ICAO code
    struct Pos pos; // Lat (deg), Lon (deg), Alt (ft), Hdg (deg)
    double gamma;   // deg – estimate straight gliding angle to the landing site
    double length;  // ft
    double width;   // ft
    double commercial; // Commercial traffic support 0:False, 1:True
    double military; // Military field 0:False, 1:True
    double headwind; // Headwind [kts]
    double crosswind; // Tailwind [kts]
    double U;       // Utility
    struct landingSite *next;   // Next landing site
} landingSite;

// Problem structure
typedef struct SearchProblem
{
    // EARTH MODEL
    struct GeoOpt GeoOpt;

    // DATASETS
    struct Aircraft *ac;
    struct Census *census;  // For population density query
    struct Census *merged;  // For holding point identification
    struct airspaceHeatmap *traffic;  // Departure/arrival air traffic heatmap based on historical ADS-B data from OpenSky Network
    struct airspaceHeatmap *heli;     // Helicopter route heatmap based on skyvector.com DC Heli Area
    struct airspaceHeatmap *prohibited; // Prohibited zone heatmap based on skyvector.com
    struct DynamicTraffic  *dynTraffic; // Dynamic air traffic based on historical ADS-B data from OpenSky Network

    // SOLVER SETTINGS
    size_t maxIter;
    enum SolverType solver;
    int exitFlag;          // 0: Success, 1: Unreachable, 2: Max. iteration is reached, 3: Initialized in prohibited area, 9: Unidentified error
    int limitedR;
    int direction_of_search; // -1: Backward, 1: Forward

    // SEARCH
    PriorityQueue *pq;
    closedList *explored;

    // STATES
    struct Pos Initial;         // Initial emergency state
    struct Pos Goal;            // Goal state
    struct Pos Touchdown;       // Touchdown state
    struct Pos InitialSearch;   // Initial search state
    struct Pos finalState;      // Final search state

    // ACTIONS
    double deltaCourse[5];      // Course angle changes [deg]
    double lmin;                // Min. path length [ft]
    double adaptiveLenghtAltitude;  // Floor altitude for adaptive segment length [ft]

    // DISCRETIZATION
    double dAltitude;           // Altitude discretization [ft]
    double dPsi;                // Heading discretization [deg]
    int dTime;                  // Time discretization [s]

    // DISTANCE - LENGTH RELATED
    double alt_to_lose;         // Altitude difference of initial and goal states [ft]
    double totalLength;         // Total path length [ft]
    double dOpt;                // Optimal path length [ft]
    double bestGlideRange;      // Best glide range [ft]
    double identRadius;         // Solution identification radius [ft]
    double du;                  // Distance-related constants for weight calculations [ft]
    double dl;                  // Distance-related constants for weight calculations [ft]

    // AIRCRAFT - FLIGHT ENVELOPE
    double v;                   // Airspeed [ft/s]
    double gammaBG;             // Best-glide flight path angle [deg]
    double gammaVFE;            // Steepest glide angle corresponding to VFE [deg]
    double gammaOpt;            // Optimal glide angle [deg]
    
    // ENVIRONMENTAL
    double windSpeed;           // Wind speed [kts]
    double windDirection;       // Wind direction [deg]

    // COST COMPUTATION
    double airspaceRiskTimeStep;    // Time step for airspace risk integration [s]
    double groundRiskTimeStep;      // Time step for ground risk integration [s]
    double airspace_dx;             // Spatial step for airspace risk integration [ft]
    double ground_dx;               // Spatial step for overflown population risk integration [ft]
    double maxPopulation;       // Max. population density within census data. Used for normalization [persons/m^2]
    double minPopulation;       // Min. population density within census data. Used for normalization [persons/m^2]
    double crossoverAlt;        // Crossover altitude / Ceiling altitude for population density consideration. [ft]
    double *goalNormal;         // Goal state normal vector pointing the goal course angle
    double courseCostMax;       // Max. course angle cost for normalization
    double courseCostMin;       // Min. course angle cost for normalization
    double populationCostNormalizer; // Overflown population cost normaliation parameter
    double goalAirspaceCost;    // Airspace cost of the goal state
    double initAirspaceCost;    // Airspace cost of the initial state
    double w_gp;                // Overflown population risk weight. Default: (1 - w_gp)
    double w_ga;                // Airspace risk weight 
    double asrisk_w_hmax;       // Airspace risk weight h_max
    double asrisk_w_hmin;       // Airspace risk weight h_min

    // Separation related
    double RNP_TER;         // Terminal area RNP requirement [NM]
    double RNP_FIN;         // Final approach RNP requirement [NM]
    double IAF_Alt;         // Initial approach fix altitude [NM] 
    double ERR_ADSB_P;      // ADS-B lateral position error [NM]
    double ERR_ADSB_V;      // ADS-B velocity error [m/s]
    double DELAY_ADSB;      // Uncompensated ADS-B delay [s] requirement by FAA § 91.227 
    double ERR_BARO_A5000;  // Barometric pressure max. error above 5000 ft MSL requirement [ft]
    double ERR_BARO_B5000;  // Barometric pressure max. error below 5000 ft MSL requirement [ft]
    double VERT_CLEARANCE;  // Vertical clearance requirement [ft]
    double HORZ_CLEARANCE;  // Horizontal clearance requirement [NM]

    // HOLDING PATTERN
    struct Pos *HoldingPoint;
    struct Pos *holdInbound;
    struct Pos *holdOutbound;
    double maxObstacleHeight;
    double minHoldArea;         // Minimum holding area (based on the turn radius: pi*R^2)
    double minHoldOutboundAltitude; // Minimum holding outbound altitude
    double *noFlyZoneVerticesLon; // No fly zone longitude coordinates
    double *noFlyZoneVerticesLat; // No fly zone latitude coordinates
    
    // DUBINS AND SEARCH PATHS
    struct DubinsPath *path2Hold;
    struct DubinsPath *finalTurn;   
    struct DubinsOpt *DubinsOpt;
    struct Path *path;
    struct Action *pathActions;
    double overflownPop;
    double airspaceOccup;
    double dGamma;                  // Flight path angle modification (Deviation from the optimal glide angle) [deg]
    double *altitudeModification;    // Altitude changes based on the flight path angle modification [ft]
    double flightTime_gp;
    double flightTime_ga;
    
    // LANDING SITES
    int siteIndex;                  // Landing site ordered index
    struct landingSite *nearest;    // Nearest landing sites
    struct landingSite *bestSite;       // The best landing site based on a utility function
    struct landingSite **rankedSites;
    int numRanked;

    // RUNTIME
    double * searchRuntime;         // Search runtime array (length numRanked) [ms]
    double totalSearchRuntime;      // Total search runtime [ms]
    double totalRiskRuntime;        // Total risk runtime of search only [ms]
    double totalHoldPlanRuntime;    // Total holding pattern identification and 3D Dubins computation runtime [ms]
    double holdPlanRiskRuntime;     // Total risk runtime of holding point identification only [ms]
    
} SearchProblem;

void setUpProblem(SearchProblem *problem, const char* cfgdir);
void readSolverType(SearchProblem *problem, const char* cfgdir);
void loadHeatmaps(SearchProblem *problem);
void loadInitialState(SearchProblem *problem, const char* cfgdir);
void loadGoalState(SearchProblem *problem, const char* cfgdir);
void loadTouchdownState(SearchProblem *problem, const char* cfgdir);
void loadAircraftParams(SearchProblem *problem, const char* cfgdir);
void loadEarthModelParams(SearchProblem *problem, const char* cfgdir);
char *loadRiskComputationParams(SearchProblem *problem, const char* cfgdir);
void loadSearchParams(SearchProblem *problem, const char* cfgdir);
void loadDynamicTrafficParams(SearchProblem *problem, const char* cfgdir);
double *getOptimalGamma(double *course, SearchProblem *problem);
double getOptimalGammaTurn(Node *node, SearchProblem *problem);
double getViableGamma(Node *node, SearchProblem *problem);
int readGammaStraight(struct Aircraft *ac);
int readGammaTurn(Aircraft *ac);
void editCFG(char *origCfg, char *caseFolder, struct Pos *initial, struct Pos *goal, struct Pos *touchdown);
int runEmergencyPlanning(SearchProblem *problem, char * casefolder, char * configfile);
void copyCFGfile(const char *sourcePath, const char *destinationPath);
landingSite* readLandingSites(const char *filename);
int selectLandingSite(SearchProblem *problem);
int cmp(const void *a, const void *b);
void initDistanceParams (SearchProblem * problem, int siteIndex);
struct DubinsPath *minRiskDubins(struct DubinsPath *dubins, SearchProblem * problem, double * totalRuntime, double * riskRuntime);

double groundRisk(SearchProblem *problem, struct DubinsPath *dubins, char *filename, int type);
double _airspaceRisk(SearchProblem *problem, struct DubinsPath *dubins, char *filename, int type);
double airspaceRisk(SearchProblem *problem, struct DubinsPath *dubins, char *filename, int type);

void freeHeatmap(struct airspaceHeatmap * hm);
void freeProblem(SearchProblem *problem);
void freePath(struct Path *path);
void freeRankedSites(struct landingSite **rankedSites);
void freeDynamicTraffic(struct DynamicTraffic *t);
void printPointerAddress(SearchProblem *problem);

#ifdef __cplusplus
}
#endif
