/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   #
#    This file is part of Gradient-guided Search for Assured Contingency Landing    #
#    Management.                                                                    #
#                                                                                   #
#    Gradient-guided Search for Assured Contingency Landing Management is free      #
#    software: you can redistribute it and/or modify it under the terms of the GNU  #
#    General Public License as published by the Free Software Foundation, either    #
#    version 3 of the License, or (at your option) any later version.               #
#                                                                                   #
#    This program is distributed in the hope that it will be useful,                #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                   #
#    GNU General Public License for more details.                                   #
#                                                                                   #
#    You should have received a copy of the GNU General Public License              #
#    along with this program. If not, see <https://www.gnu.org/licenses/>.          #
#                                                                                   #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#pragma once

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "math_utils.h"
#include "geo.h"

typedef struct FSSample {   // Flight series data sample
    struct Pos pos;    // Sample position and heading: Lat [deg], lon [deg], alt [ft], track [hdg]
    float   gs;       // Groundspeed [kts]
    int64_t ts;       // Timestamp [unix]
    bool valid;         // Flag indicating sample validity
} FSSample;

typedef struct FlightSeries {
    FSSample * sample;  // Flight series
    uint8_t type;
} FlightSeries;

typedef struct DynamicTraffic {
    char * datafile;        // Binary data file directory
    char label[10];         // Data label (e.g., IMC_LOW)
    uint32_t nflights;      // Total number of flights
    uint32_t nsteps;        // Total number of time steps
    uint32_t median;        // Median simultaneously airborne traffic
    uint32_t p90;           // .9 quantile simultaneously airborne traffic
    uint32_t peak;          // Peak simultaneously airborne traffic
    FlightSeries * flights; // Length total number of flight within data file
} DynamicTraffic;

typedef struct airspaceHeatmap {
    int LAT_SIZE;
    int LON_SIZE;
    int ALT_SIZE;
    size_t HM_SIZE;
    float * heatmap;
    float * lat_grid;
    float * lon_grid;
    float * alt_grid;
    int grid_on_cells; // 0 = point-centered (helicost), 1 = cell-centered (traffic_density)
    char * gridfile;
    char * datafile;
} airspaceHeatmap;

struct surface {
    struct PosXYZ xyz[4];      // First window
    struct PosXYZ centroid;    // Surface centroid coordinates
    struct PosXYZ normal;      // Surface normal vector coordinates
    double *k;           // The constant in the surface equation (e.g. ax + by + cz + 'k' = 0)
};

struct segment {
    struct surface surfaces[6];     // Surfaces of a segment
    struct segment *next;           // Next segment
};

struct airCorridor {
    char ident[50];
    struct segment *segment;
	struct airCorridor *next;
};

typedef struct SearchProblem SearchProblem;
typedef struct Node Node;

/* Functions to exploit air traffic density heatmap*/
void readHeatMap(airspaceHeatmap *trafficHeatmap);
float getAirTrafficDensity(struct Pos *pos, airspaceHeatmap *trafficHeatmap);

/* Analytical geometry to find the minimum distance in-between a 3D point and a set of polyhedron */
int readAirCorridors(struct airCorridor *corridors, char * filedir);

/**
    @brief Vector operations to find the relative position of a three-dimensional ECEF coordinate to air corridors
    isInsidePolyhedron : Checks if a given ECEF coordinate is inside any of the air traffic corridors
*/
int isInsidePolyhedron(struct airCorridor *corridors, struct PosXYZ pos);
int doesOverlap(struct PosXYZ pos, struct surface *surface);

/*
    Calculates angles between point-to-vertex vectors and edge vectors
*/
void point2vertexVec_to_edgeVector_angles(struct PosXYZ *pos, struct surface *surface, double angles[8]);

/*
    Performs exhaustive search to find the minimum distance
    from a given coordinate to air traffic corridors
*/
double minimum_distance(struct Pos *pos, struct airCorridor *corridors);

/*
    Air corridor occupation cost based on point-to-polyhedron distance
    and buffer distance
*/
double airCorridorOccupationCost(struct Pos *pos, struct airCorridor *corridors);

void loadDynamicTraffic(SearchProblem * problem);

double groundspeed(struct Pos * state, SearchProblem * problem);
int traversaldeltat(struct Node * node, SearchProblem * problem);
void VCPA(struct Pos * state, int aircraftID, int timeIDX, double * activeDeltav, int * activeCount, SearchProblem *problem);
double CPACost(double deltah, double deltav, SearchProblem * problem);
double JCPA(struct Pos * state, int aircraftID, int timeIDX, SearchProblem *problem);
void loadDummyIntruder(SearchProblem *problem);
