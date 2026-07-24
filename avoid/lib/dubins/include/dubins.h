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

#ifdef __cplusplus
extern "C" {
#endif

#include "traj.h"
#include "math_utils.h"
#include <limits.h>

typedef struct SearchProblem SearchProblem;

struct DubinsPath;

enum DubinsType{
    OPT, 
    RSR, 
    RSL, 
    LSR, 
    LSL, 
    RLR1, 
    RLR2,
    LRL1,
    LRL2
};

struct DubinsPath{
    struct Traj traj[4];
    struct Pos orbit1;
    struct Pos orbit2;
    enum DubinsType type;
    double dgamma;
    double hdist;           // 2D path length [NM]
    int size;               // 1: Conventional Dubins Path, 2: S-turn Dubins Path
    double gp;              // Overflown population risk
    double ga;              // Airspace occupation risk
    double risk;            // Total ground and airspace risks
    bool feasible;          // False: Infeasible flight path angle
};

struct DubinsOpt{
    int verbose;
    struct TrajOpt trajopt;
};


int Dubins(struct DubinsPath *dubins, const struct DubinsOpt *opt);

int shortestDubins(struct DubinsPath * dubins, const struct DubinsOpt * opt);

void getIntermediateWaypoint(struct Pos *interWaypoint0,
                            struct Pos *interWaypoint,
                            double straightLength,
                            double theta,
                            int extendTo,
                            struct GeoOpt *GeoOpt);

struct DubinsPath *computeSturnDubins(struct DubinsPath *dubins,
                                    double dAltitude,
                                    int extendTo,
                                    SearchProblem * problem);

struct Pos *getDubinsCoordinates(struct DubinsPath *dubins,
                                double interval, 
                                int *sampleSize,
                                struct GeoOpt *GeoOpt,
                                char * folderName);

struct Pos *getDubinsCoordinatesWithFixedTimeStep(struct DubinsPath *dubins,
                                                struct Pos * initialSample,
                                                int *sampleSize,
                                                const struct GeoOpt *GeoOpt,
                                                const double groundspeed);

struct DubinsPath *getBestSturnPath(struct DubinsPath *dubins,
                                    double dAltitude,
                                    SearchProblem *problem,
                                    double *riskRuntime);

int reduceSlope(struct DubinsPath *dubins, double dAltitudeTarget, SearchProblem *problem);

void moveWaypoint(struct DubinsPath *dubins, SearchProblem *problem);

void evaluateDubinsRisk(struct DubinsPath *dubins,
                        SearchProblem * problem);

struct DubinsPath *getMinRiskDubinsPath(struct DubinsPath *dubins,
                                        SearchProblem *problem,
                                        double *totalRuntime,
                                        double *riskRuntime);

int getDubinsWithType(struct DubinsPath * dubins, int type, const struct DubinsOpt * opt);

void updateOptimalGammaTurn(struct DubinsPath * dubins, SearchProblem * problem);

#ifdef __cplusplus
}
#endif