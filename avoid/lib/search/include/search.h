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

#include "aclm.h"
#include "geo.h"
#include "node.h"
#include "priorityQ.h"
#include "conversions.h"
#include "census.h"
#include "math_utils.h"
#include "dubins.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

typedef struct Node Node;
typedef struct closedList closedList;
typedef struct PriorityQueue PriorityQueue;
typedef struct Aircraft Aircraft;

// Action stucture
typedef struct Action {
    double deltaCourse; // deg
    double length;      // ft
    double gamma;       // deg
    double gamma_turn;  // deg
} Action;

void initializeSearch(SearchProblem * problem);
Action *getActions(struct Node * node, SearchProblem *problem);
struct Pos *result(Node *node, Action *action, SearchProblem *problem);
double *distanceCost(struct Pos *state, double length, SearchProblem *problem);
double courseCost(struct Pos *state, SearchProblem *problem);
double populationCost(struct Pos *state, Action *action, double length, SearchProblem *problem);
double airspaceCost(struct Pos *state, Action *action, double length, SearchProblem *problem);
double maxDensityAhead(struct Pos *state, SearchProblem *problem);
double dynamicAirspaceCost(struct Node * parent, struct Node * child, SearchProblem *problem);
bool finalTurnDynamicTrafficCost(const Node *node, SearchProblem *problem);
double stateCost(Node *node, SearchProblem * problem);
void courseCostNormalization(SearchProblem *problem);
void populationCostNormalizer(SearchProblem *problem);
bool dubinsFinalTurn(struct Node *node, SearchProblem * problem);
bool goalTest(Node *node, SearchProblem * problem);
bool modifyGamma(double remaining_distance, double * gammaArray, Node * node, SearchProblem *problem);
bool isInHexagon(struct Pos *state1, struct Pos *state2, double circumradius, SearchProblem * problem);
bool doesCluster(struct Node * node, SearchProblem * problem);
bool doesClusterFrontier(struct Node * node, SearchProblem * problem);
int runSearch(SearchProblem * problem);
int writeResults(SearchProblem *problem, struct DubinsPath *bestDubins, char *folderName);
void searchIntegrationCoordinates(SearchProblem * problem, double interval, struct Pos **coordinates, int *numSamples);
void searchIntegrationCoordinatesWithFixedTimeStep(SearchProblem *problem,
                                                   double dt_sample,
                                                   struct Pos **coordinates,
                                                   double **gs,
                                                   int *numSamples);
void computeCorrectedAltitudes(const struct Path * path, int depth, const Aircraft * ac, double *out_alt);
double altitudeWeight(struct Pos *state, SearchProblem *problem);

#ifdef __cplusplus
}
#endif

