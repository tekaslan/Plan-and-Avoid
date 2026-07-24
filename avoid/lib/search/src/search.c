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


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   %
#    Main Gradient-guided Search Functions                                          %
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : H. Emre Tekaslan (tekaslan@vt.edu)                                   %
#    Date    : January 2026                                                         %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "../include/search.h"

/*
    Runs gradient-guided search
*/
int runSearch(SearchProblem * problem){

    // Initialize variables
    size_t counter = 0;
    Node *node;
    Node **children;
    problem->exitFlag = 9;

    // Start the timer
    clock_t begin, end;
    begin = clock();

    // While the frontier is not empty and max. iteration has not been reached yet
    while (!pqIsEmpty(problem->pq) && counter <= problem->maxIter){  
    // while (!pqIsEmpty(problem->pq) && counter <= 15){


        // printf("Iter: %zu\n", counter);

        // Pop a state
        node = pqPop(problem->pq);
        // printf("Popped-> lat: %.3f lon: %.3f alt: %.1f hdg: %.1f, f: %.4f\n", node->state.lat, node->state.lon, node->state.alt, node->state.hdg, node->f);

        // Append the state to the closed list
        clInsert(node, problem);
        counter++;

        // Goal test
        if (goalTest(node, problem)) {

            // Terminate the timer
            end = clock();
            problem->totalSearchRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;

            // Extract the solution
            printf("- Solved.");
            problem->exitFlag = 0;
            problem->path = solution(node, &node->depth, problem);
            problem->path->depth = node->depth;
            problem->finalState = node->state;
            problem->pathActions = solutionActions(node, &node->depth, problem);
            problem->totalLength = node->parent->l*FT_2_NM + problem->finalTurn->hdist;

            return EXIT_SUCCESS;
        } else if (counter == problem->maxIter) {
            end = clock();
            problem->totalSearchRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;
            problem->exitFlag = 2;
            return EXIT_FAILURE;
        }

        // Get children states
        children = expandNode(node, problem);
        
        // Add children to the frontier
        for (int i = 0; i < 5; i++) {
            if (children[i]->f < INFINITY && !doesCluster(children[i], problem) && !doesClusterFrontier(children[i], problem)){
                pqInsert(children[i], problem->pq);
            } else {
                destroyNode(children[i]);
            }
        }
    }

    // Empty frontier flag
    if (pqIsEmpty(problem->pq)) problem->exitFlag = -3;
    end = clock();
    problem->totalSearchRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;

    return EXIT_FAILURE;
}

/*
    Initializes the gradient-guided search
*/
void initializeSearch(SearchProblem * problem) {

    // Modified flight path angle
    problem->dGamma = 0;

    // Distance related parameters
    problem->dOpt = problem->alt_to_lose / tan(fabs(problem->ac->gammaOptTurn)*DEG_2_RAD);
    problem->dl = 0.5*NM_2_FT;
    double dmin0 = minDistance(&problem->Initial, &problem->Goal, &problem->GeoOpt, problem);
    problem->du = 0.5*(problem->dOpt + dmin0);
    
    // Calculate cost normalizers
    populationCostNormalizer(problem);
    courseCostNormalization(problem);

    // Create a frontier with the priority queue
    size_t capacity = 50000;
    problem->pq = createPriorityQueue(capacity, 1);

    // Create a closed list
    problem->explored = createClosedList(capacity);

    // Make the initial search state a node
    Node *node = (Node *) malloc(sizeof(Node));
    Node * parent_node = NULL;
    double * gammaArray = getOptimalGamma(&problem->Initial.hdg, problem);
    Action action = {.deltaCourse=0, .length=0, .gamma=gammaArray[2], .gamma_turn=problem->ac->gammaBGturn};
    free(gammaArray); gammaArray = NULL;
    node = createNode(&problem->InitialSearch, parent_node, &action);
    node->f = stateCost(node, problem);
    node->t = 0;
    
    // If the initial state is in the prohibited area, return error
    double prohibitedAreaInitialization = getAirTrafficDensity(&node->state, problem->prohibited);
    if (prohibitedAreaInitialization > 0) {
        problem->exitFlag = -2;
        return;
    }

    // Airspace cost of the initial and goal states
    problem->initAirspaceCost = getAirTrafficDensity(&problem->InitialSearch, problem->traffic);
    problem->goalAirspaceCost = getAirTrafficDensity(&problem->Goal, problem->traffic);
    
    // Check if the initial state is the goal state
    if (goalTest(node, problem)) {
        printf("Solved.\n");
    }
        
    // Insert the initial state into the frontier
    pqInsert(node, problem->pq);

    // Compute the search space size
    size(&problem->InitialSearch, &problem->Goal, problem->dOpt, problem->lmin, problem->dAltitude, problem->dPsi);

    // Reset timers
    problem->totalRiskRuntime = 0;
    problem->totalSearchRuntime = 0;
}


/*
    Returns a set of feasible actions for a given state
*/
Action *getActions(struct Node * node, SearchProblem *problem)
{
    // Get distance and bearing to the goal
    double dist;
    double course;
    geo_dist(&node->state, &problem->Goal, &dist, &course, &problem->GeoOpt);

    Action *action_list = (Action *) malloc(5*sizeof(Action));
    if (action_list == NULL) {
        fprintf(stderr, "Memory allocation for action_list failed.\n");
        return NULL;  // Handle memory allocation failure
    }
    
    for (int i = 0; i < 5; i++) {

        // Set course change
        action_list[i].deltaCourse = problem->deltaCourse[i];

        // Set segment length
        action_list[i].length = problem->lmin;
        if (node->state.alt > problem->adaptiveLenghtAltitude) {
            if (node->ga > 0 || maxDensityAhead(&node->state, problem)) action_list[i].length = problem->lmin;
            // else action_list[i].length = problem->lmin + 0.1*(node->state.alt - problem->adaptiveLenghtAltitude); // Version 2
            else action_list[i].length = problem->lmin + 0.5*(node->state.alt - problem->adaptiveLenghtAltitude); // Version 3
        }

        // New course
        double newCourse = wrapTo360(node->state.hdg + action_list[i].deltaCourse);

        // Get the straight segment flight path angle
        double *gammaArray = getOptimalGamma(&newCourse, problem);
        if (problem->limitedR) action_list[i].gamma = gammaArray[0];
        else action_list[i].gamma = gammaArray[2];
        free(gammaArray); gammaArray = NULL;

        // Get the flight path angle for turning
        {
            if (action_list[i].deltaCourse != 0.0) {
                // Make a scratch copy of the parent node
                Node tmp_node = *node;

                // Build a scratch Action that represents this candidate child action
                Action tmpA = (Action){0};
                if (node->action) {
                    // Start from parent's action so any fields getOptimalGammaTurn may rely on are present
                    tmpA = *node->action;
                }
                // Overwrite with the candidate child’s parameters
                tmpA.deltaCourse = action_list[i].deltaCourse;
                tmpA.length      = action_list[i].length;
                tmpA.gamma       = action_list[i].gamma;
                tmpA.gamma_turn  = 0.0;

                // Point the scratch node at the scratch action
                tmp_node.action = &tmpA;

                // Compute the turn gamma for this candidate
                action_list[i].gamma_turn = getOptimalGammaTurn(&tmp_node, problem);
            } else {
                action_list[i].gamma_turn = 0.0;
            }
        }
    }
    return action_list;
}

/*
    Returns child state given a parent state, action,
    and aircraft parameters (turn radius, glide angles).

    Per-segment altitude accounting:
      - include exit half of previous turn
      - include entry half of current turn 
*/
struct Pos *result(Node *node, Action *action, SearchProblem *problem) {

    // Allocate new state
    struct Pos *newState = (struct Pos *) malloc(sizeof(struct Pos));
    if (newState == NULL) {
        fprintf(stderr, "Memory allocation for newState failed.\n");
        return NULL;
    }

    // Geodetic propagation distance (NM)
    double distance_nm = action->length * FT_2_NM;

    // New heading
    double newCourse = wrapTo360(node->state.hdg + action->deltaCourse);
    newState->hdg = newCourse;

    // New position
    geo_npos(&node->state, newState, &distance_nm, &newCourse, &problem->GeoOpt);

    if (node->depth == 0) {
        const double h0         = node->state.alt;
        const double gamma_str  = action->gamma;                // deg
        newState->alt = h0 - action->length * tan(gamma_str  * DEG_2_RAD);
    } else {

        // --- Altitude propagation with "defer-then-correct" entry half and exit half ---
        const double h0         = node->state.alt;                 // altitude stored at node (with prior straight assumption)
        const double R          = problem->ac->turnRadius;         // [ft]
        const double gamma_str  = action->gamma;                   // [deg] straight glide for THIS leg
        const double gamma_turn = action->gamma_turn;       // [deg] glide in turns

        // Half of THIS action's heading change [rad]
        const double half = 0.5 * fabs(action->deltaCourse) * DEG_2_RAD;

        // Entry half that actually occurred at the END of the previous leg:
        const double s_entry = R * half;       // arc length [ft]
        const double d_entry = R * tan(half);  // chord eaten from straight [ft]

        // What straight γ did you assume on the previous leg?
        const double gamma_prev_str = (node->action ? node->action->gamma : gamma_str);

        // 1) Revert prior straight assumption over the entry chord and
        // 2) Apply the actual entry-arc sink at gamma_turn.
        // This gives the true altitude at the node.
        const double h_parent_corr = h0
            + d_entry * tan(gamma_prev_str * DEG_2_RAD)   // add back straight sink you had over-counted
            - s_entry * tan(gamma_turn * DEG_2_RAD);  // apply correct turn sink

        // 3) Now fly THIS leg: start with the EXIT half of the same turn, then straight
        const double s_exit = R * half;       // exit arc length [ft]
        const double d_exit = R * tan(half);  // exit chord eaten from this leg's straight [ft]

        double l_str = action->length - d_exit;     // remaining straight in this leg

        newState->alt = h_parent_corr
                    - s_exit * tan(gamma_turn * DEG_2_RAD)
                    - l_str  * tan(gamma_str  * DEG_2_RAD);
    }

    return newState;
}

/*
    Searches the maximum glide range for the maximum
    and minimum course cost values that will be used
    for normalization
*/
void courseCostNormalization(SearchProblem *problem) 
{   

    // Define the resolution and create arrays
    size_t grid_size = 100; 
    double LAT[grid_size], LON[grid_size];
    double DOT[grid_size][grid_size], DMIN[grid_size][grid_size], HDG_COST_MAP[grid_size][grid_size];
    struct Pos pos_geod;
    struct PosXYZ pos_ned;
    pos_geod.alt = problem->Goal.alt;
    
    // Compute the bounding box
    struct Pos *N = (struct Pos *) malloc(sizeof(struct Pos));
    struct Pos *S = (struct Pos *) malloc(sizeof(struct Pos));
    struct Pos *E = (struct Pos *) malloc(sizeof(struct Pos));
    struct Pos *W = (struct Pos *) malloc(sizeof(struct Pos));

    double const maxRange = 1.1*problem->bestGlideRange*FT_2_NM;
    double const radials[4] = {0.0, 90.0, 180.0, 270.0};
    geo_npos(&problem->InitialSearch, N, &maxRange, &radials[0], &problem->GeoOpt);
    geo_npos(&problem->InitialSearch, S, &maxRange, &radials[2], &problem->GeoOpt);
    geo_npos(&problem->InitialSearch, E, &maxRange, &radials[1], &problem->GeoOpt);
    geo_npos(&problem->InitialSearch, W, &maxRange, &radials[3], &problem->GeoOpt);

    // Generate the grid of latitude and longitude
    double dmin, course;
    for (size_t i = 0; i < grid_size; i++) {
        pos_geod.lat = S->lat + (N->lat - S->lat) * i / (grid_size - 1);
        LAT[i] = pos_geod.lat;
        for (size_t j = 0; j < grid_size; j++){
            pos_geod.lon = W->lon + (E->lon - W->lon) * j / (grid_size - 1);
            LON[j] = pos_geod.lon;

            // Convert geodetic to NED
            geo_lla2ned(&problem->Goal, &pos_geod, &pos_ned);

            DOT[i][j] = fmax(0, pos_ned.Y * problem->goalNormal[0] + pos_ned.X * problem->goalNormal[1]);

            geo_dist(&pos_geod, &problem->Goal, &dmin, &course, &problem->GeoOpt);
            DMIN[i][j] = dmin*NM_2_FT;
        }
    }

    FILE * file;
    char filename[50];
    strcpy(filename, "tmp/courseCostMapLAT.csv");
    file = fopen(filename, "w");

    FILE * file2;
    strcpy(filename, "tmp/courseCostMapLON.csv");
    file2 = fopen(filename, "w");

    FILE * file3;
    strcpy(filename, "tmp/courseCostMap.csv");
    file3 = fopen(filename, "w");
    
    // Compute the course angle cost map
    for (size_t i = 0; i < grid_size; i++) {
        for (size_t j = 0; j < grid_size; j++) {
            HDG_COST_MAP[i][j] = problem->dOpt * DOT[i][j] / DMIN[i][j];
            if (j < grid_size-1)
            {
                fprintf(file, "%.6f,", LAT[i]);
                fprintf(file2, "%.6f,", LON[j]);
                fprintf(file3, "%.6f,", HDG_COST_MAP[i][j]);
            }
            else
            {
                fprintf(file, "%.6f\n", LAT[i]);
                fprintf(file2, "%.6f\n", LON[j]);
                fprintf(file3, "%.6f\n", HDG_COST_MAP[i][j]);
            }  
        }
    }
    fclose(file);
    fclose(file2);
    fclose(file3);

    // Find the max and min values for normalization
    double hdg_cost_max = HDG_COST_MAP[0][0];
    double hdg_cost_min = HDG_COST_MAP[0][0];
    for (size_t i = 0; i < grid_size; i++) {
        for (size_t j = 0; j < grid_size; j++) {
            if (HDG_COST_MAP[i][j] > hdg_cost_max) {
                hdg_cost_max = HDG_COST_MAP[i][j];
            }
            if (HDG_COST_MAP[i][j] < hdg_cost_min) {
                hdg_cost_min = HDG_COST_MAP[i][j];
            }
        }
    }

    // Store the max and min values in the problem struct
    problem->courseCostMax = hdg_cost_max;
    problem->courseCostMin = hdg_cost_min;

    // Free allocated memory
    free(N); N = NULL;
    free(S); S = NULL;
    free(E); E = NULL;
    free(W); W = NULL;
}

/*
    Calculates the population cost normalizer
*/
void populationCostNormalizer(SearchProblem *problem)
{
    problem->populationCostNormalizer = (4.0/27.0) * problem->dOpt * tan(problem->gammaOpt*DEG_2_RAD) / problem->InitialSearch.alt;
}

/*
    Returns the distance cost as a 2D array
*/
double *distanceCost(struct Pos *state, double length, SearchProblem *problem)
{

    double dmin = minDistance(state, &problem->Goal, &problem->GeoOpt, problem);

    double w;
    if (dmin >= problem->du) {
        w = 0;
    } else if ((dmin > problem->dl) & (dmin < problem->du)) {
        w = fabs(problem->du - dmin)/(problem->du - problem->dl);
    } else {
        w = 1;
    }

    double *g_d = (double *) malloc(2*sizeof(double));
    g_d[0] = w*fabs(length + dmin - problem->dOpt)/problem->dOpt;
    g_d[1] = (1-w)*min(1, dmin/problem->dOpt);

    return g_d;
}

/*
    Returns the course angle cost of a given state
*/
double courseCost(struct Pos *state, SearchProblem *problem)
{

    // Compute altitude weight
    // If 3000 ft above the goal altitude,
    // do NOT apply course angle cost to let search to freely explore the space
    double hmax = problem->Goal.alt + 3000; 
    if (hmax - state->alt < 0) return 0;
    double w = (hmax - state->alt)/(hmax);

    // Straight line distance from the state to the goal
    double dmin;
    double course;
    geo_dist(state, &problem->Goal, &dmin, &course, &problem->GeoOpt);
    dmin *= NM_2_FT;

    // Coordinate transformation
    struct PosXYZ state_ned;
    geo_lla2ned(&problem->Goal, state, &state_ned);

    // Dot product
    double v[2] = {state_ned.Y, state_ned.X};
    double dotProd = problem->goalNormal[0]*v[0] + problem->goalNormal[1]*v[1];
    double rho = max(0, dotProd);

    // Normalized course cost
    double g_chi = (rho - problem->courseCostMin)/(problem->courseCostMax - problem->courseCostMin);
    double ftol = 1.1;
    if (g_chi > ftol) {
        printf("g_chi: %f\n", g_chi);
        perror("Course cost is out of bounds.\n");
    } else if (g_chi > 1) {
        g_chi = 1;
    }

    // Distance weight
    g_chi *= problem->dOpt/dmin;

    // Altitude weight
    g_chi *= w;
    
    return g_chi;
}

/*
    Returns the overflown population cost of a segment
    given an initial state and an action using trapezoidal
    integration
*/
double populationCost(struct Pos *state, Action *action, double length, SearchProblem *problem)
{

    if (state->alt > problem->crossoverAlt) return 0;

    // Extract variables
    double lat = state->lat;
    double lon = state->lon;
    double alt = state->alt;
    double hdg = state->hdg;

    // Initialize variables
    int n = (int) floor(action->length/problem->ground_dx);
    double latArray[n];
    double lonArray[n];
    double altArray[n];
    double lengthArray[n];
    double c[n];
    double w1[n];
    double w2[n];
    double direction = wrapTo360(hdg + action->deltaCourse);

    // Generate sampled path
    struct Pos tmp_pos;
    tmp_pos.lat = lat;
    tmp_pos.lon = lon;
    tmp_pos.alt = alt;
    struct Pos new_pos;
    double stepLength = problem->ground_dx*FT_2_NM;
    double len;
    for (int i = 0; i < n; ++i) {
        len = i*stepLength;
        geo_npos(&tmp_pos, &new_pos, &len, &direction, &problem->GeoOpt);
        latArray[i] = new_pos.lat;
        lonArray[i] = new_pos.lon;
        altArray[i] = tmp_pos.alt - i*problem->ground_dx * tan(action->gamma*DEG_2_RAD);
        lengthArray[i] =  length + i*problem->ground_dx;
    }

    // Calculate population cost
    double density;
    double I = 0;
    for (int i = 0; i < n; ++i) {

        // Sample population density
        density = getPopulationDensityRTree(problem->census, latArray[i], lonArray[i], problem->maxPopulation);
        c[i] = (density - problem->minPopulation)/(problem->maxPopulation - problem->minPopulation);

        // Assert the sampled density is within the closed interval [0,1]
        if (c[i] < 0 || c[i] > 1) {
            printf("Population sample: %.4f\n", c[i]);
            exit(EXIT_FAILURE);
        } else if (isnan(c[i])) {
            c[i] = 0;
        }

        // Distance weight w1
        w1[i] = pow(1 - lengthArray[i] / problem->dOpt, 2);

        // Assert the weight is within the closed interval [0,1]
        if (w1[i] < 0 || w1[i] > 1) {
            printf("Population weight w1: %.4f\n", w1[i]);
            exit(EXIT_FAILURE);
        }

        // Altitude weight w2
        if (problem->InitialSearch.alt >= problem->crossoverAlt) {
            w2[i] = (altArray[i] >= problem->crossoverAlt) ? 0 : (1 - altArray[i] / problem->crossoverAlt);
        } else {
            w2[i] = 1 - altArray[i] / problem->InitialSearch.alt;
        }

        // Assert the weight is within the closed interval [0,1]
        if (w2[i] < 0 || w2[i] > 1) {
            printf("Population weight w2: %.4f\n", w2[i]);
            exit(EXIT_FAILURE);
        }

        // Trapezoidal integration
        if (i > 0) {
            I += 0.5 * (w1[i - 1] * w2[i - 1] * c[i - 1] + w1[i] * w2[i] * c[i]) * problem->groundRiskTimeStep;
        }
    }

    // Time averaging
    double T = action->length / (problem->ac->airspeed);
    if (T <= 0) return 0;
    else return I / (problem->populationCostNormalizer * T);
}

/*
    Altitude-dependent weight for airspace cost
*/
double altitudeWeight(struct Pos *state, SearchProblem *problem){

    // Get altitude bounds as a function of goal state airspace cost for weighting
    double w_altmax, w_altmin;
    w_altmax = problem->asrisk_w_hmax;
    w_altmin = problem->asrisk_w_hmin;
    
    // Remaining altitude to the goal
    double dh = state->alt - problem->Goal.alt;

    // Altitude weight
    double w_alt;
    if (dh >= w_altmax) return w_alt = 1;
    else if (dh >= w_altmin && dh < w_altmax) {w_alt = (dh - w_altmin)/(w_altmax - w_altmin); return pow(w_alt,2);}
    else return 0;
}

/*
    Airspace corridor occupation cost
*/
double airspaceCost(struct Pos *state, Action *action, double length, SearchProblem *problem) {

    // Get altitude weight first
    double walt = 1;
    if (problem->goalAirspaceCost > 0){
        walt = altitudeWeight(state, problem);
        if (walt == 0) return 0;
    }

    // Extract variables
    double lat = state->lat;
    double lon = state->lon;
    double alt = state->alt;
    double hdg = state->hdg;

    // Initialize variables
    int n = (int)floor(action->length/problem->airspace_dx);
    double latArray[n];
    double lonArray[n];
    double altArray[n];
    double lengthArray[n];
    double c[n];
    double direction = wrapTo360(hdg + action->deltaCourse);

    // Generate sampled path
    struct Pos tmp_pos;
    tmp_pos.lat = lat;
    tmp_pos.lon = lon;
    tmp_pos.alt = alt;
    struct Pos new_pos;
    double stepLength = problem->airspace_dx*FT_2_NM;
    double len;
    for (int i = 0; i < n; ++i) {
        len = i*stepLength;
        geo_npos(&tmp_pos, &new_pos, &len, &direction, &problem->GeoOpt);
        latArray[i] = new_pos.lat;
        lonArray[i] = new_pos.lon;
        altArray[i] = tmp_pos.alt - i*problem->airspace_dx * tan(action->gamma*DEG_2_RAD);
        lengthArray[i] =  length + i*problem->airspace_dx;
    }

    // Calculate airspace cost
    double I = 0;
    for (int i = 0; i < n; ++i) {
        // Sample airspace density
        new_pos.lat = latArray[i];
        new_pos.lon = lonArray[i];
        new_pos.alt = altArray[i];
        double traffic_cost = getAirTrafficDensity(&new_pos, problem->traffic);
        double heli_cost = getAirTrafficDensity(&new_pos, problem->heli);
        double prohibitedArea_cost = getAirTrafficDensity(&new_pos, problem->prohibited);

        if (prohibitedArea_cost > 0.5) {
            return INFINITY;
        } else {
            c[i] = 0.5*traffic_cost + 0.25*heli_cost + 0.25*prohibitedArea_cost;
        }

        // Assert the sampled density is within the closed interval [0,1]
        if (c[i] < 0) {
            printf("airspaceCost - Airspace density sample: %.4f\n", c[i]);
            exit(EXIT_FAILURE);
        } else if (isnan(c[i])) {
            c[i] = 0;
        }

        // Trapezoidal integration
        if (i > 0) {
            I += 0.5 * (c[i - 1] + c[i]) * problem->airspaceRiskTimeStep;
        }
    }

    // Time averaging
    double T = action->length / (problem->ac->airspeed);
    if (T <=  0) return 0;
    else return I * walt;
}

/*
    Estimates the average airspace density ahead in a cone
*/
double maxDensityAhead(struct Pos *state, SearchProblem *problem)
{
    double walt = 1;
    if (problem->goalAirspaceCost > 0){
        walt = altitudeWeight(state, problem);
        if (walt == 0) return 0;
    }

    // Extract position variables
    double lat = state->lat;
    double lon = state->lon;
    double alt = state->alt;
    double hdg = state->hdg;

    int n_rays = 5;           // Number of directions to sample in the cone
    int n_points = 1;         // Number of points along each ray
    int n_level = 3;          // Number of rays in the vertical plane
    double cone_angle = 60;   // Total cone spread in degrees
    double stepLength = 10*problem->lmin*FT_2_NM;      // [NM]
    double max_density = 0;
    double sum_density = 0;

    double start_angle = wrapTo360(hdg - cone_angle/2.0);
    double angle_step = cone_angle / (n_rays - 1);

    struct Pos tmp, new_pos;
    tmp.lat = lat;
    tmp.lon = lon;
    tmp.alt = alt;

    for (int i = 0; i < n_rays; i++) {
        double dir = wrapTo360(start_angle + i * angle_step);
        for (int j = 1; j <= n_points; j++) {
            double len = j * stepLength;
            for (int k = 0; k < n_level; k++) {
                double gamma = problem->gammaBG + k*(10 - problem->gammaBG)/n_level;
                geo_npos(&tmp, &new_pos, &len, &dir, &problem->GeoOpt);
                new_pos.alt = alt - j * stepLength * NM_2_FT * tan(gamma * DEG_2_RAD);

                double traffic = getAirTrafficDensity(&new_pos, problem->traffic);
                double heli = getAirTrafficDensity(&new_pos, problem->heli);
                double prohibited = getAirTrafficDensity(&new_pos, problem->prohibited);

                if (prohibited > 0.5) return INFINITY;

                double cost = 0.5*traffic + 0.25*heli + 0.25*prohibited;
                sum_density += cost;
                if (cost > max_density) max_density = cost;
            }
        }
    }
    
    return sum_density/(n_rays * n_points * n_level)*walt;
}

// /*
//     Dynamic airspace closes point of approach cost
//     integration over time for search segment
// */
// double dynamicAirspaceCost(struct Node * parent, struct Node * child, SearchProblem *problem){

//     // Get action traversal time
//     double gs = groundspeed(&child->state, problem);
//     double deltat = child->action->length/(gs * KTS_2_FTS);

//     // Time step guard
//     int t0 = (int)lround(parent->t);
//     int N  = (int)ceil(deltat);

//     int maxN = problem->dynTraffic->nsteps - t0;
//     if (maxN <= 0) return 0.0;
//     if (N > maxN) N = maxN;

//     // Extract variables
//     double lat = parent->state.lat;
//     double lon = parent->state.lon;
//     double alt = parent->state.alt;
//     double hdg = parent->state.hdg;

//     // Initialize arrays
//     double latArray[N];
//     double lonArray[N];
//     double altArray[N];
//     double c[N];
//     double direction = wrapTo360(hdg + child->action->deltaCourse);

//     // Rate of change of altitude
//     double altitudeRate = (parent->state.alt - child->state.alt)/deltat;

//     // Generate sampled path
//     struct Pos tmp_pos;
//     tmp_pos.lat = lat;
//     tmp_pos.lon = lon;
//     tmp_pos.alt = alt;
//     struct Pos new_pos;
//     double dt = 1;
//     double stepLength = (gs * KTS_2_FTS) * dt *FT_2_NM;
//     double len;
//     for (int i = 0; i < N; ++i){
//         len = i*stepLength;
//         geo_npos(&tmp_pos, &new_pos, &len, &direction, &problem->GeoOpt);
//         latArray[i] = new_pos.lat;
//         lonArray[i] = new_pos.lon;
//         altArray[i] = alt - altitudeRate * i;
//     }

//     // Initialize time index
//     int timeIDX = (int) round(parent->t);

//     // Integrate CPA along the action
//     double I = 0;
//     for (int i = 0; i < N; i++){

//         // Sample ownship state
//         new_pos.lat = latArray[i];
//         new_pos.lon = lonArray[i];
//         new_pos.alt = altArray[i];

//         // Compute the absolute maximum CPA cost considering nearby traffic
//         double Jmax = 0;
//         for (uint32_t aircraftID = 0; aircraftID < problem->dynTraffic->nflights; aircraftID++){

//             // Compute CPA for air traffic with ID: aircraftID
//             double J = JCPA(&new_pos, aircraftID, timeIDX, problem);

//             // Check if it's higher, if so assign it to the maximum cost
//             if (J > Jmax) Jmax = J;
//         }

//         // Store the cost
//         c[i] = Jmax;

//         // Update time index
//         timeIDX++;

//         // Trapezoidal integration
//         if (i > 0){
//             I += 0.5 * (c[i - 1] + c[i]) * dt;
//         }
//     }
//     return I;
// }

// /*
//     Dynamic airspace closes point of approach cost
//     integration over time for final turn Dubins
// */
// bool finalTurnDynamicTrafficCost(const Node *node, SearchProblem *problem){
//     if (!node || !problem || !problem->finalTurn || !problem->dynTraffic) {
//         return 0.0;
//     }

//     // Sample the final Dubins turn at dt=1s
//     int N = 0;
//     struct Pos *pts = getDubinsCoordinatesWithFixedTimeStep(problem->finalTurn,
//                                                             NULL,
//                                                             &N,
//                                                             &problem->GeoOpt,
//                                                             NULL,
//                                                             problem);
//     if (!pts || N < 2) {
//         if (pts) free(pts);
//         return 0.0;
//     }

//     // Absolute time index start
//     int t0 = (int) lround(node->t);

//     // Cap integration horizon to available dynamic traffic
//     int maxN = problem->dynTraffic->nsteps - t0;
//     if (maxN <= 1) {
//         free(pts);
//         return 0.0;
//     }
//     if (N > maxN) N = maxN;

//     // Integrate max CPA cost along the Dubins samples
//     const double dt = 1.0;
//     double I = 0.0;

//     // c[0]
//     double Jmax_prev = 0.0; 
//     {
//         int timeIDX = t0; // for i=0
//         struct Pos own = pts[0];

//         for (uint32_t aircraftID = 0; aircraftID < problem->dynTraffic->nflights; aircraftID++) {
//             double J = JCPA(&own, aircraftID, timeIDX, problem);
//             if (J > Jmax_prev) Jmax_prev = J;
//         }
//     }

//     for (int i = 1; i < N; i++) {
//         int timeIDX = t0 + i;
//         struct Pos own = pts[i];
//         double Jmax = 0.0;

//         for (uint32_t aircraftID = 0; aircraftID < problem->dynTraffic->nflights; aircraftID++) {
//             double J = JCPA(&own, aircraftID, timeIDX, problem);
//             if (J > Jmax) Jmax = J;
//         }

//         // Trapezoidal integration in time
//         I += 0.5 * (Jmax_prev + Jmax) * dt;
//         Jmax_prev = Jmax;
//     }
//     free(pts);

//     // Compute max. allowed cost
//     double maxAllowedCPA = 0.5 * (N - 1) * dt; // Average instantaneous risk over the final turn must be ≤ 0.5

//     // printf("maxAllowedCPA: %.4f, I: %.4f\n", maxAllowedCPA, I);

//     return I <= maxAllowedCPA;
// }


/*
    Returns the cost of state
*/
double stateCost(Node *node, SearchProblem * problem){

    // If altitude is lower than the goal altitude, return infinite cost
    if (node->state.alt < problem->Goal.alt) return INFINITY;

    // Remaining optimum-glide range
    double dbg = (node->state.alt - problem->Goal.alt)/tan(problem->gammaOpt*DEG_2_RAD);

    // Distance to the goal state
    double dist;
    double course;
    geo_dist(&node->state, &problem->Goal, &dist, &course, &problem->GeoOpt);
    dist *= NM_2_FT;
    dist -= problem->identRadius;

    // If unreachable or below the goal altitude
    if ((dbg - dist) < 0) return INFINITY;

    // Cost terms
    double *gd  = distanceCost(&node->state, node->l, problem);
    node->gd1 = gd[0];
    node->gd2 = gd[1];
    free(gd); gd = NULL;

    // Course angle cost
    node->gchi = courseCost(&node->state, problem);

    // Overflown population cost
    node->gp = 0;
    if (problem->w_gp > 0 && node->action->length != 0){
        double l = node->l - node->action->length;
        clock_t begin = clock();
        node->gp = populationCost(&node->parent->state, node->action, l, problem);
        clock_t end = clock();
        problem->totalRiskRuntime += (double) (end - begin)*1000 / CLOCKS_PER_SEC;
    }
    
    // Total cost
    if (problem->solver == SEARCH) {
        if (node->gchi >= .5) {
            return 0.2*node->gd1 + 0.3*node->gd2 + 0.2*node->gchi + 0.3*node->gp;
        }
        else if ((node->gchi > 1e-4) && (node->gchi < .5)) {
            return 0.2*node->gd1 + 0.3*node->gd2 + 0.1*node->gchi + 0.4*node->gp;
        }
        else {
            return 0.2*node->gd1 + 0.4*node->gd2 + 0.0*node->gchi + 0.4*node->gp;
        }
    } else if (problem->solver == SEARCH_AIRSPACE){

        // Airspace costs
        node->ga = 0;
        if (node->action->length != 0) {
            double l = node->l - node->action->length;
            clock_t begin = clock();
            node->ga = airspaceCost(&node->parent->state, node->action, l, problem);
            clock_t end = clock();
            problem->totalRiskRuntime += (double) (end - begin)*1000 / CLOCKS_PER_SEC;
        }

        // Cumulative risk
        if (node->parent) {
            node->cumul_gp = node->parent->cumul_gp + node->gp;
            node->cumul_ga = node->parent->cumul_ga + node->ga;
        } else {
            node->cumul_gp = node->gp;
            node->cumul_ga = node->ga;
        }

        // Cumulative path cost g(n)
        node->g = node->cumul_ga;

        // Heuristic cost h(n)
        double ga_heur = maxDensityAhead(&node->state, problem);
        if (problem->w_gp && node->gp) node->h = 0.1*node->gd1 + 0.1*node->gd2 + 0.05*node->gchi + 0.25*ga_heur + 0.5*node->gp;
        else node->h = 0.2*node->gd1 + 0.2*node->gd2 + 0.1*node->gchi + 0.5*ga_heur;

        double w = 1;
        if (problem->initAirspaceCost && (node->ga > 0) && (problem->InitialSearch.alt - node->state.alt < 500)) w = 0;
    
        // Total cost f(n)
        return node->g + w*node->h;

    } else if (problem->solver == SEARCH_AIRSPACE_DYN){

        // Heuristic cost h(n)
        // Airspace corridor occupation costs
        double prohibitedArea_cost  = getAirTrafficDensity(&node->state, problem->prohibited);
        if (prohibitedArea_cost > 0.5) return INFINITY;
        double corridor_cost        = getAirTrafficDensity(&node->state, problem->traffic);
        double heli_cost            = getAirTrafficDensity(&node->state, problem->heli);
        
        // Weighted corridor occupation cost
        double ha = 0.5*corridor_cost + 0.25*heli_cost + 0.25*prohibitedArea_cost;

        // Look ahead heuristc
        double lookahead = maxDensityAhead(&node->state, problem);

        if (problem->w_gp && node->gp) node->h = 0.1*node->gd1 + 0.1*node->gd2 + 0.05*node->gchi + 0.1*lookahead + 0.15*ha + 0.5*node->gp;
        else node->h = 0.2*node->gd1 + 0.2*node->gd2 + 0.1*node->gchi + 0.25*lookahead + 0.25*ha;

        // Cumulative risk
        // Dynamic airspace cost
        node->ga = node->parent ? dynamicAirspaceCost(node->parent, node, problem) : 0.0;
        if (node->parent) {
            node->cumul_gp = node->parent->cumul_gp + node->gp;
            node->cumul_ga = node->parent->cumul_ga + node->ga;
        } else {
            node->cumul_gp = node->gp;
            node->cumul_ga = node->ga;
        }

        // Cumulative path cost g(n)
        node->g = node->cumul_ga;

        // printf("node->g: %.2f, node->h: %.2f\n", node->g, node->h);

        // Total cost f(n)
        // printf("node->g: %.4f, node->h: %.4f, node->f: %.4f\n", node->g, node->h, node->g + node->h);
        double wh = 0.01;
        return node->g + wh * node->h;
    }

    return NAN;
}

bool dubinsFinalTurn(struct Node *node, SearchProblem * problem){
    // Get the candidate solution
    problem->path = solution(node, &node->depth, problem);
    problem->path->depth = node->depth;
    problem->pathActions = solutionActions(node, &node->depth, problem);

    double *alt_corr = (double *) malloc(problem->path->depth * sizeof(double));
    computeCorrectedAltitudes(problem->path,
                              problem->path->depth,
                              problem->ac,
                              alt_corr);

    // Define the initial and final positions for Dubins, and initialize them
    struct Pos init = {.lat = node->state.lat, .lon = node->state.lon, .alt = alt_corr[problem->path->depth-1]};
    struct Pos final = {.lat = problem->Goal.lat, .lon = problem->Goal.lon, .alt = problem->Goal.alt};

    // Define new initial and final positions for iteration
    struct Pos new_init = {.lat = node->state.lat, .lon = node->state.lon, .alt = init.alt};
    struct Pos new_final = {.lat = problem->Goal.lat, .lon = problem->Goal.lon, .alt = problem->Goal.alt};

    // Initial great-circle distance between the initial and final states
    double dist, dist2, course;
    geo_dist(&init, &final, &dist, &course, &problem->GeoOpt);

    // If the distance is less than 3 turn radius, find new initial and final positions
    double step = 10*FT_2_NM;
    double *gammaArray = getOptimalGamma(&node->state.hdg, problem);
    double gamma = gammaArray[2];
    free(gammaArray); gammaArray = NULL;

    gammaArray = getOptimalGamma(&problem->Goal.hdg, problem);
    double gammaFinal = gammaArray[2];
    free(gammaArray); gammaArray = NULL;

    double extendInitDirection = wrapTo360(node->state.hdg - 180);
    double extendFinalDirection = problem->Goal.hdg;
    double minDist = 2*problem->ac->turnRadius * FT_2_NM;
    double deaten = problem->ac->turnRadius * tan(DEG_2_RAD * fabs(node->action->deltaCourse)/2);
    double allowed_backward_traversal = node->action->length - deaten;

    while (dist <= minDist) {

        // Find the distance from initial Dubins position to the given state
        geo_dist(&init, &node->state, &dist2, &course, &problem->GeoOpt);

        // If this distance is less than allowed_backward_traversal, pull it back to increase the distance in-between the initial and final Dubins positions
        if (dist2 <= allowed_backward_traversal*FT_2_NM) {
            geo_npos(&init, &new_init, &step, &extendInitDirection, &problem->GeoOpt);
            new_init.alt += step*NM_2_FT*tan(gamma*DEG_2_RAD);
            init = new_init;
        }

        // Extend the final Dubins position toward the touchdown
        geo_npos(&final, &new_final, &step, &extendFinalDirection, &problem->GeoOpt);
        new_final.alt = problem->Touchdown.alt + (step*NM_2_FT)*tan(gammaFinal*DEG_2_RAD);
        final = new_final;

        // Update the distance in-between the initial and final Dubins positions
        geo_dist(&init, &final, &dist, &course, &problem->GeoOpt);
    }

    // Set Dubins' initial and final positions
    problem->finalTurn->traj[0].wpt.pos = init;
    problem->finalTurn->traj[0].wpt.hdg = node->state.hdg;
    problem->finalTurn->traj[3].wpt.pos = final;
    problem->finalTurn->traj[3].wpt.hdg = problem->Touchdown.hdg;

    // Turn radius
    problem->finalTurn->traj[0].wpt.rad = problem->ac->turnRadius * FT_2_NM;
    problem->finalTurn->traj[1].wpt.rad = 0.0;
    problem->finalTurn->traj[2].wpt.rad = problem->ac->turnRadius * FT_2_NM;

    // Straight flight path angle
    gammaArray = getOptimalGamma(&course, problem);
    problem->finalTurn->traj[1].wpt.gam = -gammaArray[2];
    free(gammaArray); gammaArray = NULL;

    // Turn flight path angle
    problem->finalTurn->traj[0].wpt.gam = -problem->ac->gammaOptTurn;
    problem->finalTurn->traj[2].wpt.gam = -problem->ac->gammaOptTurn;

    // Compute the initial approximate path
    shortestDubins(problem->finalTurn, problem->DubinsOpt);

    // Set the optimal turn flight path angles
    // given the heading changes of the approximate path
    updateOptimalGammaTurn(problem->finalTurn, problem);
    
    // Compute the shortest path again with the updated turn flight path angles
    int type = problem->finalTurn->type;
    shortestDubins(problem->finalTurn, problem->DubinsOpt);
    problem->finalTurn->type = type;

    // Assert final Dubins altitude
    geo_dist(&problem->Touchdown, &problem->finalTurn->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);
    double dh = problem->finalTurn->traj[3].wpt.pos.alt - problem->Touchdown.alt;
    double fingam = atan(dh / (dist*NM_2_FT)) * RAD_2_DEG;

    gammaArray = getOptimalGamma(&problem->Goal.hdg, problem);
    if (fingam < gammaArray[0] || fingam > gammaArray[1]){
        free(gammaArray); gammaArray = NULL;
        return EXIT_FAILURE;
    }

    free(gammaArray); gammaArray = NULL;

    return EXIT_SUCCESS;
}

/*
    Tests the given node for solution
*/
bool goalTest(Node *node, SearchProblem * problem){

    // Check course angle cost - it must be zero
    if (node->gchi > 0) {
        return false;
    }

    // Check great-circle distance to the goal state
    double dist;
    double course;
    geo_dist(&node->state, &problem->Goal, &dist, &course, &problem->GeoOpt);
    dist *= NM_2_FT;

    if ((dist > problem->identRadius) || (dist < problem->ac->turnRadius)) {
        return false;
    }

    // Check altitude
    double *gammaArray = getOptimalGamma(&course, problem);
    double hmax = problem->Goal.alt + problem->identRadius*tan(DEG_2_RAD*gammaArray[1]);

    if (node->state.alt > hmax) {
        free(gammaArray); gammaArray = NULL;
        return false;
    }

    // Compute the final turn
    if (dubinsFinalTurn(node, problem)) return false;

    // Check final turn loss of well-clear
    // if (!finalTurnDynamicTrafficCost(node, problem)) return false;

    // Get the feasible flight path angle range for the final approach
    gammaArray = getOptimalGamma(&problem->Goal.hdg, problem);

    // Check the flight path angle of the remaining traversal (i.e., the final approach)
    geo_dist(&problem->finalTurn->traj[3].wpt.pos, &problem->Touchdown, &dist, &course, &problem->GeoOpt);
    double dh = problem->finalTurn->traj[3].wpt.pos.alt - problem->Touchdown.alt;
    double gamma_rem = atan(dh/(dist*NM_2_FT)) * RAD_2_DEG;

    if ((gamma_rem >= gammaArray[0]) && (gamma_rem <= gammaArray[1])) {
        free(gammaArray); gammaArray = NULL;
        problem->dGamma = 0;
        return true;
    } else {
        free(gammaArray); gammaArray = NULL;
        return false;
    }

    free(gammaArray); gammaArray = NULL;
    problem->dGamma = 0;
    return false;
}

/*
    Modifies the flight path angle of the search path
    Returns 0 if feasible, 1 otherwise.
*/
bool modifyGamma(double remaining_distance, double *gammaArray, Node * node, SearchProblem *problem){

    // Altitude that can be lost in a given straight line distance
    double dh_feasible = remaining_distance*tan(gammaArray[2]*DEG_2_RAD);

    // Optimum altitude of the state
    double hopt = problem->Touchdown.alt + dh_feasible;

    // How much altitude gain is necessary
    double dh = problem->finalTurn->traj[3].wpt.pos.alt - hopt;
    if (dh > 0) return EXIT_FAILURE;

    // How much flight path angle we must change to make current the altitude optimum altitude
    double dGamma = atan(dh/node->l)*RAD_2_DEG;

    double new_gamma;
    Node *current_node = node;

    // Traverse through the parent nodes and check the flight path angle
    while (current_node != NULL && current_node->parent) {
        new_gamma = current_node->action->gamma + dGamma;
        if (new_gamma < gammaArray[0] || new_gamma > gammaArray[1]) {

            return EXIT_FAILURE;
        }
        current_node = current_node->parent;
    }

    problem->dGamma = dGamma;

    return EXIT_SUCCESS;
}


/*
    Finds whether a given state (state1) is
    inside a hexagon with centroid state (state2)
*/
bool isInHexagon(struct Pos *state1, struct Pos *state2, double circumradius, SearchProblem * problem){
    // Get the distance between states
    double dist;
    double course;
    geo_dist(state2, state1, &dist, &course, &problem->GeoOpt);
    dist *= NM_2_FT;

    // If the distance is greater than the circumradius of hexagon, state lies outside of the hexagon.
    if (dist > circumradius) {
        return false;
    }
    // If the distance is less than the radius of the inscribed circle, state lies inside of the hexagon.
    else if (dist <= circumradius*sqrt(3)/2) {
        return true;
    } else { // If neither of above is true, check the relative position of the state with respect to the hexagon
    
        int a = floor(course/60);
        double alpha = course - a*60;
        double x = circumradius*sin(60*DEG_2_RAD) / sin((120 - alpha)*DEG_2_RAD);
        if (dist <= x) {
            return true;
        } else {
            return false;
        }
    }
}


/*
    Checks if a given state creates a cluster
    with any of the states in a given closed list.
*/
bool doesCluster(struct Node * node, SearchProblem *problem){

    // Find explored nodes with similar heading and altitude
    double dPsi, dh;

    for (size_t i = 0; i < problem->explored->size; i++) {

        // Check time first!
        int t_cand = (int) lround(node->t);
        int t_rep  = (int) lround(problem->explored->nodes[i].t);
        if (abs(t_cand - t_rep) > problem->dTime) continue;         // It is a unique state.

        // Check geometry
        dPsi = fabs(wrapTo180(node->state.hdg - problem->explored->nodes[i].state.hdg));
        dh   = fabs(node->state.alt - problem->explored->nodes[i].state.alt);

        if ((dPsi <= problem->dPsi) && (dh <= 0.5*problem->dAltitude)) {

            // Check if it is in cells
            if (problem->explored->nodes[i].action->length > 1e-5) {
                if (isInHexagon(&node->state, &problem->explored->nodes[i].state, 0.5 * problem->explored->nodes[i].action->length, problem)) {
                    return true;
                }
            } else {
                if (isInHexagon(&node->state, &problem->explored->nodes[i].state, 0.5 * problem->lmin, problem)) {
                    return true;
                }
            }
        }
    }

    return false;
}



/*
    Checks if a given state creates cluster
    with any of the states in a given closed list.
*/
bool doesClusterFrontier(struct Node * node, SearchProblem * problem){

    // Find explored nodes with similar heading and altitude
    double dPsi, dh;

    for (size_t i = 0; i < problem->pq->size; i++){

        // Check time first!
        int t_cand = (int) lround(node->t);
        int t_rep  = (int) lround(problem->pq->heap[i].node->t);
        if (abs(t_cand - t_rep) > problem->dTime) continue;

        // Check geometry
        dPsi = fabs(wrapTo180(node->state.hdg - problem->pq->heap[i].node->state.hdg));
        dh   = fabs(node->state.alt - problem->pq->heap[i].node->state.alt);

        if ((dPsi <= problem->dPsi) && (dh <= 0.5*problem->dAltitude)) {

            // Check if it is in cells
            if (problem->pq->heap[i].node->action) {
                if (isInHexagon(&node->state, &problem->pq->heap[i].node->state, 0.5 * problem->pq->heap[i].node->action->length, problem)) {
                    return true;
                }
            } else {
                if (isInHexagon(&node->state, &problem->pq->heap[i].node->state, 0.5 * problem->lmin, problem)) {
                    return true;
                }
            }
        }
    }
    return false;
}

/*
    Returns coordinate samples for 
    risk integration
*/
void searchIntegrationCoordinates(SearchProblem * problem, double interval, struct Pos **coordinates, int *numSamples){
    if (!problem->path) {
        fprintf(stderr, "ERROR: pathOverflownPopulation called with NULL problem->path\n");
        return;
    }

    // Sample integration coordinates
    double *distances = (double *) malloc((problem->path->depth + 1) * sizeof(double));
    if (!distances) {
        fprintf(stderr, "malloc failed for distances - pathOverflownPopulation - search.c\n");
        exit(EXIT_FAILURE);
    }
    cumulativeDistances(problem, distances);

    double *sampled_lat, *sampled_lon, *sampled_alt;
    sampleEquidistantPoints(problem, distances, interval, &sampled_lat, &sampled_lon, &sampled_alt, numSamples);

    // Allocate memory for sampling coordinates
    *coordinates = (struct Pos *) malloc((*numSamples) * sizeof(struct Pos));
    if (!(*coordinates)) {
        fprintf(stderr, "malloc failed for coordinates - searchIntegrationCoordinates\n");
        exit(EXIT_FAILURE);
    }
    
    // Store sampling coordinates
    for (int i = 0; i < (*numSamples); i++){
        (*coordinates)[i].lat = sampled_lat[i];
        (*coordinates)[i].lon = sampled_lon[i];
        (*coordinates)[i].alt = sampled_alt[i];
    }

    // Free memory
    free(distances); distances = NULL;
    free(sampled_lat); sampled_lat = NULL;
    free(sampled_lon); sampled_lon = NULL;
    free(sampled_alt); sampled_alt = NULL;
}

/*
    Returns coordinate samples for risk integration
    sampled at fixed time step (seconds) using groundspeed
*/
/*
    Returns coordinate samples for risk integration
    sampled at fixed time step (seconds) using groundspeed

    Outputs:
      *coordinates: length *numSamples
      *gs:          length *numSamples (knots, unless you want ft/s)
*/
void searchIntegrationCoordinatesWithFixedTimeStep(SearchProblem *problem,
                                                   double dt_sample,
                                                   struct Pos **coordinates,
                                                   double **gs,
                                                   int *numSamples){
    if (!problem || !problem->path || problem->path->depth < 1) {
        fprintf(stderr, "ERROR: searchIntegrationCoordinatesWithFixedTimeStep: invalid path\n");
        *coordinates = NULL;
        if (gs) *gs = NULL;
        *numSamples = 0;
        return;
    }

    int N = problem->path->depth;

    // --- cumulative distance along path (size N+1) ---
    double *distances = (double *)malloc((N + 1) * sizeof(double));
    if (!distances) { perror("malloc distances"); exit(EXIT_FAILURE); }
    cumulativeDistances(problem, distances); // distances in ft (per your comment)

    // --- precompute groundspeed at each node (knots) ---
    double *gs_node = (double *)malloc(N * sizeof(double));
    if (!gs_node) { perror("malloc gs_node"); exit(EXIT_FAILURE); }

    for (int i = 0; i < N; i++) {
        gs_node[i] = groundspeed(&problem->path->nodes[i].state, problem); // knots
        if (gs_node[i] < 0.0) gs_node[i] = 0.0;
    }

    // --- build cumulative time along path (size N) ---
    double *cumTime = (double *)malloc(N * sizeof(double));
    if (!cumTime) { perror("malloc cumTime"); exit(EXIT_FAILURE); }
    cumTime[0] = 0.0;

    for (int i = 0; i < N - 1; i++) {
        double ds = distances[i + 1] - distances[i]; // ft
        if (ds < 0) ds = 0;

        // use average GS over the segment for time accumulation (more consistent)
        double gs_avg_kt = 0.5 * (gs_node[i] + gs_node[i + 1]);
        double v_fts = gs_avg_kt * KTS_2_FTS; // ft/s

        if (v_fts <= 1e-6) v_fts = 1e-6;
        double dt_seg = ds / v_fts;
        cumTime[i + 1] = cumTime[i] + dt_seg;
    }

    double T_end = cumTime[N - 1];
    if (dt_sample <= 0) dt_sample = 1.0;

    // Include both endpoints: t = 0, dt, ..., T_end
    int M = (int) floor(T_end / dt_sample) + 1;

    // If last sample isn't exactly T_end, add one more sample at T_end
    if (fabs((M - 1) * dt_sample - T_end) > 1e-9) M += 1;

    *coordinates = (struct Pos *)malloc(M * sizeof(struct Pos));
    if (!(*coordinates)) { perror("malloc coordinates"); exit(EXIT_FAILURE); }

    if (gs) {
        *gs = (double *)malloc(M * sizeof(double));
        if (!(*gs)) { perror("malloc gs"); exit(EXIT_FAILURE); }
    }

    *numSamples = M;

    // --- sample by time and interpolate within segments ---
    int seg = 0;

    for (int k = 0; k < M; k++) {
        double tq = k * dt_sample;
        if (k == M - 1) tq = T_end;

        while (seg < N - 2 && cumTime[seg + 1] < tq) seg++;

        double t0 = cumTime[seg];
        double t1 = cumTime[seg + 1];
        double denom = (t1 - t0);

        double alpha = 0.0;
        if (denom > 1e-12) {
            alpha = (tq - t0) / denom;
            if (alpha < 0.0) alpha = 0.0;
            if (alpha > 1.0) alpha = 1.0;
        }

        struct Pos *p0 = &problem->path->nodes[seg].state;
        struct Pos *p1 = &problem->path->nodes[seg + 1].state;

        (*coordinates)[k].lat = p0->lat + alpha * (p1->lat - p0->lat);
        (*coordinates)[k].lon = p0->lon + alpha * (p1->lon - p0->lon);
        (*coordinates)[k].alt = p0->alt + alpha * (p1->alt - p0->alt);
        (*coordinates)[k].hdg = p0->hdg + alpha * (p1->hdg - p0->hdg);

        // Return GS at sample time (linear between nodes; units: knots)
        if (gs && *gs) {
            (*gs)[k] = gs_node[seg] + alpha * (gs_node[seg + 1] - gs_node[seg]);
        }
    }

    free(distances); distances = NULL;
    free(cumTime); cumTime = NULL;
    free(gs_node); gs_node = NULL;
}



// Recompute physically-correct altitudes for a finalized path,
// charging entry half-turn to the end of the previous leg,
// exit half-turn at the start of the current leg,
// and straight in between. Does NOT modify node->state.alt.
// Writes altitudes into out_alt[i] for node i (0..depth-1).
void computeCorrectedAltitudes(const struct Path *path,
                                int depth,
                                const Aircraft *ac,
                                double *out_alt) {

    const double R = ac->turnRadius;         // [ft]
    double d_entry, s_entry, d_exit, s_exit, str_l_eff;

    out_alt[0] = path->nodes[0].state.alt;

    for (int i = 0; i < depth - 1; ++i) {

        if (i == 0){

            // Segment initial exit turn
            double half_chi = fabs(path->nodes[1].action->deltaCourse)/2 * DEG_2_RAD; // Half turn angle
            d_exit = R * tan(half_chi); // Entry eaten distance
            s_exit  = R * half_chi;      // Entry arc length
            str_l_eff = path->nodes[1].action->length - d_exit;

            double gamma_str = path->nodes[i].action->gamma * DEG_2_RAD;
            double gamma_turn =  path->nodes[i].action->gamma_turn* DEG_2_RAD;

            out_alt[i+1] = out_alt[i] - s_exit*tan(gamma_turn) - str_l_eff*tan(gamma_str);

        } else {

            // Previous segment secondary entry turn
            double half_chi = fabs(path->nodes[i+1].action->deltaCourse)/2 * DEG_2_RAD; // Half turn angle
            d_entry = R * tan(half_chi); // Entry eaten distance
            d_exit = d_entry;
            s_entry  = R * half_chi;      // Entry arc length
            s_exit = s_entry;
            
            double gamma_str = path->nodes[i].action->gamma * DEG_2_RAD;
            double gamma_turn = path->nodes[i+1].action->gamma_turn * DEG_2_RAD;
            
            // Make correction on the altitude of the current state
            out_alt[i] = out_alt[i] + d_entry*tan(gamma_str) - s_entry*tan(gamma_turn);

            // Now fly the next segment
            gamma_str = path->nodes[i+1].action->gamma * DEG_2_RAD;

            // Effectice straight length
            str_l_eff = path->nodes[i+1].action->length - d_exit;

            // Propagate the altitude till the next waypoint
            out_alt[i+1] = out_alt[i] - s_exit*tan(gamma_turn) - str_l_eff*tan(gamma_str);
        }
    }
}

// /*
//     Writes search path waypoints, explored states,
//     and Dubins waypoints into csv files.
// */
// int writeResults(SearchProblem *problem, struct DubinsPath *bestDubins, char *folderName)
// {   

//     // Write to a file
//     FILE * file;
//     char directory[200];

//     if (problem->exitFlag == 0 || problem->exitFlag == 1){
        
//         sprintf(directory, "%s/path.csv", folderName);

//         // PATH WAYPOINTS
//         file = fopen(directory, "w");
//         if (!file)
//         {
//             perror("Can't open the path csv file...");
//             return EXIT_FAILURE;
//         }

//         // Write path waypoints
//         // Allocate corrected altitude buffer
//         double *alt_corr = (double *) malloc(problem->path->depth * sizeof(double));
//         computeCorrectedAltitudes(
//             problem->path,
//             problem->path->depth,
//             problem->ac,
//             alt_corr
//         );

//         // When writing path.csv, use alt_corr[i] instead of problem->path->nodes[i].state.alt
//         for (size_t i = 0; i < problem->path->depth; i++) {
//             fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f\n",
//                 problem->path->nodes[i].state.lat,
//                 problem->path->nodes[i].state.lon,
//                 alt_corr[i],
//                 problem->path->nodes[i].state.hdg,
//                 problem->path->nodes[i].t);
//         }
//         free(alt_corr); alt_corr = NULL;
//         fclose(file);

//         // FINAL TURN DUBINS
//         sprintf(directory, "%s/finalturn_sim.csv", folderName);
//         file = fopen(directory, "w");
//         if (!file) {
//             perror("Can't open the final turn csv file...");
//             return EXIT_FAILURE;
//         }
//         for (int i = 0; i < 4; i++){
//             fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %d\n", problem->finalTurn->traj[i].wpt.pos.lat,
//                                                                 problem->finalTurn->traj[i].wpt.pos.lon,
//                                                                 problem->finalTurn->traj[i].wpt.pos.alt,
//                                                                 problem->finalTurn->traj[i].wpt.hdg,
//                                                                 problem->finalTurn->traj[i].wpt.gam,
//                                                                 problem->finalTurn->type);
//         }
//         fprintf(file,"%.6f, %.6f\n", problem->finalTurn->orbit1.lat, problem->finalTurn->orbit1.lon);
//         fprintf(file,"%.6f, %.6f\n", problem->finalTurn->orbit2.lat, problem->finalTurn->orbit2.lon);
//         fclose(file);

//         sprintf(directory, "%s/finalturn.csv", folderName);
//         file = fopen(directory, "w");
//         if (!file) {
//             perror("Can't open the final turn csv file...");
//             return EXIT_FAILURE;
//         }
//         int sampleSize;
//         // struct Pos *coordinates = getDubinsCoordinates(problem->finalTurn, interval, &sampleSize, &problem->GeoOpt, directory);
//         struct Pos *coordinates = getDubinsCoordinatesWithFixedTimeStep(problem->finalTurn, NULL, &sampleSize, &problem->GeoOpt, directory, problem);
//         free(coordinates);
//         fclose(file);

//         // PATH ACTIONS
//         sprintf(directory, "%s/actions.csv", folderName);
//         file = fopen(directory, "w");
//         if (!file) {
//             perror("Can't open the actions csv file...");
//             return EXIT_FAILURE;
//         }

//         // Write
//         for (size_t i = 0; i < problem->path->depth; i++) {
//             fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", problem->pathActions[i].deltaCourse,
//                                                 problem->pathActions[i].length,
//                                                 problem->pathActions[i].gamma,
//                                                 problem->pathActions[i].gamma_turn);
//         }
//         fclose(file);
//     }

//     // Write best/fallback Dubins coordinates
//     if (problem->exitFlag >= 0) {

//         // Open file
//         sprintf(directory, "%s/best_dubins.csv", folderName);
//         file = fopen(directory, "w");

//         // S-Turn
//         if (bestDubins->size == 2) {
//             // Get coordinates
//             int sampleSize1, sampleSize2;
//             struct Pos *coordinates1 = getDubinsCoordinatesWithFixedTimeStep(&bestDubins[0], NULL, &sampleSize1, &problem->GeoOpt, NULL, problem);
            
//             // Write
//             for (int i = 0; i < sampleSize1; i++) {
//                 fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f\n", coordinates1[i].lat, coordinates1[i].lon, coordinates1[i].alt, coordinates1[i].hdg, coordinates1[i].t);
//             }
            
//             // Get coordinates
//             struct Pos *coordinates2 = getDubinsCoordinatesWithFixedTimeStep(&bestDubins[1], &coordinates1[sampleSize1-1], &sampleSize2, &problem->GeoOpt, NULL, problem);

//             // Write
//             for (int i = 1; i < sampleSize2; i++){
//                 fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f\n", coordinates2[i].lat, coordinates2[i].lon, coordinates2[i].alt, coordinates2[i].hdg, coordinates1[sampleSize1-1].t + coordinates2[i].t);
//             }
            
//             free(coordinates1); coordinates1 = NULL;
//             free(coordinates2); coordinates2 = NULL;

//         } else {

//             int sampleSize;
//             struct Pos *coordinates1 = getDubinsCoordinatesWithFixedTimeStep(bestDubins, NULL, &sampleSize, &problem->GeoOpt, NULL, problem);

//             // Write
//             for (int i = 0; i < sampleSize; i++) {
//                 fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f\n", coordinates1[i].lat,  coordinates1[i].lon,  coordinates1[i].alt,  coordinates1[i].hdg, coordinates1[i].t);
//             }
//             free(coordinates1); coordinates1 = NULL;
//         }
//         fclose(file);
//     }

//     // Write best/fallback Dubins for simulation
//     if (problem->exitFlag >= 0) {

//         // Open file
//         sprintf(directory, "%s/best_dubins_sim.csv", folderName);
//         file = fopen(directory, "w");

//         sprintf(directory, "%s/best_dubins_sim_orbit.csv", folderName);
//         FILE * file2 = fopen(directory, "w");

//         // S-Turn
//         if (bestDubins->size == 2) {
            
//             // Write waypoint coordinates
//             for (int j = 0; j < 2; j++){
//                 for (int i = 0; i < 4; i++){
//                     fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %d\n", bestDubins[j].traj[i].wpt.pos.lat,
//                                                                         bestDubins[j].traj[i].wpt.pos.lon,
//                                                                         bestDubins[j].traj[i].wpt.pos.alt,
//                                                                         bestDubins[j].traj[i].wpt.hdg,
//                                                                         bestDubins[j].traj[i].wpt.gam,
//                                                                         bestDubins[j].type);
//                     }

//                     // Write orbit coordinates
//                     fprintf(file2,"%.6f, %.6f\n", bestDubins[j].orbit1.lat, bestDubins[j].orbit1.lon);
//                     fprintf(file2,"%.6f, %.6f\n", bestDubins[j].orbit2.lat, bestDubins[j].orbit2.lon);
//             }
//             fclose(file);
//             fclose(file2);
            
//         } else {

//             // Write waypoint coordinates
//             for (int i = 0; i < 4; i++){
//                 fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %d\n", bestDubins->traj[i].wpt.pos.lat,
//                                                                     bestDubins->traj[i].wpt.pos.lon,
//                                                                     bestDubins->traj[i].wpt.pos.alt,
//                                                                     bestDubins->traj[i].wpt.hdg,
//                                                                     bestDubins->traj[i].wpt.gam,
//                                                                     bestDubins->type);
//                 }

//             // Write orbit coordinates
//             fprintf(file2,"%.6f, %.6f\n", bestDubins->orbit1.lat, bestDubins->orbit1.lon);
//             fprintf(file2,"%.6f, %.6f\n", bestDubins->orbit2.lat, bestDubins->orbit2.lon);
            
//             fclose(file);
//             fclose(file2);
//         }
        
//     }
    
    
//     // EXPLORED SET
//     if (problem->exitFlag >= 0) {
//         sprintf(directory, "%s/explored.csv", folderName);
//         file = fopen(directory, "w");
//         if (!file) {
//             perror("Can't open the explored csv file...");
//             return EXIT_FAILURE;
//         }

//         // Write
//         for (size_t i = 0; i < problem->explored->size; i++) {
//             fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", problem->explored->nodes[i].state.lat,
//                                                     problem->explored->nodes[i].state.lon,
//                                                     problem->explored->nodes[i].state.alt,
//                                                     problem->explored->nodes[i].state.hdg,
//                                                     problem->explored->nodes[i].g,
//                                                     problem->explored->nodes[i].h,
//                                                     problem->explored->nodes[i].f);
//         }
//         fclose(file);

//         // EXPLORED SET ACTIONS
//         sprintf(directory, "%s/explored_actions.csv", folderName);
//         file = fopen(directory, "w");
//         if (!file)
//         {
//             perror("Can't open the explored actions csv file...");
//             return EXIT_FAILURE;
//         }

//         // Write
//         for (size_t i = 0; i < problem->explored->size; i++) {
//             fprintf(file, "%.6f, %.6f, %.6f\n", problem->explored->nodes[i].action->deltaCourse,
//                                                     problem->explored->nodes[i].action->length,
//                                                     problem->explored->nodes[i].action->gamma);
//         }
//         fclose(file);
//     }

//     return EXIT_SUCCESS;
// }
