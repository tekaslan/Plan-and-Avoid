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
#    Dubins Path Solver                                                             %
#    Airspace, Air Traffic, and Ground-risk Aware                                   %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : Pedro Di Donato & H. Emre Tekaslan (tekaslan@vt.edu)                 %
#    Date    : July 2026                                                            %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#                      https://scholar.google.com/citations?user=UCxHXTgAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "../include/dubins.h"

static bool dubins_optimal(struct DubinsPath * const dub,
                          struct DubinsOpt * opt);

int Dubins(struct DubinsPath * const dubins,
           struct DubinsOpt * const opt)

{ 
    int i = 0;          /* Counter */
    struct Pos cnt1 = {.lat =0.0, .lon = 0.0, .alt = 0.0}; /* Center Circle 1 */
    struct Pos cnt2 = {.lat =0.0, .lon = 0.0, .alt = 0.0}; /* Center Circle 2 */
    struct Traj* p = NULL;         /* Auxiliar pointer for traj */
    double hdgcc = 0.0;     /* HDG from center 1 to center 2 */
    double distcc = 0.0;    /* Distance form center 1 to center 2 */
    double alpha= 0.0;      /* Auxiliar angle */
    double hdgaux = 0.0;    /* Auxiliar angle */
    double distaux = 0.0;   /* Auxiliar distance */
    double rad[3] = {0.0,0.0,0.0};
    
    /* Initialize Trajectory Array */
    for(i = 0; i < 3; i++)
        rad[i] = dubins->traj[i].wpt.rad;

    /* Optimum Trajectory */
    if(dubins->type == OPT){       
        if(dubins_optimal(dubins, opt))    return EXIT_FAILURE;
        else                               return EXIT_SUCCESS;
    }
        
    /* Initial Position */
    p = dubins->traj;

    /* Find both circle centers */
    if((dubins->type < 3) || (dubins->type == 5) || (dubins->type == 6)){
        hdgaux = (dubins->traj[0].wpt.hdg + 90.0);
        p->wpt.rad = rad[0];
    }
    else{
        hdgaux = (dubins->traj[0].wpt.hdg - 90.0);
        p->wpt.rad = -rad[0];
    }
    geo_npos(&(dubins->traj[0].wpt.pos),&cnt1,&(rad[0]),&hdgaux,&(opt->trajopt.geoopt));
    dubins->orbit1 = cnt1;

    if((dubins->type == 1) || (dubins->type == 3) || 
       (dubins->type == 5) || (dubins->type == 6)){
        hdgaux = (dubins->traj[3].wpt.hdg + 90.0);
        }
    else{
        hdgaux = (dubins->traj[3].wpt.hdg - 90.0);
        }
    geo_npos(&(dubins->traj[3].wpt.pos),&cnt2,&(rad[2]),&hdgaux,&(opt->trajopt.geoopt));

    /* Find Distance and heading between circle centers */
    geo_dist(&cnt1,&cnt2,&distcc,&hdgcc,&(opt->trajopt.geoopt));
    dubins->orbit2 = cnt2;

    /* For Turn-Turn-Turn Case - Find if the solution is possible */
    if((dubins->type >= 5) && (dubins->type <= 8) &&
         ((rad[0] + 2*(rad[1]) + rad[2] < distcc) || 
            (rad[0] + distcc < rad[2]) || 
            (rad[2] + distcc < rad[0]))){
        p->next->wpt.pos.lat = -91;
        return EXIT_FAILURE; 
    }

    /* Find alpha */
    switch (dubins->type){
        case 1:
            alpha = hdgcc - RAD_2_DEG*acos((rad[0]-rad[2])/distcc); break;
        case 2:
            alpha = hdgcc - RAD_2_DEG*acos((rad[0]+rad[2])/distcc); break;
        case 3:
            alpha = hdgcc + RAD_2_DEG*acos((rad[0]+rad[2])/distcc); break;
        case 4:
            alpha = hdgcc + RAD_2_DEG*acos((rad[0]-rad[2])/distcc); break;
        case 5:
            alpha = hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        case 6:
            alpha = hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        case 7:
            alpha = hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        case 8:
            alpha = hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[2])*(rad[1] + rad[2]) +
                            (rad[0] + rad[1])*(rad[0] + rad[1]) + 
                            (distcc)*(distcc))/(2*(distcc)*(rad[0] + rad[1]))); break;      
        default: break;
    }
    
    /* Checking if alpha in NaN (no solution possible) */
    if(alpha != alpha){
        p->next->wpt.pos.lat = -91;
        return EXIT_FAILURE;
    }

    /* Waypoint 2 */
    p = p->next;
    geo_npos(&cnt1, &p->wpt.pos, &(rad[0]), &alpha, &(opt->trajopt.geoopt));
    if((dubins->type < 3) || (dubins->type == 5) || (dubins->type == 6))
        p->wpt.hdg = fmod(alpha + 90.0,360.0);  
    else    
        p->wpt.hdg = fmod(alpha + 270.0,360.0); 

    /* Heading tolerance - Avoid numerical errors */
    if(fabs(p->wpt.hdg-dubins->traj[0].wpt.hdg) < 1e-9){
        p->wpt.hdg = dubins->traj[0].wpt.hdg;
    }

    if(fabs(p->wpt.hdg-dubins->traj[3].wpt.hdg) < 1e-9){
        p->wpt.hdg = dubins->traj[3].wpt.hdg;
    } 

    /* Radius = 0 for CSC Case (only working for now) */
    if((dubins->type < 5)||(p->wpt.hdg == dubins->traj[0].wpt.hdg))    
        p->wpt.rad = 0;
    else{ 
        if((dubins->type == 5) || (dubins->type == 6))
            p->wpt.rad = -rad[1];
        else
            p->wpt.rad = rad[1];
    }

    /* Waypoint 3 */
    switch (dubins->type){
        case 1:
            distaux = distcc*sin(DEG_2_RAD*(hdgcc - alpha));
            (p->next)->wpt.rad = rad[2];
            break;
        case 2:
            distaux = distcc*sin(DEG_2_RAD*(hdgcc - alpha));
            (p->next)->wpt.rad = -rad[2];
            break;
        case 3:
            distaux = distcc*sin(DEG_2_RAD*(alpha - hdgcc));
            (p->next)->wpt.rad = rad[2];
            break;
        case 4:
            distaux = distcc*sin(DEG_2_RAD*(alpha - hdgcc));
            (p->next)->wpt.rad = -rad[2];
            break;
        case 5:
            alpha = 180.0 + hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            
            (p->next)->wpt.rad = rad[2];
            break;   
        case 6:
            alpha = 180.0 + hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            
            (p->next)->wpt.rad = rad[2];
            break;   
        case 7:
            alpha = 180.0 + hdgcc + RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            (p->next)->wpt.rad = -rad[2];
            break;  
        case 8:
            alpha = 180.0 + hdgcc - RAD_2_DEG*acos((-(rad[1] + rad[0])*(rad[1] + rad[0]) +
                                     (rad[2] + rad[1])*(rad[2] + rad[1]) +
                                     (distcc)*(distcc))/(2*(distcc)*(rad[1] + rad[2]))); 
            (p->next)->wpt.rad = -rad[2];
            break;  
        default: break; 
    }
    
    if(dubins->type < 5){
        geo_npos(&p->wpt.pos,&((p->next)->wpt.pos),&distaux,&p->wpt.hdg,&(opt->trajopt.geoopt));
        (p->next)->wpt.hdg = p->wpt.hdg;
    }
    else{
        geo_npos(&cnt2,&((p->next)->wpt.pos),&(rad[2]),&alpha,&(opt->trajopt.geoopt));
        if((dubins->type == 5) || (dubins->type == 6))
            p->next->wpt.hdg = fmod(alpha + 90.0,360.0);
        else
            p->next->wpt.hdg = fmod(alpha + 270.0,360.0);
        if(fabs(p->next->wpt.hdg - p->wpt.hdg) < 1e-9){
            p->next->wpt.hdg = p->wpt.hdg;
        }
    }

    p = p->next;

    /* Waypoint 4 */
    p = p->next;
    p->wpt.rad = 0;

    /* Get Final Distance */
    dubins->hdist = Traj_HDist(dubins->traj,0,&(opt->trajopt));

    return EXIT_SUCCESS;
}

/*
    Returns the shortest Turn-Straight-Turn Dubins path
*/
int shortestDubins(struct DubinsPath * dubins, struct DubinsOpt * const opt)
{

    // Allocate memory for tmp_path inside the loop
    struct DubinsPath *tmp_path = (struct DubinsPath *) malloc(4*sizeof(struct DubinsPath));
    double bestLength = 9999;
    int bestIndex = -1;
    for (int type = 0; type < 4; type++) {

        Traj_InitArray(tmp_path[type].traj, 4);
        Traj_CopyAll(dubins->traj, tmp_path[type].traj);

        tmp_path[type].type = type + 1;
        tmp_path[type].hdist = 0;
        
        // Perform calculations
        Dubins(&tmp_path[type], opt);
        traj_calctraj_angdist(&(tmp_path[type].traj[0]), 0, &(opt->trajopt));
        Traj_Calc3D(&(tmp_path[type].traj[0]), 0, &(opt->trajopt));

        // Check if the path is found
        if (tmp_path[type].traj[1].wpt.pos.lat < -90) {
            continue;
        }

        // Total horizontal length
        tmp_path[type].hdist = 0;
        for (int i = 0; i < 3; i++) {
            tmp_path[type].hdist += tmp_path[type].traj[i].hdist;
        }

        // Compare
        if (tmp_path[type].hdist < bestLength) {
            bestIndex = type;
            bestLength = tmp_path[type].hdist;
        }
    }

    if (bestIndex == -1) {
        free(tmp_path); tmp_path = NULL;
        return EXIT_FAILURE;
    }

    Traj_CopyAll(tmp_path[bestIndex].traj, dubins->traj);
    dubins->hdist = tmp_path[bestIndex].hdist;
    dubins->orbit1 = tmp_path[bestIndex].orbit1;
    dubins->orbit2 = tmp_path[bestIndex].orbit2;
    dubins->type = tmp_path[bestIndex].type;

    free(tmp_path); tmp_path = NULL;

    return EXIT_SUCCESS;
}

/*
    Computes a new intermediate waypoint for S-Turn paths,
    referencing an initial waypoint
*/
void getIntermediateWaypoint(struct Pos *interWaypoint0,
                            struct Pos *interWaypoint,
                            double straightLength,
                            double theta,
                            int extendTo,
                            struct GeoOpt *GeoOpt)
{
    // Compute the vertical extension distance
    double x = 0.5*straightLength * sin((90-0.5*theta)*DEG_2_RAD)/sin(0.5*theta*DEG_2_RAD);

    // Extension direction
    double course = interWaypoint0->hdg + 90*extendTo;
    course = wrapTo360(course);

    // New intermediate waypoint coordinates
    geo_npos(interWaypoint0, interWaypoint, &x, &course, GeoOpt);

    // New intermediate waypoint hdg
    interWaypoint->hdg = interWaypoint0->hdg;
}

/*
    Returns an S-Turn Dubins path
*/
struct DubinsPath *computeSturnDubins(struct DubinsPath *dubins,
                                    double dAltitude,
                                    int extendTo,
                                    SearchProblem * problem)
{

    // Set the initial intermediate waypoint along the straight segment
    struct Pos interWaypoint0, interWaypoint;
    interWaypoint0.lat = 0.5*(dubins->traj[1].wpt.pos.lat + dubins->traj[2].wpt.pos.lat);
    interWaypoint0.lon = 0.5*(dubins->traj[1].wpt.pos.lon + dubins->traj[2].wpt.pos.lon);
    interWaypoint0.alt = 0.5*(dubins->traj[1].wpt.pos.alt + dubins->traj[2].wpt.pos.alt);
    interWaypoint0.hdg = dubins->traj[1].wpt.hdg;

    // Get the initial straight segment length
    double straightLength, course;
    geo_dist(&dubins->traj[1].wpt.pos, &dubins->traj[2].wpt.pos, &straightLength, &course, &problem->GeoOpt);

    // Initialize the actual altitude loss
    double dh = dubins->traj[0].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Create two Dubins path structures to store paths
    struct DubinsPath *sTurnPath = (struct DubinsPath *) malloc(2*sizeof(struct DubinsPath));
    for (int i = 0; i < 2; i++) {
        Traj_InitArray(sTurnPath[i].traj, 4);   // Allocates memory for trajectory structure
        sTurnPath[i].size = 2;                  // Sets path structure size to 2, indicating the structure holds two Dubins paths (S-Turn)
    }

    // Iterate over theta
    double theta = 60; // [deg]
    int maxIter = 1000;
    int counter = 0;
    double step = 1e-3;
    double ftol = 3;
    while ((fabs(dh - dAltitude) > ftol) && (theta < 180) && (counter < maxIter)) {

        // Copy the original initial and final waypoints to the Sturn path structures
        Traj_CopyAll(dubins->traj, sTurnPath[0].traj);
        Traj_CopyAll(dubins->traj, sTurnPath[1].traj);
        
        // Get an intermediate waypoint
        getIntermediateWaypoint(&interWaypoint0, &interWaypoint, straightLength, theta, extendTo, &problem->GeoOpt);

        // Update the Sturn path structures
        sTurnPath[0].traj[3].wpt.pos = interWaypoint;
        sTurnPath[0].traj[3].wpt.hdg = interWaypoint.hdg;

        // Initial estimate
        shortestDubins(&sTurnPath[0], problem->DubinsOpt);

        // Update the turn flight path angles given heading changes and wind conditions
        updateOptimalGammaTurn(&sTurnPath[0], problem);
        shortestDubins(&sTurnPath[0], problem->DubinsOpt);

        sTurnPath[1].traj[0].wpt.pos = interWaypoint;
        sTurnPath[1].traj[0].wpt.hdg = interWaypoint.hdg;
        sTurnPath[1].traj[0].wpt.pos.alt = sTurnPath[0].traj[3].wpt.pos.alt;
        shortestDubins(&sTurnPath[1], problem->DubinsOpt);
        updateOptimalGammaTurn(&sTurnPath[1], problem);
        shortestDubins(&sTurnPath[1], problem->DubinsOpt);

        // Check altitude loss
        dh = sTurnPath[0].traj[0].wpt.pos.alt - sTurnPath[1].traj[3].wpt.pos.alt;

        // Update theta and iteration count
        theta -= (dAltitude - dh)*step;
        counter++;

        // Console output
        // printf("Altitude error: %.2f ft, New angle: %.1f deg\n", (dAltitude - dh), theta);
    }

    return sTurnPath;
}

void updateOptimalGammaTurn(struct DubinsPath * dubins, SearchProblem * problem){
    
    // First orbit optimal flight path angle
    double dPsi1 = dubins->traj[0].dpsi;
    Node tmp_node;
    tmp_node.state.hdg = dubins->traj[0].wpt.hdg;
    tmp_node.action = calloc(1, sizeof(Action));
    if (fabs(dPsi1) <= 90 + 1e-1){
        tmp_node.action->deltaCourse = dPsi1;
        dubins->traj[0].wpt.gam = -getOptimalGammaTurn(&tmp_node, problem);
    } else {
        int nturns = (int) fabs(dPsi1)/90;
        double rem = fmod(dPsi1, 90.0);
        double weightedSum = 0;
        int sign = dPsi1/fabs(dPsi1);
        for (int i = 0; i < nturns; i++){
            tmp_node.action->deltaCourse = 90*sign;
            double gamma = -getOptimalGammaTurn(&tmp_node, problem);
            weightedSum += gamma*90;
            tmp_node.state.hdg = wrapTo360(tmp_node.state.hdg + 90*sign);
        }
        tmp_node.action->deltaCourse = rem;
        double gamma = -getOptimalGammaTurn(&tmp_node, problem);
        weightedSum += gamma*fabs(rem);
        dubins->traj[0].wpt.gam = weightedSum/fabs(dPsi1);
    }

    // Second orbit optimal flight path angle
    double dPsi2 = dubins->traj[2].dpsi;
    tmp_node.state.hdg = dubins->traj[1].wpt.hdg;
    if (fabs(dPsi2) <= 90 + 1e-1){
        tmp_node.action->deltaCourse = dPsi2;
        dubins->traj[2].wpt.gam = -getOptimalGammaTurn(&tmp_node, problem);
    } else {
        int nturns = (int) fabs(dPsi2)/90;
        double rem = fmod(dPsi2, 90.0);
        double weightedSum = 0;
        int sign = dPsi2/fabs(dPsi2);
        for (int i = 0; i < nturns; i++){
            tmp_node.action->deltaCourse = 90*sign;
            double gamma = -getOptimalGammaTurn(&tmp_node, problem);
            weightedSum += gamma*90;
            tmp_node.state.hdg = wrapTo360(tmp_node.state.hdg + 90*sign);
        }
        tmp_node.action->deltaCourse = rem;
        double gamma = -getOptimalGammaTurn(&tmp_node, problem);
        weightedSum += gamma*fabs(rem);
        dubins->traj[2].wpt.gam = weightedSum/fabs(dPsi2);
    }
    free(tmp_node.action); tmp_node.action=NULL;
}
 
/*
   Return Dubins path coordinates and
   writes them to a CSV file.
*/
struct Pos *getDubinsCoordinates(struct DubinsPath *dubins,
                                 double interval, 
                                 int *sampleSize,
                                 struct GeoOpt *GeoOpt,
                                 char * folderName)
    {

    // Get number of samples for each segment
    // First orbit
    double orbitLengthO1 = fabs(dubins->traj[0].dpsi*DEG_2_RAD) * fabs(dubins->traj[0].wpt.rad) * NM_2_FT;
    int nSampleO1 = max(2, orbitLengthO1/interval);

    // Straight
    double dist, bearing;
    geo_dist(&dubins->traj[1].wpt.pos, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    int nSampleS = max(2, (dist*NM_2_FT)/interval);

    // Second orbit
    double orbitLengthO2 = fabs(dubins->traj[2].dpsi*DEG_2_RAD) * fabs(dubins->traj[2].wpt.rad) * NM_2_FT;
    int nSampleO2 = max(2, orbitLengthO2/interval);

    // Total number of samples
    int ntot = (nSampleO1 + nSampleS + nSampleO2) - 2;
    *sampleSize = ntot;
    struct Pos *coordinates = (struct Pos *) malloc(ntot*sizeof(struct Pos));
    
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit1, &dubins->traj[0].wpt.pos, &dist, &bearing, GeoOpt);
    
    double dPsi;
    double hdg;
    double newhdg = dubins->traj[0].wpt.hdg;
    double R = fabs(dubins->traj[0].wpt.rad);
    double alt;
    double dh = dubins->traj[0].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    for (int i = 0; i < nSampleO1; i++) {
        dPsi = i * dubins->traj[0].dpsi/(nSampleO1 - 1);
        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[0].wpt.hdg + dPsi);
        alt = dubins->traj[0].wpt.pos.alt - i*dh/(nSampleO1 - 1);

        geo_npos(&dubins->orbit1, &coordinates[i], &R, &hdg, GeoOpt);
        coordinates[i].alt = alt;
        coordinates[i].hdg = newhdg;
    }

    // STRAIGHT SEGMENT
    int idx;
    double dlat = dubins->traj[2].wpt.pos.lat - dubins->traj[1].wpt.pos.lat;
    double dlon = dubins->traj[2].wpt.pos.lon - dubins->traj[1].wpt.pos.lon;
    double dalt = dubins->traj[2].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    for (int i = 1; i < nSampleS - 1; i++)
    {
        idx = i + nSampleO1 - 1;

        coordinates[idx].lat = dubins->traj[1].wpt.pos.lat + i*dlat/(nSampleS - 1);
        coordinates[idx].lon = dubins->traj[1].wpt.pos.lon + i*dlon/(nSampleS - 1);
        coordinates[idx].alt = dubins->traj[1].wpt.pos.alt + i*dalt/(nSampleS - 1);
        coordinates[idx].hdg = newhdg;
    }

    // SECOND ORBIT
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit2, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    
    dPsi = 0;
    hdg = 0;
    newhdg = 0;
    dh = dubins->traj[2].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;
    for (int i = 0; i < nSampleO2; i++) {
        idx = i + nSampleO1 + nSampleS - 2;
        dPsi = i * dubins->traj[2].dpsi/(nSampleO2 - 1);
        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[2].wpt.hdg + dPsi);
        geo_npos(&dubins->orbit2, &coordinates[idx], &R, &hdg, GeoOpt);
        coordinates[idx].alt = dubins->traj[2].wpt.pos.alt - i*dh/(nSampleO2 - 1);
        coordinates[idx].hdg = newhdg;
    }

    if (folderName) {
        // Create a file
        FILE * file;

        // Open the file
        file = fopen(folderName, "w");
        if (!file)
        {
            perror("Can't open the dubins csv file...");
            return coordinates;
        }

        // Write coordinates
        for (int i = 0; i < ntot; i++)
        {
            fprintf(file, "%.6f, %.6f, %.6f, %.6f\n", coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, coordinates[i].hdg);
        }
        fclose(file);
    }

    return coordinates;
}

/*
   Return Dubins path coordinates and
   writes them to a CSV file.
*/
struct Pos *getDubinsCoordinatesWithFixedTimeStep_(struct DubinsPath *dubins,
                                                int *sampleSize,
                                                struct GeoOpt *GeoOpt,
                                                char * folderName,
                                                SearchProblem * problem){

    // Get ground speeds
    double gs1 = groundspeed(&dubins->traj[0].wpt.pos, problem);
    double gs2 = groundspeed(&dubins->traj[1].wpt.pos, problem);
    double gs3 = groundspeed(&dubins->traj[3].wpt.pos, problem);
                                                    
    // Sampling frequency
    double dt = 1;

    // Average ground speed orbit 1
    double gsO1 = .5 * (gs1 + gs2) * KTS_2_FTS;

    // Sampling interval
    double interval = gsO1 * dt;

    // Get number of samples for each segment
    // First orbit
    double orbitLengthO1 = fabs(dubins->traj[0].dpsi*DEG_2_RAD) * fabs(dubins->traj[0].wpt.rad) * NM_2_FT;
    int steps = lround(orbitLengthO1 / interval);
    int nSampleO1 = steps + 1;
    
    // Straight
    interval = (gs2 * KTS_2_FTS) * dt;
    double dist, bearing;
    geo_dist(&dubins->traj[1].wpt.pos, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    steps = lround((dist*NM_2_FT) / interval);
    int nSampleS = steps + 1;

    // Offset the starting point of the orbit 2 given the remaining time within the straight segment
    double gsO2 = .5 * (gs2 + gs3) * KTS_2_FTS;
    double orbitLengthO2 = fabs(dubins->traj[2].dpsi*DEG_2_RAD) * fabs(dubins->traj[2].wpt.rad) * NM_2_FT ;

    // Second orbit
    interval = gsO2 * dt;
    steps = lround(orbitLengthO2 / interval);
    int nSampleO2 = steps + 1;

    // Total number of samples
    int ntot = (nSampleO1 + nSampleS + nSampleO2) - 2;
    *sampleSize = ntot;
    struct Pos *coordinates = (struct Pos *) malloc(ntot*sizeof(struct Pos));

    // START SAMPLING
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit1, &dubins->traj[0].wpt.pos, &dist, &bearing, GeoOpt);
    double dPsi;
    double hdg;
    double newhdg = dubins->traj[0].wpt.hdg;
    double R = fabs(dubins->traj[0].wpt.rad);
    double alt;
    double dh = dubins->traj[0].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    coordinates[0].t = 0;
    for (int i = 0; i < nSampleO1; i++) {
        
        dPsi = i * dubins->traj[0].dpsi/(nSampleO1 - 1);
        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[0].wpt.hdg + dPsi);
        alt = dubins->traj[0].wpt.pos.alt - i*dh/(nSampleO1 - 1);

        geo_npos(&dubins->orbit1, &coordinates[i], &R, &hdg, GeoOpt);
        coordinates[i].alt = alt;
        coordinates[i].hdg = newhdg;

        if (i > 0) coordinates[i].t = coordinates[i-1].t + dt;
    }

    // STRAIGHT SEGMENT
    // Sample
    int idx;
    double dlat = dubins->traj[2].wpt.pos.lat - dubins->traj[1].wpt.pos.lat;
    double dlon = dubins->traj[2].wpt.pos.lon - dubins->traj[1].wpt.pos.lon;
    double dalt = dubins->traj[2].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    for (int i = 1; i < nSampleS - 1; i++){

        idx = i + nSampleO1 - 1;

        coordinates[idx].lat = dubins->traj[1].wpt.pos.lat + i*dlat/(nSampleS - 1);
        coordinates[idx].lon = dubins->traj[1].wpt.pos.lon + i*dlon/(nSampleS - 1);
        coordinates[idx].alt = dubins->traj[1].wpt.pos.alt + i*dalt/(nSampleS - 1);
        coordinates[idx].hdg = newhdg;
        coordinates[idx].t = coordinates[idx-1].t + dt;
    }

    // SECOND ORBIT
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit2, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    
    dPsi = 0;
    hdg = 0;
    newhdg = 0;
    dh = dubins->traj[2].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;
    for (int i = 0; i < nSampleO2; i++) {
        idx = i + nSampleO1 + nSampleS - 2;
        dPsi = i * dubins->traj[2].dpsi/(nSampleO2 - 1);
        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[2].wpt.hdg + dPsi);
        geo_npos(&dubins->orbit2, &coordinates[idx], &R, &hdg, GeoOpt);
        coordinates[idx].alt = dubins->traj[2].wpt.pos.alt - i*dh/(nSampleO2 - 1);
        coordinates[idx].hdg = newhdg;
        coordinates[idx].t = coordinates[idx-1].t + dt;
    }

    if (folderName) {
        // Create a file
        FILE * file;

        // Open the file
        file = fopen(folderName, "w");
        if (!file) {
            perror("Can't open the dubins csv file...");
            return coordinates;
        }

        // Write coordinates
        for (int i = 0; i < ntot; i++) {
            fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f\n", coordinates[i].lat, coordinates[i].lon, coordinates[i].alt, coordinates[i].hdg, coordinates[i].t);
        }
        fclose(file);
    }

    return coordinates;
}

/*
   Return Dubins path coordinates and
   writes them to a CSV file.
*/
struct Pos *getDubinsCoordinatesWithFixedTimeStep(struct DubinsPath *dubins,
                                                struct Pos * initialSample,
                                                int *sampleSize,
                                                struct GeoOpt *GeoOpt,
                                                char * folderName,
                                                SearchProblem * problem){

    // Get ground speeds
    double gs1 = groundspeed(&dubins->traj[0].wpt.pos, problem);
    double gs2 = groundspeed(&dubins->traj[1].wpt.pos, problem);
    double gs3 = groundspeed(&dubins->traj[3].wpt.pos, problem);
                                                    
    // Sampling frequency
    double dt = 1;

    // Average ground speed orbit 1
    double gsO1 = .5 * (gs1 + gs2) * KTS_2_FTS; // ft/s

    // Sampling interval
    double interval = gsO1 * dt;               // ft

    // Get number of samples for each segment
    // First orbit
    double dt0 = dt;
    if (initialSample){
        // Start sampling first orbit
        // First step duration so that samples land on the next integer second
        double tfrac = fmod(initialSample->t, dt);
        if (tfrac < 0.0) tfrac += dt;
        dt0 = (tfrac > 1e-12) ? (dt - tfrac) : dt;
        if (dt0 < 1e-6) dt0 = dt;
    }
    double R1 = fabs(dubins->traj[0].wpt.rad); // NM (magnitude; sign handled by dpsi)
    double orbitLengthO1 = fabs(dubins->traj[0].dpsi*DEG_2_RAD) * R1 * NM_2_FT; // ft

    // Use floor so remaining time is in [0,dt)
    int steps = (int)floor(orbitLengthO1 / interval);
    int nSampleO1 = steps + 1;

    // How much traversal will be left in the first orbit
    double remainingO1Length = orbitLengthO1 - steps * gsO1 * dt; // ft
    double reaminingO1Time   = (gsO1 > 1e-12) ? (remainingO1Length / gsO1) : 0.0; // s

    // Straight
    interval = (gs2 * KTS_2_FTS) * dt; // ft per 1s

    double dist, bearing;

    // Second orbit
    double gsO2 = .5 * (gs2 + gs3) * KTS_2_FTS; // ft/s
    double R2 = fabs(dubins->traj[2].wpt.rad);  // NM (magnitude; sign handled by dpsi)
    interval = gsO2 * dt;                       // ft per 1s
    double orbitLengthO2 = fabs(dubins->traj[2].dpsi*DEG_2_RAD) * R2 * NM_2_FT; // ft

    // Total number of samples
    int ntot = 0;
    struct Pos *coordinates = NULL;

    // Get the bearing between the initial state and first orbit center
    if (initialSample){
        geo_dist(&dubins->orbit1, initialSample, &dist, &bearing, GeoOpt);
    } else {
        geo_dist(&dubins->orbit1, &dubins->traj[0].wpt.pos, &dist, &bearing, GeoOpt);
    }
    
    double dPsi;
    double hdg;
    double newhdg = dubins->traj[0].wpt.hdg;
    double R = R1;
    double alt;
    double dh = dubins->traj[0].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;

    // Precompute orbit1 per-second heading increment (deg) from kinematics
    double R1ft = R1 * NM_2_FT;
    double sign1 = (dubins->traj[0].dpsi >= 0.0) ? 1.0 : -1.0;
    double dPsiStepO1 = (R1ft > 1e-9) ? (sign1 * (gsO1 * dt / R1ft) * RAD_2_DEG) : 0.0;

    // Total orbit1 time for altitude interpolation
    double T1 = (gsO1 > 1e-12) ? (orbitLengthO1 / gsO1) : 0.0;

    // Temporary buffer for orbit1 samples so we can seed straight from the true last sample
    struct Pos *orbit1_tmp = (struct Pos *) malloc(nSampleO1*sizeof(struct Pos));
    if (!orbit1_tmp) return NULL;

    // Seed first sample
    if (initialSample){
        orbit1_tmp[0] = *initialSample;
        newhdg = initialSample->hdg;
    } else {
        orbit1_tmp[0].t = 0.0;
    }

    // Time-marched orbit1 so the first step can be dt0 (for S-turn stitching)
    dPsi = 0.0;
    for (int i = (initialSample ? 1 : 0); i < nSampleO1; i++) {

        double dt_step = dt;
        if (initialSample && i == 1) dt_step = dt0;

        if (i > 0) orbit1_tmp[i].t = orbit1_tmp[i-1].t + dt_step;

        // Advance along orbit by dt_step
        dPsi += (R1ft > 1e-9) ? (sign1 * (gsO1 * dt_step / R1ft) * RAD_2_DEG) : 0.0;
        if (fabs(dPsi) > fabs(dubins->traj[0].dpsi)) dPsi = dubins->traj[0].dpsi;

        hdg =  wrapTo360(bearing + dPsi);
        newhdg = wrapTo360(dubins->traj[0].wpt.hdg + dPsi);

        // altitude over full orbit duration (continuous) using elapsed time since start
        double frac = (T1 > 1e-12) ? ((orbit1_tmp[i].t - orbit1_tmp[0].t) / T1) : 1.0;
        if (frac < 0.0) frac = 0.0;
        if (frac > 1.0) frac = 1.0;
        alt = dubins->traj[0].wpt.pos.alt - frac*dh;

        geo_npos(&dubins->orbit1, &orbit1_tmp[i], &R, &hdg, GeoOpt);
        orbit1_tmp[i].alt = alt;
        orbit1_tmp[i].hdg = newhdg;
    }

    // Find the starting point of the straight segment given the remaining time within the first orbit
    // IMPORTANT: within the next 1s tick, first finish orbit1, then fly straight for the leftover time.
    struct Pos end_pos_o1;
    {
        double hdg_end = wrapTo360(bearing + dubins->traj[0].dpsi); // center-to-point bearing at orbit end
        geo_npos(&dubins->orbit1, &end_pos_o1, &R1, &hdg_end, GeoOpt);
        end_pos_o1.alt = dubins->traj[1].wpt.pos.alt;               // set explicitly (geo_npos may clobber)
        end_pos_o1.hdg = dubins->traj[1].wpt.hdg;
    }

    struct Pos init_pos_str;
    double d = (gs2 * KTS_2_FTS) * (dt - reaminingO1Time) * FT_2_NM; // NM traveled on straight within that 1s tick
    double straight_hdg = dubins->traj[1].wpt.hdg;                   // tangent track heading (deg)
    geo_npos(&end_pos_o1, &init_pos_str, &d, &straight_hdg, GeoOpt);

    // Altitude at init_pos_str: linear along full straight traj[1] -> traj[2]
    double dist_full_nm, bearing_full;
    geo_dist(&dubins->traj[1].wpt.pos, &dubins->traj[2].wpt.pos, &dist_full_nm, &bearing_full, GeoOpt);
    double dalt_full = dubins->traj[2].wpt.pos.alt - dubins->traj[1].wpt.pos.alt;
    double frac_init_str = (dist_full_nm > 1e-12) ? (d / dist_full_nm) : 0.0;
    if (frac_init_str < 0.0) frac_init_str = 0.0;
    if (frac_init_str > 1.0) frac_init_str = 1.0;
    init_pos_str.alt = dubins->traj[1].wpt.pos.alt + frac_init_str * dalt_full; // IMPORTANT: set alt
    init_pos_str.hdg = straight_hdg;

    // Now we can compute straight samples count from init_pos_str to traj[2]
    geo_dist(&init_pos_str, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    double straightLength = dist * NM_2_FT;                          // ft (remaining straight after transition)
    double gsS = (gs2 * KTS_2_FTS);                                  // ft/s
    steps = (int)floor(straightLength / (gsS * dt));
    int nSampleS = steps + 1;

    // How much traversal will be left in the straight
    double remainingSLength = straightLength - steps * gsS * dt;     // ft
    double remainingSTime   = (gsS > 1e-12) ? (remainingSLength / gsS) : 0.0;

    // Second orbit: within the next 1s tick, first finish straight, then enter orbit2 for leftover time
    double extraO2Arc = gsO2 * (dt - remainingSTime);                // ft traveled on orbit2 within that tick
    if (extraO2Arc < 0.0) extraO2Arc = 0.0;
    if (extraO2Arc > orbitLengthO2) extraO2Arc = orbitLengthO2;

    // Remaining orbit2 length after consuming extraO2Arc
    double orbit2Remaining = orbitLengthO2 - extraO2Arc;
    steps = (int)floor(orbit2Remaining / (gsO2 * dt));
    int nSampleO2 = steps + 1;

    // Total number of samples
    ntot = nSampleO1 + nSampleS + nSampleO2;
    *sampleSize = ntot;

    coordinates = (struct Pos *) malloc(ntot*sizeof(struct Pos));
    if (!coordinates) { free(orbit1_tmp); return NULL; }

    // Copy orbit1 tmp into output
    for (int i = 0; i < nSampleO1; i++) coordinates[i] = orbit1_tmp[i];
    free(orbit1_tmp); orbit1_tmp = NULL;

    // STRAIGHT SEGMENT
    // Sample
    int idx;

    // Seed first straight sample as the transition point one second after the last orbit1 sample
    idx = nSampleO1;
    coordinates[idx] = init_pos_str;
    coordinates[idx].hdg = straight_hdg;
    coordinates[idx].t = coordinates[idx-1].t + dt;

    // Step along the straight each second from init_pos_str toward traj[2]
    for (int i = 1; i < nSampleS; i++){

        idx = nSampleO1 + i;

        double step_nm = (gs2 * KTS_2_FTS) * dt * FT_2_NM;
        geo_npos(&coordinates[idx-1], &coordinates[idx], &step_nm, &bearing, GeoOpt);

        // IMPORTANT: geo_npos can overwrite alt -> restore it after geo_npos
        double dist_i_nm, bear_i;
        geo_dist(&dubins->traj[1].wpt.pos, &coordinates[idx], &dist_i_nm, &bear_i, GeoOpt);
        double frac = (dist_full_nm > 1e-12) ? (dist_i_nm / dist_full_nm) : 1.0;
        if (frac < 0.0) frac = 0.0;
        if (frac > 1.0) frac = 1.0;

        coordinates[idx].alt = dubins->traj[1].wpt.pos.alt + frac * dalt_full;
        coordinates[idx].hdg = straight_hdg;
        coordinates[idx].t = coordinates[idx-1].t + dt;
    }

    // SECOND ORBIT
    // Get the bearing between the initial state and first orbit center
    geo_dist(&dubins->orbit2, &dubins->traj[2].wpt.pos, &dist, &bearing, GeoOpt);
    
    dPsi = 0;
    hdg = 0;
    newhdg = 0;
    dh = dubins->traj[2].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Total orbit2 time for altitude interpolation
    double T3 = (gsO2 > 1e-12) ? (orbitLengthO2 / gsO2) : 0.0;
    double t0_o2 = (gsO2 > 1e-12) ? (extraO2Arc / gsO2) : 0.0;

    // Compute the initial angular offset on orbit2 corresponding to extraO2Arc
    double R2ft = R2 * NM_2_FT;
    double sign2 = (dubins->traj[2].dpsi >= 0.0) ? 1.0 : -1.0;
    double dPsi0_deg = (R2ft > 1e-9) ? (sign2 * (extraO2Arc / R2ft) * RAD_2_DEG) : 0.0;

    // Seed first orbit2 sample as the transition point one second after the last straight sample
    idx = nSampleO1 + nSampleS;
    hdg = wrapTo360(bearing + dPsi0_deg);
    newhdg = wrapTo360(dubins->traj[2].wpt.hdg + dPsi0_deg);
    geo_npos(&dubins->orbit2, &coordinates[idx], &R2, &hdg, GeoOpt);

    // Set alt AFTER geo_npos
    {
        double frac0 = (T3 > 1e-12) ? (t0_o2 / T3) : 1.0;
        if (frac0 < 0.0) frac0 = 0.0;
        if (frac0 > 1.0) frac0 = 1.0;
        coordinates[idx].alt = dubins->traj[2].wpt.pos.alt - frac0*dh;
    }

    coordinates[idx].hdg = newhdg;
    coordinates[idx].t = coordinates[idx-1].t + dt;

    // Advance orbit2 by 1-second kinematic increments thereafter
    double dPsiStepO2 = (R2ft > 1e-9) ? (sign2 * (gsO2 * dt / R2ft) * RAD_2_DEG) : 0.0;

    for (int i = 1; i < nSampleO2; i++) {
        idx = nSampleO1 + nSampleS + i;

        double dPsi_i = dPsi0_deg + i * dPsiStepO2;
        if (fabs(dPsi_i) > fabs(dubins->traj[2].dpsi)) dPsi_i = dubins->traj[2].dpsi;

        hdg =  wrapTo360(bearing + dPsi_i);
        newhdg = wrapTo360(dubins->traj[2].wpt.hdg + dPsi_i);

        geo_npos(&dubins->orbit2, &coordinates[idx], &R2, &hdg, GeoOpt);

        // Set alt AFTER geo_npos
        {
            double ti = t0_o2 + i*dt;
            double frac = (T3 > 1e-12) ? (ti / T3) : 1.0;
            if (frac < 0.0) frac = 0.0;
            if (frac > 1.0) frac = 1.0;
            coordinates[idx].alt = dubins->traj[2].wpt.pos.alt - frac*dh;
        }

        coordinates[idx].hdg = newhdg;
        coordinates[idx].t   = coordinates[idx-1].t + dt;
    }

    if (folderName) {
        // Create a file
        FILE * file;

        // Open the file
        file = fopen(folderName, "w");
        if (!file) {
            perror("Can't open the dubins csv file...");
            return coordinates;
        }

        // Write samples
        fprintf(file,"Lat [deg], Lon [deg], Alt [ft], Hdg [deg], Groundspeed [kts], Timestamp [s]\n");
        for (int i = 0; i < ntot; i++) {

            double gs_sample_kts;

            if (i < nSampleO1) {
                gs_sample_kts = gsO1 / KTS_2_FTS;   // convert ft/s → kts
            }
            else if (i < nSampleO1 + nSampleS) {
                gs_sample_kts = gs2;                // already in kts
            }
            else {
                gs_sample_kts = gsO2 / KTS_2_FTS;   // convert ft/s → kts
            }

            fprintf(file, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n",
                    coordinates[i].lat,
                    coordinates[i].lon,
                    coordinates[i].alt,
                    coordinates[i].hdg,
                    gs_sample_kts,
                    coordinates[i].t);
        }

        fclose(file);
    }

    return coordinates;
}


/*
    Generates left- and right-extended S-turn Dubins paths for a given Dubins path.
    Returns the path with the minimum overflown population risk.
*/
struct DubinsPath *getBestSturnPath(struct DubinsPath *dubins,
                                    double dAltitude,
                                    SearchProblem *problem,
                                    double *riskRuntime)
{
    int extendTo[2] = {-1, 1};

    // Allocate array for storing both left and right S-turn options
    struct DubinsPath **tmp_sturn = (struct DubinsPath **) calloc(2, sizeof(struct DubinsPath *));
    if (!tmp_sturn) {
        fprintf(stderr, "Memory allocation failed for tmp_sturn\n");
        exit(EXIT_FAILURE);
    }

    // Allocate result array: two connected paths
    struct DubinsPath *bestSturn = (struct DubinsPath *) malloc(2 * sizeof(struct DubinsPath));
    if (!bestSturn) {
        fprintf(stderr, "Memory allocation failed for bestSturn\n");
        free(tmp_sturn);
        exit(EXIT_FAILURE);
    }

    double minRisk = INFINITY;
    int bestIdx = -1;
    int feasible = 0;   // Feasibility flag
    double htol = 5;    // Altitude error tolerance [ft]

    for (int i = 0; i < 2; i++) {
        tmp_sturn[i] = computeSturnDubins(dubins, dAltitude, extendTo[i], problem);
        if (!tmp_sturn[i]) continue;

        double alt_change = tmp_sturn[i][0].traj[0].wpt.pos.alt - tmp_sturn[i][1].traj[3].wpt.pos.alt;
        if (fabs(alt_change - dAltitude) > htol) continue;

        feasible = 1;

        clock_t begin = clock();
        evaluateDubinsRisk(tmp_sturn[i], problem);
        // printf("tmp_sturn[%d].risk: %f\n", i, tmp_sturn[i]->risk);
        clock_t end = clock();
        *riskRuntime += (double)(end - begin) * 1000 / CLOCKS_PER_SEC;

        if (tmp_sturn[i][0].risk < minRisk) {
            minRisk = tmp_sturn[i][0].risk;
            bestIdx = i;
        }
    }

    if (!feasible || bestIdx == -1) {
        for (int i = 0; i < 2; i++) {
            Traj_InitArray(bestSturn[i].traj, 4);
            bestSturn[i].gp = NAN;
            bestSturn[i].ga = NAN;
            bestSturn[i].risk = NAN;
            bestSturn[i].feasible = false;
        }
    } else {
        for (int i = 0; i < 2; i++) {
            Traj_InitArray(bestSturn[i].traj, 4);
            Traj_CopyAll(tmp_sturn[bestIdx][i].traj, bestSturn[i].traj);
            bestSturn[i].hdist = tmp_sturn[bestIdx][i].hdist;
            bestSturn[i].orbit1 = tmp_sturn[bestIdx][i].orbit1;
            bestSturn[i].orbit2 = tmp_sturn[bestIdx][i].orbit2;
            bestSturn[i].type = tmp_sturn[bestIdx][i].type;
            bestSturn[i].gp = tmp_sturn[bestIdx][i].gp;
            bestSturn[i].ga = tmp_sturn[bestIdx][i].ga;
            bestSturn[i].risk = tmp_sturn[bestIdx][i].risk;
            bestSturn[i].feasible = true;
        }
    }

    for (int i = 0; i < 2; i++) {
        if (tmp_sturn[i]) {
            free(tmp_sturn[i]);
            tmp_sturn[i] = NULL;
        }
    }
    free(tmp_sturn);

    return bestSturn;
}


/*
    Computes the ground and airspace risks associated 
    with the given Dubins path and coordinate
*/
void evaluateDubinsRisk(struct DubinsPath *dubins,
                        SearchProblem * problem)
{

    // Compute ground risk
    if (problem->w_gp){
        dubins->gp = groundRisk(problem, dubins, NULL, 1);
    } else dubins->gp = -1;

    // Compute airspace risk
    dubins->ga = airspaceRisk(problem, dubins, NULL, 1);

    // Total risk
    if (dubins->size < 2) dubins->risk = problem->w_gp*dubins->gp + problem->w_ga*dubins->ga;
    else dubins[0].risk = problem->w_gp*dubins->gp + problem->w_ga*dubins->ga;
    
}

/*
    Iterates over Dubins-based paths (RSR, RSL, LSR, LSL).
    Computes S-turn paths in case altitude dissipation is needed.
    Returns the one with the minimum risk.
*/
struct DubinsPath *getMinRiskDubinsPath(struct DubinsPath *dubins,
                                        SearchProblem *problem,
                                        double *totalRuntime,
                                        double *riskRuntime)
{

    // Check if a Dubins path is feasible given the initial and final waypoints
    double dist, course;
    geo_dist(&dubins->traj[0].wpt.pos, &dubins->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);
    
    while (dist*NM_2_FT < 2*problem->ac->turnRadius){
        // Extend the goal state
        moveWaypoint(dubins, problem);
        geo_dist(&dubins->traj[0].wpt.pos, &dubins->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);
    }
   
    // Allocate memory for a temp. path structure (Four for all types: RSR, RSL, LSR, LSL)
    struct DubinsPath *tmp_path = (struct DubinsPath *) malloc(4*sizeof(struct DubinsPath));

    // Allocate memory for the best path (i.e., min. risk Dubin path)
    struct DubinsPath *best =  malloc(sizeof(struct DubinsPath));
    best->dgamma = 0;   // Initialize modification on gamma
    best->size = 0; // 0 = Invalid, 1 = Dubins, 2 = S-turn

    // Target altitude loss
    double targetdAltitude = dubins->traj[0].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Iterate over all feasible paths
    double min_gp = INFINITY, min_ga = INFINITY; // Initialize minimum risks
    double minRisk = INFINITY;  // Total risk
    int bestIndex;
    double riskRuntimeSturn = 0;

    // Start timer
    clock_t begin = clock();
    for (int type = 0; type < 4; type++){

        // printf(" *** Working on Case %d\n", type+1);

        // Initialize tmp_path and compute Dubins path
        Traj_InitArray(tmp_path[type].traj, 4);
        Traj_CopyAll(dubins->traj, tmp_path[type].traj);

        // Initialize the path variables
        tmp_path[type].type = type + 1;
        tmp_path[type].hdist = 0;
        tmp_path[type].size = 0;
        tmp_path[type].gp = INFINITY;
        tmp_path[type].ga = INFINITY;
        tmp_path[type].risk = INFINITY;
        
        // Perform Dubins calculations
        Dubins(&tmp_path[type], problem->DubinsOpt);
        traj_calctraj_angdist(&(tmp_path[type].traj[0]), 0, &(problem->DubinsOpt->trajopt));
        Traj_Calc3D(&(tmp_path[type].traj[0]), 0, &(problem->DubinsOpt->trajopt));

        // If the case cannot be solved, assign NAN to the horizontal path length.
        if (isnan(tmp_path[type].traj[3].wpt.pos.alt)) {tmp_path[type].hdist = NAN; continue;};

        // // Update the optimal turn flight path angles
        updateOptimalGammaTurn(&tmp_path[type], problem);
        double gamma1 = tmp_path[type].traj[0].wpt.gam;
        double gamma2 = tmp_path[type].traj[2].wpt.gam;

        // Clean initialization
        Traj_CopyAll(dubins->traj, tmp_path[type].traj);
        tmp_path[type].traj[0].wpt.gam = gamma1;
        tmp_path[type].traj[2].wpt.gam = gamma2;

        // Initialize the path variables
        tmp_path[type].type = type + 1;
        tmp_path[type].hdist = 0;
        tmp_path[type].size = 0;
        tmp_path[type].gp = INFINITY;
        tmp_path[type].ga = INFINITY;
        tmp_path[type].risk = INFINITY;

        // Recompute the Dubins path
        // printf("Recomputing...\n");
        tmp_path[type].hdist = 0;
        Dubins(&tmp_path[type], problem->DubinsOpt);
        traj_calctraj_angdist(&(tmp_path[type].traj[0]), 0, &(problem->DubinsOpt->trajopt));
        Traj_Calc3D(&(tmp_path[type].traj[0]), 0, &(problem->DubinsOpt->trajopt));

        // Check altitude loss
        double dh = tmp_path[type].traj[0].wpt.pos.alt - tmp_path[type].traj[3].wpt.pos.alt;

        // If the actual altitude loss is less than the target, compute S-turn paths
        double htol = -3;
        if (dh - targetdAltitude < htol) {
            struct DubinsPath *sturn = getBestSturnPath(&tmp_path[type], targetdAltitude, problem, riskRuntime);

            if (sturn[0].risk <= minRisk) {

                // Update the best risk
                min_gp = sturn[0].gp;
                min_ga = sturn[0].ga;
                minRisk = sturn[0].risk;

                // Allocate memory for the best path
                if (best) {
                    best = (struct DubinsPath *) realloc(best, 2*sizeof(struct DubinsPath));
                } else {
                    best = (struct DubinsPath *) malloc(2*sizeof(struct DubinsPath));
                }

                // Copy the best path
                for (int i = 0; i < 2; i++) {
                    Traj_InitArray(best[i].traj, 4);
                    Traj_CopyAll(sturn[i].traj, best[i].traj);
                    best[i].hdist = sturn[i].hdist;
                    best[i].orbit1 = sturn[i].orbit1;
                    best[i].orbit2 = sturn[i].orbit2;
                    best[i].type = sturn[i].type;
                    best[i].dgamma = 0;
                    best[i].size = 2;
                    best[i].ga = min_ga;
                    best[i].gp = min_gp;
                    best[i].risk = minRisk;
                }
                free(sturn); sturn = NULL;
            }
        } else {
        
            // Reduce the path slope, if the altitude loss is greater than the target.
            if (dh - targetdAltitude > -htol) {
                // If slope can be reduced, compute the path risk and compare it. Continue otherwise.
                if (reduceSlope(&tmp_path[type], targetdAltitude, problem)) {
                    continue;
                }
            }

            // Get path coordinates
            if (!isnan(tmp_path[type].hdist)) {

                clock_t beginRisk = clock();
                evaluateDubinsRisk(&tmp_path[type], problem);
                clock_t endRisk = clock();
                *riskRuntime += (double) (endRisk - beginRisk)*1000 / CLOCKS_PER_SEC;
                
            } else {
                tmp_path[type].feasible = false;
                tmp_path[type].gp = NAN;
                tmp_path[type].ga = NAN;
                tmp_path[type].risk = NAN;
                continue;
            }

            // Compare
            if (tmp_path[type].risk <= minRisk) {

                // Update the best risks
                min_gp = tmp_path[type].gp;
                min_ga = tmp_path[type].ga;
                minRisk = tmp_path[type].risk;

                // Allocate memory
                if (best != NULL) {
                    best = (struct DubinsPath *) realloc(best, sizeof(struct DubinsPath));
                } else {
                    best = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));
                }

                // Hard copy the best path
                Traj_InitArray(best->traj, 4);
                Traj_CopyAll(tmp_path[type].traj, best->traj);
                best->hdist = tmp_path[type].hdist;
                best->orbit1 = tmp_path[type].orbit1;
                best->orbit2 = tmp_path[type].orbit2;
                best->type = tmp_path[type].type;
                best->dgamma = tmp_path[type].dgamma;
                best->size = 1;
                best->ga = tmp_path[type].ga;
                best->gp = tmp_path[type].gp;
                best->risk = tmp_path[type].risk;
            }
        }
    }

    // If none of the dubins paths is feasible, return NAN risk
    if (best->size == 0) {
        best->ga = NAN;
        best->gp = NAN;
        best->risk = NAN;
    }

    clock_t end = clock();
    *totalRuntime = (double) (end - begin)*1000 / CLOCKS_PER_SEC;

    // Free allocated memory
    free(tmp_path); tmp_path = NULL;

    return best;
}

/*
    Reduces the gliding angle of the straight segment of
    a given Dubins Path.
*/
int reduceSlope(struct DubinsPath *dubins, double dAltitudeTarget, SearchProblem *problem)
{
    // Actual altitude loss
    double dAltitudeActual = dubins->traj[0].wpt.pos.alt - dubins->traj[3].wpt.pos.alt;

    // Deficit altitude
    double deficit = dAltitudeActual - dAltitudeTarget;

    // How much altitude the path loses along the straight segment
    double dAltStraight = dubins->traj[1].hdist*NM_2_FT * tan(fabs(dubins->traj[1].wpt.gam)*DEG_2_RAD);

    // How much altitude the path must lose throughout the straight segment to reach the goal altitude
    double dh = dAltStraight - deficit;

    // Get the new gamma
    double gamma = atan(dh/(dubins->traj[1].hdist*NM_2_FT)) * RAD_2_DEG;

    // Get feasible flight path angle range
    double *gammaArray = getOptimalGamma(&dubins->traj[1].wpt.hdg, problem);

    // Check if the new gamma is feasible
    if ((gamma >= gammaArray[0]) && (gamma <= gammaArray[1]))
    {
        dubins->traj[1].wpt.gam = -gamma;
        dubins->dgamma = gamma - gammaArray[2];

        free(gammaArray); gammaArray = NULL;
        dubins->feasible = true;
        return EXIT_SUCCESS;
    }

    dubins->feasible = false;
    free(gammaArray); gammaArray = NULL;
    return EXIT_FAILURE;
}

/*
    Moves the goal state to make Dubins path feasible
    in case of close initial and final waypoints
*/
void moveWaypoint(struct DubinsPath *dubins, SearchProblem *problem)
{   
    // Find the new coordinate
    double dist = 0.1*problem->ac->turnRadius*FT_2_NM;
    double course = wrapTo360(dubins->traj[3].wpt.hdg - 180);
    double initalt = dubins->traj[3].wpt.pos.alt;
    geo_npos(&dubins->traj[3].wpt.pos, &dubins->traj[3].wpt.pos, &dist, &course, &problem->GeoOpt);

    // Find the new altitude
    double *gammaArray = getOptimalGamma(&dubins->traj[3].wpt.hdg, problem);
    dubins->traj[3].wpt.pos.alt = initalt + (dist*NM_2_FT)*tan(gammaArray[2]*DEG_2_RAD);
    free(gammaArray);
}

/*
    Returns the shortest Turn-Straight-Turn Dubins path
*/
int getDubinsWithType(struct DubinsPath * const dubins, int type, struct DubinsOpt * const opt) {

    // Allocate memory for tmp_path inside the loop
    struct DubinsPath *tmp_path = (struct DubinsPath *) malloc(sizeof(struct DubinsPath));

    // Compute the path
    Traj_InitArray(tmp_path->traj, 4);
    Traj_CopyAll(dubins->traj, tmp_path->traj);

    tmp_path->type = type + 1;
    tmp_path->hdist = 0;
    
    // Perform calculations
    Dubins(tmp_path, opt);
    traj_calctraj_angdist(&(tmp_path->traj[0]), 0, &(opt->trajopt));
    Traj_Calc3D(&(tmp_path->traj[0]), 0, &(opt->trajopt));

    Traj_CopyAll(tmp_path->traj, dubins->traj);
    dubins->hdist = tmp_path->hdist;
    dubins->orbit1 = tmp_path->orbit1;
    dubins->orbit2 = tmp_path->orbit2;
    dubins->type = tmp_path->type;
    free(tmp_path); tmp_path = NULL;

    return EXIT_SUCCESS;
}

#include "../../dubins/src/dubins_optimal.inc"
