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
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : Pedro Di Donato & H. Emre Tekaslan (tekaslan@vt.edu)                 %
#    Date    : April 2025                                                           %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#                      https://scholar.google.com/citations?user=UCxHXTgAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
#include "dubins.h"

static bool dubins_optimal(struct DubinsPath * const dub,
                          const struct DubinsOpt * opt);

int Dubins(struct DubinsPath * const dubins,
           const struct DubinsOpt * opt)

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
int shortestDubins(struct DubinsPath * dubins, const struct DubinsOpt * const opt)
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
struct Pos *getDubinsCoordinatesWithFixedTimeStep(struct DubinsPath *dubins,
                                                struct Pos * initialSample,
                                                int *sampleSize,
                                                const struct GeoOpt *GeoOpt,
                                                const double groundspeed){

    if (!isfinite(groundspeed) || groundspeed <= 1e-6) {
        *sampleSize = 0;
        return NULL;
    }
    // Get ground speeds
    double gs1 = groundspeed;
    double gs2 = groundspeed;
    double gs3 = groundspeed;
                                                    
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
    if (!isfinite(R1) || R1 <= 0.0) {
        *sampleSize = 0;
        return NULL;
    }
    double orbitLengthO1 = fabs(dubins->traj[0].dpsi*DEG_2_RAD) * R1 * NM_2_FT; // ft
    if (!isfinite(orbitLengthO1) || orbitLengthO1 < 0.0) {
        *sampleSize = 0;
        return NULL;
    }

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
    if (!isfinite(R2) || R2 <= 0.0) {
        *sampleSize = 0;
        return NULL;
    }
    interval = gsO2 * dt;                       // ft per 1s
    double orbitLengthO2 = fabs(dubins->traj[2].dpsi*DEG_2_RAD) * R2 * NM_2_FT; // ft
    if (!isfinite(orbitLengthO2) || orbitLengthO2 < 0.0) {
        *sampleSize = 0;
        return NULL;
    }

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
        orbit1_tmp[i].alt = alt;            // IMPORTANT: overwrite alt after geo_npos
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
    if (nSampleO1 <= 0 || nSampleS <= 0 || nSampleO2 <= 0) {
        *sampleSize = 0;
        return NULL;
    }

    if (nSampleO1 > INT_MAX - nSampleS) {
        *sampleSize = 0;
        return NULL;
    }
    if (nSampleO1 + nSampleS > INT_MAX - nSampleO2) {
        *sampleSize = 0;
        return NULL;
    }

    ntot = nSampleO1 + nSampleS + nSampleO2;

    if (ntot <= 0 || ntot > 10000000) {
        *sampleSize = 0;
        return NULL;
    }
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

    // IMPORTANT: set alt AFTER geo_npos
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

        // IMPORTANT: set alt AFTER geo_npos
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

    return coordinates;
}

/*
    Returns the shortest Turn-Straight-Turn Dubins path
*/
int getDubinsWithType(struct DubinsPath * dubins, int type, const struct DubinsOpt * opt) {

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

#ifdef __cplusplus
}
#endif
