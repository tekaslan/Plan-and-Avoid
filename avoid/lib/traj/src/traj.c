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
#    Trajectory Handling Functions                                                  %
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Author  : Pedro Di Donato                                                      %
#    Date    : 2017                                                                 %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=UCxHXTgAAAAJ&hl=en %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
#include "../include/traj.h"

int Traj_Copy(struct Traj *p, struct Traj *c){

    if(p && c){
        c->wpt   = p->wpt;
        c->wpt.rad   = p->wpt.rad;
        c->wpt.gam   = p->wpt.gam;
        c->wpt.phi   = p->wpt.phi;
        c->wpt.V     = p->wpt.V;
        c->wpt.CL    = p->wpt.CL;
        c->dpsi  = p->dpsi;
        c->hdist = p->hdist;
        c->nturn  = p->nturn;
        return EXIT_SUCCESS;
    }
    else{
        return EXIT_FAILURE;
    }

}

int Traj_CopyAll(struct Traj *p, struct Traj *c){

    int i = 0;
    for(; ; p = p->next, c = c->next, i++)
        if(Traj_Copy(p,c))  return i;
            
}

int Traj_InitArray(struct Traj *p, int npoints){
    
    int i = 0;
    for(; i < (npoints-1); i++, p = p->next){
        memset(p,0,sizeof(struct Traj));
        p->wpt.hdg = NAN;
        p->wpt.pos.alt = NAN;
        p->wpt.gam = NAN;
        p->wpt.phi = NAN;
        p->wpt.V = NAN;
        p->wpt.CL = NAN;
        p->hdist = NAN;
        p->next = (p+1);
    }
    memset(p,0,sizeof(struct Traj));
    p->wpt.hdg = NAN;
    p->wpt.gam = NAN;
    p->wpt.phi = NAN;
    p->wpt.V = NAN;
    p->wpt.CL = NAN;
    p->hdist = NAN;

    return EXIT_SUCCESS;
}

double Traj_HDist(struct Traj* p, int npoints, const struct TrajOpt *opt){

    int i = 0;
    if(!npoints) npoints = pow(2,8);
    double dist = 0.0;

    for(;((p->next != NULL) && (i < npoints)); p = p->next, i++){
        if(isnan(p->hdist))
            traj_calctraj_angdist(p,npoints,opt);
        dist += p->hdist;
    }
    
    return dist;
}

int traj_calctraj_angdist(struct Traj* p, int npoints, const struct TrajOpt *opt){

    int i = 0;
    double phi = 0.0;
    if(!npoints) npoints = pow(2,8);

    for(; ((p->next != NULL) && (i < npoints)); p = p->next, i++){
        if (fabs(p->wpt.rad) < ZERO_TOLERANCE)
        {
            if (geo_dist(&p->wpt.pos,&(p->next)->wpt.pos,&(p->hdist),&phi,&(opt->geoopt)) != 0)
                return -1;
        }
        else{
            p->dpsi = fmod((p->next)->wpt.hdg - p->wpt.hdg,360.0);
            if (p->dpsi*p->wpt.rad < 0){
                if (p->dpsi < 0)
                    p->dpsi += 360.0;
                else
                    p->dpsi -= 360.0;
            }
            p->hdist = fabs(p->wpt.rad*(p->dpsi + p->nturn*360.0)*DEG_2_RAD);
        }
    }

    return i;
}

/**
 Implements the function presented in \cite beard2012small_chp11
 for computing the angular distance between two angles. 

 If CCW the angular distance returned is negative. That is different that 
 the reference to be compatible with legacy code.

 Based on the Figure 11.8 of the reference (reproduced here):
 \image html ang_distance.png
 \image latex ang_distance.pdf 

 CW Case:
 \f[
    \Delta\vartheta = 360 + \vartheta_2 - \vartheta_1 \bmod 360
 \f]
 
 CCW Case:
 \f[
    \Delta\vartheta = -\left(360 + \vartheta_1 - \vartheta_2 \bmod 360\right)
 \f]
*/  
int traj_calcangdist(const double * const ang1, const double * const ang2, 
                     const double * const direction, 
                     double * const angdist)
{
    if(*direction > 0){   /* CW */
        *angdist = fmod(360.0 + *ang2 - *ang1, 360.0);
        return 0;
    } else if (*direction < 0){ /* CCW */
        *angdist = -fmod(360.0 - *ang2 + *ang1, 360.0);
        return 0;
    }
    else return 1; /* Error */
}
 

int traj_semiplan(const double * const ang1, const double * const ang2)
{
    double direction = 1.0;
    double dist_cw = 0.0;
    traj_calcangdist(ang1, ang2, &direction, &dist_cw);
    
    if(dist_cw > 90.0 && dist_cw < 270.0)   return 1;
    else                                    return 0;
}

int Traj_Calc3D(struct Traj* p, int npoints, const struct TrajOpt *opt){

    int i = 0;
    if(!npoints) npoints = pow(2,8);

    for(; ((p->next != NULL) && (i < npoints)); p = p->next, i++){
        if(isnan(p->hdist))
            traj_calctraj_angdist(p,npoints,opt);

        /*
            H. Emre Tekaslan: This section is commented out as it enforces initial
            altitude values to be NaN. Altitudes of waypoints must be computed
            regardless of initial altitude values.
        */
        // ==================================================================================================
        // if(isnan(p->next->wpt.pos.alt)){
            // p->next->wpt.pos.alt = p->wpt.pos.alt + p->hdist*tan(p->wpt.gam*DEG_2_RAD)*NM_2_FT;
        // }
        // else if(isnan(p->wpt.gam)){
        //     p->wpt.gam = atan2(p->next->wpt.pos.alt - p->wpt.pos.alt, p->hdist * NM_2_FT) * RAD_2_DEG;
        // }
        // ==================================================================================================

        // Compute the altitude
        p->next->wpt.pos.alt = p->wpt.pos.alt + p->hdist*tan(p->wpt.gam*DEG_2_RAD)*NM_2_FT;
    }

    return i;
}

struct Pos* Traj_ToPlot(struct Traj* traj, int *npoints, struct TrajOpt *opt){

    int i = 0, j = 0, k = 0;
    double dpsip = 0.1;
    double hdg;
    struct Traj * p = traj;
    struct Pos *pos = NULL;
    char bflag = 0;

    *npoints = 0;

    for(; p != NULL; p = p->next){
        /* Code to Avoid Infinity Loop on printing */
        if((p == traj->next) && (bflag == 3)) break;
        if((p == traj) || (p == traj->next)) bflag++; 


        (*npoints)++;
        if(fabs(p->dpsi) > 0.0)
            (*npoints) += floor(fabs(p->dpsi/dpsip));
    }
    
    pos = (struct Pos *) malloc(*npoints*sizeof(struct Pos));
    bflag = 0;

    for(i = 0, p = traj; p != NULL; p = p->next, i++){

        /* Code to Avoid Infinity Loop on printing */
        if((p == traj->next) && (bflag == 3)) break;
        if((p == traj) || (p == traj->next)) bflag++; 

        pos[i] = p->wpt.pos;
        if(p->dpsi > 0) dpsip = 0.1;
        else            dpsip = -0.1;
        if(fabs(p->dpsi) > 0.0){
            hdg = p->wpt.hdg;
            k = floor(fabs(p->dpsi/dpsip)) + 1;
            for(j = 1; j < k; j++){
                i++;
                geo_npost(&(pos[i-1]),&(pos[i]),&(p->wpt.rad),&(hdg),&dpsip,&(opt->geoopt));
                pos[i].alt = pos[i-1].alt + tan(p->wpt.gam*DEG_2_RAD)*fabs(p->wpt.rad * dpsip * DEG_2_RAD)*NM_2_FT;
                hdg += dpsip;
            }
        }
    }

    return pos;
}

int Traj_Print(FILE *output, struct Traj *p, int npoints){

    int i = 0;
    if(!npoints) npoints = pow(2,8);

    fprintf(output, "#Waypoint, Lat (deg), Lon (deg), Alt (ft), Hdg (deg), "
                    "Rad (NM), Gam (deg), V (KTAS), "
                    "Delta Psi(deg), Hdist (NM), Turns #\n");

    for(; ((p != NULL) && (i < npoints)); p = p->next, i++){

        fprintf(output, "%d,%lf,%lf,%lf,%lf,"
                        "%lf,%lf,%lf,"
                        "%lf,%lf,%d\n",
                        i, p->wpt.pos.lat, p->wpt.pos.lon, p->wpt.pos.alt, p->wpt.hdg,
                        p->wpt.rad, p->wpt.gam, p->wpt.V,
                        p->dpsi, p->hdist, p->nturn);
                        
    }

    return i;

}

double Traj_LayerTime(struct Traj* traj, double layer_min, double layer_max, struct TrajOpt *opt){

    struct Traj *p = traj;
    double time = 0.0;
    int npoints = 0;

    for(; p->next != NULL; p = p->next){

        /* Check if Segment Distance was Calculated */    
        if(p->hdist == 0.0)
            traj_calctraj_angdist(p,npoints,opt);
 
        /* Initial Segment WayPoint in the Layer */
        if((p->wpt.pos.alt > layer_min) && (p->wpt.pos.alt < layer_max)){
            /* Final Segment Waypoint also in the Layer */ 
            if((p->next->wpt.pos.alt >= layer_min) && (p->next->wpt.pos.alt <= layer_max)){
                time += p->hdist*cos(p->wpt.gam*DEG_2_RAD)*NM_2_M/p->wpt.V;
            }
            /* Final Segment Waypoint not in the Layer */
            else{
                /* Final Segment Waypoint Below the Layer */
                if(p->next->wpt.pos.alt < layer_min){
                    time += p->hdist*cos(p->wpt.gam*DEG_2_RAD)*NM_2_M/p->wpt.V*
                            (layer_min - p->wpt.pos.alt)/(p->next->wpt.pos.alt - p->wpt.pos.alt);
                }
                /* Final Segment Waypoint Above the Layer */
                else{
                    time += p->hdist*cos(p->wpt.gam*DEG_2_RAD)*NM_2_M/p->wpt.V*
                            (layer_max - p->wpt.pos.alt)/(p->next->wpt.pos.alt - p->wpt.pos.alt);
                }
            }
        }
        /* Initial Segment WayPoint NOT in the Layer while Final Segment Waypoint in the Layer*/
        else{
            if((p->next->wpt.pos.alt >= layer_min) && (p->next->wpt.pos.alt <= layer_max)){
                /* Final Segment Waypoint Below the Layer */
                if(p->next->wpt.pos.alt < layer_min){
                    time += p->hdist*cos(p->wpt.gam*DEG_2_RAD)*NM_2_M/p->wpt.V*
                            (layer_min - p->next->wpt.pos.alt)/(p->wpt.pos.alt - p->next->wpt.pos.alt);
                }
                /* Final Segment Waypoint Above the Layer */
                else{
                    time += p->hdist*cos(p->wpt.gam*DEG_2_RAD)*NM_2_M/p->wpt.V*
                            (layer_max - p->next->wpt.pos.alt)/(p->wpt.pos.alt - p->next->wpt.pos.alt);
                }
            }
        }
    }

    return time;
}


void traj_file(struct Traj* traj, struct TrajOpt *opt){

    double t = 0.0;
    struct Traj * p = traj;
    int i = 0;
    FILE * output = fopen("traj.dat","w");

    for(p = traj; p->next; p = p->next){
        double trel = 0.0;
        double hdg = 0.0;
        struct Pos pos = p->wpt.pos;

        /* Avoid empty segments */
        if(p->hdist < 1e-8) continue;

        for(i = 0; i <=10; i++){
    
            trel = i*p->hdist* 3600.0 / (10*p->wpt.V*cos(p->wpt.gam*DEG_2_RAD));
            if(!i && p != traj) trel += 0.000001;
            double dist = p->wpt.V*cos(p->wpt.gam*DEG_2_RAD)*trel/3600.0;
            if(!p->wpt.rad){
                hdg = p->wpt.hdg;
                geo_npos(&(p->wpt.pos),&pos,&dist,&(p->wpt.hdg),&(opt->geoopt));
            }
            else{
                double dpsi = (dist / p->wpt.rad) * RAD_2_DEG;
                hdg = p->wpt.hdg + dpsi;
                geo_npost(&(p->wpt.pos),&pos,&(p->wpt.rad),&(p->wpt.hdg),&dpsi,&(opt->geoopt));
            }        
            pos.alt = p->wpt.pos.alt + tan(p->wpt.gam*DEG_2_RAD)*fabs(dist)*NM_2_FT;

            fprintf(output,"%lf,"
                           "%lf,%lf,%lf,"
                           "%lf,%lf,%lf,%lf,"
                            "%lf,%lf\n",
                            t + trel,
                            pos.lat*NM_2_M, pos.lon*NM_2_M, -pos.alt*FT_2_M,
                            p->wpt.V*KT_2_MS, p->wpt.gam*DEG_2_RAD,
                            cos(hdg*DEG_2_RAD), sin(hdg*DEG_2_RAD),
                            p->wpt.CL, p->wpt.phi*DEG_2_RAD);
        }   

        t += trel;
    }

    fprintf(output,"%lf,"
                   "%lf,%lf,%lf,"
                   "%lf,%lf,%lf,%lf,"
                    "%lf,%lf\n",
                    t + 0.000001 ,
                    p->wpt.pos.lat*NM_2_M, p->wpt.pos.lon*NM_2_M, -p->wpt.pos.alt*FT_2_M,
                    p->wpt.V*KT_2_MS, p->wpt.gam*DEG_2_RAD,
                    cos(p->wpt.hdg*DEG_2_RAD), sin(p->wpt.hdg*DEG_2_RAD),
                    p->wpt.CL, p->wpt.phi*DEG_2_RAD);

    fclose(output);                                        
    
}
            
