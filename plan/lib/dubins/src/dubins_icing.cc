/**
 * @file
 *
 * @brief Contains the function to calculate the Dubins path
 */

#include "../include/dubins.h"

int Dubins_IcingDescent(Traj * traj, double layermin, DubinsOpt *opt){

    DubinsPath dubins;
    Traj final;
    Traj *p = NULL;
    int i = 0;
    double alt_error = 99.0;
    double dist = 1.0;

    /* Initialize Dubins Array */
    dubins.type = OPT;
    Traj_InitArray(dubins.traj,4);

    /* Initial position for Dubins */
    if(traj->wpt.alt > layermin){
        dist = (traj->wpt.alt - layermin)/tan(traj->gam*DEG_2_RAD)*FT_2_NM;
    }      
    else{
        return EXIT_FAILURE;
    }

    Geo_Npos(&(traj->wpt),&(dubins.traj[0].wpt), &(dist), 
             &(traj->wpt.hdg), &(opt->trajopt.geoopt));
    dubins.traj[0].wpt.hdg = traj->wpt.hdg; 
    dubins.traj[0].wpt.alt = layermin; 
    
    /* Final Position is Given */
    final = *(traj->next->next->next->next->next);
   
    /* Final Segment must the a straight line */
    dubins.traj[3].wpt.hdg = final.wpt.hdg;    
    
    /* Radius and Gammas are also given*/
    for(i = 0, p = traj->next; i < 4; i++, p = p->next){
        dubins.traj[i].rad = p->rad;
        dubins.traj[i].gam = p->gam;
    }
    
    dubins.traj[3].hdist = 0.0;    
    for(i = 0; (alt_error > 1) && (i < 100); i++){
        
        /* Calculate last position */ 
        dist = -dubins.traj[3].hdist;
        Geo_Npos(&(final.wpt),&(dubins.traj[3].wpt), &(dist), 
                 &(final.wpt.hdg), &(opt->trajopt.geoopt));
        
        /* Calculate Dubins */
        Dubins(&dubins,opt);

        /* Calculate Third Dimension of Dubins */
        traj_calctraj_angdist(&(dubins.traj[0]),4,&(opt->trajopt));
        Traj_Calc3D(&(dubins.traj[0]),4,&(opt->trajopt));

        /* Calculate Altitude Error */
        alt_error = dubins.traj[3].wpt.alt - 
                    (final.wpt.alt + dist*tan(dubins.traj[3].gam*DEG_2_RAD)*NM_2_FT);
      
        if(alt_error > (-2.0*M_PI*dubins.traj[2].rad*tan(dubins.traj[2].gam*DEG_2_RAD)*NM_2_FT) + 500) {
            dubins.traj[2].turn += 1;
            continue;
        }
 
        /* DEBUG STATEMENTS 
        printf("final dubins alt = %f, desir alt = %f, final alt = %f, error alt = %f, dist = %f\n", dubins.traj[3].wpt.alt, 
                (final.wpt.alt + dist*tan(dubins.traj[3].gam*DEG_2_RAD)*NM_2_FT), final.wpt.alt,
                alt_error, dist); */
        
        /* Update Hill Climbing */
        if((dist = 0.0) && (alt_error < 0))
            return EXIT_FAILURE;
        else
            dubins.traj[3].hdist += alt_error*FT_2_NM;
    }
     
    for(i = 0, p = traj->next; i < 4; i++, p = p->next){
        p->wpt = dubins.traj[i].wpt;
        p->rad = dubins.traj[i].rad;
        p->gam = dubins.traj[i].gam;
        p->dpsi = dubins.traj[i].dpsi;
        p->hdist = dubins.traj[i].hdist;
        p->turn = dubins.traj[i].turn;
    }
    *p = final;

    return EXIT_SUCCESS;
}
