#include "../include/dubins.h"
#include <time.h>
// #include "gmaps.h"

int Test_Dubins(){

    struct DubinsPath dubins;
    struct DubinsOpt opt;
    struct Traj *p;
    struct Pos *pos = NULL;
    int npoints;

    Traj_InitArray(dubins.traj,4);

    /* Initial Position */
    dubins.traj[0].wpt.pos.lat = 38.84150395;
    dubins.traj[0].wpt.pos.lon = -77.29862850;
    dubins.traj[0].wpt.pos.alt = 5217.60;
    dubins.traj[0].wpt.hdg = 200.1;

    /* Final Position */
    dubins.traj[3].wpt.pos.lat = 38.73137684;
    dubins.traj[3].wpt.pos.lon = -77.20341034;
    dubins.traj[3].wpt.pos.alt = 516.33;
    dubins.traj[3].wpt.hdg = 133.0;

    /* Radius */
    dubins.traj[0].wpt.rad = 0.2469;
    dubins.traj[1].wpt.rad = 0.0;
    dubins.traj[2].wpt.rad = 0.2469;

    /* Flight Path Angle */
    dubins.traj[0].wpt.gam = -4.55;
    dubins.traj[1].wpt.gam = -4.4;
    dubins.traj[2].wpt.gam = -4.55;

    /* Dubins Type Optimal */
    dubins.type = OPT;

    /* Options */
    opt.verbose = 0;
    opt.trajopt.verbose = 0;
    opt.trajopt.geoopt.verbose = 0;
    opt.trajopt.geoopt.model = WGS84;

    /* Calculate Dubins Path */
    double total_time = 0;
    clock_t begin = clock();
    Dubins(&dubins,&opt);
    clock_t end = clock();
    double time = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Dubins Elapsed time: %lf \n",time*1000);
    traj_calctraj_angdist(&(dubins.traj[0]),0,&(opt.trajopt));
    Traj_Calc3D(&(dubins.traj[0]),0,&(opt.trajopt));
    pos = Traj_ToPlot(dubins.traj, &npoints, &(opt.trajopt));

    /* Plot Trajectory */ 
    // struct Gmap gmap;
    // if(gmaps_initialize(&gmap, "data/fxb.csv")) return EXIT_FAILURE;
    // strcpy(gmap.outfile,"results/dubins.png");
    // if(gmaps_polyline(&gmap,pos,npoints))     return EXIT_FAILURE; 
 
    /* Print Output */
    printf("\n\t ### Test Dubins ###\n\n");
    for(p = dubins.traj; p; p = p->next){
        printf("lat: %f, lon: %f, alt: %.1f, hdg: %.0f, rad: %f,"
                "dpsi = %f, hdist = %f\n",
                p->wpt.pos.lat, p->wpt.pos.lon, p->wpt.pos.alt, p->wpt.hdg, 
                p->wpt.rad, p->dpsi, p->hdist);
    }
    printf("\n\n");
    
    return EXIT_SUCCESS;
}

int main(){
    
    Test_Dubins();

    return EXIT_SUCCESS;
}
