#include "../include/traj.h"

int main (int argc, char * argv[]){

    FILE *input = NULL;
    FILE *output = NULL;
    struct Traj* traj = NULL;
    int n = 0, i = 0;
    char turn, ns, ew;
    int lond, lonm, latd, latm;
    double lons, lats;
    double chk;
    const struct TrajOpt opt = {.geoopt = {.model = WGS84, .verbose = 0}, .verbose = 0};
    char name[10] = {0};
    
    if(!(input = fopen(argv[1],"r"))){
        fprintf(stderr,"Impossible to open input file: %s", argv[1]);
        return EXIT_FAILURE;
    }

    if(!(output = fopen(argv[2],"w"))){
        fprintf(stderr,"Impossible to open output file: %s", argv[2]);
        return EXIT_FAILURE;
    }

    fscanf(input,"%d\n", &n);

    traj = (struct Traj *) malloc(n*sizeof(struct Traj));
    Traj_InitArray(traj, n);

    for(i = 0; i < n; i++){
        fscanf(input, "%c %d %d %lf%c %d %d %lf%c %lf %lf",
              &turn, 
              &latd, &latm, &lats, &ns,
              &lond, &lonm, &lons, &ew,
              &(traj[i].wpt.pos.alt), &(traj[i].wpt.V));
        fgets(name, 10, input);
        traj[i].wpt.pos.lat = latd + ((double) latm) / 60.0 + lats / 3600.0;
        if(ns == 'S'){
            traj[i].wpt.pos.lat *= -1.0;
        }

        traj[i].wpt.pos.lon = lond + ((double) lonm) / 60.0 + lons / 3600.0;
        if(ew == 'W'){
            traj[i].wpt.pos.lon *= -1.0;
        }

        if(i > 0){
            if(turn == 'S'){
                if(!isnan(traj[i-1].wpt.hdg)){
                    chk = traj[i-1].wpt.hdg;
                    geo_dist(&(traj[i-1].wpt.pos), &(traj[i].wpt.pos), &(traj[i-1].hdist), &(traj[i-1].wpt.hdg), &(opt.geoopt));
                    if(fabs(traj[i-1].wpt.hdg - chk) > 0.1){
                        fprintf(stderr,"Heading Discrepance between end of turn and beginning of straight segments: %d-%d\n",
                                        i, i+1);
                    }
                }
                else{
                    geo_dist(&(traj[i-1].wpt.pos), &(traj[i].wpt.pos), &(traj[i-1].hdist), &(traj[i-1].wpt.hdg), &(opt.geoopt));
                }                    
                traj[i-1].wpt.rad = 0.0;
            }
            else{
                if(isnan(traj[i-1].wpt.hdg)){
                    traj[i-1].wpt.hdg = traj[i-2].wpt.hdg;
                }
                geo_distt(&(traj[i-1].wpt.pos), &(traj[i].wpt.pos), &(traj[i-1].wpt.rad), &(traj[i-1].wpt.hdg), &(opt.geoopt));
            }
            
        }
        if(i == n-1){
            traj[i].wpt.hdg = traj[i-1].wpt.hdg;
        }
        
    }


    Traj_Calc3D(traj, 0, &opt);
    Traj_Print(stdout, traj, 0);
    

    return EXIT_SUCCESS;
}
