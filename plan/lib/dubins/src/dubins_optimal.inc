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
#    Dubins Path Handling Functions                                                 %
#                                                                                   %
#    Author  : Pedro Di Donato                					    %
#    Date    : 2017                                                                 %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=UCxHXTgAAAAJ&hl=en %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

static void dubins_evalandcheck(struct DubinsPath * dub, 
                                struct DubinsPath * best,
                                enum DubinsType type,
                                struct DubinsOpt *opt)
{
    struct DubinsPath new;
    Traj_InitArray(new.traj,4);
    Traj_CopyAll(dub->traj, new.traj);
    new.type = type;

    if(Dubins(&new, opt)) return;

    if (new.hdist < best->hdist){ 
        Traj_CopyAll(new.traj, best->traj);
        best->hdist = new.hdist;
        Dubins(&new, opt);
    }

    return;
}
    
static bool dubins_short(const struct Traj * const traj,
                        const double * const course,
                        const double * const dist)
{

    const double alpha = traj[0].wpt.hdg - *course;
    const double beta = traj[3].wpt.hdg - *course;
    double rad = traj[0].wpt.rad;

    if(traj[2].wpt.rad > rad)   rad = traj[2].wpt.rad;

    if(*dist/rad < sqrt(4 - (fabs(cos(alpha)) + fabs(cos(beta))) *
                           (fabs(cos(alpha)) + fabs(cos(beta)))) + 
                           fabs(sin(alpha)) + fabs(sin(beta))){
        return true;
    }
    else{
        return false;
    }
}

static uint_fast8_t dubins_quadrant(const double * const hdg, 
                                    const double * const course)
{

    const double aux = *hdg - *course;
    if(((aux >= -90.0) && (aux <= 0))|| (aux >= 270.0))     return 1u;
    else if((aux >= -180) || (aux >= 180.0))                return 2u;
    else if((aux >= -270) || (aux >= 90.0))                 return 3u;
    else                                                    return 4u;

}
  
static bool dubins_optimal(struct DubinsPath * const dub,
                          struct DubinsOpt * opt)
{

    double course = 0.0;
    double dist = 0.0;
    geo_dist(&(dub->traj[0].wpt.pos),&(dub->traj[3].wpt.pos),&dist,&course,
             &(opt->trajopt.geoopt));

    const uint_fast8_t alpha = dubins_quadrant(&(dub->traj[0].wpt.hdg), 
                                               &course);
    const uint_fast8_t beta  = dubins_quadrant(&(dub->traj[3].wpt.hdg),
                                               &course);

    struct DubinsPath best;
    Traj_InitArray(best.traj,4);  
    best.hdist = 1e99;
    
    switch(alpha){
    case 1u:
        switch(beta){
        case 1u:  
            dubins_evalandcheck(dub,&best,RSL,opt); 
            break;
        case 2u:
            dubins_evalandcheck(dub,&best,RSR,opt); 
            dubins_evalandcheck(dub,&best,RSL,opt); 
            break;
        case 3u:   
            dubins_evalandcheck(dub,&best,RSR,opt); 
            dubins_evalandcheck(dub,&best,LSR,opt);
            break;
        case 4u:   
            dubins_evalandcheck(dub,&best,RSR,opt); 
            dubins_evalandcheck(dub,&best,RSL,opt); 
            dubins_evalandcheck(dub,&best,LSR,opt); 
            break;
        default:   return true;
        } 
        break;
    case 2u:
        switch(beta){
        case 1u:  
            dubins_evalandcheck(dub,&best,LSL,opt); 
            dubins_evalandcheck(dub,&best,RSL,opt); 
            break;
        case 2u:
            dubins_evalandcheck(dub,&best,LSL,opt); 
            dubins_evalandcheck(dub,&best,RSL,opt); 
            dubins_evalandcheck(dub,&best,RSR,opt); 
            break;
        case 3u:   
            dubins_evalandcheck(dub,&best,RSR,opt); 
            break;
        case 4u:   
            dubins_evalandcheck(dub,&best,RSR,opt); 
            dubins_evalandcheck(dub,&best,RSL,opt); 
            break;
        default:   return true;
        } 
        break;
    case 3u:
        switch(beta){
        case 1u:  
            dubins_evalandcheck(dub,&best,LSL,opt); 
            dubins_evalandcheck(dub,&best,LSR,opt); 
            break;
        case 2u:
            dubins_evalandcheck(dub,&best,LSL,opt);
            break;
        case 3u:   
            dubins_evalandcheck(dub,&best,RSR,opt); 
            dubins_evalandcheck(dub,&best,LSR,opt); 
            dubins_evalandcheck(dub,&best,LSL,opt); 
            break;
        case 4u:   
            dubins_evalandcheck(dub,&best,RSR,opt); 
            dubins_evalandcheck(dub,&best,LSR,opt); 
            break;
        default:   return true;
        } 
        break;
    case 4u:
        switch(beta){
        case 1u:  
            dubins_evalandcheck(dub,&best,RSL,opt); 
            dubins_evalandcheck(dub,&best,LSR,opt); 
            dubins_evalandcheck(dub,&best,LSL,opt); 
            break;
        case 2u:
            dubins_evalandcheck(dub,&best,LSL,opt); 
            dubins_evalandcheck(dub,&best,RSL,opt); 
            break;
        case 3u:   
            dubins_evalandcheck(dub,&best,LSR,opt); 
            dubins_evalandcheck(dub,&best,LSL,opt); 
            break;
        case 4u:   
            dubins_evalandcheck(dub,&best,LSR,opt); 
            break;
        default:   return true;
        } 
        break;
    default: return true;
    }

    if(dubins_short(dub->traj, &course, &dist)){
        dubins_evalandcheck(dub,&best,LRL1,opt); 
        dubins_evalandcheck(dub,&best,LRL1,opt); 
        dubins_evalandcheck(dub,&best,RLR2,opt); 
        dubins_evalandcheck(dub,&best,RLR2,opt); 
    }

    best.hdist = Traj_HDist(best.traj,0,&(opt->trajopt));

    Traj_CopyAll(best.traj, dub->traj);
    dub->hdist = best.hdist;

    return false;
 
}
