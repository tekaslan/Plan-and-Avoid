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
#ifndef TRAJ_HEADER
#define TRAJ_HEADER

#include "../../geo/include/geo.h"
#include <string.h>

/** Tolerance for checking zero values */
#define ZERO_TOLERANCE 1.0E-6

struct Waypoint;
struct Traj;
struct TrajOpt;

/** Waypoint definition */
struct Waypoint{
    /**< Position struct */
    struct Pos pos;     
    /**< Heading from this waypoint to the next [deg] */
    double hdg;         
    /**< Trajectory radius [nm] from this waypoint to the next, 0 if straight */
    double rad;         
    /**< Flight path angle from this waypoint to the next [deg] */
    double gam;         
    /**< Bank angle from this waypoint to the next [deg] */
    double phi;         
    /**< True airspeed from this point to the next [deg] */
    double V;           
    /**< CL */
    double CL;           
};

struct Traj{
    /**< Waypoint */
    struct Waypoint wpt;    
    /**< Heading changing from this waypoint to the next */
    double dpsi;            
    /**< Distance in the horizontal plane from this waypoint to the next */
    double hdist;
    /**< Number of complete 360deg turns to be executed from this waypoint to the next */
    int nturn;
    /**< Pointer to next waypoint */
    struct Traj * next;
};

struct TrajOpt{
    int verbose;
    struct GeoOpt geoopt;
};


int Traj_Copy(struct Traj *p, struct Traj *c);
int Traj_CopyAll(struct Traj *p, struct Traj *c);
int traj_calctraj_angdist(struct Traj* p, int npoints, const struct TrajOpt *opt);
double Traj_HDist(struct Traj* p, int npoints, const struct TrajOpt *opt);
int Traj_Calc3D(struct Traj* p, int npoints, const struct TrajOpt *opt);
int Traj_Print(FILE *output, struct Traj *p, int npoints);
int Traj_InitArray(struct Traj *p, int npoints);
struct Pos* Traj_ToPlot(struct Traj* traj, int *npoints, struct TrajOpt *opt);
double Traj_LayerTime(struct Traj* traj, double layer_min, double layer_max, struct TrajOpt *opt);


/**
 * @brief Computes the angular distance between two angles according 
 *        to the desired direction.
 * @param[in] ang1 - Initial angle [deg] - Between 0 and 360.
 * @param[in] ang2 - Final angle [deg] - Between 0 and 360.
 * @param[in] direction - Positive to CW, Negative for CCW.
 * @param[out] angdist - Angular distance [deg] - Between -360 and 360.
 *                       Negative if CCW.
 * @retval 0 Success.
 * @retval 1 Fail. Direction is equal to zero.
*/
int traj_calcangdist(const double * const ang1, const double * const ang2, 
                     const double * const direction, 
                      double * const angdist);


int traj_semiplan(const double * const ang1, const double * const ang2);


void traj_file(struct Traj* traj, struct TrajOpt *opt);

#endif
