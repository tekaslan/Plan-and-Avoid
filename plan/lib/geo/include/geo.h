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

/**
 * @file
 *
 * @brief Geo Main Header File
 *
 * @author Pedro Donato & H. Emre Tekaslan
 * 
 * @version June 3, 2024
 */

#ifndef GEO_HEADER
#define GEO_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <inttypes.h>
#include <proj.h>

#include "conversions.h"
#include "geo_structs.h"

#ifndef VPRINT
/** Print statements for verbose mode */
#define VPRINT(...) if(opt->verbose) fprintf(stdout,__VA_ARGS__)
#endif


// CONSTANTS
// WGS84
#define earth_semimajor (6378137)
#define earth_flattening (1 / 298.2572235630)
#define earth_semiminor (earth_semimajor*(1 - earth_flattening))

// Earth's Polar Radius [m]
// #define earth_semiminor (6356752)

// Earth's Eccentricity
#define earth_eccentricity (0.0818197909921145)

/** 
 * @brief Gives the 2D distance and heading between two positions with 
 *        different models
 * @param[in] pos1 - From position
 * @param[in] pos2 - To position
 * @param[out] dist - Distance between two points [NM]
 * @param[in] course - Course from pos1 to pos2 [deg]
 * @param[in] opt - Geo options, used for get Earth Model.
 * @retval 0 Success.
 * @retval 1 Invalid Position 1.
 * @retval 2 Invalid Position 2.
*/
uint_fast8_t geo_dist(const struct Pos* const pos1,const struct Pos* const pos2,
             double * const dist, double * const course, 
             const struct GeoOpt * const opt);

/** 
 * @brief Calculates the next 2D position given a initial position, distance and 
 *        course usign different models.
 * @param[in] pos1 - From position
 * @param[out] pos2 - To position
 * @param[in] dist - Distance between two points [NM]
 * @param[in] course - Course from pos1 to pos2 [deg]
 * @param[in] opt - Geo options, used for get Earth Model.
 * @retval 0 Success
 * @retval 1 Invalid position 1. 
*/
uint_fast8_t geo_npos(const struct Pos* const pos1, struct Pos* const pos2, 
              const double * const dist, const double * const course, 
              const struct GeoOpt *opt);
/**
 * @brief Normalize heading to between 0 and 360 as aviation standard.
 * @param[in,out] hdg Heading angle to be normalized.
*/
void geo_normhdg(double * const hdg);

/** 
 * @brief Calculates the next 2D position given a initial position, turn radius 
 *        initial heading and heading change usign different models.
 * @param[in] pos1 - From position
 * @param[out] pos2 - To position
 * @param[in] rad - Turn radius [NM]
 * @param[in] p1hdg - Initial heading [deg]
 * @param[in] dpsi - Heading change [deg]
 * @param[in] opt - Geo options, used for get Earth Model.
 * @retval 0 Success
 * @retval 1 Invalid position 1. 
*/
uint_fast8_t geo_npost( const struct Pos* const pos1, struct Pos* const pos2, 
                        const double *rad, const double *p1hdg, 
                        const double *dpsi, 
                        const struct GeoOpt * const opt);

/** 
 * @brief Calculates the turn radius between two points given their position
 *        and first point heading. 
 *        
 * @param[in] pos1 - From position
 * @param[out] pos2 - To position
 * @param[in] rad - Turn radius [NM]
 * @param[in] p1hdg - Initial heading [deg]
 * @param[in] opt - Geo options, used for get Earth Model.
 * @retval 0 Success
 * @retval 1 Invalid position 1. 
 * @retval 2 Invalid position 2. 
*/
uint_fast8_t geo_distt( const struct Pos* const pos1, 
                        const struct Pos* const pos2, 
                        double * const rad, 
                        const double * const p1hdg,
                        const struct GeoOpt * const opt);

/**
 * @brief Convert an angle from double to Angle struct.
 * @param[out] angd Angle struct.
 * @param[in] ang Angle in double format.
*/
void geo_toang(struct Angle * const angd, const double * const ang);

/**
 * @brief Convert an angle from double to Angle struct.
 * @param[in] angd Angle struct.
 * @param[out] ang Angle in double format.
*/
void geo_fromang(const struct Angle * const angd, double * const ang);

/**
 * @brief Calculates the next 3D position in a straight flight.
 * @param[in] pos1 - From position
 * @param[out] pos2 - To position
 * @param[in] dist - Distance between two points [NM]
 * @param[in] hdg - Heading from pos1 to pos2 [deg]
 * @param[in] gam - Flight path angle from pos1 to pos2 [deg]
 * @param[in] opt - Geo options, used for get Earth Model.
 * @return EXIT_SUCCESS always
*/
uint_fast8_t geo_npos3D(const struct Pos* pos1, struct Pos* pos2,
                const double *dist, const double *hdg, const double *gam, 
                const struct GeoOpt * opt);

/**
 * @brief Calculates Earth's curvature for a given latitude.
 * @param[in] pos1 - Latitude, longitude, altitude coordinate
 * @return N - Earth's curvature in meters
*/
double geo_earths_curvature(struct Pos *pos1);

/**
 * @brief Converts LLA to ECEF coordinate
 * @param[in] pos1 - LLA coordinate
 * @param[out] pos2 - ECEF coordinate
*/
void geo_lla2ecef(struct Pos *pos1, struct PosXYZ *pos2);

/**
 * @brief Converts LLA to NED coordinate
 * @param[in] pos1 - LLA coordinate
 * @param[out] pos2 - NED coordinate
*/
void geo_lla2ned(struct Pos *ref_lla, struct Pos *pos_lla, struct PosXYZ *pos_ned);

#endif
