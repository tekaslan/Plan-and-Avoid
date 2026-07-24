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
 * @brief Struc ts used in Geo Package
 * @author Pedro Donato
*/

#ifndef GEO_STRUCTS
#define GEO_STRUCTS

#include <stdint.h>
#include <stdbool.h>


struct GeoOpt;
struct Pos;
struct PosECEF;
struct Angle;

/**
 * @brief Earth Model Enumerator. Used for the Geodesy Functions and also for 
 *        the Pos struct.
*/
enum Earth {
    WGS84,   /**< WGS-84 ellipsoid model. Used by the GPS system. */
	BESSEL,  /**< Bessel ellipsoid model. Used to test Vincenty functions. */
    /** International ellipsoid model. Used to test Vincenty functions. */
	INTERNATIONAL,
    FLAT,     /**< Flat Earth. Use cartesian coordinates */
    SPHERE
};

/**
 * @brief Options for Geo package
*/
struct GeoOpt{
    /** Earth model (ellipsoid or flat Earth)*/
    enum Earth model;
    /** Verbose mode */
    bool  verbose;
};

#ifndef POSITION_STRUCT
#define POSITION_STRUCT
/**
 *  @struct Pos
 *  @brief Contains the 3 dim state vec (lat,lon,alt).
 * 
 *  @details 
 *  ## General Info##
 *  Three states are the coordinates in Geodetic Reference Frame (WGS-84)
 *  or the x,y,z coordinate in the flat Earth model. This is defined by the
 *  value used with the Earth enumerator
 *
 *  ##Legacy Code##
 *  Old AFP had a State struct very similar, but more limited since it would
 *  only work with degrees and assume North and West coordinates 
 */
struct Pos {
    double lat; /**<Geodetic latitude in deg OR y pos in NM in flat Earth model*/
    double lon; /**<Geodetic longitude in deg OR x pos in NM in flat Earth model*/
    double alt; /**<Geodetic altitude in feet OR z pos in ft in flat Earth Model*/
    double hdg;
    double t; // Timestamp [s]
};
#endif

struct PosXYZ
{
    double X;   // ECEF x-position [NM]
    double Y;   // ECEF y-position [NM]
    double Z;   // ECEF z-position [NM]
};

/** 
 * @brief Struct to manage lat/lon in angles.
*/
struct Angle{
    int_fast16_t deg; /**< Angle degrees */
    int_fast8_t min;  /**< Angle minutes */
    double sec;       /**< Angle seconds */
};

#endif

