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
 * @brief Unit Conversions
 * @author Pedro Donato 
*/

// #ifndef UNIT_CONVERSIONS
// #define UNIT_CONVERSIONS

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/* ORIGINAL VALUES */
/** @brief Degrees to radians Ref. Table 6 of \cite nist2006si */
#define DEG_2_RAD (M_PI/180.0)

/** @brief Nautical miles to meters - Ref. Table 8 of \cite nist2006si
 *         and Table 3-3 of \cite icao2010ap5} */
#define NM_2_M (1852.0)

/** @brief Foot to meter - Ref. Table 3-3 of \cite icao2010ap5 */
#define FT_2_M (0.3048)

/** @brief Knots to m/s - Ref. Table 3-3 of \cite icao2010ap5 */
#define KT_2_MS (0.514444)

/** @brief Miles per hour to Knots*/
#define MPH_2_KT (1/1.151)

/** @brief Meters per second to miles per hour*/
#define MS_2_MPH (2.23694)

/** @brief Feet per second to knots*/
#define FTS_2_KTS (0.5925)


/* SIMPLE INVERSIONS */
/** @brief Meters to nautical miles - Simple inversion */ 
#define M_2_NM (1.0/NM_2_M)

/** @brief Meter to foot - Simple inversion */
#define M_2_FT (1.0/FT_2_M)

/** @brief Knots to m/s - Simple inversion */
#define MS_2_KT (1.0/KT_2_MS)

/** @brief Radians to degrees - Simple inversion */
#define RAD_2_DEG (1.0/DEG_2_RAD)

/** @brief Meters per second to miles per hour*/
#define MPH_2_MS (1/MS_2_MPH)

/** @brief Knots to feet per second*/
#define KTS_2_FTS (1/FTS_2_KTS)

/* COMBINATIONS */

/** @brief Foot to nautical mile - Deduced from above */
#define FT_2_NM  (FT_2_M*M_2_NM)

/** @brief Nautical mile to foot - Deduced from above */
#define NM_2_FT  (NM_2_M*M_2_FT)

/** @brief Feet to kilometer - Deduced from above */
#define FT_2_KM (FT_2_M/1000.0)

// #endif
