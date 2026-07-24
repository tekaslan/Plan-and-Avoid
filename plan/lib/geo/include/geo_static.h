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
 * @brief Static functions used by geo.c
 * @author Pedro Donato (pdonato@umich.edu) 
*/

/**
 * @brief Functions to Compute Inverse Geodetics per \cite vincenty1975direct .
*/
static void geo_distvincenty(  const struct Pos* const pos1, 
                        const struct Pos* const pos2,                                                            
                        double * const dist, double * const course,
                        const struct GeoOpt * const opt);

/**
 * @brief Functions to Compute Direct Geodetics per \cite vincenty1975direct.
*/
static void geo_nposvincenty(  const struct Pos* const pos1, 
                        struct Pos* const pos2,                                                            
                        const double * const dist, const double * const course,
                        const struct GeoOpt * const opt);

/**
 * @brief Functions to Compute Inverse Geodetics assuming flat earth.
*/
static void geo_distflat(const struct Pos* const pos1, 
                         const struct Pos* const pos2,                                                            
                         double * const dist, double * const course);

/**
 * @brief Functions to Compute Direct Geodetics assuming flat earth.
*/
static void geo_nposflat(  const struct Pos* const pos1, 
                           struct Pos* const pos2,                                                            
                           const double * const dist, 
                           const double * const course);

/**
 * @brief Check is lat and lon values are valid.
*/
static bool geo_checklatlon(const struct Pos* const pos,
                            const struct GeoOpt * const opt);
