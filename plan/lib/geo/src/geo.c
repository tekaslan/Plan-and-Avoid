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
#    Geodesic Calculations                                                          %
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Modifications  : Pedro Di Donato & H. Emre Tekaslan (tekaslan@vt.edu)          %
#    Date    : June 2024                                                            %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#                      https://scholar.google.com/citations?user=UCxHXTgAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "geo.h"
#include "geo_static.h"

/* PRIVATE FUNCTIONS */
#include "vincenty.inc"
#include "flat.inc"
#include "geo_check.inc"
#include "sphere.inc"

uint_fast8_t 
geo_dist(const struct Pos* const pos1, const struct Pos* const pos2, 
             double * const dist, double * const hdg, 
             const struct GeoOpt * const opt){

    /* Check for inputs inside range */
    if(geo_checklatlon(pos1,opt)) return 1;
    if(geo_checklatlon(pos2,opt)) return 2;

    /* Check for same inputs */
    if((pos1->lat == pos2->lat) && (pos1->lon == pos2->lon)){
        *dist = 0.0;
        *hdg = 0.0;
        return 0;
    }
    
    /* Vincenty inverse formula implementation */
    if(opt->model != FLAT){
        if (opt->model == SPHERE){
            geo_distsphere(pos1, pos2, dist, hdg);
        }
        else{
            geo_distvincenty(pos1, pos2, dist, hdg, opt);
        }
    }else{
        geo_distflat(pos1, pos2, dist, hdg);
    }

    /* Normalize Heading from 0 to 360 */
    geo_normhdg(hdg);

    return 0;
}

uint_fast8_t geo_npos(const struct Pos* const pos1, struct Pos* const pos2,
              const double * const dist, const double * const hdg,
              const struct GeoOpt * opt)
{
    /* Check for inputs inside range */
    if(geo_checklatlon(pos1,opt)) return 1;

    /* Check for zero distance */
    if(*dist == 0.0){
        pos2->lat = pos1->lat;
        pos2->lon = pos1->lon;
        return 0;
    }
 
    /* Compute new position */
    if(opt->model != FLAT){
        if (opt->model == SPHERE){
            geo_npossphere(pos1, pos2, dist, hdg);
        }
        else{
            geo_nposvincenty(pos1, pos2, dist, hdg, opt);
        }
 
    }else{                                                                                                                           
        geo_nposflat(pos1, pos2, dist, hdg);                                                                                         
    }

    return 0;
}


void geo_normhdg(double *const hdg){

    while(*hdg >= 360.0){
        *hdg -= 360.0;
    }
    while(*hdg < 0.0){
        *hdg += 360.0;
    }

}

uint_fast8_t geo_npost( const struct Pos* const pos1, struct Pos* const pos2, 
                        const double *rad, const double *p1hdg, 
                        const double *dpsi, 
                        const struct GeoOpt * const opt)
{

    double hdg = *p1hdg + *dpsi/2.0;
    double dist = fabs(*rad) * sqrt(2-2*cos((*dpsi)*DEG_2_RAD));
    
    /* Check for inputs inside range */
    if(geo_checklatlon(pos1,opt)) return 1;

    /* Check for zero distance */
    if(*dpsi == 0.0){
        pos2->lat = pos1->lat;
        pos2->lon = pos1->lon;
        return 0;
    }

    /* Compute new position */ 
    geo_npos(pos1,pos2,&dist,&hdg,opt);

    return 0;
}

uint_fast8_t geo_distt( const struct Pos* const pos1, 
                        const struct Pos* const pos2, 
                        double * const rad,
                        const double * const p1hdg,
                        const struct GeoOpt * const opt)
{
    
    double dist = 0.0;
    double psi_rel = 0.0;
    double psi_p = 0.0;

    /* Check for inputs inside range */
    if(geo_checklatlon(pos1,opt)) return 1;
    if(geo_checklatlon(pos2,opt)) return 2;

    /* Check for same inputs */
    if((pos1->lat == pos2->lat) && (pos1->lon == pos2->lon)){
        *rad = 0.0;
        return 0;
    }
 
    /* Compute Distance */ 
    geo_dist(pos1, pos2, &dist, &psi_rel, opt);
    psi_p = 90.0 + *p1hdg - psi_rel;
    *rad = dist / (2.0 * cos(psi_p * DEG_2_RAD));

    return 0;
}

/**
	##General Info##
	Calls the function Geo_Npos() and calculates the new height using the 
	flight path angle \f$ \gamma\f$ value.

	New heading is equal to the old heading.

	##Legacy Code Comments##
	There was not a dedicated function for this in the legacy code. It was
	implemented inside other functions using the same idea.

	##Modification History##
	### July 31, 2014###
	\li Function name changed from WGS84_Npos_3D to Geo_Npos_3D
	\li Inclusion of the Earth model flag
*/
uint_fast8_t geo_npos3D(const struct Pos* pos1, struct Pos* pos2,
                const double *dist, const double *hdg, const double *gam, 
                const struct GeoOpt * opt)
{

	if(geo_npos(pos1,pos2,dist,hdg,opt)) return 1;
	// pos2->alt = pos1->alt + (*dist)*NM_2_FT*tan((*gam)*DEG_2_RAD);
    pos2->alt = pos1->alt + (*dist)*tan((*gam)*DEG_2_RAD);

	return 0;
}

/*
    Computes Earth's curvature in meters at given latitude
*/
double geo_earths_curvature(struct Pos *pos1)
{   
    double a = earth_semimajor;
    double b = earth_semiminor;
    double lat_radians = pos1->lat * DEG_2_RAD;
    double numerator = a * a;
    double denominator = sqrt(a * a * cos(lat_radians) * cos(lat_radians) + b * b * sin(lat_radians) * sin(lat_radians));

    return numerator / denominator;
}

/*
    Converts latitude, longitude, altitude coordinates to 
    Earth-Centered-Earth-Fixed coordinate frame.
    pos1: LLA coordinates
    pos2: ECEF coordinates in meters
*/
void geo_lla2ecef(struct Pos *pos1, struct PosXYZ *pos2)
{
    // Get Earth's curvature
    double N = geo_earths_curvature(pos1);

    // Convert coordinates
    pos2->X = (N + pos1->alt*FT_2_M)*cos(pos1->lat*DEG_2_RAD)*cos(pos1->lon*DEG_2_RAD);
    pos2->Y = (N + pos1->alt*FT_2_M)*cos(pos1->lat*DEG_2_RAD)*sin(pos1->lon*DEG_2_RAD);
    pos2->Z = (N * pow((earth_semiminor / earth_semimajor),2) + pos1->alt*FT_2_M) * sin(pos1->lat*DEG_2_RAD);
}

/*
    Converts LLA coordinates to North-East-Down system
    ref_lla: Reference LLA coordinates
    pos_lla: Query LLA coordinates
*/
void geo_lla2ned(struct Pos *ref_lla, struct Pos *pos_lla, struct PosXYZ *pos_ned)
{   

    // Find ECEF coordinates of the reference
    struct PosXYZ ref_ecef;
    struct PosXYZ pos_ecef;
    geo_lla2ecef(ref_lla, &ref_ecef);
    geo_lla2ecef(pos_lla, &pos_ecef);

    double sin_phi = sin(ref_lla->lat * DEG_2_RAD);
    double cos_phi = cos(ref_lla->lat * DEG_2_RAD);
    double sin_lambda = sin(ref_lla->lon * DEG_2_RAD);
    double cos_lambda = cos(ref_lla->lon * DEG_2_RAD);

    // Rotation Matrix for ECEF to ENU transformation
    double R[3][3] = {{-sin_lambda, cos_lambda, 0},   // North
                      {-sin_phi*cos_lambda, -sin_phi*sin_lambda, cos_phi},    // East
                      {cos_phi * cos_lambda, cos_phi * sin_lambda, sin_phi}}; // Down

    // Difference vector
    double v[3] = {(pos_ecef.X - ref_ecef.X),
                   (pos_ecef.Y - ref_ecef.Y),
                   (pos_ecef.Z - ref_ecef.Z)};

    // double t = cos_lambda * v[0] + sin_lambda * v[1];
    // double East = -sin_lambda * v[0] + cos_lambda * v[1];
    // double Up = cos_phi * t + sin_phi * v[2];
    // double North = -sin_phi * t + cos_phi * v[2];

    pos_ned->Y = (R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2]); // North
    pos_ned->X = (R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2]); // East
    pos_ned->Z = (R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]); // Down
}
