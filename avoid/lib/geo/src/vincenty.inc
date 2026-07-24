/**
 * @file
 * @brief Functions to Compute Direct and Inverse Geodetics as per 
 *        \cite vincenty1975direct
*/

#ifndef VINCENTY_INC
#define VINCENTY_INC

#include "../include/conversions.h"

/**
 * @brief Earth Ellipsoid Model Struct 
*/
struct Geo_Model{
    double a; /**< Major semi-axis of ellipsoid */
    double f; /**< Flattening */
    double b; /**< Minor semi-axis of ellipsoid */
};

/**
 * @brief Earth Ellipsoid Models Ref. \cite soler1989important
*/
static const struct Geo_Model models[3] = {
    {6378137.0, 1.0/298.257223563, (1.0-1.0/298.257223563)*6378137.0}, 
    {6377397.155, 1.0/299.1528128, (1.0-1.0/299.1528128)*6377397.155}, 
    {6378388.0, 1.0/297.0, (1.0-1.0/297.0)*6378388.0}};

/**
 * @brief Tolerance for Vincenty Method iterations. 
 *        Ref. \cite vincenty1975direct 
*/
static const double e = 1e-12;

//change the comments (in-out)
/**
 * @brief Functions to Compute Inverse Geodetics per \cite vincenty1975direct .
 * @param[in] pos1 Point 1.
 * @param[in] pos2 Point 2 which will be computed.
 * @param[out] dist Distance between points 1 and 2.
 * @param[out] course Course between points 1 and 2.
 * @param[in] opt GeoOpt struct to get the ellipsoid model.
*/
static void geo_distvincenty(  const struct Pos* const pos1, 
                        const struct Pos* const pos2,                                                            
                        double * const dist, double * const course,
                        const struct GeoOpt * const opt )
{

    /* Variables for calculation follows Vincenty */
    double U1 = atan((1.0-models[opt->model].f)*tan(pos1->lat*DEG_2_RAD));
    double U2 = atan((1.0-models[opt->model].f)*tan(pos2->lat*DEG_2_RAD));
    double L = (pos2->lon - pos1->lon)*DEG_2_RAD;

    double lambda = L;
    double lambda_old = lambda - 10.0;

    double cos2a = 0.0, sinsig = 0.0, cossigm2 = 0.0, cossig = 0.0, sig = 0.0;

    while(fabs(lambda - lambda_old) > e){

        lambda_old = lambda;

        sinsig = sqrt((cos(U2)*sin(lambda))*(cos(U2)*sin(lambda)) + 
                (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda))*
                (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda)));

        cossig = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
        sig = atan2(sinsig,cossig);
        double sina = cos(U1)*cos(U2)*sin(lambda)/sinsig;
        cos2a = 1.0 - sina*sina;
        cossigm2 = cossig - 2.0*sin(U1)*sin(U2)/cos2a;

        double C = models[opt->model].f/16.0*cos2a*
                   (4.0+models[opt->model].f*(4.0-3.0*cos2a));

        lambda = L + (1.0-C)*models[opt->model].f*sina*
            (sig + C*sinsig*(cossigm2+C*cossig*(-1.0+2.0*cossigm2*cossigm2)));
    }

    double u2 = cos2a*(models[opt->model].a*models[opt->model].a-
                models[opt->model].b*models[opt->model].b)/
                (models[opt->model].b*models[opt->model].b);

    double A = 1.0+u2/16384.0*(4096.0+u2*(-768.0+u2*(320.0-175.0*u2)));    
    double B = u2/1024.0*(256.0+u2*(-128.0+u2*(74.0-47.0*u2)));

    double dsig = B*sinsig*(cossigm2+
                            0.250*B*(cossig*(-1.0 + 2.0*cossigm2*cossigm2)
                            -1.0/6.0*B*cossigm2*(-3+4*sinsig*sinsig)*
                            (-3+4*cossigm2*cossigm2)));

    /* Output distance in Nautical Miles */
    *dist = models[opt->model].b*A*(sig-dsig) * M_2_NM;

    /* Output heading in degrees */
    *course = RAD_2_DEG*atan2(cos(U2)*sin(lambda),
            (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda)));

}

/**
 * @brief Functions to Compute Direct Geodetics per \cite vincenty1975direct.
 * @param[in] pos1 Point 1.
 * @param[out] pos2 Point 2 which will be computed.
 * @param[in] dist Distcne between points 1 and 2.
 * @param[in] course Couser between points 1 and 2.
 * @param[in] opt GeoOpt struct to get the ellipsoid model.
*/

static void geo_nposvincenty(  const struct Pos* const pos1, 
                        struct Pos* const pos2,                                                            
                        const double * const dist, const double * const course,
                        const struct GeoOpt * const opt)
{

    double a1 = *course*DEG_2_RAD;
    double s = (*dist)*NM_2_M;

    double tanU1 = (1.0-models[opt->model].f)*tan(pos1->lat*DEG_2_RAD);
    double cosU1 = 1.0/sqrt(1.0 + (tanU1*tanU1));
    double sinU1 = tanU1*cosU1;
    double sig1 = atan2(tanU1,cos(a1));
    double sina = cosU1*sin(a1);
    double cos2a = 1.0 - sina*sina;
    double u2 = cos2a*(models[opt->model].a*models[opt->model].a
                -models[opt->model].b*models[opt->model].b)/
                (models[opt->model].b*models[opt->model].b);

    double A = 1.0+u2/16384.0*(4096.0+u2*(-768.0+u2*(320.0-175.0*u2)));    
    double B = u2/1024.0*(256.0+u2*(-128.0+u2*(74.0-47.0*u2)));

    double sig = s/(models[opt->model].b*A);
    double sigo = 2.0*M_PI;
    double sigm2 = 0.0;
    while(fabs(sig-sigo) > e){
        sigm2 = 2.0*sig1+sig;
        double dsig = B*sin(sig)*(cos(sigm2)
                    +0.250*B*(cos(sig)*(-1.0 + 2.0*cos(sigm2)*cos(sigm2))
                    -1.0/6.0*B*cos(sigm2)*(-3.0+4.0*sin(sig)*sin(sig))*
                    (-3.0+4.0*cos(sigm2)*cos(sigm2))));
        sigo = sig;
        sig = s/(models[opt->model].b*A) + dsig;
    }

    /* Latitude Output in Degrees */
    pos2->lat = RAD_2_DEG*atan2((sinU1*cos(sig) + cosU1*sin(sig)*cos(a1)),
                (1-models[opt->model].f)*
                sqrt(sina*sina+(sinU1*sin(sig)-cosU1*cos(sig)*cos(a1))*
                (sinU1*sin(sig)-cosU1*cos(sig)*cos(a1))));

    double lambda = atan2(sin(sig)*sin(a1),
                          (cosU1*cos(sig) - sinU1*sin(sig)*cos(a1)));

    double C = models[opt->model].f/
               16.0*cos2a*(4.0+models[opt->model].f*(4.0-3.0*cos2a));

    double L = lambda 
               -(1.0-C)*models[opt->model].f*sina*(sig+C*sin(sig)*(cos(sigm2) 
               + C*cos(sig)*(-1.0+2.0*cos(sigm2)*cos(sigm2))));

    /* Longitude Output in Degrees and normalized */
    pos2->lon = pos1->lon+L*RAD_2_DEG;

    while(pos2->lon > 180.0){
        pos2->lon -= 360.0;
    }
    
    while(pos2->lon < -180.0){
        pos2->lon += 360.0;
    }

}

#endif