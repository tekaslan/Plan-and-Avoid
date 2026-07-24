/**
 * @file
 * @brief Sphere Earth Assumption, Haverine Formula
*/

/**
 * @brief Functions to Compute Sphere Coordinates in Haversine Formula.
 * @param[in] pos1 Point 1.
 * @param[in] pos2 Point 2 which will be computed.
 * @param[out] dist Distance between points 1 and 2.
 * @param[out] course Course between points 1 and 2.
 */
static void geo_distsphere(  const struct Pos* const pos1,
                        const struct Pos* const pos2,                                                            
                        double * const dist, double * const course)
{
    /* radius of Earth */
    double r = 6378137.0;
    double U1=pos1->lon * DEG_2_RAD;
    double P1=pos1->lat * DEG_2_RAD;
    double U2=pos2->lon * DEG_2_RAD;
    double P2=pos2->lat * DEG_2_RAD;
    
    /* Output distance in Nautical Miles */
    *dist = M_2_NM * 2 * r * asin(sqrt(pow(sin((U2 - U1) / 2) , 2.0) + cos(U1) * cos(U2) * pow(sin((P2 - P1) / 2) , 2.0)));

    /* Output heading in degrees */
    *course = RAD_2_DEG * atan2(sin(U2 - U1) * cos(P1), cos(P1) * sin(P2) - sin(P1) * cos(P2) * cos(U2 - U1));
}

/**
 * @brief Functions to Compute Sphere Coordinates in Haversine Formula.
 * @param[in] pos1 Point 1.
 * @param[out] pos2 Point 2 which will be computed.
 * @param[in] dist Distcne between points 1 and 2.
 * @param[in] course Course between points 1 and 2.
*/

static void geo_npossphere(  const struct Pos* const pos1,
                        struct Pos* const pos2,                                                            
                        const double * const dist, const double * const course)
{
    double r = 6378137.0;
    double theta = *dist/r;
    double U1=pos1->lon * DEG_2_RAD;
    double P1=pos1->lat * DEG_2_RAD;
    
    /* Latitude Output in Degrees */
    pos2->lat = RAD_2_DEG * asin(sin(P1) * cos(theta) + cos(P1) * sin(theta) * cos(*course * DEG_2_RAD));
    
    /* Longitude Output in Degrees and normalized */
    pos2->lon = RAD_2_DEG * (U1 + atan2(sin(*course * DEG_2_RAD) * sin(theta) * cos(P1), cos(theta)-sin(P1) * sin(pos2->lat * DEG_2_RAD)));
}