/**
 * @file 
 * @brief Functions for Flat Earth Assumption
*/

/**
 * @brief Functions to Compute Inverse Geodetics assuming flat earth.
 * @param[in] pos1 Point 1.
 * @param[in] pos2 Point 2 which will be computed.
 * @param[out] dist Distance between points 1 and 2.
 * @param[in] course Course between points 1 and 2.
*/
static void geo_distflat(const struct Pos* const pos1, 
                         const struct Pos* const pos2,                                                            
                         double * const dist, double * const course)
{
    *dist = sqrt(pow((pos1->lat-pos2->lat),2.0) +
                 pow((pos1->lon-pos2->lon),2.0));

    *course = RAD_2_DEG*atan2(pos2->lon - pos1->lon, pos2->lat - pos1->lat);

}

/**
 * @brief Functions to Compute Direct Geodetics assuming flat earth.
 * @param[in] pos1 Point 1.
 * @param[out] pos2 Point 2 which will be computed.
 * @param[in] dist Distance between points 1 and 2.
 * @param[in] course Course between points 1 and 2.
*/
static void geo_nposflat(  const struct Pos* const pos1, 
                           struct Pos* const pos2,                                                            
                           const double * const dist,
                           const double * const course)
{
    pos2->lat = pos1->lat + (*dist)*cos(*course*DEG_2_RAD);    
    pos2->lon = pos1->lon + (*dist)*sin(*course*DEG_2_RAD);    
}