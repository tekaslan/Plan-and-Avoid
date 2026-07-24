/**
 * @file
 * @brief Check functions for geo inputs.
*/

/**
 * @brief Check is lat and lon values are valid.
 * @param[in] pos Position struct to be checked.
 * @param[in] opt GeoOpt struct to get the model.
 * @retval TRUE Invalid position.
 * @retval FALSE Valid position.
*/
static bool geo_checklatlon(const struct Pos* const pos, 
                            const struct GeoOpt * const opt)
{
   if(opt->model != FLAT){
        if((pos->lat > 90.0) || (pos->lat < -90.0) || 
           (pos->lon > 180.0) || (pos->lon < -180.0)){
            // fprintf(stderr,"Geo ERROR - Invalid Position. \n");
            return 1;
        }
    }

    return 0;
}
