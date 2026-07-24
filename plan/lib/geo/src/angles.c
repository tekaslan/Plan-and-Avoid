/**
 * @file
 * @brief Functions to deal with Angle struct.
 * @author Pedro Donato
*/

#include "../include/geo.h"

/**
 * @brief Tolerance used to avoid values of 60 on seconds.
*/
static const double tol = 1e-12;

void geo_toang(struct Angle * const angd, const double * const ang)
{
    angd->deg = (int_fast16_t) *ang;
    if(angd->deg){
        angd->min = (int_fast8_t) fabs(fmod(*ang * 60.0, 60.0));
        angd->sec = fabs(fmod(*ang * 3600.0, 60.0));
    }
    else{
        angd->min = (int_fast8_t) fmod(*ang * 60.0, 60.0);
        if(angd->min)
            angd->sec = fabs(fmod(*ang * 3600.0, 60.0));
        else
            angd->sec = fmod(*ang * 3600.0, 60.0);
    }
        

    if(fabs(angd->sec - 60.0) < tol){
        angd->min++;
        angd->sec -= 60.0;
        if(angd->sec < 0) angd->sec = 0.0;
    }
        
}

void geo_fromang(const struct Angle * const angd, double * const ang)
{
    if(angd->deg > 0){
        *ang = angd->deg + ((double) angd->min)/60.0 + angd->sec/3600.0;
    }else if(angd->deg < 0){ 
        *ang = angd->deg - ((double) angd->min)/60.0 - angd->sec/3600.0;
    }else if(angd->min > 0){
        *ang = ((double) angd->min)/60.0 + angd->sec/3600.0;
    }else if(angd->min < 0){
        *ang = ((double) angd->min)/60.0 - angd->sec/3600.0;
    }else{
        *ang = angd->sec;
    }

}
