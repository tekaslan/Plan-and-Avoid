/**
 * @file
 */

#ifndef ISAHDR
#define ISAHDR

#include "conversions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** 
 * @brief Calculates the geopotential altitude from geometric altitude according to US Atmosphere 1976
 * @param[in]	*z	Geometric altitude in feet
 * @return 	Geodetic altitude in kilometers
 */
double ISA_H(const double *z);

/**
 * @brief Calculates the temperature in K associatedd with a altitude (until 32km)
 * @param[in] alt 	Altitude in feet
 * @todo Consider the temperature variation on that day
 */	
double ISA_T(const double *alt);

/**
	@brief Calculates the pressure in Pa associatedd with a altitude (until 32km)
	@param[in] alt 	Altitude in feet
	@return Pressure in Pa
 */
double ISA_P(const double *alt);

/**
	@brief Calculates the density in km/m^3 associatedd with a altitude (until 32km)
	@param[in] alt 	Altitude in feet
	@return Air density in km/m^3
 */
double ISA_Rho(const double *alt);
#endif
