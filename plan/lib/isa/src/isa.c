/**
 * @file
 * @brief Implements the US Standard Atmosphere 1976.
 * 
 */

#include "../include/isa.h"

static const double r0 = 6.356766e3;
static const double g0 = 9.80665;
static const double T0 = 288.15;
static const double P0 = 1.013250e5;
static const double R = 8.31432e3;
static const double M0 = 28.9644;

static const double Hb[7]  = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0};
static const double Lmb[7] = {-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0};

/** 
 * @details Simply implements equation 18 from \cite NASA1976
 * with the initial converion from ft to km from the input.
 * Also uses the value of r0 from page 8 of the same reference.
 *
 * \f[
 *    H = \frac{r_0Z}{r_0+Z} = \frac{6.356766.10^6Z}{6.356766.10^6 + Z} 
 * \f]
 * 
 * Note: important to note that the value of r0 of table 2 is wrong,
 * it is written 6.356766e6 km and should be m.
 */
double ISA_H(const double * const z){
	return (*z)*FT_2_KM;
	/*return r0*(*z)*FT_2_KM / (r0 + (*z)*FT_2_KM); */
}

double ISA_Z(const double * const h){
    return r0 *(*h)*FT_2_KM / (r0 - (*h)*FT_2_KM);

}

/**
 * @details Implements Equation 23 from \cite NASA1976
 * using the data from table 4 of the same reference. Observe that according to
 * the same reference below 86km the molecular temperature (\f$T_M\f$) is equal to the
 * kinetic temperature (\f$T\f$).  
 * \f[
 * 		T = T_M = T_{M,b}+L_{M,b}\left(H-H_b\right)
 * \f]
 * where \f$T_{M,0}=288.15K\f$.
 *
 * The table below presents the values of the parameters used. Only the first
 * three layers are modeled.
 *
 * <CENTER> 
 * | Layer ID | Geopotential Height | Temperature Gradient |
 * |:--:|:--:|:--:|
 * | \f$b\f$ | \f$H_b\f$ (km') | \f$L_{M,b}\f$ (K/km') |
 * | 0 | 0 | -6.5 |
 * | 1 | 11 | 0.0 |
 * | 2 | 20 | 1.0 |
 * </CENTER>	
 *
 * The value of \f$H\f$ is obtained calling ISA_H
 */
double ISA_T(const double * const alt){

    double Tmb = T0;
    int b = 0;

    for(; b < 6; b++){
        if(ISA_H(alt) > Hb[b+1]){
            Tmb += Lmb[b] * (Hb[b+1] - Hb[b]);
        }
        else{
            return Tmb + Lmb[b] * (ISA_H(alt) - Hb[b]);
        }
    }

    fprintf(stderr,"ISA - ERROR too high altitude! \n");
    return 0.0;
}

/**
 * @details Implements equations 33a and 33b of \cite NASA1976US:
 *
 * \f[
 *   P =
 *   P_b\left[\frac{T_{M,b}}{T_{M,b}+L_{M,b}\left(H-H_b\right)}\right]^{\frac{g_0'M_0}{R^*L_{M,b}}} 
 * \f]
 * \f[
 * 	 P = P_be^{\frac{-g_0'M_0\left(H-H_b\right)}{R^*L_{M,b}}}
 * \f]
 * with the constants given in the same reference (most of them in table 2):
 * 
 * \li \f$P_0 = 1.013250.10^5N/m^2\f$
 * \li \f$g_0' = 9.80665m^2/(s^2m')\f$
 * \li \f$M_0 = 28.96644kg/kmol\f$ - This value is presented in page 9 of the \cite NASA1976 
 * \li \f$R^* = 8.31432.10^{3}N.m/(kmol.K)\f$ - The value is table 2 has the
 * wrong exponent. The correct value is presented in page 3 of \cite NASA1976  
 * 
 * The value of \f$H\f$ is obtained calling ISA_H
 */
double ISA_P(const double * const alt)
{

    double Pb = P0;
    double Tmb = T0;
    double aux;
    int b = 0;

    for(; b < 6; b++){

        if(ISA_H(alt) > Hb[b+1]){
            Tmb += Lmb[b] * (Hb[b+1] - Hb[b]);
            if(Lmb[b] != 0.0){
                Pb *= pow((Tmb / (Tmb + Lmb[b]*(Hb[b+1] - Hb[b]))),
                          (g0*M0 / (R*Lmb[b]/1000.0)));
            }else{
                Pb *= exp(-g0*M0*(Hb[b+1] - Hb[b]) / (R * Tmb));
            }
        }else{
             if(Lmb[b] != 0.0){
                return Pb * pow((Tmb / (Tmb + Lmb[b]*(ISA_H(alt) - Hb[b]))),
                                (g0*M0 / (R*Lmb[b]/1000.0)));
            }else{
                aux = Pb * exp(-g0*M0*(ISA_H(alt) - Hb[b]) / (R * Tmb));
                return aux;
            }
        }
    }
            
    fprintf(stderr,"ISA - ERROR too high altitude! \n");
    return 0.0;
}

/**
 * @details Implements equation 42 of \cite NASA1976:
 * \f[
 *    \rho = \frac{PM_0}{R^*T_M}
 * \f]
 * where:
 * \li \f$ P \f$ is computed using ISA_P
 * \li \f$ T_M\f$ is computed using ISA_M
 * \li \f$M_0 = 28.96644kg/kmol\f$ - This value is presented in page 9 of the
 * \cite NASA1976 
 * \li \f$R^* = 8.31432.10^{3}N.m/(kmol.K)\f$ - The value is table 2 has the
 * wrong exponent. The correct value is presented in page 3 of \cite NASA1976  
 */
double ISA_Rho(const double * const alt){
	return ISA_P(alt)*M0/(R*ISA_T(alt));
}
