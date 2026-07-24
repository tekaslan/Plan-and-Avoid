/**
 * @file
 * @brief Main program to test Vincenty formulas \cite vincenty1975direct. Basic 
 *        idea is to reproduce Table II of the reference. Small deviations are 
 *        observed but considered ok.
 * @author Pedro Donato
*/

#include "../include/geo.h"

/**
 * @brief Struct corresponding to one test case according to  
 *        \cite vincenty1975direct table II (page 92).
*/
struct GeoTest{
    struct Angle phi1;  /**< Point 1 Latitude */
    struct Angle phi2;  /**< Point 2 Latitude */
    struct Angle a1;    /**< Course from point 1 to 2 */ 
    struct Angle a2;    /**< Course from point 2 to 1 */
    struct Angle L;     /**< Longitude difference between points 1 and 2 */
    double s;           /**< Distance between poinst 1 and 2 */
};

/** 
 * @brief Computes the values of the test cases of \cite vincenty1975direct
 *        presented in Table II of the reference.
 * @param[in] geo Struct with inputs for the test.
 * @param[in] opt Geo option struct. 
*/
static void geo_test(const struct GeoTest * const geo, 
                     const struct GeoOpt * const opt);

/**
 * @brief Prints Angle struct in angle format. 
 * @param[in] ang Angle struct.
*/  
static void geo_printangle(const struct Angle * const ang);

/**
 * @brief Main program to test Vincenty formulas \cite vincenty1975direct. Basic 
 *        idea is to reproduce Table II of the reference. Small deviations are 
 *        observed but considered ok.
*/
int main(){

    struct GeoTest geo;
    struct GeoOpt opt = {BESSEL, 0};
   
    fprintf(stdout,"\t\t TEST FUNCTION FOR VINCENTY FUNCTIONS\n");
    
    /* Example A of Vincenty */
    fprintf(stdout,"\n\t\t\tEXAMPLE A\n");
    geo.phi1 = (struct Angle){55, 45, 0.0};
    geo.a1   = (struct Angle){96, 36, 8.7996};
    geo.s    = 14110526.170*M_2_NM;

    geo.phi2 = (struct Angle){-33, 26, 0.0};
    geo.L    = (struct Angle){108, 13, 0.0};
    geo.a2   = (struct Angle){137, 52, 22.01454};

    geo_test(&geo, &opt);

    /* Example B of Vincenty */
    fprintf(stdout,"\n\t\t\tEXAMPLE B\n");
    opt.model = INTERNATIONAL;
    geo.phi1 = (struct Angle){37, 19, 54.95367};
    geo.a1   = (struct Angle){95, 27, 59.63089};
    geo.s    = 4085966.703*M_2_NM;

    geo.phi2 = (struct Angle){26, 7, 42.83946};
    geo.L    = (struct Angle){41, 28, 35.50729};
    geo.a2   = (struct Angle){118, 5, 58.96161};

    geo_test(&geo, &opt);

    /* Example C of Vincenty */
    fprintf(stdout,"\n\t\t\tEXAMPLE C\n");
    opt.model = INTERNATIONAL;
    geo.phi1 = (struct Angle){35, 16, 11.24862};
    geo.a1   = (struct Angle){15, 44, 23.74850};
    geo.s    = 8084823.839*M_2_NM;

    geo.phi2 = (struct Angle){67, 22, 14.77638};
    geo.L    = (struct Angle){137, 47, 28.31435};
    geo.a2   = (struct Angle){144, 55, 39.92147};

    geo_test(&geo, &opt);

    /* Example D of Vincenty */
    fprintf(stdout,"\n\t\t\tEXAMPLE D\n");
    opt.model = INTERNATIONAL;
    geo.phi1 = (struct Angle){1, 0, 0.0};
    geo.a1   = (struct Angle){89, 0, 0.0};
    geo.s    = 19960000.000*M_2_NM;

    geo.phi2 = (struct Angle){0, -59, 53.83076};
    geo.L    = (struct Angle){179, 17, 48.02997};
    geo.a2   = (struct Angle){91, 0, 06.11733};

    geo_test(&geo, &opt);

    /* Example E of Vincenty */
    fprintf(stdout,"\n\t\t\tEXAMPLE E\n");
    opt.model = INTERNATIONAL;
    geo.phi1 = (struct Angle){1, 0, 0.0};
    geo.a1   = (struct Angle){4, 59, 59.99995};
    geo.s    = 19780006.558*M_2_NM;

    geo.phi2 = (struct Angle){1, 1, 15.18952};
    geo.L    = (struct Angle){179, 46, 17.84244};
    geo.a2   = (struct Angle){174, 59, 59.88481};

    geo_test(&geo, &opt);

    return EXIT_SUCCESS;
}


static void geo_test(const struct GeoTest * const geo, 
                     const struct GeoOpt * const opt)
{
    struct Pos A; 
    struct Pos B, Bcomp;
    double alpha1, alpha1comp;
    double alpha2, alpha2comp;
    double s;
    
    geo_fromang(&(geo->phi1),&(A.lat));
    geo_fromang(&(geo->a1),&alpha1);
    geo_npos(&A, &Bcomp, &(geo->s), &alpha1, opt);

    geo_fromang(&(geo->phi2),&(B.lat));
    geo_fromang(&(geo->L),&(B.lon));
    geo_fromang(&(geo->a2),&alpha2);
    geo_dist(&A, &B, &s, &alpha1comp, opt);

    /* First line */
    geo_printangle(&(geo->phi1));
    geo_printangle(&(geo->phi2));
    printf("%.1lf\t",1.0e5*fmod((Bcomp.lat-B.lat)*3600.0,60.0));
    printf("%.1lf\n",1.0e5*fmod((alpha1comp-alpha1)*3600.0,60.0));

    /* Second Line */
    geo_printangle(&(geo->a1));
    geo_printangle(&(geo->L));
    printf("%.1lf\t",1.0e5*fmod((Bcomp.lon-B.lon)*3600.0,60.0));
    geo_dist(&B, &A, &s, &alpha2comp, opt);
    alpha2comp -= 180.0;
    if(alpha2comp < 0) alpha2comp += 360.0;
    printf("%.1lf\n",1.0e5*fmod((alpha2comp-alpha2)*3600.0,60.0));
    
    /* Third Line */
    printf("%012.3lf     \t",geo->s*NM_2_M);
    geo_printangle(&(geo->a2));  
    geo_dist(&Bcomp, &A, &s, &alpha2comp, opt);
    alpha2comp -= 180.0;
    if(alpha2comp < 0) alpha2comp += 360.0;
    printf("%.1lf\t",1.0e5*fmod((alpha2comp-alpha2)*3600.0,60.0));
    geo_dist(&A, &B, &s, &alpha1comp, opt);
    printf("%.1lf\n",1.0e3*(s-geo->s)*NM_2_M);

}


static void geo_printangle(const struct Angle * const ang){

    fprintf(stdout, "%03"PRIdFAST16"o%02"PRIdFAST8"'%.5lf''\t",
                    ang->deg, ang->min, ang->sec);
}
