/**
 * @brief Example function to use Geo library.
 * @author Pedro Donato 
*/

#include "../include/geo.h"

/**
 * @brief Main function to show how the library works. 
 *
 * @details Show example of computing distances and courses; other point given
 * distance and course; and example of wrong inputs. It will output an ERROR.
*/
int main()
{
    
    /* Initial Position Struct */
    struct Pos A = {0};

    /* Final Position Struct */
    struct Pos B = {0};
    
    /* Distance variable */
    double dist = 0.0;

    /* Course variable */
    double course = 0.0;

    /* Option struct - use WGS84 which is the used for GPS*/
    struct GeoOpt opt = {0};
    opt.model = WGS84;
    opt.verbose = 0;  

    /* Position A is FXB */
    A.lat = 42.293496;
    A.lon = -83.711954; 

    /* Position B is Big House */
    B.lat = 42.265836;
    B.lon = -83.748663;

    /* Compute distance and heading */
    geo_dist(&A, &B, &dist, &course, &opt);

    /* Print Output */
    fprintf(stdout, "\n#   Geo Library Example \n#\n");
    fprintf(stdout, "# Distance from FXB to Big House : %lf [m]\n",
                    dist * NM_2_M); 
    fprintf(stdout, "# Course from FXB to Big House : %lf [deg]\n",
                    course); 

    /* Compute opposite course */
    geo_dist(&B, &A, &dist, &course, &opt);
    fprintf(stdout, "# Course from Big House to FXB : %lf [deg]\n#\n",
                    course); 

    /* Compute a new position */
    dist = 1000.0 * M_2_NM;
    course = 45.0;
    geo_npos(&A, &B, &dist, &course, &opt);
    fprintf(stdout, "# New latitude 1000m North-East from FXB: %lf [deg]\n",
                    B.lat);
    fprintf(stdout, "# New longitude 1000m North-East from FXB: %lf [deg]\n#\n",
                    B.lon);

    /* Check for same position */
    geo_dist(&A, &A, &dist, &course, &opt);
    fprintf(stdout, "# Distance from FXB to FXB : %lf [m]\n",
                    dist * NM_2_M); 
    fprintf(stdout, "# Course from FXB to FXB : %lf [deg]\n\n",
                    course); 

    /* Check for invalid latitude */
    B.lat = 95.0;
    geo_dist(&A, &B, &dist, &course, &opt);
    fprintf(stdout, "# Distance from FXB to Invalid Point : %lf [m]\n",
                    dist * NM_2_M); 
    fprintf(stdout, "# Course from FXB to Invalid Point : %lf [deg]\n",
                    course); 

    return EXIT_SUCCESS;
}
