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

#ifndef MATH_HEADER
#define MATH_HEADER

#include <stdlib.h>
#include <math.h>
#include "geo.h"
#include "search.h"

struct surface;
typedef struct SearchProblem SearchProblem;


/*
    Returns size of an array
*/
int arr_length(double arr[]);

/*
    Flips surface vertices left to right
*/
void fliplr(struct surface *surface);

/*
    Finds the minimum value in an array 
*/  
double minimum(double *arr);

/*
    Returns dot product of given ECEF position vectors
*/
double dot(struct PosXYZ *pos1, struct PosXYZ *pos2);

/*
    Calculates 2-norm of a given vector
*/
double norm(struct PosXYZ *vec);

/*
    Cross products of ECEF position vectors pos1 and pos2 
    written in ECEF position vector cross_prod
*/
void cross(struct PosXYZ *pos1, struct PosXYZ *pos2, struct PosXYZ *cross_prod);

/*
    Normalizes a given vector
*/
void normalize(struct PosXYZ *v);

/*
    Creates a vector between given two points and assigns to the third input
*/
void createVector(struct PosXYZ *p1, struct PosXYZ *p2, struct PosXYZ *p3);

/*
    Calculates the angle between two vectors
*/
double angle_between_vectors(struct PosXYZ *vec1, struct PosXYZ *vec2);

/*
    Calculates the centroid coordinates of a given surface
*/
void surfaceCentroid(struct surface *surface);

/*
    Calculates the unit normal vector of a given surface
*/
void surfaceNormal(struct surface *surface);

/*
    Calculates the coefficient of a surface equation 'k' in the form of ax + by + cz + k = 0.
*/
void surfaceEquation(struct surface *surface);

/*
    Calculates the distance between a point and a surface
*/
double distance_point2surface(struct PosXYZ *pos, struct surface *surface);

/*
    Calculates the distance between a point and an edge defined with two points
*/
double distance_point2edge(struct PosXYZ *p, struct PosXYZ *edge_point1, struct PosXYZ *edge_point2);

/*
    Calculates the distance from a point to all vertices of a surfaces and returns the minimum distance
*/
double distance_point2vertex(struct PosXYZ *pos, struct surface *surface);

/*
    Wraps the given angle to 360 and 180 degrees
*/
double wrapTo360(double angle);
double wrapTo180(double angle);

/*
    Returns an upper bound to the discrete search space size
*/
int size(struct Pos *init, struct Pos *goal, double dOpt, double segmentLength, double dAltitude, double dPsi);


double minDistance(struct Pos *state1, struct Pos *state2, double turnRadius, struct GeoOpt *GeoOpt, SearchProblem * problem);
double min(double a, double b);
double max(double a, double b);
double *normalVector(struct Pos *state, struct GeoOpt *GeoOpt);
void cumulativeDistances(SearchProblem * problem, double *distances);
double linearInterp1d(double x0, double x1, double y0, double y1, double x);
void sampleEquidistantPoints(SearchProblem * problem, double *distances, double interval,
                            double **sampled_lat, double **sampled_lon, double **sampled_alt,
                            int *num_samples);

/*
    Returns a random float between min and max
*/
double random_double(double min, double max);

#endif
