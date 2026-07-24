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

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                                   %
#    Mathematical and Geometrical Utility Functions                                 %
#    Airspace and Ground-risk Aware                                                 %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : H. Emre Tekaslan (tekaslan@vt.edu)                                   %
#    Date    : April 2025                                                           %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "math_utils.h"


/*
    Returns size of an array
*/
int arr_length(double arr[]) {
    int i;
    int count = 0;
    for(i=0; arr[i]!='\0'; i++)
    {
        count++;
    }
    return count;
}

/*
    Flips surface vertices left to right
*/
void fliplr(struct surface *surface) {

    // Temporary array
    struct surface tmp_surface;

    // Assign values to temp array
    int j = 0;
    for (int i = 4; i > 0; i-- && j++)
    {
        tmp_surface.xyz[i-1] = surface->xyz[j];
    }

    // Assign temp array to the original array
    *surface = tmp_surface;
}

/*
    Finds the minimum value in an array 
*/  
double minimum(double *arr) { 
    // Get the size of array
    int size = arr_length(arr);

    // Assume the first element is the minimum 
    double min = arr[0]; 
  
    // Loop through the array to find the minimum 
    for (int i = 1; i < size; i++) { 
        if (arr[i] < min) { 
            // Update minimum if a smaller element is found 
            min = arr[i]; 
        } 
    } 
  
    return min; 
} 

/*
    Returns dot product of given ECEF position vectors
*/
double dot(struct PosXYZ *pos1, struct PosXYZ *pos2) {
    return (pos1->X*pos2->X + pos1->Y*pos2->Y + pos1->Z*pos2->Z);
}


/*
    Calculates 2-norm of a given vector
*/
double norm(struct PosXYZ *vec) {   
    // Vector's inner product with itself yields its norm squared.
    double n = sqrt(dot(vec,vec));

    return n;
}


/*
    Cross products of ECEF position vectors pos1 and pos2 
    written in ECEF position vector cross_prod
*/
void cross(struct PosXYZ *pos1, struct PosXYZ *pos2, struct PosXYZ *cross_prod) {

    // Cross product
    cross_prod->X = (pos1->Y*pos2->Z - pos1->Z*pos2->Y);
    cross_prod->Y = (pos1->Z*pos2->X - pos1->X*pos2->Z);
    cross_prod->Z = (pos1->X*pos2->Y - pos1->Y*pos2->X);

}

/*
    Normalizes a given vector
*/
void normalize(struct PosXYZ *v){
    double magnitude = norm(v);
    if (fabs(magnitude) < 1e-10) {
        // Optional: Warn about near-zero length vector.
        return;
    }
    v->X /= magnitude;
    v->Y /= magnitude;
    v->Z /= magnitude;
}

/*
    Creates a vector between given two points and assigns to the third input
*/
void createVector(struct PosXYZ *p1, struct PosXYZ *p2, struct PosXYZ *p3){
    // Create the vector
    p3->X = p1->X - p2->X;
    p3->Y = p1->Y - p2->Y;
    p3->Z = p1->Z - p2->Z;
}

/*
    Calculates the angle between two vectors
*/
double angle_between_vectors(struct PosXYZ *vec1, struct PosXYZ *vec2) {
    double ang = acos(dot(vec1, vec2)/(norm(vec1)*norm(vec2)));

    return ang;
}

/*
    Calculates the centroid coordinates of a given surface
*/
void surfaceCentroid(struct surface *surface){   
    // Define intermediate variables
    double X = 0, Y = 0, Z = 0;

    // Sum all coordinates
    for (int i = 0; i < 4; i++)
    {
        X = X + surface->xyz[i].X;
        Y = Y + surface->xyz[i].Y;
        Z = Z + surface->xyz[i].Z;
    }

    // Assign centroid to the surface
    surface->centroid.X = X/4;
    surface->centroid.Y = Y/4;
    surface->centroid.Z = Z/4;
}

/*
    Calculates the unit normal vector of a given surface
*/
void surfaceNormal(struct surface *surface){

    if (surface == NULL) {
        printf("surface is NULL!\n");
        return;
    }

    // Pre-allocate edge vectors
    struct PosXYZ v1;
    struct PosXYZ v2;

    // Define edge vectors
    v1.X = surface->xyz[1].X - surface->xyz[0].X;
    v1.Y = surface->xyz[1].Y - surface->xyz[0].Y;
    v1.Z = surface->xyz[1].Z - surface->xyz[0].Z;

    v2.X = surface->xyz[3].X - surface->xyz[0].X;
    v2.Y = surface->xyz[3].Y - surface->xyz[0].Y;
    v2.Z = surface->xyz[3].Z - surface->xyz[0].Z;

    // Calculate the cross product of two edge vectors which is the surface normal
    cross(&v1, &v2, &surface->normal);

    // Normalize
    double norm_val = sqrt(surface->normal.X * surface->normal.X +
        surface->normal.Y * surface->normal.Y +
        surface->normal.Z * surface->normal.Z);

    if (fabs(norm_val) < 1e-10) {
        fprintf(stderr, "Warning: Degenerate surface encountered; normal set to zero.\n");
        surface->normal.X = surface->normal.Y = surface->normal.Z = 0.0;
    } else {
        surface->normal.X /= norm_val;
        surface->normal.Y /= norm_val;
        surface->normal.Z /= norm_val;
    }
}

/*
    Calculates the coefficient of a surface equation 'k' in the form of ax + by + cz + k = 0.
*/
void surfaceEquation(struct surface *surface) {
    // Allocate memory
    surface->k = (double *) malloc(sizeof(double));

    // Calculate the surface equation constant
    *surface->k = -dot(&surface->normal, &surface->xyz[0]);
}

/*
    Calculates the distance between a point and a surface
*/
double distance_point2surface(struct PosXYZ *pos, struct surface *surface) {   
    // Get the surface equation constant
    surfaceEquation(surface);

    // Point-to-surface distance
    double d = fabs(dot(pos, &surface->normal) + *surface->k);

    return d;
}

/*
    Calculates the distance between a point and an edge defined between two points
*/
double distance_point2edge(struct PosXYZ *p, struct PosXYZ *edge_point1, struct PosXYZ *edge_point2) {   
    // Get the edge vector
    struct PosXYZ *v = (struct PosXYZ *) malloc(sizeof(struct PosXYZ));
    createVector(edge_point2, edge_point1, v);

    // Get point-to-vertex vector
    struct PosXYZ *p2v = (struct PosXYZ *) malloc(sizeof(struct PosXYZ));
    createVector(p, edge_point1, p2v);

    // Get the cross product of the edge vector and the point-to-vertex vector
    struct PosXYZ cross_prod_vec;
    cross(v, p2v, &cross_prod_vec);

    // Calculate the distance
    double d = norm(&cross_prod_vec)/norm(v);

    free(v); v = NULL;
    free(p2v); p2v = NULL;

    return d;
}

/*
    Calculates the distance from a point to all vertices of a surfaces and returns the minimum distance
*/
double distance_point2vertex(struct PosXYZ *pos, struct surface *surface) {
    // Pre-allocate distance array
    double d_new = 999;
    double d_min = 999;

    // Calculates distances
    struct PosXYZ v;
    for (int i = 0; i < 4; i++)
    {
        createVector(pos, &surface->xyz[i], &v);
        d_new = norm(&v);
        if (d_new < d_min)
        {
            d_min = d_new;
        }
    }
    return d_min;
}

/*
    Wraps the given angle to 360
*/
double wrapTo360(double angle) {
    while(angle >= 360.0){
        angle -= 360.0;
    }
    while(angle < 0.0){
        angle += 360.0;
    }
    return angle;
}

/*
    Wraps the given angle to 180
*/
double wrapTo180(double angle) {
    // Wrap the angle to the range [-180, 180)
    angle = fmod(angle + 180.0, 360.0);
    if (angle < 0) angle += 360.0; // Ensure positive result from fmod()
    return angle - 180.0;
}

/*
    Returns an upper bound to the discrete search space size
*/
int size(struct Pos *init, struct Pos *goal, double dOpt, double segmentLength, double dAltitude, double dPsi) {
    // Vertical discretization size
    int Nv = ceil((init->alt - goal->alt)/dAltitude);
    double dh = (init->alt - goal->alt)/(Nv-1);

    // Number of hexagons
    int Nhex = 0;

    // Initialize
    double alt = goal->alt;
    double circumradius;
    int N;
    int Nsum;
    for (int i = 0; i < Nv; i++) {
        if (alt > 2000) {
            circumradius = 0.5*(segmentLength + 0.1*(alt - 2000));
        } else {
            circumradius = 0.5*segmentLength;
        }

        // Number of layers of hexagons within the footprint
        N = (int) ceil((2*dOpt/sqrt(3) - 0.5*sqrt(3)*circumradius)/(sqrt(3)*circumradius));

        Nsum = 0;
        for (int j = 1; j < N; j++) {
            Nsum = Nsum + 6*j;
        }

        // Total number of hexagons
        Nhex += 1 + Nsum;

        // Update altitude
        alt += dh;

    }

    // Upper bound of the number states
    int NS = Nhex * 360/dPsi;

    return NS;
}

/*
    Returns the shortest path distance from a state to the other
    considering the optimal Dubins length
*/
double minDistance(struct Pos *state1, struct Pos *state2, double turnRadius, struct GeoOpt *GeoOpt, SearchProblem * problem)
{

    // Get great circle distance between states
    double bearing;
    double distance;
    geo_dist(state1, state2, &distance, &bearing, GeoOpt);

    // // Add turning distance
    // double dPsi = fabs(DEG_2_RAD*(wrapTo180(bearing - state1->hdg)));
    // double dTurn = turnRadius * dPsi;

    // double dmin = distance*NM_2_FT + dTurn;     // Overestimate to the minimum distance

    if (distance*NM_2_FT < 3*problem->ac->turnRadius){
        return distance*NM_2_FT;
    } else {
        
        // Initialize Dubins structure
        struct DubinsPath dubins;
        Traj_InitArray(dubins.traj, 4);
        dubins.traj[0].wpt.pos.lat = state1->lat;
        dubins.traj[0].wpt.pos.lon = state1->lon;
        dubins.traj[0].wpt.pos.alt = state1->alt;
        dubins.traj[0].wpt.hdg = state1->hdg;

        dubins.traj[3].wpt.pos.lat = state2->lat;
        dubins.traj[3].wpt.pos.lon = state2->lon;
        dubins.traj[3].wpt.pos.alt = state2->alt;
        dubins.traj[3].wpt.hdg = state2->hdg;

        dubins.traj[0].wpt.gam = -problem->ac->gammaOptTurn;
        dubins.traj[1].wpt.gam = -problem->gammaOpt;
        dubins.traj[2].wpt.gam = -problem->ac->gammaOptTurn;

        dubins.traj[0].wpt.rad = problem->ac->turnRadius*FT_2_NM;
        dubins.traj[1].wpt.rad = 0;
        dubins.traj[2].wpt.rad = problem->ac->turnRadius*FT_2_NM;

        shortestDubins(&dubins, problem->DubinsOpt);

        double dmin = dubins.hdist*NM_2_FT;
        // double dmin = dubins.hdist*NM_2_FT - problem->identRadius;

        return dmin;
    }
}

/*
    Returns the minimum of two given values
*/
double min(double a, double b) {
    return (a <= b) ? a : b;
}

/*
    Returns the maximum of two given values
*/
double max(double a, double b) {
    return (a >= b) ? a : b;
}

/*
    Calculates the unit normal vector
    for a given state
*/
double *normalVector(struct Pos *state, struct GeoOpt *GeoOpt) {
    struct Pos *tmpPos = (struct Pos *) malloc(sizeof(struct Pos));
    double const d = 3000*FT_2_NM;

    geo_npos(state, tmpPos, &d, &state->hdg, GeoOpt);
    tmpPos->alt = state->alt;

    struct PosXYZ tmp_ned;
    geo_lla2ned(state, tmpPos, &tmp_ned);

    double *normal = (double *) malloc(2*sizeof(normal));
    normal[1] = tmp_ned.X;
    normal[0] = tmp_ned.Y;

    // Normalize to a unit vector
    double magnitude = sqrt(pow(normal[0],2) + pow(normal[1], 2));
    normal[0] /= magnitude;
    normal[1] /= magnitude;

    free(tmpPos); tmpPos = NULL;

    return normal;
}

// Function to compute cumulative 3D distances along the path
void cumulativeDistances(SearchProblem * problem, double *distances) {
    double dist, course, horizontal_distance, vertical_distance;
    distances[0] = 0.0;
    for (int i = 1; i < problem->path->depth; i++) {
        geo_dist(&problem->path->nodes[i-1].state, &problem->path->nodes[i].state, &dist, &course, &problem->GeoOpt);
        horizontal_distance = dist*NM_2_FT;
        vertical_distance = problem->path->nodes[i-1].state.alt - problem->path->nodes[i].state.alt;
        distances[i] = distances[i - 1] + sqrt(horizontal_distance * horizontal_distance + vertical_distance * vertical_distance);
    }

    geo_dist(&problem->path->nodes[problem->path->depth-1].state, &problem->Goal, &dist, &course, &problem->GeoOpt);
    horizontal_distance = dist*NM_2_FT;
    vertical_distance = problem->path->nodes[problem->path->depth-1].state.alt - problem->Goal.alt;
    distances[problem->path->depth] = distances[problem->path->depth - 1] + sqrt(horizontal_distance * horizontal_distance + vertical_distance * vertical_distance);

}

// Linear interpolation function
double linearInterp1d(double x0, double x1, double y0, double y1, double x) {
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

// Function to sample equidistant points along the path
void sampleEquidistantPoints(SearchProblem * problem, double *distances, double interval,
                            double **sampled_lat, double **sampled_lon, double **sampled_alt,
                            int *num_samples) {

    int n = problem->path->depth;
    // double max_distance = distances[n - 1];
    double max_distance = distances[n];
    double dist2sN      = distances[n-1];
    *num_samples = (int) (max_distance / interval) + 1;
    int num_samples2sN = (int) (dist2sN / interval) + 1;

    *sampled_lat = (double *) malloc((*num_samples) * sizeof(double));
    *sampled_lon = (double *) malloc((*num_samples) * sizeof(double));
    *sampled_alt = (double *) malloc((*num_samples) * sizeof(double));

    for (int i = 0; i < *num_samples; i++) {
        double target_dist = i * interval;

        // Find the segment where the target distance falls
        int j = 0;
        while (j < n - 1 && distances[j + 1] < target_dist) {
            j++;
        }

        if (j < n - 1) {
            (*sampled_lat)[i] = linearInterp1d(distances[j], distances[j + 1], problem->path->nodes[j].state.lat, problem->path->nodes[j+1].state.lat, target_dist);
            (*sampled_lon)[i] = linearInterp1d(distances[j], distances[j + 1], problem->path->nodes[j].state.lon, problem->path->nodes[j+1].state.lon, target_dist);
            (*sampled_alt)[i] = linearInterp1d(distances[j], distances[j + 1], problem->path->nodes[j].state.alt, problem->path->nodes[j+1].state.alt, target_dist);
        } else if (j == n - 1){
            (*sampled_lat)[i] = linearInterp1d(distances[j], distances[j + 1], problem->path->nodes[j].state.lat, problem->Goal.lat, target_dist);
            (*sampled_lon)[i] = linearInterp1d(distances[j], distances[j + 1], problem->path->nodes[j].state.lon, problem->Goal.lon, target_dist);
            (*sampled_alt)[i] = linearInterp1d(distances[j], distances[j + 1], problem->path->nodes[j].state.alt, problem->Goal.alt, target_dist);
        } else {
            (*sampled_lat)[i] = problem->Goal.lat;
            (*sampled_lon)[i] = problem->Goal.lon;
            (*sampled_alt)[i] = problem->Goal.alt;
        }
    }
}

// Returns a random float between min and max
double random_double(double min, double max) {
    return min + ((double) rand() / (double) RAND_MAX) * (max - min);
}
