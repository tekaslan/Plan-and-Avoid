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

#ifndef CENSUS_H
#define CENSUS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <shapefil.h>
#include "../include/rtree_wrapper.h"

struct RTreeHandle;

typedef struct Census {
    SHPHandle shapefile;      // Handle for shapefile
    DBFHandle dbf;            // Handle for DBF file
    int pop_field_index;      // Index of the population field
    int area_field_index;     // Index of the area field
    RTreeHandle *rtree;       // Handle for the R-tree
    SHPObject **polygons;     // Array storing polygon objects
    int numPolygons;          // Number of polygons stored
} Census;

// Initializes Census by reading a given shapefile
struct Census *initCensus(const char *shapefile_path);

// Returns the population density at a given coordinate (linear search)
double getPopulationDensityLinearSearch(struct Census *census, double lat, double lon, double max_density);
bool isPointInPolygon(double px, double py, double *x_coords, double *y_coords, int num_vertices);

// Frees memory allocated for Census structure
void freeCensus(struct Census *census);

// Returns the bounding box of the Census dataset
double *censusBoundingBox(struct Census *census);

// Finds the polygon that contains the given (x, y) coordinate
// SHPObject *findPolygon(struct Census *census, double x, double y);
int findPolygon(struct Census* census, double x, double y);

// Checks whether a point is inside a given polygon
int pointInPolygon(SHPObject* polygon, double qx, double qy);

double getPopulationDensityRTree(struct Census *census, double lat, double lon, double max_density);

#endif // CENSUS_H
