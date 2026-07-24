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
#    Census Data Handling Functions                                                 %
#    Airspace, Air Traffic and Ground-risk Aware                                    %
#    Aircraft Contingency Landing Planner                                           %
#    Using Gradient-guided 4D Discrete Search                                       %
#    and 3D Dubins Solver                                                           %
#                                                                                   %
#    Autonomous Aerospace Systems Laboratory (A2Sys)                                %
#    Kevin T. Crofton Aerospace and Ocean Engineering Department                    %
#                                                                                   %
#    Author  : H. Emre Tekaslan (tekaslan@vt.edu)                                   %
#    Date    : July 2026                                                            %
#                                                                                   %
#    Google Scholar  : https://scholar.google.com/citations?user=uKn-WSIAAAAJ&hl=en %
#    LinkedIn        : https://www.linkedin.com/in/tekaslan/                        %
#                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "../include/census.h"
#include "../include/rtree_wrapper.h"

/*
    Initializes census structure
    Reads the given shape file
    Creates R*-Tree
*/
struct Census* initCensus(const char * shapefile_path) {

    struct Census * census = (struct Census *) malloc(sizeof(struct Census));
    if (!census) {
        perror("Failed to allocate memory for Census");
        exit(EXIT_FAILURE);
    }

    census->shapefile = SHPOpen(shapefile_path, "rb");
    if (!census->shapefile) {
        fprintf(stderr, "Error: Unable to open shapefile at %s\n", shapefile_path);
        free(census);
        exit(EXIT_FAILURE);
    }

    census->dbf = DBFOpen(shapefile_path, "rb");
    if (!census->dbf) {
        fprintf(stderr, "Error: Unable to open DBF file at %s\n", shapefile_path);
        SHPClose(census->shapefile);
        free(census);
        exit(EXIT_FAILURE);
    }

    // Get field indices
    census->pop_field_index = DBFGetFieldIndex(census->dbf, "POP20");
    census->area_field_index = DBFGetFieldIndex(census->dbf, "ALAND20");
    if (census->pop_field_index == -1 || census->area_field_index == -1) {
        fprintf(stderr, "Error: Required fields not found in DBF file.\n");
        DBFClose(census->dbf);
        SHPClose(census->shapefile);
        free(census);
        exit(EXIT_FAILURE);
    }

    rtree_init(); // Initialize the R-tree

    int numEntities;
    SHPGetInfo(census->shapefile, &numEntities, NULL, NULL, NULL);
    census->numPolygons = numEntities;
    census->polygons = (SHPObject**)malloc(numEntities * sizeof(SHPObject*));
    for (int i = 0; i < numEntities; i++) {
        SHPObject* shape = SHPReadObject(census->shapefile, i);
        if (!shape) continue;
        census->polygons[i] = shape;

        // Insert polygon bounding box into clR-tree
        rtree_insert(shape->dfXMin, shape->dfYMin, shape->dfXMax, shape->dfYMax, i);
    }

    return census;
}

/*
    Ray-casting (Point-in-Polygon)
*/
bool isPointInPolygon(double px, double py, double *x_coords, double *y_coords, int num_vertices)
{
    bool inside = false;

    for (int i = 0, j = num_vertices - 1; i < num_vertices; j = i++)
    {
        // Get the vertices of the polygon edge
        double xi = x_coords[i], yi = y_coords[i];
        double xj = x_coords[j], yj = y_coords[j];

        // Checks if a right ray-cast from the point intersects with the edge
        bool intersect = ((yi > py) != (yj > py)) &&
                         (px < (xj - xi) * (py - yi) / (yj - yi) + xi);

        // Change the flag everytime horizontal ray intesects with an edge
        if (intersect)
            inside = !inside;
    }
    return inside;
}

/*
    Returns population density
*/
double getPopulationDensityLinearSearch(struct Census *census, double lat, double lon, double max_density)
{
    SHPObject *shape = NULL;
    double density = 0.0;
    int num_shapes, shape_type;
    double min_bound[4], max_bound[4];

    // Get shapefile info
    SHPGetInfo(census->shapefile, &num_shapes, &shape_type, min_bound, max_bound);

    for (int entity = 0; entity < num_shapes; entity++) {
        shape = SHPReadObject(census->shapefile, entity);
        if (!shape) continue;

        // Use the Ray-Casting algorithm to check if the point is inside the polygon
        if (isPointInPolygon(lon, lat, shape->padfX, shape->padfY, shape->nVertices))
        {
            double population = DBFReadDoubleAttribute(census->dbf, entity, census->pop_field_index);
            double area = DBFReadDoubleAttribute(census->dbf, entity, census->area_field_index);

            if (area > 0) {
                density = population / area;
            }

            SHPDestroyObject(shape);
            break;
        }

        SHPDestroyObject(shape);
    }

    // Replace outliers and cap density
    if (isnan(density) || density < 0) {
        return 0.0;
    } else {
        return fmin(density, max_density);
    }
}

/*
    Frees up the memory
*/
void freeCensus(struct Census *census)
{
    if (!census) return;

    if (census->dbf) {
        DBFClose(census->dbf);
        census->dbf = NULL;
    }

    if (census->shapefile) {
        SHPClose(census->shapefile);
        census->shapefile = NULL;
    }

    for (int i = 0; i < census->numPolygons; i++) {
        if (census->polygons[i]) {
            SHPDestroyObject(census->polygons[i]);
            }
        }
    free(census->polygons); census->polygons = NULL;

    rtree_destroy();

    free(census);
}

// Ray-casting (Point-in-Polygon)
int pointInPolygon(SHPObject* polygon, double qx, double qy) {
    int intersections = 0;
    for (int i = 0; i < polygon->nParts; i++) {
        int startIdx = polygon->panPartStart[i];
        int endIdx = (i == polygon->nParts - 1) ? polygon->nVertices : polygon->panPartStart[i + 1];

        for (int j = startIdx, k = endIdx - 1; j < endIdx; k = j++) {
            double x1 = polygon->padfX[j], y1 = polygon->padfY[j];
            double x2 = polygon->padfX[k], y2 = polygon->padfY[k];

            if (((y1 > qy) != (y2 > qy)) &&
                (qx < (x2 - x1) * (qy - y1) / (y2 - y1) + x1)) {
                intersections++;
            }
        }
    }
    return (intersections % 2) == 1;
}

// Find Polygon Containing the Query Point
int findPolygon(struct Census* census, double x, double y) {
    int ids[100]; // Buffer for results
    int count = rtree_search(x, y, x, y, ids, 100);

    for (int i = 0; i < count; i++) {
        int id = ids[i];
        if (id >= 0 && id < census->numPolygons) {
            SHPObject* polygon = census->polygons[id];
            if (polygon && pointInPolygon(polygon, x, y)) {
                return id;
            }
        }
    }

    return -1;
}

double getPopulationDensityRTree(struct Census *census, double lat, double lon, double max_density)
{
    int id = findPolygon(census, lon, lat);
    if (id != -1) {
        // Read population and area attributes
        double population = DBFReadDoubleAttribute(census->dbf, id, census->pop_field_index);
        double area = DBFReadDoubleAttribute(census->dbf, id, census->area_field_index);

        // Calculate population density
        double density = (area > 0) ? (population / area) : 0.0;
        density = (density < 0 || isnan(density)) ? 0.0 : fmin(density, max_density);
        return density;
    } else {
        return 0; 
    }
}

/*
    Returns census data bounding box coordinates
*/
double *censusBoundingBox(struct Census *census)
{
    SHPObject *shape = NULL;
    int num_shapes, shape_type;
    double min_bound[4], max_bound[4];

    // Get shapefile info
    SHPGetInfo(census->shapefile, &num_shapes, &shape_type, min_bound, max_bound);

    double *bounds = (double *) malloc(4*sizeof(double));
    bounds[0] = min_bound[1];   // North
    bounds[1] = max_bound[1];   // South
    bounds[2] = min_bound[0];   // West
    bounds[3] = max_bound[0];   // East

    return bounds;
}