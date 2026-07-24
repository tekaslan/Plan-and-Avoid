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
#    Air Traffic Handling Functions                                                 %
#    Airspace, Air Traffic, and Ground-risk Aware                                   %
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

#include "airtraffic.h"
#define _USE_MATH_DEFINES


/*
    Reads dynamic air traffic data
    Returns DynamicTraffic structure
*/
void loadDynamicTraffic(SearchProblem * problem){

    // Open file
    char dir[300];
    sprintf(dir, "lib/airtraffic/data/%s", problem->dynTraffic->datafile);
    FILE * file = fopen(dir, "rb");
    if (!file){
        perror("Cannot open dynamic traffic data file.\n");
        exit(EXIT_FAILURE);
    }

    // Check the file has the correct formatting
    char magic[8];
    if (fread(magic, 1, 8, file) != 8) {fprintf(stderr,"Bad data format.\n"); exit(EXIT_FAILURE); }
    if (memcmp(magic, "ATD1BIN", 7) != 0) {  // compare first 7 chars
        fprintf(stderr, "Not a DynamicTraffic binary. Make sure the data formatting is accurate.\n");
        exit(EXIT_FAILURE);
    }

    // Read the traffic data label (e.g., IMC_LOW, VFR_MED, etc.)
    if (fread(problem->dynTraffic->label, 1, 10, file) != 10) {fprintf(stderr,"Bad file.\n"); exit(EXIT_FAILURE); }
    problem->dynTraffic->label[9] = '\0'; // Force termination for printing

    // Read number of flights and number of time steps within the data file
    if (fread(&problem->dynTraffic->nflights, sizeof(uint32_t), 1, file) != 1 ||
        fread(&problem->dynTraffic->nsteps,   sizeof(uint32_t), 1, file) != 1) {
        fprintf(stderr, "Failed to read dimensions.\n");
        exit(EXIT_FAILURE);
    }

    // Read data statistics
    if (fread(&problem->dynTraffic->median, sizeof(uint32_t), 1, file) != 1 ||
        fread(&problem->dynTraffic->p90,    sizeof(uint32_t), 1, file) != 1 ||
        fread(&problem->dynTraffic->peak,   sizeof(uint32_t), 1, file) != 1) {
        fprintf(stderr, "Failed to read data statistics.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for flights
    problem->dynTraffic->flights = calloc(problem->dynTraffic->nflights, sizeof(FlightSeries));
    if (!problem->dynTraffic->flights) {fprintf(stderr, "FlightSeries memory cannot be allocated.\n"); exit(EXIT_FAILURE);}

    // Allocate memory for flight data
    for (int i = 0; i < problem->dynTraffic->nflights; i++) {
        problem->dynTraffic->flights[i].sample = calloc(problem->dynTraffic->nsteps, sizeof(FSSample));
    }

    // Read aircraft type
    uint32_t meta_bytes = 0;
    if (fread(&meta_bytes, sizeof(uint32_t), 1, file) != 1) {
        perror("meta_bytes read failed.\n");
        exit(EXIT_FAILURE);
    }

    uint8_t *type = calloc(problem->dynTraffic->nflights, sizeof(uint8_t));
    if (!type) {
        perror("Type allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    if (meta_bytes > 0) {
        if (meta_bytes < problem->dynTraffic->nflights) {
            fprintf(stderr, "meta_bytes < nflights (corrupt file).\n");
            exit(EXIT_FAILURE);
        }

        if (fread(type, sizeof(uint8_t), problem->dynTraffic->nflights, file)
            != problem->dynTraffic->nflights) {
            perror("Aircraft type[] read failed");
            exit(EXIT_FAILURE);
        }
    }

    // Now read the flight data
    for (uint32_t i = 0; i < problem->dynTraffic->nflights; i++) {

        float *lat   = calloc(problem->dynTraffic->nsteps, sizeof(float));
        float *lon   = calloc(problem->dynTraffic->nsteps, sizeof(float));
        float *alt   = calloc(problem->dynTraffic->nsteps, sizeof(float));
        float *track = calloc(problem->dynTraffic->nsteps, sizeof(float));
        float *gs   = calloc(problem->dynTraffic->nsteps, sizeof(float));
        int64_t *ts = calloc(problem->dynTraffic->nsteps, sizeof(int64_t));

        fread(lat,   sizeof(float), problem->dynTraffic->nsteps, file);
        fread(lon,   sizeof(float), problem->dynTraffic->nsteps, file);
        fread(alt,   sizeof(float), problem->dynTraffic->nsteps, file);
        fread(track, sizeof(float), problem->dynTraffic->nsteps, file);
        fread(gs, sizeof(float), problem->dynTraffic->nsteps, file);
        fread(ts, sizeof(int64_t), problem->dynTraffic->nsteps, file);

        // Assign type
        problem->dynTraffic->flights[i].type = type[i];

        for (uint32_t k = 0; k < problem->dynTraffic->nsteps; k++) {
            problem->dynTraffic->flights[i].sample[k].pos.lat   = lat[k];
            problem->dynTraffic->flights[i].sample[k].pos.lon   = lon[k];
            problem->dynTraffic->flights[i].sample[k].pos.alt   = alt[k];
            problem->dynTraffic->flights[i].sample[k].pos.hdg   = track[k];
            problem->dynTraffic->flights[i].sample[k].gs = gs[k];
            problem->dynTraffic->flights[i].sample[k].ts = ts[k];
            if (isnan(problem->dynTraffic->flights[i].sample[k].pos.lat)){
                problem->dynTraffic->flights[i].sample[k].valid = false;
            } else {
                problem->dynTraffic->flights[i].sample[k].valid = true;
            }
        }
        free(lat); free(lon); free(alt); free(track); free(gs); free(ts);
    }
    
    if (problem->dynTraffic == NULL) {
        fprintf(stderr, "Dynamic traffic data are unavailable.\n");
        return;
    }

    printf(
        "Imported dynamic traffic summary:\n"
        "  Label:          %s\n"
        "  Flights:        %d\n"
        "  Time steps:     %d\n"
        "  Median traffic: %d\n"
        "  90th percentile: %d\n"
        "  Peak traffic:   %d\n",
        problem->dynTraffic->label,
        problem->dynTraffic->nflights,
        problem->dynTraffic->nsteps,
        problem->dynTraffic->median,
        problem->dynTraffic->p90,
        problem->dynTraffic->peak);

    // Close the file
    fclose(file);
}

/*
    Reads air traffic heatmap dataset and stores it
    into trafficHeatmap structure
*/
void readHeatMap(airspaceHeatmap *trafficHeatmap) {

    // Allocate memory and read grid coordinates
    trafficHeatmap->lat_grid = calloc((trafficHeatmap->LAT_SIZE) , sizeof(float));
    trafficHeatmap->lon_grid = calloc((trafficHeatmap->LON_SIZE) , sizeof(float));
    trafficHeatmap->alt_grid = calloc((trafficHeatmap->ALT_SIZE) , sizeof(float));

    if (!trafficHeatmap->lat_grid || !trafficHeatmap->lon_grid || !trafficHeatmap->alt_grid) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = fopen(trafficHeatmap->gridfile, "rb");
    if (!fp) {
        printf("%s\n", trafficHeatmap->gridfile);
        perror("Failed to open grid file\n");
        exit(EXIT_FAILURE);
    }

    size_t n_read;
    n_read = fread(trafficHeatmap->lat_grid, sizeof(float), trafficHeatmap->LAT_SIZE, fp);
    if (n_read != trafficHeatmap->LAT_SIZE) {
        fprintf(stderr, "fread lat_grid failed: expected %d, got %zu\n", trafficHeatmap->LAT_SIZE, n_read);
        exit(EXIT_FAILURE);
    }

    n_read = fread(trafficHeatmap->lon_grid, sizeof(float), trafficHeatmap->LON_SIZE, fp);
    if (n_read != trafficHeatmap->LON_SIZE) {
        fprintf(stderr, "fread lon_grid failed: expected %d, got %zu\n", trafficHeatmap->LON_SIZE, n_read);
        exit(EXIT_FAILURE);
    }

    n_read = fread(trafficHeatmap->alt_grid, sizeof(float), trafficHeatmap->ALT_SIZE, fp);
    if (n_read != trafficHeatmap->ALT_SIZE) {
        fprintf(stderr, "fread alt_grid failed: expected %d, got %zu\n", trafficHeatmap->ALT_SIZE, n_read);
        exit(EXIT_FAILURE);
    }

    fclose(fp);

    // Open air traffic data file
    fp = fopen(trafficHeatmap->datafile, "rb");
    if (!fp) {
        printf("%s\n", trafficHeatmap->datafile);
        perror("Failed to open heatmap data file.\n");
        exit(EXIT_FAILURE);
    }

    // Heatmap data size
    if (trafficHeatmap->grid_on_cells) trafficHeatmap->HM_SIZE = (trafficHeatmap->LAT_SIZE-1) * (trafficHeatmap->LON_SIZE-1) * (trafficHeatmap->ALT_SIZE-1);
    else trafficHeatmap->HM_SIZE = trafficHeatmap->LAT_SIZE * trafficHeatmap->LON_SIZE * trafficHeatmap->ALT_SIZE;

    // Allocate memory
    trafficHeatmap->heatmap = (float *) calloc(trafficHeatmap->HM_SIZE, sizeof(float));
    
    // Read air traffic heatmap data
    n_read = fread(trafficHeatmap->heatmap, sizeof(float), trafficHeatmap->HM_SIZE, fp);
    if (n_read != trafficHeatmap->HM_SIZE) {
        fprintf(stderr, "fread heatmap failed: expected %zu, got %zu\n", trafficHeatmap->HM_SIZE, n_read);
        exit(EXIT_FAILURE);
    }

    fclose(fp);

}

/*
    Returns air traffic density of a given position
    using trafficHeatmap data structure
*/
float getAirTrafficDensity(struct Pos *pos, airspaceHeatmap *trafficHeatmap)
{
    if (pos->lat < trafficHeatmap->lat_grid[0] || pos->lat > trafficHeatmap->lat_grid[trafficHeatmap->LAT_SIZE - 1]) {
        return 0;
    }
    if (pos->lon < trafficHeatmap->lon_grid[0] || pos->lon > trafficHeatmap->lon_grid[trafficHeatmap->LON_SIZE - 1]) {
        return 0;
    }
    if (pos->alt < trafficHeatmap->alt_grid[0] || pos->alt > trafficHeatmap->alt_grid[trafficHeatmap->ALT_SIZE - 1]) {
        return 0;
    }

    int i = 0, j = 0, k = 0;

    for (i = 0; i < trafficHeatmap->LAT_SIZE - 1; i++) {
        if (trafficHeatmap->lat_grid[i] <= pos->lat && pos->lat <= trafficHeatmap->lat_grid[i + 1])
            break;
    }

    for (j = 0; j < trafficHeatmap->LON_SIZE - 1; j++) {
        if (trafficHeatmap->lon_grid[j] <= pos->lon && pos->lon <= trafficHeatmap->lon_grid[j + 1])
            break;
    }

    for (k = 0; k < trafficHeatmap->ALT_SIZE - 1; k++) {
        if (trafficHeatmap->alt_grid[k] <= pos->alt && pos->alt <= trafficHeatmap->alt_grid[k + 1])
            break;
    }

    int max_i = (trafficHeatmap->grid_on_cells) ? trafficHeatmap->LAT_SIZE - 2 : trafficHeatmap->LAT_SIZE - 1;
    int max_j = (trafficHeatmap->grid_on_cells) ? trafficHeatmap->LON_SIZE - 2 : trafficHeatmap->LON_SIZE - 1;
    int max_k = (trafficHeatmap->grid_on_cells) ? trafficHeatmap->ALT_SIZE - 2 : trafficHeatmap->ALT_SIZE - 1;

    if (i > max_i) i = max_i;
    if (j > max_j) j = max_j;
    if (k > max_k) k = max_k;

    if (i < 0 || j < 0 || k < 0){
        printf("Invalid heatmap indices... %d %d %d\n", i,j,k);
    }

    int stride_lat = (trafficHeatmap->grid_on_cells) ? (trafficHeatmap->LAT_SIZE - 1) : trafficHeatmap->LAT_SIZE;
    int stride_lon = (trafficHeatmap->grid_on_cells) ? (trafficHeatmap->LON_SIZE - 1) : trafficHeatmap->LON_SIZE;


    int idx = i + stride_lat * j + stride_lat * stride_lon * k;

    return trafficHeatmap->heatmap[idx];
}


/*
    Reads air traffic corridor dataset and stores it
    into airCorridor structure
*/
int readAirCorridors(struct airCorridor *corridors, char * filedir)
{   
    printf("* * * Air Traffic Corridors * * *\n");

    // Open air traffic corridors data file
    FILE* datafile;
    char directory[200];

    // Set the directory
    strcpy(directory, filedir);

    // Open the data file
    datafile = fopen(directory,"r");
    if (!datafile){
		perror("Error");
		fprintf(stderr,"Failure opening file %s\n",directory);
		return EXIT_FAILURE;
	}

    // Read all lines
    char line[500];
    char *aux;
    char *token;
    struct airCorridor *tmp_corridor = corridors;
    tmp_corridor->segment = (struct segment *) malloc(sizeof(struct segment));
    struct segment *tmp_segment = (struct segment *) malloc(sizeof(struct segment));

    // Read a line
    fgets(line,500,datafile);
    int corridorCount = 0;
    int segmentCount  = 0;
    while (!feof(datafile))
    {

        // Check if the line contains air corridor identifier
        aux = strstr(line,"Identifier");
        while (aux)
        {
            // Update corridor count
            corridorCount++;

            // Copy the identifier into the corridor structure
            token = &aux[strcspn(aux,"H")];
            strcpy(tmp_corridor->ident, token);

            // Skip a line
            fgets(line,500,datafile);

            // Define auxilary strings to identify segment lines
            aux = strstr(line,"Segment");

            tmp_segment = tmp_corridor->segment;
            while (aux)
            {
                
                // Update segment count
                segmentCount++;

                // Read and assign vertex coordinates to segment surfaces
                // Iterate over surfaces
                for (int i = 0; i < 6; i++)            
                {

                    // Get to coordinates' line
                    fgets(line,500,datafile); fgets(line,500,datafile);

                    // Iterate over vertices
                    for (int j = 0; j < 4; j++)         
                    {   

                        // Assign coordinates of i-th surface's j-th vertex
                        // token = strtok(aux," ");
                        token = strtok(line," ");
                        if (!token) {
                            fprintf(stderr, "Error parsing coordinate on line: %s\n", line);
                            exit(EXIT_FAILURE);
                        }
                        tmp_segment->surfaces[i].xyz[j].X = atof(token);

                        token = strtok(NULL," ");
                        if (!token) {
                            fprintf(stderr, "Error parsing coordinate on line: %s\n", line);
                            exit(EXIT_FAILURE);
                        }
                        tmp_segment->surfaces[i].xyz[j].Y = atof(token);

                        token = strtok(NULL," ");
                        if (!token) {
                            fprintf(stderr, "Error parsing coordinate on line: %s\n", line);
                            exit(EXIT_FAILURE);
                        }
                        tmp_segment->surfaces[i].xyz[j].Z = atof(token);
                        
                        // Skip to the next vertex
                        fgets(line,500,datafile);
                    }
                }

                // Skip a line
                fgets(line,500,datafile);

                // Check if there is more segment to handle. If so, allocate memory for the next one.
                aux = strstr(line,"Segment");
                if (aux)
                {
                    tmp_segment->next = (struct segment *) malloc(sizeof(struct segment));
                    tmp_segment = tmp_segment->next;
                }
            }

            // Check if the line contains air corridor identifier
            aux = strstr(line,"Identifier");
            if (aux)
            {
                tmp_corridor->next = (struct airCorridor *) malloc(sizeof(struct airCorridor));
                tmp_corridor = tmp_corridor->next;
                tmp_corridor->next = NULL;
                tmp_corridor->segment = (struct segment *) malloc(sizeof(struct segment));
                tmp_segment->next = NULL;
            }
        }
    };
    tmp_corridor->next = NULL;
    tmp_segment->next = NULL;
    
    printf("Database is stored.\n\n");

    // Close the data file
    fclose(datafile);

    return EXIT_SUCCESS;
}

/*
    Checks if a given ECEF coordinate is inside any of the air traffic corridors
*/
int isInsidePolyhedron(struct airCorridor *corridors, struct PosXYZ pos)
{   

    // Initialize flag indicating if the point is inside any of corridors (0: Inside, 1:Not inside)
    int flag = 1;

    // Define a temporary corridor and segment structure to iterate over them
    struct airCorridor *tmp_corridor;
    tmp_corridor = corridors;
    struct segment *tmp_segment;

    // Define surface iteration index
    int i;
    int segmentCount;

    // Pre-allocate dot product vector for each surface of a segment
    double dot_prod[6];

    // Pre-allocate point-to-centroid vector
    struct PosXYZ point2centroidVector;

    // Iterate over all corridors as long as flag is zero
    while (tmp_corridor && flag) {
        // Reset the segment pointer for the current corridor
        tmp_segment = tmp_corridor->segment;
        segmentCount = 0;

        // Iterate over all segments as long as flag is zero
        while (tmp_segment && flag) {

            segmentCount++;

            // Initialize iteration index
            i = 0;
            for (; i < 6; i++) {
                // Get surface centroid
                surfaceCentroid(&tmp_segment->surfaces[i]);

                // Get surface normal
                surfaceNormal(&tmp_segment->surfaces[i]);

                // The second surface is defined counterclockwise, so normal vector points toward outside of the segment
                // Multiply by -1, so that the normal vector points inward.
                if (i == 1){
                    tmp_segment->surfaces[i].normal.X = -tmp_segment->surfaces[i].normal.X;
                    tmp_segment->surfaces[i].normal.Y = -tmp_segment->surfaces[i].normal.Y;
                    tmp_segment->surfaces[i].normal.Z = -tmp_segment->surfaces[i].normal.Z;
                }

                // Point to centroid vector
                createVector(&pos, &tmp_segment->surfaces[i].centroid, &point2centroidVector);

                // Calculate the dot products of surface normal and point-to-centroid vectors
                dot_prod[i] = dot(&tmp_segment->surfaces[i].normal, &point2centroidVector);

                // Ideally, if one of the dot products is positive, break the loop as the point is outside.
                // Due to truncation and rounding errors, instead of zero, a small positive value (1e-3) is used.
                if (dot_prod[i] > 1e-6){
                    break;
                }
            }

            // If index i = 6, it means that the dot product loop is completed, and
            // all dot products are non-negative; therefore, the point is inside a polyhedron.
            if (i == 6) {
                flag = 0;
            }

            // Skip to the next segment
            tmp_segment = tmp_segment->next;
        }

        // Skip to the next corridor
        tmp_corridor = tmp_corridor->next;
    }

    return flag;
}

/*
    Checks if a point overlaps with an air corridor surface along the surface normal
*/
int doesOverlap(struct PosXYZ pos, struct surface *surface)
{
    // Initialize flag indicating if the point and surface overlap (0:Overlap, 1:Doesn't overlap)
    int flag = 1;

    // Define surface edge vectors and edge normals
    struct PosXYZ v[4];
    struct PosXYZ n[4];

    // Define point-to-vertex vectors
    struct PosXYZ p2v[2];

    // Define dot product vector
    double dot_prod[4];

    // Get surface normal
    surfaceNormal(surface);

    // Get surface edge vectors and edge normals
    for (int j = 0; j < 4; j++)
    {   
        // Edge vectors
        if (j == 0)
        {
            createVector(&surface->xyz[0], &surface->xyz[3], &v[0]);
        }
        else
        {
            createVector(&surface->xyz[j], &surface->xyz[j-1], &v[j]);
        }

        // Edge normals
        cross(&v[j], &surface->normal, &n[j]);
    }

        // Get point-to-vertex vectors - Only first and third vertices are enough as each share two different edges.
        createVector(&pos, &surface->xyz[0], &p2v[0]);
        createVector(&pos, &surface->xyz[2], &p2v[1]);

        // Calculate dot products of point-to-vertex vectors and edge normals
        dot_prod[0] = dot(&p2v[0], &n[0]);
        dot_prod[1] = dot(&p2v[0], &n[1]);
        dot_prod[2] = dot(&p2v[1], &n[2]);
        dot_prod[3] = dot(&p2v[1], &n[3]);

        if (dot_prod[0] < 1e-3 && dot_prod[1] < 1e-3 && dot_prod[2] < 1e-3 && dot_prod[3] < 1e-3)
        {
            flag = 0;
        }
        
    return flag;
}

/*
    Calculates angles between point-to-vertex vectors and edge vectors
*/
void point2vertexVec_to_edgeVector_angles(struct PosXYZ *pos, struct surface *surface, double angles[8])
{
    // Create point-to-vertex vectors
    struct PosXYZ a;
    struct PosXYZ b;
    struct PosXYZ c;
    struct PosXYZ d;
    createVector(pos, &surface->xyz[0], &a);
    createVector(pos, &surface->xyz[3], &b);
    createVector(pos, &surface->xyz[1], &c);
    createVector(pos, &surface->xyz[2], &d);

    // Create edge vectors
    struct PosXYZ v1;
    struct PosXYZ v2;
    struct PosXYZ v3;
    struct PosXYZ v4;
    createVector(&surface->xyz[0], &surface->xyz[3], &v1);
    createVector(&surface->xyz[1], &surface->xyz[0], &v2);
    createVector(&surface->xyz[2], &surface->xyz[1], &v3);
    createVector(&surface->xyz[3], &surface->xyz[2], &v4);

    // Calculate angles between point-to-vertex vectors and edge vectors
    angles[0] = angle_between_vectors(&a, &v1);
    angles[1] = angle_between_vectors(&b, &v1);
    angles[2] = angle_between_vectors(&a, &v2);
    angles[3] = angle_between_vectors(&c, &v2);
    angles[4] = angle_between_vectors(&c, &v3);
    angles[5] = angle_between_vectors(&d, &v3);
    angles[6] = angle_between_vectors(&d, &v4);
    angles[7] = angle_between_vectors(&b, &v4);
}

/*
    Performs exhaustive search to find the minimum distance
    from a given coordinate to air traffic corridors
*/
double minimum_distance(struct Pos *pos, struct airCorridor *corridors)
{

    // Transform LLA coordinates to the ECEF frame
    struct PosXYZ PosXYZ;

    // LLA to ECEF transformation for X and Y coordinates. Keep the altitude same.
    double N = geo_earths_curvature(pos);
    PosXYZ.X = (N + pos->alt*FT_2_M)*cos(pos->lat*DEG_2_RAD)*cos(pos->lon*DEG_2_RAD)*M_2_NM;
    PosXYZ.Y = (N + pos->alt*FT_2_M)*cos(pos->lat*DEG_2_RAD)*sin(pos->lon*DEG_2_RAD)*M_2_NM;
    PosXYZ.Z = pos->alt*FT_2_NM;

    // Define intermediate variables
    int flag_isinside = 1;              // Flag indicating if point is inside a corridor (0:Inside, 1:Outside)
    int flag_overlaps = 1;              // Flag indicating if point overlaps with a surface along its normal vector direction (0: Overlaps, 1: Doesn't overlap)
    double d_min_corridors = 999;       // Initial min. distance from point to air corridors. This distance is what we are aiming to find.
    double d_min_segment = 999;         // Initial min. distance from point to a corridor segment
    double d_min_surface = 999;         // Initial min. distance from point to a segment surface
    double d_min_surface_new = 999;     // Initial running min. distance from point to a segment surface
    double d_min_edges = 999;           // Initial min. distance from point to a surface edge
    double d_min_edges_new = 999;       // Initial running min. distance from point to a surface edge
    double angles[8];
    int corridor_count = 1;
    int segment_count;

    // Check if the point is inside a corridor
    flag_isinside = isInsidePolyhedron(corridors, PosXYZ);

    // If so, set the minimum distance to zero.
    if (!flag_isinside) {
        d_min_corridors = 0;
    }
    // If not, search for the minimum distance
    else {   
        // Assign temporary structures to iterate over corridors and segments
        struct airCorridor *tmp_corridor;
        tmp_corridor = corridors;
        struct segment *tmp_segment;
        
        // Iterate over all corridors
        // Iterate over all corridors
        while (tmp_corridor) {
            // Get the segment
            tmp_segment = tmp_corridor->segment;

            // Iterate over all segments
            segment_count = 1;
            while (tmp_segment) {

                // Iterate over all surfaces
                for (int i = 0; i < 6; i++) {
                    // Flip left-to-right the second surface to make the vertex orientation correct
                    if (i == 1) {
                        fliplr(&tmp_segment->surfaces[i]);
                    }

                    // Check if the point overlaps with the surface
                    flag_overlaps = doesOverlap(PosXYZ, &tmp_segment->surfaces[i]);

                    // If so, get the min. distance which is the distance from a point to surface
                    if (!flag_overlaps) {
                        d_min_surface_new = distance_point2surface(&PosXYZ, &tmp_segment->surfaces[i]);

                        // Search for the minimum distance to the surface
                        if (d_min_surface_new < d_min_surface) {
                            d_min_surface = d_min_surface_new;
                        }
                    }
                    // If not, search for the min. distance
                    else {
                        // Get angles between point-to-vertex vectors and edge vectors
                        point2vertexVec_to_edgeVector_angles(&PosXYZ, &tmp_segment->surfaces[i], angles);

                        // Initialize array to store the distances to edges
                        d_min_edges = 999;

                        // Check if the point overlaps with any of the edges of the surface
                        for (int j = 0; j < 4; j++) {
                            // The angle condition holds true if the point overlaps with j-th edge
                            if (angles[2*j] <= M_PI/2 && angles[2*j+1] >= M_PI/2){   
                                if (j == 0) {
                                    d_min_edges_new = distance_point2edge(&PosXYZ, &tmp_segment->surfaces[i].xyz[3], &tmp_segment->surfaces[i].xyz[0]);
                                } else {
                                    d_min_edges_new = distance_point2edge(&PosXYZ, &tmp_segment->surfaces[i].xyz[j-1], &tmp_segment->surfaces[i].xyz[j]);
                                }

                                // Search for the minimum distance to edges
                                if (d_min_edges_new < d_min_edges) {
                                    d_min_edges = d_min_edges_new;
                                }
                            }
                        }

                        // Assign the minimum distance to edges as the minimum distance to that particular surface
                        if (d_min_edges != 999) {
                            d_min_surface_new = d_min_edges;
                        } else {
                            d_min_surface_new = distance_point2vertex(&PosXYZ, &tmp_segment->surfaces[i]);
                        }
                    }
                    // Search for the minimum distance to the surface
                    if (d_min_surface_new < d_min_surface) {
                        d_min_surface = d_min_surface_new;
                    }
                }

                // Search for the minimum distance to the segment
                if (d_min_surface < d_min_segment) {
                    d_min_segment = d_min_surface;
                }

                // Update the minimum distance to surface
                d_min_surface = 999;

                // Move to the next segment if it exists
                tmp_segment = tmp_segment->next;
                segment_count++;
            }

            // Search for the minimum distance to the corridors
            if (d_min_segment < d_min_corridors) {
                d_min_corridors = d_min_segment;
            }

            // Update the minimum distance to segment
            d_min_segment = 999;

            // Move to the next corridor if it exists
            tmp_corridor = tmp_corridor->next;
            corridor_count++;
        }

    }

    // Return the minimum distance from the point to all corridors
    return d_min_corridors;
}

/*
    Air corridor occupation cost based on point-to-polyhedron distance
    and buffer distance
*/
double airCorridorOccupationCost(struct Pos *pos, struct airCorridor *corridors) {
    
    // Point's min. distance to air traffic corridors
    double d = minimum_distance(pos, corridors);

    // Buffer distance
    double max_d = 500*FT_2_NM;
    if (d > max_d) d = max_d;
    return 1.0 - d / max_d;
}

/*
    Computes ground speed along track given a state
*/
double groundspeed(struct Pos * state, SearchProblem * problem){

    // Relative wind direction
    double angDiff = wrapTo180(problem->windDirection - state->hdg); 
    
    return problem->ac->airspeed*FTS_2_KTS - problem->windSpeed*cos(angDiff*DEG_2_RAD);
}

/*
    Returns segment traversal time (int)
*/
int traversaldeltat(struct Node * node, SearchProblem * problem){

    // Ground speed
    double gs = groundspeed(&(node->state), problem);

    // Segment flight duration rounded to the ceiling integer
    double dt = node->action->length/(gs * KTS_2_FTS);
    
    return ceil(dt);
}


/*
    Computes indices and separation of active vertical points
*/
void VCPA(struct Pos * state, int aircraftID, int timeIDX, double * activeDeltav, int * activeCount, SearchProblem *problem){

    *activeCount = 0;

    // Search for the closest point along the intruder's path vertically
    int totalActive = 0;

    // Get a sample, if not valid skip
    FSSample s = problem->dynTraffic->flights[aircraftID].sample[timeIDX];
    if (!s.valid) return;

    // Get the altitude difference
    double dv = fabs(state->alt - s.pos.alt);

    // If the difference is less than a threshold, mark it active collision point
    const double alt_threshold = 5000.0;  // FAA Advisory Circular AC No: 20-165B, Airworthiness Approval of Automatic Dependent Surveillance - Broadcast OUT Systems
    double err_baro_own      = (state->alt  >= alt_threshold) ? problem->ERR_BARO_A5000 : problem->ERR_BARO_B5000;
    double err_baro_intruder = (s.pos.alt >= alt_threshold) ? problem->ERR_BARO_A5000 : problem->ERR_BARO_B5000;
    dv -= err_baro_own + err_baro_intruder;
    dv = fmax(0.0, dv);

    if (dv <= (problem->VERT_CLEARANCE)){

        // Collision-active intruder count
        totalActive++;

        // Write altitude separation of the active point
        *activeDeltav = dv;
    }

    // Overwrite total count
    *activeCount = totalActive;
}


/*
    Computes separation violation cost
    deltah: Geofence horizontal separation
    deltav: Geofence vertical separation
*/
double CPACost(double deltah, double deltav, SearchProblem * problem){

    if ((deltah >= problem->HORZ_CLEARANCE) || (deltav >= problem->VERT_CLEARANCE)) return 0;

    // Compute the horizontal margin
    double mh = fabs(fmin(0, deltah - problem->HORZ_CLEARANCE));

    // Compute the vertical margin
    double mv = fabs(fmin(0, deltav - problem->VERT_CLEARANCE));

    // Horizontal penalty
    double Jh;
    double kh = 1;
    Jh = 1 - 2 / (1 + exp(kh * mh));

    // Vertical penalty
    double Jv;
    double kv = 0.005;
    Jv = 1 - 2 / (1 + exp(kv * mv));

    return Jh*Jv;
}

/*
    Returns maximum separation penalty
    given ownship state, intruder ID, and local time index
*/
double JCPA(struct Pos * state, int aircraftID, int timeIDX, SearchProblem *problem){

    // Check aircraft type
    uint8_t type = problem->dynTraffic->flights[aircraftID].type;
    if (type == 0){
        printf("Type 0 aircraft detected.\n");
        return 0;
    }

    // Check vertical separation first
    int activeCount;  // Loss of well-clear count
    double deltav;    // Vertical separation
    VCPA(state, aircraftID, timeIDX, &deltav, &activeCount, problem);
    
    // If no vertically active point, all well-clear.
    if (activeCount == 0) return 0;

    // Compute CPA
    // Get horizontal separation—C.G. to C.G. distance
    double deltah, course;
    FSSample sample = problem->dynTraffic->flights[aircraftID].sample[timeIDX];
    geo_dist(state, &sample.pos, &deltah, &course, &problem->GeoOpt);

    // Find geofence to geofence distance
    if (sample.pos.alt > problem->IAF_Alt) {
        deltah = fmax(0.0, deltah - problem->RNP_TER - problem->ERR_ADSB_P);
    } else {
        deltah = fmax(0.0, deltah - problem->RNP_FIN - problem->ERR_ADSB_P);
    }

    // Compute separation cost
    double J = CPACost(deltah, deltav, problem);

    return J;
}