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

#ifndef RTREE_WRAPPER_H
#define RTREE_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void *RTreeHandle;

// Initializes the R-tree
void rtree_init();

// Inserts a bounding box with an associated ID into the R-tree
void rtree_insert(double minX, double minY, double maxX, double maxY, int id);

// Searches the R-tree for bounding boxes that intersect the given region
// Returns the number of found elements and stores the results in the provided array
int rtree_search(double minX, double minY, double maxX, double maxY, int *results, int max_results);

// Destroys the R-tree and frees memory
void rtree_destroy();

#ifdef __cplusplus
}
#endif

#endif // RTREE_WRAPPER_H
