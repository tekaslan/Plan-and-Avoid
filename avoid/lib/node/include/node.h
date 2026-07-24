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

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "search.h"
#include "geo.h"
#include <stddef.h>

typedef struct SearchProblem SearchProblem;

// Node structure
typedef struct Node 
{
    struct Pos state;
    struct Node *parent;
    struct Action *action;
    double l;       // Path length so far
    double gd1;     // Optimal gliding cost
    double gd2;     // Min. remaining traversal cost
    double gchi;    // Course angle cost
    double gp;      // Local overflown population cost (ground risk)
    double ga;      // Local airspace cost
    double cumul_ga; // Cumulative airspace cost
    double cumul_gp; // Local overflown population cost
    double g;       // Cumulative cost
    double h;       // Heuristic cost
    double f;       // Total cost
    double t;       // Timestamp [s]
    int depth;      // Node depth
} Node;

// Closed list
typedef struct closedList
{
    struct Node * nodes;
    size_t size;
    size_t capacity;
} closedList;


// Path
typedef struct Path
{
    struct Node * nodes;
    size_t depth;
} Path; 

// Function prototypes
struct Action *createAction(double deltaCourse, double length, double gamma);
Node *createNode(struct Pos *state, Node *parent, struct Action *action);
void destroyNode(Node *node);
Node **expandNode(Node *node, SearchProblem *problem);
Node *childNode(Node *node, SearchProblem *problem, struct Action *action);
Path *solution(Node *node, int *solutionLength, SearchProblem *problem);
struct Action *solutionActions(Node *node, int *solutionLength, SearchProblem *problem);
Node **path(Node *node, int *pathLength);
closedList *createClosedList(size_t capacity);
void clInsert(Node *node, SearchProblem *problem);
void destroyClosedList(closedList *clist);

#ifdef __cplusplus
}
#endif