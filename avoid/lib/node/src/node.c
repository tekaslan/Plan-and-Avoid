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
#    Search Node Handling Functions                                                 %
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

#include "node.h"
#include <stdlib.h>

// Create a new action
Action *createAction(double deltaCourse, double length, double gamma)
{
    Action *action = (Action *) malloc(sizeof(Action));
    action->deltaCourse = deltaCourse;
    action->length = length;
    action->gamma = gamma;

    return action;
};

// Create a new node
Node *createNode(struct Pos *state, Node *parent, struct Action *action)
{
    Node *newNode = (Node *)malloc(sizeof(Node));

    if (newNode == NULL) {
        perror("Memory allocation failed");
        exit(1);
    }

    newNode->state = *state;

    if (parent){
        newNode->parent = parent;
    } else {
        newNode->parent = NULL;
    }
    
    if (parent) {
        newNode->l = parent->l + action->length;
    } else {
        newNode->l = action->length;
    }
    newNode->depth = parent ? parent->depth + 1 : 0; // Correct depth calculation

    if (action == NULL) {
        newNode->action = NULL; 
    } else {
        newNode->action = (struct Action *)malloc(sizeof(struct Action));
        if (newNode->action == NULL) {
            perror("Memory allocation failed for action in node creation.");
            free(newNode); newNode = NULL;
            exit(1);
        }

        // Copy the action data (important!):
        memcpy(newNode->action, action, sizeof(struct Action));
    }

    return newNode;
}


/*
    Explored list
*/
closedList *createClosedList(size_t capacity)
{   
    closedList * explored = (closedList *) malloc(sizeof(closedList));
    explored->nodes = (Node *) malloc(capacity * sizeof(Node));
    explored->capacity = capacity;
    explored->size = 0;

    return explored;
}

// Insert an item into the closed list
void clInsert(Node *node, SearchProblem *problem)
{
    if (problem->explored->size == problem->explored->capacity) {
        problem->explored->capacity *= 2;
        problem->explored->nodes = (Node *) realloc(problem->explored->nodes, sizeof(Node) * problem->explored->capacity);
    }

    problem->explored->nodes[problem->explored->size] = *node;
    problem->explored->size++;
}

// Destroy a node
void destroyNode(Node *node) {
    if (node == NULL) {
        return;
    }

    if (node->action != NULL) {
        free(node->action);
        node->action = NULL;
    }

    free(node);
    node = NULL;
}

// Destroy a tree
void destroyTree(Node *node) {
    if (node == NULL) return;

    Node *current = node;
    while (current != NULL) {
        Node *next = current->parent;
        destroyNode(current);
        current = next;
    }
}

// Expand a node to generate child nodes
Node **expandNode(Node *node, SearchProblem *problem)
{

    // Get the array of actions from the problem's actions function
    Action *actions = getActions(node, problem);

    int numChildren = 5;

    // Allocate memory for child nodes
    Node **children = malloc(numChildren * sizeof(*children));
    if (!children) {
        perror("Failed to allocate memory for child nodes");
        exit(EXIT_FAILURE);
    }

    // Create child nodes for each action
    for (int i = 0; i < numChildren; i++) {

        Action *actionCopy = (Action *) malloc(sizeof(Action));
        if (!actionCopy) {
            perror("Failed to allocate memory for actionCopy");

            // Free previously allocated children and actions
            for (int j = 0; j < i; j++) {
                destroyNode(children[j]);
            }
            free(children); children = NULL;
            free(actions); actions = NULL;
            exit(EXIT_FAILURE);
        }
        *actionCopy = actions[i];

        children[i] = childNode(node, problem, actionCopy);
        free(actionCopy); actionCopy = NULL;

        // Check if childNode returned NULL (allocation failure):
        if (children[i] == NULL) {
            perror("Failed to create child node");
            // Free previously allocated children and actions
            for (int j = 0; j <= i; j++) { // Free up to and including the current child
                if (children[j] != NULL) { // Check for NULL before freeing
                    destroyNode(children[j]); 
                }
            }
            free(children); children = NULL;
            free(actions); actions = NULL;
            exit(EXIT_FAILURE);
        }
    }
    free(actions); actions = NULL;

    return children;
}

// Generate a child node from an action
Node *childNode(Node *node, SearchProblem *problem, Action *action)
{
    // Expand parent state to a child state
    struct Pos *nextState;
    nextState = result(node, action, problem);

    // Compute the cost of the child state
    Node *newNode = createNode(nextState, node, action);

    // Timestamp (ignores groundspeed transients)
    double gs = groundspeed(nextState, problem);
    newNode->t = node->t + newNode->action->length/(gs * KT_2_MS * M_2_FT);

    // Cost
    if (nextState->alt > 0){
        newNode->f = stateCost(newNode, problem);
    } else {
        newNode->f = INFINITY;
    }

    free(nextState); nextState = NULL;

    return newNode;
}

// Extract the solution path
Path *solution(Node *node, int *solutionLength, SearchProblem *problem)
{
    int pathLength = 0;
    Node **pathNodes = path(node, &pathLength);

    Path *waypoints = (Path *) malloc(sizeof(Path));
    if (!waypoints) {
        perror("Failed to allocate memory for Path");
        exit(EXIT_FAILURE);
    }

    waypoints->nodes = (Node *) malloc(sizeof(Node) * pathLength);
    if (!waypoints->nodes) {
        perror("Failed to allocate memory for Path nodes");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < pathLength; i++) {
        // Copy state
        waypoints->nodes[i].state.lat = pathNodes[i]->state.lat;
        waypoints->nodes[i].state.lon = pathNodes[i]->state.lon;
        waypoints->nodes[i].state.hdg = pathNodes[i]->state.hdg;
        
        if (problem->altitudeModification && i < pathLength - 1 ) {
            waypoints->nodes[i].state.alt = pathNodes[i]->state.alt + problem->altitudeModification[i];
        } else {
            waypoints->nodes[i].state.alt = pathNodes[i]->state.alt;
        }

        // Copy other scalar members
        waypoints->nodes[i].l = pathNodes[i]->l;
        waypoints->nodes[i].depth = pathNodes[i]->depth;
        waypoints->nodes[i].gd1 = pathNodes[i]->gd1;
        waypoints->nodes[i].gd2 = pathNodes[i]->gd2;
        waypoints->nodes[i].gchi = pathNodes[i]->gchi;
        waypoints->nodes[i].gp = pathNodes[i]->gp;
        waypoints->nodes[i].ga = pathNodes[i]->ga;
        waypoints->nodes[i].cumul_ga = pathNodes[i]->cumul_ga;
        waypoints->nodes[i].cumul_gp = pathNodes[i]->cumul_gp;
        waypoints->nodes[i].g = pathNodes[i]->g;
        waypoints->nodes[i].t = pathNodes[i]->t;

        // Set parent to NULL to avoid accidental double free
        waypoints->nodes[i].parent = NULL;

        // Deep copy action
        if (pathNodes[i]->action) {
            waypoints->nodes[i].action = malloc(sizeof(Action));
            if (!waypoints->nodes[i].action) {
                perror("Failed to allocate memory for Action");
                exit(EXIT_FAILURE);
            }
            memcpy(waypoints->nodes[i].action, pathNodes[i]->action, sizeof(Action));
        } else {
            waypoints->nodes[i].action = NULL;
        }
    }

    waypoints->depth = pathLength;
    *solutionLength = pathLength;
    free(pathNodes); pathNodes = NULL;

    return waypoints;
}

// Extract the solution path
Action *solutionActions(Node *node, int *solutionLength, SearchProblem *problem){
    int pathLength = 0;
    Node **pathNodes = path(node, &pathLength);

    Action *actions = (Action *) malloc((pathLength-1)*sizeof(Action));

    for (int i = 1; i < pathLength; i++){
        actions[i-1].deltaCourse  = pathNodes[i-1]->action->deltaCourse;
        actions[i-1].length       = pathNodes[i-1]->action->length;
        actions[i-1].gamma        = pathNodes[i-1]->action->gamma + problem->dGamma;
        actions[i-1].gamma_turn   = pathNodes[i-1]->action->gamma_turn;
    }
    *solutionLength = pathLength - 1;
    free(pathNodes); pathNodes = NULL;

    return actions;
}

// Get the full path from root to the node
Node **path(Node *node, int *pathLength) {
    if (node == NULL) {
        *pathLength = 0;
        return NULL;
    }

    Node *temp = node;

    // Count path length
    while (temp != NULL) {
        (*pathLength)++;
        temp = temp->parent;
    }

    Node **pathNodes = (Node **)malloc(sizeof(Node *) * (*pathLength));
    if (pathNodes == NULL) {
        perror("Memory allocation failed"); // Use perror for better error messages
        exit(1);
    }

    // Fill path in correct order
    temp = node;
    for (int i = (*pathLength) - 1 ; i >= 0; i--) {
        pathNodes[i] = temp;
        if (temp->parent != NULL && (uintptr_t)temp->parent < 0x1000) { // Invalid pointer check
            fprintf(stderr, "Invalid parent pointer detected\n");
            free(pathNodes); pathNodes = NULL;
            exit(1);
        }
        temp = temp->parent;
    }
    return pathNodes;
}

void destroyClosedList(closedList *clist) {
    if (!clist) return;
    for (size_t i = 0; i < clist->size; i++) {
        if (clist->nodes[i].action) free(clist->nodes[i].action);
        // No need to free parent; it's shared and freed elsewhere
    }
    free(clist->nodes);
    free(clist);
}