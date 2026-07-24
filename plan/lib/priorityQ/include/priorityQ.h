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

#ifndef PRIORITY_QUEUE_H
#define PRIORITY_QUEUE_H

#include "../../search/include/search.h"
#include "../../node/include/node.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct Node Node;

typedef struct PriorityQueueNode {
    Node *node;              // The item stored in the priority queue
    double priority;         // Priority of the item
} PriorityQueueNode;

typedef struct PriorityQueue {
    PriorityQueueNode *heap; // Array to store the heap
    size_t size;             // Current size of the heap
    size_t capacity;         // Capacity of the heap
    int isMinHeap;           // 1 for min-heap, 0 for max-heap
} PriorityQueue;

// Function prototypes
PriorityQueue *createPriorityQueue(size_t capacity, int isMinHeap);
void destroyPriorityQueue(PriorityQueue *pq);
void pqInsert(Node *node, PriorityQueue *pq);
Node *pqPop(PriorityQueue *pq);
int pqIsEmpty(PriorityQueue *pq);
int pqContains(PriorityQueue *pq, void *key, int (*cmp)(void *, void *));
void pqRemove(PriorityQueue *pq, void *key, int (*cmp)(void *, void *));
void pqUpdatePriority(PriorityQueue *pq, void *key, double newPriority, int (*cmp)(void *, void *));

#endif // PRIORITY_QUEUE_H
