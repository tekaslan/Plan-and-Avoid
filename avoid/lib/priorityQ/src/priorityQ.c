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
#    Priority Queue Handling Functions                                              %
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
#include "priorityQ.h"
#include "node.h"
#include <stdlib.h>
#include <stdio.h>

// Helper function to swap two nodes in the heap
static void swap(PriorityQueueNode *a, PriorityQueueNode *b) {
    PriorityQueueNode temp = *a;
    *a = *b;
    *b = temp;
}

// Helper function to compare priorities
static int compare(double a, double b, int isMinHeap) {
    return isMinHeap ? (a < b) : (a > b);
}

// Heapify up, it is a binary tree
static void heapifyUp(PriorityQueue *pq, size_t index) {

    // Do this untill root
    while (index > 0) {

        // Get the parent index
        size_t parent = (index - 1) / 2;

        // Only compare it with its parent, if lower swap parent and child, and compare it again with new parent
        if (compare(pq->heap[index].priority, pq->heap[parent].priority, pq->isMinHeap)) {
            swap(&pq->heap[index], &pq->heap[parent]);
            index = parent;
        } else {
            break;
        }
    }
}

// Heapify down
static void heapifyDown(PriorityQueue *pq, size_t index) {
    size_t left, right, smallestOrLargest, size = pq->size;

    while (index < size) {
        left = 2 * index + 1;
        right = 2 * index + 2;
        smallestOrLargest = index;

        if (left < size && compare(pq->heap[left].priority, pq->heap[smallestOrLargest].priority, pq->isMinHeap)) {
            smallestOrLargest = left;
        }

        if (right < size && compare(pq->heap[right].priority, pq->heap[smallestOrLargest].priority, pq->isMinHeap)) {
            smallestOrLargest = right;
        }

        if (smallestOrLargest == index) {
            break;
        }

        swap(&pq->heap[index], &pq->heap[smallestOrLargest]);
        index = smallestOrLargest;
    }
}

// Create a priority queue
PriorityQueue *createPriorityQueue(size_t capacity, int isMinHeap) {
    PriorityQueue *pq = (PriorityQueue *) malloc(sizeof(PriorityQueue));
    pq->heap = (PriorityQueueNode *) malloc(sizeof(PriorityQueueNode) * capacity);
    pq->size = 0;
    pq->capacity = capacity;
    pq->isMinHeap = isMinHeap;
    return pq;
}

// Destroy a priority queue
void destroyPriorityQueue(PriorityQueue *pq) {
    if (!pq) return;

    if (pq->heap) {
        for (size_t i = 0; i < pq->size; ++i) {
            if (pq->heap[i].node) free(pq->heap[i].node);
        }
        free(pq->heap);
        pq->heap = NULL;
    }
    free(pq);
}


// Insert an item into the priority queue
void pqInsert(Node *node, PriorityQueue *pq)
{
    if (pq->size == pq->capacity) {
        pq->capacity *= 2;
        pq->heap = (PriorityQueueNode *)realloc(pq->heap, sizeof(PriorityQueueNode) * pq->capacity);
    }

    pq->heap[pq->size].node = node;
    pq->heap[pq->size].priority = pq->isMinHeap ? node->f : -node->f;
    heapifyUp(pq, pq->size);
    pq->size++;
}

// Pop the highest-priority item from the queue
Node *pqPop(PriorityQueue *pq)
{
    if (pq->size == 0) {
        return NULL;
    }

    // This is the node to be popped
    Node *node = pq->heap[0].node;

    // 0 index is empty now, we cannot fill it with an intermediate node as it violates binary tree structure
    // Therefore, replace the last node into the 0-index.
    pq->heap[0] = pq->heap[--pq->size];

    // Now, make sure priority order is maintained by heapifying down
    heapifyDown(pq, 0);

    return node;
}

// Check if the priority queue is empty
int pqIsEmpty(PriorityQueue *pq) {
    return pq->size == 0;
}

// Check if the queue contains a specific key
int pqContains(PriorityQueue *pq, void *key, int (*cmp)(void *, void *)) {
    for (size_t i = 0; i < pq->size; i++) {
        if (cmp(pq->heap[i].node, key) == 0) {
            return 1;
        }
    }
    return 0;
}

// Remove a specific key from the queue
void pqRemove(PriorityQueue *pq, void *key, int (*cmp)(void *, void *)) {
    for (size_t i = 0; i < pq->size; i++) {
        if (cmp(pq->heap[i].node, key) == 0) {
            pq->heap[i] = pq->heap[--pq->size];
            heapifyDown(pq, i);
            return;
        }
    }
}

// Update the priority of a specific key
void pqUpdatePriority(PriorityQueue *pq, void *key, double newPriority, int (*cmp)(void *, void *)) {
    for (size_t i = 0; i < pq->size; i++) {
        if (cmp(pq->heap[i].node, key) == 0) {
            pq->heap[i].priority = pq->isMinHeap ? newPriority : -newPriority;
            heapifyUp(pq, i);
            heapifyDown(pq, i);
            return;
        }
    }
}
