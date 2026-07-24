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
#    R*-Tree C++ Wrapper Functions                                                  %
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

#include "rtree_wrapper.h"
#include <spatialindex/SpatialIndex.h>
#include <vector>
#include <iostream>

class RTreeIndex {
public:
    SpatialIndex::IStorageManager *storage;
    SpatialIndex::ISpatialIndex *rtree;

    RTreeIndex() {
        storage = SpatialIndex::StorageManager::createNewMemoryStorageManager();
        SpatialIndex::id_type indexIdentifier = 0;
        rtree = SpatialIndex::RTree::createNewRTree(*storage, 0.7, 50, 100, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
        // rtree = SpatialIndex::RTree::createNewRTree(*storage, 0.8, 20, 200, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
    }

    ~RTreeIndex() {
        delete rtree;
        delete storage;
    }

    void insert(double minX, double minY, double maxX, double maxY, int id) {
        double low[2] = { minX, minY };
        double high[2] = { maxX, maxY };
        rtree->insertData(0, nullptr, SpatialIndex::Region(low, high, 2), id);

    }

    void search(double minX, double minY, double maxX, double maxY, std::vector<int> &results) {
        MyVisitor visitor(results);
        double low[2] = { minX, minY };
        double high[2] = { maxX, maxY };
        rtree->intersectsWithQuery(SpatialIndex::Region(low, high, 2), visitor);
    }

private:
    class MyVisitor : public SpatialIndex::IVisitor {
    public:
        std::vector<int> &ids;
        MyVisitor(std::vector<int> &r) : ids(r) {}

        void visitNode(const SpatialIndex::INode &/*node*/) {}
        void visitData(const SpatialIndex::IData &data) {
            ids.push_back(data.getIdentifier());
        }
        void visitData(std::vector<const SpatialIndex::IData *> &/*v*/) {}
    };
};

// Global R-tree instance
static RTreeIndex *g_rtree = nullptr;

extern "C" {

void rtree_init() {
    if (!g_rtree) g_rtree = new RTreeIndex();
}

void rtree_insert(double minX, double minY, double maxX, double maxY, int id) {
    if (g_rtree) g_rtree->insert(minX, minY, maxX, maxY, id);
}

int rtree_search(double minX, double minY, double maxX, double maxY, int *results, int max_results) {
    if (!g_rtree) return 0;
    std::vector<int> found_ids;
    g_rtree->search(minX, minY, maxX, maxY, found_ids);
    int count = std::min((int)found_ids.size(), max_results);
    for (int i = 0; i < count; i++) {
        results[i] = found_ids[i];
    }
    return count;
}

void rtree_destroy() {
    if (g_rtree) {
        delete g_rtree;
        g_rtree = nullptr;
    }
}

}