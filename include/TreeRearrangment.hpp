/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * Tree Search Heuristic Library (TshLib) is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tree Search Heuristic Library (TshLib) is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file TreeRearrangment.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 10 2017
 * @version 1.0
 * @maintainer Lorenzo Gatti
 * @maintainer Massimo Maiolo
 * @email lg@lorenzogatti.me
 * @email massimo.maiolo@zhaw.ch
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit: http://www.lorenzogatti.me
 */
#ifndef TSHLIB_TREEREARRANGEMENT_HPP
#define TSHLIB_TREEREARRANGEMENT_HPP

//#include "PhyTree.hpp"
#include "Utree.hpp"
#include "Utilities.hpp"
#include "Move.hpp"

namespace tshlib {

    class TreeRearrangment {
    protected:
        VirtualNode *trSourceNode_;                /* Starting node from which starting the tree exploration */

        int trSearchRadius_min;                    /* Radius of the node search (for NNI must set it to 3) */
        int trSearchRadius_max;                    /* Radius of the node search (for NNI must set it to 3) */
        std::vector<Move *>  trMoveSet;            /* Vector containing the pre-computed moves */
        TreeSearchHeuristics trStrategy;           /* Strategy used to test (apply/revert) the candidate moves  */
        TreeRearrangmentOperations trTreeCoverage;

    protected:
        /* Strategy used to define the candidate moves -- max coverage of the tree space */
        bool trPreserveBranchLenghts_;             /* Switch to preserve branch lentghs in case the move is applied (i.e NNI vs SPR) */

        Utree *UTree_;


    private:
        std::string trUID_;                        /* Tree-rearrangment ID. Useful in case of parallel independent executions */
        int trCandidateMovesFound_;                /* Number of moves found during the defineMoves call */
        bool trInitialized_;

    public:

        TreeRearrangment();

        void initialize();

        void initTreeRearrangment(VirtualNode *node_source, int radius, bool preserve_blengths);

        void initTreeRearrangment(Utree *ref_tree, int min_radius, int max_radius, bool preserve_blengths, VirtualNode *node_source);

        ~TreeRearrangment();

        VirtualNode *getSourceNode() const { return trSourceNode_ ?: nullptr; }

        void setSourceNode(VirtualNode *mset_sourcenode) { TreeRearrangment::trSourceNode_ = mset_sourcenode; }

        int getMinRadius() const { return trSearchRadius_min; }

        void setMinRadius(int mset_min_radius) { TreeRearrangment::trSearchRadius_min = mset_min_radius; }

        int getMaxRadius() const { return trSearchRadius_max; }

        void setMaxRadius(int mset_max_radius) { TreeRearrangment::trSearchRadius_max = mset_max_radius; }

        Utree *getTree() const { return UTree_; }

        void setTree(Utree *tree) { TreeRearrangment::UTree_ = tree; }

        /*!
         * @brief Perform a complete node search and fill the vector containing the candidate moves.
         * @param includeSelf bool ?
         */
        void defineMoves(bool includeSelf, bool allowDuplicatedMoves = true, TreeSearchHeuristics moveSchema = TreeSearchHeuristics::swap);

        const std::vector<VirtualNode *> updatePathBetweenNodes(unsigned long moveID, std::vector<VirtualNode *> inPath);

        bool applyMove(unsigned long moveID);

        void commitMove(int moveID);

        Move *getMove(unsigned long moveID);

        bool revertMove(unsigned long moveID);

        void printMoves();

        unsigned long getNumberOfMoves();

        Move *selectBestMove(double value);

        void storeMove(Move *inMove);

        void setTreeTopology(Utree *inTree);

        void displayRearrangmentStatus(int idMove, bool printTree);

        std::string getStrategy() const {
            std::string rtToken;
            switch (trStrategy) {
                case TreeSearchHeuristics::swap:
                    rtToken = "Swap";
                    break;
                case TreeSearchHeuristics::phyml:
                    rtToken = "phyML";
                    break;
                case TreeSearchHeuristics::mixed:
                    rtToken = "Mixed(Swap+phyML)";
                    break;
                case TreeSearchHeuristics::nosearch:
                    rtToken = "no-search";
                    break;
            }
            return rtToken;
        }

        TreeRearrangmentOperations getRearrangmentCoverage() const {
            return trTreeCoverage;
        }

        std::string getRearrangmentCoverageDescription() const {
            std::string rtToken;
            switch (trTreeCoverage) {
                case TreeRearrangmentOperations::classic_NNI:
                    rtToken = "NNI-like";
                    break;
                case TreeRearrangmentOperations::classic_SPR:
                    rtToken = "SPR-like";
                    break;
                case TreeRearrangmentOperations::classic_TBR:
                    rtToken = "TBR-like";
                    break;
                case TreeRearrangmentOperations::classic_Mixed:
                    rtToken = "Complete";
                    break;
            }
            return rtToken;
        }

        void setTreeCoverage(TreeRearrangmentOperations trTreeCoverage) {
            TreeRearrangment::trTreeCoverage = trTreeCoverage;
        }

        void setStrategy(TreeSearchHeuristics trStrategy) {
            TreeRearrangment::trStrategy = trStrategy;
        }

    protected:

        /*!
         * @brief Recursive function to retrieve all the nodes within a fixed radius from a starting node
         * @param node_source   VirtualNode pointer to the starting node
         * @param radius_min    int Radius of the search (NNI = 1, SPR > 1)
         * @param radius_max    int Radius of the search (NNI = 1, SPR > 1)
         * @param includeSelf bool ?
         */
        void getNodesInRadiusDown(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, bool includeSelf, MoveDirections direction, bool allowDuplicatedMoves,
                                  TreeSearchHeuristics moveSchema = TreeSearchHeuristics::swap);


        /*!
         * @brief Recursive function to retrieve all the nodes within a fixed radius from a starting node
         * @param node_source   PhyTree Pointer to the starting node
         * @param radius_min    int Radius of the search (NNI = 3, SPR > 3)
         * @param radius_max    int Radius of the search (NNI = 1, SPR > 1)
         * @param traverse_direction     NodePosition it indicates the direction of the node wrt parent node
         */
        void getNodesInRadiusUp(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, NodePosition traverse_direction, bool allowDuplicatedMoves,
                                TreeSearchHeuristics moveSchema = TreeSearchHeuristics::swap);

        /*!
         * @brief Append candidate move to the mset_moves vector
         * @param move Move Pointer to the candidate move object
         */
        void addMove(Move *move, bool allowDuplicatedMoves = true, TreeSearchHeuristics moveSchema = TreeSearchHeuristics::swap);
    };


}

namespace treesearchheuristics {
    using namespace tshlib;

    void performTestTreeSearch(Utree *input_tree, TreeSearchHeuristics tsh_strategy);

}


#endif //TSHLIB_TREEREARRANGEMENT_HPP
