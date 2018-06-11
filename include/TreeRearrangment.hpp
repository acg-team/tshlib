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

namespace tshlib {

    enum class TreeSearchStopCondition {
        iterations, convergence
    };

    enum class TreeSearchHeuristics {
        classic_NNI, classic_SPR, classic_Mixed, particle_swarm, hillclimbing, greedy, swap, phyml, mixed, nosearch
    };

    enum class TreeRearrangmentOperations {
        classic_NNI, classic_SPR, classic_TBR, classic_Mixed
    };

    class Move {

    private:

    protected:
        VirtualNode *move_targetnode;   /* Pointer to the target node found during the node search */
        VirtualNode *move_sourcenode;   /* Pointer to the source node  */

    public:
        int move_id;                            /* Move ID - Useful in case of parallel independent executions*/
        std::string move_name;                  /* Move Name - Unused */
        int move_radius;                        /* Move Radius */
        MoveDirections move_direction;          /* Move Direction for applying a rotation to the VirtualNode pointers */
        double move_lk;                         /* Likelihood of the move if applied */
        bool move_applied;                      /* Indicator is set to true if the move is applied to the tree */
        std::string move_class;                 /* String indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */
        MoveType move_type;                     /* Integer indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */

        /*!
         * @brief Standard constructor
         */
        Move();

        /*!
        * @brief Standard deconstructor
        */
        ~Move();

        Move(const Move &inMove) {

            this->move_id = inMove.move_id;
            this->move_name = "copy_" + inMove.move_name;
            this->move_radius = inMove.move_radius;
            this->move_lk = inMove.move_lk;
            this->move_applied = inMove.move_applied;
            this->move_class = inMove.move_class;
            this->move_type = inMove.move_type;
            this->move_direction = inMove.move_direction;
            this->move_targetnode = inMove.move_targetnode;
            this->move_sourcenode = inMove.move_sourcenode;

        }

        Move &operator=(const Move &inMove) {
            this->move_id = inMove.move_id;
            this->move_name = "copy_" + inMove.move_name;
            this->move_radius = inMove.move_radius;
            this->move_lk = inMove.move_lk;
            this->move_applied = inMove.move_applied;
            this->move_class = inMove.move_class;
            this->move_type = inMove.move_type;
            this->move_direction = inMove.move_direction;
            this->move_targetnode = inMove.move_targetnode;
            this->move_sourcenode = inMove.move_sourcenode;
        };


        std::string getMoveDirection() const {
            std::string returnString;
            switch (move_direction) {
                case MoveDirections::left:
                    returnString = "left";
                    break;
                case MoveDirections::up:
                    returnString = "up";
                    break;
                case MoveDirections::up_left:
                    returnString = "up-left";
                    break;
                case MoveDirections::up_right:
                    returnString = "up-right";
                    break;
                case MoveDirections::right :
                    returnString = "right";
                    break;
                case MoveDirections::undef :
                    returnString = "undef";
                    break;

            }

            return returnString;
        }


        /*!
         * @brief Reset the protected move_targetnode field
         */
        void deleteTargetNode();

        /*!
         * @brief Returns the target node pointer
         * @return VirtualNode pointer of the target node
         */
        VirtualNode *getTargetNode();

/*!
        * @brief Returns the source node pointer
        * @return VirtualNode pointer of the source node
        */
        VirtualNode *getSourceNode();

        /*!
         * @brief Set the protected move_targetnode field
         * @param target_node PhyTree Pointer to the target node
         */
        void setTargetNode(VirtualNode *target_node);

        void setSourceNode(VirtualNode *source_node);

        void setMoveClass(int Value);

        MoveType getMoveType() const;

        void setRadius(int radius);

        void setDirection(MoveDirections direction);

        double getLikelihood() { return move_lk; };

        void recomputeLikelihood();

        void initMove();
    };


    class TreeRearrangment {
    private:

        Utree *tree;
        std::string mset_id;                /* Tree-rearrangment ID. Useful in case of parallel independent executions */
        int mset_min_radius;                /* Radius of the node search (for NNI must set it to 3) */
        int mset_max_radius;                /* Radius of the node search (for NNI must set it to 3) */
        bool mset_preserve_blenghts;        /* Switch to preserve branch lentghs in case the move is applied (i.e NNI vs SPR) */
        std::vector<Move *> mset_moves;     /* Vector containing the pre-computed moves */
        int mset_foundmoves;                /* Number of moves found during the defineMoves call */

        VirtualNode *mset_sourcenode;       /* Starting node from which starting the tree exploration */

    public:
        std::string mset_strategy;          /* Description of the node search strategy */

        TreeRearrangment();

        void initTreeRearrangment(VirtualNode *node_source, int radius, bool preserve_blengths);

        void initTreeRearrangment(Utree *ref_tree, int min_radius, int max_radius, bool preserve_blengths, VirtualNode *node_source);

        ~TreeRearrangment();

        VirtualNode *getSourceNode() const { return mset_sourcenode ?: nullptr; }

        void setSourceNode(VirtualNode *mset_sourcenode) { TreeRearrangment::mset_sourcenode = mset_sourcenode; }

        int getMinRadius() const { return mset_min_radius; }

        void setMinRadius(int mset_min_radius) { TreeRearrangment::mset_min_radius = mset_min_radius; }

        int getMaxRadius() const { return mset_max_radius; }

        void setMaxRadius(int mset_max_radius) { TreeRearrangment::mset_max_radius = mset_max_radius; }

        Utree *getTree() const { return tree; }

        void setTree(Utree *tree) { TreeRearrangment::tree = tree; }

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
