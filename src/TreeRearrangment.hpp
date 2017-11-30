/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * tshexe is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tshexe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with likpip. If not, see <http://www.gnu.org/licenses/>.
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

#include "PhyTree.hpp"
#include "Utree.hpp"
#include "Utilities.hpp"

enum class TreeSearchHeuristics{classic_NNI, classic_SPR, classic_Mixed, particle_swarm};

class Move {

private:

protected:
    VirtualNode *move_targetnode;   /* Pointer to the target node found during the node search */

public:
    int move_id;                    /* Move ID - Useful in case of parallel independent executions*/
    std::string move_name;          /* Move Name - Unused */
    int move_radius;                /* Move Radius */
    MoveDirections move_direction;  /* Move Direction for applying a rotation to the VirtualNode pointers */
    double move_lk;                 /* Likelihood of the move if applied */
    bool move_applied;              /* Indicator is set to true if the move is applied to the tree */
    std::string move_class;         /* String indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */
    MoveType move_type;             /* Integer indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */


    /*!
     * @brief Standard constructor
     */
    Move();

    /*!
    * @brief Standard deconstructor
    */
    ~Move();

    /*!
     * @brief Reset the protected move_targetnode field
     */
    void deleteTargetNode();

    /*!
     * @brief Returns the target node pointer
     * @return PhyTree pointer of the node
     */
    VirtualNode *getTargetNode();

    /*!
     * @brief Set the protected move_targetnode field
     * @param target_node PhyTree Pointer to the target node
     */
    void setTargetNode(VirtualNode *target_node);

    void setMoveClass(int Value);

    void setRadius(int radius);

    void setDirection(MoveDirections direction);

    void recomputeLikelihood();

    void initMove();
};


class TreeRearrangment {
private:

    std::string mset_id;                /* Tree-rearrangment ID. Useful in case of parallel independent executions */
    int mset_min_radius;                /* Radius of the node search (for NNI must set it to 3) */
    int mset_max_radius;                /* Radius of the node search (for NNI must set it to 3) */
    bool mset_preserve_blenghts;        /* Switch to preserve branch lentghs in case the move is applied (i.e NNI vs SPR) */
    std::vector<Move *> mset_moves;     /* Vector containing the pre-computed moves */

    VirtualNode *mset_sourcenode;       /* Starting node from which starting the tree exploration */

public:
    std::string mset_strategy;          /* Description of the node search strategy */

    TreeRearrangment();

    virtual void initTreeRearrangment(VirtualNode *node_source, int radius, bool preserve_blengths);

    virtual void initTreeRearrangment(VirtualNode *node_source, int min_radius, int max_radius, bool preserve_blengths);

    virtual ~TreeRearrangment();

    VirtualNode *getSourceNode();

    /*!
     * @brief Perform a complete node search and fill the vector containing the candidate moves.
     * @param includeSelf bool ?
     */
    virtual void defineMoves(bool includeSelf);

    virtual bool applyMove(unsigned long moveID);

    Move *getMove(unsigned long moveID);

    virtual bool revertMove(unsigned long moveID);

    virtual void printMoves();

    unsigned long getNumberOfMoves();

    virtual void selectBestMove(unsigned long moveID);

protected:

    /*!
     * @brief Recursive function to retrieve all the nodes within a fixed radius from a starting node
     * @param node_source   VirtualNode pointer to the starting node
     * @param radius_min    int Radius of the search (NNI = 1, SPR > 1)
     * @param radius_max    int Radius of the search (NNI = 1, SPR > 1)
     * @param includeSelf bool ?
     */
    virtual void getNodesInRadiusDown(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, bool includeSelf, MoveDirections direction);


    /*!
     * @brief Recursive function to retrieve all the nodes within a fixed radius from a starting node
     * @param node_source   PhyTree Pointer to the starting node
     * @param radius_min    int Radius of the search (NNI = 3, SPR > 3)
     * @param radius_max    int Radius of the search (NNI = 1, SPR > 1)
     * @param traverse_direction     NodePosition it indicates the direction of the node wrt parent node
     */
    void getNodesInRadiusUp(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, NodePosition traverse_direction);

    /*!
     * @brief Append candidate move to the mset_moves vector
     * @param move Move Pointer to the candidate move object
     */
    void addMove(Move *move);

};

namespace treesearchheuristics{

    void testTSH(Utree *input_tree, TreeSearchHeuristics tsh_strategy);

}





#endif //TSHLIB_TREEREARRANGEMENT_HPP
