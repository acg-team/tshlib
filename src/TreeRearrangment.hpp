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

struct move_info {
    PhyTree *node1;
    PhyTree *node2;
    int ID;
    double lk;
};

struct utree_move_info {
    node *node1;
    node *node2;
    int ID;
    double lk;
};

class Move {

private:

public:
    /*!
     * @brief Standard constructor
     */
    Move();

    /*!
    * @brief Standard deconstructor
    */
    ~Move();

    /*!
     * @brief Set the protected move_targetnode field
     * @param target_node PhyTree Pointer to the target node
     */
    void setTargetNode(PhyTree *target_node);

    /*!
     * @brief Reset the protected move_targetnode field
     */
    void deleteTargetNode();

    /*!
     * @brief Returns the target node pointer
     * @return PhyTree pointer of the node
     */
    PhyTree *getTargetNode();

protected:

    int move_id;                    /* Move ID - Useful in case of parallel independent executions*/
    std::string move_name;          /* Move Name - Unused */
    std::string move_desc;          /* Move Desc - Unused */
    PhyTree *move_targetnode;       /* Pointer to the target node found during the node search */
    double move_lk;                 /* Likelihood of the move if applied */
    bool move_applied;              /* Indicator is set to true if the move is applied to the tree */
    std::string move_class;         /* Move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */

};

class TreeRearrangment {

private:

    /*!
     * @brief Recursive function to retrieve all the nodes within a fixed radius from a starting node
     * @param node_source PhyTree Pointer to the starting node
     * @param radius int Radius of the search (NNI = 1, SPR > 1)
     * @param save bool ?
     */
    void getNodesInRadius(PhyTree *node_source, int radius, bool save);

    /*!
     * @brief Recursive function to retrieve all the nodes within a fixed radius from a starting node
     * @param node_source   PhyTree Pointer to the starting node
     * @param radius        int Radius of the search (NNI = 1, SPR > 1)
     * @param direction     int ?
     */
    void getNodesInRadiusUp(PhyTree *node_source, int radius, int direction);

    /*!
     * @brief Append candidate move to the mset_moves vector
     * @param move Move Pointer to the candidate move object
     */
    void addMove(Move *move);


public:
    PhyTree *mset_sourcenode;       /* Starting node from which starting the tree exploration */
    std::string mset_id;            /* Tree-rearrangment ID. Useful in case of parallel independent executions */
    int mset_radius;                /* Radius of the node search (for NNI must set it to 1) */
    bool mset_preserve_blenghts;    /* Switch to preserve branch lentghs in case the move is applied (i.e NNI vs SPR) */
    std::vector<Move *> mset_moves; /* Vector containing the pre-computed moves */
    std::string mset_strategy;      /* Description of the node search strategy */


    TreeRearrangment(PhyTree *node_source, int radius, bool preserve_blengths);

    ~TreeRearrangment();

    /*!
     * @brief Perform a complete node search and fill the vector containing the candidate moves.
     * @param saveMove bool ?
     */
    void fillListMoves(bool saveMove);

protected:


};

void nodes_within_radius(PhyTree *start_node, PhyTree *node, int radius, std::vector<move_info> &list_nodes);

void nodes_within_radius_up(PhyTree *start_node, PhyTree *node, int radius, int direction,std::vector<move_info> &list_nodes);

void get_list_nodes_within_radius(PhyTree *node, int radius, std::vector<move_info> &list_nodes_left, std::vector<move_info> &list_nodes_right,std::vector<move_info> &list_nodes_up);

std::vector<PhyTree *> fill_with_nodes(PhyTree *n);

std::vector<PhyTree *> get_unique(std::vector<PhyTree *> &list_nodes_n1, std::vector<PhyTree *> &list_nodes_n2);

std::vector<PhyTree *> get_path_from_nodes(PhyTree *n1, PhyTree *n2);

#endif //TSHLIB_TREEREARRANGEMENT_HPP
