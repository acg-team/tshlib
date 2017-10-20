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


class Move {

private:

public:

    Move();

    ~Move();

    void setTargetNode(PhyTree *target_node);

    void deleteTargetNode();

protected:

    int move_id;
    std::string move_name;
    std::string move_desc;
    PhyTree *move_targetnode;
    double move_lk;
    bool move_applied;
    std::string move_class;

};

class TreeRearrangment {

private:
    PhyTree *mset_sourcenode;
    std::string mset_id;
    int mset_radius;
    bool mset_preserve_blenghts;
    std::vector<Move *> mset_moves;
    std::string mset_strategy;

    void getNodesInRadius(PhyTree *node_source, int radius, bool save);

    void getNodesInRadiusUp(PhyTree *node_source, int radius, int direction);

    void addMove(Move *move);


public:
    TreeRearrangment(PhyTree *node_source, int radius, bool preserve_blengths);

    ~TreeRearrangment();

    void fillListMoves(bool saveMove);

protected:

};


void
nodes_within_radius(PhyTree *start_node, PhyTree *node, int radius, bool save, std::vector<move_info> &list_nodes);

void nodes_within_radius_up(PhyTree *start_node, PhyTree *node, int radius, int direction,
                            std::vector<move_info> &list_nodes);

void get_list_nodes_within_radius(PhyTree *node, int radius, std::vector<move_info> &list_nodes);

std::vector<PhyTree *> fill_with_nodes(PhyTree *n);

std::vector<PhyTree *> get_unique(std::vector<PhyTree *> &list_nodes_n1, std::vector<PhyTree *> &list_nodes_n2);

std::vector<PhyTree *> get_path_from_nodes(PhyTree *n1, PhyTree *n2);

#endif //TSHLIB_TREEREARRANGEMENT_HPP
