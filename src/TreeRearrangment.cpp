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
 * @file TreeRearrangment.cpp
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
 */
#include "TreeRearrangment.hpp"


//===================================================================================================================
//double recompute_lk(PhyTree *tree,double k){
//    return k;
//}
//===================================================================================================================
void
nodes_within_radius(PhyTree *start_node, PhyTree *node, int radius, bool save, std::vector<move_info> &list_nodes) {

    if (!save) {
        save = true;
    } else {
        move_info m;
        m.node1 = start_node;
        m.node2 = node;
        list_nodes.push_back(m);
    }

    if (radius <= 0) {
        return;
    }

    if (!node->isLeaf()) {
        radius--;
        nodes_within_radius(start_node, node->get_left_child(), radius, save, list_nodes);
        nodes_within_radius(start_node, node->get_right_child(), radius, save, list_nodes);
    }

}

//===================================================================================================================
void nodes_within_radius_up(PhyTree *start_node, PhyTree *node, int radius, int direction,
                            std::vector<move_info> &list_nodes) {
    unsigned int idx;

    //TODO: check binary tree condition!

    move_info m;
    m.node1 = start_node;
    m.node2 = node;
    list_nodes.push_back(m);

    if (radius <= 0) {
        return;
    }

    radius--;
    if (direction == 0) {
        if (node->getParent() != NULL) {
            idx = node->indexOf();
            nodes_within_radius_up(start_node, node->getParent(), radius, idx, list_nodes);
        }
        nodes_within_radius(start_node, node->get_right_child(), radius, true, list_nodes);
    } else if (direction == 1) {
        if (node->getParent() != NULL) {
            idx = node->indexOf();
            nodes_within_radius_up(start_node, node->getParent(), radius, idx, list_nodes);
        }
        nodes_within_radius(start_node, node->get_left_child(), radius, true, list_nodes);
    }

}

//===================================================================================================================
void get_list_nodes_within_radius(PhyTree *node, int radius, std::vector<move_info> &list_nodes) {
    bool save;

    save = false;

    nodes_within_radius(node, node, radius, save, list_nodes);

    if (node->getParent() != NULL) {
        nodes_within_radius_up(node, node->getParent(), radius, node->indexOf(), list_nodes);
    }

}

//===================================================================================================================
std::vector<PhyTree *> fill_with_nodes(PhyTree *n) {
    std::vector<PhyTree *> list_nodes_n;

    list_nodes_n.push_back(n);
    while (n->getParent() != NULL) {
        n = n->getParent();
        list_nodes_n.push_back(n);
    }

    return list_nodes_n;
}

//===================================================================================================================
std::vector<PhyTree *> get_unique(std::vector<PhyTree *> &list_nodes_n1, std::vector<PhyTree *> &list_nodes_n2) {
    std::vector<PhyTree *> list_nodes;
    PhyTree *n1;
    PhyTree *n2;

    while (list_nodes_n1.size() > 0 && list_nodes_n2.size() > 0) {
        n1 = list_nodes_n1.at(list_nodes_n1.size() - 1);
        n2 = list_nodes_n2.at(list_nodes_n2.size() - 1);
        if (n1 == n2) {
            list_nodes.push_back(n1);
            list_nodes_n1.pop_back();
            list_nodes_n2.pop_back();
        } else {
            break;
        }
    }

    while (list_nodes_n1.size() > 0) {
        n2 = list_nodes_n1.at(list_nodes_n1.size() - 1);
        list_nodes.push_back(n1);
        list_nodes_n1.pop_back();
    }

    while (list_nodes_n2.size() > 0) {
        n2 = list_nodes_n2.at(list_nodes_n2.size() - 1);
        list_nodes.push_back(n2);
        list_nodes_n2.pop_back();
    }

    std::reverse(list_nodes.begin(), list_nodes.end());

    return list_nodes;
}


//===================================================================================================================
std::vector<PhyTree *> get_path_from_nodes(PhyTree *n1, PhyTree *n2) {
    std::vector<PhyTree *> list_nodes_n0;
    std::vector<PhyTree *> list_nodes_n1;
    std::vector<PhyTree *> list_nodes_n2;

    // add nodes from n1 to root
    list_nodes_n1 = fill_with_nodes(n1);

    // add nodes from n2 to root
    list_nodes_n2 = fill_with_nodes(n2);

    list_nodes_n0 = get_unique(list_nodes_n1, list_nodes_n2);

    return list_nodes_n0;
}

