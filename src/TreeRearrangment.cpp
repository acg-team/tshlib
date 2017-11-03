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
#include <numeric>
#include "TreeRearrangment.hpp"


TreeRearrangment::~TreeRearrangment() = default;


TreeRearrangment::TreeRearrangment(VirtualNode *node_source, int radius, bool preserve_blengths) {

    this->mset_sourcenode = node_source;
    this->mset_id = node_source->vnode_name + ":" + std::to_string(radius);
    this->mset_min_radius = radius;
    this->mset_max_radius = radius;
    this->mset_preserve_blenghts = preserve_blengths;
    this->mset_strategy = "undefined";

}

TreeRearrangment::TreeRearrangment(VirtualNode *node_source, int min_radius, int max_radius, bool preserve_blengths) {

    this->mset_sourcenode = node_source;
    this->mset_id = node_source->vnode_name + ":" + std::to_string(min_radius) + "-" + std::to_string(max_radius);
    this->mset_min_radius = min_radius;
    this->mset_max_radius = max_radius;
    this->mset_preserve_blenghts = preserve_blengths;
    this->mset_strategy = "mixed (NNI+SPR+TBR)";

}


void TreeRearrangment::getNodesInRadius(VirtualNode *node_source, int radius, bool includeSelf) {

    VirtualNode *node = node_source;

    auto moveInstance = new Move;

    if (!includeSelf) {
        includeSelf = true;

    } else {
        if (radius == 0) {
            moveInstance->setTargetNode(node);
            moveInstance->setMoveClass(this->mset_cur_radius);
            this->addMove(moveInstance);
        }
    }

    if (radius < 0) {
        return;
    }

    if (!node->isTerminalNode()) {
        radius--;
        this->getNodesInRadius(node->getNodeLeft(), radius, includeSelf);
        this->getNodesInRadius(node->getNodeRight(), radius, includeSelf);
    }


}


void TreeRearrangment::fillListMoves(bool includeSelf) {

    // Flag the nodes ac

    // Generate an array of integers representing the whole range of radius to use
    std::vector<int> x((unsigned long) (this->mset_max_radius - this->mset_min_radius) + 1);
    std::iota(std::begin(x), std::end(x), this->mset_min_radius);

    // For each radius in the array perform a search in the node space
    // (WARN: this loop is recursive = first left, second right, then up )
    for (auto &radius : x) {
        this->mset_cur_radius = radius;
        this->getNodesInRadius(this->mset_sourcenode, radius, includeSelf);

        if (nullptr != this->mset_sourcenode->getNodeUp()) {
            this->getNodesInRadiusUp(this->mset_sourcenode->getNodeUp(),
                                     radius - 1,
                                     this->mset_sourcenode->indexOf());
        }
    }
}


void TreeRearrangment::addMove(Move *move) {
    this->mset_moves.emplace_back(move);

}


void TreeRearrangment::getNodesInRadiusUp(VirtualNode *node_source, int radius, int direction) {

    auto vnode = new VirtualNode;
    auto moveInstance = new Move;
    unsigned int idx;

    vnode = node_source;

    //TODO: check binary tree condition!
    if (radius == 0) {
        moveInstance->setTargetNode(vnode);
        moveInstance->setMoveClass(this->mset_cur_radius);
        this->addMove(moveInstance);
    }

    if (radius > 0) {
        if (direction == 0) {
            radius--;
            if (vnode->getNodeUp() != nullptr) {

                idx = vnode->indexOf();
                this->getNodesInRadiusUp(vnode->getNodeUp(), radius, idx);
            }

            this->getNodesInRadius(vnode->getNodeRight(), radius, true);

        } else if (direction == 1) {
            radius--;
            if (vnode->getNodeUp() != nullptr) {

                idx = vnode->indexOf();
                this->getNodesInRadiusUp(vnode->getNodeUp(), radius, idx);
            }

            this->getNodesInRadius(vnode->getNodeLeft(), radius, true);

        } else if (direction == 2) {
            // If the direction is 2, we are moving across the pseudoroot, therefore no need to decrease the radious
            if (vnode->getNodeUp() != nullptr) {

                // Get the nodes in the radius from the node we reached after crossing the pseudoroot
                this->getNodesInRadius(vnode, radius, false);

            }

        }
    }
}


void TreeRearrangment::printListMoves() {

    std::cout << "[set " << this->mset_id << "] " << this->mset_strategy << " strategy" << std::endl;
    std::cout << "[class]\t(P\t; Q)" << std::endl;
    for (auto &nmove: this->mset_moves) {
        std::cout << "[" << nmove->move_class << "]\t(" << this->mset_sourcenode->vnode_name << "\t; " << nmove->getTargetNode()->vnode_name << ")" << std::endl;
    }

}


unsigned long TreeRearrangment::getNumberOfMoves() {

    return this->mset_moves.size();
}


bool TreeRearrangment::applyMove(unsigned long moveID) {
    bool status = false;
    try {


        status = this->mset_sourcenode->swapNode(this->mset_moves.at(moveID)->getTargetNode());

    } catch (const std::exception &e) {

        std::cout << "a standard exception was caught, with message '"
                  << e.what() << "'" << std::endl;

    }

    return status;

}


bool TreeRearrangment::revertMove(unsigned long moveID) {
    bool status = false;
    try {

        status = this->mset_moves.at(moveID)->getTargetNode()->swapNode(this->mset_sourcenode);


    } catch (const std::exception &e) {

        std::cout << "a standard exception was caught, with message '"
                  << e.what() << "'" << std::endl;

    }

    return status;
}


Move *TreeRearrangment::getMove(unsigned long moveID) {

    return this->mset_moves.at(moveID) ?: nullptr;
}


VirtualNode *TreeRearrangment::getSourceNode() {

    return this->mset_sourcenode ?: nullptr;
}


Move::~Move() = default;


Move::Move() {

    this->move_classes[0] = "Undef";
    this->move_classes[3] = "NNI";
    this->move_classes[4] = "SPR";
    this->move_classes[99] = "TBR";


    this->move_id = 0;
    this->move_name = "undefined";
    this->move_desc = "undefined";
    this->move_class = "undefined";

    this->move_lk = -INFINITY;

    this->move_applied = false;

}


void Move::setTargetNode(VirtualNode *target_node) {

    this->move_targetnode = target_node;

}


void Move::deleteTargetNode() {

    this->move_targetnode = nullptr;

}


VirtualNode *Move::getTargetNode() {

    return this->move_targetnode;

}


void Move::setMoveClass(int Value) {
    std::string class_name;
    if (Value >= 4) {
        class_name = this->move_classes[4];
    } else {
        class_name = this->move_classes[Value];
    }

    this->move_class = class_name;
}





void nodes_within_radius(PhyTree *start_node, PhyTree *node, int radius, std::vector<move_info> &list_nodes) {

//    if (!save) {
//          save = true;
//    } else {
    move_info m;
    m.node1 = start_node;
    m.node2 = node;
    list_nodes.push_back(m);
//    }

    if (radius <= 0) {
        return;
    }

    if (!node->isLeaf()) {
        radius--;
        nodes_within_radius(start_node, node->get_left_child(), radius, list_nodes);
        nodes_within_radius(start_node, node->get_right_child(), radius, list_nodes);
    }

}

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
        nodes_within_radius(start_node, node->get_right_child(), radius, list_nodes);
    } else if (direction == 1) {
        if (node->getParent() != NULL) {
            idx = node->indexOf();
            nodes_within_radius_up(start_node, node->getParent(), radius, idx, list_nodes);
        }
        nodes_within_radius(start_node, node->get_left_child(), radius, list_nodes);
    }

}

void
get_list_nodes_within_radius(PhyTree *node, int radius, std::vector<move_info> &list_nodes_left, std::vector<move_info> &list_nodes_right, std::vector<move_info> &list_nodes_up) {
    //bool save;

    //save = false;

    //nodes_within_radius(node, node, radius, save, list_nodes_down);
    //radius--;
    nodes_within_radius(node, node->get_left_child(), radius, list_nodes_left);
    nodes_within_radius(node, node->get_right_child(), radius, list_nodes_right);

    if (node->getParent() != NULL) {
        nodes_within_radius_up(node, node->getParent(), radius, node->indexOf(), list_nodes_up);
    }

}

std::vector<PhyTree *> fill_with_nodes(PhyTree *n) {
    std::vector<PhyTree *> list_nodes_n;

    list_nodes_n.push_back(n);
    while (n->getParent() != NULL) {
        n = n->getParent();
        list_nodes_n.push_back(n);
    }

    return list_nodes_n;
}

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

