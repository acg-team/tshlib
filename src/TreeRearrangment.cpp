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
#include <iomanip>
#include <chrono>
#include <glog/logging.h>
#include <gflags/gflags.h>
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


void TreeRearrangment::getNodesInRadiusDown(VirtualNode *node_source, int radius, MoveDirections direction, bool includeSelf) {

    VirtualNode *node = node_source;

    auto moveInstance = new Move;

    if (!includeSelf) {
        includeSelf = true;

    } else {
        if (radius == 0) {
            moveInstance->setDirection(direction);
            moveInstance->setRadius(this->mset_cur_radius);
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
        // Direction is required to know where to look at in the tree when the move is performed
        this->getNodesInRadiusDown(node->getNodeLeft(), radius, direction, includeSelf);
        this->getNodesInRadiusDown(node->getNodeRight(), radius, direction, includeSelf);
    }


}


void TreeRearrangment::defineMoves(bool includeSelf) {

    // Flag the nodes according to their position on the tree (left or right or above the source node -- p node).

    // Generate an array of integers representing the whole range of radius to use
    std::vector<int> x((unsigned long) (this->mset_max_radius - this->mset_min_radius) + 1);
    std::iota(std::begin(x), std::end(x), this->mset_min_radius);

    // For each radius in the array perform a search in the node space
    // (WARN: this loop is recursive = first left, second right, then up )
    for (auto &radius : x) {
        this->mset_cur_radius = radius;
        if (!this->mset_sourcenode->isTerminalNode()) {
            this->getNodesInRadiusDown(this->mset_sourcenode->getNodeLeft(), radius - 1, MoveDirections::left, includeSelf);
            this->getNodesInRadiusDown(this->mset_sourcenode->getNodeRight(), radius - 1, MoveDirections::right, includeSelf);
        }
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


void TreeRearrangment::getNodesInRadiusUp(VirtualNode *node_source, int radius, NodePosition traverse_direction) {

    auto vnode = new VirtualNode;
    auto moveInstance = new Move;
    NodePosition idx;
    MoveDirections moving_direction;

    vnode = node_source;

    //TODO: check binary tree condition!
    if (radius == 0) {

        switch (traverse_direction) {
            case NodePosition::left:
                moving_direction = MoveDirections::up_left;
                break;
            case NodePosition::right:
                moving_direction = MoveDirections::up_right;
                break;
            case NodePosition::up:
                moving_direction = MoveDirections::up;
                break;
            default:
                moving_direction = MoveDirections::undef;

        }
        moveInstance->setDirection(moving_direction);
        moveInstance->setRadius(this->mset_cur_radius);
        moveInstance->setTargetNode(vnode);
        moveInstance->setMoveClass(this->mset_cur_radius);
        this->addMove(moveInstance);
    }

    if (radius > 0) {
        if (traverse_direction == NodePosition::left) {
            radius--;
            if (vnode->getNodeUp() != nullptr) {

                idx = vnode->indexOf();
                this->getNodesInRadiusUp(vnode->getNodeUp(), radius, idx);
            }

            this->getNodesInRadiusDown(vnode->getNodeRight(), radius, MoveDirections::up, true);

        } else if (traverse_direction == NodePosition::right) {
            radius--;
            if (vnode->getNodeUp() != nullptr) {

                idx = vnode->indexOf();
                this->getNodesInRadiusUp(vnode->getNodeUp(), radius, idx);
            }

            this->getNodesInRadiusDown(vnode->getNodeLeft(), radius, MoveDirections::up, true);

        } else if (traverse_direction == NodePosition::up) {
            // If the traverse_direction is 2, we are moving across the pseudoroot, therefore no need to decrease the radious
            if (vnode->getNodeUp() != nullptr) {

                // Get the nodes in the radius from the node we reached after crossing the pseudoroot
                this->getNodesInRadiusDown(vnode, radius, MoveDirections::up, false);

            }

        }
    }
}


void TreeRearrangment::printMoves() {

    VLOG(2) << "[set " << this->mset_id << "] " << this->mset_strategy << " strategy" << std::endl;
    VLOG(2) << "[class]\t(P\t; Q)" << std::endl;
    for (auto &nmove: this->mset_moves) {
        VLOG(2) << "[" << nmove->move_class << "]\t" << nmove->move_radius
                  << "\t(" << this->mset_sourcenode->vnode_name << "\t; "
                  << nmove->getTargetNode()->vnode_name << ")\t"
                  << static_cast<int>(nmove->move_direction) << std::endl;
    }

}


unsigned long TreeRearrangment::getNumberOfMoves() {

    return this->mset_moves.size();
}


bool TreeRearrangment::applyMove(unsigned long moveID) {

    VirtualNode *pnode = this->mset_sourcenode;
    VirtualNode *qnode = this->mset_moves.at(moveID)->getTargetNode();
    bool revertRotations = false;
    // Swap pnode with qnode according to the direction found during the move configuration
    // If the swap is performed correctly then the function returns true otherwise false
    return pnode->swapNode(qnode, this->mset_moves.at(moveID)->move_direction, revertRotations);
}


bool TreeRearrangment::revertMove(unsigned long moveID) {

    VirtualNode *pnode = this->mset_moves.at(moveID)->getTargetNode();
    VirtualNode *qnode = this->mset_sourcenode;
    bool revertRotations = true;
    // Swap pnode with qnode according to the direction found during the move configuration
    // If the swap is performed correctly then the function returns true otherwise false
    return pnode->swapNode(qnode, MoveDirections::up, revertRotations);

}


Move *TreeRearrangment::getMove(unsigned long moveID) {

    return this->mset_moves.at(moveID) ?: nullptr;
}


VirtualNode *TreeRearrangment::getSourceNode() {

    return this->mset_sourcenode ?: nullptr;
}

void TreeRearrangment::selectBestMove(unsigned long moveID) {

}


Move::~Move() = default;


Move::Move() {

    this->move_id = 0;
    this->move_name = "undef";
    this->move_type = MoveType::undef;
    this->move_class = "undef";
    this->move_direction = MoveDirections::undef;

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

    if (Value >= 4 and Value < 10) {
        this->move_type = MoveType::SPR;
        this->move_class = "SPR";
    } else if (Value >= 10) {
        this->move_type = MoveType::TBR;
        this->move_class = "TBR";
    } else if (Value == 3) {
        this->move_type = MoveType::NNI;
        this->move_class = "NNI";
    } else {
        this->move_type = MoveType::undef;
        this->move_class = "undef";
    }
}


void Move::setRadius(int radius) {

    this->move_radius = radius;
}


void Move::setDirection(MoveDirections direction) {

    this->move_direction = direction;

}

void Move::recomputeLikelihood() {

}


void ::treesearchheuristics::testTSH(Utree *input_tree, TreeSearchHeuristics tsh_strategy) {

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    int min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
    int max_radius = input_tree->getMaxNodeDistance(); // Hard coded max value for a small tree (this ensures the complete q-node search)

    unsigned long total_exec_moves = 0;

    // Print node description with neighbors
    for (auto &vnode:input_tree->listVNodes) {
        VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

        // Initialise a new rearrangement list
        auto rearrangmentList = TreeRearrangment(vnode, min_radius, max_radius, true);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        rearrangmentList.defineMoves(false);

        // Print the list of moves for the current P node (source node)
        rearrangmentList.printMoves();

        VLOG(1) << "[tsh] Strategy " << rearrangmentList.mset_strategy << std::endl;
        VLOG(1) << "[utree rearrangment] Found " << rearrangmentList.getNumberOfMoves() << " possible moves for node " << vnode->vnode_name << std::endl;


        // For each potential move computed before, apply it to the tree topology, print the resulting newick tree, and revert it.
        for (unsigned long i = 0; i < rearrangmentList.getNumberOfMoves(); i++) {
            bool status;

            // Apply the move
            status = rearrangmentList.applyMove(i);

            //utree->saveTreeOnFile("../data/test.txt");

            if (status) {
                VLOG(2) << "[apply  move]\t" << rearrangmentList.getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                        << " | (" << rearrangmentList.getSourceNode()->vnode_name << "->" << rearrangmentList.getMove(i)->getTargetNode()->vnode_name << ")"
                        << "\t[" << rearrangmentList.getMove(i)->move_radius << "] | "
                        << input_tree->printTreeNewick(true) << std::endl;
                //utree->_testReachingPseudoRoot();
            }

            // Revert the move, and return to the original tree
            status = rearrangmentList.revertMove(i);
            //utree->saveTreeOnFile("../data/test.txt");
            if (status) {
                VLOG(2) << "[revert move]\t" << rearrangmentList.getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                        << " | (" << rearrangmentList.getMove(i)->getTargetNode()->vnode_name << "->" << rearrangmentList.getSourceNode()->vnode_name << ")"
                        << "\t[" << rearrangmentList.getMove(i)->move_radius << "] | "
                        << input_tree->printTreeNewick(true) << std::endl;
                //utree->_testReachingPseudoRoot();
            }
            total_exec_moves += rearrangmentList.getNumberOfMoves();
        }

    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    VLOG(1) << "Moves applied and reverted: " << total_exec_moves << std::endl;
    VLOG(1) << "Elapsed time: " << duration << " microseconds" << std::endl;
    VLOG(1) << "*** " << (double) duration/total_exec_moves << " microseconds/move *** " << std::endl;


}
