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
#include <limits>
#include <iomanip>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <glog/logging.h>
#include <map>
#include "TreeRearrangment.hpp"


#ifdef TSHLIB_BENCHMARK
#include <iostream>
#include <fstream>
#endif

namespace tshlib {
    TreeRearrangment::TreeRearrangment() = default;

    TreeRearrangment::~TreeRearrangment() {

        for (std::vector<Move *>::reverse_iterator i = trMoveSet.rbegin(); i < trMoveSet.rend(); i++) {
            Move *move = *i;
            delete move;
        }
        std::vector<Move *>().swap(trMoveSet);


    };

    void TreeRearrangment::getNodesInRadiusDown(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, bool includeSelf, MoveDirections direction, bool allowDuplicatedMoves) {

        VirtualNode *node = node_source;

        if (!includeSelf) {
            includeSelf = true;

        } else {
            //if (radius_max == 0 ) {
            if (radius_curr <= (radius_max - radius_min) && radius_curr >= 0) {
                auto moveInstance = new Move;
                moveInstance->initMove();
                moveInstance->setSourceNode(trSourceNode_);
                moveInstance->setDirection(direction);
                moveInstance->setRadius(radius_max - radius_curr);
                moveInstance->setTargetNode(node);
                moveInstance->setClass(trStrategy);
                addMove(moveInstance, allowDuplicatedMoves);
            }
        }

        if (radius_curr < 0) {
            return;
        }

        if (!node->isTerminalNode()) {
            radius_curr--;
            // Direction is required to know where to look at in the tree when the move is performed
            getNodesInRadiusDown(node->getNodeLeft(), radius_min, radius_curr, radius_max, includeSelf, direction, allowDuplicatedMoves);
            getNodesInRadiusDown(node->getNodeRight(), radius_min, radius_curr, radius_max, includeSelf, direction, allowDuplicatedMoves);
        }


    }

    void TreeRearrangment::getNodesInRadiusUp(VirtualNode *node_source, int radius_min, int radius_curr, int radius_max, NodePosition traverse_direction, bool allowDuplicatedMoves) {

        VirtualNode *vnode;

        NodePosition idx;
        MoveDirections moving_direction;


        vnode = node_source;

        //TODO: check binary tree condition!
        if (radius_curr <= (radius_max - radius_min) && radius_curr >= 0) {

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
            auto moveInstance = new Move;
            moveInstance->initMove();
            moveInstance->setSourceNode(trSourceNode_);
            moveInstance->setDirection(moving_direction);
            moveInstance->setRadius(radius_max - radius_curr);
            moveInstance->setTargetNode(vnode);
            moveInstance->setClass(trStrategy);
            addMove(moveInstance, allowDuplicatedMoves);
        }

        if (radius_curr > 0) {
            if (traverse_direction == NodePosition::left) {
                radius_curr--;
                if (vnode->getNodeUp() != nullptr) {

                    idx = vnode->indexOf();
                    getNodesInRadiusUp(vnode->getNodeUp(), radius_min, radius_curr, radius_max, idx, allowDuplicatedMoves);
                }

                getNodesInRadiusDown(vnode->getNodeRight(), radius_min, radius_curr, radius_max, true, MoveDirections::up, allowDuplicatedMoves);

            } else if (traverse_direction == NodePosition::right) {
                radius_curr--;
                if (vnode->getNodeUp() != nullptr) {

                    idx = vnode->indexOf();
                    getNodesInRadiusUp(vnode->getNodeUp(), radius_min, radius_curr, radius_max, idx, allowDuplicatedMoves);
                }

                getNodesInRadiusDown(vnode->getNodeLeft(), radius_min, radius_curr, radius_max, true, MoveDirections::up, allowDuplicatedMoves);

            } else if (traverse_direction == NodePosition::up) {

                // If the traverse_direction is 2, we are moving across the pseudoroot, therefore no need to decrease the radious
                if (vnode->getNodeUp() != nullptr) {

                    // Get the nodes in the radius from the node we reached after crossing the pseudoroot
                    getNodesInRadiusDown(vnode, radius_min, radius_curr - 1, radius_max - 1, false, MoveDirections::up, allowDuplicatedMoves);

                }

            }
        }
    }

    void TreeRearrangment::defineMoves(bool includeSelf, bool allowDuplicatedMoves) {
        // Flag the nodes according to their position on the tree (left or right or above the source node -- p node).
        // For each node within the radius extremities, define a move and add it to TreeRearrangment.
        // Start from the children of the current starting node (if any)

        trCandidateMovesFound_ = 0;

        if (!trSourceNode_->isTerminalNode()) {
            getNodesInRadiusDown(trSourceNode_->getNodeLeft(), trSearchRadius_min, trSearchRadius_max - 1, trSearchRadius_max, includeSelf, MoveDirections::down_left, allowDuplicatedMoves);

            getNodesInRadiusDown(trSourceNode_->getNodeRight(), trSearchRadius_min, trSearchRadius_max - 1, trSearchRadius_max, includeSelf, MoveDirections::down_right, allowDuplicatedMoves);
        }
        // If the node is a leaf, then go up
        if (nullptr != trSourceNode_->getNodeUp()) {
            getNodesInRadiusUp(trSourceNode_->getNodeUp(), trSearchRadius_min, trSearchRadius_max - 1, trSearchRadius_max, trSourceNode_->indexOf(), allowDuplicatedMoves);
        }


        VLOG(1) << "[TSH Cycle]  Found " << trCandidateMovesFound_ << " candidate moves for node " << trSourceNode_->getNodeName();
    }

    void TreeRearrangment::addMove(Move *move, bool allowDuplicatedMoves) {

        bool storeMove = true;

        if (!allowDuplicatedMoves) {
            for (auto &query:trMoveSet) {

                if (query->getTargetNode() == move->getTargetNode() && query->getSourceNode() == move->getSourceNode()) {
                    storeMove = false;
                }

                if (move->getTargetNode() == query->getSourceNode() && move->getSourceNode() == query->getTargetNode()) {
                    storeMove = false;
                }

            }
        }

        if (storeMove) {
            move->moveUID_ = (int) trMoveSet.size();
            trMoveSet.push_back(move);
            trCandidateMovesFound_++;
        }

    }

    void TreeRearrangment::printMoves() {

        VLOG(2) << "[set " << trUID_ << "] " << getStrategy() << " strategy" << std::endl;
        VLOG(2) << "[class]\t(P\t; Q)" << std::endl;
        for (auto &nmove: trMoveSet) {
            VLOG(2) << "[" << nmove->moveClassDescription_ << "]\t" << nmove->moveRadius_
                    << "\t(" << trSourceNode_->vnode_name << "\t; "
                    << nmove->getTargetNode()->vnode_name << ")\t"
                    << static_cast<int>(nmove->moveDirection_) << std::endl;
        }

    }

    unsigned long TreeRearrangment::getNumberOfMoves() {

        return trMoveSet.size();
    }

    bool TreeRearrangment::applyMove(unsigned long moveID) {

        bool outcomeExecutionMove = false;

        VirtualNode *pnode = trMoveSet.at(moveID)->getSourceNode();
        VirtualNode *qnode = trMoveSet.at(moveID)->getTargetNode();

        VLOG(2) << "[tshlib::applyMove] Debug: [MOVE#" << trMoveSet.at(moveID)->getUID() << "] direction: " << trMoveSet.at(moveID)->getDirection();
        VLOG(2) << "[tshlib::applyMove] Debug: [MOVE#" << trMoveSet.at(moveID)->getUID() << "] S: " << pnode->getNodeName() << " T: " << qnode->getNodeName();


        bool revertRotations;

        switch(trMoveSet.at(moveID)->getType()){

            case MoveType::VFNNI:

            case MoveType::FNNI:

            case MoveType::NNI:

                revertRotations = false;
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove =  pnode->swapNode(qnode, trMoveSet.at(moveID)->moveDirection_, revertRotations);

                break;

            case MoveType::SPR:
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove = _applySPR(trMoveSet.at(moveID));

                break;

            case MoveType::TBR:

                LOG(WARNING) << "[tshlib::TreeRearrangment::applyMove] Move [" << moveID << "] is TBR. Not implemented method to apply it.";

                break;



            case MoveType::undef:

                LOG(FATAL) << "[tshlib::TreeRearrangment::applyMove] Something went wrong during the application of the move [" << moveID << "]. It looks like its type is undefined!";

                break;
        }

        return outcomeExecutionMove;

    }

    bool TreeRearrangment::revertMove(unsigned long moveID) {

        bool outcomeExecutionMove = false;

        VirtualNode *pnode = trMoveSet.at(moveID)->getTargetNode();
        VirtualNode *qnode = trMoveSet.at(moveID)->getSourceNode();

        //bool revertRotations = true;
        //return pnode->swapNode(qnode, MoveDirections::up, revertRotations);

        bool revertRotations;

        switch(trMoveSet.at(moveID)->getType()){

            case MoveType::VFNNI:

            case MoveType::FNNI:

            case MoveType::NNI:
                revertRotations = true;
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove =  pnode->swapNode(qnode, MoveDirections::up, revertRotations);

                break;

            case MoveType::SPR:
                // Swap pnode with qnode according to the direction found during the move configuration
                // If the swap is performed correctly then the function returns true otherwise false
                outcomeExecutionMove = _revertSPR(trMoveSet.at(moveID));

                break;

            case MoveType::TBR:

                LOG(WARNING) << "[tshlib::TreeRearrangment::applyMove] Move [" << moveID << "] is TBR. Not implemented method to revert it.";

                break;



            case MoveType::undef:

                LOG(FATAL) << "[tshlib::TreeRearrangment::applyMove] Something went wrong during the application of the move [" << moveID << "]. It looks like its type is undefined!";

                break;
        }

        return outcomeExecutionMove;


    }

    Move *TreeRearrangment::getMove(unsigned long moveID) {

        return trMoveSet.at(moveID) ?: nullptr;
    }

    Move *TreeRearrangment::selectBestMove(double value = -std::numeric_limits<double>::infinity()) {

        Move *selectedMove = nullptr;
        for (auto &move:trMoveSet) {
            if (move->moveScore_ > value) {
                selectedMove = move;
                value = move->moveScore_;
            }
        }

        return selectedMove;

    }

    void TreeRearrangment::commitMove(int moveID) {

        // apply move
        applyMove(moveID);

        // reset node rotations
        for (auto &node:UTree_->listVNodes) {
            node->vnode_rotated = NodeRotation::undef;

        }

    }

    void TreeRearrangment::storeMove(Move *inMove) {

        inMove->moveUID_ = (int) trMoveSet.size();
        trMoveSet.push_back(inMove);

    }

    void TreeRearrangment::setTreeTopology(Utree *inTree) {
        UTree_ = inTree;
    }

    void TreeRearrangment::displayRearrangmentStatus(int idMove, bool printTree) {

        std::string start_col_line, end_col_line;

        // ------------------------------------
        // Some abbellishments for the console output
        if (getMove(idMove)->moveScore_ > 0) {
            start_col_line = "\033[1;34m";
            end_col_line = "\033[0m";
        } else {
            start_col_line = "";
            end_col_line = "";

        }

        // ------------------------------------
        // Move exection details
        if (printTree) {
            VLOG(2) << "[test  move]\t" << getMove(idMove)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << idMove
                    << " | " << start_col_line << getMove(idMove)->moveScore_ << end_col_line << "\t"
                    << " | (" << getMove(idMove)->getSourceNode()->vnode_name << "->" << getMove(idMove)->getTargetNode()->vnode_name << ")"
                    << "\t[" << getMove(idMove)->moveRadius_ << "] | " << getTree()->printTreeNewick(true);
        } else {
            VLOG(2) << "[test  move]\t" << getMove(idMove)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << idMove
                    << " | " << start_col_line << getMove(idMove)->moveScore_ << end_col_line << "\t"
                    << " | (" << getMove(idMove)->getSourceNode()->vnode_name << "->" << getMove(idMove)->getTargetNode()->vnode_name << ")"
                    << "\t[" << getMove(idMove)->moveRadius_ << "]";
        }

    }

    const std::vector<VirtualNode *> TreeRearrangment::updatePathBetweenNodes(unsigned long moveID, std::vector<VirtualNode *> inPath) {

        Move *move = getMove(moveID);

        if (move->moveDirection_ != MoveDirections::up) {

            std::vector<VirtualNode *> tmpVector_B, tmpVector_C, updatedNodesInPath;
            tmpVector_B = inPath;

            // Remove the first element of the array since it is not a likelihood component (source node)
            tmpVector_B.erase(tmpVector_B.begin());

            // Find the position of the target node (it is not necessarely the end of the vector)
            std::ptrdiff_t pos;

            if (move->moveDirection_ == MoveDirections::up_left || move->moveDirection_ == MoveDirections::up_right) {
                pos = std::distance(tmpVector_B.begin(), std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getTargetNode()));
            } else {
                pos = std::distance(tmpVector_B.begin(), std::find(tmpVector_B.begin(), tmpVector_B.end(), move->getSourceNode()));
            }

            if (pos <= tmpVector_B.size()) {

                // Copy the reference to the pointers of the node starting from the target node to the end of the vector (root)
                for (std::ptrdiff_t i = pos; i < tmpVector_B.size(); i++) {
                    tmpVector_C.push_back(tmpVector_B.at(i));
                }

                // Revert the order of the elements in the vector
                std::reverse(std::begin(tmpVector_C), std::end(tmpVector_C));

                // Add the beginning of the vector B to vector C
                for (std::ptrdiff_t i = 0; i < pos; i++) {
                    tmpVector_C.push_back(tmpVector_B.at(i));
                }

                // Return the newly ordered vector
                updatedNodesInPath = tmpVector_C;
            }

            return updatedNodesInPath;

        } else {
            return inPath;
        }
    }

    void TreeRearrangment::initialize() {

        // Define the radius for pruning and regrafting the input tree.
        switch (trTreeCoverage) {

            case tshlib::TreeRearrangmentOperations::classic_NNI:
                trSearchRadius_min = 3;
                trSearchRadius_max = 3;
                break;

            case tshlib::TreeRearrangmentOperations::classic_SPR:
                trSearchRadius_min = 4;
                trSearchRadius_max = UTree_->getMaxNodeDistance() / 2;
                break;

            case tshlib::TreeRearrangmentOperations::classic_TBR:
                trSearchRadius_min = 5;
                trSearchRadius_max = UTree_->getMaxNodeDistance() / 2;
                break;

            case tshlib::TreeRearrangmentOperations::classic_Mixed:
                trSearchRadius_min = 3;  // Minimum radius for an NNI move is 3 nodes
                trSearchRadius_max = UTree_->getMaxNodeDistance(); // Full tree traversing from any nodeInterface of the tree
                break;

        }

        VLOG(1) << "[tsh] Defined coverage radius as ["<<trSearchRadius_min<<";"<<trSearchRadius_max<<"] on a max of ["<<UTree_->getMaxNodeDistance()<<"]";

        trInitialized_=true;
    }

    bool TreeRearrangment::_applySPR(Move *move) {

        bool execstatus = false;

        VirtualNode *sourceNode = move->getSourceNode();
        VirtualNode *targetNode = move->getTargetNode();


        if (sourceNode == targetNode) {

            LOG(FATAL) << "[tshlib::_applySPR] The source and the target nodes must be different!";

        } else if (targetNode == nullptr) {

            LOG(FATAL) << "[tshlib::_applySPR] The target node is empty";


        } else {

            // References
            VirtualNode *moveStepChild = nullptr;
            VirtualNode *moveStepParent = nullptr;
            VirtualNode *parentTarget = targetNode->getNodeUp();
            VirtualNode *parentSource = sourceNode->getNodeUp();

            // Assign flags to nodes

            switch (move->getMoveDirection()) {
                case MoveDirections::up:
                case MoveDirections::up_left:
                case MoveDirections::up_right:

                    // Assign flags to nodes
                    moveStepChild = sourceNode->getSiblingNode();
                    moveStepParent = parentSource->getNodeUp();
                    // Disconnections
                    parentSource->disconnectNode();
                    moveStepChild->disconnectNode();
                    targetNode->disconnectNode();

                    break;

                case MoveDirections::down_left:
                case MoveDirections::down_right:

                    // Rotations //TODO: Check for rotation direction (it could be the way around)
                    if(move->getMoveDirection() == MoveDirections::down_left){
                        sourceNode->rotateCounterClockwise();
                    }else{
                        sourceNode->rotateClockwise();
                    }
                    // Resolve parent-source node
                    parentSource = sourceNode->getNodeUp();

                    // Assign flags to nodes
                    moveStepChild = parentSource->getNodeRight();
                    moveStepParent = parentSource->getNodeLeft();
                    // Disconnections
                    moveStepChild->disconnectNode();
                    moveStepParent->disconnectNode();
                    targetNode->disconnectNode();
                    // Re-solve
                    VirtualNodeUtils::rotateNodeClockwise(parentSource);

                    break;

                default:
                    LOG(ERROR) << "[thslib::_applySPR] Move #" << move->getUID() << " has an undefined direction.";
                    break;
            }

            // Re-Connections
            moveStepParent->connectNode(moveStepChild);
            parentTarget->connectNode(parentSource);
            parentSource->connectNode(targetNode);

            VLOG(2) << "[tshlib::_applySPR] Debug: [MOVE#" << move->getUID() << "] pS: " << parentSource->getNodeName() << " pT: " << parentTarget->getNodeName();
            VLOG(2) << "[tshlib::_applySPR] Debug: [MOVE#" << move->getUID() << "] stepParent: " << moveStepParent->getNodeName() << " stepChild: " << moveStepChild->getNodeName();


            // Set StepParent and StepChild in move description
            move->setStepChildNode(moveStepChild);
            move->setStepParentNode(moveStepParent);

            execstatus = true;

        }
        return execstatus;
    }

    bool TreeRearrangment::_revertSPR(Move *move) {

        bool execstatus = false;

        VirtualNode *sourceNode = move->getSourceNode();
        VirtualNode *parentSourceNode = move->getSourceNode()->getNodeUp();

        VirtualNode *targetNode =  move->getTargetNode();
        VirtualNode *parentTargetNode =  move->getSourceNode()->getNodeUp()->getNodeUp();

        VirtualNode *moveStepParent = move->getStepParentNode();
        VirtualNode *moveStepChild = move->getStepChildNode();


        if (sourceNode == targetNode) {

            LOG(FATAL) << "[tshlib::sprNode_revert] The source and the target nodes must be different!";

        } else if (targetNode == nullptr) {

            LOG(FATAL) << "[tshlib::sprNode_revert] The target node is empty";


        } else {

            // 1. Disconnections
            moveStepChild->disconnectNode();
            parentSourceNode->disconnectNode();
            targetNode->disconnectNode();

            // 2. Revert rotations on the parentSourceNode (resolve)
            bool rotatedcase = false;
            switch(sourceNode->getNodeRotation()){
                case NodeRotation::counterclockwise:
                    sourceNode->rotateClockwise(true);
                    VirtualNodeUtils::rotateNodeClockwise(parentSourceNode);
                    rotatedcase = true;
                    break;
                case NodeRotation::clockwise:
                    sourceNode->rotateCounterClockwise(true);
                    VirtualNodeUtils::rotateNodeCounterClockwise(parentSourceNode);
                    rotatedcase = true;
                    break;
                case NodeRotation::undef:
                    break;
            }

            // 4. Reconnections
            parentSourceNode->connectNode(moveStepChild);

            if(rotatedcase){
                parentSourceNode->connectNode(moveStepParent);
            }else{
                moveStepParent->connectNode(parentSourceNode);
            }

            parentTargetNode->connectNode(targetNode);

            execstatus = true;
        }
        return execstatus;
    }

}


//    void TreeRearrangment::initTreeRearrangment(VirtualNode *node_source, int radius, bool preserve_blengths) {

//        mset_sourcenode = node_source;
//        mset_id = node_source->vnode_name + ":" + std::to_string(radius);
//        mset_min_radius = radius;
//        mset_max_radius = radius;
//        mset_preserve_blenghts = preserve_blengths;
//        mset_strategy = "undefined";

//    }


//    void TreeRearrangment::initTreeRearrangment(Utree *ref_tree, int min_radius, int max_radius, bool preserve_blengths, VirtualNode *node_source) {

//        tree = ref_tree;
//        mset_sourcenode = node_source;
//        mset_id = node_source->vnode_name + ":" + std::to_string(min_radius) + "-" + std::to_string(max_radius);
//        mset_min_radius = min_radius;
//        mset_max_radius = max_radius;
//        mset_preserve_blenghts = preserve_blengths;
//
//        if (min_radius == 3 && max_radius == 3) {
//            mset_strategy = "standard NNI";
//        }
//
//        if (min_radius == 3 && max_radius > 3) {
//            mset_strategy = "mixed (NNI+SPR+TBR)";
//        }
//
//        if (min_radius == 4 && (max_radius > 4 && max_radius < 10)) {
//            mset_strategy = "standard SPR";
//        }


//    }



//void ::treesearchheuristics::performTestTreeSearch(Utree *input_tree, TreeSearchHeuristics tsh_strategy) {
//
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//    int min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
//    int max_radius = input_tree->getMaxNodeDistance(); // Hard coded max value for a small tree (this ensures the complete q-node search)
//
//    unsigned long total_exec_moves = 0;
//
//    // Print node description with neighbors
//    for (auto &vnode:input_tree->listVNodes) {
//        //    VirtualNode *vnode = input_tree->listVNodes.at(200);
//
//        VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;
//
//        // Initialise a new rearrangement list
//        auto rearrangmentList = new TreeRearrangment;
//
//        rearrangmentList->initTreeRearrangment(nullptr, min_radius, max_radius, true, vnode);
//
//        // Get all the target nodes with distance == radius from the source node
//        // excluding the starting node.
//        rearrangmentList->defineMoves(false, false);
//
//        // Print the list of moves for the current P node (source node)
//        rearrangmentList->printMoves();
//
//        VLOG(1) << "[tsh] Strategy " << rearrangmentList->getStrategy() << std::endl;
//        VLOG(1) << "[utree rearrangment] Found " << rearrangmentList->getNumberOfMoves() << " possible moves for node " << vnode->vnode_name << std::endl;
//
//        //
//
//        // For each potential move computed before, apply it to the tree topology, print the resulting newick tree, and revert it.
//        for (unsigned long i = 0; i < rearrangmentList->getNumberOfMoves(); i++) {
//            bool status;
//
//            // Apply the move
//            status = rearrangmentList->applyMove(i);
//
//            //utree->saveTreeOnFile("../data/test.txt");
//
//            if (status) {
//                VLOG(2) << "[apply  move]\t" << rearrangmentList->getMove(i)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << i
//                        << " | (" << rearrangmentList->getSourceNode()->vnode_name << "->" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << ")"
//                        << "\t[" << rearrangmentList->getMove(i)->moveRadius_ << "] | "
//                        << input_tree->printTreeNewick(true) << std::endl;
//                //utree->_testReachingPseudoRoot();
//            }
//
//            // Revert the move, and return to the original tree
//            status = rearrangmentList->revertMove(i);
//            //utree->saveTreeOnFile("../data/test.txt");
//            if (status) {
//                VLOG(2) << "[revert move]\t" << rearrangmentList->getMove(i)->moveClassDescription_ << "." << std::setfill('0') << std::setw(3) << i
//                        << " | (" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << "->" << rearrangmentList->getSourceNode()->vnode_name << ")"
//                        << "\t[" << rearrangmentList->getMove(i)->moveRadius_ << "] | "
//                        << input_tree->printTreeNewick(true) << std::endl;
//                //utree->_testReachingPseudoRoot();
//            }
//
//        }
//
//        total_exec_moves += rearrangmentList->getNumberOfMoves() * 2;
//        delete rearrangmentList;
//    }
//
//    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
//
//    VLOG(2) << "Moves applied and reverted: " << total_exec_moves << std::endl;
//    VLOG(2) << "Elapsed time: " << duration << " microseconds" << std::endl;
//    VLOG(2) << "*** " << (double) duration / total_exec_moves << " microseconds/move *** " << std::endl;
//
//
//#ifdef TSHLIB_BENCHMARK
//    std::ofstream myfile;
//    myfile.open ("tshlib_statistics.csv", std::ios_base::app);
//    myfile << input_tree->listVNodes.size()<< "," <<total_exec_moves << "," << duration << "," << (double) duration/total_exec_moves << std::endl;
//    myfile.close();
//#endif
//
//}
//


