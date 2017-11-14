/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * tshlib is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tshlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with likpip. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Utree.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 26 10 2017
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

#include <string>
#include <random>
#include <fstream>

#include "Utree.hpp"


void ::UtreeUtils::_traverseTree(Utree *in_tree, VirtualNode *target, PhyTree *source) {

    for (unsigned long i = 0; i < source->n_children(); i++) {
        auto ichild = new VirtualNode();
        if (!source->get_children().at(i)->isLeaf()) {

            ichild->vnode_id = source->get_children().at(i)->getNodeID();
            ichild->vnode_name = source->get_children().at(i)->getName();
            ichild->vnode_branchlength = source->get_children().at(i)->getBranchLength();

            target->connectNode(ichild);
            in_tree->addMember(ichild);
            _traverseTree(in_tree, ichild, source->get_children().at(i));

        } else {

            ichild->vnode_id = source->get_children().at(i)->getNodeID();
            ichild->vnode_name = source->get_children().at(i)->getName();
            ichild->vnode_branchlength = source->get_children().at(i)->getBranchLength();
            // Set all the other directions to null
            ichild->setNodeLeft(nullptr);
            ichild->setNodeRight(nullptr);

            // Set the LEAF flag to true
            ichild->vnode_leaf = true;
            in_tree->addMember(ichild);
            target->connectNode(ichild);

        }
    }

}

void ::UtreeUtils::convertUtree(PhyTree *in_tree, Utree *out_tree) {

    // For each node descending the root, create either a new VirtualNode
    for (unsigned long i = 0; i < in_tree->n_children(); i++) {

        auto ichild = new VirtualNode;

        ichild->vnode_id = in_tree->get_children().at(i)->getNodeID();
        ichild->vnode_name = in_tree->get_children().at(i)->getName();
        ichild->vnode_branchlength = in_tree->get_children().at(i)->getBranchLength();

        // If the node in PhyTree is a leaf, skip the recursion
        if (in_tree->get_children().at(i)->isLeaf()) {

            ichild->vnode_leaf = true;

        } else {

            _traverseTree(out_tree, ichild, in_tree->get_children().at(i));

        }
        // Add this node as starting point of the tree
        out_tree->addMember(ichild, true);

    }

    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    out_tree->startVNodes.at(0)->setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->setNodeUp(out_tree->startVNodes.at(0));

}


VirtualNode::VirtualNode() {
    // Initialise a new VirtualNode completely disconnected from the tree
    this->vnode_right = nullptr;
    this->vnode_left = nullptr;
    this->vnode_up = nullptr;
    this->vnode_branchlength = 0;
    this->vnode_leaf = false;
    this->vnode_rotated = NodeRotation::undef;

};

void VirtualNode::connectNode(VirtualNode *inVNode) {

    // Connect this node to the parent node
    inVNode->setNodeUp(this);
    //inNode->_oneway_connectNode(this);
    this->_oneway_connectNode(inVNode);
}

bool VirtualNode::isTerminalNode() {

    return this->vnode_leaf;

}

void VirtualNode::setNodeRight(VirtualNode *inVNode) {

    this->vnode_right = inVNode;

}

void VirtualNode::setNodeLeft(VirtualNode *inVNode) {

    this->vnode_left = inVNode;

}

void VirtualNode::setNodeUp(VirtualNode *inVNode) {

    this->vnode_up = inVNode;

}

VirtualNode *VirtualNode::getNodeUp() {

    return this->vnode_up ?: nullptr;
}

void VirtualNode::rotateClockwise() {

    _recursive_cw_rotation(this, false);

}

void VirtualNode::rotateClockwise(bool revertRotations) {

    _recursive_cw_rotation(this, revertRotations);

}

void VirtualNode::rotateCounterClockwise() {

    _recursive_ccw_rotation(this, false);
}

void VirtualNode::rotateCounterClockwise(bool revertRotations) {

    _recursive_ccw_rotation(this, revertRotations);
}

VirtualNode *VirtualNode::getNodeLeft() {
    return this->vnode_left ?: nullptr;
}

VirtualNode *VirtualNode::getNodeRight() {
    return this->vnode_right ?: nullptr;
}

void Utree::addMember(VirtualNode *inVNode, bool isStartNode) {

    this->listVNodes.push_back(inVNode);
    if (isStartNode) {
        this->startVNodes.push_back(inVNode);
    }

}

std::vector<VirtualNode *> Utree::findPseudoRoot(VirtualNode *inVNode, bool fixPseudoRootOnNextSubtree) {

    // This method returns a vector which contains all the nodes in the path
    // from the start node to the pseudoroot
    std::vector<VirtualNode *> path2root;

    auto *CurrentNode = new VirtualNode;

    // Add the first element of the path to the root
    path2root.push_back(inVNode);

    CurrentNode = inVNode;

    // Add all the other nodes traversing the tree in post-order
    do {
        CurrentNode = CurrentNode->getNodeUp();
        path2root.push_back(CurrentNode);

    } while (CurrentNode != CurrentNode->getNodeUp()->getNodeUp());

    if (fixPseudoRootOnNextSubtree) {

        path2root.push_back(CurrentNode->getNodeUp());

    }

    return path2root;
}

std::string Utree::printTreeNewick(bool showInternalNodeNames) {
    std::string s, terminator;

    this->_updateStartNodes();

    s += "(";
    for (unsigned long i = 0; i < this->startVNodes.size(); i++) {
        if (i == this->startVNodes.size() - 1) {
            terminator = ");";
        } else {
            terminator = ",";
        }
        s += _recursiveFormatNewick(this->startVNodes.at(i), showInternalNodeNames) + terminator;
    }

    return s;
}

std::string Utree::_recursiveFormatNewick(VirtualNode *vnode, bool showInternalNodeNames) {

    std::stringstream newick;

    if (vnode->isTerminalNode()) {

        newick << vnode->vnode_name << ":" << vnode->vnode_branchlength;

    } else {

        newick << "(";
        newick << _recursiveFormatNewick(vnode->getNodeLeft(), showInternalNodeNames);
        newick << ",";
        newick << _recursiveFormatNewick(vnode->getNodeRight(), showInternalNodeNames);
        newick << ")";

        if (showInternalNodeNames) {
            newick << vnode->vnode_name;
        }

        newick << ":" << vnode->vnode_branchlength;
    }

    return newick.str();
}

Utree::~Utree() = default;

Utree::Utree() = default;

void Utree::_updateStartNodes() {

    // Find the pseudoroot node on the first side of the tree
    bool fixPseudoRootOnNextSubtree = true;
    std::vector<VirtualNode *> sideA = this->findPseudoRoot(this->listVNodes.at(0), fixPseudoRootOnNextSubtree);
    VirtualNode *sideA_node = sideA.back();

    // Find the pseudoroot node on the opposite side of the tree
    fixPseudoRootOnNextSubtree = false;
    std::vector<VirtualNode *> sideB = this->findPseudoRoot(this->listVNodes.at(0), fixPseudoRootOnNextSubtree);
    VirtualNode *sideB_node = sideB.back();

    // Reset the original utree attribute
    this->startVNodes.clear();

    // Update the utree attribute
    this->startVNodes.push_back(sideA_node);
    this->startVNodes.push_back(sideB_node);

}

void Utree::_testReachingPseudoRoot() {

    std::string strpath;
    std::vector<VirtualNode *> vnodes2root;

    for (auto &i : this->listVNodes) {

        vnodes2root = this->findPseudoRoot(i, true);

        for (auto &tnode: vnodes2root) {
            strpath += tnode->vnode_name;

            if (tnode->vnode_rotated != NodeRotation::undef) {

                strpath += "\033[1;34m(" + std::to_string(static_cast<int>(tnode->vnode_rotated)) + ")\033[0m";

            } else {

                strpath += "(" + std::to_string(static_cast<int>(tnode->vnode_rotated)) + ")";
            }

            strpath += ">";
        }
        strpath.pop_back();

        std::cout << strpath << std::endl;
        vnodes2root.clear();
        strpath.clear();

    }

}

void Utree::saveTreeOnFile(std::string outfilepath) {

    std::ofstream outfile;

    outfile.open(outfilepath, std::ios_base::app);
    outfile << this->printTreeNewick(true) << std::endl;

}

std::string VirtualNode::printNeighbours() {

    std::stringstream description;

    if (!this->isTerminalNode()) {

        description << this->vnode_name << " (^" << this->getNodeUp()->vnode_name << ";";
        description << "<" << this->getNodeLeft()->vnode_name << ";";
        description << this->getNodeRight()->vnode_name << ">)";

    } else {

        description << this->vnode_name << " (^" << this->getNodeUp()->vnode_name << "; ";
        description << "<-;->)";
    }

    return description.str();
}

NodePosition VirtualNode::indexOf() {
    NodePosition node_position;
    VirtualNode *parent = this->getNodeUp();
    if (parent->getNodeLeft() == this) {

        node_position = NodePosition::left;

    } else if (parent->getNodeRight() == this) {

        node_position = NodePosition::right;

    } else if (parent->getNodeUp() == this) {

        node_position = NodePosition::up;

    } else {
        //perror("The pointer index of the current node is unknown");
        node_position = NodePosition::undef;
    }
    return node_position;
}

bool VirtualNode::swapNode(VirtualNode *targetNode, MoveDirections move_direction, bool revertRotations) {

    bool execstatus = false;

    if (this == targetNode) {

        perror("The source and the target nodes must be different!");
        exit(EXIT_FAILURE);

    } else if (targetNode == nullptr) {

        perror("The target node is empty");


    } else {

        VirtualNode *pnode, *qnode, *pnode_parent, *qnode_parent;

        pnode = this;
        qnode = targetNode;
        pnode_parent = this->getNodeUp();
        qnode_parent = targetNode->getNodeUp();

/*
        std::cout << pnode->printNeighbours() << std::endl;
        std::cout << qnode->printNeighbours() << std::endl;
        std::cout << pnode_parent->printNeighbours() << std::endl;
        std::cout << qnode_parent->printNeighbours() << std::endl;
*/
        bool invertQParents = false;
        bool invertPParents = false;

        switch (move_direction) {

            case MoveDirections::left:
                pnode->rotateCounterClockwise();
                pnode_parent = pnode->getNodeUp();

                pnode->disconnectNode();
                qnode->disconnectNode();

                qnode_parent->connectNode(pnode);
                pnode_parent->connectNode(qnode);

                break;

            case MoveDirections::right:

                pnode->rotateClockwise();
                pnode_parent = pnode->getNodeUp();

                pnode->disconnectNode();
                qnode->disconnectNode();

                qnode_parent->connectNode(pnode);
                pnode_parent->connectNode(qnode);

                break;

            case MoveDirections::up_left:

                qnode->rotateCounterClockwise();
                qnode_parent = qnode->getNodeUp();
                qnode->disconnectNode();

                pnode->disconnectNode();
                pnode_parent->connectNode(qnode);

                qnode_parent->connectNode(pnode);

                break;

            case MoveDirections::up_right:

                qnode->rotateClockwise();
                qnode_parent = qnode->getNodeUp();
                qnode->disconnectNode();

                pnode->disconnectNode();
                pnode_parent->connectNode(qnode);

                qnode_parent->connectNode(pnode);

                break;

            case MoveDirections::up:

                // Prune nodes
                pnode->disconnectNode();
                qnode->disconnectNode();

                // Save value of rotation
                if (qnode->vnode_rotated != NodeRotation::undef) {
                    invertQParents = true;
                }

                if (pnode->vnode_rotated != NodeRotation::undef) {
                    invertPParents = true;
                }

                // Reset rotation
                pnode->resetNodeDirections(revertRotations);
                qnode->resetNodeDirections(revertRotations);

                // Regraft nodes
                if (invertQParents) {

                    qnode->connectNode(pnode_parent);

                } else {

                    pnode_parent->connectNode(qnode);
                }

                if (invertPParents) {

                    pnode->connectNode(qnode_parent);

                } else {

                    qnode_parent->connectNode(pnode);

                }


                break;

            default:

                std::cout << "I cannot move the node" << std::endl;
                exit(EXIT_FAILURE);
        }

        execstatus = true;
    }

    return execstatus;
}

void VirtualNode::disconnectNode() {

    // 1. Get the location where this node is connected on the parent node and disconnect
    switch (this->indexOf()) {
        case NodePosition::left:
            // The node is connected on the left side
            this->getNodeUp()->setNodeLeft(nullptr);
            break;
        case NodePosition::right:
            // The node is connected on the right side
            this->getNodeUp()->setNodeRight(nullptr);
            break;
        case NodePosition::up:
            // The node is connect on the up side (pseudoroot)
            this->getNodeUp()->setNodeUp(nullptr);
            break;
        default:
            break;
    }

    // 2. Disconnect the parent node
    this->setNodeUp(nullptr);
}

void VirtualNode::_oneway_connectNode(VirtualNode *inVNode) {

    // Check which direction is still available (either left or right)
    if (!this->getNodeLeft()) {
        this->setNodeLeft(inVNode);
        return;
    } else if (!this->getNodeRight()) {
        this->setNodeRight(inVNode);
        return;
    } else if (!this->getNodeUp()) {
        this->setNodeUp(inVNode);
        return;
    } else {
        perror("No direction available in the VirtualNode to add a new child");
    }

}

bool VirtualNode::isParent(VirtualNode *inVNode) {

    bool status = false;

    VirtualNode *CurrentNode = inVNode;

    do {
        CurrentNode = CurrentNode->getNodeUp();
        if (CurrentNode == this) {
            status = true;
            break;
        };
    } while (CurrentNode != CurrentNode->getNodeUp()->getNodeUp());

    return status;
}

void VirtualNode::resetNodeDirections(bool revertRotations) {

    if (this->vnode_rotated != NodeRotation::undef) {
        switch (this->vnode_rotated) {
            case NodeRotation::clockwise:
                this->rotateCounterClockwise(revertRotations);
                break;
            case NodeRotation::counterclockwise:
                this->rotateClockwise(revertRotations);
                break;
            default:
                break;
        }
    }

}

void VirtualNode::_recursive_cw_rotation(VirtualNode *vnode, bool revertRotations) {

    bool continueTraverse = true;
    NodePosition curr_node_position = NodePosition::undef;

    if (revertRotations and vnode->vnode_rotated == NodeRotation::undef) { return; }

    if (revertRotations) {

        continueTraverse = true;

    } else {
        if (vnode->getNodeUp() != nullptr) {
            if (vnode->getNodeUp()->getNodeUp() == vnode) {
                continueTraverse = false;
            }
        }
        curr_node_position = vnode->indexOf();
    }

    // rotate the pointers only if the node is not a leaf
    if (!vnode->isTerminalNode()) {
        auto curr_vn_up = new VirtualNode;
        auto curr_vn_left = new VirtualNode;
        auto curr_vn_right = new VirtualNode;
        auto n_up = new VirtualNode;


        curr_vn_up = vnode->getNodeUp();
        curr_vn_left = vnode->getNodeLeft();
        curr_vn_right = vnode->getNodeRight();

        vnode->setNodeUp(curr_vn_right);
        vnode->setNodeLeft(curr_vn_up);
        vnode->setNodeRight(curr_vn_left);

        // Store the information about the rotation
        if (vnode->vnode_rotated != NodeRotation::undef) {
            vnode->vnode_rotated = NodeRotation::undef;
        } else {
            vnode->vnode_rotated = NodeRotation::clockwise;
        }

        if (revertRotations) {
            n_up = vnode->getNodeUp();

            switch (n_up->vnode_rotated) {
                case NodeRotation::clockwise:
                    curr_node_position = NodePosition::left;
                    break;
                case NodeRotation::counterclockwise:
                    curr_node_position = NodePosition::right;
                    break;
                default:
                    curr_node_position = NodePosition::undef;
            }

        } else {
            n_up = curr_vn_up;
        }

        if (continueTraverse) {
            if (n_up != nullptr) {
                switch (curr_node_position) {
                    case NodePosition::left:
                        _recursive_ccw_rotation(n_up, revertRotations);
                        break;
                    case NodePosition::right:
                        _recursive_cw_rotation(n_up, revertRotations);
                        break;
                    case NodePosition::up:
                        _recursive_cw_rotation(n_up, revertRotations);
                    default:
                        return;
                }
            }


        } else {

            //std::cout << "This node is the pseudoroot" << std::endl;

        }
    }

}

void VirtualNode::_recursive_ccw_rotation(VirtualNode *vnode, bool revertRotations) {

    bool continueTraverse = true;
    NodePosition curr_node_position = NodePosition::undef;

    if (revertRotations and vnode->vnode_rotated == NodeRotation::undef) { return; }

    if (revertRotations) {

        continueTraverse = true;

    } else {
        if (vnode->getNodeUp() != nullptr) {
            if (vnode->getNodeUp()->getNodeUp() == vnode) {

                continueTraverse = false;
            }
        }
        curr_node_position = vnode->indexOf();
    }


    // rotate the pointers only if the node is not a leaf
    if (!vnode->isTerminalNode()) {

        auto curr_vn_up = new VirtualNode;
        auto curr_vn_left = new VirtualNode;
        auto curr_vn_right = new VirtualNode;

        auto n_up = new VirtualNode;

        curr_vn_up = vnode->getNodeUp();
        curr_vn_left = vnode->getNodeLeft();
        curr_vn_right = vnode->getNodeRight();

        // Revert pointers
        vnode->setNodeUp(curr_vn_left);
        vnode->setNodeLeft(curr_vn_right);
        vnode->setNodeRight(curr_vn_up);

        // Store the information about the rotation
        if (vnode->vnode_rotated != NodeRotation::undef) {
            vnode->vnode_rotated = NodeRotation::undef;
        } else {
            vnode->vnode_rotated = NodeRotation::counterclockwise;
        }

        if (revertRotations) {
            n_up = vnode->getNodeUp();

            switch (n_up->vnode_rotated) {
                case NodeRotation::clockwise: //clockwise
                    curr_node_position = NodePosition::left; // force to rotate counterclockwise
                    break;
                case NodeRotation::counterclockwise: //counterclockwise
                    curr_node_position = NodePosition::right; // force to rotate clockwise
                    break;
                default:
                    curr_node_position = NodePosition::undef;
            }

        } else {
            n_up = curr_vn_up;
        }

        if (continueTraverse) {
            if (n_up != nullptr) {
                switch (curr_node_position) {
                    case NodePosition::left:
                        _recursive_ccw_rotation(n_up, revertRotations);
                        break;
                    case NodePosition::right:
                        _recursive_cw_rotation(n_up, revertRotations);
                        break;
                    case NodePosition::up:
                        _recursive_ccw_rotation(n_up, revertRotations);
                    default:
                        return;
                }
            }

        } else {

            //std::cout << "This node is the pseudoroot" << std::endl;

        }

    }

}

VirtualNode::~VirtualNode() = default;

