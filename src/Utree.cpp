//
// Created by Lorenzo Gatti on 26/10/17.
//

#include <string>
#include <random>

#include "Utree.hpp"
#include "TreeRearrangment.hpp"


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
    for (int i = 0; i < in_tree->n_children(); i++) {

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


    /*std::default_random_engine generator;
    std::uniform_int_distribution<unsigned long> distribution(0, out_tree->listVNodes.size() - 1);
    unsigned long randomIndex = distribution(generator);


    for (unsigned long i = 0; i < out_tree->listVNodes.size(); i++) {
        if (i != randomIndex) {
            out_tree->listVNodes.at(i)->setNodeUp(out_tree->listVNodes.at(randomIndex));
        } else {
            if (randomIndex == 0) {
                randomIndex = randomIndex + 1;
            } else {
                randomIndex = randomIndex - 1;
            }
            out_tree->listVNodes.at(i)->setNodeUp(out_tree->listVNodes.at(randomIndex));
        }
    }*/

}


VirtualNode::VirtualNode() {
    // Initialise a new VirtualNode completely disconnected from the tree
    this->vnode_right = nullptr;
    this->vnode_left = nullptr;
    this->vnode_up = nullptr;
    this->vnode_branchlength = 0;
    this->vnode_leaf = false;
    this->vnode_rotated = 0;

};

void VirtualNode::connectNode(VirtualNode *inNode) {

    // Connect this node to the parent node
    inNode->setNodeUp(this);
    //inNode->_oneway_connectNode(this);
    this->_oneway_connectNode(inNode);
}

bool VirtualNode::isTerminalNode() {

    return this->vnode_leaf;

}

void VirtualNode::setNodeRight(VirtualNode *inNode) {

    this->vnode_right = inNode;

}

void VirtualNode::setNodeLeft(VirtualNode *inNode) {

    this->vnode_left = inNode;

}

void VirtualNode::setNodeUp(VirtualNode *inNode) {

    this->vnode_up = inNode;

}

VirtualNode *VirtualNode::getNodeUp() {

    return this->vnode_up ?: nullptr;
}

void VirtualNode::RotateClockwise() {

    _recursive_cw_rotation(this, false);

}

void VirtualNode::RotateClockwise(bool revertRotations) {

    _recursive_cw_rotation(this, revertRotations);

}

void VirtualNode::RotateCounterClockwise() {

    _recursive_ccw_rotation(this, false);
}

void VirtualNode::RotateCounterClockwise(bool revertRotations) {

    _recursive_ccw_rotation(this, revertRotations);
}

VirtualNode *VirtualNode::getNodeLeft() {
    return this->vnode_left ?: nullptr;
}

VirtualNode *VirtualNode::getNodeRight() {
    return this->vnode_right ?: nullptr;
}

void Utree::addMember(VirtualNode *iNode, bool isStartNode) {

    this->listVNodes.push_back(iNode);
    if (isStartNode) {
        this->startVNodes.push_back(iNode);
    }

}

std::vector<VirtualNode *> Utree::findPseudoRoot(VirtualNode *iNode) {

    // This method returns a vector which contains all the nodes in the path
    // from the start node to the pseudoroot
    std::vector<VirtualNode *> path2root;

    auto *CurrentNode = new VirtualNode;

    // Add the first element of the path to the root
    path2root.push_back(iNode);

    CurrentNode = iNode;

    // Add all the other nodes traversing the tree in post-order
    do {
        CurrentNode = CurrentNode->getNodeUp();
        path2root.push_back(CurrentNode);

    } while (CurrentNode != CurrentNode->getNodeUp()->getNodeUp());

    if (this->fixPseudoRootOnNextSubtree) {

        path2root.push_back(CurrentNode->getNodeUp());

    }

    return path2root;
}

std::string Utree::printTreeNewick(bool showInternalNodeNames) {
    std::string s, terminator;

    this->_updateStartNodes();

    s += "(";
    for (int i = 0; i < this->startVNodes.size(); i++) {
        if (i == this->startVNodes.size() - 1) {
            terminator = ");";
        } else {
            terminator = ",";
        }
        s += _recursiveFormatNewick(this->startVNodes.at(i), showInternalNodeNames) + terminator;
    }

    return s;
}

std::string Utree::_recursiveFormatNewick(VirtualNode *n, bool showInternalNodeNames) {

    std::stringstream newick;

    if (n->isTerminalNode()) {
        //std::stringstream newick;
        newick << n->vnode_name << ":" << n->vnode_branchlength;
        //return newick.str();
    } else {
        //std::stringstream newick;
        newick << "(";
        newick << _recursiveFormatNewick(n->getNodeLeft(), showInternalNodeNames);
        newick << ",";
        newick << _recursiveFormatNewick(n->getNodeRight(), showInternalNodeNames);
        newick << ")";
        if (showInternalNodeNames) {
            newick << n->vnode_name;
        }
        newick << ":" << n->vnode_branchlength;
        //return newick.str();
    }
    return newick.str();

}

Utree::~Utree() = default;

Utree::Utree() {

    this->fixPseudoRootOnNextSubtree = true;

}

void Utree::_updateStartNodes() {

    bool currValue = this->fixPseudoRootOnNextSubtree;
    this->fixPseudoRootOnNextSubtree = true;
    std::vector<VirtualNode *> sideA = this->findPseudoRoot(this->listVNodes.at(0));
    VirtualNode *sideA_node = sideA.back();

    this->fixPseudoRootOnNextSubtree = false;
    std::vector<VirtualNode *> sideB = this->findPseudoRoot(this->listVNodes.at(0));
    VirtualNode *sideB_node = sideB.back();

    this->fixPseudoRootOnNextSubtree = currValue;

    this->startVNodes.clear();

    this->startVNodes.push_back(sideA_node);
    this->startVNodes.push_back(sideB_node);

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

int VirtualNode::indexOf() {

    VirtualNode *parent = this->getNodeUp();
    if (parent->getNodeLeft() == this) {
        return 0;
    } else if (parent->getNodeRight() == this) {
        return 1;
    } else {
        //perror("The pointer index of the current node is unknown");
        return 2;
    }
}

bool VirtualNode::swapNode(VirtualNode *targetNode, MoveDirections move_direction) {

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


/*        std::cout << pnode->printNeighbours() << std::endl;
        std::cout <<qnode->printNeighbours() << std::endl;
        std::cout <<pnode_parent->printNeighbours() << std::endl;
        std::cout <<qnode_parent->printNeighbours() << std::endl;*/

        bool invertQParents = false;
        bool invertPParents = false;

        switch (move_direction) {

            case MoveDirections::left:
                pnode->RotateCounterClockwise();
                pnode_parent = pnode->getNodeUp();

                pnode->disconnectNode();
                qnode->disconnectNode();

                qnode_parent->connectNode(pnode);
                pnode_parent->connectNode(qnode);

                break;

            case MoveDirections::right:

                pnode->RotateClockwise();
                pnode_parent = pnode->getNodeUp();

                pnode->disconnectNode();
                qnode->disconnectNode();

                qnode_parent->connectNode(pnode);
                pnode_parent->connectNode(qnode);


                break;

            case MoveDirections::up:


                // Prune nodes
                pnode->disconnectNode();
                qnode->disconnectNode();

                // Save value of rotation
                if (qnode->vnode_rotated != 0) {
                    invertQParents = true;
                }

                if (pnode->vnode_rotated != 0) {
                    invertPParents = true;
                }

                // Reset rotation
                pnode->ResetNodeDirections();
                qnode->ResetNodeDirections();

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

                //pnode_parent->connectNode(qnode);
                //qnode_parent->connectNode(pnode);

                break;

            default:

                std::cout << "I cannot move the node" << std::endl;
                exit(EXIT_FAILURE);
        }
/*
        // If both nodes are in the pathway to the root
        if (pnode->isParent(qnode)) {

            pnode->RotateClockwise();
            pnode_parent = pnode->getNodeUp();

            pnode->disconnectNode();
            qnode->disconnectNode();

            qnode_parent->connectNode(pnode);
            pnode_parent->connectNode(qnode);

        } else if (qnode->isParent(this)) {

            qnode->RotateClockwise();
            //qnode->RotateCounterClockwise();
            qnode_parent = qnode->getNodeUp();

            qnode->disconnectNode();
            pnode->disconnectNode();

            pnode_parent->connectNode(qnode);
            qnode_parent->connectNode(pnode);

            // If the nodes are not in the pathway to the root
        } else if (!pnode->isParent(qnode) and !qnode->isParent(this)) {
            bool invertQParents = false;
            bool invertPParents = false;

            // Prune nodes
            pnode->disconnectNode();
            qnode->disconnectNode();

            // Save value of rotation
            if (qnode->vnode_rotated != 0) {
                invertQParents = true;
            }

            if (pnode->vnode_rotated != 0) {
                invertPParents = true;
            }

            // Reset rotation
            pnode->ResetNodeDirections();
            qnode->ResetNodeDirections();

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


        }*/

        execstatus = true;
    }

    return execstatus;
}

void VirtualNode::disconnectNode() {

    // 1. Get the location where this node is connected on the parent node and disconnect
    switch (this->indexOf()) {
        case 0:
            // The node is connected on the left side
            this->getNodeUp()->setNodeLeft(nullptr);
            break;
        case 1:
            // The node is connected on the right side
            this->getNodeUp()->setNodeRight(nullptr);
            break;
        case 2:
            // The node is connect on the up side (pseudoroot)
            this->getNodeUp()->setNodeUp(nullptr);
            break;
        default:
            break;
    }

    // 2. Disconnect the parent node
    this->setNodeUp(nullptr);
}

void VirtualNode::_oneway_connectNode(VirtualNode *inNode) {

    // Check which direction is still available (either left or right)
    if (!this->getNodeLeft()) {
        this->setNodeLeft(inNode);
        return;
    } else if (!this->getNodeRight()) {
        this->setNodeRight(inNode);
        return;
    } else if (!this->getNodeUp()) {
        this->setNodeUp(inNode);
        return;
    } else {
        perror("No direction available in the VirtualNode to add a new child");
    }

}

bool VirtualNode::isParent(VirtualNode *inNode) {

    bool status = false;

    VirtualNode *CurrentNode = inNode;

    do {
        CurrentNode = CurrentNode->getNodeUp();
        if (CurrentNode == this) {
            status = true;
            break;
        };
    } while (CurrentNode != CurrentNode->getNodeUp()->getNodeUp());

    return status;
}

void VirtualNode::ResetNodeDirections() {

    if (this->vnode_rotated != 0) {
        switch (this->vnode_rotated) {
            case 1:
                this->RotateCounterClockwise(true);
                break;
            case 2:
                this->RotateClockwise(true);
                break;
            default:
                break;
        }
    }

}

void VirtualNode::_recursive_cw_rotation(VirtualNode *n, bool revertRotations) {

    bool continueTraverse = true;
    if (n->getNodeUp() != nullptr) {
        if (n->getNodeUp()->getNodeUp() == n) {
            continueTraverse = false;
        }
    }

    if (revertRotations and n->vnode_rotated == 0) { return; }

    // rotate the pointers only if the node is not a leaf
    if (!n->isTerminalNode()) {
        auto curr_vn_up = new VirtualNode;
        auto curr_vn_left = new VirtualNode;
        auto curr_vn_right = new VirtualNode;
        auto n_left = new VirtualNode;
        auto n_right = new VirtualNode;

        curr_vn_up = n->getNodeUp();
        curr_vn_left = n->getNodeLeft();
        curr_vn_right = n->getNodeRight();

        //n->setNodeUp(curr_vn_left);
        //n->setNodeLeft(curr_vn_right);
        //n->setNodeRight(curr_vn_up);

        n->setNodeUp(curr_vn_right);
        n->setNodeLeft(curr_vn_up);
        n->setNodeRight(curr_vn_left);

        // Store the information about the rotation
        if (n->vnode_rotated != 0) {
            n->vnode_rotated = 0;
        } else {
            n->vnode_rotated = 1;
        }

        if (revertRotations) {
            n_left = curr_vn_left;
            n_right = curr_vn_right;
        } else {
            n_left = n->getNodeLeft();
            n_right = n->getNodeRight();
        }

        if (continueTraverse) {

            // Move to the newly set left side
            if (n->getNodeLeft() != nullptr) {
                _recursive_cw_rotation(n->getNodeLeft(), revertRotations);
            }
            // Move to the newly set right side
            if (n->getNodeRight() != nullptr) {
                _recursive_cw_rotation(n->getNodeRight(), revertRotations);
            }

        } else {

            //std::cout << "This node is the pseudoroot" << std::endl;

        }
    }

}

void VirtualNode::_recursive_ccw_rotation(VirtualNode *n, bool revertRotations) {
    bool continueTraverse = true;
    if (n->getNodeUp() != nullptr) {
        if (n->getNodeUp()->getNodeUp() == n) {
            continueTraverse = false;
        }
    }

    if (revertRotations and n->vnode_rotated == 0) { return; }

    // rotate the pointers only if the node is not a leaf
    if (!n->isTerminalNode()) {

        auto curr_vn_up = new VirtualNode;
        auto curr_vn_left = new VirtualNode;
        auto curr_vn_right = new VirtualNode;
        auto n_left = new VirtualNode;
        auto n_right = new VirtualNode;

        curr_vn_up = n->getNodeUp();
        curr_vn_left = n->getNodeLeft();
        curr_vn_right = n->getNodeRight();

        // Revert pointers
        //n->setNodeUp(curr_vn_right);
        //n->setNodeLeft(curr_vn_up);
        //n->setNodeRight(curr_vn_left);
        n->setNodeUp(curr_vn_left);
        n->setNodeLeft(curr_vn_right);
        n->setNodeRight(curr_vn_up);

        // Store the information about the rotation
        if (n->vnode_rotated != 0) {
            n->vnode_rotated = 0;
        } else {
            n->vnode_rotated = 2;
        }

        if (revertRotations) {
            n_left = curr_vn_left;
            n_right = curr_vn_right;
        } else {
            n_left = n->getNodeLeft();
            n_right = n->getNodeRight();
        }

        if (continueTraverse) {

            // Move to the newly set left side
            if (n_left != nullptr) {
                _recursive_ccw_rotation(n_left, revertRotations);
            }
            // Move to the newly set right side
            if (n_right != nullptr) {
                _recursive_ccw_rotation(n_right, revertRotations);
            }


        } else {

            //std::cout << "This node is the pseudoroot" << std::endl;

        }

    }

}

VirtualNode::~VirtualNode() = default;
