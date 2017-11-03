//
// Created by Lorenzo Gatti on 26/10/17.
//

#include <string>
#include <random>

#include "Utree.hpp"


void ::UtreeUtils::_traverseTree(VirtualNode *target, PhyTree *source) {

    for (unsigned long i = 0; i < source->n_children(); i++) {
        auto ichild = new VirtualNode();
        if (!source->get_children().at(i)->isLeaf()) {

            ichild->vnode_id = source->get_children().at(i)->getNodeID();
            ichild->vnode_name = source->get_children().at(i)->getName();
            ichild->vnode_branchlength = source->get_children().at(i)->getBranchLength();

            target->addMember(ichild);
            _traverseTree(ichild, source->get_children().at(i));

        } else {

            ichild->vnode_id = source->get_children().at(i)->getNodeID();
            ichild->vnode_name = source->get_children().at(i)->getName();
            ichild->vnode_branchlength = source->get_children().at(i)->getBranchLength();
            // Set all the other directions to null
            ichild->setNodeLeft(nullptr);
            ichild->setNodeRight(nullptr);

            // Set the LEAF flag to true
            ichild->vnode_leaf = true;
            target->addMember(ichild);

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

            _traverseTree(ichild, in_tree->get_children().at(i));
        }

        out_tree->addMember(ichild);
    }


    // Collapse multiforcating trees to star tree pointing to the same pseudoroot
    // Pick root node at random within the node-vector

    std::default_random_engine generator;
    std::uniform_int_distribution<unsigned long> distribution(0, out_tree->topology.size() - 1);
    unsigned long randomIndex = distribution(generator);


    for (unsigned long i = 0; i < out_tree->topology.size(); i++) {
        if (i != randomIndex) {
            out_tree->topology.at(i)->setNodeUp(out_tree->topology.at(randomIndex));
        } else {
            if (randomIndex == 0) {
                randomIndex = randomIndex + 1;
            } else {
                randomIndex = randomIndex - 1;
            }
            out_tree->topology.at(i)->setNodeUp(out_tree->topology.at(randomIndex));
        }
    }

}


VirtualNode::VirtualNode() {
    // Initialise a new VirtualNode completely disconnected from the tree
    this->vnode_right = nullptr;
    this->vnode_left = nullptr;
    this->vnode_up = nullptr;
    this->vnode_branchlength = 0;
    this->vnode_leaf = false;

};

void VirtualNode::addMember(VirtualNode *inNode) {

    // Connect this node to the parent node
    inNode->setNodeUp(this);

    // Check which direction is still available (either left or right)
    if (!this->vnode_left) {
        this->setNodeLeft(inNode);
        return;
    } else if (!this->vnode_right) {
        this->setNodeRight(inNode);
        return;
    } else {
        perror("No direction available in the VirtualNode to add a new child");
    }


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

    auto *curr_vn_up = new VirtualNode;
    auto *curr_vn_left = new VirtualNode;
    auto *curr_vn_right = new VirtualNode;

    this->setNodeUp(curr_vn_left);
    this->setNodeLeft(curr_vn_right);
    this->setNodeRight(curr_vn_up);

}

void VirtualNode::RotateCounterClockwise() {

    auto *curr_vn_up = new VirtualNode;
    auto *curr_vn_left = new VirtualNode;
    auto *curr_vn_right = new VirtualNode;

    this->setNodeUp(curr_vn_right);
    this->setNodeLeft(curr_vn_up);
    this->setNodeRight(curr_vn_left);

}

VirtualNode *VirtualNode::getNodeLeft() {
    return this->vnode_left ?: nullptr;
}

VirtualNode *VirtualNode::getNodeRight() {
    return this->vnode_right ?: nullptr;
}


void Utree::addMember(VirtualNode *iNode) {

    this->topology.push_back(iNode);


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

    if (this->fixPseudoRootOnLeaf) {

        path2root.push_back(CurrentNode->getNodeUp());

    }

    return path2root;
}

std::string Utree::printTreeNewick() {
    std::string s, terminator;
    s += "(";
    for (int i = 0; i < this->topology.size(); i++) {
        if (i == this->topology.size() - 1) {
            terminator = ");";

        } else {
            terminator = ",";
        }
        s += _recursiveFormatNewick(this->topology.at(i)) + terminator;

    }

    return s;

}

std::string Utree::_recursiveFormatNewick(VirtualNode *n) {

    if (n->isTerminalNode()) {
        std::stringstream newick;
        newick << n->vnode_name << ":" << n->vnode_branchlength;
        return newick.str();
    } else {
        std::stringstream newick;
        newick << "(";
        newick << _recursiveFormatNewick(n->getNodeLeft());
        newick << ",";
        newick << _recursiveFormatNewick(n->getNodeRight());
        newick << ")";
        newick << ":" << n->vnode_branchlength;
        return newick.str();
    }


}

Utree::~Utree() {

}

Utree::Utree() {

    this->fixPseudoRootOnLeaf = true;

}






