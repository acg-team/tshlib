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
#include <iomanip>
#include <glog/logging.h>
#include <map>

#include "Utree.hpp"
#include "Alignment.hpp"

/*
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
            ichild->_setNodeLeft(nullptr);
            ichild->_setNodeRight(nullptr);

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

    out_tree->startVNodes.at(0)->_setNodeUp(out_tree->startVNodes.at(1));
    out_tree->startVNodes.at(1)->_setNodeUp(out_tree->startVNodes.at(0));

}
*/

using namespace tshlib;
VirtualNode *UtreeUtils::getPseudoRoot(VirtualNode *vn){

    VirtualNode *vtemp=vn;
    while(!vtemp->isRootNode()){
        vtemp=vtemp->getNodeUp();
    }

    return vtemp;
}

namespace tshlib {

    VirtualNode::VirtualNode() {
        // Initialise a new VirtualNode completely disconnected from the tree
        this->vnode_right = nullptr;
        this->vnode_left = nullptr;
        this->vnode_up = nullptr;
        this->vnode_branchlength = 0;
        this->vnode_leaf = false;
        this->vnode_rotated = NodeRotation::undef;
        //this->vnode_character = NULL;
        //this->vnode_seqid = -1;

    };


    void VirtualNode::connectNode(VirtualNode *inVNode) {

        // Connect this node to the parent node
        inVNode->_setNodeUp(this);
        //inNode->_oneway_connectNode(this);
        this->_oneway_connectNode(inVNode);
    }


    bool VirtualNode::isTerminalNode() {

        return this->vnode_leaf;
    }


    bool VirtualNode::isRootNode() {

        return this->getNodeUp() == nullptr;
    }


    void VirtualNode::setNodeName(const std::string s) {

        this->vnode_name = s;
    }


    void VirtualNode::_setNodeRight(VirtualNode *inVNode) {

        this->vnode_right = inVNode;

    }


    void VirtualNode::_setNodeLeft(VirtualNode *inVNode) {

        this->vnode_left = inVNode;

    }


    void VirtualNode::_setNodeUp(VirtualNode *inVNode) {

        this->vnode_up = inVNode;

    }


    const std::string VirtualNode::getNodeName() {

        return this->vnode_name;
    }


    double VirtualNode::getIota() {

        return this->vnode_iota;
    }


    void VirtualNode::setIota(double iota) {

        this->vnode_iota = iota;
    }


    double VirtualNode::getBeta() {

        return this->vnode_beta;
    }


    void VirtualNode::setBeta(double beta) {

        this->vnode_beta = beta;
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


    void VirtualNode::_traverseVirtualNodeTree() {

        std::cout << "[Name] " << this->vnode_name << std::endl;

        if (this->vnode_up) {
            std::cout << "[Up name] " << this->getNodeUp()->vnode_name << std::endl;
        } else {
            std::cout << "[Up name] " << " ' ' " << std::endl;
        }

        std::cout << "[Node branch length] " << this->vnode_branchlength << std::endl;
        std::cout << "[Node iota] " << this->vnode_iota << std::endl;
        std::cout << "[Node beta] " << this->vnode_beta << std::endl;
        //std::cout << "[Node Pr] " << this->vnode_Pr.rows() << " x " << this->vnode_Pr.cols() << std::endl;

        //if(this->vnode_character!=NULL){
        //    std::cout<<"[Node char] "<<this->vnode_character<<std::endl;
        //} else{
        //    std::cout<<"[Node char] '' "<<std::endl;
        //}

        //std::cout<<"[Node setA_flag] "<<this->vnode_setA_backup<<std::endl;

        if (this->isTerminalNode()) {
            std::cout << std::endl;
        } else {

            std::cout << "[Left name] " << this->getNodeLeft()->vnode_name << std::endl;
            std::cout << "[Right name] " << this->getNodeRight()->vnode_name << std::endl;
            std::cout << std::endl;

            this->getNodeLeft()->_traverseVirtualNodeTree();
            this->getNodeRight()->_traverseVirtualNodeTree();
        }


    }


    double VirtualNode::computeTotalTreeLength() {

        if (this->isTerminalNode()) {
            return this->vnode_branchlength;
        } else {
            return this->vnode_branchlength +
                   this->getNodeLeft()->computeTotalTreeLength() +
                   this->getNodeRight()->computeTotalTreeLength();
        }
    }

/*
    void VirtualNode::recombineFv() {

        // For each internal node
        if (!this->isTerminalNode()) {

            // For each column of the alignment
            for (unsigned int k = 0; k < this->vnode_Fv_operative.size(); k++) {

                Eigen::VectorXd &fvL = this->getNodeLeft()->vnode_Fv_operative.at(k);
                Eigen::VectorXd &fvR = this->getNodeRight()->vnode_Fv_operative.at(k);

                this->vnode_Fv_operative.at(k) = this->getPr() * (fvL).cwiseProduct(fvR);

            }

            // Before moving to the next node, recompute the empty column fv quantities
            Eigen::VectorXd &fvE_L = this->getNodeLeft()->vnode_Fv_empty_operative;
            Eigen::VectorXd &fvE_R = this->getNodeRight()->vnode_Fv_empty_operative;

            this->vnode_Fv_empty_operative = this->getPr() * (fvE_L).cwiseProduct(fvE_R);


        }

    }


    void VirtualNode::revertFv() {

        this->vnode_Fv_operative.clear();
    }


    void VirtualNode::keepFv() {

        this->vnode_Fv_best.clear();
        this->vnode_Fv_best = this->vnode_Fv_operative;
        this->vnode_Fv_empty_best = this->vnode_Fv_empty_operative;

    }
*/

    void VirtualNode::clearChildren() {
        this->vnode_left = nullptr;
        this->vnode_right = nullptr;
    }

/*
    void VirtualNode::printAncestralFlagOnFile(FILE *fid) {

        //fprintf(fid,"%s %d\n",this->vnode_name.c_str(),this->vnode_setA_backup);

        if (this->isTerminalNode()) {

        } else {
            this->getNodeLeft()->printAncestralFlagOnFile(fid);
            this->getNodeRight()->printAncestralFlagOnFile(fid);
        }

    }
*/

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

        VirtualNode *CurrentNode;

        // Add the first element of the path to the root
        //path2root.push_back(inVNode);

        CurrentNode = inVNode;

        // Add all the other nodes traversing the tree in post-order
        while(CurrentNode!= nullptr){

            path2root.push_back(CurrentNode);

            if(CurrentNode == CurrentNode->getNodeUp()->getNodeUp()){break;}


            if(CurrentNode->getNodeUp() != nullptr){
                CurrentNode = CurrentNode->getNodeUp();
            }else{
                CurrentNode = nullptr;
            }

        }
        /*
        do {
            if(CurrentNode->getNodeUp() != nullptr){
                CurrentNode = CurrentNode->getNodeUp();
                path2root.push_back(CurrentNode);
            }else{
                break;
            }

        } while (CurrentNode->getNodeUp() != nullptr && CurrentNode != CurrentNode->getNodeUp()->getNodeUp());
        */
        if (fixPseudoRootOnNextSubtree) {

            if(CurrentNode->getNodeUp() != nullptr) {
                path2root.push_back(CurrentNode->getNodeUp());
            }

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


    Utree::~Utree() {

        for (std::vector<VirtualNode *>::reverse_iterator i = this->listVNodes.rbegin(); i < this->listVNodes.rend(); i++) {
            VirtualNode *vnode = *i;
            delete vnode;
        }

        delete this->rootnode;

    }


    double Utree::computeTotalTreeLength() {
        double t_length;

        t_length = 0.0;
        for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
            t_length += this->listVNodes.at(i)->vnode_branchlength;
        }

        return t_length;
    }


    void Utree::setIota(double tau, double mu) {
        VirtualNode *vn;
        double T;

        if (fabs(mu) < 1e-8) {
            perror("ERROR in set_iota: mu too small");
        }

        T = tau + 1 / mu;

        if (fabs(T) < 1e-8) {
            perror("ERROR in set_iota: T too small");
        }

        for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
            vn = this->listVNodes.at(i);

            //if (vn->isRootNode()) {
            //    vn->vnode_iota=(1/mu)/T;
            //} else {
            vn->vnode_iota = vn->vnode_branchlength / T;
            //}
        }

    }


    void Utree::setBeta(double tau, double mu) {
        VirtualNode *vn;

        if (fabs(mu) < 1e-8) {
            perror("ERROR : mu too small");
        }

        for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
            vn = this->listVNodes.at(i);

            //if (vn->isRootNode()) {
            //    vn->vnode_beta=1.0;
            //} else {
            if (fabs(vn->vnode_branchlength) < 1e-8) {
                perror("ERROR : branch_length too small");
            }
            vn->vnode_beta = (1 - exp(-mu * vn->vnode_branchlength)) / (mu * vn->vnode_branchlength);
            //}

        }

    }


    void Utree::_printUtree() {

        VirtualNode *vn;

        for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
            vn = this->listVNodes.at(i);

            std::cout << "[Node name] " << vn->vnode_name << std::endl;

            std::cout << "[Up child node name] " << vn->getNodeUp()->vnode_name << std::endl;

            if (!vn->isTerminalNode()) {
                std::cout << "[L. child node name] " << vn->getNodeLeft()->vnode_name << std::endl;
                std::cout << "[R. child node name] " << vn->getNodeRight()->vnode_name << std::endl;
            }
            std::cout << "[Node branch length] " << vn->vnode_branchlength << std::endl;
            std::cout << "[Node iota] " << vn->vnode_iota << std::endl;
            std::cout << "[Node beta] " << vn->vnode_beta << std::endl;
            //std::cout << "[Node Pr] " << vn->vnode_Pr.rows() << " x " << vn->vnode_Pr.cols() << std::endl;

//        if(vn->vnode_character!=NULL){
//            std::cout<<"[Node char] "<<vn->vnode_character<<std::endl;
//        } else{
//            std::cout<<"[Node char] '' "<<std::endl;
//        }

            //std::cout<<"[Node setA_flag] "<<vn->vnode_setA_backup<<std::endl;

            std::cout << std::endl;
        }

    }

/*
    void Utree::clearFv() {

        for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
            this->listVNodes.at(i)->vnode_Fv_backup.clear();
        }
    }

*/
    Utree::Utree() {

        // Added virtual root to utree
        auto root = new VirtualNode;
        root->vnode_name = "root";

        this->rootnode = root;

    };


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

            VLOG(1) << strpath << std::endl;
            vnodes2root.clear();
            strpath.clear();

        }

    }


    void Utree::saveTreeOnFile(std::string outfilepath) {

        std::ofstream outfile;

        outfile.open(outfilepath, std::ios_base::app);
        outfile << this->printTreeNewick(true) << std::endl;

    }


    int Utree::getMaxNodeDistance() {

        long int max_distance = 0;

        for (auto &vnode:this->listVNodes) {

            if (vnode->isTerminalNode()) {

                if (max_distance < this->findPseudoRoot(vnode).size()) {
                    max_distance = this->findPseudoRoot(vnode).size();
                };

            }
        }
        return max_distance * 2;
    }


    std::string VirtualNode::printNeighbours() {

        std::stringstream description;

        if (!this->isTerminalNode()) {

            if (this->getNodeUp()) {
                description << this->vnode_name << " (^" << this->getNodeUp()->vnode_name << ";";
            } else {
                description << this->vnode_name << " (^NULL;";
            }

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
                this->getNodeUp()->_setNodeLeft(nullptr);
                break;
            case NodePosition::right:
                // The node is connected on the right side
                this->getNodeUp()->_setNodeRight(nullptr);
                break;
            case NodePosition::up:
                // The node is connect on the up side (pseudoroot)
                this->getNodeUp()->_setNodeUp(nullptr);
                break;
            default:
                break;
        }

        // 2. Disconnect the parent node
        this->_setNodeUp(nullptr);
    }


    void VirtualNode::_oneway_connectNode(VirtualNode *inVNode) {

        // Check which direction is still available (either left or right)
        if (!this->getNodeLeft()) {
            this->_setNodeLeft(inVNode);
            return;
        } else if (!this->getNodeRight()) {
            this->_setNodeRight(inVNode);
            return;
        } else if (!this->getNodeUp()) {
            this->_setNodeUp(inVNode);
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

/*
    void VirtualNode::_recursiveSetDescCount(Alignment &MSA, bool isReference, int colnum) {

        if (this->isTerminalNode()) {
            //this->vnode_descCount_backup = (this->vnode_character == '-' ? 0 : 1);
            if (isReference) {
                this->vnode_descCount_backup.at(colnum) = (MSA.align_dataset.at(this->vnode_seqid)->seq_data.at(colnum) == '-' ? 0 : 1);
            } else {
                this->vnode_descCount_operative.at(colnum) = (MSA.align_dataset.at(this->vnode_seqid)->seq_data.at(colnum) == '-' ? 0 : 1);
            }
        } else {
            this->getNodeLeft()->_recursiveSetDescCount(MSA, isReference, colnum);
            this->getNodeRight()->_recursiveSetDescCount(MSA, isReference, colnum);
            if (isReference) {
                this->vnode_descCount_backup.at(colnum) = this->getNodeLeft()->vnode_descCount_backup.at(colnum) +
                                                          this->getNodeRight()->vnode_descCount_backup.at(colnum);
            } else {
                this->vnode_descCount_operative.at(colnum) = this->getNodeLeft()->vnode_descCount_operative.at(colnum) +
                                                             this->getNodeRight()->vnode_descCount_operative.at(colnum);
            }
        }

    }


    void VirtualNode::_recursiveSetAncestralFlag(Alignment &MSA, int colnum, bool isReference) {
        int descCount;

        if (this->isTerminalNode()) {
            if (isReference) {
                descCount = this->vnode_descCount_backup.at(colnum);
                //this->setSetA(descCount == num_gaps);
                this->vnode_setA_backup.at(colnum) = (descCount == MSA.align_num_characters.at(colnum));
            } else {
                descCount = this->vnode_descCount_operative.at(colnum);
                //this->setSetA(descCount == num_gaps);
                this->vnode_setA_operative.at(colnum) = (descCount == MSA.align_num_characters.at(colnum));
            }
        } else {
            this->getNodeLeft()->_recursiveSetAncestralFlag(MSA, colnum, isReference);
            this->getNodeRight()->_recursiveSetAncestralFlag(MSA, colnum, isReference);
            if (isReference) {
                descCount = this->vnode_descCount_backup.at(colnum);
                this->vnode_setA_backup.at(colnum) = (descCount == MSA.align_num_characters.at(colnum));
            } else {
                descCount = this->vnode_descCount_operative.at(colnum);
                this->vnode_setA_operative.at(colnum) = (descCount == MSA.align_num_characters.at(colnum));
            }
            //descCount = this->vnode_descCount_backup;
            //this->setSetA(descCount == num_gaps);
        }

    }


    void VirtualNode::setAncestralFlag(Alignment &MSA, int colnum, bool isReference) {

        this->_recursiveSetDescCount(MSA, isReference, colnum);

        this->_recursiveSetAncestralFlag(MSA, isReference, false);

    }


    void Utree::setLeafState(std::string s) {

        for (auto &node:this->listVNodes) {
            if (node->isTerminalNode()) {
                node->setLeafCharacter(s.at(node->vnode_seqid));
            }
        }

    }

    void Utree::prepareSetADesCountOnNodes(int numcol, int lengthAlphabet) {

        for (auto &node:this->listVNodes) {

            node->initialiseLikelihoodComponents(numcol, lengthAlphabet);

        }

    }
*/
    void Utree::printAllNodesNeighbors() {

        for (auto &node:this->listVNodes) {
            VLOG(2) << node->printNeighbours();
        }

    }

    void Utree::addVirtualRootNode() {


        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {
            // Update starting nodes in the startVNodes array
            this->_updateStartNodes();

            // Clear children nodes in the virtual root
            this->rootnode->clearChildren();

            // Disconnect the bidirectional connection in the pseudoroot nodes
            this->startVNodes.at(0)->disconnectNode();

            // Reconnect nodes to virtual root node
            this->rootnode->connectNode(this->startVNodes.at(0));
            this->rootnode->connectNode(this->startVNodes.at(1));
        }
    }

    void Utree::removeVirtualRootNode() {

        if (rootnode->getNodeRight() != nullptr && rootnode->getNodeLeft() != nullptr) {
            // Get reference children nodes from virtual root
            VirtualNode *left = this->rootnode->getNodeLeft();
            VirtualNode *right = this->rootnode->getNodeRight();

            // Disconnect nodes from virtual root node
            left->disconnectNode();
            right->disconnectNode();

            // Establish bidirectional connection
            left->_setNodeUp(right);
            right->_setNodeUp(left);

            // Clear residual connection with virtual root node
            this->rootnode->clearChildren();

        }
    }

    std::vector<VirtualNode *> Utree::computePathBetweenNodes(VirtualNode *vnode_1, VirtualNode *vnode_2) {

        std::vector<VirtualNode *> list_vnode_to_root;

        std::vector<VirtualNode *> path2root_1 = this->findPseudoRoot(vnode_1, false);
        std::vector<VirtualNode *> path2root_2 = this->findPseudoRoot(vnode_2, false);

        std::reverse(path2root_1.begin(), path2root_1.end());
        std::reverse(path2root_2.begin(), path2root_2.end());

        // Test map

        std::map<VirtualNode *,int> tmpMap;

        for(int i=0; i<path2root_1.size(); i++){
            tmpMap.insert(std::pair<VirtualNode *,int>(path2root_1.at(i), i));
        }

        for(int i=0; i<path2root_2.size(); i++) {

            if (tmpMap.find(path2root_2.at(i)) == tmpMap.end()){

                tmpMap.insert(std::pair<VirtualNode *,int>(path2root_2.at(i), tmpMap.size()));

            }
        }

        auto newvector = std::vector<VirtualNode *>(tmpMap.size());
        // Create a map iterator and point to beginning of map
        //std::map<VirtualNode *,int>::iterator it = tmpMap.begin();

        // Iterate over the map using c++11 range based for loop
        for (std::pair<VirtualNode *,int> element : tmpMap) {

            // Accessing KEY from element
            VirtualNode *node = element.first;

            // Accessing VALUE from element.
            int order = element.second;

            newvector[order] = node;

        }

        std::reverse(newvector.begin(), newvector.end());


        //path2root_1.insert( path2root_1.end(), path2root_2.begin(), path2root_2.end() );
        //auto last = std::unique(path2root_1.begin(), path2root_1.end());
        // v now holds {1 2 3 4 5 6 7 x x x x x x}, where 'x' is indeterminate
        //path2root_1.erase(last, path2root_1.end());
        //list_vnode_to_root = path2root_1;

        //list_vnode_to_root = this->_unique(path2root_1, path2root_2);


        //delete(path2root_1);
        //delete(path2root_2);
        return newvector;
    }

    std::vector<VirtualNode *> Utree::_unique(std::vector<VirtualNode *> &list_nodes_n1, std::vector<VirtualNode *> &list_nodes_n2) {

        std::vector<VirtualNode *> list_nodes;
        VirtualNode *n1;
        VirtualNode *n2;

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
            n1 = list_nodes_n1.at(list_nodes_n1.size() - 1);
            list_nodes.push_back(n1);
            list_nodes_n1.pop_back();
        }

        while (list_nodes_n2.size() > 0) {
            n2 = list_nodes_n2.at(list_nodes_n2.size() - 1);
            list_nodes.push_back(n2);
            list_nodes_n2.pop_back();
        }

        return list_nodes;
    }

    std::vector<VirtualNode *> Utree::getPostOrderNodeList() {

        bool removeRoot = false;
        std::vector<VirtualNode *> rlist;
        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {
            addVirtualRootNode();
            removeRoot = true;
        }
        _getPostOrderNodeList(rlist, rootnode);


        if(removeRoot) {
            removeVirtualRootNode();
        }

        return rlist;
    }

    std::vector<VirtualNode *> Utree::getPostOrderNodeList(VirtualNode *startNode) {

        bool removeRoot = false;
        std::vector<VirtualNode *> rlist;

        if (rootnode->getNodeRight() == nullptr && rootnode->getNodeLeft() == nullptr) {
            addVirtualRootNode();
            removeRoot = true;
        }
        _getPostOrderNodeList(rlist, startNode);

        if(removeRoot) {
            removeVirtualRootNode();
        }
        return rlist;
    }


    void Utree::_getPostOrderNodeList(std::vector<VirtualNode *> &rlist, VirtualNode *node ){


        if (!node->isTerminalNode()) {

            _getPostOrderNodeList(rlist, node->getNodeLeft());
            _getPostOrderNodeList(rlist, node->getNodeRight());

        }

        rlist.push_back(node);

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
            VirtualNode *curr_vn_up;
            VirtualNode *curr_vn_left;
            VirtualNode *curr_vn_right;
            VirtualNode *n_up;


            curr_vn_up = vnode->getNodeUp();
            curr_vn_left = vnode->getNodeLeft();
            curr_vn_right = vnode->getNodeRight();

            vnode->_setNodeUp(curr_vn_right);
            vnode->_setNodeLeft(curr_vn_up);
            vnode->_setNodeRight(curr_vn_left);

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

            VirtualNode *curr_vn_up;
            VirtualNode *curr_vn_left;
            VirtualNode *curr_vn_right;

            VirtualNode *n_up;

            curr_vn_up = vnode->getNodeUp();
            curr_vn_left = vnode->getNodeLeft();
            curr_vn_right = vnode->getNodeRight();

            // Revert pointers
            vnode->_setNodeUp(curr_vn_left);
            vnode->_setNodeLeft(curr_vn_right);
            vnode->_setNodeRight(curr_vn_up);

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
/*
    void VirtualNode::initialiseLikelihoodComponents(int numcol, int lengthAlphabet) {

        this->vnode_descCount_backup.resize(numcol);
        this->vnode_descCount_operative.resize(numcol);
        this->vnode_descCount_best.resize(numcol);

        this->vnode_setA_backup.resize(numcol);
        this->vnode_setA_operative.resize(numcol);
        this->vnode_setA_best.resize(numcol);

        //this->vnode_Fv_empty_backup.resize(lengthAlphabet);
        //this->vnode_Fv_empty_operative.resize(lengthAlphabet);
        //this->vnode_Fv_empty_best.resize(lengthAlphabet);
        this->vnode_Fv_terminal.resize(numcol);
        this->vnode_Fv_backup.resize(numcol);
        this->vnode_Fv_operative.resize(numcol);
        this->vnode_Fv_partial_operative.resize(numcol);
        this->vnode_Fv_best.resize(numcol);

    }
*/
    bool VirtualNode::isPseudoRootNode() {
        return this->getNodeUp()->getNodeUp() == this;
    }
/*
    void VirtualNode::_printFV() {

        __print2Dmat(this, this->vnode_Fv_operative);

        std::vector<Eigen::VectorXd> temp;
        temp.push_back(this->vnode_Fv_empty_operative);
        __print2Dmat(this,temp );



    }

    void VirtualNode::__print2Dmat(VirtualNode *node, std::vector<Eigen::VectorXd> input){
        std::string line;
        std::ostringstream sout;

        double rows = input.at(0).size();
        double cols = input.size();

        // Initialization of the 2D vector
        std::vector<std::vector<double> > array2D;
        array2D.resize(rows);
        for (int i = 0; i < rows; ++i) {
            array2D[i].resize(cols);
        }

        for (int c = 0; c < input.size(); c++) {
            for (int a = 0; a < input.at(c).size(); a++) {
                array2D[a][c] = input.at(c)(a);
            }
        }
        VLOG(3) << "-----------------------------------------------------";
        VLOG(3) << "Node: " << node->vnode_name;


        for (int i = 0; i < rows; ++i) {
            line = "";
            for (int c = 0; c < cols; c++) {

                if (i == 0) {
                    sout << std::setfill('0') << std::setw(3) << c << "      | ";
                }

                line += std::to_string(array2D[i][c]);
                line += " | ";

            }
            if (i == 0) {
                VLOG(3) << sout.str();
            }
            VLOG(3) << line;
        }

    }

*/
    VirtualNode::~VirtualNode() = default;

}