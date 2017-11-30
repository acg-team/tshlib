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
#include "Alignment.hpp"

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

std::vector<VirtualNode *> UtreeUtils::fill_with_nodes(VirtualNode *n) {
    std::vector<VirtualNode *> list_nodes_n;

    VirtualNode *tmp;

    list_nodes_n.push_back(n);
    tmp=n;
    while (tmp->getNodeUp() != NULL) {
        tmp = tmp->getNodeUp();
        list_nodes_n.push_back(tmp);
    }

    return list_nodes_n;
}

std::vector<VirtualNode *> UtreeUtils::get_unique(std::vector<VirtualNode *> &list_nodes_n1, std::vector<VirtualNode *> &list_nodes_n2) {
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

    std::reverse(list_nodes.begin(), list_nodes.end());

    return list_nodes;
}

std::vector<VirtualNode *> UtreeUtils::get_path_from_nodes(VirtualNode *vn1, VirtualNode *vn2) {
    std::vector<VirtualNode *> list_nodes_n0;
    std::vector<VirtualNode *> list_nodes_n1;
    std::vector<VirtualNode *> list_nodes_n2;

    // add nodes from n1 to root
    list_nodes_n1 = fill_with_nodes(vn1);

    // add nodes from n2 to root
    list_nodes_n2 = fill_with_nodes(vn2);

    list_nodes_n0 = get_unique(list_nodes_n1, list_nodes_n2);

    return list_nodes_n0;
}

void UtreeUtils::recombineAllFv(std::vector<VirtualNode *> list_vnode_to_root){
    VirtualNode *vn;

    for(unsigned int k=0;k<list_vnode_to_root.size();k++){
        vn=list_vnode_to_root.at(k);
        vn->recombineFv();
    }

}

void UtreeUtils::revertAllFv(std::vector<VirtualNode *> list_vnode_to_root){
    VirtualNode *vn;

    for(unsigned int k=0;k<list_vnode_to_root.size();k++){
        vn=list_vnode_to_root.at(k);
        vn->revertFv();
    }

}

void UtreeUtils::keepAllFv(std::vector<VirtualNode *> list_vnode_to_root){
    VirtualNode *vn;

    for(unsigned int k=0;k<list_vnode_to_root.size();k++){
        vn=list_vnode_to_root.at(k);
        vn->keepFv();
    }

}

void ::UtreeUtils::associateNode2Alignment(Alignment *inMSA, Utree *inTree) {

    for(auto &node:inTree->listVNodes){

        if(node->isTerminalNode()){

            for(int i=0; i<inMSA->align_dataset.size();i++){

                if(inMSA->align_dataset.at(i)->seq_name.compare(node->vnode_name)==0){

                    node->vnode_seqid = i;
                    break;
                }

            }

        }

    }


}

VirtualNode *UtreeUtils::getPseudoRoot(VirtualNode *vn){

    VirtualNode *vtemp=vn;
    while(!vtemp->isRootNode()){
        vtemp=vtemp->getNodeUp();
    }

    return vtemp;
}

VirtualNode::VirtualNode() {
    // Initialise a new VirtualNode completely disconnected from the tree
    this->vnode_right = nullptr;
    this->vnode_left = nullptr;
    this->vnode_up = nullptr;
    this->vnode_branchlength = 0;
    this->vnode_leaf = false;
    this->vnode_rotated = NodeRotation::undef;
    this->vnode_character = NULL;
    this->vnode_seqid = -1;

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

bool VirtualNode::isRootNode() {

    return this->getNodeUp()->getNodeUp()==this;

}

void VirtualNode::setNodeName(const std::string s){
    this->vnode_name=s;
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

void VirtualNode:: setLeafCharacter(char ch) {

    //std::cout<<"I am node:"<<this->vnode_name<<" and I am setting: "<<ch<<std::endl;

    this->vnode_character = ch;
}

void VirtualNode::setSetA(bool b) {
    this->vnode_setA = b;
}

std::string VirtualNode::getNodeName(){
    return this->vnode_name;
}

void VirtualNode::setMSAFv(Eigen::VectorXd &fv) {
    this->vnode_Fv.push_back(fv);
}

bool VirtualNode::getSetA() {
    return this->vnode_setA;
}

double VirtualNode::getIota() {
    return this->vnode_iota;
}

void VirtualNode::setIota(double iota){
    this->vnode_iota=iota;
}


double VirtualNode::getBeta() {
    return this->vnode_beta;
}

void VirtualNode::setBeta(double beta){
    this->vnode_beta=beta;
}

const Eigen::MatrixXd & VirtualNode::getPr() {

    return this->vnode_Pr;

}

char VirtualNode::getLeafCharacter(){
    return this->vnode_character;
}

void VirtualNode::setChild(VirtualNode *vn){

    if(this->vnode_left!=NULL && this->vnode_right!=NULL){
        perror("ERROR, more than 2 children is not allowed");
    }

    if(this->vnode_left==NULL){
        this->vnode_left=vn;
        vn->setNodeUp(this);
    }else{
        this->vnode_right=vn;
        vn->setNodeUp(this);
    }
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

void VirtualNode::_traverseVirtualNodeTree(){

    std::cout<<"[Name] "<<this->vnode_name<<std::endl;

    if(this->vnode_up!=NULL){
        std::cout<<"[Up name] "<<this->getNodeUp()->vnode_name<<std::endl;
    }else{
        std::cout<<"[Up name] "<<" ' ' "<<std::endl;
    }

    std::cout<<"[Node branch length] "<<this->vnode_branchlength<<std::endl;
    std::cout<<"[Node iota] "<<this->vnode_iota<<std::endl;
    std::cout<<"[Node beta] "<<this->vnode_beta<<std::endl;
    std::cout<<"[Node Pr] "<<this->vnode_Pr.rows() << " x " << this->vnode_Pr.cols()<<std::endl;

    if(this->vnode_character!=NULL){
        std::cout<<"[Node char] "<<this->vnode_character<<std::endl;
    } else{
        std::cout<<"[Node char] '' "<<std::endl;
    }

    std::cout<<"[Node setA_flag] "<<this->vnode_setA<<std::endl;

    if(this->isTerminalNode()){
        std::cout<<std::endl;
    }else{

        std::cout<<"[Left name] "<<this->getNodeLeft()->vnode_name<<std::endl;
        std::cout<<"[Right name] "<<this->getNodeRight()->vnode_name<<std::endl;
        std::cout<<std::endl;

        this->getNodeLeft()->_traverseVirtualNodeTree();
        this->getNodeRight()->_traverseVirtualNodeTree();
    }


}

void VirtualNode::setNodeParent(VirtualNode *vn){
    //if(this->getNodeUp()!=NULL){
    //    perror("ERROR: root node already assigned");
    //}else{
        this->vnode_up=vn;
    //}
}

double VirtualNode::computeTotalTreeLength() {

    if (this->isTerminalNode()) {
        return this->vnode_branchlength;
    } else {
        return this->vnode_branchlength+
               this->getNodeLeft()->computeTotalTreeLength()+
               this->getNodeRight()->computeTotalTreeLength();
    }
}

void VirtualNode::setAllIotas(double tau, double mu){
    double T;

    if (fabs(mu) < 1e-8) {
        perror("ERROR in set_iota: mu too small");
    }

    T = tau + 1 / mu;

    if (fabs(T) < 1e-8) {
        perror("ERROR in set_iota: T too small");
    }

    if (this->vnode_up==NULL) {
            this->vnode_iota=(1/mu)/T;
    } else {
        this->getNodeLeft()->setAllIotas(tau,mu);
        this->getNodeRight()->setAllIotas(tau,mu);
        this->vnode_iota = this->vnode_branchlength / T;
    }

}

void VirtualNode::setAllBetas(double mu){

    if (fabs(mu) < 1e-8) {
        perror("ERROR in set_iota: mu too small");
    }

    if (this->vnode_up==NULL) {
        this->vnode_beta=1.0;
    } else {
        this->getNodeLeft()->setAllBetas(mu);
        this->getNodeRight()->setAllBetas(mu);
        if (fabs(this->vnode_branchlength) < 1e-8) {
            perror("ERROR : branch_length too small");
        }
        this->vnode_beta= (1.0 - exp(-mu * this->vnode_branchlength)) / (mu * this->vnode_branchlength);
    }

}

void VirtualNode::recombineFv(){
    unsigned int lL,lR;

    if(this->isTerminalNode()){
        return;
    }

    lL=this->getNodeLeft()->vnode_Fv.size();
    lR=this->getNodeRight()->vnode_Fv.size();

    if(lL!=lR){
        perror("ERROR: left fv size different from right fv size");
        exit(EXIT_FAILURE);
    }

    this->vnode_Fv_temp.clear();
    for(unsigned int k=0;k<lL;k++){
        Eigen::VectorXd &fvL=this->getNodeLeft()->vnode_Fv.at(k);
        Eigen::VectorXd &fvR=this->getNodeRight()->vnode_Fv.at(k);
        Eigen::VectorXd fv0=fvL.cwiseProduct(fvR);
        this->vnode_Fv_temp.push_back(fv0);
    }

}

void VirtualNode::revertFv(){

    this->vnode_Fv_temp.clear();

}

void VirtualNode::keepFv(){

    this->vnode_Fv_best.clear();
    this->vnode_Fv_best=this->vnode_Fv_temp;

}

void VirtualNode::clearChildren(){
    this->vnode_left= nullptr;
    this->vnode_right= nullptr;
}

void VirtualNode::printAncestralFlagOnFile(FILE *fid){

    fprintf(fid,"%s %d\n",this->vnode_name.c_str(),this->vnode_setA);

    if(this->isTerminalNode()){

    }else{
        this->getNodeLeft()->printAncestralFlagOnFile(fid);
        this->getNodeRight()->printAncestralFlagOnFile(fid);
    }

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


double Utree::computeTotalTreeLength() {
    double t_length;

    t_length=0.0;
    for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
        t_length+=this->listVNodes.at(i)->vnode_branchlength;
    }

    return t_length;
}

void Utree::setIota(double tau, double mu){
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
        vn=this->listVNodes.at(i);

        //if (vn->isRootNode()) {
        //    vn->vnode_iota=(1/mu)/T;
        //} else {
            vn->vnode_iota = vn->vnode_branchlength / T;
        //}
    }

}

void Utree::setBeta(double tau, double mu){
    VirtualNode *vn;

    if (fabs(mu) < 1e-8) {
        perror("ERROR : mu too small");
    }

    for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
        vn=this->listVNodes.at(i);

        //if (vn->isRootNode()) {
        //    vn->vnode_beta=1.0;
        //} else {
            if (fabs(vn->vnode_branchlength) < 1e-8) {
                perror("ERROR : branch_length too small");
            }
            vn->vnode_beta= (1 - exp(-mu * vn->vnode_branchlength)) / (mu * vn->vnode_branchlength);
        //}

    }

}

void Utree::_printUtree(){

    VirtualNode *vn;

    for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
        vn=this->listVNodes.at(i);

        std::cout<<"[Node name] "<<vn->vnode_name<<std::endl;

        std::cout << "[Up child node name] " << vn->getNodeUp()->vnode_name<< std::endl;

        if(!vn->isTerminalNode()) {
            std::cout << "[L. child node name] " << vn->getNodeLeft()->vnode_name<< std::endl;
            std::cout << "[R. child node name] " << vn->getNodeRight()->vnode_name<< std::endl;
        }
        std::cout<<"[Node branch length] "<<vn->vnode_branchlength<<std::endl;
        std::cout<<"[Node iota] "<<vn->vnode_iota<<std::endl;
        std::cout<<"[Node beta] "<<vn->vnode_beta<<std::endl;
        std::cout<<"[Node Pr] "<<vn->vnode_Pr.rows() << " x " << vn->vnode_Pr.cols()<<std::endl;

        if(vn->vnode_character!=NULL){
            std::cout<<"[Node char] "<<vn->vnode_character<<std::endl;
        } else{
            std::cout<<"[Node char] '' "<<std::endl;
        }

        std::cout<<"[Node setA_flag] "<<vn->vnode_setA<<std::endl;

        std::cout<<std::endl;
    }

}

void Utree::setPr(int extended_alphabet_size){

    VirtualNode *vn;

    for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
        vn=this->listVNodes.at(i);

        if (vn->vnode_Pr.rows() * vn->vnode_Pr.cols() != 0) {
            vn->vnode_Pr.resize(0, 0);
        }

        //if (!vn->isRootNode()) {
            vn->vnode_Pr.resize(extended_alphabet_size,extended_alphabet_size);

            for (int i = 0; i < extended_alphabet_size; i++) {
                for (int j = 0; j < extended_alphabet_size; j++) {
                    vn->vnode_Pr(i,j) = 0.1 * vn->vnode_branchlength; //fake
                }
            }

        //}
    }

}

void Utree::clearFv() {

    for (unsigned long i = 0; i < this->listVNodes.size(); i++) {
        this->listVNodes.at(i)->vnode_Fv.clear();
    }
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

void VirtualNode::_recursiveSetDescCount() {

    if (this->isTerminalNode()) {
        this->vnode_descCount = (this->vnode_character == '-' ? 0 : 1);
    } else {
        this->getNodeLeft()->_recursiveSetDescCount();
        this->getNodeRight()->_recursiveSetDescCount();
        this->vnode_descCount = this->getNodeLeft()->vnode_descCount +
                                this->getNodeRight()->vnode_descCount;
    }

}

void VirtualNode::_recursiveSetAncestralFlag(std::string &MSA_col, int num_gaps) {
    int descCount;

    if (this->isTerminalNode()) {
        descCount = this->vnode_descCount;
        this->setSetA(descCount == num_gaps);
    } else {
        this->getNodeLeft()->_recursiveSetAncestralFlag(MSA_col,num_gaps);
        this->getNodeRight()->_recursiveSetAncestralFlag(MSA_col,num_gaps);

        descCount = this->vnode_descCount;
        this->setSetA(descCount == num_gaps);
    }

}

void VirtualNode::setAncestralFlag(std::string MSA_col){

    int num_gaps;

    num_gaps=AlignUtils::countNumGapsInMSAColumn(MSA_col);

    //this->_updateStartNodes();

    //for (unsigned long i = 0; i < this->startVNodes.size(); i++) {
    this->_recursiveSetDescCount();
    //_recursiveSetDescCount(this->startVNodes.at(1));
    //}
    //for (unsigned long i = 0; i < this->startVNodes.size(); i++) {
    this->_recursiveSetAncestralFlag(MSA_col,num_gaps);
    //_recursiveSetAncestralFlag(this->startVNodes.at(1),MSA_col,num_gaps);
    //}

}

void Utree::setLeafState(std::string s){

    for (auto &node:this->listVNodes) {
        if(node->isTerminalNode()) {
            node->setLeafCharacter(s.at(node->vnode_seqid));
        }
    }

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

