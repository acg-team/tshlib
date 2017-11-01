//
// Created by Lorenzo Gatti on 26/10/17.
//

#ifndef TSHLIB_UTREE_HPP
#define TSHLIB_UTREE_HPP

#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/src/Core/IO.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "PhyTree.hpp"


/*
 * VirtualNode contains the virtual directions for traversing the tree. Each VirtualNode represents
 * either an internal node or a leaf.
 *
 *         [VN]
 *          ||
 *          ||
 *  +----------------+
 *  |       up       |
 *  |     /    \     |  Node of a real phylogenetic tree
 *  | left  -  right |
 *  +----------------+
 *    //          \\
 *   //            \\
 * [VN]            [VN]
 */

class VirtualNode {
private:

public:

    int vnode_id;                                   /* Node ID - Useful in case of parallel independent executions */
    std::string vnode_name;                         /* Node Name - Useful in case of parallel independent executions */
    double vnode_branchlength;                      /* Branch length connecting the node the parent node */
    double vnode_iota;                              /* PIP Iota value computed on the branch connecting the node to the parent */
    double vnode_beta;                              /* PIP Beta value computed on the branch connecting the node to the parent */
    Eigen::MatrixXd vnode_Pr;                       /* Pr matrix computed recursively */
    double vnode_tau;                               /* Tau value up to this node */
    double vnode_nu;                                /* Nu value associated to this node */
    std::vector<Eigen::VectorXd> vnode_Fv;          /* Fv matrix computed recursively */
    bool vnode_setA;                                /* Flag: Include node in computing the set A -- might be not necessary */
    int vnode_descCount;                            /*  ? */
    char vnode_character;                           /* Character associated to this node (only if terminal node) */
    int vnode_depth;                                /* Depth level of the node in the tree */
    bool vnode_leaf;                                /* Flag: terminal node in the tree topology */

    VirtualNode();

    ~VirtualNode();

    virtual void addMember(VirtualNode *inNode);

    virtual void RotateClockwise();

    virtual void RotateCounterClockwise();

    void getMemberNeighbors(int radius);


    bool isTerminalNode();


protected:
    VirtualNode *vnode_up;                          /* NodeUp - This is the pointer to the VirtualNode above */
    VirtualNode *vnode_left;                        /* NodeLeft  -  This is the pointer to the VirtualNode on the leftside */
    VirtualNode *vnode_right;                       /* NodeRight - This is the pointer to the VirtualNode on the rightside */

    virtual void setNodeRight(VirtualNode inNode);

    virtual void setNodeLeft(VirtualNode inNode);

    virtual void setNodeUp(VirtualNode inNode);

    VirtualNode getNodeLeft();

    VirtualNode getNodeRight();

    VirtualNode getNodeUp();

};

template<class NodeType>
class VirtualLeaf : public VirtualNode {

};

template<class NodeType>
class VirtualInternalNode : public VirtualNode {
public:

};

class Utree {
public:

    VirtualNode *utree_start_node;

    Utree();

    ~Utree();

    void setPseudoRoot();

    void removePseudoRoot();

    void FindPseudoRoot();

    void printTree();

    void printTreeNewick();


protected:

private:


};

void createUtree(PhyTree *in_tree, Utree *out_tree);


#endif //TSHLIB_UTREE_HPP
