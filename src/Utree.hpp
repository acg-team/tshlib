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
#include "Utilities.hpp"


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
    bool vnode_leaf;                                /* Flag: terminal node in the tree listVNodes */
    int vnode_move_direction;                       /* Int: This attribute is used to perform the correct rotation of the p-node w.r.t q-node. */
    int vnode_rotated;                              /* Flag: if node was rotaded during a tree rearrangement move */

    /*!
     *  Standard constructor
     */
    VirtualNode();

    /*!
     * Virtual deconstructor
     */
    virtual ~VirtualNode();

    /*!
     * @brief This function connects the current node to another one. It automatically performs a bidirectional connection
     * @param inVNode Target node to apply the connection to
     */
    virtual void connectNode(VirtualNode *inVNode);

    /*!
     * @brief This function disconnect the current node from any other one above it. The function is bidirectional.
     */
    virtual void disconnectNode();

    /*!
     * @brief This function perfomrms a one-step clockwise rotation of the virtualnode pointers
     */
    virtual void rotateClockwise();

    virtual void rotateClockwise(bool revertRotations);

    /*!
     * @brief This function perfomrms a one-step counter-clockwise rotation of the virtualnode pointers
     */
    virtual void rotateCounterClockwise();

    virtual void rotateCounterClockwise(bool revertRotations);

    /*!
     * @brief This function resets any rotation previously performed on the node
     */
    virtual void resetNodeDirections();

    //void getMemberNeighbors(int radius);

    /*!
     * @brief The function prints the neighborhood of the node in the format <^nodeUp;nodeLeft;nodeRight >
     * @return std::string with the node neighborhood
     */
    std::string printNeighbours();


    virtual void setNodeRight(VirtualNode *inVNode);

    virtual void setNodeLeft(VirtualNode *inVNode);

    virtual void setNodeUp(VirtualNode *inVNode);

    VirtualNode *getNodeUp();

    VirtualNode *getNodeLeft();

    VirtualNode *getNodeRight();

    /*!
     * @brief This function returns true if the node is terminal
     * @return boolean value (true or false)
     */
    bool isTerminalNode();

    /*!
     * @brief This function return the index of node as seen from the parent immediate above it.
     * @return
     */
    int indexOf();

    /*!
     * @brief This function swaps the current node and another one passed in the argument. The swap takes in consideration the topology of the tree.
     * @param targetNode
     * @return boolean value if the execution was performed correctly.
     */
    bool swapNode(VirtualNode *targetNode, MoveDirections move_direction);

    /*!
     * @brief The function checks if the current node is a parent of another node using a recursive structure
     * @param inVNode VirtualNode pointer
     * @return False if the node passed is not parent of the current one, True otherwise
     */
    bool isParent(VirtualNode *inVNode);



protected:
    VirtualNode *vnode_up;                          /* NodeUp - This is the pointer to the VirtualNode above */
    VirtualNode *vnode_left;                        /* NodeLeft  -  This is the pointer to the VirtualNode on the leftside */
    VirtualNode *vnode_right;                       /* NodeRight - This is the pointer to the VirtualNode on the rightside */

    virtual void _oneway_connectNode(VirtualNode *inVNode);

    virtual void _recursive_cw_rotation(VirtualNode *vnode, bool revertRotations);

    virtual void _recursive_ccw_rotation(VirtualNode *vnode, bool revertRotations);

};


class Utree {
public:

    std::vector<VirtualNode *> listVNodes;
    std::vector<VirtualNode *> startVNodes;

    Utree();

    virtual ~Utree();

    /*!
     * @brief This function finds the pseudoroot traversing the tree from starting node until a bidirectional connection is found
     * @param inVNode starting node
     * @return std::vector of VirtualNode pointers from the starting point until the pseudoroot
     */
    std::vector<VirtualNode *> findPseudoRoot(VirtualNode *inVNode, bool fixPseudoRootOnNextSubtree = false);

    virtual void addMember(VirtualNode *inVNode, bool isStartNode = false);

    std::string printTreeNewick(bool showInternalNodeNames);

    virtual void saveTreeOnFile(std::string outfilepath);

    virtual void _testReachingPseudoRoot();

protected:
    std::string _recursiveFormatNewick(VirtualNode *vnode, bool showInternalNodeNames);

    virtual void _updateStartNodes();

private:

};


namespace UtreeUtils {

    void _traverseTree(Utree *in_tree, VirtualNode *target, PhyTree *source);

    void convertUtree(PhyTree *in_tree, Utree *out_tree);
}

#endif //TSHLIB_UTREE_HPP
