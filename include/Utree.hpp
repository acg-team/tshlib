/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
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
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Utree.hpp
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

#ifndef TSHLIB_UTREE_HPP
#define TSHLIB_UTREE_HPP

#include <string>
#include <vector>


#include "Utilities.hpp"

//#include "Alignment.hpp"

namespace tshlib {

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
        int vnode_depth;                                /* Depth level of the node in the tree */
        bool vnode_leaf;                                /* Flag: terminal node in the tree listVNodes */
        int vnode_move_direction;                       /* Int: This attribute is used to perform the correct rotation of the p-node w.r.t q-node. */
        NodeRotation vnode_rotated;

        NodeRotation getVnode_rotated() const;

        /* Flag: if node was rotaded during a tree rearrangement move */
        int vnode_seqid;                                /* Needed to associate the seqID to the nodeID

        /*!
         *  Standard constructor
         */
        VirtualNode();

        int getVnode_id() const;

        void setVnode_id(int vnode_id);

        VirtualNode(const VirtualNode &inNode);

        /*!
         * Virtual deconstructor
         */
         ~VirtualNode();

        /*!
         * @brief This function connects the current node to another one. It automatically performs a bidirectional connection
         * @param inVNode Target node to apply the connection to
         */
        void connectNode(VirtualNode *inVNode);


        void setBranchLength(double blength);


        /*!
         * @brief This function disconnect the current node from any other one above it. The function is bidirectional.
         */
        void disconnectNode();

        /*!
         * @brief This function perfomrms a one-step clockwise rotation of the virtualnode pointers
         */
        void rotateClockwise();

        void rotateClockwise(bool revertRotations);

        /*!
         * @brief This function perfomrms a one-step counter-clockwise rotation of the virtualnode pointers
         */
        void rotateCounterClockwise();

        void rotateCounterClockwise(bool revertRotations);

        /*!
         * @brief This function resets any rotation previously performed on the node
         */
        void resetNodeDirections(bool revertRotations);

        //void getMemberNeighbors(int radius);

        /*!
         * @brief The function prints the neighborhood of the node in the format <^nodeUp;nodeLeft;nodeRight >
         * @return std::string with the node neighborhood
         */
        std::string printNeighbours();

        void setNodeName(std::string s);

        const std::string getNodeName();

        VirtualNode *getNodeUp();

        VirtualNode *getNodeLeft();

        VirtualNode *getNodeRight();

        void _traverseVirtualNodeTree();

        double computeTotalTreeLength();

        VirtualNode *getSiblingNode();

        int getNodeLevel();

        void setNodeLevel(int level);

        void clearChildren();

        /*!
         * @brief This function returns true if the node is terminal
         * @return boolean value (true or false)
         */
        bool isTerminalNode();

        /*!
         * @brief This function returns true if the node is root
         * @return boolean value (true or false)
         */
        bool isRootNode();

        /*!
         * @brief This function returns true if the node is pseudoroot
         * @return boolean value (true or false)
         */
        bool isPseudoRootNode();

        /*!
         * @brief This function return the index of node as seen from the parent immediate above it.
         * @return
         */
        NodePosition indexOf();

        /*!
         * @brief This function swaps the current node and another one passed in the argument. The swap takes in consideration the topology of the tree.
         * @param targetNode
         * @return boolean value if the execution was performed correctly.
         */
        bool swapNode(VirtualNode *targetNode, MoveDirections move_direction, bool revertRotations);

        /*!
         * @brief The function checks if the current node is a parent of another node using a recursive structure
         * @param inVNode VirtualNode pointer
         * @return False if the node passed is not parent of the current one, True otherwise
         * @deprecated
         */
        bool isParent(VirtualNode *inVNode);

        void _setNodeRight(VirtualNode *inVNode);

        void _setNodeLeft(VirtualNode *inVNode);

        void _setNodeUp(VirtualNode *inVNode);

        void _bidirectionalUpwardConnection(VirtualNode *inNode);

        NodeRotation getNodeRotation() const {
            return vnode_rotated;
        }

    protected:

        VirtualNode *vnode_up;                          /* NodeUp - This is the pointer to the VirtualNode above */
        VirtualNode *vnode_left;                        /* NodeLeft  -  This is the pointer to the VirtualNode on the leftside */
        VirtualNode *vnode_right;                       /* NodeRight - This is the pointer to the VirtualNode on the rightside */


        void _oneway_connectNode(VirtualNode *inVNode);

        void _recursive_cw_rotation(VirtualNode *vnode, bool revertRotations);

        void _recursive_ccw_rotation(VirtualNode *vnode, bool revertRotations);


    };


    class Utree {
    private:
        bool initialized_treeDepth;

    public:

        std::vector<VirtualNode *> listVNodes;
        std::vector<VirtualNode *> startVNodes;
        VirtualNode *rootnode;

        Utree();

        Utree(const Utree &rhs) { /* copy construction from rhs*/ }

        Utree &operator=(const Utree &rhs) {};

        ~Utree();

        /*!
         * @brief This function finds the pseudoroot traversing the tree from starting node until a bidirectional connection is found
         * @param inVNode starting node
         * @return std::vector of VirtualNode pointers from the starting point until the pseudoroot
         */
        std::vector<VirtualNode *> findPseudoRoot(VirtualNode *inVNode, bool fixPseudoRootOnNextSubtree = false);

        /*!
         * @brief This function adds a complete VirtualNode to the listVNodes attribute of the class Utree
         * @param inVNode       Pointer to the VirtualNode to add
         * @param isStartNode
         */
        void addMember(VirtualNode *inVNode, bool isStartNode = false);

        /*!
         * @brief
         * @param showInternalNodeNames
         * @return
         */
        std::string printTreeNewick(bool showInternalNodeNames);

        /*!
         * @brief
         * @return
         */
        int getMaxNodeDistance();

        /*!
         * @brief
         * @param outfilepath
         */
        void saveTreeOnFile(std::string outfilepath);

        /*!
         * @brief
         */
        void printAllNodesNeighbors();

        /*!
         * @brief This function create a virtual root on the utree object. It breaks the link between two pseudoroot virtualnodes
         */
        void addVirtualRootNode();

        /*!
         * @brief This function remove the virtual root node added with the companion function addVirtualRootNode()
         */
        void removeVirtualRootNode();

        void _testReachingPseudoRoot();

        double computeTotalTreeLength();

        void _printUtree();

        std::vector<VirtualNode *> computePathBetweenNodes(VirtualNode *vnode_1, VirtualNode *vnode_2);

        std::vector<VirtualNode *> _unique(std::vector<VirtualNode *> &list_nodes_n1, std::vector<VirtualNode *> &list_nodes_n2);

        std::vector<VirtualNode *> getPostOrderNodeList();

        std::vector<VirtualNode *> getPostOrderNodeList(VirtualNode *startNode);

        void _getPostOrderNodeList(std::vector<VirtualNode *> &rlist, VirtualNode *node);

        void computeTreeDepth();

        int getTreeDepthAtNode(VirtualNode *vnode);

        void computeNodeDepth(VirtualNode *vnode);

    protected:
        std::string _recursiveFormatNewick(VirtualNode *vnode, bool showInternalNodeNames);


        void _updateStartNodes();

    private:

    };



    namespace VirtualNodeUtils{

        void rotateNodeClockwise(VirtualNode *vnode);

        void rotateNodeCounterClockwise(VirtualNode *vnode);
    }


    namespace UtreeUtils {
        using namespace tshlib;

        //void _traverseTree(Utree *in_tree, VirtualNode *target, PhyTree *source);

        //void convertUtree(PhyTree *in_tree, Utree *out_tree);

        /*!
         * @brief This function has been ported in Utree::_unique
         * @param n
         * @return
         */
        //std::vector<VirtualNode *> get_unique(std::vector<VirtualNode *> &list_nodes_n1, std::vector<VirtualNode *> &list_nodes_n2);


        //std::vector<VirtualNode *> fill_with_nodes(VirtualNode *n);

    /*!
         * @brief This function has been ported in Utree::computePathBetweenNodes
         * @param n
         * @return
         */
        //std::vector<VirtualNode *> get_path_from_nodes(VirtualNode *vn1, VirtualNode *vn2);

        //void associateNode2Alignment(Alignment *inMSA, Utree *inTree);

        VirtualNode *getPseudoRoot(VirtualNode *vn);

    }
}
#endif //TSHLIB_UTREE_HPP
