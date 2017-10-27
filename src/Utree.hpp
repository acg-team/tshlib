//
// Created by Lorenzo Gatti on 26/10/17.
//

#ifndef TSHEXE_UTREE_H
#define TSHEXE_UTREE_H

#include <string>
#include "PhyTree.hpp"
#include "TreeRearrangment.hpp"


/*
 *          r
 *       /    \
 *   left  -  right
 */

class VirtualNode {
private:
public:
    VirtualNode *vnode_root;
    VirtualNode *vnode_left;
    VirtualNode *vnode_right;
    int vnode_id;
    std::string vnode_name;
    double vnode_branchlength;
    double vnode_iota;
    double vnode_beta;
    Eigen::MatrixXd vnode_Pr;
    double vnode_tau;
    double vnode_nu;
    std::vector<Eigen::VectorXd> vnode_MSA_fv;
    bool vnode_setA;
    int vnode_descCount;
    char vnode_character;
    int vnode_depth;

protected:
};

class Utree {
public:

    Utree();

    ~Utree();

    void setPseudoRoot();

    void removePseudoRoot();

    void addMember();

    void removeMember();

    void swapMembers();

    void getMemberNeighbors(int radius);

    void printTree();

protected:
private:
    VirtualNode utree_start_node;
    std::vector<VirtualNode *> members;
};
/*

std::string utree_formatNewickR(node *n, bool is_root);

std::string utree_formatNewick(node *utree_pseudo_root);

void print_utree(node *n);

void print_utree_rec(node *n);

void print_node_neighbours(node *n);

void utree_nodes_within_radius(node *start_node, node *new_node, int radius, std::vector<utree_move_info> &list_nodes);

void utree_get_list_nodes_within_radius(node *n,
                                        int radius,
                                        std::vector<utree_move_info> &list_nodes_left,
                                        std::vector<utree_move_info> &list_nodes_right,
                                        std::vector<utree_move_info> &list_nodes_up);

void copy_vector(std::vector<node *> &dest, std::vector<node *> &source);
*/


#endif //TSHEXE_UTREE_H
