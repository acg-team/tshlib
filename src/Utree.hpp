//
// Created by Lorenzo Gatti on 26/10/17.
//

#ifndef TSHEXE_UTREE_H
#define TSHEXE_UTREE_H

#include <string>
#include "PhyTree.hpp"
#include "TreeRearrangment.hpp"


class Utree {

};


std::string utree_formatNewickR(node *n, bool is_root);

std::string utree_formatNewick(node *utree_pseudo_root);

void print_utree(node *n);

void print_utree_rec(node *n);

void print_node_neighbours(node *n);

void utree_nodes_within_radius(node *start_node, node *new_node, int radius,
                               std::vector<utree_move_info> &list_nodes);

void utree_get_list_nodes_within_radius(node *n,
                                        int radius,
                                        std::vector<utree_move_info> &list_nodes_left,
                                        std::vector<utree_move_info> &list_nodes_right,
                                        std::vector<utree_move_info> &list_nodes_up);

void copy_vector(std::vector<node *> &dest, std::vector<node *> &source);


#endif //TSHEXE_UTREE_H
