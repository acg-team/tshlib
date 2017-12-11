/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti
 *******************************************************************************
 *
 * This file is part of tshexe
 *
 * tshexe is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tshexe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with likpip. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file deprecated.cpp
 * @author Lorenzo Gatti
 * @date 11 11 2017
 * @version 1.0
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit: http://www.lorenzogatti.me
 */
#include <TreeRearrangment.hpp>
#include <Utree.hpp>
#include <loguru.hpp>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "deprecated.hpp"

std::__1::string utree_formatNewickR(node *n, bool is_root) {

    if (n->next == nullptr) {
        return n->data->getName();
    } else {
        std::__1::stringstream newick;
        if (is_root) {
            newick << "(";
            newick << utree_formatNewickR(n->back, false) << ":" << n->back->data->getBranchLength();
            newick << ",";
            newick << utree_formatNewickR(n->next->back, false) << ":" << n->next->back->data->getBranchLength();
            newick << ",";
            newick << utree_formatNewickR(n->next->next->back, false) << ":" << n->next->next->back->data->getBranchLength();
            newick << ")";
        } else {
            newick << "(";
            newick << utree_formatNewickR(n->next->back, false) << ":" << n->next->back->data->getBranchLength();
            newick << ",";
            newick << utree_formatNewickR(n->next->next->back, false) << ":" << n->next->next->back->data->getBranchLength();
            newick << ")";
        }

        return newick.str();
    }

}

std::__1::string utree_formatNewick(node *utree_pseudo_root) {
    std::__1::string s;

    if (utree_pseudo_root->next == nullptr) {
        return nullptr;
    }

    s = utree_formatNewickR(utree_pseudo_root, true) + ";";

    return s;
}

void print_node_neighbours(node *n) {

    std::__1::string descnode;
    descnode += n->data->getName() + " ";

    if (n->next != nullptr) {
        descnode += "(^" + n->back->data->getName() + ";";
        descnode += "<" + n->next->back->data->getName() + ";";
        descnode += n->next->next->back->data->getName() + ">)";

    } else {
        descnode += "(^" + n->back->data->getName() + "; <-;->)";
    }

    LOG_S(INFO) << descnode;
}

void print_utree_rec(node *n) {

    print_node_neighbours(n);

    if (n->next != nullptr) {
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }

}

/*!
 * @brief This function prints the unrooted tree created by lists of struct node
 * @deprecated This function is not more in use
 * @param n
 */
void print_utree(node *n) {

    print_node_neighbours(n);

    if (n->next != nullptr) {
        print_utree_rec(n->back);
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }
}

void utree_nodes_within_radius(node *start_node, node *new_node, int radius, std::__1::vector<utree_move_info> &list_nodes) {

    utree_move_info m;
    m.node1 = start_node;
    if (new_node->next != nullptr) {
        new_node = new_node->next;
    }
    m.node2 = new_node;
    list_nodes.push_back(m);

    if (radius <= 0) {
        return;
    }

    if (new_node->next != nullptr) {
        radius--;
        utree_nodes_within_radius(start_node, new_node->back, radius, list_nodes);
        utree_nodes_within_radius(start_node, new_node->next->back, radius, list_nodes);
    }

}

/*!
 * @brief This function returns the list of nodes withtin a defined radius starting from a pnode
 * @deprecated This function is not more in use
 * @param n
 * @param radius
 * @param list_nodes_left
 * @param list_nodes_right
 * @param list_nodes_up
 */
void utree_get_list_nodes_within_radius(node *n,
                                        int radius,
                                        std::__1::vector<utree_move_info> &list_nodes_left,
                                        std::__1::vector<utree_move_info> &list_nodes_right,
                                        std::__1::vector<utree_move_info> &list_nodes_up) {

    if (n->next != nullptr) {
        utree_nodes_within_radius(n, n->back, radius, list_nodes_up);
        utree_nodes_within_radius(n, n->next->back, radius, list_nodes_left);
        utree_nodes_within_radius(n, n->next->next->back, radius, list_nodes_right);
    }

}

void copy_vector(std::__1::vector<node *> &dest, std::__1::vector<node *> &source) {

    for (auto m : source) {
        node *n = new node;
        n->next = m->next;
        n->back = m->back;
        n->data = m->data;
        n->ID = m->ID;
        dest.push_back(n);
    }

}

/*!
 * @brief This function performs NNI and SPR moves using utree (node) structures
 * @deprecated
 * @param tree
 * @param utree
 * @param source
 * @param target
 * @param file_tree_idx
 */

void SPR_move(PhyTree *tree, std::__1::vector<node *> &utree, node *source, node *target, int file_tree_idx) {
    node *p_child_1;
    node *p_child_2;
    node *q_child;
    bool valid_move;
    FILE *fid;
    char tree_filename[80];
    std::__1::string ss;

    LOG_S(INFO) << "node_1:";
    print_node_neighbours(source);
    LOG_S(INFO) << "node_2:";
    print_node_neighbours(target);
    //LOG_S(INFO)<<"";


    p_child_1 = source->next->back;
    p_child_2 = source->next->next->back;
    q_child = target->back;

    valid_move = tree->utree_swap(source, target);

    if (valid_move) {

        LOG_S(INFO) << "-------------";
        LOG_S(INFO) << utree_formatNewick(utree.at(0));
        LOG_S(INFO) << "-------------";


        //---------------------------------------------------------
        //file_tree_idx++;
        sprintf(tree_filename, "%s_%d.nwk", "../data/out/tree", file_tree_idx);
        LOG_S(INFO) << tree_filename;
        fid = fopen(tree_filename, "w");
        ss = utree_formatNewick(utree.at(0));
        fprintf(fid, "%s", ss.c_str());
        fclose(fid);
        //---------------------------------------------------------


        source->next->back = p_child_1;
        p_child_1->back = source->next;
        source->next->next->back = p_child_2;
        p_child_2->back = source->next->next;
        target->back = q_child;
        q_child->back = target;

        LOG_S(INFO) << "****************";
        LOG_S(INFO) << utree_formatNewick(utree.at(0));
        LOG_S(INFO) << "****************";

    } else {
        LOG_S(INFO) << "I am skipping this...";
    }

}


void nodes_within_radius(PhyTree *start_node, PhyTree *node, int radius, std::vector<move_info> &list_nodes) {

//    if (!save) {
//          save = true;
//    } else {
    move_info m;
    m.node1 = start_node;
    m.node2 = node;
    list_nodes.push_back(m);
//    }

    if (radius <= 0) {
        return;
    }

    if (!node->isLeaf()) {
        radius--;
        nodes_within_radius(start_node, node->get_left_child(), radius, list_nodes);
        nodes_within_radius(start_node, node->get_right_child(), radius, list_nodes);
    }

}

void nodes_within_radius_up(PhyTree *start_node, PhyTree *node, int radius, int direction,
                            std::vector<move_info> &list_nodes) {
    unsigned int idx;

    //TODO: check binary tree condition!

    move_info m;
    m.node1 = start_node;
    m.node2 = node;
    list_nodes.push_back(m);

    if (radius <= 0) {
        return;
    }

    radius--;
    if (direction == 0) {
        if (node->getParent() != NULL) {
            idx = node->indexOf();
            nodes_within_radius_up(start_node, node->getParent(), radius, idx, list_nodes);
        }
        nodes_within_radius(start_node, node->get_right_child(), radius, list_nodes);
    } else if (direction == 1) {
        if (node->getParent() != NULL) {
            idx = node->indexOf();
            nodes_within_radius_up(start_node, node->getParent(), radius, idx, list_nodes);
        }
        nodes_within_radius(start_node, node->get_left_child(), radius, list_nodes);
    }

}

void
get_list_nodes_within_radius(PhyTree *node, int radius, std::vector<move_info> &list_nodes_left, std::vector<move_info> &list_nodes_right, std::vector<move_info> &list_nodes_up) {
    //bool save;

    //save = false;

    //nodes_within_radius(node, node, radius, save, list_nodes_down);
    //radius--;
    nodes_within_radius(node, node->get_left_child(), radius, list_nodes_left);
    nodes_within_radius(node, node->get_right_child(), radius, list_nodes_right);

    if (node->getParent() != NULL) {
        nodes_within_radius_up(node, node->getParent(), radius, node->indexOf(), list_nodes_up);
    }

}

std::vector<PhyTree *> fill_with_nodes(PhyTree *n) {
    std::vector<PhyTree *> list_nodes_n;

    list_nodes_n.push_back(n);
    while (n->getParent() != NULL) {
        n = n->getParent();
        list_nodes_n.push_back(n);
    }

    return list_nodes_n;
}

std::vector<PhyTree *> get_unique(std::vector<PhyTree *> &list_nodes_n1, std::vector<PhyTree *> &list_nodes_n2) {
    std::vector<PhyTree *> list_nodes;
    PhyTree *n1;
    PhyTree *n2;

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
        n2 = list_nodes_n1.at(list_nodes_n1.size() - 1);
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

std::vector<PhyTree *> get_path_from_nodes(PhyTree *n1, PhyTree *n2) {
    std::vector<PhyTree *> list_nodes_n0;
    std::vector<PhyTree *> list_nodes_n1;
    std::vector<PhyTree *> list_nodes_n2;

    // add nodes from n1 to root
    list_nodes_n1 = fill_with_nodes(n1);

    // add nodes from n2 to root
    list_nodes_n2 = fill_with_nodes(n2);

    list_nodes_n0 = get_unique(list_nodes_n1, list_nodes_n2);

    return list_nodes_n0;
}


std::string create_col_MSA(std::vector<std::pair<std::string, std::string>> &MSA, int index) {
    std::string colMSA;

    for (unsigned int i = 0; i < MSA.size(); i++) {
        colMSA.append(MSA.at(i).second, index, 1);
    }

    return colMSA;
}


