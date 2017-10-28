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
 * @file main.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 10 2017
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
 * @see For more information visit:
 */
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#define LOGURU_IMPLEMENTATION 1
#define LOGURU_WITH_STREAMS 1

#include <loguru.hpp>

#include "PhyTree.hpp"
#include "TreeRearrangment.hpp"
#include "Alignment.hpp"
#include "Likelihood.hpp"
#include "newick.hpp"


std::string utree_formatNewickR(node *n, bool is_root) {

    if (n->next == NULL) {
        return n->data->getName();
    } else {
        std::stringstream newick;
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

std::string utree_formatNewick(node *utree_pseudo_root) {
    std::string s;

    if (utree_pseudo_root->next == NULL) {
        return NULL;
    }

    s = utree_formatNewickR(utree_pseudo_root, true) + ";";

    return s;
}

void print_node_neighbours(node *n) {

    std::string descnode;
    descnode += n->data->getName() + " ";

    if (n->next != NULL) {
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

    if (n->next != NULL) {
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }

}

void print_utree(node *n) {

    print_node_neighbours(n);

    if (n->next != NULL) {
        print_utree_rec(n->back);
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }
}

void utree_nodes_within_radius(node *start_node, node *new_node, int radius, std::vector<utree_move_info> &list_nodes) {

    utree_move_info m;
    m.node1 = start_node;
    if (new_node->next != NULL) {
        new_node = new_node->next;
    }
    m.node2 = new_node;
    list_nodes.push_back(m);

    if (radius <= 0) {
        return;
    }

    if (new_node->next != NULL) {
        radius--;
        utree_nodes_within_radius(start_node, new_node->back, radius, list_nodes);
        utree_nodes_within_radius(start_node, new_node->next->back, radius, list_nodes);
    }

}

void utree_get_list_nodes_within_radius(node *n,
                                        int radius,
                                        std::vector<utree_move_info> &list_nodes_left,
                                        std::vector<utree_move_info> &list_nodes_right,
                                        std::vector<utree_move_info> &list_nodes_up) {

    if (n->next != NULL) {
        utree_nodes_within_radius(n, n->back, radius, list_nodes_up);
        utree_nodes_within_radius(n, n->next->back, radius, list_nodes_left);
        utree_nodes_within_radius(n, n->next->next->back, radius, list_nodes_right);
    }

}

void copy_vector(std::vector<node *> &dest, std::vector<node *> &source) {

    for (unsigned int i = 0; i < source.size(); i++) {
        node *n = new node;
        node *m = source.at(i);
        n->next = m->next;
        n->back = m->back;
        n->data = m->data;
        n->ID = m->ID;
        dest.push_back(n);
    }

}


void SPR_move(PhyTree *tree, std::vector<node *> &utree, node *source, node *target, int file_tree_idx) {
    node *p_child_1;
    node *p_child_2;
    node *q_child;
    bool valid_move;
    FILE *fid;
    char tree_filename[80];
    std::string ss;

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


int main(int argc, char **argv) {
    loguru::init(argc, argv);
    LOG_S(INFO) << "test";
    std::string tree_file = argv[1];
    PhyTree *tree = nullptr;
    double mu;
    double lambda;
    double tau;
    double nu;

    mu = 0.1;
    lambda = 0.2;

    //------------------------------------------------------------------------------------------------------------------
    // INIT TREE

    // tree filename
    std::ifstream tree_str(tree_file.c_str());

    // read newick file
    tree = newick_parser::parse_newick(&tree_str);

    // set name of internal nodes
    tree->set_missing_node_name("V");

    // compute total tree length
    tau = tree->computeLength();

    // compute the normalizing Poisson intensity
    nu = compute_nu(tau, lambda, mu);

    // set insertion probability to each node
    tree->set_iota(tau, mu);

    // set survival probability to each node
    tree->set_beta(tau, mu);

    // print newick tree
    LOG_S(INFO) << tree->formatNewick();
    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA

    //Alignment ;
    auto *alignment = new Alignment;

    std::vector<std::pair<std::string, std::string> > MSA;

    //TODO: Create a class MSA (+ weight per each column as in codonPhyML)
    std::string seq1_label = "A";
    std::string seq2_label = "B";
    std::string seq3_label = "C";
    std::string seq4_label = "D";
    std::string seq5_label = "E";

    std::string seq1_DNA = "A-C";
    std::string seq2_DNA = "-TC";
    std::string seq3_DNA = "-TC";
    std::string seq4_DNA = "-CC";
    std::string seq5_DNA = "-CG";


    MSA.emplace_back(seq1_label, seq1_DNA);
    MSA.emplace_back(seq2_label, seq2_DNA);
    MSA.emplace_back(seq3_label, seq3_DNA);
    MSA.emplace_back(seq4_label, seq4_DNA);
    MSA.emplace_back(seq5_label, seq5_DNA);

    alignment->addSequence(seq1_label, seq1_DNA);
    alignment->addSequence(seq2_label, seq2_DNA);
    alignment->addSequence(seq3_label, seq3_DNA);
    alignment->addSequence(seq4_label, seq4_DNA);
    alignment->addSequence(seq5_label, seq5_DNA);

    LOG_S(INFO) << alignment->align_dataset.size();

    //----------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION

    int is_DNA_AA_Codon;
    int extended_alphabet_size;
    int num_leaves;
    unsigned long MSA_len;
    double log_col_lk;
    double logLK;
    Eigen::VectorXd pi;
    double p0;

    num_leaves = MSA.size();

    // 1:DNA, 2:AA, 3:Codon
    is_DNA_AA_Codon = 1;

    // DNA alphabet
    extended_alphabet_size = 5;

    // set "pseudo" probability matrix
    tree->tmp_initPr(extended_alphabet_size); //TODO: pass Q from codonPhyML (?)

    // set Pi, steady state frequencies
    pi = Eigen::VectorXd::Zero(extended_alphabet_size);
    pi[0] = 0.25;
    pi[1] = 0.25;
    pi[2] = 0.25;
    pi[3] = 0.25;
    pi[4] = 0.0;

    // get MSA length
    MSA_len = static_cast<unsigned long>(alignment->getAlignmentSize());

    logLK = 0.0;
    // compute log_col_lk
    for (int i = 0; i < MSA_len; i++) {

        // extract MSA column
        std::string s = create_col_MSA(MSA, i);

        // assign char at the leaves
        tree->set_leaf_state(s);

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        tree->set_ancestral_flag(s);

        tree->clear_fv();

        // compute column likelihood
        //TODO: Add weight per column
        log_col_lk = compute_col_lk(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);

        LOG_S(INFO) << "log_col_lk=" << log_col_lk;

        logLK += log_col_lk;
    }

    // compute empty column likelihood
    p0 = compute_log_lk_empty_col(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);

    LOG_S(INFO) << "p0=" << p0;

    logLK += phi(MSA_len, nu, p0);
    //------------------------------------------------------------------------------------------------------------------
    // BUILD UNROOTED TREE

    std::vector<node *> utree;
    node *utree_pseudo_root;

    utree = tree->create_unrooted_tree(tree, num_leaves);

    utree_pseudo_root = utree.at(0);

    //------------------------------------------------------------------------------------------------------------------
    // GET ALL NODES WITHIN RADIUS
    // TODO: implement high-level heuristic to decide which TS to use according to node level

    int radius;

    radius = 3;

    //----------------------------------------------------
    std::vector<utree_move_info> utree_nni_spr_stack_left;
    std::vector<utree_move_info> utree_nni_spr_stack_right;
    std::vector<utree_move_info> utree_nni_spr_stack_up;
    node *u_node;

    u_node = utree_pseudo_root;

    //---------------------------------------------------------------
    u_node = utree.at(0);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(4);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(7);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(3);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(10);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(11);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(12);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);

    u_node = utree.at(13);
    LOG_S(INFO) << "UNODE:"; //<<u_node->data->getName()<<"";
    print_node_neighbours(u_node);
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);
    //---------------------------------------------------------------
    LOG_S(INFO) << "size list_left:" << utree_nni_spr_stack_left.size();
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        LOG_S(INFO) << "list[" << i << "]=(" << (utree_nni_spr_stack_left.at(i)).node1->data->getName() << ";"
                    << (utree_nni_spr_stack_left.at(i)).node2->data->getName() << ")";
    }

    LOG_S(INFO) << "size list_right:" << utree_nni_spr_stack_right.size();
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        LOG_S(INFO) << "list[" << i << "]=(" << (utree_nni_spr_stack_right.at(i)).node1->data->getName() << ";"
                    << (utree_nni_spr_stack_right.at(i)).node2->data->getName() << ")";
    }

    LOG_S(INFO) << "size list_up:" << utree_nni_spr_stack_up.size();
    for (unsigned int i = 0; i < utree_nni_spr_stack_up.size(); i++) {
        LOG_S(INFO) << "list[" << i << "]=(" << (utree_nni_spr_stack_up.at(i)).node1->data->getName() << ";"
                    << (utree_nni_spr_stack_up.at(i)).node2->data->getName() << ")";
    }

    //------------------------------------------------------------------------------------------------------------------
    // PERFORM SPR MOVES and RECOMPUTE logLK

    int max_idx;
    double max_val;
    std::vector<PhyTree *> p;
    bool valid_move;
    utree_move_info un;

    node *p_child_1;
    node *p_child_2;
    node *q_child;

    FILE *fid;
    int file_tree_idx;
    char tree_filename[80];
    std::string ss;

    node *source;
    node *target;

    //---------------------------------------------------------
    file_tree_idx = 0;
    sprintf(tree_filename, "%s_%d.nwk", "../data/out/tree", file_tree_idx);
    fid = fopen(tree_filename, "w");
    ss = utree_formatNewick(utree.at(0));
    fprintf(fid, "%s", ss.c_str());
    fclose(fid);
    //---------------------------------------------------------

    max_idx = -1;
    max_val = -INFINITY;
    //---------------------------------------------------------------------
    LOG_S(INFO) << "processing left";
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        un = utree_nni_spr_stack_left.at(i);

        source = un.node1;
        target = un.node2;

        file_tree_idx++;
        SPR_move(tree, utree, source, target, file_tree_idx);
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        un = utree_nni_spr_stack_left.at(i);

        source = un.node1;
        if (un.node2->next != NULL) {
            target = un.node2->next;

            file_tree_idx++;
            SPR_move(tree, utree, source, target, file_tree_idx);
        }
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        un = utree_nni_spr_stack_left.at(i);

        source = un.node1->next->next;
        target = un.node2;

        file_tree_idx++;
        SPR_move(tree, utree, source, target, file_tree_idx);
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        un = utree_nni_spr_stack_left.at(i);

        source = un.node1->next->next;
        if (un.node2->next != NULL) {
            target = un.node2->next;

            file_tree_idx++;
            SPR_move(tree, utree, source, target, file_tree_idx);
        }
    }
    //---------------------------------------------------------------------
    LOG_S(INFO) << "processing right";
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        un = utree_nni_spr_stack_right.at(i);

        source = un.node1;
        target = un.node2;

        file_tree_idx++;
        SPR_move(tree, utree, source, target, file_tree_idx);
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        un = utree_nni_spr_stack_right.at(i);

        source = un.node1;
        if (un.node2->next != NULL) {
            target = un.node2->next;

            file_tree_idx++;
            SPR_move(tree, utree, source, target, file_tree_idx);
        }
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        un = utree_nni_spr_stack_right.at(i);

        source = un.node1->next;
        target = un.node2;

        file_tree_idx++;
        SPR_move(tree, utree, source, target, file_tree_idx);
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        un = utree_nni_spr_stack_right.at(i);

        source = un.node1->next;
        if (un.node2->next != NULL) {
            target = un.node2->next;

            file_tree_idx++;
            SPR_move(tree, utree, source, target, file_tree_idx);
        }
    }
    //---------------------------------------------------------------------
    LOG_S(INFO) << "processing up";
    for (unsigned int i = 0; i < utree_nni_spr_stack_up.size(); i++) {
        un = utree_nni_spr_stack_up.at(i);

        source = un.node1->next;
        target = un.node2;

        file_tree_idx++;
        SPR_move(tree, utree, source, target, file_tree_idx);
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_up.size(); i++) {
        un = utree_nni_spr_stack_up.at(i);

        source = un.node1->next;
        if (un.node2->next != NULL) {
            target = un.node2->next;

            file_tree_idx++;
            SPR_move(tree, utree, source, target, file_tree_idx);
        }
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_up.size(); i++) {
        un = utree_nni_spr_stack_up.at(i);

        source = un.node1->next->next;
        target = un.node2;

        file_tree_idx++;
        SPR_move(tree, utree, source, target, file_tree_idx);
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_up.size(); i++) {
        un = utree_nni_spr_stack_up.at(i);

        source = un.node1->next->next;
        if (un.node2->next != NULL) {
            target = un.node2->next;

            file_tree_idx++;
            SPR_move(tree, utree, source, target, file_tree_idx);
        }
    }


    LOG_S(INFO) << "max_val:" << max_val << " at index: " << max_idx;

    utree_nni_spr_stack_left.empty();
    utree_nni_spr_stack_right.empty();
    utree_nni_spr_stack_up.empty();

    exit(0);
}
