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
#include <iomanip>

#define LOGURU_IMPLEMENTATION 1
#define LOGURU_WITH_STREAMS 1

#include <loguru.hpp>
#include <Utree.hpp>

#include "PhyTree.hpp"
#include "TreeRearrangment.hpp"
#include "Alignment.hpp"
#include "Likelihood.hpp"
#include "newick.hpp"


std::string utree_formatNewickR(node *n, bool is_root) {

    if (n->next == nullptr) {
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

    if (utree_pseudo_root->next == nullptr) {
        return nullptr;
    }

    s = utree_formatNewickR(utree_pseudo_root, true) + ";";

    return s;
}


void print_node_neighbours(node *n) {

    std::string descnode;
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


void print_utree(node *n) {

    print_node_neighbours(n);

    if (n->next != nullptr) {
        print_utree_rec(n->back);
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }
}


void utree_nodes_within_radius(node *start_node, node *new_node, int radius, std::vector<utree_move_info> &list_nodes) {

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


void utree_get_list_nodes_within_radius(node *n,
                                        int radius,
                                        std::vector<utree_move_info> &list_nodes_left,
                                        std::vector<utree_move_info> &list_nodes_right,
                                        std::vector<utree_move_info> &list_nodes_up) {

    if (n->next != nullptr) {
        utree_nodes_within_radius(n, n->back, radius, list_nodes_up);
        utree_nodes_within_radius(n, n->next->back, radius, list_nodes_left);
        utree_nodes_within_radius(n, n->next->next->back, radius, list_nodes_right);
    }

}


void copy_vector(std::vector<node *> &dest, std::vector<node *> &source) {

    for (auto m : source) {
        node *n = new node;
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


namespace fasta_parser {

    std::vector<std::pair<std::string, std::string>> getBuffer(std::string in_fastafile) {
        std::ifstream in_fasta(in_fastafile.c_str());
        std::string line;
        std::string line2;
        std::vector<std::pair<std::string, std::string>> sequences;
        while (in_fasta.is_open()) {
            if (!std::getline(in_fasta, line)) {
                in_fasta.close();
                break;
            }
            if (line.find('>') != std::string::npos) {

                LOG_S(DEBUG2) << "Found header: " << line;

                // parse next line since it must be the sequence
                std::getline(in_fasta, line2);
                //std::pair<std::string, std::string> p;

                sequences.emplace_back(line, line2);


            }
        }

        return sequences;
    }


    Alignment *createAlignment(std::string in_fastafile) {

        auto *alignment = new Alignment;
        std::vector<std::pair<std::string, std::string>> msa;
        msa = fasta_parser::getBuffer(in_fastafile);

        for (int i = 0; i < msa.size(); i++) {
            alignment->addSequence(msa.at(i).first, msa.at(i).second);
        }

        return alignment;

    };


}



int main(int argc, char **argv) {

    /* LOGURU LOGGING ENGINE
     * Possible Loggin Levels
     * FATAL   -3
     * ERROR   -2
     * WARNING -1
     * INFO    +0
     * DEBUG1  +1
     * DEBUG1  +2
     */
    loguru::init(argc, argv);

    LOG_S(DEBUG1) << "test";
    std::string tree_file = argv[1];
    std::string msa_file = argv[2];
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
    nu = LKFunc::compute_nu(tau, lambda, mu);

    // set insertion probability to each node
    tree->set_iota(tau, mu);

    // set survival probability to each node
    tree->set_beta(tau, mu);

    // print newick tree
    LOG_S(DEBUG1) << "[Initial Tree Topology] " << tree->formatNewick();
    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA FROM FILE

    auto *alignment = new Alignment;
    alignment = fasta_parser::createAlignment(msa_file);

    //----------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION

    int is_DNA_AA_Codon;
    int extended_alphabet_size;
    unsigned long num_leaves;
    unsigned long MSA_len;
    double log_col_lk;
    double logLK;
    Eigen::VectorXd pi;
    double p0;


    num_leaves = alignment->align_dataset.size();
    LOG_S(DEBUG1) << "[Sequences in MSA] Leaves: " << num_leaves;

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

    //------------------------------------------------------------------------------------------------------------------
    // COMPUTE LK GIVEN TREE TOPOLOGY AND MSA
    logLK = 0.0;
    // compute log_col_lk
    for (int i = 0; i < MSA_len; i++) {

        // extract MSA column
        std::string s = alignment->extractColumn(i);
        LOG_S(DEBUG1) << "[Extracted column] (" << i << ") = " << s;

        // assign char at the leaves
        tree->set_leaf_state(s);

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        tree->set_ancestral_flag(s);

        // Initialise FV matrices at each node
        tree->clear_fv();

        // Compute column likelihood
        //TODO: Add weight per column
        log_col_lk = LKFunc::compute_col_lk(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);

        LOG_S(DEBUG1) << "[Initial LK] P(c" << i << ") = " << log_col_lk;

        logLK += log_col_lk;
    }

    // compute empty column likelihood
    p0 = LKFunc::compute_log_lk_empty_col(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);
    LOG_S(DEBUG1) << "[Initial LK] p0 = " << p0;

    logLK += LKFunc::phi(MSA_len, nu, p0);
    LOG_S(DEBUG1) << "[Initial LK] LK = " << logLK;
    //------------------------------------------------------------------------------------------------------------------
    // BUILD UNROOTED TREE

    auto real_utree = new Utree;

    std::vector<node *> utree;
    node *utree_pseudo_root;

    UtreeUtils::convertUtree(tree, real_utree);
    LOG_S(DEBUG1) << "[Initial Utree listVNodes] " << real_utree->printTreeNewick(true);

    //------------------------------------------------------------------------------------------------------------------
    // TEST Function findPseudoRoot either fixing the root on a leaf or on the middle node

    std::string strpath;
    std::vector<VirtualNode *> path2root;
    //VirtualNode *startNode = real_utree->listVNodes.at(0)->getNodeLeft()->getNodeLeft();
    VirtualNode *startNode = real_utree->listVNodes.at(0);
    // ------------------------------------
    real_utree->fixPseudoRootOnNextSubtree = false;
    path2root = real_utree->findPseudoRoot(startNode);


    for (auto &node: path2root) {
        strpath += "->" + node->vnode_name;
    }

    LOG_S(DEBUG2) << "[Path back to root] (Root on middle node) " << strpath;
    path2root.clear();
    strpath.clear();

    // ------------------------------------
    real_utree->fixPseudoRootOnNextSubtree = true;
    path2root = real_utree->findPseudoRoot(startNode);


    for (auto &tnode: path2root) {
        strpath += "->" + tnode->vnode_name;
    }
    LOG_S(DEBUG2) << "[Path back to root] (Root on next subtree) " << strpath;
    path2root.clear();
    strpath.clear();

    //------------------------------------------------------------------------------------------------------------------
    // GET ALL NODES WITHIN RADIUS

    int min_radius = 3;
    int max_radius = 6;

    // Print node description with neighbors
    for (auto &vnode:real_utree->listVNodes) {
        LOG_S(DEBUG2) << "[utree neighbours] " << vnode->printNeighbours();

        // Initialise a new rearrangement list
        auto ts_list = TreeRearrangment(vnode, min_radius, max_radius, true);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        ts_list.fillListMoves(false);

        // Print the list
        ts_list.printListMoves();

        LOG_S(DEBUG1) << "[tsh] Strategy " << ts_list.mset_strategy;
        LOG_S(DEBUG2) << "[utree rearrangment] Found " << ts_list.getNumberOfMoves() << " possible moves for node " << vnode->vnode_name;

        for (unsigned long i = 0; i < ts_list.getNumberOfMoves(); i++) {
            bool status;


            try {
                status = ts_list.applyMove(i);
                if (status) {

                    LOG_S(DEBUG2) << "[apply move]\t" << ts_list.getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                                  << " | (" << ts_list.getSourceNode()->vnode_name << "->" << ts_list.getMove(i)->getTargetNode()->vnode_name << ")\t| "
                                  << real_utree->printTreeNewick(true);


                }
            } catch (const std::exception &e) {

                std::cout << "a standard exception was caught, with message '"
                          << e.what() << "'" << std::endl;

            }

            status = ts_list.revertMove(i);
            if (status) {

                LOG_S(DEBUG2) << "[revert move]\t" << ts_list.getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                              << " | (" << ts_list.getMove(i)->getTargetNode()->vnode_name << "->" << ts_list.getSourceNode()->vnode_name << ")\t| "
                              << real_utree->printTreeNewick(true);
            }
        }

    }

    exit(0);
}
