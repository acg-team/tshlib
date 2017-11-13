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
#include <TreeRearrangment.hpp>
#include <Alignment.hpp>
#include <Likelihood.hpp>
#include <newick.hpp>


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
     * Possible logging Levels
     * FATAL   -3
     * ERROR   -2
     * WARNING -1
     * INFO    +0
     * DEBUG1  +1
     * DEBUG1  +2
     */
    loguru::init(argc, argv);
    LOG_S(DEBUG1) << "test";
    //------------------------------------------------------------------------------------------------------------------
    std::string tree_file = argv[1];
    std::string msa_file = argv[2];
    //------------------------------------------------------------------------------------------------------------------
    PhyTree *tree = nullptr;
    double mu;
    double lambda;
    double tau;
    double nu;
    //------------------------------------------------------------------------------------------------------------------
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

    Alignment *alignment = fasta_parser::createAlignment(msa_file);

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

    auto utree = new Utree;

    //std::vector<node *> utree;
    //node *utree_pseudo_root;

    UtreeUtils::convertUtree(tree, utree);
    std::cout << "[Initial Utree listVNodes] " << utree->printTreeNewick(true) << std::endl;

    // Save tree to file
    utree->saveTreeOnFile("../data/test.txt");

    //------------------------------------------------------------------------------------------------------------------
    // TEST Function findPseudoRoot either fixing the root on a subtree or on the middle node

    std::string strpath;
    std::vector<VirtualNode *> path2root;
    //VirtualNode *startNode = utree->listVNodes.at(0)->getNodeLeft()->getNodeLeft();
    VirtualNode *startNode = utree->listVNodes.at(0);
    // ------------------------------------
    bool fixPseudoRootOnNextSubtree = false;
    path2root = utree->findPseudoRoot(startNode, fixPseudoRootOnNextSubtree);

    for (auto &node: path2root) {
        strpath += "->" + node->vnode_name;
    }

    LOG_S(DEBUG2) << "[Path back to root] (Root on middle node) " << strpath;
    path2root.clear();
    strpath.clear();

    // ------------------------------------
    fixPseudoRootOnNextSubtree = true;
    path2root = utree->findPseudoRoot(startNode, fixPseudoRootOnNextSubtree);


    for (auto &tnode: path2root) {
        strpath += "->" + tnode->vnode_name;
    }
    LOG_S(DEBUG2) << "[Path back to root] (Root on next subtree) " << strpath;
    path2root.clear();
    strpath.clear();

    //------------------------------------------------------------------------------------------------------------------
    // Get all the nodes between the radius boundaries

    int min_radius = 3;
    int max_radius = 6;

    // Print node description with neighbors
    for (auto &vnode:utree->listVNodes) {
        LOG_S(DEBUG2) << "[utree neighbours] " << vnode->printNeighbours();

        // Initialise a new rearrangement list
        auto rearrangmentList = TreeRearrangment(vnode, min_radius, max_radius, true);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        rearrangmentList.defineMoves(false);

        // Print the list of moves for the current P node (source node)
        rearrangmentList.printMoves();

        LOG_S(DEBUG1) << "[tsh] Strategy " << rearrangmentList.mset_strategy;
        LOG_S(DEBUG2) << "[utree rearrangment] Found " << rearrangmentList.getNumberOfMoves() << " possible moves for node " << vnode->vnode_name;

        /*
        // For each potential move computed before, apply it to the tree topology, print the resulting newick tree, and revert it.
        for (unsigned long i = 0; i < rearrangmentList.getNumberOfMoves(); i++) {
            bool status;

            // Apply the move
            status = rearrangmentList.applyMove(i);
            utree->saveTreeOnFile("../data/test.txt");

            if (status) {
                LOG_S(DEBUG2) << "[apply move]\t" << rearrangmentList.getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                              << " | (" << rearrangmentList.getSourceNode()->vnode_name << "->" << rearrangmentList.getMove(i)->getTargetNode()->vnode_name << ")\t| "
                              << utree->printTreeNewick(true);
                utree->_testReachingPseudoRoot();
            }


            // Revert the move, and return to the original tree
            status = rearrangmentList.revertMove(i);

            if (status) {
                LOG_S(DEBUG2) << "[revert move]\t" << rearrangmentList.getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                              << " | (" << rearrangmentList.getMove(i)->getTargetNode()->vnode_name << "->" << rearrangmentList.getSourceNode()->vnode_name << ")\t| "
                              << utree->printTreeNewick(true);
                utree->_testReachingPseudoRoot();
            }


        }
        */
    }

    exit(0);
}
