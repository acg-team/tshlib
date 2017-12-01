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
#include <glog/logging.h>
#include <gflags/gflags.h>
//#define LOGURU_IMPLEMENTATION 1
//#define LOGURU_WITH_STREAMS 1

//#include <loguru.hpp>
#include <Utree.hpp>
#include <TreeRearrangment.hpp>
#include <Likelihood.hpp>
#include <newick.hpp>

void testSetAinRootPath(unsigned long MSA_len, Alignment *alignment, Utree *utree, std::vector<VirtualNode *> &list_vnode_to_root);

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

                //std::cout << "Found header: " << line;

                // parse next line since it must be the sequence
                std::getline(in_fasta, line2);
                //std::pair<std::string, std::string> p;

                sequences.emplace_back(line.substr(1,line.size()), line2);


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
    // Initialize Google's logging library.
    FLAGS_alsologtostderr = true;
    google::InitGoogleLogging("THSLIB Main File");
    LOG(INFO) << "TSHLIB Initialising" << std::endl;
    //------------------------------------------------------------------------------------------------------------------
    std::string tree_file = argv[1];
    std::string msa_file = argv[2];

    //------------------------------------------------------------------------------------------------------------------
    PhyTree *tree = nullptr;
    double mu;
    double lambda;
    double tau;
    double nu;
    int is_DNA_AA_Codon;
    int extended_alphabet_size;
    unsigned long num_leaves;
    unsigned long MSA_len;
    double log_col_lk;
    double logLK;
    Eigen::VectorXd pi;
    double p0;

    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA FROM FILE
    Alignment *alignment = fasta_parser::createAlignment(msa_file);

    num_leaves = alignment->align_dataset.size();
    LOG(INFO) << "[Sequences in MSA] Leaves: " << num_leaves << std::endl;

    // 1:DNA, 2:AA, 3:Codon
    is_DNA_AA_Codon = 1;

    // DNA alphabet
    extended_alphabet_size = 5;

    mu = 0.1;
    lambda = 0.2;

    //------------------------------------------------------------------------------------------------------------------
    // INIT ROOTED TREE

    // tree filename
    std::ifstream tree_str(tree_file.c_str());

    // read newick file
    tree = newick_parser::parse_newick(&tree_str);

    // set name of internal nodes
    tree->set_missing_node_name("V"); //TODO: check whether V is unique

    //tree->printOnlyName();

    //------------------------------------------------------------------------------------------------------------------
    // BUILD UNROOTED TREE
    auto utree = new Utree;
    UtreeUtils::convertUtree(tree, utree);
    LOG(INFO) << "[Initial utree] " << utree->printTreeNewick(true) << std::endl;
    utree->prepareSetADesCountOnNodes((int) alignment->getAlignmentSize());
    UtreeUtils::associateNode2Alignment(alignment,utree);

    // Add the root
    auto *root = new VirtualNode;
    root->prepareSetA_DescCount((int) alignment->getAlignmentSize());
    //    double T = tau + 1 / mu;
    //    double iota = (1/mu)/T;
    //    double beta = 1.0;

    root->setNodeParent(nullptr);
    root->setNodeName("Root");
    root->clearChildren();

    VirtualNode *pseudo_root1;
    VirtualNode *pseudo_root2;
    pseudo_root1=UtreeUtils::getPseudoRoot(utree->startVNodes.at(0));
    pseudo_root2=pseudo_root1->getNodeUp();
    root->setChild(pseudo_root1);
    root->setChild(pseudo_root2);


    //---------------------------------------------------------
    // compute total tree length

    tau = tree->computeLength();
    VLOG(2) << "[PhyTree length] " << tau;
    tau=utree->computeTotalTreeLength();
    VLOG(2) << "[Utree length] " << tau;
    tau=root->computeTotalTreeLength();
    VLOG(2) << "[Tree length (from root)] " << tau;

    // compute the normalizing Poisson intensity
    nu = LKFunc::compute_nu(tau, lambda, mu);

    // set insertion probability to each node
    tree->set_iota(tau, mu);
    utree->setIota(tau,mu);
    root->setAllIotas(tau,mu);

    // set survival probability to each node
    tree->set_beta(tau, mu);
    utree->setBeta(tau,mu);
    root->setAllBetas(mu);

    //tree->print_local_var();
    //std::cout<<std::endl;
    //utree->_printUtree();

    // set "pseudo" probability matrix
    tree->tmp_initPr(extended_alphabet_size); //TODO: pass Q from codonPhyML (?)
    utree->setPr(extended_alphabet_size);

    //root->_traverseVirtualNodeTree();

    // print newick tree
    LOG(INFO) << "[Initial Tree Topology] " << tree->formatNewick() << std::endl;


    //----------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION
    num_leaves = alignment->align_dataset.size();
    LOG(INFO) << "[Sequences in MSA] Leaves: " << num_leaves << std::endl;

    // 1:DNA, 2:AA, 3:Codon
    is_DNA_AA_Codon = 1;

    // DNA alphabet
    extended_alphabet_size = 5;

    // set "pseudo" probability matrix
    tree->tmp_initPr(extended_alphabet_size); //TODO: pass Q from codonPhyML (?)

    // set Pi, steady state frequencies
    pi = Eigen::VectorXd::Zero(extended_alphabet_size);
    pi << 0.25,0.25,0.25,0.25,0.25;
//    pi[0] = 0.25;
//    pi[1] = 0.25;
//    pi[2] = 0.25;
//    pi[3] = 0.25;
//    pi[4] = 0.0;

    // get MSA length
    MSA_len = static_cast<int>(alignment->getAlignmentSize());

    // COMPUTE LK GIVEN TREE TOPOLOGY AND MSA
    logLK = 0.0;

    double log_col_lk0;
    bool isReferenceRun = true;
    // compute log_col_lk
    for (int i = 0; i < MSA_len; i++) {

        // extract MSA column
        std::string s = alignment->extractColumn(i);
        VLOG(2) << "[Extracted column] (" << i << ") = " << s << std::endl;

        // assign char at the leaves
        tree->set_leaf_state(s);
        utree->setLeafState(s);

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        tree->set_ancestral_flag(s);
        root->setAncestralFlag(s, i, isReferenceRun);

        //root->_traverseVirtualNodeTree();

        // Initialise FV matrices at each node
        tree->clear_fv();
        utree->clearFv();

        // Compute column likelihood
        //TODO: Add weight per column
        log_col_lk0 = LKFunc::compute_col_lk(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);
        log_col_lk = LKFunc::compute_col_lk(root, pi, is_DNA_AA_Codon, extended_alphabet_size, i);

        VLOG(2) << "[Initial LK] P(c" << i << ") = " << log_col_lk;

        VLOG(2) << "[diff lk] " << abs(log_col_lk0-log_col_lk);

        logLK += log_col_lk;

    }

    // compute empty column likelihood
    p0 = LKFunc::compute_log_lk_empty_col(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);
    VLOG(2) << "[Initial LK] p0 = " << p0 << std::endl;

    p0 = LKFunc::compute_log_lk_empty_col(root, pi, is_DNA_AA_Codon, extended_alphabet_size);
    VLOG(2) << "[Initial LK] p0 = " << p0 << std::endl;

    logLK += LKFunc::phi(MSA_len, nu, p0);
    VLOG(1) << "[Initial LK] LK = " << logLK << std::endl;

    //------------------------------------------------------------------------------------------------------------------
    // Remove the root
    VirtualNode *vnL;
    VirtualNode *vnR;
    vnL=root->getNodeLeft();
    vnR=root->getNodeRight();
    vnL->setNodeParent(vnR);
    vnR->setNodeParent(vnL);
    //delete root;
    //utree->_printUtree();
    //----------------------------------------------

    // Save tree to file
    //utree->saveTreeOnFile("../data/test.txt");
    utree->printAllNodesNeighbors();

    //------------------------------------------------------------------------------------------------------------------
    // TEST
    // findPseudoRoot either fixing the root on a subtree or on the middle node

    std::string strpath;
    std::vector<VirtualNode *> path2root;

    VirtualNode *startNode = utree->listVNodes.at(0);
    // ------------------------------------
    // Retrieve the root (on the middle root node) starting from a "random" point on the tree
    bool fixPseudoRootOnNextSubtree = false;
    path2root = utree->findPseudoRoot(startNode, fixPseudoRootOnNextSubtree);

    for (auto &node: path2root) strpath += "->" + node->vnode_name;

    VLOG(2) << "[Path back to root] (Root on middle node) " << strpath << std::endl;
    path2root.clear();
    strpath.clear();

    // ------------------------------------
    // Retrieve the root (on the subtree immediately after the pseudo root ndoe) starting from a "random" point on the tree
    fixPseudoRootOnNextSubtree = true;
    path2root = utree->findPseudoRoot(startNode, fixPseudoRootOnNextSubtree);

    for (auto &tnode: path2root) strpath += "->" + tnode->vnode_name;

    VLOG(2) << "[Path back to root] (Root on next subtree) " << strpath << std::endl;
    path2root.clear();
    strpath.clear();

    //------------------------------------------------------------------------------------------------------------------
    // DEFINE, APPLY & REVERT TREE REARRANGEMENTS
    // Get all the nodes between the radius boundaries and for each of them build the move list

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    unsigned long total_exec_moves = 0;

    int min_radius = 3;  // Minimum radius for an NNI move is 3 nodes
    int max_radius = utree->getMaxNodeDistance(); // Hard coded max value for a small tree (this ensures the complete q-node search)

    std::vector<VirtualNode *> list_vnode_to_root;

    bool is_the_best_move = false;
    int ID_best_move;

    // Print node description with neighbors
    for (auto &vnode:utree->listVNodes) {
        VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

        // Initialise a new rearrangement list
        auto rearrangmentList = new TreeRearrangment;

        rearrangmentList->initTreeRearrangment(vnode, min_radius, max_radius, true);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        rearrangmentList->defineMoves(false);

        // Print the list of moves for the current P node (source node)
        //rearrangmentList.printMoves();

        VLOG(1) << "[tsh] Strategy " << rearrangmentList->mset_strategy << std::endl;
        VLOG(1) << "[utree rearrangment] Found " << rearrangmentList->getNumberOfMoves() << " possible moves for node " << vnode->vnode_name << std::endl;


        // For each potential move computed before, apply it to the tree topology, print the resulting newick tree, and revert it.
        for (unsigned long i = 0; i < rearrangmentList->getNumberOfMoves(); i++) {
            bool status;

//            //Add the root
//            root->setNodeName("Root");
//            root->clearChildren();
//            pseudo_root1 = UtreeUtils::getPseudoRoot(utree->startVNodes.at(0));
//            pseudo_root2 = pseudo_root1->getNodeUp();
//            root->setChild(pseudo_root1);
//            root->setChild(pseudo_root2);

            VirtualNode *source;
            VirtualNode *target;
            list_vnode_to_root.clear();

            source = rearrangmentList->getSourceNode();
            target = rearrangmentList->getMove(i)->getTargetNode();

            std::vector<VirtualNode *> path2root_1 = utree->findPseudoRoot(source, false);
            std::vector<VirtualNode *> path2root_2 = utree->findPseudoRoot(target, false);

            list_vnode_to_root = UtreeUtils::get_unique(path2root_1,path2root_2);


            // Apply the move
            status = rearrangmentList->applyMove(i);

            if (status) {
                VLOG(2) << "[apply  move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                        << " | (" << rearrangmentList->getSourceNode()->vnode_name << "->" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << ")"
                        << "\t[" << rearrangmentList->getMove(i)->move_radius << "] | "
                        << utree->printTreeNewick(true) << std::endl;
                //utree->_testReachingPseudoRoot();
            }

            //utree->saveTreeOnFile("../data/test.txt");

            if(status) {
                //utree->printAllNodesNeighbors();


                testSetAinRootPath(MSA_len, alignment, utree, list_vnode_to_root);

                utree->addRootNode();

                //UtreeUtils::recombineAllFv(list_vnode_to_root);

                //TODO: here goes the smart lk recomputation
                ID_best_move = i; // index (ID) of the best move
                is_the_best_move=true;

                utree->removeRootNode();


            }



            // Revert the move, and return to the original tree
            status = rearrangmentList->revertMove(i);
            //utree->saveTreeOnFile("../data/test.txt");

            if(status){

                if(!is_the_best_move){
                    //UtreeUtils::revertAllFv(list_vnode_to_root); // clear not necessary
                }else{
                    //UtreeUtils::keepAllFv(list_vnode_to_root);
                }

            }

            if (status) {
                VLOG(2) << "[revert move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                              << " | (" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << "->" << rearrangmentList->getSourceNode()->vnode_name << ")"
                              << "\t[" << rearrangmentList->getMove(i)->move_radius << "] | "
                              << utree->printTreeNewick(true) << std::endl;
                //utree->_testReachingPseudoRoot();
            }
            total_exec_moves += rearrangmentList->getNumberOfMoves()*2;
        }
        delete rearrangmentList;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    VLOG(0) << "Moves applied and reverted: " << total_exec_moves << std::endl;
    VLOG(0) << "Elapsed time: " << duration << " microseconds" << std::endl;
    VLOG(0) << "*** " << (double) duration/total_exec_moves << " microseconds/move *** " << std::endl;

    //treesearchheuristics::testTSH(utree, TreeSearchHeuristics::classic_Mixed);


    exit(0);
}

void testSetAinRootPath(unsigned long MSA_len, Alignment *alignment, Utree *utree, std::vector<VirtualNode *> &list_vnode_to_root) {
    for (int n = 0; n < utree->listVNodes.size(); n++) {
        VirtualNode *tempnode = utree->listVNodes.at(n);

        for (int msa_col = 0; msa_col < MSA_len; msa_col++) {

            std::__1::string s = alignment->extractColumn(msa_col);
            utree->setLeafState(s);

            utree->findPseudoRoot(tempnode, false).back()->setAncestralFlag(s, msa_col, false);
            utree->findPseudoRoot(tempnode, true).back()->setAncestralFlag(s, msa_col, false);
            //utree->findPseudoRoot()
            // root->setAncestralFlag(s, msa_col, false);

            bool temp = tempnode->vnode_setA_temp.at(msa_col);
            bool ref = tempnode->vnode_setA.at(msa_col);

            if (ref == false && temp == true) {
                if (find(list_vnode_to_root.begin(), list_vnode_to_root.end(), tempnode) != list_vnode_to_root.end()) {
                } else {

                    //VLOG(1) << "[Not found] Root <" << root->getNodeLeft()->vnode_name << ";" << root->getNodeRight()->vnode_name << ">" ;
                    LOG(FATAL) << "[Not found] This node> " << tempnode->vnode_name << "@" << msa_col + 1 << "temp: " << temp << " ref:  " << ref;
                }
            }
        }
    }
}
