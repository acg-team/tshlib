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
#include <Optimization_Brent.hpp>

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

                sequences.emplace_back(line.substr(1, line.size()), line2);


            }
        }

        return sequences;
    }


    Alignment *createAlignment(std::string in_fastafile, AlignmentAlphabet alphabet) {

        auto *alignment = new Alignment;
        std::vector<std::pair<std::string, std::string>> msa;
        msa = fasta_parser::getBuffer(in_fastafile);

        for (int i = 0; i < msa.size(); i++) {
            alignment->addSequence(msa.at(i).first, msa.at(i).second);
        }

        alignment->getAlignmentSize();

        alignment->align_num_characters.resize((unsigned long) alignment->align_length);
        alignment->alphabet = alphabet;

        switch(alphabet){
            case AlignmentAlphabet::dna:
                alignment->align_alphabetsize = 5;
                break;
            case AlignmentAlphabet::aa:
                alignment->align_alphabetsize = 21;
                break;
            case AlignmentAlphabet::codon:
                alignment->align_alphabetsize = 65;
                break;
            default:
                LOG(FATAL) << "Alphabet not implemented";

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
    double logLK;
    Eigen::VectorXd pi;
    Eigen::MatrixXd Q;

    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA FROM FILE
    Alignment *alignment = fasta_parser::createAlignment(msa_file, AlignmentAlphabet::dna);

    num_leaves = alignment->align_dataset.size();
    LOG(INFO) << "[Sequences in MSA] Leaves: " << num_leaves << std::endl;
    alignment->countNumberCharactersinColumn();
    // 1:DNA, 2:AA, 3:Codon


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
    LOG(INFO) << "[Initial Utree Topology] " << utree->printTreeNewick(true) << std::endl;
    utree->prepareSetADesCountOnNodes((int) alignment->getAlignmentSize(), extended_alphabet_size);
    UtreeUtils::associateNode2Alignment(alignment, utree);

    // TODO: refactoring of prepareSetADescCountonNodes + initialiseLikelihoodComponents (one function to rule them all)
    // Add the root
    utree->addVirtualRootNode();
    utree->rootnode->initialiseLikelihoodComponents((int) alignment->getAlignmentSize(), extended_alphabet_size);

    //---------------------------------------------------------
    // compute total tree length

    //tau = tree->computeLength();
    //VLOG(2) << "[PhyTree length] " << tau;
    tau = utree->computeTotalTreeLength();
    VLOG(2) << "[Utree length] " << tau;
    //tau = root->computeTotalTreeLength();
    tau = utree->rootnode->computeTotalTreeLength();
    VLOG(2) << "[Tree length (from root)] " << tau;


    // set Pi, steady state frequencies
    pi = Eigen::VectorXd::Zero(extended_alphabet_size);
    pi << 0.25, 0.25, 0.25, 0.25, 0.25;

    // Fill Q matrix as for JC69
    Q = Eigen::MatrixXd::Zero(extended_alphabet_size,extended_alphabet_size);
    for(int r = 0; r<Q.rows()-1; r++){
        for(int c=0; c<Q.cols()-1; c++ ){

            if(r==c){
                Q(r,c) =  -3.0/4.0-mu;
            }else{
                Q(r,c) =  1.0/4.0;
            }

        }
        Q(r,extended_alphabet_size-1) = mu;
    }



    auto likelihood = new Likelihood();

    likelihood->Init(utree, pi, Q, mu, lambda);

    // compute the normalizing Poisson intensity
    //nu = likelihood->setNu(tau, lambda, mu);

/*
    auto pip_model = new PIP;
    pip_model->lambda = lambda;
    pip_model->mu = mu;
*/


    //root->_traverseVirtualNodeTree();

    // print newick tree
    //LOG(INFO) << "[Initial PhyTree Topology] " << tree->formatNewick() << std::endl;


    //----------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION

    // DNA alphabet
    // TODO: extended_alphabet_size to be deleted
    extended_alphabet_size = 5;

    // COMPUTE LK GIVEN TREE TOPOLOGY AND MSA
    logLK = 0.0;

    bool isReferenceRun = true;


    std::vector<VirtualNode *> allnodes_postorder;
    likelihood->fillNodeListComplete_bottomUp(allnodes_postorder, utree->rootnode);

    // set survival probability to each node
    likelihood->setAllIotas(allnodes_postorder);
    likelihood->setAllBetas(allnodes_postorder);

    likelihood->setInsertionHistories(allnodes_postorder,*alignment);
    // set "pseudo" probability matrix
    likelihood->setPr(utree, extended_alphabet_size);


    //TODO: Add weight per column
    likelihood->computeFV(allnodes_postorder, *alignment);
    //likelihood->unloadParametersOperative();

    // Initialise likelihood components on the tree
    //for (int i = 0; i < alignment->getAlignmentSize(); i++) {

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        //utree->rootnode->setAncestralFlag(*alignment, i, isReferenceRun);

        //root->_traverseVirtualNodeTree();

        // Initialise FV matrices at each node
        //tree->clear_fv();
        //utree->clearFv();

        // Compute column likelihood
        //likelihood->computeLikelihoodComponents(utree->rootnode, likelihood->pi, is_DNA_AA_Codon, extended_alphabet_size, i, *alignment);

    //}

    // compute empty column likelihood
    //likelihood->computeLikelihoodComponents_EmptyColumn(utree->rootnode, likelihood->pi, is_DNA_AA_Codon, extended_alphabet_size);

    //p0 = LKFunc::compute_log_lk_empty_col(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);
    //VLOG(2) << "[Initial LK] p0 = " << p0 << std::endl;

    //logLK += likelihood->phi((int) alignment->align_length, nu, p0);
    //VLOG(1) << "[Initial LK] LK = " << logLK << std::endl;

    //----------------------------------------------
    //likelihood->loadParametersOperative();

    logLK = likelihood->computePartialLK(allnodes_postorder, *alignment, likelihood->pi);
    likelihood->saveLikelihoodComponents();

    VLOG(2) << "[Initial LK] Full Tree from partial lk routines " << logLK;

    //------------------------------------------------------------------------------------------------------------------
    // Remove the root
    utree->removeVirtualRootNode();

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

    bool computeMoveLikelihood = true;
    std::vector<VirtualNode *> list_vnode_to_root;

    // Print node description with neighbors
    for (auto &vnode:utree->listVNodes) {
        VLOG(2) << "[utree neighbours] " << vnode->printNeighbours() << std::endl;

        // Initialise a new rearrangement list
        auto rearrangmentList = new TreeRearrangment;

        rearrangmentList->initTreeRearrangment(utree, min_radius, max_radius, true, vnode);

        // Get all the target nodes with distance == radius from the source node
        // excluding the starting node.
        rearrangmentList->defineMoves(false);

        // Print the list of moves for the current P node (source node)
        //rearrangmentList.printMoves();

        VLOG(1) << "[tsh] Strategy " << rearrangmentList->mset_strategy << std::endl;
        VLOG(1) << "[utree rearrangment] Found " << rearrangmentList->getNumberOfMoves() << " possible moves for node " << vnode->vnode_name << std::endl;

        std::string start_col_line, end_col_line;

        // For each potential move computed before, apply it to the tree topology, print the resulting newick tree, and revert it.
        for (unsigned long i = 0; i < rearrangmentList->getNumberOfMoves(); i++) {
            bool status;
            logLK = 0;

            // Prepare the list of nodes involved in the move
            // TODO: This list belongs specifically to the tree-rearrangement definition -> move to TreeRearrangment class
            list_vnode_to_root.clear();
            list_vnode_to_root = utree->computePathBetweenNodes(rearrangmentList->getSourceNode(), rearrangmentList->getMove(i)->getTargetNode());
            list_vnode_to_root.push_back(utree->rootnode);
            //std::reverse(list_vnode_to_root.begin(), list_vnode_to_root.end());

            // Apply the move
            status = rearrangmentList->applyMove(i);

            // Print root reachability from every node (includes rotations)
            //utree->_testReachingPseudoRoot();

            // Print tree on file
            //utree->saveTreeOnFile("../data/test.txt");
            bool isLKImproved = false;

            if (computeMoveLikelihood) {

                //utree->printAllNodesNeighbors();
                //testSetAinRootPath(MSA_len, alignment, utree, list_vnode_to_root)

                // Add the root
                utree->addVirtualRootNode();

                // Compute the full likelihood from the list of nodes involved in the rearrangment
                likelihood->recombineAllFv(list_vnode_to_root);
                likelihood->setInsertionHistories(list_vnode_to_root,*alignment);
                logLK = likelihood->computePartialLK(list_vnode_to_root, *alignment, pi);


                //VLOG(2) << "[Tree LK] Before Brent: " << logLK;
                double max_lenght = pi[1] * 1.1;
                double min_lenght = pi[1] * 0.9;

                //isLKImproved = Generic_Brent_Lk(&likelihood->pi[1], min_lenght, max_lenght, SMALL, BRENT_ITMAX,
                //                                LKFunc::LKcore , *likelihood, list_vnode_to_root, *alignment, logLK);

                rearrangmentList->getMove(i)->move_lk = logLK;
                //VLOG(2) << "[Tree LK] Ater Brent: " << logLK;

                utree->removeVirtualRootNode();

            }

            if(rearrangmentList->getMove(i)->move_lk>0){
                start_col_line = "\033[1;34m";
                end_col_line = "\033[0m";
            }else{
                start_col_line = "";
                end_col_line = "";

            }

            VLOG(2) << "[apply  move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                    << " | (" << isLKImproved <<") " << start_col_line<< rearrangmentList->getMove(i)->move_lk<<end_col_line << "\t"
                    << " | (" << rearrangmentList->getSourceNode()->vnode_name << "->" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << ")"
                    << "\t[" << rearrangmentList->getMove(i)->move_radius << "] | "
                    << utree->printTreeNewick(true) << std::endl;



            // Revert the move, and return to the original tree
            status = rearrangmentList->revertMove(i);

            // Print tree on file
            //utree->saveTreeOnFile("../data/test.txt");

            if (computeMoveLikelihood) {


                // Add the root
                utree->addVirtualRootNode();
                //if (!rearrangmentList->getMove(i)->move_applied) {
                //likelihood->revertAllFv(list_vnode_to_root); // clear not necessary
                // Compute the full likelihood from the list of nodes involved in the rearrangment
                likelihood->recombineAllFv(list_vnode_to_root);
                logLK = likelihood->computePartialLK(list_vnode_to_root, *alignment, pi);
                rearrangmentList->getMove(i)->move_lk = logLK;

                utree->removeVirtualRootNode();
                //} else {
                    //UtreeUtils::keepAllFv(list_vnode_to_root);
                //}

            }
            if(rearrangmentList->getMove(i)->move_lk>0){
                start_col_line = "\033[1;34m";
                end_col_line = "\033[0m";
            }else{
                start_col_line = "";
                end_col_line = "";

            }

            VLOG(2) << "[revert move]\t" << rearrangmentList->getMove(i)->move_class << "." << std::setfill('0') << std::setw(3) << i
                    << " | (" << isLKImproved <<") " << start_col_line<< rearrangmentList->getMove(i)->move_lk<<end_col_line << "\t"
                    << " | (" << rearrangmentList->getMove(i)->getTargetNode()->vnode_name << "->" << rearrangmentList->getSourceNode()->vnode_name << ")"
                    << "\t[" << rearrangmentList->getMove(i)->move_radius << "] | "
                    << utree->printTreeNewick(true) << std::endl;

            // Print root reachability from every node (includes rotations)
            //utree->_testReachingPseudoRoot();

            // Count moves performed
            total_exec_moves += rearrangmentList->getNumberOfMoves() * 2;
        }

        // Clean memory
        delete rearrangmentList;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    VLOG(0) << "Moves applied and reverted: " << total_exec_moves << std::endl;
    VLOG(0) << "Elapsed time: " << duration << " microseconds" << std::endl;
    VLOG(0) << "*** " << (double) duration / total_exec_moves << " microseconds/move *** " << std::endl;

    //treesearchheuristics::testTSH(utree, TreeSearchHeuristics::classic_Mixed);


    exit(0);
}

void testSetAinRootPath(unsigned long MSA_len, Alignment *alignment, Utree *utree, std::vector<VirtualNode *> &list_vnode_to_root) {
    for (int n = 0; n < utree->listVNodes.size(); n++) {
        VirtualNode *tempnode = utree->listVNodes.at(n);

        for (int msa_col = 0; msa_col < MSA_len; msa_col++) {

            std::__1::string s = alignment->extractColumn(msa_col);
            utree->setLeafState(s);

            utree->findPseudoRoot(tempnode, false).back()->setAncestralFlag(*alignment, msa_col, false);
            utree->findPseudoRoot(tempnode, true).back()->setAncestralFlag(*alignment, msa_col, false);
            //utree->findPseudoRoot()
            // root->setAncestralFlag(s, msa_col, false);

            bool temp = tempnode->vnode_setA_operative.at(msa_col);
            bool ref = tempnode->vnode_setA_backup.at(msa_col);

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
