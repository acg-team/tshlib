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
#include <type_traits>

#include "PhyTree.hpp"
#include "TreeRearrangment.hpp"
#include "Alignment.hpp"
#include "Likelihood.hpp"
#include "newick.hpp"
#include "Alphabet.h"




//===================================================================================================================
int main(int argc, char **argv) {


    //PhyTree *t1;
    //PhyTree *t2;
    //std::string tree_file="/home/max/PIP_C++/NNI_SPR/tree_5_leaves_r_bl.nwk";
    std::string tree_file = argv[1];
    PhyTree *tree = nullptr;
    double mu;
    double lambda;
    double tau;
    double nu;

    mu = 0.1;
    lambda = 0.2;


    //----------------------------------------------------------
    // INIT TREE

    // tree filename
    std::ifstream tree_str(tree_file.c_str());

    // read newick file
    tree = newick_parser::parse_newick(&tree_str);

    // set name of internal nodes
    tree->set_missing_node_name("V");

    // compute total tree length
    tau = tree->computeLength();

//    std::cout<<tau<<"\n";

    // compute the normalizing Poisson intensity
    nu = compute_nu(tau, lambda, mu);

    // set insertion probability to each node
    tree->set_iota(tau, mu);

    // set survival probability to each node
    tree->set_beta(tau, mu);

    std::cout << tree->formatNewick() << "\n\n";
    tree->print();
    std::cout << "\n";
    //----------------------------------------------------------
    // LOAD MSA

    //Alignment ;
    auto *alignment = new Alignment;


    std::vector<std::pair<std::string, std::string> > MSA;

/*    std::string seq1_label="A";
    std::string seq1_DNA="ACGT";
    std::string seq2_label="B";
    std::string seq2_DNA="TAGC";
    std::string seq3_label="C";
    std::string seq3_DNA="GCAT";
    std::string seq4_label="D";
    std::string seq4_DNA="CGTA";
    std::string seq5_label="E";
    std::string seq5_DNA="GTCA";
*/
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

    std::cout << alignment->align_dataset.size() << std::endl;
    //----------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION
    int is_DNA_AA_Codon;
    int extended_alphabet_size;

    is_DNA_AA_Codon = 1; // 1:DNA, 2:AA, 3:Codon
    extended_alphabet_size = 5; // DNA alphabet

    // set "pseudo" probability matrix
    tree->tmp_initPr(extended_alphabet_size); //TODO: pass Q from codonPhyML (?)

    unsigned long MSA_len;
    double lk;
    Eigen::VectorXd pi;

    // set Pi, steady state frequencies
    pi = Eigen::VectorXd::Zero(extended_alphabet_size);
    pi[0] = 0.25;
    pi[1] = 0.25;
    pi[2] = 0.25;
    pi[3] = 0.25;
    pi[4] = 0.0;

    // get MSA length
    MSA_len = MSA.at(0).second.size();

    //std::cout<<"MSA_len="<<MSA_len<<"\n";

    double LK = 0;

    // compute lk
    for (int i = 0; i < MSA_len; i++) {

        // extract MSA column
        std::string s = create_col_MSA(MSA, i);
        //std::cout<<"col["<<i<<"]="<<s<<"\n";

        // assign char at the leaves
        tree->set_leaf_state(s);
        //print_leaf_state(tree);

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        tree->set_ancestral_flag(s);
        //print_descCount(tree);
        //print_ancestral_flag(tree);

        tree->clear_fv();

        // compute column likelihood
        //TODO: Add weight per column
        lk = compute_col_lk(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);

        std::cout << "col_lk=" << lk << "\n";

        LK += lk;
    }

    double p0;
    compute_lk_empty_col(*tree, p0, pi, is_DNA_AA_Codon, extended_alphabet_size);
    p0 = log(p0);

    std::cout << "p0=" << p0 << "\n";

    LK += phi(MSA_len, nu, p0);
    //----------------------------------------------------------
    // GET ALL NODES WITHIN RADIUS TODO: implement high-level heuristic to decide which TS to use according to node level

    int radius;
    PhyTree *node;
    std::vector<move_info> nni_spr_stack;

    node = tree->get_left_child();
    radius = 3;

    auto ts_stack = new TreeRearrangment(node, radius, true);
    ts_stack->fillListMoves(false);
    std::cout << "--------------------------------------------------------------------\n";
    std::cout << "size TreeRearrangement list:" << ts_stack->mset_moves.size() << "\n";

    for (unsigned int i = 0; i < ts_stack->mset_moves.size(); i++) {


        std::cout << "list[" << i << "]=(" << ts_stack->mset_sourcenode->getName() << ";"
                  << ts_stack->mset_moves.at(i)->getTargetNode()->getName() << ")\n";
    }

    std::cout << "--------------------------------------------------------------------\n";
    get_list_nodes_within_radius(node, radius, nni_spr_stack);

    std::cout << "size move_info list:" << nni_spr_stack.size() << "\n";

    for (unsigned int i = 0; i < nni_spr_stack.size(); i++) {
        std::cout << "list[" << i << "]=(" << (nni_spr_stack.at(i)).node1->getName() << ";"
                  << (nni_spr_stack.at(i)).node2->getName() << ")\n";
    }


    //----------------------------------------------------------
    // PERFORM SPR MOVES and RECOMPUTE LK

    //move_info m;
    move_info n;
    int max_idx;
    double max_val;
    std::vector<PhyTree *> p;
    bool valid_move;

    max_val = -INFINITY;
    for (unsigned int i = 0; i < nni_spr_stack.size(); i++) {

        // perform SPR move
        std::cout << "\n\nPerform SPR move\n";
        n = nni_spr_stack.at(i);

        //std::cout<<"ID: "<<n.ID<<"\n";
        std::cout << "n.t1=" << n.node1->getName() << " : n.t2=" << n.node2->getName() << "\n";
        valid_move = tree->swap2(n.node1, n.node2);

        if (valid_move) {
            // print newick
            std::cout << "after SPR move\n";
            std::cout << tree->formatNewick() << "\n";


            // get all nodes in the SPR path
            p = get_path_from_nodes(n.node1, n.node2);

            // update all fv values
            // TODO: Refactor as method on node
            update_fv_values(p, extended_alphabet_size);

            //TODO recompute the sum

            // store index of max
            if (n.lk > max_val) {
                max_val = n.lk;
                max_idx = i;
            }


            // rollback SPR move
            std::cout << "Perform SPR move rollback\n";
            tree->swap2(n.node1, n.node2);

            // print newick
            std::cout << "after rollback\n";
            std::cout << tree->formatNewick() << "\n";

            p.clear();
        }
    }

    std::cout << "max_val:" << max_val << " at index: " << max_idx << "\n";

    //	nni_spr_stack.pop_back();
    nni_spr_stack.empty();

    return 0;
}

//----------------------------------------------------------
//	for(int i=0;i<10;i++){
//		Eigen::VectorXd fv;
//		fv=Eigen::VectorXd::Zero(5);
//		fv[0]=i+1;
//		tree->append_MSA_fv(fv);
//	}
//
//	for(int i=0;i<10;i++){
//		Eigen::VectorXd fv=tree->get_MSA_fv(i);
//		std::cout<<"fv:\n"<<fv<<"\n";
//	}
//----------------------------------------------------------
//	std::vector<PhyTree *> p;
//
////	t1=tree->get_right_child();
////	t2=tree->get_left_child()->get_left_child()->get_left_child();
//	t1=tree->get_left_child()->get_left_child()->get_left_child()->get_right_child();
//	t2=tree->get_left_child()->get_left_child()->get_left_child()->get_left_child();
//
//	p=get_path_from_nodes(t1,t2);
//
//	for(unsigned int i=0;i<p.size();i++){
//		std::cout<<"p["<<i<<"]="<<(p.at(i))->getName()<<"\n";
//	}
//	update_fv_values(p); //TODO fill all the MSA_fv vectors
//----------------------------------------------------------
//	PhyTree *t1;
//
//	t1=tree->get_right_child()->copy();
//	t1->null_parent();
//	double bl=t1->getBranchLength();
//
//	std::cout<<t1->getName()<<"\n";
//	std::cout<<tree->n_children()<<"\n";
//
//	tree->deleteChild(1);
//
//	std::cout<<tree->n_children()<<"\n";
//
//	tree->addChild(t1,bl,0);
//
//	tree->print();
//	std::cout<<"\n";
//
//	std::cout<<tree->formatNewick()<<"\n\n";
//----------------------------------------------------------
//	PhyTree *t1;
//	PhyTree *t2;
//
//	t1=tree->get_right_child();
//	t2=tree->get_left_child()->get_left_child()->get_left_child();
//
////	tree->swap(tree,1,tree->get_left_child(),1);
//	tree->swap2(t1,t2);
//
//	std::cout<<tree->formatNewick()<<"\n\n";
//----------------------------------------------------------
//	std::stack<move_info> nni_spr_stack; // STL Stack object
//
//	for(int i=0;i<10;i++){
//		move_info m;
//		m.ID=i;
//		nni_spr_stack.push(m);
//	}
//
//	for(int i=0;i<10;i++){
//		move_info mm = nni_spr_stack.top();
//		std::cout<<"ID: "<<mm.ID<<"\n";
//		nni_spr_stack.pop();
//	}
//----------------------------------------------------------





