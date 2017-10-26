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

#define LOGURU_IMPLEMENTATION 1

#include <loguru.hpp>

#include "PhyTree.hpp"
#include "Utree.hpp"
#include "TreeRearrangment.hpp"
#include "Alignment.hpp"
#include "Likelihood.hpp"
#include "newick.hpp"
#include "Alphabet.h"

int main(int argc, char **argv) {

    loguru::init(argc, argv);
    //LOG_SCOPE_F(INFO, "Will indent all log messages within this scope.");


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
    LOG_F(INFO, "Parsing tree file %s ", argv[1]);

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
    std::cout << tree->formatNewick() << "\n\n";

    //------------------------------------------------------------------------------------------------------------------
    // LOAD MSA
    auto *alignment = new Alignment;

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
    unsigned long num_leaves;
    unsigned long MSA_len;
    double log_col_lk;
    double logLK;
    Eigen::VectorXd pi;
    double p0;

    num_leaves = alignment->align_dataset.size();
    //num_leaves = MSA.size();

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
    //MSA_len = MSA.at(0).second.size();
    MSA_len = static_cast<unsigned long>(alignment->getAlignmentSize());

    logLK = 0.0;
    // compute log_col_lk
    for (int i = 0; i < MSA_len; i++) {

        // extract MSA column

        std::string s = alignment->extractColumn(i);
        //std::string s = create_col_MSA(MSA, i);

        // assign char at the leaves
        tree->set_leaf_state(s);

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        tree->set_ancestral_flag(s);

        tree->clear_fv();

        // compute column likelihood
        //TODO: Add weight per column
        log_col_lk = compute_col_lk(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);

        std::cout << "log_col_lk=" << log_col_lk << "\n";

        logLK += log_col_lk;
    }

    // compute empty column likelihood
    p0 = compute_log_lk_empty_col(*tree, pi, is_DNA_AA_Codon, extended_alphabet_size);

    std::cout << "p0=" << p0 << std::endl;

    logLK += phi(MSA_len, nu, p0);
    //------------------------------------------------------------------------------------------------------------------
    // BUILD UNROOTED TREE

    std::vector<node *> utree;
    node *utree_pseudo_root;

    utree = tree->create_unrooted_tree(tree, num_leaves);

    utree_pseudo_root = utree.at(0);

    /*
    for(unsigned int i=0;i<utree.size();i++){
        node *n=utree.at(i);
        std::cout<<"\n"<<i<<"\n";
        std::cout<<n->data->getName()<<"\n";
        if(n->next!=NULL){
            std::cout<<"next: "<<n->next->data->getName()<<"\n";
        }
        if(n->back!=NULL){
            std::cout<<"back: "<<n->back->data->getName()<<"\n";
        }
    }
    */
    //------------------------------------------------------------------------------------------------------------------
    // GET ALL NODES WITHIN RADIUS
    // TODO: implement high-level heuristic to decide which TS to use according to node level

    int radius;

    radius = 3;

    //----------------------------------------------------
//    PhyTree *t_node;
//    std::vector<move_info> nni_spr_stack_left;
//    std::vector<move_info> nni_spr_stack_right;
//    std::vector<move_info> nni_spr_stack_up;
//    t_node = tree->get_left_child()->get_left_child();
    //get_list_nodes_within_radius(t_node, radius, nni_spr_stack_left,nni_spr_stack_right,nni_spr_stack_up);
    //----------------------------------------------------
    std::vector<utree_move_info> utree_nni_spr_stack_left;
    std::vector<utree_move_info> utree_nni_spr_stack_right;
    std::vector<utree_move_info> utree_nni_spr_stack_up;
    node *u_node;

    u_node = utree_pseudo_root;
    //utree_get_list_nodes_within_radius(u_node, radius,utree_nni_spr_stack_left,utree_nni_spr_stack_right,utree_nni_spr_stack_up);

    //for(unsigned int k=0;k<utree.size();k++) {


    //---------------------------------------------------------------
//    u_node=utree.at(0);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node,
//                                       radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_right,
//                                       utree_nni_spr_stack_up);

    u_node = utree.at(5);
    std::cout << "UNODE:" << u_node->data->getName() << std::endl;
    utree_get_list_nodes_within_radius(u_node, radius,
                                       utree_nni_spr_stack_left,
                                       utree_nni_spr_stack_right,
                                       utree_nni_spr_stack_up);
//    u_node=utree.at(6);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node, radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left);
//    u_node=utree.at(9);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node, radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left);
//    u_node=utree.at(10);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node, radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left);
//    u_node=utree.at(11);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node, radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left);
//    u_node=utree.at(12);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node, radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left);
//    u_node=utree.at(13);
//    std::cout<<"UNODE:"<<u_node->data->getName()<<"\n";
//    utree_get_list_nodes_within_radius(u_node, radius,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left,
//                                       utree_nni_spr_stack_left);
    //---------------------------------------------------------------
    std::cout << "size list_left:" << utree_nni_spr_stack_left.size() << std::endl;
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        std::cout << "list[" << i << "]=(" << (utree_nni_spr_stack_left.at(i)).node1->data->getName() << ";"
                  << (utree_nni_spr_stack_left.at(i)).node2->data->getName() << ")\n";
    }

    std::cout << "size list_right:" << utree_nni_spr_stack_right.size() << std::endl;
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        std::cout << "list[" << i << "]=(" << (utree_nni_spr_stack_right.at(i)).node1->data->getName() << ";"
                  << (utree_nni_spr_stack_right.at(i)).node2->data->getName() << ")\n";
    }

    std::cout << "size list_up:" << utree_nni_spr_stack_up.size() << std::endl;
    for (unsigned int i = 0; i < utree_nni_spr_stack_up.size(); i++) {
        std::cout << "list[" << i << "]=(" << (utree_nni_spr_stack_up.at(i)).node1->data->getName() << ";"
                  << (utree_nni_spr_stack_up.at(i)).node2->data->getName() << ")\n";
    }
    //------------------------------------------------------------------------------------------------------------------
    // PERFORM SPR MOVES and RECOMPUTE logLK

    //move_info m;
    //move_info n;
    int max_idx;
    double max_val;
    std::vector<PhyTree *> p;
    bool valid_move;
    utree_move_info un;
    //std::vector<node *> utree_copy;

    node *p_child_1;
    node *p_child_2;
    node *q_child;

    /*
    std::cout<<"UTREE:\n";
    print_utree(utree.at(0));
    std::cout<<"\n\n";
    */

    node *source;

    max_idx = -1;
    max_val = -INFINITY;
    //---------------------------------------------------------------------
    std::cout << "processing left" << std::endl;
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        un = utree_nni_spr_stack_left.at(i);

        source = un.node1;

        std::cout << "node_1:" << std::endl;
        print_node_neighbours(source);
        std::cout << "node_2:" << std::endl;
        print_node_neighbours(un.node2);
        std::cout << std::endl;


        p_child_1 = source->next->back;
        p_child_2 = source->next->next->back;
        q_child = un.node2->back;

        valid_move = tree->utree_swap(source, un.node2);

        if (valid_move) {

            std::cout << "-------------" << std::endl;
            std::cout << utree_formatNewick(utree.at(0)) << std::endl;
            std::cout << "-------------" << std::endl;

            source->next->back = p_child_1;
            p_child_1->back = source->next;
            source->next->next->back = p_child_2;
            p_child_2->back = source->next->next;
            un.node2->back = q_child;
            q_child->back = un.node2;

            p.clear();
        } else {
            std::cout << "I am skipping this..." << std::endl;
        }
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
        un = utree_nni_spr_stack_left.at(i);

        source = un.node1->next->next;

        p_child_1 = source->next->back;
        p_child_2 = source->next->next->back;
        q_child = un.node2->back;

        valid_move = tree->utree_swap(source, un.node2);

        if (valid_move) {

            std::cout << "-------------" << std::endl;
            std::cout << utree_formatNewick(utree.at(0)) << std::endl;
            std::cout << "-------------" << std::endl;

            source->next->back = p_child_1;
            p_child_1->back = source->next;
            source->next->next->back = p_child_2;
            p_child_2->back = source->next->next;
            un.node2->back = q_child;
            q_child->back = un.node2;

            p.clear();
        } else {
            std::cout << "I am skipping this..." << std::endl;
        }
    }
    //---------------------------------------------------------------------
    std::cout << "processing right" << std::endl;
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        un = utree_nni_spr_stack_right.at(i);

        source = un.node1;

        std::cout << "node_1:" << std::endl;
        print_node_neighbours(source);
        std::cout << "node_2:" << std::endl;
        print_node_neighbours(un.node2);
        std::cout << std::endl;


        p_child_1 = source->next->back;
        p_child_2 = source->next->next->back;
        q_child = un.node2->back;

        valid_move = tree->utree_swap(source, un.node2);

        if (valid_move) {

            std::cout << "-------------" << std::endl;
            std::cout << utree_formatNewick(utree.at(0)) << std::endl;
            std::cout << "-------------" << std::endl;

            source->next->back = p_child_1;
            p_child_1->back = source->next;
            source->next->next->back = p_child_2;
            p_child_2->back = source->next->next;
            un.node2->back = q_child;
            q_child->back = un.node2;

            p.clear();
        } else {
            std::cout << "I am skipping this..." << std::endl;
        }
    }
    for (unsigned int i = 0; i < utree_nni_spr_stack_right.size(); i++) {
        un = utree_nni_spr_stack_right.at(i);

        source = un.node1->next->next;

        p_child_1 = source->next->back;
        p_child_2 = source->next->next->back;
        q_child = un.node2->back;

        valid_move = tree->utree_swap(source, un.node2);

        if (valid_move) {

            std::cout << "-------------" << std::endl;
            std::cout << utree_formatNewick(utree.at(0)) << std::endl;
            std::cout << "-------------" << std::endl;

            source->next->back = p_child_1;
            p_child_1->back = source->next;
            source->next->next->back = p_child_2;
            p_child_2->back = source->next->next;
            un.node2->back = q_child;
            q_child->back = un.node2;

            p.clear();
        } else {
            std::cout << "I am skipping this..." << std::endl;
        }
    }
    //---------------------------------------------------------------------


//    for (unsigned int i = 0; i < utree_nni_spr_stack_left.size(); i++) {
//        // perform SPR move
//
//        //std::cout << "\n\nPerform SPR move\n";
//        un = utree_nni_spr_stack_left.at(i);
//
//        /*
//        if(un.node1->next!=NULL) {
//            un.node1 = un.node1->next->next;
//        }
//        if(un.node2->next!=NULL) {
//            un.node2 = un.node2->next;
//        }
//        */
//
////        std::cout << "node_1=" << un.node1->data->getName();
////        if(un.node1->next!=NULL) {
////            std::cout << "(" << un.node1->back->data->getName() << ";";
////            std::cout << un.node1->next->back->data->getName()<< ";";
////            std::cout << un.node1->next->next->back->data->getName()<<")";
////        }
////        std::cout<<"\n";
////
////        std::cout << "node_2=" << un.node2->data->getName();
////        if(un.node2->next!=NULL) {
////            std::cout << "(" << un.node2->back->data->getName() << ";";
////            std::cout << un.node2->next->back->data->getName()<< ";";
////            std::cout << un.node2->next->next->back->data->getName()<<")";
////        }
////        std::cout<<"\n";
//
//
////        std::cout << "node_1:\n";
////        print_node_neighbours(un.node1);
////        std::cout << "node_2:\n";
////        print_node_neighbours(un.node2);
////        std::cout<<"\n";
//
//
//
//        //utree_copy=utree;
//        //std::memcpy(&utree_copy,&utree,sizeof(node *)*(utree.size()));
////        typedef typename std::vector<node *>::const_iterator iter;
////        iter it1=utree.begin();
////        iter it2=utree.end();
////        utree_copy=std::copy(it1,it2,node *);
//
//        //utree_pseudo_root->ID=1234;
//
//        //copy_vector(utree_copy,utree);
//
////        node *n1=utree.at(0);
////        node *n2=utree_copy.at(0);
////
////        std::cout<<n1->ID<<";"<<n2->ID<<"\n";
////        utree_pseudo_root->ID=4321;
////        std::cout<<n1->ID<<";"<<n2->ID<<"\n";
//
//        source=un.node1;
//
//
//        p_child_1=source->next->back;
//        p_child_2=source->next->next->back;
//        q_child=un.node2->back;
//
//        valid_move = tree->utree_swap(source,un.node2);
//
//
//        if (valid_move) {
//
//            // print newick
//            //std::cout << "after SPR move\n";
//            //std::cout << tree->formatNewick() << "\n";
//            std::cout<<"-------------\n";
//            //print_utree(utree.at(0));
//            std::cout<<utree_formatNewick(utree.at(0))<<"\n";
//            std::cout<<"-------------\n";
//
//
//
//
//
//            // get all nodes in the SPR path
//            //p = get_path_from_nodes(n.node1, n.node2);
//
//            // update all fv values
//            // TODO: Refactor as method on node
//            //update_fv_values(p, extended_alphabet_size);
//
//            //TODO recompute the sum
//
//            // store index of max
////            if (n.lk > max_val) {
////                max_val = n.lk;
////                max_idx = i;
////            }
//
//
//
//            // rollback SPR move
//            //std::cout << "Perform rollback\n";
//
//
//            source->next->back=p_child_1;
//            p_child_1->back=source->next;
//            source->next->next->back=p_child_2;
//            p_child_2->back=source->next->next;
//            un.node2->back=q_child;
//            q_child->back=un.node2;
//
////            un.node2=un.node2->next->next;
////
////
////            std::cout << "node_1:\n";
////            print_node_neighbours(un.node1);
////            std::cout << "node_2:\n";
////            print_node_neighbours(un.node2);
////            std::cout<<"\n";
////
////            tree->utree_swap(un.node1,un.node2);
//
////            for(unsigned int j=utree.size()-1;j>=0;j--){
////                node* r_node = utree.at(j);
////                delete r_node;
////            }
////            utree.clear();
//
//            //utree = utree_copy;
//
//            //std::cout<<utree.at(0)->data->getName()<<"\n";
//
////            // print newick
//            /*
//            std::cout << "after rollback\n";
//            std::cout << tree->formatNewick() << "\n";
//            std::cout<<"-------------\n";
//            print_utree(utree.at(0));
//            std::cout<<"-------------\n";
//            */
//
//            p.clear();
//        }
//    }


    std::cout << "max_val:" << max_val << " at index: " << max_idx << std::endl;

    //	nni_spr_stack.pop_back();
    utree_nni_spr_stack_left.empty();
    utree_nni_spr_stack_right.empty();
    utree_nni_spr_stack_up.empty();


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




/*
    node = tree;
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_right_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child()->get_right_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child()->get_left_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child()->get_left_child()->get_right_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child()->get_left_child()->get_left_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child()->get_left_child()->get_left_child()->get_left_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
    node = tree->get_left_child()->get_left_child()->get_left_child()->get_right_child();
    radius = 4;
    get_list_nodes_within_radius(node, radius, nni_spr_stack);
     */
//----------------------------------------------------
/*
move_info m;
m.node1 = tree->get_left_child()->get_left_child()->get_left_child()->get_left_child();
m.node2 = tree->get_right_child();
nni_spr_stack.push_back(m);
*/
//----------------------------------------------------
/*
std::cout << "size list_left:" << nni_spr_stack_left.size() << "\n";
for (unsigned int i = 0; i < nni_spr_stack_left.size(); i++) {
    std::cout << "list[" << i << "]=(" << (nni_spr_stack_left.at(i)).node1->getName() << ";"
              << (nni_spr_stack_left.at(i)).node2->getName() << ")\n";
}

std::cout << "size list_right:" << nni_spr_stack_right.size() << "\n";
for (unsigned int i = 0; i < nni_spr_stack_right.size(); i++) {
    std::cout << "list[" << i << "]=(" << (nni_spr_stack_right.at(i)).node1->getName() << ";"
              << (nni_spr_stack_right.at(i)).node2->getName() << ")\n";
}

std::cout << "size list_up:" << nni_spr_stack_up.size() << "\n";
for (unsigned int i = 0; i < nni_spr_stack_up.size(); i++) {
    std::cout << "list[" << i << "]=(" << (nni_spr_stack_up.at(i)).node1->getName() << ";"
              << (nni_spr_stack_up.at(i)).node2->getName() << ")\n";
}
 */



//    PhyTree *ctree;
//
//    max_idx = -1;
//    max_val = -INFINITY;
//    for (unsigned int i = 0; i < nni_spr_stack_left.size(); i++) {
//
//    // perform SPR move
//    std::cout << "\n\nPerform SPR move\n";
//    n = nni_spr_stack_left.at(i);
//
//        for(unsigned int child=0;child<2;child++) {
//
//        ctree = tree->copy();
//
//        //std::cout<<"ID: "<<n.ID<<"\n";
//        std::cout << "n.t1=" << n.node1->getName() << " : n.t2=" << n.node2->getName() <<" position :"<<child<< "\n";
//
//        valid_move = tree->swap_left(n.node1, n.node2,child);
//
//        if (valid_move) {
//            // print newick
//            std::cout << "after SPR move\n";
//            std::cout << tree->formatNewick() << "\n";
//
//            /*
//
//            // get all nodes in the SPR path
//            p = get_path_from_nodes(n.node1, n.node2);
//
//            // update all fv values
//            // TODO: Refactor as method on node
//            update_fv_values(p, extended_alphabet_size);
//
//            //TODO recompute the sum
//
//            // store index of max
//            if (n.log_col_lk > max_val) {
//                max_val = n.log_col_lk;
//                max_idx = i;
//            }
//
//            */
//
//            // rollback SPR move
//            std::cout << "Perform SPR move rollback\n";
//            //tree->swap2(n.node1, n.node2);
//
//            delete tree;
//            tree = ctree;
//
//            // print newick
//            std::cout << "after rollback\n";
//            std::cout << tree->formatNewick() << "\n";
//
//            p.clear();
//            }
//        }
//    }

