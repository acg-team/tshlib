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
 * @author Lorenzo Gatti & Massimo Maiolo
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
#include "PhyTree.hpp"
#include "newick.hpp"
#include "nni_spr.hpp"

//=======================================================================================================
//DP-PIP
std::string create_col_MSA(std::vector<std::pair<std::string,std::string>> &MSA,int index){
    std::string colMSA;

    for(unsigned int i=0;i<MSA.size();i++){
        colMSA.append(MSA.at(i).second,index,1);
    }

    return colMSA;
}
//=======================================================================================================
//DP-PIP
Eigen::VectorXd go_down(PhyTree &tree,int is_DNA_AA_Codon,int dim_alphabet){
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;
    char ch;

    if(tree.isLeaf()){

        fv=Eigen::VectorXd::Zero(dim_alphabet+1);
        int idx;

        if(is_DNA_AA_Codon==1){
            ch=tree.get_leaf_character();
            idx=mytable[(int)ch];
        }else if(is_DNA_AA_Codon==2){
            ch=tree.get_leaf_character();
            idx=mytableAA[(int)ch];
        }else{
            perror("go_down not implemented for codon model yet\n");
            exit(EXIT_FAILURE);
        }
        idx=idx<0?dim_alphabet:idx;
        fv[idx]=1.0;
    }else{

        fvL=go_down(tree[0],is_DNA_AA_Codon,dim_alphabet);
        fvR=go_down(tree[1],is_DNA_AA_Codon,dim_alphabet);

        fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

    }

    return fv;
}
//=======================================================================================================
//DP-PIP
/*
double compute_col_lk(PhyTree &tree,
                      std::string &MSA_col,
                      Eigen::VectorXd &pi,
                      int is_DNA_AA_Codon){


    double pr;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;
    Eigen::VectorXd fv;
    double fv0;
    int dim_alphabet;

    if(is_DNA_AA_Codon==1){
        dim_alphabet=4;
    }else if(is_DNA_AA_Codon==2){
        dim_alphabet=20;
    }else if(is_DNA_AA_Codon==3){
        dim_alphabet=61;
    }else{
        perror("ERROR: alphabet not recognized\n");
        exit(EXIT_FAILURE);
    }

    fvL=go_down(tree[0],is_DNA_AA_Codon,dim_alphabet);
    fvR=go_down(tree[1],is_DNA_AA_Codon,dim_alphabet);

    fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

    fv0=fv.dot(pi);

    pr=tree.get_iota()*tree.get_beta()*fv0;

    pr=log(pr);

    return pr;
}
*/
//=======================================================================================================
//DP-PIP
Eigen::VectorXd compute_lk_empty_col(PhyTree &node,double &lk,Eigen::VectorXd &pi,int is_DNA_AA_Codon,int dim_extended_alphabet){
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;

    //std::cout<<"node="<<node.getName()<<" fv:"<<fv<<"\n";

    if(node.isLeaf()){
        fv=Eigen::VectorXd::Zero(dim_extended_alphabet);

        fv[dim_extended_alphabet-1]=1.0;

        lk+=node.get_iota()*(1-node.get_beta()+node.get_beta()*(fv.dot(pi)));

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_setA()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";*/

        return fv;
    }else{

        fvL=compute_lk_empty_col(node[0],lk,pi,is_DNA_AA_Codon,dim_extended_alphabet);
        fvR=compute_lk_empty_col(node[1],lk,pi,is_DNA_AA_Codon,dim_extended_alphabet);

        fv=(node.get_left_child()->get_Pr()*fvL).cwiseProduct(node.get_right_child()->get_Pr()*fvR);

        lk+=node.get_iota()*(1-node.get_beta()+node.get_beta()*(fv.dot(pi)));

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_setA()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";
*/
        return fv;
    }

}
//=======================================================================================================
//DP-PIP
Eigen::VectorXd compute_lk_recursive(PhyTree &node,double &lk,Eigen::VectorXd &pi,int is_DNA_AA_Codon,int dim_extended_alphabet){
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;

    //std::cout<<"node="<<node.getName()<<" fv:"<<fv<<"\n";

    if(node.isLeaf()){
        fv=Eigen::VectorXd::Zero(dim_extended_alphabet);

        //fv[0]=1.0;
        int idx;

        if(is_DNA_AA_Codon==1){
            idx=mytable[(int)node.get_leaf_character()];
        }else if(is_DNA_AA_Codon==2){
            idx=mytableAA[(int)node.get_leaf_character()];
        }else{
            perror("not implemented for codon model yet\n");
            exit(EXIT_FAILURE);
        }

        idx=idx<0?dim_extended_alphabet-1:idx;
        fv[idx]=1.0;

        if(node.get_setA()){
            lk+=node.get_iota()*node.get_beta()*(fv.dot(pi));
        }

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_setA()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";*/

        node.set_MSA_fv(fv);

        return fv;
    }else{

        fvL=compute_lk_recursive(node[0],lk,pi,is_DNA_AA_Codon,dim_extended_alphabet);
        fvR=compute_lk_recursive(node[1],lk,pi,is_DNA_AA_Codon,dim_extended_alphabet);

        fv=(node.get_left_child()->get_Pr()*fvL).cwiseProduct(node.get_right_child()->get_Pr()*fvR);

        if(node.get_setA()){
            lk+=node.get_iota()*node.get_beta()*(fv.dot(pi));
        }

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_setA()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";
*/

        node.set_MSA_fv(fv);

        return fv;
    }

}
//=======================================================================================================
//DP-PIP
double compute_col_lk(  PhyTree &tree,
                        Eigen::VectorXd &pi,
                        int is_DNA_AA_Codon,
                        int alphabet_size){

    double lk;

    compute_lk_recursive(tree,lk,pi,is_DNA_AA_Codon,alphabet_size);

    return log(lk);
}
//===================================================================================================================
//double recompute_lk(PhyTree *tree,double k){
//    return k;
//}
//===================================================================================================================
void nodes_within_radius(PhyTree *start_node,PhyTree *node,int radius,bool save,std::vector<move_info> &list_nodes){

    if(!save){
        save=true;
    }else{
        move_info m;
        m.node1=start_node;
        m.node2=node;
        list_nodes.push_back(m);
    }

    if(radius<=0){
        return;
    }

    if(!node->isLeaf()){
        radius --;
        nodes_within_radius(start_node,node->get_left_child(),radius,save,list_nodes);
        nodes_within_radius(start_node,node->get_right_child(),radius,save,list_nodes);
    }

}
//===================================================================================================================
void nodes_within_radius_up(PhyTree *start_node,PhyTree *node,int radius,int direction,std::vector<move_info> &list_nodes){
    index_t idx;

    //TODO: check binary tree condition!

    move_info m;
    m.node1=start_node;
    m.node2=node;
    list_nodes.push_back(m);

    if(radius<=0){
        return;
    }

    radius --;
    if(direction==0){
        if(node->getParent()!=NULL){
            idx=node->indexOf();
            nodes_within_radius_up(start_node,node->getParent(),radius,idx,list_nodes);
        }
        nodes_within_radius(start_node,node->get_right_child(),radius,true,list_nodes);
    }else if(direction==1){
        if(node->getParent()!=NULL){
            idx=node->indexOf();
            nodes_within_radius_up(start_node,node->getParent(),radius,idx,list_nodes);
        }
        nodes_within_radius(start_node,node->get_left_child(),radius,true,list_nodes);
    }

}
//===================================================================================================================
void get_list_nodes_within_radius(PhyTree *node,int radius,std::vector<move_info> &list_nodes){
    bool save;

    save=false;

    nodes_within_radius(node,node,radius,save,list_nodes);

    if(node->getParent()!=NULL){
        nodes_within_radius_up(node,node->getParent(),radius,node->indexOf(),list_nodes);
    }

}
//===================================================================================================================
std::vector<PhyTree *> fill_with_nodes(PhyTree *n){
    std::vector<PhyTree *> list_nodes_n;

    list_nodes_n.push_back(n);
    while(n->getParent()!=NULL){
        n=n->getParent();
        list_nodes_n.push_back(n);
    }

    return list_nodes_n;
}
//===================================================================================================================
std::vector<PhyTree *> get_unique(std::vector<PhyTree *> &list_nodes_n1,std::vector<PhyTree *> &list_nodes_n2){
    std::vector<PhyTree *> list_nodes;
    PhyTree *n1;
    PhyTree *n2;

    while(list_nodes_n1.size()>0 && list_nodes_n2.size()>0){
        n1=list_nodes_n1.at(list_nodes_n1.size()-1);
        n2=list_nodes_n2.at(list_nodes_n2.size()-1);
        if(n1==n2){
            list_nodes.push_back(n1);
            list_nodes_n1.pop_back();
            list_nodes_n2.pop_back();
        }else{
            break;
        }
    }

    while(list_nodes_n1.size()>0){
        n2=list_nodes_n1.at(list_nodes_n1.size()-1);
        list_nodes.push_back(n1);
        list_nodes_n1.pop_back();
    }

    while(list_nodes_n2.size()>0){
        n2=list_nodes_n2.at(list_nodes_n2.size()-1);
        list_nodes.push_back(n2);
        list_nodes_n2.pop_back();
    }

    std::reverse(list_nodes.begin(), list_nodes.end());

    return list_nodes;
}
//===================================================================================================================
double compute_nu(double tau,double lambda,double mu){

    if(fabs(mu)<1e-8){
        perror("ERROR in compute_nu: mu too small");
    }

    return lambda*(tau+1/mu);
}
//===================================================================================================================
std::vector<PhyTree *> get_path_from_nodes(PhyTree *n1,PhyTree *n2){
    std::vector<PhyTree *> list_nodes_n0;
    std::vector<PhyTree *> list_nodes_n1;
    std::vector<PhyTree *> list_nodes_n2;

    // add nodes from n1 to root
    list_nodes_n1=fill_with_nodes(n1);

    // add nodes from n2 to root
    list_nodes_n2=fill_with_nodes(n2);

    list_nodes_n0=get_unique(list_nodes_n1,list_nodes_n2);

    return list_nodes_n0;
}
//===================================================================================================================
double phi(int m,double nu,double p0){
    double p;
    double log_factorial_m;

    log_factorial_m=0;
    for(int i=1;i<=m;i++){
        log_factorial_m+=log(i);
    }

    p=-log_factorial_m+m*log(nu)+(nu*(p0-1));

    return p;
}
//===================================================================================================================
int main(int argc, char** argv)
{
    PhyTree *t1;
    PhyTree *t2;
    std::string tree_file="/home/max/PIP_C++/NNI_SPR/tree_5_leaves_r_bl.nwk";
    PhyTree* tree = NULL;
    double mu;
    double lambda;
    double tau;
    double nu;

    mu=0.1;
    lambda=0.2;

    //----------------------------------------------------------
    // INIT TREE

    // tree filename
    std::ifstream tree_str(tree_file.c_str());

    // read newick file
    tree = newick_parser::parse_newick(&tree_str);

    // set name of internal nodes
    tree->set_missing_node_name("V");

    // compute total tree length
    tau=tree->computeLength();

//    std::cout<<tau<<"\n";

    // compute the normalizing Poisson intensity
    nu=compute_nu(tau,lambda,mu);

    // set insertion probability to each node
    tree->set_iota(tau,mu);

    // set survival probability to each node
    tree->set_beta(tau,mu);

	std::cout<<tree->formatNewick()<<"\n\n";
	tree->print();
	std::cout<<"\n";
    //----------------------------------------------------------
    // LOAD MSA

    std::vector< std::pair<std::string,std::string> > MSA;

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

    std::string seq1_label="A";
    std::string seq2_label="B";
    std::string seq3_label="C";
    std::string seq4_label="D";
    std::string seq5_label="E";

    std::string seq1_DNA="A";
    std::string seq2_DNA="-";
    std::string seq3_DNA="-";
    std::string seq4_DNA="-";
    std::string seq5_DNA="-";

    MSA.push_back(std::make_pair(seq1_label,seq1_DNA));
    MSA.push_back(std::make_pair(seq2_label,seq2_DNA));
    MSA.push_back(std::make_pair(seq3_label,seq3_DNA));
    MSA.push_back(std::make_pair(seq4_label,seq4_DNA));
    MSA.push_back(std::make_pair(seq5_label,seq5_DNA));
    //----------------------------------------------------------
    // INITIAL LIKELIHOOD COMPUTATION
    int is_DNA_AA_Codon;
    int extended_alphabet_size;

    is_DNA_AA_Codon=1; // 1:DNA, 2:AA, 3:Codon
    extended_alphabet_size=5; // DNA alphabet

    // set "pseudo" probability matrix
    tree->tmp_initPr(extended_alphabet_size);

    int MSA_len;
    double lk;
    Eigen::VectorXd pi;

    // set Pi, steady state frequencies
    pi=Eigen::VectorXd::Zero(extended_alphabet_size);
    pi[0]=0.25;
    pi[1]=0.25;
    pi[2]=0.25;
    pi[3]=0.25;
    pi[4]=0.0;

    // get MSA length
    MSA_len=MSA.at(0).second.size();

    //std::cout<<"MSA_len="<<MSA_len<<"\n";

    double LK=0;

    // compute lk
    for(int i=0;i<MSA_len;i++){

        // extract MSA column
        std::string s=create_col_MSA(MSA,i);
        //std::cout<<"col["<<i<<"]="<<s<<"\n";

        // assign char at the leaves
        set_leaf_state(tree,s);
        //print_leaf_state(tree);

        // set ancestral flag (1=plausible insertion location, 0=not plausible insertion location)
        set_ancestral_flag(tree,s);
        //print_descCount(tree);
        //print_ancestral_flag(tree);

        tree->clear_fv();

        // compute column likelihood
        lk=compute_col_lk(*tree,pi,is_DNA_AA_Codon,extended_alphabet_size);

        std::cout<<"col_lk="<<lk<<"\n";

        LK+=lk;
    }

    double p0;
    compute_lk_empty_col(*tree,p0,pi,is_DNA_AA_Codon,extended_alphabet_size);
    p0=log(p0);

    std::cout<<"p0="<<p0<<"\n";

    LK+=phi(MSA_len,nu,p0);
    //----------------------------------------------------------
    // GET ALL NODES WITHIN RADIUS

    int radius;
    PhyTree* node;
    std::vector<move_info> nni_spr_stack;

    node=tree->get_left_child();
    radius=3;

    get_list_nodes_within_radius(node,radius,nni_spr_stack);

    std::cout<<"size list:"<<nni_spr_stack.size()<<"\n";

    for(unsigned int i=0;i<nni_spr_stack.size();i++){
        std::cout<<"list["<<i<<"]=("<<(nni_spr_stack.at(i)).node1->getName()<<";"<<(nni_spr_stack.at(i)).node2->getName()<<")\n";
    }
    //----------------------------------------------------------
    // PERFORM SPR MOVES and RECOMPUTE LK

    move_info m;
    move_info n;
    int max_idx;
    double max_val;
    std::vector<PhyTree *> p;
    bool valid_move;

    max_val=-INFINITY;
    for(unsigned int i=0;i<nni_spr_stack.size();i++){

        // perform SPR move
        std::cout<<"\n\nPerform SPR move\n";
        n= nni_spr_stack.at(i);

        //std::cout<<"ID: "<<n.ID<<"\n";
        std::cout<<"n.t1="<<n.node1->getName()<<" : n.t2="<<n.node2->getName()<<"\n";
        valid_move=tree->swap2(n.node1,n.node2);

        if(valid_move){
            // print newick
            std::cout << "after SPR move\n";
            std::cout << tree->formatNewick() << "\n";

            /*
            // get all nodes in the SPR path
            p = get_path_from_nodes(n.node1, n.node2);

            // update all fv values
            update_fv_values(p, extended_alphabet_size); //TODO recompute the sum

            // store index of max
            if (n.lk > max_val) {
                max_val = n.lk;
                max_idx = i;
            }
             */

            // rollback SPR move
            std::cout << "Perform SPR move rollback\n";
            tree->swap2(n.node1, n.node2);

            // print newick
            std::cout << "after rollback\n";
            std::cout << tree->formatNewick() << "\n";

            p.clear();
        }
    }

    std::cout<<"max_val:"<<max_val<<" at index: "<<max_idx<<"\n";

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





