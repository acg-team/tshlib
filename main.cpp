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
#include "nni_spr.h"

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
//template <class ALPHABET>
//double compute_lk_down(PhyTree &tree,std::string &s,Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,int ugglyflag){
//
//	double pr;
//	int idx;
//	Eigen::VectorXd fvL;
//	Eigen::VectorXd fvR;
//	Eigen::VectorXd fv;
//	double fv0;
//
//	if(tree.isLeaf()){
//
//		idx=0;
//		fv=go_down<ALPHABET>(tree,s,idx,ugglyflag);
//		fv0=fv.dot(pi);
//		pr=tree.get_iota()*tree.get_beta()*fv0;
//
//		return pr;
//
//	}else{
//
//		idx=0;
//		fv=go_down<ALPHABET>(tree,s,idx,ugglyflag);
//		fv0=fv.dot(pi);
//		pr=tree.get_iota()*tree.get_beta()*fv0;
//
//		bool flagL=true;
//		bool flagR=true;
//		idx=0;
//		allgaps<ALPHABET>(tree[0],s,idx,flagL);
//		int ixx=idx;
//		allgaps<ALPHABET>(tree[1],s,idx,flagR);
//
//		int len;
//		if(flagR){
//			std::string sL;//=stringFromSequence(s);
//			len=ixx;
//			sL=s.substr(0,len);
//			return pr + compute_lk_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
//		}
//
//		if(flagL){
//			std::string sR;//=stringFromSequence(s);
//			sR=s.substr(ixx);
//			return pr + compute_lk_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
//		}
//
//	}
//
//	return pr;
//}
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
//=======================================================================================================
//DP-PIP
Eigen::VectorXd compute_lk_recursive(PhyTree &node,double &lk,Eigen::VectorXd &pi,int is_DNA_AA_Codon,int dim_alphabet){
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;
    std::cout<<"node="<<node.getName()<<" fv:"<<fv<<"\n";

    if(node.isLeaf()){
        fv=Eigen::VectorXd::Zero(dim_alphabet+1);

        fv[0]=1.0;

        if(node.get_setA()){
            lk+=node.get_iota()*node.get_beta()*(fv.dot(pi));
        }

        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_setA()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";

        return fv;
    }else{

        fvL=compute_lk_recursive(node[0],lk,pi,is_DNA_AA_Codon,dim_alphabet);
        fvR=compute_lk_recursive(node[1],lk,pi,is_DNA_AA_Codon,dim_alphabet);

        fv=(node.get_left_child()->get_Pr()*fvL).cwiseProduct(node.get_right_child()->get_Pr()*fvR);

        if(node.get_setA()){
            lk+=node.get_iota()*node.get_beta()*(fv.dot(pi));
        }

        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_setA()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";
        return fv;
    }

}
//=======================================================================================================
//DP-PIP
double compute_col_lk_prova(PhyTree &tree,
                            std::string &MSA_col,
                            Eigen::VectorXd &pi,
                            int is_DNA_AA_Codon){


    int dim_alphabet;
    double lk;
//	Eigen::VectorXd fv;

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

    compute_lk_recursive(tree,lk,pi,is_DNA_AA_Codon,dim_alphabet);

    return lk;
}
//===================================================================================================================
double recompute_lk(PhyTree *tree,double k){
    return k;
}
//===================================================================================================================
void nodes_within_radius(PhyTree *node,int radius,bool save,std::vector<PhyTree *> &list_nodes){

    if(!save){
        save=true;
    }else{
        list_nodes.push_back(node);
//		std::cout<<"nodes_within_radius: saving node = "<<node->getName()<<" ; radius "<<radius<<"\n";
    }

    if(radius<=0){
        return;
    }

    if(!node->isLeaf()){
        radius --;
        nodes_within_radius(node->get_left_child(),radius,save,list_nodes);
        nodes_within_radius(node->get_right_child(),radius,save,list_nodes);
    }

}
//===================================================================================================================
void nodes_within_radius_up(PhyTree *node,int radius,int direction,std::vector<PhyTree *> &list_nodes){
    index_t idx;

    //TODO: check binary tree condition!

    list_nodes.push_back(node);

    if(radius<=0){
        return;
    }

    radius --;
    if(direction==0){
        if(node->getParent()!=NULL){
            idx=node->indexOf();
            nodes_within_radius_up(node->getParent(),radius,idx,list_nodes);
        }
        nodes_within_radius(node->get_right_child(),radius,true,list_nodes);
    }else if(direction==1){
        if(node->getParent()!=NULL){
            idx=node->indexOf();
            nodes_within_radius_up(node->getParent(),radius,idx,list_nodes);
        }
        nodes_within_radius(node->get_left_child(),radius,true,list_nodes);
    }

}
//===================================================================================================================
void get_list_nodes_within_radius(PhyTree *node,int radius,std::vector<PhyTree *> &list_nodes){
    bool save;

    save=false;

    nodes_within_radius(node,radius,save,list_nodes);

    if(node->getParent()!=NULL){
        nodes_within_radius_up(node->getParent(),radius,node->indexOf(),list_nodes);
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

    std::ifstream tree_str(tree_file.c_str());

    tree = newick_parser::parse_newick(&tree_str);

    tree->set_missing_node_name("V");

    tau=tree->computeLength();

    nu=compute_nu(tau,lambda,mu);

    tree->set_tau(tau);

    tree->set_nu(nu);

    tree->set_iota(tau,mu);

    tree->set_beta(tau,mu);

//	std::cout<<tree->formatNewick()<<"\n\n";
//	tree->print();
//	std::cout<<"\n";
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
    std::vector< std::pair<std::string,std::string> > MSA;

    std::string seq1_label="A";
    std::string seq1_DNA="ACGT";
    std::string seq2_label="B";
    std::string seq2_DNA="TAGC";
    std::string seq3_label="C";
    std::string seq3_DNA="GCAT";
    std::string seq4_label="D";
    std::string seq4_DNA="CGTA";
    std::string seq5_label="E";
    std::string seq5_DNA="GTCA";

    MSA.push_back(std::make_pair(seq1_label,seq1_DNA));
    MSA.push_back(std::make_pair(seq2_label,seq2_DNA));
    MSA.push_back(std::make_pair(seq3_label,seq3_DNA));
    MSA.push_back(std::make_pair(seq4_label,seq4_DNA));
    MSA.push_back(std::make_pair(seq5_label,seq5_DNA));

    int alphabet_size=5;

    tree->tmp_initPr(alphabet_size);

    int MSA_len;
    double lk;
    Eigen::VectorXd pi;
    int is_DNA_AA_Codon;

    is_DNA_AA_Codon=1; // 1:DNA, 2:AA, 3:Codon

    pi=Eigen::VectorXd::Zero(alphabet_size);
    pi[0]=0.25;
    pi[1]=0.25;
    pi[2]=0.25;
    pi[3]=0.25;
    pi[4]=0.0;

    MSA_len=MSA.at(0).second.size();

    std::cout<<"MSA_len="<<MSA_len<<"\n";

    for(int i=0;i<MSA_len;i++){

        std::string s=create_col_MSA(MSA,i);

        set_ancestral_flag(tree,s);

        set_leaf_state(tree,s);

        std::cout<<"col["<<i<<"]="<<s<<"\n";

//		lk=compute_col_lk(*tree,s,pi,is_DNA_AA_Codon);

        lk=compute_col_lk_prova(*tree,s,pi,is_DNA_AA_Codon);

        std::cout<<"col_lk="<<lk<<"\n";
    }


    exit(EXIT_SUCCESS);
    //----------------------------------------------------------
    int radius;
    PhyTree* node;
    std::vector<PhyTree *> list_nodes;

    node=tree->get_left_child();
    radius=3;
    get_list_nodes_within_radius(node,radius,list_nodes);


    std::cout<<"size list:"<<list_nodes.size()<<"\n";

    for(unsigned int i=0;i<list_nodes.size();i++){
        std::cout<<"list["<<i<<"]="<<(list_nodes.at(i))->getName()<<"\n";
    }
    //----------------------------------------------------------
    move_info m;
    move_info n;
    std::vector<move_info> nni_spr_stack; // STL Stack object

    t1=tree->get_right_child();
    t2=tree->get_left_child()->get_left_child()->get_left_child();
    m.ID=1;
    m.node1=t1;
    m.node2=t2;
    nni_spr_stack.push_back(m);

    t1=tree->get_left_child();
    t2=tree->get_left_child()->get_left_child();
    m.ID=2;
    m.node1=t1;
    m.node2=t2;
    nni_spr_stack.push_back(m);


    int max_idx;
    double max_val;

    max_val=-INFINITY;
    for(unsigned int i=0;i<nni_spr_stack.size();i++){

        // perform SPR move
        std::cout<<"Perform SPR move\n";
        n= nni_spr_stack.at(i);
        std::cout<<"ID: "<<n.ID<<"\n";
        std::cout<<"n.t1="<<n.node1->getName()<<" : n.t2="<<n.node2->getName()<<"\n";
        tree->swap2(n.node1,n.node2);

        std::cout<<tree->get_right_child()->getName()<<" : ";
        std::cout<<tree->get_right_child()->getParent()->getName()<<" \n";

        // print newick
        std::cout<<"after SPR move\n";
        std::cout<<tree->formatNewick()<<"\n";

        // compute new lk
        n.lk=recompute_lk(tree,i*10);
        // store index max
        if(n.lk>max_val){
            max_val=n.lk;
            max_idx=i;
        }

        // rollback SPR move
        std::cout<<"Perform SPR move rollback\n";
        tree->swap2(n.node1,n.node2);

        // print newick
        std::cout<<"after rollback\n";
        std::cout<<tree->formatNewick()<<"\n";
    }

    std::cout<<"max_val:"<<max_val<<" at index: "<<max_idx<<"\n";

    //	nni_spr_stack.pop_back();
    nni_spr_stack.empty();
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

    return 0;
}





/*
//===================================================================================================================
//===================================================================================================================
int main(int argc, char **argv) {

    //std::string tree_file = "/home/max/PIP_C++/NNI_SPR/tree_5_leaves_r_bl.nwk";
    std::string tree_file = argv[1];
    PhyTree *tree = nullptr;

    std::ifstream tree_str(tree_file.c_str());
    tree = newick_parser::parse_newick(&tree_str);

    tree->set_missing_node_name("V");

    std::cout << tree->formatNewick() << "\n\n";

    tree->print();
    std::cout << "\n";

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
    PhyTree *t1;
    PhyTree *t2;

    t1 = tree->get_right_child();
    t2 = tree->get_left_child()->get_right_child();

//	tree->swap(tree,1,tree->get_left_child(),1);
    tree->swap2(t1, t2);

    std::cout << tree->formatNewick() << "\n\n";
    //----------------------------------------------------------

    return 0;
}*/





