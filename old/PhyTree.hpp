/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
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
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file PhyTree.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 10 2017
 * @version 2.0.2
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

#ifndef TSHLIB_PHYTREE_HPP
#define TSHLIB_PHYTREE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <sstream>
#include <cmath>
//=====================
//DP-PIP
#include <Eigen/Core>
#include <Eigen/src/Core/IO.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
//=====================

class PhyTree;

typedef struct node {
    node *next;
    node *back;
    PhyTree *data;
    int ID;
};




class PhyTree {

private:
    std::string name;
    int node_id;
    double branch_length;
    double branch_support;
    PhyTree *parent;
    std::vector<PhyTree *> children;
    //================================================================
    //DP-PIP
    double iota;
    double beta;
    Eigen::MatrixXd Pr;
    double tau;
    double nu;
    std::vector<Eigen::VectorXd> MSA_fv;
    bool setA;
    int descCount;
    char character;


    void print_prefix(std::string prefix) const {

        //std::cout << prefix << "(" << this->branch_length << ") " << this->name << std::endl;

        //====================================================================================
        //DP-PIP
        std::cout << prefix << "(" << this->branch_length << ") " << this->name << " " << " iota: " << this->iota
                  << " beta: " << this->beta << std::endl;
        //====================================================================================

        for (std::vector<PhyTree *>::const_iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->print_prefix(prefix + "  ");
        }
    }

    std::string formatNewickR() const {
        if (this->isLeaf()) {
            return this->getName();
        } else {
            std::stringstream newick;
            newick << "(";
            std::vector<PhyTree *>::const_iterator i = this->children.begin();
            newick << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
            for (++i; i < this->children.end(); ++i) {
                newick << "," << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
            }
            newick << ")";
            return newick.str();
        }
    }

public:

    PhyTree(std::string name = "") {
        this->parent = NULL;
        this->branch_length = 0;
        this->branch_support = 1;
        this->name = name;
        this->node_id = 0;

        //==============================
        //DP-PIP
        this->iota = 0;
        this->beta = 0;
        this->Pr.resize(0, 0);
        this->tau = 0;
        this->nu = 0;
        //==============================

    }

    ~PhyTree() {
        assert(this->parent == NULL);
        for (std::vector<PhyTree *>::reverse_iterator i = this->children.rbegin(); i < this->children.rend(); ++i) {
            PhyTree *child = *i;
            child->parent = NULL;
            child->branch_length = 0;

            delete child;
        }
    }

    PhyTree *copy() {
        PhyTree *out = new PhyTree();
        out->branch_length = this->branch_length;
        out->branch_support = this->branch_support;
        out->name = this->name;

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            out->addChild((*i)->copy(), (*i)->branch_length, (*i)->branch_support);
        }
        return out;
    }

    void addChild(PhyTree *child, double branch_length = 0, double branch_support = 1) {
        assert(child != this);
        assert(child->parent == NULL);

        this->children.push_back(child);
        child->parent = this;
        child->branch_length = branch_length;
        child->branch_support = branch_support;
    }

    unsigned int indexOf() {
        PhyTree *parent = this->parent;
        assert(parent != NULL);

        for (unsigned int i = 0; i < parent->children.size(); ++i) {
            if (parent->children[i] == this) {
                return i;
            }
        }

        assert(false);
        return -1;
    }

    void pluck() {
        assert(this->parent != NULL);

        unsigned int index = this->indexOf();
        std::vector<PhyTree *>::iterator iter = this->parent->children.begin() + index;

        this->parent->children.erase(iter);
        this->parent = NULL;

        this->branch_length = 0;
        this->branch_support = 1;
    }

    PhyTree *pluckChild(unsigned int index) {
        std::vector<PhyTree *>::iterator iter = this->children.begin() + index;
        PhyTree *child = *iter;

        this->children.erase(iter);
        child->parent = NULL;
        child->branch_length = 0;
        this->branch_support = 1;

        return child;
    }

    void deleteChild(unsigned int index) {
        PhyTree *child = this->pluckChild(index);
        delete child;
    }

    unsigned int countLeaves() {
        if (this->isLeaf()) {
            return 1;
        } else {
            unsigned int leaves = 0;
            for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
                leaves += (*i)->countLeaves();
            }
            return leaves;
        }
    }

    double computeLength() {
        if (this->isLeaf()) {
            return 0;
        } else {
            double length = 0;
            for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
                length += (*i)->branch_length + (*i)->computeLength();
            }
            return length;
        }
    }

    std::string getName() const {
        return this->name;
    }

    int getNodeID() {
        return this->node_id;
    }

    PhyTree *getParent() {
        return this->parent;
    }

    const PhyTree *getParent() const {
        return this->parent;
    }

    double getBranchLength() const {
        return this->branch_length;
    }

    double getBranchSupport() const {
        return this->branch_support;
    }

    PhyTree &operator[](int i) {
        return *this->children[i];
    }

    const PhyTree &operator[](int i) const {
        return *this->children[i];
    }

    std::vector<PhyTree *>::iterator firstChild() {
        return this->children.begin();
    }

    std::vector<PhyTree *>::iterator lastChild() {
        return this->children.end();
    }

    std::vector<PhyTree *>::const_iterator firstChild() const {
        return this->children.begin();
    }

    std::vector<PhyTree *>::const_iterator lastChild() const {
        return this->children.end();
    }

    unsigned int n_children() const {
        return this->children.size();
    }

    bool isLeaf() const {
        return this->children.empty();
    }

    void print() const {
        this->print_prefix("");
    }

    std::string formatNewick() const {
        return this->formatNewickR() + ";";
    }

    //DP-PIP
    std::vector<PhyTree *> get_children() {

        return this->children;
    }

    //DP-PIP
    PhyTree *get_left_child() {

        if (this->n_children() > 0) {
            return this->children[0];
        } else {
            return NULL;
        }

    }

    //DP-PIP
    PhyTree *get_right_child() {

        if (this->n_children() > 0) {
            return this->children[1];
        } else {
            return NULL;
        }

    }

    //DP-PIP
    double get_iota() const {
        return this->iota;
    }

    //DP-PIP
    double get_beta() const {
        return this->beta;
    }

    //DP-PIP
    static bool empty(Eigen::VectorXd &v) {

        if (v.rows() * v.cols() == 0) {
            return true;
        } else {
            return false;
        }
    }

    //DP-PIP
    void setName(std::string name) {
        this->name = name;
    }

    void setNodeID(int id) {
        if (id) {
            this->node_id = id;
        }
    }

    //DP-PIP
    void set_missing_node_name(std::string s) {

        static int num = 0;

        if (this->isLeaf()) {

        } else {

            if (this->name.empty()) {
                std::string name(s);
                name.append(std::to_string(num++));
                this->setName(name);
            }

            //std::cout<<"name: "<<this->name<<"\n";

        }

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_missing_node_name(s);
        }

    }

    //DP-PIP
    const Eigen::MatrixXd &get_Pr() {

        return this->Pr;

    }

    //DP-PIP
    void set_iota(double tau, double mu) {

        if (fabs(mu) < 1e-8) {
            perror("ERROR in set_iota: mu too small");
        }

        double T = tau + 1 / mu;

        if (fabs(T) < 1e-8) {
            perror("ERROR in set_iota: T too small");
        }

        if (this->parent == NULL) {
            this->iota = (1 / mu) / T;
        } else {
            this->iota = this->branch_length / T;
        }

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_iota(tau, mu);
        }

    }

    //DP-PIP
    double get_nu_local() const {
        return this->nu;
    }

    //DP-PIP
    double get_tau_local() const {
        return this->tau;
    }

    //DP-PIP
    void compute_nu_local(double tau, double lambda, double mu) {
        this->nu = lambda * (tau + 1 / mu);
    }

    //DP-PIP
    void compute_tau_local() {

        if (this->isLeaf()) {
            this->tau = this->branch_length;
        } else {
            this->get_left_child()->compute_tau_local();
            this->get_right_child()->compute_tau_local();

            this->tau = this->get_left_child()->tau + this->get_right_child()->tau;
        }

    }

    //DP-PIP
    void set_iota_local_down(double tau, double mu) {

        this->iota = this->branch_length / (tau + 1 / mu);

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_iota_local_down(tau, mu);
        }

    }

    //DP-PIP
    void set_iota_local(double tau, double mu) {

        this->iota = (1 / mu) / (tau + 1 / mu);

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_iota_local_down(tau, mu);
        }


    }

    //DP-PIP
    void set_beta(double tau, double mu) {

        if (fabs(mu) < 1e-8) {
//			error("ERROR in set_beta: mu too small");
        }

        if (this->parent == NULL) {
            this->beta = 1;
        } else {

            if (fabs(this->branch_length) < 1e-8) {
//				error("ERROR in set_beta: branch_length too small");
            }

            this->beta = (1 - exp(-mu * this->branch_length)) / (mu * this->branch_length);
        }

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_beta(tau, mu);
        }

    }

    //DP-PIP
    void set_beta_down(double tau, double mu) {

        if (fabs(mu) < 1e-8) {
//			error("ERROR in set_beta: mu too small");
        }

        if (fabs(this->branch_length) < 1e-8) {
//			error("ERROR in set_beta: branch_length too small");
        }

        this->beta = (1 - exp(-mu * this->branch_length)) / (mu * this->branch_length);

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_beta_down(tau, mu);
        }

    }

    //DP-PIP
    void set_beta_local(double tau, double mu) {

        this->beta = 1.0;

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_beta_down(tau, mu);
        }

    }

    //DP-PIP
    void set_tau(double tau) {

        this->tau = tau;

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_tau(tau);
        }

    }

    //DP-PIP
    void set_nu(double nu) {

        this->nu = nu;

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_nu(nu);
        }

    }

    //DP-PIP
    void print_local_var() {

        std::cout << "----------------------\n";
        std::cout << "Name: " << this->name << "\n";
        if(this->n_children()==2){
            std::cout << "L. child name: " << this->get_left_child()->getName() << "\n";
            std::cout << "R. child name: " << this->get_right_child()->getName() << "\n";
        }
        std::cout << "tau: " << this->tau << "\n";
        std::cout << "nu: " << this->nu << "\n";
        std::cout << "iota: " << this->iota << "\n";
        std::cout << "beta: " << this->beta << "\n";

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->print_local_var();
        }
    }

    //DP-PIP
    void print_br() {

        std::cout << "Name: " << this->name << " bl: " << this->branch_length << std::endl;

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->print_br();
        }
    }

    //DP-PIP
    void printOnlyName() {

        /*
        if(this->isLeaf()){
            std::cout<<"Leaf name: "<<this->name<<std::endl;
        }else{
            std::cout<<"Node name: "<<this->name<<std::endl;
            this->children[0]->printOnlyName();
            this->children[1]->printOnlyName();
        }
        */

        std::cout << "node name: " << this->name << "\n";
        if (this->n_children() == 2) {
            std::cout << "      L child: " << this->get_left_child()->getName() << "\n";
            std::cout << "      R child: " << this->get_right_child()->getName() << "\n";
        }


        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->printOnlyName();
        }


    }

    void null_parent() {
        this->parent = NULL;
    }

    bool swap_up(PhyTree *p, PhyTree *q, unsigned int child) {}

    bool swap_right(PhyTree *p, PhyTree *q, unsigned int child) {}

    bool swap_left(PhyTree *p, PhyTree *q,unsigned int child) {

        if(p->isLeaf() && q->isLeaf()){
            return false;
        }

        PhyTree *p0;
        PhyTree *p1;

        p0=p->children[0];
        p1=p->children[1];

        if(p0->isLeaf() && p1->isLeaf()){
            return false;
        }

        if(p0->isLeaf()){
            PhyTree *tmp;
            tmp=p0;
            p0=p1;
            p1=tmp;
        }

        p->children[0]=q->children[0];
        q->children[0]->parent=p;

        p->children[1]=q;
        q->parent=p;



//        node *p1;
//        node *p2;
//        node *p3;
//        node *p4;
//        node *p5;
//        node *p6;
//        node *p7;
//        node *p8;
//        node *p9;
//
//        //-----------------
//        p1->next=p2;
//        p1->back=NULL;
//        p1->data=p;
//        //-----------------
//        p2->next=p3;
//        p2->back=p4;
//        p2->data=p;
//        //-----------------
//        p3->next=p1;
//        p3->back=p7;
//        p3->data=p;
//        //-----------------
//        p4->next=p5;
//        p4->back=p2;
//        p4->data=p->children[0];
//        //-----------------
//        p5->next=p6;
//        p5->back=p8;
//        p5->data=p->children[0];
//        //-----------------
//        p6->next=p4;
//        p6->back=p9;
//        p6->data=p->children[0];
//        //-----------------
//        p7->next=NULL;
//        p7->back=p3;
//        p7->data=p->children[1];
//        //-----------------
//        p8->next=NULL;
//        p8->back=p5;
//        p8->data=p->children[0]->children[0];
//        //-----------------
//        p9->next=NULL;
//        p9->back=p6;
//        p9->data=p->children[0]->children[1];
//        //-----------------
//
//        p4->back=p7;
//        p7->back=p4;
//
//        p2->back=p8;
//        p8->back=p2;
//
//        p3->back=p5;
//        p5->back=p3;


        return true;
    }

    bool swap2(PhyTree *p, PhyTree *q) {
        PhyTree *parent_p;        // parent of p
        PhyTree *parent_q;        // parent of t2q
        PhyTree *parent_parent_p; //parent of parent of p
        PhyTree *sister_p; // t1's sister
        unsigned int index_p;
        unsigned int index_sister_p;
        unsigned int index_parent_p;
        unsigned int index_q;

        if(p->parent==NULL){
            return false;
        }
        if(q->parent==NULL){
            return false;
        }
        if(p->parent->parent==NULL){
            return false;
        }

//        if(p->isLeaf() || q->isLeaf()){
//            return false;
//        }

        parent_p=p->parent;
        index_p=p->indexOf();

        index_sister_p=index_p==0?1:0;
        sister_p=parent_p->children[index_sister_p];

        parent_q=q->parent;
        index_q=q->indexOf();

        parent_parent_p=parent_p->parent;
        index_parent_p=parent_p->indexOf();

        if(parent_q==p){
            PhyTree *sister_q;
            PhyTree *child_q;
            unsigned int index_sister_q;
            index_sister_q=index_q==0?1:0;
            sister_q=parent_q->children[index_sister_q];
            child_q=q->children[index_sister_q];
            parent_q->children[index_sister_q]=child_q;
            child_q->parent=p;
            q->children[index_sister_q]=sister_q;
            sister_q->parent=q;
            return true;
        }
        if(parent_p==q){
            PhyTree *child_p;
            unsigned int index_child_p;
            index_child_p=index_sister_p; //by convention
            child_p=p->children[index_child_p];
            parent_p->children[index_sister_p]=child_p;
            child_p->parent=parent_p;
            p->children[index_child_p]=sister_p;
            sister_p->parent=p;
            return true;
        }


        //remove parent_p node and attach p's sister to its grand parent
        parent_parent_p->children[index_parent_p]=sister_p;
        sister_p->parent=parent_parent_p;

        //add parent of p node as q's parent
        parent_q->children[index_q]=parent_p;
        parent_p->parent=parent_q;
        parent_p->children[index_sister_p]=q;
        q->parent=parent_p;

        return true;
    }

    void swap(PhyTree *t1, PhyTree *t2) {
        PhyTree *pt1;        // parent of t1
        PhyTree *pt2;        // parent of t2
        unsigned int index1;    // left/right child index for t1
        unsigned int index2;    // left/right child index for t2

        index1 = t1->indexOf();
        pt1 = t1->parent;

        index2 = t2->indexOf();
        pt2 = t2->parent;

        if (t1 == pt2) {
            std::cout << "CASO1\n";
            pt1->children[index1] = t2;
            t2->parent = pt1;
            t1->children[index2] = t2->children[index2];
            t2->children[index2] = t1;
            t1->parent = t2;
            t1->children[index2]->parent = t1;
        } else if (t2 == pt1) {
            std::cout << "CASO2\n";
            pt2->children[index2] = t1;
            t1->parent = pt2;
            t2->children[index1] = t1->children[index1];
            t1->children[index1] = t2;
            t2->parent = t1;
            t2->children[index1]->parent = t2;
        } else {
            std::cout << "CASO3\n";
            /*
            pt1->children[index1] = t2;
            pt2->children[index2] = t1;
            t1->parent = pt2;
            t2->parent = pt1;
            */
            //PhyTree *ct1;
            PhyTree *ct2;
            //ct1=t1->children[index1];
            ct2 = t2->children[index2];

            pt1->children[index1] = t2;
            t2->parent = pt1;

            t1->parent = t2;
            t2->children[index2] = t1;

            pt2->children[index2] = ct2;
            ct2->parent = pt2;
        }
    }

    void append_MSA_fv(Eigen::VectorXd fv) {
        this->MSA_fv.push_back(fv);
    }

    Eigen::VectorXd &get_MSA_fv(int idx) {
        return this->MSA_fv.at(idx);
    }

    void set_MSA_fv(Eigen::VectorXd &fv) {
        this->MSA_fv.push_back(fv);
    }

    unsigned int get_MSA_fv_size() {
        return this->MSA_fv.size();
    }

    void update_MSA_fv(int dim_alphabet) {
        unsigned int l1, l2;
        PhyTree *n1;
        PhyTree *n2;
        Eigen::VectorXd fv0, fv1, fv2;


        fv0 = Eigen::VectorXd::Zero(dim_alphabet);

        n1 = this->children[0];
        n2 = this->children[1];

        l1 = n1->get_MSA_fv_size();
        l2 = n2->get_MSA_fv_size();

        if (l1 != l2) {
            perror("diff size\n");
            exit(EXIT_FAILURE);
        }

        this->MSA_fv.clear();

        for (unsigned int k = 0; k < l1; k++) {
            fv1 = n1->get_MSA_fv(k);
            fv2 = n2->get_MSA_fv(k);
            //TODO: use cwiseproduct
            for (int i = 0; i < dim_alphabet; i++) {
                fv0[i] = fv1[i] * fv2[i];
            }
            this->MSA_fv.push_back(fv0);
        }

    }

    void set_setA(bool b) {
        this->setA = b;
    }

    bool get_InsertionHistories() {
        return this->setA;
    }

    void set_descCount(int b) {
        this->descCount = b;
    }

    int get_descCount() {
        return this->descCount;
    }

    void set_leaf_character(char ch) {
        this->character = ch;
    }

    char get_leaf_character() {
        return this->character;
    }

    void clear_fv() {

        this->MSA_fv.clear();

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->clear_fv();
        }
    }

    void tmp_initPr(int dim) {

        if (this->Pr.rows() * this->Pr.cols() != 0) {
            this->Pr.resize(0, 0);
        }

        if (this->parent != NULL) {
            this->Pr.resize(dim, dim);
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    this->Pr(i, j) = 0.1 * this->getBranchLength();
                }
            }
//			std::cout<<"NODE="<<this->getName()<<"\n";
//			std::cout<<this->Pr<<"\n";
        }

        for (std::vector<PhyTree *>::iterator i = this->children.begin(); i < this->children.end(); ++i) {
            (*i)->tmp_initPr(dim);
        }

    }

    void set_leaf_state(std::string &MSA_col) {
        int idx;
        idx = 0;
        set_leaf_state(MSA_col, idx);

    }

    void set_leaf_state(std::string &MSA_col, int &idx) {

        if (this->isLeaf()) {
            this->set_leaf_character(MSA_col[idx]);
            idx++;
        } else {
            this->children[0]->set_leaf_state(MSA_col, idx);
            this->children[1]->set_leaf_state(MSA_col, idx);
        }

    }

    void set_ancestral_flag(std::string &MSA_col, int num_gaps) {
        int descCount;
        if (this->isLeaf()) {
            descCount = this->descCount;
            this->set_setA(descCount == num_gaps);
        } else {
            this->children[0]->set_ancestral_flag(MSA_col, num_gaps);
            this->children[1]->set_ancestral_flag(MSA_col, num_gaps);
            descCount = this->descCount;
            this->set_setA(descCount == num_gaps);
        }

    }

    void set_ancestral_flag(std::string &MSA_col) {
        int num_gaps;

        num_gaps = 0;
        // TODO: wrapping up this function separately
        for (int i = 0; i < MSA_col.size(); i++) {
            num_gaps += (MSA_col[i] != '-');
        }

        //std::cout<<"num non gaps="<<num_gaps<<"\n";

//        set_descCount(node);
        this->set_descCount();
        this->set_ancestral_flag(MSA_col, num_gaps);
    }

    void set_descCount() {

        if (this->isLeaf()) {
            this->descCount = this->get_leaf_character() == '-' ? 0 : 1;
            //std::cout<<"set_descCount:"<<node->getName()<<": "<<node->get_descCount()<<"\n";
        } else {
            this->children[0]->set_descCount();
            this->children[1]->set_descCount();
            this->descCount = this->children[0]->descCount + this->children[1]->descCount;
            //std::cout<<"set_descCount:"<<node->getName()<<": "<<node->get_descCount()<<"\n";
        }

    }

    void create_unrooted_tree(std::vector<node *> &utree,PhyTree *tree,node *parent){


        if(tree->isLeaf()){

            node *n = new node;

            parent->back=n;
            n->next=NULL;
            n->back=parent;
            n->data=tree;

            utree.push_back(n);
        }else{

            node *node_1= new node;
            node *node_2= new node;
            node *node_3= new node;

            parent->back=node_1;

            node_1->next=node_2;
            node_1->back=parent;
            node_1->data=tree;

            node_2->next=node_3;
            node_2->data=tree;

            node_3->next=node_1;
            node_3->data=tree;

            utree.push_back(node_1);
            utree.push_back(node_2);
            utree.push_back(node_3);

            create_unrooted_tree(utree,tree->children[0],node_2);
            create_unrooted_tree(utree,tree->children[1],node_3);
        }

    }

    std::vector<node *> create_unrooted_tree(PhyTree *tree,int num_leaves){
        std::vector<node *> utree;

        if(tree->n_children()!=2){
            perror("not binary tree\n");
            exit(EXIT_FAILURE);
        }

        if(num_leaves<3){
            perror("min 3 leaves\n");
            exit(EXIT_FAILURE);
        }

        if(tree->children[0]->n_children()!=2){
            PhyTree *tmp;
            tmp=tree->children[0];
            tree->children[0]=tree->children[1];
            tree->children[1]=tmp;
        }

        node *pseudo_root_1 = new node;
        node *pseudo_root_2 = new node;
        node *pseudo_root_3 = new node;

        pseudo_root_1->next=pseudo_root_2;
        pseudo_root_2->next=pseudo_root_3;
        pseudo_root_3->next=pseudo_root_1;

        pseudo_root_1->data = tree->children[0];
        pseudo_root_2->data = tree->children[0];
        pseudo_root_3->data=tree->children[0];

        utree.push_back(pseudo_root_1);
        utree.push_back(pseudo_root_2);
        utree.push_back(pseudo_root_3);

        create_unrooted_tree(utree,tree->children[1],pseudo_root_1);
        create_unrooted_tree(utree,tree->children[0]->children[0],pseudo_root_2);
        create_unrooted_tree(utree,tree->children[0]->children[1],pseudo_root_3);

        return utree;
    }

    bool utree_swap(node *p, node *q) {

        if(p==q){
            perror("same pointer!");
            exit(EXIT_FAILURE);
        }

        if(p->next==NULL || q->next==NULL){
            return false;
        }

        node *n1;
        node *n2;

        n1=p->next->back;
        n2=p->next->next->back;

        n1->back=n2;
        n2->back=n1;

        p->next->back=q->back;
        q->back->back=p->next;

        p->next->next->back=q;
        q->back=p->next->next;

        return true;

    }

};


PhyTree *midpointRoot(PhyTree *root);

std::vector<std::string> get_tree_order_ancestral(const PhyTree *tree);

std::vector<std::string> get_tree_order(const PhyTree *tree);

void update_fv_values(std::vector<PhyTree *> &p, int alphabet_size);

void set_ancestral_flag(PhyTree *node, std::string &MSA_col, int num_gaps);

void set_ancestral_flag(PhyTree *node, std::string &MSA_col);

void set_leaf_state(PhyTree *node, std::string &MSA_col, int &idx);

void set_leaf_state(PhyTree *node, std::string &MSA_col);

void print_ancestral_flag(PhyTree *node);

void print_leaf_state(PhyTree *node);

void print_descCount(PhyTree *tree);


#endif /* TSHLIB_PHYTREE_HPP */
