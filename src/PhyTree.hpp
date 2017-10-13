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
 * @file PhyTree.hpp
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

#ifndef TSHLIB_PHYTREE_HPP
#define TSHLIB_PHYTREE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <sstream>
#include <cmath>
#include "../main.hpp"
//=====================
//DP-PIP
#include <Eigen/Core>
#include <Eigen/src/Core/IO.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
//=====================

class PhyTree {

private:
    std::string name;
    double branch_length;
    double branch_support;
    PhyTree *parent;
    std::vector<PhyTree*> children;
    //================================================================
    //DP-PIP
    double iota;
    double beta;
    Eigen::MatrixXd Pr;
    double tau;
    double nu;
    std::vector<Eigen::VectorXd> MSA_fv;
    bool setA;
    char character;
    //================================================================

    void print_prefix(std::string prefix) const {

        //std::cout << prefix << "(" << this->branch_length << ") " << this->name << std::endl;

        //====================================================================================
        //DP-PIP
        std::cout << prefix << "(" << this->branch_length << ") " << this->name << " " <<" iota: " << this->iota <<" beta: "<<this->beta<<std::endl;
        //====================================================================================

        for(std::vector<PhyTree*>::const_iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->print_prefix(prefix+"  ");
        }
    }

//	void fixDistancesR() {
////		if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
////			if(std::isnan(this->branch_length)) this->branch_length = cmdlineopts.max_dist;
////			this->branch_length = std::min(std::max(cmdlineopts.min_dist,this->branch_length),cmdlineopts.max_dist);
////		} else {
//			if(std::isnan(this->branch_length)) this->branch_length = cmdlineopts.max_pdist;
//			this->branch_length = std::min(std::max(cmdlineopts.min_pdist,this->branch_length),cmdlineopts.max_pdist);
////		}
//		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
//			(*i)->fixDistancesR();
//		}
//	}

    std::string formatNewickR() const {
        if(this->isLeaf()) {
            return this->getName();
        } else {
            std::stringstream newick;
            newick << "(";
            std::vector<PhyTree*>::const_iterator i=this->children.begin();
            newick << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
            for(++i; i < this->children.end(); ++i) {
                newick << "," << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
            }
            newick << ")";
            return newick.str();
        }
    }

public:

    //================================================================
    //DP-PIP
    //ugly public
    //double **fv_array;
    //================================================================

    PhyTree(std::string name="") {
        this->parent = NULL;
        this->branch_length = 0;
        this->branch_support = 1;
        this->name = name;

        //==============================
        //DP-PIP
        this->iota = 0;
        this->beta = 0;
        this->Pr.resize(0,0);
        this->tau=0;
        this->nu=0;
        //==============================

    }

    ~PhyTree() {
        assert(this->parent == NULL);
        for(std::vector<PhyTree*>::reverse_iterator i=this->children.rbegin(); i < this->children.rend(); ++i) {
            PhyTree *child = *i;
            child->parent = NULL;
            child->branch_length = 0;

            delete child;
        }
    }

    PhyTree* copy() {
        PhyTree* out = new PhyTree();
        out->branch_length = this->branch_length;
        out->branch_support = this->branch_support;
        out->name = this->name;

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            out->addChild((*i)->copy(),(*i)->branch_length,(*i)->branch_support);
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

    index_t indexOf() {
        PhyTree *parent = this->parent;
        assert(parent != NULL);

        for(index_t i=0; i < parent->children.size(); ++i) {
            if(parent->children[i] == this) {
                return i;
            }
        }

        assert(false);
        return -1;
    }

    void pluck() {
        assert(this->parent != NULL);

        index_t index = this->indexOf();
        std::vector<PhyTree*>::iterator iter = this->parent->children.begin()+index;

        this->parent->children.erase(iter);
        this->parent = NULL;

        this->branch_length = 0;
        this->branch_support = 1;
    }

    PhyTree* pluckChild(index_t index) {
        std::vector<PhyTree*>::iterator iter = this->children.begin()+index;
        PhyTree *child = *iter;

        this->children.erase(iter);
        child->parent = NULL;
        child->branch_length = 0;
        this->branch_support = 1;

        return child;
    }

    void deleteChild(index_t index) {
        PhyTree *child = this->pluckChild(index);
        delete child;
    }

//	void fixDistances() {
//		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
//			(*i)->fixDistancesR();
//		}
//	}

    index_t countLeaves() {
        if(this->isLeaf()) {
            return 1;
        } else {
            index_t leaves = 0;
            for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
                leaves += (*i)->countLeaves();
            }
            return leaves;
        }
    }

    double computeLength() {
        if(this->isLeaf()) {
            return 0;
        } else {
            double length = 0;
            for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
                length += (*i)->branch_length + (*i)->computeLength();
            }
            return length;
        }
    }

    std::string getName() const {
        return this->name;
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

    PhyTree& operator[](int i) {
        return *this->children[i];
    }

    const PhyTree& operator[](int i) const {
        return *this->children[i];
    }

    std::vector<PhyTree*>::iterator firstChild() {
        return this->children.begin();
    }

    std::vector<PhyTree*>::iterator lastChild() {
        return this->children.end();
    }

    std::vector<PhyTree*>::const_iterator firstChild() const {
        return this->children.begin();
    }

    std::vector<PhyTree*>::const_iterator lastChild() const {
        return this->children.end();
    }

    index_t n_children() const {
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

    //=======================================================================================================
    //DP-PIP
    std::vector<PhyTree*> get_children(){

        return this->children;
    }
    //=======================================================================================================
    //DP-PIP
    PhyTree* get_left_child(){

        if(this->n_children()>0){
            return this->children[0];
        }else{
            return NULL;
        }

    }
    //=======================================================================================================
    //DP-PIP
    PhyTree* get_right_child(){

        if(this->n_children()>0){
            return this->children[1];
        }else{
            return NULL;
        }

    }
    //=======================================================================================================
    //DP-PIP
    double get_iota() const {
        return this->iota;
    }
    //=======================================================================================================
    //DP-PIP
    double get_beta() const {
        return this->beta;
    }
    //=======================================================================================================
    //DP-PIP
    //unsigned int get_Pc_size() const {
    //	return this->Pc.size();
    //}
    //=======================================================================================================
    //DP-PIP
    static bool empty(Eigen::VectorXd &v){

        if(v.rows() * v.cols() == 0){
            return true;
        }else{
            return false;
        }
    }
    //=======================================================================================================
    /*
    //DP-PIP
    double get_Pc(int index){

        if(!empty(this->Pc)){
            return (double)this->Pc(index);
        }else{
            return 0.0;
        }

    }
    */
    //=======================================================================================================
    //DP-PIP
    void setName(std::string name){
        this->name=name;
    }
    //=======================================================================================================
    //DP-PIP
    void set_missing_node_name(std::string s){

        static int num=0;

        if(this->isLeaf()){

        }else{

            if(this->name.empty()){
                std::string name(s);
                name.append(std::to_string(num++));
                this->setName(name);
            }

            //std::cout<<"name: "<<this->name<<"\n";

        }

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_missing_node_name(s);
        }

    }
    //=======================================================================================================
    /*
    //DP-PIP
    void set_Pc(int index,double val){

        if(!empty(this->Pc)){
            this->Pc[index]=val;
        }

    }
    */
    //=======================================================================================================
    /*
    //DP-PIP
    void setZero_Pc(int len){

        this->Pc.setZero(len);

    }
    */
    //=======================================================================================================
    //DP-PIP
    const Eigen::MatrixXd& get_Pr(){

        return this->Pr;

    }
    //=======================================================================================================
    //DP-PIP
    void set_iota(double tau,double mu){

        if(fabs(mu)<1e-8){
//			error("ERROR in set_iota: mu too small");
        }

        double T=tau+1/mu;

        if(fabs(T)<1e-8){
//			error("ERROR in set_iota: T too small");
        }

        if(this->parent==NULL){
            this->iota=(1/mu)/T;
        }else{
            this->iota=this->branch_length/T;
        }

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_iota(tau,mu);
        }

    }
    //=======================================================================================================
    //DP-PIP
    double get_nu_local() const {
        return this->nu;
    }
    //=======================================================================================================
    //DP-PIP
    double get_tau_local() const {
        return this->tau;
    }
    //=======================================================================================================
    //DP-PIP
    void compute_nu_local(double tau,double lambda,double mu){
        this->nu=lambda*(tau+1/mu);
    }
    //=======================================================================================================
    //DP-PIP
    void compute_tau_local(){

        if(this->isLeaf()) {
            this->tau=this->branch_length;
        } else {
            this->get_left_child()->compute_tau_local();
            this->get_right_child()->compute_tau_local();

            this->tau=this->get_left_child()->tau+this->get_right_child()->tau;
        }

    }
    //=======================================================================================================
    //DP-PIP
    void set_iota_local_down(double tau,double mu){

        this->iota=this->branch_length/(tau+1/mu);

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_iota_local_down(tau,mu);
        }

    }
    //=======================================================================================================
    //DP-PIP
    void set_iota_local(double tau,double mu){

        this->iota=(1/mu)/(tau+1/mu);

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_iota_local_down(tau,mu);
        }


    }
    //=======================================================================================================
    //DP-PIP
    void set_beta(double tau,double mu){

        if(fabs(mu)<1e-8){
//			error("ERROR in set_beta: mu too small");
        }

        if(this->parent==NULL){
            this->beta=1;
        }else{

            if(fabs(this->branch_length)<1e-8){
//				error("ERROR in set_beta: branch_length too small");
            }

            this->beta=(1-exp(-mu*this->branch_length))/(mu*this->branch_length);
        }

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_beta(tau,mu);
        }

    }
    //=======================================================================================================
    //DP-PIP
    void set_beta_down(double tau,double mu){

        if(fabs(mu)<1e-8){
//			error("ERROR in set_beta: mu too small");
        }

        if(fabs(this->branch_length)<1e-8){
//			error("ERROR in set_beta: branch_length too small");
        }

        this->beta=(1-exp(-mu*this->branch_length))/(mu*this->branch_length);

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_beta_down(tau,mu);
        }

    }
    //=======================================================================================================
    //DP-PIP
    void set_beta_local(double tau,double mu){

        this->beta=1.0;

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_beta_down(tau,mu);
        }

    }
    //=======================================================================================================
    //DP-PIP
    void set_tau(double tau){

        this->tau=tau;

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_tau(tau);
        }

    }
    //=======================================================================================================
    //DP-PIP
    void set_nu(double nu){

        this->nu=nu;

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->set_nu(nu);
        }

    }
    //=======================================================================================================
    //DP-PIP
    void print_local_var(){

        std::cout<<"----------------------\n";
        std::cout<<"Name: "<<this->name<<"\n";
        std::cout<<"tau: "<<this->tau<<"\n";
        std::cout<<"nu: "<<this->nu<<"\n";
        std::cout<<"iota: "<<this->iota<<"\n";
        std::cout<<"beta: "<<this->beta<<"\n";

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->print_local_var();
        }
    }
    //=======================================================================================================
    //DP-PIP
    void print_br(){

        std::cout<<"Name: "<<this->name<<" bl: "<<this->branch_length<<std::endl;

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->print_br();
        }
    }
    //=======================================================================================================
    //DP-PIP
//	void print_Pr(){
//
//		std::cout<<"node name: "<<this->name<<std::endl;
//		std::cout<<"Pr, row: "<<this->Pr.rows()<<" columns: "<<this->Pr.cols()<<" size: "<<this->Pr.size()<<"\n";
//		Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
//		std::cout << this->Pr.format(HeavyFmt) << std::endl;
//
//		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
//			(*i)->print_Pr();
//		}
//	}
    //=======================================================================================================
    //DP-PIP
    void printOnlyName(){

        /*
        if(this->isLeaf()){
            std::cout<<"Leaf name: "<<this->name<<std::endl;
        }else{
            std::cout<<"Node name: "<<this->name<<std::endl;
            this->children[0]->printOnlyName();
            this->children[1]->printOnlyName();
        }
        */

        std::cout<<"node name: "<<this->name<<"\n";
        if(this->n_children()==2){
            std::cout<<"      L child: "<<this->get_left_child()->getName()<<"\n";
            std::cout<<"      R child: "<<this->get_right_child()->getName()<<"\n";
        }



        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->printOnlyName();
        }


    }
    //=======================================================================================================
    //DP-PIP
//	template <class ALPHABET>
//	void initPrPIP(const ModelFactory<ALPHABET> *model_factory,double gamma_rate){
//
//		// TODO: check parsedistance
////		Model<ALPHABET> model = model_factory->getPIPModel(this->branch_length);
//
//
////		double tmp=this->branch_length*gamma_rate;
////
////		printf("initPrPIP %lf\n",tmp);
//
//		//@gamma_distribution
//		Model<ALPHABET> model = model_factory->getPIPModel(this->branch_length*gamma_rate);
//
//		if(this->Pr.rows()*this->Pr.cols()!=0){
//			this->Pr.resize(0,0);
//		}
//		if(this->parent!=NULL){
//			this->Pr=model.P_PIP;
//		}
//
////		std::cout<<"\n";
////		std::cout<<"Pr at ("<<this->getName()<<") with bl: "<<this->branch_length<<"\n";
////		std::cout<<this->Pr<<std::endl;
////		std::cout<<"\n";
//
//
//		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
//			(*i)->initPrPIP<ALPHABET>(model_factory,gamma_rate);
//		}
//
//	}
    //=======================================================================================================
    /*
    //DP-PIP
    template <class ALPHABET>
    void init_fv_map(const std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &sequences){

        if(this->isLeaf()){
            const std::string name = this->getName();

            typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>> ::const_iterator vect_iter;
            vect_iter it = std::find_if(sequences.begin(),sequences.end(),CompareFirst<ALPHABET>(name));

            if (it == sequences.end()){
                error("ERROR sequence name doesn't match any tree leaves");
                exit(EXIT_FAILURE);
            }

            std::string s=stringFromSequence<ALPHABET>(it->second);
            sequence_t<ALPHABET> seq=it->second;

            for(unsigned int i=0;i<seq.length();i++){
                Eigen::VectorXd x=Eigen::VectorXd::Zero(ALPHABET::DIM+1);
                x[seq[i].value()]=1.0;
                this->fv_map[s.substr(i,1)]=x;
                //this->fv_map[seq.substr(i,1)]=x;
            }

            Eigen::VectorXd x=Eigen::VectorXd::Zero(ALPHABET::DIM+1);
            x[ALPHABET::DIM]=1.0;

            this->fv_map[GAP_STR]=x;


        }else{

            this->children[0]->init_fv_map<ALPHABET>(sequences);
            this->children[1]->init_fv_map<ALPHABET>(sequences);

        }
    }
    */
    //=======================================================================================================
    //DP-PIP
    /*
    template <class ALPHABET>
    void update_fv_map(std::map<std::string,Eigen::VectorXd> &temp_fv_map,std::vector<std::pair<std::string,sequence_t<ALPHABET> > > MSA){

        typedef typename std::map<std::string,Eigen::VectorXd>::iterator MapIterator;

        unsigned int len=MSA.at(0).second.length();

        for(unsigned int i=0;i<len;i++){
            sequence_t<ALPHABET> seq=create_col_MSA(MSA,i);
            std::string s=stringFromSequence<ALPHABET>(seq);
            MapIterator it=temp_fv_map.find(s);
            if(it == temp_fv_map.end()){
                std::cout<<"ERRORE in update_fv_map: "<<s<<"\n";
            }
            this->fv_map[s]=it->second;
        }

        sequence_t<ALPHABET> seq_gap=create_col_MSA_gap(MSA);
        std::string s_gap=stringFromSequence<ALPHABET>(seq_gap);
        MapIterator it=temp_fv_map.find(s_gap);
        if(it == temp_fv_map.end()){
            std::cout<<"ERRORE in update_fv_map: "<<s_gap<<"\n";
        }
        this->fv_map[s_gap]=it->second;

    }
    */
    //=======================================================================================================
    /*
    //DP-PIP
    template <class ALPHABET>
    void compute_fv(Eigen::VectorXd &fv,sequence_t<ALPHABET> sL,sequence_t<ALPHABET> sR){

        typedef typename std::map<std::string,Eigen::VectorXd>::iterator MapIterator;

        MapIterator itL=this->children[0]->fv_map.find(stringFromSequence<ALPHABET>(sL));
        if(itL == this->children[0]->fv_map.end()){
            std::cout<<"ERRORE in new_get_fv: "<<stringFromSequence<ALPHABET>(sL);
            exit(EXIT_FAILURE);
        }
        Eigen::VectorXd &fvL=itL->second;

        MapIterator itR=this->children[1]->fv_map.find(stringFromSequence<ALPHABET>(sR));
        if(itR == this->children[1]->fv_map.end()){
            std::cout<<"ERRORE in new_get_fv: "<<stringFromSequence<ALPHABET>(sR);
            exit(EXIT_FAILURE);
        }
        Eigen::VectorXd &fvR=itR->second;

        fv=(this->children[0]->Pr*fvL).cwiseProduct(this->children[1]->Pr*fvR);

    }
    */
    //=======================================================================================================
    void null_parent(){
        this->parent=NULL;
    }
    //=======================================================================================================
    void swap(PhyTree *t1,index_t index1,PhyTree *t2,index_t index2){
        PhyTree *t0;

        t0=t1->children[index1];
        t1->children[index1]=t2->children[index2];
        t2->children[index2]=t0;

        t1->children[index1]->parent=t1;
        t2->children[index2]->parent=t2;

    }
    //=======================================================================================================
    void swap2(PhyTree *t1,PhyTree *t2){
        PhyTree *pt1; 		// parent of t1
        PhyTree *pt2; 		// parent of t2
        index_t index1; 	// left/right child index for t1
        index_t index2; 	// left/right child index for t2

        index1=t1->indexOf();
        pt1=t1->parent;

        index2=t2->indexOf();
        pt2=t2->parent;

        pt1->children[index1]=t2;
        pt2->children[index2]=t1;

        t1->parent=pt2;
        t2->parent=pt1;
    }
    //=======================================================================================================
    void append_MSA_fv(Eigen::VectorXd fv){
        this->MSA_fv.push_back(fv);
    }
    //=======================================================================================================
    Eigen::VectorXd &get_MSA_fv(int idx){
        return this->MSA_fv.at(idx);
    }
    //=======================================================================================================
    void set_MSA_fv(Eigen::VectorXd &fv,int idx){
        this->MSA_fv.at(idx)=fv;
    }
    //=======================================================================================================
    unsigned int get_MSA_fv_size(){
        return this->MSA_fv.size();
    }
    //=======================================================================================================
    void update_MSA_fv(int dim_alphabet){
        unsigned int l1,l2;
        PhyTree *n1;
        PhyTree *n2;
        Eigen::VectorXd fv0,fv1,fv2;


        fv0=Eigen::VectorXd::Zero(dim_alphabet);

        n1=this->children[0];
        n2=this->children[1];

        l1=n1->get_MSA_fv_size();
        l2=n2->get_MSA_fv_size();

        if(l1!=l2){
            perror("diff size\n");
            exit(EXIT_FAILURE);
        }

        this->MSA_fv.clear();

        for(unsigned int k=0;k<l1;k++){
            fv1=n1->get_MSA_fv(k);
            fv2=n2->get_MSA_fv(k);
            for(int i=0;i<dim_alphabet;i++){
                fv0[i]=fv1[i]*fv2[i];
            }
            this->MSA_fv.push_back(fv0);
        }

    }
    //=======================================================================================================
    void set_setA(bool b){
        this->setA=b;
    }
    //=======================================================================================================
    bool get_setA(){
        return this->setA;
    }
    //=======================================================================================================
    void set_leaf_character(char ch){
        this->character=ch;
    }
    //=======================================================================================================
    char get_leaf_character(){
        return this->character;
    }
    //=======================================================================================================
    void tmp_initPr(int dim){

        if(this->Pr.rows()*this->Pr.cols()!=0){
            this->Pr.resize(0,0);
        }

        if(this->parent!=NULL){
            this->Pr.resize(dim,dim);
            for(int i=0;i<dim;i++){
                for(int j=0;j<dim;j++){
                    this->Pr(i,j)=0.1;
                }
            }
//			std::cout<<"NODE="<<this->getName()<<"\n";
//			std::cout<<this->Pr<<"\n";
        }

        for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
            (*i)->tmp_initPr(dim);
        }

    }
    //=======================================================================================================
};


PhyTree* midpointRoot(PhyTree *root);

std::vector<std::string> get_tree_order_ancestral(const PhyTree *tree);
std::vector<std::string> get_tree_order(const PhyTree *tree);

void update_fv_values(std::vector<PhyTree *> &p,int alphabet_size);
bool set_ancestral_flag(PhyTree *node,std::string &MSA_col,int &idx);
void set_ancestral_flag(PhyTree *node,std::string &MSA_col);
void set_leaf_state(PhyTree *node,std::string &MSA_col,int &idx);
void set_leaf_state(PhyTree *node,std::string &MSA_col);


#endif /* TSHLIB_PHYTREE_HPP */
