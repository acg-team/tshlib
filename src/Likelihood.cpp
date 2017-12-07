/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti and Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti and Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * tshexe is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tshexe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with likpip. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Likelihood.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 10 2017
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
#include <fstream>
#include <random>
#include <glog/logging.h>
#include <gflags/gflags_declare.h>
#include <glog/vlog_is_on.h>
#include "Likelihood.hpp"

namespace LKFunc {
    Eigen::VectorXd
    compute_lk_empty_col(PhyTree &node, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {
        Eigen::VectorXd fv;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;

        if (node.isLeaf()) {
            fv = Eigen::VectorXd::Zero(dim_extended_alphabet);

            fv[dim_extended_alphabet - 1] = 1.0;

            lk += node.get_iota() * (1 - node.get_beta() + node.get_beta() * (fv.dot(pi)));

        } else {

            fvL = compute_lk_empty_col(node[0], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);
            fvR = compute_lk_empty_col(node[1], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);

            fv = (node.get_left_child()->get_Pr() * fvL).cwiseProduct(node.get_right_child()->get_Pr() * fvR);

            lk += node.get_iota() * (1 - node.get_beta() + node.get_beta() * (fv.dot(pi)));

        }
        return fv;
    }




    double compute_log_lk_empty_col(PhyTree &node, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {
        double p0;

        compute_lk_empty_col(node, p0, pi, is_DNA_AA_Codon, dim_extended_alphabet);

        return log(p0);
    }



    Eigen::VectorXd
    compute_lk_recursive(PhyTree &node, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {
        Eigen::VectorXd fv;
        Eigen::VectorXd fvL;
        Eigen::VectorXd fvR;

        if (node.isLeaf()) {
            fv = Eigen::VectorXd::Zero(dim_extended_alphabet);

            int idx;
            // TODO: Maybe this can be antoher function independt.
            if (is_DNA_AA_Codon == 1) {
                idx = mytable[(int) node.get_leaf_character()];
            } else if (is_DNA_AA_Codon == 2) {
                idx = mytableAA[(int) node.get_leaf_character()];
            } else {
                perror("not implemented for codon model yet\n");
                exit(EXIT_FAILURE);
            }

            idx = idx < 0 ? dim_extended_alphabet - 1 : idx;  //TODO: check the alphabet size and the gap index.
            fv[idx] = 1.0;

            if (node.get_InsertionHistories()) {
                //std::cout<<"SETA1 "<<node.getName()<<std::endl;
                lk += node.get_iota() * node.get_beta() * (fv.dot(pi));
            }

            node.set_MSA_fv(fv);

        } else {

            fvL = compute_lk_recursive(node[0], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);
            fvR = compute_lk_recursive(node[1], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);

            fv = (node.get_left_child()->get_Pr() * fvL).cwiseProduct(node.get_right_child()->get_Pr() * fvR);

            if (node.get_InsertionHistories()) {
                //std::cout<<"SETA2 "<<node.getName()<<std::endl;
                lk += node.get_iota() * node.get_beta() * (fv.dot(pi));
            }


            node.set_MSA_fv(fv);


        }
        return fv;
    }



    double compute_col_lk(PhyTree &tree, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int alphabet_size) {

        double lk = 0.0;

        compute_lk_recursive(tree, lk, pi, is_DNA_AA_Codon, alphabet_size);

        return log(lk);
    }



}


double Likelihood::compute_log_lk_empty_col(VirtualNode *root, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {
    //       VirtualNode *pseudo_root;
    double p0;

    //     pseudo_root=utree.startVNodes.at(0);

    compute_lk_empty_col(root, p0, pi, is_DNA_AA_Codon, dim_extended_alphabet);

//        pseudo_root=utree.startVNodes.at(1);

    //      compute_lk_empty_col(pseudo_root, p0, pi, is_DNA_AA_Codon, dim_extended_alphabet);

    return log(p0);
}

void
Likelihood::compute_lk_empty_col(VirtualNode *vnode, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {

    if (vnode->isTerminalNode()) {
        vnode->vnode_Fv_empty = Eigen::VectorXd::Zero(dim_extended_alphabet);

        vnode->vnode_Fv_empty[dim_extended_alphabet - 1] = 1.0;

        lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (vnode->vnode_Fv_empty.dot(pi)));

    } else {

        compute_lk_empty_col(vnode->getNodeLeft(), lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);
        compute_lk_empty_col(vnode->getNodeRight(), lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);

        vnode->vnode_Fv_empty = (vnode->getNodeLeft()->getPr() * vnode->getNodeLeft()->vnode_Fv_empty).cwiseProduct(vnode->getNodeRight()->getPr() * vnode->getNodeRight()
                ->vnode_Fv_empty);

        lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (vnode->vnode_Fv_empty.dot(pi)));

    }
}

Eigen::VectorXd
Likelihood::compute_lk_recursive(VirtualNode *vnode, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet, int colnum) {
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;

    if (vnode->isTerminalNode()) {
        fv = Eigen::VectorXd::Zero(dim_extended_alphabet);

        int idx;
        // TODO: Maybe this can be antoher function independt.
        if (is_DNA_AA_Codon == 1) {
            idx = mytable[(int) vnode->getLeafCharacter()];
        } else if (is_DNA_AA_Codon == 2) {
            idx = mytableAA[(int) vnode->getLeafCharacter()];
        } else {
            perror("not implemented for codon model yet\n");
            exit(EXIT_FAILURE);
        }

        idx = idx < 0 ? dim_extended_alphabet - 1 : idx;  //TODO: check the alphabet size and the gap index.
        fv[idx] = 1.0;

        if (vnode->getSetA(colnum)) {
            //std::cout<<"SETA1* "<<vnode->vnode_name<<std::endl;
            lk += vnode->getIota() * vnode->getBeta() * (fv.dot(pi));
        }
        vnode->vnode_Fv.push_back(fv);
        //vnode->setMSAFv(fv);

    } else {

        fvL = compute_lk_recursive(vnode->getNodeLeft(), lk, pi, is_DNA_AA_Codon, dim_extended_alphabet, colnum);
        fvR = compute_lk_recursive(vnode->getNodeRight(), lk, pi, is_DNA_AA_Codon, dim_extended_alphabet, colnum);

        fv = (vnode->getNodeLeft()->getPr() * fvL).cwiseProduct(vnode->getNodeRight()->getPr() * fvR);

        if (vnode->getSetA(colnum)) {
            //std::cout<<"SETA2* "<<vnode->vnode_name<<std::endl;
            lk += vnode->getIota() * vnode->getBeta() * (fv.dot(pi));
        }

        vnode->vnode_Fv.push_back(fv);
        //vnode->setMSAFv(fv);


    }
    return fv;
}

double Likelihood::compute_col_lk(VirtualNode *root, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int alphabet_size, int colnum) {

    double lk = 0.0;

    compute_lk_recursive(root, lk, pi, is_DNA_AA_Codon, alphabet_size, colnum);

    return log(lk);
}

double Likelihood::compute_nu(double tau, double lambda, double mu) {

    if (fabs(mu) < 1e-8) {
        perror("ERROR in compute_nu: mu too small");
    }

    return lambda * (tau + 1 / mu);
}


double Likelihood::phi(int m, double nu, double p0) {
    double p;
    double log_factorial_m;

    log_factorial_m = 0;
    for (int i = 1; i <= m; i++) {
        log_factorial_m += log(i);
    }

    p = -log_factorial_m + m * log(nu) + (nu * (p0 - 1));

    return p;
}

/*
std::vector<VirtualNode *> UtreeUtils::get_path_from_nodes(VirtualNode *vn1, VirtualNode *vn2) {
    std::vector<VirtualNode *> list_nodes_n0;
    std::vector<VirtualNode *> list_nodes_n1;
    std::vector<VirtualNode *> list_nodes_n2;


    // add nodes from n1 to root
    list_nodes_n1 = fill_with_nodes(vn1);


    // add nodes from n2 to root
    list_nodes_n2 = fill_with_nodes(vn2);

    list_nodes_n0 = get_unique(list_nodes_n1, list_nodes_n2);

    return list_nodes_n0;
}
*/

void Likelihood::recombineAllFv(std::vector<VirtualNode *> list_vnode_to_root){
    VirtualNode *vn;

    for(unsigned int k=0;k<list_vnode_to_root.size();k++){
        vn=list_vnode_to_root.at(k);
        vn->recombineFv();
    }

}

void Likelihood::revertAllFv(std::vector<VirtualNode *> list_vnode_to_root){
    VirtualNode *vn;

    for(unsigned int k=0;k<list_vnode_to_root.size();k++){
        vn=list_vnode_to_root.at(k);
        vn->revertFv();
    }

}

void Likelihood::keepAllFv(std::vector<VirtualNode *> list_vnode_to_root){
    VirtualNode *vn;

    for(unsigned int k=0;k<list_vnode_to_root.size();k++){
        vn=list_vnode_to_root.at(k);
        vn->keepFv();
    }

}

double Likelihood::computePartialLK(std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment, Eigen::VectorXd &pi) {

    double lk = 0;
    for(int alignment_column=0; alignment_column<alignment.getAlignmentSize(); alignment_column++) {
        double lk_col = 0;
        for (auto &vnode:list_vnode_to_root) {
            int currNodeDescCount;
            if (!vnode->isTerminalNode()) {

                currNodeDescCount = vnode->getNodeLeft()->vnode_descCount.at(alignment_column) + vnode->getNodeRight()->vnode_descCount.at(alignment_column);

                if (currNodeDescCount != alignment.align_num_characters.at(alignment_column)) {

                    lk_col += vnode->getIota() * vnode->getBeta() * (vnode->vnode_Fv_temp.at(alignment_column).dot(pi));

                }
            }else{

                lk_col += vnode->getIota() * vnode->getBeta() * (vnode->vnode_Fv_temp.at(alignment_column).dot(pi));
            }

        }

        VLOG(2) << "[tree lk] column: " << alignment_column << " LK: " << lk_col;
        lk = lk+log(lk_col);
    }

    return lk;


}

double Likelihood::computeLogLkEmptyColumnBothSides(VirtualNode *source, VirtualNode *target, Eigen::VectorXd &pi, int m, double nu, int dim_extended_alphabet) {

    double lk;
    double lk_sideA = computeLkEmptyColumn(source, pi, dim_extended_alphabet);
    double lk_sideB = computeLkEmptyColumn(target, pi, dim_extended_alphabet);

    lk = lk_sideA + lk_sideB;
    // TODO: getnodeup is fermi uno prima del nodo root.
    lk = phi(m, nu, lk);

    return lk;
}

double Likelihood::computeLkEmptyColumn(VirtualNode *vnode, Eigen::VectorXd &pi, int dim_extended_alphabet ) {

    Eigen::VectorXd fv;
    double lk = 0;

    if(vnode->getNodeUp()){
        if (vnode->isTerminalNode()) {
            fv = ::Eigen::DenseBase<::Eigen::Matrix<double, -1, 1, 0, -1, 1>>::Zero(dim_extended_alphabet);

            fv[dim_extended_alphabet - 1] = 1.0;

            lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (fv.dot(pi)));

        } else {

            fv = (vnode->getNodeLeft()->getPr() *
                    vnode->getNodeLeft()->vnode_Fv_empty_temp).cwiseProduct(vnode->getNodeRight()->getPr() *
                                                                                                     vnode->getNodeRight()->vnode_Fv_empty_temp);

            lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (fv.dot(pi)));

        }


        computeLkEmptyColumn(vnode->getNodeUp(), pi, dim_extended_alphabet);
    }



    return lk;


}

void Likelihood::recombineEmptyFv(VirtualNode *vnode, Eigen::VectorXd &pi, int dim_extended_alphabet) {

    if(vnode->getNodeUp()){
        if (vnode->isTerminalNode()) {

            vnode->vnode_Fv_empty_temp = ::Eigen::DenseBase<::Eigen::Matrix<double, -1, 1, 0, -1, 1>>::Zero(dim_extended_alphabet);
            vnode->vnode_Fv_empty_temp[dim_extended_alphabet - 1] = 1.0;

        }else{
            // TODO: source target should have to initialise first the children and only after themselves.
            vnode->vnode_Fv_empty_temp =  vnode->getNodeLeft()->vnode_Fv_empty_temp.cwiseProduct(vnode->getNodeRight()->vnode_Fv_empty_temp);

        }
        recombineEmptyFv(vnode->getNodeUp(), pi, dim_extended_alphabet);

    }


}

void Likelihood::recombineAllEmptyFv(VirtualNode *source, VirtualNode *target, Eigen::VectorXd &pi, int dim_extended_alphabet) {

    // TODO: avoid double passing on common nodes
    if(!source->isTerminalNode()){

        recombineEmptyFv(source->getNodeLeft(), pi, dim_extended_alphabet);
        recombineEmptyFv(source->getNodeRight(), pi, dim_extended_alphabet);

    }else{

        recombineEmptyFv(source, pi, dim_extended_alphabet);

    }

    if(!target->isTerminalNode()){

        recombineEmptyFv(target->getNodeLeft(), pi, dim_extended_alphabet);
        recombineEmptyFv(target->getNodeRight(), pi, dim_extended_alphabet);

    }else{

        recombineEmptyFv(target, pi, dim_extended_alphabet);

    }


}

void Likelihood::Init(Utree *tree, Eigen::VectorXd &pi, Eigen::MatrixXd &Q) {

    this->tree = tree;
    this->pi = pi;
    this->Q = Q;

}

Likelihood::Likelihood() = default;
Likelihood::~Likelihood() = default;





double LKFunc::LKcore(Likelihood &lk, std::vector<VirtualNode *> &list_node, Alignment &alignment){
    return lk.computePartialLK(list_node, alignment, lk.pi);

}
