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
#include <iomanip>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include "Likelihood.hpp"

namespace LKFunc {
    /*
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


*/
    void ExtendNodeListonSetA(VirtualNode *qnode, std::vector<VirtualNode *> &list_node, int i);
}

/*
double Likelihood::computeLikelihoodComponents_EmptyColumn(VirtualNode *root, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {

    _computeLikelihoodComponentsEmptyColumn_recursive(root, pi, is_DNA_AA_Codon, dim_extended_alphabet);

}

void
Likelihood::_computeLikelihoodComponentsEmptyColumn_recursive(VirtualNode *vnode, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {

    if (vnode->isTerminalNode()) {
        vnode->vnode_Fv_empty_backup = Eigen::VectorXd::Zero(dim_extended_alphabet);

        vnode->vnode_Fv_empty_backup[dim_extended_alphabet - 1] = 1.0;

        //lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (vnode->vnode_Fv_empty_backup.dot(pi)));

    } else {

        _computeLikelihoodComponentsEmptyColumn_recursive(vnode->getNodeLeft(), pi, is_DNA_AA_Codon, dim_extended_alphabet);
        _computeLikelihoodComponentsEmptyColumn_recursive(vnode->getNodeRight(), pi, is_DNA_AA_Codon, dim_extended_alphabet);

        vnode->vnode_Fv_empty_backup = (vnode->getNodeLeft()->getPr() * vnode->getNodeLeft()->vnode_Fv_empty_backup).cwiseProduct(vnode->getNodeRight()->getPr() * vnode->getNodeRight()
                ->vnode_Fv_empty_backup);

        //lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (vnode->vnode_Fv_empty_backup.dot(pi)));

    }
}

Eigen::VectorXd
Likelihood::_computeLikelihoodComponents_recursive(VirtualNode *vnode, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet, int colnum, Alignment &MSA) {
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;

    if (vnode->isTerminalNode()) {
        fv = Eigen::VectorXd::Zero(dim_extended_alphabet);

        int idx;

        // TODO: Maybe this can be antoher function independt.
        if (is_DNA_AA_Codon == 1) {
            idx = mytable[(int) MSA.align_dataset.at(vnode->vnode_seqid)->seq_data.at(colnum)];
        } else if (is_DNA_AA_Codon == 2) {
            idx = mytableAA[(int) MSA.align_dataset.at(vnode->vnode_seqid)->seq_data.at(colnum)];
        } else {
            perror("not implemented for codon model yet\n");
            exit(EXIT_FAILURE);
        }

        idx = idx < 0 ? dim_extended_alphabet - 1 : idx;  //TODO: check the alphabet size and the gap index.
        fv[idx] = 1.0;

        //if (vnode->getSetA(colnum)) {
            //std::cout<<"SETA1* "<<vnode->vnode_name<<std::endl;
        //    lk += vnode->getIota() * vnode->getBeta() * (fv.dot(pi));
        //}
        vnode->vnode_Fv_backup.push_back(fv);
        vnode->vnode_Fv_empty_backup = Eigen::VectorXd::Zero(dim_extended_alphabet);
        vnode->vnode_Fv_empty_backup(dim_extended_alphabet-1) = 1;
        //vnode->setMSAFv(fv);

    } else {

        fvL = _computeLikelihoodComponents_recursive(vnode->getNodeLeft(), pi, is_DNA_AA_Codon, dim_extended_alphabet, colnum, MSA);
        fvR = _computeLikelihoodComponents_recursive(vnode->getNodeRight(), pi, is_DNA_AA_Codon, dim_extended_alphabet, colnum, MSA);

        fv = (vnode->getNodeLeft()->getPr() * fvL).cwiseProduct(vnode->getNodeRight()->getPr() * fvR);

        //if (vnode->getSetA(colnum)) {
            //std::cout<<"SETA2* "<<vnode->vnode_name<<std::endl;
        //    lk += vnode->getIota() * vnode->getBeta() * (fv.dot(pi));
        //}

        vnode->vnode_Fv_backup.push_back(fv);
        //vnode->setMSAFv(fv);


    }
    return fv;
}

double Likelihood::computeLikelihoodComponents(VirtualNode *root, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int alphabet_size, int colnum, Alignment &MSA) {
    _computeLikelihoodComponents_recursive(root, pi, is_DNA_AA_Codon, alphabet_size, colnum, MSA);
}
*/
void Likelihood::setNu() {

    if (fabs(this->mu) < 1e-8) {
        perror("ERROR in setNu: mu too small");
    }

    this->nu = this->lambda * (this->tau + 1 / this->mu);
    //return lambda * (tau + 1 / mu);
}


double Likelihood::phi(int m, double p0) {
    double p;
    double log_factorial_m;

    log_factorial_m = 0;
    for (int i = 1; i <= m; i++) {
        log_factorial_m += log(i);
    }


    p = -log_factorial_m + m * log(this->nu) + (this->nu * (p0 - 1));

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

void Likelihood::recombineAllFv(std::vector<VirtualNode *> list_vnode_to_root) {

    for (auto &vnode:list_vnode_to_root) {
        //std::cout << "[RecombineFV]" << vnode->printNeighbours() << std::endl;
        vnode->recombineFv();
        vnode->_printFV();
    }

}

void Likelihood::revertAllFv(std::vector<VirtualNode *> list_vnode_to_root) {

    for (auto &vnode:list_vnode_to_root) {
        vnode->revertFv();
    }

}

void Likelihood::keepAllFv(std::vector<VirtualNode *> list_vnode_to_root) {

    for (auto &vnode:list_vnode_to_root) {
        vnode->keepFv();
    }


}



double Likelihood::computePartialLK_WholeAlignment(std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment) {
    double lk_log = 0;

    for (int alignment_column = 0; alignment_column < alignment.getAlignmentSize(); alignment_column++) {
        lk_log+=log(computePartialLK(list_vnode_to_root,  alignment, alignment_column));
    }

    double lk_empty = computePartialEmptyLK(list_vnode_to_root,  alignment);

    // compute PHi
    double log_phi_value = phi(alignment.getAlignmentSize(), lk_empty);
    lk_log += log_phi_value;
    return lk_log;
}

double Likelihood::computePartialEmptyLK(std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment) {
    double lk_empty = 0;

    for (auto &vnode:list_vnode_to_root) {

        if (vnode->isTerminalNode()) {
            lk_empty += vnode->getIota() * (1 - vnode->getBeta());

        } else {
            Eigen::VectorXd &Left = vnode->getNodeLeft()->vnode_Fv_empty_operative;
            Eigen::VectorXd &Right = vnode->getNodeRight()->vnode_Fv_empty_operative;
            lk_empty += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * ((Left).cwiseProduct(Right)).dot(this->pi));
        }

    }

    return lk_empty;
}



double Likelihood::computePartialLK(std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment , int alignment_column) {

    double lk_col = 0;

    for (auto &vnode:list_vnode_to_root) {

        if (vnode->vnode_setA_operative.at(alignment_column)) {

            VLOG(3) << "Likelhood for setA of  " << vnode->vnode_name;

            if(vnode->isTerminalNode()){

                lk_col += vnode->getIota() * vnode->getBeta() * vnode->vnode_Fv_terminal.at(alignment_column);

            }else{
                Eigen::VectorXd &Left = vnode->getNodeLeft()->vnode_Fv_operative.at(alignment_column);
                Eigen::VectorXd &Right = vnode->getNodeRight()->vnode_Fv_operative.at(alignment_column);
                lk_col += vnode->getIota() * vnode->getBeta() * ((Left).cwiseProduct(Right)).dot(this->pi);

            }

        }


    }

    return lk_col;
}

/*

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
                    vnode->getNodeLeft()->vnode_Fv_empty_operative).cwiseProduct(vnode->getNodeRight()->getPr() *
                                                                                                     vnode->getNodeRight()->vnode_Fv_empty_operative);

            lk += vnode->getIota() * (1 - vnode->getBeta() + vnode->getBeta() * (fv.dot(pi)));

        }


        computeLkEmptyColumn(vnode->getNodeUp(), pi, dim_extended_alphabet);
    }



    return lk;


}

void Likelihood::recombineEmptyFv(VirtualNode *vnode, Eigen::VectorXd &pi, int dim_extended_alphabet) {

    if(vnode->getNodeUp()){
        if (vnode->isTerminalNode()) {

            vnode->vnode_Fv_empty_operative = ::Eigen::DenseBase<::Eigen::Matrix<double, -1, 1, 0, -1, 1>>::Zero(dim_extended_alphabet);
            vnode->vnode_Fv_empty_operative[dim_extended_alphabet - 1] = 1.0;

        }else{
            // TODO: source target should have to initialise first the children and only after themselves.
            vnode->vnode_Fv_empty_operative =  vnode->getNodeLeft()->vnode_Fv_empty_operative.cwiseProduct(vnode->getNodeRight()->vnode_Fv_empty_operative);

        }
        //switch(vnode->indexOf()
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
*/

void Likelihood::Init(Utree *inUtree, Eigen::VectorXd &valPi, Eigen::MatrixXd &valQ, double valMu, double valLambda) {

    this->tree = inUtree;
    this->pi = valPi;
    this->Q = valQ;

    this->tau = inUtree->computeTotalTreeLength();
    this->mu = valMu;
    this->lambda = valLambda;
    this->setNu();

}


void Likelihood::computePr(std::vector<VirtualNode *> &listNodes, int extended_alphabet_size) {

    // TODO: the following 4 lines have to be computed only when Q changes!!!
    Eigen::EigenSolver<MatrixExtended> solver(this->Q);
    this->sigma = solver.eigenvalues().real();
    this->V = solver.eigenvectors().real();
    this->Vi = this->V.inverse();

    for (auto &vnode:listNodes) {

        // Clear the matrix containing the probability matrix
        vnode->vnode_Pr.resize(0, 0);

        // Recompute the Q matrix with the exponential equal to the branch lenght
        vnode->vnode_Pr.resize(extended_alphabet_size, extended_alphabet_size);

        vnode->vnode_Pr = this->V * (this->sigma * vnode->vnode_branchlength).array().exp().matrix().asDiagonal() * this->Vi;

    }
}

void Likelihood::compileNodeList_postorder(std::vector<VirtualNode *> &nodelist, VirtualNode *vnode) {

    if (!vnode->isTerminalNode()) {

        compileNodeList_postorder(nodelist, vnode->getNodeLeft());
        compileNodeList_postorder(nodelist, vnode->getNodeRight());
        nodelist.push_back(vnode);

    } else {

        nodelist.push_back(vnode);
    }


}

void Likelihood::loadLikelihoodComponents_Operative() {

    for (auto &vnode:this->tree->listVNodes) {
        std::swap(vnode->vnode_Fv_backup, vnode->vnode_Fv_operative);
        std::swap(vnode->vnode_Fv_empty_backup, vnode->vnode_Fv_empty_operative);
        std::swap(vnode->vnode_descCount_backup, vnode->vnode_descCount_operative);

    }
    if (this->tree->rootnode) {
        std::swap(this->tree->rootnode->vnode_Fv_backup, this->tree->rootnode->vnode_Fv_operative);
        std::swap(this->tree->rootnode->vnode_Fv_empty_backup, this->tree->rootnode->vnode_Fv_empty_operative);
        std::swap(this->tree->rootnode->vnode_descCount_backup, this->tree->rootnode->vnode_descCount_operative);

    }

}

void Likelihood::unloadLikelihoodComponents_Operative() {

    for (auto &vnode:this->tree->listVNodes) {
        std::swap(vnode->vnode_Fv_operative, vnode->vnode_Fv_backup);
        std::swap(vnode->vnode_Fv_empty_operative, vnode->vnode_Fv_empty_backup);
        std::swap(vnode->vnode_descCount_operative, vnode->vnode_descCount_backup);
    }
    if (this->tree->rootnode) {
        std::swap(this->tree->rootnode->vnode_Fv_operative, this->tree->rootnode->vnode_Fv_backup);
        std::swap(this->tree->rootnode->vnode_Fv_empty_operative, this->tree->rootnode->vnode_Fv_empty_backup);
        std::swap(this->tree->rootnode->vnode_descCount_operative, this->tree->rootnode->vnode_descCount_backup);

    }
}

void Likelihood::setInsertionHistories(std::vector<VirtualNode *> &listNodes, Alignment &MSA) {
    //std::cout<<"\t";
    //for (auto &vnode:listNodes) {

    //    std::cout << "\t"<< vnode->vnode_name;
    //}

    //std::cout << std::endl;

    for (int i = 0; i < MSA.getAlignmentSize(); i++) {
        //std::cout << "["<<std::setfill('0') << std::setw(2)<<i<< "]\t";

        for (auto &vnode:listNodes) {

            if (vnode->isTerminalNode()) {

                vnode->vnode_descCount_operative.at(i) = (MSA.align_dataset.at(vnode->vnode_seqid)->seq_data.at(i) == '-' ? 0 : 1);

            } else {

                vnode->vnode_descCount_operative.at(i) = vnode->getNodeLeft()->vnode_descCount_operative.at(i) +
                                                         vnode->getNodeRight()->vnode_descCount_operative.at(i);

            }

            vnode->vnode_setA_operative.at(i) = (vnode->vnode_descCount_operative.at(i) == MSA.align_num_characters.at(i));
            //  std::cout << vnode->vnode_setA_operative.at(i) << "\t";
        }
        //std::cout << std::endl;
    }

}

void Likelihood::computeFV(std::vector<VirtualNode *> &listNodes, Alignment &MSA) {

    for (auto &vnode:listNodes) {
        for (int i = 0; i < MSA.getAlignmentSize(); i++) {
            if (vnode->isTerminalNode()) {
                Eigen::VectorXd fv = Eigen::VectorXd::Zero(MSA.align_alphabetsize);

                int idx;

                // TODO: Maybe this can be antoher function independt.
                if (MSA.alphabet == AlignmentAlphabet::dna) {
                    idx = mytable[(int) MSA.align_dataset.at(vnode->vnode_seqid)->seq_data.at(i)];
                } else if (MSA.alphabet == AlignmentAlphabet::aa) {
                    idx = mytableAA[(int) MSA.align_dataset.at(vnode->vnode_seqid)->seq_data.at(i)];
                } else {
                    perror("not implemented for codon model yet\n");
                    exit(EXIT_FAILURE);
                }

                //TODO: check the alphabet size and the gap index.
                idx = idx < 0 ? MSA.align_alphabetsize - 1 : idx;
                fv[idx] = 1.0;


                vnode->vnode_Fv_terminal.at(i) = fv.dot(this->pi);
                vnode->vnode_Fv_operative.at(i) = vnode->getPr() * fv;


            } else {

                Eigen::VectorXd &fvL = vnode->getNodeLeft()->vnode_Fv_operative.at(i);
                Eigen::VectorXd &fvR = vnode->getNodeRight()->vnode_Fv_operative.at(i);

                if (!vnode->isRootNode()) {

                    vnode->vnode_Fv_operative.at(i) = vnode->getPr() * ((fvL).cwiseProduct(fvR));

                } else {

                    vnode->vnode_Fv_operative.at(i) = (fvL).cwiseProduct(fvR);
                }

            }

        }

        if (vnode->isTerminalNode()) {

            Eigen::VectorXd fv = Eigen::VectorXd::Zero(MSA.align_alphabetsize);
            fv(MSA.align_alphabetsize - 1) = 1.0;

            vnode->vnode_Fv_empty_operative = vnode->getPr() * fv;

        } else {

            Eigen::VectorXd &fvE_L = vnode->getNodeLeft()->vnode_Fv_empty_operative;
            Eigen::VectorXd &fvE_R = vnode->getNodeRight()->vnode_Fv_empty_operative;

            if (!vnode->isRootNode()) {

                vnode->vnode_Fv_empty_operative = vnode->getPr() * ((fvE_L).cwiseProduct(fvE_R));

            } else {

                vnode->vnode_Fv_empty_operative = (fvE_L).cwiseProduct(fvE_R);
            }

        }

        vnode->_printFV();
    }

}

void Likelihood::saveLikelihoodComponents() {


    for (auto &vnode:this->tree->listVNodes) {

//        std::copy(vnode->vnode_Fv_operative.begin(),
//                  vnode->vnode_Fv_operative.end(),
//                  std::back_inserter(vnode->vnode_Fv_backup));

        vnode->vnode_Fv_backup = vnode->vnode_Fv_operative;

        vnode->vnode_Fv_empty_backup = vnode->vnode_Fv_empty_operative;

        vnode->vnode_descCount_backup = vnode->vnode_descCount_operative;

//        std::copy(vnode->vnode_descCount_operative.begin(),
//                  vnode->vnode_descCount_operative.end(),
//                  std::back_inserter(vnode->vnode_descCount_backup));

    }
    if (this->tree->rootnode) {

//        std::copy(this->tree->rootnode->vnode_Fv_operative.begin(),
//                  this->tree->rootnode->vnode_Fv_operative.end(),
//                  std::back_inserter(this->tree->rootnode->vnode_Fv_backup));


        this->tree->rootnode->vnode_Fv_backup = this->tree->rootnode->vnode_Fv_operative;
        this->tree->rootnode->vnode_Fv_empty_backup = this->tree->rootnode->vnode_Fv_empty_operative;
        this->tree->rootnode->vnode_descCount_backup = this->tree->rootnode->vnode_descCount_operative;
//        std::copy(this->tree->rootnode->vnode_descCount_operative.begin(),
//                  this->tree->rootnode->vnode_descCount_operative.end(),
//                  std::back_inserter(this->tree->rootnode->vnode_descCount_backup));


    }


}


void Likelihood::setAllIotas(std::vector<VirtualNode *> &listNodes) {
    double T;

    if (fabs(this->mu) < 1e-8) {
        perror("ERROR in set_iota: mu too small");
    }

    T = this->tau + 1 / this->mu;

    if (fabs(T) < 1e-8) {
        perror("ERROR in set_iota: T too small");
    }


    for (auto &vnode:listNodes) {

        if (vnode->isRootNode()) {

            vnode->vnode_iota = (1 / this->mu) / T;

        } else {

            vnode->vnode_iota = vnode->vnode_branchlength / T;
        }

    }
}


void Likelihood::setAllBetas(std::vector<VirtualNode *> &listNodes) {

    if (fabs(this->mu) < 1e-8) {
        perror("ERROR in set_iota: mu too small");
    }


    for (auto &vnode:listNodes) {

        if (vnode->isRootNode()) {

            vnode->vnode_beta = 1.0;

        } else {

            vnode->vnode_beta = (1.0 - exp(-this->mu * vnode->vnode_branchlength)) / (this->mu * vnode->vnode_branchlength);
        }

    }


}

void Likelihood::optimiseLambda(int m, double p0) {


    this->lambda = exp(phi(m - 1, p0)) / exp(phi(m, p0)) - 1;

}



Likelihood::Likelihood() = default;

Likelihood::~Likelihood() = default;


double LKFunc::LKcore(Likelihood &lk, std::vector<VirtualNode *> &list_node, Alignment &alignment) {
    return lk.computePartialLK_WholeAlignment(list_node, alignment);

}

double LKFunc::LKRearrangment(Likelihood &lk, std::vector<VirtualNode *> &list_node_complete, Alignment &alignment) {

    std::vector<VirtualNode *> temp_list;

    double lk_log = 0;

    double lk_empty = lk.computePartialEmptyLK(list_node_complete, alignment);


    for(int i=0; i<alignment.getAlignmentSize(); i++){
        ExtendNodeListonSetA(list_node_complete.back(), temp_list, i);

        lk_log += log(lk.computePartialLK(temp_list, alignment, i));

        temp_list.clear();
    }

    // compute PHi
    double log_phi_value = lk.phi(alignment.getAlignmentSize(), lk_empty);
    lk_log += log_phi_value;
    return lk_log;

}

void LKFunc::ExtendNodeListonSetA(VirtualNode *qnode, std::vector<VirtualNode *> &list_node, int i) {
    VirtualNode *temp = qnode;

    list_node.push_back(temp);

    do{
            if (temp->isTerminalNode()){
                break;
            }

            if(temp->getNodeLeft()->vnode_setA_operative.at(i)){

                temp = temp->getNodeLeft();

            }else if(temp->getNodeRight()->vnode_setA_operative.at(i)){

                temp = temp->getNodeRight();

            }else{

                break;
            }

            list_node.push_back(temp);


    }while(temp->vnode_setA_operative.at(i));
}

