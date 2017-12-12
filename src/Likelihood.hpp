/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
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
 * @file Likelihood.hpp
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
#ifndef TSHEXE_LIKELIHOOD_HPP
#define TSHEXE_LIKELIHOOD_HPP

#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "PhyTree.hpp"
#include "Utree.hpp"
#include "Alignment.hpp"


const char mytable[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, 2, -1, 1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, 0, 0, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, 1,
                           -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0,
                           -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

const char mytableAA[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, 0, 20, 1, 2, 3, 4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14,
                             15, 16, 20, 17, 18, 20, 19, 20, -1, -1, -1, -1, -1, -1, 0, 20, 1, 2, 3,
                             4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14, 15, 16, 20, 17, 18, 20,
                             19, 20, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                             -1, -1, -1, -1, -1, -1, -1, -1, -1};


class Likelihood {

    typedef Eigen::Matrix<double,  5, 5> MatrixExtended;
    typedef Eigen::Matrix<double,  5, 1> VectorExtended;

private:


public:
    //PhyTree *link_node;
    Utree *tree;
    Eigen::VectorXd pi;
    MatrixExtended Q;
    double mu;
    double tau;
    double lambda;
    double nu;
    MatrixExtended V;
    MatrixExtended Vi;
    VectorExtended sigma;

    Likelihood();

     ~Likelihood();


    void Init(Utree *inUtree, Eigen::VectorXd &valPi, Eigen::MatrixXd &valQ, double valMu, double valLambda);

    void computeFV(std::vector<VirtualNode *> &listNodes, Alignment &MSA);

    /*!
     * @brief compute the normalizing Poisson intensity
     * @param tau
     * @param lambda
     * @param mu
     * @return
     */
    void setNu();

    double phi(int m, double p0);

    void compileNodeList_postorder(std::vector<VirtualNode *> &nodelist, VirtualNode *vnode);

    void recombineAllFv(std::vector<VirtualNode *> list_vnode_to_root);

    void revertAllFv(std::vector<VirtualNode *> list_vnode_to_root);

    void keepAllFv(std::vector<VirtualNode *> list_vnode_to_root);

    double computePartialLK(std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment);

    void computePr(std::vector<VirtualNode *> &listNodes, int extended_alphabet_size);

    void loadLikelihoodComponents_Operative();

    void unloadLikelihoodComponents_Operative();

    void setInsertionHistories(std::vector<VirtualNode *> &listNodes, Alignment &MSA);

    void saveLikelihoodComponents();

    void setAllIotas(std::vector<VirtualNode *> &listNodes);

    void setAllBetas(std::vector<VirtualNode *> &listNodes);

//    double computeLogLkEmptyColumnBothSides(VirtualNode *source, VirtualNode *target, Eigen::VectorXd &pi, int m, double nu, int dim_extended_alphabet);
//    double computeLkEmptyColumn(VirtualNode *vnode, Eigen::VectorXd &pi, int dim_extended_alphabet );
//    void recombineEmptyFv(VirtualNode *vnode, Eigen::VectorXd &pi, int dim_extended_alphabet);
//    void recombineAllEmptyFv(VirtualNode *source, VirtualNode *target, Eigen::VectorXd &pi, int dim_extended_alphabet );
//    void _computeLikelihoodComponentsEmptyColumn_recursive(VirtualNode *vnode, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet);
//    double computeLikelihoodComponents_EmptyColumn(VirtualNode *root, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet);
//    Eigen::VectorXd _computeLikelihoodComponents_recursive(VirtualNode *vn, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet, int colnum, Alignment &MSA);
//    double computeLikelihoodComponents(VirtualNode *root, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int alphabet_size, int colnum, Alignment &MSA);

protected:

};



namespace LKFunc {

    double LKcore(Likelihood &lk, std::vector<VirtualNode *> &list_node, Alignment &alignment);

    // Eigen::VectorXd compute_lk_empty_col(PhyTree &node, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet);

    // double compute_log_lk_empty_col(PhyTree &node, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet);

    // Eigen::VectorXd compute_lk_recursive(PhyTree &node, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet);

    // double compute_col_lk(PhyTree &tree, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int alphabet_size);

}
#endif //TSHEXE_LIKELIHOOD_HPP
