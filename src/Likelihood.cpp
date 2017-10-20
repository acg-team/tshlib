/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti
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
#include "Likelihood.hpp"

//=======================================================================================================
//DP-PIP
Eigen::VectorXd
compute_lk_empty_col(PhyTree &node, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;

    //std::cout<<"node="<<node.getName()<<" fv:"<<fv<<"\n";

    if (node.isLeaf()) {
        fv = Eigen::VectorXd::Zero(dim_extended_alphabet);

        fv[dim_extended_alphabet - 1] = 1.0;

        lk += node.get_iota() * (1 - node.get_beta() + node.get_beta() * (fv.dot(pi)));

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_InsertionHistories()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";*/

        return fv;
    } else {

        fvL = compute_lk_empty_col(node[0], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);
        fvR = compute_lk_empty_col(node[1], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);

        fv = (node.get_left_child()->get_Pr() * fvL).cwiseProduct(node.get_right_child()->get_Pr() * fvR);

        lk += node.get_iota() * (1 - node.get_beta() + node.get_beta() * (fv.dot(pi)));

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_InsertionHistories()<<"\n";
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
Eigen::VectorXd
compute_lk_recursive(PhyTree &node, double &lk, Eigen::VectorXd &pi, int is_DNA_AA_Codon, int dim_extended_alphabet) {
    Eigen::VectorXd fv;
    Eigen::VectorXd fvL;
    Eigen::VectorXd fvR;

    //std::cout<<"node="<<node.getName()<<" fv:"<<fv<<"\n";

    if (node.isLeaf()) {
        fv = Eigen::VectorXd::Zero(dim_extended_alphabet);

        //fv[0]=1.0;
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
            lk += node.get_iota() * node.get_beta() * (fv.dot(pi));
        }

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_InsertionHistories()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";*/

        node.set_MSA_fv(fv);

    } else {

        fvL = compute_lk_recursive(node[0], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);
        fvR = compute_lk_recursive(node[1], lk, pi, is_DNA_AA_Codon, dim_extended_alphabet);

        fv = (node.get_left_child()->get_Pr() * fvL).cwiseProduct(node.get_right_child()->get_Pr() * fvR);

        if (node.get_InsertionHistories()) {
            lk += node.get_iota() * node.get_beta() * (fv.dot(pi));
        }

/*        std::cout<<"NODE:"<<node.getName()<<"\n";
        std::cout<<"setA="<<node.get_InsertionHistories()<<"\n";
        std::cout<<"iota="<<node.get_iota()<<"\n";
        std::cout<<"beta="<<node.get_beta()<<"\n";
        std::cout<<"fv:\n";
        std::cout<<fv<<"\n";
        std::cout<<"lk:\n";
        std::cout<<lk<<"\n\n";
*/

        node.set_MSA_fv(fv);


    }
    return fv;
}

//=======================================================================================================
//DP-PIP
double compute_col_lk(PhyTree &tree,
                      Eigen::VectorXd &pi,
                      int is_DNA_AA_Codon,
                      int alphabet_size) {

    double lk;

    compute_lk_recursive(tree, lk, pi, is_DNA_AA_Codon, alphabet_size);

    return log(lk);
}


//===================================================================================================================
double compute_nu(double tau, double lambda, double mu) {

    if (fabs(mu) < 1e-8) {
        perror("ERROR in compute_nu: mu too small");
    }

    return lambda * (tau + 1 / mu);
}

//===================================================================================================================
double phi(int m, double nu, double p0) {
    double p;
    double log_factorial_m;

    log_factorial_m = 0;
    for (int i = 1; i <= m; i++) {
        log_factorial_m += log(i);
    }

    p = -log_factorial_m + m * log(nu) + (nu * (p0 - 1));

    return p;
}