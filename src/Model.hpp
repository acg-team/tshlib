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
 * License along with tshlib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Model.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 12 12 2017
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
 * @see For more information visit: http://www.lorenzogatti.me
 */
#ifndef TSHEXE_MODEL_HPP
#define TSHEXE_MODEL_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "Utree.hpp"


class Model {
protected:
    Utree *tree;
    Alignment *alignment;

public:
    double lk;

    Model();
    virtual ~Model();

/*
    virtual double computeLikelihood();
    virtual void setTree();
    virtual Utree* getTree();
    virtual void setAlignment();
    virtual Alignment *getAlignment();
*/

};

class modelPIP : public Model{

public:

    double lambda;
    double mu;
    Eigen::VectorXd pi;
    Eigen::MatrixXd Q;
    double nu;
    double tau;
    double phi;
    Eigen::MatrixXd V;
    Eigen::MatrixXd Vi;

    modelPIP();
    virtual ~modelPIP();

    /*virtual double computeLikelihood();*/

};

class Parameters{

public:
    Parameters();
    virtual ~Parameters();

};


class parametersPIP: public Parameters{

    enum parameters{lamda, mu, tau, branch, qmatrix, pi, undef};

public:
    Eigen::VectorXd Fv_empty;
    Eigen::MatrixXd Fv;
    Eigen::MatrixXd Pr;

    parametersPIP();
    virtual ~parametersPIP();

};




#endif //TSHEXE_MODEL_HPP
