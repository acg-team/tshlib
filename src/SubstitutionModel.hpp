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

/*
template <class parType>
class Parameter{

public:
    std::string par_name;
    parType par_value;
    bool par_optimisable;
    bool par_nodetype;

    Parameter(const std::string &par_name, parType par_value, bool par_optimisable, bool par_nodetype) : par_name(par_name),
                                                                                                         par_value(par_value),
                                                                                                         par_optimisable(par_optimisable),
                                                                                                         par_nodetype(par_nodetype) {}

    virtual ~Parameter() = default;

};

template <typename parType>
class SubstitutionModel {
protected:
    Utree *tree;
    Alignment *alignment;
    
    Parameter<Eigen::VectorXd> pi;
    Parameter<Eigen::VectorXd> Q;
    Parameter<Eigen::MatrixXd> V;
    Parameter<Eigen::MatrixXd> Vi;

    std::vector<std::unique_ptr<Parameter<parType>>> submod_parameters;
    //std::vector<Parameter *> submod_parameters;

public:

    double lk;

    SubstitutionModel();
    virtual ~SubstitutionModel();




    virtual double computeLikelihood();
    virtual void setTree();
    virtual Utree* getTree();
    virtual void setAlignment();
    virtual Alignment *getAlignment();


};

class Model_JC69 : public SubstitutionModel{

public:

    Model_JC69(double valScalingFactor);

    virtual ~Model_JC69();

};


template <class RevSubModel>
class Model_PIP : public SubstitutionModel, public RevSubModel{

    enum class param{lamda, mu, tau, branch, qmatrix, pi, undef};

public:

    Parameter<double> lambda;
    Parameter<double> mu;
    Parameter<double> nu;
    Parameter<double> tau;
    Parameter<double> phi;

    Model_PIP();
    virtual ~Model_PIP();

};


class NodeParameters_PIP{

public:
    Eigen::VectorXd Fv_empty;
    Eigen::MatrixXd Fv;
    Eigen::MatrixXd Pr;


    NodeParameters_PIP();
    virtual ~NodeParameters_PIP();

};


*/

#endif //TSHEXE_MODEL_HPP
