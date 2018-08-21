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
 * License along with tshlib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Model.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 12 12 2017
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
 * @see For more information visit: http://www.lorenzogatti.me
 */
#ifndef TSHEXE_MODEL_HPP
#define TSHEXE_MODEL_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <map>

#include "Utree.hpp"

namespace tshlib {

    class Parameter {
    public:
        std::string name;
        bool optimisable;

        explicit Parameter();

        explicit Parameter(bool optimisable);

        virtual ~Parameter();

        virtual void Optimise();

        virtual void UpdateList();
    };


    class topology : public Parameter {
        //TODO: include all the methods from tree_rearrangment
    public:
        Utree *value;

        topology(bool optimisable, Utree *tree);

        ~topology() override = default;

        void Optimise() override {};

        void UpdateList() override {};

    };

    class msa : public Parameter {

    public:
        Alignment *value;

        msa(bool optimisable, Alignment *alignment);

        ~msa() override = default;

        void Optimise() override {};

        void UpdateList() override {};

    };


    class Vi : public Parameter {

    public:
        Eigen::MatrixXd value;

        Vi(bool optimisable, const Eigen::MatrixXd &Vi);

        ~Vi() override = default;

        void Optimise() override {};

        void UpdateList() override {};

    };

    class pi : public Parameter {
    public:
        Eigen::VectorXd value;

        pi(bool optimisable, const Eigen::VectorXd &pi);

        ~pi() override = default;

        void Optimise() override {};

        void UpdateList() override {};


    };

    class Q : public Parameter {
    public:
        Eigen::MatrixXd value;

        Q(bool optimisable, const Eigen::MatrixXd &Qmatrix);

        ~Q() override = default;

        void Optimise() override {};

        void UpdateList() override {};

    };

    class V : public Parameter {

    public:
        Eigen::MatrixXd value;

        V(bool optimisable, const Eigen::MatrixXd &V);

        ~V() override = default;

        void Optimise() override {};

        void UpdateList() override {};
    };

    class ModelFactory {
    public:
        std::vector<Parameter *> params;
        std::map<std::string, Parameter *> parameters;

    };


}

#endif //TSHEXE_MODEL_HPP
