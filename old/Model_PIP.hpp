/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti
 *******************************************************************************
 *
 * This file is part of Tree Search Heuristic Library (TshLib)
 *
 * Tree Search Heuristic Library (TshLib) is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tree Search Heuristic Library (TshLib) is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Tree Search Heuristic Library (TshLib). If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Model_PIP.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 12 2017
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
#ifndef TSHEXE_MODEL_PIP_HPP
#define TSHEXE_MODEL_PIP_HPP

#include <glog/logging.h>
#include "ModelFactory.hpp"

namespace tshlib {
    class pip_Nu : public Parameter {
    public:
        double value;

        pip_Nu(bool optimisable, double value);

        ~pip_Nu() override = default;

    };

    class pip_Tau : public Parameter {
    public:
        double value;

        pip_Tau(bool optimisable, double value);

        ~pip_Tau() override = default;
    };

    class pip_Phi : public Parameter {
    public:
        double value;

        pip_Phi(bool optimisable, double value);

        ~pip_Phi() override = default;
    };

    class pip_Lambda : public Parameter {
    public:
        double value;

        pip_Lambda(bool optimisable, double value);

        ~pip_Lambda() override = default;

        void Optimise() override {};

        void UpdateList() override {};
    };

    class pip_Mu : public Parameter {

    public:
        double value;

        pip_Mu(bool optimisable, double value);

        ~pip_Mu() override = default;

        void Optimise() override {};

        void UpdateList() override {};
    };

    class PIP : public ModelFactory {
    public:
        // TODO: Check for better ways to extend a derived class using template
        PIP(ModelFactory *model, double lambda_value, double mu_value) {

            params = model->params;

            // Extend model parameters
            parameters = model->parameters;

            // Initialise model parameters
            auto lambda = new pip_Lambda(true, lambda_value);
            auto mu = new pip_Mu(true, mu_value);

            Q *qmatrix = dynamic_cast<Q *>(parameters["Q"]);

            long rows = qmatrix->value.rows();
            long cols = qmatrix->value.cols();

            (qmatrix->value).conservativeResize(rows + 1, cols + 1);

            rows = qmatrix->value.rows();
            cols = qmatrix->value.cols();

            // Fill Q matrix as for JC69
            for (int r = 0; r < rows - 1; r++) {
                for (int c = 0; c < cols - 1; c++) {

                    if (r == c) {
                        qmatrix->value(r, c) = qmatrix->value(r, c) - mu->value;
                    } else {
                        qmatrix->value(r, c) = qmatrix->value(r, c);
                    }

                }
                qmatrix->value(r, cols - 1) = mu->value;
            }

            // Store new parameters
            parameters[lambda->name] = lambda;
            parameters[mu->name] = mu;

            params.push_back(lambda);
            params.push_back(mu);


        };

        void setLambda();

        void setMu();

        void computeNu();

        void computeTau();

        void computePhi();

    };


}
#endif //TSHEXE_MODEL_PIP_HPP
