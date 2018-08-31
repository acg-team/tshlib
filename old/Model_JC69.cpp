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
 * @file Model_JC69.cpp
 * @author Lorenzo Gatti
 * @date 18 12 2017
 * @version 3.0.1
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
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
#include "Model_JC69.hpp"

namespace tshlib {
    jc69_Alpha::jc69_Alpha(bool optimisable, double value) : Parameter(optimisable), value(value) {name="jc69_Alpha";}

    JC69::JC69(double aplha_value) {

        // Initialise parameters
        auto alpha = new jc69_Alpha(true, aplha_value);

        Eigen::MatrixXd Qmatrix = Eigen::MatrixXd::Zero(4,4);

        // Fill Q matrix as for JC69
        for(int r = 0; r<Qmatrix.rows(); r++){
            for(int c=0; c<Qmatrix.cols(); c++ ){

                if(r==c){
                    Qmatrix(r,c) =  -3*alpha->value/4.0;
                }else{
                    Qmatrix(r,c) =  alpha->value/4.0;
                }

            }
        }
        auto Qval = new Q(true, Qmatrix);

        parameters[alpha->name] = alpha;
        parameters[Qval->name] = Qval;

        // Save parameters in the model parameter vector
        params.push_back(alpha);
        params.push_back(Qval);
    }

    JC69::~JC69() = default;
}