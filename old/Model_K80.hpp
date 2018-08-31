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
 * @file Model_K80.hpp
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
#ifndef TSHEXE_MODEL_K80_HPP
#define TSHEXE_MODEL_K80_HPP

#include "ModelFactory.hpp"

namespace tshlib {
    class k80_kappa : public Parameter {

    public:
        double value;

        k80_kappa(bool optimisable, double value);

        ~k80_kappa() override = default;

        void Optimise() override {};

        void UpdateList() override {};

    };


    class K80 : public ModelFactory {

    public:

        K80(double value_kappa);

        virtual ~K80();

    };
}
#endif //TSHEXE_MODEL_K80_HPP
