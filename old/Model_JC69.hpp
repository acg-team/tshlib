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
 * @file Model_JC69.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 12 2017
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
#ifndef TSHEXE_MODEL_JC69_HPP
#define TSHEXE_MODEL_JC69_HPP
#include "ModelFactory.hpp"

namespace tshlib {
    class jc69_Alpha : public Parameter {
    public:
        double value;

        jc69_Alpha(bool optimisable, double value);

        ~jc69_Alpha() override = default;

        void Optimise() override {};

        void UpdateList() override {};
    };

    class JC69 : public ModelFactory {
    public:
        explicit JC69(double alpha_value);

        virtual ~JC69();

    };
}

#endif //TSHEXE_MODEL_JC69_HPP
