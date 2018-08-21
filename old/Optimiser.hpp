/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of Tree Search Heuristic Library (TshLib)
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
 * @file Optimiser.hpp
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
#ifndef TSHEXE_OPTIMISER_HPP
#define TSHEXE_OPTIMISER_HPP


#include "Utree.hpp"
#include "Likelihood.hpp"


//class Optimiser{
//public:
//    ModelFactory *model;
//
//    Optimiser(ModelFactory *model);
//
//    virtual ~Optimiser();
//
//
//};


namespace tshlib {
    template<class Model>
    class Optimiser {
        enum class LimitType {
            score, iterations
        };

    private:

        Optimiser::LimitType limit_type;
        double limit_value;
        std::vector<VirtualNode *> nodeList;

    public:

        Model *model;

        Optimiser(Model *model);

        virtual ~Optimiser();

        virtual void computeSteps();

        virtual void executeOptimisation();


    protected:


    };

}
#endif //TSHEXE_OPTIMISER_HPP
