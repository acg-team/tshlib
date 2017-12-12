/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshexe
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
#ifndef TSHEXE_OPTIMISER_HPP
#define TSHEXE_OPTIMISER_HPP


#include "Utree.hpp"
#include "Likelihood.hpp"

template <class modelType>
class Optimiser {
    enum LimitType {score, iterations};

private:

    Optimiser::LimitType limit_type;
    double limit_value;
    std::vector<VirtualNode *> nodeList;
    double initScore;
    double currScore;

public:

    Utree *getTree() const;

    Alignment *getAlignment() const;

    void setAlignment(Alignment *alignment);

    Optimiser();
    Optimiser(Utree *tree, Alignment *alignment, modelType *model, Optimiser::LimitType optlimit, double vallimit);
    Optimiser(Utree *tree, Alignment *alignment, modelType *model, double initialScore);

    virtual void computeSteps();

    virtual void executeOptimisation();

    virtual ~Optimiser();


protected:

    Utree *tree;
    Alignment *alignment;
    modelType *model;
    void *(objFunction);

};


#endif //TSHEXE_OPTIMISER_HPP
