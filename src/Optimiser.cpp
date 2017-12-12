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
 * @file Optimiser.cpp
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
#include "Optimiser.hpp"



template<class modelType>
Optimiser<modelType>::Optimiser(Utree *tree, Alignment *alignment, modelType *model, Optimiser::LimitType optlimit, double vallimit) {

    this->alignment = alignment;
    this->tree = tree;
    this->model = model;
    this->currScore = 0;
    this->initScore = 0;
    this->limit_type = optlimit;
    this->limit_value = vallimit;
    this->nodeList.clear();

}

template<class modelType>
Optimiser<modelType>::Optimiser(Utree *tree, Alignment *alignment, modelType *model, double initialScore) {

    this->alignment = alignment;
    this->tree = tree;
    this->model = model;
    this->currScore = 0;
    this->initScore = initialScore;
    this->limit_type = Optimiser::score;
    this->limit_value = 0;
    this->nodeList.clear();


}

template<class modelType>
Utree *Optimiser<modelType>::getTree() const {
    return tree;
}

template<class modelType>
Alignment *Optimiser<modelType>::getAlignment() const {
    return alignment;
}

template<class modelType>
void Optimiser<modelType>::setAlignment(Alignment *alignment) {
    Optimiser::alignment = alignment;
}

template<class modelType>
void Optimiser<modelType>::computeSteps() {

/*    switch(this->model->parameters){

        case modelType::model

            break;


    }*/

}

template<class modelType>
void Optimiser<modelType>::executeOptimisation() {

}

template<class modelType>
Optimiser<modelType>::Optimiser() = default;

template<class modelType>
Optimiser<modelType>::~Optimiser() = default;



