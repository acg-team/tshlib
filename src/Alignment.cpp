/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * tshexe is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tshexe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with likpip. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Alignment.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 10 2017
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
 * @see For more information visit:
 */

#include <string>
#include <vector>
#include "Alignment.hpp"
#include "Alphabet.h"

//=======================================================================================================

//DP-PIP
std::string create_col_MSA(std::vector<std::pair<std::string, std::string>> &MSA, int index) {
    std::string colMSA;

    for (unsigned int i = 0; i < MSA.size(); i++) {
        colMSA.append(MSA.at(i).second, index, 1);
    }

    return colMSA;
}
//=======================================================================================================



//================================

std::string Alignment::extractAlignmentColumn(int index) {
    return std::string();
}


Alignment::Alignment(int dim) {

    this->align_alphabetsize = dim;
}


Alignment::Alignment() {

    this->align_alphabetsize = 0;
}

void Alignment::addSequence(std::string label, std::string data) {

    this->align_dataset.emplace_back(new Sequence(label, data));
}

void Alignment::addSequence(std::string label, std::string data, std::vector<unsigned int> weight) {

    this->align_dataset.emplace_back(new Sequence(label, data, weight));

}

//================================

Alignment::~Alignment() = default;
//================================



Sequence::Sequence(std::string label, std::string data) {

    this->seq_name = label;
    this->seq_data = data;
    this->seq_compressed = false;

    for (std::vector<int>::size_type i = 0; i != this->seq_data.size(); i++) {
        this->seq_weight.push_back(1);
    }

}

//================================
Sequence::Sequence(std::string label, std::string data, std::vector<unsigned int> weight) {

    this->seq_name = label;
    this->seq_data = data;
    this->seq_compressed = true;
    this->seq_weight = weight;

}

//================================
Sequence::~Sequence() = default;

