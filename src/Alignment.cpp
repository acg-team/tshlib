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


Alignment::Alignment(bool compressed) {

    this->align_compressed = compressed;
    this->align_alphabetsize = 0;
}


Alignment::Alignment() {
    this->align_compressed = false;
    this->align_alphabetsize = 0;
}


Alignment::~Alignment() = default;


void Alignment::addWeight(std::vector<unsigned int> column_weight) {
    if (this->align_compressed) {
        this->align_weight = column_weight;
        this->align_compressed = true;
    }
}


std::string Alignment::extractColumn(int index) {

    std::string column;

    for (std::vector<int>::size_type i = 0; i != this->align_dataset.size(); i++) {
        column.push_back(this->align_dataset.at(i)->seq_data[index]);

    }

    return column;
}


void Alignment::addSequence(std::string label, std::string data) {

    this->align_dataset.emplace_back(new Sequence(label, data));

    // Add automatic weight if the alignment is not compressed
    if (!this->align_compressed && this->align_weight.size() != data.size()) {
        this->align_weight.clear();
        for (std::vector<int>::size_type i = 0; i != data.size(); i++) {
            this->align_weight.emplace_back(1);
        }
    }


}


long int Alignment::getAlignmentSize() {

    long length = 0;

    for (std::vector<int>::size_type i = 0; i != this->align_dataset.size(); i++) {

        if (length < this->align_dataset[i]->seq_data.size()) {

            length = this->align_dataset[i]->seq_data.size();
        }


    }
    return length;
}


std::string create_col_MSA(std::vector<std::pair<std::string, std::string>> &MSA, int index) {
    std::string colMSA;

    for (unsigned int i = 0; i < MSA.size(); i++) {
        colMSA.append(MSA.at(i).second, index, 1);
    }

    return colMSA;
}
