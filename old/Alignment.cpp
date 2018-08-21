/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
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
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Alignment.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 10 2017
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
 * @see For more information visit:
 */

#include <string>
#include <vector>
#include "Alignment.hpp"
namespace tshlib {

    Alignment::Alignment(bool compressed) {
        this->align_alphabetsize = 0;
        this->align_compressed = compressed;
        this->align_alphabetsize = 0;
    }

    Alignment::Alignment() {
        this->align_compressed = false;
        this->align_alphabetsize = 0;

    }

    Alignment::~Alignment() = default;


    void Alignment::addWeight(std::vector<int> column_weight) {
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

        this->align_length = length;

        return length;
    }

    void Alignment::countNumberCharactersinColumn() {

        int numCharacters;

        for (int c = 0; c < this->getAlignmentSize(); c++) {
            numCharacters = 0;

            for (int i = 0; i < this->align_dataset.size(); i++) {

                char ch = this->align_dataset.at(i)->seq_data.at(c);
                if (ch != '-') {
                    numCharacters++;
                }

                //num_gaps += (int)(this->align_dataset.at(c)->seq_data.at(i) != '-');
            }

            this->align_num_characters.at(c) = numCharacters;

        }


    }


    Alignment_AA::Alignment_AA() {

        this->align_alphabetsize = DIM;
        this->alphabet = AlignmentAlphabet::aa;

    }

    Alignment_DNA::Alignment_DNA() {

        this->align_alphabetsize = DIM;
        this->alphabet = AlignmentAlphabet::dna;

    }

    Alignment_Codon::Alignment_Codon() {

        this->align_alphabetsize = DIM;
        this->alphabet = AlignmentAlphabet::codon;

    }

}