/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti
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
 * @file Alignment.hpp
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
#ifndef TSHLIB_ALIGNMENT_HPP
#define TSHLIB_ALIGNMENT_HPP

#include <string>
#include <vector>
#include "Alphabet.h"

class Sequence {
private:
public:
    std::string seq_name;
    std::string seq_data;
    std::vector<unsigned int> seq_weight;
    bool seq_compressed;

    Sequence(std::string label, std::string data);

    Sequence(std::string label, std::string data, std::vector<unsigned int> weight);

    ~Sequence();

protected:

};


class Alignment {

private:

public:

    int align_alphabetsize;
    std::vector<Sequence *> align_dataset;

    Alignment();

    Alignment(int dim);

    void addSequence(std::string label, std::string data);

    void addSequence(std::string label, std::string data, std::vector<unsigned int> weight);

    std::string extractAlignmentColumn(int index);

    ~Alignment();

protected:


};



std::string create_col_MSA(std::vector<std::pair<std::string, std::string>> &MSA, int index);

#endif //TSHLIB_ALIGNMENT_HPP
