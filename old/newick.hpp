/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
 *******************************************************************************
 *
 * This file is part of tshlib
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
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file newick.hpp
 * @author Lorenzo Gatti & Massimo Maiolo
 * @date 11 10 2017
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

#ifndef TSHLIB_NEWICK_HPP
#define TSHLIB_NEWICK_HPP

#include "PhyTree.hpp"
#include <iostream>
#include <string>

namespace newick_parser {

    class LexerException : public std::exception {
    private:
        std::string error;
    public:
        LexerException(std::string error) throw() {
            this->error = error;
        }

        ~LexerException() throw() {
        }

        virtual const char *what() const throw() {
            return error.c_str();
        }
    };

    class ParserException : public std::exception {
    private:
        std::string error;
    public:
        ParserException(std::string error) throw() {
            this->error = error;
        }

        ~ParserException() throw() {
        }

        virtual const char *what() const throw() {
            return error.c_str();
        }
    };

    PhyTree *parse_newick(std::istream *in) throw(ParserException, LexerException);

}
#endif /* TSHLIB_NEWICK_HPP */
