/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti
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
 * License along with tshlib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Alphabet.hpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 19 12 2017
 * @version 3.0.1
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

#ifndef TSHLIB_ALPHABET_HPP
#define TSHLIB_ALPHABET_HPP

#include <string>
#include <cassert>
#include <iostream>

namespace tshlib {
#define GAP_CHAR '-'

//=======================================================
//DP-PIP
#define MATCH  '1'
#define GAPX   '2'
#define GAPY   '3'
#define UNKNOW '4'
// for gamma distribution
#define MODELPAREPS   1.E-4
//=======================================================


    class Alphabet {

    protected:
        char data;

        Alphabet() { data = -1; };

    public:
        enum {
            DIM = 0
        };

        static const Alphabet GAP;
        static const Alphabet X;
        static const Alphabet stripStart;
        static const Alphabet stripEnd;

        //=======================================================
        //DP-PIP
        static const Alphabet match;
        static const Alphabet gapX;
        static const Alphabet gapY;
        static const Alphabet unknow;
        //=======================================================


        int value() const;

        std::string asString();

        char asChar() const;

        bool isGap() const;

        bool isValid() const;

        bool isUnknown() const;

        bool operator==(const Alphabet &other) const;

        bool operator!=(const Alphabet &other) const;

        bool operator<(const Alphabet &other) const;
    };

    class AA : Alphabet {
    public:
        explicit AA(char c);

        explicit AA(int i);

        AA() {};

        enum {
            DIM = 20
        };

        static const AA GAP;
        static const AA X;
        static const AA stripStart;
        static const AA stripEnd;


        //=======================================================
        //DP-PIP
        static const AA match;
        static const AA gapX;
        static const AA gapY;
        static const AA unknow;
        //=======================================================


        int value() const;

        std::string asString() const;

        char asChar() const;

        bool isGap() const { return this->data == this->GAP.data; }

        bool isValid() const { return this->value() >= 0 && this->value() < this->DIM; }

        bool isUnknown() const { return this->data == this->X.data; }

        bool operator==(const AA &other) const { return this->data == other.data; }

        bool operator!=(const AA &other) const { return this->data != other.data; }

        bool operator<(const AA &other) const { return this->data < other.data; }
    };

    class Codon : Alphabet {
    public:
        explicit Codon(char c1, char c2, char c3);

        explicit Codon(int i);

        Codon() {};

        enum {
            DIM = 61
        };

        static const Codon GAP;
        static const Codon X;
        static const Codon stripStart;
        static const Codon stripEnd;


        //=======================================================
        //DP-PIP
        static const Codon match;
        static const Codon gapX;
        static const Codon gapY;
        static const Codon unknow;
        //=======================================================


        int value() const;

        std::string asString() const;

        char asChar() const;

        operator AA() const { return AA(this->asChar()); }

        bool isGap() const { return this->data == this->GAP.data; }

        bool isValid() const { return this->value() >= 0 && this->value() < this->DIM; }

        bool isUnknown() const { return this->data == this->X.data; }

        bool operator==(const Codon &other) const { return this->data == other.data; }

        bool operator!=(const Codon &other) const { return this->data != other.data; }

        bool operator<(const Codon &other) const { return this->data < other.data; }
    };

    class DNA : Alphabet {
    public:
        explicit DNA(char c);

        explicit DNA(int i);

        DNA() {};

        enum {
            DIM = 4
        };

        static const DNA GAP;
        static const DNA X;
        static const DNA stripStart;
        static const DNA stripEnd;


        //=======================================================
        //DP-PIP
        static const DNA match;
        static const DNA gapX;
        static const DNA gapY;
        static const DNA unknow;
        //=======================================================


        int value() const;

        std::string asString() const;

        char asChar() const;

        bool isGap() const { return this->data == this->GAP.data; }

        bool isValid() const { return this->value() >= 0 && this->value() < this->DIM; }

        bool isUnknown() const { return this->data == this->X.data; }

        bool operator==(const DNA &other) const { return this->data == other.data; }

        bool operator!=(const DNA &other) const { return this->data != other.data; }

        bool operator<(const DNA &other) const { return this->data < other.data; }
    };

}
#endif // TSHLIB_ALPHABET_HPP
