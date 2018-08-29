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
 * @file Alphabet.cpp
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

#include "Alphabet.hpp"

namespace tshlib {
    const char aa_translation_table[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, 0, 20, 1, 2, 3, 4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14,
                                            15, 16, 20, 17, 18, 20, 19, 20, -1, -1, -1, -1, -1, -1, 0, 20, 1, 2, 3,
                                            4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14, 15, 16, 20, 17, 18, 20,
                                            19, 20, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1};

    const char dna_translation_table[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, 2, -1, 1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, 0, 0, -1, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, 1,
                                             -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0,
                                             -1, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    const char dna_inv_translation_table[5] = {'T', 'C', 'A', 'G', 'X'};

    const char aa_inv_translation_table[21] = {'A', 'C', 'D', 'E', 'F', 'G', 'H',
                                               'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};

    const char *codon_inv_translation_table[62] = {"TTT", "TTC", "TTA",
                                                   "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TGT", "TGC", "TGG",
                                                   "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC",
                                                   "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG",
                                                   "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC",
                                                   "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
                                                   "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG", "XXX"};

    const char codon_from_product_table[64] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1,
                                               -1, 10, 11, -1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                                               26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
                                               44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60};

    const char codon_inv2_translation_table[63] = {
            "FfUuSs5$YyCcWLl1|PpOoHhQqRrBbIi!MTt7+NnKkZz:;Vvw^Aa4@DdEeGg6&X"};

    const char codon_inv3_translation_table[63] = {
            "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGGX"};

    const char codon_translation_table[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, 31, -1, -1, 7, -1, 60, -1, -1, -1, -1, 36, -1,
                                               -1, -1, -1, -1, 15, -1, -1, 51, 6, 59, 35, -1, -1, 43, 44, -1, -1, -1,
                                               -1, 52, 49, 27, 10, 53, 55, 0, 57, 21, 29, -1, 39, 13, 32, 37, 19, 17,
                                               23, 25, 4, 33, 2, 45, 12, 61, 8, 41, -1, -1, -1, 48, -1, -1, 50, 28,
                                               11, 54, 56, 1, 58, 22, 30, -1, 40, 14, -1, 38, 20, 18, 24, 26, 5, 34,
                                               3, 46, 47, -1, 9, 42, -1, 16, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    const AA AA::X = AA('X');
    const DNA DNA::X = DNA('X');
    const Codon Codon::X = Codon('X', 'X', 'X');

    const AA AA::GAP = AA(GAP_CHAR);
    const DNA DNA::GAP = DNA(GAP_CHAR);
    const Codon Codon::GAP = Codon(GAP_CHAR, GAP_CHAR, GAP_CHAR);

    const AA AA::stripStart = AA('M');
    const DNA DNA::stripStart = DNA(GAP_CHAR);
    const Codon Codon::stripStart = Codon('A', 'T', 'G');

    const AA AA::stripEnd = AA(GAP_CHAR);
    const DNA DNA::stripEnd = DNA(GAP_CHAR);
    const Codon Codon::stripEnd = Codon('X', 'X', 'X');


//=======================================================
//DP-PIP

    const AA AA::match = AA(MATCH);
    const DNA DNA::match = DNA(MATCH);
    const Codon Codon::match = Codon(MATCH, MATCH, MATCH);

    const AA AA::gapX = AA(GAPX);
    const DNA DNA::gapX = DNA(GAPX);
    const Codon Codon::gapX = Codon(GAPX, GAPX, GAPX);

    const AA AA::gapY = AA(GAPY);
    const DNA DNA::gapY = DNA(GAPY);
    const Codon Codon::gapY = Codon(GAPY, GAPY, GAPY);

    const AA AA::unknow = AA(UNKNOW);
    const DNA DNA::unknow = DNA(UNKNOW);
    const Codon Codon::unknow = Codon(UNKNOW, UNKNOW, UNKNOW);
//=======================================================


    AA::AA(char c) {
        if (c == '_' || c == '-' || c == '.' || c == ' ')
            c = GAP_CHAR;
        this->data = c;
    }

    AA::AA(int i) {
        if (i < 0 || i > AA::DIM) {
            this->data = '?';
        } else {
            this->data = aa_inv_translation_table[i];
        }
    }

    int AA::value() const {
        return aa_translation_table[(int) this->data];
    }

    std::string AA::asString() const {
        return std::string() + this->data;
    }

    char AA::asChar() const {
        return this->data;
    }

    Codon::Codon(char c1, char c2, char c3) {
        int c;

        if (c1 == '_' || c2 == '_' || c3 == '_' ||
            c1 == '-' || c2 == '-' || c3 == '-' ||
            c1 == ' ' || c2 == ' ' || c3 == ' ' ||
            c1 == '.' || c2 == '.' || c3 == '.')
            goto gap_codon;

        c = dna_translation_table[(int) (unsigned char) c3];
        if (c < 0) goto invalid_codon;
        if (c >= 4) goto unknown_codon;

        c += 4 * dna_translation_table[(int) (unsigned char) c2];
        if (c < 0) goto invalid_codon;
        if (c >= 16) goto unknown_codon;

        c += 16 * dna_translation_table[(int) (unsigned char) c1];
        if (c < 0) goto invalid_codon;
        if (c >= 64) goto unknown_codon;

        this->data = codon_from_product_table[c];
        return;

        invalid_codon:
        this->data = -1;
        return;

        unknown_codon:
        this->data = this->DIM;
        return;

        gap_codon:
        this->data = this->DIM + 1;
        return;
    }

    Codon::Codon(int i) {
        if (i < 0 || i > Codon::DIM) {
            this->data = -1;
        } else {
            this->data = i;
        }
    }

    int Codon::value() const {
        if (this->isGap()) return -1;
        return this->data;
    }

    std::string Codon::asString() const {
        if (this->isGap()) {
            char gap[4] = {GAP_CHAR, GAP_CHAR, GAP_CHAR, '\0'};
            return std::string((char *) gap);
        } else if (!this->isValid()) {
            return std::string("XXX");
        }
        return std::string(codon_inv_translation_table[(int) this->data]);
    }

    char Codon::asChar() const {
        if (this->isGap()) {
            return GAP_CHAR;
        } else if (!this->isValid()) {
            return 'X';
        }
        return codon_inv3_translation_table[(int) this->data];
    }

    DNA::DNA(char c) {
        if (c == '_' || c == '-' || c == '.' || c == ' ')
            c = GAP_CHAR;
        this->data = c;
    }

    DNA::DNA(int i) {
        if (i < 0 || i >= (int) sizeof(dna_inv_translation_table)) {
            this->data = '?';
        } else {
            this->data = dna_inv_translation_table[i];
        }
    }

    int DNA::value() const {
        return dna_translation_table[(int) this->data];
    }

    std::string DNA::asString() const {
        return std::string() + this->data;
    }

    char DNA::asChar() const {
        return this->data;
    }

}