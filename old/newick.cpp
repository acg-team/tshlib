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
 * @file newick.cpp
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


#include <iostream>
#include <sstream>
#include "newick.hpp"

using namespace std;

namespace newick_parser {

    class TokenBuffer {

    private:
        istream *in;
        string last;
        bool valid;

        static bool issep(char c) {
            return (c == ',' || c == ':' || c == '(' || c == ')' || c == ';');
        }

        string next_token() throw(LexerException) {
            string token = "";
            while (*this->in) {
                char c = this->in->get();
                if (isspace(c)) {
                    continue;
                } else if (issep(c)) {
                    return token + c;
                } else {
                    token += c;
                    while (*this->in && (!isspace(this->in->peek()) && !issep(this->in->peek()))) {
                        char c = this->in->get();
                        token += c;
                    }
                    return token;
                }
            }

            throw LexerException("Unexpected EOF or I/O error");
        }

    public:
        TokenBuffer(istream *in) {
            this->in = in;
            valid = false;
        }

        string peek() throw(LexerException) {
            if (!valid) {
                last = next_token();
            }
            valid = true;
            return last;
        }

        string next() throw(LexerException) {
            if (!valid) {
                last = next_token();
            }
            valid = false;
            return last;
        }
    };

    static double parse_double(const string &str) throw() {
        istringstream ss(str.c_str());
        double out = 0;

        ss >> out;

        return out;
    }

    static PhyTree *parse_tree(TokenBuffer &buffer) throw(ParserException, LexerException) {
        PhyTree *t = new PhyTree();

        string tok = buffer.next();

        if (tok != "(")
            throw ParserException("Unexpected token: '" + tok + "', expected: '('");

        do {
            PhyTree *child = NULL;
            tok = buffer.peek();

            if (tok == "(")
                child = parse_tree(buffer);
            else
                child = new PhyTree(buffer.next());

            tok = buffer.next();
            if (tok != ":") {
                double support = parse_double(tok);
                (void) support; // ignore branch support

                tok = buffer.next();
                if (tok != ":")
                    throw ParserException("Unexpected token: '" + tok + "', expected: ':'");
            }

            tok = buffer.next();
            double dist = parse_double(tok);

            t->addChild(child, dist);

            tok = buffer.peek();
            if (tok == ")") {
                tok = buffer.next();
                break;
            }

            tok = buffer.next();
            if (tok != ",")
                throw ParserException("Unexpected token: '" + tok + "', expected: ','");
        } while (true);

        return t;
    }

    PhyTree *parse_newick(istream *in) throw(ParserException, LexerException) { TokenBuffer buffer(in);
        PhyTree *t = parse_tree(buffer);

        string tok = buffer.next();
        if (tok != ";") {
            if (tok != ":") {
                (void) parse_double(tok);
                tok = buffer.next();
            }
            if (tok != ":")
                throw ParserException("Unexpected token: " + tok);
            tok = buffer.next();
            (void) parse_double(tok);
            tok = buffer.next();
        }

        if (tok != ";")
            throw ParserException("Unexpected token: " + tok);

        return t;
    }

}
