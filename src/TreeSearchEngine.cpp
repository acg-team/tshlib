/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti
 *******************************************************************************
 *
 * This file is part of tshexe
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
 * @file TreeSearchEngine.cpp
 * @author Lorenzo Gatti
 * @date 11 10 2017
 * @version 1.0
 * @maintainer Lorenzo Gatti
 * @email lg@lorenzogatti.me
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

#include <string>
#include "TreeSearchEngine.hpp"


// constructor of TreeSearchEngine,
TreeSearchEngine::TreeSearchEngine(std::string name) {
    this->test_field = name;
}

//copy constructor for making a new copy of a TreeSearchEngine
TreeSearchEngine::TreeSearchEngine(const TreeSearchEngine &copy_from) {

}

//copy assignment for assigning a value from one TreeSearchEngine to another
TreeSearchEngine &TreeSearchEngine::operator=(const TreeSearchEngine &copy_from) {
}

// destructor, just an example
TreeSearchEngine::~TreeSearchEngine() {
    //delete[] this->test_field ;
}

void TreeSearchEngine::setString(std::string input = "") {
    this->test_field = input;
}

int TreeSearchEngine::getLength() {
    return (int) this->test_field.length();
}