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
 * @file TreeSearchEngine.hpp
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
#ifndef TSHEXE_TREESEARCHENGINE_HPP
#define TSHEXE_TREESEARCHENGINE_HPP

#include <string>

using namespace std;


class TreeSearchEngine {

private:
    string test_field;

public:

    TreeSearchEngine(string name);

    TreeSearchEngine(const TreeSearchEngine &copy_from);

    TreeSearchEngine &operator=(const TreeSearchEngine &copy_from);

    ~TreeSearchEngine();

    void setString(string input);

    int getLength();

};

#endif //TSHEXE_TREESEARCHENGINE_HPP