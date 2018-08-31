/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * Tree Search Heuristic Library (TshLib) is a free software: you can redistribute
 * it and/or modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tree Search Heuristic Library (TshLib) is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file TreeSearchEngine.hpp
 * @author Lorenzo Gatti
 * @date 11 10 2017
 * @version 3.0.1
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
#ifndef TSHLIB_TREESEARCH_HPP
#define TSHLIB_TREESEARCH_HPP

#include <string>

using namespace std;

namespace tshlib {

    class TreeSearch {

    private:
        string test_field;

    public:

        TreeSearch(string name);

        TreeSearch(const TreeSearch &copy_from);

        TreeSearch &operator=(const TreeSearch &copy_from);

        ~TreeSearch();

        void setString(string input);

        int getLength();

    };
}
#endif //TSHLIB_TREESEARCH_HPP
