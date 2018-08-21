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
 * @file tshlib_utree_test.cpp
 * @author Lorenzo Gatti
 * @date 27 04 2018
 * @version 2.0.2
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


#include "gtest/gtest.h"
#include <Utree.hpp>
#include "testData.hpp"

using namespace tshlib;




TEST(TreeDepth, Leaves_and_Internal_0) {
    auto t = new testData;
    Utree *tree = t->getTree();
    tree->computeTreeDepth();

    EXPECT_EQ(tree->getTreeDepthAtNode(tree->listVNodes.at(1)), 4);

}

TEST(TreeDepth, Leaves_and_Internal_1) {
    auto t = new testData;
    Utree *tree = t->getTree();
    tree->computeTreeDepth();

    EXPECT_EQ(tree->getTreeDepthAtNode(tree->listVNodes.at(3)), 2);

}

TEST(TreeDepth, Leaves_and_Internal_2) {
    auto t = new testData;
    Utree *tree = t->getTree();
    tree->computeTreeDepth();

    EXPECT_EQ(tree->getTreeDepthAtNode(tree->listVNodes.at(4)), 2);

}


