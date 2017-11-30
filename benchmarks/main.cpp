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
 * License along with tshexe. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file main.cpp
 * @author Lorenzo Gatti
 * @date 28 11 2017
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
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <glog/logging.h>
#include <gflags/gflags.h>

#include <Utree.hpp>
#include <TreeRearrangment.hpp>
#include <Alignment.hpp>
#include <Likelihood.hpp>
#include <newick.hpp>

#include "main.hpp"

int main(int argc, char **argv) {


    // Initialize Google's logging library.
    FLAGS_alsologtostderr = true;
    google::InitGoogleLogging("THSLIB Main File");

    //------------------------------------------------------------------------------------------------------------------
    std::string tree_file = argv[1];
    PhyTree *tree = nullptr;
    // INIT TREE

    // tree filename
    std::ifstream tree_str(tree_file.c_str());
    LOG(INFO) << "TSHLIB: " << tree_file << std::endl;

    // read newick file
    tree = newick_parser::parse_newick(&tree_str);

    // set name of internal nodes
    tree->set_missing_node_name("V");
    VLOG(1) << "[Initial Tree Topology] " << tree->formatNewick() << std::endl;

    //------------------------------------------------------------------------------------------------------------------
    // BUILD UNROOTED TREE

    auto utree = new Utree;
    UtreeUtils::convertUtree(tree, utree);
    delete tree;

    VLOG(1) << "[Initial utree] " << utree->printTreeNewick(true) << std::endl;

    //------------------------------------------------------------------------------------------------------------------
    // Execute tree rearrangements

    treesearchheuristics::testTSH(utree, TreeSearchHeuristics::classic_Mixed);

    //------------------------------------------------------------------------------------------------------------------
    // Clear up


    delete utree;

    exit(0);
}