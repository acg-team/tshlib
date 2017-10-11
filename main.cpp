/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti & Massimo Maiolo
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
 * License along with likpip. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file main.cpp
 *  * @date 11 10 2017
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
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "PhyTree.hpp"
#include "newick.hpp"

//===================================================================================================================
//===================================================================================================================
int main(int argc, char **argv) {

    //std::string tree_file = "/home/max/PIP_C++/NNI_SPR/tree_5_leaves_r_bl.nwk";
    std::string tree_file = argv[1];
    PhyTree *tree = nullptr;

    std::ifstream tree_str(tree_file.c_str());
    tree = newick_parser::parse_newick(&tree_str);

    tree->set_missing_node_name("V");

    std::cout << tree->formatNewick() << "\n\n";

    tree->print();
    std::cout << "\n";

    //----------------------------------------------------------
//	PhyTree *t1;
//
//	t1=tree->get_right_child()->copy();
//	t1->null_parent();
//	double bl=t1->getBranchLength();
//
//	std::cout<<t1->getName()<<"\n";
//	std::cout<<tree->n_children()<<"\n";
//
//	tree->deleteChild(1);
//
//	std::cout<<tree->n_children()<<"\n";
//
//	tree->addChild(t1,bl,0);
//
//	tree->print();
//	std::cout<<"\n";
//
//	std::cout<<tree->formatNewick()<<"\n\n";
    //----------------------------------------------------------
    PhyTree *t1;
    PhyTree *t2;

    t1 = tree->get_right_child();
    t2 = tree->get_left_child()->get_right_child();

//	tree->swap(tree,1,tree->get_left_child(),1);
    tree->swap2(t1, t2);

    std::cout << tree->formatNewick() << "\n\n";
    //----------------------------------------------------------

    return 0;
}





