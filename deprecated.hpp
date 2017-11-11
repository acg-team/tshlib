/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
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
 * @file deprecated.h
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 11 11 2017
 * @version 1.0
 * @maintainer Lorenzo Gatti
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
#ifndef TSHEXE_DEPRECATED_H
#define TSHEXE_DEPRECATED_H

struct utree_move_info {
    node *node1;
    node *node2;
    int ID;
    double lk;
};

struct move_info {
    PhyTree *node1;
    PhyTree *node2;
    int ID;
    double lk;
};


void nodes_within_radius(PhyTree *start_node, PhyTree *node, int radius, std::vector <move_info> &list_nodes);

void nodes_within_radius_up(PhyTree *start_node, PhyTree *node, int radius, int direction, std::vector <move_info> &list_nodes);

void get_list_nodes_within_radius(PhyTree *node, int radius, std::vector <move_info> &list_nodes_left, std::vector <move_info> &list_nodes_right,
                                  std::vector <move_info> &list_nodes_up);

std::vector<PhyTree *> fill_with_nodes(PhyTree *n);

std::vector<PhyTree *> get_unique(std::vector < PhyTree * > &list_nodes_n1, std::vector < PhyTree * > &list_nodes_n2);

std::vector<PhyTree *> get_path_from_nodes(PhyTree *n1, PhyTree *n2);

std::string create_col_MSA(std::vector <std::pair<std::string, std::string>> &MSA, int index);


#endif //TSHEXE_DEPRECATED_H
