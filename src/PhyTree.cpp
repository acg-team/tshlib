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
 * @file PhyTree.cpp
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

#include <assert.h>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include "PhyTree.hpp"


//=======================================================================================================
static void maxDistPairR(const PhyTree *root, std::vector<double> &distances, std::vector<const PhyTree *> &leaves,
                         const PhyTree *&max, double &max_dist) {
    if (!root->isLeaf()) {
        if (root->n_children() != 2) {
//			error("multifurcations not supported");
        }
        std::vector<double> distances2;
        std::vector<const PhyTree *> leaves2;
        maxDistPairR(&(*root)[0], distances, leaves, max, max_dist);
        maxDistPairR(&(*root)[1], distances2, leaves2, max, max_dist);

        assert(distances.size() == leaves.size());
        assert(distances2.size() == leaves2.size());

        for (unsigned int i = 0; i < distances.size(); ++i) {
            distances[i] += (*root)[0].getBranchLength();
        }
        for (unsigned int j = 0; j < distances2.size(); ++j) {
            distances2[j] += (*root)[1].getBranchLength();
        }

        for (unsigned int i = 0; i < distances.size(); ++i) {
            for (unsigned int j = 0; j < distances2.size(); ++j) {
                if (distances[i] + distances2[j] > max_dist) {
                    max_dist = distances[i] + distances2[j];
                    if (distances[i] > distances2[j]) max = leaves[i];
                    else max = leaves2[j];
                }
            }
        }

        distances.insert(distances.end(), distances2.begin(), distances2.end());
        leaves.insert(leaves.end(), leaves2.begin(), leaves2.end());
    } else {
        distances.push_back(0);
        leaves.push_back(root);
    }
}

//=======================================================================================================
static double maxDistPair(const PhyTree *root, const PhyTree *&max) {
    std::vector<double> distances;
    std::vector<const PhyTree *> leaves;
    double max_dist = -INFINITY;

    maxDistPairR(root, distances, leaves, max, max_dist);

    return max_dist;
}

//=======================================================================================================
PhyTree *midpointRoot(PhyTree *root) {
    // find the two leaves with max distance (or at least one of them plus the distance to the other)
    const PhyTree *max = NULL;
    double dist = maxDistPair(root, max);

    assert(max != NULL);
    PhyTree *current = (PhyTree *) max;

    // find the edge where to place the root
    dist /= 2;
    while (current != root && dist - current->getBranchLength() > 0) {
        dist -= current->getBranchLength();
        PhyTree *parent = current->getParent();
        current = parent;
    }
    if (current == root) {
#ifdef DEBUG
        std::cerr << "tree already midpoint rooted (error " << dist << ")" << std::endl;
#endif
        return root;
    }

    PhyTree *new_root = new PhyTree("new_root");
    double current_dist = current->getBranchLength() - dist;
    double current_support = current->getBranchSupport();
    PhyTree *parent = current->getParent();

    current->pluck();
    new_root->addChild(current, dist, current_support);
    current = new_root;

    while (parent != root) {
        double new_dist = parent->getBranchLength();
        double new_support = parent->getBranchSupport();
        PhyTree *new_parent = parent->getParent();

        parent->pluck();
        current->addChild(parent, current_dist, current_support);

        current = parent;
        parent = new_parent;
        current_dist = new_dist;
        current_support = new_support;
    }

    assert(root->n_children() == 1);

    current_dist += (*root)[0].getBranchLength();
    current_support = std::max(current_support, (*root)[0].getBranchSupport());
    PhyTree *other = root->pluckChild(0);
    current->addChild(other, current_dist, current_support);

    delete root;

    return new_root;
}

//=======================================================================================================
static std::string list_to_name(const std::vector<std::string> &leaves) {
    std::vector<std::string> sorted = leaves;
    std::sort(sorted.begin(), sorted.end());

    std::stringstream ss;
    bool first = true;
    ss << "(";
    for (std::vector<std::string>::const_iterator it = sorted.begin(); it != sorted.end(); ++it) {
        const std::string &name = *it;
        if (name[0] != '(') {
            if (!first) {
                ss << ",";
            } else {
                first = false;
            }
            ss << name;
        }
    }
    ss << ")";
    return ss.str();
}

//=======================================================================================================
std::vector<std::string> get_tree_order_ancestral(const PhyTree *tree) {
    std::vector<std::string> order;

    if (tree->isLeaf()) {
        order.push_back(tree->getName());
    } else {
        for (std::vector<PhyTree *>::const_iterator i = tree->firstChild(); i < tree->lastChild(); ++i) {
            std::vector<std::string> subtree_order;
            subtree_order = get_tree_order_ancestral(*i);

            unsigned int pos = order.size();
            order.insert(order.begin() + pos, subtree_order.begin(), subtree_order.end());

            if (i != tree->firstChild()) {
                std::string anc_name = list_to_name(order);
                order.insert(order.begin() + pos, anc_name);
            }
        }
    }

    return order;
}

//=======================================================================================================
static void get_tree_order_rec(const PhyTree *tree, std::vector<std::string> &order) {
    if (tree->isLeaf()) {
        order.push_back(tree->getName());
    } else {
        for (std::vector<PhyTree *>::const_iterator i = tree->firstChild(); i < tree->lastChild(); ++i) {
            get_tree_order_rec(*i, order);
        }
    }
}

//=======================================================================================================
std::vector<std::string> get_tree_order(const PhyTree *tree) {
//	if(cmdlineopts.ancestral_flag) {
    return get_tree_order_ancestral(tree);
//	}

    std::vector<std::string> order;
    get_tree_order_rec(tree, order);
    return order;
}

//=======================================================================================================
void update_fv_values(std::vector<PhyTree *> &p, int alphabet_size) {
    PhyTree *n;

    for (unsigned int k = 0; k < p.size(); k++) {
        n = p.at(k);
        n->update_MSA_fv(alphabet_size);
    }

}

//=======================================================================================================
bool set_ancestral_flag(PhyTree *node, std::string &MSA_col, int &idx) {
    bool b0, b1, b2;

    if (node->getParent() == NULL) {
        node->set_setA(true);
        return true;
    } else if (node->isLeaf()) {
        b0 = (MSA_col[idx] == '-');
        idx++;
        node->set_setA(b0);
        std::cout << "I am leaf :" << node->getName() << " MSA-col[" << idx << "]=" << MSA_col[idx] << " is gap? ="
                  << b0 << " setA=" << node->get_setA() << "\n";
        return b0;
    } else {
        b1 = set_ancestral_flag(node->get_left_child(), MSA_col, idx);
        b2 = set_ancestral_flag(node->get_right_child(), MSA_col, idx);
        b0 = (b1 && b2);
        node->set_setA(b0);
        std::cout << "I am node :" << node->getName() << " child all gap?=" << b0 << " setA=" << node->get_setA()
                  << "\n";
        return b0;
    }

}

//=======================================================================================================
void set_ancestral_flag(PhyTree *node, std::string &MSA_col) {
    int idx;

    idx = 0;
    set_ancestral_flag(node, MSA_col, idx);

}

//=======================================================================================================
void set_leaf_state(PhyTree *node, std::string &MSA_col, int &idx) {

    if (node->isLeaf()) {
        node->set_leaf_character(MSA_col[idx]);
        idx++;
    } else {
        set_leaf_state(node->get_left_child(), MSA_col, idx);
        set_leaf_state(node->get_right_child(), MSA_col, idx);
    }

}

//=======================================================================================================
void set_leaf_state(PhyTree *node, std::string &MSA_col) {
    int idx;

    idx = 0;
    set_leaf_state(node, MSA_col, idx);

}
//=======================================================================================================





