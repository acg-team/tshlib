//
// Created by Lorenzo Gatti on 26/10/17.
//

#include <string>
#include "Utree.hpp"

std::string utree_formatNewickR(node *n, bool is_root) {

    if (n->next == NULL) {
        return n->data->getName();
    } else {
        std::stringstream newick;
        if (is_root) {
            newick << "(";
            newick << utree_formatNewickR(n->back, false) << ":" << n->back->data->getBranchLength();
            newick << ",";
            newick << utree_formatNewickR(n->next->back, false) << ":" << n->next->back->data->getBranchLength();
            newick << ",";
            newick << utree_formatNewickR(n->next->next->back, false) << ":"
                   << n->next->next->back->data->getBranchLength();
            newick << ")";
        } else {
            newick << "(";
            newick << utree_formatNewickR(n->next->back, false) << ":" << n->next->back->data->getBranchLength();
            newick << ",";
            newick << utree_formatNewickR(n->next->next->back, false) << ":"
                   << n->next->next->back->data->getBranchLength();
            newick << ")";
        }

        return newick.str();
    }

}

std::string utree_formatNewick(node *utree_pseudo_root) {
    std::string s;

    if (utree_pseudo_root->next == NULL) {
        return NULL;
    }

    /*
    s=utree_formatNewickR(utree_pseudo_root->back)+
            utree_formatNewickR(utree_pseudo_root->next->back)+
            utree_formatNewickR(utree_pseudo_root->next->next->back)+ ";";
    */

    s = utree_formatNewickR(utree_pseudo_root, true) + ";";

    return s;
}

void print_node_neighbours(node *n) {

    std::cout << n->data->getName() << " ";

    if (n->next != NULL) {
        std::cout << "(";
        std::cout << "^" << n->back->data->getName() << ";";
        std::cout << "<" << n->next->back->data->getName() << ";";
        std::cout << n->next->next->back->data->getName() << ">";
        std::cout << ")";
    } else {
        std::cout << "(";
        std::cout << "^" << n->back->data->getName() << ";";
        std::cout << "<" << "-" << ";";
        std::cout << "-" << ">";
        std::cout << ")";
    }

    std::cout << "\n";

}

void print_utree_rec(node *n) {

    print_node_neighbours(n);

    if (n->next != NULL) {
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }

}

void print_utree(node *n) {

    print_node_neighbours(n);

    if (n->next != NULL) {
        print_utree_rec(n->back);
        print_utree_rec(n->next->back);
        print_utree_rec(n->next->next->back);
    }

}

void utree_nodes_within_radius(node *start_node, node *new_node, int radius,
                               std::vector<utree_move_info> &list_nodes) {

    /*
    utree_move_info m;
    m.node1 = start_node;
    m.node2 = new_node;
    list_nodes.push_back(m);

    if (radius <= 0) {
        return;
    }

    if (new_node->next!=NULL){
        radius--;
        utree_nodes_within_radius(start_node,new_node->next->back,radius,list_nodes);
        utree_nodes_within_radius(start_node,new_node->next->next->back,radius,list_nodes);
    }
    */

    utree_move_info m;
    m.node1 = start_node;
    if (new_node->next != NULL) {
        new_node = new_node->next;
    }
    m.node2 = new_node;
    list_nodes.push_back(m);

    if (radius <= 0) {
        return;
    }

    if (new_node->next != NULL) {
        radius--;
        utree_nodes_within_radius(start_node, new_node->back, radius, list_nodes);
        utree_nodes_within_radius(start_node, new_node->next->back, radius, list_nodes);
    }

}

void utree_get_list_nodes_within_radius(node *n,
                                        int radius,
                                        std::vector<utree_move_info> &list_nodes_left,
                                        std::vector<utree_move_info> &list_nodes_right,
                                        std::vector<utree_move_info> &list_nodes_up) {

    if (n->next != NULL) {
        utree_nodes_within_radius(n, n->back, radius, list_nodes_up);
        utree_nodes_within_radius(n, n->next->back, radius, list_nodes_left);
        utree_nodes_within_radius(n, n->next->next->back, radius, list_nodes_right);
    }

}

void copy_vector(std::vector<node *> &dest, std::vector<node *> &source) {

    for (unsigned int i = 0; i < source.size(); i++) {
        node *n = new node;
        node *m = source.at(i);
        n->next = m->next;
        n->back = m->back;
        n->data = m->data;
        n->ID = m->ID;
        dest.push_back(n);
    }

}
