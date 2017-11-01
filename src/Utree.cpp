//
// Created by Lorenzo Gatti on 26/10/17.
//

#include <string>
#include "Utree.hpp"
#include "PhyTree.hpp"


/*

void create_unrooted_tree(std::vector<node *> &utree,PhyTree *tree,node *parent){


    if(tree->isLeaf()){

        node *n = new node;

        parent->back=n;
        n->next=NULL;
        n->back=parent;
        n->data=tree;

        utree.push_back(n);
    }else{

        node *node_1= new node;
        node *node_2= new node;
        node *node_3= new node;

        parent->back=node_1;

        node_1->next=node_2;
        node_1->back=parent;
        node_1->data=tree;

        node_2->next=node_3;
        node_2->data=tree;

        node_3->next=node_1;
        node_3->data=tree;

        utree.push_back(node_1);
        utree.push_back(node_2);
        utree.push_back(node_3);

        create_unrooted_tree(utree,tree->children[0],node_2);
        create_unrooted_tree(utree,tree->children[1],node_3);
    }

}
*/

void cascade(VirtualNode *parent, VirtualNode *child) {

    if (!parent->isTerminalNode()) {

        auto ichild = new VirtualNode();
        cascade(child, ichild);

    }

}

void createUtree(PhyTree *in_tree, Utree *out_tree) {


    // For each node descending the root, create either a new VirtualInternalNode or a VirtualLeaf
    for (int i = 0; i < in_tree->get_children().at(0)->n_children(); i++) {

        auto inode = new VirtualNode;

        cascade(out_tree->utree_start_node, inode);


        out_tree->utree_start_node->addMember(inode);
    }
}


VirtualNode::VirtualNode() {

}

void VirtualNode::addMember(VirtualNode *inNode) {

    this->vnode_up;
    this->vnode_right;
    this->vnode_left;

}

bool VirtualNode::isTerminalNode() {
    return this->vnode_leaf;
}
