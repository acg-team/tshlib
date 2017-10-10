//===================================================================================================================
//===================================================================================================================
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <limits>
#include "PhyTree.h"
#include "newick.h"
#include "main.h"
//===================================================================================================================
//===================================================================================================================
int main(int argc, char** argv)
{

	std::string tree_file="/home/max/PIP_C++/NNI_SPR/tree_5_leaves_r_bl.nwk";
	PhyTree* tree = NULL;

	std::ifstream tree_str(tree_file.c_str());
	tree = newick_parser::parse_newick(&tree_str);

	tree->set_missing_node_name("V");

	std::cout<<tree->formatNewick()<<"\n\n";

	tree->print();
	std::cout<<"\n";

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

	t1=tree->get_right_child();
	t2=tree->get_left_child()->get_right_child();

//	tree->swap(tree,1,tree->get_left_child(),1);
	tree->swap2(t1,t2);

	std::cout<<tree->formatNewick()<<"\n\n";
	//----------------------------------------------------------

	return 0;
}





