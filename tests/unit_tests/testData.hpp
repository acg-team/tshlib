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
 * @file testData.hpp
 * @author Lorenzo Gatti
 * @date 07 06 2018
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
#ifndef TSHLIB_TESTDATA_HPP
#define TSHLIB_TESTDATA_HPP

#include <Utree.hpp>
using namespace tshlib;

class testData{
public:

    testData() {}
    virtual ~testData() = default;

    Utree *getTree() {
        auto *testTree = new Utree();
        auto *nodeA = new VirtualNode();
        nodeA->setBranchLength(0.1);
        nodeA->setNodeName("A");
        nodeA->setVnode_id(1);
        testTree->addMember(nodeA);
        auto *nodeB = new VirtualNode();
        nodeB->setBranchLength(0.2);
        nodeB->setNodeName("B");
        nodeB->setVnode_id(2);
        testTree->addMember(nodeB);
        auto *nodeC = new VirtualNode();
        nodeC->setBranchLength(0.3);
        nodeC->setNodeName("C");
        nodeC->setVnode_id(3);
        testTree->addMember(nodeC);
        auto *nodeD = new VirtualNode();
        nodeD->setBranchLength(0.4);
        nodeD->setNodeName("D");
        nodeD->setVnode_id(4);
        testTree->addMember(nodeD);
        auto *nodeE = new VirtualNode();
        nodeE->setBranchLength(0.5);
        nodeE->setNodeName("E");
        nodeE->setVnode_id(5);
        testTree->addMember(nodeE);
        auto *nodeF = new VirtualNode();
        nodeF->setBranchLength(0.6);
        nodeF->setNodeName("F");
        nodeF->setVnode_id(6);
        testTree->addMember(nodeF);
        auto *nodeG = new VirtualNode();
        nodeG->setBranchLength(0.7);
        nodeG->setNodeName("G");
        nodeG->setVnode_id(7);
        testTree->addMember(nodeG);
        auto *nodeH = new VirtualNode();
        nodeH->setBranchLength(0.8);
        nodeH->setNodeName("H");
        nodeH->setVnode_id(8);
        testTree->addMember(nodeH);
        auto *nodeV1 = new VirtualNode();
        nodeV1->setBranchLength(0.8);
        nodeV1->setNodeName("V1");
        nodeV1->setVnode_id(9);
        testTree->addMember(nodeV1);
        auto *nodeV2 = new VirtualNode();
        nodeV2->setBranchLength(0.7);
        nodeV2->setNodeName("V2");
        nodeV2->setVnode_id(10);
        testTree->addMember(nodeV2);
        auto *nodeV3 = new VirtualNode();
        nodeV3->setBranchLength(0.6);
        nodeV3->setNodeName("V3");
        nodeV3->setVnode_id(11);
        testTree->addMember(nodeV3);
        auto *nodeV4 = new VirtualNode();
        nodeV4->setBranchLength(0.5);
        nodeV4->setNodeName("V4");
        nodeV4->setVnode_id(12);
        testTree->addMember(nodeV4);
        auto *nodeV5 = new VirtualNode();
        nodeV5->setBranchLength(0.4);
        nodeV5->setNodeName("V5");
        nodeV5->setVnode_id(13);
        testTree->addMember(nodeV5);
        auto *nodeV6 = new VirtualNode();
        nodeV6->setBranchLength(0.3);
        nodeV6->setNodeName("V6");
        nodeV6->setVnode_id(14);
        testTree->addMember(nodeV6);

        nodeV1->connectNode(nodeA);
        nodeV1->connectNode(nodeB);
        nodeV2->connectNode(nodeC);
        nodeV2->connectNode(nodeV1);
        nodeV3->connectNode(nodeD);
        nodeV3->connectNode(nodeV2);
        nodeV6->connectNode(nodeG);
        nodeV6->connectNode(nodeH);
        nodeV5->connectNode(nodeF);
        nodeV5->connectNode(nodeV6);
        nodeV4->connectNode(nodeE);
        nodeV4->connectNode(nodeV5);
        nodeV4->connectNode(nodeV3);

        return testTree;
    }

};


#endif //TSHLIB_TESTDATA_HPP
