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
 * @file Move.hpp
 * @author Lorenzo Gatti
 * @date 11 06 2018
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
#ifndef TSHLIB_MOVE_HPP
#define TSHLIB_MOVE_HPP

#include "Utilities.hpp"
#include <Utree.hpp>

namespace tshlib {
    class Move {

    private:

    protected:
        VirtualNode *moveTargetNode_;       /* Pointer to the target node found during the node search */
        VirtualNode *moveSourceNode_;       /* Pointer to the source node  */
        VirtualNode *moveStepParentNode_;   /* Pointer to the node that acts as step parent during the SPR move */
        VirtualNode *moveStepChildNode_;    /* Pointer to the node that acts as step parent during the SPR move */
        bool onPseudoRoot_ ;                /* If the SPR node is infact a pseudoroot */


    public:
        int moveUID_;                           /* Move UID - Useful in case of parallel independent executions*/
        std::string moveName_;                   /* Move Name - Unused */
        int moveRadius_;                        /* Move Radius */
        MoveDirections moveDirection_;          /* Move Direction for applying a rotation to the VirtualNode pointers */
        double moveScore_;                      /* Likelihood of the move if applied */
        bool moveApplied_;                      /* Indicator is set to true if the move is applied to the tree */
        std::string moveClassDescription_;                 /* String indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */
        MoveType moveType_;                     /* Integer indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */
        TreeSearchHeuristics moveStrategy_;     /* Store the strategy used to generate this candidate.

        /*!
         * @brief Standard constructor
         */
        Move();

        /*!
        * @brief Standard deconstructor
        */
        ~Move();

        Move(const Move &inMove) {

            moveUID_ = inMove.moveUID_;
            moveName_ = "copy_" + inMove.moveName_;
            moveRadius_ = inMove.moveRadius_;
            moveScore_ = inMove.moveScore_;
            moveApplied_ = inMove.moveApplied_;
            moveClassDescription_ = inMove.moveClassDescription_;
            moveType_ = inMove.moveType_;
            moveDirection_ = inMove.moveDirection_;
            moveTargetNode_ = inMove.moveTargetNode_;
            moveSourceNode_ = inMove.moveSourceNode_;
            moveStepParentNode_ = inMove.moveStepParentNode_;

        }


        Move &operator=(const Move &inMove) {
            moveUID_ = inMove.moveUID_;
            moveName_ = "copy_" + inMove.moveName_;
            moveRadius_ = inMove.moveRadius_;
            moveScore_ = inMove.moveScore_;
            moveApplied_ = inMove.moveApplied_;
            moveClassDescription_ = inMove.moveClassDescription_;
            moveType_ = inMove.moveType_;
            moveDirection_ = inMove.moveDirection_;
            moveTargetNode_ = inMove.moveTargetNode_;
            moveSourceNode_ = inMove.moveSourceNode_;
            moveStepParentNode_ = inMove.moveStepParentNode_;

        };

        void initMove();


        std::string getDirection() const {
            std::string rtToken;
            switch (moveDirection_) {
                case MoveDirections::down_left:
                    rtToken = "left";
                    break;
                case MoveDirections::up:
                    rtToken = "up";
                    break;
                case MoveDirections::up_left:
                    rtToken = "up-left";
                    break;
                case MoveDirections::up_right:
                    rtToken = "up-right";
                    break;
                case MoveDirections::down_right :
                    rtToken = "right";
                    break;
                case MoveDirections::undef :
                    rtToken = "undef";
                    break;

            }

            return rtToken;
        }


        void setClass(TreeSearchHeuristics tsStrategy);


        std::string getClass() const {

            std::string rtToken;
            switch (moveType_) {

                case MoveType::VFNNI:
                    rtToken = "vfNNI";
                    break;

                case MoveType::FNNI:
                    rtToken = "fNNI";
                    break;

                case MoveType::NNI:
                    rtToken = "NNI";
                    break;

                case MoveType::SPR:
                    rtToken = "SPR";
                    break;

                case MoveType::TBR:
                    rtToken = "TBR";
                    break;

                case MoveType::undef:
                    rtToken = "undef";
                    break;
            }

            return rtToken;
        }

        /*!
         * @brief Reset the protected move_targetnode field
         */
        void deleteTargetNode();

        /*!
         * @brief Returns the target node pointer
         * @return VirtualNode pointer of the target node
         */
        VirtualNode *getTargetNode(){
            return moveTargetNode_;
        };

/*!
        * @brief Returns the source node pointer
        * @return VirtualNode pointer of the source node
        */
        VirtualNode *getSourceNode(){
            return moveSourceNode_;
        };

        /*!
         * @brief Set the protected move_targetnode field
         * @param target_node PhyTree Pointer to the target node
         */
        void setTargetNode(VirtualNode *target_node){
            moveTargetNode_ = target_node;
        };

        void setSourceNode(VirtualNode *source_node){
            moveSourceNode_ = source_node;
        };

        void setStepParentNode(VirtualNode *stepParent){
            moveStepParentNode_ = stepParent;
        };

        VirtualNode *getStepParentNode(){
            return moveStepParentNode_;
        };

        void setStepChildNode(VirtualNode *stepChild){
            moveStepChildNode_ = stepChild;
        };

        VirtualNode *getStepChildNode(){
            return moveStepChildNode_;
        };


        MoveType getType() const{
            return moveType_;
        };


        void setRadius(int radius){
            Move::moveRadius_ = radius;
        };

        int getRadius() const {
            return moveRadius_;
        }

        void setDirection(MoveDirections direction){
            Move::moveDirection_ = direction;
        };

        int getUID() const {
            return moveUID_;
        }

        void setUID(int moveUID_) {
            Move::moveUID_ = moveUID_;
        }

        const std::string &getName() const {
            return moveName_;
        }

        void setName(const std::string &moveName_) {
            Move::moveName_ = moveName_;
        }

        double getScore() const {
            return moveScore_;
        }

        void setScore(double moveScore_) {
            Move::moveScore_ = moveScore_;
        }

        bool isOverPseudoRoot_() const {
            return onPseudoRoot_;
        }

        void setOnPseudoRoot_(bool onPseudoRoot_) {
            Move::onPseudoRoot_ = onPseudoRoot_;
        }

        MoveDirections getMoveDirection() const {
            return moveDirection_;
        }


        TreeSearchHeuristics getMoveStrategy() const {
            return moveStrategy_;
        }

    };


}

#endif //TSHLIB_MOVE_HPP
