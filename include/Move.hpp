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
#ifndef TSHLIB_MOVE_HPP
#define TSHLIB_MOVE_HPP

namespace tshlib {
    class Move {

    private:

    protected:
        VirtualNode *move_targetnode;   /* Pointer to the target node found during the node search */
        VirtualNode *move_sourcenode;   /* Pointer to the source node  */

    public:
        int move_id;                            /* Move ID - Useful in case of parallel independent executions*/
        std::string move_name;                  /* Move Name - Unused */
        int move_radius;                        /* Move Radius */
        MoveDirections move_direction;          /* Move Direction for applying a rotation to the VirtualNode pointers */
        double move_lk;                         /* Likelihood of the move if applied */
        bool move_applied;                      /* Indicator is set to true if the move is applied to the tree */
        std::string move_class;                 /* String indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */
        MoveType move_type;                     /* Integer indicating the move class (i.e. NNI,SPR,TBR) - Usefull in case of mixed tree-search strategies */

        /*!
         * @brief Standard constructor
         */
        Move();

        /*!
        * @brief Standard deconstructor
        */
        ~Move();

        Move(const Move &inMove) {

            move_id = inMove.move_id;
            move_name = "copy_" + inMove.move_name;
            move_radius = inMove.move_radius;
            move_lk = inMove.move_lk;
            move_applied = inMove.move_applied;
            move_class = inMove.move_class;
            move_type = inMove.move_type;
            move_direction = inMove.move_direction;
            move_targetnode = inMove.move_targetnode;
            move_sourcenode = inMove.move_sourcenode;

        }

        Move &operator=(const Move &inMove) {
            move_id = inMove.move_id;
            move_name = "copy_" + inMove.move_name;
            move_radius = inMove.move_radius;
            move_lk = inMove.move_lk;
            move_applied = inMove.move_applied;
            move_class = inMove.move_class;
            move_type = inMove.move_type;
            move_direction = inMove.move_direction;
            move_targetnode = inMove.move_targetnode;
            move_sourcenode = inMove.move_sourcenode;
        };


        std::string getMoveDirection() const {
            std::string rtToken;
            switch (move_direction) {
                case MoveDirections::left:
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
                case MoveDirections::right :
                    rtToken = "right";
                    break;
                case MoveDirections::undef :
                    rtToken = "undef";
                    break;

            }

            return rtToken;
        }

        std::string getMoveClass() const {

            std::string rtToken;
            switch (move_type) {

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
        VirtualNode *getTargetNode();

/*!
        * @brief Returns the source node pointer
        * @return VirtualNode pointer of the source node
        */
        VirtualNode *getSourceNode();

        /*!
         * @brief Set the protected move_targetnode field
         * @param target_node PhyTree Pointer to the target node
         */
        void setTargetNode(VirtualNode *target_node);

        void setSourceNode(VirtualNode *source_node);

        void setMoveClass(int Value);

        MoveType getMoveType() const;

        void setRadius(int radius);

        void setDirection(MoveDirections direction);

        double getLikelihood() { return move_lk; };

        void initMove();
    };
}

#endif //TSHLIB_MOVE_HPP
