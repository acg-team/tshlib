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
 * @file Move.cpp
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

#include <numeric>
#include <limits>
#include <iomanip>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <include/Utilities.hpp>
#include <include/Utree.hpp>

#include "Move.hpp"


namespace tshlib {
    Move::~Move() = default;

    Move::Move() = default;


    void Move::initMove() {
        moveUID_ = 0;
        moveName_ = "undef";
        moveType_ = MoveType::undef;
        moveClass_ = "undef";
        moveDirection_ = MoveDirections::undef;

        moveScore_ = -std::numeric_limits<double>::infinity();

        moveApplied_ = false;
    }



    void Move::deleteTargetNode() {

        moveTargetNode_ = nullptr;

    }




    void Move::setClass(int Value) {

        if (Value >= 4 and Value < 10) {
            moveType_ = MoveType::SPR;
        } else if (Value >= 10) {
            moveType_ = MoveType::TBR;
        } else if (Value == 3) {
            moveType_ = MoveType::NNI;
        } else {
            moveType_ = MoveType::undef;
        }

        moveClass_ = getClass();
    }







}