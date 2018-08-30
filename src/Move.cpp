/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti
 *******************************************************************************
 *
 * This file is part of tshlib
 *
 * Tree Search Heuristic Library (TshLib) is a free software: you can redistribute
 * it and/or modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tree Search Heuristic Library (TshLib) is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with TshLib. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Move.cpp
 * @author Lorenzo Gatti
 * @date 11 06 2018
 * @version 3.0.1
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
#include <Utilities.hpp>


#include "Move.hpp"


namespace tshlib {
    Move::~Move() = default;

    Move::Move() = default;


    void Move::initMove() {
        moveUID_ = 0;
        moveName_ = "undef";
        moveType_ = MoveType::undef;
        moveClassDescription_ = "undef";
        moveDirection_ = MoveDirections::undef;

        moveScore_ = -std::numeric_limits<double>::infinity();

        moveApplied_ = false;
    }


    void Move::deleteTargetNode() {

        moveTargetNode_ = 0;

    }


    void Move::setClass(TreeSearchHeuristics tsStrategy, bool _location__overpseudoroot) {

        // Set tree search strategy associated to this move
        moveStrategy_ = tsStrategy;

        // Get radius of the current move
        int radius = getRadius();

        if (radius == 3) {
            switch (tsStrategy) {
                case TreeSearchHeuristics::swap:
                case TreeSearchHeuristics::phyml:
                case TreeSearchHeuristics::mixed:
                    moveType_ = MoveType::NNI;
                    break;
                case TreeSearchHeuristics::nosearch:
                    moveType_ = MoveType::undef;
                    break;
            }
        } else if (radius == 4) {
            switch (tsStrategy) {
                case TreeSearchHeuristics::swap:
                    moveType_ = MoveType::FNNI;
                    break;
                case TreeSearchHeuristics::phyml:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::mixed:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::nosearch:
                    moveType_ = MoveType::undef;
                    break;
            }
        } else if (radius > 4) {
            switch (tsStrategy) {
                case TreeSearchHeuristics::swap:
                    moveType_ = MoveType::VFNNI;
                    break;
                case TreeSearchHeuristics::phyml:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::mixed:
                    moveType_ = MoveType::SPR;
                    break;
                case TreeSearchHeuristics::nosearch:
                    moveType_ = MoveType::undef;
                    break;
            }
        } else {

            moveType_ = MoveType::undef;
        }

        // If either the source or the target node define are on the pseudoroot, and the tree-search strategy is [PhyML], then the
        // movetype is handled as TBR
        //if ((moveTargetNode_->isPseudoRootNode() || moveSourceNode_->isPseudoRootNode()) && tsStrategy == TreeSearchHeuristics::phyml)
        if (_location__overpseudoroot && tsStrategy == TreeSearchHeuristics::phyml)
            moveType_ = MoveType::TBR;


        moveClassDescription_ = getClass();
    }


}