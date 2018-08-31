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
 * @file Utilities.hpp
 * @author Lorenzo Gatti
 * @date 10 11 2017
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
 */

#ifndef TSHLIB_UTILITIES_H
#define TSHLIB_UTILITIES_H
namespace tshlib {
    enum class NodeRotation {
        undef = 0, clockwise = 1, counterclockwise = 2
    };

    enum class NodePosition {
        left, right, up, undef
    };

    enum class MoveDirections {
        down_left, down_right, up, up_right, up_left, undef
    };

    enum class MoveType {
        NNI = 1, SPR = 2, TBR = 3, FNNI = 4, VFNNI = 5, undef = 0
    };

    enum class TreeSearchStopCondition {
        iterations, convergence
    };

    enum class TreeSearchHeuristics {
        swap, phyml, mixed, nosearch
    };

    enum class StartingNodeHeuristics {
        particle_swarm, hillclimbing, greedy, undef
    };

    enum class TreeRearrangmentOperations {
        classic_NNI, classic_SPR, classic_TBR, classic_Mixed
    };

}
#endif //TSHLIB_UTILITIES_H
