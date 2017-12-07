/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti
 *
 *
 * Copyright (C) 2015-2017 by Lorenzo Gatti
 *******************************************************************************
 *
 * This file is part of tshexe
 *
 * tshexe is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tshexe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with tshexe. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file Optimization_Brent.hpp
 * @author Lorenzo Gatti
 * @date 06 12 2017
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
#ifndef TSHEXE_OPTIMIZATION_BRENT_HPP
#define TSHEXE_OPTIMIZATION_BRENT_HPP

#include <cfloat>
#include "Utree.hpp"
#include "Likelihood.hpp"

#define  BRENT_ITMAX        10000
#define  BRENT_CGOLD    0.3819660
#define  BRENT_ZEPS        1.e-10
#define SHFT(a, b, c, d)                (a)=(b);(b)=(c);(c)=(d);
#define SMALL DBL_MIN
#define  UNLIKELY          -1.e10

int sgn(double d);

double Generic_Brent_Lk(double *param, double ax, double cx, double tol,
                        int n_iter_max, double (*obj_func)(Likelihood &, std::vector<VirtualNode *> &, Alignment &),
                        Likelihood &likelihood, std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment, double start_LK);

#endif //TSHEXE_OPTIMIZATION_BRENT_HPP
