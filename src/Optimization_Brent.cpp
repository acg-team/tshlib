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
 * @file Optimization_Brent.cpp
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
#include <cmath>
#include "Optimization_Brent.hpp"
#include "Likelihood.hpp"

namespace tshlib {
    int sgn(double d) {
        double eps = 1e-10;
        return d < -eps ? -1 : d > eps;
    }

    bool Generic_Brent_Lk(double *param, double ax, double cx, double tol,
                          int n_iter_max, double (*obj_func)(Likelihood &, std::vector<VirtualNode *> &, Alignment &),
                          Likelihood &likelihood, std::vector<VirtualNode *> &list_vnode_to_root, Alignment &alignment, double &start_LK) {
        int iter;
        double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
        double e = 0.0;
        double old_lnL, init_lnL, cur_lnL;
        double bx = *param;

        d = 0.0;
        a = ((ax < cx) ? ax : cx);
        b = ((ax > cx) ? ax : cx);
        x = w = v = bx;
        old_lnL = UNLIKELY;
        (*param) = fabs(bx);
        fw = fv = fx = fu = -start_LK;
        init_lnL = fw;

        old_lnL = -start_LK;

        for (iter = 1; iter <= BRENT_ITMAX; iter++) {
            xm = 0.5 * (a + b);
            tol2 = 2.0 * (tol1 = tol * fabs(x) + BRENT_ZEPS);

            //cur_lnL = (stree) ? (stree->tree->c_lnL) : (tree->c_lnL);
            cur_lnL = old_lnL;


            if (fabs(e) > tol1) {
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) p = -p;
                q = fabs(q);
                etemp = e;
                e = d;
                if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                    d = BRENT_CGOLD * (e = (x >= xm ? a - x : b - x));
                } else {
                    d = p / q;
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2) d = sgn(xm - x);
                }
            } else {
                d = BRENT_CGOLD * (e = (x >= xm ? a - x : b - x));
            }

            u = (fabs(d) >= tol1 ? x + d : x + sgn(d));
            (*param) = fabs(u);

            //old_lnL = (stree) ? (stree->tree->c_lnL) : (tree->c_lnL);
            //old_lnL = 0;
            fu = -(*obj_func)(likelihood, list_vnode_to_root, alignment);

            cur_lnL = fu;
            // Did I reach convergence or max iteration? If yes, stop it!
            if ((fabs(cur_lnL - old_lnL) < tol) || (iter > n_iter_max - 1)) {
                (*param) = x;
                //return (*obj_func)(branch, tree);
                //return (stree) ? (stree->tree->c_lnL) : (tree->c_lnL);
                start_LK = -cur_lnL;
                return true;
            }

            if (fu <= fx) {
                if (u >= x) a = x; else b = x;
                SHFT(v, w, x, u)
                SHFT(fv, fw, fx, fu)
            } else {
                if (u < x) a = u; else b = u;
                if (fu < fw || fabs(w - x) < SMALL) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                } else if (fu < fv || fabs(v - x) < SMALL || fabs(v - w) < SMALL) {
                    v = u;
                    fv = fu;
                }
            }

            old_lnL = fu;
        }
        // TODO: it should return a boolean
        return (false);
        /* Not Reached ??  *param=x;   */
        /* Not Reached ??  return fx; */
    }
}
