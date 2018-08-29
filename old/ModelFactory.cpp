/*******************************************************************************
 * Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
 *
 *
 * Copyright (C) 2015-2018 by Lorenzo Gatti & Massimo Maiolo
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
 * @file Model.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 12 12 2017
 * @version 3.0.1
 * @maintainer Lorenzo Gatti
 * @maintainer Massimo Maiolo
 * @email lg@lorenzogatti.me
 * @email massimo.maiolo@zhaw.ch
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
#include "ModelFactory.hpp"

namespace tshlib {

    void Parameter::Optimise() {}

    void Parameter::UpdateList() {}

    Parameter::Parameter() = default;

    Parameter::Parameter(bool optimisable) : optimisable(optimisable) {}

    Parameter::~Parameter() = default;

    topology::topology(bool optimisable, Utree *tree) : Parameter(optimisable), value(tree) { name = "topology"; }

    msa::msa(bool optimisable, Alignment *alignment) : Parameter(optimisable), value(alignment) { name = "msa"; }

    Vi::Vi(bool optimisable, const Eigen::MatrixXd &Vi) : Parameter(optimisable), value(Vi) { name = "vi"; }

    pi::pi(bool optimisable, const Eigen::VectorXd &pi) : Parameter(optimisable), value(pi) { name = "pi"; }

    Q::Q(bool optimisable, const Eigen::MatrixXd &Q) : Parameter(optimisable), value(Q) { name = "Q"; }

    V::V(bool optimisable, const Eigen::MatrixXd &V) : Parameter(optimisable), value(V) { name = "V"; }

}