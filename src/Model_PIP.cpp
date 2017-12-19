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
 * @file Model_PIP.cpp
 * @author Lorenzo Gatti
 * @author Massimo Maiolo
 * @date 18 12 2017
 * @version 1.0
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
#include "Model_PIP.hpp"

pip_Lambda::pip_Lambda(bool optimisable, double value) : Parameter(optimisable), value(value) { name = "pip_lambda";}

pip_Mu::pip_Mu(bool optimisable, double value) : Parameter(optimisable), value(value) {name = "pip_mu";}

pip_Phi::pip_Phi(bool optimisable, double value) : Parameter(optimisable), value(value) {name = "pip_phi";}

pip_Tau::pip_Tau(bool optimisable, double value) : Parameter(optimisable), value(value) {name = "pip_Tau";}

pip_Nu::pip_Nu(bool optimisable, double value) : Parameter(optimisable), value(value) {name = "pip_Nu";}


