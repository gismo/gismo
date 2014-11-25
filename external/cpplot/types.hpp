/**
 * Copyright 2011, 2012 Jonatan Olofsson
 *
 * This file is part of cpplot.
 *
 * cpplot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cpplot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cpplot.  If not, see <http://www.gnu.org/licenses/>.
 */

/****************************************************************************
License: Gnu Public license (GPL) v3
* Author: Jonatan Olofsson (jonatan.olofsson@gmail.com)
* Version: 0.1
* Based on
Author: Yuichi Katori (yuichi.katori@gmail.com)
Project:MATPLOT++ (MATLAB-like plotting tool in C++).
Version:0.3.13
****************************************************************************/


/**
 * In this file, the main types of cpplot are defined
 */

#ifndef _CPPLOT_TYPES_HPP_
#define _CPPLOT_TYPES_HPP_

#include <boost/shared_ptr.hpp>

namespace cpplot {
    class axes_t_t;
    /**
     * The axes_t is a shared pointer that will keep the
     * axes object it's pointing to alive and available.
     */
    typedef boost::shared_ptr<axes_t_t> axes_t;
    typedef std::map<int, axes_t> axess_t;

    class layer_t_t;
    /**
     * The figure_t is a shared pointer that will keep the
     * layer object it's pointing to alive and available.
     */
    typedef boost::shared_ptr<layer_t_t> layer_t;
    typedef std::map<std::string, layer_t> layers_t;

    class figure_t_t;
    /**
     * The figure_t is a shared pointer that will keep the
     * figure object it's pointing to alive and available.
     */
    typedef boost::shared_ptr<figure_t_t> figure_t;
    typedef std::map<int, figure_t> figures_t;
}

#include "drawing.hpp"

#endif
