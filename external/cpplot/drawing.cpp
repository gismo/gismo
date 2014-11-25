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

#include "cpplot_common.hpp"

namespace cpplot {
    // coordinate transform  ///
    // axes coordination
    float drawing_t_t::ctx(const double x) {
        if(ca->XScale == math::linear_scale) {//linear_scale
            return ( x - ca->XLim[0] ) / (ca->XLim[1] - ca->XLim[0] );
        } else {//log
            return ( std::log10(x) - std::log10(ca->XLim[0]) )/( std::log10(ca->XLim[1]) - std::log10(ca->XLim[0]) );
        }
    }
    float drawing_t_t::cty(const double y) {
        if(ca->YScale == math::linear_scale) {//linear_scale
            return ( y - ca->YLim[0] ) / ( ca->YLim[1] - ca->YLim[0] );
        } else {//log
            return ( std::log10(y) - std::log10(ca->YLim[0]) ) / ( std::log10(ca->YLim[1]) - std::log10(ca->YLim[0]) );
        }
    }
    float drawing_t_t::ct3x(const double x) {
        return -1 + 2*( x - ca->XLim[0] ) / ( ca->XLim[1] - ca->XLim[0] );
    }
    float drawing_t_t::ct3y(const double y) {
        return -1 + 2*( y - ca->YLim[0] ) / ( ca->YLim[1] - ca->YLim[0] );
    }
    float drawing_t_t::ct3z(const double z) {
        return -1 + 2*( z - ca->ZLim[0] ) / ( ca->ZLim[1] - ca->ZLim[0] );
    }
}
