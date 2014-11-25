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

#include "math.hpp"

namespace cpplot {
    namespace math {
        dmat peaks(const int n) {
            float x1 = 1, y1 = 0;
            float x2 = -1, y2 = 1;
            float x3 = -0.5,y3 = -1;
            float sr1, sr2, sr3;
            float sigma = 0.4;
            float a1 = 1, a2 = 0.5, a3 = 0.3;
            double x, y;
            //vector< vector< double > > Z(n,n);
            dmat Z(n,dvec(n));
            for(int i = 0; i < n; ++i) {
                for(int j = 0; j < n; ++j) {
                    x = -2.0+4.0*j/(n-1);
                    y = -2.0+4.0*i/(n-1);
                    sr1 = (x-x1)*(x-x1)+(y-y1)*(y-y1);
                    sr2 = (x-x2)*(x-x2)+(y-y2)*(y-y2);
                    sr3 = (x-x3)*(x-x3)+(y-y3)*(y-y3);
                    Z[i][j] = a1/(1+sr1/sigma) +a2/(1+sr2/sigma) +a3/(1+sr3/sigma) ;
                }
            }
            return Z;
        }
    }
}
