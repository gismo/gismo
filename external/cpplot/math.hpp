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

#ifndef _LIB_CPPLOT_MATH_
#define _LIB_CPPLOT_MATH_

#include <vector>
#include <algorithm>
#include <limits>

namespace cpplot {
    typedef std::vector<double> dvec; ///< General vector of doubles
    typedef std::vector< std::vector<double> > dmat; ///< General matrix of doubles
    typedef std::vector< std::vector<float> > tcvec; ///< Matrix of floats
    typedef std::vector< std::vector< std::vector<float> > > tcmat; ///< 4D tensor of floats

    namespace math {
        enum scale { linear_scale, logarithmic_scale}; ///< Indicates the scale of an axis

        /// Get the maximum valued element from vector
        template<typename T>
        T max(const std::vector<T>& x) {
            return *std::max_element(x.begin(), x.end());
        }

        /// Get the minimum valued element from vector
        template<typename T>
        T min(const std::vector<T>& x) {
            return *std::min_element(x.begin(), x.end());
        }

        /// Get the maximum valued element from a matrix
        template<typename T>
        T max(const std::vector<std::vector<T> >& x) {
            double res = std::numeric_limits<double>::min();
            for(unsigned int i = 0; i < x.size(); ++i) {
                res = std::max(res, max(x[i]));
            }
            return res;
        }

        /// Get the minimum valued element from matrix
        template<typename T>
        T  min(const std::vector<std::vector<T> >& x) {
            double res = std::numeric_limits<double>::max();
            for(unsigned int i = 0; i < x.size(); ++i) {
                res = std::min(res, min(x[i]));
            }
            return res;
        }

        /**
         * Generate n equidistant numbers between [min, max]
         */
        template<typename T>
        T linspace(const double min, const double max, unsigned int n) {
            if(n<1) { n = 1; }
            T a(n);
            for(unsigned int i = 0; i < n; ++i) {
                a[i] = min + (max-min)*i/(n-1);
            }
            return a;
        }

        dmat peaks(const int n); ///<
    }
}

#endif
