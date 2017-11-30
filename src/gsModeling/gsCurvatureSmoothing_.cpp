/** @file gsCurvatureSmoothing_.cpp

    @brief Computes a closed B-spline curve with a smaller number of curvature
    extrema compared to a given closed B-spline curve  i.e. some kind of
    smoothing the curvature of the curve. This smoothting can be done with the
    help of two methods - total variation and Hadenfeld's algorithm (see Jan
    Hadenfeld, Iteratives Glätten von B-Spline Kurven und B-Spline Flächen,
    Shaker Verlag, PhD Thesis)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl
*/

#include <gsModeling/gsCurvatureSmoothing.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsCurvatureSmoothing<real_t>;

} // namespace gismo
