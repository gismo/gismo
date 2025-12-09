/** @file gsKnotVector_var.cpp

    @brief Explicit instantiations for gsKnotVector with autodiff::var type

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsCore/gsTemplateTools.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsKnotVector.hpp>

#include <autodiff/reverse/var.hpp>

namespace gismo
{
    // Explicit template instantiation for autodiff::var (reverse-mode AD)
    typedef autodiff::reverse::detail::Variable<double> var_t;
    
    TEMPLATE_INST gsKnotVector<var_t>;
}
