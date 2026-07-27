/** @file gsNurbsXmlRegistration.h

    @brief Registers the gsNurbs XML readers/writers in gsXmlRegistry.

    The only file to touch when gsNurbs gains a serializable type.
    Priorities reproduce the historical dispatch order of
    gsIO/gsXmlUtils.hpp exactly (see modularization report S3).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <gsCore/gsCurve.h>
#include <gsCore/gsSurface.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorNurbs.h>

namespace gismo {
namespace internal {

template <class T>
void gsNurbsXmlRegisterTypes()
{
    typedef gsGeometry<T> G;
    typedef gsCurve<T>    C;
    typedef gsSurface<T>  S;

    gsXmlRegisterAs<G, gsBSpline<T> >               ::enroll(100);
    gsXmlRegisterAs<G, gsNurbs<T> >                 ::enroll(110);
    gsXmlRegisterAs<G, gsTensorBSpline<2,T> >       ::enroll(120);
    gsXmlRegisterAs<G, gsTensorBSpline<3,T> >       ::enroll(130);
    gsXmlRegisterAs<G, gsTensorBSpline<4,T> >       ::enroll(140);
    gsXmlRegisterAs<G, gsTensorBSpline<1,T> >       ::enrollGet(); // historically read-only
    gsXmlRegisterAs<G, gsTensorNurbs<2,T> >         ::enroll(150);
    gsXmlRegisterAs<G, gsTensorNurbs<3,T> >         ::enroll(160);
    gsXmlRegisterAs<G, gsTensorNurbs<4,T> >         ::enroll(170);

    gsXmlRegisterAs<C, gsBSpline<T> >               ::enroll(100);
    gsXmlRegisterAs<C, gsNurbs<T> >                 ::enroll(110);

    gsXmlRegisterAs<S, gsTensorBSpline<2,T> >       ::enroll(100);
    gsXmlRegisterAs<S, gsTensorNurbs<2,T> >         ::enroll(110);
}

#ifndef GISMO_BUILD_LIB
namespace {
struct gsNurbsXmlReg
{ gsNurbsXmlReg() { gsNurbsXmlRegisterTypes<real_t>(); } };
static const gsNurbsXmlReg s_gsNurbsXmlReg;
}
#endif

} // namespace internal
} // namespace gismo
