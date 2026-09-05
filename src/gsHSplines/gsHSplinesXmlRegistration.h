/** @file gsHSplinesXmlRegistration.h

    @brief Registers the gsHSplines XML readers/writers in gsXmlRegistry.

    The only file to touch when gsHSplines gains a serializable type.
    Priorities reproduce the historical dispatch order of
    gsIO/gsXmlUtils.hpp exactly, including its asymmetries (gsHBSpline<1>
    was writable but never readable through the abstract bases).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <gsCore/gsCurve.h>
#include <gsCore/gsSurface.h>
#include <gsHSplines/gsTHBSpline.h> // also provides the gsHBSpline alias
#include <gsHSplines/gsRationalTHBSpline.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsRationalTHBSplineBasis.h>

namespace gismo {
namespace internal {

template <class T>
void gsHSplinesXmlRegisterTypes()
{
    typedef gsGeometry<T> G;
    typedef gsCurve<T>    C;
    typedef gsSurface<T>  S;

    gsXmlRegisterAs<G, gsTHBSpline<1,T> >         ::enroll(180);
    gsXmlRegisterAs<G, gsTHBSpline<2,T> >         ::enroll(190);
    gsXmlRegisterAs<G, gsTHBSpline<3,T> >         ::enroll(200);
    // (210 is gsTrimSurface, registered by gsModeling)
    gsXmlRegisterAs<G, gsHBSpline<1,T> >          ::enrollPut(220); // write-only historically
    gsXmlRegisterAs<G, gsHBSpline<2,T> >          ::enroll(230);
    gsXmlRegisterAs<G, gsHBSpline<3,T> >          ::enroll(240);
    gsXmlRegisterAs<G, gsRationalTHBSpline<1,T> > ::enroll(250);
    gsXmlRegisterAs<G, gsRationalTHBSpline<2,T> > ::enroll(260);
    gsXmlRegisterAs<G, gsRationalTHBSpline<3,T> > ::enroll(270);

    gsXmlRegisterAs<C, gsHBSpline<1,T> >          ::enrollPut(120); // write-only historically

    gsXmlRegisterAs<S, gsTHBSpline<2,T> >         ::enroll(120);
    gsXmlRegisterAs<S, gsHBSpline<2,T> >          ::enroll(130);

    typedef gsBasis<T> B;
    gsXmlRegisterAs<B, gsHTensorBasis<1,T> >        ::enrollPut(180);
    gsXmlRegisterAs<B, gsHTensorBasis<2,T> >        ::enrollPut(190);
    gsXmlRegisterAs<B, gsHTensorBasis<3,T> >        ::enrollPut(200);
    gsXmlRegisterAs<B, gsHTensorBasis<4,T> >        ::enrollPut(210);
    gsXmlRegisterAs<B, gsTHBSplineBasis<1,T> >      ::enrollGet();
    gsXmlRegisterAs<B, gsTHBSplineBasis<2,T> >      ::enrollGet();
    gsXmlRegisterAs<B, gsTHBSplineBasis<3,T> >      ::enroll(220);
    gsXmlRegisterAs<B, gsTHBSplineBasis<4,T> >      ::enrollGet();
    gsXmlRegisterAs<B, gsHBSplineBasis<1,T> >       ::enrollGet();
    gsXmlRegisterAs<B, gsHBSplineBasis<2,T> >       ::enrollGet();
    gsXmlRegisterAs<B, gsHBSplineBasis<3,T> >       ::enrollGet();
    gsXmlRegisterAs<B, gsHBSplineBasis<4,T> >       ::enrollGet();
    gsXmlRegisterAs<B, gsRationalTHBSplineBasis<1,T> >::enrollGet();
    gsXmlRegisterAs<B, gsRationalTHBSplineBasis<2,T> >::enrollGet();
    gsXmlRegisterAs<B, gsRationalTHBSplineBasis<3,T> >::enroll(230);
}

#ifndef GISMO_BUILD_LIB
namespace {
struct gsHSplinesXmlReg
{ gsHSplinesXmlReg() { gsHSplinesXmlRegisterTypes<real_t>(); } };
static const gsHSplinesXmlReg s_gsHSplinesXmlReg;
}
#endif

} // namespace internal
} // namespace gismo
