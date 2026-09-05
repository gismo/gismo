/** @file gsWriteParasolid.hpp

    @brief Provides implementation of gsWriteParasolid function.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mokris
*/

#include <gsParasolid/gsPKUtils.h>

#include <gsParasolid/gsFrustrum.h>

#include <gsCore/gsLinearAlgebra.h>
#include <gsDomain/gsKnotVector.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsBSpline.h>

namespace gismo {

namespace extensions {

bool validMultiplicities(const std::vector<index_t>& mult,
                         const int deg)
{
    for (size_t i = 1; i != mult.size() - 1; i++)
    {
        if (mult[i] == deg + 1)
        {
            std::cout <<
                "Only B-Splines with (inner) multiplicity less (<) "
                "than degree + 1 are supported. \n"
                "Parasolid restriction."
                      << std::endl;
            return false;
        }
    }
    return true;
}

template<class T>
bool createPK_BSURF(const gsTensorBSpline< 2, T> & bsp,
                    PK_BSURF_t & bsurf,
                    bool closed_u,
                    bool closed_v)
{
    typedef typename gsKnotVector<T>::mult_t mult_t;
    std::vector<mult_t> mult;

    // Gismo and Parasolid store the coefficients in different order.
    // Therefore we create the BSURF with u and v swapped and then transpose it.
    for (index_t dim = 0; dim != 2; dim++)
    {
        const int deg = bsp.basis().degree(dim);
        mult = bsp.basis().knots(dim).multiplicities();

        if (!validMultiplicities(mult, deg))
        {
            return false;
        }
    }

    // Translate to parasolid standard form, ie fill up parasolid
    // spline data record
    PK_BSURF_sf_t sform;   // B-spline data holder (standard form)

    // Degrees
    sform.u_degree      = bsp.basis().degree(1);
    sform.v_degree      = bsp.basis().degree(0);

    // Knots in u-direction
    std::vector<T> gknot0 = bsp.basis().knots(1).unique();
    std::vector<mult_t> gmult0 = bsp.basis().knots(1).multiplicities();
    sform.n_u_knots     = gknot0.size();
    sform.u_knot        = gknot0.data();
    sform.u_knot_mult   = gmult0.data();

    // Knots in v-direction
    std::vector<T> gknot1 = bsp.basis().knots(0).unique();
    std::vector<mult_t> gmult1 = bsp.basis().knots(0).multiplicities();
    sform.n_v_knots     = gknot1.size();
    sform.v_knot        = gknot1.data();
    sform.v_knot_mult   = gmult1.data();

    // Control points
    sform.n_u_vertices  = bsp.basis().size(1);
    sform.n_v_vertices  = bsp.basis().size(0);
    gsMatrix<T> coefs   = bsp.coefs();
    const int n = bsp.geoDim();
    if ( n < 3 )
    {
        coefs.conservativeResize(gsEigen::NoChange, 3);
        coefs.rightCols(3-n).setZero();
    }
    coefs.transposeInPlace();
    coefs.resize(3*bsp.basis().size(), 1);
    sform.vertex_dim    = 3; // always 3 for surfaces
    sform.vertex        = coefs.data();

    // Attributes
    sform.is_rational   = PK_LOGICAL_false;
    sform.form          = PK_BSURF_form_unset_c;
    sform.u_knot_type   = PK_knot_unset_c;
    sform.v_knot_type   = PK_knot_unset_c;
    sform.is_u_periodic = PK_LOGICAL_false;
    sform.is_v_periodic = PK_LOGICAL_false;
    sform.is_u_closed   = PK_LOGICAL_false;
    sform.is_v_closed   = PK_LOGICAL_false;

    if (closed_u)
    {
        sform.is_v_closed = PK_LOGICAL_true;
    }

    if (closed_v)
    {
        sform.is_u_closed = PK_LOGICAL_true;
    }

    sform.self_intersecting = PK_self_intersect_unset_c;
    sform.convexity         = PK_convexity_unset_c;

    // Create parasolid surface with the previous spline data
    PK_ERROR_code_t err = PK_BSURF_create(&sform, &bsurf);
    PARASOLID_ERROR(PK_BSURF_create, err);

    // Transposition (new on 2019-02-26).
    PK_BSURF_reparameterise_o_t options;
    PK_BSURF_reparameterise_o_m(options);
    options.transpose = PK_LOGICAL_true;

    err = PK_BSURF_reparameterise(bsurf,&options);
    PARASOLID_ERROR(PK_BSURF_reparameterise, err);

    return true;
}


template<class T>
bool createPK_BCURVE( const gsBSpline<T>& curve,
                      PK_BCURVE_t& bcurve)
{
    typedef typename gsKnotVector<T>::mult_t mult_t;
    
    PK_BCURVE_sf_t sform; // B-curve data holder (standard form)

    // Degree
    sform.degree = curve.degree();

    // Knots
    std::vector<T> knots = curve.basis().knots().unique();
    std::vector<mult_t> mult = curve.basis().knots().multiplicities();
    sform.n_knots = knots.size();
    sform.knot = knots.data();
    sform.knot_mult = mult.data();


    // Control points
    sform.n_vertices = curve.basis().size();
    gsMatrix<T> coefs = curve.coefs();
    const int n = curve.geoDim();
    if (n < 3)
    {
        coefs.conservativeResize(gsEigen::NoChange, 3);
        coefs.rightCols(3 - n).setZero();
    }
    coefs.transposeInPlace();
    coefs.resize(3 * curve.basis().size(), 1);
    sform.vertex_dim = 3;
    sform.vertex = coefs.data();

    // Attributes
    sform.is_rational = PK_LOGICAL_false;
    sform.form = PK_BCURVE_form_unset_c;
    sform.knot_type = PK_knot_unset_c;
    sform.is_periodic = PK_LOGICAL_false;
    sform.is_closed = PK_LOGICAL_false;
    sform.self_intersecting = PK_self_intersect_unset_c;

    PK_ERROR_code_t err = PK_BCURVE_create(&sform, &bcurve);
    PARASOLID_ERROR(PK_BCURVE_create, err);

    return true;
}


template<class T>
bool createPK_GEOM( const gsGeometry<T> & ggeo,
                    PK_GEOM_t & pgeo)
{
    // Identify input gismo geometry
    if ( const gsTensorBSpline<2,T> * tbsp =
         dynamic_cast<const gsTensorBSpline<2,T> *>(&ggeo) )
    {
        return createPK_BSURF(*tbsp, pgeo);
    }
// the following lines produce warnings if called from multipatch version of gsWriteParasolid,
// because it already assumes that the geometries are surfaces
    else if ( const gsBSpline<>* bspl =
              dynamic_cast< const gsBSpline<>* >(&ggeo) )
    {
        return createPK_BCURVE(*bspl, pgeo);
    }
    else
    {
        gsInfo << "Cannot write "<<ggeo<<" to parasolid file.\n";
        return false;
    }
}

} // namespace extensions

} // namespace gismo
