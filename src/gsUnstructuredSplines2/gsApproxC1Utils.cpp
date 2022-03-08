/** @file gsApproxC1Utils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsUnstructuredSplines2/gsApproxC1Utils.h>


namespace gismo
{
    void createGluingDataSpace(const gsGeometry<real_t> & patch, const gsBasis<real_t> & basis, index_t dir,
                               gsBSplineBasis<real_t> & result, index_t p_tilde, index_t r_tilde)
    {
        const gsBSplineBasis<real_t> basis_1 = dynamic_cast<const gsBSplineBasis<real_t> &>(basis.component(dir));

        if (p_tilde == -1)
            p_tilde = basis_1.degree() - 1;
        if (r_tilde == -1)
            r_tilde = p_tilde - 1;

        gsKnotVector<real_t> kv_gluingData(basis_1.knots().unique(), p_tilde, r_tilde);

        // Geometry check
        const gsBSplineBasis<real_t> basis_geo = dynamic_cast<const gsBSplineBasis<real_t> &>(patch.basis().component(dir));
        gsKnotVector<real_t> kv_geo = dynamic_cast<const gsKnotVector<real_t> &>(basis_geo.knots());
        index_t p_geo = kv_geo.degree();

        std::vector<real_t> unique_geo = kv_geo.unique();
        for (size_t i = 1; i < unique_geo.size()-1; i++) // First and last are ignored
        {
            real_t knot_geo = unique_geo.at(i);
            index_t r_geo = p_geo - kv_geo.multiplicities().at(i) - 1;
            //index_t p_basis = kv_gluingData.degree();
            //index_t r_basis = p_basis - kv_gluingData.multiplicities().at(kv_gluingData.uFind(knot_geo).uIndex());
            if (r_tilde > r_geo)
                kv_gluingData.insert(knot_geo, r_tilde-r_geo);
            else if (r_tilde < r_geo && r_tilde != 1)
                kv_gluingData.remove(knot_geo, r_geo-r_tilde);
        }
        result = gsBSplineBasis<real_t>(kv_gluingData); // S(\tilde{p},\tilde{r},h)
    }

    void createPlusSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        result = gsBSplineBasis<real_t>(basis_1);

        if (m != 1)
            result.elevateContinuity(1);

        // Geometry check
        const gsBSplineBasis<real_t> basis_geo = dynamic_cast<const gsBSplineBasis<real_t> &>(patch.basis().component(dir));
        gsKnotVector<real_t> kv_geo = dynamic_cast<const gsKnotVector<real_t> &>(basis_geo.knots());
        index_t p_geo = kv_geo.degree();

        gsKnotVector<real_t> kv_basis = dynamic_cast<const gsKnotVector<real_t> &>(result.knots());
        index_t p_basis = kv_basis.degree();
        std::vector<real_t> unique_geo = kv_geo.unique();
        for (size_t i = 1; i < unique_geo.size()-1; i++) // First and last are ignored
        {
            real_t knot_geo = unique_geo.at(i);
            index_t r_geo = p_geo - kv_geo.multiplicities().at(i);
            index_t r_basis = p_basis - kv_basis.multiplicities().at(kv_basis.uFind(knot_geo).uIndex());
            if (r_basis > r_geo)
                kv_basis.insert(knot_geo, r_basis-r_geo);
            else if (r_basis < r_geo)
                kv_basis.remove(knot_geo, r_geo-r_basis);
        }

        result = gsBSplineBasis<real_t>(kv_basis);
    }

    void createMinusSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & result)
    {
        gsBSplineBasis<real_t> basis_1 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(dir));

        index_t m, p;
        p = basis_1.degree();
        m = basis_1.knots().multiplicityIndex(p+1);

        result = gsBSplineBasis<real_t>(basis_1);

        result.degreeDecrease(1);
        if (m != 1)
            result.elevateContinuity(1);


        // Geometry check
        const gsBSplineBasis<real_t> basis_geo = dynamic_cast<const gsBSplineBasis<real_t> &>(patch.basis().component(dir));
        gsKnotVector<real_t> kv_geo = dynamic_cast<const gsKnotVector<real_t> &>(basis_geo.knots());
        index_t p_geo = kv_geo.degree();

        gsKnotVector<real_t> kv_basis = dynamic_cast<const gsKnotVector<real_t> &>(result.knots());
        index_t p_basis = kv_basis.degree();
        std::vector<real_t> unique_geo = kv_geo.unique();
        for (size_t i = 1; i < unique_geo.size()-1; i++) // First and last are ignored
        {
            real_t knot_geo = unique_geo.at(i);
            index_t r_geo = p_geo - kv_geo.multiplicities().at(i) - 1;
            index_t r_basis = p_basis - kv_basis.multiplicities().at(kv_basis.uFind(knot_geo).uIndex());
            if (r_basis > r_geo)
                kv_basis.insert(knot_geo, r_basis-r_geo);
            else if (r_basis < r_geo)
                kv_basis.remove(knot_geo, r_geo-r_basis);
        }

        result = gsBSplineBasis<real_t>(kv_basis);
    }

    void createEdgeSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & basis_plus,
                        gsBSplineBasis<real_t> & basis_minus, gsBSplineBasis<real_t> & basis_gluingData,
                        gsTensorBSplineBasis<2, real_t> & result)
    {
        gsTensorBSplineBasis<2, real_t> basis_edge = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis);

        basis_edge.component(dir).setDegreePreservingMultiplicity(basis_plus.degree()+basis_gluingData.degree()-1);

        index_t r_plus, r_minus, r_edge, r_gD, r;
        r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
        r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
        r_gD = basis_gluingData.degree() - basis_gluingData.knots().multiplicityIndex(basis_gluingData.degree()+1);
        r_edge = basis_edge.degree(dir) - basis_edge.knots(dir).multiplicityIndex(basis_edge.degree(dir)+1);

        r = math::min(r_gD, math::min(r_plus-1, r_minus));
        if (r_edge > r)
            basis_edge.component(dir).reduceContinuity(r_edge - r);
        else if (r_edge < r)
            basis_edge.component(dir).elevateContinuity(r - r_edge);

        result = gsTensorBSplineBasis<2, real_t>(basis_edge);
    }

    void createEdgeSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, index_t dir, gsBSplineBasis<real_t> & basis_plus,
                        gsBSplineBasis<real_t> & basis_minus, gsTensorBSplineBasis<2, real_t> & result)
    {
        gsTensorBSplineBasis<2, real_t> basis_edge = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis);

        gsKnotVector<real_t> basis_edge_knots1 = basis_edge.knots(dir);
        gsKnotVector<real_t> basis_plus_knots = basis_plus.knots();
        gsKnotVector<real_t> basis_minus_knots = basis_minus.knots();
        
        basis_edge.component(dir).setDegreePreservingMultiplicity(basis_plus.degree());

        index_t r_plus, r_minus, r_edge, r;
        r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
        r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
        r_edge = basis_edge.degree(dir) - basis_edge.knots(dir).multiplicityIndex(basis_edge.degree(dir)+1);

        r = math::min(r_plus-1, r_minus);
        if (r_edge > r)
            basis_edge.component(dir).reduceContinuity(r_edge - r);
        else if (r_edge < r)
            basis_edge.component(dir).elevateContinuity(r - r_edge);

        // Geometry check
        const gsBSplineBasis<real_t> basis_geo = dynamic_cast<const gsBSplineBasis<real_t> &>(patch.basis().component(dir));
        gsKnotVector<real_t> kv_geo = dynamic_cast<const gsKnotVector<real_t> &>(basis_geo.knots());
        index_t p_geo = kv_geo.degree();

        gsKnotVector<real_t> & kv_basis = dynamic_cast<gsKnotVector<real_t> &>(basis_edge.knots(dir));
        index_t p_basis = kv_basis.degree();
        std::vector<real_t> unique_geo = kv_geo.unique();
        for (size_t i = 1; i < unique_geo.size()-1; i++) // First and last are ignored
        {
            real_t knot_geo = unique_geo.at(i);
            index_t r_geo = p_geo - kv_geo.multiplicities().at(i) - 1;
            index_t r_basis = p_basis - kv_basis.multiplicities().at(kv_basis.uFind(knot_geo).uIndex());
            if (r_basis > r_geo)
                kv_basis.insert(knot_geo, r_basis-r_geo);
            else if (r_basis < r_geo)
                kv_basis.remove(knot_geo, r_geo-r_basis);
        }
        result = gsTensorBSplineBasis<2, real_t>(basis_edge);
    }

    // Corner Vertex
    void createVertexSpace(const gsGeometry<real_t> & patch, gsBasis<real_t> & basis, bool isInterface_1,
                           bool isInterface_2, gsTensorBSplineBasis<2, real_t> & result, index_t p_tilde, index_t r_tilde)
    {
        gsTensorBSplineBasis<2, real_t> basis_vertex = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(basis);

        for (index_t dir = 0; dir < basis.domainDim(); dir ++)
        {
            if (dir == 0 ? isInterface_1 : isInterface_2) // If edge is an interface
            {
                gsBSplineBasis<real_t> basis_gluingData, basis_plus, basis_minus;
                createGluingDataSpace(patch, basis, dir, basis_gluingData, p_tilde, r_tilde);
                createPlusSpace(patch, basis, dir, basis_plus);
                createMinusSpace(patch, basis, dir, basis_minus);

                basis_vertex.component(dir).setDegreePreservingMultiplicity(basis_plus.degree()+basis_gluingData.degree()-1);

                index_t r_plus, r_minus, r_edge, r_gD, r;
                r_plus = basis_plus.degree() - basis_plus.knots().multiplicityIndex(basis_plus.degree()+1); // p+1, since c++ starts at 0
                r_minus = basis_minus.degree() - basis_minus.knots().multiplicityIndex(basis_minus.degree()+1);
                r_gD = basis_gluingData.degree() - basis_gluingData.knots().multiplicityIndex(basis_gluingData.degree()+1);
                r_edge = basis_vertex.degree(dir) - basis_vertex.knots(dir).multiplicityIndex(basis_vertex.degree(dir)+1);

                r = math::min(r_gD, math::min(r_plus, r_minus));
                if (r_edge > r)
                    basis_vertex.component(dir).reduceContinuity(r_edge - r);
                else if (r_edge < r)
                    basis_vertex.component(dir).elevateContinuity(r - r_edge);
            }
            else
            {
                index_t r_12, p_12;
                p_12 = basis.degree(dir);
                r_12 = p_12 - basis_vertex.component(dir).knots().multiplicityIndex(p_12+1);
                if (r_12 == p_12 - 1) // == basis_vertex_1.degree(1)
                    basis_vertex.component(dir).reduceContinuity(1); // In the case for the max. smoothness
            }
        }

        result = gsTensorBSplineBasis<2, real_t>(basis_vertex);
    }
}



