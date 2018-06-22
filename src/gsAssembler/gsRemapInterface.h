/** @file gsRemapInterface.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Seiler, R. Schneckenleitner
    Created on: 2018-06-12
*/


#pragma once

#include <gsIO/gsOptionList.h>
#include <gsAssembler/gsAssembler.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsAffineFunction.h>
#include <gsAssembler/gsQuadRule.h>



namespace gismo {


template <class T>
class gsRemapInterface : public gsFunction<T>
{
public:
    /// Shared pointer for gsAffineFunction
    typedef memory::shared_ptr< gsRemapInterface > Ptr;

    /// Unique pointer for gsAffineFunction
    typedef memory::unique_ptr< gsRemapInterface > uPtr;


    gsRemapInterface() {};
    gsRemapInterface(const gsMultiPatch<T> & mp): m_domain(mp)
    {
        GISMO_ASSERT(m_domain.nPatches() == 2, "Domain must consist of two patches!");
        m_interface = m_domain.interfaces()[0]; // Assume that the 2 patches match
    }

    ~gsRemapInterface() {};

    // Member for constructing the reparametrization
    // Fills in m_reparamInterfaceMap
    void constructReparam();

    // Memeber function for constructing the breakpoints
    // Fills in m_breakpoints
    void constructBreaks();

    // Memeber function for evaluating the reparameterization at the quadrature nodes
    void mapQuadNodes(); // eval_into() ???

    void breakpoint_iterator();

    // Helper to compute the closest point to lti on the other patch via Newton's method
    gsMatrix<T> closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & R, const gsMatrix<T> & lti);

    gsMatrix<T> giveBreakpoints() { return m_breakpoints; }
    gsBSpline<T> giveInterfaceMap() { return m_reparamInterfaceMap; }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual int domainDim() const { return m_domain.geoDim(); }

private:
    gsMultiPatch<T> m_domain;
    gsMatrix<T> m_breakpoints;
    gsBSpline<T> m_reparamInterfaceMap;
    //gsAffineFunction<T> m_affineTrans;
    boundaryInterface m_interface;

    gsMatrix<T> m_quadNodes1, m_quadNodes2;

    // Member to enrich a matrix of 1D points to a matrix of m_domain.geoDim() points
    void enrichToVector(const int boundarySide, const gsMatrix<T> & intervals, gsMatrix<T> & pts);

}; // End gsRemapInterface






template<class T>
void gsRemapInterface<T>::constructBreaks() {
    // computes break points per element
    const gsAffineFunction <T> interfaceMap(m_domain.getMapForInterface(m_interface));

    const int patch1 = m_interface.first().patch; // Numeration correct when using a domain with more than two patches?
    const int patch2 = m_interface.second().patch;
    const gsBasis <T> &B1 = m_domain.basis(patch1);// (!) unknown 0
    const gsBasis <T> &B2 = m_domain.basis(patch2);

/*
    const int bSize1 = B1.numElements(m_interface.first().side()); // for hierarchical splines
    const int bSize2 = B2.numElements(m_interface.second().side());
    const int ratio = bSize1 / bSize2;
    GISMO_ASSERT(bSize1 >= bSize2 && bSize1 % bSize2 == 0,
                 "DG assumes nested interfaces. Got bSize1=" <<
                                                             bSize1 << ", bSize2=" << bSize2 << ".");
*/

    // Initialize geometry evaluators
    //typename gsGeometry<T>::Evaluator geoEval1(m_domain[patch1].evaluator(evFlags));
    //typename gsGeometry<T>::Evaluator geoEval2(m_domain[patch2].evaluator(evFlags));

    // Initialize domain element iterators
    // iterates along knot spans in parametric domain
    // this alone is not sufficient to find shifted knot spans
    typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( m_interface.first() .side() );
    typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( m_interface.second().side() );



    // Compute knots in physical domain by evaluating left and right geometry maps at the knot values --------------------------------------
    int numelP1 = domIt1->numElements();
    int numelP2 = domIt2->numElements();
    gsMatrix <T> physicalKnotsP1(m_domain.geoDim(), numelP1 + 1), physicalKnotsP2(m_domain.geoDim(), numelP2 + 1), dummy;

    domIt1->reset();
    domIt2->reset();

    for (int i = 0; i <= numelP1; i++)
    {
        m_domain.patch(patch1).eval_into(domIt1->lowerCorner().row(1), dummy);
        physicalKnotsP1.col(i) = dummy;
        domIt1->next();
    }

    for (int i = 0; i <= numelP2; i++)
    {
        m_domain.patch(patch2).eval_into(domIt2->lowerCorner().row(1), dummy);
        physicalKnotsP2.col(i) = dummy;
        domIt2->next();
    }

    m_domain.patch(patch1).eval_into(domIt1->upperCorner().row(1), dummy);
    physicalKnotsP1.block(0, numelP1, 2, 1) = dummy;

    m_domain.patch(patch2).eval_into(domIt2->upperCorner().row(1), dummy);
    physicalKnotsP2.block(0, numelP2, 2, 1) = dummy;

    gsMatrix<T> physicalBreaks(m_domain.geoDim(), numelP1+numelP2+2);

    physicalBreaks << physicalKnotsP1, physicalKnotsP2;



    //  merge the two physical knot matrices
    gsSortedVector<T> parameterBreaks; //TODO: provide sorting

    domIt1->reset();
    //for (; domIt1->good(); domIt1->next()) {
    //    count++;

        // Determine "corners", i.e. boundaries for intergration for the current element -----
        // this merges knot span boundaries in the physical domain,
        // boundaries in the parametric domain are computed later

        //std::cout << "Element Nummer: " << currentElement << std::endl;
    //    std::vector<gsMatrix<T> > corners;
    //    corners.push_back(Rvalues(1, currentElement - 1));

        // Sorting in physical space probably not needed
        /*for (int i = 0; i < numelP2; i++) {
            if (Lvalues(1, i) > Rvalues(1, currentElement - 1) && Lvalues(1, i) < Rvalues(1, currentElement))
                corners.push_back(Lvalues(1, i));
        }
        */


    //}

    //Determine fixed coordinate of patch1
    // fixedDir ==  0 corresponds to fixed u and 1 corresponds to a fixed v
    index_t fixedDir = m_interface.first().side().direction();

    m_breakpoints = gsMatrix<T> (m_domain.geoDim(), physicalBreaks.cols() - 2 );
    gsMatrix<T> G2_parametric_LC;

    for (int i = 0;  i < physicalBreaks.cols() - 2; i++)
    {
        // computes knot span boundaries in the parametric domain, i.e. the preimages of
        // the boundaries in the physical domain
        m_domain[patch1].invertPoints(physicalBreaks.col(i),G2_parametric_LC);
        if(fixedDir)
            parameterBreaks.push_sorted_unique(G2_parametric_LC(0,0)); // sort w.r.t. u direction
        else
            parameterBreaks.push_sorted_unique(G2_parametric_LC(0,1)); // sort w.r.t. v direction
    }

    for (int i = 0;  i < parameterBreaks.size(); i++)
    {
        if(fixedDir)
            m_breakpoints.col(i) << parameterBreaks[i], G2_parametric_LC(0,1);
        else
            m_breakpoints.col(i) << G2_parametric_LC(0,0), parameterBreaks[i];
    }

}



template<class T>
void gsRemapInterface<T>::breakpoint_iterator()
{
    const int patch1 = m_interface.first().patch; // Numeration correct when using a domain with more than two patches?
    const int patch2 = m_interface.second().patch;
    index_t fixedDir = m_interface.first().side().direction();

    // --- MAIN PART ---------------------------------------------------------------

    for (int i = 0;  i < m_breakpoints.cols() - 1; i++)
    {

        // computes knot span boundaries in the parametric domain, i.e. the preimages of
        // the boundaries in the physical domain
        gsMatrix<T> G2_parametric_UC, G2_parametric_LC, lowerCorner(m_domain.geoDim(), 1), upperCorner(m_domain.geoDim(), 1),
                    quNodes1, quNodes2;
        gsVector<T> quWeights;

        //m_domain[patch1].invertPoints(m_breakpoints.col(i+1),G2_parametric_UC);
        //m_domain[patch1].invertPoints(m_breakpoints.col(i),G2_parametric_LC);

        upperCorner.setZero();
        lowerCorner.setZero();

        upperCorner = m_breakpoints.col(i+1);
        lowerCorner = m_breakpoints.col(i);

        //gsQuadRule<T> QuRule = gsQuadrature::get(basis1, options, side1.direction()); // Quadrature rule

        //QuadRule.mapTo( lowerCorner, upperCorner, quNodes1, quWeights);
        //eval_into(quNodes1,quNodes2); // compute the quadrature nodes on the other box

        //gsMatrix<T> transformed_quadNodes(m_domain.geoDim(),quNodes1.cols());
        // TODO: transformed_quadNodes needs to map from R to R 2
        // we need to reparameterize the nonzero / non-one values of quNodes 2, and pass them to the right
        // row of quNodes2
        //m_reparamInterfaceMap.eval_into(quNodes2.row(1), transformed_quadNodes);


        // Perform required evaluations on the quadrature nodes
        /*
        if ( count % ratio == 0 ) // next master element ?, for hierarchical splines
            domIt2->next();
         */

    }

    //++currentElement;
}

template<class T>
void gsRemapInterface<T>::constructReparam()
{
    const int patch1 = m_interface.first().patch; // Numeration correct when using a domain with more than two patches?
    const int patch2 = m_interface.second().patch;
    const patchSide side = m_interface.first();
    const int numIntervals = 11;
    //const gsAffineFunction<T> affineTans(m_domain.getMapForInterface(m_interface));

    index_t fixedDir = m_interface.first().side().direction();

    gsGeometry<> & GEO_L_ref =  m_domain.patch(m_interface.first().patch);
    gsGeometry<> & GEO_R_ref =  m_domain.patch(m_interface.second().patch);
    gsGeometry<>::uPtr GEO_L_ptr = GEO_L_ref.boundary(m_interface.first());
    gsGeometry<>::uPtr GEO_R_ptr = GEO_R_ref.boundary(m_interface.second());

    // compute interface values
    gsMultiBasis<T> mb(m_domain);

    std::vector<boxCorner> corners;

    // Assume tensor structure
    gsTensorBSplineBasis<2, T>* tb = dynamic_cast<gsTensorBSplineBasis<2, T>* >(&mb.basis(patch1));


    T firstKnot, lastKnot;
    if(fixedDir) // v is fixed
    {
        firstKnot = tb->knots(0).first();
        lastKnot = tb->knots(0).last();
    }
    else // u is fixed
    {
        //gsInfo << "I am here\n";
        firstKnot = tb->knots(1).first();
        lastKnot = tb->knots(1).last();
    }

    gsVector<T> upper(1);
    upper << firstKnot;
    gsVector<T> lower(1);
    lower << lastKnot;

    gsMatrix<T> t_vals = uniformPointGrid(lower, upper, numIntervals); // uniformly distributed samples between first and last knot on the parameter domain

    gsMatrix<T> samples_left, samples_right;
    gsMatrix<T> find_start_value;

    // Get the corresponding edges
    //Edge 1, {(u,v) : u = 0}
    //Edge 2, {(u,v) : u = 1}
    //Edge 3, {(u,v) : v = 0}
    //Edge 4, {(u,v) : v = 1}
    int firstParameterEdge = m_interface.first().index(); // = 1
    int secondParameterEdge = m_interface.second().index(); // = 2

    gsInfo << "left boundary: " << m_interface.first().index() << "\n";
    gsInfo << "right boundary: " << m_interface.second().index() << "\n";

    //gsMatrix<T> vals2dPatch1(t_vals.rows()+1, t_vals.cols()), vals2dPatch2(t_vals.rows()+1, t_vals.cols());
    gsMatrix<T> vals2dPatch1, vals2dPatch2;
    enrichToVector(firstParameterEdge, t_vals, vals2dPatch1);
    enrichToVector(secondParameterEdge, t_vals, vals2dPatch2);


    GEO_L_ref.eval_into(vals2dPatch1,samples_left);
    GEO_R_ref.eval_into(vals2dPatch2,samples_right);

    //gsInfo << "vals2dPatch1:\n" << vals2dPatch1 << "\n vals2dPatch2:\n" << vals2dPatch2 << std::endl;
    //std::cout << "samples left:\n" << samples_left << "\n samples right:\n" << samples_right << std::endl;

    gsMatrix<> B(numIntervals, m_domain.geoDim());

    for (int i=0; i<t_vals.cols();i++)
    {
        // find a suitable start value for the Newton iteration
        find_start_value = (samples_right.colwise()) - samples_left.col(i);
        //find_start_value.noalias() = find_start_value.cwiseAbs();

        size_t row;
        size_t col;

        size_t *row_ptr = &row;
        size_t *col_ptr = &col;

        //gsMatrix<T> smallest = find_start_value.colwise().sum();
        //smallest.minCoeff(row_ptr,col_ptr);
        find_start_value.colwise().squaredNorm().minCoeff(row_ptr,col_ptr);

        //gsMatrix<T> b_null = t_vals(*row_ptr, *col_ptr);
        gsMatrix<T> b_null = samples_right.col(*col_ptr);
        // pass on L if mapping b of "L(b)" is to be found
        // pass on R, if mapping a of "R(a)" is to be found
        gsMatrix<T> b = closestPoint(b_null, GEO_R_ref, samples_left.col(i));
        B.row(i) = b.transpose();


    }

    // the coefficients to fit
    std::cout << "B:\n" << B << std::endl;

    // if one wants to find the coefficients to represent the identity
    // check the error
    // assume that the left map is the identity
    // make sure that "/zuAgnes5Matching.xml" is loaded for checking the error!!
    gsMatrix<T> eval_orig, eval_fit, B2, id;
    //gsFunctionExpr<> p("x", "y", 2);

    gsKnotVector<> KV(0,1,0,4);
    gsCurveFitting<> fit(t_vals.transpose(),B,KV);

    fit.compute();
    const gsBSpline<> & fitted_curve = fit.curve();
    std::cout << "Hi, I'm the resulting curve: \n" << fitted_curve << std::endl;

    gsMatrix<> eval_points = uniformPointGrid(lower, upper, 1000); // to check the error

    //p.eval_into(eval_points, eval_orig);
    fitted_curve.eval_into(eval_points, eval_fit);

    enrichToVector(firstParameterEdge, eval_points, id);

    GEO_L_ref.eval_into(id,B2);

    double error = 0;

    for (int i = 0; i < eval_points.cols(); i++)
        error += (B2.col(i) - eval_fit.col(i)).squaredNorm();

    error = std::sqrt(error);

    std::cout << "Error: " << error << std::endl;

}

template<class T>
gsMatrix<T> gsRemapInterface<T>::closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & intMap, const gsMatrix<T> & lti)
{
    double stopping_criterion = std::numeric_limits<T>::max();
    double eps = 0.001;
    gsMatrix<> b_neu = b_null;

    gsMatrix<> evalGeomap, diffPatch1Patch2, jacPatch2, deriv2Geomap;

    // iterate until stopping criterion is satisfied
    while (stopping_criterion > eps)
    {
        gsMatrix<> b_alt = b_neu;

        // compute value of function to minimize
        intMap.eval_into(b_alt,evalGeomap);
        intMap.jacobian_into(b_alt,jacPatch2);
        intMap.deriv2_into(b_alt, deriv2Geomap);

        diffPatch1Patch2 = (lti - evalGeomap);

        jacPatch2.transposeInPlace();

        //gsInfo << "derivative: " << jacPatch2 << "\n";

        // generate the objective g
        gsMatrix<T> g(m_domain.geoDim(), 1);
        g = jacPatch2 * diffPatch1Patch2;

        // compute the jacobian of g
        // be carefull, is very error-prone
        gsMatrix<T> g_deriv = gsMatrix<T>::Zero(m_domain.geoDim(), m_domain.geoDim());
        g_deriv(0,0) = deriv2Geomap(0,0) * diffPatch1Patch2(0,0) - jacPatch2(0,0) * jacPatch2(0,0) + deriv2Geomap(3,0) * diffPatch1Patch2(1,0) - jacPatch2(0,1)*jacPatch2(0,1);
        g_deriv(0,1) = deriv2Geomap(2,0) * diffPatch1Patch2(0,0) - jacPatch2(0,0) * jacPatch2(1,0) + deriv2Geomap(5,0) * diffPatch1Patch2(1,0) - jacPatch2(0,1)*jacPatch2(1,1);
        g_deriv(1,0) = deriv2Geomap(2,0) * diffPatch1Patch2(0,0) - jacPatch2(0,1) * jacPatch2(0,0) + deriv2Geomap(5,0) * diffPatch1Patch2(1,0) - jacPatch2(1,1)*jacPatch2(0,1);
        g_deriv(1,1) = deriv2Geomap(1,0) * diffPatch1Patch2(0,0) - jacPatch2(0,1) * jacPatch2(0,1) + deriv2Geomap(4,0) * diffPatch1Patch2(1,0) - jacPatch2(1,1)*jacPatch2(1,1);

        // compute the next Newton step
        b_neu = b_alt + g_deriv.fullPivHouseholderQr().solve(g);

        stopping_criterion = (b_alt - b_neu).squaredNorm();
    }

    return b_neu;

}


    template<class T>
void gsRemapInterface<T>::mapQuadNodes()
{
    gsQuadRule<T> QuRule ; // Quadrature rule
    gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights
    // Evaluation flags for the Geometry map
    unsigned evFlags(0);

}

template <typename T>
void gsRemapInterface<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    const int patch1 = m_interface.first().patch; // Numeration correct when using a domain with more than two patches?
    const int patch2 = m_interface.second().patch;
    const gsBasis <T> &B1 = m_domain.basis(patch1);
    const gsBasis <T> &B2 = m_domain.basis(patch2);

    // TODO: Maybe we can do better using the reparametrization, inversion is very expensive in general
    gsMatrix<T> physNodes1;
    // Get the nodes on the physical patch 1
    B1.evalFunc_into(u, m_domain[patch1].coefs(), physNodes1);

    // Invert the nodes on the physical domain for patch 2
    m_domain[patch2].invertPoints(physNodes1, result);
}

// Function to enhance a sequence of 1D points in an interval to 2D points in the parameter domain
// A matrix pts with the already correct dimensions is expected
// boundarySide is the boundary side according to gismo's boundary::side
// intervals is the number of samples to enrich to more dimensions
// pts is the matrix populated with the corresponding points for more dimensions
// only works for 2d at the moment!!
// TODO: Generalize for arbitrary dimensions
template<class T>
void gsRemapInterface<T>::enrichToVector(const int boundarySide, const gsMatrix<T> & intervals, gsMatrix <T> &pts)
{
    pts.resize(m_domain.geoDim(), intervals.cols());

    switch (boundarySide)
    {
        case 1 :
            // u = 0
            pts.row(0) = gsMatrix<T>::Zero(1, intervals.cols());
            pts.row(1) = intervals;

            break;
        case 2 :
            //u = 1;
            pts.row(0) = gsMatrix<T>::Ones(1, intervals.cols());
            pts.row(1) = intervals;
            break;
        case 3 :
            //v = 0;
            pts.row(0) = intervals;
            pts.row(1) = gsMatrix<T>::Zero(1, intervals.cols());
            break;
        case 4 :
            //v = 1;
            pts.row(0) = intervals;
            pts.row(1) = gsMatrix<T>::Ones(1, intervals.cols());
            break;
    }
}

} // End namespace::gismo
