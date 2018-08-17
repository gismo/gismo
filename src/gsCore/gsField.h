/** @file gsField.h

    @brief Provides declaration of the Field class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsUtils/gsPointGrid.h>

namespace gismo
{

/**
 * \brief A scalar of vector field defined on a m_parametric geometry.
 *
 * A gsField is, generally speaking, some mathematical function that is defined on a domain of interest
 * (the name "field" is motivated by, e.g., "scalar field" or "vector field").
 *
 * The gsField combines the following:\n
 * - <b>Geometric information</b> on the domain:\n
 *
 * The domain can be represented as one single patch or as a
 * collection of multiple patches (a.k.a. subdomains).\n This
 * information is stored in a member of the gsMultiPatch class.\n
 *
 * - The <b>function</b> defined on the domain:\n
 *
 * For each patch (a.k.a. subdomain), the gsField contains a member
 * of class gsFunction (which represents the "local field", so to
 * say).  On this, the operations of gsFunction can be carried out
 * (e.g., function evaluation or computation of derivatives).\n
 * Remark: The collection of patch-wise gsFunction is stored in the
 * private member gsField::m_fields.
 *
 * Note that the geometry representation of a single patch can be
 * extracted by calling the member function gsField::patch.
 *
 * The "local field" on a single patch can be extracted by calling gsField::function.
 *
 * \ingroup Core
 */
template<class T>
class gsField
{
private:
    /*
     * Begin of: Stuff from former gsNorms.hpp (not in stable anymore)
     */
    static T computeL2Distance(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numEvals)
    {
        GISMO_ASSERT(u.targetDim() == v.targetDim(), "Functions need to have same target dimension");

        const int d = geo.parDim();
        assert(d == geo.geoDim());

        // compute the tensor Gauss rule
        gsMatrix<T> nodes;
        gsVector<T> weights;
        gsMatrix<T> range = geo.basis().support();

        // Number of nodes of the underlying Gauss rule to use
        const index_t nodesPerInterval = 1;
        const int nodesPerElement = math::ipow(nodesPerInterval, d);
        const int numElements = (numEvals + nodesPerElement - 1) / nodesPerElement;
        std::vector<std::vector<T> > intervals;
        uniformIntervals<T>(range.col(0), range.col(1), intervals, numElements);

        // perform the quadrature
        gsGaussRule<T> quRule(gsVector<index_t>::Constant(d, nodesPerInterval));
        const int numPts = quRule.numNodes();

        gsTensorDomainIterator<T> domIt(intervals);
        gsMatrix<T> geo_pts, geo_jac, u_val, v_val;
        T sum = 0.0;

        for (; domIt.good(); domIt.next())
        {
            // Map the Quadrature rule to the element
            quRule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), nodes, weights);

            // only compute the geometry points if either function is not parametrized
            geo_pts = (!isParametrized_u || !isParametrized_v) ?
                      geo.eval(nodes) : gsMatrix<T>();
            geo_jac = geo.jacobian(nodes);

            // evaluate u and v
            u_val = isParametrized_u ? u.eval(nodes) : u.eval(geo_pts);
            v_val = isParametrized_v ? v.eval(nodes) : v.eval(geo_pts);

            for (index_t k = 0; k < numPts; ++k)
            {
                const T funcDet = math::abs(geo_jac.block(0, k * d, d, d).determinant());
                const gsVector<T> diff = u_val.col(k) - v_val.col(k);
                sum += weights[k] * funcDet * diff.dot(diff);
            }
        }

        return math::sqrt(sum);
    }

    static T igaL2DistanceOnElt( const gsGeometry<T> & geo,
                          const gsGeometry<T> & func,
                          const gsFunction<T> & v,
                          const bool & v_isParam,
                          const typename gsBasis<T>::domainIter & domIt,
                          const gsQuadRule<T> & quRule)
    {
        gsMapData<T> mdGeo, mdFunc;
        mdGeo.flags = NEED_VALUE | NEED_MEASURE;
        mdFunc.flags = NEED_VALUE;

        gsVector<T> quWeightsGeo, quWeightsFunc;
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), mdGeo.points, quWeightsGeo);
        quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), mdFunc.points, quWeightsFunc);

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        //geoEval->evaluateAt(quNodes);
        //funcEval->evaluateAt(quNodes);
        geo.computeMap(mdGeo);
        func.computeMap(mdFunc);
        const gsMatrix<T> & func_vals = mdFunc.values[0];

        // Evaluate function v
        gsMatrix<T> v_val = v_isParam ? v.eval(mdGeo.points)
                                      : v.eval(mdGeo.values[0]);

        // perform the quadrature
        T sum(0.0);

        for (index_t k = 0; k < quRule.numNodes(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeightsGeo[k] * math::abs(mdGeo.measure(k));
            const gsVector<T> diff = func_vals.col(k) - v_val.col(k);
            sum += weight * diff.dot(diff);
        }
        return sum;
    }

    static T igaL2Distance(const gsGeometry<T>& patch,
                    const gsGeometry<T>& func,
                    const gsFunction<T>& v,
                    bool v_isParam)
    {
        //typename gsGeometryEvaluator<T>::uPtr geoEval ( patch.evaluator(NEED_VALUE | NEED_MEASURE ));
        // assuming real-valued function
        //typename gsGeometryEvaluator<T>::uPtr funcEval ( func.evaluator(NEED_VALUE) );

        // degree of the underlying Gauss rule to use
        gsVector<int> numNodes(func.parDim());
        for (index_t i = 0; i < numNodes.size(); ++i)
            numNodes[i] = func.basis().degree(i) + 1;
        //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( func.basis() );

        T sum(0);
        typename gsBasis<T>::domainIter domIt = func.basis().makeDomainIterator();
        gsGaussRule<T> quRule(numNodes);
        for (; domIt->good(); domIt->next())
        {
            sum += igaL2DistanceOnElt(patch, func, v, v_isParam, domIt, quRule);
        }
        return math::sqrt(sum);
    }

    // if u is given as a field, the information is split into
    // "function"- and "geometry"-information and passed to igaL2Distance.
    static T igaFieldL2Distance(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam)
    {
        T dist(0);

        for (int i = 0; i < u.nPieces(); ++i)
        {
            // extract the "function"-part of the gsField
            const gsGeometry<T> & func = static_cast<const gsGeometry<T> &>( u.function(i));
            // call igaL2Distance( patch, func, v, v_isParam)
            const T curDist = igaL2Distance(u.patch(i), func, v, v_isParam);
            dist += curDist * curDist;
        }

        return math::sqrt(dist);
    }

    static T igaFieldL2Distance(const gsField<T>& u, const gsFunction<T>& v, const gsMultiBasis<T>& B, bool v_isParam)
    {
        T dist(0);

        for (int i = 0; i < u.nPieces(); ++i)
        {
            // extract the "function"-part of the gsField
            const gsFunction<T> & func = u.function(i);
            // call igaL2Distance( patch, func, v, v_isParam)
            const T curDist = igaL2Distance(u.patch(i), func, v, B[i], v_isParam);
            dist += curDist * curDist;
        }

        return math::sqrt(dist);
    }

    static T computeL2Distance(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numEvals)
    {
        T dist = T();

        for (int i = 0; i < u.nPieces(); ++i)
        {
            T curDist = computeL2Distance(u.patch(i), u.function(i), u.isParametrized(), v, isParametrized_v, numEvals);
            dist += curDist * curDist;
        }

        return math::sqrt(dist);
    }

    static T computeL2Distance(const gsField<T>& u, const gsField<T>& v, int numEvals)
    {
        T dist = T();
        GISMO_ASSERT(u.nPieces() == v.nPieces(), "Fields not matching: " << u.nPieces() << " != " << v.nPieces());

        for (int i = 0; i < u.nPieces(); ++i)
        {
            T curDist = computeL2Distance(u.patch(i), u.function(i), u.isParametrized(),
                                          v.function(i), v.isParametrized(), numEvals);
            dist += curDist * curDist;
        }

        return math::sqrt(dist);
    }
    /*
     * End of: Stuff from former gsNorms.hpp (not in stable anymore)
     */


public:
    typedef memory::shared_ptr< gsField >  Ptr;// todo: remove
    typedef memory::unique_ptr< gsField >  uPtr;// todo: remove

    gsField(): m_patches(NULL) { }

    gsField( const gsFunctionSet<T> & mp, 
             typename gsFunctionSet<T>::Ptr fs, 
             const bool isparam)
    : m_patches(&mp), m_fields(fs), m_parametric(isparam)
    { }

    gsField( const gsGeometry<T> & sp, const gsFunctionSet<T> & pf, const bool isparam = false) 
    : m_patches(&sp), m_fields(memory::make_shared_not_owned(&pf)), m_parametric(isparam)
    { }

    gsField( const gsGeometry<T> & sp, const gsGeometry<T> & pf) 
    : m_patches(&sp), m_fields(memory::make_shared_not_owned(&pf)), m_parametric(true)
    { }

    gsField( const gsMultiPatch<T> & mp, const gsFunctionSet<T> & f, const bool isparam = false) 
    : m_patches(&mp), m_fields(memory::make_shared_not_owned(&f)), m_parametric(isparam)
    { }

    gsField( const gsMultiPatch<T> & mp, const gsMultiPatch<T> & f) 
    : m_patches(&mp), m_fields(memory::make_shared_not_owned(&f)), m_parametric(true)
    { }

public:
    
// TO DO:

// EVAL_physical_position
// Need to solve x = m_geometry(u) for u

    // Return a point in  the physical domain at parameter value u
    /**
     * @brief Maps points \a u from the parameter domain to the physical domain.
     *
     * @param[in] u Evaluation points as gsMatrix of size <em>d</em> x <em>n</em>.\n
     * \a d denotes the dimension of the parameter domain (i.e., d = parDim()).\n
     * \a n denotes the number of evaluation points.\n
     * Each column of \a u corresponds to one evaluation point.
     * @param[in] i Index of the considered patch/subdomain.
     * @returns uPtr The <em>j</em>-th column of \a uPtr corresponds
     * to the image of the point \a u_j (which is defined by the \a j-th column of the input parameter \a u).
     */
    gsMatrix<T> point(const gsMatrix<T>& u, int i = 0) const
    {
        return m_patches->piece(i).eval(u);
    }

    // Return the value of the Field at parameter value u
    // TO DO: rename to evalParam()
    /**
     * @brief Evaluation of the field at points \a u.
     *
     * @param[in] u Evaluation points as gsMatrix of size <em>d</em> x <em>n</em>.\n
     * \a d denotes the dimension of the parameter domain (i.e., d = parDim()).\n
     * \a n denotes the number of evaluation points.\n
     * Each column of \a u corresponds to one evaluation point.
     * @param[in] i Index of the considered patch/subdomain.
     * @returns uPtr The <em>j</em>-th column of \a uPtr corresponds
     * to the value of the field at the point \a u_j (which is defined by the \a j-th column of the input parameter \a u).
     */
    gsMatrix<T> value(const gsMatrix<T>& u, int i = 0)  const
    {
        return m_parametric
            ? m_fields->piece(i).eval(u)
            : m_fields->piece(i).eval( point(u, i) );
    }

    // Return the value of the Field at physical value u 
    // TO DO: rename to evalPhys()
    gsMatrix<T> pvalue(const gsMatrix<T>& u, int i)  const
    {
        GISMO_ASSERT(!m_parametric, "Cannot compute physical value");
        return ( m_fields->piece(i).eval(u) ); 
    }

    /// Computes the L2-distance between the two fields, on the physical domain
    T distanceL2(gsField<T> const & field, int numEvals= 1000) const 
    {
        return computeL2Distance(*this, field, numEvals);
    }

    /// Computes the L2-distance between the field and a function \a func on the physical domain
    T distanceL2(gsFunction<T> const & func, 
                 bool isFunc_param = false,
                 int numEvals=1000) const
    {
        if (m_parametric) // isogeometric field
            return igaFieldL2Distance(*this, func, isFunc_param);
        else
            return computeL2Distance(*this, func, isFunc_param, numEvals);
    }

    /// Computes the L2-distance between the field and a function \a
    /// func on the physical domain, using mesh from B
    T distanceL2(gsFunction<T> const & func,
                 gsMultiBasis<T> const & B,
                 bool isFunc_param = false,
                 int numEvals=1000) const
    {
        if (m_parametric) // isogeometric field
            return igaFieldL2Distance(*this, func, B, isFunc_param);
        else
            return computeL2Distance(*this, func, isFunc_param, numEvals);
    }

    /// Computes the H1-distance between the field and a function \a
    /// func on the physical domain
    T distanceH1(gsFunction<T> const & func, 
                 bool isFunc_param = false,
                 int = 1000) const
    {
        if ( m_parametric ) // isogeometric field
            return igaFieldH1Distance(*this, func, isFunc_param);
        else
        {
            gsWarn <<"H1 seminorm not implemented.\n";
            return -1;
        }
    }

    /// Computes the H1-distance between the field and a function \a
    /// func on the physical domain, using mesh from B
    T distanceH1(gsFunction<T> const & func,
                 gsMultiBasis<T> const & B,
                 bool isFunc_param = false,
                 int = 1000) const
    {
        if ( m_parametric ) // isogeometric field
            return igaFieldH1Distance(*this, func, B,isFunc_param);
        else
        {
            gsWarn <<"H1 seminorm not implemented.\n";
            return -1;
        }
    }

    /// Computes the DG-distance between the field and a function \a
    /// func on the physical domain
    T distanceDG(gsFunction<T> const & func, 
                 bool isFunc_param = false,
                 int = 1000) const
    {
        if ( m_parametric ) // isogeometric field
            return igaFieldDGDistance(*this, func, isFunc_param);
        else
        {
            gsWarn <<"DG norm not implemented.\n";
            return -1;
        }
    }
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { 
        os << ( m_parametric ? "Parameterized f" : "F") 
           << "unction field.\n Defined on " << m_patches;
        return os; 
    }
    
    /// \brief Returns the dimension of the parameter domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 2).
    int parDim() const { return m_patches->domainDim(); }

    /// \brief Returns the dimension of the physical domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 3).
    int geoDim() const { return m_patches->targetDim(); }

    /// \brief Returns the dimension of the physical domain
    /// (e.g., if the domain is a surface in three-dimensional space, it returns 3).
    int dim() const { return m_fields->targetDim(); }

    /// Returns the number of patches.
    GISMO_DEPRECATED int nPatches()  const { return m_patches->nPieces(); }

    /// Returns the number of pieces.
    int nPieces()  const { return m_patches->nPieces(); }

    const gsGeometry<T> & geometry() const 
    {
        GISMO_ASSERT(dynamic_cast<const gsGeometry<T>*>(m_patches),
                     "No geometry in field. The domain is"<< *m_patches);
        return *static_cast<const gsGeometry<T>*>(m_patches);
    }

    /// Returns gsMultiPatch containing the geometric information on the domain.
    const gsMultiPatch<T> & patches() const    
    { 
        GISMO_ASSERT(dynamic_cast<const gsMultiPatch<T>*>(m_patches),
                     "No patches in field. The field domain is "<< *m_patches);
        return *static_cast<const gsMultiPatch<T>*>(m_patches);
    }

    /// Returns the gsGeometry of patch \a i.
    const gsGeometry<T> & patch(int i=0) const 
    {
        GISMO_ASSERT( i<nPieces(),
                      "gsField: Invalid patch index.");
        GISMO_ASSERT(dynamic_cast<const gsGeometry<T>*>(&m_patches->piece(i)),
                     "No geometry in field. The domain is"<< m_patches->piece(i));
        return static_cast<const gsGeometry<T>&>(m_patches->piece(i));
    }

    /// Returns the gsFunction of patch \a i.
    const gsFunction<T> & function(int i=0) const  
    { 
        GISMO_ASSERT(dynamic_cast<const gsFunction<T>*>(&m_patches->piece(i)),
                     "No function in field. The domain is"<< m_patches->piece(i));
        return static_cast<const gsFunction<T>&>(m_fields->piece(i)); 
    }

    /// Attempts to return an Isogeometric function for patch i
    const gsGeometry<T> & igaFunction(int i=0) const
    {
        GISMO_ASSERT(i<m_fields->nPieces(),
                     "gsField: Invalid patch index.");
        GISMO_ASSERT(m_parametric,
                     "Cannot get an IGA function from non-parametric field.");
        GISMO_ASSERT(dynamic_cast<const gsGeometry<T>*>(&m_fields->piece(i)),
                     "Cannot return an igaFunction from a function of type: "<< m_fields->piece(i) );
        return static_cast<const gsGeometry<T> &>(m_fields->piece(i));
    }

    /// True if the field function is defined on parametric
    /// coordinates (i.e. the same coordinate system as the patches)
    bool isParametric() const { return m_parametric; }
    
    /// True if the field function is parametrized by a set of basis
    /// functions in parametric coordinates (i.e. the same coordinate
    /// system as the patches)
    bool isParametrized() const 
    { return m_parametric && dynamic_cast<const gsGeometry<T>*>(&m_fields->piece(0));}

    /** \brief Returns the coefficient vector (if it exists)
        corresponding to the function field for patch \a i.
    
    Returns the coefficients of the field corresponding to the \a i-th
    patch. This is only possible in the case when the field is defined
    in terms of basis functions (ie. it derives from gsGeometry).
    
    */
    const gsMatrix<T> & coefficientVector(int i=0) const
    {
        return igaFunction(i).coefs();
    }

// Data members
private:

    /// The isogeometric field is defined on this multipatch domain
    const gsFunctionSet<T> * m_patches;

    // If there are many patches, one field per patch

    /// \brief Vector containing "local fields" for each patch/subdomain.
    ///
    /// For each patch/subdomain, the "local field" is represented by
    /// a gsFunction. This local field can be accessed with
    /// gsField::function.
    typename gsFunctionSet<T>::Ptr m_fields;

    /**
     * @brief \a True iff this is an isogeometric field.
     *
     * If \a m_parametric is \a true, the evaluation points for calling gsField::value have to be placed in the
     * \a parameter domain.
     *
     * If \a m_parametric is \a false, then the evaluation points are in the \a physical domain.
     * This applies to, e.g., given exact solutions which are defined on the physical domain.
     */
    bool m_parametric;// True iff this is an Isogeometric field, living on parameter domain

}; // class gsField


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsField<T>& b)
{return b.print(os); }


} // namespace gismo
