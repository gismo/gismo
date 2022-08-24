/** @file gsExprEvaluator.h

    @brief Generic expressions evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsIO/gsParaviewCollection.h>
#include<gsAssembler/gsQuadrature.h>

namespace gismo
{

/**
   \brief Generic evaluator of isogeometric expressions

   The expressions may be scalar ot vector-valued. Computed quatities
   can be global or element.wise.
*/
template<class T>
class gsExprEvaluator
{
private:
    typename gsExprHelper<T>::Ptr m_exprdata;

    expr::gsFeElement<T> m_element;

private:
    std::vector<T> m_elWise;
    T              m_value;

    gsOptionList m_options;

public:
    typedef std::vector< boundaryInterface > intContainer;
    typedef std::vector< patchSide > bContainer;

    typedef typename gsExprHelper<T>::element     element;
    typedef typename gsExprHelper<T>::geometryMap geometryMap;
    typedef typename gsExprHelper<T>::variable    variable;

public:

    gsExprEvaluator() : m_exprdata(gsExprHelper<T>::make()),
    m_options(defaultOptions()) { }

    gsExprEvaluator(typename gsExprHelper<T>::Ptr env)
    : m_exprdata(env), m_value(-1), m_options(defaultOptions())
    { }

    //gsExprEvaluator(typename gsExprHelper<T> env)

    gsExprEvaluator(const gsExprAssembler<T> & o)
    : m_exprdata(o.exprData()), m_options(defaultOptions())
    { }

    gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
        opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
        opt.addInt ("plot.npts", "Number of sampling points for plotting", 3000 );
        opt.addSwitch("plot.elements", "Include the element mesh in plot (when applicable)", false);
        //opt.addSwitch("plot.cnet", "Include the control net in plot (when applicable)", false);
        return opt;
    }

    gsOptionList & options() {return m_options;}

public:

    /// Returns the last computed value
    T value() const { return m_value; }

    /// The number of lastly computed values
    size_t nValues() const { return m_elWise.size(); }

    /// Returns an std::vector containing the last computed values per element.
    const std::vector<T> & elementwise() const { return m_elWise; }

    /// Returns a vector containing the last computed values per element.
    gsAsConstVector<T> allValues() const { return gsAsConstVector<T>(m_elWise); }

    /// Returns the last computed values per element, resized as a matrix
    gsAsConstMatrix<T> allValues(index_t nR, index_t nC) const
    { return gsAsConstMatrix<T>(m_elWise, nR, nC); }

    /// \brief Sets the domain of integration.
    /// \warning Must be called before any computation is requested
    void setIntegrationElements(const gsMultiBasis<T> & mesh)
    { m_exprdata->setMultiBasis(mesh); }

    /// Registers \a mp as an isogeometric geometry map and return a handle to it
    geometryMap getMap(const gsMultiPatch<T> & mp) //conv->tmp->error
    { return m_exprdata->getMap(mp); }

    /// Registers \a g as an isogeometric geometry map and return a handle to it
    geometryMap getMap(const gsFunction<T> & gm)
    { return m_exprdata->getMap(gm); }

    /// Registers \a func as a variable and returns a handle to it
    variable getVariable(const gsFunctionSet<T> & func, index_t dim = 1)
    { return m_exprdata->getVar(func, dim); }

    /// Registers \a func as a variable defined on \a G and returns a handle to it
    variable getVariable(const gsFunctionSet<T> & func, geometryMap G)
    { return m_exprdata->getVar(func, G); }

    /// Returns a handle to an isogeometric element
    element getElement() const { return m_element; }

    /// Calculates the square root of the lastly computed quantities (eg. integrals)
    void calcSqrt()
    {
        gsAsVector<T> en(m_elWise);
        en.array() = en.array().sqrt();
        m_value = math::sqrt(m_value);
    }

    /// Calculates the \a p-th root of the lastly computed quantities (eg. integrals)
    void calcRoot(const index_t p)
    {
        gsAsVector<T> en(m_elWise);
        en.array() = en.array().pow(static_cast<T>(1)/p);
        m_value = math::pow(m_value, static_cast<T>(1)/p );
    }

public:

    /// Calculates the integral of the expression \a expr on the whole integration domain
    template<class E>
    T integral(const expr::_expr<E> & expr)
    { return compute_impl<E,false,plus_op>(expr); }

    /// Calculates the integral of the expression \a expr on each element
    template<class E>
    T integralElWise(const expr::_expr<E> & expr)
    { return compute_impl<E,true,plus_op>(expr); }

    // overloads for scalars
    T integral(const T & val) { return integral<T>(val); }
    T integralElWise(const T & val) { return integralElWise<T>(val); }

    /// Calculates the integral of the expression \a expr on the
    /// boundary of the integration domain
    template<class E> // note: integralBdrElWise not offered
    T integralBdr(const expr::_expr<E> & expr)
    { return computeBdr_impl<E,plus_op>(expr,
      m_exprdata->multiBasis().topology().boundaries()); }

    /// Calculates the integral of the expression \a expr on the
    /// boundaries contained in \a bdrlist
    template<class E> // note: integralBdrElWise not offered
    T integralBdr(const expr::_expr<E> & expr, const bContainer & bdrlist)
    { return computeBdr_impl<E,plus_op>(expr,bdrlist); }

    /// Calculates the integral of the expression \a expr on the
    /// interfaces of the (multi-basis) integration domain
    template<class E> // note: elementwise integral not offered
    T integralInterface(const expr::_expr<E> & expr)
    { return computeInterface_impl<E,plus_op>(expr,
      m_exprdata->multiBasis().topology().interfaces()); }

    /// Calculates the integral of the expression \a expr on the
    /// interfaces \a iFaces of the integration domain
    template<class E> // note: elementwise integral not offered
    T integralInterface(const expr::_expr<E> & expr, const intContainer & iFaces)
    { return computeInterface_impl<E,plus_op>(expr, iFaces); }

    /// Calculates the maximum value of the expression \a expr by
    /// sampling over a finite number of points
    template<class E>
    T max(const expr::_expr<E> & expr)
    { return compute_impl<E,false,max_op>(expr); }

    /// Calculates the maximum value of the expression \a expr by
    /// on each element by sampling over a finite number of points
    template<class E>
    T maxElWise(const expr::_expr<E> & expr)
    { return compute_impl<E,true,max_op>(expr); }

    /// Calculates the maximum of the expression \a expr on the
    /// interfaces of the (multi-basis) integration domain
    template<class E> // note: elementwise integral not offered
    T maxInterface(const expr::_expr<E> & expr)
    { return computeInterface_impl<E,max_op>(expr, m_exprdata->multiBasis().topology().interfaces()); }

    /// Calculates the minimum value of the expression \a expr by
    /// sampling over a finite number of points
    template<class E>
    T min(const expr::_expr<E> & expr)
    { return compute_impl<E,false,min_op>(expr); }

    /// Calculates the minimum value of the expression \a expr
    /// on each element by sampling over a finite number of points
    template<class E>
    T minElWise(const expr::_expr<E> & expr)
    { return compute_impl<E,true,min_op>(expr); }

    /// Computes values of the expression \a expr
    /// at the grid points \a git of patch \a patchId
#ifdef __DOXYGEN__
    template<class E> void
#else
    template<class E, int mode, short_t d>
    typename util::enable_if<E::ScalarValued,void>::type
#endif
    eval(const expr::_expr<E> & expr,
         gsGridIterator<T,mode,d> & git,
         const index_t patchInd = 0);

    template<class E, int mode, short_t d>
    typename util::enable_if<!E::ScalarValued,void>::type
    eval(const expr::_expr<E> & expr,
              gsGridIterator<T,mode,d> & git,
              const index_t patchInd = 0);

    /// Computes value of the expression \a testExpr
    /// at the point \a pt of patch \a patchId
    template<class E>
#ifdef __DOXYGEN__
    gsAsConstMatrix<T>
#else
    typename util::enable_if<E::ScalarValued,gsAsConstMatrix<T> >::type
#endif
    eval(const expr::_expr<E> & testExpr, const gsVector<T> & pt,
         const index_t patchInd = 0);

    template<class E>
    typename util::enable_if<!E::ScalarValued,gsAsConstMatrix<T> >::type
    eval(const expr::_expr<E> & testExpr, const gsVector<T> & pt,
         const index_t patchInd = 0);

    /// Computes value of the expression \a expr at the point \a pt of
    /// patch \a patchId, and displays the result
    template<class E> void
    testEval(const expr::_expr<E> & expr,
             const gsVector<T> & pt, const index_t patchInd = 0)
    {
        // expr.printDetail(gsInfo);
        gsInfo << "Result:\n"<< eval(expr,pt,patchInd) <<"\n";
    }

    void info() const { m_exprdata->print(); }

    // Interpolates the expression \a expr over the isogeometric domain \a G
    template<class E> void interpolate(const expr::_expr<E> &)
    {
        GISMO_NO_IMPLEMENTATION
        // for all patches
        //   get anchors of patch
        //   evaluate expr
        //   solve system
        //   add to multipatch
        // return multipatch
    }

    ///\brief Creates a paraview file named \a fn containing valies of the
    //( expression \a expr over the isogeometric domain \a G.
    ///
    /// Plotting properties are controlled by entries in the options
    template<class E>
    void writeParaview(const expr::_expr<E> & expr,
                       geometryMap G, std::string const & fn)
    { writeParaview_impl<E,true>(expr,G,fn); }

    ///\brief Creates a paraview file named \a fn containing valies of the
    //( expression \a expr over the parametric domain.
    ///
    /// Plotting properties are controlled by entires in the options
    template<class E>
    void writeParaview(const expr::_expr<E> & expr,
                       std::string const & fn)
    { writeParaview_impl<E,false>(expr,m_exprdata->getMap(),fn); }


private:

    template<class E, bool gmap>
    void writeParaview_impl(const expr::_expr<E> & expr,
                            geometryMap G, std::string const & fn);

    template<class E, bool storeElWise, class _op>
    T compute_impl(const expr::_expr<E> & expr);

    template<class E, class _op>
    T computeBdr_impl(const expr::_expr<E> & expr, const bContainer & bdrlist);

    template<class E, class _op>
    T computeInterface_impl(const expr::_expr<E> & expr, const intContainer & iFaces);

    template<class E>
    void computeGrid_impl(const expr::_expr<E> & expr, const index_t patchInd);

    struct plus_op
    {
        static inline T init() { return 0; }
        static inline void acc(const T contrib, const T w, T & res) { res += w * contrib; }
    };
    struct min_op
    {
        static inline T init() { return math::limits::max(); }
        static inline void acc (const T contrib, const T, T & res)
        {
            res = math::min(contrib, res);
        }
    };
    struct max_op
    {
        static inline T init() { return math::limits::min(); }
        static inline void acc (const T contrib, const T, T & res)
        {
            res = math::max(contrib, res);
        }
    };

};

template<class T>
template<class E, bool storeElWise, class _op>
T gsExprEvaluator<T>::compute_impl(const expr::_expr<E> & expr)
{
    // GISMO_ASSERT( expr.isScalar(), // precompute
    //               "Expecting scalar expression instead of "
    //               <<expr.cols()<<" x "<<expr.rows() );
    //expr.print(gsInfo); // precompute

    gsQuadRule<T> QuRule;  // Quadrature rule
    gsVector<T> quWeights; // quadrature weights

    // initialize flags
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
    m_exprdata->setFlags(expr, SAME_ELEMENT, SAME_ELEMENT);

    // Computed value
    T elVal;
    m_value = _op::init();
    m_elWise.clear();
    if ( storeElWise )
        m_elWise.reserve(m_exprdata->multiBasis().totalElements());
    for (unsigned patchInd=0; patchInd < m_exprdata->multiBasis().nBases(); ++patchInd)
    {
        // Quadrature rule
        QuRule =  gsQuadrature::get(m_exprdata->multiBasis().basis(patchInd), m_options);
        //gsDebugVar(QuRule.numNodes());

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().piece(patchInd).makeDomainIterator();
        m_element.set(*domIt);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(patchInd);

            // Compute on element
            elVal = _op::init();
            for (index_t k = 0; k != quWeights.rows(); ++k) // loop over quadrature nodes
                _op::acc(expr.val().eval(k), quWeights[k], elVal);

            //gsDebugVar(elVal);
            _op::acc(elVal, 1, m_value);
            if ( storeElWise )
                m_elWise.push_back( elVal );
        }
    }

    return m_value;
}

template<class T>
template<class E, class _op>
T gsExprEvaluator<T>::computeBdr_impl(const expr::_expr<E> & expr,
                                      const bContainer & bdrlist)
{
    // GISMO_ASSERT( expr.isScalar(),
    //               "Expecting scalar expression instead of "
    //               <<expr.cols()<<" x "<<expr.rows() );

    //expr.print(gsInfo);

    gsQuadRule<T> QuRule;  // Quadrature rule
    gsVector<T> quWeights; // quadrature weights

    // initialize flags
    m_exprdata->setFlags(expr, SAME_ELEMENT, SAME_ELEMENT);

    // Computed value
    T elVal;
    m_value = _op::init();
    m_elWise.clear();

    for (typename gsBoxTopology::const_biterator bit = //!! not multipatch!
             bdrlist.begin(); bit != bdrlist.end(); ++bit)
    {
        // Quadrature rule
        QuRule = gsQuadrature::get(m_exprdata->multiBasis().basis(bit->patch), m_options,bit->direction());
        m_exprdata->mapData.side = bit->side();

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().piece(bit->patch).makeDomainIterator(bit->side());
        m_element.set(*domIt);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(bit->patch);

            // Compute on element
            elVal = _op::init();
            for (index_t k = 0; k != quWeights.rows(); ++k) // loop over quadrature nodes
                _op::acc(expr.val().eval(k), quWeights[k], elVal);

            _op::acc(elVal, 1, m_value);
            //if ( storeElWise ) m_elWise.push_back( elVal );
        }
    }

    return m_value;
}

template<class T>
template<class E, class _op>
T gsExprEvaluator<T>::computeInterface_impl(const expr::_expr<E> & expr, const intContainer & iFaces)
{
    // GISMO_ASSERT( expr.isScalar(),
    //               "Expecting scalar expression instead of "
    //               <<expr.cols()<<" x "<<expr.rows() );

    //expr.print(gsInfo);

    gsQuadRule<T> QuRule;  // Quadrature rule
    gsVector<T> quWeights; // quadrature weights

    // initialize flags
    m_exprdata->setFlags(expr, SAME_ELEMENT, SAME_ELEMENT);

    // Computed value
    T elVal;
    m_value = _op::init();
    m_elWise.clear();

    for (typename gsBoxTopology::const_iiterator iit = //!! not multipatch!
             iFaces.begin(); iit != iFaces.end(); ++iit)
    {
        const boundaryInterface & iFace = *iit;
        const index_t patch1 = iFace.first().patch;
        // const index_t patch2 = iFace.second().patch; //!
        // Quadrature rule
        QuRule = gsQuadrature::get(m_exprdata->multiBasis().basis(patch1),
                                   m_options, iFace.first().side().direction());

        m_exprdata->mapData.side = iFace.first().side();

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(patch1).makeDomainIterator(iFace.first().side());
        m_element.set(*domIt);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(patch1);

            // Compute on element
            elVal = _op::init();
            for (index_t k = 0; k != quWeights.rows(); ++k) // loop over quadrature nodes
                _op::acc(expr.val().eval(k), quWeights[k], elVal);

            _op::acc(elVal, 1, m_value);
            //if ( storeElWise ) m_elWise.push_back( elVal );
        }
    }

    return m_value;
}


template<class T>
template<class E, int mode, short_t d>
typename util::enable_if<E::ScalarValued,void>::type
gsExprEvaluator<T>::eval(const expr::_expr<E> & expr,
                         gsGridIterator<T,mode,d> & git,
                         const index_t patchInd)
{ // to remove
    // bug: fails due to gsFeVariable::rows() before evaluation
    // GISMO_ASSERT( expr.isScalar(), "Expecting scalar");

    m_exprdata->initFlags();
    expr.setFlag();
    m_elWise.clear();
    m_elWise.reserve(git.numPoints());

    for( git.reset(); git; ++git )
    {
        m_exprdata->points() = *git;
        m_exprdata->precompute(patchInd);
        m_elWise.push_back( expr.val().eval(0) );

        // equivalent:
        //m_elWise.push_back( m_exprdata->eval(expr).value() );

    }
    m_value = m_elWise.back(); // not used
}

template<class T>
template<class E, int mode, short_t d>
typename util::enable_if<!E::ScalarValued,void>::type
gsExprEvaluator<T>::eval(const expr::_expr<E> & expr,
                         gsGridIterator<T,mode,d> & git,
                         const index_t patchInd)
{
    // bug: fails due to gsFeVariable::rows() before evaluation
    // GISMO_ASSERT( expr.isScalar(), "Expecting scalar");

    m_exprdata->setFlags(expr);
    m_elWise.clear();
    m_elWise.reserve(git.numPoints());
    gsMatrix<T> tmp;
    for( git.reset(); git; ++git )
    {
        m_exprdata->points() = *git; // ..
        m_exprdata->precompute(patchInd);
        // m_exprdata->eval_into(expr, tmp);
        tmp = expr.eval(0);
        m_elWise.insert(m_elWise.end(), tmp.data(), tmp.data()+tmp.size());
    }
    m_value = 0; // not used
}


template<class T>
template<class E>
typename util::enable_if<E::ScalarValued,gsAsConstMatrix<T> >::type
gsExprEvaluator<T>::eval(const expr::_expr<E> & expr, const gsVector<T> & pt,
                         const index_t patchInd)
{
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE);
    m_exprdata->setFlags(expr);
    m_elWise.clear();
    m_exprdata->points() = pt;
    m_exprdata->precompute(patchInd);

    // expr.printDetail(gsInfo); //

    m_value = expr.val().eval(0);
    return gsAsConstMatrix<T>(&m_value,1,1);
}

template<class T>
template<class E>
typename util::enable_if<!E::ScalarValued,gsAsConstMatrix<T> >::type
gsExprEvaluator<T>::eval(const expr::_expr<E> & expr, const gsVector<T> & pt,
                         const index_t patchInd)
{
    m_exprdata->setFlags(expr);
    m_exprdata->points() = pt;
    m_exprdata->precompute(patchInd);

    // expr.printDetail(gsInfo); //after precompute

    gsMatrix<T> tmp = expr.eval(0);
    // const index_t r = expr.rows();
    // const index_t c = expr.cols();
    // gsInfo <<"tmp - "<< tmp.dim() <<" rc - "<< r <<", "<<c<<"\n";
    const index_t r = tmp.rows();
    const index_t c = tmp.cols();
    m_elWise.resize(r*c);
    gsAsMatrix<T>(m_elWise, r, c) = tmp; //expr.eval(0);
    return gsAsConstMatrix<T>(m_elWise, r, c);
}

template<class T>
template<class E, bool gmap>
void gsExprEvaluator<T>::writeParaview_impl(const expr::_expr<E> & expr,
                                            geometryMap G,
                                            std::string const & fn)
    {
        m_exprdata->setFlags(expr);

        //if false, embed topology ?
        const index_t n = m_exprdata->multiBasis().nBases();
        gsParaviewCollection collection(fn);
        std::string fileName;

        gsMatrix<T> pts, vals, ab;

        const bool mesh = m_options.askSwitch("plot.elements");

        for ( index_t i=0; i != n; ++i )
        {
            fileName = fn + util::to_string(i);
            unsigned nPts = m_options.askInt("plot.npts", 3000);
            ab = m_exprdata->multiBasis().piece(i).support();
            gsGridIterator<T,CUBE> pt(ab, nPts);
            eval(expr, pt, i);
            nPts = pt.numPoints();
            vals = allValues(m_elWise.size()/nPts, nPts);

            if (gmap) // Forward the points ?
            {
                eval(G, pt, i);
                pts = allValues(m_elWise.size()/nPts, nPts);
            }

            gsWriteParaviewTPgrid( gmap ? pts : pt.toMatrix(), // parameters
                                  vals,
                                  pt.numPointsCwise(), fileName );
            collection.addPart(fileName, ".vts");

            if ( mesh )
            {
                fileName+= "_mesh";
                gsMesh<T> msh(m_exprdata->multiBasis().basis(i), 2);
                static_cast<const gsGeometry<T>&>(G.source().piece(i)).evaluateMesh(msh);
                gsWriteParaview(msh, fileName, false);
                collection.addPart(fileName, ".vtp");
            }
        }
        collection.save();
    }

} //namespace gismo
