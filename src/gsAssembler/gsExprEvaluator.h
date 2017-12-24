/** @file gsExprEvaluator.h

    @brief Generic expressions evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsAssembler/gsExprHelper.h>
#include<gsAssembler/gsExprAssembler.h>

#include<gsIO/gsParaviewCollection.h>

namespace gismo
{

/**
   Generic evaluator of (scalar) element-wise isogeometric expressions
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

    typedef typename gsExprHelper<T>::element     element;
    typedef typename gsExprHelper<T>::geometryMap geometryMap;
    typedef typename gsExprHelper<T>::variable    variable;

public:

    gsExprEvaluator() : m_exprdata(gsExprHelper<T>::New()),
    m_options(defaultOptions()) { }

    gsExprEvaluator(const gsExprAssembler<T> & o)
    : m_exprdata(o.exprData()), m_options(defaultOptions()) //++
    { }

    gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
        opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
        return opt;
    }

    gsOptionList & options() {return m_options;}
    
public:
    
    T value() const { return m_value; }

    gsAsConstVector<T> allValues() const { return gsAsConstVector<T>(m_elWise); }

    gsAsConstMatrix<T> allValues(index_t nR, index_t nC) const
    { return gsAsConstMatrix<T>(m_elWise, nR, nC); }
    
    void setIntegrationElements(const gsMultiBasis<T> & mesh)
    { m_exprdata->setMultiBasis(mesh); }

    geometryMap setMap(const gsMultiPatch<T> & mp) //conv->tmp->error
    { return m_exprdata->setMap(mp); }

    geometryMap setMap(const gsGeometry<T> & mp)
    { return m_exprdata->setMap(mp); }

    variable setVariable(const gsFunctionSet<T> & func, index_t dim = 1)
    { return m_exprdata->setVar(func, dim); }

    variable setVariable(const gsFunctionSet<T> & func, geometryMap G)
    { return m_exprdata->setVar(func, G); }

    element getElement() const { return m_element; }

    void calcSqrt()
    {
        gsAsVector<T> en(m_elWise);
        en.array() = en.array().sqrt();
        m_value = math::sqrt(m_value);
    }
    
    void calcRoot(const index_t p)
    {
        gsAsVector<T> en(m_elWise);
        en.array() = en.array().pow(static_cast<T>(1)/p);
        m_value = math::pow(m_value, static_cast<T>(1)/p );
    }

public:
    
    template<class E> 
    T integral(const expr::_expr<E> & expr)
    { return compute_impl<E,false,plus_op>(expr); }

    // specialize (overload) for scalar
    T integral(const T & val) { return integral<T>(val); }

    template<class E> 
    T integralElWise(const expr::_expr<E> & expr)
    { return compute_impl<E,true,plus_op>(expr); }

    template<class E> // note: integralBdrElWise not offered
    T integralBdr(const expr::_expr<E> & expr)
    { return computeBdr_impl<E,plus_op>(expr); }

	template<class E> // note: elementwise integral not offered
    T integralInterface(const expr::_expr<E> & expr)
    { return computeInterface_impl<E,plus_op>(expr, m_exprdata.multiBasis().topology().interfaces()); }

	template<class E> // note: elementwise integral not offered
    T integralInterface(const expr::_expr<E> & expr, const intContainer & iFaces)
    { return computeInterface_impl<E,plus_op>(expr, iFaces); }
    
    template<class E> 
    T max(const expr::_expr<E> & expr)
    { return compute_impl<E,false,max_op>(expr); }

    template<class E> 
    T maxElWise(const expr::_expr<E> & expr)
    { return compute_impl<E,true,max_op>(expr); }

    template<class E> 
    T min(const expr::_expr<E> & expr)
    { return compute_impl<E,false,min_op>(expr); }

    template<class E> 
    T minElWise(const expr::_expr<E> & expr)
    { return compute_impl<E,true,min_op>(expr); }

    template<class E, int mode, int d>
    typename util::enable_if<E::ScalarValued,void>::type
    eval(const expr::_expr<E> & expr,
         gsGridIterator<T,mode,d> & git,
         const index_t patchInd = 0);

    template<class E, int mode, int d>
    typename util::enable_if<!E::ScalarValued,void>::type
    eval(const expr::_expr<E> & expr,
              gsGridIterator<T,mode,d> & git,
              const index_t patchInd = 0);

    template<class E>
    typename util::enable_if<E::ScalarValued,gsAsConstMatrix<T> >::type
    eval(const expr::_expr<E> & testExpr, const gsVector<T> & pt,
         const index_t patchInd = 0);

    template<class E>
    typename util::enable_if<!E::ScalarValued,gsAsConstMatrix<T> >::type
    eval(const expr::_expr<E> & testExpr, const gsVector<T> & pt,
         const index_t patchInd = 0);
    
    template<class E> void
    testEval(const expr::_expr<E> & expr,
             const gsVector<T> & pt, const index_t patchInd = 0)
    {
//        expr.printDetail(gsInfo);
        gsInfo << "Result:\n"<< eval(expr,pt,patchInd) <<"\n";
    }
    
    template<class E>
    void interpolate(const expr::_expr<E> & expr)
    {
        // for all patches
        //   get anchors of patch
        //   evaluate expr
        //   solve system
        //   add to multipatch
        // return multipatch
    }

    template<class E>
    void writeParaview(const expr::_expr<E> & expr,
                       geometryMap G,
                       std::string const & fn, 
                       unsigned nPts = 3000, bool mesh = true)
    {
        //embed topology
        const index_t n = m_exprdata->multiBasis().nBases();
        gsParaviewCollection collection(fn);
        std::string fileName;

        gsMatrix<T> pts, vals, ab;
        
        for ( index_t i=0; i != n; ++i )
        {
            fileName = fn + util::to_string(i);

            ab = m_exprdata->multiBasis().piece(i).support();
            gsGridIterator<T,CUBE> pt(ab, nPts);
            eval(expr, pt, i);
            nPts = pt.numPoints();
            vals = allValues(m_elWise.size()/nPts, nPts);

            // Forward the points
            eval(G, pt, i);
            pts = allValues(m_elWise.size()/nPts, nPts); // give ?
            
            gsWriteParaviewTPgrid(pts, //pt.toMatrix(), // parameters
                                  vals,
                                  pt.numPointsCwise(), fileName );
            collection.addPart(fileName, ".vts");

            if ( mesh ) 
            {
                fileName+= "_mesh";
                gsMesh<T> msh(m_exprdata->multiBasis().basis(i), 2);
                static_cast<const gsGeometry<>&>(G.source().piece(i)).evaluateMesh(msh);
                gsWriteParaview(msh, fileName, false);
                collection.addPart(fileName, ".vtp");
            }
        }
        collection.save();
    }

        
private:
    
    template<class E, bool storeElWise, class _op>
    T compute_impl(const expr::_expr<E> & expr);

    template<class E, class _op>
    T computeBdr_impl(const expr::_expr<E> & expr);

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
        static inline void acc (const T contrib, const T w, T & res)
        {
            GISMO_UNUSED(w);
            res = math::min(contrib, res);
        }
    };
    struct max_op
    {
        static inline T init() { return math::limits::min(); }
        static inline void acc (const T contrib, const T w, T & res)
        {
            GISMO_UNUSED(w);
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
        
    for (unsigned patchInd=0; patchInd < m_exprdata->multiBasis().nBases(); ++patchInd)
    {
        // Quadrature rule
        QuRule = gsGaussRule<T>(m_exprdata->multiBasis().basis(patchInd), m_options);

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

            _op::acc(elVal, 1, m_value);
            if ( storeElWise )
                m_elWise.push_back( elVal );
        }
    }
        
    return m_value;
}    

template<class T>
template<class E, class _op>
T gsExprEvaluator<T>::computeBdr_impl(const expr::_expr<E> & expr)
{
    GISMO_ASSERT( expr.isScalar(),
                  "Expecting scalar expression instead of "
                  <<expr.cols()<<" x "<<expr.rows() );
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
             m_exprdata->multiBasis().topology().bBegin(); bit != m_exprdata->multiBasis().topology().bEnd(); ++bit)
    {
        // Quadrature rule
        QuRule = gsGaussRule<T>(m_exprdata->multiBasis().basis(bit->patch), m_options,
                                bit->direction());
        
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
    GISMO_ASSERT( expr.isScalar(),
                  "Expecting scalar expression instead of "
                  <<expr.cols()<<" x "<<expr.rows() );
	
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
		const int patch1 = iFace.first().patch;
		const int patch2 = iFace.second().patch;
        // Quadrature rule
        QuRule = gsGaussRule<T>(m_exprdata->multiBasis().basis(patch1), m_options,
                                iFace.first().side().direction());

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
template<class E, int mode, int d>
typename util::enable_if<E::ScalarValued,void>::type
gsExprEvaluator<T>::eval(const expr::_expr<E> & expr,
                         gsGridIterator<T,mode,d> & git,
                         const index_t patchInd)
{ // to remove
    // bug: fails due to gsFeVariable::rows() before evaluation
    // GISMO_ASSERT( expr.isScalar(), "Expecting scalar"); 

    m_exprdata->setFlags(expr);
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
template<class E, int mode, int d>
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

    expr.printDetail(gsInfo); //
    
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
    
    expr.printDetail(gsInfo); //after precompute

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
        
} //namespace gismo
