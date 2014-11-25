/*
* gsBasisEvaluator.h created on 24.07.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/**
    \brief The ValueTransformationType enum
**/

enum ValueTransformationType
{
    NO_TRANSFORMATION,
    INVERSE_COMPOSITION,
    DIV_CONFORMING,
    CURL_CONFORMING // not implemented
};

template <typename T> class gsBasisEvaluator;



template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator ( const gsBasis<T> &basis, unsigned flags=0, const gsGeometry<T> *geo=NULL, ValueTransformationType geoTrans=INVERSE_COMPOSITION );

template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator ( const std::vector<gsBasis<T> *> &basis, unsigned flags=0, const gsGeometry<T> *geo=NULL, ValueTransformationType geoTrans=INVERSE_COMPOSITION );



template <typename T>
class gsBasisEvaluator
{
public:
    gsBasisEvaluator(unsigned flags)
        : m_flags(flags)
    {}
    virtual ~gsBasisEvaluator() {}

            void       addFlags    (unsigned newFlags) {this->setFlags(m_flags|newFlags);}
    virtual void       setFlags    (unsigned newFlags) {m_flags = newFlags;}
            unsigned   getFlags    () const            {return m_flags;}
            unsigned   getGeoFlags () const            {return m_geo_flags;}


    virtual       void          evaluateAt ( const gsMatrix<T> &points)  = 0;
    virtual       void          evaluateAt ( const gsMatrix<T> &points, const gsGeometryEvaluator<T> &geoEval)  = 0;

                   unsigned    getParDim   () const {return m_parDim;}
                   unsigned    getTarDim   () const {return m_tarDim;}
                   unsigned    getSpaceDim () const {return m_spaceDim;}

    const          gsMatrix<unsigned>       & actives     ()          const {return m_actives;}

    const          gsMatrix<T>              & values      ()          const {return m_values;}
    const typename gsMatrix<T>::constColumn   value       (index_t k) const {return m_values.col(k);}


    const          gsMatrix<T>              & derivs      ()          const {return m_derivs;}
    const typename gsMatrix<T>::constColumn   deriv       (index_t k) const {return m_derivs.col(k);}

    const          gsMatrix<T>              & jacobians   ()          const {return m_jacobians;}
    const typename gsMatrix<T>::constColumns  jacobian    (index_t k) const
    {
        GISMO_ASSERT(this->m_flags & NEED_JACOBIAN, "Jacobians not computed");
        return m_jacobians.middleCols( k*m_parDim , m_parDim);
    }
    const          gsMatrix<T>              & derivs2     ()          const {return m_2ndDers;}
    const typename gsMatrix<T>::constColumn & deriv2      (index_t k) const
    {
        GISMO_ASSERT(this->m_flags & NEED_2ND_DER, "2nd derivative not computed");
        return m_2ndDers.middleCol( k);
    }
    const typename gsMatrix<T>::constColumn   div         (index_t k) const {return m_divs.col(k);}
    const typename gsMatrix<T>::constColumn   curl        (index_t k) const {return m_curls.col(k);}
    const typename gsMatrix<T>::constColumn   laplacian   (index_t k) const {return m_laps.col(k);}

    const          gsMatrix<T>              & divs        ()          const {return m_divs;}
    const          gsMatrix<T>              & curls       ()          const {return m_curls;}
    const          gsMatrix<T>              & laplacians  ()          const {return m_laps;}
protected:
    unsigned           m_flags;
    unsigned           m_geo_flags;
    unsigned           m_parDim;
    unsigned           m_tarDim;
    unsigned           m_spaceDim;

    gsMatrix<unsigned> m_actives;

    gsMatrix<T>        m_values;

    gsMatrix<T>        m_derivs;
    gsMatrix<T>        m_jacobians;
    gsMatrix<T>        m_2ndDers;

    gsMatrix<T>        m_divs;
    gsMatrix<T>        m_curls;
    gsMatrix<T>        m_laps;
};


/**
    \brief actual implementation of the basis evaluators

    \tparam T the coefficient type
    \tparam ParDim the dimension of the parametric domain, usually f$[0,1]^{ParDim}f$
    \tparam TarDim the dimension of the image of the functions in our discrete space
    \tparam geoTransformType the type of mapping between the values and derivatives between the
            parametric domain and the physical domain: for instance composition withj the inverse or
            div conformorming transformation

**/
template <typename T, int ParDim, int TarDim, typename geoTransformType>
class gsGenericBasisEvaluator
    : public gsBasisEvaluator<T>
{
public:
    gsGenericBasisEvaluator(gsBasis<T> *(&basis)[TarDim],unsigned flags);
public:
    virtual void setFlags   ( unsigned newFlags);
    virtual void evaluateAt ( const gsMatrix<T> &points);
    virtual void evaluateAt ( const gsMatrix<T> &points, const gsGeometryEvaluator<T> &geoEval);
protected:
    int                         m_max_deriv;
    gsBasis<T>                 *m_basis[TarDim];
    gsMatrix<T>                 m_basis_vals[TarDim];
    unsigned                    m_active_shift[TarDim];

    using  gsBasisEvaluator<T>::m_flags;
    using  gsBasisEvaluator<T>::m_geo_flags;
    using  gsBasisEvaluator<T>::m_parDim;
    using  gsBasisEvaluator<T>::m_tarDim;
    using  gsBasisEvaluator<T>::m_spaceDim;

    using  gsBasisEvaluator<T>::m_actives;

    using  gsBasisEvaluator<T>::m_values;

    using  gsBasisEvaluator<T>::m_derivs;
    using  gsBasisEvaluator<T>::m_jacobians;
    using  gsBasisEvaluator<T>::m_2ndDers;

    using  gsBasisEvaluator<T>::m_divs;
    using  gsBasisEvaluator<T>::m_curls;
    using  gsBasisEvaluator<T>::m_laps;
};



} // namespace gismo

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBasisEvaluator.hpp)
#endif
