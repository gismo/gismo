/** @file gsFieldCreator.h

    @brief Provides declaration of gsFieldCreator functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsField.h>

namespace gismo
{


/**
   @brief Generates a field with value the absolute difference (error)
   between and isogeometric function and a function defined on the
   physical domain

   For a multipatch geometry, use the static method of gismo::gsFieldCreator

   \ingroup Core
*/
template<class T>
class gsAbsError : public gsFunction<T> 
{
public:
    /// Shared pointer for gsAbsError
    typedef memory::shared_ptr< gsAbsError > Ptr;

    /// Unique pointer for gsAbsError
    typedef memory::unique_ptr< gsAbsError > uPtr;

    gsAbsError( gsFunction<T> const & f1, 
                gsGeometry<T> const & geo, 
                gsFunction<T> const & f2,
                bool _f2param = false)
    : m_geo(geo), m_f1(f1), m_f2(f2), f2param(_f2param)
    { 
            
    }

    GISMO_CLONE_FUNCTION(gsAbsError)

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_geo.eval_into(u, tmp);

        if(!f2param)
            result.noalias() = ( m_f1.eval(u) - m_f2.eval(tmp) ).cwiseAbs();
        else
            result.noalias() = ( m_f1.eval(u) - m_f2.eval(u) ).cwiseAbs(); // f2 parametric function

    }

    short_t domainDim() const { return m_f1.domainDim(); }
    short_t targetDim() const { return m_f1.targetDim(); }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Absolute error.\n"; return os; };
private:
    const gsGeometry<T> & m_geo;

    const gsFunction<T> & m_f1;
    const gsFunction<T> & m_f2;

    bool f2param;

private:
//    gsAbsError() { }
};


/**
   @brief Generates a field with value being the gradient of an isogeometric function

   For a multipatch geometry, use the static method of gismo::gsFieldCreator

   \ingroup Core
*/
template<class T>
class gsGradientField : public gsFunction<T> 
{
public:
    /// Shared pointer for gsGradientField
    typedef memory::shared_ptr< gsGradientField > Ptr;

    /// Unique pointer for gsGradientField
    typedef memory::unique_ptr< gsGradientField > uPtr;

    gsGradientField( gsGeometry<T> const & geo, gsFunction<T> const & f)
    : m_geo(geo), m_f(f)
    { }

    GISMO_CLONE_FUNCTION(gsGradientField)

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_geo.eval_into(u, tmp);
        m_f.deriv_into(tmp, result);
    }

    short_t domainDim() const { return m_geo.parDim(); }
    short_t targetDim() const { return m_geo.geoDim(); }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Gradient field of "; return os; };
private:
    const gsGeometry<T> & m_geo;
    const gsFunction<T> & m_f  ;
};

/**
   @brief Generates a field with value the Jacobian determinant of a geometry

   For a multipatch geometry, use the static method of gismo::gsFieldCreator

   \ingroup Core
*/
template<class T>
class gsJacDetField : public gsFunction<T> 
{
public:
    /// Shared pointer for gsJacDetField
    typedef memory::shared_ptr< gsJacDetField > Ptr;

    /// Unique pointer for gsJacDetField
    typedef memory::unique_ptr< gsJacDetField > uPtr;

    gsJacDetField( gsGeometry<T> const & geo)
    : m_geo(geo), m_dim(m_geo.parDim())
    { 
        GISMO_ENSURE(m_dim == m_geo.geoDim(), "Not extended to surface case yet." );
    }

    GISMO_CLONE_FUNCTION(gsJacDetField)

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1, u.cols() );
        gsMatrix<T> tmp;
        for (index_t k=0; k!= u.cols(); ++k)
            result(0,k) = m_geo.deriv(u.col(k)).reshape(m_dim,m_dim).determinant();
    }

    short_t domainDim() const { return m_geo.parDim(); }
    short_t targetDim() const { return 1; }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Jacobian determinant field of "; return os; };
private:
    const gsGeometry<T> & m_geo;
    const int m_dim;
};

/**
   @brief Generates the normal field of a geometry

   For a multipatch geometry, use the static method of gismo::gsFieldCreator

   \ingroup Core
*/
template<class T>
class gsNormalField : public gsFunction<T> 
{
public:
    /// Shared pointer for gsNormalField
    typedef memory::shared_ptr< gsNormalField > Ptr;

    /// Unique pointer for gsNormalField
    typedef memory::unique_ptr< gsNormalField > uPtr;

    gsNormalField( gsGeometry<T> const & geo) : m_geo(geo)
    { 
            
    }

    GISMO_CLONE_FUNCTION(gsNormalField)

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        const short_t ParDim = m_geo.parDim();
        result.resize(ParDim+1, u.cols()) ;

        gsMatrix<T> Jk;

        for( index_t j=0; j < u.cols(); ++j )
        {
            Jk = m_geo.jacobian( u.col(j) );
            T alt_sgn(1.0);
            gsMatrix<T> mm(ParDim,ParDim);
            for (int i = 0; i <= ParDim; ++i) // for all components of the normal vector
            {
                Jk.rowMinor(i, mm);
                result(i,j) = alt_sgn * mm.determinant();
                alt_sgn = -alt_sgn;
            }
        }
    }

    short_t domainDim() const { return m_geo.parDim(); }
    short_t targetDim() const { return m_geo.geoDim();}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "NormalField"; return os; };
private:
    const gsGeometry<T> & m_geo;

private:
//    gsNormalField() { }
};

/**
   @brief Generates a field that attaches the parameter values on each
   physical point

   For a multipatch geometry, use the static method of gismo::gsFieldCreator

   \ingroup Core
*/
template<class T>
class gsParamField : public gsFunction<T> 
{
public:
    /// Shared pointer for gsParamField
    typedef memory::shared_ptr< gsParamField > Ptr;

    /// Unique pointer for gsParamField
    typedef memory::unique_ptr< gsParamField > uPtr;

    gsParamField(gsGeometry<T> const & geo)
    : m_geo(geo)
    { 
            
    }

    GISMO_CLONE_FUNCTION(gsParamField)

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { result = u; }

    short_t domainDim() const { return m_geo.domainDim(); }
    short_t targetDim() const { return m_geo.domainDim(); }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Parameter field.\n"; return os; };

private:
    const gsGeometry<T> & m_geo;
};

/**
   @brief Generates a field that indicates the boundary sides on the geometry

   For a multipatch geometry, use the static method of gismo::gsFieldCreator

   In Paraview, choose blot as a color map and Wireframe as representation

   \ingroup Core
*/
template<class T>
class gsBoundaryField : public gsFunction<T> 
{
public:
    /// Shared pointer for gsBoundaryField
    typedef memory::shared_ptr< gsBoundaryField > Ptr;

    /// Unique pointer for gsBoundaryField
    typedef memory::unique_ptr< gsBoundaryField > uPtr;

    explicit gsBoundaryField(gsGeometry<T> const & geo_)
    : geo(geo_), m_supp(geo.support())
    { }

    GISMO_CLONE_FUNCTION(gsBoundaryField)

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { 
        result.setZero(1, u.cols() );

        const int d = geo.parDim();
        for (boxSide c=boxSide::getFirst(d); c<boxSide::getEnd(d); ++c)
        {
            const index_t dir = c.direction();
            const T par = m_supp(dir, c.parameter() );

            for (index_t v = 0; v != u.cols(); ++v) // for all columns of u
            {
                if ( math::abs( u(dir,v) - par ) < 1e-2 )
                    result(0,v) = static_cast<T>(c);
            }
        }
    }

    short_t domainDim() const { return geo.domainDim(); }
    short_t targetDim() const { return 1; }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Boundary side indicator field.\n"; return os; };

private:
    const gsGeometry<T> & geo;
    gsMatrix<T>           m_supp;
};



/**
   @brief Class that creates standard fields on a given parametric
   (multipatch) geometry.

   \ingroup Core
*/
template<class T>
struct gsFieldCreator
{

    static gsField<T> absError(gsField<T> const & field, gsFunction<T> const & f, bool fparam = false)
    {
        const gsMultiPatch<T> & mp = field.patches();
    
        gsPiecewiseFunction<T> * nFields = new gsPiecewiseFunction<T>(mp.nPatches());
        
        for (size_t k=0; k< mp.nPatches(); ++k)
            nFields->addPiecePointer( new gsAbsError<T>(field.function(k), mp.patch(k), f, fparam) );
        
        return gsField<T>(mp, typename gsPiecewiseFunction<T>::Ptr(nFields), true );
    }
    
    
    static gsField<T> gradient(gsMultiPatch<T> const & mp, gsFunction<T> const & f)
    {
        gsPiecewiseFunction<T> * nFields = new gsPiecewiseFunction<T>(mp.nPatches());
        
        
        for (size_t k=0; k< mp.nPatches(); ++k)
            nFields->addPiecePointer( new gsGradientField<T>(mp.patch(k), f) );
        
        return gsField<T>(mp, typename gsPiecewiseFunction<T>::Ptr(nFields), true );
    }

    static gsField<T> normal(gsMultiPatch<T> const & mp)
    {
        gsPiecewiseFunction<T> * nFields = new gsPiecewiseFunction<T>(mp.nPatches());
        
        
        for (size_t k=0; k< mp.nPatches(); ++k)
            nFields->addPiecePointer( new gsNormalField<T>(mp.patch(k)) );
        
        return gsField<T>(mp, typename gsPiecewiseFunction<T>::Ptr(nFields), true );
    }


    static gsField<T> jacDet(gsMultiPatch<T> const & mp)
    {
        gsPiecewiseFunction<T> * nFields = new gsPiecewiseFunction<T>(mp.nPatches());
        
        
        for (size_t k=0; k< mp.nPatches(); ++k)
            nFields->addPiecePointer( new gsJacDetField<T>(mp.patch(k)) );

        return gsField<T>(mp, typename gsPiecewiseFunction<T>::Ptr(nFields), true );
    }

    static gsField<T> parameters(gsMultiPatch<T> const & mp)
    {
        gsPiecewiseFunction<T> * nFields = new gsPiecewiseFunction<T>(mp.nPatches());
        
        
        for (size_t k=0; k< mp.nPatches(); ++k)
            nFields->addPiecePointer( new gsParamField<T>(mp.patch(k)) );

        return gsField<T>(mp, typename gsPiecewiseFunction<T>::Ptr(nFields), true );
    }

    static gsField<T> boundarySides(gsMultiPatch<T> const & mp)
    {
        gsPiecewiseFunction<T> * nFields = new gsPiecewiseFunction<T>(mp.nPatches());
        for (size_t k=0; k< mp.nPatches(); ++k)
            nFields->addPiecePointer( new gsBoundaryField<T>(mp.patch(k)) );
        
        return gsField<T>(mp, typename gsPiecewiseFunction<T>::Ptr(nFields), true );
    }

}; // struct gsFieldCreator



} // namespace gismo
