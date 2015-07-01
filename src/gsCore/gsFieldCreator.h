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
    gsAbsError( gsFunction<T> const & f1, 
                gsGeometry<T> const & geo, 
                gsFunction<T> const & f2 )
    : m_geo(geo), m_f1(f1), m_f2(f2)
    { 
            
    }

    gsAbsError * clone() const
    { return new gsAbsError(*this); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_geo.eval_into(u, tmp);
        result.noalias() = ( *m_f1.eval(u) - *m_f2.eval(tmp) ).cwiseAbs();
    }

    int domainDim() const { return m_f1.domainDim(); }
    int targetDim() const { return m_f1.targetDim(); }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Absolute error.\n"; return os; };
private:
    const gsGeometry<T> & m_geo;

    const gsFunction<T> & m_f1;
    const gsFunction<T> & m_f2;

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
    gsGradientField( gsGeometry<T> const & geo, gsFunction<T> const & f)
    : m_geo(geo), m_f(f)
    { 
            
    }

    gsGradientField * clone() const
    { return new gsGradientField(*this); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> tmp;
        m_geo.eval_into(u, tmp);
        m_f.newderiv_into(tmp, result);
    }

    int domainDim() const { return m_geo.parDim(); }
    int targetDim() const { return m_geo.geoDim(); }

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
    gsJacDetField( gsGeometry<T> const & geo)
    : m_geo(geo), m_dim(m_geo.parDim())
    { 
        GISMO_ENSURE(m_dim == m_geo.geoDim(), "Not extended to surface case yet." );
    }

    gsJacDetField * clone() const
    { return new gsJacDetField(*this); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize(1, u.cols() );
        gsMatrix<T> tmp;
        for (index_t k=0; k!= u.cols(); ++k)
            result(0,k) = m_geo.deriv(u.col(k))->reshape(m_dim,m_dim).determinant();
    }

    int domainDim() const { return m_geo.parDim(); }
    int targetDim() const { return 1; }

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
    gsNormalField( gsGeometry<T> const & geo) : m_geo(geo)
    { 
            
    }

    gsNormalField * clone() const
    { return new gsNormalField(*this); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        const int ParDim = m_geo.parDim();
        result.resize(ParDim+1, u.cols()) ;

        gsMatrix<T> Jk;

        for( index_t j=0; j < u.cols(); ++j )
        {
            Jk = m_geo.jac( u.col(j) );
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

    int domainDim() const { return m_geo.parDim(); }
    int targetDim() const { return m_geo.geoDim();}

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
    gsParamField(gsGeometry<T> const & geo)
    : m_geo(geo)
    { 
            
    }

    gsParamField * clone() const
    { return new gsParamField(*this); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { result = u; }

    int domainDim() const { return m_geo.domainDim(); }
    int targetDim() const { return m_geo.domainDim(); }

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
    gsBoundaryField(gsGeometry<T> const & geo_)
    : geo(geo_), m_supp(geo.support())
    { }

    gsBoundaryField * clone() const
    { return new gsBoundaryField(*this); }

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

    int domainDim() const { return geo.domainDim(); }
    int targetDim() const { return 1; }

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

    static typename gsField<T>::uPtr absError(gsField<T> const & field, gsFunction<T> const & f)
    {
        const gsMultiPatch<T> & mp = field.patches();
    
        std::vector<gsFunction<T> *> nFields;
        nFields.reserve( mp.nPatches() );
        
        for (unsigned k=0; k< mp.nPatches(); ++k)
            nFields.push_back( new gsAbsError<T>(field.function(k), mp.patch(k), f) );
        
        return typename gsField<T>::uPtr( new gsField<T>(mp, nFields, true ) );
    }
    
    
    static typename gsField<T>::uPtr gradient(gsMultiPatch<T> const & mp, gsFunction<T> const & f)
    {
        std::vector<gsFunction<T> *> nFields;
        nFields.reserve( mp.nPatches() );
        
        for (unsigned k=0; k< mp.nPatches(); ++k)
            nFields.push_back( new gsGradientField<T>(mp.patch(k), f) );
        
        return typename gsField<T>::uPtr( new gsField<T>(mp, nFields, true ) );
    }

    static typename gsField<T>::uPtr normal(gsMultiPatch<T> const & mp)
    {
        std::vector<gsFunction<T> *> nFields;
        nFields.reserve( mp.nPatches() );
        
        for (unsigned k=0; k< mp.nPatches(); ++k)
            nFields.push_back( new gsNormalField<T>(mp.patch(k)) );
        
        return typename gsField<T>::uPtr( new gsField<T>(mp, nFields, true ) );
    }


    static typename gsField<T>::uPtr jacDet(gsMultiPatch<T> const & mp)
    {
        std::vector<gsFunction<T> *> nFields;
        nFields.reserve( mp.nPatches() );
        
        for (unsigned k=0; k< mp.nPatches(); ++k)
            nFields.push_back( new gsJacDetField<T>(mp.patch(k)) );
        
        return typename gsField<T>::uPtr( new gsField<T>(mp, nFields, true ) );
    }

    static typename gsField<T>::uPtr parameters(gsMultiPatch<T> const & mp)
    {
        std::vector<gsFunction<T> *> nFields;
        nFields.reserve( mp.nPatches() );
        
        for (unsigned k=0; k< mp.nPatches(); ++k)
            nFields.push_back( new gsParamField<T>(mp.patch(k)) );
        
        return typename gsField<T>::uPtr( new gsField<T>(mp, nFields, true ) );
    }

    static typename gsField<T>::uPtr boundarySides(gsMultiPatch<T> const & mp)
    {
        std::vector<gsFunction<T> *> nFields;
        nFields.reserve( mp.nPatches() );
        
        for (unsigned k=0; k< mp.nPatches(); ++k)
            nFields.push_back( new gsBoundaryField<T>(mp.patch(k)) );
        
        return typename gsField<T>::uPtr( new gsField<T>(mp, nFields, true ) );
    }

}; // struct gsFieldCreator



} // namespace gismo
