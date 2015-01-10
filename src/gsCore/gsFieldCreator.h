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
   @brief AbsError

*/

template<class T>
class gsAbsError : public gismo::gsFunction<T> 
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


template<class T>
class gsGradientField : public gismo::gsFunction<T> 
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
        m_f.deriv_into(tmp, result);
        result.resize( m_geo.geoDim(), u.cols() );
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

template<class T>
class gsJacDetField : public gismo::gsFunction<T> 
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


template<class T>
class gsNormalField : public gismo::gsFunction<T> 
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
   @brief Class that creates standard fields on a given parametric
   (multipatch) geometry.
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

}; // struct gsFieldCreator



}; // namespace gismo
