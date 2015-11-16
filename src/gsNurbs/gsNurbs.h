/** @file gsNurbs.h

    @brief Represents a NURBS curve/function with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>

#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBoehm.h>


namespace gismo
{

/** \brief
    A NURBS function of one argument, with arbitrary target dimension.

    This is the geometry type associated with gsNurbsBasis.

    \tparam T coefficient type

    \ingroup geometry
    \ingroup Nurbs
*/

template<class T, class KnotVectorType>
class gsNurbs :  public gsGeoTraits<1,T>::GeometryBase
{

public: 
    typedef typename gsGeoTraits<1,T>::GeometryBase Base;

    /// Coefficient type
    typedef T Scalar_t;

    typedef gsRationalBasis<gsBSplineBasis<T> > RBasis;
    typedef gsNurbsBasis<T,KnotVectorType> Basis;

    /// Shared pointer for gsNurbs
    typedef memory::shared_ptr< gsNurbs<T> > Ptr;

public:

    /// Default empty constructor
    gsNurbs() : Base() { }
    
    /// Construct NURBS by NURBS basis functions and coefficient matrix
    gsNurbs( const Basis & basis, const gsMatrix<T> & coefs ) :
    Base(basis, coefs) { }

    /// Construct NURBS by NURBS basis functions and coefficient matrix
    gsNurbs( const Basis & basis, gsMovable< gsMatrix<T> > coefs ) :
    Base( basis, coefs) { }

    /// Construct B-Spline by a knot vector, degree, weights and coefficient matrix
    gsNurbs( const gsKnotVector<T>& KV, const gsMatrix<T> & w, const gsMatrix<T> & coefs )
    {
        this->m_coefs = coefs;
        this->m_basis = new Basis(KV, w);
    }
    
    /// Construct B-Spline by a knot vector, degree, weights and coefficient matrix
    gsNurbs( const gsKnotVector<T>& KV, gsMovable< gsMatrix<T> > w, gsMovable< gsMatrix<T> > coefs )
    {
        this->m_coefs = coefs;
        this->m_basis = new Basis(KV, w);
    }
    
    /// Construct B-Spline by a knot vector, degree and projective coefficient matrix
    gsNurbs( const gsKnotVector<T>& KV, const gsMatrix<T> * pcoefs ) :
    Base( new Basis(KV, & pcoefs->rightCols(1) ), pcoefs) 
    {
        // TO DO: divide pcoefs by the weights
    }
    
    /// Clone function. Used to make a copy of the (derived) geometry
    virtual gsNurbs * clone() const
    { return new gsNurbs(*this); }

    GISMO_BASIS_ACCESSORS

public:

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "NURBS curve "<< "of degree "<< this->basis().degree()
         << " over knots "<< this->basis().knots() <<",\n";
        os << "weights: ["<< this->weights().transpose()<< " ]\n ";
        os << "with control points "<< this->m_coefs << ".\n"; 
        return os; 
    };


//////////////////////////////////////////////////
// Additional members for univariate B-Splines
//////////////////////////////////////////////////

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return this->basis().knots().first(); };

    /// Returns the starting value of the domain of the basis
    T domainEnd() const { return this->basis().knots().last(); };

    /// Access to i-th weight
    T & weight(int i) const { return this->basis().weight(i); }

    /// Returns the weights of the rational basis
    const gsMatrix<T> & weights() const { return this->basis().weights(); }
    gsMatrix<T> & weights()       { return this->basis().weights(); }

    void isCompatible( gsGeometry<T> * other )
    {
        // compatible curves: same degree, same first/last p+1 knots
    };
    void makeCompatible( gsGeometry<T> * other )
    {
        // compatible curves: same degree, same first/last p+1 knots
    };

    void merge( gsGeometry<T> * other )
    {
        //check for BSpline
        gsNurbs<T> *  bother = static_cast<gsNurbs<T> *>( other );

        //check degree
        if ( this->basis().degree() != bother->basis().degree() )
        {std::cout<<"gsNurbs: Cannot merge with different degree curve"<<"\n"; return;}
            
        //check geometric dimension

        //check that it touches *this curve

        // merge knot vector
        this->basis().knots().merge( bother->basis().knots() ) ;

        // merge coefficients
        int n= this->coefsSize();
        this->m_coefs.conservativeResize( n + bother->coefsSize() -1, Eigen::NoChange ) ;
            
        this->m_coefs.block( n,0,bother->coefsSize()-1,bother->geoDim() ) =
            bother->m_coefs.block( 1,0,bother->coefsSize()-1,bother->geoDim() ) ;

        // merge Weights
        this->weights().conservativeResize( n + bother->coefsSize() -1, Eigen::NoChange ) ;
            
        this->weights().block( n,0,bother->coefsSize()-1,1 ) =
            bother->weights().block( 1,0,bother->coefsSize()-1, 1 ) ;
    };
    
    /// Insert the given new knot (multiplicity \a i) without changing the curve.
    void insertKnot( T knot, int i = 1 )
    {
        assert( i>0);
        // TO DO: There is also Oslo Algorithm and others

        std::cout <<"before \n"<< this->m_coefs.transpose() <<"\n";

        gsMatrix<T> tmp( this->m_coefs.rows(), this->m_coefs.cols()+1 );
        tmp.leftCols(this->m_coefs.cols()) =  
            this->weights().asDiagonal() * this->m_coefs;
        tmp.rightCols(1) = this->weights();

        gsBoehm( this->basis().knots(), tmp, knot, i); 

        int l = tmp.cols() -1;
        this->basis().setWeights( tmp.rightCols( 1 ) );
        for ( index_t k = 0; k< this->m_coefs.rows(); ++k)
            for ( index_t j = 0; j< this->m_coefs.cols(); ++j)
                tmp(k,j) /= tmp(k,l);           
        //tmp.col(j).cwiseQuotient( *this->basis().m_weights );
            
        this->m_coefs =  tmp.leftCols( this->m_coefs.cols() );
    }

    
    void elevate( )
    {
        // uses knot insertion
        //gsMatrix<T> * c = this->coefs();
    };    


    //void toProjective() { m_weights=w; } ;
    

protected:
    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...


// Data members
private:

    bool projective;

}; // class gsNurbs


} // namespace gismo
