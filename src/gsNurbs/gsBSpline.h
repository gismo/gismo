/** @file gsBSpline.h

    @brief Represents a B-spline curve/function with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBSplineSolver.h>

#include <gsUtils/gsPointGrid.h>


namespace gismo
{
    
/** \brief
    A B-spline function of one argument, with arbitrary target dimension.

    This is the geometry type associated with gsBSplineBasis.
    
    \tparam T coefficient type

    \ingroup geometry
    \ingroup Nurbs
*/
    
    
template<class T>
class gsBSpline : public gsGeoTraits<1,T>::GeometryBase
{
public: 
    typedef gsKnotVector<T> KnotVectorType;

    typedef typename gsGeoTraits<1,T>::GeometryBase Base;

    typedef gsBSplineBasis<T> Basis;

    /// Shared pointer for gsBSpline
    typedef memory::shared_ptr< gsBSpline > Ptr;

    /// Unique pointer for gsBSpline
    typedef memory::unique_ptr< gsBSpline > uPtr;
    
public:
    
    /// Default empty constructor.
    gsBSpline() { }

    // enable swap
    using gsGeometry<T>::swap;
    
    /// Construct B-Spline by basis and coefficient matrix.
    gsBSpline( const Basis & basis, gsMatrix<T> coefs )
    : Base( basis, give(coefs) )
    {
        if( this->basis().isPeriodic() )
            this->basis().expandCoefs(m_coefs);
    }
    
    /// Construct B-Spline by a knot vector and coefficient matrix.
    gsBSpline(KnotVectorType KV, gsMatrix<T> coefs, bool periodic = false )
    {
        this->m_basis = new Basis(give(KV));
        m_coefs.swap(coefs);
            
        if( periodic )
        {
            const index_t sz = this->basis().size();

            this->basis().setPeriodic();

            if ( m_coefs.rows() == sz )
            {
                this->basis().trimCoefs(m_coefs);
            }
            else if  ( m_coefs.rows() == this->basis().size() )
            {
                this->basis().expandCoefs(m_coefs);
            }
            else
            {
                GISMO_ERROR("Wrong number of coefficients for periodic basis.");
            }
        }
        else // non-periodic
        {
            if( this->m_coefs.rows() + KV.degree() + 1 != static_cast<int>( KV.size() ) )
                gsWarn << "gsBSpline Warning: #Knots="<< KV.size()<< ", #coefs="<< this->m_coefs.rows() <<"\n";
        }
    }
    

    /// @brief Construct a B-spline from an interval and knot vector specification.
    /// \param u0 starting parameter
    /// \param u1 end parameter
    /// \param interior number of interior knots
    /// \param degree degree of the spline space
    /// \param coefs coefficients of the spline space
    /// \param mult_interior multiplicity at the interior knots
    /// \param periodic specifies whether the B-spline is periodic
    ///
    /// \ingroup Nurbs
    gsBSpline(T u0, T u1, unsigned interior, int degree,
              gsMatrix<T> coefs, unsigned mult_interior=1, bool periodic = false)
    {
        this->m_basis = new Basis(u0, u1, interior, degree, mult_interior, periodic );
        if( periodic )
        {
            GISMO_ASSERT( this->basis().numCrossingFunctions() == 1,
                          "Knot vector with clamped knots should lead to m_periodic == 1.");

            this->m_coefs = this->basis().perCoefs(coefs);
            GISMO_ASSERT( this->m_coefs.rows() >= this->basis().size(),
                          "You should give at least as many coefficients as the size of the basis after converting to periodic.");
            if( this->m_coefs.rows() > this->basis().size() )
                gsWarn << "Some of the last coefficients are going to be lost during conversion to periodic basis.\n";
        }
        else // non-periodic
        {
            this->m_coefs.swap(coefs);
            GISMO_ASSERT( this->m_coefs.rows() == this->basis().size(),
                          "Number of coefficients does not match the size of the basis.");
        }
    }

    GISMO_CLONE_FUNCTION(gsBSpline)
    
    GISMO_BASIS_ACCESSORS    
    
public:
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "BSpline curve "<< "of degree "<< this->basis().degree()<< ", "<<  this->basis().knots() <<".\n";
        os << "with control points "<< this->m_coefs.row(0)<<" ... "<<this->m_coefs.bottomRows(1) << ".\n";
        if( this->basis().isPeriodic() )
            os << "Periodic with overlay " << this->basis().numCrossingFunctions() << ".\n";
        return os;
    }
    
    /// Returns the starting value of the domain of the basis
    T domainStart() const { return this->basis().knots().first(); }
    
    /// Returns the end value of the domain of the basis
    T domainEnd() const { return this->basis().knots().last(); }
    
    /// Returns a reference to the knot vector
    KnotVectorType & knots(const int i = 0) 
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT( i==0, "Requested knots of invalid direction "<< i );
        return this->basis().knots();
    }

    /// Returns a (const )reference to the knot vector
    const KnotVectorType & knots(const int i = 0) const 
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT( i==0, "Requested knots of invalid direction "<< i );
        return this->basis().knots(); 
    }

    /// Returns the degree of the B-spline
    short_t degree(short_t i = 0) const
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT( i==0, "Requested knots of invalid direction "<< i );
        return this->basis().degree(); 
    }
    
    // compatible curves: same degree, same first/last p+1 knots
    void isCompatible( gsGeometry<T> * )
    { GISMO_NO_IMPLEMENTATION }
    
    // compatible curves: same degree, same first/last p+1 knots
    void makeCompatible( gsGeometry<T> * )
    { GISMO_NO_IMPLEMENTATION }


    /// Merge other B-spline into this one.
    void merge( gsGeometry<T> * otherG );

    /// Insert the given new knot (multiplicity \a i) without changing the curve.
    void insertKnot( T knot, int i = 1)
    {
        if (i==0) return;
        //if ( i==1)
        //single knot insertion: Boehm's algorithm
        //else
        //knot with multiplicity:   Oslo algorithm
        if( this->basis().isPeriodic() )
        {
            int borderKnotMult = this->basis().borderKnotMult();
            KnotVectorType & knots = this->knots();
            unsigned deg = this->basis().degree();


            GISMO_ASSERT( knot != knots[deg] && knot != knots[knots.size() - deg - 1],
                          "You are trying to increase the multiplicity of the p+1st knot but the code is not ready for that.\n");

            // If we would be inserting to "passive" regions, we rather insert the knot into the mirrored part.
            // Adjustment of the mirrored knots is then desirable.
            if( knot < knots[deg - borderKnotMult + 1] )
            {
                knot += this->basis()._activeLength();
            }
            else if( knot > knots[knots.size() - deg + borderKnotMult - 2] )
            {
                knot -= this->basis()._activeLength();
            }
            // If necessary, we update the mirrored part of the knot vector.
            if((knot < knots[2*deg + 1 - borderKnotMult]) || (knot >= knots[knots.size() - 2*deg - 2 + borderKnotMult]))
                this->basis().enforceOuterKnotsPeriodic();

            // We copy some of the control points to pretend for a while that the basis is not periodic.
            //gsMatrix<T> trueCoefs = this->basis().perCoefs( this->coefs() );
            gsBoehm( this->basis().knots(), this->coefs(), knot, i );
            //this->coefs() = trueCoefs;
            //this->coefs().conservativeResize( this->basis().size(), this->coefs().cols() );
        }
        else // non-periodic
            gsBoehm( this->basis().knots(), this->coefs() , knot, i);
    }

    /// Insert the given new knots in the range \a [\em inBegin .. \em inEend )
    /// without changing the curve.
    /// \param inBegin iterator pointing to the first knot to be inserted
    /// \param inEnd iterator pointing to one position after the last
    /// knot to be inserted
    template <class It>
    void insertKnots( It inBegin, It inEnd)
    {
        if( this->basis().isPeriodic() )
        {
            // We assume we got valid (i.e., non-NULL) iterators; I don't think we have a reasonable way to test it in GISMO_ASSERT.

            GISMO_ASSERT( (*inBegin > *(this->knots().begin()))
                          && (*(inEnd-1) < *(this->knots().end()-1)),
                          "Please, ask me to insert knots inside the knot interval only." );

            std::vector<T> newKnots(inBegin,inEnd);

            // It can happen that user would ask us to insert knots outside the active range.
            // Should it happen, we shift the values and then sort the knot vector to remain non-decreasing.

            T activeLength = this->basis()._activeLength();
            T blue1 = *(this->basis().knots().begin()  + this->degree() );
            T blue2 = *(this->basis().knots().end()    - this->degree() );
            typename std::vector<T>::iterator begin = newKnots.begin();
            typename std::vector<T>::iterator end   = newKnots.end();
            for( typename std::vector<T>::iterator it = begin; it != end; ++it )
            {
                if( *it < blue1 )
                    *it += activeLength;
                else if( *it > blue2 )
                    *it -= activeLength;
                else if( (*it == blue1) || (*it == blue2) )
                {
                    // I would like to erase it and go on but for that the iterator is not enough.
                    // ANGELOS, wouldn't it be better to put a GISMO_ASSERT here instead? Or throw an exception?
                    gsWarn << "Interrupting insertKnots(): Cannot insert value equal to the p+1st knot from the beginning or end into a periodic curve.\n";
                    return;
                }
            }
            /* // This way it would be nicer but it does not cover all the dubious cases.
            for( It it = begin; *it < blue1; ++it )
                *it += activeLength;
            for( It it = end-1; *it > blue2; --it )
                *it -= activeLength;*/

            std::sort( begin, end );

            // We copy some of the control points to pretend for a while that the basis is not periodic.
            //gsMatrix<T> trueCoefs = this->basis().perCoefs( this->coefs() );
            gsBoehmRefine( this->basis().knots(), this->coefs(), this->degree(), begin, end );
            //this->basis().enforceOuterKnotsPeriodic();
            //this->coefs() = trueCoefs; // Isn't this memory-leaking?
            // Periodic basis is restored using the so-called Student's Algorithm: forget as much as possible =).
            //this->coefs().conservativeResize( this->basis().size(), this->coefs().cols() );
        }
        else // non-periodic
            gsBoehmRefine( this->basis().knots(), this->coefs(), this->degree(), inBegin, inEnd);
    }

    // Look at gsGeometry class for a description
    void degreeElevate(short_t const i = 1, short_t const dir = -1);

    /// @brief Returns true iff the point p is contained (approximately) on
    /// the curve, with the given tolerance.
    // Under Construction..( TO DO )
    bool contains( gsMatrix<T> const & p, T const & tol = 1e-6 )
        {
            assert( p.cols()==1 );
            gsBSplineSolver<T> slv;
            std::vector<T> roots;
            short_t dim = this->geoDim();
            gsMatrix<T> ev, tmp(1,1);
            int i(1);

            slv.allRoots( *this, roots, 0, p(0,0), 1e-6, 100);
            for ( typename std::vector<T>::const_iterator it=roots.begin();
                  it != roots.end(); ++it)
            {
                tmp(0,0)= *it;
                this->eval_into(tmp,ev);
                for (i=1; i!=dim; ++i )
                    if ( math::abs( ev(i,0) - p(i,0) ) > tol ) 
                        break;
                if (i==dim) 
                    return true;
            }
            return false;
        }

    /// Sample \a npoints uniformly distributed (in parameter domain) points on the curve.
    gsMatrix<T> sample(int npoints = 50) const
    {      
        gsMatrix<T> images;
        gsMatrix<T> interval = this->parameterRange();
        const gsMatrix<T> pts = gsPointGrid(interval, npoints );
        this->eval_into( pts, images );
        return images;
    }

    /// Tries to convert the curve into periodic.
    void setPeriodic(bool flag = true)
    {
        this->basis().setPeriodic(flag);
        this->m_coefs = this->basis().perCoefs( this->m_coefs );
    }

    inline int numCoefs() const { return this->m_coefs.rows() - this->basis().numCrossingFunctions(); }

    inline bool isPeriodic() const { return this->basis().isPeriodic(); }


    /// Returns true if this curve is closed
    bool isClosed(T tol = 1e-10) const 
    { 
        return ( this->basis().isPeriodic() || 
                 (m_coefs.row(0) - m_coefs.row(m_coefs.rows()-1)).squaredNorm()<tol );  
    }

    /// \brief Return true if point \a u is on the curve with
    /// tolerance \a tol
    bool isOn(gsMatrix<T> const &u, T tol = 1e-3) const;

    /// \brief Return true if point \a u is an endpoint (corner) of
    /// the curve with tolerance \a tol
    bool isPatchCorner(gsMatrix<T> const &v, T tol = 1e-3) const;

    /// \brief returns the tensor-index \a curr of the corner control
    /// point \a v, or an invalid index if the corner is not found
    /// within the tolerance \a tol
    void findCorner(const gsMatrix<T>   & v, 
                    gsVector<index_t,1> & curr,
                    T tol = 1e-3);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the starting point of the curve. Assumes that \a v is
    /// either the starting or the ending point of the curve
    void setOriginCorner(gsMatrix<T> const &v);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the ending point of the curve. Assumes that \a v is
    /// either the starting or the ending point of the curve
    void setFurthestCorner(gsMatrix<T> const &v);

    void swapDirections(const unsigned i, const unsigned j);
    
protected:
    
    using Base::m_coefs;
    
    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...

}; // class gsBSpline


/*
// Product of spline functions
template<class T>
gsBSpline<T> operator*(const gsBSpline<T> & lhs, const gsBSpline<T> & rhs) 
{
    // Note: quick-fix implementation. Better algorithms  exist
    gsKnotVector<T> kv = lhs.knots().knotUnion( lhs.knots() );
    kv.degreeElevate( lhs.degree() + rhs.degree() - kv.degree() );

    gsMatrix<T> pts;
    kv.greville_into(pts);
    const gsMatrix<T> ev = (lhs.eval(pts)).array() * (rhs.eval(pts)).array();
    
    // fixme: avoid temporaries here
    return *(static_cast<gsBSpline<T>*>(gsBSplineBasis<T>(kv).interpolateData(ev,pts)));
}
*/

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBSpline.hpp)
#endif
