/** @file gsTensorNurbs.h

    @brief Represents a tensor-product NURBS patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

#include <gsTensor/gsTensorTools.h> // todo: move to hpp

namespace gismo
{

/** \brief 
    A tensor product Non-Uniform Rational B-spline function
    (NURBS) of parametric dimension \em d, with arbitrary target
    dimension.

    This is the geometry type associated with gsTensorNurbsBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector the NURBS bases use

    \ingroup geometry
    \ingroup Nurbs
*/

template<short_t d, class T>
class gsTensorNurbs : public gsGeoTraits<d,T>::GeometryBase
{

public: 
    typedef gsKnotVector<T> KnotVectorType;

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef T Scalar_t;
    
    typedef gsTensorBSplineBasis<d,T> TBasis;      // underlying tensor basis
    
    /// Family type
    typedef gsBSplineBasis<T>  Family_t;
    
    // rational version of tensor basis (basis for this geometry)
    typedef gsTensorNurbsBasis<d,T>   Basis;

    /// Associated boundary geometry type
    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::RatGeometry BoundaryGeometryType;
    //typedef gsTensorNurbs<static_cast<short_t>(d-1),T> BoundaryGeometryType;

    /// Associated boundary basis type
    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::RatBasis BoundaryBasisType;

    /// Shared pointer for gsTensorNurbs
    typedef memory::shared_ptr< gsTensorNurbs > Ptr;

    /// Unique pointer for gsTensorNurbs
    typedef memory::unique_ptr< gsTensorNurbs > uPtr;

public:

    /// Default empty constructor
    gsTensorNurbs() : Base() { }

    gsTensorNurbs( const Basis & basis,  gsMatrix<T> coefs ) :
    Base( basis, give(coefs) ) { }

    /// Construct 2D tensor NURBS by knot vectors, degrees and coefficient matrix
    /// All weights are set to be 1.
    gsTensorNurbs( gsKnotVector<T> const& KV1, gsKnotVector<T> const & KV2,
                   gsMatrix<T> tcoefs)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "<< d
                     <<"D NURBS using 2 knot-vectors.");

        gsBSplineBasis<T>    * Bu    = new gsBSplineBasis<T>(KV1);
        gsBSplineBasis<T>    * Bv    = new gsBSplineBasis<T>(KV2);

        TBasis   *tbasis = new TBasis(Bu,Bv) ;//d==2
      
        this->m_basis = new Basis(tbasis) ;
        this->m_coefs.swap(tcoefs);
        GISMO_ASSERT(tbasis->size()== m_coefs.rows(), 
                     "Coefficient matrix for the NURBS does not have "
                     "the expected number of control points (rows)." );
    }

    /// Construct 2D tensor NURBS by knot vectors, degrees, weights and coefficient matrix
    gsTensorNurbs( gsKnotVector<T> const& KV1, gsKnotVector<T> const & KV2,
                   gsMatrix<T> tcoefs, gsMatrix<T> wgts)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "<< d
                     <<"D NURBS using 2 knot-vectors.");

        gsBSplineBasis<T>    * Bu    = new gsBSplineBasis<T>(KV1);
        gsBSplineBasis<T>    * Bv    = new gsBSplineBasis<T>(KV2);

        TBasis   *tbasis = new TBasis(Bu,Bv) ;//d==2
      
        GISMO_ASSERT(tbasis->size()== tcoefs.rows(), 
                     "Coefficient matrix for the NURBS does not have "
                     "the expected number of control points (rows)." );

        this->m_basis = new Basis(tbasis , give(wgts)) ;
        this->m_coefs.swap(tcoefs);
    }

    /// Construct 3D tensor NURBS by knot vectors, degrees and coefficient matrix
    /// \a tcoefs, \a wgts become empty after the constructor is called
    gsTensorNurbs( gsKnotVector<T> const & KV1, 
                   gsKnotVector<T> const & KV2, 
                   gsKnotVector<T> const & KV3,
                   gsMatrix<T>  tcoefs, 
                   gsMatrix<T> wgts )
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "<< d
                     <<"D NURBS using 3 knot-vectors.");
      
        gsBSplineBasis<T> * Bu= new gsBSplineBasis<T>(KV1);
        gsBSplineBasis<T> * Bv= new gsBSplineBasis<T>(KV2);
        gsBSplineBasis<T> * Bw= new gsBSplineBasis<T>(KV3);
        TBasis *tbasis = new TBasis(Bu,Bv,Bw) ;//d==3
      
        GISMO_ASSERT(tbasis->size()== tcoefs.rows(), 
                     "Coefficient matrix for the NURBS does not have "
                     "the expected number of control points (rows)." );

        this->m_basis = new Basis(tbasis, give(wgts));
        this->m_coefs.swap(tcoefs);
    }

    /// Construct 3D tensor NURBS by knot vectors, degrees and coefficient matrix
    /// All weights are set to be 1.
    gsTensorNurbs( gsKnotVector<T> const & KV1, 
                   gsKnotVector<T> const & KV2, 
                   gsKnotVector<T> const & KV3,
                   gsMatrix<T> tcoefs)
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "<< d
                     <<"D NURBS using 3 knot-vectors.");
      
        gsBSplineBasis<T> * Bu= new gsBSplineBasis<T>(KV1);
        gsBSplineBasis<T> * Bv= new gsBSplineBasis<T>(KV2);
        gsBSplineBasis<T> * Bw= new gsBSplineBasis<T>(KV3);
        TBasis *tbasis = new TBasis(Bu,Bv,Bw) ;//d==3
      
        GISMO_ASSERT(tbasis->size()== tcoefs.rows(), 
                     "Coefficient matrix for the NURBS does not have "
                     "the expected number of control points (rows)." );

        this->m_basis = new Basis(tbasis) ;
        this->m_coefs.swap(tcoefs);
    }

    /// Construct 4D tensor NURBS by knot vectors, degrees and coefficient matrix
    /// \a tcoefs, \a wgts become empty after the constructor is called
    gsTensorNurbs( gsKnotVector<T> const & KV1,
                   gsKnotVector<T> const & KV2,
                   gsKnotVector<T> const & KV3,
                   gsKnotVector<T> const & KV4,
                   gsMatrix<T>  tcoefs,
                   gsMatrix<T> wgts )
    {
        GISMO_ASSERT(d==4, "Wrong dimension: tried to make a "<< d
                                                              <<"D NURBS using 4 knot-vectors.");

        std::vector< gsBSplineBasis<T>* > cbases;
        cbases.reserve(4);
        cbases.push_back(new gsBSplineBasis<T>(KV1) );
        cbases.push_back(new gsBSplineBasis<T>(KV2) );
        cbases.push_back(new gsBSplineBasis<T>(KV3) );
        cbases.push_back(new gsBSplineBasis<T>(KV4) );
        TBasis * tbasis = gsBSplineTraits<d,T>::Basis::New(cbases); //d==4

        GISMO_ASSERT(tbasis->size()== tcoefs.rows(),
                     "Coefficient matrix for the NURBS does not have "
                         "the expected number of control points (rows)." );

        this->m_basis = new Basis(tbasis, give(wgts));
        this->m_coefs.swap(tcoefs);
    }

    GISMO_BASIS_ACCESSORS
    
    public:

// ***********************************************
// Virtual member functions required by the base class
// ***********************************************

    GISMO_CLONE_FUNCTION(gsTensorNurbs)

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "Tensor-NURBS geometry "<< "R^"<< this->parDim() << 
            " --> R^"<< this->geoDim()<< ", #control pnts= "<< this->coefsSize() <<": "
         << this->coef(0) <<" ... "<< this->coef(this->coefsSize()-1); 
        os << "\nweights: "
           << this->basis().weights().at(0) <<" ... "
           << this->basis().weights().at(this->coefsSize()-1)
           <<"\n" ;
        return os; }

// ***********************************************
// Additional members for tensor NURBS
// ***********************************************

    /// Returns a reference to the knot vector \a i
    const KnotVectorType & knots(const int i) const 
    { return this->basis().source().knots(i); } 

    KnotVectorType & knots(const int i) 
    { return this->basis().source().knots(i); } 

    /// Inserts knot \a knot at direction \a dir, \a i times
    void insertKnot( T knot, int dir, int i = 1)
    {
        GISMO_ASSERT( i>0, "multiplicity must be at least 1");
        GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
                      "Invalid basis component "<< dir <<" requested for degree elevation" );
        
        const index_t n = this->m_coefs.cols();
        gsTensorBSplineBasis<d,T> & tbs = this->basis().source();

        gsVector<index_t,d> sz;
        tbs.size_cwise(sz);
        
        swapTensorDirection(0, dir, sz, m_coefs  );
        std::swap(sz[0],sz[dir]);
        swapTensorDirection(0, dir, sz, weights());
        const index_t nc = sz.template tail<static_cast<short_t>(d-1)>().prod();
        m_coefs  .resize( sz[0], n * nc );
        weights().resize( sz[0], nc     );
        
        gsBoehm(tbs.knots(dir), weights(), knot, i, false);
        gsBoehm(tbs.knots(dir), m_coefs  , knot, i, true );
        sz[0] = m_coefs.rows();

        const index_t ncoef = sz.prod();
        m_coefs  .resize(ncoef, n );
        weights().resize(ncoef, 1 );
        swapTensorDirection(0, dir, sz, m_coefs  );
        std::swap(sz[0],sz[dir]);
        swapTensorDirection(0, dir, sz, weights());
    }

    /// Access to i-th weight
    T & weight(int i) const { return this->basis().weight(i); }

    /// Returns the NURBS weights
    gsMatrix<T> & weights() const { return this->basis().weights(); }

    /// Returns the NURBS weights as non-const reference
    gsMatrix<T> & weights() { return this->basis().weights(); }
    
    /// Returns the degree of the basis wrt direction i 
    short_t degree(unsigned i) const
    { return this->basis().source().component(i).degree(); }

/// Toggle orientation wrt coordinate k
    void reverse(unsigned k)
    { 
        GISMO_ASSERT(d==2, "only 2D for now");

        gsVector<int> str(d); 
        gsVector<int> sz (d); 

        gsTensorBSplineBasis<d,T> & tbsbasis = this->basis().source();

        sz[0]  = tbsbasis.component(k ).size();
        sz[1]  = tbsbasis.component(!k).size();
        str[0] = tbsbasis.source().stride( k );
        str[1] = tbsbasis.source().stride(!k );
    
        gsMatrix<T> & w  = tbsbasis.weights();

        for  ( int i=0; i< sz[0]; i++ )
            for  ( int j=0; j< sz[1]/2; j++ )
            {
                this->m_coefs.row(i*str[0] + j*str[1] ).swap(
                    this->m_coefs.row(i*str[0] + (sz[1]-j-1)*str[1] )
                    );

                w.row(i*str[0] + j*str[1] ).swap(
                    w.row(i*str[0] + (sz[1]-j-1)*str[1] )
                    );
            }
        tbsbasis.component(k).reverse();
    }

    /// Constucts an isoparametric slice of this tensor Nurbs curve by fixing
    /// \a par in direction \a dir_fixed. The resulting tensor Nurbs has
    /// one less dimension and is given back in \a result.
    void slice(index_t dir_fixed, T par, BoundaryGeometryType & result) const
    {
        GISMO_ASSERT(d-1>=0,"d must be greater or equal than 1");
        GISMO_ASSERT(dir_fixed>=0 && static_cast<unsigned>(dir_fixed)<d,"cannot fix a dir greater than dim or smaller than 0");
        // construct the d-1 basis
        boxSide side(dir_fixed,0);
        //typename gsTensorNurbsBasis<static_cast<short_t>(d-1),T>::uPtr tbasis = this->basis().boundaryBasis(side);
        typename BoundaryBasisType::uPtr tbasis = this->basis().boundaryBasis(side);

        if(d==1)
        {
            gsMatrix<T> val(1,1),point;
            val(0,0)=par;
            this->eval_into(val,point);
            //result = gsTensorNurbs<static_cast<short_t>(d-1), T>(*tbasis, point );
            result = BoundaryGeometryType(*tbasis, point );
        }
        else
        {
            const int mult   = this->basis().knots(dir_fixed).multiplicity(par);
            const int degree = this->basis().degree(dir_fixed);

            gsMatrix<T> coefs;
            if( mult>=degree )
            {
                // no knot insertion needed, just extract the right coefficients
                const gsKnotVector<T>& knots = this->basis().knots(dir_fixed);
                const index_t index = (knots.iFind(par) - knots.begin()) - this->basis().degree(dir_fixed);
                gsVector<index_t,d> sizes;
                this->basis().size_cwise(sizes);
                constructCoefsForSlice<d, T>(dir_fixed, index, this->coefs(), sizes, coefs);
            }
            else
            {
                // clone the basis and inserting up to degree knots at par
                gsTensorNurbs<d,T>* clone = this->clone().release();

                gsVector<index_t,d> intStrides;
                this->basis().stride_cwise(intStrides);
                gsTensorBoehm(
                    clone->basis().knots(dir_fixed),clone->coefs(),par,dir_fixed,
                    intStrides.template cast<unsigned>(), degree-mult,true);

                // extract right ceofficients
                const gsKnotVector<T>& knots = clone->basis().knots(dir_fixed);
                const index_t index = (knots.iFind(par) - knots.begin()) - clone->basis().degree(dir_fixed);
                gsVector<index_t,d> sizes;
                clone->basis().size_cwise(sizes);
                constructCoefsForSlice<d, T>(dir_fixed, index, clone->coefs(), sizes, coefs);
                delete clone;
            }

            // construct the object
            //result = gsTensorBSpline<static_cast<short_t>(d-1),T>(*tbasis, give(coefs) );
            //result = BoundaryGeometry(*tbasis, give(coefs) );
            result = BoundaryGeometryType (*tbasis, coefs );
        }
    }

    void swapDirections(const unsigned i, const unsigned j)
    {
        gsVector<index_t,d> sz;
        this->basis().size_cwise(sz);
        swapTensorDirection(i, j, sz, m_coefs);
        this->basis().swapDirections(i,j);
    }



protected:
    // todo: check function: check the coefficient number, degree, knot vector ...

    using gsGeometry<T>::m_coefs;
    using gsGeometry<T>::m_basis;

}; // class gsTensorNurbs


// ***********************************************
// ***********************************************


} // namespace gismo
