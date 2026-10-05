/** @file gsTensorBSpline.hpp

    @brief Provides implementation of a tensor-product B-spline patch
    of arbitrary dimension

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConstantFunction.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSplineAlgorithms.h>
#include <gsNurbs/gsBoehm.h>

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

#include <gsTensor/gsTensorTools.h>

namespace gismo
{

template<short_t d, class T>
void constructCoefsForSlice(index_t dir_fixed, index_t index,
                            const gsMatrix<T> & fullCoefs,
                            const gsVector<index_t,d> & sizes,
                            gsMatrix<T>& result)
{
    gsVector<index_t,d> lowerCorner, upperCorner;
    lowerCorner.setZero();
    upperCorner = sizes;
    lowerCorner[dir_fixed] = index;
    upperCorner[dir_fixed] = index + 1;
    // to do: gsMatrix<index_t> ind = gsTensorBasis::coefSlice(dim_fixed, index) ?

    // Collect the boundary coefficients
    result.resize( sizes.prod() / sizes[dir_fixed], fullCoefs.cols() );
    gsVector<index_t,d> str, cur = lowerCorner;
    tensorStrides(sizes,str);
    index_t r = 0;
    do {
        result.row(r++) = fullCoefs.row( cur.dot(str) );
    } while ( nextLexicographic(cur, lowerCorner, upperCorner) );
}

template<short_t d, class T>
gsTensorBSpline<d,T>::gsTensorBSpline(gsMatrix<T> const & corner,
                                      KnotVectorType KV1, KnotVectorType KV2)
{
    GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "<< d<<"D tensor B-spline using 2 knot-vectors.");

    std::vector<Family_t*> cbases;
    const index_t n1 = KV1.size() - KV1.degree() - 1;
    const index_t n2 = KV2.size() - KV2.degree() - 1;

    cbases.push_back(new gsBSplineBasis<T>(give(KV1)) );
    cbases.push_back(new gsBSplineBasis<T>(give(KV2)) );
    Basis * tbasis = Basis::New(cbases); //d==2


    GISMO_ASSERT( (corner.rows()==4) && (corner.cols()==3),
                  "gsTensorBSpline: Please make sure that the size of *corner* is 4-by-3");

    gsMatrix<T> pcp (n1*n2, 3);
    // set up CPs on boundary first. The inner CPs on each boundary curve are
    // uniformly linear dependent on the two corner CPs
    index_t j=0; // boundary v=0
    for (index_t i=0; i<=n1-1; i++)
    {
        for (index_t xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(0,xi) + i/((T)(n1-1))*( corner(1,xi) - corner(0,xi) );
        }
    }
    j=n2-1; // boundary v=1
    for (index_t i=0; i<=n1-1; i++)
    {
        for (index_t xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(3,xi) + i/((T)(n1-1))*( corner(2,xi) - corner(3,xi) );
        }
    }
    index_t i=0; // boundary u=0;
    for (j=0; j<=n2-1; j++)
    {
        for (index_t xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(0,xi) + j/((T)(n2-1))*( corner(3,xi) - corner(0,xi) );
        }
    }
    i=n1-1; // boundary u=1;
    for (j=0; j<=n2-1; j++)
    {
        for (index_t xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(1,xi) + j/((T)(n2-1))*( corner(2,xi) - corner(1,xi) );
        }
    }
    // uniformly linear dependent in horizontal direction
    for (j=1; j<=n2-2; j++)
    {
        for (i=1; i<=n1-2; i++)
        {
            for (index_t xi=0; xi<=2; xi++) //specification of x or y or z
            {
                pcp(i+j*n1,xi)=pcp(0+j*n1,xi) + i/((T)(n1-1))*( pcp(n1-1+j*n1,xi)-pcp(0+j*n1,xi) );
            }
        }
    }

    this->m_basis = tbasis;
    this->m_coefs.swap( pcp );
}

// todo: move to hpp
template<short_t d, class T>
void gsTensorBSpline<d,T>::slice(index_t dir_fixed,T par,
                                 BoundaryGeometryType & result) const
{
    GISMO_ASSERT(d-1>=0,"d must be greater or equal than 1");
    GISMO_ASSERT(dir_fixed>=0 && static_cast<unsigned>(dir_fixed)<d,"cannot fix a dir greater than dim or smaller than 0");
    // construct the d-1 basis
    boxSide side(dir_fixed,0);
    typename BoundaryBasisType::uPtr tbasis = this->basis().boundaryBasis(side);

    if(d==1)
    {
        gsMatrix<T> val(1,1),point;
        val(0,0)=par;
        this->eval_into(val,point);
        result = BoundaryGeometryType(*tbasis, point );
    }
    else
    {
        const index_t mult   = this->basis().knots(dir_fixed).multiplicity(par);
        const short_t degree = this->basis().degree(dir_fixed);
        index_t index;
        gsMatrix<T> coefs;
        if( mult>=degree )
        {
            // no knot insertion needed, just extract the right coefficients
            const gsKnotVector<T>& knots = this->basis().knots(dir_fixed);
            if (par == *knots.domainEnd())
                index = knots.size()-this->basis().degree(dir_fixed)-2;
            else
                index = (knots.iFind(par) - knots.begin()) - this->basis().degree(dir_fixed);
            gsVector<index_t,d> sizes;
            this->basis().size_cwise(sizes);
            constructCoefsForSlice<d, T>(dir_fixed, index, this->coefs(), sizes, coefs);
        }
        else
        {
            // clone the basis and inserting up to degree knots at par
            gsTensorBSpline<d,T>* clone = this->clone().release();

            gsVector<index_t,d> intStrides;
            this->basis().stride_cwise(intStrides);
            gsTensorBoehm(
                clone->basis().knots(dir_fixed),clone->coefs(),par,dir_fixed,
                intStrides.template cast<size_t>(), degree-mult,true);

            // extract right ceofficients
            const gsKnotVector<T>& knots = clone->basis().knots(dir_fixed);
            if (par == *knots.domainEnd())
                index = knots.size()-clone->basis().degree(dir_fixed)-2;
            else
                index = (knots.iFind(par) - knots.begin()) - clone->basis().degree(dir_fixed);
            gsVector<index_t,d> sizes;
            clone->basis().size_cwise(sizes);
            constructCoefsForSlice<d, T>(dir_fixed, index, clone->coefs(), sizes, coefs);
            delete clone;
        }

        // construct the object
        //result = gsTensorBSpline<static_cast<short_t>(d-1),T>(*tbasis, give(coefs) );
        //result = BoundaryGeometry(*tbasis, give(coefs) );
        result = BoundaryGeometryType(*tbasis, coefs );
    }
}

template<short_t d, class T>
void gsTensorBSpline<d,T>::reverse(unsigned k)
{
    gsTensorBSplineBasis<d,T> & tbsbasis = this->basis();
    gsVector<index_t,d> sz;
    tbsbasis.size_cwise(sz);
    flipTensorVector(k, sz, m_coefs);
    tbsbasis.component(k).reverse();
}


template<short_t d, class T>
void gsTensorBSpline<d,T>::swapDirections(const unsigned i, const unsigned j)
{
    gsVector<index_t,d> sz;
    this->basis().size_cwise(sz);
    swapTensorDirection(i, j, sz, m_coefs);
    this->basis().swapDirections(i,j);
}

template<short_t d, class T>
void gsTensorBSpline<d,T>::toggleOrientation()
{
    swapDirections(0,1);
}


template<short_t d, class T>
bool gsTensorBSpline<d,T>::isPatchCorner(gsMatrix<T> const &v, T tol) const
{
    gsVector<index_t,d> str(d), vupp(d), curr = gsVector<index_t,d>::Zero(d);
    this->basis().stride_cwise(str);
    this->basis().size_cwise(vupp);
    vupp.array() -= 1;

    do // loop over all vertices
    {
        if ( (v - m_coefs.row(curr.dot(str))).squaredNorm() < tol )
            return true;
    }
    while ( nextCubeVertex(curr, vupp) );

    return false;
}

template<short_t d, class T>
void gsTensorBSpline<d,T>::findCorner(const gsMatrix<T> & v,
                                      gsVector<index_t,d> & curr,
                                      T tol)
{
    gsVector<index_t,d> str , // Tensor strides
        vupp; // Furthest corner

    this->basis().stride_cwise(str);
    this->basis().size_cwise(vupp);
    vupp.array() -= 1;

    curr.setZero();
    do // loop over all vertices
    {
        if ( (v - m_coefs.row(curr.dot(str))).squaredNorm() < tol )
            return;
    }
    while ( nextCubeVertex(curr, vupp) );

    // Corner not found, Invalidate the result
    vupp.array() += 1;
    curr.swap(vupp);
    gsWarn<<"Point "<< v <<" is not an corner of the patch. (Call isPatchCorner() first!).\n";
}

template<short_t d, class T>
void gsTensorBSpline<d,T>::setOriginCorner(gsMatrix<T> const &v)
{
    gsVector<index_t,d> curr;
    findCorner(v, curr);
    if ( curr[0] == this->basis().size(0) )
        return;
    for(unsigned k = 0; k!=d; ++k)
        if ( curr[k] != 0 )
            this->reverse(k);
}

template<short_t d, class T>
void gsTensorBSpline<d,T>::setFurthestCorner(gsMatrix<T> const &v)
{
    gsVector<index_t,d> curr;
    findCorner(v, curr);
    if ( curr[0] == this->basis().size(0) )
        return;
    for(unsigned k = 0; k!=d; ++k)
        if ( curr[k] == 0 )
            this->reverse(k);
}


template<short_t d, class T>
void gsTensorBSpline<d,T>::degreeElevate(short_t const i, short_t const dir)
{
    if (dir == -1)
    {
        for (short_t j = 0; j < d; ++j)
            degreeElevate(i, j);
        return;
    }

    if (knots(dir).numLeftGhosts() != 0 || knots(dir).numRightGhosts() != 0)
    {
        gsGeometry<T>::degreeElevate(i, dir);
        return;
    }

    GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
                  "Invalid basis component "<< dir <<" requested for degree elevation" );

    const index_t n = this->m_coefs.cols();

    gsVector<index_t,d> sz;
    this->basis().size_cwise(sz);

    swapTensorDirection(0, dir, sz, this->m_coefs);
    this->m_coefs.resize( sz[0], n * sz.template tail<static_cast<short_t>(d-1)>().prod() );

    bspline::degreeElevateBSpline(this->basis().component(dir), this->m_coefs, i);
    sz[0] = this->m_coefs.rows();

    this->m_coefs.resize( sz.prod(), n );
    swapTensorDirection(0, dir, sz, this->m_coefs);
}

template<short_t d, class T>
void gsTensorBSpline<d,T>::insertKnot( T knot, int dir, int i)
{
    GISMO_ASSERT( i>0, "multiplicity must be at least 1");


    GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
                  "Invalid basis component "<< dir <<" requested for degree elevation" );

    const index_t n = this->m_coefs.cols();

    gsVector<index_t,d> sz;
    this->basis().size_cwise(sz);

    swapTensorDirection(0, dir, sz, this->m_coefs);
    this->m_coefs.resize( sz[0], n * sz.template tail<static_cast<short_t>(d-1)>().prod() );

    gsBoehm( this->basis().component(dir).knots(), this->coefs() , knot, i);
    sz[0] = this->m_coefs.rows();

    this->m_coefs.resize( sz.prod(), n );
    swapTensorDirection(0, dir, sz, this->m_coefs);
}


template<short_t d, class T>
index_t gsTensorBSpline<d,T>::removeKnot( T knot, short_t dir, short_t i, T tol)
{
    GISMO_ASSERT( i > 0, "Must remove at least once.");
    GISMO_ASSERT( dir >= 0 && dir < d,
                  "Invalid direction " << dir );

    gsVector<index_t,d> sz;
    this->basis().size_cwise(sz);

    gsVector<index_t,d> intStrides;
    this->basis().stride_cwise(intStrides);

    return gsTensorKnotRemove<T>(
        this->basis().component(dir).knots(),
        this->coefs(),
        knot,
        dir,
        intStrides.template cast<unsigned>(),
        i,
        true,
        tol);
}


template<short_t d, class T>
typename gsGeometry<T>::uPtr gsTensorBSpline<d,T>::localRep(const gsMatrix<T> & u) const
{
    std::vector<KnotVectorType> kv(d); // the local knot-vectors
    gsVector<index_t,d> cfirst, clast; // tensor-indices of local coefficients
    index_t sz = 1; // number of control points in the local representation

    // Fill in the data defined above
    for(unsigned i = 0; i!=d; ++i)
    {
        const short_t deg = degree(i);
        typename KnotVectorType::const_iterator span = knots(i).iFind(u(i,0));

        sz       *= deg + 1;
        clast[i]  = span - knots(i).begin();
        cfirst[i] = clast[i] - deg;
        kv[i]     = KnotVectorType(deg, span - deg, span + deg + 2);
    }

    // Collect the local coefficients
    const gsMatrix<T> & allCoefs = this->coefs();
    gsMatrix<T> coefs(sz, allCoefs.cols() );
    gsVector<index_t,d> str, cur = cfirst;
    basis().stride_cwise(str);
    index_t r = 0;
    do {
        coefs.row(r++) = allCoefs.row( cur.dot(str) );
    } while ( nextCubePoint(cur, cfirst, clast) );

    // All set, return the local representation
    return Basis(kv).makeGeometry(give(coefs));
}


template<short_t d, class T>
std::ostream & gsTensorBSpline<d,T>::print(std::ostream &os) const
{
    os << "Tensor BSpline geometry "<< "R^"<< d <<
        " --> R^"<< this->geoDim()
       << ", #control pnts= "<< this->coefsSize();
    if ( m_coefs.size() )
        os << ": "<< this->coef(0) <<" ... "<< this->coef(this->coefsSize()-1);
    if ( m_basis )
        os<<"\nBasis:\n" << this->basis() ;
    return os;
}

template<short_t d, class T>
std::vector<gsGeometry<T>* > gsTensorBSpline<d,T>::uniformSplit(index_t dir) const
{
    // 1. insert p+1 in all directions
    // 2. recover 2^d patches
    GISMO_ASSERT( (dir > -2) && (dir < static_cast<index_t>(d)),
                  "Invalid basis component "<< dir <<" requested for geometry splitting" );
    std::vector<gsGeometry<T>* > result_temp, result;
    gsVector<T> midpoints;
    if(dir==-1)
    {
        result.reserve(math::exp2(d));
        midpoints.setZero(d);

        for(unsigned i=0; i<d;++i)
            midpoints(i)= (basis().knots(i).sbegin().value() + (--basis().knots(i).send()).value())/ (T)(2);

        for(unsigned i=0; i<d;++i)
        {
            result_temp.clear();

            //one could uniform the if-statement and the for-loop by setting result[0] = this,
            //however, the const prevents this.
            if(result.size()==0)
            {
                gsTensorBSpline<d,T>* left = new gsTensorBSpline<d,T>();
                gsTensorBSpline<d,T>* right = new gsTensorBSpline<d,T>();
                this->splitAt(i,midpoints(i),*left,*right);
                result_temp.push_back(left);
                result_temp.push_back(right);
            }
            for(size_t j=0; j<result.size();j++)
            {
                gsTensorBSpline<d,T>* left = new gsTensorBSpline<d,T>();
                gsTensorBSpline<d,T>* right = new gsTensorBSpline<d,T>();
                static_cast<gsTensorBSpline<d,T>*>(result[j])->splitAt(i,midpoints(i),*left,*right);

                result_temp.push_back(left);
                result_temp.push_back(right);
            }


            freeAll(result);
            result = result_temp;
        }
    }
    else
    {
        result.reserve(2);
        T xi =  (basis().knots(dir).sbegin().value() + (--basis().knots(dir).send()).value())/(T)(2);
        gsTensorBSpline<d,T>* left = new gsTensorBSpline<d,T>();
        gsTensorBSpline<d,T>* right = new gsTensorBSpline<d,T>();

        splitAt(dir,xi,*left,*right);

        result.push_back(left);
        result.push_back(right);
    }
    return result;

}


template<short_t d, class T>
void gsTensorBSpline<d,T>::splitAt( index_t dir,T xi, gsTensorBSpline<d,T>& left,  gsTensorBSpline<d,T>& right) const
{
    GISMO_ASSERT( (dir >= 0) && (dir < static_cast<index_t>(d)),
                  "Invalid basis component "<< dir <<" requested for geometry splitting" );

    GISMO_ASSERT(basis().knots(dir).sbegin().value()<xi && xi< (--basis().knots(dir).send()).value() , "splitting point "<<xi<<" not in the knotvector");

    //First make a copy of the actual object, to allow const
    gsTensorBSpline<d,T> copy(*this);

    // Extract a reference to the knots, the basis and coefs of the copy
    KnotVectorType & knots = copy.basis().knots(dir);
    gsTensorBSplineBasis<d,T> & base = copy.basis();

    // some constants
    const int p = base.degree(dir);                      // degree
    const index_t mult = p + 1 - knots.multiplicity(xi); // multiplicity

    //insert the knot, such that its multiplicity is p+1
    if (mult>0)
        copy.insertKnot(xi, dir, mult);

    //swap the direction dir with 0, to be able to extract the coefs.
    copy.swapDirections(0,dir);

    gsMatrix<T> & coefs = copy.coefs();
    const index_t tDim  = coefs.cols();

    //some more constants
    gsVector<index_t,d> sizes;                    // number of coefs in each dir
    base.size_cwise(sizes);
    const index_t sz = sizes.prod();          // total number of coefs

    //find the number of coefs left from xi (in direction 0)
    const index_t nL = knots.uFind(xi).firstAppearance();
    index_t nR = base.size(0) - nL;

    //Split the coefficients
    gsMatrix<T> coefL, coefR;
    coefL.setZero(sizes.tail(d-1).prod()*(nL), tDim);
    coefR.setZero(sz-coefL.rows(), tDim);

    index_t kL,kR,i;
    i=kL=kR=0;
    while(i<sz)
    {
        coefL.block(kL,0,nL, tDim) = coefs.block(i,0,nL, tDim);
        coefR.block(kR,0,nR, tDim) = coefs.block(i+nL,0,nR, tDim);

        kL+=nL;
        kR+=nR;

        i+= nL + nR;
    }

    //build up the new geometries
    //build the knot vector for direction 0 (swapped!)
    typename KnotVectorType::iterator it = knots.iFind(xi);
    typename KnotVectorType::knotContainer matL(knots.begin(),++it);
    it-=p+1; // move the iterator to the beginning of the inserted knots
    typename KnotVectorType::knotContainer matR(it, knots.end());
    KnotVectorType knotsL(give(matL),p);
    KnotVectorType knotsR(give(matR),p);

    // rescale the splitted knot vector (not mandatory)
    // knotsL.affineTransformTo(0,1);
    // knotsR.affineTransformTo(0,1);

    //collect the other directions
    std::vector<KnotVectorType> KVL, KVR;
    KVL.push_back(knotsL);
    KVR.push_back(knotsR);
    for(i=1; i<static_cast<index_t>(d);++i)
    {
        KVL.push_back(base.knots(i));
        KVR.push_back(base.knots(i));
    }

    //finally the two new geometries
    left  = gsTensorBSpline<d,T>(Basis(give(KVL)), give(coefL));
    left.swapDirections(0,dir);
    right = gsTensorBSpline<d,T>(Basis(give(KVR)), give(coefR));
    right.swapDirections(0,dir);
}


template<short_t d, class T>
std::vector<gsGeometry<T>* >
gsTensorBSpline<d,T>::splitAtMult(index_t minMult, index_t dir) const
{
    GISMO_ASSERT( (dir >= -1) && (dir < static_cast<index_t>(d)),
                  "Invalid basis component "<< dir <<" requested for splitting" );
    std::vector<gsGeometry<T>* > result;

    if (-1==dir)
    {
        std::vector<gsGeometry<T>* > tmpi, tmp;
        result = this->splitAtMult(minMult,0);
        for(short_t i=1; i<d;++i)
        {
            tmp.swap(result);
            result.clear();
            for(size_t j=0; j!=tmp.size();++j)
            {
                tmpi = static_cast<gsTensorBSpline<d,T>*>(tmp[j])
                    ->splitAtMult(minMult,i);
                delete tmp[j];
                result.insert( result.end(), tmpi.begin(), tmpi.end() );
            }
        }
    }
    else
    {
        gsTensorBSpline<d,T> * tmp = new gsTensorBSpline<d,T>(*this);
        //iterate over knots
        for (typename KnotVectorType::uiterator it = knots(dir).ubegin()+1;
             it!=knots(dir).uend()-1; ++it)
        {
            if (it.multiplicity()>=minMult)
            {
                gsTensorBSpline<d,T> * o = new gsTensorBSpline<d,T>();
                tmp->splitAt(dir,*it,*o,*tmp);
                result.push_back(o);
            }
        }
        result.push_back(tmp);
    }
    return result;
}


template<short_t d, class T>
typename gsGeometry<T>::uPtr
gsTensorBSpline<d,T>:: iface(const boundaryInterface & bi,
                             const gsGeometry<T> & other) const
{
    // Grab boundary control point indices in matching configuration
    gsMatrix<index_t> bdr0, bdr1;
    this->basis().matchWith(bi, other.basis(), bdr0, bdr1);

    //from here: Assume linear curves, merge control points (todo: add option for this)
    index_t b[2];
    b[0]=b[1]=0;
    std::list<std::pair<const gsMatrix<T> *,index_t> > cv;//patch,cp-index
    //maybe: check if both ifaces are identical using a flag...

    cv.push_back( std::make_pair(&this->coefs(), bdr0.at(b[0]++) ) );
    do {
        T dist0=(cv.back().first->row(cv.back().second)-this->coef(bdr0.at(b[0]))).squaredNorm();
        if ( 0 == dist0 ) { b[0]++; continue; } //skip double point
        T dist1=(cv.back().first->row(cv.back().second)-other.coef(bdr1.at(b[1]))).squaredNorm();
        if ( 0 == dist1 ) { b[1]++; continue; } //skip double point

        // gsDebugVar(dist0);
        // gsDebugVar(dist1);
        // gsDebugVar( this->coef(bdr0.at(b[0]+1)) );
        // gsDebugVar( other.coef(bdr1.at(b[1]))   );
        if (dist0>dist1)
            cv.push_back( std::make_pair(&other.coefs(), bdr1.at(b[1]++) ) );
        else
            cv.push_back( std::make_pair(&this->coefs(), bdr0.at(b[0]++) ) );

    } while ( b[0]!=bdr0.size() && b[1]!=bdr1.size() );

    //gsDebugVar(cv.size());
    while ( b[0]<bdr0.size() )
        cv.push_back( std::make_pair(&this->coefs(), bdr0.at(b[0]++) ) );

    //gsDebugVar(cv.size());
    while ( b[1]<bdr1.size() )
        cv.push_back( std::make_pair(&other.coefs(), bdr1.at(b[1]++) ) );

    //gsDebugVar(cv.size());

    // temporary fix: the last point is always doubled.
    cv.pop_back();

    // Construct interface geometry using cv and uniform knots (polyline)
    gsMatrix<T> cf(cv.size(),this->geoDim());
    index_t c = 0;
    for(typename std::list<std::pair<const gsMatrix<T> *,index_t> >::iterator
            it = cv.begin(); it!=cv.end(); ++it)
        cf.row(c++) = it->first->row(it->second);
    gsKnotVector<T> kv(0,1,c-2,2,1);
    gsBSplineBasis<T> bs(kv);
    return bs.makeGeometry(cf);
}


// ============================================================================
//  Calculus methods
// ============================================================================

// ----------------------------------------------------------------------------
//  toBezier
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::toBezier() const
{
    gsTensorBSpline<d,T> result(*this);
    for (short_t k = 0; k < d; ++k)
    {
        const short_t p = result.degree(k);
        const KnotVectorType & kv0 = result.knots(k);
        const T first = kv0.first();
        const T last  = kv0.last();

        // Collect unique interior knots BEFORE any insertion (avoids invalidation)
        std::vector<T> interior;
        interior.reserve(static_cast<size_t>(kv0.uSize()));
        for (size_t i = 0; i < kv0.uSize(); ++i)
        {
            const T xi = kv0.uValue(i);
            if (xi > first && xi < last)
                interior.push_back(xi);
        }

        for (const T xi : interior)
        {
            const int mult = result.knots(k).multiplicity(xi);
            const int need = p + 1 - mult;
            if (need > 0)
                result.insertKnot(xi, k, need);
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
//  multiply / squared / cubed  (Bernstein product formula)
// ----------------------------------------------------------------------------

namespace {
/// Binomial coefficient C(n,k) as a compile-time-safe integer computation
inline index_t binom(index_t n, index_t k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    index_t result = 1;
    for (index_t i = 0; i < k; ++i)
    {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}
} // anonymous namespace

template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::multiply(const gsTensorBSpline<d,T> & A,
                                                    const gsTensorBSpline<d,T> & B,
                                                    bool keepBezier)
{
    GISMO_ASSERT(A.targetDim() == 1 && B.targetDim() == 1,
                 "multiply() requires scalar (targetDim==1) splines, got "
                 << A.targetDim() << " and " << B.targetDim() << ".");

    // 1. Breakpoints of the product and the continuity it inherits.
    //
    //    Where A is C^a and B is C^b the product AB is C^{min(a,b)}, so the
    //    minimal space is S(pA+pB, C^{min(a,b)}), i.e. interior multiplicity
    //    pA+pB-min(a,b).  This has to be read off the ORIGINAL knot vectors:
    //    the Bézier extraction in step 2 sets every interior multiplicity to
    //    degree+1 and destroys the continuity information.
    std::vector< std::vector<T> >   breaks(d);
    std::vector< std::vector<index_t> > targetMult(d);
    for (short_t k = 0; k < d; ++k)
    {
        const KnotVectorType & kvA = A.knots(k);
        const KnotVectorType & kvB = B.knots(k);
        GISMO_ASSERT(kvA.first() == kvB.first() && kvA.last() == kvB.last(),
                     "multiply(): the two factors live on different parametric "
                     "domains in direction " << k << ".");

        const short_t pA = A.degree(k), pB = B.degree(k);

        std::vector<T> & u = breaks[k];
        for (size_t i = 0; i < kvA.uSize(); ++i)
        {
            const T xi = kvA.uValue(i);
            if (xi > kvA.first() && xi < kvA.last()) u.push_back(xi);
        }
        for (size_t i = 0; i < kvB.uSize(); ++i)
        {
            const T xi = kvB.uValue(i);
            if (xi > kvB.first() && xi < kvB.last()) u.push_back(xi);
        }
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), u.end()), u.end());

        targetMult[k].reserve(u.size());
        for (size_t i = 0; i < u.size(); ++i)
        {
            const index_t mA = kvA.has(u[i]) ? static_cast<index_t>(kvA.multiplicity(u[i])) : 0;
            const index_t mB = kvB.has(u[i]) ? static_cast<index_t>(kvB.multiplicity(u[i])) : 0;
            // A factor with no knot at u[i] is a polynomial across it, hence
            // C^inf; capping its continuity at its degree instead would
            // understate the compression whenever the two degrees differ.
            const index_t cA = (mA == 0) ? pB : pA - mA;
            const index_t cB = (mB == 0) ? pA : pB - mB;
            targetMult[k].push_back(pA + pB - std::min(cA, cB));
        }
    }

    // 2. Bézier extraction of both factors on the common breakpoint set, so
    //    that the two coefficient grids have the same element structure.
    gsTensorBSpline<d,T> bA(A), bB(B);
    for (short_t k = 0; k < d; ++k)
    {
        const index_t pA = bA.degree(k), pB = bB.degree(k);
        for (size_t i = 0; i < breaks[k].size(); ++i)
        {
            const T xi = breaks[k][i];
            const index_t mA = bA.knots(k).has(xi) ? static_cast<index_t>(bA.knots(k).multiplicity(xi)) : 0;
            if (pA + 1 - mA > 0) bA.insertKnot(xi, k, pA + 1 - mA);
            const index_t mB = bB.knots(k).has(xi) ? static_cast<index_t>(bB.knots(k).multiplicity(xi)) : 0;
            if (pB + 1 - mB > 0) bB.insertKnot(xi, k, pB + 1 - mB);
        }
    }

    gsVector<index_t,d> degA, degB, degR;
    for (short_t k = 0; k < d; ++k)
    {
        degA[k] = static_cast<index_t>(bA.degree(k));
        degB[k] = static_cast<index_t>(bB.degree(k));
        degR[k] = degA[k] + degB[k];
    }

    // 3. Result basis: every breakpoint at multiplicity degR+1 (Bézier).
    std::vector<KnotVectorType> kvs_result(d);
    for (short_t k = 0; k < d; ++k)
    {
        const KnotVectorType & kv = bA.knots(k);
        std::vector<T> knots_result;
        knots_result.reserve(kv.uSize() * (degR[k] + 1));
        for (size_t i = 0; i < kv.uSize(); ++i)
        {
            const T xi = kv.uValue(i);
            for (index_t m = 0; m < degR[k] + 1; ++m)
                knots_result.push_back(xi);
        }
        kvs_result[k] = KnotVectorType(knots_result, static_cast<short_t>(degR[k]));
    }
    std::vector<Family_t*> cbases;
    cbases.reserve(d);
    for (short_t k = 0; k < d; ++k)
        cbases.push_back(new gsBSplineBasis<T>(kvs_result[k]));
    Basis * result_basis = Basis::New(cbases);

    // 4. Sizes and strides of the two Bézier grids and of the result.
    gsVector<index_t,d> szA, szB;
    bA.basis().size_cwise(szA);
    bB.basis().size_cwise(szB);

    gsVector<index_t,d> n_elems;
    for (short_t k = 0; k < d; ++k)
    {
        n_elems[k] = szA[k] / (degA[k] + 1);
        GISMO_ASSERT(n_elems[k] == szB[k] / (degB[k] + 1),
                     "multiply(): element counts differ in direction " << k
                     << " after Bézier extraction (" << n_elems[k] << " vs "
                     << szB[k] / (degB[k] + 1) << ").");
    }

    gsVector<index_t,d> strA, strB, res_sz, strR;
    strA[0] = 1; strB[0] = 1; strR[0] = 1;
    for (short_t k = 0; k < d; ++k)
        res_sz[k] = n_elems[k] * (degR[k] + 1);
    for (short_t k = 1; k < d; ++k)
    {
        strA[k] = strA[k-1] * szA[k-1];
        strB[k] = strB[k-1] * szB[k-1];
        strR[k] = strR[k-1] * res_sz[k-1];
    }

    gsMatrix<T> res_coefs(res_sz.prod(), 1);
    res_coefs.setZero();

    const gsMatrix<T> & ca = bA.coefs();
    const gsMatrix<T> & cb = bB.coefs();

    // 5. Element-by-element Bernstein convolution (see the class doc).
    gsVector<index_t,d> elem_idx;
    elem_idx.setZero();
    do
    {
        gsVector<index_t,d> k_idx;
        k_idx.setZero();
        do
        {
            T val = T(0);
            gsVector<index_t,d> i_idx;
            i_idx.setZero();
            do
            {
                bool valid = true;
                gsVector<index_t,d> j_idx;
                T weight = T(1);
                for (short_t kd = 0; kd < d; ++kd)
                {
                    j_idx[kd] = k_idx[kd] - i_idx[kd];
                    if (j_idx[kd] < 0 || j_idx[kd] > degB[kd])
                    { valid = false; break; }
                    weight *= static_cast<T>(binom(degA[kd], i_idx[kd]))
                            * static_cast<T>(binom(degB[kd], j_idx[kd]))
                            / static_cast<T>(binom(degR[kd], k_idx[kd]));
                }
                if (valid)
                {
                    index_t fi = 0, fj = 0;
                    for (short_t kd = 0; kd < d; ++kd)
                    {
                        fi += (elem_idx[kd] * (degA[kd]+1) + i_idx[kd]) * strA[kd];
                        fj += (elem_idx[kd] * (degB[kd]+1) + j_idx[kd]) * strB[kd];
                    }
                    val += weight * ca(fi, 0) * cb(fj, 0);
                }

                short_t carry = 0;
                do
                {
                    ++i_idx[carry];
                    if (i_idx[carry] > degA[carry]) { i_idx[carry] = 0; ++carry; }
                    else break;
                } while (carry < d);
                if (carry == d) break;
            } while (true);

            index_t fr = 0;
            for (short_t kd = 0; kd < d; ++kd)
                fr += (elem_idx[kd] * (degR[kd]+1) + k_idx[kd]) * strR[kd];
            res_coefs(fr, 0) = val;

            short_t carry = 0;
            do
            {
                ++k_idx[carry];
                if (k_idx[carry] > degR[carry]) { k_idx[carry] = 0; ++carry; }
                else break;
            } while (carry < d);
            if (carry == d) break;
        } while (true);

        short_t carry = 0;
        do
        {
            ++elem_idx[carry];
            if (elem_idx[carry] >= n_elems[carry]) { elem_idx[carry] = 0; ++carry; }
            else break;
        } while (carry < d);
        if (carry == d) break;
     } while (true);

    // The geometry constructor clones the basis, so the one built above is
    // ours to release.
    gsTensorBSpline<d,T> result(*result_basis, give(res_coefs));
    delete result_basis;

    // 6. Compress back to the minimal space determined in step 1.
    if (!keepBezier)
    {
        for (short_t k = 0; k < d; ++k)
            for (size_t i = 0; i < breaks[k].size(); ++i)
            {
                const T xi = breaks[k][i];
                index_t mult = static_cast<index_t>(result.knots(k).multiplicity(xi));
                while (mult > targetMult[k][i])
                {
                    if (result.removeKnot(xi, k, 1) == 0) break;
                    --mult;
                }
            }
    }
    return result;
}

template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::squared(bool keepBezier) const
{
    GISMO_ASSERT(this->targetDim() == 1,
                 "squared() is only implemented for scalar (targetDim==1) splines.");
    return multiply(*this, *this, keepBezier);
}

template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::cubed(bool keepBezier) const
{
    GISMO_ASSERT(this->targetDim() == 1,
                 "cubed() is only implemented for scalar (targetDim==1) splines.");
    // The square is taken in its minimal space rather than in Bézier form:
    // multiply() derives the continuity of the product from the knot
    // multiplicities of its arguments, and a Bézier square would advertise
    // C^{-1} and so block all compression of c³.
    return multiply(this->squared(false), *this, keepBezier);
}

// ----------------------------------------------------------------------------
//  grad(dir) — partial derivative in direction dir
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::grad(short_t dir) const
{
    GISMO_ASSERT(dir >= 0 && static_cast<unsigned>(dir) < d,
                 "grad: invalid direction " << dir);
    const short_t p = this->degree(dir);
    GISMO_ASSERT(p >= 1, "grad: cannot differentiate a degree-0 spline.");

    const index_t nc = this->m_coefs.cols();  // number of components

    gsVector<index_t,d> sz;
    this->basis().size_cwise(sz);             // sz[k] = number of B-spline DOFs in dir k

    // Swap target direction to position 0 (same pattern as insertKnot / degreeElevate)
    gsMatrix<T> coefs = this->m_coefs;
    swapTensorDirection(0, dir, sz, coefs);
    // After swap: sz[0] holds the size in the direction we want to differentiate.
    // Flatten all other directions into columns.
    const index_t n_other = sz.template tail<static_cast<short_t>(d-1)>().prod();
    coefs.resize(sz[0], nc * n_other);

    // Apply B-spline derivative formula row-by-row (each row = one "fiber")
    //   d_i = p * (c_{i+1} - c_i) / (t_{i+p+1} - t_{i+1})
    const index_t n_dir = sz[0];
    const KnotVectorType & kv = this->knots(dir);

    gsMatrix<T> dcoefs(n_dir - 1, nc * n_other);
    for (index_t i = 0; i < n_dir - 1; ++i)
    {
        const T denom = kv[i + p + 1] - kv[i + 1];
        if (denom == T(0))
            dcoefs.row(i).setZero();
        else
            dcoefs.row(i) = (T(p) / denom) * (coefs.row(i + 1) - coefs.row(i));
    }

    // Unflatten back to tensor layout with updated size in the derivative direction
    sz[0] = n_dir - 1;
    dcoefs.resize(sz.prod(), nc);
    swapTensorDirection(0, dir, sz, dcoefs);

    // Build result knot vector in direction dir: drop first and last knot
    // (standard formula: KV of derivative basis = KV[1:-1])
    std::vector<KnotVectorType> kvs(d);
    for (short_t k = 0; k < d; ++k)
    {
        if (k == dir)
        {
            const KnotVectorType & kv_full = this->knots(dir);
            const index_t sz_full = kv_full.size();
            // New knot vector: skip first and last knot, degree = p-1
            std::vector<T> kv_new;
            kv_new.reserve(static_cast<size_t>(sz_full - 2));
            for (index_t i = 1; i < sz_full - 1; ++i)
                kv_new.push_back(kv_full[i]);
            kvs[k] = KnotVectorType(kv_new, static_cast<short_t>(p - 1));
        }
        else
            kvs[k] = this->knots(k);
    }

    std::vector<Family_t*> cbases;
    cbases.reserve(d);
    for (short_t k = 0; k < d; ++k)
        cbases.push_back(new gsBSplineBasis<T>(kvs[k]));
    Basis * dbasis = Basis::New(cbases);
    return gsTensorBSpline<d,T>(*dbasis, give(dcoefs));
}

// ----------------------------------------------------------------------------
//  grad() — all partial derivatives
// ----------------------------------------------------------------------------
template<short_t d, class T>
std::vector< gsTensorBSpline<d,T> > gsTensorBSpline<d,T>::grad() const
{
    std::vector< gsTensorBSpline<d,T> > result;
    result.reserve(d);
    for (short_t k = 0; k < d; ++k)
        result.push_back(this->grad(k));
    return result;
}

// ----------------------------------------------------------------------------
//  makeCompatible / linearCombination — common-space arithmetic
// ----------------------------------------------------------------------------
template<short_t d, class T>
void gsTensorBSpline<d,T>::makeCompatible(std::vector< gsTensorBSpline<d,T> > & splines)
{
    if (splines.size() < 2) return;

    for (short_t k = 0; k < d; ++k)
    {
        // 1. Common degree.  Elevation has to precede the knot union because
        //    degreeElevate raises every interior multiplicity by the same
        //    amount in order to preserve continuity, so a union taken before
        //    elevation would be undone by it.
        short_t pmax = 0;
        for (size_t s = 0; s != splines.size(); ++s)
            if (splines[s].degree(k) > pmax) pmax = splines[s].degree(k);
        for (size_t s = 0; s != splines.size(); ++s)
        {
            const short_t p = splines[s].degree(k);
            if (p < pmax) splines[s].degreeElevate(pmax - p, k);
        }

        // 2. Union of the interior knots, with the maximal multiplicity.
        const T first = splines[0].knots(k).first();
        const T last  = splines[0].knots(k).last();
        std::vector<T> vals;
        for (size_t s = 0; s != splines.size(); ++s)
        {
            const KnotVectorType & kv = splines[s].knots(k);
            GISMO_ASSERT(kv.first() == first && kv.last() == last,
                         "makeCompatible(): spline " << s << " lives on a "
                         "different parametric domain in direction " << k << ".");
            for (size_t i = 0; i < kv.uSize(); ++i)
            {
                const T xi = kv.uValue(i);
                if (xi > first && xi < last) vals.push_back(xi);
            }
        }
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());

        for (size_t i = 0; i != vals.size(); ++i)
        {
            const T xi = vals[i];
            index_t mmax = 0;
            for (size_t s = 0; s != splines.size(); ++s)
            {
                const KnotVectorType & kv = splines[s].knots(k);
                const index_t m = kv.has(xi) ? static_cast<index_t>(kv.multiplicity(xi)) : 0;
                if (m > mmax) mmax = m;
            }
            for (size_t s = 0; s != splines.size(); ++s)
            {
                const KnotVectorType & kv = splines[s].knots(k);
                const index_t m = kv.has(xi) ? static_cast<index_t>(kv.multiplicity(xi)) : 0;
                if (mmax > m) splines[s].insertKnot(xi, k, mmax - m);
            }
        }
    }
}

template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::linearCombination(
                                    T a, const gsTensorBSpline<d,T> & A,
                                    T b, const gsTensorBSpline<d,T> & B)
{
    GISMO_ASSERT(A.targetDim() == B.targetDim(),
                 "linearCombination(): target dimensions differ ("
                 << A.targetDim() << " vs " << B.targetDim() << ").");

    std::vector< gsTensorBSpline<d,T> > v;
    v.reserve(2);
    v.push_back(A);
    v.push_back(B);
    makeCompatible(v);

    GISMO_ASSERT(v[0].coefs().rows() == v[1].coefs().rows(),
                 "linearCombination(): coefficient size mismatch after "
                 "makeCompatible() (" << v[0].coefs().rows() << " vs "
                 << v[1].coefs().rows() << ").");

    v[0].coefs() = a * v[0].coefs() + b * v[1].coefs();
    return v[0];
}

// ----------------------------------------------------------------------------
//  div() — divergence (for vector-valued splines with targetDim == d)
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::div() const
{
    GISMO_ASSERT(this->targetDim() == static_cast<short_t>(d),
                 "div() requires targetDim == d (got targetDim="
                 << this->targetDim() << ", d=" << d << ").");

    // For each component k, extract column k of coefs → scalar spline,
    // differentiate in direction k, then sum all terms in a common space.
    //
    // After grad(k), term k has degree p-1 in direction k and p elsewhere.
    // Degree-elevate direction k by 1 so every term is already at degree
    // (p,...,p) before makeCompatible(); this also fixes the d==1 case, where
    // there is nothing for makeCompatible() to elevate against.
    std::vector< gsTensorBSpline<d,T> > terms;
    terms.reserve(d);
    for (short_t k = 0; k < d; ++k)
    {
        gsMatrix<T> col_k = this->m_coefs.col(k);
        gsTensorBSpline<d,T> comp_k(
            static_cast<const Basis&>(this->basis()), col_k);
        gsTensorBSpline<d,T> dk = comp_k.grad(k);
        dk.degreeElevate(1, k);
        terms.push_back(give(dk));
    }

    makeCompatible(terms);

    gsTensorBSpline<d,T> result = terms[0];
    for (size_t k = 1; k != terms.size(); ++k)
    {
        GISMO_ASSERT(result.m_coefs.rows() == terms[k].m_coefs.rows(),
                     "div(): coefficient size mismatch after knot equalization ("
                     << result.m_coefs.rows() << " vs " << terms[k].m_coefs.rows() << ")");
        result.m_coefs += terms[k].m_coefs;
    }
    return result;
}

// ----------------------------------------------------------------------------
//  lapl() — Laplacian for scalar splines (targetDim == 1)
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::lapl(bool keepBezier) const
{
    GISMO_ASSERT(this->targetDim() == 1,
                 "lapl() is only implemented for scalar (targetDim==1) splines.");
    for (short_t k = 0; k < d; ++k)
        GISMO_ASSERT(this->degree(k) >= 2,
                     "lapl(): degree must be >= 2 in all directions (degree("
                     << k << ")=" << this->degree(k) << ").");

    // D2_k = ∂²c/∂x_k² has degree p_k-2 and knot vector kv_k[2:-2] in direction
    // k, and is unchanged in the other directions.  Every term is elevated back
    // to degree p in direction k, then all terms are summed in the common space
    // produced by makeCompatible().
    //
    // keepBezier == false (default): that common space is S(p, C^{p-3}) —
    //   two differentiations in direction k raise the interior multiplicity
    //   there by 2, and the union propagates that to the other terms.
    //
    // keepBezier == true: each unique interior knot of the original basis is
    //   additionally inserted twice in every direction j != k, which drives the
    //   union to Bézier (C^{-1}) form.
    std::vector< gsTensorBSpline<d,T> > terms;
    terms.reserve(d);

    for (short_t k = 0; k < d; ++k)
    {
        gsTensorBSpline<d,T> D2 = this->grad(k).grad(k);
        D2.degreeElevate(2, k);

        if (keepBezier)
        {
            for (short_t j = 0; j < d; ++j)
            {
                if (j == k) continue;
                const KnotVectorType & orig_kv = this->knots(j);
                const T first = orig_kv.first();
                const T last  = orig_kv.last();
                for (size_t i = 0; i < orig_kv.uSize(); ++i)
                {
                    const T xi = orig_kv.uValue(i);
                    if (xi <= first || xi >= last) continue;
                    D2.insertKnot(xi, j, 2);
                }
            }
        }
        terms.push_back(give(D2));
    }

    makeCompatible(terms);

    gsTensorBSpline<d,T> result = terms[0];
    for (size_t k = 1; k != terms.size(); ++k)
    {
        GISMO_ASSERT(result.m_coefs.rows() == terms[k].m_coefs.rows(),
                     "lapl(): coefficient size mismatch between directions ("
                     << result.m_coefs.rows() << " vs " << terms[k].m_coefs.rows() << ")");
        result.m_coefs += terms[k].m_coefs;
    }
    return result;
}

// ----------------------------------------------------------------------------
//  hess() — Hessian of a scalar spline, stored as a targetDim == d*d spline
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::hess() const
{
    GISMO_ASSERT(this->targetDim() == 1,
                 "hess() is only implemented for scalar (targetDim==1) splines.");
    for (short_t k = 0; k < d; ++k)
        GISMO_ASSERT(this->degree(k) >= 2,
                     "hess(): degree must be >= 2 in all directions (degree("
                     << k << ")=" << this->degree(k) << ").");

    // H_ij = ∂²c/∂x_i∂x_j lives in a different space for each (i,j): degree
    // p-2 in direction i when i==j, p-1 in directions i and j otherwise, with
    // the interior multiplicities of the input unchanged.  makeCompatible()
    // elevates every entry back to degree p — which raises the multiplicity in
    // a differentiated direction by the number of differentiations there — and
    // then takes the union, so the common space is degree p with interior
    // multiplicity m+2, i.e. two orders less continuous than the input.
    // That is exactly the space lapl() sums in, as it must be: the trace of
    // the Hessian is the Laplacian.
    std::vector< gsTensorBSpline<d,T> > terms;
    terms.reserve(d*d);
    for (short_t i = 0; i < d; ++i)
        for (short_t j = 0; j < d; ++j)
            terms.push_back(this->grad(i).grad(j));

    makeCompatible(terms);

    const index_t n = terms[0].coefs().rows();
    gsMatrix<T> coefs(n, d*d);
    for (size_t e = 0; e != terms.size(); ++e)
    {
        GISMO_ASSERT(terms[e].coefs().rows() == n,
                     "hess(): coefficient size mismatch after makeCompatible() ("
                     << terms[e].coefs().rows() << " vs " << n << ").");
        coefs.col(e) = terms[e].coefs().col(0);
    }

    return gsTensorBSpline<d,T>(terms[0].basis(), give(coefs));
}

namespace internal
{

/// @brief Get a Tensor BSpline from XML data
///
/// \ingroup Nurbs
template<short_t d, class T>
class gsXml< gsTensorBSpline<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorBSpline<TMPLA2(d,T)>);
    static std::string tag ()  { return "Geometry"; }
    static std::string type () { return "TensorBSpline" +  to_string(d); }

    static gsTensorBSpline<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTensorBSpline<d,T> >( node );
    }

    static gsXmlNode * put (const gsTensorBSpline<d,T> & obj,
                            gsXmlTree & data)
    {
        return putGeometryToXml(obj,data);
    }
    GSXML_GET_INTO(gsTensorBSpline<TMPLA2(d,T)>);
};



}// namespace internal

} // namespace gismo
