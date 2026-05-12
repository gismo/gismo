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
    const int n1 = KV1.size() - KV1.degree() - 1;
    const int n2 = KV2.size() - KV2.degree() - 1;

    cbases.push_back(new gsBSplineBasis<T>(give(KV1)) );
    cbases.push_back(new gsBSplineBasis<T>(give(KV2)) );
    Basis * tbasis = Basis::New(cbases); //d==2


    GISMO_ASSERT( (corner.rows()==4) && (corner.cols()==3),
                  "gsTensorBSpline: Please make sure that the size of *corner* is 4-by-3");

    gsMatrix<T> pcp (n1*n2, 3);
    // set up CPs on boundary first. The inner CPs on each boundary curve are
    // uniformly linear dependent on the two corner CPs
    int j=0; // boundary v=0
    for (int i=0; i<=n1-1; i++)
    {
        for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(0,xi) + i/((T)(n1-1))*( corner(1,xi) - corner(0,xi) );
        }
    }
    j=n2-1; // boundary v=1
    for (int i=0; i<=n1-1; i++)
    {
        for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(3,xi) + i/((T)(n1-1))*( corner(2,xi) - corner(3,xi) );
        }
    }
    int i=0; // boundary u=0;
    for (j=0; j<=n2-1; j++)
    {
        for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(0,xi) + j/((T)(n2-1))*( corner(3,xi) - corner(0,xi) );
        }
    }
    i=n1-1; // boundary u=1;
    for (j=0; j<=n2-1; j++)
    {
        for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
        {
            pcp(i+j*n1,xi)=corner(1,xi) + j/((T)(n2-1))*( corner(2,xi) - corner(1,xi) );
        }
    }
    // uniformly linear dependent in horizontal direction
    for (j=1; j<=n2-2; j++)
    {
        for (i=1; i<=n1-2; i++)
        {
            for (unsigned int xi=0; xi<=2; xi++) //specification of x or y or z
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
        const int mult   = this->basis().knots(dir_fixed).multiplicity(par);
        const int degree = this->basis().degree(dir_fixed);
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
                intStrides.template cast<unsigned>(), degree-mult,true);

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
index_t gsTensorBSpline<d,T>::removeKnot( T knot, short_t dir, short_t i)
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
        true);
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
        for (index_t i = 0; i < kv0.uSize(); ++i)
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
//  squared  (Bernstein product formula)
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
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::squared(bool keepBezier) const
{
    GISMO_ASSERT(this->targetDim() == 1,
                 "squared() is only implemented for scalar (targetDim==1) splines.");

    // 1. Convert to Bézier form (C^{-1})
    gsTensorBSpline<d,T> bez = this->toBezier();

    // degrees in the Bézier form
    gsVector<index_t,d> deg;
    for (short_t k = 0; k < d; ++k)
        deg[k] = static_cast<index_t>(bez.degree(k));

    // 2. Build the result knot vectors: each unique Bézier knot gets
    //    multiplicity 2*p_k+1, giving degree 2*p_k with C^{-1} continuity.
    std::vector<KnotVectorType> kvs_result(d);
    for (short_t k = 0; k < d; ++k)
    {
        const KnotVectorType & kv = bez.knots(k);
        const index_t p  = deg[k];
        const index_t p2 = 2 * p;
        std::vector<T> knots_result;
        for (index_t i = 0; i < kv.uSize(); ++i)
        {
            const T xi = kv.uValue(i);
            for (index_t m = 0; m < p2 + 1; ++m)
                knots_result.push_back(xi);
        }
        kvs_result[k] = KnotVectorType(knots_result, static_cast<short_t>(p2));
    }

    // 3. Assemble result basis
    std::vector<Family_t*> cbases;
    cbases.reserve(d);
    for (short_t k = 0; k < d; ++k)
        cbases.push_back(new gsBSplineBasis<T>(kvs_result[k]));
    Basis * result_basis = Basis::New(cbases);

    // 4. Compute result coefficients via the Bernstein product formula.
    //
    //    In direction k, c_i^2 has coefficient index running from 0 to 2p_k.
    //    The formula for c^2 in the Bézier basis (using the product rule for
    //    Bernstein polynomials) is:
    //
    //      q_{k_0,...,k_{d-1}} = sum_{i+j=k (per dim)}
    //          prod_dim [ C(p_dim, i_dim) * C(p_dim, j_dim) / C(2p_dim, k_dim) ]
    //          * b_{i_0,...} * b_{j_0,...}
    //
    //    We loop over elements (Bézier patches), extract the (p+1)^d local
    //    coefficients, apply the formula, and write the (2p+1)^d output coefs.

    // sizes of Bézier grid per direction
    gsVector<index_t,d> bez_sz;
    bez.basis().size_cwise(bez_sz);

    // number of elements per direction
    gsVector<index_t,d> n_elems;
    for (short_t k = 0; k < d; ++k)
        n_elems[k] = bez_sz[k] / (deg[k] + 1); // = number of Bézier patches

    // strides for Bézier coef flat index
    gsVector<index_t,d> bez_stride;
    bez_stride[0] = 1;
    for (short_t k = 1; k < d; ++k)
        bez_stride[k] = bez_stride[k-1] * bez_sz[k-1];

    // result sizes and strides
    gsVector<index_t,d> res_sz;
    for (short_t k = 0; k < d; ++k)
        res_sz[k] = n_elems[k] * (2*deg[k] + 1);

    gsVector<index_t,d> res_stride;
    res_stride[0] = 1;
    for (short_t k = 1; k < d; ++k)
        res_stride[k] = res_stride[k-1] * res_sz[k-1];

    gsMatrix<T> res_coefs(res_sz.prod(), 1);
    res_coefs.setZero();

    const gsMatrix<T> & bc = bez.coefs();  // Bézier coefs (n_bez x 1)

    // Iterate over all elements (multi-index of patches)
    gsVector<index_t,d> elem_idx;
    elem_idx.setZero();
    do
    {
        // For each element, iterate over result Bernstein index k ∈ [0,2p]^d
        gsVector<index_t,d> k_idx;
        k_idx.setZero();
        do
        {
            T val = T(0);
            // Accumulate: sum over i where j = k - i, 0 <= i,j <= p (per dim)
            // Loop over i_idx in [0,p]^d
            gsVector<index_t,d> i_idx;
            i_idx.setZero();
            // We need j_idx = k_idx - i_idx valid (all entries in [0,p])
            do
            {
                bool valid = true;
                gsVector<index_t,d> j_idx;
                T weight = T(1);
                for (short_t kd = 0; kd < d; ++kd)
                {
                    j_idx[kd] = k_idx[kd] - i_idx[kd];
                    if (j_idx[kd] < 0 || j_idx[kd] > deg[kd])
                    { valid = false; break; }
                    // Bernstein product weight: C(p,i)*C(p,j)/C(2p,k)
                    weight *= static_cast<T>(binom(deg[kd], i_idx[kd]))
                            * static_cast<T>(binom(deg[kd], j_idx[kd]))
                            / static_cast<T>(binom(2*deg[kd], k_idx[kd]));
                }
                if (valid)
                {
                    // Flat index in Bézier coef array for i and j
                    index_t fi = 0, fj = 0;
                    for (short_t kd = 0; kd < d; ++kd)
                    {
                        fi += (elem_idx[kd] * (deg[kd]+1) + i_idx[kd]) * bez_stride[kd];
                        fj += (elem_idx[kd] * (deg[kd]+1) + j_idx[kd]) * bez_stride[kd];
                    }
                    val += weight * bc(fi, 0) * bc(fj, 0);
                }

                // Increment i_idx
                short_t carry = 0;
                do
                {
                    ++i_idx[carry];
                    if (i_idx[carry] > deg[carry]) { i_idx[carry] = 0; ++carry; }
                    else break;
                } while (carry < d);
                if (carry == d) break;
            } while (true);

            // Flat index in result coef array
            index_t fr = 0;
            for (short_t kd = 0; kd < d; ++kd)
                fr += (elem_idx[kd] * (2*deg[kd]+1) + k_idx[kd]) * res_stride[kd];
            res_coefs(fr, 0) = val;

            // Increment k_idx
            short_t carry = 0;
            do
            {
                ++k_idx[carry];
                if (k_idx[carry] > 2*deg[carry]) { k_idx[carry] = 0; ++carry; }
                else break;
            } while (carry < d);
            if (carry == d) break;
        } while (true);

        // Increment elem_idx
        short_t carry = 0;
        do
        {
            ++elem_idx[carry];
            if (elem_idx[carry] >= n_elems[carry]) { elem_idx[carry] = 0; ++carry; }
            else break;
        } while (carry < d);
        if (carry == d) break;
     } while (true);

    gsTensorBSpline<d,T> result(*result_basis, give(res_coefs));

    if (!keepBezier)
    {
        // Remove interior knots to reach minimal space S(2p, C^{p-1}).
        // In direction k: degree is 2*p_k, target continuity C^{p_k-1}
        //   ⟹ interior multiplicity = 2*p_k - (p_k - 1) = p_k + 1
        for (short_t k = 0; k < d; ++k)
        {
            const int targetMult = static_cast<int>(deg[k]) + 1; // p+1
            std::vector<T> interior;
            {
                const KnotVectorType & kv0 = result.knots(k);
                const T first = kv0.first(), last = kv0.last();
                for (index_t i = 0; i < kv0.uSize(); ++i)
                {
                    const T xi = kv0.uValue(i);
                    if (xi > first && xi < last) interior.push_back(xi);
                }
            }
            for (const T xi : interior)
            {
                int mult = result.knots(k).multiplicity(xi);
                while (mult > targetMult)
                {
                    if (result.removeKnot(xi, k, 1) == 0) break;
                    --mult;
                }
            }
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
//  cubed  (Bernstein triple-product via squared × self)
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::cubed(bool keepBezier) const
{
    GISMO_ASSERT(this->targetDim() == 1,
                 "cubed() is only implemented for scalar (targetDim==1) splines.");

    // c³ = c² × c.  Compute c² in full Bézier (we need to multiply again)
    // and self in Bézier, then do the Bernstein product.
    gsTensorBSpline<d,T> c2 = this->squared(true);  // Bézier of c²
    gsTensorBSpline<d,T> c1 = this->toBezier();     // Bézier of c

    // degrees of c² (= 2*p_k) and c (= p_k)
    gsVector<index_t,d> deg2, deg1;
    for (short_t k = 0; k < d; ++k)
    {
        deg2[k] = static_cast<index_t>(c2.degree(k));  // = 2*p_k
        deg1[k] = static_cast<index_t>(c1.degree(k));  // = p_k
    }
    // Original degrees
    gsVector<index_t,d> deg;
    for (short_t k = 0; k < d; ++k)
        deg[k] = deg1[k];

    // Build result knot vectors: degree 3*p_k, Bézier (mult 3*p_k+1)
    std::vector<KnotVectorType> kvs_result(d);
    for (short_t k = 0; k < d; ++k)
    {
        const KnotVectorType & kv = c1.knots(k);
        const index_t p3 = 3 * deg[k];
        std::vector<T> knots_result;
        for (index_t i = 0; i < kv.uSize(); ++i)
        {
            const T xi = kv.uValue(i);
            for (index_t m = 0; m < p3 + 1; ++m)
                knots_result.push_back(xi);
        }
        kvs_result[k] = KnotVectorType(knots_result, static_cast<short_t>(p3));
    }
    std::vector<Family_t*> cbases3;
    cbases3.reserve(d);
    for (short_t k = 0; k < d; ++k)
        cbases3.push_back(new gsBSplineBasis<T>(kvs_result[k]));
    Basis * result_basis3 = Basis::New(cbases3);

    // sizes
    gsVector<index_t,d> c2_sz, c1_sz;
    c2.basis().size_cwise(c2_sz);
    c1.basis().size_cwise(c1_sz);

    // number of elements per direction (same for c2 and c1 since same knot structure)
    gsVector<index_t,d> n_elems;
    for (short_t k = 0; k < d; ++k)
        n_elems[k] = c1_sz[k] / (deg1[k] + 1);

    gsVector<index_t,d> c2_stride, c1_stride;
    c2_stride[0] = 1; c1_stride[0] = 1;
    for (short_t k = 1; k < d; ++k)
    {
        c2_stride[k] = c2_stride[k-1] * c2_sz[k-1];
        c1_stride[k] = c1_stride[k-1] * c1_sz[k-1];
    }

    gsVector<index_t,d> res_sz;
    for (short_t k = 0; k < d; ++k)
        res_sz[k] = n_elems[k] * (3*deg[k] + 1);
    gsVector<index_t,d> res_stride;
    res_stride[0] = 1;
    for (short_t k = 1; k < d; ++k)
        res_stride[k] = res_stride[k-1] * res_sz[k-1];

    gsMatrix<T> res_coefs3(res_sz.prod(), 1);
    res_coefs3.setZero();

    const gsMatrix<T> & bc2 = c2.coefs();
    const gsMatrix<T> & bc1 = c1.coefs();

    // Bernstein product: c³_{k} = sum_{i+j=k} C(2p,i)*C(p,j)/C(3p,k) * c²_i * c_j
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
                    if (j_idx[kd] < 0 || j_idx[kd] > deg1[kd])
                    { valid = false; break; }
                    weight *= static_cast<T>(binom(deg2[kd], i_idx[kd]))
                            * static_cast<T>(binom(deg1[kd], j_idx[kd]))
                            / static_cast<T>(binom(deg2[kd]+deg1[kd], k_idx[kd]));
                }
                if (valid)
                {
                    index_t fi = 0, fj = 0;
                    for (short_t kd = 0; kd < d; ++kd)
                    {
                        fi += (elem_idx[kd] * (deg2[kd]+1) + i_idx[kd]) * c2_stride[kd];
                        fj += (elem_idx[kd] * (deg1[kd]+1) + j_idx[kd]) * c1_stride[kd];
                    }
                    val += weight * bc2(fi, 0) * bc1(fj, 0);
                }
                short_t carry = 0;
                do
                {
                    ++i_idx[carry];
                    if (i_idx[carry] > deg2[carry]) { i_idx[carry] = 0; ++carry; }
                    else break;
                } while (carry < d);
                if (carry == d) break;
            } while (true);

            index_t fr = 0;
            for (short_t kd = 0; kd < d; ++kd)
                fr += (elem_idx[kd] * (3*deg[kd]+1) + k_idx[kd]) * res_stride[kd];
            res_coefs3(fr, 0) = val;

            short_t carry = 0;
            do
            {
                ++k_idx[carry];
                if (k_idx[carry] > 3*deg[carry]) { k_idx[carry] = 0; ++carry; }
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

    gsTensorBSpline<d,T> result3(*result_basis3, give(res_coefs3));

    if (!keepBezier)
    {
        // Remove interior knots to reach minimal space S(3p, C^{p-1}).
        // Degree is 3*p_k, continuity C^{p_k-1}
        //   ⟹ interior multiplicity = 3*p_k - (p_k - 1) = 2*p_k + 1
        for (short_t k = 0; k < d; ++k)
        {
            const int targetMult = 2 * static_cast<int>(deg[k]) + 1;
            std::vector<T> interior;
            {
                const KnotVectorType & kv0 = result3.knots(k);
                const T first = kv0.first(), last = kv0.last();
                for (index_t i = 0; i < kv0.uSize(); ++i)
                {
                    const T xi = kv0.uValue(i);
                    if (xi > first && xi < last) interior.push_back(xi);
                }
            }
            for (const T xi : interior)
            {
                int mult = result3.knots(k).multiplicity(xi);
                while (mult > targetMult)
                {
                    if (result3.removeKnot(xi, k, 1) == 0) break;
                    --mult;
                }
            }
        }
    }
    return result3;
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
        GISMO_ASSERT(math::abs(denom) > T(0),
                     "grad: zero denominator at i=" << i <<
                     " (kv[" << i+p+1 << "]=" << kv[i+p+1] <<
                     ", kv[" << i+1 << "]=" << kv[i+1] << ")");
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

    // Compute all d terms first.
    // After grad(k), term k has degree p-1 in direction k and p elsewhere.
    // Degree-elevate direction k by 1 so every term lives in degree (p,...,p).
    std::vector< gsTensorBSpline<d,T> > terms;
    terms.reserve(d);
    for (short_t k = 0; k < d; ++k)
    {
        // Extract scalar component k
        gsMatrix<T> col_k = this->m_coefs.col(k);
        gsTensorBSpline<d,T> comp_k(
            static_cast<const Basis&>(this->basis()), col_k);
        gsTensorBSpline<d,T> dk = comp_k.grad(k);
        dk.degreeElevate(1, k);   // lift dir k back to degree p
        terms.push_back(dk);
    }

    // Bring all terms to a common knot vector by knot insertion.
    // The common knot vector per direction is the union of all terms' knot vectors.
    // We do this by inserting missing knots from every term into a running result.
    gsTensorBSpline<d,T> result = terms[0];
    for (short_t k = 1; k < static_cast<short_t>(d); ++k)
    {
        // Insert into result any knots present in terms[k] but missing/lower-mult in result
        for (short_t dir = 0; dir < d; ++dir)
        {
            const KnotVectorType & kv_src = terms[k].knots(dir);
            const T first = kv_src.first();
            const T last  = kv_src.last();
            for (index_t i = 0; i < kv_src.uSize(); ++i)
            {
                const T xi = kv_src.uValue(i);
                if (xi <= first || xi >= last) continue;
                const int m_src = kv_src.multiplicity(xi);
                const int m_res = result.knots(dir).has(xi)
                    ? result.knots(dir).multiplicity(xi) : 0;
                if (m_src > m_res)
                    result.insertKnot(xi, dir, m_src - m_res);
            }
        }
        // Also bring terms[k] up to result's knot structure
        for (short_t dir = 0; dir < d; ++dir)
        {
            const KnotVectorType & kv_res = result.knots(dir);
            const T first = kv_res.first();
            const T last  = kv_res.last();
            for (index_t i = 0; i < kv_res.uSize(); ++i)
            {
                const T xi = kv_res.uValue(i);
                if (xi <= first || xi >= last) continue;
                const int m_res = kv_res.multiplicity(xi);
                const int m_src = terms[k].knots(dir).has(xi)
                    ? terms[k].knots(dir).multiplicity(xi) : 0;
                if (m_res > m_src)
                    terms[k].insertKnot(xi, dir, m_res - m_src);
            }
        }
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

    // We compute D2_k = ∂²c/∂x_k² for each direction k, then bring all terms
    // to the same B-spline space and sum coefficients.
    //
    // After two differentiations in direction k, D2_k has:
    //   - degree p_k - 2 in direction k
    //   - degree p_j in all other directions j
    //   - knot vector kv_k[2:-2] in direction k (drop 2 from each end)
    //   - knot vector kv_j unchanged in other directions
    //
    // keepBezier == false (default, minimal space S(p, C^{p-3})):
    //   Bring all D2_k terms to the same knot structure by knot union only —
    //   no degree elevation.  Direction k has degree p-2, continuity C^{p-3};
    //   other directions j stay at degree p.  Degree-elevate direction k from
    //   p-2 to p-2... wait, we want the *sum* in a common space.
    //
    //   The natural common space for all D2_k is:
    //     - degree p in every direction
    //     - interior multiplicity = orig_mult + 2 in the differentiated direction
    //       (because two differentiations reduce continuity by 2: mult increases by 2)
    //     - other directions: use knot *union* without extra insertions.
    //   We achieve this by:
    //     - For D2_k: degree-elevate direction k by 2 (p-2 → p).
    //     - For D2_k in other directions j: insert each unique interior knot of
    //       *only* the D2_j term's direction-j knot vector (which has +2 mult)
    //       into direction j of D2_k, bringing it up to the union of all terms.
    //   This gives the minimal space S(p, C^{p-3}).
    //
    // keepBezier == true (old over-elevated behaviour):
    //   Additionally insert each unique interior knot of the original kv_j twice
    //   in every other direction j ≠ k, then the result is a Bézier (C^{-1}).

    gsTensorBSpline<d,T> result;
    bool init = false;

    for (short_t k = 0; k < d; ++k)
    {
        // ∂²c/∂x_k²
        gsTensorBSpline<d,T> D2 = this->grad(k).grad(k);

        // Degree-elevate direction k from p_k-2 back to p_k
        D2.degreeElevate(2, k);

        if (keepBezier)
        {
            // Old behaviour: insert each unique interior knot of the
            // original basis twice in every other direction j.
            for (short_t j = 0; j < d; ++j)
            {
                if (j == k) continue;
                const KnotVectorType & orig_kv = this->knots(j);
                const T first = orig_kv.first();
                const T last  = orig_kv.last();
                for (index_t i = 0; i < orig_kv.uSize(); ++i)
                {
                    const T xi = orig_kv.uValue(i);
                    if (xi <= first || xi >= last) continue;
                    D2.insertKnot(xi, j, 2);
                }
            }
        }
        // (When keepBezier==false we do NOT do those extra insertions.
        //  The knot-union step below handles alignment between D2_k terms.)

        if (!init)
        {
            result = give(D2);
            init   = true;
        }
        else
        {
            // Bring result and D2 to the same knot structure (union per direction).
            for (short_t dir = 0; dir < d; ++dir)
            {
                // Insert into result knots present in D2 but missing/lower-mult
                {
                    const KnotVectorType & kv_src = D2.knots(dir);
                    const T first = kv_src.first(), last = kv_src.last();
                    for (index_t i = 0; i < kv_src.uSize(); ++i)
                    {
                        const T xi = kv_src.uValue(i);
                        if (xi <= first || xi >= last) continue;
                        const int ms = kv_src.multiplicity(xi);
                        const int mr = result.knots(dir).has(xi)
                                     ? result.knots(dir).multiplicity(xi) : 0;
                        if (ms > mr) result.insertKnot(xi, dir, ms - mr);
                    }
                }
                // Insert into D2 knots present in result but missing/lower-mult
                {
                    const KnotVectorType & kv_res = result.knots(dir);
                    const T first = kv_res.first(), last = kv_res.last();
                    for (index_t i = 0; i < kv_res.uSize(); ++i)
                    {
                        const T xi = kv_res.uValue(i);
                        if (xi <= first || xi >= last) continue;
                        const int mr = kv_res.multiplicity(xi);
                        const int ms = D2.knots(dir).has(xi)
                                     ? D2.knots(dir).multiplicity(xi) : 0;
                        if (mr > ms) D2.insertKnot(xi, dir, mr - ms);
                    }
                }
            }
            GISMO_ASSERT(result.m_coefs.rows() == D2.m_coefs.rows(),
                         "lapl(): coefficient size mismatch between directions ("
                         << result.m_coefs.rows() << " vs " << D2.m_coefs.rows() << ")");
            result.m_coefs += D2.m_coefs;
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
//  hess() — placeholder
// ----------------------------------------------------------------------------
template<short_t d, class T>
gsTensorBSpline<d,T> gsTensorBSpline<d,T>::hess() const
{
    GISMO_NO_IMPLEMENTATION;
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
