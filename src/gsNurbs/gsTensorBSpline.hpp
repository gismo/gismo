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
typename gsGeometry<T>::uPtr gsTensorBSpline<d,T>::localRep(const gsMatrix<T> & u) const
{
    std::vector<KnotVectorType> kv(d); // the local knot-vectors
    gsVector<index_t,d> cfirst, clast; // tensor-indices of local coefficients
    index_t sz = 1; // number of control points in the local representation

    // Fill in the data defined above
    for(unsigned i = 0; i!=d; ++i)
    {
        const int deg = degree(i);
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
        return result;
    }

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
};



}// namespace internal

} // namespace gismo
