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

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBoehm.h>

#include <gsUtils/gsMultiIndexIterators.h>

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

#include <gsTensor/gsTensorTools.h>

#include <gsCore/gsConstantFunction.h>

namespace gismo
{

template<unsigned d, class T, class KnotVectorType>
gsTensorBSpline<d,T,KnotVectorType>::gsTensorBSpline(gsMatrix<T> const & corner, KnotVectorType const& KV1, KnotVectorType const & KV2)
{
  assert(d==2);

  gsBSplineBasis<T,KnotVectorType> * Bu= new gsBSplineBasis<T,KnotVectorType>(KV1);
  gsBSplineBasis<T,KnotVectorType> * Bv= new gsBSplineBasis<T,KnotVectorType>(KV2);
  Basis *tbasis = new Basis(Bu,Bv) ;//d==2

  int n1 = KV1.size() - KV1.degree() - 1;
  int n2 = KV2.size() - KV2.degree() - 1;

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
template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::slice(index_t dir_fixed,T par,
                                                BoundaryGeometryType & result) const
{
    GISMO_ASSERT(d-1>=0,"d must be greater or equal than 1");
    GISMO_ASSERT(dir_fixed>=0 && static_cast<unsigned>(dir_fixed)<d,"cannot fix a dir greater than dim or smaller than 0");
    // construct the d-1 basis
    boxSide side(dir_fixed,0);
    BoundaryBasisType *tbasis = this->basis().boundaryBasis(side) ;

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

        gsMatrix<T> coefs;
        if( mult>=degree )
        {
            // no knot insertion needed, just extract the right coefficients
            constructCoefsForSlice(dir_fixed,par,*this,coefs);
        }
        else
        {
            // clone the basis and inserting upto degree knots at par
            gsTensorBSpline<d,T,KnotVectorType>* clone = this->clone();

            gsVector<index_t,d> intStrides;
            this->basis().stride_cwise(intStrides);
            gsTensorBoehm<T,KnotVectorType,gsMatrix<T> >(
                        clone->basis().knots(dir_fixed),clone->coefs(),par,dir_fixed,
                        intStrides.template cast<unsigned>(), degree-mult,true);

            // extract right ceofficients
            constructCoefsForSlice(dir_fixed,par,*clone,coefs);
            delete clone;
        }

        // construct the object
        //result = gsTensorBSpline<d-1,T>(*tbasis, give(coefs) );
        //result = BoundaryGeometry(*tbasis, give(coefs) );
        result = BoundaryGeometryType(*tbasis, coefs );
    }
    delete tbasis;
}

template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::reverse(unsigned k)
{ 
    gsTensorBSplineBasis<d,T,KnotVectorType> & tbsbasis = this->basis();
    gsVector<int,d> sz;
    tbsbasis.size_cwise(sz);
    flipTensorVector(k, sz, m_coefs);
    tbsbasis.component(k).reverse();
}


template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::
swapDirections(const unsigned i, const unsigned j)
{
    gsVector<int,d> sz;
    this->basis().size_cwise(sz);
    swapTensorDirection(i, j, sz, m_coefs);
    this->basis().swapDirections(i,j);
}

template<unsigned d, class T, class KnotVectorType>
bool gsTensorBSpline<d,T,KnotVectorType>::isPatchCorner(gsMatrix<T> const &v, T tol) const
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

template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::findCorner(const gsMatrix<T> & v, 
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

template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::setOriginCorner(gsMatrix<T> const &v)
{
    gsVector<index_t,d> curr;
    findCorner(v, curr);
    if ( curr[0] == this->basis().size(0) )
        return;
    for(unsigned k = 0; k!=d; ++k)
        if ( curr[k] != 0 )
            this->reverse(k);
}

template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::setFurthestCorner(gsMatrix<T> const &v)
{
    gsVector<index_t,d> curr;
    findCorner(v, curr);
    if ( curr[0] == this->basis().size(0) )
        return;
    for(unsigned k = 0; k!=d; ++k)
        if ( curr[k] == 0 )
            this->reverse(k);
}


template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::degreeElevate(int const i, int const dir)
{
    if (dir == -1)
    {
        for (unsigned j = 0; j < d; ++j)
            degreeElevate(i, j);
        return;
    }
    
    GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
                  "Invalid basis component "<< dir <<" requested for degree elevation" );

    const index_t n = this->m_coefs.cols();
    
    gsVector<index_t,d> sz;
    this->basis().size_cwise(sz);
    
    swapTensorDirection(0, dir, sz, this->m_coefs);
    this->m_coefs.resize( sz[0], n * sz.template tail<d-1>().prod() );

    bspline::degreeElevateBSpline(this->basis().component(dir), this->m_coefs, i);
    sz[0] = this->m_coefs.rows();

    this->m_coefs.resize( sz.prod(), n );
    swapTensorDirection(0, dir, sz, this->m_coefs);
}


template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::constructCoefsForSlice(unsigned dir_fixed,T par,const gsTensorBSpline<d,T,KnotVectorType>& geo,gsMatrix<T>& result) const
{
    const gsTensorBSplineBasis<d,T,KnotVectorType>& base = geo.basis();
    const gsMatrix<T>& fullCoefs=geo.coefs();
    // pick the right coefficients and store them in coefs
    const unsigned degree = base.degree(dir_fixed);
    const KnotVectorType& knots = base.knots(dir_fixed);
    const int index = knots.findspan(par)-degree;
    gsVector<index_t,d> sizes,lowerCorner,upperCorner;
    base.size_cwise( sizes );
    lowerCorner.setZero();
    lowerCorner(dir_fixed)=index;
    upperCorner=sizes;
    upperCorner(dir_fixed)=index+1;

    // to do: gsMatrix<index_t> ind = gsTensorBasis::coefSlice(dim_fixed, index) ?

    gsTensorGridIterator<index_t> gridIter(sizes);
    gsTensorGridIterator<index_t> * iter = gridIter.makeSubGridIterator(lowerCorner,upperCorner);
    index_t size=1;
    for(unsigned i = 0;i<d;++i)
        if(dir_fixed!=i)
            size*=sizes(i);
    result.resize(size,fullCoefs.cols());
    index_t i=0;
    for(iter->first();iter->good();iter->next())
    {
        result.row(i)=fullCoefs.row(iter->flatIndex());
        ++i;
    }

    delete iter;
}


template<unsigned d, class T, class KnotVectorType>
std::ostream & gsTensorBSpline<d,T,KnotVectorType>::print(std::ostream &os) const
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


namespace internal
{

/// @brief Get a Tensor BSpline from XML data
///
/// \ingroup Nurbs
template<unsigned d, class T, class KnotVectorType>
class gsXml< gsTensorBSpline<d,T, KnotVectorType> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorBSpline<TMPLA3(d,T,KnotVectorType)>);
    static std::string tag ()  { return "Geometry"; }
    static std::string type () { return "TensorBSpline" +  to_string(d); }

    static gsTensorBSpline<d,T,KnotVectorType> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTensorBSpline<d,T,KnotVectorType> >( node );
    }
    
    static gsXmlNode * put (const gsTensorBSpline<d,T,KnotVectorType> & obj,
                            gsXmlTree & data)
    {
        return putGeometryToXml(obj,data);
    }
};



}// namespace internal

} // namespace gismo
