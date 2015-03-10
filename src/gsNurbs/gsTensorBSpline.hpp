
#pragma once 

#include <gsIO/gsXmlUtils.h>

#include <gsNurbs/gsBSplineBasis.h>

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
void gsTensorBSpline<d,T,KnotVectorType>::slice(index_t dir_fixed,T par,gsTensorBSpline<d-1,T>& result) const
{
    typedef gsTensorBSplineBasis<d-1,T,KnotVectorType> newBasis;
    GISMO_ASSERT(d-1>0,"cannot take iso slice of a curve");
    // construct the d-1 basis
    std::vector<gsBSplineBasis<T,KnotVectorType>* > bases;
    for(unsigned i=0;i<d;++i)
        if( static_cast<unsigned>(i) != dir_fixed)
            bases.push_back(new gsBSplineBasis<T,KnotVectorType>(this->basis().knots(i)));
    newBasis *tbasis = new newBasis(bases) ;

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
        gsTensorBSpline<d,T>* clone = this->clone();

        gsVector<unsigned> strides;
        gsVector<int> intStrides;
        this->basis().stride_cwise(intStrides);
        strides=intStrides.cast<unsigned>();
        gsTensorBoehm<T,KnotVectorType,gsMatrix<T> >(
                    clone->basis().knots(dir_fixed),clone->coefs(),par,dir_fixed,
                    strides,degree-mult,true);

        // extract right ceofficients
        constructCoefsForSlice(dir_fixed,par,*clone,coefs);
        delete clone;
    }

    // construct the object
    result = gsTensorBSpline<d-1,T>(*tbasis, give(coefs) );
    delete tbasis;
}

template<unsigned d, class T, class KnotVectorType>
void gsTensorBSpline<d,T,KnotVectorType>::constructCoefsForSlice(index_t dir_fixed,T par,const gsTensorBSpline<d,T>& geo,gsMatrix<T>& result) const
{
    const gsTensorBSplineBasis<d,T,KnotVectorType>& base = geo.basis();
    const gsMatrix<T>& fullCoefs=geo.coefs();
    // pick the right coefficients and store them in coefs
    const unsigned degree = base.degree(dir_fixed);
    const KnotVectorType& knots = base.knots(dir_fixed);
    const int index = knots.findspan(par)-degree;
    gsVector<index_t,d> sizes,lowerCorner,upperCorner;

    gsVector<unsigned,d> ssizes;
    base.size_cwise( ssizes );
    sizes=ssizes.template cast<index_t>();
    lowerCorner.setZero();
    lowerCorner(dir_fixed)=index;
    upperCorner=sizes;
    upperCorner(dir_fixed)=index+1;

    // to do: gsMatrix<index_t> ind = gsTensorBasis::coefSlice(dim_fixed, index) ?

    gsTensorGridIterator<index_t> gridIter(sizes);
    gsTensorGridIterator<index_t> *iter = gridIter.makeSubGridIterator(lowerCorner,upperCorner);
    index_t size=1;
    for(index_t i = 0;i<d;++i)
        if(dir_fixed!=i)
            size*=sizes(i);
    result.resize(size,fullCoefs.cols());
    index_t i=0;
    for(iter->first();iter->good();iter->next())
    {
        result.row(i)=fullCoefs.row(iter->flatIndex());
        ++i;
    }
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
