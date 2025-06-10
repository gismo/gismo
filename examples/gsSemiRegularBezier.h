/** @file gsSemiRegularBezier.h

    @brief A semiregular Bezier basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, E. Skënderaj
*/

#pragma once

namespace gismo
{

template<short_t d, class T = real_t>
class gsSemiRegularBezier
{ };

template<class T>
class gsSemiRegularBezier<2,T> : public gsBasis<T>
{
    typedef gsBSplineBasis<T> bsbasis;
    
    bsbasis m_u;
    std::vector<bsbasis> m_v;
    std::vector<index_t> m_offset;

public:
    typedef memory::shared_ptr< gsSemiRegularBezier > Ptr;
    typedef memory::unique_ptr< gsSemiRegularBezier > uPtr;
    typedef T Scalar_t;
    static const bool IsRational = false;
    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

    gsSemiRegularBezier(int deg = 1);

public:
    // these are virtual functions that we must implement
    gsSemiRegularBezier * clone_impl() const
    { return new gsSemiRegularBezier(*this); }

    std::ostream &print(std::ostream &os) const;
    short_t domainDim() const { return 2; }
    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    {
        gsGenericGeometry<2,T> * g = new gsGenericGeometry<2,T>(*this, give(coefs) );
        return memory::unique_ptr<gsGeometry<T> >(g);
    }

public:
    index_t size() const;
    virtual void anchors_into(gsMatrix<T>& result) const;
    virtual void anchor_into(index_t i, gsMatrix<T>& result) const;
    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;
    virtual void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;
    virtual void numActive_into(const gsMatrix<T> & u, gsVector<index_t>& result) const;
    virtual bool isActive(const index_t i, const gsVector<T> & u) const;
    virtual void activeCoefs_into(const gsVector<T> & u, const gsMatrix<T> & coefs,
                                  gsMatrix<T>& result) const;
    virtual gsMatrix<index_t> allBoundary( ) const;
    virtual gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset) const;
    gsMatrix<index_t> boundary(boxSide const & s) const
    { return this->boundaryOffset(s,0); }
    virtual index_t functionAtCorner(boxCorner const & c) const;
    virtual gsMatrix<T> support() const;
    virtual gsMatrix<T> support(const index_t & i) const;
    gsMatrix<T> supportInterval(index_t dir) const;
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;
    virtual void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    virtual void derivSingle_into(index_t i,
                                  const gsMatrix<T> & u,
                                  gsMatrix<T>& result ) const;
    virtual void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    virtual void deriv2Single_into(index_t i,
                                   const gsMatrix<T> & u,
                                   gsMatrix<T>& result ) const;
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> > & result,
                                  bool sameElement = false) const;
    GISMO_DEPRECATED
    virtual void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u,
                                        int n, gsMatrix<T>& result) const;
    virtual void evalDerSingle_into(index_t i, const
                                    gsMatrix<T> & u, int n,
                                    gsMatrix<T>& result) const;
    virtual typename gsBasis<T>::uPtr create() const;
    virtual const gsSemiRegularBezier & source () const { return *this; }
    virtual gsSemiRegularBezier & source () { return *this; }

    virtual size_t numElements(boxSide const & s = 0) const;
    virtual size_t elementIndex(const gsVector<T> & u ) const;
    virtual gsMatrix<T> elementInSupportOf(index_t j) const;
    virtual void refine(gsMatrix<T> const & boxes, int refExt = 0);
    virtual void unrefine(gsMatrix<T> const & boxes, int refExt = 0);
    virtual std::vector<index_t> asElements(gsMatrix<T> const & boxes, int refExt = 0) const;
    virtual std::vector<index_t> asElementsUnrefine(gsMatrix<T> const & boxes, int refExt = 0) const;
    virtual void refineElements(std::vector<index_t> const & boxes);
    virtual void unrefineElements(std::vector<index_t> const & boxes);
    virtual void refineElements_withCoefs(gsMatrix<T> & coefs,std::vector<index_t> const & boxes);
    virtual void unrefineElements_withCoefs(gsMatrix<T> & coefs,std::vector<index_t> const & boxes);
    virtual void uniformRefine(int numKnots = 1, int mul=1, short_t dir=-1);
    virtual void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1, short_t const dir = -1);
    virtual void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor>& transfer,
                                            int numKnots = 1, int mul = 1);
    virtual void uniformCoarsen(int numKnots = 1);
    virtual void uniformCoarsen_withCoefs(gsMatrix<T>& coefs, int numKnots = 1);
    virtual void uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor>& transfer,
                                             int numKnots = 1);
    virtual void degreeElevate(short_t const & i = 1, short_t const dir = -1);
    virtual void degreeReduce(short_t const & i = 1, short_t const dir = -1);
    virtual void degreeIncrease(short_t const & i = 1, short_t const dir = -1);
    virtual void degreeDecrease(short_t const & i = 1, short_t const dir = -1);
    void setDegree(short_t const& i);
    void setDegreePreservingMultiplicity(short_t const& i);
    virtual void elevateContinuity(int const & i = 1);
    virtual void reduceContinuity(int const & i = 1);

    virtual short_t maxDegree() const;
    virtual short_t minDegree() const;
    virtual short_t totalDegree() const;
    virtual short_t degree(short_t i) const;
    memory::unique_ptr<gsGeometry<T> > interpolateData(gsMatrix<T> const& vals,
                                    gsMatrix<T> const& pts ) const;
    virtual void reverse();
    virtual void matchWith(const boundaryInterface & bi, const gsSemiRegularBezier<2,T> & other,
                           gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther, index_t offset = 0) const;
protected:
};

template<class T>
gsSemiRegularBezier<2,T>::gsSemiRegularBezier(int deg)
{
    const int numV = deg + 1; // number of rows

    m_v.reserve(numV);
    m_offset.resize(numV + 1);
    m_offset[0] = 0;

    // Fixed U-direction Bézier basis: degree deg, one segment
    gsKnotVector<T> knotU(0.0, 1.0, 0.0, deg + 1);
    m_u = bsbasis(knotU);

    // Build V-direction Bézier bases per row: degrees 0 to deg
    for (int i = 0; i < numV; ++i)
    {
        int vDeg = i; // degree of the ith row in V-direction

        std::vector<T> knotsV;
        // start and end multiplicities = degree + 1 ⇒ Bézier
        for (int j = 0; j < vDeg + 1; ++j)
            knotsV.push_back(0.0);
        for (int j = 0; j < vDeg + 1; ++j)
            knotsV.push_back(1.0);

        gsKnotVector<T> kv(knotsV, vDeg);
        m_v.push_back(bsbasis(kv));

        m_offset[i + 1] = m_offset[i] + m_v.back().size(); // cumulative basis function count
    }
}



template<class T>
index_t gsSemiRegularBezier<2,T>::size() const
{
    index_t sz = 0;
    for(auto & bs : m_v)
        sz += bs.size();
    return sz;
}

template<class T>
std::ostream & gsSemiRegularBezier<2,T>::print(std::ostream &os) const
{
    gsInfo<<"Semi-regular basis  with "<<size()<<" functions.\n";
    for(size_t i = 0; i!=m_v.size(); ++i)
        os<< "i="<<i<<": " << m_v[i].knots()<<"\n";
    return os;
}

template<class T> inline
void gsSemiRegularBezier<2,T>::anchors_into(gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T> inline
void gsSemiRegularBezier<2,T>::anchor_into(index_t, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::connectivity(const gsMatrix<T> &, gsMesh<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
{
    const index_t nact = (m_u.degree()+1);
    result.resize( nact*nact, u.cols() );
    gsMatrix<index_t> actV, actU;
    
    for (index_t i = 0; i < u.cols(); ++i) // for all columns of u
    {     
        actU = m_u.active( u.row(0).col(i) );
        for (index_t a = 0; a!=actU.size(); ++a )
        {
            actV = m_v[actU.at(a)].active( u.row(1).col(i) );
            actV.array() += m_offset[ actU.at(a) ];
            result.col(i).segment(a*nact, nact) = actV;
            
        }
    }
}



template <class T>
bool gsSemiRegularBezier<2,T>::isActive(const index_t, const gsVector<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::numActive_into(const gsMatrix<T> &, gsVector<index_t>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::activeCoefs_into(const gsVector<T> &, const gsMatrix<T> &,
                                  gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<index_t>
gsSemiRegularBezier<2,T>::allBoundary() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<index_t>
gsSemiRegularBezier<2,T>::boundaryOffset(boxSide const &,index_t) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
index_t
gsSemiRegularBezier<2,T>::functionAtCorner(boxCorner const &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsSemiRegularBezier<2,T>::support() const
{
    gsMatrix<T> box(2,2);
    box.row(0) = m_u.support();
    box.row(1) = m_v.front().support();
    return box;
}

template<class T>
gsMatrix<T> gsSemiRegularBezier<2,T>::support(const index_t & i) const
{
    index_t a = std::upper_bound(m_offset.begin(), m_offset.end(),i) - m_offset.begin() - 1;
    gsMatrix<T> box(2,2);
    box.row(0) = m_u.support(a);
    box.row(1) = m_v[a].support(i-m_offset[a]);
    return box;
}

template<class T>
gsMatrix<T> gsSemiRegularBezier<2,T>::supportInterval(index_t dir) const
{ return support().row(dir); }

template<class T>
void gsSemiRegularBezier<2,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    const index_t nact = (m_u.degree()+1);
    result.resize( nact*nact, u.cols() );
    gsMatrix<index_t> actU;
    gsMatrix<T> valU, valV;
    
    for (index_t i = 0; i < u.cols(); ++i) // for all columns of u
    {
        actU = m_u.active( u.row(0).col(i) );
        valU = m_u.eval  ( u.row(0).col(i) );
        for (index_t a = 0; a!=actU.size(); ++a )
        {
            valV = m_v[actU.at(a)].eval( u.row(1).col(i) );
            result.col(i).segment(a*nact, nact) = valU.at(a) * valV;
        }
    }
}

template<class T>
void gsSemiRegularBezier<2,T>::evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    index_t a = std::upper_bound(m_offset.begin(), m_offset.end(),i) - m_offset.begin() - 1;
    result = m_u.evalSingle(a,u.row(0)).array() * 
        m_v[a].evalSingle(i-m_offset[a],u.row(1)).array();
}

template<class T>
void gsSemiRegularBezier<2,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    const index_t nact = (m_u.degree()+1);
    const index_t nder = 2;
    result.resize( nact*nact*nder, u.cols() );
    gsMatrix<index_t> actU;
    std::vector<gsMatrix<T> > valU, valV;

    for (index_t i = 0; i < u.cols(); ++i) // for all columns of u
    {
        actU = m_u.active( u.row(0).col(i) );
        valU = m_u.evalAllDers(u.row(0).col(i), 1);
        for (index_t a = 0; a!=actU.size(); ++a )
        {
            valV = m_v[actU.at(a)].evalAllDers( u.row(1).col(i), 1);
            for (index_t k = 0; k!=nact; ++k)
            {
                result((a*nact+k)*nder  ,i) = valU[1].at(a) * valV[0].at(k);
                result((a*nact+k)*nder+1,i) = valU[0].at(a) * valV[1].at(k);
            }
        }
    }
}

template<class T>
void gsSemiRegularBezier<2,T>::derivSingle_into(index_t,
                                  const gsMatrix<T> &,
                                  gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    const index_t nact = (m_u.degree()+1);
    const index_t nder = 3;
    result.resize( nact*nact*nder, u.cols() );
    gsMatrix<index_t> actU;
    std::vector<gsMatrix<T> > valU, valV;

    for (index_t i = 0; i < u.cols(); ++i) // for all columns of u
    {
        actU = m_u.active( u.row(0).col(i) );
        valU = m_u.evalAllDers(u.row(0).col(i), 2);
        for (index_t a = 0; a!=actU.size(); ++a )
        {
            valV = m_v[actU.at(a)].evalAllDers( u.row(1).col(i), 2);
            for (index_t k = 0; k!=nact; ++k)
            {
                result((a*nact+k)*nder,i  ) = valU[2].at(a) * valV[0].at(k);
                result((a*nact+k)*nder+1,i) = valU[0].at(a) * valV[2].at(k);
                result((a*nact+k)*nder+2,i) = valU[1].at(a) * valV[1].at(k);
            }
        }
    }
}

template<class T>
void gsSemiRegularBezier<2,T>::deriv2Single_into(index_t,
                                   const gsMatrix<T> &,
                                   gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template <typename T> //__attribute__ ((fallthrough)) // todo
void gsSemiRegularBezier<2,T>::evalAllDers_into(const gsMatrix<T> & u, const int n,
                                  std::vector<gsMatrix<T> > & result,
                                  bool sameElement) const
{
    GISMO_UNUSED(sameElement);
    //gsWarn << "generic evalAllDers called from "<<typeid(*this).name()<< "\n";
    result.resize(n+1);

    switch(n)
    {
    case 0:
        eval_into(u, result[0]);
        break;
    case 1:
        eval_into (u, result[0]);
        deriv_into(u, result[1]);
        break;
    case 2:
        eval_into  (u, result[0]);
        deriv_into (u, result[1]);
        deriv2_into(u, result[2]);
        break;
    default:
        GISMO_ERROR("evalAllDers implemented for order up to 2<"<<n ); //<< " for "<<*this);
        break;
    }
}

template<class T>
void gsSemiRegularBezier<2,T>::evalAllDersSingle_into(index_t, const gsMatrix<T> &,
                                        int, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::evalDerSingle_into(index_t, const
                                    gsMatrix<T> &, int,
                                    gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
typename gsBasis<T>::uPtr gsSemiRegularBezier<2,T>::create() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsSemiRegularBezier<2,T>::numElements(boxSide const &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsSemiRegularBezier<2,T>::elementIndex(const gsVector<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsSemiRegularBezier<2,T>::elementInSupportOf(index_t) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
std::vector<index_t> gsSemiRegularBezier<2,T>::asElements(gsMatrix<T> const &, int) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
std::vector<index_t> gsSemiRegularBezier<2,T>::asElementsUnrefine(gsMatrix<T> const &, int) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::refine(gsMatrix<T> const &, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::unrefine(gsMatrix<T> const &, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::refineElements(std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::unrefineElements(std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::refineElements_withCoefs(gsMatrix<T> &,std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::unrefineElements_withCoefs(gsMatrix<T> &,std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::uniformRefine(int, int, short_t)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::uniformRefine_withCoefs(gsMatrix<T>& , int , int , short_t const )
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> &,
                                            int, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::uniformCoarsen(int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::uniformCoarsen_withCoefs(gsMatrix<T>& , int )
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> &,
                                            int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::degreeElevate(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::degreeReduce(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::degreeIncrease(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::degreeDecrease(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::elevateContinuity(int const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::reduceContinuity(int const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsSemiRegularBezier<2,T>::maxDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsSemiRegularBezier<2,T>::minDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsSemiRegularBezier<2,T>::totalDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsSemiRegularBezier<2,T>::degree(short_t k) const
{ return 0==k ? m_u.degree() : m_v.front().degree(); }

template<class T>
void gsSemiRegularBezier<2,T>::reverse()
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBezier<2,T>::matchWith(const boundaryInterface &, const gsSemiRegularBezier<2,T> &,
               gsMatrix<index_t> &, gsMatrix<index_t> &, index_t) const
{ GISMO_NO_IMPLEMENTATION }

}; // namespace gismo