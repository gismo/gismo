/** @file gsSemiRegularBasis.h

    @brief A semiregular B-spline basis

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
class gsSemiRegularBasis
{ };

template<class T>
class gsSemiRegularBasis<2,T> : public gsBasis<T>
{
    typedef gsBSplineBasis<T> bsbasis;
    
    bsbasis m_u;
    std::vector<bsbasis> m_v;
    std::vector<index_t> m_offset;
    std::vector<index_t> qLevels; // q-levels for each v-basis
    index_t m_deg;

public:
    typedef memory::shared_ptr< gsSemiRegularBasis > Ptr;
    typedef memory::unique_ptr< gsSemiRegularBasis > uPtr;
    typedef T Scalar_t;
    static const bool IsRational = false;
    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

    gsSemiRegularBasis(int deg = 1, int numLayers = 1);



public:
    // these are virtual functions that we must implement
    gsSemiRegularBasis * clone_impl() const
    { return new gsSemiRegularBasis(*this); }

    std::ostream &print(std::ostream &os) const;
    short_t domainDim() const { return 2; }
    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    {
        gsGenericGeometry<2,T> * g = new gsGenericGeometry<2,T>(*this, give(coefs) );
        return memory::unique_ptr<gsGeometry<T> >(g);
    }

public:
    index_t size() const;
    // Accessors (getters)
    const gsBSplineBasis<T> & getUBasis() const;
    const std::vector<gsBSplineBasis<T>> & getVBases() const;
    const gsBSplineBasis<T> & getVBasisAt(index_t i) const;
    const gsVector<index_t> getqLevels() const;
    
;


    // Mutators (setters)
    void setUBasis(const gsBSplineBasis<T> & basis);
    void setVBasisAt(index_t i, const gsBSplineBasis<T> & basis);
    void decrease_by_one();
    void setqLevels(const gsVector<index_t>& newLevels);
    void updateVBases();
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
    virtual const gsSemiRegularBasis & source () const { return *this; }
    virtual gsSemiRegularBasis & source () { return *this; }

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
    virtual void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                           gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther, index_t offset = 0) const;

    virtual memory::shared_ptr<gsDomain<T> > domain() const
    {
    
    static gsTensorBSplineBasis<2> basis2dtp(m_u.knots(), m_u.knots());
    return basis2dtp.domain();
    }

    void swapDirections(unsigned a, unsigned b)
    {
        std::swap(m_dirMap[a], m_dirMap[b]);
        std::swap(m_flip[a],   m_flip[b]);
    }

    void reverseDirection(unsigned d)
    {
        m_flip[d] = !m_flip[d];
    }

protected:
std::array<unsigned,2> m_dirMap {{0,1}};   // logical (u,v) → canonical (U,V)
std::array<bool,2>     m_flip   {{false,false}}; // flip logical axis?

};

template<class T>
gsSemiRegularBasis<2,T>::gsSemiRegularBasis(int deg, int numLayers)
{
    m_deg = deg;
    const int numV = (1 << numLayers) + deg; // Total vertical basis rows
    m_v.reserve(numV);
    m_offset.resize(numV + 1);
    m_offset[0] = 0;

    // Construct shared U-direction knot vector
    std::vector<T> knotsU;
    int numInteriorU = (1 << numLayers) - 1; // 2^l - 1
    for (int j = 0; j < deg + 1; ++j)
        knotsU.push_back(0.0);
    for (int j = 1; j <= numInteriorU; ++j)
        knotsU.push_back(static_cast<T>(j) / (1 << numLayers));
    for (int j = 0; j < deg + 1; ++j)
        knotsU.push_back(1.0);

    // Set shared basis in U direction
    gsKnotVector<T> kvU(knotsU, deg);
    m_u = bsbasis(kvU);

    // Build V-direction basis functions per row
    //reserve q to be a vector with numV elements where each element is = numLayers
    // This means that the knot vector for all vBases is the same as knotsU
    qLevels.resize(numV);
    for (int i = 0; i < numV; ++i)
    {
        qLevels[i] = numLayers; //constant means global knot vector for all vBases same as knotsU
    }
    for (int i = 0; i < numV; ++i)
    {
        unsigned q = qLevels[i]; //constant means global knot vector for all vBases same as knotsU


        std::vector<T> knotsV;
        int numInteriorV = (1 << q) - 1;
        
        for (int j = 0; j < deg + 1; ++j)
            knotsV.push_back(0.0);
        for (int j = 1; j <= numInteriorV; ++j)
            knotsV.push_back(static_cast<T>(j) / (1 << q));
        for (int j = 0; j < deg + 1; ++j)
            knotsV.push_back(1.0);

        gsKnotVector<T> kv(knotsV, deg);
        m_v.push_back(bsbasis(kv));
        m_offset[i + 1] = m_offset[i] + m_v.back().size();
    }
}


template<class T>
const gsBSplineBasis<T> & gsSemiRegularBasis<2,T>::getUBasis() const
{
    return m_u;
}

template<class T>
const std::vector<gsBSplineBasis<T>> & gsSemiRegularBasis<2,T>::getVBases() const
{
    return m_v;
}

template<class T>
const gsBSplineBasis<T> & gsSemiRegularBasis<2,T>::getVBasisAt(index_t i) const
{
    GISMO_ASSERT(i < m_v.size(), "Index out of bounds in getVBasisAt");
    return m_v[i];
}

template<class T>
void gsSemiRegularBasis<2,T>::setUBasis(const gsBSplineBasis<T> & basis)
{
    m_u = basis;
}

template<class T>
void gsSemiRegularBasis<2,T>::setVBasisAt(index_t i, const gsBSplineBasis<T> & basis)
{
    GISMO_ASSERT(i < m_v.size(), "Index out of bounds in setVBasisAt");
    m_v[i] = basis;
}

//getter for qLevels of a semi-regular basis object
template <typename T>
const gsVector<index_t> gsSemiRegularBasis<2,T>::getqLevels() const
{
    gsVector<index_t> v(qLevels.size());
    for (size_t i = 0; i < qLevels.size(); ++i)
        v[i] = qLevels[i];
    return v;
}


template <typename T>
void gsSemiRegularBasis<2,T>::decrease_by_one()
{
    for (size_t i = 0; i < qLevels.size(); ++i)
    {
        if (qLevels[i] > 0)
        {
            qLevels[i] -= 1;
            break; // stop after decreasing one element
        }
    }
}

template <typename T>
void gsSemiRegularBasis<2,T>::setqLevels(const gsVector<index_t>& newLevels)
{
    qLevels.resize(newLevels.size());
    for (index_t i = 0; i < newLevels.size(); ++i)
        qLevels[i] = newLevels[i];
}


template<class T>
void gsSemiRegularBasis<2,T>::updateVBases()
{
    m_v.clear();
    m_offset[0] = 0;
    int numV = qLevels.size();

    for (int i = 0; i < numV; ++i)
    {
        unsigned q = qLevels[i];
        std::vector<T> knotsV;
        int numInteriorV = (1 << q) - 1;

        for (int j = 0; j < m_deg + 1; ++j)
            knotsV.push_back(0.0);

        for (int j = 1; j <= numInteriorV; ++j)
            knotsV.push_back(static_cast<T>(j) / (1 << q));

        for (int j = 0; j < m_deg + 1; ++j)
            knotsV.push_back(1.0);

        gsKnotVector<T> kv(knotsV, m_deg);
        m_v.push_back(bsbasis(kv));

        m_offset[i + 1] = m_offset[i] + m_v.back().size();
    }
}













template<class T>
index_t gsSemiRegularBasis<2,T>::size() const
{
    index_t sz = 0;
    for(auto & bs : m_v)
        sz += bs.size();
    return sz;
}

template<class T>
std::ostream & gsSemiRegularBasis<2,T>::print(std::ostream &os) const
{
    gsInfo<<"Adaptive semi-regular basis  with "<<size()<<" functions.\n";
    for(size_t i = 0; i!=m_v.size(); ++i)
        os<< "i="<<i<<": " << m_v[i].knots()<<"\n";
    return os;
}






template<class T> inline
void gsSemiRegularBasis<2,T>::anchors_into(gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T> inline
void gsSemiRegularBasis<2,T>::anchor_into(index_t, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::connectivity(const gsMatrix<T> &, gsMesh<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
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
bool gsSemiRegularBasis<2,T>::isActive(const index_t, const gsVector<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::numActive_into(const gsMatrix<T> & u,
                                             gsVector<index_t>& result) const
{
    result.resize(u.cols());
    gsMatrix<index_t> actives;

    for (index_t i = 0; i < u.cols(); ++i)
    {
        // Extract active basis indices at this parametric point
        this->active_into(u.col(i), actives);

        // The number of active functions is just the number of rows
        result(i) = actives.rows();
    }
}


template<class T>
void gsSemiRegularBasis<2,T>::activeCoefs_into(const gsVector<T> &, const gsMatrix<T> &,
                                  gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<index_t>
gsSemiRegularBasis<2,T>::allBoundary() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<index_t>
gsSemiRegularBasis<2,T>::boundaryOffset(const boxSide& s, index_t offset) const
{

    // Map logical side (from outside callers) to this basis' canonical axes
    const unsigned ldir = s.direction();          // logical axis: 0(u) or 1(v)
    const bool     lmax = s.parameter();          // logical 0/1-side

    const unsigned cdir = m_dirMap[ldir];         // canonical axis
    const bool     cmax = m_flip[ldir] ? !lmax : lmax;

    const boxSide  canonSide(cdir, cmax);

    // From here on, use canonSide instead of s:
    const short_t  dir   = canonSide.direction();
    const bool     isMax = canonSide.parameter();


    const index_t nRows = static_cast<index_t>(m_v.size());
    GISMO_ASSERT(nRows > 0, "Semi-regular basis has no rows.");

    if (dir == 1)
    {
        // - v-constant edge (angular-edge): return the whole row j in increasing-u order
        const index_t j = isMax ? (nRows - 1 - offset) : offset;
        GISMO_ASSERT(j < nRows, "Offset exceeds number of v-rows: j=" << j << ", nRows=" << nRows);

        const index_t rowLen = m_v[j].size();
        const index_t base   = m_offset[j];

        gsMatrix<index_t> res(rowLen, 1);
        for (index_t i = 0; i < rowLen; ++i)
            res(i,0) = base + i;      // global index of (j,i)

        return res;
    }
    else
    {
        // -u-constant edge (radius-edge): one entry per row, ordered by increasing v (row index)
        // Ensure all rows are big enough for this offset
        for (index_t j = 0; j < nRows; ++j)
        {
            const index_t rowLen = m_v[j].size();
            GISMO_ASSERT(offset < rowLen,
                "Row " << j << " too short for requested offset " << offset
                       << " (rowLen=" << rowLen << ")");
        }

        gsMatrix<index_t> res(nRows, 1);
        for (index_t j = 0; j < nRows; ++j)
        {
            const index_t rowLen = m_v[j].size();
            const index_t i      = isMax ? (rowLen - 1 - offset) : offset;
            res(j,0)             = m_offset[j] + i;  // global index of (j,i)
        }
        return res;
    }
}

template<class T>
index_t
gsSemiRegularBasis<2,T>::functionAtCorner(boxCorner const &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsSemiRegularBasis<2,T>::support() const
{
    gsMatrix<T> box(2,2);
    box.row(0) = m_u.support();
    box.row(1) = m_v.front().support();
    return box;
}

template<class T>
gsMatrix<T> gsSemiRegularBasis<2,T>::support(const index_t & i) const
{
    index_t a = std::upper_bound(m_offset.begin(), m_offset.end(),i) - m_offset.begin() - 1;
    gsMatrix<T> box(2,2);
    box.row(0) = m_u.support(a);
    box.row(1) = m_v[a].support(i-m_offset[a]);
    return box;
}

template<class T>
gsMatrix<T> gsSemiRegularBasis<2,T>::supportInterval(index_t dir) const
{ return support().row(dir); }

template<class T>
void gsSemiRegularBasis<2,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
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
void gsSemiRegularBasis<2,T>::evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    index_t a = std::upper_bound(m_offset.begin(), m_offset.end(),i) - m_offset.begin() - 1;
    result = m_u.evalSingle(a,u.row(0)).array() * 
        m_v[a].evalSingle(i-m_offset[a],u.row(1)).array();
}

template<class T>
void gsSemiRegularBasis<2,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
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
void gsSemiRegularBasis<2,T>::derivSingle_into(index_t,
                                  const gsMatrix<T> &,
                                  gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
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
void gsSemiRegularBasis<2,T>::deriv2Single_into(index_t,
                                   const gsMatrix<T> &,
                                   gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template <typename T> //__attribute__ ((fallthrough)) // todo
void gsSemiRegularBasis<2,T>::evalAllDers_into(const gsMatrix<T> & u, const int n,
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
void gsSemiRegularBasis<2,T>::evalAllDersSingle_into(index_t, const gsMatrix<T> &,
                                        int, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::evalDerSingle_into(index_t, const
                                    gsMatrix<T> &, int,
                                    gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
typename gsBasis<T>::uPtr gsSemiRegularBasis<2,T>::create() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsSemiRegularBasis<2,T>::numElements(boxSide const &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsSemiRegularBasis<2,T>::elementIndex(const gsVector<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsSemiRegularBasis<2,T>::elementInSupportOf(index_t) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
std::vector<index_t> gsSemiRegularBasis<2,T>::asElements(gsMatrix<T> const &, int) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
std::vector<index_t> gsSemiRegularBasis<2,T>::asElementsUnrefine(gsMatrix<T> const &, int) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::refine(gsMatrix<T> const &, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::unrefine(gsMatrix<T> const &, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::refineElements(std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::unrefineElements(std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::refineElements_withCoefs(gsMatrix<T> &,std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::unrefineElements_withCoefs(gsMatrix<T> &,std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::uniformRefine(int, int, short_t)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::uniformRefine_withCoefs(gsMatrix<T>& , int , int , short_t const )
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> &,
                                            int, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::uniformCoarsen(int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::uniformCoarsen_withCoefs(gsMatrix<T>& , int )
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> &,
                                            int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::degreeElevate(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::degreeReduce(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::degreeIncrease(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::degreeDecrease(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::elevateContinuity(int const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::reduceContinuity(int const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsSemiRegularBasis<2,T>::maxDegree() const
{
    return m_u.degree();
}

template<class T>
short_t gsSemiRegularBasis<2,T>::minDegree() const
{
    return m_u.degree();
}

template<class T>
short_t gsSemiRegularBasis<2,T>::totalDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsSemiRegularBasis<2,T>::degree(short_t k) const
{ return 0==k ? m_u.degree() : m_v.front().degree(); }

template<class T>
void gsSemiRegularBasis<2,T>::reverse()
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsSemiRegularBasis<2,T>::matchWith(const boundaryInterface & bi,
                                        const gsBasis<T>       & other,
                                        gsMatrix<index_t>      & bndThis,
                                        gsMatrix<index_t>      & bndOther,
                                        index_t offset) const
{
    // 1. Query both boundaries
    const boxSide sideThis  = bi.first().side();
    const boxSide sideOther = bi.second().side();

    bndThis  = this->boundaryOffset(sideThis, offset);
    bndOther = other.boundaryOffset(sideOther, offset);

    // 2. Sizes must match
    GISMO_ASSERT(bndThis.rows() == bndOther.rows(),
        "Boundary mismatch: " 
        << bndThis.rows() << " vs " << bndOther.rows());

    if (bndThis.rows() <= 1)
        return; // 1 dof – trivial

    // 3. Get orientation data from topology
    const index_t s1 = bi.first().direction();           // normal
    const gsVector<bool> & dirOr = bi.dirOrientation();  // true=same, false=flipped

    // Tangential direction = 1 - s1
    const index_t t1 = 1 - s1;

    // 4. If flipped along tangential direction, reverse order
    if (!dirOr[t1])
    {
        for (index_t i = 0, j = bndOther.rows() - 1; i < j; ++i, --j)
            bndOther.row(i).swap(bndOther.row(j));
    }
}



}; // namespace gismo


