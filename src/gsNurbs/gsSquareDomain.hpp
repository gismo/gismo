/** @file gsSquareDomain.hpp

    @brief Implementation of the gsSquareDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsCore/gsFuncData.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsBoxTopology.h>
#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>
#include <gsUtils/gsPointGrid.h>
#include <algorithm>
#include <limits>

namespace gismo
{

template <class T>
gsSquareDomain<T>::gsSquareDomain(const gsGeometry<T> & domain, bool slide)
{
    m_domain = memory::make_unique(domain.clone().release());
    // "Slide" must be in m_options before applyOptions() runs _initMapper,
    // which reads it.
    m_options.addSwitch("Slide", "Boundary controls slide along the boundary", slide);
    this->applyOptions();
}

template <class T>
gsSquareDomain<T>::gsSquareDomain(const gsBasis<T> & basis, bool slide)
{
    gsMatrix<T> coefs = basis.anchors();
    coefs.transposeInPlace();
    m_domain = basis.makeGeometry(give(coefs));
    m_options.addSwitch("Slide", "Boundary controls slide along the boundary", slide);
    this->applyOptions();
}

// Copy constructor
template <class T>
gsSquareDomain<T>::gsSquareDomain(const gsSquareDomain<T> & other)
{
    m_domain = memory::make_unique(other.domain().clone().release());
    m_mapper = other.mapper();
    m_interfaces = other.m_interfaces;
    m_indices = other.m_indices;
    m_options = other.m_options;
}

// Copy assignment
template <class T>
gsSquareDomain<T> & gsSquareDomain<T>::operator=(const gsSquareDomain<T> & other)
{
    if (this != &other)
    {
        m_domain = memory::make_unique(other.domain().clone().release());
        m_mapper = other.mapper();
        m_interfaces = other.m_interfaces;
        m_indices = other.m_indices;
        m_options = other.m_options;
    }
    return *this;
}

template <class T>
void gsSquareDomain<T>::addInterface(const boundaryInterface & bi)
{
    m_interfaces.push_back(bi);
}

template <class T>
gsOptionList & gsSquareDomain<T>::options()
{
    return m_options;
}

template <class T>
void gsSquareDomain<T>::applyOptions()
{
    this->_initMapper(*m_domain,m_mapper);
    this->_initIndices(*m_domain,m_mapper,m_indices);
}

template <class T>
const gsGeometry<T> & gsSquareDomain<T>::domain() const
{
    return *m_domain;
}

template <class T>
const gsDofMapper & gsSquareDomain<T>::mapper() const
{
    return m_mapper;
}

template <class T>
gsMatrix<T> gsSquareDomain<T>::support() const
{
    return m_domain->support();
}

template <class T>
short_t gsSquareDomain<T>::maxDegree() const
{
    return m_domain->basis().maxDegree();
}

template <class T>
short_t gsSquareDomain<T>::domainDim() const
{
    return m_domain->domainDim();
}

template <class T>
short_t gsSquareDomain<T>::targetDim() const
{
    return m_domain->targetDim();
}

template <class T>
void gsSquareDomain<T>::eval_into(const gsMatrix<T> & param, gsMatrix<T> & result) const
{
    m_domain->eval_into(param,result);
}

template <class T>
void gsSquareDomain<T>::deriv_into(const gsMatrix<T> & param, gsMatrix<T> & result) const
{
    m_domain->deriv_into(param,result);
}

template <class T>
void gsSquareDomain<T>::compute(const gsMatrix<T> & param, gsFuncData<T> & out) const
{
    // Single-pass evaluation -- see the note in the header.
    m_domain->compute(param,out);
}

template <class T>
void gsSquareDomain<T>::deriv2_into(const gsMatrix<T> & param, gsMatrix<T> & result) const
{
    m_domain->deriv2_into(param,result);
}

template <class T>
void gsSquareDomain<T>::setControls(const gsVector<T> & controls)
{
    GISMO_ENSURE((size_t)controls.rows()==m_indices.size(),"Wrong size of controls vector");
    index_t ii;
    for (index_t d = 0; d!=m_domain->domainDim(); d++)
        for (index_t k = 0; k!=m_domain->coefs().rows(); k++)
        {
            ii = m_mapper.index(k,0,d);
            if (m_mapper.is_free_index(ii))
                m_domain->coefs()(k,d) = controls[ii];
        }
}

template <class T>
gsVector<T> gsSquareDomain<T>::getControls() const
{
    gsVector<T> controls(m_indices.size());
    index_t ii;
    for (index_t d = 0; d!=m_domain->domainDim(); d++)
        for (index_t k = 0; k!=m_domain->coefs().rows(); k++)
        {
            ii = m_mapper.index(k,0,d);
            if (m_mapper.is_free_index(ii))
                controls[ii] = m_domain->coefs()(k,d);
        }
    return controls;
}

// template <class T>
// const T & gsSquareDomain<T>::control(index_t i) const
// {
//     return m_domain->coefs()(m_indices[i].first,m_indices[i].second);
// }

// template <class T>
// T & gsSquareDomain<T>::control(index_t i)
// {
//     return m_domain->coefs()(m_indices[i].first,m_indices[i].second);
// }

template <class T>
size_t gsSquareDomain<T>::nControls() const
{
    return m_indices.size();
}

template <class T>
void gsSquareDomain<T>::control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const
{
    gsMatrix<T> tmp;

    result.resize(targetDim()*nControls(), points.cols());
    result.setZero();
    for (index_t p = 0; p!=points.cols(); p++)
    {
        gsAsMatrix<T> res = result.reshapeCol(p,nControls(),targetDim());
        for (index_t k = 0; k!=m_domain->coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain->targetDim(); d++)
                if (m_mapper.is_free(k,0,d))
                {
                    // @hverhelst: Why don't we use eval_into and put the actives in the right place?
                    m_domain->basis().evalSingle_into(k,points.col(p),tmp); // evaluate basis function k
                    res(m_mapper.index(k,0,d),d) = tmp(0,0); // tmp is a single value (1 point, 1 basis function)
                }
    }
}

template <class T>
void gsSquareDomain<T>::detFromJacobian_into(const gsMatrix<T> & jacobians, short_t dd,
                                             gsMatrix<T> & result)
{
    GISMO_ENSURE(dd >= 1 && jacobians.rows() == dd*dd,
                 "detFromJacobian_into: jacobians must be (dd*dd) x nPoints, got "
                 << jacobians.rows() << " rows for dd = " << dd);
    const index_t numPt = jacobians.cols();
    result.resize(1, numPt);
    if (dd == 2)
    {
        for (index_t q = 0; q != numPt; ++q)
            result(0, q) = jacobians(0, q) * jacobians(3, q) - jacobians(1, q) * jacobians(2, q);
    }
    else if (dd == 3)
    {
        // Rule of Sarrus / cofactor expansion, mirroring the adj(J_c)^T
        // cofactors at gsAdaptiveParametrization.hpp:1130-1141, spelled out
        // explicitly to avoid the general reshapeCol().determinant() path's
        // LU factorization for this common case.
        for (index_t q = 0; q != numPt; ++q)
        {
            const T J00 = jacobians(0, q), J01 = jacobians(1, q), J02 = jacobians(2, q);
            const T J10 = jacobians(3, q), J11 = jacobians(4, q), J12 = jacobians(5, q);
            const T J20 = jacobians(6, q), J21 = jacobians(7, q), J22 = jacobians(8, q);
            result(0, q) = J00 * (J11 * J22 - J12 * J21)
                          - J01 * (J10 * J22 - J12 * J20)
                          + J02 * (J10 * J21 - J11 * J20);
        }
    }
    else
    {
        gsMatrix<T> Jc(dd, dd);
        for (index_t q = 0; q != numPt; ++q)
        {
            Jc = jacobians.reshapeCol(q, dd, dd);
            result(0, q) = Jc.determinant();
        }
    }
}

template <class T>
void gsSquareDomain<T>::detJacobian_into(const gsMatrix<T> & points, gsMatrix<T> & result) const
{
    GISMO_ENSURE(domainDim() == targetDim(),
                 "detJacobian_into requires a square Jacobian (domainDim == targetDim)");
    gsMatrix<T> dsig;
    this->deriv_into(points, dsig);
    detFromJacobian_into(dsig, domainDim(), result);
}

template <class T>
void gsSquareDomain<T>::detJacobianDeriv_into(const gsMatrix<T> & points, gsMatrix<T> & result,
                                              gsMatrix<T> * jacobian) const
{
    const index_t nc = nControls();
    const index_t dd = domainDim();
    GISMO_ENSURE(dd == targetDim(), "detJacobianDeriv_into requires a square Jacobian (domainDim == targetDim)");

    const index_t numPt = points.cols();
    result.resize(nc, numPt);
    result.setZero();
    if (jacobian)
        jacobian->resize(dd * dd, numPt);

    // Single basis-level pass (active_into + deriv_into, combined via compute())
    // instead of separately evaluating a geometry-level deriv_into (itself an
    // active_into+deriv_into pair) plus basis().deriv_into plus
    // basis().active_into -- all of those would walk the same points over
    // the same basis redundantly.
    gsFuncData<T> tmp(NEED_DERIV | NEED_ACTIVE);
    m_domain->basis().compute(points, tmp);
    const gsMatrix<T> & coefsFull = m_domain->coefs();

    gsMatrix<T> coefM, adjJT(dd, dd);
    for (index_t p = 0; p != numPt; ++p)
    {
        const typename gsMatrix<index_t>::constColumn active = tmp.active(p);
        const index_t nActive = active.rows();
        const typename gsFuncData<T>::matrixView derivs = tmp.deriv(p); // dd x nActive

        coefM.resize(nActive, dd);
        for (index_t loc = 0; loc != nActive; ++loc)
            coefM.row(loc) = coefsFull.row(active(loc));

        // dsigma_dxi(j,c) = d(sigma_c)/dxi_j -- exactly the component-major
        // (row c*dd+j) layout detJacobianDeriv_into's jacobian
        // out-param documents, so it is written out directly with no
        // transpose.
        const gsMatrix<T> dsigma_dxi = derivs * coefM;
        if (jacobian)
            jacobian->reshapeCol(p, dd, dd) = dsigma_dxi;

        // Jsigma(c,j) = d(sigma_c)/dxi_j
        gsMatrix<T> Jsigma = dsigma_dxi.transpose();

        // Cofactor matrix transpose (== det(J)*J^{-1,T} by Cramer's rule, but
        // computed directly): avoids the explicit inverse, whose internal
        // division by det(J) is exactly the quantity the fold barrier drives
        // to zero -- catastrophic cancellation right where this matters.
        if (dd == 2)
        {
            adjJT(0, 0) =  Jsigma(1, 1); adjJT(0, 1) = -Jsigma(1, 0);
            adjJT(1, 0) = -Jsigma(0, 1); adjJT(1, 1) =  Jsigma(0, 0);
        }
        else
        {
            const T detJ = Jsigma.determinant();
            adjJT = detJ * Jsigma.inverse().transpose();
        }

        for (index_t loc = 0; loc != nActive; ++loc)
        {
            const index_t k = active(loc);
            for (index_t d = 0; d != dd; ++d)
            {
                if (!m_mapper.is_free(k, 0, d))
                    continue;
                index_t ii = m_mapper.index(k, 0, d);
                T val = T(0);
                for (index_t j = 0; j != dd; ++j)
                    val += adjJT(d, j) * derivs(j, loc);
                // matchDof (see _initMapper / addInterface) can map two active coefficients
                // to the same free index, both active at a seam point: += is required here,
                // since = would silently drop whichever contribution is written first.
                result(ii, p) += val;
            }
        }
    }
}

template <class T>
void gsSquareDomain<T>::perturb(T factor)
{
    gsVector<T> controls = getControls();
    controls += factor * gsVector<T>::Random(m_indices.size());
    setControls(controls);
}

template <class T>
T gsSquareDomain<T>::minJacobian(index_t nPerElement) const
{
    GISMO_ENSURE(domainDim() == targetDim(),
                 "minJacobian requires a square Jacobian (domainDim == targetDim)");
    GISMO_ENSURE(nPerElement >= 2, "minJacobian: nPerElement must be >= 2 "
                 "(got " << nPerElement << "); 0 samples nothing and returns the "
                 "numeric_limits sentinel, 1 degenerates to each element's lower "
                 "corner");

    const gsBasis<T> & sb = m_domain->basis();
    gsVector<unsigned> np(domainDim());
    np.setConstant(nPerElement);
    gsMatrix<T> pts, det;
    T mn = std::numeric_limits<T>::max();
    for (auto & elem : sb.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        this->detJacobian_into(pts, det);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, det(0, q));
    }
    return mn;
}

template <class T>
template <short_t d>
typename gsGeometry<T>::uPtr gsSquareDomain<T>::_detJacobianSplineImpl(bool keepBezier) const
{
    typedef gsTensorBSpline<d,T> Spline;

    const gsTensorBSplineBasis<d,T> * sbasis =
        dynamic_cast<const gsTensorBSplineBasis<d,T> *>(&m_domain->basis());
    GISMO_ENSURE(sbasis != nullptr,
                 "detJacobianSpline requires sigma's basis to be a gsTensorBSplineBasis<d,T>");

    const gsMatrix<T> & coefsFull = m_domain->coefs();

    // J[i][j] = d(sigma_i)/d(xi_j), each a scalar spline on sigma's own basis.
    std::vector< std::vector<Spline> > J(d);
    for (short_t i = 0; i != d; ++i)
    {
        gsMatrix<T> ci = coefsFull.col(i);
        Spline comp(*sbasis, ci);
        for (short_t j = 0; j != d; ++j)
            J[i].push_back(comp.grad(j));
    }

    // Leibniz expansion det(J) = sum_perm sign(perm) * prod_i J[i][perm[i]]
    // over the d! permutations of (0..d-1). Permutations are generated in
    // lexicographic order via std::next_permutation; each one is signed
    // independently by its inversion-count parity (the number of pairwise
    // swaps needed to sort it back to the identity), since consecutive
    // permutations from std::next_permutation are not in general a single
    // transposition apart. d==2 reduces to the two-term form
    // J00*J11 - J01*J10; d==3 reproduces the six-term rule of Sarrus.
    std::vector<short_t> perm(d);
    for (short_t i = 0; i != d; ++i) perm[i] = i;

    Spline det;
    bool first = true;
    do
    {
        index_t inversions = 0;
        for (short_t i = 0; i != d; ++i)
            for (short_t j = i+1; j != d; ++j)
                if (perm[i] > perm[j]) ++inversions;
        const T sign = (inversions % 2 == 0) ? T(1) : T(-1);

        Spline term = J[0][perm[0]];
        for (short_t i = 1; i != d; ++i)
            term = Spline::multiply(term, J[i][perm[i]], keepBezier);

        if (first) { det = give(term); first = false; }
        else       det = Spline::linearCombination(T(1), det, sign, term);
    } while (std::next_permutation(perm.begin(), perm.end()));

    return typename gsGeometry<T>::uPtr(new Spline(give(det)));
}

template <class T>
typename gsGeometry<T>::uPtr gsSquareDomain<T>::detJacobianSpline(bool keepBezier) const
{
    GISMO_ENSURE(domainDim() == targetDim(),
                 "detJacobianSpline requires a square Jacobian (domainDim == targetDim)");
    switch (domainDim())
    {
        case 2: return _detJacobianSplineImpl<2>(keepBezier);
        case 3: return _detJacobianSplineImpl<3>(keepBezier);
        default: GISMO_ERROR("detJacobianSpline: unsupported domainDim() " << domainDim());
    }
}

template <class T>
T gsSquareDomain<T>::minDetJCoefficient() const
{
    return this->detJacobianSpline()->coefs().minCoeff();
}

template <class T>
void gsSquareDomain<T>::_initMapper(const gsGeometry<T> & domain, gsDofMapper & mapper) const
{
    mapper = gsDofMapper(domain.basis(),domain.domainDim());

    gsBoxTopology topology(domain.domainDim(),1);
    topology.addAutoBoundaries();
    for (typename gsBoxTopology::biterator it = topology.bBegin(); it != topology.bEnd(); ++it)
    {
        gsMatrix<index_t> boundary = domain.basis().boundary(*it);
        if (m_options.askSwitch("Slide",false))
            mapper.markBoundary(0,boundary,it->direction());
        else
            for (index_t d = 0; d!=domain.domainDim(); d++)
                mapper.markBoundary(0,boundary,d);

    }

    gsMatrix<index_t> int10, int20, int11, int21;
    for (typename std::vector<boundaryInterface>::const_iterator ifc = m_interfaces.cbegin();
                                                                 ifc!= m_interfaces.cend(); ++ifc)
    {
        GISMO_ENSURE(ifc->first().patch == 0 && ifc->second().patch == 0,"Interfaces are only supported for single patch domains");
        GISMO_ENSURE(ifc->first().direction() == ifc->second().direction(),"Interfaces are only supported for same direction");

        int10 = domain.basis().boundaryOffset(ifc->first().side(), 0);
        int20 = domain.basis().boundaryOffset(ifc->second().side(), 0);
        int11 = domain.basis().boundaryOffset(ifc->first().side(), 1);
        int21 = domain.basis().boundaryOffset(ifc->second().side(), 1);
        GISMO_ASSERT(int10.rows() == int20.rows() && int10.rows() == int11.rows() && int10.rows() == int21.rows(),"Boundary offsets have different sizes");
        index_t sz = int10.rows();
        for (index_t d = 0; d!=domain.domainDim(); d++)
        {
            if (d!=ifc->first().direction())
            {
                for ( index_t k=0; k<sz; ++k)
                {
                    mapper.matchDof( ifc->first().patch, int10(k,0), ifc->second().patch, int20(k,0), d );
                    mapper.matchDof( ifc->first().patch, int10(k,0), ifc->second().patch, int11(k,0), d );
                    mapper.matchDof( ifc->first().patch, int20(k,0), ifc->second().patch, int21(k,0), d );
                }
            }
        }
    }

    mapper.finalize();
}

template <class T>
void gsSquareDomain<T>::_initIndices(const gsGeometry<T> & domain, const gsDofMapper & mapper, std::vector<std::pair<index_t,index_t>> & indices) const
{
    indices.resize(mapper.freeSize());
    // std::vector<index_t> i(mapper.freeSize());
    // std::vector<index_t> j(mapper.freeSize());
    index_t ii;
    for (index_t k = 0; k!=domain.coefs().rows(); k++)
        for (index_t d = 0; d!=domain.domainDim(); d++)
        {
            ii = mapper.index(k,0,d);
            if (mapper.is_free_index(ii))
                indices[mapper.index(k,0,d)] = std::make_pair(k,d);
        }
}


namespace internal
{

/// @brief Get a gsSquareDomain from XML data
template<class T>
class gsXml< gsSquareDomain<T> >
{
private:
    gsXml() { }
public:
    // GSXML_COMMON_FUNCTIONS(gsSquareDomain<TMPLA2(T)>);
    static std::string tag () { return gsXml<gsGeometry<T>>::tag(); }
    static std::string type () { return gsXml<gsGeometry<T>>::type(); }

    // static gsGeometry<T> * get (gsXmlNode * node)
    // {
    //     return getGeometryFromXml< gsGeometry<T> >( node );
    // }

    static gsXmlNode * put (const gsSquareDomain<T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj.domain(),data);
    }
};

} // internal

} // namespace gismo

