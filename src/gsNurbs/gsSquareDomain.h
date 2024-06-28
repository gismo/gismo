/** @file gsMultiBasis.h

    @brief Provides declaration of MultiBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsDofMapper.h>

namespace gismo
{

template <short_t DIM, class T>
class gsSquareDomain : public gsFunction<T>
{
    using Base = gsFunction<T> ;

public:
    // default constructor
    // gsSquareDomain()

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  domain  The domain
     */
    gsSquareDomain(const gsTensorBSpline<DIM,T> & domain)
    {
        m_domain = domain;
        this->_initMapper(m_domain,m_mapper);
        this->_initIndices(m_domain,m_mapper,m_indices);
    }

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  basis  The basis
     * @param[in]  domain  The domain
     */
    gsSquareDomain(const gsTensorBSplineBasis<DIM,T> & basis)
    {
        gsMatrix<T> coefs = basis.anchors();
        coefs.transposeInPlace();
        m_domain = gsTensorBSpline<DIM,T>(basis,coefs);
        this->_initMapper(m_domain,m_mapper);
        this->_initIndices(m_domain,m_mapper,m_indices);
    }


    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  numElevation  The number elevation
     * @param[in]  numRefine     The number refine
     */
    gsSquareDomain(index_t numElevation = 0, index_t numRefine = 0)
    {
        m_domain = *gsNurbsCreator<T>::BSplineSquare();
        m_domain.degreeElevate(numElevation);
        index_t numKts = pow(2, numRefine) - 1;
        m_domain.uniformRefine(numKts);

        this->_initMapper(m_domain,m_mapper);
        this->_initIndices(m_domain,m_mapper,m_indices);
    }

    const gsTensorBSpline<DIM,T> & domain() const
    {
        return m_domain;
    }

    const gsDofMapper & mapper() const
    {
        return m_mapper;
    }

    gsMatrix<T> support() const override
    {
        return m_domain.support();
    }

    short_t maxDegree() const
    {
        return m_domain.basis().maxDegree();
    }

    short_t domainDim() const override
    {
        return m_domain.domainDim();
    }

    short_t targetDim() const override
    {
        return m_domain.domainDim();
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        m_domain.eval_into(u,result);
    }

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        m_domain.deriv_into(u,result);
    }

    /// Returns the controls of the function
    // DO NOT WORK YET
        //   gsAsConstVector<T> controls() const override { return gsAsConstVector<T>(m_domain.coefs(),m_indices); };
        //   gsAsVector<T>      controls()       override { return gsAsVector<T>     (m_domain.coefs(),m_indices); };  

    /// Returns the \a i th control of the function
    // const typename gsMatrix<T>::CoeffReturnType & control(index_t i) const override { return gsAsConstVector<T>(m_parameters.data(),m_parameters.size())(i);}
    const T & control(index_t i) const override { return m_domain.coefs()(m_indices[i].first,m_indices[i].second);}
          T & control(index_t i)       override { return m_domain.coefs()(m_indices[i].first,m_indices[i].second);}


    // const gsVector<T> & parameters() const { return m_parameters; };
    //       gsVector<T> & parameters()       { return m_parameters; };

    /// Returns the number of controls of the function
    size_t nControls() const override
    {
        return m_mapper.freeSize();
    }

    /// Returns the control derivative
    virtual void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override
    {
        gsMatrix<T> tmp;

        result.resize(targetDim()*nControls(), points.cols());
        result.setZero();
        for (index_t p = 0; p!=points.cols(); p++)
        {
            gsAsMatrix<T> res = result.reshapeCol(p,nControls(),targetDim());
            for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
                for (index_t d = 0; d!=m_domain.targetDim(); d++)
                    if (m_mapper.is_free(k,0,d))
                    {
                        m_domain.basis().evalSingle_into(k,points.col(p),tmp); // evaluate basis function k
                        res(m_mapper.index(k,0,d),d) = tmp(0,0); // tmp is a single value (1 point, 1 basis function)
                    }
        }
    }

private:
    void _initMapper(const gsTensorBSpline<DIM,T> & domain, gsDofMapper & mapper) const
    {
        // Mapper storing control points
        mapper = gsDofMapper(domain.basis(),domain.targetDim());

        gsBoxTopology topology(DIM,1);
        topology.addAutoBoundaries();
        // gsMatrix<index_t> boundary = domain.basis().allBoundary();
        gsDebugVar(topology);
        for (typename gsBoxTopology::biterator it = topology.bBegin(); it != topology.bEnd(); ++it)
        {
            gsMatrix<index_t> boundary = domain.basis().boundary(*it);
            gsDebugVar(boundary);
            // for (index_t a = 0; a!=boundary.rows(); a++)
                // for (index_t d = 0; d!=domain.targetDim(); d++)
                mapper.markBoundary(0,boundary,it->direction());
        }
        mapper.finalize();
    }

    void _initIndices(const gsTensorBSpline<DIM,T> & domain, const gsDofMapper & mapper, std::vector<std::pair<index_t,index_t>> & indices) const
    {
        indices.resize(mapper.freeSize());
        // std::vector<index_t> i(mapper.freeSize());
        // std::vector<index_t> j(mapper.freeSize());
        for (index_t k = 0; k!=domain.coefs().rows(); k++)
            for (index_t d = 0; d!=domain.targetDim(); d++)
                if (mapper.is_free(k,0,d))
                    indices[mapper.index(k,0,d)] = std::make_pair(k,d);
    }

protected:
    gsTensorBSpline<DIM,T> m_domain;
    gsDofMapper m_mapper;
    std::vector<std::pair<index_t,index_t>> m_indices;
};

} // namespace gismo
