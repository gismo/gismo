/** @file gsC1SplineBase.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

/*
    TO DO
*/

#pragma once

#include<gsIO/gsOptionList.h>
#include<gsC1Basis/gsC1Basis.h>

namespace gismo
{

template<short_t d,class T>
class gsC1SplineBase
{

private:
    typedef gsC1Basis<d,T> Basis;
    typedef typename std::vector<Basis> C1BasisContainer;

    /// Shared pointer for gsC1SplineBase
    typedef memory::shared_ptr< gsC1SplineBase > Ptr;

    /// Unique pointer for gsC1SplineBase
    typedef memory::unique_ptr< gsC1SplineBase > uPtr;

public:

    virtual ~gsC1SplineBase() {};

    gsC1SplineBase(gsMultiPatch<T> & mp, gsMultiBasis<T> & mb)
    : m_patches(mp), m_multiBasis(mb)
    {

    }

public:

    virtual gsOptionList & options() { return m_options; }
    virtual void defaultOptions() { };

    virtual void setOptions(gsOptionList opt) {m_options.update(opt, gsOptionList::addIfUnknown); };

    virtual void init() = 0;
    virtual void compute() = 0;

    void uniformRefine()
    {
        index_t p = m_multiBasis.minCwiseDegree();
        index_t r = m_options.getInt("discreteRegularity");

        m_multiBasis.uniformRefine(1,p-r);
    }

    void getMultiBasis(gsMultiBasis<T> & multiBasis_result)
    {
        multiBasis_result.clear();

        std::vector<gsBasis<T> *> basis_temp = std::vector<gsBasis<T> *>(m_patches.nPatches());
        for (size_t np = 0; np < m_patches.nPatches(); np++) {
            gsC1Basis<2, real_t>::uPtr basis = gsC1Basis<d, T>::make(m_bases[np]);
            basis_temp[np] = static_cast<gsBasis<> *>(basis.release());
        }

        multiBasis_result = gsMultiBasis<>(basis_temp, m_patches.topology());
    };

    gsSparseMatrix<T> & getSystem() { return m_matrix; };
    void setSystem(gsSparseMatrix<T> & system) { m_matrix = system; };

    T getMinMeshSize()
    {
        T meshSize = 1.0;
        for (size_t np = 0; np < m_patches.nPatches(); np++)
            if (m_multiBasis.basis(np).getMinCellLength() < meshSize)
                meshSize = m_multiBasis.basis(np).getMinCellLength();
        return meshSize;
    }

private:
    // Helper functions

protected:
    // Data members

    // Put here the members of the shared functions
    /// Multipatch
    gsMultiPatch<T> & m_patches;
    gsMultiBasis<T> m_multiBasis;

    /// Optionlist
    gsOptionList m_options;

    /// C1 Basis
    C1BasisContainer m_bases;

    /// System matrix
    gsSparseMatrix<T> m_matrix;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1SplineBase.hpp)
#endif
