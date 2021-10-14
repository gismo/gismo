/** @file gsContainerBasisBase.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

/*
    TO DO
*/

#pragma once

#include<gsIO/gsOptionList.h>
#include<gsUnstructuredSplines/gsContainerBasis.h>

namespace gismo
{

    template<short_t d,class T>
    class gsContainerBasisBase
    {

    private:
        typedef gsContainerBasis<d,T> basis;
        typedef typename std::vector<basis> basisContainer;

        /// Shared pointer for gsContainerBasisBase
        typedef memory::shared_ptr< gsContainerBasisBase > Ptr;

        /// Unique pointer for gsContainerBasisBase
        typedef memory::unique_ptr< gsContainerBasisBase > uPtr;

    public:

        virtual ~gsContainerBasisBase() {};

        gsContainerBasisBase(gsMultiPatch<T> & mp, gsMultiBasis<T> & mb)
                : m_patches(mp), m_multiBasis(mb)
        {

        }

    public:

        virtual gsOptionList & options() { return m_options; }
        virtual void defaultOptions() { };

        virtual void setOptions(gsOptionList opt) {m_options.update(opt, gsOptionList::addIfUnknown); };

        virtual void init() = 0;
        virtual void compute() = 0;
        

        void getMultiBasis(gsMultiBasis<T> & multiBasis_result)
        {
            multiBasis_result.clear();

            std::vector<gsBasis<T> *> basis_temp = std::vector<gsBasis<T> *>(m_patches.nPatches());
            for (size_t np = 0; np < m_patches.nPatches(); np++) {
                gsContainerBasis<2, real_t>::uPtr basis = gsContainerBasis<d, T>::make(m_bases[np]);
                basis_temp[np] = static_cast<gsBasis<> *>(basis.release());
            }

            multiBasis_result = gsMultiBasis<>(basis_temp, m_patches.topology());
        };

        gsSparseMatrix<T> & getSystem() { return m_matrix; };
        void setSystem(gsSparseMatrix<T> & system) { m_matrix = system; };

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
        basisContainer m_bases;

        /// System matrix
        gsSparseMatrix<T> m_matrix;
    };

}