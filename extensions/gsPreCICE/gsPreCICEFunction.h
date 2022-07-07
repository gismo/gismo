/** @file gsPreCICEFunction.h

    @brief Provides declaration of ConstantFunction class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsPreCICE/gsPreCICE.h>

namespace gismo
{

/**
    @brief Class defining a globally constant function

    \tparam T value type

    \ingroup function
    \ingroup Core
*/

template <class T>
class gsPreCICEFunction : public gsFunction<T>
{
public:
    typedef gsGeometry<T> Base;

    /// Shared pointer for gsPreCICEFunction
    typedef memory::shared_ptr< gsPreCICEFunction > Ptr;

    /// Unique pointer for gsPreCICEFunction
    typedef memory::unique_ptr< gsPreCICEFunction > uPtr;

    gsPreCICEFunction() { }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^{\text{dim(val)}} \f$
    gsPreCICEFunction(        gsPreCICE<T> *    interface,
                        const index_t &         meshID,
                        const index_t &         dataID,
                        const gsMultiPatch<T> & patches,
                        const bool parametric = false)
    :
    m_interface(interface),
    m_meshID(meshID),
    m_dataID(dataID),
    m_patches(patches),
    m_parametric(parametric),
    m_patchID(0),
    m_domainDim(m_patches.domainDim())
    {
    }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^{\text{dim(val)}} \f$
    static uPtr make(   const gsPreCICE<T> *    interface,
                        const index_t &         meshID,
                        const index_t &         dataID,
                        const gsMultiPatch<T> & patches,
                        const bool parametric = false)
    { return uPtr(new gsPreCICEFunction(interface, meshID, dataID, patches, parametric)); }

    GISMO_CLONE_FUNCTION(gsPreCICEFunction)

    const gsPreCICEFunction<T> & piece(const index_t) const
    {
        return *this;
    }

    // Documentation in gsFunction class
    virtual short_t domainDim() const
    { return m_domainDim; }

    // Documentation in gsFunction class
    virtual short_t targetDim() const
    { return 1; }

    // Documentation in gsFunction class
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);

        gsMatrix<T> coords(m_patches.targetDim(),u.cols());
        if (m_parametric)
            m_patches.patch(m_patchID).eval_into(u,coords);
        else
            coords = u;

        m_interface->readBlockScalarData(m_meshID,m_dataID,coords,result);
    }

    // Documentation in gsFunction class
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    // Documentation in gsFunction class
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> > & result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice

        // GISMO_NO_IMPLEMENTATION;

        gsMatrix<T> tmp;
        this->eval_into(u,tmp);
        result.push_back(tmp);
    }

    // Documentation in gsFunction class
    virtual std::ostream &print(std::ostream &os) const
    {
        return os;
    }


private:

    gsPreCICE<T> * m_interface;
    gsMultiPatch<T> m_patches;
    bool m_parametric;
    index_t m_meshID, m_dataID;
    index_t m_patchID;
    index_t m_domainDim;

};

}
