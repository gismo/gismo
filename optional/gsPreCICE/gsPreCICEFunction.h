/** @file gsPreCICEFunction.h

    @brief Provides declaration of ConstantFunction class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-...)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsPreCICE/gsPreCICE.h>

namespace gismo
{

/**
 * @brief      Class defining a gsFunction that reads from the precice::SolverInterface
 *
 * @tparam     T     Number format
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

    /// Default constructor
    gsPreCICEFunction() { }

    /**
     * @brief      Constructs a new instance of the gsPreCICEFunction
     *
     * @param      interface   The precice::SolverInterface (see \a gsPreCICE)
     * @param[in]  meshName      The ID of the mesh on which the data is located
     * @param[in]  dataName      The ID of the data
     * @param[in]  patches     The geometry
     * @param[in]  parametric  Specifies whether the data is defined on the parametric domain or not
     */
    gsPreCICEFunction(        gsPreCICE<T> *    interface,
                        const std::string &     meshName,
                        const std::string &     dataName,
                        const gsMultiPatch<T> & patches,
                        const index_t &         targetDim,
                        const bool parametric = false)
    :
    m_interface(interface),
    m_meshName(meshName),
    m_dataName(dataName),
    m_patches(patches),
    m_parametric(parametric),
    m_patchID(0),
    m_domainDim(m_patches.domainDim()),
    m_targetDim(targetDim)
    {
    }

    /// Constructs a function pointer
    static uPtr make(   const gsPreCICE<T> *    interface,
                        const std::string &     meshName,
                        const std::string &     dataName,
                        const gsMultiPatch<T> & patches,
                        const index_t &         targetDim,
                        const bool parametric = false)
    { return uPtr(new gsPreCICEFunction(interface, meshName, dataName, patches, targetDim, parametric)); }

    GISMO_CLONE_FUNCTION(gsPreCICEFunction)

    /// Access a piece
    const gsPreCICEFunction<T> & piece(const index_t) const
    {
        return *this;
    }

    /// See \a gsFunction
    virtual short_t domainDim() const
    { return m_domainDim; }

    /// Gives the targetDomain, currently only scalar functions (todo)
    virtual short_t targetDim() const
    { return m_targetDim; }

    /// See \a gsFunction
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> coords;
        this->_getCoords(u,coords);
        m_interface->readData(m_meshName,m_dataName,coords,result);
    }

    /// See \a gsFunction
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \a gsFunction
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \a gsFunction
    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> > & result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        gsMatrix<T> tmp;
        this->eval_into(u,tmp);
        result.push_back(tmp);
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsPreCICEFunction, defined using gsPreCICE (precice::SolverInterface):\n";
        os << m_interface;
        return os;
    }

protected:
    void _getCoords(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);

        // Does not work
        // if (m_parametric)
        //     m_interface->readBlockScalarData(m_meshName,m_dataName,m_patches.patch(m_patchID).eval(u),result);
        // if (m_parametric)
        //     m_interface->readBlockScalarData(m_meshName,m_dataName,u,result);

        result.resize(m_patches.targetDim(),u.cols());
        if (m_parametric)
            m_patches.patch(m_patchID).eval_into(u,result);
        else
            result = u;
    }

protected:

    gsPreCICE<T> * m_interface;
    std::string m_meshName;
    std::string m_dataName;
    gsMultiPatch<T> m_patches;
    bool m_parametric;
    index_t m_patchID;
    index_t m_domainDim, m_targetDim;

};

}
