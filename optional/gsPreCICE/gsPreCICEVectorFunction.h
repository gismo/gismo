/** @file gsPreCICEVectorFunction.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-...)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEFunction.h>

namespace gismo
{

/**
 * @brief      Class defining a gsFunction that reads from the precice::SolverInterface
 *
 * @tparam     T     Number format
 */
template <class T>
class gsPreCICEVectorFunction : public gsPreCICEFunction<T>
{
    using Base = gsPreCICEFunction<T>;
public:

    /// Shared pointer for gsPreCICEVectorFunction
    typedef memory::shared_ptr< gsPreCICEVectorFunction > Ptr;

    /// Unique pointer for gsPreCICEVectorFunction
    typedef memory::unique_ptr< gsPreCICEVectorFunction > uPtr;

    /// Default constructor
    gsPreCICEVectorFunction() { }

    /**
     * @brief      Constructs a new instance of the gsPreCICEVectorFunction
     *
     * @param      interface   The precice::SolverInterface (see \a gsPreCICE)
     * @param[in]  meshID      The ID of the mesh on which the data is located
     * @param[in]  dataID      The ID of the data
     * @param[in]  patches     The geometry
     * @param[in]  parametric  Specifies whether the data is defined on the parametric domain or not
     */
    gsPreCICEVectorFunction(        gsPreCICE<T> *    interface,
                        const index_t &         meshID,
                        const index_t &         dataID,
                        const gsMultiPatch<T> & patches,
                        const index_t &         targetDim=1,
                        const bool parametric = false)
    :
    Base(interface,meshID,dataID,patches,parametric)
    {
        m_targetDim = targetDim;
    }

    /// Constructs a function pointer
    static uPtr make(   const gsPreCICE<T> *    interface,
                        const index_t &         meshID,
                        const index_t &         dataID,
                        const gsMultiPatch<T> & patches,
                        const index_t &         targetDim=1,
                        const bool parametric = false)
    { return uPtr(new gsPreCICEVectorFunction(interface, meshID, dataID, patches, targetDim, parametric)); }

    GISMO_CLONE_FUNCTION(gsPreCICEVectorFunction)

    using Base::_getCoords;

    /// See \a gsFunction
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> coords;
        gsMatrix<T> tmp;
        _getCoords(u,coords);

        m_interface->readBlockVectorData(m_meshID,m_dataID,coords,tmp);
        result = tmp;
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsPreCICEVectorFunction, defined using gsPreCICE (precice::SolverInterface):\n";
        os << m_interface;
        return os;
    }

private:

    using Base::m_interface;
    using Base::m_meshID;
    using Base::m_dataID;
    using Base::m_targetDim;

};

}
