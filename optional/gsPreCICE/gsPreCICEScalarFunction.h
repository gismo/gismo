/** @file gsPreCICEScalarFunction.h

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

namespace gismo
{

/**
 * @brief      Class defining a gsFunction that reads from the precice::SolverInterface
 *
 * @tparam     T     Number format
 */
template <class T>
class gsPreCICEScalarFunction : public gsPreCICEFunction<T>
{
    using Base = gsPreCICEFunction<T>;
public:

    /// Shared pointer for gsPreCICEScalarFunction
    typedef memory::shared_ptr< gsPreCICEScalarFunction > Ptr;

    /// Unique pointer for gsPreCICEScalarFunction
    typedef memory::unique_ptr< gsPreCICEScalarFunction > uPtr;

    /// Default constructor
    gsPreCICEScalarFunction() { }

    /**
     * @brief      Constructs a new instance of the gsPreCICEScalarFunction
     *
     * @param      interface   The precice::SolverInterface (see \a gsPreCICE)
     * @param[in]  meshID      The ID of the mesh on which the data is located
     * @param[in]  dataID      The ID of the data
     * @param[in]  patches     The geometry
     * @param[in]  parametric  Specifies whether the data is defined on the parametric domain or not
     */
    gsPreCICEScalarFunction(        gsPreCICE<T> *    interface,
                        const index_t &         meshID,
                        const index_t &         dataID,
                        const gsMultiPatch<T> & patches,
                        const bool parametric = false)
    :
    Base(interface,meshID,dataID,patches,parametric){}

    /// Constructs a function pointer
    static uPtr make(   const gsPreCICE<T> *    interface,
                        const index_t &         meshID,
                        const index_t &         dataID,
                        const gsMultiPatch<T> & patches,
                        const bool parametric = false)
    { return uPtr(new gsPreCICEScalarFunction(interface, meshID, dataID, patches, 1, parametric)); }

    GISMO_CLONE_FUNCTION(gsPreCICEScalarFunction)

    /// Gives the targetDomain, currently only scalar functions (todo)
    virtual short_t targetDim() const
    { return 1; }

    using Base::_getCoords;

    /// See \a gsFunction
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> coords;
        _getCoords(u,coords);
        m_interface->readBlockScalarData(m_meshID,m_dataID,coords,result);
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsPreCICEScalarFunction, defined using gsPreCICE (precice::SolverInterface):\n";
        os << m_interface;
        return os;
    }

private:

    using Base::m_interface;
    using Base::m_meshID;
    using Base::m_dataID;

};

}
