/** @file gsLaplacePde.h

    @brief Describes a Laplace equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsPde/gsPde.h>
#include <gsCore/gsPiecewiseFunction.h>

namespace gismo
{

/** @brief
    The Laplace equation.

    This class describes a Laplace equation. This class exists mainly
    for compatibility reasons.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsLaplacePde : public gsPde<T>
{
protected:
    using gsPde<T>::m_domain;

public:
    gsLaplacePde( ) { }

    gsLaplacePde(const gsMultiPatch<T> & domain)
    {
        m_domain = domain;
        this->m_unknownDim.setOnes(1);
    }

    gsLaplacePde(const gsMultiPatch<T>         & domain,
                 const gsBoundaryConditions<T> & bc)
    : gsPde<T>(domain,bc)
    {
        this->m_unknownDim.setOnes(1);
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
	    os<<"Laplace equation  -\u0394u = 0\n";
	    return os; 
	}

}; // class gsLaplacePde

} // namespace gismo
