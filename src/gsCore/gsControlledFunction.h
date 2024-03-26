/** @file gsControlledFunction.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
*/

//! [Include namespace]
#include <gsCore/gsFunction.h>

#  define return_t  auto

namespace gismo
{

template<class T>
class gsControlledFunction : public gsFunction<T>
{
public:
    gsControlledFunction()
    {}

    virtual short_t domainDim() const = 0;

    const gsVector<T> & parameters() const { return m_parameters; };
          gsVector<T> & parameters()       { return m_parameters; };

    size_t parameterSize() const {return m_parameters.size();}

    T & parameter(index_t i) {return m_parameters[i];}

    // void precomputeDerivatives();

    // void parameterDerivative(index_t i)
    // {
    //     return m_derivates[i];
    // }


protected:
    gsVector<T> m_parameters;

};

}// End namespace