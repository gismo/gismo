/** @file gsAdaptiveRefUtils.h

    @brief Provides class for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once


#include <iostream>

namespace gismo
{


/**
 * @brief      This class provides a function that returns a constant error on each element
 * 
 * The elements of the provided basis are extracted and the errors from the error container
 * are associated to each element. The function can be passed to a Paraview export.
 *
 * @tparam     T     { description }
 */
template<typename T>
class gsElementErrorPlotter : public gsFunction<T>
{
public:
    /**
     * @brief      Constructs a function to plot the error of elements
     *
     * @param[in]  basis   The basis
     * @param[in]  errors  The errors per element
     */
    gsElementErrorPlotter(const gsBasis<T>& basis, const std::vector<T>& errors ) 
    : m_basis(basis),m_errors(errors)
    { }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
    {
        // Initialize domain element iterator -- using unknown 0
        res.setZero(1,u.cols());
        gsMatrix<T> supp(m_basis.domainDim(),2);
        for(index_t i=0; i<u.cols();++i)
        {
            int iter =0;
            // Start iteration over elements

            typename gsBasis<T>::domainIter domIt = m_basis.makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                bool flag = true;
                supp.col(0) = domIt->lowerCorner();
                supp.col(1) = domIt->upperCorner();

                // if (const gsComposedBasis<T> * cb = dynamic_cast<const gsComposedBasis<T> *>(&m_basis))
                // {
                    
                //     gsMatrix<T> pointsInv = gsMatrix<T>::Zero(supp.rows(),supp.cols());
                //     pointsInv.colwise() = domIt->centerPoint(); // Initialize with the centerpoint
                //     cb->composition()->invertPoints(supp,pointsInv,1e-4,true);
                //     GISMO_ASSERT(!(pointsInv.array()==std::numeric_limits<T>::infinity()).any(),"Point inversion failed");
                //     pointsInv.swap(supp);
                // }

                for(int d=0; d<domainDim();++d )
                {
                    if(supp.col(0)(d)> u(d,i) || u(d,i) > supp.col(1)(d)) // replace with .array().any()
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                     res(0,i) = m_errors.at(iter);
                     break;
                }
                iter++;
            }
        }
    }

    short_t domainDim() const { return m_basis.dim();}

private:
    const gsBasis<T>& m_basis;
    const std::vector<T>& m_errors;
};

} // namespace gismo