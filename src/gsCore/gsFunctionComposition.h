/** @file gsFunctionComposition.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

//! [Include namespace]
#include <gsCore/gsFunction.h>

#  define return_t  auto

namespace gismo
{

template<class T>
class gsFunctionComposition : public gsFunction<T>
{
public:
    gsFunctionComposition(const std::vector<gsFunction<T> *> functions)
    :
    m_functions(functions)
    {
        for (size_t l = 0; l!=m_functions.size()-1; l++)
            GISMO_ENSURE(m_functions[l+1]->domainDim()==m_functions[l]->targetDim(),
                "Domain dimension of function "<<l+1<<
                " should be equal to the target dimension of function "<<l<<
                ", but functions[l+1]->domainDim() = "<<m_functions[l+1]->domainDim()<<
                " and functions[l]->targetDim() = )"<<m_functions[l]->targetDim());
    }

    short_t domainDim() const { return m_functions.front()->domainDim(); }
    short_t targetDim() const { return m_functions.back()->targetDim(); }

    gsMatrix<T> support() const { return m_functions.front()->support(); }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsMatrix<T> coord = u;
        for (size_t l = 0; l!=m_functions.size(); l++)
        {
            m_functions[l]->eval_into(coord,result);
            coord = result;
        }
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        index_t domainDim, targetDim;
        gsMatrix<T> coord = u, newcoord, deriv, tmp, tmpresult;

        m_functions[0]->deriv_into(coord,tmpresult);
        domainDim = m_functions[0]->domainDim();
        targetDim = m_functions[0]->targetDim();

        for (index_t l = 1; l!=m_functions.size(); l++)
        {
            // Compute the new coord for the next function
            m_functions[l-1]->eval_into(coord,newcoord);
            coord = newcoord;

            // evaluate the derivatives on coord
            // The derivatives are structured as follows (each col is a point of u):
            // [[df1/dx1
            //   df1/dx2
            //   ...
            //   df2/dx1
            //   df2/dx2
            //   ...
            //   dfn/dx1
            //   dfn/dx2
            //   ...]]
            m_functions[l]->deriv_into(coord,deriv);
            tmp.resize(m_functions[l]->targetDim()*domainDim,u.cols());
            for (index_t k = 0; k!=u.cols(); k++)
            {
                gsAsMatrix<T,Dynamic,Dynamic> resultMat = tmpresult.reshapeCol(k,domainDim,targetDim);
                gsAsMatrix<T,Dynamic,Dynamic> derivMat = deriv.reshapeCol(k,m_functions[l]->domainDim(),m_functions[l]->targetDim());
                // The product has size:
                // (domainDim x targetDim) x (m_functions[l]->domainDim(),m_functions[l]->targetDim())
                //  =
                // (domainDim x m_functions[l]->targetDim())
                gsAsMatrix<T,Dynamic,Dynamic> tmpMat = tmp.reshapeCol(k,domainDim,m_functions[l]->targetDim());
                tmpMat = resultMat*derivMat;

            }
            targetDim = m_functions[l]->targetDim();
            tmpresult = tmp;
        }
        result = tmpresult;
    }

    std::ostream &print(std::ostream &os) const
    {
        os <<"Composite function:\n";
        for (size_t f = 0; f!=m_functions.size(); f++)
        {
            os << "* Function "<<f
               << " ( R^" << m_functions[f]->domainDim() << " --> R^" << m_functions[f]->targetDim() << "):\n"
               << *m_functions[f]<<"\n";
        }
        return os;
    }

protected:
    const std::vector<gsFunction<T> *> m_functions;
};
}
