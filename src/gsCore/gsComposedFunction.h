/** @file gsComposedFunction.h

    @brief Provides the implementation of a function composed by another function

    Given a parametric domain (xi,eta), a composition (u,v) = C(xi,eta),
    the function f(xi,eta) is evaluated as f(C(xi,eta)).
    The derivatives are defined with respect to xi, eta

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst
        S. Imperatore
*/

#pragma once

//! [Include namespace]
#include <gsCore/gsFunction.h>

namespace gismo
{

template<class T>
class gsComposedFunction : public gsFunction<T>
{
public:

    typedef memory::shared_ptr< gsComposedFunction > Ptr;
    typedef memory::unique_ptr< gsComposedFunction > uPtr;

    typedef gsFunction<T>   CompositionT;
    typedef gsFunction<T>   FunctionT;

    GISMO_CLONE_FUNCTION(gsComposedFunction)

public:

    /// @brief Empty constructor
    gsComposedFunction();

    /**
     * @brief Construct a composed function from pointers
     *
     * @param[in] composition   the composition
     * @param[in] basis         the basis to be composed
     */
    gsComposedFunction( const CompositionT  * composition,
                        const FunctionT * function);

    /**
     * @brief Construct a composed function from references
     *
     * @param[in] composition   the composition
     * @param[in] basis         the basis to be composed
     */
    gsComposedFunction( const CompositionT  & composition,
                        const FunctionT & function);

    /**
     * @brief Construct a composed function from unique pointers
     *
     * @param[in] composition   the composition
     * @param[in] basis         the basis to be composed
     */
    gsComposedFunction( typename CompositionT::Ptr composition,
                        typename FunctionT::Ptr function);

    /// Return the composition
    const CompositionT & composition() const;
    /// Return the function
    const FunctionT & function() const;

    /// See \ref gsFunction for more documentation
    short_t domainDim() const override;

    /// See \ref gsFunction for more documentation
    short_t targetDim() const override;

    /// See \ref gsFunction for more documentation
    gsMatrix<T> support() const override;

    // void evalAllDers_into(const gsMatrix<T> & u, int n,
    //                         std::vector<gsMatrix<T> >& result,
    //                         bool sameElement) const
    // {
    //     gsMatrix<T> coords = m_composition->eval(u);
    //     this->_applyBounds(coords);
    //     m_basis->evalAllDers_into(coords,n,result,sameElement);
    // }



    /// See \ref gsFunction for more documentation
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    /// See \ref gsFunction for more documentation
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    /// See \ref gsFunction for more documentation
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    std::ostream &print(std::ostream &os) const override;

protected:

    typename CompositionT::Ptr   m_composition;
    typename FunctionT::Ptr  m_function;

};


/*
//  Implementation for arbitrary number of functions

template<class T>
class gsComposedFunction : public gsFunction<T>
{
public:

    // gsComposedFunction(const gsFunction<T> & composition, const gsFunction<T> & function)
    // :
    // gsComposedFunction(std::vector<gsFunction<T>*>(composition.clone().release(),
    //                                                function.clone().release()))
    // {}

    gsComposedFunction(const gsFunction<T> * composition, const gsFunction<T> * function)
    {
        gsDebugVar(composition);
        gsDebugVar(*composition);
        m_functions.push_back(memory::make_shared_not_owned(composition));
        m_functions.push_back(memory::make_shared_not_owned(function));
        for (size_t l = 0; l!=m_functions.size()-1; l++)
            GISMO_ENSURE(m_functions[l+1]->domainDim()==m_functions[l]->targetDim(),
                "Domain dimension of function "<<l+1<<
                " should be equal to the target dimension of function "<<l<<
                ", but functions[l+1]->domainDim() = "<<m_functions[l+1]->domainDim()<<
                " and functions[l]->targetDim() = )"<<m_functions[l]->targetDim());
    }

    gsComposedFunction()
    {}

    // gsComposedFunction(std::vector<const gsFunction<T> *> functions)
    // :
    // m_functions(functions)
    // {
    //     for (size_t l = 0; l!=m_functions.size()-1; l++)
    //         GISMO_ENSURE(m_functions[l+1]->domainDim()==m_functions[l]->targetDim(),
    //             "Domain dimension of function "<<l+1<<
    //             " should be equal to the target dimension of function "<<l<<
    //             ", but functions[l+1]->domainDim() = "<<m_functions[l+1]->domainDim()<<
    //             " and functions[l]->targetDim() = )"<<m_functions[l]->targetDim());
    // }

    // ~gsComposedFunction()
    // {
    // }

public:

    // /// Move constructor
    // gsComposedFunction( gsComposedFunction&& other )
    // {

    // }

    // /// Move assignment operator
    // gsComposedFunction& operator= ( gsComposedFunction&& other )
    // {
    //     freeAll(m_functions);
    //     m_functions = give(other.m_functions);
    //     return *this;
    // }

public:

    short_t domainDim() const { return m_functions.front()->domainDim(); }
    short_t targetDim() const { return m_functions.back()->targetDim(); }

    gsMatrix<T> support() const { return m_functions.front()->support(); }

    // void evalAllDers_into(const gsMatrix<T> & u, int n,
    //                         std::vector<gsMatrix<T> >& result,
    //                         bool sameElement) const
    // {
    //     gsMatrix<T> coords = m_composition->eval(u);
    //     this->_applyBounds(coords);
    //     m_basis->evalAllDers_into(coords,n,result,sameElement);
    // }



    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        gsDebugVar(m_functions.front());
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

        for (size_t l = 1; l!=m_functions.size(); l++)
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

    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    std::ostream &print(std::ostream &os) const
    {
        os <<"Composite function:\n";
        for (size_t f = 0; f!=m_functions.size(); f++)
        {
            os << "* Function "<<f
               << " ( R^" << m_functions[f]->domainDim() << " --> R^" << m_functions[f]->targetDim() << "):\n"
               << *m_functions[f]<<"\n"
               << "(address: "<<m_functions[f]<<")\n";
        }
        return os;
    }

    index_t numCompositions() const { return m_functions.size()-1; };

    // const gsFunction<T> & composition(const index_t i) const { return *m_functions[i]; }
    //       gsFunction<T> & composition(const index_t i)       { return *m_functions[i]; }

    const gsFunction<T> * composition(const index_t i) const { return  m_functions[i]; }

protected:
    std::vector<typename gsFunction<T>::Ptr> m_functions;
};
*/
} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsComposedFunction.hpp)
#endif
