/** @file gsComposedFunction.h

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

namespace gismo
{

template<class T>
class gsComposedFunction : public gsFunction<T>
{
public:
    typedef gsBasis<T>      BasisT;
    typedef gsFunction<T>   CompositionT;

    GISMO_OVERRIDE_CLONE_FUNCTION(gsComposedFunction)

public:

    gsComposedFunction(const gsFunction<T> & composition, const gsFunction<T> & function)
    :
    m_composition(&composition),
    m_function(&function)
    {
        GISMO_ENSURE(m_function->domainDim()==m_composition->targetDim(),
            "Domain dimension of the function "<<
            " should be equal to the target dimension of the composition "<<
            ", but basis.domainDim() = "<<m_function->domainDim()<<
            " and composition.targetDim() = )"<<m_composition->targetDim());
    }

    gsComposedFunction()
    {}

    short_t domainDim() const { return m_composition->domainDim(); }
    short_t targetDim() const { return m_function->targetDim(); }

    gsMatrix<T> support() const { return m_composition->support(); }

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
        gsMatrix<T> coords = m_composition->eval(u);
        this->_applyBounds(coords);
        m_function->eval_into(coords,result);
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        index_t domainDim, targetDim;
        domainDim = m_composition->domainDim();
        targetDim = m_composition->targetDim();

        gsFuncData<T> fd(NEED_VALUE | NEED_DERIV);
        m_composition->compute(u,fd);

        gsMatrix<T> coord, deriv, tmp, compderiv;
        coord = fd.values[0];
        compderiv = fd.values[1];

        this->_applyBounds(coord);
        m_function->deriv_into(coord,deriv);
        result.resize(m_function->targetDim()*domainDim,u.cols());
        for (index_t k = 0; k!=u.cols(); k++)
        {
            gsAsMatrix<T,Dynamic,Dynamic> compderivMat = compderiv.reshapeCol(k,domainDim,targetDim);
            gsAsMatrix<T,Dynamic,Dynamic> derivMat = deriv.reshapeCol(k,m_function->domainDim(),m_function->targetDim());
            // The product has size:
            // (domainDim x targetDim) x (m_function->domainDim(),m_function->targetDim())
            //  =
            // (domainDim x m_function->targetDim())
            gsAsMatrix<T,Dynamic,Dynamic> resultMat = result.reshapeCol(k,domainDim,m_function->targetDim());
            resultMat = compderivMat*derivMat;
        }
    }

    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    std::ostream &print(std::ostream &os) const
    {
        os <<"Composite function:\n";
        os << "* Composition ( R^" << m_composition->domainDim() << " --> R^" << m_composition->targetDim() << "):\n"
            << *m_composition<<"\n"
            << "(address: "<<m_composition<<")\n";
        os << "* Function ( R^" << m_function->domainDim() << " --> R^" << m_function->targetDim() << "):\n"
            << *m_function<<"\n"
            << "(address: "<<m_function<<")\n";
        return os;
    }

private:
    void _applyBounds(gsMatrix<T> & coords) const
    {
        // for (index_t k=0; k!=coords.cols(); k++)
        // {
        //     gsDebugVar(coords.col(k));
        //     gsDebugVar(m_composition->support().col(0));
        //     gsDebugVar(m_composition->support().col(1));
        //     coords.col(k) = coords.col(k).cwiseMax(m_composition->support().col(0));
        //     coords.col(k) = coords.col(k).cwiseMin(m_composition->support().col(1));
        // }
    }

protected:

    const CompositionT * m_composition;
    const gsFunction<T>* m_function;

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
}
