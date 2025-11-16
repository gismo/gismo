/** @file gsFeElement.h

    @brief Defines an element as an expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

#include <gsExpressions/gsFeSpace.h>

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for a finite element solution, representing a function given by a vector of
 *        coefficients in a gsFeSpace.
 *
 *        Typically it used for accessing the solution of a boundary-value
 *        problem.
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class gsFeSolution : public _expr<gsFeSolution<T> >
{
protected:
    const gsFeSpace<T> _u;
    gsMatrix<T> * _Sv; ///< Pointer to a coefficient vector
    bool m_isAcross; ///< true when this expression is evaluated across an interface

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    bool isAcross() const { return m_isAcross; }

    gsFeSolution right() const
    {
        gsFeSolution ac(*this);
        ac.m_isAcross = true;
        return ac;
    }

    gsFeSolution left() const { return gsFeSolution(*this); }

    explicit gsFeSolution(const gsFeSpace<T> & u) : _u(u), _Sv(NULL) { }

    gsFeSolution(const gsFeSpace<T> & u, gsMatrix<T> & Sv) : _u(u), _Sv(&Sv) { }

    const gsFeSpace<T> & space() const {return _u;};

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(check(), "Invalid state in gsFeSolution");
        const gsDofMapper & map = _u.mapper();
        auto & act = _u.data().actives.col(1 == _u.data().actives.cols() ? 0:k );
        res.setZero(_u.dim(), 1);
        for (index_t c = 0; c!=_u.dim(); c++) // for all components
        {
            for (index_t i = 0; i!=_u.data().actives.rows(); ++i)
            {
                const index_t ii = map.index( act[i], _u.data().patchId, c);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.at(c) += _Sv->at(ii) * _u.data().values[0](i,k);
                else
                    res.at(c) += _u.data().values[0](i,k) *
                        _u.fixedPart().at( map.global_to_bindex(ii) );
            }
        }
        return res;
    }

    // Performs validity checks for the solution object
    bool check() const
    {
        if ( _Sv->size()!=_u.mapper().freeSize() )
        {
            gsWarn<< "The solution vector has wrong dimensions: "
                  <<_Sv->size()<<" != "<<_u.mapper().freeSize() <<". ";
            return false;
        }
        if((size_t)_u.source().size()*dim()!=_u.mapper().mapSize())
        {
            gsWarn<< "The solution space is inconsistent: "
                  <<_u.source().size()*dim()<<" != "<<_u.mapper().mapSize()<<". ";
            return false;
        }
        return true;
    }

    //template<class U>
    //linearComb(U & ie){ sum up ie[_u] times the _Sv  }
    // ie.eval(k), _u.data().actives(), fixedPart() - see lapl_expr

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t rows() const {return _u.dim(); }

    static index_t cols() {return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_VALUE | NEED_ACTIVE;
    }

    void print(std::ostream &os) const { os << "s"; }

public:
    index_t dim() const { return _u.dim();}

    index_t parDim() const
    { return _u.source().domainDim(); }

    //gsDofMapper & mapper() {return _u.mapper();}
    const gsDofMapper & mapper() const {return _u.mapper();}

    inline const gsMatrix<T> & fixedPart() const {return _u.fixedPart();}
    gsMatrix<T> & fixedPart() {return _u.fixedPart();}

    //gsFuncData<T> & data() {return _u.data();}
    const gsFuncData<T> & data() const {return _u.data();}

    void setSolutionVector(gsMatrix<T>& solVector)
    { _Sv = & solVector; }

    /// @brief Sets all coefficients of the solution vector that belong to
    ///    patch \a p , and refer to the specified \a component equal to \a value.
    /// @param component The index of the component to be set.
    /// @param value The value that the coefficients will be set to.
    /// @param patch The index of the patch whose coefficients will be set. By default all patches are affected,
    void setComponent(index_t component, real_t value, index_t patch=-1)
    {
        gsMatrix<T> & solVector = *_Sv;
        const gsDofMapper & mapper = _u.mapper();

        index_t patchStart, patchEnd;
        if (patch==-1){
            patchStart = 0;
            patchEnd   = _u.mapper().numPatches();
        }
        else{
            patchStart = patch;
            patchEnd   = patch + 1;
        }

        for (index_t p=patchStart; p!=patchEnd; ++p)
        {
            for (size_t i = 0; i != mapper.patchSize(p, component); ++i)
            {
                const index_t ii = mapper.index(i, p, component);
                if ( mapper.is_free_index(ii) ) // DoF value is in the solVector
                    solVector.at(ii) = value;
            }
        }
    }

    const gsMatrix<T> & coefs() const { return *_Sv; }
    //gsMatrix<T> & coefs() { return *_Sv; } // wd4702 ?

    //const gsMatrix<T> & coefs(component, patch) const { return *_Sv; }

    /// val: perturbation value, j: global index, p: patch
    void perturbLocal(T val, index_t j, index_t p = 0)
    {
        GISMO_UNUSED(p);
        // GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        //if (_u.mapper().is_free_index(j) )
        //{
            GISMO_ASSERT(j<_Sv->size(), "Solution vector is not initialized/allocated, sz="<<_Sv->size() );
            _Sv->at(j) += val;
            //}
        //else
        //    _u.fixedPart().at( _u.mapper().global_to_bindex(j) ) += val;
    }

    /// Extract the coefficients of piece \a p
    void extract(gsMatrix<T> & result, const index_t p = 0) const
    { _u.getCoeffs(*_Sv, result, p); }

    /// Extracts ALL the coefficients in a solution vector; including
    /// coupled and boundary DoFs
    void extractFull(gsMatrix<T> & result) const
    {
        index_t offset;
        const index_t dim = _u.dim();
        const size_t totalSz = _u.mapper().mapSize();
        result.resize(totalSz, 1);
        for (size_t p=0; p!=_u.mapper().numPatches(); ++p)
        {
            offset = _u.mapper().offset(p);
            // Reconstruct solution coefficients on patch p

            for (index_t c = 0; c!=dim; c++) // for all components
            {
                const index_t sz  = _u.mapper().patchSize(p,c);

                // loop over all basis functions (even the eliminated ones)
                for (index_t i = 0; i < sz; ++i)
                {
                    //gsDebugVar(i);
                    const int ii = _u.mapper().index(i, p, c);
                    //gsDebugVar(ii);
                    if ( _u.mapper().is_free_index(ii) ) // DoF value is in the solVector
                    {
                        result(i+offset,0) = _Sv->at(ii);
                    }
                    else // eliminated DoF: fill with Dirichlet data
                        result(i+offset,0) =  _u.fixedPart().at( _u.mapper().global_to_bindex(ii) );
                }
                offset += sz;
            }
        }
    }

    /// Extract this variable as a multipatch object
    void extract(gsMultiPatch<T> & result) const
    {
        result.clear();

        if( const gsMultiBasis<T>* basis = dynamic_cast<const gsMultiBasis<T>* >(&_u.source()) )
            for (size_t i = 0; i != basis->nBases(); ++i)
            {
                memory::unique_ptr<gsGeometry<T> > p(this->extractPiece(i));
                result.addPatch(*p);
            }
    }

    /// Extract this variable as a gsMappedSpline object
    void extract(gsMappedSpline<2,T> & result) const
    {
        if( const gsMappedBasis<2,T>* basis = dynamic_cast<const gsMappedBasis<2,T>* >(&_u.source()) )
        {
            gsMatrix<T> coefs;
            this->extractFull(coefs);
            coefs.resize(coefs.rows()/_u.dim(),_u.dim());
            result.init(*basis,coefs);
        }
    }

    /// Extract the piece \a p as a gsGeometry pointer
    memory::unique_ptr<gsGeometry<T> > extractPiece(const index_t p) const
    {
        if ( const gsBasis<T> * b = dynamic_cast<const gsBasis<T>*>(&_u.source().piece(p)) )
        {
            gsMatrix<T> cf;
            extract(cf, p);
            return b->makeGeometry(give(cf));
        }
        GISMO_ERROR("gsFeSolution: Extraction error");
    }

    // insert g-coefficients to the solution vector
    void insert(const gsGeometry<T> & g, const index_t p = 0) const
    {
        insert(g.coefs(), p);
    }

    // insert g-coefficients to the solution vector
    void insert(const gsMatrix<T> & cf, const index_t p = 0) const
    {
        gsMatrix<T> & sol = *_Sv;
        //gsMatrix<T> & fixedPart = _u.fixedPart();
        const gsDofMapper & mapper = _u.mapper();
        for (index_t c = 0; c!=_u.dim(); c++) // for all components
        {
            for (index_t i = 0; i != cf.rows(); ++i)
            {
                const index_t ii = mapper.index(i, p, c);
                if ( mapper.is_free_index(ii) ) // DoF value is in the solVector
                    sol.at(ii) = cf(i, c);
                /*
                  else
                  {
                  fixedPart.row(m_sd->mapper.global_to_bindex(ii)) = cf.row(i);
                  }
                */
            }
        }
    }
};

}// namespace expr
}// namespace gismo