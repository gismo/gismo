/** @file gsFeSpace.h

    @brief Defines a space expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for a geometry map
 * @ingroup Expressions
 * @tparam T The scalar type
 */
template<class T>
class gsGeometryMap : public _expr<gsGeometryMap<T> >
{
    const gsFunctionSet<T> * m_fs; ///< Evaluation source for this geometry map
    const gsMapData<T>     * m_fd; ///< Temporary variable storing flags and evaluation data
    //index_t d, n;

    bool m_isAcross; ///< true when the patch evaluated is across an interface

public:
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    bool isAcross() const { return m_isAcross; }

    gsGeometryMap right() const
    {
        gsGeometryMap ac;
        ac.m_fs = m_fs;
        ac.m_isAcross = true;
        return ac;
    }

    gsGeometryMap left() const
    {
        gsGeometryMap ac;
        ac.m_fs = m_fs;
        ac.m_isAcross = false;
        return ac;
    }

    /// Returns the function source
    const gsFunctionSet<T> & source() const {return *m_fs;}

    /// Returns the function data
    const gsMapData<T> & data() const
    {
        GISMO_ASSERT(NULL!=m_fd, "gsGeometryMap: invalid data "<< m_fs <<","<<m_fd);
        return *m_fd;
    }

    index_t targetDim() const { return m_fs->targetDim();}
    index_t domainDim() const { return m_fs->domainDim();}

    /// Copy the coefficients of another gsGeometryMap to this one, if they are compatible.
    void copyCoefs( const gsGeometryMap<T> & other) const
    {
        GISMO_ASSERT( dynamic_cast<const gsMultiPatch<T>*>( this->m_fs ), "error");
        const gsMultiPatch<T> & thisMP  = static_cast<const gsMultiPatch<T>&>(*this->m_fs );
        GISMO_ASSERT( dynamic_cast<const gsMultiPatch<T>*>( other.m_fs ), "error");
        const gsMultiPatch<T> & otherMP = static_cast<const gsMultiPatch<T>&>(*other.m_fs );
        GISMO_ASSERT( (thisMP.domainDim()==otherMP.domainDim())&&
                      (thisMP.geoDim()==otherMP.geoDim())&&
                      (thisMP.coefsSize() == otherMP.coefsSize())&&
                      (thisMP.nPatches()==otherMP.nPatches()),
                "The geometryMaps are not compatible!");

        // For every patch of the MultiPatch
        for ( index_t p=0; p < thisMP.nPatches(); p++ )
        {
            // Copy coeffs of the other MultiPatch
            thisMP.patch(p).coefs() = otherMP.patch(p).coefs();
        }

    }   //end copyCoffs

    void deformBy( const gsFeSolution<T> & displacement) const
    {
        const index_t dim = m_fs->domainDim();

        const gsMultiBasis<T> & mb = static_cast<const gsMultiBasis<T>&>(displacement.space().source());
        const gsMultiPatch<T> & mp = static_cast<const gsMultiPatch<T>&>(*this->m_fs );
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<T>*>(&displacement.space().source()), "error");
        GISMO_ASSERT( dynamic_cast<const gsMultiPatch<T>*>( this->m_fs), "error");

        // For every patch of the MultiPatch
        for ( size_t p=0; p < mp.nPatches(); p++ )
        {
            // Get the patch's coefficients
            gsMatrix<T> &result = mp.patch(p).coefs();

            // Number of basis functions of patch with index p
            const index_t sz  = mb[p].size();

            // For all components
            for (index_t c = 0; c!=dim; c++)
            {
                // loop over all basis functions (even the eliminated ones)
                for (index_t i = 0; i < sz; ++i)
                {
                    const int ii = displacement.mapper().index(i, p, c);
                    if ( displacement.mapper().is_free_index(ii) ) // DoF value is in the defVector
                    {
                        result(i,c) += displacement.coefs().at(ii);
                    }
                    else
                    {
                        result(i,c) += displacement.fixedPart().at( displacement.mapper().global_to_bindex(ii));
                    }
                }
            }
        }
    } // end deformBy

public:
    typedef T Scalar;

    friend class gismo::gsExprHelper<Scalar>;

    void print(std::ostream &os) const { os << "G"; }

    auto eval(const index_t k) const -> decltype(m_fd->values[0].col(k))
    { return m_fd->values[0].col(k); }

protected:

    gsGeometryMap() : m_fs(NULL), m_fd(NULL), m_isAcross(false) { }

    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}
    void setData(const gsMapData<Scalar> & val) { m_fd = &val;}

public:

    index_t rows() const { return m_fs->targetDim(); }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(*this);
        m_fd->flags |= NEED_VALUE;
    }
};

// Traits for gsGeometryMap
template <typename T>  struct expr_traits<gsGeometryMap<T> >
{
    typedef T Scalar;
    typedef const gsGeometryMap<T> Nested_t; // nesting without ref!
};

}// namespace expr
}// namespace gismo