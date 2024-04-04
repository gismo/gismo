/** @file gsMappedBasis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiPatch.h>

#include <gsMSplines/gsMappedSingleBasis.h>
#include <gsMSplines/gsWeightMapper.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsCore/gsBasisFun.h>

namespace gismo
{

// forward declarations of the mapper classes
template<short_t d,class T> class gsMappedSingleBasis;

template<short_t d,class T>
class gsMappedBasis  : public gsFunctionSet<T>
{
    friend class gsMappedSingleBasis<d,T>;

private:
    typedef std::vector<std::pair<int,int> >::iterator step_iter;
    typedef gsBasis<T> BasisType;
    typedef typename std::vector<BasisType *>::const_iterator ConstBasisIter;
    typedef typename std::vector<BasisType *>::iterator BasisIter;
    typedef typename std::vector<gsMatrix<T> *>::const_iterator ConstMatrixPtrIter;
    typedef typename gsSparseMatrix<T,0,index_t>::InnerIterator InIterMat;
    typedef typename std::vector<index_t> IndexContainer;
    typedef typename std::vector<index_t>::const_iterator ConstIndexIter;
    typedef typename std::vector<T> WeightContainer;
    typedef typename std::vector<T>::const_iterator ConstWeightIter;
    typedef gsEigen::PermutationMatrix<Dynamic,Dynamic,index_t> gsPermutationMatrix;

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

public:
    /// Shared pointer for gsMappedBasis
    typedef memory::shared_ptr< gsMappedBasis > Ptr;

    /// Unique pointer for gsMappedBasis
    typedef memory::unique_ptr< gsMappedBasis > uPtr;

    gsMappedBasis() : m_mapper(nullptr)
    { }

    gsMappedBasis(gsMultiBasis<T> const & mb, const gsSparseMatrix<T> & m)
    : m_mapper(nullptr)
    {
        init(mb,m);
    }

    gsMappedBasis( const gsMappedBasis& other );

    virtual ~gsMappedBasis();

    void init(gsMultiBasis<T> const & mb, const gsSparseMatrix<T> & m)
    {
        GISMO_ASSERT(mb.domainDim()==d, "Error in dimensions");
        m_topol  = mb.topology();

        delete m_mapper;
        m_mapper = new gsWeightMapper<T>(m);

        freeAll(m_bases);
        m_bases.reserve(mb.nBases());
        cloneAll(mb.patchBases(),m_bases);
        m_sb.clear();
        m_sb.reserve(m_bases.size());
        for ( size_t q=0; q!=m_bases.size(); ++q )
            m_sb.push_back( gsMappedSingleBasis<d,T>(this,q) );

        m_mapper->optimize(gsWeightMapper<T>::optSourceToTarget);
    }

    index_t nPieces() const {return m_topol.nBoxes();}

public:
    //////////////////////////////////////////////////
    // getters for the private fields
    //////////////////////////////////////////////////

    /// getter for (const) m_bases[i]
    BasisType const & getBase(index_t i) const
    { return *m_bases[i]; }

    /// getter for m_bases[i]
    BasisType & getBase(index_t i)
    { return *m_bases[i]; }

     /// getter for m_bases
    const std::vector<BasisType*> getBases() const;

     /// getter for m_bases
    std::vector<BasisType*> getBasesCopy() const
    {
        std::vector<BasisType*> bases;
        cloneAll(m_bases,bases);
        return bases;
    }

    /// getter for m_mapper
    gsWeightMapper<T> const & getMapper() const
    { return *m_mapper; }

    /// getter for m_mapper
    gsWeightMapper<T> * getMapPointer() const
    { return m_mapper; }

    /// getter for m_topol
    gsBoxTopology const & getTopol() const
    { return m_topol; }

public:
    //////////////////////////////////////////////////
    // general functions for interacting with this class
    //////////////////////////////////////////////////

    /// Clone function. Used to make a copy of a derived basis
    GISMO_CLONE_FUNCTION(gsMappedBasis)

    // Prints the object to the stream
    //std::ostream & print(std::ostream & os) const { }

    short_t domainDim() const
    {
        GISMO_ASSERT(m_bases.size()>0,"there should be at least one basis provided.");
        return m_bases[0]->domainDim();
    }

    /// Returns the approximate number of (global) basis functions in the patch with given \a index
    /// all the sizes summed up should give the total number of (global) basis functions
    /// the default argument -1 will give the total number of basis functions for all patches
    index_t size(const index_t index) const;
    index_t size() const { return m_mapper->getNrOfTargets(); }

    index_t globalSize() const { return m_mapper->getNrOfTargets(); }

    /// Returns the number of local basis functions in the basis
    index_t localSize() const { return m_mapper->getNrOfSources(); }

    /// Returns the number of local basis functions of the patch with given \a index in the basis
    index_t localSize(const index_t index) const
    { return m_bases[index]->size(); }

    /// Returns the maximal polynomial degree of the patches.
    short_t maxDegree() const;

    /// Returns the polynomial degree of direction \a i of \a patch
    short_t degree(const index_t patch,const short_t i) const
    { return m_bases[patch]->degree(i); }

    /// returns the amount of patches of the multi patch
    size_t nPatches() const
    { return m_bases.size(); }

private:
    /// helper function for boundary and innerBoundaries
    void addLocalIndicesOfPatchSide(const patchSide& ps,index_t offset,std::vector<index_t>& locals) const;

public:
    void reorderDofs(const gsPermutationMatrix& permMatrix)
    {
        (*m_mapper)*=permMatrix;
        m_mapper->optimize(gsWeightMapper<T>::optSourceToTarget);
    }

    void reorderDofs_withCoef(const gsPermutationMatrix& permMatrix,gsMatrix<T>& coefs)
    {
        reorderDofs(permMatrix);
        coefs*=permMatrix;
    }

    /// gets all indices of global basis functions on the boundary
    /// upto a given offset
    void boundary(std::vector<index_t> & indices,index_t offset = 0) const;

    /// gets all indices of global basis functions on the inner boundary
    /// upto a given offset
    void innerBoundaries(std::vector<index_t> & indices,index_t offset = 0) const;

    /// gives back the gsMappedSingleBasis object set to the patch i, which
    /// ressembles the composite basis on one patch
    gsMappedSingleBasis<d,T> & getMappedSingleBasis(const index_t i) const
    {
        // if (m_sb.empty())
        // {
        //     m_sb.reserve(m_bases.size());
        //     for (size_t q = 0; q!=m_bases.size(); ++q)
        //         m_sb.push_back( gsMappedSingleBasis<d,T>(this,q) );
        // }
        return m_sb[i];
    }

    const gsMappedSingleBasis<d,T> & piece(const index_t k) const { return getMappedSingleBasis(k); }
    //const gsFunctionSet & piece(const index_t k) const { return m_sb[k]; }

    /// gives back the domain iterator of the boundary side \a s of a given \a patch
    domainIter makeDomainIterator(const index_t patch,const boxSide & s) const
    { return m_bases[patch]->makeDomainIterator(s); }

    /** exports the patch \a i of this geometry (with coefs) to a Geometry object
     *  of the underlying basis. The ownership of this Geometry will go to the
     *  caller of this function.
     *
     * \param[in] localCoef : the coefficients to the local basis functions
     */
    gsGeometry<T>* exportPatch(const index_t i,gsMatrix<T> const & localCoef) const;

    /** exports this geometry (with coefs) to a gsMultiPatch object.
     *
     * \param[in] localCoef : the coefficients to the local basis functions
     */
    gsMultiPatch<T> exportToPatches(gsMatrix<T> const & localCoef) const;

private:
    // Avoid warnings for hidden overloads w.r.t gsFunctionSet
    void active_into(const gsMatrix<T> & u,gsMatrix<index_t>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    void eval_into(const gsMatrix<T> & u,gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    void deriv_into(const gsMatrix<T> & u,gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    void deriv2_into(const gsMatrix<T> & u,gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION; }
    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> >& result ) const
    { GISMO_NO_IMPLEMENTATION; }

public:
    //////////////////////////////////////////////////
    // functions for evaluating and derivatives
    //////////////////////////////////////////////////

    /** \brief Returns the indices of active (non-zero) basis functions of \a patch
     * at points <em>u</em>, as a list of indices, in <em>result</em>.
     *
     * \param[in] u  gsMatrix containing evaluation points. Each column represents one evaluation point.
     * \param[out]  result For every column \a i of \a u, a column containing the indices of the
     *   active basis functions at evaluation point <em>u</em>.col(<em>i</em>).
     */
    void active_into(const index_t patch, const gsMatrix<T> & u,
                     gsMatrix<index_t>& result) const; //global BF active on patch at point

    /// Returns the number of active (nonzero) basis functions at points \a u in \a result.
    virtual void numActive_into(const index_t patch,const gsMatrix<T> & u, gsVector<index_t>& result) const
    {
        GISMO_UNUSED(patch); GISMO_UNUSED(u); GISMO_UNUSED(result);
        GISMO_NO_IMPLEMENTATION
    }

    /// @name Evaluation functions
    /// @{

    /// \brief Evaluates nonzero basis functions of \a patch at point \a u into \a result.
    ///
    /// Let...\n
    /// \a d denote the dimension of the parameter domain.\n
    /// \a k denote the number of active (i.e., non-zero) basis functions (see active_into()).
    /// \a n denote the number of evaluation points.\n
    ///
    /// The \a n <b>evaluation points \a u</b> are given in a gsMatrix of size <em>d</em> x <em>n</em>.
    /// Each column of \a u represents one evaluation point.\n
    /// \n
    /// The gsMatrix <b>\a result</b> contains the computed function values in the following form:\n
    /// Column \a j of \a result corresponds to one evaluation point (specified by the <em>j</em>-th column of \a u).
    /// The column contains the values of all active functions "above" each other.\n
    ///
    /// For example, for scalar basis functions \a Bi : (x,y,z)-> R, a colum represents\n
    /// (B1, B2, ... , Bn)^T,\n
    /// where the order the basis functions \a Bi is as returned by active() and active_into().
    ///
    /// \param[in] u Evaluation points given as gsMatrix of size <em>d</em> x <em>n</em>.
    /// See above for details.
    /// \param[in,out] result gsMatrix of size <em>k</em> x <em>n</em>.
    /// See above for details.
    ///
    void eval_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    void deriv_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    void deriv2_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    gsPiecewiseFunction<T> basisFunction(index_t global_BF)
    {
        const size_t np = nPatches();
        gsPiecewiseFunction<T> bf(np);
        for( size_t i = 0; i!=np; ++i)
            bf.addPiece( getMappedSingleBasis(i).function(global_BF) );
        return bf;
    }

    /// Evaluate the \a global_BF-th basis function on \a patch at points \a u into \a result.
    void evalSingle_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    void derivSingle_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    void deriv2Single_into(const index_t patch, const index_t global_BF, const gsMatrix<T> & u, gsMatrix<T>& result ) const;

    /// Evaluate a single basis function \a i at points \a u.
    gsMatrix<T> evalSingle(const index_t patch, const index_t global_BF, const gsMatrix<T> & u) const
    {
        gsMatrix<T> result;
        this->evalSingle_into(patch, global_BF, u, result);
        return result;
    }

    /// Evaluate a single basis function \a i derivative at points \a u.
    gsMatrix<T> derivSingle(const index_t patch, const index_t global_BF, const gsMatrix<T> & u) const
    {
        gsMatrix<T> result;
        this->derivSingle_into(patch, global_BF, u, result);
        return result;
    }

    /// Evaluate the second derivative of a single basis function \a i at points \a u.
    gsMatrix<T> deriv2Single(const index_t patch, const index_t global_BF, const gsMatrix<T> & u) const
    {
        gsMatrix<T> result;
        this->deriv2Single_into(patch, global_BF, u, result);
        return result;
    }

    /// @brief Evaluate the nonzero basis functions of \a patch and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const index_t patch, const gsMatrix<T> & u,
                          const index_t n, std::vector<gsMatrix<T> >& result ) const;

    /// @brief Evaluate the basis function \a global_BF at \a patch and its derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDersSingle_into(const index_t patch,const index_t global_BF, const gsMatrix<T> & u,const index_t n,gsMatrix<T> & result ) const;

    /// @}

    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

public:
    //////////////////////////////////////////////////
    // functions for converting global to local coefs and back
    //////////////////////////////////////////////////

    /** converts the global coefficients to the local ones
     *
     * \param[in] globalCoefs : the coefficients to the global basis functions
     * \param[out] localCoefs : the coefficients to the local basis functions
     */
    void global_coef_to_local_coef(gsMatrix<T> const & globalCoefs,gsMatrix<T> & localCoefs) const
    { m_mapper->mapToSourceCoefs(globalCoefs,localCoefs); }

    /** converts the local coefficients to the global ones
     *
     * \param[in] localCoefs : the coefficients to the local basis functions
     * \param[out] globalCoefs : the coefficients to the global basis functions
     */
    void local_coef_to_global_coef(gsMatrix<T> const & localCoefs,gsMatrix<T> & globalCoefs) const
    {  m_mapper->mapToTargetCoefs(localCoefs,globalCoefs); }

    /// See gsMappedBasis::_getPatchIndex
    index_t getPatch(index_t const localIndex)
    {
        return _getPatch(localIndex);
    }
    index_t getPatchIndex(index_t const localIndex)
    {
        return _getPatchIndex(localIndex);
    }

protected:
    //////////////////////////////////////////////////
    // functions for working with Indexes
    //////////////////////////////////////////////////

    /** gets the number of the patch, where the local basis function with the
     *  given index is on.
     *
     * \param[in] localIndex : the accumulated index of the basis function
     */
    index_t _getPatch(index_t const localIndex) const;

    /** gets the patchIndex of the basis function specified by its localIndex.
     *  This is the index of this basis function in the basis of this patch.
     *
     * \param[in] localIndex : the accumulated index of the basis function
     */
    index_t _getPatchIndex(const index_t localIndex) const;

    /** gets the local index of the basis function of the patch with the
     *  given number.
     *
     * \param[in] patch : the number of the patch
     * \param[in] patchIndex : the index of the basis function on the patch
     */
    index_t _getLocalIndex(index_t const patch,index_t const patchIndex) const
    { return _getFirstLocalIndex(patch)+patchIndex; }

    /** gets the local index of the first basis function of the patch with the
     *  given number.
     *
     * \param[in] patch : the number of the patch
     */
    index_t _getFirstLocalIndex(index_t const patch) const;

    /** gets the local index of the last basis function of the patch with the
     *  given number.
     *
     * \param[in] patch : the number of the patch
     */
    index_t _getLastLocalIndex(index_t const patch) const
    { return _getFirstLocalIndex(patch)+m_bases[patch]->size()-1; }

    // Data members
protected:
    /// Topology, specifying the relation (connections) between the patches
    gsBoxTopology m_topol;

    /// Vector of local bases
    std::vector<BasisType *> m_bases;

    /// Map between the local basis functions and the newly created ones
    gsWeightMapper<T> * m_mapper;
    // gsSparseMatrix<T> r:C, c:B

    /// Underlying bases per patch
    mutable std::vector<gsMappedSingleBasis<d,T> > m_sb;

    // Make gsMultiBasis a member instead of m_bases and m_topol?
};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsMappedBasis
   */
  // void pybind11_init_gsMappedBasis1(pybind11::module &m);
  void pybind11_init_gsMappedBasis2(pybind11::module &m);
  // void pybind11_init_gsMappedBasis3(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMappedBasis.hpp)
#endif
