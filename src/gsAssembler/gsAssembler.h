
// assemble mass stiffness moments
 
#pragma once

#include <ostream>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsDofMapper.h>
#include <gsPde/gsPde.h>
//#include <gismo.h>


namespace gismo
{

/// @brief A linear system consisting of sparse matrix and right hand side.
///
/// This is a very simple class which only holds pointers to a sparse matrix
/// and a right-hand side vector. It does no processing on its own.
template <class T>
class gsSparseSystem
{
public:
    /// Empty constructor.
    gsSparseSystem()
    { m_matrix = 0; m_rhs = 0; }

    /// Constructor from given matrix and right-hand side
    gsSparseSystem(gsSparseMatrix<T>* m, gsVector<T>* rhs)
    { m_matrix = m; m_rhs = rhs; }

    /// Return the matrix.
    gsSparseMatrix<T>* matrix() const  { return m_matrix; }

    /// Return the right-hand side.
    gsVector<T>*       rhs   () const  { return m_rhs; }

    /// Return the number of degrees of freedom in the system.
    index_t size()      { return this->matrix()->cols(); }

    /// Delete the matrix and right-hand side.
    void free()
    {
        delete m_matrix; m_matrix = 0;
        delete m_rhs;    m_rhs = 0;
    }

    /// Print the linear system to a stream.
    friend std::ostream &operator<<(std::ostream &os, const gsSparseSystem<T> & system)
    { 
        os << "gsSparseSystem:\n";
        if ( system.m_matrix )
            os << "* Left hand side:\n" << *system.m_matrix ;    
        if ( system.m_rhs )
            os << "Right hand side (transposed):\n"  << system.m_rhs->transpose() ;
        os<<"\n";
        return os;
    }

private:
    gsSparseMatrix<T> *m_matrix;
    gsVector<T> *m_rhs;
};



/** @brief
    Abstract base class for assemblers.
    
    An assembler takes a single-patch geometry, a discretization basis on that
    geometry, and a PDE and assembles the patch-local stiffness matrix and
    right-hand side for the PDE using the given basis.
    It can also apply various boundary conditions.

    An assembler is initalized by passing it a gsGeometry describing the
    patch. A problem is then assembled by calling gsAssembler::assemble().
 */

template<class T>
class gsAssembler
{
public:

    /// Default empty constructor
    gsAssembler() : m_geometry(NULL) { }
    
    /// Constructor using a geometry
    gsAssembler(const gsGeometry<T> & geom)
    { 
        m_geometry = &geom;
    }
    
    virtual ~gsAssembler() { } //destructor
    
public:

    /// Get the assembler's geometry.
    const gsGeometry<T> & geometry() const    { return *m_geometry; }

    /// Set the assembler's geometry.
    virtual void setGeometry(const gsGeometry<T> & geom)
    { 
        m_geometry = &geom;
    } 

    /** @brief
        Assemble a sparse linear system for the given PDE on the assembler's geometry.

        \param basis    the discretization basis
        \param mapper   the gsDofMapper which describes interior, interface, and boundary dofs
        \param ddof     values for eliminated Dirichlet dofs
        \param pde      the partial differential equation to be assembled
        \param patchIndex index of the patch currently being assembled
        \return the sparse linear system that was assembled
     */
    virtual gsSparseSystem<T>
    assemble( const gsBasis<T>& basis, const gsDofMapper& mapper, 
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0) = 0;

    /** @brief
        Assemble the mass matrix for the given discretization basis on the assembler's geometry.

        \param basis    the discretization basis
        \param mapper   the gsDofMapper which describes interior, interface, and boundary dofs
        \param patchIndex index of the patch currently being assembled
        \return         the mass matrix as a sparse matrix
     */
    virtual gsSparseMatrix<T> *
    assembleMass( const gsBasis<T>& basis, const gsDofMapper& mapper, int patchIndex=0 )
    { GISMO_NO_IMPLEMENTATION }

    /// \brief Returns the stiffness matrix (for Laplace/Poisson problem) for a given basis.
    ///
    /// \warning If the problem is symmetric, only the lower part is filled!
    ///
    /// To get a full view use
    ///
    /// gsSparseMatrix<T> Afull;
    ///  Afull = * safe(A.template selfadjointView<Lower>();
    ///
    /// or
    ///
    ///    gsSparseMatrix<> Afull;
    ///    Afull = A.selfadjointView<Lower>();
    ///
    ///
    ///
    /// Or to copy the matrix to the complete matrix use
    ///
    /// gsSparsematrix<T> Acomplete = A.selfadjointView<Lower>();
    ///
    ///
    virtual gsSparseMatrix<T> * stiffness( const gsBasis<T>& B )
    { GISMO_NO_IMPLEMENTATION }


    gsSparseMatrix<T> * stiffnessFull( const gsBasis<T>& B )
    {
        gsSparseMatrix<T> * L = new gsSparseMatrix<T>;
        *L = safe( stiffness(B) )->template selfadjointView<Lower>();
        L->makeCompressed();
        return L;
    }

    /// \brief Returns the mass matrix for a given basis.
    ///
    /// \warning If the problem is symmetric, only the lower part is filled!
    virtual gsSparseMatrix<T> * massMatrix( const gsBasis<T>& B )
    { GISMO_NO_IMPLEMENTATION }

    gsSparseMatrix<T> * massMatrixFull( const gsBasis<T>& B )
    {
        gsSparseMatrix<T> * L = new gsSparseMatrix<T>;
        *L = safe( massMatrix(B) )->template selfadjointView<Lower>();
        L->makeCompressed();
        return L;
    }

    virtual gsVector<T> * moments( const gsBasis<T>& B, gsFunction<T> const& f )
    { GISMO_NO_IMPLEMENTATION }

    /// Compute Boundary moments of basis B with respect to function f
    /// along side s of geometry
    virtual gsVector<T> * boundaryMoments( const gsBasis<T>& B,
                                           gsFunction<T> const& f,
                                           boundary::side const& s )
    { GISMO_NO_IMPLEMENTATION }

    /// Add contribution of boundary condition \a bc to the system \a system
    /// \param B        discretization basis for the patch
    /// \param bc       the boundary condition to apply
    /// \param mapper   contains the global Dof map for the \a system
    /// \param system   the linear system where the contribution shall be added
    virtual void applyBoundary( const gsBasis<T> & B,
                        const boundary_condition<T> & bc,
                        const gsDofMapper& mapper,
                        gsSparseSystem<T> & system )
    { GISMO_NO_IMPLEMENTATION }

    /// Add contribution of interface \a bi to the system \a system
    /// \param B1       discretization basis for the first patch
    /// \param B2       discretization basis for the second patch
    /// \param geo2     the adjacent patch to m_geometry
    /// \param bi       description of the interface which joins the two patches
    /// \param mapper   contains the global Dof map for the \a system
    /// \param system   the linear system where the contribution shall be added
    virtual void applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                          const gsGeometry<T> & geo2,
                          const boundaryInterface & bi,
                          const gsDofMapper& mapper,
                          gsSparseSystem<T> & system )
    { GISMO_NO_IMPLEMENTATION }


    // generic helper functions for assembling routines


    /// Add contributions from local to global stiffness matrix, identity DOF mapping
    static void localToGlobal(const gsMatrix<T>& localStiffness, 
                              const gsMatrix<unsigned>& localDofs,
                              gsSparseMatrix<T>& K, 
                              bool symmetric);

    /// Add contributions from local stiffness matrix/load vector to
    /// global stiffness matrix/load vector, eliminating Dirichlet BCs
    static void localToGlobal_withBC(const gsMatrix<T>& localStiffness, 
                                     const gsVector<T>& localRhs, 
                                     const gsMatrix<T>& dirbc, 
                                     const gsDofMapper& mapper, 
                                     const gsVector<int>& loc2glob, 
                                     gsSparseMatrix<T>& K, gsVector<T>& f, 
                                     bool symmetric);

    /// Add contributions from local to global stiffness matrix, eliminating Dirichlet BCs
    static void localToGlobal_withBC(const gsMatrix<T>& localStiffness,
                                     const gsDofMapper& mapper,
                                     const gsVector<int>& loc2glob,
                                     gsSparseMatrix<T>& K,
                                     bool symmetric);

    /// Allocate a sparse matrix and reserve enough space in it for assembling using the given matrix
    static gsSparseMatrix<T> *initMatrix(int nDofs, const gsBasis<T>& basis, bool symmetric);

    /// Allocate sparse matrix and right-hand side and initialize them
    static gsSparseSystem<T> initLinearSystem(int nDofs, const gsBasis<T>& basis, bool symmetric);

    T getMu(const gsBasis<T>& b)
    {
        // TODO: basis should provide a meshSize() method
        const T h = math::pow( (T) b.size(), -1.0 / b.dim() );
        const T bdeg = (T)b.degree(0);
        return ( (bdeg+b.dim())* (bdeg+1) * 2.0 / h );
        //return ( 2.0 / (h * h) );
    }

// Data members
protected:

    const gsGeometry<T> * m_geometry;

}; // class gsAssembler


//////////////////////////////////////////////////
//////////////////////////////////////////////////



template <class T>
void gsAssembler<T>::localToGlobal(const gsMatrix<T>& localStiffness, 
                                   const gsMatrix<unsigned>& localDofs,
                                   gsSparseMatrix<T>& K, 
                                   bool symmetric)
{
    const int numActive = localDofs.rows();

    for (index_t i = 0; i < numActive; ++i)
    {
        const int ii = localDofs(i,0);
        for (index_t j = 0; j < numActive; ++j)
        {
            const int jj = localDofs(j,0);
            // if matrix is symmetric, store only lower triangular part
            if (!symmetric || jj <= ii)
                K.coeffRef(ii, jj) += localStiffness(i, j);
        }
    }
}


template <class T>
void gsAssembler<T>::localToGlobal_withBC(const gsMatrix<T>& localStiffness,
                                          const gsVector<T>& localRhs,
                                          const gsMatrix<T>& dirbc,
                                          const gsDofMapper& mapper,
                                          const gsVector<int>& loc2glob,
                                          gsSparseMatrix<T>& K,
                                          gsVector<T>& f,
                                          bool symmetric)
{
    const int numActive = loc2glob.size();

    for (index_t i=0; i < numActive; ++i)
    {
        const int ii = loc2glob[i];
        if ( mapper.is_free_index(ii) )
        {
            f[ii] += localRhs[i];

            for (index_t j=0; j < numActive; ++j)
            {
                const int jj = loc2glob[j];
                if ( mapper.is_free_index(jj) )
                {
                    // if matrix is symmetric, store only lower triangular part
                    if (!symmetric || jj <= ii)
                        K.coeffRef(ii, jj) += localStiffness(i, j);
                }
                else if ( mapper.is_boundary_index(jj) )        // Dirichlet boundary condition?
                {
                    f[ii] -= dirbc( mapper.global_to_bindex(jj) ) * localStiffness(i, j);
                }
            }
        }
    }
}

/// Add contributions from local to global matrix, eliminating Dirichlet BCs
template <class T>
void gsAssembler<T>::localToGlobal_withBC(const gsMatrix<T>& localMatrix,
                                          const gsDofMapper& mapper,
                                          const gsVector<int>& loc2glob,
                                          gsSparseMatrix<T>& K,
                                          bool symmetric)
{
    const int numActive = loc2glob.size();

    for (index_t i=0; i < numActive; ++i)
    {
        const int ii = loc2glob[i];
        if ( !mapper.is_free_index(ii) )
            continue;

        for (index_t j=0; j < numActive; ++j)
        {
            const int jj = loc2glob[j];
            if ( !mapper.is_free_index(jj) )
                continue;

            // if matrix is symmetric, store only lower triangular part
            if (!symmetric || jj <= ii)
                K.coeffRef(ii, jj) += localMatrix(i, j);
        }
    }
}


template <class T>
gsSparseMatrix<T> * gsAssembler<T>::initMatrix(int nDofs, const gsBasis<T>& basis, bool symmetric)
{
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(nDofs, nDofs);

    // estimate max nz per row (only valid for tensor product bases right now)
    unsigned nzRowsPerCol = 1;
    for (int i = 0; i < basis.dim(); ++i)
        nzRowsPerCol *= 2 * basis.degree(i) + 1;    // diagonal entry plus p interactions on each side

    // Unfortunately, this optimization doesn't work in general because Neumann/interface
    // dofs are not ordered the same way as interior dofs, thus they can have more nzs.
    //if (symmetric)
    //    nzRowsPerCol = (nzRowsPerCol + 1) / 2;
    if (nDofs > 0)
    K->reserve( gsVector<int>::Constant(nDofs, nzRowsPerCol) );
    return K;
}

template <class T>
gsSparseSystem<T> gsAssembler<T>::initLinearSystem(int nDofs, const gsBasis<T>& basis, bool symmetric)
{
    gsSparseMatrix<T> * K = initMatrix(nDofs, basis, symmetric);
    gsVector<T> * b = new gsVector<T>;
    b->setZero( nDofs );
    return gsSparseSystem<T>( K, b );
}


}; // namespace gismo

