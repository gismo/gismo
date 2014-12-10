
// assemble mass stiffness moments

#pragma once

#include <ostream>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsConstantFunction.h>
#include <gsPde/gsPde.h>
#include <gsCore/gsDomain.h>
#include <gsCore/gsDofMapper.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

namespace gismo
{




/** @brief
A merger of gsAssembler and gsGaussAssembler
 */

template<class T>
class gsAssemblerUtils
{
private:

    /// Default empty constructor
    gsAssemblerUtils()  { }

    /// Constructor using a geometry

    virtual ~gsAssemblerUtils() { } //destructor

public:


    
    //JS2 add documentation here


    static gsVector<int> getNumIntNodesFor(const gsBasis<T>& b)
    {
        gsVector<int> numNodes( b.dim() );
        for (int i = 0; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }


    static gsVector<int> getNumIntNodesForSide(const gsBasis<T>& b, int dir)
    {
        gsVector<int> numNodes ( b.dim() );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = b.degree(i) + 1;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }

    static gsVector<int> getNumIntNodesForInterface(const gsBasis<T>& b1, const gsBasis<T>& b2,
                           const boundaryInterface & bi, bool left = true)
    {
        // assumes matching orientation
        gsVector<int> numNodes ( b1.dim() );
        const int dir = ( left ? direction( bi.first().side() ) : direction( bi.second().side() ) );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

    static gsVector<int> getNumIntNodesForCoupled(const gsBasis<T>& b1, const gsBasis<T>& b2)
    {
        gsVector<int> numNodes ( b1.dim() );
        for (int i = 0; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }


    static T getMu(const gsBasis<T>& b)
    {
        // TODO: basis should provide a meshSize() method
        const T h = math::pow( (T) b.size(), -1.0 / b.dim() );
        const T bdeg = (T)b.degree(0);
        return ( (bdeg+b.dim())* (bdeg+1) * 2.0 / h );
        //return ( 2.0 / (h * h) );
    }


    static void localToGlobal(const gsMatrix<T>& localStiffness,
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
    
    /// Add contributions from local stiffness matrix/load vector to
    /// global stiffness matrix/load matrix, eliminating Dirichlet BCs
    static void localToGlobal_withBC(const gsMatrix<T>& localStiffness,
                                     const gsMatrix<T>& localRhs,
                                     const gsMatrix<T>& dirbc,
                                     const gsDofMapper& mapper,
                                     const gsVector<int>& loc2glob,
                                     gsSparseMatrix<T>& K,
                                     gsMatrix<T>& f,
                                     bool symmetric)
{
    const int numActive = loc2glob.size();
    for (index_t i=0; i < numActive; ++i)
    {
        const int ii = loc2glob[i];
        if ( mapper.is_free_index(ii) )
        {
            f.row(ii) += localRhs.row(i);
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
                    f.row(ii) -= dirbc.row( mapper.global_to_bindex(jj) ) * localStiffness(i, j);
                }

            }
        }
    }
}


}; // class gsAssembler


} // namespace gismo

