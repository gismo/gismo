// gauss rule: input: #nodes
// To assemble planar (2d, 3d) problems

#pragma once

#include <iostream>
#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsPoissonPde.h>
#include <gsPde/gsSurfacePoissonPde.h>
#include <gsPde/gsStokesPde.h>
#include <gsPde/gsConvDiffRePde.h>
#include <gsPde/gsEulerBernoulliBeamPde.h>

namespace gismo
{    

/** @brief
    Implementation of an assembler using Gauss quadrature.
*/
    
template<class T = real_t>
class gsGaussAssembler : public gsAssembler<T>
{
public:
    /// Default empty constructor
    gsGaussAssembler() : gsAssembler<T>() { }
    
    /// Construct by a provided geometry
    gsGaussAssembler(const gsGeometry<T> & geom) : gsAssembler<T>(geom)
    {
        this->setGeometry( geom );
    }
    
    ~gsGaussAssembler()                 //destructor
    { }

    virtual void setGeometry(const gsGeometry<T> & geom)
    { 
      this->m_geometry = &geom;
      // Caution: integration points should be set wrt the degree of
      // the discretization basis, not the geometry map
    }

public:

    static gsVector<int> getNumIntNodesFor(const gsBasis<T>& b)
    {
        gsVector<int> numNodes( b.dim() );
        for (int i = 0; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }

    /// Maps Gauss nodes from side1 to side2
    static void mapGaussNodes(const gsMatrix<T> & nodes1, 
                              const  boundary::side & side1,
                              const  boundary::side & side2,
                              T fixedParam,
                              gsMatrix<T> & nodes2 )
    {
        const int dir1 = direction(side1);
        const int dir2 = direction(side2);
        const int d    = nodes1.rows();
        //const int par2 = parameter(side2);
        nodes2.resize( d, nodes1.cols() );

        if ( dir1 == dir2 )
        {
            nodes2.row(dir2).setConstant(fixedParam);
            //nodes2.row( !dir2 )  = nodes1.row( !dir2 );
            nodes2.topRows(dir2)        = nodes1.topRows( dir2 );
            nodes2.bottomRows(d-dir2-1) = nodes1.bottomRows( d-dir2-1 );
        }
        else
        {
            GISMO_ASSERT( nodes1.rows() == 2, "Implemented for 2D");
            nodes2.row(dir2).setConstant(fixedParam);
            nodes2.row( !dir2 )  = nodes1.row( dir2 );
        }

        /*
        for (int i = 0; i < dir1; ++i)
            nodes2.row(i)  = nodes1.row(i);
        nodes2.row(dir2).array() = fixedParam;
        for (int i = dir+1; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        */

    }

    virtual gsSparseSystem<T>
    assemble( const gsBasis<T>& basis, const gsDofMapper& mapper,
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0)
    {
        const gsPoissonPde<T> *poisson = dynamic_cast<const gsPoissonPde<T>*>(&pde);
        if (poisson)
            return assemblePoisson(basis, mapper, ddof, *poisson, patchIndex);
        
        const gsSurfacePoissonPde<T> *surfacepoisson = 
            dynamic_cast<const gsSurfacePoissonPde<T>*>(&pde);
        if (surfacepoisson)
            return assembleSurfacePoisson(basis, mapper, ddof, *surfacepoisson, patchIndex);
        
        const gsConvDiffRePde<T> *cdr = dynamic_cast<const gsConvDiffRePde<T>*>(&pde);
        if (cdr)
            return assembleCDR(basis, mapper, ddof, *cdr, patchIndex);

        const gsEulerBernoulliBeamPde<T> *ebbeam = dynamic_cast<const gsEulerBernoulliBeamPde<T>*>(&pde);
        if (ebbeam)
            return assembleBeam(basis, mapper, ddof, *ebbeam, patchIndex);

        const gsStokesPde<T> * sst = dynamic_cast<const gsStokesPde<T>*>(&pde);
        if (sst)
            return assembleStokes(basis, mapper, ddof, *sst, patchIndex);

        GISMO_ERROR("Unknown PDE type in assemble()");
    }

    virtual gsSparseMatrix<T> *
    assembleMass( const gsBasis<T>& basis, const gsDofMapper& mapper, int patchIndex=0 );

    /// Assembler for single patch Poisson equation
    gsSparseSystem<T> assemblePoisson( const gsBasis<T>& basis, const gsDofMapper& mapper, 
                                       const gsMatrix<T> & ddof, const gsPoissonPde<T> & pde, 
                                       int patchIndex=0);
    
    /// Assembler for single patch Surface Poisson equation
    gsSparseSystem<T> assembleSurfacePoisson( const gsBasis<T>& basis, const gsDofMapper& mapper, 
                                       const gsMatrix<T> & ddof, const gsSurfacePoissonPde<T> & pde, 
                                       int patchIndex=0);

    /// Assemble convection-diffusion-reaction problem
    gsSparseSystem<T> assembleCDR( const gsBasis<T>& basis, const gsDofMapper& mapper,
                                   const gsMatrix<T> & ddof, const gsConvDiffRePde<T> & pde, 
                                   int patchIndex=0);

    /// Assemble Euler-Bernoulli beam equation
    gsSparseSystem<T> assembleBeam( const gsBasis<T>& basis, const gsDofMapper& mapper,
                                    const gsMatrix<T> & ddof, const gsEulerBernoulliBeamPde<T> & pde,
                                    int patchIndex=0 );

    /// Assembler for Stationary Stokes equations
    gsSparseSystem<T> assembleStokes( const gsBasis<T>& basis_u,
                                      //const gsBasis<T>& basis_p,
                                      const gsDofMapper& mapper,
                                      const gsMatrix<T> & ddof, const gsStokesPde<T> & pde, 
                                      int patchIndex=0);
/*
    /// Assembler a PDE equation on a single patch -- under testing
    gsSparseSystem<T> assemblePde( const gsBasis<T>& basis , const gsDofMapper& mapper, 
                                   const gsMatrix<T> & ddof, const gsPoissonPde<T> & pde, 
                                   int patchIndex=0);
*/
    void applyBoundary( const gsBasis<T>   & B,
                        const boundary_condition<T> & bc,
                        const gsDofMapper& mapper,
                        gsSparseSystem<T> & system );
    
    void applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                  const gsGeometry<T> & geo2,
                  const boundaryInterface & bi,
                  const gsDofMapper& mapper,
                  gsSparseSystem<T> & system );


    gsSparseMatrix<T> * stiffness( const gsBasis<T>& B);

    gsSparseMatrix<T> * massMatrix( const gsBasis<T>& B );
    
    gsVector<T> * moments( const gsBasis<T>& B, gsFunction<T> const & f ) ;
    
    /// Compute Boundary moments of basis B with respect to function f
    /// along side s of geometry
    gsVector<T> * boundaryMoments( const gsBasis<T>& B,
                                   gsFunction<T> const& f,
                                   boundary::side const& s );
    
    static gsMatrix<T> * innerProduct ( const gsBasis<T>& B1, const gsBasis<T>& B2);
    static gsMatrix<T> * innerProduct1( const gsBasis<T>& B1, const gsBasis<T>& B2);
    static gsMatrix<T> * innerProduct2( const gsBasis<T>& B1, const gsBasis<T>& B2);

    /// Compute the determinants needed for computing a surface
    /// integral at all columns of x
    static void computeSurfaceDets_into( const gsGeometry<T>& surf, 
                                         const gsMatrix<T>& x, gsMatrix<T>& result );

private:

    /// Add contribution of Nitsche Dirichlet boundary to matrix K
    /// \param B is a boundary basis, \param f is the Dirichlet function
    void boundaryNitsche( const gsBasis<T>   & B,
                          const int patch, 
                          const boundary::side s,
                          const gsFunction<T> & f,
                          const gsDofMapper& mapper,
                          gsSparseSystem<T> & system );

    /// Add contribution of Neumann boundary condition to he \a system
    /// \param B is a boundary basis, \param f is the Dirichlet function    
    void boundaryNeumann( const gsBasis<T>   & B,
                          const int patch, 
                          const boundary::side s,
                          const gsFunction<T> & f,
                          const gsDofMapper& mapper,
                          gsSparseSystem<T> & system );

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
                                                    const boundary::side side1, 
                                                    const boundary::side side2,
                                                    bool left = true)
    {
        // assumes matching orientation
        gsVector<int> numNodes ( b1.dim() );
        const int dir = ( left ? direction( side1 ) : direction( side2 ) );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = math::max( b1.degree(i), b2.degree(i) ) + 1 ;
            //numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b1.dim(); ++i)
            numNodes[i] = math::max( b1.degree(i), b2.degree(i) ) + 1 ;
            //numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

    static gsVector<int> getNumIntNodesForCoupled(const gsBasis<T>& b1, const gsBasis<T>& b2)
    {
        gsVector<int> numNodes ( b1.dim() );
        for (int i = 0; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

    T minElementVol(const gsBasis<T>& b);

}; // class gsGaussAssembler
    

// function B=adj(A)
// n=size(A,1); d=1:n-1;
// B=zeros(n); c=2*mod(size(A,1),2)-1;
// AA=[A,A;A,A]';
// for j=1:n
//     for k=1:n
//         B(j,k)=c^(j+k)*det(AA(j+d,k+d));
//     end
// end


//////////////////////////////////////////////////
//////////////////////////////////////////////////

} // namespace gismo



#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsGaussAssembler.hpp)
#endif
