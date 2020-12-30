/** @file gsAssembler.hpp

    @brief Provides generic assembler routines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals and load vector

#include <gsCore/gsFuncData.h>


namespace gismo
{

template<class T>
gsOptionList gsAssembler<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("DirichletStrategy", "Method for enforcement of Dirichlet BCs [11..14]", 11 );
    opt.addInt("DirichletValues"  , "Method for computation of Dirichlet DoF values [100..103]", 101);
    opt.addInt("InterfaceStrategy", "Method of treatment of patch interfaces [0..3]", 1  );
    opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    opt.addReal("bdA", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 2.0  );
    opt.addInt ("bdB", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 1    );
    opt.addReal("bdO", "Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]", 0.333);
    return opt;
}

template<class T>
void gsAssembler<T>::refresh()
{
    gsWarn <<" gsAssembler<T>::Refresh is an empty call\n";
}

template<class T>
void gsAssembler<T>::assemble()
{GISMO_NO_IMPLEMENTATION}

template<class T>
void gsAssembler<T>::assemble(const gsMultiPatch<T> &)
{GISMO_NO_IMPLEMENTATION}

template<class T>
gsAssembler<T> * gsAssembler<T>::create() const
{GISMO_NO_IMPLEMENTATION}

template<class T>
gsAssembler<T> * gsAssembler<T>::clone() const
{GISMO_NO_IMPLEMENTATION}

template<class T>
bool gsAssembler<T>::check()
{
    const gsBoundaryConditions<T> & m_bConditions = m_pde_ptr->bc();

    // Check if boundary conditions are OK
    const index_t np = m_bases.front().nBases();
    for (typename gsBoundaryConditions<T>::const_iterator it =
         m_bConditions.dirichletBegin() ; it != m_bConditions.dirichletEnd(); ++it )
    {
        GISMO_ENSURE( it->ps.patch < np && it->ps.patch >= 0,
                      "Problem: a Dirichlet boundary condition is set "
                      "on a patch id which does not exist.");
    }

    for (typename gsBoundaryConditions<T>::const_iterator it =
         m_bConditions.neumannBegin() ; it != m_bConditions.neumannEnd(); ++it )
    {
        GISMO_ENSURE( it->ps.patch < np && it->ps.patch >= 0,
                      "Problem: a Neumann boundary condition is set "
                      "on a patch id which does not exist.");
    }

    //TODO: add check m_bases[i].nBases() == pde.domain().nPatches().

    if ( m_pde_ptr->domain().nPatches() == 0)
        gsWarn<< "No domain given ! \n";

    // /*
    const dirichlet::strategy dirStr = static_cast<dirichlet::strategy>(m_options.getInt("DirichletStrategy"));
    if ( 0 == m_pde_ptr->bc().size() && dirStr!=dirichlet::none && static_cast<dirichlet::values>(dirStr)==dirichlet::homogeneous )
        gsWarn<< "No boundary conditions given ! \n";
    //*/

    return true;
}

template <class T>
void gsAssembler<T>::scalarProblemGalerkinRefresh()
{
    // Check for coherency
    GISMO_ASSERT(this->check(), "Incoherent data in Assembler");

    GISMO_ASSERT(1==m_bases.size(), "Expecting a single discrete space "
                                    "for standard scalar Galerkin");

    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper mapper;
    m_bases.front().getMapper(
        static_cast<dirichlet::strategy>(m_options.getInt("DirichletStrategy")),
        static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
        this->pde().bc(), mapper, 0);

    if ( 0 == mapper.freeSize() ) // Are there any interior dofs ?
        gsWarn << " No internal DOFs, zero sized system.\n";

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(mapper);//1,1
}

template<class T>
void gsAssembler<T>::penalizeDirichletDofs(short_t unk)
{
    GISMO_ASSERT( m_options.getInt("DirichletStrategy")
                  == dirichlet::penalize, "Incorrect options");

    static const T PP = 1e9; // magic number

    const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(unk)];
    const gsDofMapper     & mapper = m_system.colMapper(unk);
    // Note: dofs in the system, however dof values need to be computed
    const gsDofMapper & bmap = mbasis.getMapper(dirichlet::elimination,
                       static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                                                m_pde_ptr->bc(), unk) ;

    GISMO_ENSURE( m_ddof[unk].rows() == mapper.boundarySize() &&
                  m_ddof[unk].cols() == m_pde_ptr->numRhs(),
                  "The Dirichlet DoFs were not computed.");

    // BCs
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = m_pde_ptr->bc().dirichletBegin();
          it != m_pde_ptr->bc().dirichletEnd(); ++it )
    {
        const gsBasis<T> & basis = mbasis[it->patch()];

        gsMatrix<index_t> bnd = basis.boundary(it->side() );
        for (index_t k=0; k!= bnd.size(); ++k)
        {
            // free dof position
            const index_t ii = mapper.index ( bnd(k) , it->patch() );
            // boundary dof position
            const index_t bb = bmap  .bindex( bnd(k) , it->patch() );

            m_system.matrix()(ii,ii) = PP;
            m_system.rhs().row(ii)   = PP * m_ddof[unk].row(bb);
        }
    }

    // Corner values
    for ( typename gsBoundaryConditions<T>::const_citerator
          it = m_pde_ptr->bc().cornerBegin();
          it != m_pde_ptr->bc().cornerEnd(); ++it )
    {
        const index_t i  = mbasis[it->patch].functionAtCorner(it->corner);
        const index_t ii = mapper.bindex( i , it->patch );
        m_system.matrix()(ii,ii)       = PP;
        m_system.rhs().row(ii).setConstant(PP * it->value);
    }
}

template<class T>
void gsAssembler<T>::setFixedDofs(const gsMatrix<T> & coefMatrix, short_t unk, size_t patch)
{
    GISMO_ASSERT( m_options.getInt("DirichletValues") == dirichlet::user, "Incorrect options");

    const dirichlet::strategy dirStr = static_cast<dirichlet::strategy>(m_options.getInt("DirichletStrategy"));
    const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(unk)];
    const gsDofMapper & mapper =
        dirichlet::elimination == dirStr ? m_system.colMapper(unk)
        : mbasis.getMapper(dirichlet::elimination,
                           static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                           m_pde_ptr->bc(), unk) ;

    GISMO_ASSERT(m_ddof[unk].rows()==mapper.boundarySize() &&
                 m_ddof[unk].cols() == m_pde_ptr->numRhs(),
                 "Fixed DoFs were not initialized");

    // for every side with a Dirichlet BC
    for ( typename gsBoundaryConditions<T>::const_iterator
          it =  m_pde_ptr->bc().dirichletBegin();
          it != m_pde_ptr->bc().dirichletEnd()  ; ++it )
    {
        const index_t k = it->patch();
        if ( k == static_cast<index_t>(patch) )
        {
            // Get indices in the patch on this boundary
            const gsMatrix<index_t> boundary =
                    mbasis[k].boundary(it->side());

            //gsInfo <<"Setting the value for: "<< boundary.transpose() <<"\n";

            for (index_t i=0; i!= boundary.size(); ++i)
            {
                // Note: boundary.at(i) is the patch-local index of a
                // control point on the patch
                const index_t ii  = mapper.bindex( boundary.at(i) , k );

                m_ddof[unk].row(ii) = coefMatrix.row(boundary.at(i));
            }
        }
    }
}

template<class T>
void gsAssembler<T>::setFixedDofVector(gsMatrix<T> vals, short_t unk)
{
   if(m_ddof.size()==0)
       m_ddof.resize(m_system.numColBlocks());
    m_ddof[unk].swap(vals);
    // Assuming that the DoFs are already set by the user
    GISMO_ENSURE( m_ddof[unk].rows() == m_system.colMapper(unk).boundarySize()
                  , //&& m_ddof[unk].cols() == m_pde_ptr->numRhs(),
                  "The Dirichlet DoFs were not provided correctly.");
}


template<class T>
void gsAssembler<T>::computeDirichletDofs(short_t unk)
{
    //if ddof-size is not set
    //fixme: discuss if this is really the right place for this.

    // correction: m_ddof should have an array of ddofs for every unknown,
    // not for every block in the sparse system;
    // number of columns in each array is a number of components in the corresponding unknown;
    // sparse system gets this info during construction;
    // m_system.numUnknown() provides the number of unknown;
    // m_system.unkSize(i) provides the number of components in the unknown i
    if(m_ddof.size()==0)
        m_ddof.resize(m_system.numUnknowns());

    if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche)
        return; // Nothing to compute

    const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(unk)];
    const gsDofMapper & mapper =
            dirichlet::elimination == m_options.getInt("DirichletStrategy") ?
        m_system.colMapper(unk) :
        mbasis.getMapper(dirichlet::elimination,
                         static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                         m_pde_ptr->bc(), unk);

    //gsDebugVar(m_options.getInt("DirichletAAAStrategy"));
    //gsDebugVar(m_options.getInt("DirichletValues"));

    switch ( m_options.getInt("DirichletValues") )
    {
    case dirichlet::homogeneous:
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros
        m_ddof[unk].setZero(mapper.boundarySize(), m_system.unkSize(unk) * m_pde_ptr->numRhs() );
        break;
    case dirichlet::interpolation:
        computeDirichletDofsIntpl(mapper, mbasis,unk);
        break;
    case dirichlet::l2Projection:
        computeDirichletDofsL2Proj(mapper, mbasis,unk);
        break;
    case dirichlet::user :
         // Assuming that the DoFs are already set by the user
        GISMO_ENSURE( m_ddof[unk].size() == mapper.boundarySize()*m_system.unkSize(unk)* m_pde_ptr->numRhs(), "The Dirichlet DoFs are not set.");
        m_ddof[unk].resize(mapper.boundarySize(), m_system.unkSize(unk)* m_pde_ptr->numRhs());
            break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values.");
    }

    // Corner values
    for ( typename gsBoundaryConditions<T>::const_citerator
          it = m_pde_ptr->bc().cornerBegin();
          it != m_pde_ptr->bc().cornerEnd(); ++it )
    {
        if(it->unknown == unk)
        {
            const index_t i  = mbasis[it->patch].functionAtCorner(it->corner);
            const index_t ii = mapper.bindex( i , it->patch );
            m_ddof[unk].row(ii).setConstant(it->value);
        }
        else
            continue;

    }
}


// SKleiss: Note that this implementation is not useable for (T)HB-Splines!
//
// 1. Computation of the Dirichlet values explicitly uses gsPointGrid(rr)
// Computing a grid of evaluation points does not make sense for any locally
// refined basis.
// Also, "component(i)" is used.
// I'm afraid this makes sense ONLY FOR TENSOR-PRODUCT bases.
//
// 2. As of now (16.May 2014), the boundaryBasis of (T)HB-spline basis is not
// implemented, as far as I know.
//
// 3. gsInterpolate uses the anchors of the boundary basis.
// With truncated hierarchical B-splines, the use of classical anchors does
// not work, because functions might be truncated to zero at these points.
template<class T> //
void gsAssembler<T>::computeDirichletDofsIntpl(const gsDofMapper & mapper,
                                               const gsMultiBasis<T> & mbasis,
                                               const short_t unk_)
{
    m_ddof[unk_].resize(mapper.boundarySize(), m_system.unkSize(unk_) * m_pde_ptr->numRhs() );
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = m_pde_ptr->bc().dirichletBegin();
          it != m_pde_ptr->bc().dirichletEnd(); ++it )
    {
        //const int unk = it->unknown();
        const index_t k   = it->patch();
        if(it->unknown()!=unk_)
            continue;
        const gsBasis<T> & basis = mbasis[k];

        // Get dofs on this boundary
        const gsMatrix<index_t> boundary = basis.boundary(it->side());

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const index_t ii= mapper.bindex( boundary.at(i) , k );
                m_ddof[unk_].row(ii).setZero();
            }
            continue;
        }

        // Get the side information
        short_t dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( this->patches().parDim() );

        for ( short_t i=0; i < this->patches().parDim(); ++i)
        {
            if ( i==dir )
            {
                gsVector<T> b(1);
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {
                rr.push_back( basis.component(i).anchors().transpose() );
            }
        }

        GISMO_ASSERT(it->function()->targetDim() == m_system.unkSize(unk_) * m_pde_ptr->numRhs(),
                     "Given Dirichlet boundary function does not match problem dimension."
                     <<it->function()->targetDim()<<" != "<<m_system.unkSize(unk_) << " * " << m_system.rhs().cols()<<"\n");

        // Compute dirichlet values
        gsMatrix<T> fpts;
        if ( it->parametric() )
            fpts = it->function()->eval( gsPointGrid<T>( rr ) );
        else
            fpts = it->function()->eval( m_pde_ptr->domain()[it->patch()].eval(  gsPointGrid<T>( rr ) ) );

        // Interpolate dirichlet boundary
        typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            const index_t ii = mapper.bindex( boundary.at(l) , it->patch() );
            m_ddof[unk_].row(ii) = dVals.row(l);
        }
    }
}

template<class T>
void gsAssembler<T>::computeDirichletDofsL2Proj(const gsDofMapper & mapper,
                                                const gsMultiBasis<T> & ,
                                                const short_t unk_)
{
    m_ddof[unk_].resize( mapper.boundarySize(), m_system.unkSize(unk_)* m_pde_ptr->numRhs());

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( mapper.boundarySize(), m_system.unkSize(unk_)* m_pde_ptr->numRhs() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<index_t> globIdxAct;
    gsMatrix<T> basisVals;

    gsMapData<T> md(NEED_MEASURE);

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          iter = m_pde_ptr->bc().dirichletBegin();
          iter != m_pde_ptr->bc().dirichletEnd(); ++iter )
    {
        if (iter->isHomogeneous() )
            continue;

        GISMO_ASSERT(iter->function()->targetDim() == m_system.unkSize(unk_)* m_pde_ptr->numRhs(),
                     "Given Dirichlet boundary function does not match problem dimension."
                     <<iter->function()->targetDim()<<" != "<<m_system.unkSize(unk_)<<"x"<<m_system.rhs().cols()<<"\n");

        const short_t unk = iter->unknown();
        if(unk!=unk_)
            continue;
        const index_t patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            md.points, quWeights);

            //geoEval->evaluateAt( md.points );
            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            basis.eval_into( md.points, basisVals);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(md.points.col(0), globIdxAct );
            mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index( globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            // Do the actual assembly:
            for( index_t k=0; k < md.points.cols(); k++ )
            {
                const T weight_k = quWeights[k] * md.measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const unsigned ii = mapper.global_to_bindex( globIdxAct( i ));

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const unsigned j = eltBdryFcts[j0];
                        const unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j

                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_ddof[unk_] = solver.compute( globProjMat ).solve ( globProjRhs );

} // computeDirichletDofsL2Proj

template<class T>
void gsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                       gsMultiPatch<T>& result, short_t unk) const
{
    // we might need to get a result even without having the system ..
    //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    // fixme: based on \a m_options and \a unk choose the right dof mapper
    const gsDofMapper & mapper = m_system.colMapper(unk);

    result.clear(); // result is cleared first

    /*
    GISMO_ASSERT(solVector.rows() == m_dofs,
                 "The provided solution vector does not match the system."
                 " Expected: "<<mapper.freeSize()<<", Got:"<<solVector.rows() );
    */
    // is the situation whtn solVector has more than one columns important?
    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[unk].cols() );

    // to do: test unknown_dim == dim

    gsMatrix<T> coeffs;
    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const size_t sz  = m_bases[m_system.colBasis(unk)][p].size();
        coeffs.resize(sz, dim);

        for (size_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof[unk].row( mapper.bindex(i, p) ).head(dim);
            }
        }

        result.addPatch( m_bases[m_system.colBasis(unk)][p].makeGeometry( give(coeffs) ) );
    }

    // AM: result topology ?
}

template<class T>
void gsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                       gsMultiPatch<T>& result,
                                       const gsVector<index_t> & unknowns) const
{
    // we might need to get a result even without having the system ..
    //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    GISMO_ASSERT(solVector.cols()==1, "Vector valued output only works for single rhs");
    index_t idx;

    const short_t dim = unknowns.rows();

    // fixme: based on \a m_options and \a unk choose the right dof mapper
    std::vector<gsDofMapper> mappers(dim);
    for(short_t unk = 0; unk<dim;++unk)
        mappers[unk] = m_system.colMapper(unknowns[unk]);

    result.clear(); // result is cleared first

    gsMatrix<T> coeffs;
    gsVector<index_t> basisIndices(dim);
    for(short_t unk = 0; unk<dim;++unk)
        basisIndices[unk] = m_system.colBasis(unknowns[unk]);

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        const size_t sz  = m_bases[basisIndices[0]][p].size(); //must be equal for all unk
        coeffs.resize(sz, dim);

        for(short_t unk = 0; unk<dim;++unk)
        {
            // Reconstruct solution coefficients on patch p
            for (size_t i = 0; i < sz; ++i)
            {
                if ( mappers[unk].is_free(i, p) ) // DoF value is in the solVector
                {
                    m_system.mapToGlobalColIndex(i,p,idx,unknowns[unk]);
                    coeffs(i,unk) = solVector(idx,0);
                    // coeffs(i,unk) = solVector(mappers[unk].index(i, p),0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,unk) = m_ddof[unknowns[unk]](mappers[unk].bindex(i, p),0);
                }
            }
        }
        result.addPatch( m_bases[basisIndices[0]][p].makeGeometry( give(coeffs) ) );
    }

    // AM: result topology ?
}


template<class T>
gsField<T> gsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                                short_t unk) const
{
    typename gsMultiPatch<T>::Ptr result(new gsMultiPatch<T>);
    constructSolution(solVector,*result, unk);
    return gsField<T>(m_pde_ptr->domain(), result, true);
}


//This silently assumes the same basis for all components
template<class T>
void gsAssembler<T>::updateSolution(const gsMatrix<T>& solVector,
                                    gsMultiPatch<T>& result, T theta) const
{
    // GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");
    index_t idx;

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        // Update solution coefficients on patch p
        const size_t sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < m_system.numColBlocks(); ++j)
        {
            const gsDofMapper & mapper = m_system.colMapper(j);
            for (size_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    m_system.mapToGlobalColIndex(i,p,idx,j);
                    coeffs(i,j) += theta * solVector(idx,0);
                }
            }
        }
    }
}

}// namespace gismo
