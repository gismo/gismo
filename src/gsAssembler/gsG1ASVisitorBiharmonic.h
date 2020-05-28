/** @file gsG1ASVisitorBiharmonic.h

    @brief Visitor for a simple Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

#pragma once

namespace gismo
{

/** \brief Visitor for the biharmonic equation.
 *
 * Assembles the bilinear terms
 * \f[ (\Delta u,\Delta v)_\Omega \text{ and } (f,v)_\Omega \f]
 * For \f[ u = g \quad on \quad \partial \Omega \f],
 *
 */

template <class T>
class gsG1ASVisitorBiharmonic
{
public:

    gsG1ASVisitorBiharmonic(const gsPde<T> & pde)
    {
        rhs_ptr = static_cast<const gsBiharmonicPde<T>&>(pde).rhs() ;
    }

    /** \brief Constructor for gsG1ASVisitorBiharmonic.
     *
     * \param[in] rhs Given right-hand-side function/source term that, for
     */
    gsG1ASVisitorBiharmonic(const gsFunction<T> & rhs) :
        rhs_ptr(&rhs)
    {
        GISMO_ASSERT( rhs.targetDim() == 1 ,"Not yet tested for multiple right-hand-sides");
    }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;

    }

    void initialize(const gsBasis<T> & basis,
                    const index_t ,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule)
    {
        // Setup Quadrature
        rule = gsQuadrature::get(basis, options);

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
                         const gsGeometry<T>    & geo,
                         gsMatrix<T>            & quNodes)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);
        numActive = actives.rows();

        //deriv2_into()
        //col(point) = B1_xx B2_yy B1_zz B_xy B1_xz B1_xy B2_xx ...

        // Evaluate basis functions on element
        basis.evalAllDers_into(md.points, 2, basisData);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geo.computeMap(md);

        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into(md.values[0], rhsVals); // Dim: 1 X NumPts

        if(md.dim.first +1 == md.dim.second)
        {
            gsMatrix<T> geoMapDeriv1 = geo.deriv(md.points); // First derivative of the geometric mapping with respect to the parameter coordinates
            gsMatrix<T> geoMapDeriv2 = geo.deriv2(md.points); // Second derivative of the geometric mapping with respect to the parameter coordinates

//            FIRST FUNDAMENTAL FORM: G = J^T * J
//
//            G = | G11   G12|
//                | G21   G22|
//
//            INVERSE OF THE FIRST FUNDAMENTAL FORM
//
//                      1    | G22  -G12|      1
//            G^-1 = ------- |          | = ------- G* ^-1
//                    det(G) | -G21  G11|    det(G)

            // First fundamental form
            gsMatrix<T> G11 = ( geoMapDeriv1.row(0).cwiseProduct(geoMapDeriv1.row(0)) +
                                geoMapDeriv1.row(2).cwiseProduct(geoMapDeriv1.row(2)) +
                                geoMapDeriv1.row(4).cwiseProduct(geoMapDeriv1.row(4)));

//          G12 = G21
            gsMatrix<T> G12 = ( geoMapDeriv1.row(0).cwiseProduct(geoMapDeriv1.row(1)) +
                                geoMapDeriv1.row(2).cwiseProduct(geoMapDeriv1.row(3)) +
                                geoMapDeriv1.row(4).cwiseProduct(geoMapDeriv1.row(5)));

            gsMatrix<T> G22 = ( geoMapDeriv1.row(1).cwiseProduct(geoMapDeriv1.row(1)) +
                                geoMapDeriv1.row(3).cwiseProduct(geoMapDeriv1.row(3)) +
                                geoMapDeriv1.row(5).cwiseProduct(geoMapDeriv1.row(5)));


            // Derivative of the first fundamental form
            gsMatrix<T> DuG11 = 2 * ( geoMapDeriv2.row(0).cwiseProduct(geoMapDeriv1.row(0)) +
                                      geoMapDeriv2.row(3).cwiseProduct(geoMapDeriv1.row(2)) +
                                      geoMapDeriv2.row(6).cwiseProduct(geoMapDeriv1.row(4)) );

            gsMatrix<T> DvG11 = 2 * ( geoMapDeriv2.row(2).cwiseProduct(geoMapDeriv1.row(0)) +
                                      geoMapDeriv2.row(5).cwiseProduct(geoMapDeriv1.row(2)) +
                                      geoMapDeriv2.row(8).cwiseProduct(geoMapDeriv1.row(4)) );

//          DuG12 = DuG21
            gsMatrix<T> DuG12 = (     geoMapDeriv2.row(0).cwiseProduct(geoMapDeriv1.row(1)) +
                                      geoMapDeriv2.row(2).cwiseProduct(geoMapDeriv1.row(0)) +
                                      geoMapDeriv2.row(3).cwiseProduct(geoMapDeriv1.row(3)) +
                                      geoMapDeriv2.row(5).cwiseProduct(geoMapDeriv1.row(2)) +
                                      geoMapDeriv2.row(6).cwiseProduct(geoMapDeriv1.row(5)) +
                                      geoMapDeriv2.row(8).cwiseProduct(geoMapDeriv1.row(4)) );

//          DvG12 = DvG21
            gsMatrix<T> DvG21 = (     geoMapDeriv2.row(2).cwiseProduct(geoMapDeriv1.row(1)) +
                                      geoMapDeriv2.row(1).cwiseProduct(geoMapDeriv1.row(0)) +
                                      geoMapDeriv2.row(5).cwiseProduct(geoMapDeriv1.row(3)) +
                                      geoMapDeriv2.row(4).cwiseProduct(geoMapDeriv1.row(2)) +
                                      geoMapDeriv2.row(8).cwiseProduct(geoMapDeriv1.row(5)) +
                                      geoMapDeriv2.row(7).cwiseProduct(geoMapDeriv1.row(4)) );


            gsMatrix<T> DuG22 = 2 * ( geoMapDeriv2.row(2).cwiseProduct(geoMapDeriv1.row(1)) +
                                      geoMapDeriv2.row(5).cwiseProduct(geoMapDeriv1.row(3)) +
                                      geoMapDeriv2.row(8).cwiseProduct(geoMapDeriv1.row(5)) );

            gsMatrix<T> DvG22 = 2 * ( geoMapDeriv2.row(1).cwiseProduct(geoMapDeriv1.row(1)) +
                                      geoMapDeriv2.row(4).cwiseProduct(geoMapDeriv1.row(3)) +
                                      geoMapDeriv2.row(7).cwiseProduct(geoMapDeriv1.row(5)) );

            gsMatrix<T> detG = G11.cwiseProduct(G22) - G12.cwiseProduct(G12);

//          1 / sqrt( det( G ) )
            gsMatrix<T> sqrtDetG_inv;
            sqrtDetG_inv.resize(1, numActive);

//          1 / ( 2 * det( G )^( 3/2 ) )
            gsMatrix<T> sqrtDetG_inv_derivative;
            sqrtDetG_inv_derivative.resize(1, numActive);

            for(index_t k = 0; k < md.points.cols(); k++)
            {
                sqrtDetG_inv(0, k) = 1 / sqrt(detG(0, k));
                sqrtDetG_inv_derivative(0, k) = 1 / (2 * detG(0, k) * sqrt(detG(0, k)) );
            }

            gsMatrix<T> & basisGrads = basisData[1];
            gsMatrix<T> & basis2ndDerivs = basisData[2];

//          div ( sqrt( det( G ) ) * ( 1 / det( G ) * G* ^-1 * grad( u ) ) )
            surfParametricLaplace.resize(numActive, md.points.cols());

            for(index_t i = 0; i < numActive - 1; i++)
            {
//              Du( 1 / sqrt( det( G ) ) ) * G* ^-1 * grad( u ) ) +
//              Dv( 1 / sqrt( det( G ) ) ) * G* ^-1 * grad( u ) ) +
//              1 / sqrt( det( G ) ) * ( div ( G* ^-1 * grad( u ) ) )
                surfParametricLaplace.row(i) = sqrtDetG_inv_derivative.cwiseProduct(
                                               2 * G12.cwiseProduct( DuG12 ) -
                                               G22.cwiseProduct( DuG11 ) -
                                               G11.cwiseProduct( DuG22 ) ).cwiseProduct(
                                               G22.cwiseProduct( basisGrads.row( i * 2 ) ) -
                                               G12.cwiseProduct( basisGrads.row( i * 2 + 1 ) ) )
                                               +
                                               sqrtDetG_inv_derivative.cwiseProduct(
                                               2 * G12.cwiseProduct( DvG21 ) -
                                               G22.cwiseProduct( DvG11 ) -
                                               G11.cwiseProduct( DvG22 ) ).cwiseProduct(
                                               G11.cwiseProduct( basisGrads.row( i * 2 + 1 ) ) -
                                               G12.cwiseProduct( basisGrads.row( i * 2 ) ) )
                                               +
                                               sqrtDetG_inv.cwiseProduct(
                                               DuG22.cwiseProduct( basisGrads.row( i * 2 ) ) +
                                               G22.cwiseProduct( basis2ndDerivs.row( i * 3 ) ) -
                                               DuG12.cwiseProduct( basisGrads.row( i * 2 + 1 ) ) -
                                               G12.cwiseProduct( basis2ndDerivs.row( i * 3 + 2 ) ) +
                                               DvG11.cwiseProduct( basisGrads.row( i *2 + 1 ) ) +
                                               G11.cwiseProduct( basis2ndDerivs.row( i * 2 + 1 ) ) -
                                               DvG21.cwiseProduct( basisGrads.row( i * 2 ) ) -
                                               G12.cwiseProduct( basis2ndDerivs.row( i * 3 + 2 ) ) );
            }
            rhsVals = rhsVals.cwiseProduct( detG.cwiseProduct( sqrtDetG_inv ) );
        }

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
    }


    inline void assemble(gsDomainIterator<T>    & ,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData[0];
        gsMatrix<T> & basisGrads = basisData[1];
        gsMatrix<T> & basis2ndDerivs = basisData[2];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * md.measure(k);

            // Compute physical laplacian at k as a 1 x numActive matrix
            if(md.dim.first == md.dim.second)
            {
                transformLaplaceHgrad(md, k, basisGrads, basis2ndDerivs, physBasisLaplace);

                localMat.noalias() += weight * (physBasisLaplace.transpose() * physBasisLaplace);
                localRhs.noalias() += weight * ( basisVals.col(k) * rhsVals.col(k).transpose() ) ;

            }
            else
            if(md.dim.first + 1 == md.dim.second)
            {
                localMat.noalias() += weight * (surfParametricLaplace.col(k) * surfParametricLaplace.col(k).transpose());
                localRhs.noalias() += weight * ( basisVals.col(k) * rhsVals.col(k).transpose() ) ;
            }

        }

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs[0], 0, 0);
    }

    /*
    inline void localToGlobal(const gsDofMapper     & mapper,
                              const gsMatrix<T>     & eliminatedDofs,
                              const int patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);
        //const int numActive = actives.rows();

        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = actives(i);
            if ( mapper.is_free_index(ii) )
            {
                rhsMatrix.row(ii) += localRhs.row(i);

                for (index_t j=0; j < numActive; ++j)
                {
                    const int jj = actives(j);
                    if ( mapper.is_free_index(jj) )
                    {
                        sysMatrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else
                    {
                        rhsMatrix.row(ii).noalias() -= localMat(i, j) *
                            eliminatedDofs.row( mapper.global_to_bindex(jj) );
                    }
                }
            }
        }
    }
    */


protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;

protected:
    // Basis values

    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>        physBasisLaplace;
    gsMatrix<T>        surfParametricLaplace;

    gsMatrix<unsigned> actives;
    index_t numActive;


protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;
};


} // namespace gismo

