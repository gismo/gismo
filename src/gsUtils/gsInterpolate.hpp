//interpolation

# pragma once

#include <gsCore/gsGeometry.h>
#include <gsUtils/gsCollocationMatrix.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo
{


template<class T>
gsGeometry<T> * gsInterpolate( gsBasis<T> const& g, gsMatrix<T> const& pts, gsMatrix<T> const& vals )
{
  GISMO_ASSERT (g.dim()  == pts.rows() , "Wrong dimension of the points" );
  GISMO_ASSERT (g.size() == pts.cols() , "Expecting as many points as the basis functions." );
  GISMO_ASSERT (g.size() == vals.cols(), "Expecting as many values as the number of points." );

  gsSparseMatrix<T>  Cmat;
  gsCollocationMatrix_into( g , pts,Cmat );
  gsMatrix<T> x (g.size(), vals.rows());

  //Eigen::ConjugateGradient< gsSparseMatrix<T> > solver(Cmat);// for symmetric - does not apply here
  //Eigen::BiCGSTAB<gsSparseMatrix<T>, Eigen::DiagonalPreconditioner<T> > solver;
  //Eigen::BiCGSTAB<gsSparseMatrix<T>, Eigen::IdentityPreconditioner > solver;
  Eigen::BiCGSTAB< gsSparseMatrix<T>,  Eigen::IncompleteLUT<T> > solver( Cmat );

  // Solves for many right hand side  columns
  x =  solver.solve( vals.transpose() ); //toDense()
  
  // gsInfo <<"gs Interpolate error : " << solver.error() << std::"\n";
  // gsInfo <<"gs Interpolate iters : " << solver.iterations() << std::"\n";
  // gsInfo <<"intpl sol : " << x.transpose() << std::"\n";
  
  return g.makeGeometry( give(x) );
}


template<class T>
gsGeometry<T> * gsInterpolate( const gsBasis<T>& g, const gsFunction<T>& f )
{
    // Caution: not tested!
    gsMatrix<T> pts = g.anchors();
    gsMatrix<T> fpts;
    f.eval_into( pts, fpts);

    gsGeometry<T>* result = gsInterpolate<T>(g, pts, fpts );
    return result;
}


template <class T>
void gsInterpolateBoundary( const gsBasis<T> & basis,
                           const gsFunction<T> & f,
                           const gsGeometry<T> & geo,
                           const gsVector<int> & Sides,
                           gsVector<unsigned> & vecIdx,
                           gsMatrix<T> & vecCoeff,
                           bool getGlobalData)
{

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
        geo.evaluator(NEED_VALUE | NEED_MEASURE) );

    const int NB = basis.size();

    // Vector of flags indicating whether a given DOF is a
    // "true" boundary DOF or not.
    gsVector<unsigned> Flag_IsBdryGlob( NB );
    Flag_IsBdryGlob.setZero();

    gsMatrix< unsigned > act_Idx_raw;

    // Iterate through all given boundaries (Sides)
    // to get all true boundary DOF.
    std::auto_ptr< gsDomainIterator<T> > bdIt;
    for( int side_idx = 0; side_idx < Sides.size(); side_idx++)
    {
        int side = Sides[ side_idx ];
        bdIt = basis.makeDomainIterator(static_cast<boundary::side>(side));

        gsVector<int> numIntNodes( basis.dim() );
        numIntNodes.setOnes();
        bdIt->computeQuadratureRule(numIntNodes);

        // Iterate over all elements/cells of the side.
        for (; bdIt->good(); bdIt->next())
        {
            // Get all active functions. These, however, may include
            // functions which are zero at the boundary.
            act_Idx_raw = bdIt->activeFuncs;

            GISMO_ASSERT(act_Idx_raw.cols() == 1,
                         "Something wrong with bdIt->activeFuncs!");

            bdIt->evaluateBasis(0);

            // Test all these functions whether they are really nonzero at
            // the midpoint of the boundary element/cell.
            //
            // !!! WARNING !!!
            // This tests a variable of type T (maybe double or so)
            // for equality!
            for( int i=0; i < act_Idx_raw.size(); i++ )
                if( bdIt->basisValue( i, 0 ) != 0 )
                    Flag_IsBdryGlob[ act_Idx_raw( i, 0) ] = 1;
        }
    }

    // NSys will be the size of the global system which we will solve.
    int NSys = 0;
    for( int i=0; i < Flag_IsBdryGlob.size(); i++ )
        NSys += Flag_IsBdryGlob[i]; // Count the 1's in the flag-list.

    // Set up an inverse mapping
    // from the global index of the boundary DOF
    // to the place in the system to be solved.
    // "-1" corresponds to "not a boundary DOF".
    gsVector<int> Map_GlobToSys( NB );
    Map_GlobToSys.setConstant(-1);

    int ii = 0;
    for( int i=0; i < Flag_IsBdryGlob.size(); i++ )
        if( Flag_IsBdryGlob[i] == 1 )
        {
            Map_GlobToSys[i] = ii;
            ii += 1;
        }

    // Set up matrix and right-hand-side for the projection.
    gsSparseMatrix<T> gloA(NSys,NSys);
    gsMatrix<T> gloB(NSys,1);
    gloA.setZero();
    gloB.setZero();

    gsMatrix<T> locA;
    gsMatrix<T> locB;

    // Again loop over all sides in Sides.
    for( int side_idx = 0; side_idx < Sides.size(); side_idx++)
    {
        int side = Sides[ side_idx ];
        bdIt = basis.makeDomainIterator(static_cast<boundary::side>(side));

        // Manually set up the quadrature rule for the boundary domain iterator.
        gsVector<int> numIntNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numIntNodes[i] = basis.degree(i) + 1;
        int dir = direction(side);
        numIntNodes[dir] = 1;
        bdIt->computeQuadratureRule(numIntNodes);

        // Iterate over all elements/cells of the side.
        for (; bdIt->good(); bdIt->next())
        {
            // Get all active functions.
            act_Idx_raw = bdIt->activeFuncs;

            // Check whether these actives are truly active.
            std::vector<int> act_Idx_true;
            for( int i=0; i < act_Idx_raw.size(); i++)
                if( Flag_IsBdryGlob[ act_Idx_raw(i,0) ] == 1 )
                    act_Idx_true.push_back(i);

            // act_Idx_true is now a mapping
            // from set of truly active functions
            // to the local indices of all active functions.
            int locN = act_Idx_true.size();

            // Set up computation of local contributions to the global system.
            locA.resize( locN, locN );
            locB.resize( locN, 1 );
            locA.setZero();
            locB.setZero();

            // Evaluate the basis functions and the geometry mapping.
            bdIt->evaluateBasis(0);
            geoEval->evaluateAt(bdIt->quNodes);

            // Evaluate the function that should be interpolated.
            gsMatrix<T> f_Vals;
            f.eval_into( geoEval->values(), f_Vals );

            // Get the truly active basis functions.
            typename gsMatrix<T>::Block BasisFcts_raw = bdIt->basisValues();
            gsMatrix<T> BasisFcts_true( locN, BasisFcts_raw.cols() );
            for( unsigned i=0; i < act_Idx_true.size(); i++ )
                BasisFcts_true.row( i ) = BasisFcts_raw.row( act_Idx_true[i] );

            gsMatrix<T> BasisFcts_pp;
            // Finally, loop over quadrature points and compute the local contributions.
            for( int pp = 0; pp < BasisFcts_true.cols(); pp++)
            {
                const T weight = bdIt->quWeights[pp] * geoEval->measure(pp);
                BasisFcts_pp = BasisFcts_true.col(pp);

                locA.noalias() += weight * ( BasisFcts_pp * BasisFcts_pp.transpose() );
                locB.noalias() += BasisFcts_pp * ( weight * f_Vals.col(pp) ).transpose();
            }

            // Put the local contributions to the corresponding places
            // in the global system.
            int ii;
            int jj;
            for( int i=0; i < locN; i++)
            {
                // ii is the global index of the i-th truly active function.
                ii = act_Idx_raw( act_Idx_true[i], 0);
                for( int j=0; j < locN; j++)
                {
                    jj = act_Idx_raw( act_Idx_true[j], 0);
                    gloA.coeffRef( Map_GlobToSys[ii], Map_GlobToSys[jj] ) += locA(i,j);
                    // ...Map_GlobToSys[ii] is the place in the system matrix gloA
                    // corresponding to the function with global index ii.
                }
                gloB( Map_GlobToSys[ii], 0 ) += locB( i, 0 );
            }
        } // bdIt
    } // for side

    gloA.makeCompressed();

    Eigen::ConjugateGradient< gsSparseMatrix<T> > solver;
    //Eigen::Matrix<T,Dynamic,Dynamic> gloX;
    vecCoeff = solver.compute( gloA ).solve ( gloB );

    if( getGlobalData )
    {
        // If a global vector should be returned...
        // ...the vector with the indices of the DOF is simply the "identity", ...
        vecIdx.resize( NB );
        // ... and the computed coefficients will have to be put in the
        // correct positions in a global vector. The rest of the global
        // vector are set to zeros.
        gsMatrix<T> vecCoeff_glob( NB , 1);
        vecCoeff_glob.setZero();

        for( int i=0; i < NB; i++)
        {
            vecIdx[i] = i; // The "identity" vector.

            // Put the coefficients in the correct places.
            if( Flag_IsBdryGlob[i] == 1 )
                vecCoeff_glob( i ,0 ) = vecCoeff( Map_GlobToSys[ i ], 0 );
        }
        vecCoeff.resize( NB, 1 );
        vecCoeff = vecCoeff_glob;
    }
    else
    {
        // If only the "relevant" vectors of coefficients and coefficient values
        // should be returned, put all the indices of the "true" boundary DOFs
        // into the vector vecIdx.
        vecIdx.resize( NSys );
        for( int i=0; i < Flag_IsBdryGlob.size(); i++)
            if( Flag_IsBdryGlob[i] == 1 )
                vecIdx[ Map_GlobToSys[ i ] ] = i;
    }

} // gsInterpolateBoundary





};// namespace gismo





