/** @file gsPatchPreconditionersCreator.hpp

    @brief Provides preconditioners that live on the parameter domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/
#pragma once

#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsKroneckerOp.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

namespace {

template<typename T>
gsBoundaryConditions<T> boundaryConditionsForDirection( const gsBoundaryConditions<T>& bc, index_t direction )
{
    gsBoundaryConditions<T> result;

    for ( index_t i = 1; i <= 2; ++i)
    {
        patchSide global(0,i+2*direction), local(0,i);
        const boundary_condition<T>* cond = bc.getConditionFromSide(global);
        if (cond)
            result.addCondition(local,cond->type(),cond->function());
    }
    return result;
}

template<typename T>
gsBoundaryConditions<T> boundaryConditionsForComponent( const gsBoundaryConditions<T>& bc, short_t component )
{
    gsBoundaryConditions<T> result;

    for ( typename gsBoundaryConditions<T>::const_bciterator it = bc.beginAll(); it != bc.endAll(); ++it )
        for ( typename gsBoundaryConditions<T>::const_iterator cond = it->second.begin(); cond != it->second.end(); ++cond )
            if ( cond->unknown() == 0 && cond->patch() == 0
                 && ( cond->unkComponent() == -1 || cond->unkComponent() == component )
                 && !result.getConditionFromSide(cond->ps) )
                result.addCondition(cond->ps,cond->type(),cond->function());
    return result;
}

template<typename T>
void eliminateDirichlet1D(const gsBoundaryConditions<T>& bc,
                          const gsOptionList& opt, gsSparseMatrix<T> & result)
{
    dirichlet::strategy ds = (dirichlet::strategy)opt.askInt("DirichletStrategy",dirichlet::elimination);
    if (ds == dirichlet::elimination)
    {
        patchSide west(0,boundary::west), east(0,boundary::east);
        index_t i = 0;
        if (bc.getConditionFromSide( west ) && bc.getConditionFromSide( west )->type() == condition_type::dirichlet ) i += 1;
        if (bc.getConditionFromSide( east ) && bc.getConditionFromSide( east )->type() == condition_type::dirichlet ) i += 2;
        if (i%2 + i/2 >= result.rows() || i%2 + i/2 >= result.cols())
            result.resize(0,0);
        else switch ( i )
             {
             case 0: break;
             case 1: result = result.block( 1, 1, result.rows()-1, result.cols()-1 ); break;
             case 2: result = result.block( 0, 0, result.rows()-1, result.cols()-1 ); break;
             case 3: result = result.block( 1, 1, result.rows()-2, result.cols()-2 ); break;
             }
    }
    else
        GISMO_ERROR("Unknown Dirichlet strategy.");
}

template<typename T>
gsSparseMatrix<T> assembleMass(const gsBasis<T>& basis)
{
    gsExprAssembler<T> mass(1,1);
    mass.setIntegrationDomain(basis.domain());
    typename gsExprAssembler<T>::space u = mass.getSpace(basis);
    mass.initMatrix();
    mass.assemble( u * u.tr() );
    gsSparseMatrix<T> result;
    mass.matrix_into(result);
    return result;
}

template<typename T>
gsSparseMatrix<T> assembleStiffness(const gsBasis<T>& basis)
{
    gsExprAssembler<T> stiff(1,1);
    stiff.setIntegrationDomain(basis.domain());
    typename gsExprAssembler<T>::space u = stiff.getSpace(basis);
    stiff.initMatrix();
    stiff.assemble( gismo::expr::grad(u) * gismo::expr::grad(u).tr() );
    gsSparseMatrix<T> result;
    stiff.matrix_into(result);
    return result;
}


// Quadrature nodes and weights of all elements of a univariate basis, with
// the values, the derivatives and the indices of the basis functions that
// are active at the nodes (one column per node)
template<typename T>
struct univariateQuadrature
{
    gsMatrix<T>       nodes;
    gsVector<T>       weights;
    gsMatrix<T>       values;
    gsMatrix<T>       derivs;
    gsMatrix<index_t> actives;
    T                 volume; // sum of the weights
};

template<typename T>
void setupUnivariateQuadrature(const gsBasis<T>& basis, const gsOptionList& opt,
                               univariateQuadrature<T>& result)
{
    typename gsQuadRule<T>::uPtr rule = gsQuadrature::getPtr(basis, opt);
    memory::shared_ptr< gsDomain<T> > domain = basis.domain();
    typename gsBasis<T>::domainIter domIt    = domain->beginAll();
    typename gsBasis<T>::domainIter domItEnd = domain->endAll();

    std::vector<T> nodes, weights;
    gsMatrix<T> elNodes;
    gsVector<T> elWeights;
    for (; domIt < domItEnd; ++domIt )
    {
        rule->mapTo( domIt.lowerCorner(), domIt.upperCorner(), elNodes, elWeights );
        for ( index_t q=0; q<elWeights.rows(); ++q )
        {
            nodes.push_back( elNodes(0,q) );
            weights.push_back( elWeights[q] );
        }
    }

    const index_t n = (index_t)nodes.size();
    result.nodes.resize(1,n);
    result.weights.resize(n);
    for ( index_t q=0; q<n; ++q )
    {
        result.nodes(0,q)  = nodes[q];
        result.weights[q]  = weights[q];
    }
    result.volume = result.weights.sum();
    basis.eval_into  (result.nodes, result.values );
    basis.deriv_into (result.nodes, result.derivs );
    basis.active_into(result.nodes, result.actives);
}

// Assembles the univariate matrix with the entries
// sum_q weights[q] coef[q] f_i(nodes[q]) f_j(nodes[q]), where the f_i are
// the derivatives of the basis functions if stiffness is true and their
// values otherwise
template<typename T>
gsSparseMatrix<T> assembleWeighted(const univariateQuadrature<T>& qu,
                                   const gsVector<T>& coef, bool stiffness, index_t sz)
{
    const gsMatrix<T>& f = stiffness ? qu.derivs : qu.values;
    const index_t K = f.rows();
    gsSparseEntries<T> entries;
    entries.reserve( K*K*f.cols() );
    for ( index_t q=0; q<f.cols(); ++q )
    {
        const T w = qu.weights[q] * coef[q];
        for ( index_t i=0; i<K; ++i )
            for ( index_t j=0; j<K; ++j )
                entries.add( qu.actives(i,q), qu.actives(j,q), w * f(i,q) * f(j,q) );
    }
    gsSparseMatrix<T> result(sz,sz);
    result.setFrom(entries);
    result.makeCompressed();
    return result;
}

// Averages of the diagonal of |det J| J^{-1} J^{-T} over the anchors of the
// basis of the geometry
template<typename T>
gsVector<T> averageMetricDiagonal(const gsGeometry<T>& geo)
{
    const index_t d = geo.domainDim();
    gsMatrix<T> pts, jacs;
    geo.basis().anchors_into(pts);
    geo.jacobian_into(pts, jacs);
    gsVector<T> result = gsVector<T>::Zero(d);
    for ( index_t i=0; i<pts.cols(); ++i )
    {
        const gsMatrix<T> J = jacs.block(0, i*d, d, d);
        const T detJ = math::abs(J.determinant());
        const gsMatrix<T> Jinv = J.inverse();
        for ( index_t k=0; k<d; ++k )
            result[k] += detJ * Jinv.row(k).squaredNorm();
    }
    result /= (T)(pts.cols());
    return result;
}

// Sums of |det J| (J^{-1})_{kc}^2 over the tensor grid of the quadrature
// nodes, where J is the Jacobian of the geometry. For the direction k, the
// component c and the direction j, result[(k*d+c)*d+j] has one entry per
// node of the direction j: the sum over the nodes of all other directions,
// weighted with their quadrature weights. The grid is traversed line by
// line, so the memory is linear in the univariate sizes.
//
// Returns false if the determinant vanishes or changes its sign.
template<typename T>
bool assembleMetricSums(const gsGeometry<T>& geo,
                        const std::vector< univariateQuadrature<T> >& qu,
                        std::vector< gsVector<T> >& result)
{
    const index_t d = (index_t)qu.size();
    result.resize(d*d*d);
    for ( index_t i=0; i<d*d; ++i )
        for ( index_t j=0; j<d; ++j )
            result[i*d+j].setZero( qu[j].weights.rows() );

    const index_t n0 = qu[0].weights.rows();
    gsMatrix<T> pts(d,n0), jacs, J, Jinv, line(d,d);
    pts.row(0) = qu[0].nodes;
    gsVector<index_t> idx = gsVector<index_t>::Zero(d);
    T sign = 0;

    // Loop over the lines in the direction 0
    while (true)
    {
        T w = 1;
        for ( index_t l=1; l<d; ++l )
        {
            pts.row(l).setConstant( qu[l].nodes(0,idx[l]) );
            w *= qu[l].weights[idx[l]];
        }
        geo.jacobian_into(pts, jacs);

        line.setZero();
        for ( index_t i=0; i<n0; ++i )
        {
            J = jacs.block(0, i*d, d, d);
            const T detJ = J.determinant();
            if (sign == 0)
                sign = detJ > 0 ? 1 : -1;
            if (!(sign * detJ > 0))
                return false;
            Jinv = J.inverse();
            for ( index_t k=0; k<d; ++k )
                for ( index_t c=0; c<d; ++c )
                {
                    const T val = sign * detJ * Jinv(k,c) * Jinv(k,c);
                    result[(k*d+c)*d][i] += w * val;
                    line(k,c) += qu[0].weights[i] * val;
                }
        }

        for ( index_t j=1; j<d; ++j )
        {
            T wj = 1;
            for ( index_t l=1; l<d; ++l )
                if (l != j)
                    wj *= qu[l].weights[idx[l]];
            for ( index_t k=0; k<d; ++k )
                for ( index_t c=0; c<d; ++c )
                    result[(k*d+c)*d+j][idx[j]] += wj * line(k,c);
        }

        // Next line
        index_t l = 1;
        for ( ; l<d; ++l )
        {
            if (++idx[l] < qu[l].weights.rows())
                break;
            idx[l] = 0;
        }
        if (l == d)
            break;
    }
    return true;
}

// Univariate matrices of the separable approximation
//   c_k(xi) ~ a_k(xi_k) prod_{j!=k} m_j(xi_j)
// of the coefficients c_k = sum_c factors(k,c) |det J| (J^{-1})_{kc}^2,
// given the sums from assembleMetricSums: stiffness matrices weighted with
// a_k and mass matrices weighted with m_k.
//
// a_k is the mean of c_k over the hyperplanes orthogonal to the direction k.
// m_j is the geometric mean over k!=j of the means of c_k over the
// hyperplanes orthogonal to the direction j, each divided by the mean of c_k
// over the whole parameter domain. All means are weighted with the
// quadrature weights. The approximation is exact if each c_k is of the
// given form; for a constant Jacobian, a_k = c_k and m_j = 1.
template<typename T>
void assembleTensorSeparable(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    const std::vector< univariateQuadrature<T> >& qu,
    const std::vector< gsVector<T> >& sums,
    const gsMatrix<T>& factors,
    std::vector< gsSparseMatrix<T> >& local_stiff,
    std::vector< gsSparseMatrix<T> >& local_mass
    )
{
    const index_t d = basis.dim();

    // mean[k*d+j]: means of c_k over the hyperplanes orthogonal to j
    std::vector< gsVector<T> > mean(d*d);
    gsVector<T> total(d);
    for ( index_t k=0; k<d; ++k )
    {
        for ( index_t j=0; j<d; ++j )
        {
            T vol = 1;
            for ( index_t l=0; l<d; ++l )
                if (l != j)
                    vol *= qu[l].volume;
            mean[k*d+j].setZero( qu[j].weights.rows() );
            for ( index_t c=0; c<d; ++c )
                mean[k*d+j] += (factors(k,c) / vol) * sums[(k*d+c)*d+j];
        }
        total[k] = qu[k].weights.dot(mean[k*d+k]) / qu[k].volume;
    }

    local_stiff.resize(d);
    local_mass.resize(d);
    for ( index_t j=0; j<d; ++j )
    {
        const gsVector<T>& a = mean[j*d+j];
        gsVector<T> m = gsVector<T>::Ones( qu[j].weights.rows() );
        for ( index_t k=0; k<d; ++k )
            if (k != j)
                m.array() *= mean[k*d+j].array() / total[k];
        if (d > 2)
            for ( index_t q=0; q<m.rows(); ++q )
                m[q] = math::pow( m[q], (T)(1)/(T)(d-1) );

        GISMO_ASSERT ( (a.array() > 0).all() && (m.array() > 0).all(),
            "gsPatchPreconditionersCreator: the separable weights have to be positive." );

        // The i-th univariate matrix belongs to the direction d-1-i
        const index_t sz = basis.component(j).size();
        local_stiff[d-1-j] = assembleWeighted(qu[j], a, true , sz);
        local_mass [d-1-j] = assembleWeighted(qu[j], m, false, sz);
        const gsBoundaryConditions<T> bcj = boundaryConditionsForDirection(bc,j);
        eliminateDirichlet1D(bcj, opt, local_stiff[d-1-j]);
        eliminateDirichlet1D(bcj, opt, local_mass [d-1-j]);
    }
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorMass(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();
    std::vector< gsSparseMatrix<T> > result(d);
    for ( index_t i=0; i!=d; ++i )
    {
        result[i] = assembleMass(basis.component(d-1-i));
        eliminateDirichlet1D(boundaryConditionsForDirection(bc,d-1-i), opt, result[i]);
    }
    return result;
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorStiffness(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();
    std::vector< gsSparseMatrix<T> > result(d);
    for ( index_t i=0; i!=d; ++i )
    {
        result[i] = assembleStiffness(basis.component(d-1-i));
        eliminateDirichlet1D(boundaryConditionsForDirection(bc,d-1-i), opt, result[i]);
    }
    return result;
}

} // anonymous namespace

template<typename T>
gsSparseMatrix<T> gsPatchPreconditionersCreator<T>::massMatrix(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(basis, bc, opt);
    const index_t d = local_mass.size();
    gsSparseMatrix<T> result = local_mass[d-1];
    for (index_t i=d-2; i>-1; --i)
        result = local_mass[i].kron(result);
    return result;
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::massMatrixOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(basis, bc, opt);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeMatrixOp(local_mass[i].moveToPtr());

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::massMatrixInvOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_mass = assembleTensorMass(basis, bc, opt);

    std::vector<OpPtr> local_mass_op(d);
    for (index_t i=0; i<d; ++i)
        local_mass_op[i] = makeSparseCholeskySolver(local_mass[i]);

    return gsKroneckerOp<T>::make(local_mass_op);
}

template<typename T>
gsSparseMatrix<T> gsPatchPreconditionersCreator<T>::stiffnessMatrix(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T alpha,
    T beta
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

    if ( beta!=1 )
    {
        for (index_t i=0; i!=d; ++i)
            local_stiff[i] *= beta;
    }

    gsSparseMatrix<T> K = give(local_stiff[d-1]);
    gsSparseMatrix<T> M = give(local_mass [d-1]);

    for (index_t i=d-2; i>-1; --i)
    {
        K  = local_mass[i].kron(K);
        K += local_stiff[i].kron(M);
        if ( i != 0 || alpha != 0 )
            M = local_mass[i].kron(M);
    }
    if (alpha==1)
        K += M;
    else if (alpha!=0)
        K += alpha * M;
    return K;
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::stiffnessMatrixOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T alpha,
    T beta
    )
{
    const index_t d = basis.dim();

    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

    std::vector<OpUPtr> local_stiff_op(d);
    std::vector<OpPtr > local_mass_op (d);
    for (index_t i=0; i<d; ++i)
    {
        if (beta!=1)
            local_stiff[i] *= beta;
        local_stiff_op[i] = makeMatrixOp(local_stiff[i].moveToPtr());
        local_mass_op [i] = makeMatrixOp(local_mass [i].moveToPtr());
    }
    OpUPtr K = give(local_stiff_op[0]);
    OpPtr  M = give(local_mass_op [0]);

    for (index_t i=1; i<d; ++i)
    {
        K = gsSumOp<T>::make(
            gsKroneckerOp<T>::make( give(K),      local_mass_op [i]  ),
            gsKroneckerOp<T>::make( M,       give(local_stiff_op[i]) )
            );
        if ( i < d-1 || alpha != 0 )
            M = gsKroneckerOp<T>::make(M, local_mass_op[i]);
    }
    if (alpha==1)
        K = gsSumOp<T>::make( give(K), M );
    else if (alpha!=0)
        K = gsSumOp<T>::make( give(K), gsScaledOp<T>::make( M, alpha ) );

    return K;
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T alpha,
    T beta,
    T gamma
    )
{
    GISMO_ASSERT ( beta != 0, "gsPatchPreconditionersCreator<T>::fastDiagonalizationOp() does not work for beta==0." );

    return fastDiagonalization(basis, bc, opt, alpha, gsVector<T>::Constant(basis.dim(),beta), gamma);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T alpha,
    const gsVector<T>& beta
    )
{
    GISMO_ASSERT ( beta.rows() == basis.dim(),
        "gsPatchPreconditionersCreator<T>::fastDiagonalizationOp() needs one beta per direction." );
    GISMO_ASSERT ( (beta.array() > 0).all(),
        "gsPatchPreconditionersCreator<T>::fastDiagonalizationOp() needs positive betas." );

    return fastDiagonalization(basis, bc, opt, alpha, beta, 0);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(
    const gsBasis<T>& basis,
    const gsGeometry<T>& geo,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();
    GISMO_ASSERT ( geo.domainDim() == d && geo.targetDim() == d,
        "gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(): geometry and basis do not fit." );

    std::vector< univariateQuadrature<T> > qu(d);
    for ( index_t k=0; k<d; ++k )
        setupUnivariateQuadrature(basis.component(k), opt, qu[k]);

    std::vector< gsVector<T> > sums;
    if ( !assembleMetricSums(geo, qu, sums) )
    {
        // The geometry map is not regular: one average per direction
        return fastDiagonalization(basis, bc, opt, 0, averageMetricDiagonal(geo), 0);
    }

    std::vector< gsSparseMatrix<T> > local_stiff, local_mass;
    const gsMatrix<T> factors = gsMatrix<T>::Ones(d,d);
    assembleTensorSeparable(basis, bc, opt, qu, sums, factors, local_stiff, local_mass);
    return fastDiagonalization(local_stiff, local_mass, 0, 0);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T lambda,
    T mu,
    const gsGeometry<T>* geo
    )
{
    std::vector< gsBoundaryConditions<T> > bcs;
    for ( index_t c=0; c<basis.dim(); ++c )
        bcs.push_back( boundaryConditionsForComponent(bc,(short_t)c) );
    return elasticityFastDiagonalizationOp(basis, bcs, opt, lambda, mu, geo);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp(
    const gsBasis<T>& basis,
    const std::vector< gsBoundaryConditions<T> >& bc,
    const gsOptionList& opt,
    T lambda,
    T mu,
    const gsGeometry<T>* geo
    )
{
    GISMO_ASSERT ( (index_t)bc.size() == basis.dim(),
        "gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp() needs boundary conditions for each component." );
    GISMO_ASSERT ( mu > 0 && lambda + mu > 0,
        "gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp() needs mu>0 and lambda+mu>0." );

    const index_t d = basis.dim();
    typename gsBlockOp<T>::uPtr result = gsBlockOp<T>::make(d,d);

    gsVector<T> scaling = gsVector<T>::Ones(d);
    if (geo)
    {
        GISMO_ASSERT ( geo->domainDim() == d && geo->targetDim() == d,
            "gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp(): geometry and basis do not fit." );

        std::vector< univariateQuadrature<T> > qu(d);
        for ( index_t k=0; k<d; ++k )
            setupUnivariateQuadrature(basis.component(k), opt, qu[k]);

        std::vector< gsVector<T> > sums;
        if ( assembleMetricSums(*geo, qu, sums) )
        {
            // The coefficient of the direction k in the block of the component c is
            // mu (|det J| J^{-1} J^{-T})_{kk} + (lambda+mu) |det J| (J^{-1})_{kc}^2
            std::vector< gsSparseMatrix<T> > local_stiff, local_mass;
            for ( index_t c=0; c<d; ++c )
            {
                gsMatrix<T> factors = gsMatrix<T>::Constant(d,d,mu);
                factors.col(c).array() += lambda+mu;
                assembleTensorSeparable(basis, bc[c], opt, qu, sums, factors, local_stiff, local_mass);
                result->addOperator(c, c, fastDiagonalization(local_stiff, local_mass, 0, 0) );
            }
            return result;
        }

        // The geometry map is not regular: one average per direction
        scaling = averageMetricDiagonal(*geo);
    }

    for ( index_t c=0; c<d; ++c )
    {
        gsVector<T> beta = mu * scaling;
        beta[c] += (lambda+mu) * scaling[c];
        result->addOperator(c, c,
            fastDiagonalization(basis, bc[c], opt, 0, beta, 0) );
    }
    return result;
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::fastDiagonalization(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T alpha,
    const gsVector<T>& beta,
    T gamma
    )
{
    const index_t d = basis.dim();

    // Assemble univariate
    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

    // The i-th univariate matrix belongs to the direction d-1-i
    for (index_t i=0; i!=d; ++i)
        local_stiff[i] *= beta[d-1-i];

    return fastDiagonalization(local_stiff, local_mass, alpha, gamma);
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::fastDiagonalization(
    const std::vector< gsSparseMatrix<T> >& local_stiff,
    const std::vector< gsSparseMatrix<T> >& local_mass,
    T alpha,
    T gamma
    )
{
    const index_t d = (index_t)local_stiff.size();

    // Determine overall size
    index_t sz = 1;
    for ( index_t i=0; i<d; ++i )
        sz *= local_stiff[i].rows();

    // Initialize the diagonal with 1
    gsMatrix<T> diag;
    diag.setConstant(sz,1,alpha); // This is the pure-mass part!

    index_t glob = sz; // Indexing value for setting up the Kronecker product

    typedef typename gsMatrix<T>::GenSelfAdjEigenSolver EVSolver;
    typedef typename EVSolver::EigenvectorsType evMatrix;
    typedef typename EVSolver::RealVectorType evVector;
    EVSolver ges;

    std::vector<OpPtr> Qop(d);
    std::vector<OpPtr> QTop(d);
    gsMatrix<T> ev;

    T avg_term = gamma;

    // Now, setup the Q's and update the D's
    for ( index_t i=0; i<d; ++i )
    {
        // Solve generalized eigenvalue problem
        ges.compute(local_stiff[i], local_mass[i], gsEigen::ComputeEigenvectors);
        // Q^T M Q = I, or M = Q^{-T} Q^{-1}
        // Q^T K Q = D, or K = Q^{-T} D Q^{-1}

        if (gamma != (T)(0))
        {
            gsMatrix<T> etrans = ges.eigenvectors().transpose()*local_mass[i]*gsMatrix<T>::Ones(local_mass[i].rows(),1);
            GISMO_ASSERT((etrans.block(1, 0, etrans.rows()-1, 1).array() < (T)(1)/100000000).all(),
                "gsPatchPreconditionerCreator::fastDiagonalizationOp: gamma!=0 only allowed for pure Neumann.");
            avg_term *= etrans(0,0) * etrans(0,0);
        }

        // From the eigenvalues, we setup the matrix D already in an Kroneckerized way.
        const evVector & D = ges.eigenvalues();

        const index_t loc = D.rows();
        glob /= loc;
        const index_t glob2 = sz / loc / glob;

        for ( index_t l=0; l<loc; ++l )
            for ( index_t m=0; m<glob; ++m )
                for ( index_t n=0; n<glob2; ++n )
                    diag( m + l*glob + n*loc*glob, 0 ) += D(l,0);

        // Finally, we store the eigenvectors
        ev.swap(const_cast<evMatrix&>(ges.eigenvectors()));

        // These are the operators representing the eigenvectors
        typename gsMatrixOp< gsMatrix<T> >::Ptr matrOp = makeMatrixOp( ev.moveToPtr() );
        Qop [i] = matrOp;
        // Here we are safe as long as we do not want to apply QTop after Qop got destroyed.
        QTop[i] = makeMatrixOp( matrOp->matrix().transpose() );
    }

    GISMO_ASSERT( glob == 1,
        "gsPatchPreconditionerCreator::fastDiagonalizationOp: Internal error." );

    diag(0,0) += avg_term;

    for ( index_t l=0; l<sz; ++l )
        diag( l, 0 ) = (T)(1)/diag( l, 0 );

    memory::unique_ptr< gsEigen::DiagonalMatrix<T,Dynamic> > diag_mat( new gsEigen::DiagonalMatrix<T,Dynamic>( give(diag) ) );

    return gsProductOp<T>::make(
        gsKroneckerOp<T>::make(QTop),
        makeMatrixOp(give(diag_mat)),
        gsKroneckerOp<T>::make(Qop)
        );
}

namespace {

// Get the tilde basis
template<typename T>
void tildeSpaceBasis_oneside(const gsTensorBSplineBasis<1,T>& basis, bool isLeftHandSide, gsMatrix<T>& tildeBasis, gsMatrix<T>& complBasis, bool bc = 0, const bool odd = true)
{
    // bc == false: Neumann (or any other not-eliminating bc), bc == true: Dirichlet

    const index_t b = (index_t)bc;
    const index_t p = basis.degree();
    const T h = basis.knots().maxIntervalLength();

    if (p-b<=0)
    {
        tildeBasis.resize(0,0); complBasis.resize(0,0);
        return;
    }

    const T u = (isLeftHandSide ? basis.knots().first() : basis.knots().last());
    gsMatrix<T> U(1,1);
    U(0,0) = u;

    std::vector< gsMatrix<T> > allDerivs;
    basis.evalAllDers_into(U, p-1, allDerivs);

    // Collect all derivatives in matrix (rows: derivatives, columns: basis functions)
    // normalize with h^(deriv)
    // use only odd derivatives if odd = true and only even derivatives if odd = false
    // skip last (on left) or first (on right) basis function since it's always in S-tilde
    gsMatrix<T> derivs;
    // If we have Dirichlet bc, we reduce the number of rows and columns by 1. We basically
    // eliminate the first row and the first column (first basis function, 0th derivative)
    // for the right-side, we have to remove the last basis function.
    derivs.setZero(p-b, p-b);

    const index_t offset = (isLeftHandSide ? b : 1);

    index_t i_start;
    if (odd) i_start = 1;
    else i_start = 2*b;

    for (index_t i = i_start; i < p; i += 2)
        for (index_t j = 0; j < p-b; ++j)
            derivs(i-b, j) = math::pow(h, i) * allDerivs[i](j+offset);

    typename gsMatrix<T>::JacobiSVD svd = derivs.jacobiSvd(gsEigen::ComputeFullV);

    index_t n_tilde;
    if (odd) n_tilde = (p + 1) / 2 - b;
    else n_tilde = p / 2 - b;

    tildeBasis = svd.matrixV().rightCols(n_tilde);
    complBasis = svd.matrixV().leftCols(p - b - n_tilde);
}

// Compute a basis for S-tilde and one for its orthogonal complement
template<typename T>
void tildeSpaceBasis(const gsBasis<T>& basis, gsSparseMatrix<T>& B_tilde, gsSparseMatrix<T>& B_compl, const gsBoundaryConditions<T>& bc, const bool odd = true)
{
    GISMO_ASSERT( nullptr != (dynamic_cast<const gsTensorBSplineBasis<1,T>*>(&basis)),
                    "gsPatchPreconditionersCreator<T>::getTildeSpaceBasisTransformation and "
                    "gsPatchPreconditionersCreator<T>::subspaceCorrectedMassSmootherOp only work with tensor-B-spline bases." );
    const gsTensorBSplineBasis<1,T>& bbasis = static_cast<const gsTensorBSplineBasis<1,T>&>(basis);

    patchSide west(0,boundary::west), east(0,boundary::east);
    bool bwest = ( bc.getConditionFromSide( west ) && bc.getConditionFromSide( west )->type() == condition_type::dirichlet );
    bool beast = ( bc.getConditionFromSide( east ) && bc.getConditionFromSide( east )->type() == condition_type::dirichlet );

    gsMatrix<T> b_L, b_compl_L;
    gsMatrix<T> b_R, b_compl_R;

    // Contruct space with vanishing odd derivatives
    tildeSpaceBasis_oneside(bbasis, true,  b_L, b_compl_L, bwest, odd);
    tildeSpaceBasis_oneside(bbasis, false, b_R, b_compl_R, beast, odd);

    const index_t n = bbasis.size() - (index_t)bwest - (index_t)beast;
    const index_t n_L = b_L.cols();
    const index_t m_L = b_L.rows();
    const index_t n_R = b_R.cols();
    const index_t m_R = b_R.rows();
    const index_t n_c_L = b_compl_L.cols();
    const index_t m_c_L = b_compl_L.rows();
    const index_t n_c_R = b_compl_R.cols();
    const index_t m_c_R = b_compl_R.rows();
    const index_t n_I = n - n_L - n_R - n_c_L - n_c_R;

    //GISMO_ENSURE ( n_I >= 0, "tildeSpaceBasis: Too few knots for that spline degree." );
    if ( n_I <= 0 )
    {
        static bool warned = false;
        if (!warned)
        {
            gsWarn << "tildeSpaceBasis was called with too few knots for that spline degree.\n"
                   << "So, we assume that S_tilde is empty.\n";
            warned = true;
        }

        B_tilde.resize(n,0);

        gsSparseEntries<T> E_compl;
        E_compl.reserve(n);
        for (index_t i = 0; i < n; ++i)
            E_compl.add(i,i,1.0);
        B_compl.resize(n,n);
        B_compl.setFrom(E_compl);
        return;
    }

    gsSparseEntries<T> E_tilde, E_compl;

    // put b_L into upper left block of S-tilde basis
    for (index_t j = 0; j < n_L; ++j)
    {
        for (index_t i = 0; i < m_L; ++i)
            E_tilde.add(i, j, b_L(i, j));
    }
    // fill identity matrix into interior part of S-tilde basis
    for (index_t j = 0; j < n_I; ++j)
    {
        E_tilde.add(m_L + j, n_L + j, 1.0);
    }
    // put b_R into lower right block of S-tilde basis
    for (index_t j = 0; j < n_R; ++j)
    {
        for (index_t i = 0; i < m_R; ++i)
            E_tilde.add(m_L + n_I + i, n_L + n_I + j, b_R(i, j));
    }
    B_tilde.resize(n, n_L + n_I + n_R);
    B_tilde.setFrom(E_tilde);

    // put b_compl_L into upper left block of complement basis
    for (index_t j = 0; j < n_c_L; ++j)
    {
        for (index_t i = 0; i < m_c_L; ++i)
            E_compl.add(i, j, b_compl_L(i, j));
    }
    // put b_compl_R into lower right block of complement basis
    for (index_t j = 0; j < n_c_R; ++j)
    {
        for (index_t i = 0; i < m_c_R; ++i)
            E_compl.add(m_c_L + n_I + i, n_c_L + j, b_compl_R(i, j));
    }
    B_compl.resize(n, n_c_L + n_c_R);
    B_compl.setFrom(E_compl);

}

template<typename T>
void constructTildeSpaceBasis(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    std::vector< gsSparseMatrix<T> >& B_tilde,
    std::vector< gsSparseMatrix<T> >& B_l2compl,
    const bool odd = true
    )
{
    const index_t d = basis.dim();
    B_tilde.resize(d);
    B_l2compl.resize(d);
    dirichlet::strategy ds = (dirichlet::strategy)opt.askInt("DirichletStrategy",dirichlet::elimination);
    if (ds == dirichlet::elimination)
    {
        for ( index_t i=0; i<d; ++i )
            tildeSpaceBasis(basis.component(d-1-i), B_tilde[i], B_l2compl[i], boundaryConditionsForDirection(bc,d-1-i), odd);
    }
    else
        GISMO_ERROR("Unknown Dirichlet strategy.");

}

// Constructs a matrix for swapping a tensor product
// from A (x) B (x) C (x) D (x) E  to   A (x) D (x) C (x) B (x) E,
// where only the dimensions of those matrices have to be given.
// Note that also those A, B, etc. could be tensor products (here as dimension just provide the products)
// Note that also those A, B, etc. could vanish. Then just provide a 1 as dimension, i.e., a scalar.
// So, literally every thinkable swap is possible.
template<typename T>
gsSparseMatrix<T> kroneckerSwap( index_t e, index_t d, index_t c, index_t b, index_t a )
{
    const index_t sz = a*b*c*d*e;
    gsSparseMatrix<T> result(sz,sz);
    gsSparseEntries<T> entries;
    entries.reserve(sz);
    for ( index_t i=0; i<a; ++i )
        for ( index_t j=0; j<b; ++j )
            for ( index_t k=0; k<c; ++k )
                for ( index_t l=0; l<d; ++l )
                    for ( index_t m=0; m<e; ++m )
                        entries.add( i+a*(j+b*(k+c*(l+d*m))), i+a*(l+d*(k+c*(j+b*m))), 1. );

    result.setFrom(entries);
    result.makeCompressed();

    return result;
}

} // anonymous namespace

template<typename T>
std::pair< std::vector< gsSparseMatrix<T> >, std::vector< gsSparseMatrix<T> > > gsPatchPreconditionersCreator<T>::getTildeSpaceBasisTransformation(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    std::vector< gsSparseMatrix<T> > B_tilde, B_l2compl;
    constructTildeSpaceBasis( basis, bc, opt, B_tilde, B_l2compl, !opt.askSwitch( "UseVanishingEvenDerivatives", false ) );
    return std::pair< std::vector< gsSparseMatrix<T> >, std::vector< gsSparseMatrix<T> > >( B_tilde, B_l2compl );
}

template<typename T>
typename gsPatchPreconditionersCreator<T>::OpUPtr gsPatchPreconditionersCreator<T>::subspaceCorrectedMassSmootherOp(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt,
    T sigma,
    T alpha,
    T beta
    )
{
    GISMO_ASSERT ( beta != 0, "gsPatchPreconditionersCreator<T>::subspaceCorrectedMassSmootherOp() does not work for beta==0." );

    // Get some properties
    const index_t d = basis.dim();
    const T h = basis.getMinCellLength();

    // Assemble univariate
    std::vector< gsSparseMatrix<T> > local_stiff = assembleTensorStiffness(basis, bc, opt);
    std::vector< gsSparseMatrix<T> > local_mass  = assembleTensorMass(basis, bc, opt);

    if (beta != 1)
    {
        for (index_t i=0; i!=d; ++i)
            local_stiff[i] *= beta;
    }

    // Setup of basis
    std::vector< gsSparseMatrix<T> > B_tilde(d), B_l2compl(d);
    constructTildeSpaceBasis(basis, bc, opt, B_tilde, B_l2compl);

    std::vector< gsSparseMatrix<T> > M_compl(d), K_compl(d), B_compl(d);
    std::vector< typename gsLinearOperator<T>::Ptr > M_tilde_inv(d);
    for ( index_t i=0; i<d; ++i )
    {
        // Transform the complement
        typename gsLinearOperator<T>::Ptr M_inv = makeSparseCholeskySolver(local_mass[i]);
        gsMatrix<T> B_compl_dense;
        M_inv->apply( B_l2compl[i], B_compl_dense );
        B_compl[i] = B_compl_dense.sparseView();

        // Setup of matrices and corresponding solvers
        gsSparseMatrix<T> M_tilde = B_tilde[i].transpose() * local_mass[i] * B_tilde[i];
        M_tilde_inv[i] = makeSparseCholeskySolver( M_tilde );
        M_compl[i] = B_compl[i].transpose() * local_mass[i] * B_compl[i];
        K_compl[i] = B_compl[i].transpose() * local_stiff[i] * B_compl[i];
    }

    // Setup of final operator
    typename gsSumOp<T>::uPtr result = gsSumOp<T>::make();

    for ( index_t type = 0; type < (1<<d); ++ type )
    {
        std::vector< typename gsLinearOperator<T>::Ptr > correction(0);
        gsSparseMatrix<T> transfer;

        std::vector< gsSparseMatrix<T>* > transfers(d);

        index_t numberInteriors = 0;

        // Setup of transfer
        for ( index_t j = 0; j<d; ++ j )
        {
            if ( type & ( 1 << j ) )
                transfers[j] = &(B_compl[d-1-j]);
            else
            {
                transfers[j] = &(B_tilde[d-1-j]);
                ++numberInteriors;
            }

            if ( j == 0 )
                transfer = *(transfers[j]);
            else
                transfer = transfers[j]->kron(transfer);

        }

        // If the subspace is not present, ignore it.
        if ( transfer.cols() == 0 )
            continue;

        // Setup of swap, where the boundary part is shifted to the begin

        index_t left = 1, current = transfers[d-1]->cols(), right = 1;
        for ( index_t j = 0; j < d-1; ++j )
            left *= transfers[j]->cols();

        for ( index_t j = d-1; j >= 0; --j )
        {
            if ( type & ( 1 << j ) )
            {
                transfer = transfer * kroneckerSwap<T>( right, current, left, 1, 1 );
                if ( j > 0 )
                {
                    left *= current;
                    current = transfers[j-1]->cols();
                    left /= current;
                }
            }
            else
            {
                if ( j > 0 )
                {
                    right *= current;
                    current = transfers[j-1]->cols();
                    left /= current;
                }
            }

        }

        // Setup of interior correction
        for ( index_t j = d-1; j>=0; --j )
        {
            if ( ! ( type & ( 1 << j ) ) )
                correction.push_back( M_tilde_inv[d-1-j] );
        }

        // If we are in the interior, we have to do the scaling here as there is no boundary correction.
        if ( numberInteriors == d )
            correction[0] = gsScaledOp<T>::make( correction[0], (T)(1)/( alpha + beta*(T)(numberInteriors)/(sigma*h*h) ) );

        // Setup of bondary correction
        if ( numberInteriors < d )
        {
            gsSparseMatrix<T> bc_matrix;

            {
                gsSparseMatrix<T> s(0,0);
                for ( index_t k = d-1; k>=0; --k )
                {
                    if ( type & ( 1 << k ) )
                    {
                        if ( s.rows() == 0 )
                            s = M_compl[d-1-k];
                        else
                            s = M_compl[d-1-k].kron(s);
                    }
                }
                bc_matrix = ( alpha + beta*(T)(numberInteriors)/(sigma*h*h) ) * s;
            }

            for ( index_t j = d-1; j>=0; --j )
            {
                if ( type & ( 1 << j ) )
                {
                    gsSparseMatrix<T> s(0,0);
                    for ( index_t k = d-1; k>=0; --k )
                    {
                        if ( type & ( 1 << k ) )
                        {
                            gsSparseMatrix<T>* chosenMatrix;
                            if ( j == k )
                                chosenMatrix = &(K_compl[d-1-k]);
                            else
                                chosenMatrix = &(M_compl[d-1-k]);

                            if ( s.rows() == 0 )
                                s = *chosenMatrix;
                            else
                                s = s.kron(*chosenMatrix);
                        }
                    }
                    bc_matrix += s;
                }
            }

            correction.push_back(makeSparseCholeskySolver(bc_matrix));
        }

        typename gsMatrixOp< gsSparseMatrix<T> >::Ptr transOp = makeMatrixOp(transfer.moveToPtr());
        // Setup of whole operator
        // The correction is the Kronecker-product of the operators in the vector correction.
        result->addOperator(
            gsProductOp<T>::make(
                makeMatrixOp( transOp->matrix().transpose() ),
                gsKroneckerOp<T>::make( correction ),
                transOp
                )
            );
    }

    return give(result);

}

} // namespace gismo
