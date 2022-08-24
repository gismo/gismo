#pragma once

#include <gsModeling/gsFitting.h>
#include <gsIO/gsWriteParaview.h>
#include <gsSolver/gsSolverUtils.h>

#include <iostream>
#include <fstream>

namespace gismo {

template <unsigned d, class T>
class gsFittingRWF : public gsFitting<T>
{
public:
    gsFittingRWF(gsMatrix<T> const & param_values,
                 gsMatrix<T> const & points,
                 gsBasis<T> & basis)
        : gsFitting<T>(param_values, points, basis)
    {}

    /// Iterate ...
    bool nextIteration(const T alpha, const T & toll, const T & tolerance, gsGeometry<T>& lambda, const bool ref_unif, const bool condcheck, const int dreg, const bool save);

protected:
    /// TODO: docu
    void compute( const gsGeometry<T>& lambda, const bool condcheck, const int dreg);

    /// TODO: docu
    void applySmoothing(const gsGeometry<T>& lambda, gsSparseMatrix<T> & A_mat, const int dreg);

public:
    /// Difference to 1D version: result has to have the basis with the correct knots.
    void findLambda(gsGeometry<T> &result, T l_min, T l_max);

    // Only works in 1D and for ordered data points. Compare Paper 1, 3.2.Bounded Slope regularization, first example
    T getL2ApproxErrorMidpoint(const std::vector<T> & errors) const;

    // Assumption: h is uniform & the same in both directions (i.e., 2D).
    // almost midpoint rule but the length is slightly off.
    T getL2ApproxErrorMidpointUniform(const std::vector<T> & errors) const;

    /// Compute the root mean square error for a vector with point--wise errors.
    /// RMSE = sqrt(1/N sum errors^2)
    T getRMSE(const std::vector<T> & errors) const;
};

template<unsigned d, class T>
bool gsFittingRWF<d,T>::nextIteration(const T alpha, const T & toll, const T & tolerance, gsGeometry<T> & lambda, const bool ref_unif, const bool condcheck, const int dreg, const bool save)
{
    if ( this->m_pointErrors.size() != 0 )
    {
        if ( this->m_max_error > tolerance )
        {
            if (ref_unif)  // Uniform refine the fitting basis
            {
                this->m_basis->uniformRefine(1);
            }
        }
        else
        {
            gsDebug << "Tolerance reached.\n";
            return false;
        }
    }

    // We run one fitting step and compute the errors
    this->compute(lambda,condcheck, dreg);
    this->computeErrors();

    return true;
}

template<unsigned d, class T>
void gsFittingRWF<d,T>::compute( const gsGeometry<T>& lambda, const bool condcheck, const int dreg)
{
    // Wipe out previous result
    if ( this->m_result )
        delete this->m_result;

    const int num_basis=this->m_basis->size();
    const int dimension=this->m_points.cols();

    // left side matrix
    gsSparseMatrix<T> A_mat(num_basis, num_basis);
    // To optimize sparse matrix an estimation of nonzero elements per
    // column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    // right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    m_B.setZero(); // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    this->assembleSystem(A_mat, m_B);
    applySmoothing(lambda, A_mat,dreg);
    A_mat.makeCompressed();

    typename gsSparseSolver<T>:: LU solver( A_mat );

    if(condcheck)
    {
        //solver.compute(A_mat);
        Eigen::SparseMatrix<T> I(A_mat.rows(),A_mat.cols());
        I.setIdentity();
        auto A_inv = solver.solve(I);
        real_t condnum = A_mat.norm() * A_inv.norm();           // This should be * not / right?
        gsInfo << "Condition number: " << condnum << std::endl;
    }

    // Solves for many right hand side columns
    gsMatrix<T> x;
    x = solver.solve(m_B);
    // finally generate the B-spline
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();
}

template<unsigned d, class T>
void gsFittingRWF<d,T>::applySmoothing(const gsGeometry<T>& lambda, gsSparseMatrix<T> & A_mat, const int dreg)
{
    const int dim = this->m_basis->dim();
    int stride;
    if (dreg==2)
        stride = dim*(dim+1)/2;               // number of 2.derivatives
    else if (dreg==1)                         // number of 1.derivatives
        stride = dim;
    else
        stride = dim*(dim+1)*(dim+2)/(3*2);   // number of 3.derivatives

    gsVector<int> numNodes(dim);
    for ( int i = 0; i!= dim; ++i )
        numNodes[i] = this->m_basis->degree(i) + math::floor(lambda.degree(i)/2);//+1;
    gsGaussRule<T> QuRule( numNodes ); // Reference Quadrature rule
    gsMatrix<T> quNodes, localA;
    gsVector<T> quWeights;
    gsMatrix<unsigned> actives;
    gsMatrix<T> reg;
    gsMatrix<T> derAll;

    typename gsBasis<T>::domainIter domIt = this->m_basis->makeDomainIterator();

    //T maximum = 0.0;
    for (; domIt->good(); domIt->next() )
    {
        // ----------- Higher Quadrature rule for a product?
        // Map the Quadrature rule to the element and compute basis derivatives
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );
        if (dreg==2)
        {
            //this->m_basis->deriv2_into(quNodes, derAll[3]);
            this->m_basis->deriv2_into(quNodes, derAll);
        }
        else if (dreg==1)
            this->m_basis->deriv_into(quNodes,derAll);
        else
            this->m_basis->deriv3_into(quNodes, derAll);
        this->m_basis->active_into(domIt->center, actives);
        const index_t numActive = actives.rows();
        localA.setZero(numActive, numActive );

        //Quadrature for (B-Spline) lambda
        lambda.eval_into(quNodes, reg);

        // perform the quadrature
        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            const T weight = quWeights[k] * reg(k);

            for (index_t i=0; i!=numActive; ++i)
                for (index_t j=0; j!=numActive; ++j)
                {
                    T localAij = 0; // temporary variable

                    for (int s = 0; s < stride; s++)
                    {
                        // Pure first/second/third derivatives
                        if (s < dim)
                        {
                            localAij += derAll(i * stride + s, k) *
                                        derAll(j * stride + s, k);
                        }
                        // Mixed derivatives
                        else
                        {
                            // 2* dudv N_i * dudv N_j
                            if (dreg==2)
                            {
                                localAij += 2 * derAll(i * stride + s, k) *
                                                derAll(j * stride + s, k) ;
                            }
                            // 3 * du^2dv N_i * du^2dv N_j + ...
                            else if (dreg==3)
                            {
                                localAij += 3 * derAll(i * stride + s, k) *
                                                derAll(j * stride + s, k) ;
                            }
                        }
                    }

                    localA(i, j) += weight * localAij;
                }
        }

        for (index_t i=0; i!=numActive; ++i)
        {
            const int ii = actives(i,0);
            for (index_t j=0; j!=numActive; ++j)
                A_mat( ii, actives(j,0) ) += localA(i,j);
        }
    }
}



template<unsigned d, class T>
void gsFittingRWF<d,T>::findLambda(gsGeometry<T> &result, T l_min, T l_max)
{
    // Future improvements: test if it works for general d.
    // Inititalize all coefficients to l_max
    for (unsigned i=0; i<result.coefsSize(); i++)
    {
        result.coef(i,0) = l_max;
    }

    // Find coefficients which have data in their support
    gsMatrix<unsigned> actives_new;
    gsMatrix<real_t> evals_new;
    for (index_t i=0; i< this->m_param_values.cols(); i++)
    {
        result.basis().active_into(this->m_param_values.col(i), actives_new);
        result.basis().eval_into(this->m_param_values.col(i), evals_new);
        for(index_t j=0; j<actives_new.rows(); j++)
        {
            if (evals_new(j,0) != 0)
                result.coef(actives_new(j,0), 0) = l_min;
        }
    }
}

template<unsigned d, class T>
T gsFittingRWF<d,T>::getL2ApproxErrorMidpoint(const std::vector<T> & errors) const
{
    GISMO_ASSERT(d==1,"Only implemented for 1D");
    T result = 0;

    gsVector<T> paramdist (errors.size()-1);
    for (size_t j = 0; j < errors.size() - 1 ;j++)
    {
        paramdist(j) = this->m_param_values(0,j+1)-this->m_param_values(0,j);
    }

    gsVector<> localh(errors.size());
    localh[0] = 0.5* paramdist[0];
    localh[errors.size()-1] = 0.5* paramdist[errors.size()-2];
    for (size_t j=1; j<errors.size()-1; j++)
    {
        localh[j] = 0.5*paramdist[j]+0.5*paramdist[j-1];
    }

    for (size_t i=0; i < errors.size(); i++)
    {
         result += errors[i]*errors[i]*localh[i];
    }

    return sqrt(result);
}

template<unsigned d, class T>
T gsFittingRWF<d,T>::getL2ApproxErrorMidpointUniform(const std::vector<T> & errors) const
{
    GISMO_ASSERT(d==2,"Implemented only for 2D");
    T result = 0;

    for (size_t i=0; i < errors.size(); i++)
    {
         result += errors[i]*errors[i];
    }

    T h = std::abs(std::max((this->m_param_values(1,0)-this->m_param_values(0,0)),
                                 (this->m_param_values(1,1)-this->m_param_values(0,1))));

    return sqrt(result * h * h);
}

template<unsigned d, class T>
T gsFittingRWF<d,T>::getRMSE(const std::vector<T> & errors) const
{
    // TODO: gismo assert d = 2.
    // type : L^2 error if h is uniform & the same in both directions
    T result = 0;

    for (size_t i=0; i < errors.size(); i++)
    {
         result += errors[i]*errors[i];
    }

    return sqrt(result / errors.size());
}


template <unsigned d, class T>
class gsFittingRWFErrorGuided : public gsFittingRWF<d,T>
{
public:
    gsFittingRWFErrorGuided(gsMatrix<T> const & param_values,
                            gsMatrix<T> const & points,
                            gsBasis<T> & basis)
        : gsFittingRWF<d,T>(param_values, points, basis)
    {}

protected:
    /// TODO: Explain parameters.
    index_t index(index_t i, index_t m1, index_t j = 0, index_t m2 = 1,  index_t k =0) const
    {
        return k*m1*m2 + j*m1 + i;
    }

    // Note: could be moved to gsTensorBasis.
    index_t reverseindex(index_t i, index_t m1, index_t j = 0, index_t m2 = 1,  index_t k =0, index_t m3 = 1) const
    {
        return index(m1-1-i,m1,m2-1-j,m2,m3-1-k);
    }

    void DecreaseLambdaErrorGuided(gsGeometry<T>& lambda,T toll,T alpha) const;

    void ConstraintCorrection1direction(gsTensorBSpline<d,T>& lambda, T alpha, bool forward) const;
    void ConstraintCorrection(gsTensorBSpline<d,T> &lambda, T alpha);

    void writeParaviewPointsLog( gsMatrix<T>& pts, gsMatrix<T>& vals, gsVector<unsigned> &np, const std::string& fn);

public:
    bool nextIteration(const T alpha, const T & toll, const T & tolerance, gsTensorBSpline<d,T>& lambda, const bool ref_unif, const bool dec_lambda, const bool condcheck, const int dreg, const bool save);

    void writeParaviewLog(const gsTensorBSpline<d,T> & field, std::string const & fn, unsigned npts=1000, bool mesh=false);
};

template<unsigned d, class T>
bool gsFittingRWFErrorGuided<d,T>::nextIteration(const T alpha, const T & toll, const T & tolerance, gsTensorBSpline<d,T> & lambda, const bool ref_unif, const bool dec_lambda, const bool condcheck, const int dreg, const bool save)
{
    if (dec_lambda)
    {
        // Does not change lambda if m_pointErrors.size() == 0
        DecreaseLambdaErrorGuided(lambda,toll,alpha);
        ConstraintCorrection(lambda,alpha);
    }

    return gsFittingRWF<d,T>::nextIteration(alpha, toll, tolerance, lambda, ref_unif, condcheck, dreg, save);
}

template<unsigned d, class T>
void gsFittingRWFErrorGuided<d, T>::DecreaseLambdaErrorGuided(gsGeometry<T>& lambda,T toll,T alpha) const
{
    // Find all points where the error is too large, decrease the corresponding control points of lambda with support there
    gsMatrix<T> paramsWhereErrTooBig;
    for ( size_t i=0; i < this->m_pointErrors.size(); i++ )
    {
        if (this->m_pointErrors[i] > toll)
        {
            index_t oldCols = paramsWhereErrTooBig.cols();
            paramsWhereErrTooBig.conservativeResize(2, oldCols + 1);
            paramsWhereErrTooBig.col(oldCols) = this->m_param_values.col(i);
        }
    }
    std::set<index_t> coefsToDecrease;
    gsMatrix<unsigned> actives;
    for (index_t i=0; i< paramsWhereErrTooBig.cols(); i++)
    {
        lambda.basis().active_into(paramsWhereErrTooBig.col(i), actives);
        for(index_t j=0; j<actives.rows(); j++)
        {
            coefsToDecrease.insert(actives(j,0));
        }
    }

    for ( std::set<index_t>::const_iterator it=coefsToDecrease.begin() ; it!=coefsToDecrease.end(); ++it )
    {
        if ( lambda.coef(*it, 0) > 1e-11 )
        {
            lambda.coef(*it, 0) *= alpha;
        }
    }
}

template<unsigned d, class T>
void gsFittingRWFErrorGuided<d, T>::ConstraintCorrection1direction(gsTensorBSpline<d,T>& lambda, T alpha, bool forward) const
{
    const gsTensorBSplineBasis<d, T> &basis = lambda.basis();
    index_t m1 =         basis.component(0).size();
    index_t m2 = d > 1 ? basis.component(1).size() : 1;
    index_t m3 = d > 2 ? basis.component(2).size() : 1;

    // Trick: Works for d=1 and d=2 as well, because next_id_v and next_id_w are equal to this_id.
    for (index_t i=0; i<m1; i++)
    {
        for (index_t j=0; j<m2; j++)
        {
            for (index_t k=0; k<m3; k++)
            {
                // In 1D i should go only up to m1-1.
                // Similarly for j and k in 2D and 3D, respectively.
                index_t next_i = std::min(i+1,m1-1);
                index_t next_j = std::min(j+1,m2-1);
                index_t next_k = std::min(k+1,m3-1);

                // Compute the neighbouring indices.
                index_t this_id   = forward ? index(i,      m1, j,      m2, k)      : reverseindex(i     , m1, j,      m2, k,     m3);
                index_t next_id_u = forward ? index(next_i, m1, j,      m2, k)      : reverseindex(next_i, m1, j,      m2, k,     m3);
                index_t next_id_v = forward ? index(i,      m1, next_j, m2, k)      : reverseindex(i,      m1, next_j, m2, k,     m3);
                index_t next_id_w = forward ? index(i,      m1, j,      m2, next_k) : reverseindex(i,      m1, j,      m2, next_k,m3);

                // Correct the neighbouring CPs.
                lambda.coef(next_id_u, 0) = std::min(lambda.coef(next_id_u, 0), lambda.coef(this_id,0)/alpha);
                lambda.coef(next_id_v, 0) = std::min(lambda.coef(next_id_v, 0), lambda.coef(this_id,0)/alpha);
                lambda.coef(next_id_w, 0) = std::min(lambda.coef(next_id_w, 0), lambda.coef(this_id,0)/alpha);
            }
        }
    }
}

template<unsigned d, class T>
void gsFittingRWFErrorGuided<d,T>::ConstraintCorrection(gsTensorBSpline<d,T> &lambda, T alpha)
{
    GISMO_ASSERT(d<4, "Implemented only for d=1, 2 or 3.");

    ConstraintCorrection1direction(lambda, alpha, true);
    ConstraintCorrection1direction(lambda, alpha, false);
}


template<unsigned d, class T>
void gsFittingRWFErrorGuided<d,T>::writeParaviewPointsLog( gsMatrix<T>& pts, gsMatrix<T>& vals, gsVector<unsigned> &np, const std::string& fn)
{
    gsParaviewCollection collection(fn);
    const int n = 2;
    const int dom = 2;

    // Das innere aus writeParaviewLog.

    if ( 3 - dom > 0 )
    {
        np.conservativeResize(3);
        np.bottomRows(3-dom).setOnes();
    }
    else if (dom > 3)
    {
        gsWarn<< "Cannot plot 4D data.\n"; return;
    }
    if ( 3 - n > 0 )
    {
        pts.conservativeResize(3,pts.cols() );
        pts.bottomRows(3-n).setZero();
    }
    else if (n > 3)
    {
        gsWarn<< "Data is more than 3 dimensions.\n";
    }
    if ( vals.rows() == 2)
    {
        vals.conservativeResize(3,pts.cols() );
        vals.bottomRows(1).setZero(); // 3-field.dim()
    }
    gsWriteParaviewTPgrid(pts, vals, np.template cast<index_t>(), fn);
    collection.addPart(fn, ".vts");

    collection.save();
}

template<unsigned d, class T>
void gsFittingRWFErrorGuided<d,T>::writeParaviewLog(const gsTensorBSpline<d,T> & field, std::string const & fn, unsigned npts, bool mesh)
{
    gsMatrix<T> ab = field.support();
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);

    gsVector<unsigned> np = uniformSampleCount(a, b, npts);

    gsMatrix<T> eval_geo = gsPointGrid(a, b, np);
    gsMatrix<T> eval_field = field.eval(eval_geo);

    writeParaviewPointsLog(eval_geo, eval_field, np, fn);
}

} // namespace
