#ifndef GSFITTINGRWF_H
#define GSFITTINGRWF_H

#endif // GSFITTINGRWF_H
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
    void compute( const gsTensorBSpline<d,T>& lambda, const bool condcheck, const int dreg);
    void compute( const gsBSpline<T>& lambda, const bool condcheck);
    //void applySmoothing(const gsTensorBSpline<d,T>& lambda, gsSparseMatrix<T> & A_mat);
    void applySmoothing(const gsGeometry<T>& lambda, gsSparseMatrix<T> & A_mat, const int dreg);
    void applySmoothing(const gsBSpline<T>& lambda, gsSparseMatrix<T> & A_mat);

    /// Refinement with THB, lambda remains unchanged.
    bool nextIteration(const T & tolerance, const T & err_threshold, const gsTensorBSpline<d,T>& lambda);

    /// Refinement with B-splines (in 1D), lambda unchanged.
    bool nextIteration(const T alpha, const T &toll, const T & tolerance, gsBSpline<T> &lambda, const bool l_ref, const bool reg, const T l_min, const T l_max, const bool condcheck, const int strategy, const bool remove);

    /// ...
    bool nextLambdaIteration(const T alpha, const T & toll, const T & tolerance, gsTensorBSpline<d,T>& lambda, const bool ref_unif, const bool dec_lambda, const bool condcheck, const int dreg);
    bool nextLambdaIteration(const T & tolerance, const T & err_threshold, gsBSpline<T>& lambda);
    void writeParaviewLog(const gsTensorBSpline<d,T> & field, std::string const & fn, unsigned npts=1000, bool mesh=false);
    void writeParaviewPointsLog( gsMatrix<T>& pts, gsMatrix<T>& vals, gsVector<unsigned> &np, const std::string& fn);
    void writeParaviewDeriv2Log(const gsGeometry<T> & field, std::string const & fn, unsigned npts);

    void findLambda1D(gsBSplineBasis<T> basis_lambda, gsBSpline<T>& result, T l_min, T l_max ,const index_t numLRef);

    /// Difference to 1D version: result has to have the basis with the correct knots.
    void findLambda2D(gsTensorBSpline<d,T> &result, T l_min, T l_max);

    T max(gsMatrix<T> & ders) const;

    T elemLength() const;

    T getL2error(const std::vector<T> & errors) const;
    T getL2errorTensor(const std::vector<T> & errors) const;
    T getRMSE(const std::vector<T> & errors) const;
};

template<unsigned d, class T>
T gsFittingRWF<d,T>::max(gsMatrix<T> & ders) const
{
    T maximum=0;
    for (index_t i=0; i< ders.rows(); i++)
        for (index_t j=0; j< ders.cols(); j++)
            if (std::abs(ders(i,j))>maximum)
                maximum=std::abs(ders(i,j));
    return maximum;
}

template<unsigned d, class T>
T gsFittingRWF<d,T>::elemLength() const
{
    gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d, T> *> (this->m_basis);
    return basis->getMaxCellLength();
}

template<unsigned d, class T>
void gsFittingRWF<d,T>::compute( const gsTensorBSpline<d,T>& lambda, const bool condcheck, const int dreg)
{
    // Wipe out previous result
    if ( this->m_result )
        delete this->m_result;

    const int num_basis=this->m_basis->size();
    const int dimension=this->m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis, num_basis);
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    m_B.setZero(); // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    this->assembleSystem(A_mat, m_B);

    applySmoothing(lambda, A_mat, dreg);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>:: LU solver( A_mat );
    gsMatrix<T> x;
    if (condcheck)
    {
        //T cnum;
        //cnum = gsSolverUtils<T>::conditionNumber(A_mat);
        //gsInfo << "Condition number: " << cnum << std::endl;

        //solver.compute(A_mat);
        Eigen::SparseMatrix<T> I(A_mat.rows(),A_mat.cols());
        I.setIdentity();
        auto A_inv = solver.solve(I);
        real_t condnum = A_mat.norm() * A_inv.norm();           // This should be * not / right?
        gsInfo << "Condition number: " << condnum << std::endl;

        // For bash, condition number check
        T h = this->m_basis->getMaxCellLength();

        std::string filename("conditioning_"+std::to_string(this->m_basis->degree(0))+"_"+std::to_string(lambda.degree(0))+ ".csv");
        std::ofstream file_out;

        file_out.open(filename, std::ios::app);
        file_out <<  " " << h << " " << condnum << "\n";
        file_out.close();

        //x = A_inv * m_B;
    }

    /*if ( solver.preconditioner().info() != Eigen::Success )
    {
        std::cerr<<  "The preconditioner failed. Aborting.";// << std::endl;
        this->m_result = NULL;
        return;
    }
    */
    // Solves for many right hand side  columns
    //else
    //    x = solver.solve(m_B); //toDense()
    x = solver.solve(m_B); //toDense()
    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();
}

template<unsigned d, class T>
void gsFittingRWF<d,T>::compute( const gsBSpline<T>& lambda, const bool condcheck)
{
    // Wipe out previous result
    if ( this->m_result )
        delete this->m_result;

    const int num_basis=this->m_basis->size();
    const int dimension=this->m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis, num_basis);
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    m_B.setZero(); // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    this->assembleSystem(A_mat, m_B);

    int dreg =2;
    applySmoothing(lambda, A_mat,dreg);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>:: LU solver( A_mat );

    if(condcheck)
    {
        gsHTensorBasis<d, T> & basis_temp = static_cast<gsHTensorBasis<d, T> &> (*this->m_basis);
        gsBSplineBasis<T>& basis = static_cast<gsBSplineBasis<T>&> (basis_temp.tensorLevel(0));
        //gsKnotVector<T> & knots = basis.knots();
        //gsInfo << "Knots in condcheck: " << knots.asMatrix() << std::endl;
        gsInfo << "h1: " << basis.getMinCellLength() << std::endl;
        gsInfo << "h2: " << this->m_basis->getMinCellLength() << std::endl;

        solver.compute(A_mat);
        Eigen::SparseMatrix<T> I(A_mat.rows(),A_mat.cols());
        I.setIdentity();
        auto A_inv = solver.solve(I);
        real_t condnum = A_mat.norm() * A_inv.norm();           // This should be * not / right?
        gsInfo << "Condition number: " << condnum << std::endl;

        // For bash, condition number check
        T h = this->m_basis->getMinCellLength();

        std::string filename("conditioning_"+std::to_string(this->m_basis->degree(0))+"_"+std::to_string(lambda.degree())+ ".csv");
        std::ofstream file_out;

        file_out.open(filename, std::ios::app);
        file_out <<  " " << h << " " << condnum << "\n";
        file_out.close();
    }

    //if ( solver.preconditioner().info() != Eigen::Success )
    //{
    //    gsWarn<<  "The preconditioner failed. Aborting.\n";
    //    this->m_result = NULL;
    //    return;
    //}
    // Solves for many right hand side  columns
    gsMatrix<T> x;
    x = solver.solve(m_B); //toDense()
    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve
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
bool gsFittingRWF<d,T>::nextIteration(const T alpha, const T & toll, const T & tolerance, gsBSpline<T> & lambda, const bool l_ref, const bool reg, const T l_min, const T l_max, const bool condcheck, const int strategy, const bool remove)
{
    if ( this->m_pointErrors.size() != 0 )
    {

        if ( this->m_max_error > tolerance || strategy == 6 )
        {
            if (strategy == 2)
            {
                // Uniform Refinement. TODO: Change to gsHTensorBasis.
                gsTHBSplineBasis<1,T>* basis = static_cast<gsTHBSplineBasis<1,T> *> (this->m_basis);
                basis->uniformRefine();
                if (l_ref)
                {
                    lambda.basis().uniformRefine();
                    findLambda1D(lambda.basis(), lambda, l_min, l_max,0);
                }
            }

            if (strategy == 6)
            {
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
                    if ( lambda.coef(*it, 0) > 1e-11 ) // Could also try other lambda min than 1e-11?
                    {
                        lambda.coef(*it, 0) *= alpha;
                    }
                }

                // ------ New Code for putting a threshold on the maximum incline/decline of lambda, version Bert
                //index_t optk = std::ceil(lambda.basis().getMaxCellLength() * tol_slope);   // I tried different parameters, from 4 to 18
                //T optk = lambda.basis().getMaxCellLength() * tol_slope;
                T optk = 1; // for the adaptive chosen alpha we fix jumps to 1.
                index_t m = lambda.basis().size();

                for (index_t i=0; i<m-1; i++)
                {
                    lambda.coef(i+1, 0) = math::min(lambda.coef(i+1, 0), lambda.coef(i,0)/math::pow(alpha,optk));
                }
                for (index_t i=0; i<m-1; i++)
                {
                    lambda.coef(m-i-2, 0) = math::min(lambda.coef(m-i-2, 0), lambda.coef(m-i-1,0)/math::pow(alpha,optk));
                }
            }
        }
        else
        {
            gsDebug << "Tolerance reached.\n";
            return false;
        }
    }

    // Knot removal for stabilization:
    if (remove)
    {
        gsHTensorBasis<d, T> & basis_temp = static_cast<gsHTensorBasis<d, T> &> (*this->m_basis);
        gsBSplineBasis<T>& basis = static_cast<gsBSplineBasis<T>&> (basis_temp.tensorLevel(0));
        gsKnotVector<T> knots = basis.knots();

        real_t h = basis.getMinCellLength();
        real_t h_left = -1.0/2;
        real_t h_right = 1.0/2;
        real_t h_ll = h_left  + h;
        real_t h_lr = h_left  - h;
        real_t h_rl = h_right + h;
        real_t h_rr = h_right - h;

        knots.remove(h_left);
        knots.remove(h_right);
        knots.remove(h_ll);
        knots.remove(h_lr);
        knots.remove(h_rl);
        knots.remove(h_rr);

        this ->m_basis = new gsTHBSplineBasis<d, T>(gsBSplineBasis<T>(knots));
    }

    // We run one fitting step and compute the errors
    if (reg)
        this->compute(lambda, condcheck);
    else
    {
        gsFitting<T>::compute(0.0);
    }
    this->computeErrors();

    // Knot insertion after stabilization, such that uniform refinement works as expected:
    if (remove)
    {
        gsHTensorBasis<d, T>& basis_temp = static_cast<gsHTensorBasis<d, T> &> (*this->m_basis);
        gsBSplineBasis<T>& basis = static_cast<gsBSplineBasis<T>&> (basis_temp.tensorLevel(0));
        gsKnotVector<T> knots = basis.knots();

        real_t h = basis.getMinCellLength();
        real_t h_left = -1.0/2;
        real_t h_right = 1.0/2;
        real_t h_ll = h_left+h;
        real_t h_lr = h_left-h;
        real_t h_rl = h_right+h;
        real_t h_rr = h_right-h;

        knots.insert(h_left);
        knots.insert(h_right);
        knots.insert(h_ll);
        knots.insert(h_lr);
        knots.insert(h_rl);
        knots.insert(h_rr);

        this ->m_basis = new gsTHBSplineBasis<d, T>(gsBSplineBasis<T>(knots));
    }

    return true;
}

template<unsigned d, class T>
bool gsFittingRWF<d,T>::nextIteration(const T & tolerance, const  T & err_threshold, const gsTensorBSpline<d,T> & lambda)
{
    // INVARIANT
    // look at iterativeRefine

    if ( this->m_pointErrors.size() != 0 )
    {

        if ( this->m_max_error > tolerance )
        {
            // if err_treshold is -1 we refine the m_ref percent of the whole domain
            T threshold = (err_threshold >= 0) ? err_threshold : this->setRefineThreshold(this->m_pointErrors);

            std::vector<unsigned> boxes = this->getBoxes(this->m_pointErrors, threshold);
            if(boxes.size()==0)
                return false;

            gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d,T> *> (this->m_basis);
            basis->refineElements(boxes);

            gsDebug << "inserted " << boxes.size() / (2 * d + 1) << " boxes.\n";
        }
        else
        {
            gsDebug << "Tolerance reached.\n";
            return false;
        }
    }

    // We run one fitting step and compute the errors
    this->compute(lambda);
    this->computeErrors();

    return true;
}

template<unsigned d, class T>
bool gsFittingRWF<d,T>::nextLambdaIteration(const T alpha, const T & toll, const T & tolerance, gsTensorBSpline<d,T> & lambda, const bool ref_unif, const bool dec_lambda, const bool condcheck, const int dreg)
{
    if ( this->m_pointErrors.size() != 0 )
    {
        if ( this->m_max_error > tolerance )
        {
            if (ref_unif)
            {
                gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d,T> *> (this->m_basis);
                basis->uniformRefine(1);

            }
            if (dec_lambda)
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

                // Variation Bert in 2D:
                index_t m1 = lambda.basis().size(0);
                index_t m2 = lambda.basis().size(1);

                // First we go through all rows, keep v constant

                gsInfo << "Starting constraint correction in u direction" << std::endl;
                for (index_t j=0; j<m2; j++)
                {
                    for (index_t i=0; i<m1-1; i++)
                    {
                        // Adaptation for row j left to right:
                        index_t this_id = lambda.basis().index(i,j);
                        index_t next_id = lambda.basis().index(i+1,j);
                        lambda.coef(next_id, 0) = std::min(lambda.coef(next_id, 0), lambda.coef(this_id,0)/alpha);
                    }
                    for (index_t i=0; i<m1-1; i++)
                    {
                        // Adaptation for row j right to left:
                        index_t this_id = lambda.basis().index(m1-i-1,j);
                        index_t next_id = lambda.basis().index(m1-i-2,j);
                        lambda.coef(next_id, 0) = std::min(lambda.coef(next_id, 0), lambda.coef(this_id,0)/alpha);
                    }
                }

                // Second we go through all columns, keep u constant

                gsInfo << "Starting constraint correction in v direction" << std::endl;
                for (index_t i=0; i<m1; i++)
                {
                    for (index_t j=0; j<m2-1; j++)
                    {
                        // Adaptation for column i bottom to top:
                        index_t this_id = lambda.basis().index(i,j);
                        index_t next_id = lambda.basis().index(i,j+1);
                        lambda.coef(next_id, 0) = std::min(lambda.coef(next_id, 0), lambda.coef(this_id,0)/alpha);
                    }
                    for (index_t j=0; j<m2-1; j++)
                    {
                        // Adaptation for column i top to bottom:
                        index_t this_id = lambda.basis().index(i,m2-j-1);
                        index_t next_id = lambda.basis().index(i,m2-j-2);
                        lambda.coef(next_id, 0) = std::min(lambda.coef(next_id, 0), lambda.coef(this_id,0)/alpha);
                    }
                }

                gsInfo << "Done constraint correction" << std::endl;

                gsTensorBSpline<d,T> lambda_log (lambda);
                for (index_t i = 0; i<lambda_log.basis().size(); i++)
                {
                    lambda_log.coef(i)(0)=std::log10(lambda_log.coef(i)(0));
                }
                // gsInfo << "Writing lambda.xml" << std::endl;
                gsWriteParaview(lambda_log,"lambda_log",(m1+1)*(m2+1));
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
bool gsFittingRWF<d,T>::nextLambdaIteration(const T & tolerance, const  T & err_threshold, gsBSpline<T> & lambda)
{
    if ( this->m_pointErrors.size() != 0 )
    {
        gsInfo <<"Regularisation function lambda is uniformly refined." << std::endl;
        lambda.uniformRefine(1);
        if ( this->m_max_error > tolerance )
        {
            gsMatrix<T> paramsWhereErrTooBig;
            gsInfo << "There are " << this->m_pointErrors.size() << " errors, "
                   << this->m_param_values.cols() << " parameters and "
                   << paramsWhereErrTooBig.cols() << " too big ones so far." << std::endl;
            for ( size_t i=0; i < this->m_pointErrors.size(); i++ )
            {
                if (this->m_pointErrors[i] > tolerance)   // Should be err_threshold?
                {
                    index_t oldCols = paramsWhereErrTooBig.cols();
                    paramsWhereErrTooBig.conservativeResize(2, oldCols + 1);
                    paramsWhereErrTooBig.col(oldCols) = this->m_param_values.col(i);
                }
            }
            gsInfo << "There were " << paramsWhereErrTooBig.cols() << " points with the error above tolerance." << std::endl;

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

            gsInfo << coefsToDecrease.size() << " coefs should be made smaller." << std::endl;

            for ( std::set<index_t>::const_iterator it=coefsToDecrease.begin() ; it!=coefsToDecrease.end(); ++it )
            {
                lambda.coef(*it, 0) *= 0.1;
            }
            gsInfo << "Writing lambda.xml" << std::endl;
            gsWriteParaview(lambda,"lambda_"+std::to_string(lambda.coefs().rows()), 2000, false, true);
        }
        else
        {
            gsDebug << "Tolerance reached.\n";
            return false;
        }
    }

    // We run one fitting step and compute the errors
    this->compute(lambda);
    this->computeErrors();

    return true;
}


template<unsigned d, class T>
void gsFittingRWF<d,T>::writeParaviewPointsLog( gsMatrix<T>& pts, gsMatrix<T>& vals, gsVector<unsigned> &np, const std::string& fn)
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
void gsFittingRWF<d,T>::writeParaviewLog(const gsTensorBSpline<d,T> & field, std::string const & fn, unsigned npts, bool mesh)
{
    gsMatrix<T> ab = field.support();
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);

    gsVector<unsigned> np = uniformSampleCount(a, b, npts);

    gsMatrix<T> eval_geo = gsPointGrid(a, b, np);
    gsMatrix<T> eval_field = field.eval(eval_geo);

    writeParaviewPointsLog(eval_geo, eval_field, np, fn);
}

template<unsigned d, class T>
void gsFittingRWF<d,T>::writeParaviewDeriv2Log(const gsGeometry<T> &field, std::string const & fn, unsigned npts)
{
    gsMatrix<T> ab = field.support();
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);
    gsVector<unsigned> np = uniformSampleCount(a, b, npts);
    gsMatrix<T> eval_geo = gsPointGrid(a, b, np);
    gsMatrix<T> deriv2_field = field.deriv2(eval_geo);
    gsMatrix<T> coefs (deriv2_field.cols(),2);

    for (index_t j=0; j< deriv2_field.cols() ; j++)
    {
        coefs(j,0) = eval_geo(0,j);
        coefs(j,1) = deriv2_field(0,j)/math::abs(deriv2_field(0,j)) * math::log( 1 + math::abs(deriv2_field(0,j)));
    }

    std::vector<T> knots (eval_geo.cols()+2);
    knots[0] = eval_geo(0,0);
    for (index_t i=0; i<eval_geo.cols(); i++)
    {
        knots[i+1] = eval_geo(0,i);
    }
    knots[eval_geo.cols()+1] = eval_geo(0,eval_geo.cols()-1);

    gsKnotVector<T> knotvector (knots,1);
    gsBSplineBasis<> basis( knotvector );
    gsBSpline<T> result(basis,coefs);

    gsWriteParaview(result, fn, 10000, false, true);
}


template<unsigned d, class T>
void gsFittingRWF<d,T>::findLambda1D(gsBSplineBasis<T> basis_lambda, gsBSpline<T> &result, T l_min, T l_max, const index_t numLRef)
{
    // Future improvements: test if it works for general d.

    // set coefficients
    int n = basis_lambda.size();
    gsMatrix<> coefs_lambda(n,1);

    // Inititalize coefficients
    for (index_t i=0; i<n; i++)
    {
        coefs_lambda(i,0) = l_max;
    }

    // Find coefficients which have data in their support
    gsMatrix<unsigned> actives_new;
    gsMatrix<T> evals_new;
    for (index_t i=0; i< this->m_param_values.cols(); i++)
    {
        basis_lambda.active_into(this->m_param_values.col(i), actives_new);
        basis_lambda.eval_into(this->m_param_values.col(i), evals_new);
        for(index_t j=0; j<actives_new.rows(); j++)
        {
            if (evals_new(j,0) != 0)
                coefs_lambda(actives_new(j,0),0) = l_min;
        }
    }
    result = gsBSpline<T>(basis_lambda, coefs_lambda);
}

template<unsigned d, class T>
void gsFittingRWF<d,T>::findLambda2D(gsTensorBSpline<d,T> &result, T l_min, T l_max)
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
T gsFittingRWF<d,T>::getL2error(const std::vector<T> & errors) const
{
    // TODO: gismo assert d = 1.
    // type : L^2 error
    T result = 0;

    // Get local h: Only works if param are ordered
    gsVector<T> paramdist (errors.size()-1);

    for (size_t j = 0; j < errors.size() - 1 ;j++)
    {
        paramdist(j) = this->m_param_values(0,j+1)-this->m_param_values(0,j);
    }
    gsVector<> localh(errors.size());

    // ToDo: add 0.5 (dist paramdist, rand), Annahme: Rand = 1. und letzter param
    localh[0] = 0.5* paramdist[0];
    localh[errors.size()-1] = 0.5* paramdist[errors.size()-2];
    for (size_t j=1; j<errors.size()-1; j++)
    {
        localh[j] = 0.5*paramdist[j]+0.5*paramdist[j-1];
    }

    for (size_t i=0; i < errors.size(); i++)
    {
         result += errors[i]*errors[i]*localh[1];           // For error guided, Berts curve: localh[i], otherwise: localh[1] for uniform sampled data
    }

    return sqrt(result);
}

template<unsigned d, class T>
T gsFittingRWF<d,T>::getL2errorTensor(const std::vector<T> & errors) const
{
    // TODO: gismo assert d = 2.
    // type : L^2 error if h is uniform & the same in both directions
    T result = 0;

    for (size_t i=0; i < errors.size(); i++)
    {
         result += errors[i]*errors[i];
    }

    T h = std::abs(std::max((this->m_param_values(1,0)-this->m_param_values(0,0)),
                                 (this->m_param_values(1,1)-this->m_param_values(0,1))));

    gsInfo << "h: " << h << std::endl;

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

}
