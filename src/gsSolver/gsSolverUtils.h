/** @file gsSolverUtils.h

    @brief Utility class for PDE's solver related utils.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

namespace gismo
{


/** @brief
Utility class for PDE's solver related utils.

 \ingroup Solver
 */
template<class T> class gsSolverUtils
{
private:

    /// Default empty constructor
    gsSolverUtils()  { }

    ~gsSolverUtils() { } //destructor

public:

    /// \brief Finds the convergence rate for a list of errors and a list with elements sizes by the least square method.
    ///
    /// This function takes a list with errors and a list with elements sizes
    /// and find the convergence rate by fitting a straight line after taking
    /// the logarithm of the values. It uses the Least Square method to find
    /// the slope/gradient of the fitted line.
    /// \param[in] error_list list for the errors (in decreasing order)
    /// \param[in] h_list list of the elements sizes (in decreasing order)
    /// PS: This method of finding convergence rate require some initial refinements
    static T convergenceRateLS(std::vector<T> const & error_list, 
			       std::vector<T> const & h_list)
    {
        if (error_list.size() != h_list.size())
        {
            gsWarn << "The lists does not have same size!" << std::endl;
            return 0;
        }

        const index_t Nsize = error_list.size();

        // Define the matrix and vector for the least square problem
        gsMatrix<T> A(Nsize,2); 
        A.col(1).setOnes();
        A.col(0)      = gsAsConstMatrix<T>(h_list, Nsize, 1)    .array().log();
        gsVector<T> y = gsAsConstMatrix<T>(error_list, Nsize, 1).array().log();

        // Solve the least square problem
        gsVector<T> rate = A.fullPivHouseholderQr().solve(y);

        // Return the convergence rate
        return rate(0);
    }

    /// \brief Finds the spectral condition number of a small matrix.
    ///
    /// Find the condition number of a matrix by fist finding the eigenvalues
    /// of the matrix and then dividing the highest (absolute) eigenvalue by the
    /// lowest (absolute) eigenvalue. This method is computationally expensive
    /// and will only work on small matrices. The method assumes that the
    /// eigenvalues are reel.
    /// \param[in] matrix is a square matrix hows eigenvalues are found.
    /// \param[in] removeSingularity is true: ONE eigenvalue is removed to avoid inf condition number.
    static T conditionNumber(const gsMatrix<T> & matrix, bool removeSingularity = false)
    {
        GISMO_ENSURE(matrix.cols() == matrix.rows(), "Matrix must be square!");

        unsigned N_size = matrix.cols();

        if (N_size >= 10000)
            gsWarn << "Matrix has dimension " << N_size <<". It might take long to find eigenvalues...\n";

        //Compute eigenvalues
        //gsMatrix<T> eigen_values = A.eigenvalues().real();
        Eigen::EigenSolver<typename gsMatrix<T>::Base > eigen_values;
        eigen_values.compute(matrix, false);

        T tmp = std::abs(eigen_values.eigenvalues()(0,0).real());

        T eigen_low = tmp;
        T eigen_high= tmp;

        //Find lowest and highest eigenvalue
        for (unsigned k=0; k< N_size; ++k)
        {
            tmp = std::abs(eigen_values.eigenvalues()(k,0).real());

            // Remove the one zero eigenvalue from not using pressure avg.
            if(tmp < 1e-13 && removeSingularity)
            {
                removeSingularity = false;
                std::cout << "Removed the eigen value: " << tmp << std::endl;
                //In case the first eigenvalue has value zero, we use the second eigenvalue.
                if (k == 0)
                {
                    tmp = std::abs(eigen_values.eigenvalues()(k+1,0).real());
                    eigen_low = tmp;
                    eigen_high= tmp;
                }
                continue;
            }
            // Store its absolute value if lower the eigen_low or higher then eigen_high
            if( tmp < eigen_low ) {eigen_low  = tmp;}
            if( tmp > eigen_high) {eigen_high = tmp;}
        }
        //gsInfo << "Highest eigen value: " << eigen_high << " Lowest eigen value: " << eigen_low<< "\n";
        //Return condition number
        return eigen_high/eigen_low;
    }

    /// \brief Finds the spectral condition number of a small matrix.
    static T conditionNumber(const gsSparseMatrix<T> & matrix )
    {
        unsigned N_size = matrix.cols();

        if (N_size >= 10000)
            gsWarn << "Matrix has dimension " << N_size <<". It might take long to find eigenvalues...";

        // Copy to a dense matrix
        gsMatrix<T> matrixDence(matrix);

        return conditionNumber(matrixDence);
    }
};

/** @brief
class for calculating the eigenvalues of the lanczos matrix (see conjugate gradient for usage)

 \ingroup Solver
 */
template<class T>
class gsLanczosMatrix
{
    /**
     * @brief deriv derivative of characteristic polynomial
     * @param lambda evaluation point
     * @param k order of the polynomial (use n to start the recursion)
     * @return the derivative at position lambda
     */
    T deriv(T lambda , unsigned k)
    {
        if(k==0)
            return 0;
        else if(k==1)
            return -1;
        else
        {
           return (m_delta[k-1]-lambda)*deriv(lambda,k-1) - value(lambda,k-1)-m_gamma[k-2]*m_gamma[k-2]*deriv(lambda,k-2);
        }
    }

    /**
     * @brief value value of characteristic polynomial
     * @param lambda evaluation point
     * @param k order of the polynomial (use n to start the recursion)
     * @return the value at position lambda
     */
    T value(T lambda, unsigned k)
    {
        if(k==0)
            return 1;
        else if(k==1)
            return m_delta[0]-lambda;
        else
        {
           return (m_delta[k-1]-lambda)*value(lambda,k-1) - m_gamma[k-2]*m_gamma[k-2]*value(lambda,k-2);
        }
    }

    /**
     * @brief newtonIteration starts a newton iteration for searching the zeros of the charakteristic polynomial with initial value x0;
     * NOTE THIS FUNCTION NEEDS SOME OPTIMIZATION BECAUSE THE VALUES AND DERIVATIVES IN ONE ITERATION ARE CALCULATION TOO OFTEN DUE TO THE RECURSION.
     * @param x0 the initial value
     * @return the zero point (= Eigenvalue of the matrix)
     */
    T newtonIteration(T x0)
    {
        int iter =0;
        T res =1;
        T x_old = x0;
        T x_new;
        while(iter<m_maxIter && res > m_tol)
        {
            x_new = x_old - value(x_old,n)/deriv(x_old,n);
            res = std::abs(x_old - x_new);

            x_old = x_new;
            iter++;
        }
        return x_new;
    }

public:

    /**
     * @brief gsLanczosMatrix constructor for the lanczos matrix. The lanczos matrix is a symmetric tridiagonal matrix with diagonal delta and offdiagonal gamma.
     * @param _gamma  the diagonal
     * @param _delta the offdiagonal
     * @param _maxIter the number of maximal iterations of the Newton algorithm
     * @param _tol tolerace for the newton algorithm
     */
    gsLanczosMatrix(const std::vector<T> & _gamma, const std::vector<T> & _delta, int _maxIter =20, T _tol = 1.e-6): m_gamma(_gamma) , m_delta(_delta), m_maxIter(_maxIter), m_tol(_tol) { n = m_delta.size();}

    /**
     * @brief maxEigenvalue calculates the largest eigenvalue
     * @return the largest eigenvalue
     */
    T maxEigenvalue()
    {
        // x0 is rowsumNorm
        T x0 =std::abs(m_delta[0])+std::abs(m_gamma[0]);
        for(unsigned i=1;i<n-2;i++)
            if(std::abs(m_delta[i])+std::abs(m_gamma[i])+ std::abs(m_gamma[i-1])>x0)
                x0 = std::abs(m_delta[i])+std::abs(m_gamma[i])+ std::abs(m_gamma[i-1]);

        if(std::abs(m_delta[n-1])+std::abs(m_gamma[n-1])>x0)
            x0 = std::abs(m_delta[n-1])+std::abs(m_gamma[n-2]);
        return newtonIteration(x0);
    }

    /**
     * @brief minEigenvalue calculates the smallest eigenvalue
     * @return the smallest eigenvalue
     */
    T minEigenvalue()
    {
        T x0 = 0;
        return newtonIteration(x0);
    }

    /**
     * @brief matrixForm this function return the matrix form of the lanczos matrix as a gsSparseMatrix
     * @param L the lanczos matrix
     */
    void matrixForm(gsSparseMatrix<T> & L)
    {
        L.resize(n,n);
        std::vector<Eigen::Triplet<real_t> > list;
        list.reserve(3*n);

        list.push_back(Eigen::Triplet<real_t>(0,0,m_delta[0]));
        for(unsigned i = 1; i<n;i++)
        {
            list.push_back(Eigen::Triplet<real_t>(i,i-1,m_gamma[i-1]));
            list.push_back(Eigen::Triplet<real_t>(i-1,i,m_gamma[i-1]));
            list.push_back(Eigen::Triplet<real_t>(i,i,m_delta[i]));
        }
        L.setFromTriplets(list.begin(),list.end());

    }


private:
    const std::vector<T>& m_gamma;
    const std::vector<T>& m_delta;
    int m_maxIter;
    T m_tol;

    unsigned n;
};

}//namespace gismo
