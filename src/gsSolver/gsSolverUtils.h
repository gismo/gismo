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


/// \brief Utility class for PDE's solver related utils.
///
/// \ingroup Solver
template<class T>
class gsSolverUtils
{
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
            gsWarn << "The lists does not have same size!\n";
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

        // Compute eigenvalues
        //gsMatrix<T> eigen_values = A.eigenvalues().real();
        typename gsMatrix<T>::EigenSolver eigen_values;
        eigen_values.compute(matrix, false);

        T tmp = math::abs(eigen_values.eigenvalues()(0,0).real());

        T eigen_low = tmp;
        T eigen_high= tmp;

        // Find lowest and highest eigenvalue
        for (unsigned k=0; k< N_size; ++k)
        {
            tmp = math::abs(eigen_values.eigenvalues()(k,0).real());

            // Remove the one zero eigenvalue from not using pressure avg.
            if (tmp < 1e-13 && removeSingularity)
            {
                removeSingularity = false;
                gsDebug << "Removed the eigen value: " << tmp << "\n";
                // In case the first eigenvalue has value zero, we use the second eigenvalue.
                if (k == 0)
                {
                    tmp = math::abs(eigen_values.eigenvalues()(k+1,0).real());
                    eigen_low = tmp;
                    eigen_high= tmp;
                }
                continue;
            }
            // Store its absolute value if lower the eigen_low or higher then eigen_high
            if (tmp < eigen_low ) {eigen_low  = tmp;}
            if (tmp > eigen_high) {eigen_high = tmp;}
        }
        //gsDebug << "Highest eigen value: " << eigen_high << " Lowest eigen value: " << eigen_low<< "\n";
        // Return condition number
        return eigen_high/eigen_low;
    }

    /// \brief Finds the spectral condition number of a small matrix.
    static T conditionNumber(const gsSparseMatrix<T> & matrix )
    {
        unsigned N_size = matrix.cols();

        if (N_size >= 10000)
            gsWarn << "Matrix has dimension " << N_size <<". It might take long to find eigenvalues...\n" << std::flush;

        // Copy to a dense matrix
        gsMatrix<T> matrixDense(matrix);

        return conditionNumber(matrixDense);
    }

private:
    gsSolverUtils() {} // No objects of this class

};

} // namespace gismo
