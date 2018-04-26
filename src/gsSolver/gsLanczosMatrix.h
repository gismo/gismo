/** @file gsLanczosMatrix.h

    @brief Class for representing a Lanczos matrix and calculating its eigenvalues

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

namespace gismo
{

/** 
 * @brief Class for representing a Lanczos matrix and calculating its eigenvalues
 * 
 * The Lanczos matrix is a symmetric tridiagonal matrix with diagonal delta and offdiagonal gamma.
 * 
 * @ingroup Solver
 */
template<class T>
class gsLanczosMatrix
{
public:

    /**
     * @brief Constructor for the Lanczos matrix
     * The Lanczos matrix is a symmetric tridiagonal matrix with diagonal delta and offdiagonal gamma.
     * 
     * @param _gamma  the diagonal
     * @param _delta the offdiagonal
     * @param _maxIter the number of maximal iterations of the Newton algorithm
     * @param _tol tolerace for the newton algorithm
     */
    gsLanczosMatrix(const std::vector<T> & _gamma, const std::vector<T> & _delta, int _maxIter =20, T _tol = 1.e-6)
     : m_gamma(_gamma) , m_delta(_delta), m_maxIter(_maxIter), m_tol(_tol), m_n(m_delta.size()) {}

    /**
     * @brief Calculates the largest eigenvalue
     * @return the largest eigenvalue
     */
    T maxEigenvalue()
    {
         //This function might be very slow, use instead the matrix form
        if(m_n==1)
            return m_delta[0];

        // x0 is rowsumNorm
        T x0 =math::abs(m_delta[0])+math::abs(m_gamma[0]);


        for(unsigned i=1;i<m_n-2;i++)
            if(math::abs(m_delta[i])+math::abs(m_gamma[i])+ math::abs(m_gamma[i-1])>x0)
                x0 = math::abs(m_delta[i])+math::abs(m_gamma[i])+ math::abs(m_gamma[i-1]);

        if(math::abs(m_delta[m_n-1])+math::abs(m_gamma[m_n-2])>x0)
            x0 = math::abs(m_delta[m_n-1])+math::abs(m_gamma[m_n-2]);

        return newtonIteration(x0);
    }

    /**
     * @brief Calculates the smallest eigenvalue
     * @return the smallest eigenvalue
     */
    T minEigenvalue()
    {
        //This function might be very slow, use instead the matrix form
        if(m_n==1)
            return m_delta[0];
        T x0 = 0;
        return newtonIteration(x0);
    }

    /**
     * @brief This function returns the Lanczos matrix as a gsSparseMatrix
     * 
     * @param[out] L the location the matrix is to be written to
     */
    void matrixForm(gsSparseMatrix<T> & L)
    {
        L.resize(m_n,m_n);
        std::vector<Eigen::Triplet<T> > list;
        list.reserve(3*m_n);

        list.push_back(Eigen::Triplet<T>(0,0,m_delta[0]));
        for(unsigned i = 1; i<m_n;i++)
        {
            list.push_back(Eigen::Triplet<T>(i,i-1,m_gamma[i-1]));
            list.push_back(Eigen::Triplet<T>(i-1,i,m_gamma[i-1]));
            list.push_back(Eigen::Triplet<T>(i,i,m_delta[i]));
        }
        L.setFromTriplets(list.begin(),list.end());

    }

private:

    /**
     * @brief Derivative of characteristic polynomial
     * 
     * @param lambda evaluation point
     * @param k order of the polynomial (use n to start the recursion)
     * @return the derivative at position lambda
     */
    T deriv(T lambda , unsigned k)
    {
        if(k==0)
            return T(0);
        else if(k==1)
            return T(-1);
        else
        {
           return (m_delta[k-1]-lambda)*deriv(lambda,k-1) - value(lambda,k-1)-m_gamma[k-2]*m_gamma[k-2]*deriv(lambda,k-2);
        }
    }

    /**
     * @brief Value of characteristic polynomial
     * 
     * @param lambda evaluation point
     * @param k order of the polynomial (use n to start the recursion)
     * @return the value at position lambda
     */
    T value(T lambda, unsigned k)
    {
        if(k==0)
            return T(1);
        else if(k==1)
            return m_delta[0]-lambda;
        else
        {
           return (m_delta[k-1]-lambda)*value(lambda,k-1) - m_gamma[k-2]*m_gamma[k-2]*value(lambda,k-2);
        }
    }

    /**
     * @brief Newton iteration for searching the zeros of the characteristic polynomial
     * 
     * @param x0 the initial value
     * @return the zero point (= Eigenvalue of the matrix)
     */
    T newtonIteration(T x0)
    {
        int iter =0;
        T res =1;
        T x_old = x0;
        T x_new = x0;
        while(iter<m_maxIter && res > m_tol)
        {
            x_new = x_old - value(x_old,n)/deriv(x_old,n);
            res = math::abs(x_old - x_new);

            x_old = x_new;
            iter++;
        }
        return x_new;
    }

private:
    const std::vector<T>& m_gamma;
    const std::vector<T>& m_delta;
    int m_maxIter;
    T m_tol;
    size_t m_n;
};

} // namespace gismo
