/** @file gsLanczosMatrix.h

    @brief Class for representing a Lanczos matrix and calculating its eigenvalues

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, S. Takacs
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
     * @param gamma    The off diagonal (the object stores a reference to this vector)
     * @param delta    The diagonal (the object stores a reference to this vector)
     */
    gsLanczosMatrix(const std::vector<T> & gamma, const std::vector<T> & delta)
     : m_gamma(gamma), m_delta(delta), m_n(m_delta.size())
    { GISMO_ASSERT( m_delta.size() == m_gamma.size() + 1, "Size mismatch." ); }

    /**
     * @brief Calculates the largest eigenvalue
     *
     * @param maxIter  The number of maximal iterations of the Newton algorithm
     * @param tol      Tolerace for the Newton algorithm
     *
     * @return The largest eigenvalue
     */
    T maxEigenvalue(index_t maxIter = 20, T tol = 1.e-6)
    {
        if (m_n==1)
            return m_delta[0];

        // Use Gerschgorin theorem to compute upper bound for eigenvalues
        T x0 = 0;
        for (size_t i=0;i<m_n;++i)
        {
            const T tmp = math::abs(m_delta[i])
                        + ( i<m_n-1 ? math::abs(m_gamma[i])   : T(0) )
                        + ( i>0     ? math::abs(m_gamma[i-1]) : T(0) );
            if (tmp>x0) x0 = tmp;
        }

        return newtonIteration(x0, maxIter, tol);
    }

    /**
     * @brief Calculates the smallest eigenvalue
     *
     * @param maxIter  The number of maximal iterations of the Newton algorithm
     * @param tol      Tolerace for the Newton algorithm
     *
     * @return The smallest eigenvalue
     */
    T minEigenvalue(index_t maxIter = 20, T tol = 1.e-6)
    {
        if (m_n==1)
            return m_delta[0];

        T x0 = 0;
        return newtonIteration(x0, maxIter, tol);
    }

    /**
     * @brief This function returns the Lanczos matrix as \a gsSparseMatrix
     */
    gsSparseMatrix<T> matrix()
    {
        gsSparseMatrix<T> L(m_n,m_n);
        gsSparseEntries<T> list;
        list.reserve(3*m_n);

        list.add(0,0,m_delta[0]);
        for (size_t i = 1; i<m_n;i++)
        {
            list.add(i,i-1,m_gamma[i-1]);
            list.add(i-1,i,m_gamma[i-1]);
            list.add(i,i,m_delta[i]);
        }
        L.setFrom(list);
        return L;
    }

private:

    /**
     * @brief Evalutates characteristic polynomial
     *
     * @param lambda evaluation point
     * @return the value and the derivative at position lambda
     */
    std::pair<T,T> eval( T lambda )
    {
        std::vector<T> value(m_n+1);
        std::vector<T> deriv(m_n+1);

        value[0] = T(1);
        value[1] = m_delta[0]-lambda;
        deriv[0] = T(0);
        deriv[1] = T(-1);
        for (size_t k=2; k<m_n+1; ++k)
        {
            value[k] = (m_delta[k-1]-lambda) * value[k-1] - m_gamma[k-2]*m_gamma[k-2]*value[k-2];
            deriv[k] = (m_delta[k-1]-lambda) * deriv[k-1] - value[k-1] - m_gamma[k-2]*m_gamma[k-2]*deriv[k-2];
        }
        return std::pair<T,T>(value[m_n],deriv[m_n]);
    }

    /**
     * @brief Newton iteration for searching the zeros of the characteristic polynomial
     *
     * @param x0 the initial value
     * @param maxIter  The number of maximal iterations of the Newton algorithm
     * @param tol      Tolerace for the Newton algorithm
     *
     * @return the root (= eigenvalue of the matrix)
     */
    T newtonIteration(T x0, index_t maxIter, T tol)
    {
        index_t iter = 0;
        T res = 1;
        T x_old = x0;
        T x_new = x0;
        while (iter < maxIter && res > tol)
        {
            const std::pair<T,T> ev = eval(x_old);
            const T& value = ev.first;
            const T& deriv = ev.second;

            x_new = x_old - value/deriv;
            res = math::abs(x_old - x_new);

            x_old = x_new;
            iter++;
        }
        return x_new;
    }

private:
    const std::vector<T>& m_gamma;
    const std::vector<T>& m_delta;
    size_t m_n;
};

} // namespace gismo
