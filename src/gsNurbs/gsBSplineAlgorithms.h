/** @file gsBSplineAlgorithms.h

    @brief Implementation of common algorithms for B-splines

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/** @namespace gismo::bspline

    @brief
    This namespace contains implementation of B-spline related algorithms.
    
    \ingroup Nurbs
*/
namespace bspline
{
    /// Input: an iterator pointing to the the biggest knot less than \a u,
    /// evaluation point \a u, output the value of all basis functions
    /// which are active at \a u.
    /// Uses B-spline recursion for evaluation
    template <class T, typename KnotIterator, typename Derived>
    void evalBasis( T u,
                    KnotIterator knot,
                    int deg,
                    gsEigen::MatrixBase<Derived> const & result )
    {       
        STACK_ARRAY(T, left, deg  + 1);
        STACK_ARRAY(T, right, deg + 1);
        evalBasis( u, knot, deg, left, right, result);
    }

    template <class T, typename KnotIterator, typename Derived>
    void evalBasis( T u,
                    KnotIterator knot,
                    int deg,
                    T left[], T right[], 
                    gsEigen::MatrixBase<Derived> const & result )
    {
        gsEigen::MatrixBase<Derived>& res = const_cast<gsEigen::MatrixBase<Derived>&>(result);
        
        res(0,0)= (T)(1.0);  // 0-th degree function value
        
        for(int j=1; j<= deg; j++) // For all degrees ( ndu column)
        {
            left[j]  = u - *(knot+1-j);
            right[j] = *(knot+j) - u;
            T saved = (T)(0) ;
            for(int r=0; r<j ; r++) // For all (except the last)  basis functions of degree j ( ndu row)
            {
                const T temp = res(r,0) / ( right[r+1] + left[j-r] );
                res(r,0)     = saved + right[r+1] * temp ;// r-th function value of degree j
                saved = left[j-r] * temp ;
            }
            res(j,0)     = saved;
        }
    }

    /// Evaluation for degree 1 B-spline basis
    template <class T, typename KnotIterator, typename Derived>
    void evalDeg1Basis(  const T & u,
                         const KnotIterator knot,
                         gsEigen::MatrixBase<Derived> const & result )
    {
        result(0,0) =  u       - (*knot+1);
        result(1,0) =  (*knot) - u        ;
        const T dn  = (*knot)  - (*knot+1) ;
        if( dn !=  0.0 )
            result.array /= dn;
    }
    
    /// Evaluation for degree 2 B-spline basis
    template <class T, typename KnotIterator, typename Derived>
    void evalDeg2Basis(  const T & u,
                         const KnotIterator knot,
                         gsEigen::MatrixBase<Derived> const & result )
    {
        
    }
    
    /// Evaluation for degree 3 B-spline basis
    template <class T, typename KnotIterator, typename MatrixType>
    void evalDeg3Basis(  const T & u,
                         const KnotIterator knot,
                         MatrixType & result )
    {
        
    }
    

    /// Input: parameter position \a u, KnotIterator \a knot identifying the active interval,
    /// degree \a deg, Output: table \a N.
    /// Computes deBoor table used in Algorithm A2.5 in NURBS book.
    /// Used by evalDerSingle_into(...), evalAllDerSingle_into(...) and evalSingle_into(...).
    template <class T, typename KnotIterator>
    void deBoorTriangle( T u,
                         KnotIterator knot,
                         int deg,
                         T N[] )
    {
        const int table_size = deg+1;
        T saved;

        for( int j = 0; j <= deg; j++ ) // Initialize zeroth-degree functions
            if( u >= *(knot+j) && u < *(knot+j+1) )
                N[ j*table_size ] = 1;
            else
                N[ j*table_size ] = 0;
        for( int k = 1; k <= deg; k++ ) // Compute full triangular table
        {
            if( N[(k-1)] == 0 )
                saved = 0;
            else
                saved = ( (u-*(knot) )*N[ k-1 ])/( *(knot+k)- *(knot));
            for( int j = 0; j < deg-k+1; j++)
            {
                if( N[ (j+1)*table_size + (k-1) ] == 0 )
                {
                    N[ j*table_size + k ] = saved;
                    saved = 0;
                }
                else
                {
                    const T Uleft = *(knot+j+1);
                    const T Uright = *(knot+j+k+1);
                    const T temp = N[ (j+1)*table_size + (k-1) ]/(Uright-Uleft);
                    N[ j*table_size + k ] = saved + (Uright-u)*temp;
                    saved = (u-Uleft)*temp;
                }
            }
        }
    }


    /// Input: a sequence of 2*p+1 knots starting at iterator position
    /// \a knot, evaluation point \a u, output the value of the
    /// B-spline defined by \a coefs
    /// Uses DeBoor's algorithm for evaluation
    /// \todo add an overload with stride parameter for \a coefs
    template <class T, typename KnotIterator, typename MatrixType>
    void evalGeo(  const T & u,
                   const KnotIterator & knot,
                   int deg,
                   MatrixType coefs, // makes a local copy of coefs
                   MatrixType & result )
    { 
        result.resize(deg+1, 1 ) ;
        T tmp;    
        
        for ( int r=0; r!= deg; r++ ) // triangle step
            //for ( int r=s; r< deg; r++ ) // TO DO: improve using multiplicity s
            for ( int i=deg; i>r; i-- ) // recursive computation
            {
                tmp= ( u -  *(knot+i) ) / ( *(knot+i+deg-r) - *(knot+i) );
                coefs.row(i) = (T(1)-tmp)*coefs.row(i-1) + tmp* coefs.row(i) ;
            }
        result.noalias() = coefs.row(deg);
    }
    
    /// Input: a sequence of 2*p+1 knots starting at iterator position \a knot,
    /// evaluation point \a u, output the value of all basis functions
    /// which are active on \a u.
    /// Uses DeBoor's algorithm for evaluation
    template <class T, typename KnotIterator, typename MatrixType>
    void evalBasis2(  const T & u,
                      const KnotIterator & knot,
                      int deg,
                      MatrixType & result )
    {
        MatrixType pseudo_coefs;
        pseudo_coefs.setIdentity(deg+1,deg+1);
        evalGeoAlg(u, knot, deg, pseudo_coefs, result);
    }
    
    
    /// Input: a sequence of p+2 knots, evaluation point. Output: the
    /// value of the basis function supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void evalBasisSingle( ) 
    { 


    }

    /// Input: a sequence of 2*p+1 knots, evaluation point, output
    /// the derivatives evaluated at all functions supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void derivBasis( ) 
    { }

    /// Input: a sequence of p+2 knots, evaluation point. Output: the
    /// value of the derivative of basis function supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void derivBasisSingle( ) 
    { }

    /// Input: a sequence of 2*p+1 knots, evaluation point, output
    /// the value of all derivatives up to order k supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void allDersBasis( ) 
    { }

    /// Input: a sequence of p+2 knots, evaluation point, output
    /// the value of all derivatives up to order k supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void allDersSingle( ) 
    { }



namespace detail
{

/// Weights of the degree elevation of a Bezier curve from degree p to p+m.
///
/// Row i, column j of the (p+m+1) x (p+1) result is the hypergeometric weight
/// C(p,j) C(m,i-j) / C(p+m,i): the i-th control point of the elevated Bezier
/// curve is the sum over j of this weight times the j-th original control
/// point. All weights are nonnegative and every row sums to one.
///
/// Each row is generated from the ratio of consecutive weights, started from a
/// product of factors in (0,1]; no binomial coefficient or factorial is
/// formed, so there is no overflow for any m (a binomial like C(p+m,i) exceeds
/// the range of float already for moderate m).
template<class T>
gsMatrix<T> bezierElevationWeights(const short_t p, const short_t m)
{
    const index_t ph = p + m;
    gsMatrix<T> w = gsMatrix<T>::Zero(ph + 1, p + 1);
    for (index_t i = 0; i <= ph; ++i)
    {
        const index_t jlo = (std::max)(index_t(0), i - m);
        const index_t jhi = (std::min)(index_t(p), i);
        // weight at j = jlo: C(m,i)/C(ph,i) if i <= m, else C(p,i-m)/C(ph,i)
        // = C(p,ph-i)/C(ph,ph-i)
        T h = T(1);
        if (i <= m)
            for (index_t l = 0; l < i; ++l)
                h *= static_cast<T>(m - l) / static_cast<T>(ph - l);
        else
            for (index_t l = 0; l < ph - i; ++l)
                h *= static_cast<T>(p - l) / static_cast<T>(ph - l);
        w(i, jlo) = h;
        for (index_t j = jlo; j < jhi; ++j)
        {
            h *= static_cast<T>(p - j) / static_cast<T>(j + 1);
            h *= static_cast<T>(i - j) / static_cast<T>(m - i + j + 1);
            w(i, j + 1) = h;
        }
    }
    return w;
}

} // namespace detail


/// @brief Increase the degree of a 1D B-spline from degree p to degree p + m,
/// keeping the curve.
///
/// \tparam Basis_t a 1D B-spline basis (gsBSplineBasis); its knot vector is
/// replaced by the one of degree p+m, with every multiplicity raised by m.
/// \param basis the basis; elevated in place by basis.degreeElevate(m)
/// \param coefs one row per basis function and one column per coordinate; on
/// exit the coefficients of the elevated basis (basis.size() rows). All
/// operations act on whole rows, so any number of columns is efficient.
/// \param m the number of degrees to add
///
/// Precondition: the first and the last knot have multiplicity p+1 (clamped
/// ends). Interior multiplicities 1,...,p+1 are supported; multiplicity p+1
/// is a C^{-1} break.
///
/// Algorithm: Piegl, Tiller, The NURBS Book, 2nd ed., Springer 1997, Sect. 5.5,
/// Algorithm A5.9 (DegreeElevateCurve). Per knot interval [u_a,u_b]:
///  1. Bezier extraction. The control points of the polynomial piece in
///     Bezier form are obtained by inserting u_b into the local window of the
///     coefficients until it has multiplicity p (Boehm's algorithm); each
///     inserted point is a convex combination of two neighbouring points.
///  2. Bezier degree elevation. The p+m+1 Bezier points of the piece in
///     degree p+m are the convex combinations of the p+1 original ones with
///     the weights C(p,j) C(m,i-j) / C(p+m,i).
///  3. Knot removal. The elevated pieces form a spline of degree p+m in which
///     the junction u_a has been inserted up to full Bezier multiplicity
///     p+m. The elevated spline has multiplicity mu+m there (mu the original
///     multiplicity), so p-mu copies of u_a are exactly removable. They are
///     removed by Tiller's algorithm (A5.8): the p-mu+1 coefficients
///     at the junction are recovered from the Bezier points on its two
///     sides, the left half from the points of the left piece and the right
///     half from those of the right piece, so that each removal chain is
///     about (p-mu)/2 steps long. Junctions with mu = p (C^0) and mu = p+1
///     (C^{-1}) need no removal.
///
/// Why it is stable: an output coefficient is computed from the p+1
/// coefficients of one local window by steps 1 and 2 (convex combinations,
/// no cancellation) and then by at most about (p-mu)/2 steps of step 3 on
/// each side. The removal steps invert insertion steps of the elevated
/// curve; they are exact in exact arithmetic (the knots are exactly
/// removable), and their weights are ratios of knot differences, which are
/// bounded in terms of p, m and the local mesh ratio (the ratio of
/// neighbouring interval lengths). On uniform and quasi-uniform knot
/// vectors the amplification of round-off is therefore a constant, and the
/// error of every coefficient is O((p+m)^2 eps max|c|), independent of the
/// number of coefficients. On meshes whose interval lengths grow strongly
/// in one direction the removal weights grow with the mesh ratio, the
/// removal recurrence rewrites coefficients produced by earlier removals,
/// and round-off can grow along the curve.
///
/// A scheme which instead integrates the derivative coefficients
/// (the divided differences of the coefficients, which are the coefficients of
/// the derivatives up to a factor p!/(p-j)!) from the clamped start across all
/// knot intervals is not stable: the integration constants of the low
/// derivative orders are never re-anchored, so each interval adds a
/// round-off error of the size eps*|c^{(p)}|, which is ~(1/h)^p times larger
/// than eps*|c|, to an integral that is carried over O(n) intervals.
///
/// Complexity: O(K p (p+m) d) operations with K the number of knot intervals
/// and d the number of columns of \a coefs; O(p (p+m)) additional memory
/// besides the output.
template<class Basis_t>
void degreeElevateBSpline(Basis_t &basis, 
                          gsMatrix<typename Basis_t::Scalar_t> & coefs,
                          short_t m)
{
    typedef typename Basis_t::Scalar_t T;

    GISMO_ASSERT(m >= 0 && m<512, "Can only elevate degree by a positive (not enormous) amount.");
    GISMO_ASSERT(basis.size() == coefs.rows(), "Invalid coefficients");

    if (m==0) return;

    const short_t p      = basis.degree();
    const index_t d      = coefs.cols();

    // Knots of the original basis (basis.degreeElevate(m) replaces them)
    const std::vector<T> U(basis.knots().begin(), basis.knots().end());
    const index_t M = static_cast<index_t>(U.size()) - 1; // last knot index
    GISMO_ASSERT(U[0] == U[p] && U[M] == U[M - p],
                 "The knot vector must be clamped (end multiplicity p+1).");

    basis.degreeElevate(m);
    const short_t ph      = basis.degree(); // p + m
    const index_t nNew    = basis.size();

    if (0 == p)
    {
        // piecewise constants: the degree-m Bezier points of a constant are all equal to it
        GISMO_ASSERT(nNew == coefs.rows() * (m + 1), "Invalid degree 0 elevation.");
        gsMatrix<T> Q0(nNew, d);
        for (index_t k = 0; k < coefs.rows(); ++k)
            Q0.middleRows(k * (m + 1), m + 1) = coefs.row(k).replicate(m + 1, 1);
        coefs.swap(Q0);
        return;
    }

    const gsMatrix<T> bezw = detail::bezierElevationWeights<T>(p, m);

    gsMatrix<T> Q(nNew, d);              // elevated coefficients
    std::vector<T> Uh(ph + 1, U[0]);     // elevated knots, built along
    Uh.reserve(basis.knots().size());

    gsMatrix<T> bpts(p + 1, d);          // Bezier points of the current piece
    gsMatrix<T> ebpts(ph + 1, d);        // the same, elevated
    gsMatrix<T> nextbpts((std::max)(index_t(p), index_t(1)), d); // seeds of the next piece
    std::vector<T> alfs(p);              // insertion weights

    Q.row(0) = coefs.row(0);
    bpts = coefs.topRows(p + 1);

    index_t cind = 1;                    // next row of Q to be written
    index_t a = p, b = p + 1;            // last knot index of u_a and of u_b
    index_t r = -1;                      // knot copies to insert at u_b, p-mult
    T ua = U[0];

    while (b < M)
    {
        // multiplicity of the next distinct knot u_b
        index_t i = b;
        while (b < M && U[b] == U[b + 1]) ++b;
        const index_t mul = b - i + 1;
        const T ub = U[b];
        const index_t oldr = r;
        r = p - mul;
        // rows lbz..rbz of the elevated Bezier points are kept; those next to
        // the end points belong to the junction coefficients
        const index_t lbz = (oldr > 0) ? (oldr + 2) / 2
                          : ((oldr == 0 || a == p) ? 1 : 0);
        const index_t rbz = (r > 0) ? ph - (r + 1) / 2 : ph;

        // 1. Bezier extraction: insert u_b r times into the window
        if (r > 0)
        {
            const T numer = ub - ua;
            for (index_t k = p; k > mul; --k)
                alfs[k - mul - 1] = numer / (U[a + k] - ua);
            for (index_t j = 1; j <= r; ++j)
            {
                const index_t save = r - j, s = mul + j;
                for (index_t k = p; k >= s; --k)
                {
                    const T al = alfs[k - s];
                    bpts.row(k) = al * bpts.row(k) + (T(1) - al) * bpts.row(k - 1);
                }
                nextbpts.row(save) = bpts.row(p);
            }
        }

        // 2. Bezier degree elevation (only rows lbz..ph are used below)
        for (index_t k = lbz; k <= ph; ++k)
        {
            ebpts.row(k).setZero();
            for (index_t j = (std::max)(index_t(0), k - m); j <= (std::min)(index_t(p), k); ++j)
                ebpts.row(k) += bezw(k, j) * bpts.row(j);
        }

        // 3. Remove the knot u_a oldr times (Tiller's removal, A5.8)
        const index_t kind = static_cast<index_t>(Uh.size());
        if (oldr > 1)
        {
            index_t first = kind - 2, last = kind;
            const T den = ub - ua;
            const T bet = (ub - Uh[kind - 1]) / den;
            for (index_t tr = 1; tr < oldr; ++tr)
            {
                index_t ii = first, jj = last, kj = jj - kind + 1;
                while (jj - ii > tr)
                {
                    if (ii < cind)
                    {
                        const T alf = (ub - Uh[ii]) / (ua - Uh[ii]);
                        Q.row(ii) = alf * Q.row(ii) + (T(1) - alf) * Q.row(ii - 1);
                    }
                    if (jj >= lbz)
                    {
                        if (jj - tr <= kind - ph + oldr)
                        {
                            const T gam = (ub - Uh[jj - tr]) / den;
                            ebpts.row(kj) = gam * ebpts.row(kj) + (T(1) - gam) * ebpts.row(kj + 1);
                        }
                        else
                            ebpts.row(kj) = bet * ebpts.row(kj) + (T(1) - bet) * ebpts.row(kj + 1);
                    }
                    ++ii; --jj; --kj;
                }
                --first; ++last;
            }
        }

        // the knot u_a with its elevated multiplicity
        if (a != p)
            for (index_t k = 0; k < ph - oldr; ++k)
                Uh.push_back(ua);

        for (index_t j = lbz; j <= rbz; ++j)
            Q.row(cind++) = ebpts.row(j);

        if (b < M)
        {
            // set up the next piece
            for (index_t j = 0; j < r; ++j)
                bpts.row(j) = nextbpts.row(j);
            for (index_t j = (std::max)(r, index_t(0)); j <= p; ++j)
                bpts.row(j) = coefs.row(b - p + j);
            a = b; ++b; ua = ub;
        }
        else
            for (index_t k = 0; k <= ph; ++k)
                Uh.push_back(ub);
    }

    GISMO_ASSERT(cind == nNew && static_cast<index_t>(Uh.size()) == basis.knots().size(),
                 "Degree elevation produced an inconsistent number of coefficients.");
    coefs.swap(Q);
}


} // namespace bspline

}// namespace gismo
