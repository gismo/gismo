/** @file gsTensorUtils.h

    @brief Utility functions related to tensor-structured objects

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Matucci
*/

#pragma once

#include <gsTensor/gsTensorTools.h>
#include <gsMatrix/gsFiberMatrix.h>

//Note: hadamard product is just m1.array() * m2.array()
namespace gismo
{
template<typename Derived>
inline bool is_finite(const gsEigen::MatrixBase<Derived>& x)
{
   return ( (x - x).array() == (x - x).array()).all();
}

template<typename Derived>
inline bool is_nan_free(const gsEigen::MatrixBase<Derived>& x)
{
    return ((x.array() == x.array())).all();
}
}


namespace gismo
{

template<class T>
gsVector<index_t> tuckerRank(const std::vector<gsMatrix<T> > & sk)
{
    const index_t d = sk.size();
    gsVector<index_t> res(d);
    for (index_t i = 0; i!=d; ++i)
        res[i] = sk[i].cols();
    return res;
}

// Computes unf * KP(skPtr)
template<class T>
void multiplyKp(const gsMatrix<T> & unf,
                const std::vector<gsMatrix<T>*> & skPtr,
                gsMatrix<T> & result)
{
    const index_t d = skPtr.size() + 1;

    // const index_t n = unf.cols(); // many rhs !

    gsVector<index_t> nc(skPtr.size());
    for(index_t k = 0; k!=d; ++k) nc[k] = skPtr[k]->cols();
    result.resize(unf.rows(), nc);

    result = unf;
    for(index_t k = 0; k+1!=d; ++k)
    {
        //result = result * (skPtr
    }

}

// Computes KP(skPtr) * mat
template<class T>
void multiplyKp(const std::vector<gsMatrix<T> > & cwise,
                const gsMatrix<T> & mat,
                gsMatrix<T> & result)
{
    const index_t rank = cwise.front().cols() / cwise.front().rows();
    index_t sz = 1;
    for(size_t i = 0; i!= cwise.size(); ++i)
        sz *= cwise[i].rows();

    gsMatrix<T,-1,-1,ColMajor> q0, q1;

    const index_t n = mat.cols();
    result.setZero(sz, n);
    for (index_t r = 0; r != rank; ++r)
    {
        q0 = mat;

        for (size_t i = 0; i != cwise.size(); ++i)
        {
            const index_t sz_i  = cwise[i].rows();//cols
            const index_t r_i   = sz / sz_i;
            const index_t rsz_i = r * sz_i;

            q0.resize(sz_i, n * r_i );
            q1.resize(r_i , n * sz_i);

            // Multiple columns
            for ( index_t k = 0; k!=n; ++k)
                q1.middleCols(k*sz_i, sz_i).noalias() =
                    (  cwise[i].block(0, rsz_i, sz_i, sz_i)
                       * q0.middleCols(k*r_i, r_i) ).transpose();

            // Single column
            //q1.noalias() = (m_cwise[i].block(0, rsz_i, sz_i, sz_i) * q0).transpose();

            q1.swap( q0 );
        }
        q0.resize(sz, n);
        result += q0;
    }
}

// Computes KP(skPtr) * mat
template<class T>
void multiplyKpInvUpperTr(const std::vector<gsMatrix<T> > & cwise,
                const gsMatrix<T> & mat,
                   gsMatrix<T> & result)
{
    const index_t rank = cwise.front().cols() / cwise.front().rows();
    index_t sz = 1;
    for(size_t i = 0; i!= cwise.size(); ++i)
        sz *= cwise[i].rows();

    gsMatrix<T,-1,-1,ColMajor> q0, q1;

    const index_t n = mat.cols();
    result.setZero(sz, n);
    for (index_t r = 0; r != rank; ++r)
    {
        q0 = mat;

        for (size_t i = 0; i != cwise.size(); ++i)
        {
            const index_t sz_i  = cwise[i].rows();//cols
            const index_t r_i   = sz / sz_i;
            const index_t rsz_i = r * sz_i;

            q0.resize(sz_i, n * r_i );
            q1.resize(r_i , n * sz_i);

            // Multiple columns
            for ( index_t k = 0; k!=n; ++k)
                q1.middleCols(k*sz_i, sz_i).noalias() =
                    cwise[i].block(0, rsz_i, sz_i, sz_i).
                    template triangularView<gsEigen::Upper>().solve(q0.middleCols(k*r_i, r_i)).transpose();

            // Single column
            //q1.noalias() = m_cwise[i].block(0, rsz_i, sz_i, sz_i).fullPivLu().solve(q0).transpose();

            q1.swap( q0 );
        }
        q0.resize(sz, n);
        result += q0;
    }
}


// Computes skMat.transpose() * unf * KP(skPtr)
template<class T>
void tuckerCore(const gsMatrix<T> & skMat,
                const gsMatrix<T> & unf,
                const std::vector<gsMatrix<T>*> & skPtr,
                const gsVector<index_t> & perm,
                gsMatrix<T> & core)
{
    const index_t d = perm.size();
    gsVector<index_t> psz(d);
    //for(index_t k = 0; k!=d; ++k) psz[k] = rk[perm[k]];

    // Compute the core tensor
    // (r1 x n1) * (n1 x n2*n3) * (n2*n3 x r2*r3)
    //core.noalias() = skeleton.front().transpose() * unf.front() * tmp;
    //core.resize(pr, 1); for.. permute psz
    //permuteTensorVector(perm, psz, core);// note: perm is the one for unf[0]
    // note: Inv(perm)=perm
}


inline int bandWidth(const gsSparseMatrix<> & mat, index_t upto)
{
    upto = math::min(upto, mat.cols());
    int nz = 0;
    for (index_t i = 0; i != upto; ++i)
    {
        const index_t ni = mat.col(i).nonZeros();
        if (ni>nz) nz = ni;
    }
    return nz;
}

inline int bandWidth(const gsSparseMatrix<> & mat)
{ return bandWidth(mat, mat.cols()); }

template<class T>
void reserveKroneckerProduct_V0(const std::vector<gsSparseMatrix<T> > & m,
                             gsSparseMatrix<T> & result)
{
    index_t nzPerCol = 1;
    for (size_t i = 0; i != m.size(); ++i)
        nzPerCol *= bandWidth(m[i], m[i].rows());

    //result.reserve(nzPerCol); // this seems to be globally for all columns
    result.reservePerColumn(nzPerCol);
/*
    gsVector<index_t> nzPerCol(result.cols());
    for ( index_t k=0; k!=result.cols(); ++k )
        for ( size_t i=0; i!=m.size(); ++i )
            nzPerCol[k] *= m[i].col(k).nonZeros();
    result.reserve(nzPerCol);
*/
}

template<class T>
void reserveKroneckerProduct_V0(const std::vector<gsSparseMatrix<T> > & m,
                             gsFiberMatrix<T> & result)
{
    index_t nzPerCol = 1;
    for (size_t i = 0; i != m.size(); ++i)
        nzPerCol *= bandWidth(m[i], m[i].rows());

    //result.reserve(nzPerCol); // this seems to be globally for all columns
    result.reservePerColumn(nzPerCol);
    /*
        gsVector<index_t> nzPerCol(result.cols());
        for ( index_t k=0; k!=result.cols(); ++k )
            for ( size_t i=0; i!=m.size(); ++i )
                nzPerCol[k] *= m[i].col(k).nonZeros();
        result.reserve(nzPerCol);
    */
}


template<class T>
void reserveKroneckerProduct_V0(const std::vector<gsSparseMatrix<T>*> & m,
                             gsSparseMatrix<T> & result)
{
    index_t nzPerCol = 1;
    for (size_t i = 0; i != m.size(); ++i)
        nzPerCol *= bandWidth(*m[i], m[i]->rows());
    gsDebugVar(nzPerCol);
    result.reservePerColumn(nzPerCol);

}

template<class T, int D>
void reserveKroneckerProduct(const std::vector<gsFiberMatrix<T>> & m,
                           std::vector<index_t> & result)
{

    typedef typename gsFiberMatrix<T>::iterator cIter;
    const index_t n = (D==-1 ? m.size() : D);// how to do static??
    GISMO_ASSERT(n>1, "Need at least 2 matrices");

    gsVector<index_t,D> rsz(n);

    for( index_t i = 0; i!= n; ++i)
        rsz [i] = m[i].rows();// row sizes (used for columns as well)

    result.resize( rsz.prod(), 0);


    gsVector<index_t,D> c(n);
    std::vector<cIter> r(n);
    c.setZero();
    index_t t;

    //gsGridIterator<cIter,CUBE,D,true> it(r);
    index_t cc = -1;

    do // for each column
    {
        // const index_t cc = c.dot(rstr); // rows==cols
        ++cc;

        for( index_t i = 0; i!=n; ++i)
            r[i] = m[i].begin(c[i]);

        do // for each row
        {
            ++result[cc];
            // GISMO_ASSERT(tmp - tmp ==tmp - tmp, "Inf");
            // GISMO_ASSERT(tmp == tmp, "Nan");

            // start next position
            for (t = 0; t != n; ++t)
            {
                if ( ++r[t] )
                    break;
                else
                    if (t == n - 1)
                    { ++t; break;}
                    else
                        r[t] = m[t].begin(c[t]);
            }
            // end next position

        } while (t!=n);

    } while (nextLexicographic(c, rsz));//assuming square blocks
}

// Selects most important singular values
template<typename Vec>
index_t gsPickHead(Vec & cont, const typename Vec::Scalar tol)
{
    typedef typename Vec::Scalar Scalar_t;
    const Scalar_t * sv0 = cont.data();

// Error in spectral norm:
//    Scalar_t * svR = std::lower_bound( sv0 ,sv0 + cont.rows(),
//                                   (*sv0) * tol, // relative tolerance (spectral norm)
//                                   //tol,        // absolute tolerance (spectral norm)
//                                   std::greater<Scalar_t>() );

    //Pick according to the Frobenius norm
    index_t r = cont.rows(); //Maximal rank (not index)
    Scalar_t err = 0.0;
    //Scalar_t errBefore = 0.0;
    const Scalar_t tol2 = tol; //Now, tolerance is squared in gsSVDDecomposition, if we do it here we have to change the full decomposition to use squared tolerance
    //gsDebug << "Squared tolerance: " << std::scientific << tol*tol << std::fixed << "\n";
    if(tol2 !=0.0)
    {
        while (err < tol2 && r >= 1)
        {
            //errBefore = err;
            err += (*(sv0 + r - 1)) * (*(sv0 + r - 1));

            // gsDebugVar(r);
            // gsDebug << std::scientific;
            // gsDebugVar(err);
            // gsDebug << std::fixed;

            r--;
        }

        if (err >= tol2)    //If the while loop has stopped because the error became too big
                            // we have to go one step backwards,
                            // if not, the rank is 0
        {
            r++;
            //err = errBefore;
        }
    }


    gsDebug<<"Picked "<< r <<" first of "<< std::scientific << cont.head( math::min<index_t>(r+2, cont.rows())).transpose()
            << ". Squared Residual error in Frobenius norm: " << std::scientific << err <<"(tol: "<<tol<<")\n" << std::fixed;

    //gsDebug << "Picked " << r << " first singular values.\n" << "Squared Residual error in Frobenius norm: " << std::scientific << err << "\n" << std::fixed;
    return r;


        /* // better since we expect small rank ?
          Scalar_t * svR = std::find_if( sv0 ,sv0 + cont.rows(),
                                       //(*sv0) * tol, // relative tolerance
                                       std::bind1st(std::greater<Scalar_t>(),tol));
        */

        //gsInfo<<"Picked "<< svR - sv0 <<" first of "<<
        //    cont.head( math::min<index_t>((svR-sv0+2), cont.rows())).transpose() <<"..\n";
        //return svR - sv0;
}

/*
// Selects most important singular values and **updates** truncation
// error
template<typename Vec>
index_t gsPickHead(Vec & cont, const typename Vec::Scalar tol,
                   typename Vec::Scalar & truncError)
{
    typedef typename Vec::Scalar Scalar_t;
    const Scalar_t * sv0 = cont.data();
    const Scalar_t * svN = sv0 + cont.rows();
    const Scalar_t * svR = std::lower_bound( sv0 ,svN,
                                             //(*sv0) * tol, // relative tolerance
                                             tol, // absolute tolerance
                                       std::greater<Scalar_t>() );

    truncError = gsVector<Scalar_t>::Map(svR, svN - svR ).norm();
    return svR - sv0;
}
*/

//template<class T, int D>
//getUnfoldings(const gsMatrix<T> & tensor,
//              const gsVector<index_t,D> & sz,
//              std::vector<gsMatrix<T> > & unfoldings, ..


template<class T>
void gsKroneckerProductVectorSum2(const gsMatrix<T> & m1,
                                  const gsMatrix<T> & m2,
                                  gsVector<T> & result)
{
    const index_t s2 = m2.rows(), s1 = m1.rows();
    const index_t rk = m1.cols();
    result.resize(s1*s2, 1);

    for (index_t i = 0; i != s2; i++) // for all rows of m2 -
    {
        for (index_t r = 0; r != s1; r++)  // for all rows of m1
        {

            T tmp = m1(r,0)*m2(i,0);
            for (index_t t = 1; t != rk; ++t)
                tmp +=  m1(r,t)*m2(i,t);
            result(i*s1+r) = tmp;

        }
    }


}


// 2D Dense Kronecker product-sum of many matrices stacked one after the other
// (note: not needed. this is only needed for Sparse matrices (for global comp.) )
template<class T>
void gsKroneckerProductSum2(const gsMatrix<T> & m1,
                            const gsMatrix<T> & m2,
                            gsMatrix<T> & result)
{
    // Assumes square coordinate matrices
    const index_t s2 = m2.rows(),
                  s1 = m1.rows();
    const index_t rk = m1.cols() / s1;

    result.resize (s1*s2, s1*s2);

    for (index_t i = 0; i != s2; i++) // for all rows of m2 -
        for (index_t j = 0; j != s2; j++) // for all cols of m2 -
            for (index_t r = 0; r != s1; r++)  // for all rows of m1
                for (index_t c = 0; c != s1; c++) // for all cols of m1
                {
                    T tmp = m1(r,c)*m2(i,j);
                    for (index_t t = 1; t != rk; ++t)
                        tmp +=  m1(r,t*s1 + c)*m2(i,t*s2+j);
                    //result.block(i*s1, j*s1, s1, s1) = m2(i,j) * m1 ;
                    result(i*s1+r, j*s1+c) = tmp;
                }
}

// 2D Dense Kronecker product
template<class T>
void gsKroneckerProduct2(const gsMatrix<T> & m1,
                         const gsMatrix<T> & m2,
                         gsMatrix<T> & result)
{
    const index_t r  = m1.rows(), c  = m1.cols();
    const index_t r2 = m2.rows(), c2 = m2.cols();
    result.resize(r*r2, c*c2);

    for (index_t i = 0; i != r2; i++) // for all rows of m2
        for (index_t j = 0; j != c2; j++) // for all cols of m2
            result.block(i*r, j*c, r, c) = m2(i,j) * m1 ;
}

// 2D Dense Khatri-Rao product
template<class T>
void gsKhatriRaoProduct(const gsMatrix<T> & m1,
                        const gsMatrix<T> & m2,
                        gsMatrix<T> & result)
{
    const index_t r  = m1.rows(), c = m1.cols();
    const index_t r2 = m2.rows();
    GISMO_ASSERT(c==m2.cols(), "Column sizes do not match");

    result.resize(r*r2, c);

    for (index_t j = 0; j != c; j++) // for all cols of m2
        for (index_t i = 0; i != r2; i++) // for all rows of m2
            result.block(i*r, j, r, 1) = m2(i,j) * m1.col(j) ;
}

/// Dense Kronecker product of a set of matrices
/// We adopt a definition of KP ..
///
template<class T>
void gsKroneckerProduct(const std::vector<gsMatrix<T>*> & m, gsMatrix<T> & result)
{
    //typdef typename std::vector<gsMatrix<T> >::const_iteretor miter;
    const index_t n = m.size() - 1;
    if ( n==0 )
    {   // Special case of a single matrix, copy to output
        result = *m.front();
        return;
    }

    const index_t lr = m.front()->rows(),
                  lc = m.front()->cols();

    gsVector<index_t> rsz(n), csz(n), rstr(n), cstr(n);
    for( index_t i = 0; i!= n; ++i)
    {
        // sizes and strides
        rsz [i] = m[i+1]->rows();
        csz [i] = m[i+1]->cols();
        rstr[i] = lr;
        cstr[i] = lc;
        for( index_t j = i; j!=0; --j)
        {
            rstr[i] *= m[j]->rows();
            cstr[i] *= m[j]->cols();
        }
    }

    result.resize( lr*rsz.prod(), lc*csz.prod() );

    gsVector<index_t> r(n), c = gsVector<index_t>::Zero(n);
    do // for each column
    {
        const index_t cc = c.dot(cstr);
        r.setZero();
        do // for each row
        {
            T tmp = m[1]->coeff(r[0],c[0]);
            for( index_t i = 1; i!=n; ++i)
                tmp *= m[i+1]->coeff(r[i],c[i]);

            result.block( r.dot(rstr), cc, lr, lc).noalias() =
                tmp * (*m.front());

        } while (nextLexicographic(r, rsz));
    } while (nextLexicographic(c, csz));
}


template<class T>
void gsKroneckerProductStandard(const std::vector<gsMatrix<T>*> & m, gsMatrix<T> & result)
{
    std::vector<gsMatrix<T>*> m_r = m;
    std::reverse(m_r.begin(),m_r.end());
    gsKroneckerProduct(m_r,result);
}

/*
  Computes the (reversed) Kronecker product
 */
template<class T>
void gsKroneckerProduct(const std::vector<gsMatrix<T> > & m, gsMatrix<T> & result)
{
    gsKroneckerProduct(asVectorPtr(m), result);
}

template<class T> // TODO: merge with ptr version
void gsKroneckerProduct(const std::vector<gsSparseMatrix<T> > & m, gsSparseMatrix<T> & result)
{
/* // to test: is this faster ? --> also: make static dim
   gsLRSparseMatrix<T> lrm(m);
   result.resize( lrm.rows(), lrm.cols() );
   reserveKroneckerProduct_V0(m, result);

   index_t r, c;
   for(typename gsLRSparseMatrix<T>::iterator it(m); it; ++it)
   {
       const T val = it.value(r,c);
       result.insert(r, c) = val;
   }
//*/

// /*
    typedef typename gsSparseMatrix<T>::iterator cIter;

    const index_t n = m.size();
    GISMO_ASSERT(n>1, "Need at least 2 matrices");

    gsVector<index_t> rsz(n), csz(n), rstr(n), cstr(n);

    for( index_t i = 0; i!= n; ++i)
    {
        // sizes and strides
        rsz [i] = m[i].rows();
        csz [i] = m[i].cols();
    }

    tensorStrides(rsz, rstr);
    tensorStrides(csz, cstr);

    result.resize( rsz.prod(), csz.prod() );
    reserveKroneckerProduct_V0(m,result);

    gsVector<index_t> c(n);
    std::vector<cIter> r(n);
    c.setZero();
    index_t t;

    do // for each column
    {
        const index_t cc = c.dot(cstr);
        for( index_t i = 0; i!=n; ++i)
            r[i] = m[i].begin(c[i]);

        do // for each row
        {
            index_t rr = r[0].index();
            T tmp      = r[0].value();
            for( index_t i = 1; i!=n; ++i)
            {
                tmp *= r[i].value();
                rr  += rstr[i] * r[i].index();
            }

            result.insert(rr,cc) = tmp;

            // start next position
            for (t = 0; t != n; ++t)
                if ( ++r[t] )
                    break;
                else
                    if (t == n - 1)
                    { ++t; break;}
                    else
                        r[t] = m[t].begin(c[t]);
            // end next position

        } while (t!=n);

    } while (nextLexicographic(c, csz));
//*/
    result.makeCompressed();
}

template<class T>
void gsKroneckerProduct(const std::vector<gsSparseMatrix<T>*> & m,
                        gsSparseMatrix<T> & result)
{
// /*
    typedef typename gsSparseMatrix<T>::iterator cIter;

    const index_t n = m.size();
    GISMO_ASSERT(n>1, "Need at least 2 matrices");

    gsVector<index_t> rsz(n), csz(n), rstr(n), cstr(n);

    for( index_t i = 0; i!= n; ++i)
    {
        // sizes and strides
        rsz [i] = m[i]->rows();
        csz [i] = m[i]->cols();
    }

    tensorStrides(rsz, rstr);
    tensorStrides(csz, cstr);

    result.resize( rsz.prod(), csz.prod() );
    reserveKroneckerProduct_V0(m,result);

    gsVector<index_t> c(n);
    std::vector<cIter> r(n);
    c.setZero();
    index_t t;

    do // for each column
    {
        const index_t cc = c.dot(cstr);
        for( index_t i = 0; i!=n; ++i)
            r[i] = m[i]->begin(c[i]);

        do // for each row
        {
            index_t rr = r[0].index();
            T tmp      = r[0].value();
            for( index_t i = 1; i!=n; ++i)
            {
                tmp *= r[i].value();
                rr  += rstr[i] * r[i].index();
            }

            result.insert(rr,cc) = tmp;

            // start next position
            for (t = 0; t != n; ++t)
                if ( ++r[t] )
                    break;
                else
                    if (t == n - 1)
                    { ++t; break;}
                    else
                        r[t] = m[t]->begin(c[t]);
            // end next position

        } while (t!=n);

    } while (nextLexicographic(c, csz));
//*/
    result.makeCompressed();
}

template<class T>
void gsKroneckerProduct(const gsSparseMatrix<T> & m1,
                        const gsSparseMatrix<T> & m2,
                        gsSparseMatrix<T> & result)
{
    GISMO_ERROR("to do");
}

// Dense Kronecker product-sum of many matrices stacked one after the other
// generalizes the case of 2 matrices
// (note: not needed. this is only needed for Sparse matrices (for global comp.) )
template<class T>
void gsKroneckerProductSum(const std::vector<gsMatrix<T> > & m, gsMatrix<T> & result)
{
    GISMO_ERROR("No implementation");
}

template<class T, int D>
void gsKroneckerProductSum(const std::vector<gsSparseMatrix<T> > & m,
                           gsSparseMatrix<T> & result)
{
/* // to test: is this faster ? --> currently NOT: make static dim!
   index_t n = 1;
   for( size_t i = 0; i!= m.size(); ++i)
   n *= m[i].rows();// row sizes (used for columns as well)

   //gsLRSparseMatrix<T,D> lrm(m);
   result.resize(n, n);
   reserveKroneckerProduct_V0(m, result);

   index_t r, c;
   for(typename gsLRSparseMatrix<T>::iterator it(m); it; ++it)
   {
       const T val = it.value(r,c);
       result.insert(r, c) = val;
   }
//*/

///*
    // idea 1: make specific versions for 2D, 3D using for loops --> loop unrolling
    // idea 2: cache nested multiplications within the computation

    typedef typename gsSparseMatrix<T>::iterator cIter;
    const index_t n = (D==-1 ? m.size() : D);// how to do static??
    GISMO_ASSERT(n>1, "Need at least 2 matrices");

    gsVector<index_t,D> rsz(n), rstr(n);

    for( index_t i = 0; i!= n; ++i)
        rsz [i] = m[i].rows();// row sizes (used for columns as well)

    const index_t rk = m.front().cols() / rsz[0];// assuming square blocks
    const index_t rkn = n * rk;

    tensorStrides(rsz,rstr);

    result.resize( rsz.prod(), rsz.prod() );
    reserveKroneckerProduct_V0(m,result);

    gsVector<index_t,D> c(n);
    std::vector<cIter> r(rkn);
    //gsVector<cIter,D> r(rkn);
    c.setZero();
    index_t t;

    //gsGridIterator<cIter,CUBE,D,true> it(r);
    index_t cc = -1;

    do // for each column
    {
        // const index_t cc = c.dot(rstr); // rows==cols
        ++cc;

        for( index_t i = 0; i!=n; ++i)
            for( index_t k = 0; k!=rkn; k+=n)
                r[k+i] = m[i].begin((k/n)*rsz[i] + c[i]);

        do // for each row
        {
            index_t rr = r[0].index();// row index in the result
            T tmp      = r[0].value();
            for( index_t i = 1; i!=n; ++i)
            {
                tmp *= r[i].value();
                rr  += rstr[i] * r[i].index();
            }

            for( index_t k = n; k!=rkn; k+=n)
            {
                T tmpr = r[k].value();
                for( index_t i = 1; i!=n; ++i)
                    tmpr *= r[k+i].value();
                tmp += tmpr;
            }

            result.insert(rr,cc) = tmp;

            // GISMO_ASSERT(tmp - tmp ==tmp - tmp, "Inf");
            // GISMO_ASSERT(tmp == tmp, "Nan");

            // start next position
            for (t = 0; t != n; ++t)
            {
                // assuming common structure for all skeleton matrices
                for( index_t k = n; k!=rkn; k+=n) ++r[k+t];
                if ( ++r[t] )
                    break;
                else
                    if (t == n - 1)
                    { ++t; break;}
                    else
                        for( index_t k = 0; k!=rkn; k+=n)
                        r[k+t] = m[t].begin((k/n)*rsz[t] + c[t]);
            }
            // end next position

        } while (t!=n);

    } while (nextLexicographic(c, rsz));//assuming square blocks
//*/

    // Done outside *****
    //result.makeCompressed();
}

template<class T, int D>
void gsKroneckerProductSum(const std::vector<gsFiberMatrix<T> > & m,
                           gsFiberMatrix<T> & result)
{
/* // to test: is this faster ? --> currently NOT: make static dim!
   index_t n = 1;
   for( size_t i = 0; i!= m.size(); ++i)
   n *= m[i].rows();// row sizes (used for columns as well)

   //gsLRSparseMatrix<T,D> lrm(m);
   result.resize(n, n);
   reserveKroneckerProduct_V0(m, result);

   index_t r, c;
   for(typename gsLRSparseMatrix<T>::iterator it(m); it; ++it)
   {
       const T val = it.value(r,c);
       result.insert(r, c) = val;
   }
//*/

///*
    // idea 1: make specific versions for 2D, 3D using for loops --> loop unrolling
    // idea 2: cache nested multiplications within the computation

    typedef typename gsFiberMatrix<T>::iterator cIter;
    const index_t n = (D==-1 ? m.size() : D);// how to do static??
    GISMO_ASSERT(n>1, "Need at least 2 matrices");

    gsVector<index_t,D> rsz(n), rstr(n);

    for( index_t i = 0; i!= n; ++i)
        rsz [i] = m[i].rows();// row sizes (used for columns as well)

    const index_t rk = m.front().cols() / rsz[0];// assuming square blocks
    const index_t rkn = n * rk;

    tensorStrides(rsz,rstr);

    gsVector<index_t,D> c(n);
    std::vector<cIter> r(rkn);
    //gsVector<cIter,D> r(rkn);
    c.setZero();
    index_t t;

    //gsGridIterator<cIter,CUBE,D,true> it(r);
    index_t cc = -1;

    do // for each column
    {
        // const index_t cc = c.dot(rstr); // rows==cols
        ++cc;

        for( index_t i = 0; i!=n; ++i)
            for( index_t k = 0; k!=rkn; k+=n)
                r[k+i] = m[i].begin((k/n)*rsz[i] + c[i]);

        auto& fiber = result.fiber(cc);
        do // for each row
        {
            index_t rr = r[0].index();// row index in the result
            T tmp      = r[0].value();
            for( index_t i = 1; i!=n; ++i)
            {
                tmp *= r[i].value();
                rr  += rstr[i] * r[i].index();
            }

            for( index_t k = n; k!=rkn; k+=n)
            {
                T tmpr = r[k].value();
                for( index_t i = 1; i!=n; ++i)
                    tmpr *= r[k+i].value();
                tmp += tmpr;
            }

            fiber.insertBack(rr) =  tmp;

            // GISMO_ASSERT(tmp - tmp ==tmp - tmp, "Inf");
            // GISMO_ASSERT(tmp == tmp, "Nan");

            // start next position
            for (t = 0; t != n; ++t)
            {
                // assuming common structure for all skeleton matrices
                for( index_t k = n; k!=rkn; k+=n) ++r[k+t];
                if ( ++r[t] )
                    break;
                else
                    if (t == n - 1)
                    { ++t; break;}
                    else
                        for( index_t k = 0; k!=rkn; k+=n)
                        r[k+t] = m[t].begin((k/n)*rsz[t] + c[t]);
            }
            // end next position

        } while (t!=n);

    } while (nextLexicographic(c, rsz));//assuming square blocks
//*/

    // Done outside *****
    //result.makeCompressed();
}

template<class T>
void gsKhatriRaoProduct(const std::vector<gsMatrix<T>*> & m, gsMatrix<T> & result)
{
    const index_t n = m.size() - 1;
    if ( n==0 )
    {   // Special case of a single matrix, copy to output
        result = *m.front();
        return;
    }
    const index_t ncols = m.front()->cols();
    GISMO_ASSERT( ncols == m[1]->cols(), "Expecting the same number of columns");
    const index_t lr = m.front()->rows();

    gsVector<index_t> r = gsVector<index_t>::Zero(n), rsz(n), rstr(n);
    for( index_t i = 0; i!= n; ++i)
    {
        // sizes and strides
        rsz [i] = m[i+1]->rows();
        rstr[i] = lr;
        for( index_t j = i; j!=0; --j)
            rstr[i] *= m[j]->rows();
    }

    result.resize( lr*rsz.prod(), ncols );
    do// for each row
    {
        for( index_t j = 0; j!= ncols; ++j) // for each column
        {
            T tmp = m[1]->coeff(r[0], j);
            for( index_t i = 1; i!=n; ++i)
                tmp *= m[i+1]->coeff(r[i], j);

            result.block( r.dot(rstr), j, lr, 1).noalias() =
                tmp * m.front()->col(j);
        }
    } while (nextLexicographic(r, rsz));
}

template<class T, int D>
void gsKhatriRaoProductTr(const std::vector<gsMatrix<T>*> & m, gsMatrix<T> & result)
{
    static const int dm = (D==-1 ? -1 : D-1);
    const index_t n = (D==-1 ? m.size() : D) - 1;

    if ( n==0 )
    {   // Special case of a single matrix, copy (tranposed) to output
        result = m.front()->transpose();
        return;
    }
    const index_t ncols = m.front()->cols();
    GISMO_ASSERT( ncols == m[1]->cols(), "Expecting the same number of columns");
    const index_t lr = m.front()->rows();

    gsVector<index_t,dm> r = gsVector<index_t,D>::Zero(n), rsz(n), rstr(n);
    for( index_t i = 0; i!= n; ++i)
    {
        // sizes and strides
        rsz [i] = m[i+1]->rows();
        rstr[i] = lr;
        for( index_t j = i; j!=0; --j)
            rstr[i] *= m[j]->rows();
    }

    result.resize(ncols, lr*rsz.prod() );
    do// for each row
    {
        for( index_t j = 0; j!= ncols; ++j) // for each column
        {
            T tmp = m[1]->coeff(r[0], j);
            for( index_t i = 1; i!=n; ++i)
                tmp *= m[i+1]->coeff(r[i], j);

            result.block(j, r.dot(rstr), 1, lr).noalias() = // transposed
                (tmp * m.front()->col(j)).transpose();
        }
    } while (nextLexicographic(r, rsz));
}

template<class T>
void gsKhatriRaoProduct(const std::vector<gsMatrix<T> > & m, gsMatrix<T> & result)
{
    gsKhatriRaoProduct(asVectorPtr(m), result);
}

//note: not needed
//template<class T>
//void gsKhatriRaoProduct(const std::vector<gsSparseMatrix<T> > & m, gsSparseMatrix<T> & result)

template<class T>
T gsTuckerDiffNormSqr(const gsMatrix<T> & tensor,
                      const std::vector<gsMatrix<T> > & skeleton,
                      const gsMatrix<T> & core)
{

    gsMatrix<T> tmp;
    gsTucker2Full(tuckerRank(skeleton),skeleton,core, tmp);
    return (tmp-tensor).squaredNorm();

    const index_t d = skeleton.size();
    gsVector<index_t> rk(d), r(d), rstr(d), s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
    {
        sz[i] = skeleton[i].rows();
        rk[i] = skeleton[i].cols();
    }

    tensorStrides(sz,sstr);
    tensorStrides(rk,rstr);

    T result = 0.0;
    s.setZero();
    do
    {
        T cur = 0.0;
        cur = 0.0;
        r.setZero();
        do
        {
            T tmp = core( r.dot(rstr), 0);
            for (index_t i = 0; i!=d; ++i)
                tmp *= skeleton[i](s[i],r[i]);
            cur += tmp;

        } while (nextLexicographic(r, rk));

        const T diff = cur - tensor(s.dot(sstr),0);
        result += diff*diff;

    } while (nextLexicographic(s, sz));

    return result;
}

template<class T>
T gsTuckMaxDiff(const gsMatrix<T> & tensor,
                const std::vector<gsMatrix<T> > & skeleton,
                const gsMatrix<T> & core)
{
    GISMO_ASSERT(0!=core.size(), "Core is empty");

    const index_t d = skeleton.size();
    gsVector<index_t> rk(d), r(d), rstr(d), s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
    {
        sz[i] = skeleton[i].rows();
        rk[i] = skeleton[i].cols();
    }

    tensorStrides(sz,sstr);
    tensorStrides(rk,rstr);

    T result = 0.0;
    s.setZero();
    do
    {
        T cur = 0.0;
        cur = 0.0;
        r.setZero(); // gsLowRankTensor
        do
        {
            T tmp = core( r.dot(rstr), 0);
            for (index_t i = 0; i!=d; ++i)
                tmp *= skeleton[i](s[i],r[i]);
            cur += tmp;

        } while (nextLexicographic(r, rk));

        const index_t ss = s.dot(sstr);
        result = math::max(result,
                           math::abs(tensor.at(ss)-cur) );

        const T diff = cur - tensor(s.dot(sstr),0);
        result += diff*diff;

    } while (nextLexicographic(s, sz));

    return result;
}


// Tucker fitting (approximation) with fixed rank (HOOI)
template<class T>
void gsTuckerDecFixed(const gsMatrix<T> & tensor,
                      const gsVector<int> & sz,
                      const gsVector<int> & rk,// rank
                      std::vector<gsMatrix<T> > & skeleton,
                      gsMatrix<T> & core,
                      T tol = 1e-11,
                      index_t kmax = 10)
{
    const index_t d   = sz.size();
    const index_t dm1 = d-1;
    const index_t pr  = tensor.size();

    GISMO_ASSERT( (rk.array()>0).all(), "Rank "<<rk.transpose()<<" needs to be positive");
    GISMO_ASSERT( (sz.array()>=rk.array()).all(), "Rank "<<rk.transpose()<<" is too big");

    gsEigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    const unsigned svdOpts = gsEigen::ComputeThinU;

    // *** Phase I - Initial Guess taken from truncated HOSVD

    gsVector<index_t> perm(d), psz(d);
    skeleton.resize(d);

    // Store the unfoldings of the input tensor once and for all
    std::vector<gsMatrix<T> > unf(d, tensor );// initially equal to input config
    std::vector<std::vector<gsMatrix<T>*> > skPtr(d);

    for(index_t i = dm1; i!=-1; --i)// reversed to keep the perm we need
    {
        // Compute unfoldings
        psz = sz;
        // Permutation: shift to k, reverse tail
        perm = gsVector<index_t>::LinSpaced(d,0,dm1);
        std::rotate (perm.data()    , perm.data() + i, perm.data() + d);
        std::reverse(perm.data() + 1, perm.data() + d );
        permuteTensorVector(perm, psz, unf[i]);
        unf[i].resize(psz[0], pr / psz[0]);

        // Keep selections of skeletons for later use (in Kronecker products)
        skPtr[i].reserve(dm1);
        for(index_t j = 1; j!=d; ++j)
            skPtr[i].push_back( &skeleton[perm[j]] );

        // Compute initial guess
        svd.compute(unf[i], svdOpts);
        skeleton[i] = svd.matrixU().leftCols(math::min(rk[i],sz[i]) ); // nj x rj
    }

    // *** Phase II - ALS Iteration (higher-order orthogonal iteration)
    gsMatrix<T> tmp, SHT;
    kmax = math::max<index_t>(1,kmax);
    T obj = 0; // objective function: single hole tensor norm
    index_t i;
    for(i = 0; i!=kmax; ++i)
    {
        for(index_t j = dm1; j!=-1; --j)
        {
            //multiplyKP( unf[j], skPtr[j]);
            gsKroneckerProduct(skPtr[j], tmp);
            // (nj x ni*nk)  * (ni*nk x ri*rk)
            SHT.noalias() = unf[j] * tmp;
            svd.compute(SHT, svdOpts);
            skeleton[j] = svd.matrixU().leftCols(skeleton[j].cols());
        }

        // Stop if no improvement in objective
        const T newObj = SHT.squaredNorm();
        if ( math::abs(obj - newObj) <= tol )
            break;
        obj = newObj;
    }
    gsInfo<<"Tucker-ALS iterations: "<< i+1 <<"\n";
    // Compute the core tensor
    core.noalias() = skeleton.front().transpose() * SHT;
    core.resize(core.size(), 1);
    for(index_t k = 0; k!=d; ++k)
        psz[k] = skeleton[perm[k]].cols();
    permuteTensorVector(perm, psz, core);// note: perm is the one for unf[0]
}

// HOOI: Uses truncated HOSVD as starting point and does ALS iterations
// (higher-order orthogonal iteration)
template<class T>
bool gsTuckerDec(const gsMatrix<T> & tensor,
                 const gsVector<int> & sz,
                 std::vector<gsMatrix<T> > & skeleton,
                 gsMatrix<T> & core,
                 T tol = 1e-11,
                 const index_t kmax = 10)
{
    const index_t d  = sz.size();
    const index_t dm1 = d-1;
    const index_t pr  = tensor.size(); //=sz.prod();
    GISMO_ASSERT( tensor.cols()==1, "Expecting a single column" );
    GISMO_ASSERT( pr==sz.prod()   , "Invalid tensor/sizes" );

    // Check for zero tensor
    if ( tensor.array().abs().maxCoeff() <= tol )
    {
        skeleton.reserve(d);
        for(index_t i = 0; i!=d; ++i)
            skeleton[i].setZero(sz[i],1);
        core.setZero(1, 1);
        return true;
    }

    gsEigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    const unsigned svdOpts = gsEigen::ComputeThinU;

    // *** Phase I - Compute unfoldings and HOSVD
    gsVector<index_t> rk(d);// Tucker rank
    gsVector<index_t> perm(d), psz(d);
    skeleton.resize(d);

    // Store the unfoldings of the input tensor once and for all
    // This extra storage is needed only for (faster) ALS iterations
    std::vector<gsMatrix<T> > unf(d, tensor );
    std::vector<std::vector<gsMatrix<T>*> > skPtr(d);

    for(index_t j = dm1; j!=-1; --j)// reversed to keep the perm we need
    {
        // Compute unfoldings, starting from 0,1,..,d-1
        psz = sz;
        // Permutation: shift to k, reverse tail
        perm = gsVector<index_t>::LinSpaced(d,0,dm1);
        std::rotate (perm.data()    , perm.data() + j, perm.data() + d);
        std::reverse(perm.data() + 1, perm.data() + d );
        permuteTensorVector(perm, psz, unf[j]);
        unf[j].resize(psz[0], pr / psz[0]);

        // Keep permuted skeletons for later use in ALS (in Kronecker products)
        skPtr[j].reserve(dm1);
        for(index_t l = 1; l!=d; ++l)
            skPtr[j].push_back( &skeleton[perm[l]] );

        // Compute SVD of unfolding and truncate upto tolerance
        svd.compute(unf[j], svdOpts);
        rk[j] = gsPickHead(svd.singularValues(), tol);
        skeleton[j] = svd.matrixU().leftCols(rk[j]); // nj x rj
    }

    // *** Phase II - ALS iteration
    // At this point we have the HOSVD. Now we apply ALS iterations.
    // This can reduce the Tucker rank even more, and most likely to
    // the optimal numbers for our tolerance
    gsMatrix<T> tmp, SHT;
    bool improved = false;
    for(index_t i = 0; i!=kmax; ++i)
    {
        //gsInfo <<"\n Skeleton:\n" << skeleton.front() <<"\n";
        for(index_t j = dm1; j!=-1; --j)
        {
            //multiplyKP( unf[j], skPtr[j]);
            gsKroneckerProduct(skPtr[j], tmp);
            // "single hole" tensor
            SHT.noalias() = unf[j] * tmp; // (nj x ni*nk) * (ni*nk x ri*rk)
            svd.compute(SHT, svdOpts);

            // Truncate upto tolerance
            const index_t new_rj = gsPickHead(svd.singularValues(), tol);
            if ( new_rj != rk[j] )
                improved = true;
            {   // Note: maybe we need to update all skeletons once one is improved?
                rk[j] = new_rj;
                skeleton[j] = svd.matrixU().leftCols(new_rj); // nj x new_rj
            }
        }

        //if ( (i+1) % 8 == 0 &&  !improved ) break;
        if ( !improved ) break;

        //if ( improved ) gsInfo<<"Tucker Rank reduced to "<<rk.transpose()<<"\n";
    }

    if ( improved )
        gsInfo << "\nTucker: Still improving after "<<kmax
               <<" itererations, try more iterations\n";

    // Compute the core tensor
    if ( 0==kmax ) // If no ALS iterations at all
    {
        gsKroneckerProduct(skPtr.front(), tmp);
        core.noalias() = skeleton.front().transpose() * unf.front() * tmp;
    }
    else // SHT already computed
        core.noalias() = skeleton.front().transpose() * SHT;

    core.resize(core.size(), 1); for(index_t k = 0; k!=d; ++k) psz[k] = rk[perm[k]];
    permuteTensorVector(perm, psz, core);// note: can we avoid permute here (ie. perm=id) ?
    return true; // always true !
}

/* TO DO --> bottom up
template<class T>
bool gsTuckerAlsIncremental(const gsMatrix<T> & tensor,
                            const gsVector<int> & sz,
                            std::vector<gsMatrix<T> > & skeleton,
                            gsMatrix<T> & core,
                            T tol = 1e-11,
                            const index_t kmax = 10,
                            const index_t rmax = 10)
{
    const index_t d  = sz.size();
    const index_t dm1 = d-1;
    const index_t pr  = tensor.size(); //=sz.prod();
    GISMO_ASSERT( tensor.cols()==1, "Expecting a single column" );
    GISMO_ASSERT( pr==sz.prod()   , "Invalid tensor/sizes" );

    // Check for zero tensor
    if ( tensor.array().abs().maxCoeff() <= tol )
    {
        skeleton.reserve(d);
        for(index_t i = 0; i!=d; ++i)
            skeleton[i].setZero(sz[i],1);
        core.setZero(1, 1);
        return true;
    }

    Eigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    const unsigned svdOpts = Eigen::ComputeThinU;

    // *** Phase I - Compute unfoldings and HOSVD
    gsVector<index_t> rk(d);// Tucker rank
    gsVector<index_t> perm(d), psz(d);
    skeleton.resize(d);

    // Store the unfoldings of the input tensor once and for all
    // This extra storage is needed only for (faster) ALS iterations
    std::vector<gsMatrix<T> > unf(d, tensor );
    std::vector<std::vector<gsMatrix<T>*> > skPtr(d);

    for(index_t j = dm1; j!=-1; --j)// reversed to keep the perm we need
    {
        // Compute unfoldings, starting from 0,1,..,d-1
        psz = sz;
        // Permutation: shift to k, reverse tail
        perm = gsVector<index_t>::LinSpaced(d,0,dm1);
        std::rotate (perm.data()    , perm.data() + j, perm.data() + d);
        std::reverse(perm.data() + 1, perm.data() + d );
        permuteTensorVector(perm, psz, unf[j]);
        unf[j].resize(psz[0], pr / psz[0]);

        // Keep permuted skeletons for later use in ALS (in Kronecker products)
        skPtr[j].reserve(dm1);
        for(index_t l = 1; l!=d; ++l)
            skPtr[j].push_back( &skeleton[perm[l]] );

        // Compute SVD of unfolding and truncate upto tolerance
        svd.compute(unf[j], svdOpts);
        rk[j] = gsPickHead(svd.singularValues(), tol);
        skeleton[j] = svd.matrixU().leftCols(rk[j]); // nj x rj
    }

    // *** Phase II - ALS iteration
    // At this point we have the HOSVD. Now we apply ALS iterations, incrementally


// This can reduce the Tucker rank even more, and most likely to
    // the optimal numbers for our tolerance
    gsMatrix<T> tmp, SHT;
    bool improved = false;
    for(index_t i = 0; i!=kmax; ++i)
    {
        //gsInfo <<"\n Skeleton:\n" << skeleton.front() <<"\n";
        for(index_t j = dm1; j!=-1; --j)
        {
            //multiplyKP( unf[j], skPtr[j]);
            gsKroneckerProduct(skPtr[j], tmp);
            // "single hole" tensor
            SHT.noalias() = unf[j] * tmp; // (nj x ni*nk) * (ni*nk x ri*rk)
            svd.compute(SHT, svdOpts);

            // Truncate upto tolerance
            const index_t new_rj = gsPickHead(svd.singularValues(), tol);
            if ( new_rj != rk[j] )
                improved = true;
            {   // Note: maybe we need to update all skeletons once one is improved?
                rk[j] = new_rj;
                skeleton[j] = svd.matrixU().leftCols(new_rj); // nj x new_rj
            }
        }

        if ( (i+1) % 8 == 0 &&  !improved ) break;

        if ( improved ) gsInfo<<"Tucker Rank reduced to "<<rk.transpose()<<"\n";
    }

    if ( improved )
        gsInfo << "\nTucker: Still improving after "<<kmax
               <<" itererations, try more iterations\n";

    // Compute the core tensor
    if ( 0==kmax ) // If no ALS iterations at all
    {
        gsKroneckerProduct(skPtr.front(), tmp);
        core.noalias() = skeleton.front().transpose() * unf.front() * tmp;
    }
    else // SHT already computed
        core.noalias() = skeleton.front().transpose() * SHT;

    core.resize(core.size(), 1); for(index_t k = 0; k!=d; ++k) psz[k] = rk[perm[k]];
    permuteTensorVector(perm, psz, core);// note: can we avoid permute here (ie. perm=id) ?
    return true; // always true !
}
//*/


template<class T>
bool gsHoSvd(const gsMatrix<T> & tensor,
                 const gsVector<int> & sz,
                 std::vector<gsMatrix<T> > & skeleton,
                 gsMatrix<T> & core,
                 T tol = 1e-11)
{
    return gsTuckerDec(tensor, sz, skeleton, core, tol, 0);// run Tucker without ALS
    // Note: one optimization in the code below is only 1 unfolding is stored

    GISMO_ASSERT( tensor.cols()==1, "Expecting a single column" );
    GISMO_ASSERT( sz.prod()==tensor.size(), "Invalid tensor/sizes" );

    const index_t d  = sz.size();

    if ( tensor.squaredNorm() <= tol )
    {
        // Zero tensor as rank-1, for compatibility
        skeleton.resize(d);
        for(index_t i = 0; i!=d; ++i)
            skeleton[i].setZero(sz[i],1);
        core.setZero(1, 1);
        return true;
    }

    const index_t dm1 = d-1;

    gsEigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    const unsigned svdOpts = gsEigen::ComputeThinU;

    gsVector<index_t> rk(d); // Tucker rank
    gsVector<index_t> perm(d), psz(d);
    gsMatrix<T> tmp, unf;

    skeleton.resize(d);
    for(index_t i = dm1; i!=-1; --i)// reversed to keep the perm we need
    {
        // Compute permuted unfoldings, starting from 0,1,..,d-1
        unf = tensor;
        psz = sz;
        // Permutation: shift to k, reverse tail
        perm = gsVector<index_t>::LinSpaced(d,0,dm1);
        std::rotate (perm.data()    , perm.data() + i, perm.data() + d);
        std::reverse(perm.data() + 1, perm.data() + d );
        permuteTensorVector(perm, psz, unf);
        unf.resize(psz[0], tensor.size() / psz[0]);

        // Compute SVD and truncate upto tolerance
        svd.compute(unf, svdOpts);
        rk[i] = gsPickHead(svd.singularValues(), tol);
        skeleton[i] = svd.matrixU().leftCols(rk[i]);
    }

    // Compute the core tensor
    std::vector<gsMatrix<T>*> skPtr; skPtr.reserve(dm1);
    for(index_t j = 1; j!=d; ++j)
        skPtr.push_back( &skeleton[perm[j]] );
    gsKroneckerProduct(skPtr, tmp);
    core.noalias() = skeleton.front().transpose() * unf * tmp; // todo: avoid tmp
    core.resize(core.size(), 1); for(index_t k = 0; k!=d; ++k) psz[k] = rk[perm[k]];
    permuteTensorVector(perm, psz, core);// note: Inv(perm)=perm
    return true;
}

// shorthand for: compute Tucker, convert fibers to CP.
// manybe to be removed
template<class T>
bool gsCpFromTucker(const gsMatrix<T> & tensor,
                    const gsVector<int> & sz,
                    std::vector<gsMatrix<T> > & skeleton,
                    T tol = 1e-11,
                    const index_t kmax = 10)
{

    gsMatrix<T>  tucCore;
    std::vector<gsMatrix<T> > tucSk;
    // Compute Tucker tensor
    const bool converged = gsTuckerDec(tensor, sz, tucSk, tucCore, tol, kmax);

    // Convert to CP
    gsTuckerFibers2Cp(tucCore, tucSk, skeleton, tol );

    return converged;
}

template<class T>
void gsTucker2Full(const gsVector<index_t> & rk,// rank ---> remove
                   const std::vector<gsMatrix<T> > & skeleton,
                   const gsMatrix<T> & core,
                   gsMatrix<T> & result)
{
    const index_t d = rk.size();
    //const gsVector<index_t> rk = tuckerRank(skeleton); // to do ?
    gsVector<index_t> r(d), rstr(d), s(d), sstr(d), sz(d);

    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz,sstr);
    tensorStrides(rk,rstr);

    result.resize(sz.prod(),1);

    s.setZero();
    do
    {
        T & cur = result.at( s.dot(sstr) );
        cur = 0.0;
        r.setZero();
        do
        {
            T tmp = core( r.dot(rstr) );
            for (index_t i = 0; i!=d; ++i)
                tmp *= skeleton[i](s[i],r[i]);
            cur += tmp;

        } while (nextLexicographic(r, rk));
    } while (nextLexicographic(s, sz));

}

// a matrix-based computation of the full tensor
// This is less efficient, probably, but demonstrates n-mode multiplication
template<class T>
void gsTucker2Fullv2(const std::vector<gsMatrix<T> > & skeleton,
                   const gsMatrix<T> & core,
                   gsMatrix<T> & result)
{
    const gsVector<index_t> rk = tuckerRank(skeleton);
    const index_t d = rk.size();
    // const index_t csz = rk.prod();
    gsVector<index_t> perm(rk), sz(d), psz(d);
    for (index_t i = 0; i!=d; ++i) // tuckerdims
        sz[i] = skeleton[i].rows();

    gsMatrix<T> tmp;
    result = core;
    for (index_t i = 0; i!=d; ++i)
    {
        tmp.swap(result);
        // Get unfolding matching the skeleton[i]
        psz = rk;
        for (index_t j = 0; j!=i; ++j) psz[j] = sz[j];
        perm = gsVector<index_t>::LinSpaced(d,0,d-1);
        std::rotate (perm.data()    , perm.data() + i, perm.data() + d);
        std::reverse(perm.data() + 1, perm.data() + d );
        permuteTensorVector(perm, psz, tmp);
        tmp.resize(psz[0], psz.prod() / psz[0]);

        // Multiply by skeleton[i]
        result.noalias() = skeleton[i] * tmp;
        result.resize(result.size(),1);

        //bring result back to original state
        psz[0] = sz[i];
        permuteTensorVector(perm, psz, result);
    }
}

// a Kronecker product-based computation of the full tensor
// This is less efficient, probably, but it demonstrates tensor algebra
template<class T>
void gsTucker2Fullv3(const std::vector<gsMatrix<T> > & skeleton,
                   const gsMatrix<T> & core,
                   gsMatrix<T> & result)
{
    const gsVector<index_t> rk = tuckerRank(skeleton);
    const index_t d = rk.size();
    const index_t csz = rk.prod();
    gsVector<index_t> perm(rk), sz(d), psz(d);
    for (index_t i = 0; i!=d; ++i) // tuckerdims
        sz[i] = skeleton[i].rows();

    gsAsConstMatrix<T> core1(core.data(), rk[0], csz / rk[0] );
    std::vector<gsMatrix<T>*> sk; sk.reserve(d-1);
    //for(index_t j = d-1; j!=0; --j) // gsKroneckerProductStandard
    for(index_t j = 1; j!=d; ++j) // gsKroneckerProduct
        sk.push_back( const_cast<gsMatrix<T>*>(&skeleton[j]) );
    gsMatrix<T> tmp;
    gsKroneckerProduct(sk, tmp);
    // (n1 x r1) * (r1 x r2*r3) * (r2*r3 x n2*n3)
    result.noalias() = skeleton.front() * core1 * tmp.transpose();
    result.resize(result.size(),1);
}




// ----------------- Canonical

// gsMatrix::pseudoinverse()

template<class T> inline
void gsPseudoInv(const gsMatrix<T> & mat, gsMatrix<T> & result)
{
    const index_t n = mat.rows();
    result = mat.colPivHouseholderQr().solve(gsMatrix<T>::Identity(n,n));

    /*
    //Make the SVD compute U and V
    const unsigned svdOpts = Eigen::ComputeThinU | Eigen::ComputeThinV;

    Eigen::JacobiSVD< typename gsMatrix<T>::Base > svd(mat, svdOpts);

    result.noalias() = svd.matrixV() * ( 1 / svd.singularValues().array() ).matrix() * svd.matrixU().transpose();
    */
}

template<class T>
T gsCanoDiffNormSqr(const gsMatrix<T> & tensor,
                    const std::vector<gsMatrix<T> > & skeleton)
{
    GISMO_ASSERT(tensor.cols()==1,"Expecting single tensor");
    gsVector<index_t> rnk = gsVector<index_t>::Constant(1,skeleton.front().cols());
    return gsCanoDiffNormSqr(tensor, rnk, skeleton).value() ;
}

template<class T>
gsVector<T> gsCanoDiffNormSqr(const gsMatrix<T> & tensor,
                              const gsVector<index_t> & ranks,//one for each column of \a tensor
                              const std::vector<gsMatrix<T> > & skeleton)
{
    const index_t d = skeleton.size();
    const index_t n = ranks.rows();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz,sstr);

    gsVector<T> result(n), cur(n);
    result.setZero();
    s.setZero();
    do
    {
        cur.setZero();

        index_t pos = 0;
        for(index_t k = 0; k!=n; ++k)
        {
            for(index_t r = 0; r!=ranks[k]; ++r)
            {
                T tmp = 1.0;
                for (index_t i = 0; i!=d; ++i)
                    tmp *= skeleton[i](s[i], pos+r);

                cur[k] += tmp;
            }
            pos += ranks[k];
        }

        cur -= tensor.row(s.dot(sstr)).transpose();
        result.array() += cur.array() * cur.array();

    } while (nextLexicographic(s, sz));

    return result;
}

template<class T>
T gsMaxDiff(const gsMatrix<T> & tensor,
            const std::vector<gsMatrix<T> > & skeleton)
{
    GISMO_ASSERT(tensor.cols()==1,"Expecting single tensor");

    const index_t d = skeleton.size();
    const index_t rk = skeleton.front().cols();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz,sstr);

    T result = 0.0;
    s.setZero();
    do
    {
        T cur = 0.0;

        for(index_t r = 0; r!=rk; ++r)
        {
            T tmp = 1.0;
            for (index_t i = 0; i!=d; ++i)
                tmp *= skeleton[i](s[i],r);
            cur += tmp;
        }

        const index_t ss = s.dot(sstr);
        result = math::max(result,
                           math::abs(tensor.at(ss)-cur) );

    } while (nextLexicographic(s, sz));

    return result;
}

template<class T>
void gsCpResidue(const gsMatrix<T> & tensor,
                 const std::vector<gsMatrix<T> > & skeleton,
                 gsMatrix<T> & result)
{
    GISMO_ASSERT(tensor.cols()==1,"Expecting single tensor");

    const index_t d  = skeleton.size();
    const index_t rk = skeleton.front().cols();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz, sstr);
    result.resizeLike(tensor);

    s.setZero();
    do
    {
        T cur = 0.0;

        for(index_t r = 0; r!=rk; ++r)
        {
            T tmp = skeleton[0](s[0],r);
            for (index_t i = 1; i!=d; ++i)
                tmp *= skeleton[i](s[i],r);
            cur += tmp;
        }

        const index_t ss = s.dot(sstr);
        result.at(ss) = tensor.at(ss) - cur;

    } while (nextLexicographic(s, sz));
}

template<class T>
gsVector<T> gsMaxDiff(const gsMatrix<T> & tensor,
                      const gsVector<index_t> & ranks,
                      const std::vector<gsMatrix<T> > & skeleton)
{
    GISMO_ASSERT(ranks.sum()==skeleton.front().cols(), "error");

    const index_t d = skeleton.size();
    const index_t n = ranks.rows();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz,sstr);

    gsVector<T> result(n), cur(n);
    result.setZero();
    s.setZero();
    do
    {
        cur.setZero();

        index_t pos = 0;
        for(index_t k = 0; k!=n; ++k)
        {
            for(index_t r = 0; r!=ranks[k]; ++r)
            {
                T tmp = 1.0;
                for (index_t i = 0; i!=d; ++i)
                    tmp *= skeleton[i](s[i], pos+r);

                cur[k] += tmp;
            }
            pos += ranks[k];
        }

        cur = (tensor.row(s.dot(sstr)).transpose()-cur).array().abs();
        result = result.cwiseMax( cur );

    } while (nextLexicographic(s, sz));

    return result;
}

// ALS-CP with fixed rank
// TODO: Common - ELS line search
template<class T, int D>
bool gsCpAls(const gsMatrix<T> & tensor,
             const gsVector<index_t,D> & sz,
             const int rk,// rank
             std::vector<gsMatrix<T> > & skeleton,
             T tol = 1e-11,
             const index_t kmax = 10)
{
    GISMO_ASSERT( rk>0, "Rank needs to be positive");
    const index_t d  = sz.size();
    const index_t pr = sz.prod();
    skeleton.resize(d);

    // 1. Store the unfoldings of the input tensor once and for all
    std::vector<gsMatrix<T> > unf(d, tensor);
    std::vector<std::vector<gsMatrix<T>*> > skPtr(d);

    gsVector<index_t,D> perm(d), psz(d);
    for(index_t i = 0; i!=d; ++i)
    {
        // Compute unfoldings
        psz = sz;
        perm = gsVector<index_t,D>::LinSpaced(d,0,d-1);
        std::rotate (perm.data()  , perm.data()+i, perm.data()+d);
        std::reverse(perm.data()+1, perm.data()+d );
        permuteTensorVector(perm, psz, unf[i]);
        unf[i].resize(psz[0], pr / psz[0]);

        // Keep permutations of skeletons for Khatri-Rao products
        skPtr[i].reserve(d-1);
        for(index_t j = 1; j!=d; ++j)
            skPtr[i].push_back( &skeleton[perm[j]] );
    }

    // *** Phase I - Initial Guess
    if ( !skeleton.front().cols() ) // empty matrix ?
    {
        // No initial point given -- compute default initial solution
        gsEigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
        const unsigned svdOpts = gsEigen::ComputeThinU;//Note: ThinU produces error if r is big
        for(index_t i = 0; i!=d; ++i)
        {
            skeleton[i].resize(sz[i], rk);
            const index_t n = math::min(rk,sz[i]);
            svd.compute(unf[i], svdOpts);
            skeleton[i].leftCols(n) = svd.matrixU().leftCols(n);
            if ( n<rk) // Initial solution of more than sz[i] skeleton vectors ??
            {
                //gsWarn<<"Rank bigger than the mode size.\n";
                //skeleton[i].rightCols(rk-n).setConstant(1e-1);
                skeleton[i].rightCols(rk-n).setConstant(1); // TODO: initialize by residual rank-1
            }
        }
    }
    else
    {
        //gsInfo<<"CP: Using given initial.\n";
        for(index_t i = 0; i!=d; ++i)
            if (skeleton[i].cols() < rk)
            {
                GISMO_ASSERT(skeleton[i].rows() == sz[i], "Invalid initial solution.");
                skeleton[i].conservativeResize(gsEigen::NoChange, rk);
                // Note: Numerical problem here!
                skeleton[i].rightCols(rk-math::min(rk,sz[i])).setConstant(1e-1);
            }
            // else: the initial point is already fine
    }

    T resNorm = gsCanoDiffNormSqr(tensor, skeleton);

    // *** Phase II - ALS Iteration
    T delta;
    std::vector<gsVector<T> > skNorm(d);
    gsMatrix<T> tmp, tmpInv;

    for(index_t i = 0; i!=kmax; ++i)
    {
        for(index_t j = 0; j!=d; ++j)
        {
            //gsKhatriRaoProduct(skPtr[j], tmp);
            //tmp.transposeInPlace();//
            gsKhatriRaoProductTr<T,D>(skPtr[j], tmp);
            gsPseudoInv(tmp, tmpInv);// to do: inline this
            //const index_t n = tmp.rows();// Eigen bug ? col-major..
            //tmpInv = tmp.transpose().colPivHouseholderQr().solve(gsMatrix<T>::Identity(n,n));

            skeleton[j].noalias() = unf[j] * tmpInv;
            skNorm[j] = skeleton[j].colwise().norm();

            // POSSIBILITY: equilibrate inside this micro-step ?

            //GISMO_ASSERT( skNorm[j].array().all() > 1e-12, "CP: Norm=0" );
        }

        // Equilibrate
        for(index_t j = 0; j!=rk; ++j)
        {
            delta = skNorm[0][j];
            for(index_t k = 1; k!=d; ++k)
                delta *= skNorm[k][j];
            delta = math::pow(delta, 1.0 / d);

            if (delta > 1e-11) // avoid division by zero (skNorm)
            {
                for(index_t k = 0; k!=d; ++k)
                    skeleton[k].col(j) *= delta / skNorm[k][j];
            }
        }
        // to do
        //for(index_t k = 0; k!=d; ++k)
        //    skeleton[k].colwise().cwiseProduct( delta / skNorm[k].array() );

        const T curResNorm = gsCanoDiffNormSqr(tensor, skeleton);

        // Stopping if the residual norm does not decrease
        if ( math::abs( resNorm - curResNorm ) < tol )  // tol2
        {
            //gsDebug<<"CP: ALS converged at iteration "<<i<<", tol: "<<tol<<".\n";
            return true;
        }

        resNorm = curResNorm;
    }
    // /*
      gsDebug << "CP: No convergence for tol="<<tol<<", iter="
      <<kmax<<", r="<<rk<<"). max-diff: "
      <<gsMaxDiff(tensor, skeleton)<<"\n";
    //*/
    return false;
}


template<class T>
gsVector<T> gsCanoDistL2Sqr(const gsMatrix<T> & tensor,
                            const gsVector<index_t> & ranks,//one for each column of \a tensor
                            const std::vector<gsMatrix<T> > & skeleton)
{
    const index_t d = skeleton.size();
    const index_t n = ranks.rows();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    // Mass matrix
    gsSparseMatrix<T> M[d]; //fun.basis().component(i)
    std::vector<gsMatrix<T> > & skM(d);
    for (index_t i = 0; i!=d; ++i)
        skM[i] = M[i] * skeleton[i]; // M * skeleton

    gsMatrix<T> tensorM; // M*tensor
    // gsLRSparseMatrix<T> Mass(M); Mass.apply(tensor, tensorM);

    tensorStrides(sz,sstr);

    gsVector<T> result(n), cur(n), curM(n);
    result.setZero();
    s.setZero();
    do
    {
        cur .setZero();
        curM.setZero();

        index_t pos = 0;
        for(index_t k = 0; k!=n; ++k)
        {
            for(index_t r = 0; r!=ranks[k]; ++r)
            {
                T tmp  = 1.0, tmpM = 1.0;
                for (index_t i = 0; i!=d; ++i)
                {
                    tmp  *= skeleton[i](s[i], pos+r);
                    tmpM *= skM     [i](s[i], pos+r);
                }

                cur [k] += tmp ;
                curM[k] += tmpM;
            }
            pos += ranks[k];
        }

        const index_t gl = s.dot(sstr);
        cur  -= tensor .row(gl).transpose(); // CP - Tensor
        curM -= tensorM.row(gl).transpose(); // M *CP - M * Tensor
        result.array() += cur.array() * curM.array();

    } while (nextLexicographic(s, sz));

    return result;
}

/*
template<class T>
bool gsL2SvdAls(const gsGeometry<T> & fun,
                const int rk,// rank
                std::vector<gsMatrix<T> > & skeleton,
                T tol = 1e-11,
                const index_t kmax = 10)
{
    GISMO_ASSERT( rk>0, "Rank needs to be positive");
    const index_t d  = fun.parDim();
    const index_t pr = fun.basis().size();
    skeleton.resize(d);

    // Mass matrix
    gsSparseMatrix<T> Mass; //fun.basis().


}
//*/

// ALS-CP Deflation algorithm
// Comon et. al.
template<class T, int D>
bool gsCpAlsDeflation(const gsMatrix<T> & tensor,
                      const gsVector<index_t,D> & sz,
                      const int rk,// rank
                      std::vector<gsMatrix<T> > & skeleton,
                      T tol = 1e-11,
                      const index_t kmax = 10)
{
    GISMO_ASSERT( rk>0, "Rank needs to be positive");
    const index_t d  = sz.size();

    std::vector<std::vector<gsMatrix<T> > > tmpSk(rk);
    gsMatrix<T> YY(tensor);
    // Initialize
    for(index_t r = 0; r!=rk; ++r)
    {
        // Best Rank-1 approximation (try Tucker as well!)
        gsCpAls(YY, sz, 1, tmpSk[r], tol, kmax);
        // Substract from tensor
        gsCpResidue(YY, tmpSk[r], YY);
    }

    // Iterate
    gsMatrix<T> residue(YY);
    T rNorm1, rNorm0 = YY.squaredNorm();
    index_t i;
    for(i = 0; i!=kmax; ++i)
    {
        for(index_t r = 0; r!=rk; ++r)
        {
            // Y = X_r + R
            gsCp2Full(tmpSk[r], YY);
            YY += residue;

            // X_r = bestrankone(Y)
            gsCpAls(YY, sz, 1, tmpSk[r], tol, kmax);

            // R = Y - X_r
            gsCpResidue(YY, tmpSk[r], residue);
        }

        // Stopping criterion
        rNorm1 = residue.squaredNorm();
        if ( math::abs(rNorm0-rNorm1) <= tol )
            break;
        rNorm0 = rNorm1;
    }

    //gsInfo <<"-->"<<i<<", Error: "<< rNorm1 <<"\n";

    skeleton.resize(d);
    for(index_t j = 0; j!=d; ++j)
    {
        skeleton[j].resize(sz[j], rk);
        for(index_t r = 0; r!=rk; ++r)
            skeleton[j].col(r) = tmpSk[r][j];
    }

    return i!=kmax;
}

// ALS-CP Deflation algorithm
// Comon et. al.
template<class T, int D>
bool gsCpAlsGreedy(const gsMatrix<T> & tensor,
                   const gsVector<index_t,D> & sz,
                   std::vector<gsMatrix<T> > & skeleton,
                   T tol = 1e-11,
                   const index_t kmax = 20,
                   const index_t rmax = 80)
{
    const index_t d  = sz.size();

    std::vector<gsMatrix<T> > tmpSk;
    gsMatrix<T> residue(tensor);
    skeleton.resize(d);
    skeleton.clear();

    T rNorm1, rNorm0 = residue.squaredNorm();
    index_t r;
    for(r = 0; r!=rmax; ++r)
    {
        // Best Rank-1 approximation (try Tucker as well!)
        gsCpAls(residue, sz, 1, tmpSk, tol, kmax);

        // update initial solution
        for ( index_t k = 0; k!=d; ++k)
        {
            // note: re-allocation might be inefficient here
            skeleton[k].conservativeResize(sz[k], r+1 );
            skeleton[k].col(r) = tmpSk[k];
        }

        // Update residue
        gsCpResidue(residue, tmpSk, residue);

        // Stopping criterion
        rNorm1 = residue.squaredNorm();
        if ( math::abs(rNorm0-rNorm1) <= tol )
            break;
        rNorm0 = rNorm1;
    }
    return r!=rmax;
}


// Same as overload, but Without output of the canonical ranks
template<int method, class T>
void gsCpDec(const gsMatrix<T> & tensor,
             const gsVector<int> & sz,
             std::vector<gsMatrix<T> > & skeleton,
             T tol = 1e-11,
             const index_t kmax = 10)
{
    gsVector<index_t> rank;
    gsCpDec(tensor, sz, skeleton, rank, tol, kmax);
}

//Sequential Rank-One Approximation
//template<int method, class T>
template<class T>
bool gsBestRankOne(const gsMatrix<T> & tensor,
                   const gsVector<int> & sz,
                   std::vector<gsMatrix<T> > & skeleton,
                   T tol = 1e-11,
                   const index_t kmax = 20)
{
    skeleton.clear();// gsCpAls accepts initial solution
    //return gsCpAlsCore(tensor, sz, skeleton, tol, kmax, 1); // why no good quality ?
    return gsCpAls(tensor, sz, 1, skeleton, tol, kmax);

    /* //UC: sub-optimal, from Common et. al.
    const index_t d  = sz.size();
    const index_t pr = sz.prod();
    skeleton.resize(d);
    Eigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    gsMatrix<T> unf(tensor);

    unf.resize(psz[0], pr / sz[0]);
    svd.compute(unf, Eigen::ComputeThinU|Eigen::ComputeThinV);
    skeleton[0] = svd.matrixU().leftCols(1);
    unf = svd.matrixV();

    for(index_t i = 1; i+1<d; ++i)
    {
        Eigen::ComputeThinV;
    */
}

/**

Input: a skeleton of certain rank
Output:  A new skeleton of rank one more
 */
template<class T>
void gsCpRankIncrement(const gsMatrix<T> & tensor,
                       const gsVector<int> & sz,
                       std::vector<gsMatrix<T> > & skeleton,
                       T tol = 1e-11,
                       const index_t kmax = 10) // maximum number of iterations
{
    GISMO_ASSERT( tensor.cols()==1, "Expecting a single column" );

    const index_t d  = sz.size();
    skeleton.resize(d); // in case it is empty
    const index_t rk = skeleton.front().cols();
    if (rk == 0 )
    {
        // Just compute rank-1 approximation and return
        gsCpAls(tensor, sz, 1, skeleton, tol, kmax);
        return;
    }

    gsMatrix<T> residue;
    gsCpResidue(tensor, skeleton, residue);

    // Get rank-1 approximation of the residue
    std::vector<gsMatrix<T> > skeleton1;
    gsCpAls(residue, sz, 1, skeleton1, tol, kmax);

    // update input skeleton
    for ( index_t k = 0; k!=d; ++k)
    {
        // note: re-allocation might be inefficient here
        skeleton[k].conservativeResize(sz[k], rk+1);
        skeleton[k].col(rk) = skeleton1[k];
    }

    // Refine skeleton by ALS
    gsCpAls(tensor, sz, rk+1, skeleton, tol, kmax);
}

// Successive rank-1 plus ALS correction
template<class T, int D>
bool gsCpAlsIncremental(const gsMatrix<T> & tensor,
                        const gsVector<index_t,D> & sz,
                        std::vector<gsMatrix<T> > & skeleton,
                        T tol = 1e-11,
                        const index_t kmax = 20,
                        const index_t rmax = 50)
{
    GISMO_ASSERT( tensor.cols()==1, "Expecting a single column" );

    const index_t d = sz.size();
    skeleton.clear();
    skeleton.resize(d);

    std::vector<gsMatrix<T> > skeleton1;
    skeleton1.reserve(d);
    gsMatrix<T> tmp, residue(tensor);

    for (index_t i = 0; i<rmax; ++i)
    {
        // best rank-1 of the residue
        skeleton1.clear();
        gsCpAls(residue, sz, 1, skeleton1, tol, kmax);

        // update initial solution
        for ( index_t k = 0; k!=d; ++k)
        {
            // note: re-allocation might be inefficient here
            skeleton[k].conservativeResize(sz[k], i+1);
            skeleton[k].col(i) = skeleton1[k];
        }

        gsCpAls(tensor, sz, i+1, skeleton, tol, kmax);
        //gsCpAlsDeflation(tensor, sz, i+1, skeleton, tol, kmax);

        gsCpResidue(tensor, skeleton, residue); // todo: residue comp.

        // to do: inside CpDecFixed ?

        const T rn = residue.squaredNorm();
        if ( rn <= tol ) // tol2
        //if ( residue.norm() < tol )
        //if ( residue.cwiseAbs().maxCoeff() < tol )
        {
            //gsInfo<<"Converged after iter: "<<i <<": res="<<residue.cwiseAbs().maxCoeff()<<"\n";
            return true;
        }
        else
            gsDebug <<"Current rank is "<< i +1 <<" and residual norm is "<< std::scientific << rn << " " << math::sqrt(rn) <<"\n";
    }

    return false;
}

// CP performed on the Tucker core of the HOSVD
template<class T>
bool gsCpAlsCore(const gsMatrix<T> & tensor,
                 const gsVector<int> & sz,
                 std::vector<gsMatrix<T> > & skeleton,
                 T tol = 1e-11,
                 const index_t kmax = 20,
                 const index_t rmax = 50)
{
    gsMatrix<T>  tucCore;
    std::vector<gsMatrix<T> > tucSk;
    // Compute HOSVD
    bool converged = gsHoSvd(tensor, sz, tucSk, tucCore, tol); // ERR: truncate all unfoldings (rel. spectral norm)

    gsDebug<<"HoSVD OK !\n";

    if (!converged)
        gsWarn<<"Problem in HOSVD!\n";

    // todo: get maximal canonical rank ( r1*..rd / max(ri)

    // Run CP on the core tensor
    std::vector<gsMatrix<T> > skCore;
    converged = gsCpAlsIncremental<T>(tucCore, tuckerRank(tucSk), skCore, tol, kmax, rmax);
    //gsDebugVar(tuckerRank(tucSk).transpose());

    // Compute skeleton
    const index_t d = skCore.size();
    skeleton.resize(d);
    for(index_t k = 0; k!=d; ++k)
        skeleton[k].noalias() = tucSk[k] * skCore[k];

    return converged;
}

/**
   Multiple tensors of the same structure

   method 0 : ALS with incremental rank-1 steps
   method 1 : converted from Tucker
   method 2 : ALS truncated at a specified error threshold
 */
template<int method, class T>
bool gsCpDec(const gsMatrix<T> & tensor,
             const gsVector<index_t> & sz,
             std::vector<gsMatrix<T> > & skeleton,
             gsVector<index_t> & ranks, // output the canonical ranks
             T tol = 1e-11,
             const index_t kmax = 10,
             const index_t rmax = 30) // maximum rank (per component)
{
    const index_t d = sz.size();
    ranks.resize(tensor.cols());
    std::vector<gsMatrix<T> > curSk;
    bool converged = true;

    //The squared norm is used inside the functions, so
    const T tol2 = tol*tol;
    //tol *= tol;

    // TO DO: gsHoSvd here, eg. with bool option, see gsCpAlsCore

    // 1. First colunm
    if ( tensor.col(0).squaredNorm() <= tol2 )
    {
        ranks[0] = 0;
        converged = true;
    }
    else
    {
        switch (method)
        {
        case 0:
            converged = gsCpAlsIncremental<T>(tensor.col(0), sz, skeleton, tol, kmax, rmax);
            break;
        case 1:
            converged = gsCpFromTucker<T>(tensor.col(0), sz, skeleton, tol, kmax);
            break;
        case 2:
            converged = gsCpAlsGreedy<T>(tensor.col(0), sz, skeleton, tol, kmax, rmax);
            break;
        case 3:
            converged = gsCpAlsCore<T>(tensor.col(0), sz, skeleton, tol2, kmax, rmax);
            break;
        default:
            GISMO_ERROR("No such method");
        }

        ranks[0] = skeleton.front().cols();
    }
    // if (!converged)// should have converged by now
    //     gsDebug << "CP("<<0<<"): Could not satisfy tolerance "<<tol
    //             <<" for rank upto "<<kmax<<"\n";

    // 2. For the rest of the columns
    index_t s = ranks[0];
    for ( index_t c = 1; c < tensor.cols(); ++c) // for all columns
    {
        if ( tensor.col(c).array().abs().maxCoeff() <= tol )
        //if ( tensor.col(c).squaredNorm() < tol )
        {
            ranks[c] = 0;
            converged = true;
            continue;
        }

        curSk.clear();//-> probably to be removed

        switch (method)
        {
        case 0:
            converged = gsCpAlsIncremental<T>(tensor.col(c), sz, curSk, tol, kmax, rmax);
            break;
        case 1:
            converged = gsCpFromTucker<T>(tensor.col(c), sz, curSk, tol, kmax);
            break;
        // case 2:
        //     converged = gsCpAlsDeflation<T>(tensor.col(c), sz, curSk, tol, kmax);
            break;
        case 3:
            converged = gsCpAlsCore<T>(tensor.col(c), sz, curSk, tol2, kmax, rmax);
            break;
        default:
            GISMO_ERROR("No such method");
        }

        ranks[c] = curSk.front().cols();
        s+=ranks[c];
         if (!converged)
             gsDebug << "CP("<<c<<")\n";

        // Add the current skeleton vectors
        for ( index_t k = 0; k!=d; ++k)
        {
            //GISMO_ASSERT( is_finite(curSk[k]), "Has infs: "<< k );
            skeleton[k].conservativeResize(gsEigen::NoChange, s);
            skeleton[k].rightCols(ranks[c]) = curSk[k];
        }

    } // end for all columns

    //gsDebugVar(ranks.transpose() );
    //gsDebugVar(ranks.sum()       );
    //gsDebugVar(skeleton.front().cols());

    GISMO_ASSERT(skeleton.front().cols() == ranks.sum(), "error");
    return converged;
}
//*/

template<class T>
void gsCp2Full(const std::vector<gsMatrix<T> > & skeleton,
                 gsVector<T> & coeffs,// coefficients
                 gsMatrix<T> & result)
{
    const index_t d = skeleton.size();
    const index_t rk = skeleton.front().cols();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz,sstr);

    result.resize(sz.prod(),1);

    s.setZero();
    do
    {
        T & cur = result.at(s.dot(sstr) );
        cur = 0.0;

        for(index_t r = 0; r!=rk; ++r)
        {
            // to do: more efficient
            T tmp = coeffs[r];
            for (index_t i = 0; i!=d; ++i)
                tmp *= skeleton[i](s[i],r);
            cur += tmp;
        }

    } while (nextLexicographic(s, sz));

}

template<class T>
void gsCp2Full(const std::vector<gsMatrix<T> > & skeleton,
                 gsMatrix<T> & result)
{
    const index_t d = skeleton.size();
    const index_t rk = skeleton.front().cols();
    gsVector<index_t> s(d), sstr(d), sz(d);
    for (index_t i = 0; i!=d; ++i)
        sz[i] = skeleton[i].rows();

    tensorStrides(sz,sstr);

    result.resize(sz.prod(),1);

    s.setZero();
    do
    {
        T & cur = result.at(s.dot(sstr) );
        cur = 0.0;

        for(index_t r = 0; r!=rk; ++r)
        {
            // to do: more efficient
            T tmp = 1.0;
            for (index_t i = 0; i!=d; ++i)
                tmp *= skeleton[i](s[i],r);
            cur += tmp;
        }

    } while (nextLexicographic(s, sz));

}


// Using Khatri-Rhao product
template<class T>
void gsCp2Fullv2(const std::vector<gsMatrix<T> > & skeleton,
                 gsMatrix<T> & result)
{

}


template<class T>
void gsTuckerFibers2Cp(const gsMatrix<T> & core,
                       const std::vector<gsMatrix<T> > & tuc,
                       std::vector<gsMatrix<T> > & can,
                       T tol = 1e-11)
{
    GISMO_ASSERT( tuc.size()==3, "Implemented for 3D only" );
    GISMO_ASSERT( core.cols()==1, "Expecting a single column" );

    static const index_t d = tuc.size();
    gsVector<index_t> rk = tuckerRank(tuc);
    index_t pr = rk.prod();

    // TODO: permute core+skeleton s.t. r_k is sorted in increasing order

    // const T cnorm = core.norm();
    // tol *= cnorm;

    // Unfold all but the 2 largest modes
    const index_t u = rk[d-2], v = rk[d-1];
    const gsAsConstMatrix<T> core0(core.data(), pr/(u*v), u*v);

    std::vector<gsMatrix<T> > CC(d);
    // todo: check if conservativeResize is better than pre-allocation
    const index_t bound = pr / rk[d-1];// upper bound (ideally rk[d-1] = max{rk})
    for (index_t i = 0; i!= d; ++i)
        CC[i].resize(rk[i], bound);
    // coefficients of the canonical decomposition
    gsVector<T> coeffs(bound);

    // SVD and options
    gsEigen::JacobiSVD< typename gsMatrix<T>::Base > svd;
    const unsigned svdOpts = gsEigen::ComputeThinU | gsEigen::ComputeThinV;

    gsMatrix<T> A;
    index_t count = 0;

    for(index_t r = 0; r!= core0.rows(); ++r) // for all fibers
    {
        A = core0.row(r);// copy (row is not cont. in memory)
        A.resize(u, v); // r_1, r_2

        svd.compute(A, svdOpts);

        // Pick most important singular values (at most v = r_{d-1})
        const index_t nsv = gsPickHead(svd.singularValues(), tol);
        coeffs.segment(count,nsv) = svd.singularValues().head(nsv);

        for (index_t i = 0; i+2 < d; ++i)
        {
            //3D  Note: this step is currently restricted to d==3, i==0
            CC[i].middleCols(count,nsv) = gsVector<T>::Unit(rk[i], r).replicate(1,nsv);
        }
        CC[d-2].middleCols(count,nsv) = svd.matrixU().leftCols(nsv);
        CC[d-1].middleCols(count,nsv) = svd.matrixV().leftCols(nsv);
        count += nsv;
    }

    // Compute the canonical vectors
    can.resize(d);
    for (index_t i = 0; i!=d; ++i)
        can[i].noalias() = tuc[i] * CC[i].leftCols(count);

    // Absorve coefficients in the first skeleton matrix
    can.front().noalias() = can.front() * coeffs.head(count).asDiagonal();
}


template<class T>
void gsSvdDecomposition(const gsMatrix<T> & tensor2, const gsVector<index_t> & sz,
                        gsVector<index_t> & crank, // Note: vector functions
                        std::vector< gsMatrix<T> > & skeleton,
                        T tol, const int maxrank = -1)
{
    const index_t d = sz.size();
    GISMO_ASSERT( d == 2 && sz.prod() == tensor2.rows(), "Invalid tensor structure");

    //gsDebugVar( sz.transpose() );

    skeleton.resize(d);

    crank.resize(tensor2.cols());

    gsMatrix<T> tmp;
    // Permutations ?
    gsVector<T> sqsv;

    gsEigen::JacobiSVD< typename gsMatrix<T>::Base > svd; // TODO gsAsConstMatrix<T>
    const unsigned svdOpts = gsEigen::ComputeThinU | gsEigen::ComputeThinV;

    //gsInfo<<"svd threshold: "<<svd.threshold() <<"\n";

    //gsDebugVar(tensor2.cols() );
    //gsDebugVar(tensor2.rows() );
    bool createdMatrix = false;
    for ( index_t i = 0; i<tensor2.cols(); ++i ) //Loop over the components(i.e. 1 for mass, 9 for stiffness)
    {
        //Check if the component is actually 0, then we don't have to do SVD
        if( tensor2.col(i).squaredNorm() <= tol*tol)
        {
            crank[i] = 0;
            continue;
        }

        gsAsConstMatrix<T> Amat( tensor2.col(i).data(), sz[0], sz[1]);
        svd.compute(Amat, svdOpts);
        //const index_t numSV = svd.singularValues().size();
        int rank = gsPickHead(svd.singularValues(), tol*tol); //can be 0 if the value is below the tolerance
        if(maxrank >= 0) rank = std::min(maxrank, rank); //If a maximal rank is prescribed, set it here

        //std::min(rmax,pickhead);
        //gsDebug << "Singular Values: \n" << std::scientific << svd.singularValues() << std::fixed << "\n";
        //gsDebug << "Rank: "<< rank <<" (full was "<<svd.singularValues().size()<<")\n";
        sqsv = svd.singularValues().topRows(rank).array().sqrt();
        // gsDebug<< "Singular values: "<< svd.singularValues().topRows(rank).transpose()
        //        <<" ("<< ( rank < numSV ? svd.singularValues()[rank] : 0)<<",..)\n";

        // special case of zero rank
        if ( 0 != rank )
        {
            // Append to matrix
            if (!createdMatrix) //If it is the first time that the rank is not zero.
            {
                skeleton[0] = svd.matrixU().leftCols(rank) * sqsv.asDiagonal();
                skeleton[1] = svd.matrixV().leftCols(rank) * sqsv.asDiagonal();
                createdMatrix = true;
            }
            else// Append to matrix
            {
                tmp.resize(skeleton[0].rows(), skeleton[0].cols() + rank);
                tmp << skeleton[0], // 2D
                        svd.matrixU().leftCols(rank) * sqsv.asDiagonal();
                tmp.swap(skeleton[0]);
                // skeleton[0].conservativeResize(gsEigen::NoChange, skeleton[0].cols()+rank);
                // skeleton[0].rightCols(rank)= svd.matrixU().leftCols(rank)*sqsv.asDiagonal();

                tmp.resize(skeleton[1].rows(), skeleton[1].cols() + rank);
                tmp << skeleton[1], // 2D
                        svd.matrixV().leftCols(rank) * sqsv.asDiagonal();
                tmp.swap(skeleton[1]);
                // skeleton[1].conservativeResize(gsEigen::NoChange, skeleton[1].cols()+rank);
                // skeleton[1].rightCols(rank)= svd.matrixV().leftCols(rank)*sqsv.asDiagonal();
            }
        }

        crank[i] = rank;

        //gsDebug<< "Tolerance:"<< tol <<"\n";
        //gsDebug<< "Low-rank Approx. error (rank "<<rank
        //       <<", out of "<<svd.singularValues().size()<<"): "
        //<< ( skeleton[0].rightCols(rank) * skeleton[1].rightCols(rank).transpose() - Amat ).norm() <<"\n";

        gsInfo << "Tolerance:"<< tol <<"\n";
        gsInfo << "Low-rank Approx. error (rank "<<rank
               <<", out of "<<svd.singularValues().size()<<"): "
        << ( skeleton[0].rightCols(rank) * skeleton[1].rightCols(rank).transpose() - Amat ).norm() <<"\n";
    }
}

template<class T>
void gsRandomTensor(const gsVector<index_t> & sz,
                    const gsVector<index_t> & rk,
                    std::vector<gsMatrix<T> > & skeleton,
                    gsMatrix<T> & core)
{
    const index_t d = rk.size();
    skeleton.resize(d);
    for (index_t i = 0; i!=d; ++i)
        skeleton[i].setRandom(sz[i], rk[i]);
    core.setRandom(rk.prod(), 1);
}

// Boris:
// 1. best rank 1. (a)
// 2. best rank 1 of residue. (b)
// 3. Best rank 2 of original, using (a) and (b) as initial solutions. (c)
// 4. Best rank 1 of residue. (d)
// 5. Best rank 3 of original, using (c) and (d) as initial solution
// etc.



/*
inline bool isEqualPosition
(
    const gsCompressedSparseFiber<2>::const_iterator & csf,
    const std::vector<gsSparseMatrix<>::iterator> & spm
)
{
    for(size_t i =0; i!=spm.size();++i)
        if ( spm[i].index() != csf.index(i) )
            return false;
    return true;
}

inline bool isLessPosition
(
    const gsCompressedSparseFiber<2>::const_iterator & csf,
    const std::vector<gsSparseMatrix<>::iterator> & spm
)
{
    size_t first1(0), first2(0), last1(d), last2(d);
    for (; first1 != last1 && first2 != last2; ++first1, ++first2)
    {
        if ( csf.index(first1) < spm[first2].index() )
            return true;
        if ( spm[first2].index() < csf.index(first1) )
            return false;
    }
    return first1 == last1 && first2 != last2;
}

template<short_t d>
bool findNextMatch(
    gsCompressedSparseFiber<d>::const_iterator & csf,
    gsCompressedSparseFiber<d>::const_iterator & csf_end,
    std::vector<gsSparseMatrix<>::iterator> & spm
    )
{
    while ( csf && spm.back() )
    {
        if ( isEqualPosition(csf,spm) )
            return true;
        if ( isLessPosition(csf,spm) )
        {
            gsVector<index_t,d> idx(d);
            for(short_t k = 0; k!=d;++k)
                idx[k] = spm[k].index();
            csf = std::lower_bound( csf, csf_end, idx);

            if ( !++csf ) break;
        }
        else // spm < csf
        {
            spm = lower_bound( spm, csf.tensorIndex() );
            //CSF.
        }
    }
    return false;
}

inline void gsductSum
(
    const std::vector<gsSparseMatrix<> > & tkl,
    const index_t n,
    const gsCompressedSparseFiber<2> & rk,
    const gsCompressedSparseFiber<2> & rl,
    const index_t r_init,
    const index_t c_init,
    gsSparseEntries<> & result
)
{
    // Type definition
    using spmat_it = gsSparseMatrix<>::iterator;
    using csf_it   = gsCompressedSparseFiber<2>::const_iterator;

    const auto dd = static_cast<short_t>(tkl.size());
    const bool symmetric =  r_init!=c_init;

    std::vector<spmat_it> tkl_it(dd);
    gsVector<index_t,2> t_low, t_upp;
    csf_it rk_it, rk_end;
    real_t value;
    for (csf_it rl_it = rl.begin(); rl_it; ++rl_it)
    {
        for (short_t d=0; d!=dd; ++d)
        {
            tkl_it[d] = tkl[d].begin(rl_it.index(k));
            t_low [d] = tkl_it[d].index();
            t_upp [d] = (tkl.end(rl_it.index(k))-1).index();
        }

        rk_it  = rk.firstInRange(t_low, t_upp);
        while (rk_it)
        {


            rk_it.nextInRange(t_low, t_upp);
        }


        //idea: make CSF from tkl_it ? !! TKL has RANK-
        // go level by level in the tree and store tkl values.

        //todo: after the columnr rl.begin(), slide the pointers rk_it, rk_end
        // (restrict range for fast search)

        //add/iterate products of tkl.. inside the leafs of rk ? (++ with tkl following)
        rk_it  = rk.lower_bound(t_low);
        rk_end = rk.upper_bound(t_upp);
        //    findNextMatchTP(rk_it,r_end,..): find tkl_it[k] inside the CSF, advance together...

        // get arrays from tkl: &tkl[k].value()
        for(;rk_it!=rk_end;++rk_it) // WRONG: we need next in range t_low,t_upp
        {
            // for sum-fact: identify branches in the tree whilst ++rk_it

            // int i = rk_it.index(k) - tkl_it[k].index();
            // value at position i: tkl_it[k][i]

            value = 1.0;
            for (short_t d=0; d!=dd; ++d)
            {
                const index_t w = rk_it.index(k) - tkl_it[k].index();
                value *= tkl_it[k][w];
            }

            result.add(r_init + rk_it.back().row(), c_init + c, value);
            if (symmetric)
                result.add(c_init + c, r_init + rk_it.back().row(), value);
        }
    }
}
*/

} // namespace gismo
