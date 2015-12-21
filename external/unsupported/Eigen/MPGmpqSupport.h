// GMP support for Eigen linear algebra library
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Angelos Mantzaflaris, 2015

#pragma once

#include <gmpxx.h>

namespace Eigen 
{  	
  template<> struct NumTraits<mpq_class>
    : GenericNumTraits<mpq_class>
  {
    enum {
      IsInteger = 0,
      IsSigned  = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 10,
      AddCost  = 10,
      MulCost  = 40
    };

    typedef mpq_class Real;
    typedef mpq_class NonInteger;
    
    inline static Real epsilon  ()              { return 0; }
    inline static Real epsilon  (const Real& x) { return 0; }
    inline static Real dummy_precision()        { return 0; }
  };

  namespace internal {

  template<> inline mpq_class random<mpq_class>()
  {
    gmp_randclass rr(gmp_randinit_default);
    return mpq_class(rr.get_z_bits(125),rr.get_z_bits(125));
  }

  template<> inline mpq_class random<mpq_class>(const mpq_class& a, const mpq_class& b)
  {
    return a + (b-a) * random<mpq_class>();
  }

  inline bool isMuchSmallerThan(const mpq_class& a, const mpq_class& b, const mpq_class& eps)
  {
    return ::abs(a) <= ::abs(b) * eps;
  }

  inline bool isApprox(const mpq_class& a, const mpq_class& b, const mpq_class& eps)
  {
      return ::abs(a-b)<eps;
  }

  inline bool isApproxOrLessThan(const mpq_class& a, const mpq_class& b, const mpq_class& eps)
  {
    return a <= b;
  }

  template<> inline long double cast<mpq_class,long double>(const mpq_class& x)
  { return x.get_d(); }

  template<> inline double cast<mpq_class,double>(const mpq_class& x)
  { return x.get_d(); }

  template<> inline long cast<mpq_class,long>(const mpq_class& x)
  { return x.get_d(); }

  template<> inline int cast<mpq_class,int>(const mpq_class& x)
  { return x.get_d(); }

  // G+Smo
  template<> inline size_t cast<mpq_class,size_t>(const mpq_class& x)
  { return x.get_d(); }

  template<> inline unsigned cast<mpq_class,unsigned>(const mpq_class& x)
  { return x.get_d(); }

    // Specialize GEBP kernel and traits
    template<>
    class gebp_traits<mpq_class, mpq_class, false, false>
    {
    public:
      typedef mpq_class ResScalar;
      enum {
        nr = 2, // must be 2 for proper packing...
        mr = 1,
        WorkSpaceFactor = nr,
        LhsProgress = 1,
        RhsProgress = 1
      };
    };

    template<typename Index, int mr, int nr, bool ConjugateLhs, bool ConjugateRhs>
    struct gebp_kernel<mpq_class,mpq_class,Index,mr,nr,ConjugateLhs,ConjugateRhs>
    {
      typedef mpq_class num_t;

      EIGEN_DONT_INLINE
      void operator()(num_t* res, Index resStride, const num_t* blockA, const num_t* blockB, 
                      Index rows, Index depth, Index cols, num_t alpha,
                      Index strideA=-1, Index strideB=-1, Index offsetA=0, 
                      Index offsetB=0, num_t* /*unpackedB*/ = 0)
      {
        num_t acc1, acc2, tmp;
        
        if(strideA==-1) strideA = depth;
        if(strideB==-1) strideB = depth;

        for(Index j=0; j<cols; j+=nr)
        {
          Index actual_nr = (std::min<Index>)(nr,cols-j);
          num_t *C1 = res + j*resStride;
          num_t *C2 = res + (j+1)*resStride;
          for(Index i=0; i<rows; i++)
          {
            num_t *B = const_cast<num_t*>(blockB) + j*strideB + offsetB*actual_nr;
            num_t *A = const_cast<num_t*>(blockA) + i*strideA + offsetA;
            acc1 = 0;
            acc2 = 0;
            for(Index k=0; k<depth; k++)
            {
              mpq_mul(tmp.__get_mp(), A[k].__get_mp(), B[0].__get_mp());
              mpq_add(acc1.__get_mp(), acc1.__get_mp(), tmp.__get_mp());
              
              if(actual_nr==2) 
              {
                mpq_mul(tmp.__get_mp(), A[k].__get_mp(), B[1].__get_mp());
                mpq_add(acc2.__get_mp(), acc2.__get_mp(), tmp.__get_mp());
              }
              
              B+=actual_nr;
            }
            
            mpq_mul(acc1.__get_mp(), acc1.__get_mp(), alpha.__get_mp());
            mpq_add(C1[i].__get_mp(), C1[i].__get_mp(), acc1.__get_mp());
            
            if(actual_nr==2) 
            {
              mpq_mul(acc2.__get_mp(), acc2.__get_mp(), alpha.__get_mp());
              mpq_add(C2[i].__get_mp(), C2[i].__get_mp(), acc2.__get_mp());
            }
          }
        }
      }
    };

  } // end namespace internal

// Complex scalar division.
template <class U, class V, class W, class Y>
std::complex<mpq_class> 
cdiv(const __gmp_expr<mpq_t, U> & xr,
     const __gmp_expr<mpq_t, V> & xi,
     const __gmp_expr<mpq_t, W> & yr,
     const __gmp_expr<mpq_t, Y> & yi)
{
  //using std::abs;
  mpq_class r,d;
  if (abs(yr) > abs(yi))
  {
      r = yi/yr;
      d = yr + r*yi;
      return std::complex<mpq_class>(mpq_class((xr + r*xi)/d), mpq_class((xi - r*xr)/d));
  }
  else
  {
      r = yr/yi;
      d = yi + r*yr;
      return std::complex<mpq_class>(mpq_class((r*xr + xi)/d), mpq_class((r*xi - xr)/d));
  }
}

}//namaspace Eigen

