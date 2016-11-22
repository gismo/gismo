// GMP support for Eigen linear algebra library
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Warning: GNU GMP C++ is licensed under GPL/LGPL.
//
// Angelos Mantzaflaris, 2015-2016


#pragma once

#include <gmpxx.h>

namespace Eigen 
{  	

/*
  // Specialize for mpq_class
  namespace numext 
  {
  EIGEN_DEVICE_FUNC
  EIGEN_ALWAYS_INLINE mpq_class mini(const mpq_class& x, const mpq_class& y)
  { return ::min(x,y); }
  
  EIGEN_DEVICE_FUNC
  EIGEN_ALWAYS_INLINE mpq_class maxi(const mpq_class& x, const mpq_class& y)
  { return ::max(x,y); }
  }
*/

// AM: to do
//template <class U>
//template<> struct NumTraits<__gmp_expr<mpq_t, U> >

  template<> struct NumTraits<mpq_class>
    : GenericNumTraits<mpq_class>
  {
    enum {
      IsInteger = 0,
      IsSigned  = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = HugeCost,
      AddCost  = HugeCost,
      MulCost  = HugeCost,
    };

    typedef mpq_class Real;
    typedef mpq_class NonInteger;
    
    // Constants
    static inline Real Pi () { return Real(static_cast<double>(EIGEN_PI)); }

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
        Vectorizable = false,
        LhsPacketSize = 1,
        RhsPacketSize = 1,
        ResPacketSize = 1,
        NumberOfRegisters = 1,
        nr = 1,
        mr = 1,
        LhsProgress = 1,
        RhsProgress = 1
      };
      typedef ResScalar LhsPacket;
      typedef ResScalar RhsPacket;
      typedef ResScalar ResPacket;

    };

    template<typename Index, typename DataMapper, bool ConjugateLhs, bool ConjugateRhs>
    struct gebp_kernel<mpq_class,mpq_class,Index,DataMapper,1,1,ConjugateLhs,ConjugateRhs>
    {
      typedef mpq_class num_t;

      EIGEN_DONT_INLINE
      void operator()(const DataMapper& res, const num_t* blockA, const num_t* blockB, 
                      Index rows, Index depth, Index cols, const num_t& alpha,
                      Index strideA=-1, Index strideB=-1, Index offsetA=0, Index offsetB=0)
      {
        if(rows==0 || cols==0 || depth==0)
          return;

        num_t  acc1(0), tmp(0);        

        if(strideA==-1) strideA = depth;
        if(strideB==-1) strideB = depth;

        for(Index i=0; i<rows; ++i)
        {
          for(Index j=0; j<cols; ++j)
          {
            const num_t *A = blockA + i*strideA + offsetA;
            const num_t *B = blockB + j*strideB + offsetB;
            acc1 = 0;
            for(Index k=0; k<depth; k++)
            {
              mpq_mul(tmp.__get_mp() , A[k].__get_mp(), B[0].__get_mp());
              mpq_add(acc1.__get_mp(), acc1.__get_mp(), tmp.__get_mp() );
            }
            
            mpq_mul(acc1.__get_mp()    , acc1.__get_mp() , alpha.__get_mp());
            mpq_add(res(i,j).__get_mp(), res(i,j).__get_mp(), acc1.__get_mp() );
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

