
#pragma once

#include <gmpxx.h>

namespace Eigen 
{  	
  template<> struct NumTraits<mpq_class>
    : GenericNumTraits<mpq_class>
  {
    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 10,
      AddCost = 10,
      MulCost = 40
    };

    typedef mpq_class Real;
    typedef mpq_class NonInteger;
    
    inline static Real epsilon  ()     {    return 0; }
    inline static Real epsilon  (const Real& x) { return 0; }

    inline static Real dummy_precision()   
    { 
        return 0;
    }
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
      typedef mpq_class mpreal;

      EIGEN_DONT_INLINE
      void operator()(mpreal* res, Index resStride, const mpreal* blockA, const mpreal* blockB, Index rows, Index depth, Index cols, mpreal alpha,
                      Index strideA=-1, Index strideB=-1, Index offsetA=0, Index offsetB=0, mpreal* /*unpackedB*/ = 0)
      {
       /*
        mpreal acc1, acc2, tmp;
        
        if(strideA==-1) strideA = depth;
        if(strideB==-1) strideB = depth;

        for(Index j=0; j<cols; j+=nr)
        {
          Index actual_nr = (std::min<Index>)(nr,cols-j);
          mpreal *C1 = res + j*resStride;
          mpreal *C2 = res + (j+1)*resStride;
          for(Index i=0; i<rows; i++)
          {
            mpreal *B = const_cast<mpreal*>(blockB) + j*strideB + offsetB*actual_nr;
            mpreal *A = const_cast<mpreal*>(blockA) + i*strideA + offsetA;
            acc1 = 0;
            acc2 = 0;
            for(Index k=0; k<depth; k++)
            {
              mpfr_mul(tmp.mpfr_ptr(), A[k].mpfr_ptr(), B[0].mpfr_ptr(), mpreal::get_default_rnd());
              mpfr_add(acc1.mpfr_ptr(), acc1.mpfr_ptr(), tmp.mpfr_ptr(),  mpreal::get_default_rnd());
              
              if(actual_nr==2) {
                mpfr_mul(tmp.mpfr_ptr(), A[k].mpfr_ptr(), B[1].mpfr_ptr(), mpreal::get_default_rnd());
                mpfr_add(acc2.mpfr_ptr(), acc2.mpfr_ptr(), tmp.mpfr_ptr(),  mpreal::get_default_rnd());
              }
              
              B+=actual_nr;
            }
            
            mpfr_mul(acc1.mpfr_ptr(), acc1.mpfr_ptr(), alpha.mpfr_ptr(), mpreal::get_default_rnd());
            mpfr_add(C1[i].mpfr_ptr(), C1[i].mpfr_ptr(), acc1.mpfr_ptr(),  mpreal::get_default_rnd());
            
            if(actual_nr==2) {
              mpfr_mul(acc2.mpfr_ptr(), acc2.mpfr_ptr(), alpha.mpfr_ptr(), mpreal::get_default_rnd());
              mpfr_add(C2[i].mpfr_ptr(), C2[i].mpfr_ptr(), acc2.mpfr_ptr(),  mpreal::get_default_rnd());
            }
          }
        }
        //*/ // commented
      }
    };

  } // end namespace internal
}

