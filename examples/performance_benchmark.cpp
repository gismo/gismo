/** @file performance_benchmark.cpp

    @brief G+Smo performance benchmark

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

//! [Implement make_vector]
template<typename T>
std::vector<T> make_vector(T value, std::size_t size)
{
  std::vector<T> v;
  for (std::size_t i=0; i<size; ++i)
    v.push_back(value);
  return v;
}
//! [Implement make_vector]

//! [Implement memory safeguard]
template<typename T>
class memory_safeguard
{
public:
  template<typename...Args>
  memory_safeguard(Args... args)
  {
    if (T::size(args...) > gsSysInfo::getMemoryInBytes())
      GISMO_ERROR("Insufficient memory");
  }
};
//! [Implement memory safeguard]

//! [Implement benchmark native C array memcopy]
/**
 * Benchmark: native C array memcopy
 */
template<typename T>
class benchmark_c_array_memcopy
{
private:
  memory_safeguard<benchmark_c_array_memcopy> _msg;
  index_t n;
  T *m_x, *m_y;

public:
  benchmark_c_array_memcopy(index_t n)
    : _msg(n), n(n), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
  }

  ~benchmark_c_array_memcopy()
  {
    delete[] m_x;
    delete[] m_y;
  }

  uint64_t operator()()
  {
#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_y[i] = m_x[i];

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_y[n-1];
    GISMO_UNUSED(tmp);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (2 * uint64_t(n) * sizeof(T));
  }

  static std::string descr()
  {
    return "Memory copy (native C array)";
  }

  static std::string label()
  {
    return "memcopyCarray";
  }

  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark native C array memcopy]

//! [Implement benchmark native C array dot-product]
/**
 * Benchmark: native C array dot-product
 */
template<typename T>
class benchmark_c_array_dotproduct
{
private:
  memory_safeguard<benchmark_c_array_dotproduct> _msg;
  const index_t n;
  T *m_x, *m_y;

public:
  benchmark_c_array_dotproduct(index_t n)
    : _msg(n), n(n), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;
  }

  ~benchmark_c_array_dotproduct()
  {
    delete[] m_x;
    delete[] m_y;
  }

  uint64_t operator()()
  {
    volatile T sum = 0.0;

#pragma omp parallel for simd reduction(+:sum) schedule(static)
    for (index_t i=0; i<n; ++i)
      sum += m_x[i] * m_y[i];

    GISMO_UNUSED(sum);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (2 * uint64_t(n) * sizeof(T));
  }

  static std::string descr()
  {
    return "Dot product (native C array)";
  }

  static std::string label()
  {
    return "dotproductCarray";
  }

  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark native C array dot-product]

//! [Implement benchmark native C array AXPY]
/**
 * Benchmark: native C array AXPY
 */
template<typename T>
class benchmark_c_array_axpy
{
private:
  memory_safeguard<benchmark_c_array_axpy> _msg;
  const index_t n;
  T *m_x, *m_y, *m_z;

public:
  benchmark_c_array_axpy(index_t n)
    : _msg(n), n(n), m_x(new T[n]), m_y(new T[n]), m_z(new T[n])
  {
#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;
  }

  ~benchmark_c_array_axpy()
  {
    delete[] m_x;
    delete[] m_y;
    delete[] m_z;
  }

  uint64_t operator()()
  {
#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_z[i] = (T)3.414 * m_x[i] + m_y[i];

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_z[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * 3 * n;
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (3 * uint64_t(n) * sizeof(T));
  }

  static std::string descr()
  {
    return "AXPY (native C array)";
  }

  static std::string label()
  {
    return "axpyCarray";
  }

  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark native C array AXPY]

//! [Implement benchmark native C array dense matrix-vector multiplication]
/**
 * Benchmark: native C array dense matrix-vector multiplication
 */
template<typename T>
class benchmark_c_array_dense_matmul
{
private:
  memory_safeguard<benchmark_c_array_dense_matmul> _msg;
  const index_t n;
  T *m_A, *m_x, *m_y;

public:
  benchmark_c_array_dense_matmul(index_t n)
    : _msg(n), n(n), m_A(new T[n*n]), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n*n; ++i)
      m_A[i] = (T)1.0;

#pragma omp parallel for simd schedule(static)
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
  }

  ~benchmark_c_array_dense_matmul()
  {
    delete[] m_A;
    delete[] m_x;
    delete[] m_y;
  }

  uint64_t operator()()
  {
#pragma omp parallel for schedule(static)
    for (index_t i=0; i<n; ++i) {
      T sum = (T)0.0;
#pragma omp simd
      for (index_t j=0; j<n; ++j) {
        sum += m_A[n*i+j] * m_x[j];
      }
      m_y[i] = sum;
    }

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_y[n-1];
    GISMO_UNUSED(tmp);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (2 * uint64_t(n) * uint64_t(n) + uint64_t(n)) * sizeof(T);
  }

  static std::string descr()
  {
    return "Dense matrix-vector multiplication (native C array)";
  }

  static std::string label()
  {
    return "densematmulCarray";
  }
  
  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark native C array dense matrix-vector multiplication]

//! [Implement benchmark eigen vector memcopy]
/**
 * Benchmark: Eigen vector memcopy
 */
template<typename T>
class benchmark_eigen_memcopy
{
private:
  memory_safeguard<benchmark_eigen_memcopy> _msg;
  const index_t n;
  gsVector<T> x,y;

public:
  benchmark_eigen_memcopy(index_t n)
    : _msg(n), n(n), x(n), y(n)
  {
    x.fill((T)1.0);
  }

  uint64_t operator()()
  {
    y.noalias() = x;

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = y[n-1];
    GISMO_UNUSED(tmp);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (2 * uint64_t(n) * sizeof(T));
  }

  static std::string descr()
  {
    return "Memory copy (gsVector)";
  }

  static std::string label()
  {
    return "memcopyEigen";
  }
  
  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark eigen vector memcopy]

//! [Implement benchmark eigen vector dot-product]
/**
 * Benchmark: Eigen vector dot-product
 */
template<typename T>
class benchmark_eigen_dotproduct
{
private:
  memory_safeguard<benchmark_eigen_dotproduct> _msg;
  const index_t n;
  gsVector<T> x, y;

public:
  benchmark_eigen_dotproduct(index_t n)
    : _msg(n), n(n), x(n), y(n)
  {
    x.fill((T)1.0);
    y.fill((T)1.0);
  }

  uint64_t operator()()
  {
    volatile T sum = y.dot(x);
    GISMO_UNUSED(sum);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (2 * uint64_t(n) * sizeof(T));
  }

  static std::string descr()
  {
    return "Dot product (gsVector)";
  }

  static std::string label()
  {
    return "dotproductEigen";
  }

  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark eigen vector dot-product]

//! [Implement benchmark eigen vector AXPY]
/**
 * Benchmark: Eigen vector AXPY
 */
template<typename T>
class benchmark_eigen_axpy
{
private:
  memory_safeguard<benchmark_eigen_axpy> _msg;
  const index_t n;
  gsVector<T> x, y, z;

public:
  benchmark_eigen_axpy(index_t n)
    : _msg(n), n(n), x(n), y(n), z(n)
  {
    x.fill((T)1.0);
    y.fill((T)1.0);
  }

  uint64_t operator()()
  {
    z.noalias() = (T)3.141*x + y;

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = z[n-1];
    GISMO_UNUSED(tmp);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (3 * uint64_t(n) * sizeof(T));
  }

  static std::string descr()
  {
    return "AXPY (gsVector)";
  }

  static std::string label()
  {
    return "axpyEigen";
  }
  
  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark eigen vector AXPY]

//! [Implement benchmark eigen dense matrix-vector multiplication]
/**
 * Benchmark: Eigen dense matrix-vector multiplication
 */
template<typename T>
class benchmark_eigen_dense_matmul
{
private:
  memory_safeguard<benchmark_eigen_dense_matmul> _msg;
  const index_t n;
  gsMatrix<T> A;
  gsVector<T> x, y;

public:
  benchmark_eigen_dense_matmul(index_t n)
    : _msg(n), n(n), A(n,n), x(n), y(n)
  {
    A.fill(1.0);
    x.fill(1.0);
  }

  uint64_t operator()()
  {
    y.noalias() = A*x;

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = y[n-1];
    GISMO_UNUSED(tmp);

    return size();
  }

  constexpr uint64_t size() const
  {
    return size(n);
  }

  static constexpr uint64_t size(index_t n)
  {
    return (2 * uint64_t(n) * uint64_t(n) + uint64_t(n)) * sizeof(T);
  }

  static std::string descr()
  {
    return "Dense matrix-vector multiplication (gsMatrix/gsVector)";
  }

  static std::string label()
  {
    return "densematmulEigen";
  }
  
  static constexpr gismo::metric metric()
  {
    return gismo::metric::bandwidth_gb_sec;
  }
};
//! [Implement benchmark eigen dense matrix-vector multiplication]

//! [Implement benchmark Poisson 2d visitor]
/**
 * Benchmark: Visitor-based Poisson 2d
 */
template<typename T>
class benchmark_poisson2d_visitor
{
private:
  memory_safeguard<benchmark_poisson2d_visitor> _msg;
  int numPatches, numRefine, degree;
  gsMultiPatch<T> geo;
  gsMultiBasis<T> bases;
  gsConstantFunction<T> f;
  gsBoundaryConditions<T> bc;
  gsPoissonAssembler<T> assembler;

public:
  template<typename... Args>
  benchmark_poisson2d_visitor(std::tuple<Args...> args)
    : benchmark_poisson2d_visitor(std::get<0>(args), std::get<1>(args), std::get<2>(args))
  {}
    
  benchmark_poisson2d_visitor(int numPatches, int numRefine=0, int degree=1)
    : _msg(numPatches, numRefine, degree),
      numPatches(numPatches), numRefine(numRefine), degree(degree),
      geo(gsNurbsCreator<>::BSplineSquareGrid(numPatches, numPatches, 1.0)),
      bases(geo), f(0.0, 0.0, 2)
  {
    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
      bases.uniformRefine();

    // k-refinement (set degree)
    for (std::size_t i = 0; i < bases.nBases(); ++ i)
      bases[i].setDegreePreservingMultiplicity(degree);

    // create assembler
    assembler = gsPoissonAssembler<T>(geo, bases, bc, f, dirichlet::nitsche, iFace::glue);
  }

  uint64_t operator()()
  {
    assembler.assemble();
    return sizeof(T) * (assembler.matrix().nonZeros() + assembler.rhs().rows());
  }

  constexpr uint64_t size() const
  {
    return size(numPatches, numRefine, degree);
  }

  static constexpr uint64_t size(index_t numPatches, index_t numRefine, index_t degree)
  {
    // Estimated memory
    // system matrix : 1.33 * ndofs * (2*p+1)^2
    // r.h.s. vector :        ndofs
    //
    // The factor 1.33 is used because Eigen shows better performance
    // if 33% more memory is allocated during the step-by-step assembly
    return sizeof(T) * ( 1.33 * math::pow(2*degree+1,2) + 1 ) *
      (/* numPatches^2 * DOFs per patch */
       math::pow(numPatches,2) * math::pow((1<<numRefine)+degree,(T)2)
       /* remove duplicate DOFs at patch interfaces (2 directions) */
       - 2 * (numPatches * (numPatches-1)) * math::pow( (1<<numRefine)+degree,(T)1)
       /* add interior points at patch corners that have been removed before */
       + math::pow(numPatches-1,2) );
  }

  static std::string descr()
  {
    return "Visitor-based Poisson 2d assembler";
  }

  static std::string label()
  {
    return "assemble2dVisitorAssembler";
  }
  
  static constexpr gismo::metric metric()
  {
    return (gismo::metric)(gismo::metric::runtime_sec + gismo::metric::speedup);
  }
};
//! [Implement benchmark Poisson 2d visitor]

//! [Implement benchmark Poisson 3d visitor]
/**
 * Benchmark: Visitor-based Poisson 3d assembler
 */
template<typename T>
class benchmark_poisson3d_visitor
{
private:
  memory_safeguard<benchmark_poisson3d_visitor> _msg;
  int numPatches, numRefine, degree;
  gsMultiPatch<T> geo;
  gsMultiBasis<T> bases;
  gsConstantFunction<T> f;
  gsBoundaryConditions<T> bc;
  gsPoissonAssembler<T> assembler;

public:
  template<typename... Args>
  benchmark_poisson3d_visitor(std::tuple<Args...> args)
    : benchmark_poisson3d_visitor(std::get<0>(args), std::get<1>(args), std::get<2>(args))
  {}
  
  benchmark_poisson3d_visitor(int numPatches, int numRefine=0, int degree=1)
    : _msg(numPatches, numRefine, degree),
      numPatches(numPatches), numRefine(numRefine), degree(degree),
      geo(gsNurbsCreator<>::BSplineCubeGrid(numPatches, numPatches, numPatches, 1.0)),
      bases(geo), f(0.0, 0.0, 0.0, 3)
  {
    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
      bases.uniformRefine();

    // k-refinement (set degree)
    for (std::size_t i = 0; i < bases.nBases(); ++ i)
      bases[i].setDegreePreservingMultiplicity(degree);

    // create assembler
    assembler = gsPoissonAssembler<T>(geo, bases, bc, f, dirichlet::nitsche, iFace::glue);
  }

  uint64_t operator()()
  {
    assembler.assemble();
    return sizeof(T) * (assembler.matrix().nonZeros() + assembler.rhs().rows());
  }

  constexpr uint64_t size() const
  {
    return size(numPatches, numRefine, degree);
  }

  static constexpr uint64_t size(index_t numPatches, index_t numRefine, index_t degree)
  {
    // Estimated memory
    // system matrix : 1.33 * ndofs * (2*p+1)^3
    // r.h.s. vector :        ndofs
    //
    // The factor 1.33 is used because Eigen shows better performance
    // if 33% more memory is allocated during the step-by-step assembly
    return sizeof(T) * 1.33 * (numPatches * ((1<<numRefine)+degree-1)+1) *
      (/* numPatches^2 * DOFs per patch */
       math::pow(numPatches,2) * math::pow((1<<numRefine)+degree,(T)2)
       /* remove duplicate DOFs at patch interfaces (2 directions) */
       - 2 * (numPatches * (numPatches-1)) * math::pow( (1<<numRefine)+degree,(T)1)
       /* add interior points at patch corners that have been removed before */
       + math::pow(numPatches-1,2) );
  }

  static std::string descr()
  {
    return "Visitor-based Poisson 3d assembler";
  }

  static std::string label()
  {
    return "assemble3dVisitorAssembler";
  }
  
  static constexpr gismo::metric metric()
  {
    return (gismo::metric)(gismo::metric::runtime_sec + gismo::metric::speedup);
  }
};
//! [Implement benchmark Poisson 3d visitor]

//! [Implement benchmark Poisson 2d expression assembler]
/**
 * Benchmark: Expression assembler-based Poisson 2d
 */
template<typename T>
class benchmark_poisson2d_expression_assembler
{
private:
  memory_safeguard<benchmark_poisson2d_expression_assembler> _msg;
  int numPatches, numRefine, degree;
  gsMultiPatch<T> geo;
  gsMultiBasis<T> bases;
  gsBoundaryConditions<T> bc;
  
  gsExprAssembler<T> A;
  typename gsExprAssembler<>::geometryMap G;
  typename gsExprAssembler<>::space u;

  gsFunctionExpr<T> f;
  expr::gsComposition<T> ff;
  
public:
  template<typename... Args>
  benchmark_poisson2d_expression_assembler(std::tuple<Args...> args)
    : benchmark_poisson2d_expression_assembler(std::get<0>(args), std::get<1>(args), std::get<2>(args))
  {}
    
  benchmark_poisson2d_expression_assembler(int numPatches, int numRefine=0, int degree=1)
    : _msg(numPatches, numRefine, degree),
      numPatches(numPatches), numRefine(numRefine), degree(degree),
      geo(gsNurbsCreator<>::BSplineSquareGrid(numPatches, numPatches, 1.0)),
      bases(geo, true), A(1,1), G(A.getMap(geo)), u(A.getSpace(bases)),
      f("0.0", 2), ff(A.getCoeff(f, G))
  {    
    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
      bases.uniformRefine();
    
    // k-refinement (set degree)
    for (std::size_t i = 0; i < bases.nBases(); ++ i)
      bases[i].setDegreePreservingMultiplicity(degree);
    
    // set the geometry map to boundary conditions
    bc.setGeoMap(geo);
    
    // setup boundary conditions
    u.setup(bc, dirichlet::l2Projection, 0);

    // set elements used for numerical integration
    A.setIntegrationElements(bases);
    
    // initialize the system
    A.initSystem();   
  }

  uint64_t operator()()
  {
    // Compute the system matrix and right-hand side
    A.assemble(
               igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
               ,
               u * ff * meas(G) //rhs vector
               );
        
    return sizeof(T) * (A.matrix().nonZeros() + A.rhs().rows());
  }

  constexpr uint64_t size() const
  {
    return size(numPatches, numRefine, degree);
  }

  static constexpr uint64_t size(index_t numPatches, index_t numRefine, index_t degree)
  {
    // Estimated memory
    // system matrix : 1.33 * ndofs * (2*p+1)^2
    // r.h.s. vector :        ndofs
    //
    // The factor 1.33 is used because Eigen shows better performance
    // if 33% more memory is allocated during the step-by-step assembly
    return sizeof(T) * ( 1.33 * math::pow(2*degree+1,2) + 1 ) *
      (/* numPatches^2 * DOFs per patch */
       math::pow(numPatches,2) * math::pow((1<<numRefine)+degree,(T)2)
       /* remove duplicate DOFs at patch interfaces (2 directions) */
       - 2 * (numPatches * (numPatches-1)) * math::pow( (1<<numRefine)+degree,(T)1)
       /* add interior points at patch corners that have been removed before */
       + math::pow(numPatches-1,2) );
  }

  static std::string descr()
  {
    return "Expression assembler-based Poisson 2d assembler";
  }

  static std::string label()
  {
    return "assemble2dExpressionAssembler";
  }
  
  static constexpr gismo::metric metric()
  {
    return (gismo::metric)(gismo::metric::runtime_sec + gismo::metric::speedup);
  }
};
//! [Implement benchmark Poisson 2d expression assembler]

//! [Implement benchmark Poisson 3d expression assembler]
/**
 * Benchmark: Expression assembler-based Poisson 3d
 */
template<typename T>
class benchmark_poisson3d_expression_assembler
{
private:
  memory_safeguard<benchmark_poisson3d_expression_assembler> _msg;
  int numPatches, numRefine, degree;
  gsMultiPatch<T> geo;
  gsMultiBasis<T> bases;
  gsBoundaryConditions<T> bc;
  
  gsExprAssembler<T> A;
  typename gsExprAssembler<>::geometryMap G;
  typename gsExprAssembler<>::space u;

  gsFunctionExpr<T> f;
  expr::gsComposition<T> ff;
  
public:
  template<typename... Args>
  benchmark_poisson3d_expression_assembler(std::tuple<Args...> args)
    : benchmark_poisson3d_expression_assembler(std::get<0>(args), std::get<1>(args), std::get<2>(args))
  {}
    
  benchmark_poisson3d_expression_assembler(int numPatches, int numRefine=0, int degree=1)
    : _msg(numPatches, numRefine, degree),
      numPatches(numPatches), numRefine(numRefine), degree(degree),
      geo(gsNurbsCreator<>::BSplineCubeGrid(numPatches, numPatches, numPatches, 1.0)),
      bases(geo, true), A(1,1), G(A.getMap(geo)), u(A.getSpace(bases)),
      f("0.0", 3), ff(A.getCoeff(f, G))
  {    
    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
      bases.uniformRefine();
    
    // k-refinement (set degree)
    for (std::size_t i = 0; i < bases.nBases(); ++ i)
      bases[i].setDegreePreservingMultiplicity(degree);
    
    // set the geometry map to boundary conditions
    bc.setGeoMap(geo);
    
    // setup boundary conditions
    u.setup(bc, dirichlet::l2Projection, 0);

    // set elements used for numerical integration
    A.setIntegrationElements(bases);
    
    // initialize the system
    A.initSystem();    
  }

  uint64_t operator()()
  {
    // Compute the system matrix and right-hand side
    A.assemble(
               igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
               ,
               u * ff * meas(G) //rhs vector
               );
        
    return sizeof(T) * (A.matrix().nonZeros() + A.rhs().rows());
  }

  constexpr uint64_t size() const
  {
    return size(numPatches, numRefine, degree);
  }

  static constexpr uint64_t size(index_t numPatches, index_t numRefine, index_t degree)
  {
    // Estimated memory
    // system matrix : 1.33 * ndofs * (2*p+1)^3
    // r.h.s. vector :        ndofs
    //
    // The factor 1.33 is used because Eigen shows better performance
    // if 33% more memory is allocated during the step-by-step assembly
    return sizeof(T) * 1.33 * (numPatches * ((1<<numRefine)+degree-1)+1) *
      (/* numPatches^2 * DOFs per patch */
       math::pow(numPatches,2) * math::pow((1<<numRefine)+degree,(T)2)
       /* remove duplicate DOFs at patch interfaces (2 directions) */
       - 2 * (numPatches * (numPatches-1)) * math::pow( (1<<numRefine)+degree,(T)1)
       /* add interior points at patch corners that have been removed before */
       + math::pow(numPatches-1,2) );
  }

  static std::string descr()
  {
    return "Expression assembler-based Poisson 3d assembler";
  }

  static std::string label()
  {
    return "assemble3dExpressionAssembler";
  }
  
  static constexpr gismo::metric metric()
  {
    return (gismo::metric)(gismo::metric::runtime_sec + gismo::metric::speedup);
  }
};
//! [Implement benchmark Poisson 3d expression assembler]

int main(int argc, char *argv[])
{
  //! [Parse command line]
  gsBenchmark benchmark;
  std::string fn;
  bool list=false, all=false;
  std::vector<index_t>  benchmarks, msizes, nruns, nthreads, patches, subdivides, vsizes;
  index_t msizemin      = 10;
  index_t nrunsmax      = 100;
  index_t nrunsmin      = 1;
  index_t patchesmax    = 128;
  index_t patchesmin    = 1;
  index_t subdividemax  = 10;
  index_t subdividemin  = 0;
  index_t vsizemin      = 100;
  real_t  patchesfactor = 2;
  real_t  msizesfactor  = 2;
  real_t  nrunsfactor   = 1.5;
  real_t  vsizesfactor  = 4;
  index_t msizemax = (index_t) math::min((real_t)std::numeric_limits<index_t>::max(),
                                         std::sqrt((real_t)(0.8) * sizeof(real_t)*gsSysInfo::getMemoryInBytes()));
  index_t vsizemax = (index_t) math::min((real_t)std::numeric_limits<index_t>::max(),
                                         (real_t)(0.8) * sizeof(real_t)*gsSysInfo::getMemoryInBytes());

  gsCmdLine cmd("G+Smo performance benchmark.");
  cmd.printVersion();

  cmd.addReal("M", "msizesfactor", "Growth factor for the sequence of msizes (only used if '-m' is not given)", msizesfactor);
  cmd.addReal("P", "patchesfactor", "Growth factor for the sequence of patches (only used if '-p' is not given)", patchesfactor);
  cmd.addReal("R", "runsfactor", "Growth factor for the sequence of runs (only used if '-r' is not given)", nrunsfactor);
  cmd.addReal("V", "vsizesfactor", "Growth factor for the sequence of vsizes (only used if '-v' is not given)", vsizesfactor);
  cmd.addInt("", "msizemax", "Maximum number of unknowns in matrix/vector benchmarks (only used if '-m' is not given)", msizemax);
  cmd.addInt("", "msizemin", "Minimum number of unknowns in matrix/vector benchmarks (only used if '-m'is not given)", msizemin);
  cmd.addInt("", "patchesmax", "Maximum number of patches in assembly benchmarks (only used if '-p' is not given)", patchesmax);
  cmd.addInt("", "patchesmin", "Minimum number of patches in assembly benchmarks (only used if '-p' is not given)", patchesmin);
  cmd.addInt("", "runsmax", "Maximum number of runs (only used if '-r' is not given)", nrunsmax);
  cmd.addInt("", "runsmin", "Mminimum number of runs (only used if '-r' is not given)", nrunsmin);
  cmd.addInt("", "subdividemax", "Maximum number of subdivisions (h-refinement) in assembly benchmarks (only used if '-r' is not given)", subdividemax);
  cmd.addInt("", "subdividemin", "Minimum number of subdivisions (h-refinement) in assembly benchmarks (only used if '-r' is not given)", subdividemin);
  cmd.addInt("", "vsizemax", "Maximum number of unknowns in vector benchmarks (only used if '-v' is not given)", vsizemax);
  cmd.addInt("", "vsizemin", "Mminimum number of unknowns in vector benchmarks (only used if '-v' is not given)", vsizemin);
  cmd.addMultiInt("b", "benchmarks", "List of benchmarks to be run", benchmarks);
  cmd.addMultiInt("m", "msizes", "Number of unknowns in matrix/vector benchmarks (auto-generated if not given)", msizes);
  cmd.addMultiInt("p", "patches", "Number of patches in assembly benchmarks (auto-generated if not given)", patches);
  cmd.addMultiInt("r", "runs", "Number of runs over which the results are averaged (auto-generated if not given)", nruns);
  cmd.addMultiInt("s", "subdivide", "Number of subdivisions (h-refinement) in assembly benchmarks (auto-generated if not given)", subdivides);
  cmd.addMultiInt("t", "threads", "Number of OpenMP threads to be used for the benchmark (auto-generated if not given)", nthreads);
  cmd.addMultiInt("v", "vsizes", "Number of unknowns in vector benchmarks (auto-generated if not given)", vsizes);
  cmd.addString("o", "output", "Name of the output file", fn);
  cmd.addSwitch("list", "List all benchmarks and exit", list);
  cmd.addSwitch("all", "Run all benchmarks", all);
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  //! [Parse command line]

  //! [List benchmarks and exit]
  if (list) {
    gsInfo << "\nThe following benchmarks are available:\n"
           << "#01: " << benchmark_c_array_memcopy<real_t>::descr() << "\n"
           << "#02: " << benchmark_eigen_memcopy<real_t>::descr() << "\n"
           << "#03: " << benchmark_c_array_dotproduct<real_t>::descr() << "\n"
           << "#04: " << benchmark_eigen_dotproduct<real_t>::descr() << "\n"
           << "#05: " << benchmark_c_array_axpy<real_t>::descr() << "\n"
           << "#06: " << benchmark_eigen_axpy<real_t>::descr() << "\n"
           << "#07: " << benchmark_c_array_dense_matmul<real_t>::descr() << "\n"
           << "#08: " << benchmark_eigen_dense_matmul<real_t>::descr() << "\n"
           << "#09: " << benchmark_poisson2d_visitor<real_t>::descr()
           <<            " with increasing number of patches" << "\n"
           << "#10: " << benchmark_poisson2d_visitor<real_t>::descr()
           <<            " with increasing number of subdivisions" << "\n"
           << "#11: " << benchmark_poisson3d_visitor<real_t>::descr()
           <<            " with increasing number of patches" << "\n"
           << "#12: " << benchmark_poisson3d_visitor<real_t>::descr()
           <<            " with increasing number of subdivisions" << "\n"
           << "#13: " << benchmark_poisson2d_expression_assembler<real_t>::descr()
           <<            " with increasing number of patches" << "\n"
           << "#14: " << benchmark_poisson2d_expression_assembler<real_t>::descr()
           <<            " with increasing number of subdivisions" << "\n"
           << "#15: " << benchmark_poisson3d_expression_assembler<real_t>::descr()
           <<            " with increasing number of patches" << "\n"
           << "#16: " << benchmark_poisson3d_expression_assembler<real_t>::descr()
           <<            " with increasing number of subdivisions" << "\n";
      
    return EXIT_SUCCESS;
  }
  //! [List benchmarks and exit]

  //! [Default configuration]
  // If empty fill with all benchmarks 1, 2, ...
  if (all) {
    benchmarks.clear();
    for(index_t i=1; i<=16; ++i)
      benchmarks.push_back(i);
  }
  
  // If empty fill with 1, 2, 4, ..., maximum number of OpenMP threads
  if (nthreads.empty()) {
    for(index_t i=1; i<=omp_get_max_threads(); i*=2)
      nthreads.push_back(i);
  }

  // If empty fill with msizemin*msizesfactor^k, k=0, 1, 2, ..., msizemax
  if (msizes.empty()) {
    for(index_t i=msizemin;;) {
      msizes.push_back(i);
      if (i<=math::min(msizemax, std::numeric_limits<index_t>::max()) / (msizesfactor*msizesfactor))
        i*=msizesfactor;
      else
        break;
    }
  }

  // If empty fill with patchesmin, ..., patchesmax
  if (patches.empty()) {
    for(index_t i=patchesmin; i<=patchesmax; i*=patchesfactor)
      patches.push_back(i);
  }

  // If empty fill with subdividemin, ..., subdividemax
  if (subdivides.empty()) {
    for(index_t i=subdividemin; i<subdividemax; ++i)
      subdivides.push_back(i);
  }
  
  // If empty fill with vsizemin*vsizesfactor^k, k=0, 1, 2, ..., vsizemax
  if (vsizes.empty()) {
    for(index_t i=vsizemin;;) {
      vsizes.push_back(i);
      if (i<=math::min(vsizemax, std::numeric_limits<index_t>::max()) / vsizesfactor)
        i*=vsizesfactor;
      else
        break;
    }
  }

  // If empty fill with nrunsmax/nrunsfactor^k, k=0, 1, 2, ..., nrunsmin
  if (nruns.empty()) {
    index_t k = nrunsmax;
    for(index_t i=0; i<(index_t)math::max(msizes.size(), patches.size(),
                                         subdivides.size(), vsizes.size()); ++i) {
      nruns.push_back(k);
      k = math::max(nrunsmin, (index_t)(k/nrunsfactor));
    }
  }

  if (nruns.size()<math::max(msizes.size(),vsizes.size()))
    GISMO_ERROR("|nruns| must have the same size as max(|msizes|,|vsizes|)");
  //! [Default configuration]

  //! [Execute benchmarks]
  for (auto bit=benchmarks.cbegin(); bit!=benchmarks.cend(); ++bit) {
    switch((index_t)(*bit)) {


    case (1): {
      // Benchmark: memcopy native C arrays
      benchmark.create<benchmark_c_array_memcopy<real_t> >
        (vsizes, nruns, nthreads);
      break;
    }
      
    case (2): {
      // Benchmark: memcopy gsVector
      benchmark.create<benchmark_eigen_memcopy<real_t> >
        (vsizes, nruns, nthreads);
      break;
    }

    case (3): {
      // Benchmark: dot-product native C array
      benchmark.create<benchmark_c_array_dotproduct<real_t> >
        (vsizes, nruns, nthreads);
      break;
    }

    case (4): {
      // Benchmark: dot-product gsVector
      benchmark.create<benchmark_eigen_dotproduct<real_t> >
        (vsizes, nruns, nthreads);
      break;
    }

    case (5): {
      // Benchmark: axpy native C array
      benchmark.create<benchmark_c_array_axpy<real_t> >
        (vsizes, nruns, nthreads);
      break;
    }

    case (6): {
      // Benchmark: axpy gsVector
      benchmark.create<benchmark_eigen_axpy<real_t> >
        (vsizes, nruns, nthreads);
      break;
    }

    case (7): {
      // Benchmark: dense matrix-vector multiplication native C array
      benchmark.create<benchmark_c_array_dense_matmul<real_t> >
        (msizes, nruns, nthreads);
      break;
    }

    case (8): {
      // Benchmark: dense matrix-vector multiplication gsMatrix/gsVector
      benchmark.create<benchmark_eigen_dense_matmul<real_t> >
        (msizes, nruns, nthreads);
      break;
    }

    case (9): {
      // Benchmark: visitor-based Poisson 2d assembler with increasing number of patches
      benchmark.create<benchmark_poisson2d_visitor<real_t> >
        (util::zip(patches,
                   make_vector((index_t)0, patches.size()),  // subdivisions : 0
                   make_vector((index_t)3, patches.size())), // degree       : 3
         nruns, nthreads, " with increasing number of patches (#subdivisions=0, degree=3)");
      break;
    }

    case (10): {
      // Benchmark: visitor-based Poisson 2d assembler with increasing number of subdivisions
      benchmark.create<benchmark_poisson2d_visitor<real_t> >
        (util::zip(make_vector((index_t)1, subdivides.size()),  // patches : 1
                   subdivides,
                   make_vector((index_t)3, subdivides.size())), // degree  : 3
         nruns, nthreads, " with increasing number of subdivisions (#patches=1, degree=3)");
      break;
    }

    case (11): {
      // Benchmark: visitor-based Poisson 3d assembler with increasing number of patches
      benchmark.create<benchmark_poisson3d_visitor<real_t> >
        (util::zip(patches,
                   make_vector((index_t)0, patches.size()),  // subdivisions : 0
                   make_vector((index_t)2, patches.size())), // degree       : 2
         nruns, nthreads, " with increasing number of patches (#subdivisions=0, degree=2)");
      break;
    }

    case (12): {
      // Benchmark: visitor-based Poisson 3d assembler with increasing number of subdivisions
      benchmark.create<benchmark_poisson3d_visitor<real_t> >
        (util::zip(make_vector((index_t)1, subdivides.size()),  // patches : 1
                   subdivides,
                   make_vector((index_t)2, subdivides.size())), // degree  : 2
         nruns, nthreads, " with increasing number of subdivisions (#patches=1, degree=2)");
      break;
    }
      
    case (13): {
      // Benchmark: expression assembler-based Poisson 2d assembler with increasing number of patches
      benchmark.create<benchmark_poisson2d_expression_assembler<real_t> >
        (util::zip(patches,
                   make_vector((index_t)0, patches.size()),  // subdivisions : 0
                   make_vector((index_t)3, patches.size())), // degree       : 3
         nruns, nthreads, " with increasing number of patches (#subdivisions=0, degree=3)");
      break;
    }

    case (14): {
      // Benchmark: expression assembler-based Poisson 2d assembler with increasing number of subdivision
      benchmark.create<benchmark_poisson2d_expression_assembler<real_t> >
        (util::zip(make_vector((index_t)1, subdivides.size()),  // patches : 1
                   subdivides,
                   make_vector((index_t)3, subdivides.size())), // degree  : 3
         nruns, nthreads, " with increasing number of subdivisions (#patches=1, degree=3)");
      break;
    }

    case (15): {
      // Benchmark: expression assembler-based Poisson 3d assembler with increasing number of patches
      benchmark.create<benchmark_poisson3d_expression_assembler<real_t> >
        (util::zip(patches,
                   make_vector((index_t)0, patches.size()),  // subdivisions : 0
                   make_vector((index_t)2, patches.size())), // degree       : 2
         nruns, nthreads, " with increasing number of patches (#subdivisions=0, degree=2)");
      break;
    }

    case (16): {
      // Benchmark: expression assembler-based Poisson 3d assembler with increasing number of subdivision
      benchmark.create<benchmark_poisson3d_expression_assembler<real_t> >
        (util::zip(make_vector((index_t)1, subdivides.size()),  // patches : 1
                   subdivides,
                   make_vector((index_t)2, subdivides.size())), // degree  : 2
         nruns, nthreads, " with increasing number of subdivisions (#patches=1, degree=2)");
      break;
    }
      
    default:
      GISMO_ERROR("Invalid benchmark");
    }

  } // benchmark loop

  { // Memory copy ratio
    auto bmA = benchmark.find(benchmark_c_array_memcopy<real_t>::label());
    auto bmB = benchmark.find(benchmark_eigen_memcopy<real_t>::label());
    
    if (bmA != std::end(benchmark.get()) && bmB != std::end(benchmark.get())) {
      auto bm = util::ratio("memcopyRatio",
                            "Memory copy (gsVector : native C array)", *bmB, *bmA);
      benchmark.get().push_back( give(bm) );
    }
  }

  { // Dot product ratio
    auto bmA = benchmark.find(benchmark_c_array_dotproduct<real_t>::label());
    auto bmB = benchmark.find(benchmark_eigen_dotproduct<real_t>::label());
    
    if (bmA != std::end(benchmark.get()) && bmB != std::end(benchmark.get())) {
      auto bm = util::ratio("dotproductRatio",
                            "Dot product (gsVector : native C array)", *bmB, *bmA);
      benchmark.get().push_back( give(bm) );
    }
  }

  { // AXPY ratio
    auto bmA = benchmark.find(benchmark_c_array_axpy<real_t>::label());
    auto bmB = benchmark.find(benchmark_eigen_axpy<real_t>::label());
    
    if (bmA != std::end(benchmark.get()) && bmB != std::end(benchmark.get())) {
      auto bm = util::ratio("axpyRatio",
                            "AXPY (gsVector : native C array)", *bmB, *bmA);
      benchmark.get().push_back( give(bm) );
    }
  }

  { // Dense matrix-vector multiplication ratio
    auto bmA = benchmark.find(benchmark_c_array_dense_matmul<real_t>::label());
    auto bmB = benchmark.find(benchmark_eigen_dense_matmul<real_t>::label());
    
    if (bmA != std::end(benchmark.get()) && bmB != std::end(benchmark.get())) {
      auto bm = util::ratio("densematmulRatio",
                            "Dense matrix-vector multiplication (gsMatrix/gsVector : native C array)",
                            *bmB, *bmA);
      benchmark.get().push_back( give(bm) );
    }
  }  
  
  if (fn.empty())
    gsInfo << benchmark << "\n";
  else if (gsFileManager::getExtension(fn) == "tex") {
    std::ofstream file;
    file.open(fn);
    benchmark.to_tikz(file);
    file.close();
  }
  else if (gsFileManager::getExtension(fn) == "xml") {
    gsFileData<> file;
    file << benchmark;
    file.save("result.xml");
  }
  else {
    GISMO_ERROR("Unsupported file extension");
  }
  //! [Execute benchmarks]
  
  return EXIT_SUCCESS;
}
