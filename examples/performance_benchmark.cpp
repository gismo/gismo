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

//! [Create benchmark macro]
#define CREATE_BENCHMARK(_benchmark, _label, _sizes, _metric)           \
  gsInfo << "=== " << _benchmark<real_t>::name() << "\n";               \
  auto bmark = benchmark.add(_label, _benchmark<real_t>::name());       \
  for (auto it=_sizes.cbegin(); it!=_sizes.cend(); ++it) {              \
    gsInfo << "... " << (*it) << std::flush;                            \
    try {                                                               \
      _benchmark<real_t> benchmark(*it);                                \
      auto results = gsBenchmark::run(nthreads, nruns, benchmark, _metric); \
      std::string meminfo;                                              \
      uint64_t memsize = benchmark.size();                              \
      if (memsize<1024)                                                 \
        meminfo = util::to_string(memsize)+" B";                        \
      else if (memsize<1024*1024)                                       \
        meminfo = util::to_string(memsize/1024)+" KB";                  \
      else if (memsize<1024*1024*1024)                                  \
        meminfo = util::to_string(memsize/(1024*1024))+" MB";           \
      else                                                              \
        meminfo = util::to_string(memsize/(1024*1024*1024))+" GB";      \
      bmark->add(_label, meminfo, results);                             \
    } catch(...) { gsInfo << "[failed!]"; }                             \
    gsInfo << "\n";                                                     \
  }
//! [Create benchmark macro]

//! [Implement memory safeguard]
template<typename T>
class memory_safeguard
{
public:
  memory_safeguard(index_t n)
  {
    if (T::size(n) > gsSysInfo::getMemoryInBytes())
      throw std::runtime_error("Insufficient memory");
  }
};
//! [Implement memory safeguard]

//! [Implement benchmarks]
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
#pragma omp parallel for simd
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
  }

  ~benchmark_c_array_memcopy()
  {
    delete[] m_x;
    delete[] m_y;
  }

  index_t operator()()
  {
#pragma omp parallel for simd
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
    return (2 * n * sizeof(T));
  }

  static std::string name()
  {
    return "Memory copy (native C array)";
  }
};

/**
 * Benchmark: native C array dot-product
 */
template<typename T>
class benchmark_c_array_dotproduct
{
private:
  memory_safeguard<benchmark_c_array_dotproduct> _msg;
  constexpr index_t n;
  T *m_x, *m_y;

public:
  benchmark_c_array_dotproduct(index_t n)
    : _msg(n), n(n), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd
    for (index_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;
  }

  ~benchmark_c_array_dotproduct()
  {
    delete[] m_x;
    delete[] m_y;
  }

  index_t operator()()
  {
    volatile T sum = 0.0;

#pragma omp parallel for simd reduction(+:sum)
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
    return (2 * n * sizeof(T));
  }

  static std::string name()
  {
    return "Dot-product (native C array)";
  }
};

/**
 * Benchmark: native C array AXPY
 */
template<typename T>
class benchmark_c_array_axpy
{
private:
  memory_safeguard<benchmark_c_array_axpy> _msg;
  constexpr index_t n;
  T *m_x, *m_y, *m_z;

public:
  benchmark_c_array_axpy(index_t n)
    : _msg(n), n(n), m_x(new T[n]), m_y(new T[n]), m_z(new T[n])
  {
#pragma omp parallel for simd
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd
    for (index_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;
  }

  ~benchmark_c_array_axpy()
  {
    delete[] m_x;
    delete[] m_y;
    delete[] m_z;
  }

  index_t operator()()
  {
#pragma omp parallel for simd
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
    return (3 * n * sizeof(T));
  }

  static std::string name()
  {
    return "AXPY (native C array)";
  }
};

/**
 * Benchmark: native C array dense matrix-vector multiplication
 */
template<typename T>
class benchmark_c_array_dense_matmul
{
private:
  memory_safeguard<benchmark_c_array_dense_matmul> _msg;
  constexpr index_t n;
  T *m_A, *m_x, *m_y;

public:
  benchmark_c_array_dense_matmul(index_t n)
    : _msg(n), n(n), m_A(new T[n*n]), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd
    for (index_t i=0; i<n*n; ++i)
      m_A[i] = (T)1.0;

#pragma omp parallel for simd
    for (index_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
  }

  ~benchmark_c_array_dense_matmul()
  {
    delete[] m_A;
    delete[] m_x;
    delete[] m_y;
  }

  index_t operator()()
  {
    for (index_t i=0; i<n; ++i) {
      T sum = (T)0.0;
#pragma omp parallel for simd reduction(+:sum)
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
    return (2 * n * n + n) * sizeof(T);
  }

  static std::string name()
  {
    return "Dense matrix-vector multiplication (native C array)";
  }
};

/**
 * Benchmark: Eigen vector memcopy
 */
template<typename T>
class benchmark_eigen_memcopy
{
private:
  memory_safeguard<benchmark_eigen_memcopy> _msg;
  constexpr index_t n;
  gsVector<T> x,y;

public:
  benchmark_eigen_memcopy(index_t n)
    : _msg(n), n(n), x(n), y(n)
  {
    x.fill((T)0.0);
  }

  index_t operator()()
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
    return (2 * n * sizeof(T));
  }

  static std::string name()
  {
    return "Memory copy (gsVector)";
  }
};

/**
 * Benchmark: Eigen vector dot-product
 */
template<typename T>
class benchmark_eigen_dotproduct
{
private:
  memory_safeguard<benchmark_eigen_dotproduct> _msg;
  constexpr index_t n;
  gsVector<T> x, y;

public:
  benchmark_eigen_dotproduct(index_t n)
    : _msg(n), n(n), x(n), y(n)
  {
    x.fill((T)0.0);
    y.fill((T)0.0);
  }

  index_t operator()()
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
    return (2 * n * sizeof(T));
  }

  static std::string name()
  {
    return "Dot-product (gsVector)";
  }
};

/**
 * Benchmark: Eigen vector AXPY
 */
template<typename T>
class benchmark_eigen_axpy
{
private:
  memory_safeguard<benchmark_eigen_axpy> _msg;
  constexpr index_t n;
  gsVector<T> x, y, z;

public:
  benchmark_eigen_axpy(index_t n)
    : _msg(n), n(n), x(n), y(n), z(n)
  {
    x.fill((T)0.0);
    y.fill((T)0.0);
  }

  index_t operator()()
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
    return (3 * n * sizeof(T));
  }

  static std::string name()
  {
    return "AXPY (gsVector)";
  }
};

/**
 * Benchmark: Eigen dense matrix-vector multiplication
 */
template<typename T>
class benchmark_eigen_dense_matmul
{
private:
  memory_safeguard<benchmark_eigen_dense_matmul> _msg;
  constexpr index_t n;
  gsMatrix<T> A;
  gsVector<T> x, y;

public:
  benchmark_eigen_dense_matmul(index_t n)
    : _msg(n), n(n), A(n,n), x(n), y(n)
  {
    A.fill(0.0);
    x.fill(0.0);
  }

  index_t operator()()
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
    return (2 * n * n + n) * sizeof(T);
  }

  static std::string name()
  {
    return "Dense matrix-vector multiplication (gsMatrix/gsVector)";
  }
};

/**
 * Benchmark: Poisson 2D
 */
template<typename T>
class benchmark_poisson2d_visitor
{
private:
  memory_safeguard<benchmark_poisson2d_visitor> _msg;
  gsMultiPatch<T> geo;
  gsMultiBasis<T> bases;
  gsConstantFunction<T> f;
  gsBoundaryConditions<T> bcInfo;
  gsPoissonAssembler<T> assembler;

public:
  benchmark_poisson2d_visitor(int npatches, int refine=0, int degree=1)
    : _msg(0), geo(gsNurbsCreator<>::BSplineSquareGrid(npatches, npatches, 1.0)),
      bases(geo), f(0.0, 0.0, 2)
  {
    // h-refine each basis
    for (int i = 0; i < refine; ++i)
      bases.uniformRefine();

    // k-refinement (set degree)
    for (std::size_t i = 0; i < bases.nBases(); ++ i)
      bases[i].setDegreePreservingMultiplicity(degree);

    assembler = gsPoissonAssembler<T>(geo, bases, bcInfo, f, dirichlet::nitsche, iFace::glue);
  }

  index_t operator()()
  {
    assembler.assemble();

    return sizeof(T) * assembler.numDofs();
  }

  constexpr uint64_t size() const
  {
    return size(0);
  }

  static constexpr uint64_t size(index_t n)
  {
    return sizeof(T);
  }

  static std::string name()
  {
    return "Visitor-based Poisson2d";
  }
};
//! [Implement benchmarks]




int main(int argc, char *argv[])
{
  //! [Parse command line]
  gsBenchmark benchmark;
  std::vector<index_t>  benchmarks, nthreads, msizes, vsizes;
  std::string fn;
  index_t nruns = 1,
    msizesmax = (index_t) std::sqrt(real_t(0.8) * sizeof(real_t)*gsSysInfo::getMemoryInBytes()),
    vsizesmax = (index_t)          (real_t(0.8) * sizeof(real_t)*gsSysInfo::getMemoryInBytes());

  gsCmdLine cmd("G+Smo performance benchmark.");
  cmd.printVersion();

  cmd.addInt("r", "runs", "Number of runs over which the results are averaged", nruns);
  cmd.addInt("M", "msizesmax", "Maximum number of unknowns in matrix/vector benchmarks (automated generation of sequence)", msizesmax);
  cmd.addInt("V", "vsizesmax", "Maximum number of unknowns in vector benchmarks (automated generation of sequence)", vsizesmax);
  cmd.addMultiInt("b", "benchmarks", "List of benchmarks to be run", benchmarks);
  cmd.addMultiInt("m", "msizes", "Number of unknowns in matrix/vector benchmarks", msizes);
  cmd.addMultiInt("t", "threads", "Number of OpenMP threads to be used for the benchmark", nthreads);
  cmd.addMultiInt("v", "vsizes", "Number of unknowns in vector benchmarks", vsizes);
  cmd.addString("o", "output", "Name of the output file", fn);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  //! [Parse command line]

  //! [Default configuration]
  // If empty fill with all benchmarks 1, ..., 5
  if (benchmarks.empty()) {
    for(index_t i=1; i<=9; ++i)
      benchmarks.push_back(i);
  }

  // If empty fill with 1, 2, 4, ..., maximum number of OpenMP threads
  if (nthreads.empty()) {
    for(index_t i=1; i<=omp_get_max_threads(); i*=2)
      nthreads.push_back(i);
  }

  // If empty fill with 10, 100, 1.000, ..., 80% of Sqrt(total system memory)
  if (msizes.empty()) {
    for(index_t i=10;;) {
      msizes.push_back(i);
      if (i<=std::min(msizesmax, std::numeric_limits<index_t>::max()) / 64)
        i*=8;
      else
        break;
    }
  }

  // If empty fill with 100, 1.000, 10.000, ... 80% of total system memory
  if (vsizes.empty()) {
    for(index_t i=100;;) {
      vsizes.push_back(i);
      if (i<=std::min(vsizesmax, std::numeric_limits<index_t>::max()) / 8)
        i*=8;
      else
        break;
    }
  }

  //! [Default configuration]

  //! [Execute benchmarks]
  for (auto bit=benchmarks.cbegin(); bit!=benchmarks.cend(); ++bit) {
    switch((index_t)(*bit)) {

    case (1): {
      // Benchmark: memcopy native C arrays
      CREATE_BENCHMARK(benchmark_c_array_memcopy, "memcopyCarray",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    case (2): {
      // Benchmark: memcopy gsVector
      CREATE_BENCHMARK(benchmark_eigen_memcopy, "memcopyEigen",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    case (3): {
      // Benchmark: dot-product native C array
      CREATE_BENCHMARK(benchmark_c_array_dotproduct, "dotproductCarray",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    case (4): {
      // Benchmark: dot-product gsVector
      CREATE_BENCHMARK(benchmark_eigen_dotproduct, "dotproductEigen",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    case (5): {
      // Benchmark: axpy native C array
      CREATE_BENCHMARK(benchmark_c_array_axpy, "axpyCarray",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    case (6): {
      // Benchmark: axpy gsVector
      CREATE_BENCHMARK(benchmark_eigen_axpy, "axpyEigen",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    case (7): {
      // Benchmark: dense matrix-vector multiplication native C array
      CREATE_BENCHMARK(benchmark_c_array_dense_matmul, "densematmulCarray",
                       msizes, metric::bandwidth_gb_sec);
      break;
    }

    case (8): {
      // Benchmark: dense matrix-vector multiplication gsMatrix/gsVector
      CREATE_BENCHMARK(benchmark_eigen_dense_matmul, "densematmulEigen",
                       msizes, metric::bandwidth_gb_sec);
      break;
    }

    case (9): {
      // Benchmark: visitor-based Poisson 2D assembly
      CREATE_BENCHMARK(benchmark_poisson2d_visitor, "assemblerVisitor",
                       vsizes, metric::bandwidth_gb_sec);
      break;
    }

    default:
      throw std::runtime_error("Invalid benchmark");
    }

  } // benchmark loop

  if (fn.empty())
    gsInfo << benchmark << "\n";
  else {
    std::ofstream file;
    file.open(fn);
    file << benchmark << "\n";
    file.close();
  }
  //! [Execute benchmarks]

  return  EXIT_SUCCESS;
}
