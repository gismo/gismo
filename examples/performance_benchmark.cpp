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
#include <gsCore/gsJITCompiler.h>

#include <array>
#include <iomanip>

using namespace gismo;
//! [Include namespace]


//! [Implement benchmarks]
/**
 * Benchmark: native C array memcopy
 */
template<typename T>
class benchmark_c_array_memcopy
{
private:
  std::size_t n;
  T *m_x, *m_y;

public:
  benchmark_c_array_memcopy(std::size_t n)
    : n(n), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
  }

  ~benchmark_c_array_memcopy()
  {
    delete[] m_x;
    delete[] m_y;
  }

  std::size_t operator()()
  {
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_y[i] = m_x[i];

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_y[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * 2 * n;
  }
};

/**
 * Benchmark: native C array dot-product
 */
template<typename T>
class benchmark_c_array_dotproduct
{
private:
  std::size_t n;
  T *m_x, *m_y;

public:
  benchmark_c_array_dotproduct(std::size_t n)
    : n(n), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;
  }

  ~benchmark_c_array_dotproduct()
  {
    delete[] m_x;
    delete[] m_y;
  }

  std::size_t operator()()
  {
    volatile T sum = 0.0;

#pragma omp parallel for simd reduction(+:sum)
    for (std::size_t i=0; i<n; ++i)
      sum += m_x[i] * m_y[i];

    GISMO_UNUSED(sum);

    return sizeof(T) * 2 * n;
  }
};

/**
 * Benchmark: native C array AXPY
 */
template<typename T>
class benchmark_c_array_axpy
{
private:
  std::size_t n;
  T *m_x, *m_y, *m_z;

public:
  benchmark_c_array_axpy(std::size_t n)
    : n(n), m_x(new T[n]), m_y(new T[n]), m_z(new T[n])
  {
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;
  }

  ~benchmark_c_array_axpy()
  {
    delete[] m_x;
    delete[] m_y;
  }

  std::size_t operator()()
  {
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_z[i] = (T)3.414 * m_x[i] + m_y[i];

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_z[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * 3 * n;
  }
};

/**
 * Benchmark: native C array dense matrix-vector multiplication
 */
template<typename T>
class benchmark_c_array_dense_matmul
{
private:
  std::size_t n;
  T *m_A, *m_x, *m_y;

public:
  benchmark_c_array_dense_matmul(std::size_t n)
    : n(n), m_A(new T[n*n]), m_x(new T[n]), m_y(new T[n])
  {
#pragma omp parallel for simd
    for (std::size_t i=0; i<n*n; ++i)
      m_A[i] = (T)1.0;

#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
  }

  ~benchmark_c_array_dense_matmul()
  {
    delete[] m_A;
    delete[] m_x;
    delete[] m_y;
  }

  std::size_t operator()()
  {
    for (std::size_t i=0; i<n; ++i) {
      T sum = (T)0.0;
#pragma omp parallel for simd reduction(+:sum)
      for (std::size_t j=0; j<n; ++j)
        sum += m_A[n*i+j] * m_x[j];
      m_y[i] = sum;
    }

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_y[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * (2*n*n + n);
  }
};

/**
 * Benchmark: Eigen vector memcopy
 */
template<typename T>
class benchmark_eigen_vector_memcopy
{
private:
  std::size_t n;
  gsVector<T> x,y;

public:
  benchmark_eigen_vector_memcopy(std::size_t n)
    : n(n), x(n), y(n)
  {
    x.fill((T)0.0);
  }

  std::size_t operator()()
  {
    y.noalias() = x;

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = y[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * 2 * n;
  }
};

/**
 * Benchmark: Eigen vector dot-product
 */
template<typename T>
class benchmark_eigen_vector_dotproduct
{
private:
  std::size_t n;
  gsVector<T> x, y;

public:
  benchmark_eigen_vector_dotproduct(std::size_t n)
    : n(n), x(n), y(n)
  {
    x.fill((T)0.0);
    y.fill((T)0.0);
  }

  std::size_t operator()()
  {
    volatile T sum = y.dot(x);
    GISMO_UNUSED(sum);

    return sizeof(T) * 2 * n;
  }
};

/**
 * Benchmark: Eigen vector AXPY
 */
template<typename T>
class benchmark_eigen_vector_axpy
{
private:
  std::size_t n;
  gsVector<T> x, y, z;

public:
  benchmark_eigen_vector_axpy(std::size_t n)
    : n(n), x(n), y(n), z(n)
  {
    x.fill((T)0.0);
    y.fill((T)0.0);
  }

  std::size_t operator()()
  {
    z.noalias() = (T)3.141*x + y;

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = z[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * 3 * n;
  }
};

/**
 * Benchmark: Eigen dense matrix-vector multiplication
 */
template<typename T>
class benchmark_eigen_vector_dense_matmul
{
private:
  std::size_t n;
  gsMatrix<T> A;
  gsVector<T> x, y;

public:
  benchmark_eigen_vector_dense_matmul(std::size_t n)
    : n(n), A(n,n), x(n), y(n)
  {
    A.fill(0.0);
    x.fill(0.0);
  }

  std::size_t operator()()
  {
    y.noalias() = A*x;

    // Needed to make sure the compiler does not eliminate this code block
    T tmp = y[n-1];
    GISMO_UNUSED(tmp);

    return sizeof(T) * (2*n*n + n);
  }
};
//! [Implement benchmarks]

int main(int argc, char *argv[])
{
  //! [Parse command line]
  gsBenchmark benchmark;
  std::vector<int> nthreads, ssizes, dsizes, vsizes;
  std::string fn;
  int nruns=1;

  gsCmdLine cmd("G+Smo performance benchmark.");
  cmd.printVersion();

  cmd.addInt("r", "runs", "Number of runs over which the results are averaged", nruns);
  cmd.addMultiInt("d", "dsizes", "Number of unknowns in dense matrix benchmarks", dsizes);
  cmd.addMultiInt("s", "ssizes", "Number of unknowns in sparse matrix benchmarks", ssizes);
  cmd.addMultiInt("t", "threads", "Number of OpenMP threads to be used for the benchmark", nthreads);
  cmd.addMultiInt("v", "vsizes", "Number of unknowns in vector benchmarks", vsizes);
  cmd.addString("o", "output", "Name of the output file", fn);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  //! [Parse command line]

  //! [Default configuration]
  // If empty fill with 1, 2, 4, ..., maximum number of OpenMP threads
  if (nthreads.empty()) {
    for(int i=1; i<=omp_get_max_threads(); i*=2)
      nthreads.push_back(i);
  }

  // If empty fill with 10, 100, 1.000, 10.000
  if (dsizes.empty()) {
    dsizes.push_back(1e1);
    dsizes.push_back(1e2);
    dsizes.push_back(1e3);
    dsizes.push_back(1e4);
  }
  
  // If empty fill with 100, 1.000, 10.000, 100.000, 1.000.000
  if (ssizes.empty()) {
    ssizes.push_back(1e2);
    ssizes.push_back(1e3);
    ssizes.push_back(1e4);
    ssizes.push_back(1e5);
    ssizes.push_back(1e6);
  }

  // If empty fill with 100, 1.000, 10.000, 100.000, 1.000.000
  if (vsizes.empty()) {
    vsizes.push_back(1e2);
    vsizes.push_back(1e3);
    vsizes.push_back(1e4);
    vsizes.push_back(1e5);
    vsizes.push_back(1e6);
    vsizes.push_back(1e7);
    vsizes.push_back(1e8);
  }
  //! [Default configuration]

  //! [Execute benchmarks]
  {
    auto bm = benchmark.add("memcopy", "memory copy");
    {
      gsInfo << "=== Native C array memcopy\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_memcopy<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("nativememcopy",
                  "native("+util::to_string(sizeof(double)*(double)*it / 1024 / 1024, 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }

    {
      gsInfo << "=== gsVector memcopy\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_eigen_vector_memcopy<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("eigenmemcopy",
                  "eigen("+util::to_string(sizeof(double)*(double)*it / 1024 / 1024, 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }

  {
    auto bm = benchmark.add("dotprod", "dot-product");
    {
      gsInfo << "=== Native C array dot-product\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_dotproduct<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("nativedotproduct",
                  "native("+util::to_string(sizeof(double)*(double)*it / 1024 / 1024, 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }

    {
      gsInfo << "=== gsVector dot-product\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_eigen_vector_dotproduct<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("eigendotproduct",
                  "eigen("+util::to_string(sizeof(double)*(double)*it / 1024 / 1024, 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }

  {
    auto bm = benchmark.add("axpy", "axpy");
    {
      gsInfo << "=== Native C array AXPY\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_axpy<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("nativeaxpy",
                  "native("+util::to_string(sizeof(double)*(double)*it / 1024 / 1024, 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }

    {
      gsInfo << "=== gsVector AXPY\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_eigen_vector_axpy<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("eigenaxpy",
                  "eigen("+util::to_string(sizeof(double)*(double)*it / 1024 / 1024, 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }

  {
    auto bm = benchmark.add("densemvmul", "Dense matrix-vector multiply");
    {
      gsInfo << "=== Native C array dense matrix-vector multiplication\n";
      for (auto it=dsizes.cbegin(); it!=dsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=dsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_dense_matmul<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("nativdensemvmul",
                  "native("+util::to_string(std::pow(sizeof(double)*(double)*it / 1024 / 1024, 2), 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }

    {
      gsInfo << "=== gsMatrix/gsVector dense matrix-vector multiplication\n";
      for (auto it=dsizes.cbegin(); it!=dsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=dsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_eigen_vector_dense_matmul<real_t> benchmark(*it);
          auto results = gsBenchmark::run(nthreads, nruns, benchmark, metric::bandwidth_gb_sec);
          bm->add("eigenmvmul",
                  "eigen("+util::to_string(std::pow(sizeof(double)*(double)*it / 1024 / 1024, 2), 0)+" MB)",
                  results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }
  
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
