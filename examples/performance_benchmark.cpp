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

#include <array>

using namespace gismo;
//! [Include namespace]

/**
   Benchmark driver
*/
template<typename T>
std::vector< std::array<double,3> >
benchmark_driver(const std::vector<int>& nthreads, int nruns, T& benchmark)
{
  gsStopwatch stopwatch;
  std::size_t nbytes;
  double bandwidth;

  std::vector< std::array<double,3> > results;
  
  for (auto it=nthreads.cbegin(); it!=nthreads.cend(); ++it) {

    omp_set_num_threads(*it);
    bandwidth = 0.0;
      
    for (int run=0; run<nruns; ++run) {
      stopwatch.restart();
      nbytes = benchmark();
      stopwatch.stop();
      
      bandwidth += 1e-9*nbytes/stopwatch.elapsed();       
    }

    results.push_back( { static_cast<double>(*it), bandwidth/(double)nruns, stopwatch.elapsed() } );
  }

  return results;
}

/**
   Benchmark: native C array memcopy
*/
template<typename T>
class benchmark_c_array_memcopy
{
private:
  std::size_t n;

public:
  benchmark_c_array_memcopy(std::size_t n)
    : n(n)
  {}
  
  std::size_t operator()()
  {
    T* m_x = new T[n];
    T* m_y = new T[n];
    
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;
    
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_y[i] = m_x[i];
    
    // Needed to make sure the compiler does not eliminate this code block
    T tmp = m_y[n-1];
    GISMO_UNUSED(tmp);
    
    delete[] m_x;
    delete[] m_y;
    
    return sizeof(T) * 3 * n;
  }
};

/**
   Benchmark: native C array dot-product
*/
template<typename T>
class benchmark_c_array_dotproduct
{
private:
  std::size_t n;

public:
  benchmark_c_array_dotproduct(std::size_t n)
    : n(n)
  {}
  
  std::size_t operator()()
  {
    T* m_x = new T[n];
    T* m_y = new T[n];
    
#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_x[i] = (T)1.0;

#pragma omp parallel for simd
    for (std::size_t i=0; i<n; ++i)
      m_y[i] = (T)1.0;

    T sum = 0.0;
    
#pragma omp parallel for simd reduction(+:sum)
    for (std::size_t i=0; i<n; ++i)
      sum += m_x[i] * m_y[i];
    
    // Needed to make sure the compiler does not eliminate this code block
    T tmp = sum;
    GISMO_UNUSED(tmp);
    
    delete[] m_x;
    delete[] m_y;
    
    return sizeof(T) * 4 * n;
  }
};

/**
   Benchmark: eigen vector memcopy
*/
template<typename T>
class benchmark_eigen_vector_memcopy
{
private:
  std::size_t n;

public:
  benchmark_eigen_vector_memcopy(std::size_t n)
    : n(n)
  {}
  
  std::size_t operator()()
  {
    gsVector<T> x(n);
    gsVector<T> y(n);
    
    x.fill((T)1.0);
    y = x;
    
    // Needed to make sure the compiler does not eliminate this code block
    T tmp = y[n-1];
    GISMO_UNUSED(tmp);
    
    return sizeof(T) * 3 * n;
  }
};

/**
   Benchmark: eigen vector dot-product
*/
template<typename T>
class benchmark_eigen_vector_dotproduct
{
private:
  std::size_t n;

public:
  benchmark_eigen_vector_dotproduct(std::size_t n)
    : n(n)
  {}
  
  std::size_t operator()()
  {
    gsVector<T> x(n);
    gsVector<T> y(n);
    
    x.fill((T)1.0);
    y.fill((T)1.0);

    T sum = x.dot(y);
    
    // Needed to make sure the compiler does not eliminate this code block
    T tmp = sum;
    GISMO_UNUSED(tmp);
    
    return sizeof(T) * 4 * n;
  }
};
  
int main(int argc, char *argv[])
{
  //! [Parse command line]
  std::vector<int> nthreads;
  std::vector<double> bandwidths;
  int n=1000000000;
  int nruns=1;
  
  gsCmdLine cmd("G+Smo performance benchmark.");
  cmd.addMultiInt("t", "threads",
                  "Number of OpenMP threads to be used for the benchmark", nthreads);
  cmd.addInt("n", "nlength",
              "Number of unknowns in vector-type benchmarks", n);
  cmd.addInt("r", "runs",
             "Number of runs over which the results are averaged", nruns);
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
  if (nthreads.empty()) {
    for(int i=1; i<=omp_get_max_threads(); i*=2)
      nthreads.push_back(i);
  }

  {
    gsInfo << "=== Native C array memcopy ===\n";
    benchmark_c_array_memcopy<real_t> benchmark(n);
    auto results = benchmark_driver(nthreads, nruns, benchmark);
    for (auto it=results.cbegin(); it!=results.cend(); ++it)     
      gsInfo << "[OMP=" << (*it)[0] << "] " << (*it)[1] << " GB/s\n";
  }

  {
    gsInfo << "== gsVector memcopy ===\n";
    benchmark_eigen_vector_memcopy<real_t> benchmark(n);
    auto results = benchmark_driver(nthreads, nruns, benchmark);
    for (auto it=results.cbegin(); it!=results.cend(); ++it)     
      gsInfo << "[OMP=" << (*it)[0] << "] " << (*it)[1] << " GB/s\n";
  }

  {
    gsInfo << "=== Native C array dot-product ===\n";
    benchmark_c_array_dotproduct<real_t> benchmark(n);
    auto results = benchmark_driver(nthreads, nruns, benchmark);
    for (auto it=results.cbegin(); it!=results.cend(); ++it)     
      gsInfo << "[OMP=" << (*it)[0] << "] " << (*it)[1] << " GB/s\n";
  }

  {
    gsInfo << "== gsVector dot-product ===\n";
    benchmark_eigen_vector_dotproduct<real_t> benchmark(n);
    auto results = benchmark_driver(nthreads, nruns, benchmark);
    for (auto it=results.cbegin(); it!=results.cend(); ++it)     
      gsInfo << "[OMP=" << (*it)[0] << "] " << (*it)[1] << " GB/s\n";
  }
  
  return  EXIT_SUCCESS;
}
