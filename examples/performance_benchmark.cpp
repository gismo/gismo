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

//! [Implement benchmark infrastructure]

/**
 * Benchmark metrics
 */
enum class benchmark_metric {
  bandwidth_kb_sec,
  bandwidth_mb_sec,
  bandwidth_gb_sec,
  bandwidth_tb_sec,  
  perf_kflop_sec,
  perf_mflop_sec,
  perf_gflop_sec,
  perf_tflop_sec,  
  runtime_sec,
};

/**
 *   Benchmark: driver function
 */
template<typename T>
std::vector< std::array<double,4> >
benchmark_driver(const std::vector<int>& nthreads, int nruns, T& benchmark, benchmark_metric metric)
{
  gsStopwatch stopwatch;
  std::size_t benchmark_result;
  double benchmark_metric, benchmark_runtime;

  std::vector< std::array<double,4> > results;

  try {
    for (auto it=nthreads.cbegin(); it!=nthreads.cend(); ++it) {

      omp_set_num_threads(*it);
      benchmark_runtime = 0.0;
      benchmark_metric = 0.0;

      for (int run=0; run<nruns; ++run) {
        stopwatch.restart();
        benchmark_result = benchmark();
        stopwatch.stop();
        benchmark_runtime += stopwatch.elapsed();
        
        switch(metric) {
        case benchmark_metric::bandwidth_kb_sec: case benchmark_metric::perf_kflop_sec:
          benchmark_metric += 1e-3*benchmark_result/stopwatch.elapsed();
          break;
        case benchmark_metric::bandwidth_mb_sec: case benchmark_metric::perf_mflop_sec:
          benchmark_metric += 1e-6*benchmark_result/stopwatch.elapsed();
          break;
        case benchmark_metric::bandwidth_gb_sec: case benchmark_metric::perf_gflop_sec:
          benchmark_metric += 1e-9*benchmark_result/stopwatch.elapsed();
          break;
        case benchmark_metric::bandwidth_tb_sec: case benchmark_metric::perf_tflop_sec:
          benchmark_metric += 1e-12*benchmark_result/stopwatch.elapsed();
          break;
        case benchmark_metric::runtime_sec:
          benchmark_metric += stopwatch.elapsed();
          break;
        default:
          throw std::runtime_error("Unsupported metric");
        }
        
      }

      results.push_back( { static_cast<double>(*it)        /* number of OpenMP threads */,
                           benchmark_runtime/(double)nruns /* averaged elapsed time in seconds */,
                           benchmark_metric/(double)nruns  /* averaged benchmark metric */,
                           (double)metric}                 /* benchmark metric */ );
    }
  } catch(...) {}

  return results;
}

/**
 * Benchmark LaTeX output
 */
class benchmark_latex
{
public:
  /**
   * Result set class
   */
  class result_set
  {
  public:
    result_set(const std::string& label,
               const std::string& title,
               const std::vector< std::array<double,4> >& results)
      : label(label),
        title(title),
        results(results)
    {
    }

    const std::string& get_label() const
    { return label; }
    
    const std::string& get_title() const
    { return title; }

    const std::vector< std::array<double,4> >& get_results() const
    { return results; }

    std::ostream &print(std::ostream &os) const
    {
      os << "\\pgfplotstableread[row sep=\\\\,col sep=&]{\n"
         << "threads & " << label << " \\\\\n";

      for (auto it=results.cbegin(); it!=results.cend(); ++it)
        os << (*it)[0] << "&" << (*it)[2] << "\\\\\n";

      os << "}\\data" << label << "\n";

      return os;
    }

  private:
    const std::string label, title;
    std::vector< std::array<double,4> > results;
  };

  /**
   * Benchmark set class
   */
  class benchmark_set
  {
  public:
    benchmark_set(const std::string& label,
                  const std::string& title)
      : id('A'),
        label(label),
        title(title)
    {}

    ~benchmark_set()
    {
      for (auto it=results.begin(); it!=results.end(); ++it)
        delete (*it);
    }

    void add_results(const std::string& label,
                     const std::string& title,
                     const std::vector< std::array<double,4> >& results)
    {
      this->results.emplace_back(new result_set(label+std::string(1,id++), title, results));
    }

    const std::string& get_label() const
    { return label; }
    
    const std::string& get_title() const
    { return title; }

    const std::vector<result_set*>& get_results() const
    { return results; }

    std::ostream &print(std::ostream &os) const
    {
      for (auto it=results.cbegin(); it!=results.cend(); ++it)
        (*it)->print(os);

      os << "\\begin{tikzpicture}\n"
         << "\\begin{axis}[\n"
         << "name=MyAxis,\n"
         << "width=\\textwidth,\n"
         << "height=.5\\textwidth,\n"
         << "legend pos=outer north east,\n"
        
         << "symbolic x coords={";
      
      for (auto it=(*results.cbegin())->get_results().cbegin();
           it!=(*results.cbegin())->get_results().cend(); ++it)
        os << (*it)[0] << (it!=(*results.cbegin())->get_results().cend()-1 ? "," : "");
      os << "},\n"
        
         << "xlabel={OpenMP threads},\n";
      
      switch((benchmark_metric)(*(*results.cbegin())->get_results().cbegin())[4]) {        
      case benchmark_metric::bandwidth_kb_sec:        
        os << "ylabel={Bandwidth in KB/s},\n";
        break;
      case benchmark_metric::bandwidth_mb_sec:        
        os << "ylabel={Bandwidth in MB/s},\n";
        break;
      case benchmark_metric::bandwidth_gb_sec:        
        os << "ylabel={Bandwidth in GB/s},\n";
        break;
      case benchmark_metric::bandwidth_tb_sec:        
        os << "ylabel={Bandwidth in TB/s},\n";
        break;
      case benchmark_metric::perf_kflop_sec:        
        os << "ylabel={Berformance in kFLOP/s},\n";
        break;
      case benchmark_metric::perf_mflop_sec:        
        os << "ylabel={Berformance in mFLOP/s},\n";
        break;
      case benchmark_metric::perf_gflop_sec:        
        os << "ylabel={Berformance in gFLOP/s},\n";
        break;
      case benchmark_metric::perf_tflop_sec:        
        os << "ylabel={Berformance in tFLOP/s},\n";
        break;
      case benchmark_metric::runtime_sec:
        os << "ylabel={Runtime in seconds},\n";
        break;
      default:
        throw std::runtime_error("Unsupported metric");
      }
      
      os << "title={" << title << "},\n"
         << "]";

      for (auto it=results.cbegin(); it!=results.cend(); ++it)
        os << "\\addplot table[x=threads,y="
           << (*it)->get_label()
           << "]{\\data"
           << (*it)->get_label()
           << "};\n";

      os << "\\legend{";
      for (auto it=results.cbegin(); it!=results.cend(); ++it)
        os << (*it)->get_title() << (it!=results.cend()-1 ? "," : "");
      os << "}\n"
        
         << "\\end{axis}\n"

         << "\\path let \\p1=(MyAxis.west), \\p2=(MyAxis.east) in "
         << "node[below right, align=left, text=black, text width=\\x2-\\x1]\n"
         << "at ($(MyAxis.south west)+(0,-30pt)$) {%\n"
         << "G+Smo " << gsCmdLine::getGismoVersion()
         << ", Eigen " << gsCmdLine::getEigenVersion()
         << " (" << gsCmdLine::getCompilerVersion()
         << ", " << gsCmdLine::getCppVersion()
         << ", " << gsCmdLine::getStdLibVersion()
         << (gsCmdLine::getExtraLibsVersion().empty()
             ? "), \n"
             : gsCmdLine::getExtraLibsVersion()+"), \n")

         << "CPU " << gsCmdLine::getCpuInfo() << ", "
         << "Memory " << gsCmdLine::getMemoryInfo()  << ", ";

      gsJITCompilerConfig jit; jit.load("config/jit.xml");
      std::string flags = jit.getFlags();
      os << "Compiler flags ";
      
      for (auto token=strtok(&flags[0], " "); token!=NULL; token=strtok(NULL, " ")) {
        if (token[0]=='-') {
          if (token[1]=='I' || token[1]=='L' || token[1]=='l' || token[1]=='W')
            continue;
          os << "\\verb!" << token << "! ";
        }
      }
      
      os << "};\n"
           << "\\end{tikzpicture}\n";

      return os;
    }

  private:
    char id;
    const std::string label,title;
    std::vector< result_set* > results;
  };

public:
  ~benchmark_latex()
    {
      for (auto it=benchmarks.begin(); it!=benchmarks.end(); ++it)
        delete (*it);
    }

  benchmark_set* add_benchmark(const std::string& label,
                               const std::string& title)
  {
    benchmarks.emplace_back(new benchmark_set(label, title));
    return benchmarks.back();
  }

  const std::vector< benchmark_set* >& get_benchmarks() const
  { return benchmarks; }

  std::ostream &print(std::ostream &os) const
  {
    os << "\\documentclass[tikz]{standalone}\n"
       << "\\usepackage{pgfplots}\n"
       << "\\usepackage{verbatim}\n"
       << "\\begin{document}\n"
       << "\\usetikzlibrary{calc}\n";

    for (auto it=benchmarks.cbegin(); it!=benchmarks.cend(); ++it)
      (*it)->print(os);

    os << "\\end{document}\n";
    return os;
  }

private:
  std::vector< benchmark_set* > benchmarks;
};

/// Print (as string) operator
std::ostream &operator<<(std::ostream &os, const benchmark_latex& obj)
{ return obj.print(os); }
//! [Implement benchmark infrastructure]


int main(int argc, char *argv[])
{
  //! [Parse command line]
  benchmark_latex latex;
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
  }
  //! [Default configuration]

  //! [Execute benchmarks]
  {
    auto bm = latex.add_benchmark("memcopy", "memory copy");
    {
      gsInfo << "=== Native C array memcopy\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_memcopy<real_t> benchmark(*it);
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("nativememcopy",
                          "native("+util::to_string((double)*it,0)+")",
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
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("eigenmemcopy",
                          "eigen("+util::to_string((double)*it,0)+")",
                          results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }

  {
    auto bm = latex.add_benchmark("dotprod", "dot-product");
    {
      gsInfo << "=== Native C array dot-product\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_dotproduct<real_t> benchmark(*it);
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("nativedotproduct",
                          "native("+util::to_string((double)*it,0)+")",
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
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("eigendotproduct",
                          "eigen("+util::to_string((double)*it,0)+")",
                          results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }

  {
    auto bm = latex.add_benchmark("axpy", "axpy");
    {
      gsInfo << "=== Native C array AXPY\n";
      for (auto it=vsizes.cbegin(); it!=vsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=vsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_axpy<real_t> benchmark(*it);
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("nativeaxpy",
                          "native("+util::to_string((double)*it,0)+")",
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
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("eigenaxpy",
                          "eigen("+util::to_string((double)*it,0)+")",
                          results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }

  {
    auto bm = latex.add_benchmark("densemvmul", "Dense matrix-vector multiply");
    {
      gsInfo << "=== Native C array dense matrix-vector multiplication\n";
      for (auto it=dsizes.cbegin(); it!=dsizes.cend(); ++it) {
        gsInfo << (*it) << (it!=dsizes.cend()-1 ? "." : "\n") << std::flush;
        try {
          benchmark_c_array_dense_matmul<real_t> benchmark(*it);
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("nativdensemvmul",
                          "native("+util::to_string((double)*it,0)+")",
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
          auto results = benchmark_driver(nthreads, nruns, benchmark, benchmark_metric::bandwidth_gb_sec);
          bm->add_results("eigenmvmul",
                          "eigen("+util::to_string((double)*it,0)+")",
                          results);
        } catch(...) { gsInfo << "failed!"; }
      }
    }
  }
  
  if (fn.empty())
    gsInfo << latex << "\n";
  else {
    std::ofstream file;
    file.open(fn);
    file << latex << "\n";
    file.close();
  }
  //! [Execute benchmarks]

  return  EXIT_SUCCESS;
}
