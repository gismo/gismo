/** @file gsBenchmark.h

    @brief Provides a generic benchmarking framework.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsSysInfo.h>
#include <gsCore/gsJITCompiler.h>
#include <gsUtils/gsStopwatch.h>

#include <array>
#include <vector>

namespace gismo
{
  
/**
 * Benchmark metrics
 */
enum metric {
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
 * Benchmark result
 */
typedef std::array<double,4> gsBenchmarkResult;
  
/**
 *   Benchmark: driver function
 */


/**
 * Benchmark class
 */
class gsBenchmark
{
public:
  /**
   * Benchmark result set class
   */
  class gsBenchmarkResultSet
  {
  public:
    gsBenchmarkResultSet(const std::string& label,
                         const std::string& title,
                         const std::vector<gsBenchmarkResult>& results)
      : label(label),
        title(title),
        results(results)
    {
    }

    const std::string& get_label() const
    { return label; }
    
    const std::string& get_title() const
    { return title; }

    const std::vector<gsBenchmarkResult>& get() const
    { return results; }

    std::ostream &print(std::ostream &os) const;

  private:
    const std::string label, title;
    std::vector<gsBenchmarkResult> results;
  };

  /**
   * Benchmark set class
   */
  class gsBenchmarkSet
  {
  public:
    gsBenchmarkSet(const std::string& label,
                   const std::string& title)
      : id('A'),
        label(label),
        title(title)
    {}

    ~gsBenchmarkSet()
    {
      for (auto it=results.begin(); it!=results.end(); ++it)
        delete (*it);
    }

    void add(const std::string& label,
             const std::string& title,
             const std::vector<gsBenchmarkResult>& results)
    {
      this->results.emplace_back(new gsBenchmarkResultSet(label+std::string(1,id++),
                                                          title, results));
    }

    const std::string& get_label() const
    { return label; }
    
    const std::string& get_title() const
    { return title; }

    const std::vector<gsBenchmarkResultSet*>& get() const
    { return results; }

    std::ostream &print(std::ostream &os) const;

  private:
    char id;
    const std::string label,title;
    std::vector<gsBenchmarkResultSet*> results;
  };

public:
  ~gsBenchmark()
    {
      for (auto it=benchmarks.begin(); it!=benchmarks.end(); ++it)
        delete (*it);
    }

  gsBenchmarkSet* add(const std::string& label,
                      const std::string& title)
  {
    benchmarks.emplace_back(new gsBenchmarkSet(label, title));
    return benchmarks.back();
  }

  const std::vector<gsBenchmarkSet*>& get() const
  { return benchmarks; }

  std::ostream &print(std::ostream &os) const;

  template<typename T>
  static std::vector<gsBenchmarkResult>
  run(const std::vector<int>& nthreads, int nruns, T& benchmark, metric metric)
  {
    gsStopwatch stopwatch;
    std::size_t benchmark_result;
    double benchmark_metric, benchmark_runtime;
    
    std::vector<gsBenchmarkResult> results;
    
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
          case metric::bandwidth_kb_sec: case metric::perf_kflop_sec:
            benchmark_metric += 1e-3*benchmark_result/stopwatch.elapsed();
            break;
          case metric::bandwidth_mb_sec: case metric::perf_mflop_sec:
            benchmark_metric += 1e-6*benchmark_result/stopwatch.elapsed();
            break;
          case metric::bandwidth_gb_sec: case metric::perf_gflop_sec:
            benchmark_metric += 1e-9*benchmark_result/stopwatch.elapsed();
            break;
          case metric::bandwidth_tb_sec: case metric::perf_tflop_sec:
            benchmark_metric += 1e-12*benchmark_result/stopwatch.elapsed();
            break;
          case metric::runtime_sec:
            benchmark_metric += stopwatch.elapsed();
            break;
          default:
            throw std::runtime_error("Unsupported metric");
          }
          
        }
        
        results.push_back(
                          { static_cast<double>(*it)        /* number of OpenMP threads */,
                            benchmark_runtime/(double)nruns /* averaged elapsed time in seconds */,
                            benchmark_metric/(double)nruns  /* averaged benchmark metric */,
                            (double)metric                 /* benchmark metric */
                          });
      }
    } catch(...) {}
    
    return results;
  }
  
private:
  std::vector<gsBenchmarkSet*> benchmarks;
};

/// Print (as string) operator
std::ostream &operator<<(std::ostream &os, const gsBenchmark& obj)
{ return obj.print(os); }

} // namespace gismo
