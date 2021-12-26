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
#include <gsUtils/gsStopwatch.h>

namespace gismo
{
/**
   @brief Enumerator that defines the benchmark metrics.
   
   These definitions are used to control the output of the benchmark framework
*/
enum metric : uint64_t {
    speedup          = 1<<0,
    bandwidth_kb_sec = 1<<1,
    bandwidth_mb_sec = 1<<2,
    bandwidth_gb_sec = 1<<3,
    bandwidth_tb_sec = 1<<4,
    perf_kflop_sec   = 1<<5,
    perf_mflop_sec   = 1<<6,
    perf_gflop_sec   = 1<<7,
    perf_tflop_sec   = 1<<8,
    runtime_sec      = 1<<9
};
  
/**
   @brief Struct that represents a single benchmark result
*/
class Result1
{
public:
  int threads;
  gismo::metric metric;
  double value;
  double runtime;   
};
  typedef std::array<real_t,4> Result;

/**
   @brief Struct that represents a collection of benchmark results for
   a single benchmark instance

   This struct can be used to hold a series of results of a single
   benchmark instance (i.e. fixed problem size and problem
   configuration) for different numbers of threads.
*/
class gsBenchmarkResultSet
{
public:
  gsBenchmarkResultSet(const std::string& label,
                       const std::string& title,
                       const std::vector<Result>& results)
    : label(label),
      title(title),
      results(results) {}
  
  const std::string& get_label() const
  { return label; }
  
  const std::string& get_title() const
  { return title; }
  
  const std::vector<Result>& get() const
  { return results; }

  std::ostream &to_tikz(std::ostream &os) const;
  std::ostream &print(std::ostream &os) const;
  
private:
  const std::string label, title;
  std::vector<Result> results;
};

/**
   @brief Struct that represents a collection of benchmark sets for a
   series of benchmark instance

   This struct can be used to hold a series of benchmark instances
   (i.e. a series of problem sizes and configurations)
*/
class gsBenchmarkSet
{
public:
  gsBenchmarkSet(const std::string& _label,
                 const std::string& _title)
    : id('A'),
      label(_label),
      title(_title)
  {}
  
  ~gsBenchmarkSet()
  {
    for (auto it=results.begin(); it!=results.end(); ++it)
      delete (*it);
  }
  
  void add(const std::string& _label,
           const std::string& _title,
           const std::vector<Result>& _results)
  {
    this->results.emplace_back(new gsBenchmarkResultSet(_label+std::string(1,id++),
                                                        _title, _results));
  }
  
  const std::string& get_label() const
  { return label; }
  
  const std::string& get_title() const
  { return title; }
  
  const std::vector<gsBenchmarkResultSet*>& get() const
  { return results; }

  std::ostream &to_tikz(std::ostream &os) const;
  std::ostream &print(std::ostream &os) const;
  
private:
  char id;
  const std::string label,title;
  std::vector<gsBenchmarkResultSet*> results;
};
  
/**
   @brief Class that collects all benchmark results
 */
class GISMO_EXPORT gsBenchmark
{
public:
  ~gsBenchmark()
    {
      for (auto it=benchmarks.begin(); it!=benchmarks.end(); ++it)
        delete (*it);
    }

  gsBenchmarkSet* add(const std::string& _label,
                      const std::string& _title)
  {
    benchmarks.emplace_back(new gsBenchmarkSet(_label, _title));
    return benchmarks.back();
  }

  const std::vector<gsBenchmarkSet*>& get() const
  { return benchmarks; }

  std::ostream &to_tikz(std::ostream &os) const;
  std::ostream &print(std::ostream &os) const;

  template<typename T>
  static std::vector<Result>
  run(const std::vector<index_t>& nthreads, index_t nruns, T& benchmark, gismo::metric metric)
{
    gsStopwatch stopwatch;
    uint64_t benchmark_result(0);
    real_t benchmark_metric, benchmark_runtime;

    std::vector<Result> results;

    try {
      for (auto it=nthreads.cbegin(); it!=nthreads.cend(); ++it) {

        omp_set_num_threads(*it);
        benchmark_runtime = 0.0;

        stopwatch.restart();

        for (index_t run=0; run<nruns; ++run) {
          benchmark_result = benchmark();
        }

        stopwatch.stop();
        benchmark_runtime = stopwatch.elapsed()/(real_t)nruns;

        switch(metric & ~gismo::metric::speedup) {
        case gismo::metric::bandwidth_kb_sec: case gismo::metric::perf_kflop_sec:
          benchmark_metric = 1e-3*benchmark_result/benchmark_runtime;
          break;
        case gismo::metric::bandwidth_mb_sec: case gismo::metric::perf_mflop_sec:
          benchmark_metric = 1e-6*benchmark_result/benchmark_runtime;
          break;
        case gismo::metric::bandwidth_gb_sec: case gismo::metric::perf_gflop_sec:
          benchmark_metric = 1e-9*benchmark_result/benchmark_runtime;
          break;
        case gismo::metric::bandwidth_tb_sec: case gismo::metric::perf_tflop_sec:
          benchmark_metric = 1e-12*benchmark_result/benchmark_runtime;
          break;
        case gismo::metric::runtime_sec:
          benchmark_metric = benchmark_runtime;
          break;
        default:
          GISMO_ERROR("Unsupported metric");
        }

        Result res;
        res[0]= static_cast<real_t>(*it); // number of OpenMP threads
        res[1]= benchmark_runtime;        // averaged elapsed time in seconds
        res[2]= benchmark_metric;         // averaged benchmark metric
        res[3]= (real_t)metric;           // benchmark metric
        results.push_back( give(res) );
      }
    } catch(...) {}

    // Convert to relative values (speedup relative to first entry)
    if (metric & gismo::metric::speedup) {
      benchmark_runtime = results.front().at(1);
      benchmark_metric  = results.front().at(2);

      for (auto &it : results) {
        it.at(1) = benchmark_runtime / it.at(1);
        it.at(2) = benchmark_metric  / it.at(2);
      }
    }

    return results;
  }

private:
  std::vector<gsBenchmarkSet*> benchmarks;
};

/// Print (as string) operator
inline std::ostream &operator<<(std::ostream &os, const gsBenchmark& obj)
{ return obj.print(os); }

} // namespace gismo
