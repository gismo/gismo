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
    speedup          = 1 << 0,
    bandwidth_kb_sec = 1 << 1,
    bandwidth_mb_sec = 1 << 2,
    bandwidth_gb_sec = 1 << 3,
    bandwidth_tb_sec = 1 << 4,
    perf_kflop_sec   = 1 << 5,
    perf_mflop_sec   = 1 << 6,
    perf_gflop_sec   = 1 << 7,
    perf_tflop_sec   = 1 << 8,
    runtime_sec      = 1 << 9
};
  
/**
   @brief Struct that represents a single benchmark result
*/
class gsBenchmarkResult
{
public:
  int           threads;
  double        runtime;
  double        value;
  gismo::metric metric;
};

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
  /// \brief Constructor
  gsBenchmarkResultSet(const std::string& label,
                       const std::string& descr,
                       const std::vector<gsBenchmarkResult>& results)
    : label(label),
      descr(descr),
      results(results) {}

  /// \brief Returns the label
  const std::string& get_label() const
  { return label; }

  /// \brief Returns the descr
  const std::string& get_descr() const
  { return descr; }

  /// \brief Returns constant reference to the results
  const std::vector<gsBenchmarkResult>& get() const
  { return results; }

  /// \brief Serializes the content to LaTeX TIKZ
  std::ostream &to_tikz(std::ostream &os) const;

  /// \brief Pretty-prints the content
  std::ostream &print(std::ostream &os) const;
  
private:
  const std::string label, descr;
  std::vector<gsBenchmarkResult> results;
};

/**
   @brief Struct that represents a collection of benchmark sets for a
   series of benchmark instances

   This struct can be used to hold a series of benchmark instances
   (i.e. a series of problem sizes and configurations)
*/
class gsBenchmarkSet
{
public:
  /// \brief Constructor
  gsBenchmarkSet(const std::string& _label,
                 const std::string& _descr)
    : id('A'),
      label(_label),
      descr(_descr)
  {}

  /// \brief Destructor
  ~gsBenchmarkSet()
  {
    for (auto it=results.begin(); it!=results.end(); ++it)
      delete (*it);
  }

  /// \brief Adds a benchmark to the benchmark set
  void add(const std::string& _label,
           const std::string& _descr,
           const std::vector<gsBenchmarkResult>& _results)
  {
    this->results.emplace_back(new gsBenchmarkResultSet(_label+std::string(1,id++),
                                                        _descr, _results));
  }

  /// \brief Returns the label
  const std::string& get_label() const
  { return label; }

  /// \brief Returns the descr
  const std::string& get_descr() const
  { return descr; }

  /// \brief Returns constant reference to the result sets
  const std::vector<gsBenchmarkResultSet*>& get() const
  { return results; }

  /// \brief Serializes the content to LaTeX TIKZ
  std::ostream &to_tikz(std::ostream &os) const;

  /// \brief Pretty-prints the content
  std::ostream &print(std::ostream &os) const;
  
private:
  char id;
  const std::string label,descr;
  std::vector<gsBenchmarkResultSet*> results;
};
  
/**
   @brief Class that collects all benchmark results
 */
class GISMO_EXPORT gsBenchmark
{
public:
  /// \brief Destructor
  ~gsBenchmark()
  {
    for (auto&& it : benchmarks)
      delete (it);
  }

  /// \brief Returns constant reference to the benchmarks
  const std::vector<gsBenchmarkSet*>& get() const
  { return benchmarks; }

  /// \brief Serializes the content to LaTeX TIKZ
  std::ostream &to_tikz(std::ostream &os) const;

  /// \brief Pretty-prints the content
  std::ostream &print(std::ostream &os) const;

  /// \brief Returns iterator to benchmark set
  const gsBenchmarkSet* find(const std::string& label) const
  {
    for (const auto& it : benchmarks)
      if (it->get_label() == label)
        return it;
    return nullptr;
  }

  /// \brief Creates a new benchmark set, adds it to the benchmark and
  /// returns a pointer to the benchmark set to the calling routine
  gsBenchmarkSet* create(const std::string& _label,
                         const std::string& _descr)
  {
    benchmarks.emplace_back(new gsBenchmarkSet(_label, _descr));
    return benchmarks.back();
  }
  
  /// \brief Creates a new benchmark set, adds it to the benchmark and
  /// returns a pointer to the benchmark set to the calling routine
  template<typename Test, typename Iterator>
  gsBenchmarkSet* create(const Iterator             & sizes,
                         const std::vector<index_t> & runs,
                         const std::vector<index_t> & threads,
                         const std::string          & extra_name="")
  {
    GISMO_ASSERT(sizes.size()==runs.size(), "Problem sizes and number of runs must have the same length");
    
    auto benchmark = create(Test::label(), Test::name()+extra_name);    
    gsInfo << "[" << benchmark->get_label() << "] " << benchmark->get_descr() << "\n";
    
    auto riter = runs.begin();
    for (auto it : sizes) {
      gsInfo << util::to_string(it) << "(" << *riter << ")"<< std::flush;
      try {
        Test test(it);
        auto results = gsBenchmark::run_test(test, Test::metric(), threads, *riter++);
        std::string meminfo;
        uint64_t memsize = test.size();
        if (memsize<1024)
          meminfo = util::to_string(memsize)+" B";
        else if (memsize<1024*1024)
          meminfo = util::to_string(memsize/1024)+" KB";
        else if (memsize<1024*1024*1024)
          meminfo = util::to_string(memsize/(1024*1024))+" MB";
        else
          meminfo = util::to_string(memsize/(1024*1024*1024))+" GB";
        benchmark->add(Test::label(), meminfo, results);
      } catch(...) { gsInfo << "[failed!]"; }
      gsInfo << "...";
    }
    gsInfo << "\n";
    return benchmark;
  }

    /// \brief Creates a new benchmark set, adds it to the benchmark and
  /// returns a pointer to the benchmark set to the calling routine
  template<typename Test, typename... T>
  gsBenchmarkSet* create(const util::zip_helper<T...>& sizes,
                         const std::vector<index_t>  & runs,
                         const std::vector<index_t>  & threads,
                         const std::string           & extra_name="")
  {
    GISMO_ASSERT(sizes.size()==runs.size(), "Problem sizes and number of runs must have the same length");
       
    auto benchmark = create(Test::label(), Test::name()+extra_name);
    gsInfo << "[" << benchmark->get_label() << "] " << benchmark->get_descr() << "\n";
    
    auto riter = runs.begin();
    for (auto it : sizes) {
      gsInfo << util::to_string(it) << "(" << *riter << ")"<< std::flush;
      try {
        Test test(it);
        auto results = gsBenchmark::run_test(test, Test::metric(), threads, *riter++);
        std::string meminfo;
        uint64_t memsize = test.size();
        if (memsize<1024)
          meminfo = util::to_string(memsize)+" B";
        else if (memsize<1024*1024)
          meminfo = util::to_string(memsize/1024)+" KB";
        else if (memsize<1024*1024*1024)
          meminfo = util::to_string(memsize/(1024*1024))+" MB";
        else
          meminfo = util::to_string(memsize/(1024*1024*1024))+" GB";
        benchmark->add(Test::label(), meminfo, results);
      } catch(...) { gsInfo << "[failed!]"; }
      gsInfo << "...";
    }
    gsInfo << "\n";
    return benchmark;
  }

private:
  /// \brief Runs the benchmark instance \a benchmark for the
  /// specified number of \a threads and \a runs and returns an \a
  /// std::vector of \a gsBenchmarkResult that represent the
  /// respective benchmark results measured in the specified \a metric
  template<typename Test>
  static std::vector<gsBenchmarkResult>
  run_test(Test& test, gismo::metric metric,
           const std::vector<index_t>& threads, index_t runs)
  {
    gsStopwatch stopwatch;
    uint64_t result(0);
    real_t value, runtime;
    
    std::vector<gsBenchmarkResult> results;

    try {
      for (auto it=threads.cbegin(); it!=threads.cend(); ++it) {

        omp_set_num_threads(*it);
        runtime = 0.0;

        stopwatch.restart();

        for (index_t run=0; run<runs; ++run) {
          result = test();
        }

        stopwatch.stop();
        runtime = stopwatch.elapsed()/(real_t)runs;

        switch(metric & ~gismo::metric::speedup) {
        case gismo::metric::bandwidth_kb_sec: case gismo::metric::perf_kflop_sec:
          value = 1e-03 * result / runtime;
          break;
        case gismo::metric::bandwidth_mb_sec: case gismo::metric::perf_mflop_sec:
          value = 1e-06 * result / runtime;
          break;
        case gismo::metric::bandwidth_gb_sec: case gismo::metric::perf_gflop_sec:
          value = 1e-09 * result / runtime;
          break;
        case gismo::metric::bandwidth_tb_sec: case gismo::metric::perf_tflop_sec:
          value = 1e-12 * result / runtime;
          break;
        case gismo::metric::runtime_sec:
          value = runtime;
          break;
        default:
          GISMO_ERROR("Unsupported metric");
        }

        gsBenchmarkResult result;
        result.threads = static_cast<int>(*it); // number of OpenMP threads
        result.runtime = runtime;               // averaged elapsed time in seconds
        result.value   = value;                 // averaged benchmark value
        result.metric  = metric;                // benchmark metric
        results.push_back( give(result) );
      }
    } catch(...) {}

    // Convert to relative values (speedup relative to first entry)
    if (metric & gismo::metric::speedup) {
      runtime = results.front().runtime;
      value   = results.front().value;

      for (auto &it : results) {
        it.runtime = runtime / it.runtime;
        it.value   = value   / it.value;
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
