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
#include <gsIO/gsXml.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{
/**
   @brief Enumerator that defines the benchmark metrics.

   These definitions are used to control the output of the benchmark framework
*/
enum metric : uint64_t {
  speedup          = 1 <<  0,
  ratio            = 1 <<  1,
  bandwidth_kb_sec = 1 <<  2,
  bandwidth_mb_sec = 1 <<  3,
  bandwidth_gb_sec = 1 <<  4,
  bandwidth_tb_sec = 1 <<  5,
  perf_kflop_sec   = 1 <<  6,
  perf_mflop_sec   = 1 <<  7,
  perf_gflop_sec   = 1 <<  8,
  perf_tflop_sec   = 1 <<  9,
  runtime_sec      = 1 << 10
};

/**
   @brief Class that represents a single benchmark result

   A \a gsBenchmarkResult object is the most atomic unit of the
   benchmark framework. It represents the result of a single run for a
   fixed problem size and configuration and a fixed number of
   threads. A series of runs for different numbers of threads is
   collected in a \a gsBenchmarkResultSet object.
*/
class gsBenchmarkResult
{
public:
  int           threads;
  double        runtime;
  double        value;
  gismo::metric metric;
};

namespace internal
{
/// @brief Get a gsBenchmarkResult from XML data
template<>
class gsXml< gsBenchmarkResult >
{
private:
  gsXml() { }
  typedef gsBenchmarkResult Object;
public:
  GSXML_COMMON_FUNCTIONS(Object);
  static std::string tag ()  { return "BenchmarkResult"; }
  static std::string type () { return "BenchmarkResult"; }

  GSXML_GET_POINTER(Object);

  static void get_into (gsXmlNode * node, Object & obj)
  {
    gsXmlNode * child;

    child = node->first_node("threads");
    if (child != NULL) obj.threads = atoi(child->value());

    child = node->first_node("runtime");
    if (child != NULL) obj.runtime = atof(child->value());

    child = node->first_node("value");
    if (child != NULL) obj.value = atof(child->value());

    child = node->first_node("metric");
    if (child != NULL) obj.metric = (gismo::metric)atol(child->value());
  }

  static gsXmlNode * put (const Object & obj, gsXmlTree & data )
  {
    gsXmlNode * result = makeNode("BenchmarkResult", data);

    result->append_node( makeNode("threads", util::to_string(obj.threads), data) );
    result->append_node( makeNode("runtime", util::to_string(obj.runtime), data) );
    result->append_node( makeNode("value",   util::to_string(obj.value),   data) );
    result->append_node( makeNode("metric",  util::to_string(obj.metric),  data) );

    return result;
  }
};
} // namespace internal

/**
   @brief Class that represents a set of benchmark results

   A \a gsBenchmarkResultSet object holds a set of benchmark results
   (\a gsBenchmarkResult) for a fixed problem size and configuration
   but for different numbers of threads.
*/
class gsBenchmarkResultSet
{
public:
  /// \brief Default constructor
  gsBenchmarkResultSet() = default;

  /// \brief Constructor
  gsBenchmarkResultSet(const std::string& label,
                       const std::string& descr,
                       const std::vector<gsBenchmarkResult>& results)
    : label(label),
      descr(descr),
      results( give(std::vector<gsBenchmarkResult>(results)) ) {}

  /// \brief Constructor
  gsBenchmarkResultSet(const std::string& label,
                       const std::string& descr,
                       std::vector<gsBenchmarkResult>&& results)
    : label(label),
      descr(descr),
      results( give(results) ) {}

  /// \brief Returns the label
  const std::string& get_label() const
  { return label; }

  /// \brief Returns the descr
  const std::string& get_descr() const
  { return descr; }

  /// \brief Returns constant reference to the results
  const std::vector<gsBenchmarkResult>& get() const
  { return results; }

  /// \brief Returns non-constant reference to the results
  std::vector<gsBenchmarkResult>& get()
  { return results; }

  /// \brief Serializes the content to LaTeX TIKZ
  std::ostream &to_tikz(std::ostream &os) const;

  /// \brief Pretty-prints the content
  std::ostream &print(std::ostream &os) const;

private:
  std::string label, descr;
  std::vector<gsBenchmarkResult> results;
};

/// Print (as string) operator
inline std::ostream &operator<<(std::ostream &os, const gsBenchmarkResultSet& obj)
{ return obj.print(os); }

namespace internal
{
/// @brief Get a gsBenchmarkResultSet from XML data
template<>
class gsXml< gsBenchmarkResultSet >
{
private:
  gsXml() { }
  typedef gsBenchmarkResultSet Object;
public:
  GSXML_COMMON_FUNCTIONS(Object);
  static std::string tag ()  { return "BenchmarkResultSet"; }
  static std::string type () { return "BenchmarkResultSet"; }

  GSXML_GET_POINTER(Object);

  static void get_into (gsXmlNode * node, Object & obj)
  {
    gsXmlNode * child;
    std::string label, descr;

    child = node->first_node("label");
    if (child != NULL) label = child->value();

    child = node->first_node("descr");
    if (child != NULL) descr = child->value();

    std::vector<gsBenchmarkResult> results;

    child = node->first_node(gsXml< gsBenchmarkResult >::tag().c_str());
    for (; child; child = child->next_sibling() ) {
      gsBenchmarkResult result;
      gsXml< gsBenchmarkResult >::get_into(child, result);
      results.push_back( give(result) );
    }

    obj = gsBenchmarkResultSet(label, descr, give(results));
  }

  static gsXmlNode * put (const Object & obj, gsXmlTree & data )
  {
    gsXmlNode * results = makeNode("BenchmarkResultSet", data);

    results->append_node( makeNode("label", obj.get_label(), data) );
    results->append_node( makeNode("descr", obj.get_descr(), data) );

    for (const auto& it : obj.get()) {
      results->append_node( gsXml< gsBenchmarkResult >::put(it, data) );
    }

    return results;
  }
};
} // namespace internal

/**
   @brief Class that represents a collection of benchmark sets for a
   series of benchmark instances

   This struct can be used to hold a series of benchmark instances
   (i.e. a series of problem sizes and configurations)<
*/
class gsBenchmarkSet
{
public:
  /// \brief Default Constructor
  gsBenchmarkSet() = default;

  /// \brief Constructor
  gsBenchmarkSet(const std::string& label,
                 const std::string& descr,
                 const std::vector<gsBenchmarkResultSet>& results)
    : label(label),
      descr(descr),
      results( give(std::vector<gsBenchmarkResultSet>(results)) ) {}

  /// \brief Constructor
  gsBenchmarkSet(const std::string& label,
                 const std::string& descr,
                 std::vector<gsBenchmarkResultSet>&& results)
    : label(label),
      descr(descr),
      results( give(results) ) {}

  /// \brief Returns the label
  const std::string& get_label() const
  { return label; }

  /// \brief Returns the descr
  const std::string& get_descr() const
  { return descr; }

  /// \brief Returns constant reference to the result sets
  const std::vector<gsBenchmarkResultSet>& get() const
  { return results; }

  /// \brief Returns non-constant reference to the result sets
  std::vector<gsBenchmarkResultSet>& get()
  { return results; }

  /// \brief Serializes the content to LaTeX TIKZ
  std::ostream &to_tikz(std::ostream &os) const;

  /// \brief Pretty-prints the content
  std::ostream &print(std::ostream &os) const;

private:
  std::string label, descr;
  std::vector<gsBenchmarkResultSet> results;
};

/// Print (as string) operator
inline std::ostream &operator<<(std::ostream &os, const gsBenchmarkSet& obj)
{ return obj.print(os); }

namespace internal
{
/// @brief Get a gsBenchmarkSet from XML data
template<>
class gsXml< gsBenchmarkSet >
{
private:
  gsXml() { }
  typedef gsBenchmarkSet Object;
public:
  GSXML_COMMON_FUNCTIONS(Object);
  static std::string tag ()  { return "BenchmarkSet"; }
  static std::string type () { return "BenchmarkSet"; }

  GSXML_GET_POINTER(Object);

  static void get_into (gsXmlNode * node, Object & obj)
  {
    gsXmlNode * child;
    std::string label, descr;

    child = node->first_node("label");
    if (child != NULL) label = child->value();

    child = node->first_node("descr");
    if (child != NULL) descr = child->value();

    std::vector<gsBenchmarkResultSet> results;

    child = node->first_node(gsXml< gsBenchmarkResultSet >::tag().c_str());
    for (; child; child = child->next_sibling() ) {
      gsBenchmarkResultSet _results;
      gsXml< gsBenchmarkResultSet >::get_into(child, _results);
      results.push_back( give(_results) );
    }

    obj = gsBenchmarkSet(label, descr, give(results) );
  }

  static gsXmlNode * put (const Object & obj, gsXmlTree & data )
  {
    gsXmlNode * results = makeNode("BenchmarkSet", data);

    results->append_node( makeNode("label", obj.get_label(), data) );
    results->append_node( makeNode("descr", obj.get_descr(), data) );

    for (const auto& it : obj.get()) {
      results->append_node( gsXml< gsBenchmarkResultSet >::put(it, data) );
    }

    return results;
  }
};
} // namespace internal

/**
   @brief Class that represents a collection of benchmarks

   This is the top-level class of the benchmark framework and the only
   one that should be used by the user directly.
 */
class GISMO_EXPORT gsBenchmark
{
public:
  /// \brief Returns constant reference to the benchmarks
  const std::vector<gsBenchmarkSet>& get() const
  { return benchmarks; }

  /// \brief Returns non-constant reference to the benchmarks
  std::vector<gsBenchmarkSet>& get()
  { return benchmarks; }

  /// \brief Serializes the content to LaTeX TIKZ
  std::ostream &to_tikz(std::ostream &os) const;

  /// \brief Pretty-prints the content
  std::ostream &print(std::ostream &os) const;

  /// \brief Returns iterator to benchmark set
  const std::vector<gsBenchmarkSet>::const_iterator find(const std::string& label) const
  {
    for (auto it = benchmarks.cbegin(); it != benchmarks.cend(); ++it)
      if (it->get_label() == label)
        return it;
    return benchmarks.cend();
  }

  /// \brief Creates a new benchmark set, adds it to the benchmark and
  /// returns a pointer to the benchmark set to the calling routine
  template<typename Test, typename Iterator>
  const gsBenchmarkSet& create(const Iterator             & sizes,
                               const std::vector<index_t> & runs,
                               const std::vector<index_t> & threads,
                               const std::string          & extra_descr="")
  {
      //GISMO_ASSERT(sizes.size()==runs.size(), "Problem sizes and number of runs must have the same length");

    gsInfo << "[" << Test::label() << "] "
           << Test::descr()+extra_descr << "\n";

    std::vector<gsBenchmarkResultSet> results;
    char id('A');

    auto riter = runs.begin();
    for (const auto& it : sizes) {
      gsInfo << util::to_string(it) << "(" << *riter << ")"<< std::flush;
      try {
        Test test(it);
        uint64_t memsize = test.size();
        std::string meminfo;
        if (memsize<1024)
          meminfo = util::to_string(memsize)+" B";
        else if (memsize<1024*1024)
          meminfo = util::to_string(memsize/1024)+" KB";
        else if (memsize<1024*1024*1024)
          meminfo = util::to_string(memsize/(1024*1024))+" MB";
        else
          meminfo = util::to_string(memsize/(1024*1024*1024))+" GB";

        results.push_back( give(gsBenchmarkResultSet(Test::label()+std::string(1,id++), meminfo,
                                                     give(gsBenchmark::run(test, Test::metric(), threads, *riter++)))) );
      } catch(...) { gsInfo << "[failed!]"; }
      gsInfo << "...";
    }
    gsInfo << "\n";

    gsBenchmarkSet benchmark(Test::label(), Test::descr()+extra_descr, give(results) );
    benchmarks.push_back( give(benchmark) );
    return benchmarks.back();
  }

private:
  /// \brief Runs the benchmark instance \a benchmark for the
  /// specified number of \a threads and \a runs and returns an \a
  /// std::vector of \a gsBenchmarkResult that represent the
  /// respective benchmark results measured in the specified \a metric
  template<typename Test>
  static std::vector<gsBenchmarkResult>
  run(Test& test, gismo::metric metric, const std::vector<index_t>& threads, index_t runs)
  {
    std::vector<gsBenchmarkResult> results;
    gsStopwatch stopwatch;
    uint64_t result(0);
    real_t value, runtime;

    try {
      for (const auto& it : threads) {

        omp_set_num_threads(it);
        runtime = 0.0;

        stopwatch.restart();

        for (index_t run=0; run<runs; ++run) {
          result = test();
        }

        stopwatch.stop();
        runtime = stopwatch.elapsed()/(real_t)runs;

        switch(metric & ~gismo::metric::speedup & ~gismo::metric::ratio) {
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
        result.threads = static_cast<int>(it); // number of OpenMP threads
        result.runtime = runtime;              // averaged elapsed time in seconds
        result.value   = value;                // averaged benchmark value
        result.metric  = metric;               // benchmark metric
        results.push_back( give(result) );
      }
    } catch(...) {}

    // Convert to relative values (speedup relative to first entry)
    if (metric & gismo::metric::speedup) {
      runtime = results.front().runtime;
      value   = results.front().value;

      for (auto& it : results) {
        it.runtime = runtime / it.runtime;
        it.value   = value   / it.value;
      }
    }

    return results;
  }

private:
  std::vector<gsBenchmarkSet> benchmarks;
};

/// Print (as string) operator
inline std::ostream &operator<<(std::ostream &os, const gsBenchmark& obj)
{ return obj.print(os); }

namespace internal
{
/// @brief Get a gsBenchmark from XML data
template<>
class gsXml< gsBenchmark >
{
private:
  gsXml() { }
  typedef gsBenchmark Object;
public:
  GSXML_COMMON_FUNCTIONS(Object);
  static std::string tag ()  { return "Benchmark"; }
  static std::string type () { return "Benchmark"; }

  GSXML_GET_POINTER(Object);

  static void get_into (gsXmlNode * node, Object & obj)
  {
    gsXmlNode * child;

    child = node->first_node(gsXml< gsBenchmarkSet >::tag().c_str());
    for (; child; child = child->next_sibling() ) {
      gsBenchmarkSet benchmark;
      gsXml< gsBenchmarkSet >::get_into(child, benchmark);
      obj.get().push_back( give(benchmark) );
    }
  }

  static gsXmlNode * put (const Object & obj, gsXmlTree & data )
  {
    gsXmlNode * results = makeNode("Benchmark", data);

    for (const auto& it : obj.get()) {
      results->append_node( gsXml< gsBenchmarkSet >::put(it, data) );
    }

    return results;
  }
};
} // namespace internal

namespace util {

  /// \brief Returns the ratio of the two given benchmark result sets
  GISMO_EXPORT gsBenchmarkResultSet ratio(const std::string& label,
                                          const std::string& descr,
                                          const gsBenchmarkResultSet objA,
                                          const gsBenchmarkResultSet objB);

  /// \brief Returns the ratio of the two given benchmark sets
  GISMO_EXPORT gsBenchmarkSet ratio(const std::string& label,
                                    const std::string& descr,
                                    const gsBenchmarkSet objA,
                                    const gsBenchmarkSet objB);
} // namespace util

} // namespace gismo
