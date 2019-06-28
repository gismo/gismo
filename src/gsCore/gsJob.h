/** @file gsJob.h

    @brief Provides declaration of job execution class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moeller
*/

#pragma once

#include <chrono>
#include <functional>
#include <future>

namespace gismo {

namespace detail {
#ifdef GISMO_WITH_PYTHON

template<typename T>
T py_cast(PyObject* obj);

template<typename T>
PyObject* py_cast(T obj);

template<std::size_t index=0, typename...Args>
int py_cast_tuple_impl(const std::tuple<Args...>& val, PyObject* obj);

template<std::size_t index, typename...Args>
int py_cast_tuple_impl(const std::tuple<Args...>& val,
                       typename std::enable_if< (index < sizeof...(Args)-1), PyObject*>::type obj)
{
  PyObject* pValue = py_cast<>(std::get<index>(val));
  if (!pValue)
    return 1;
  
  // pValue reference gets stolen here
  PyTuple_SetItem(obj, index, pValue);
  return py_cast_tuple_impl<index+1>(val, obj);
}

template<std::size_t index, typename...Args>
int py_cast_tuple_impl(const std::tuple<Args...>& val,
                       typename std::enable_if< (index >= sizeof...(Args)-1), PyObject*>::type obj)
{
  PyObject* pValue = py_cast<>(std::get<index>(val));
  if (!pValue)
    return 1;
  
  // pValue reference gets stolen here
  PyTuple_SetItem(obj, index, pValue);
  return 0;
} 

template<typename...Args>
PyObject*
py_cast_tuple(const std::tuple<Args...>& val)
{
  PyObject* obj = PyTuple_New(sizeof...(Args));
  py_cast_tuple_impl<0>(val, obj);
  
  return obj;
}
  
#else

template<typename T, typename U>
T py_cast(U)
{
  throw "This feature requires -DGISMO_WITH_PYTHON";
}

template<typename T, typename U>
T py_cast_tuple(U)
{
  throw "This feature requires -DGISMO_WITH_PYTHON";
}
  
#endif

} // namespace detail

/**
     @brief Gismo job executor enumerator

     The Gismo job executor enumerator defines the supported job
   execution units
  */
enum class gsJob_enum
{
  Python,
  CXX
};

/**
   @brief Gismo job class

   The Gismo job class implements a job
*/
template<enum gsJob_enum>
class gsJob;

/**
    @brief Gismo job class specialization for Python executor

    The Gismo job class implements a job that can be executed
   within the Gismo stream class for the Python executor
  */
template<>
class gsJob<gsJob_enum::Python>
{
private:
  // Timing information
  std::chrono::time_point<std::chrono::high_resolution_clock> _timing_start;
  std::chrono::time_point<std::chrono::high_resolution_clock> _timing_stop;

#ifdef GISMO_WITH_PYTHON

private:
  // Python parameter and return value objects
  PyObject* _param = NULL;
  PyObject* _value = NULL;

  // Python module and global/local dictionary objects
  PyObject* _module = NULL;
  PyObject* _global = NULL;
  PyObject* _local  = NULL;

  // Handle to asynchronous operation
  std::future<void> _handle;

  // Python script
  const std::string _script;

  // Names of run, wait, and query methods
  const std::string _method_run;
  const std::string _method_wait;
  const std::string _method_query;

  // Launches the job
  void launch(const std::string& method, const gsJob* depend)
  {
    // Wait for dependent job (if any)
    if (depend != NULL)
      depend->_handle.wait();

    // Ensure that the current thread is ready to call the Python C API
    // regardless of the current state of Python, or of the global interpreter
    // lock.
    PyGILState_STATE gilState;
    gilState = PyGILState_Ensure();

    // Define function in the newly created module
    _value = PyRun_String(_script.c_str(), Py_file_input, _global, _local);

    // _value would be null if the Python syntax is wrong
    if (_value == NULL) {
      if (PyErr_Occurred()) {
        PyErr_Print();
      }
      throw "An error occured: Python syntax is wrong!";
    }

    // _value is the result of the executing code,
    // chuck it away because we've only declared a function
    Py_DECREF(_value);

    // Get a pointer to the function
    PyObject* _func = PyObject_GetAttrString(_module, method.c_str());

    // Double check we have actually found it and it is callable
    if (!_func || !PyCallable_Check(_func)) {
      if (PyErr_Occurred()) {
        PyErr_Print();
      }
      throw "An error occured: Cannot find find function \"" + method + "\"!";
    }

    // Call the function, passing to it the arguments
    _timing_start = std::chrono::high_resolution_clock::now();
    _value = PyObject_CallObject(_func, NULL);
    _timing_stop = std::chrono::high_resolution_clock::now();

    PyErr_Print();

    // Decrease reference counter
    Py_DECREF(_func);

    // Release any resources previously acquired
    PyGILState_Release(gilState);
  }

public:
  /// Default constructor deleted
  gsJob() = delete;

  /// Explicit job constructor
  template <typename...Args>
  explicit gsJob(std::tuple<Args...> py_param,
                 const std::string& script,
                 const std::string& method_run,
                 const std::string& method_wait,
                 const std::string& method_query,
                 PyObject* module,
                 PyObject* global,
                 PyObject* local,
                 const gsJob* depend = NULL)
    : _script(script)
    , _method_run(method_run)
    , _method_wait(method_wait)
    , _method_query(method_query)
    , _module(module)
    , _global(global)
    , _local(local)
  {
    // Convert parameters to Python objects
    _param = detail::py_cast_tuple<>(py_param);
    
    // Start asynchronous execution
    _handle = run(depend);
  }

  /// Destructor
  ~gsJob()
  {
    this->wait();

    // Ensure that the current thread is ready to call the Python C
    // API regardless of the current state of Python, or of the global
    // interpreter lock.
    PyGILState_STATE gilState;
    gilState = PyGILState_Ensure();
    
    // Decrease reference counter
    Py_DECREF(_value);

    // Release any resources previously acquired
    PyGILState_Release(gilState);
  }

  /// (Re-)runs the job
  std::future<void> run(const gsJob* depend = NULL)
  {
    return std::async(std::launch::async, &gsJob::launch, this, _method_run, depend);
  }

  /// Waits for job to complete and returns the result
  // nlohmann::json get() noexcept
  // {
  //   this->wait();
  //   return detail::py_cast<nlohmann::json>(_value);
  // }

  /// Waits for job to complete
  void wait() noexcept
  {
    if (_method_wait.empty()) {
      _handle.wait();

    } else {
      std::future<void> __handle;
      __handle =
        std::async(std::launch::async, &gsJob::launch, this, _method_run, this);
      __handle.wait();
    }
  }

  /// Returns true if job has completed, or false if not
  bool query() noexcept
  {

    if (_method_query.empty()) {
      return _handle.valid();

    } else {
      std::future<void> __handle;
      __handle = std::async(
        std::launch::async, &gsJob::launch, this, _method_query, this);
      __handle.wait();
      return detail::py_cast<bool>(_value);
    }
  }

#else
  
  /// Default constructor
  gsJob(const std::string& py_script,
        const std::string& py_method,
        PyObject* py_module,
        PyObject* py_global,
        PyObject* py_local,
        gsJob* py_depend = NULL)
  {
    throw "This feature requires -DGISMO_WITH_PYTHON";
  }
  
  /// Waits for job to complete and returns the result
  // const nlohmann::json& get() const noexcept
  // {
  //   throw "This feature requires -DGISMO_WITH_PYTHON";
  // }

  /// Waits for job to complete
  void wait() const noexcept
  {
    throw "This feature requires -DGISMO_WITH_PYTHON";
  }

  /// Returns true if job has completed, or false if not
  bool query() const noexcept
  {
    throw "This feature requires -DGISMO_WITH_PYTHON";

    return false;
  }

#endif

  /// Returns duration of the job
  template<class Rep = double, class Period = std::ratio<1>>
  const std::chrono::duration<Rep, Period> duration() const noexcept
  {
    return std::chrono::duration_cast<std::chrono::duration<Rep, Period>>(
      _timing_stop - _timing_start);
  }
};

/**
  @brief Gismo job class specialization for C++ executor

  The Gismo job class implements a job that can be executed
  within the Gismo stream class for the C++ executor
*/
template<>
class gsJob<gsJob_enum::CXX>
{
private:
  // Timing information
  std::chrono::time_point<std::chrono::high_resolution_clock> _timing_start;
  std::chrono::time_point<std::chrono::high_resolution_clock> _timing_stop;

  // Handle to asynchronous operation
  std::future<void> _handle;

  // Callback method
  const std::function<void(void)>& _method;

  // Launches the job
  void launch(const std::function<void(void)>& method, const gsJob* depend)
  {
    // Wait for dependent job (if any)
    if (depend != NULL)
      depend->_handle.wait();

    // Run method
    _timing_start = std::chrono::high_resolution_clock::now();
    method();
    _timing_stop = std::chrono::high_resolution_clock::now();
  }

public:
  /// Default constructor deleted
  gsJob() = delete;

  /// Explicit job constructur
  explicit gsJob(const std::function<void(void)>& method,
                const gsJob* depend = NULL)
    : _method(method)
  {
    // Start asynchronous execution
    _handle = run(depend);
  }

  /// Destructor
  ~gsJob() { this->wait(); }

  /// (Re-)runs the job
  std::future<void> run(const gsJob* depend = NULL)
  {
    return std::async(std::launch::async, &gsJob::launch, this, _method, depend);
  }

  /// Waits for job to complete
  void wait() const noexcept { _handle.wait(); }

  /// Returns true if job has completed, or false if not
  bool query() const noexcept { return _handle.valid(); }

  /// Returns duration of the job
  template<class Rep = double, class Period = std::ratio<1>>
  const std::chrono::duration<Rep, Period> duration() const noexcept
  {
    return std::chrono::duration_cast<std::chrono::duration<Rep, Period>>(
      _timing_stop - _timing_start);
  }
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsJob.hpp)
#endif
