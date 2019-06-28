/** @file gsStream.h

    @brief Provides declaration of the Stream class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moeller
*/

#pragma once

#include <deque>
#include <functional>

#ifdef GISMO_WITH_PYTHON
#include <Python.h>
#endif

#include <gsCore/gsJob.h>

namespace gismo {

/**
     @brief Gismo stream class

     The Gismo stream class implements a stream execution unit with
     support of synchronous and asynchronous execution of jobs
  */
template<enum gsJob_enum>
class gsStream;

/**
   @brief Gismo stream class specialization for the Python execution
   unit

   The Gismo stream class specialization for the Python
   execution unit supports synchronous and asynchronous execution of
   Python jobs using the C++11 std::async feature
*/
template<>
class gsStream<gsJob_enum::Python>
{
#ifdef GISMO_WITH_PYTHON
private:
  // Python module and global/local dictionary objects
  PyObject* _module = NULL;
  PyObject* _global = NULL;
  PyObject* _local  = NULL;

  // gsJob queue
  std::deque<gsJob<gsJob_enum::Python>*> _stream;

  // Global counter
  static std::size_t _counter;

  // Global Python thread state
  static PyThreadState* _state;

  // Shuts down Python interpreter gracefully
  static void shutdown()
  {
    PyEval_RestoreThread(_state);
    Py_Finalize();
  }

public:
  /// Default constructor
  gsStream()
  {
    // Check if Python interpreter has been started
    if (Py_IsInitialized() == 0) {
      // Initialize the Python interpreter
      Py_Initialize();

      // Create GIL (global interpreter lock)/enable threads
      PyEval_InitThreads();

      // Release the global interpreter lock and reset the thread state to NULL,
      // returning the previous thread state (which is not NULL).
      _state = PyEval_SaveThread();

      // Register Py_Finalize to be performed on exit
      // if (std::atexit(shutdown) != 0)
      //  throw "An error occured: Py_Finalize could not be registered for being
      //  called on exit!";
    }

    // Ensure that the current thread is ready to call the Python C API
    // regardless of the current state of Python, or of the global interpreter
    // lock.
    PyGILState_STATE gilState;
    gilState = PyGILState_Ensure();

    // Create global dictionary with builtins
    _global = PyDict_New();
    PyDict_SetItemString(_global, "__builtins__", PyEval_GetBuiltins());

    // Create a new module object
    std::string module = "gismo" + std::to_string(_counter++);
    _module = PyModule_New(module.c_str());
    PyModule_AddStringConstant(_module, "__file__", "");

    // Get the local dictionary object
    _local = PyModule_GetDict(_module);

    // Release any resources previously acquired
    PyGILState_Release(gilState);
  }

  /// Destructor
  ~gsStream()
  {
    while (!_stream.empty()) {
      delete _stream.front();
      _stream.pop_front();
    }

    // Ensure that the current thread is ready to call the Python C API
    // regardless of the current state of Python, or of the global interpreter
    // lock.
    PyGILState_STATE gilState;
    gilState = PyGILState_Ensure();

    // Decrease references
    Py_DECREF(_local);
    Py_DECREF(_global);
    Py_DECREF(_module);

    // Release any resources previously acquired
    PyGILState_Release(gilState);
  }

  /// Creates a new job and appends it to the current stream
  template <typename...Args>
  gsJob<gsJob_enum::Python>* run(std::tuple<Args...> py_param,
                                 const std::string&  py_script,
                                 const std::string&  py_method_run = "run",
                                 const std::string&  py_method_wait = "wait",
                                 const std::string&  py_method_query = "query")
  {
    gsJob<gsJob_enum::Python>* qjob;
    
    if (_stream.empty())
      qjob = new gsJob<gsJob_enum::Python>(py_param,
                                           py_script,
                                           py_method_run,
                                           py_method_wait,
                                           py_method_query,
                                           _module,
                                           _global,
                                           _local,
                                           NULL);
    else
      qjob = new gsJob<gsJob_enum::Python>(py_param,
                                           py_script,
                                           py_method_run,
                                           py_method_wait,
                                           py_method_query,
                                           _module,
                                           _global,
                                           _local,
                                           _stream.back());
    
    _stream.push_back(qjob);
    return qjob;
  }
  
  /// Returns the number of jobs in the stream
  std::size_t size() const { return _stream.size(); }
  
  /// Waits until all jobs in the current stream have completed
  void wail() const { _stream.back()->wait(); }
  
  /// Returns true if all jobs in the current stream have
  /// completed, or false if not
  bool query() const { return _stream.back()->query(); }
  
#else
  
  /// Default constructor
  gsStream() { throw "This feature requires -DGISMO_WITH_PYTHON"; }

  /// Creates a new job and appends it to the current stream
  template <typename...Args>
  gsJob<gsJob_enum::Python>* run(std::tuple<Args...> py_param,
                                 const std::string&  py_script,
                                 const std::string&  py_method = "run",
                                 const std::string&  py_method_wait = "wait",
                                 const std::string&  py_method_query = "query")
  {
    throw "This feature requires -DGISMO_WITH_PYTHON";
    
    gsJob<gsJob_enum::Python>* qjob;
    return qjob;
  }

  /// Returns the number of jobs in the stream
  std::size_t size() const
  {
    throw "This feature requires -DGISMO_WITH_PYTHON";
    
    return 0;
  }

  /// Waits until all jobs in the current stream have completed
  void wail() const { throw "This feature requires -DGISMO_WITH_PYTHON"; }

  /// Returns true if all operations in the current stream have
  /// completed, or false if not
  bool query() const
  {
    throw "This feature requires -DGISMO_WITH_PYTHON";

    return false;
  }
  
#endif // GISMO_WITH_PYTHON

  // Returns iterator to the beginning
  std::deque<gsJob<gsJob_enum::Python>*>::iterator begin() noexcept
  {
    return _stream.begin();
  }

  // Returns iterator to the end
  std::deque<gsJob<gsJob_enum::Python>*>::iterator end() noexcept
  {
    return _stream.end();
  }

  // Returns constant iterator to the beginning
  std::deque<gsJob<gsJob_enum::Python>*>::const_iterator cbegin() const noexcept
  {
    return _stream.cbegin();
  }

  // Returns constant iterator to the end
  std::deque<gsJob<gsJob_enum::Python>*>::const_iterator cend() const noexcept
  {
    return _stream.cend();
  }
};

// Initialization of static class members
std::size_t gsStream<gsJob_enum::Python>::_counter = 0;
PyThreadState* gsStream<gsJob_enum::Python>::_state = NULL;

/**
   @brief Gismo stream class specialization for the C++ execution unit

   The Gismo stream class specialization for the Python
   execution unit supports synchronous and asynchronous execution of
   Python jobs using the C++11 std::async feature
*/
template<>
class gsStream<gsJob_enum::CXX>
{
private:
  // gsJob queue
  std::deque<gsJob<gsJob_enum::CXX>*> _stream;

public:
  /// Default constructor
  gsStream() {}

  /// Destructor
  ~gsStream()
  {
    while (!_stream.empty()) {
      delete _stream.front();
      _stream.pop_front();
    }
  }

  /// Creates a new job and appends it to the current stream
  gsJob<gsJob_enum::CXX>* run(std::function<void(void)> exec)
  {
    gsJob<gsJob_enum::CXX>* qjob;

    if (_stream.empty())
      qjob = new gsJob<gsJob_enum::CXX>(exec);
    else
      qjob = new gsJob<gsJob_enum::CXX>(exec, _stream.back());

    _stream.push_back(qjob);
    return qjob;
  }

  /// Returns the number of jobs in the stream
  const std::size_t size() const { return _stream.size(); }

  /// Waits until all jobs in the current stream have completed
  void wail() const { _stream.back()->wait(); }

  /// Returns true if all jobs in the current stream have completed,
  /// or false if not
  bool query() const { return _stream.back()->query(); }

  /// Returns the duration of all jobs in the current stream
  template<class Rep = double, class Period = std::ratio<1>>
  const std::chrono::duration<Rep, Period> duration() const noexcept
  {
    std::chrono::duration<Rep, Period> _duration;

    for (auto it = this->cbegin(); it != this->cbegin(); it++)
      _duration += (*it)->duration();

    return _duration;
  }

  // Returns iterator to the beginning
  std::deque<gsJob<gsJob_enum::CXX>*>::iterator begin() noexcept
  {
    return _stream.begin();
  }

  // Returns iterator to the end
  std::deque<gsJob<gsJob_enum::CXX>*>::iterator end() noexcept
  {
    return _stream.end();
  }

  // Returns constant iterator to the beginning
  std::deque<gsJob<gsJob_enum::CXX>*>::const_iterator cbegin() const noexcept
  {
    return _stream.cbegin();
  }

  // Returns constant iterator to the end
  std::deque<gsJob<gsJob_enum::CXX>*>::const_iterator cend() const noexcept
  {
    return _stream.cend();
  }
};

} // namespace Gismo
