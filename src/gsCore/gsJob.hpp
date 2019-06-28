/** @file gsJob.hpp

    @brief Provides declaration of the Job class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moeller
*/

#pragma once

namespace gismo {

namespace detail {

#ifdef GISMO_WITH_PYTHON

template<>
short
py_cast<short>(PyObject* obj)
{
  return PyLong_AsLong(obj);
}

template<>
unsigned short
py_cast<unsigned short>(PyObject* obj)
{
  return PyLong_AsUnsignedLong(obj);
}

template<>
int
py_cast<int>(PyObject* obj)
{
  return PyLong_AsLong(obj);
}

template<>
unsigned int
py_cast<unsigned int>(PyObject* obj)
{
  return PyLong_AsUnsignedLong(obj);
}

template<>
long
py_cast<long>(PyObject* obj)
{
  return PyLong_AsLong(obj);
}

template<>
unsigned long
py_cast<unsigned long>(PyObject* obj)
{
  return PyLong_AsUnsignedLong(obj);
}

template<>
long long
py_cast<long long>(PyObject* obj)
{
  return PyLong_AsLongLong(obj);
}

template<>
unsigned long long
py_cast<unsigned long long>(PyObject* obj)
{
  return PyLong_AsUnsignedLongLong(obj);
}

template<>
void*
py_cast<void*>(PyObject* obj)
{
  return PyBytes_AsVoidPtr(obj);
}
  
template<>
float
py_cast<float>(PyObject* obj)
{
  return PyFloat_AsDouble(obj);
}

template<>
double
py_cast<double>(PyObject* obj)
{
  return PyFloat_AsDouble(obj);
}

template<>
long double
py_cast<long double>(PyObject* obj)
{
  return PyFloat_AsDouble(obj);
}

template<>
std::string
py_cast<std::string>(PyObject* obj)
{
  return PyBytes_AsString(PyUnicode_AsUTF8String(obj));
}

template<>
PyObject*
py_cast<short>(short val)
{
  return PyLong_FromLong(val);
}

template<>
PyObject*
py_cast<unsigned short>(unsigned short val)
{
  return PyLong_FromUnsignedLong(val);
}

template<>
PyObject*
py_cast<int>(int val)
{
  return PyLong_FromLong(val);
}

template<>
PyObject*
py_cast<unsigned int>(unsigned int val)
{
  return PyLong_FromUnsignedLong(val);
}

template<>
PyObject*
py_cast<long>(long val)
{
  return PyLong_FromLong(val);
}

template<>
PyObject*
py_cast<unsigned long>(unsigned long val)
{
  return PyLong_FromUnsignedLong(val);
}

template<>
PyObject*
py_cast<long long>(long long val)
{
  return PyLong_FromLongLong(val);
}

template<>
PyObject*
py_cast<unsigned long long>(unsigned long long val)
{
  return PyLong_FromUnsignedLongLong(val);
}

template<>
PyObject*
py_cast<void*>(void* val)
{
  return PyLong_FromVoidPtr(val);
}

template<>
PyObject*
py_cast<float>(float val)
{
  return PyLong_FromDouble(val);
}

template<>
PyObject*
py_cast<double>(double val)
{
  return PyLong_FromDouble(val);
}

template<>
PyObject*
py_cast<long double>(long double val)
{
  return PyLong_FromDouble(val);
}

template<>
PyObject*
py_cast<std::string>(std::string val)
{
  return PyBytes_FromString(PyUnicode_FromUTF8String(val));
}
  
#endif

} // namespace detail
  
} // namespace gismo
