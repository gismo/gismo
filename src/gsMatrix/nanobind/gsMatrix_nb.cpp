/** @file gsMatrix_nb.cpp

    @brief NanoBind bindings for gsVector, gsMatrix, gsSparseMatrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_NANOBIND

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>

#include <gismo.h>

namespace nb = nanobind;

#if NB_VERSION_MAJOR >= 2
constexpr ssize_t nb_any = -1;
#else
constexpr size_t nb_any = nb::any;
#endif

template<class T>
void bind_gsVector(nb::module_& m, const char* name)
{
    using Vec = gismo::gsVector<T>;
    nb::class_<Vec>(m, name)
        .def(nb::init<>())
        .def(nb::init<index_t>())
        .def("size", [](const Vec& v) { return v.size(); })
        .def("rows", [](const Vec& v) { return v.rows(); })
        .def("setZero", [](Vec& v) { v.setZero(); return &v; }, nb::rv_policy::reference)
        .def("setOnes", [](Vec& v) { v.setOnes(); return &v; }, nb::rv_policy::reference)
        .def("setConstant", [](Vec& v, T val) { v.setConstant(val); return &v; }, nb::rv_policy::reference)
        .def("__getitem__", [](const Vec& v, index_t i) -> T {
            if (i < 0) i += v.size();
            if (i < 0 || i >= v.size()) throw nb::index_error();
            return v[i];
        })
        .def("__setitem__", [](Vec& v, index_t i, T val) {
            if (i < 0) i += v.size();
            if (i < 0 || i >= v.size()) throw nb::index_error();
            v[i] = val;
        })
        .def("__len__", [](const Vec& v) { return v.size(); })
        .def("__repr__", [](const Vec& v) {
            std::ostringstream s;
            s << "gsVector(" << v.size() << ")\n" << v.transpose();
            return s.str();
        })
        .def("tolist", [](const Vec& v) {
            nb::list l;
            for (index_t i = 0; i < v.size(); ++i)
                l.append(v[i]);
            return l;
        })
        .def("numpy", [](const Vec& v) {
            nb::object np = nb::module_::import_("numpy");
            nb::object arr = np.attr("empty")(v.size(), nb::arg("dtype") = "float64");
            auto buf = nb::cast<nb::ndarray<T, nb::shape<nb_any>>>(arr);
            std::memcpy(buf.data(), v.data(), v.size() * sizeof(T));
            return arr;
        })
        ;
}

template<class T>
void bind_gsMatrix(nb::module_& m, const char* name)
{
    using Mat = gismo::gsMatrix<T>;
    nb::class_<Mat>(m, name)
        .def(nb::init<>())
        .def(nb::init<int, int>())
        .def("rows", [](const Mat& m) { return m.rows(); })
        .def("cols", [](const Mat& m) { return m.cols(); })
        .def("size", [](const Mat& m) { return m.size(); })
        .def("dim", &Mat::dim)
        .def("setZero", [](Mat& m) { m.setZero(); return &m; }, nb::rv_policy::reference)
        .def("setOnes", [](Mat& m) { m.setOnes(); return &m; }, nb::rv_policy::reference)
        .def("setConstant", [](Mat& m, T val) { m.setConstant(val); return &m; }, nb::rv_policy::reference)
        .def("setIdentity", [](Mat& m) { m.setIdentity(); return &m; }, nb::rv_policy::reference)
        .def("setIdentity", [](Mat& m, index_t rows, index_t cols) { m.setIdentity(rows, cols); return &m; }, nb::rv_policy::reference)
        .def("resize", [](Mat& m, index_t rows, index_t cols) { m.resize(rows, cols); })
        .def("conservativeResize", [](Mat& m, index_t rows, index_t cols) { m.conservativeResize(rows, cols); })
        .def("clear", &Mat::clear)
        .def("transpose", [](const Mat& m) { return Mat(m.transpose()); })
        .def("norm", [](const Mat& m) { return m.norm(); })
        .def("squaredNorm", [](const Mat& m) { return m.squaredNorm(); })
        .def("sum", [](const Mat& m) { return m.sum(); })
        .def("minCoeff", [](const Mat& m) { return m.minCoeff(); })
        .def("maxCoeff", [](const Mat& m) { return m.maxCoeff(); })
        .def("trace", [](const Mat& m) { return m.trace(); })
        .def("determinant", [](const Mat& m) { return m.determinant(); })
        .def("inverse", [](const Mat& m) { return Mat(m.inverse()); })
        .def("__call__", [](const Mat& m, index_t i, index_t j) -> T {
            return m(i, j);
        })
        .def("set", [](Mat& m, index_t i, index_t j, T val) {
            m(i, j) = val;
        })
        .def("__repr__", [](const Mat& m) {
            std::ostringstream s;
            s << "gsMatrix(" << m.rows() << " x " << m.cols() << ")\n" << m;
            return s.str();
        })
        .def("numpy", [](const Mat& m) {
            nb::object np = nb::module_::import_("numpy");
            nb::object arr = np.attr("empty")(
                nb::make_tuple(m.rows(), m.cols()),
                nb::arg("dtype") = "float64",
                nb::arg("order") = "F");
            auto buf = nb::cast<nb::ndarray<T, nb::shape<nb_any, nb_any>>>(arr);
            std::memcpy(buf.data(), m.data(), m.rows() * m.cols() * sizeof(T));
            return arr;
        })
        ;
}

template<class T>
void bind_gsSparseMatrix(nb::module_& m, const char* name)
{
    using SpMat = gismo::gsSparseMatrix<T>;
    nb::class_<SpMat>(m, name)
        .def(nb::init<>())
        .def(nb::init<index_t, index_t>())
        .def("rows", [](const SpMat& m) { return m.rows(); })
        .def("cols", [](const SpMat& m) { return m.cols(); })
        .def("nonZeros", [](const SpMat& m) { return m.nonZeros(); })
        .def("coeff", [](const SpMat& m, index_t i, index_t j) -> T { return m.coeff(i, j); })
        .def("__repr__", [](const SpMat& m) {
            std::ostringstream s;
            s << "gsSparseMatrix(" << m.rows() << " x " << m.cols()
              << ", nnz=" << m.nonZeros() << ")";
            return s.str();
        })
        ;
}

NB_MODULE(matrix, m) {
    m.doc() = "G+Smo Matrix module";

    bind_gsVector<real_t>(m, "gsVector");
    bind_gsVector<index_t>(m, "gsVectorInt");
    bind_gsMatrix<real_t>(m, "gsMatrix");
    bind_gsMatrix<index_t>(m, "gsMatrixInt");
    bind_gsSparseMatrix<real_t>(m, "gsSparseMatrix");
}

#endif // GISMO_WITH_NANOBIND
