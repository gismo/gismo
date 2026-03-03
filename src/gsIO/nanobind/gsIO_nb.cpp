/** @file gsIO_nb.cpp

    @brief NanoBind bindings for the gsIO module

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_NANOBIND

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <gismo.h>

namespace nb = nanobind;
using namespace gismo;

NB_MODULE(gsIO, m) {
    m.doc() = "G+Smo IO module";

    nb::class_<gsOptionList>(m, "gsOptionList")
        .def(nb::init<>())
        .def("getString", &gsOptionList::getString)
        .def("getInt", &gsOptionList::getInt)
        .def("getReal", &gsOptionList::getReal)
        .def("getSwitch", &gsOptionList::getSwitch)
        .def("setString", &gsOptionList::setString)
        .def("setInt", &gsOptionList::setInt)
        .def("setReal", &gsOptionList::setReal)
        .def("setSwitch", &gsOptionList::setSwitch)
        .def("addString", &gsOptionList::addString)
        .def("addInt", &gsOptionList::addInt)
        .def("addReal", &gsOptionList::addReal)
        .def("addSwitch", &gsOptionList::addSwitch)
        .def("remove", &gsOptionList::remove)
        .def("hasGroup", &gsOptionList::hasGroup)
        .def("getGroup", &gsOptionList::getGroup)
        .def("size", &gsOptionList::size)
        .def("__repr__", [](const gsOptionList& o) {
            std::ostringstream os; os << o; return os.str();
        });

    nb::class_<gsFileData<real_t>>(m, "gsFileData")
        .def(nb::init<>())
        .def(nb::init<const std::string &>())
        .def("read", static_cast<bool (gsFileData<real_t>::*)(const std::string &, bool)>(&gsFileData<real_t>::read),
             nb::arg("fn"), nb::arg("recursive") = false)
        .def("save", static_cast<void (gsFileData<real_t>::*)(const std::string &, bool) const>(&gsFileData<real_t>::save),
             nb::arg("fname"), nb::arg("compress") = false)
        .def("hasAnyMultiPatch", [](const gsFileData<real_t>& fd) {
            return fd.template hasAny< gsMultiPatch<real_t> >();
        })
        .def("hasAnyGeometry", [](const gsFileData<real_t>& fd) {
            return fd.template hasAny< gsGeometry<real_t> >();
        })
        .def("hasAnyBasis", [](const gsFileData<real_t>& fd) {
            return fd.template hasAny< gsBasis<real_t> >();
        })
        .def("numData", &gsFileData<real_t>::numData)
        .def("getFirst", [](gsFileData<real_t>& fd) {
            gsMultiPatch<real_t> mp;
            fd.getFirst(mp);
            return mp;
        })
        .def("clear", &gsFileData<real_t>::clear)
        .def("__repr__", [](const gsFileData<real_t>& fd) {
            std::ostringstream os; os << fd; return os.str();
        });

    m.def("gsReadFile", [](const std::string& fn) {
        gsMultiPatch<real_t> mp;
        gsReadFile<>(fn, mp);
        return mp;
    }, nb::arg("filename"));

    m.def("gsWriteParaview", [](const gsMultiPatch<real_t>& mp,
                                const std::string& fn, index_t npts) {
        gsWriteParaview(mp, fn, npts);
    }, nb::arg("mp"), nb::arg("filename"), nb::arg("npts") = 1000);
}

#endif // GISMO_WITH_NANOBIND
