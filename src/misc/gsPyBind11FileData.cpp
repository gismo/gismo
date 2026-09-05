/** @file gsPyBind11FileData.cpp

    @brief Python bindings of gsFileData (moved out of gsIO/gsFileData_.cpp:
    the typed getId/add bindings reference gsNurbs/gsPde types, which the
    base IO module must not depend on - modularization stream S3, step A7).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsFileData.h>

#ifdef GISMO_WITH_PYBIND11
#include <gsContainers/gsMultiPatch.h>
#include <gsContainers/gsMultiBasis.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsPde/gsBoundaryConditions.h>
#include <gsFunctionExpr/gsFunctionExpr.h>
#include <gsIO/gsOptionList.h>
#include <gsIO/gsOptionListXml.h>
#include <gsMatrix/gsSparseMatrix.h>
#endif




namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

gsMatrix<real_t> getMatrix(std::string & filename)
{
    gsFileData<> file;
    file.read(filename);
    gsMatrix<real_t> result;
    file.getAnyFirst(result);
    return result;
}

namespace py = pybind11;
  void pybind11_init_gsFileData(py::module &m) {

    using Class = gsFileData<real_t>;
    py::class_<Class> fd(m, "gsFileData");

    fd.def(py::init<>())
      .def(py::init<const std::string&>())
      
      .def("read", &Class::read, py::arg("filename"), py::arg("recursive")=false)
      .def("clear", &Class::clear)
      .def("numData", &Class::numData)
      .def("save",           &Class::save,           py::arg("fname")="dump", py::arg("compress")=false)
      .def("saveCompressed", &Class::saveCompressed, py::arg("fname")="dump")
      .def("dump",           &Class::dump,           py::arg("fname")="dump")

      .def("addComment", &Class::addComment)
      .def("lastPath", &Class::lastPath)
      .def("setFloatPrecision", &Class::setFloatPrecision)
      .def("getFloatPrecision", &Class::getFloatPrecision)

      // .def("getId", static_cast<const gsBasis<real_t> & (Class::*)(const size_t) const > (&Class::getId))
      .def("getId", static_cast<void (Class::*)(const int &, gsMultiPatch<real_t> &           ) const > (&Class::getId<gsMultiPatch<real_t>>), "Gets a gsMultiPatch by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsGeometry<real_t> &             ) const > (&Class::getId<gsGeometry<real_t>>  ), "Gets a gsGeometry by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsBSpline<real_t> &              ) const > (&Class::getId<gsBSpline<real_t>>   ), "Gets a gsBSpline by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsTensorBSpline<2,real_t> &      ) const > (&Class::getId<gsTensorBSpline<2,real_t>>), "Gets a gsTensorBSpline by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsTensorBSpline<3,real_t> &      ) const > (&Class::getId<gsTensorBSpline<3,real_t>>), "Gets a gsTensorBSpline by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsTensorBSpline<4,real_t> &      ) const > (&Class::getId<gsTensorBSpline<4,real_t>>), "Gets a gsTensorBSpline by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsBSplineBasis<real_t>  &        ) const > (&Class::getId<gsBSplineBasis<real_t>>), "Gets a gsBSplineBasis by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsTensorBSplineBasis<2,real_t>  &) const > (&Class::getId<gsTensorBSplineBasis<2,real_t>>), "Gets a gsTensorBSplineBasis by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsTensorBSplineBasis<3,real_t>  &) const > (&Class::getId<gsTensorBSplineBasis<3,real_t>>), "Gets a gsTensorBSplineBasis by id")
      .def("getId", static_cast<void (Class::*)(const int &, gsTensorBSplineBasis<4,real_t>  &) const > (&Class::getId<gsTensorBSplineBasis<4,real_t>>), "Gets a gsTensorBSplineBasis by id")

      .def("getId", static_cast<memory::unique_ptr<gsMultiPatch<real_t>>            (Class::*)(const int &) const > (&Class::getId<gsMultiPatch<real_t>>), "Gets a gsMultiPatch by id")
      .def("getId", static_cast<memory::unique_ptr<gsGeometry<real_t>>              (Class::*)(const int &) const > (&Class::getId<gsGeometry<real_t>>  ), "Gets a gsGeometry by id")
      .def("getId", static_cast<memory::unique_ptr<gsBSpline<real_t>>               (Class::*)(const int &) const > (&Class::getId<gsBSpline<real_t>>   ), "Gets a gsBSpline by id")
      .def("getId", static_cast<memory::unique_ptr<gsTensorBSpline<2,real_t>>       (Class::*)(const int &) const > (&Class::getId<gsTensorBSpline<2,real_t>>), "Gets a gsTensorBSpline by id")
      .def("getId", static_cast<memory::unique_ptr<gsTensorBSpline<3,real_t>>       (Class::*)(const int &) const > (&Class::getId<gsTensorBSpline<3,real_t>>), "Gets a gsTensorBSpline by id")
      .def("getId", static_cast<memory::unique_ptr<gsTensorBSpline<4,real_t>>       (Class::*)(const int &) const > (&Class::getId<gsTensorBSpline<4,real_t>>), "Gets a gsTensorBSpline by id")
      .def("getId", static_cast<memory::unique_ptr<gsBSplineBasis<real_t>>          (Class::*)(const int &) const > (&Class::getId<gsBSplineBasis<real_t>>), "Gets a gsBSplineBasis by id")
      .def("getId", static_cast<memory::unique_ptr<gsTensorBSplineBasis<2,real_t>>  (Class::*)(const int &) const > (&Class::getId<gsTensorBSplineBasis<2,real_t>>), "Gets a gsTensorBSplineBasis by id")
      .def("getId", static_cast<memory::unique_ptr<gsTensorBSplineBasis<3,real_t>>  (Class::*)(const int &) const > (&Class::getId<gsTensorBSplineBasis<3,real_t>>), "Gets a gsTensorBSplineBasis by id")
      .def("getId", static_cast<memory::unique_ptr<gsTensorBSplineBasis<4,real_t>>  (Class::*)(const int &) const > (&Class::getId<gsTensorBSplineBasis<4,real_t>>), "Gets a gsTensorBSplineBasis by id")

      .def("add", static_cast<void (Class::*)(const gsMultiPatch<real_t> &, int) > (&Class::add<gsMultiPatch<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsMultiPatch to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsMultiBasis<real_t> &, int) > (&Class::add<gsMultiBasis<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsMultiBasis to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsBSpline<real_t> &, int) > (&Class::add<gsBSpline<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsBSpline to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsTensorBSpline<2, real_t> &, int) > (&Class::add<gsTensorBSpline<2, real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsTensorBSpline to the filedata.")

      .def("add", static_cast<void (Class::*)(const gsBoundaryConditions<real_t> &, int) > (&Class::add<gsBoundaryConditions<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsBoundaryConditions to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsFunctionExpr<real_t> &, int) > (&Class::add<gsFunctionExpr<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsFunctionExpr to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsOptionList &, int) > (&Class::add<gsOptionList>), py::arg("object"), py::arg("id")=-1, "Add gsOptionList to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsMatrix<real_t> &, int) > (&Class::add<gsMatrix<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsMatrix to the filedata.")
      .def("add", static_cast<void (Class::*)(const gsSparseMatrix<real_t> &, int) > (&Class::add<gsSparseMatrix<real_t>>), py::arg("object"), py::arg("id")=-1, "Add gsSparseMatrix to the filedata.")

      //.def("addSparse", static_cast<void (Class::*)(const gsSparseMatrix<real_t> &) > (&Class::add<gsSparseMatrix<real_t>>), "Add gsSparseMatrix to the filedata.")
      
       .def("getAnyFirst", static_cast<bool (Class::*)(gsMultiPatch<real_t> &) const > (&Class::getAnyFirst<gsMultiPatch<real_t>>), "Get gsMultiPatch to the filedata.")
       .def("getAll", static_cast<std::vector< memory::unique_ptr<gsGeometry<real_t>> > (Class::*)() const > (&Class::getAll<gsGeometry<real_t>>), "Get gsGeometry to the filedata.")

        // getLabel overloads (Phase 0a: read labeled sections from XML)
        // Single-arg form: auto-dispatches based on XML tag
        .def("getLabel",
            [](const Class& self, const std::string& label) -> py::object {
                std::string tag = self.getLabelTag(label);
                if (tag.empty())
                    throw py::value_error("No object with label '" + label + "' found.");
                if (tag == internal::gsXml<gsOptionList>::tag()) {
                    auto obj = self.getLabel<gsOptionList>(label);
                    return py::cast(obj.release(), py::return_value_policy::take_ownership);
                }
                if (tag == internal::gsXml<gsMultiPatch<real_t>>::tag()) {
                    auto obj = self.getLabel<gsMultiPatch<real_t>>(label);
                    return py::cast(obj.release(), py::return_value_policy::take_ownership);
                }
                if (tag == internal::gsXml<gsBoundaryConditions<real_t>>::tag()) {
                    auto obj = self.getLabel<gsBoundaryConditions<real_t>>(label);
                    return py::cast(obj.release(), py::return_value_policy::take_ownership);
                }
                if (tag == internal::gsXml<gsMatrix<real_t>>::tag()) {
                    auto obj = self.getLabel<gsMatrix<real_t>>(label);
                    return py::cast(obj.release(), py::return_value_policy::take_ownership);
                }
                if (tag == internal::gsXml<gsSparseMatrix<real_t>>::tag()) {
                    auto obj = self.getLabel<gsSparseMatrix<real_t>>(label);
                    return py::cast(obj.release(), py::return_value_policy::take_ownership);
                }
                throw py::value_error(
                    "getLabel: unrecognized tag '" + tag + "' for label '" + label + "'. "
                    "Use the two-argument form fd.getLabel(label, result) to load this type explicitly.");
            },
            py::arg("label"),
            "Get object by label, auto-detecting type "
            "(OptionList, MultiPatch, boundaryConditions, Matrix, SparseMatrix). "
            "Raises ValueError if the label is not found or the type is not supported by the single-argument form.")

        .def("getLabel", static_cast<void (Class::*)(const std::string &, gsMultiPatch<real_t> &) const > (&Class::getLabel<gsMultiPatch<real_t>>), py::arg("label"), py::arg("result"), "Get gsMultiPatch by label (output param)")
       .def("getLabel", static_cast<void (Class::*)(const std::string &, gsOptionList &) const > (&Class::getLabel<gsOptionList>), py::arg("label"), py::arg("result"), "Get gsOptionList by label (output param)")
       .def("getLabel", static_cast<void (Class::*)(const std::string &, gsBoundaryConditions<real_t> &) const > (&Class::getLabel<gsBoundaryConditions<real_t>>), py::arg("label"), py::arg("result"), "Get gsBoundaryConditions by label (output param)")
       .def("getLabel", static_cast<void (Class::*)(const std::string &, gsMatrix<real_t> &) const > (&Class::getLabel<gsMatrix<real_t>>), py::arg("label"), py::arg("result"), "Get gsMatrix by label (output param)")
       .def("getLabel", static_cast<void (Class::*)(const std::string &, gsSparseMatrix<real_t> &) const > (&Class::getLabel<gsSparseMatrix<real_t>>), py::arg("label"), py::arg("result"), "Get gsSparseMatrix by label (output param)")
       
       // hasLabel (Phase 0a: check if a label exists)
       .def("hasLabel", &Class::hasLabel, py::arg("label"), "Check if label exists")
       
       // addWithLabel overloads (Phase 0a: write labeled sections to XML)
       .def("addWithLabel", static_cast<void (Class::*)(const gsMultiPatch<real_t> &, std::string) > (&Class::addWithLabel<gsMultiPatch<real_t>>), py::arg("object"), py::arg("label"), "Add gsMultiPatch with label")
       .def("addWithLabel", static_cast<void (Class::*)(const gsOptionList &, std::string) > (&Class::addWithLabel<gsOptionList>), py::arg("object"), py::arg("label"), "Add gsOptionList with label")
       .def("addWithLabel", static_cast<void (Class::*)(const gsMatrix<real_t> &, std::string) > (&Class::addWithLabel<gsMatrix<real_t>>), py::arg("object"), py::arg("label"), "Add gsMatrix with label")
       .def("addWithLabel", static_cast<void (Class::*)(const gsSparseMatrix<real_t> &, std::string) > (&Class::addWithLabel<gsSparseMatrix<real_t>>), py::arg("object"), py::arg("label"), "Add gsSparseMatrix with label")
       .def("addWithLabel", static_cast<void (Class::*)(const gsBoundaryConditions<real_t> &, std::string) > (&Class::addWithLabel<gsBoundaryConditions<real_t>>), py::arg("object"), py::arg("label"), "Add gsBoundaryConditions with label")
       
       // getId for additional types (Phase 0a)
       .def("getId", static_cast<void (Class::*)(const int &, gsBoundaryConditions<real_t> &) const > (&Class::getId<gsBoundaryConditions<real_t>>), py::arg("id"), py::arg("result"), "Gets gsBoundaryConditions by id")
       .def("getId", static_cast<void (Class::*)(const int &, gsOptionList &) const > (&Class::getId<gsOptionList>), py::arg("id"), py::arg("result"), "Gets gsOptionList by id")
       
       // getFirst for gsMultiPatch (Phase 0a: alias for getAnyFirst)
       .def("getFirst", static_cast<bool (Class::*)(gsMultiPatch<real_t> &) const > (&Class::getFirst<gsMultiPatch<real_t>>), py::arg("result"), "Get first gsMultiPatch")

       // Work around to obtain the matrix from Filedata. Standard way is not working!
       .def("getMatrix", &getMatrix, "Get any first gsMatrix.")

      .def("bufferSize", &Class::bufferSize)
      .def("print", &Class::print)
      .def("contents", &Class::contents)
      .def("numTags", &Class::numTags)
      

      // .def("getId", (void (gsFileData<real_t>::*)(const int&, Object&)) &Class::getId<Object>)

      .def("__repr__",
           [](const gsFileData<> &obj) {
             std::stringstream os;
             os << obj;
             return os.str();
           })
      ;    
  }
#endif

} // namespace gismo
