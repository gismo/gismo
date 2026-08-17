
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiBasis.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsMultiBasis<real_t>;

TEMPLATE_INST bool gsMultiBasis<real_t>::
repairInterfaceFindElements<2>(const boundaryInterface & bi,
                               std::vector<index_t> & refEltsFirst,
                               std::vector<index_t> & refEltsSecond );

TEMPLATE_INST bool gsMultiBasis<real_t>::
repairInterfaceFindElements<3>(const boundaryInterface & bi,
                               std::vector<index_t> & refEltsFirst,
                               std::vector<index_t> & refEltsSecond );

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsMultiBasis(py::module &m)
  {
    using Base = gsFunctionSet<real_t>;
    using Class = gsMultiBasis<real_t>;
    py::class_<Class,Base>(m, "gsMultiBasis")

      // Constructors
      .def(py::init<>())
      .def(py::init<const gsMultiPatch<real_t> &>(),
            py::arg("numeratorOnly") = false) //default arguments

      // Member functions
      .def("domainDim", &Class::domainDim, "Returns the domain dimension of the multipatch")
      .def("targetDim", &Class::targetDim, "Returns the target dimension of the multipatch")
      .def("nBases", &Class::nBases, "Returns the number of patches stored in the multipatch")
      .def("basis", static_cast<const gsBasis<real_t> & (Class::*)(const size_t) const > (&Class::basis), "Gets a const reference to basis with index i",py::return_value_policy::reference)
      .def("basis", static_cast<      gsBasis<real_t> & (Class::*)(const size_t)       > (&Class::basis) , "Gets a const reference to basis with index i",py::return_value_policy::reference)
      // Note: Bindings with unique pointers are not possible https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html
      // .def("addPatch", static_cast<void (Class::*)(typename gsBasis<real_t>::uPtr)> (&Class::addBasis), "Adds a patch")
      .def("addBasis", static_cast<void (Class::*)(         gsBasis<real_t> *    )> (&Class::addBasis), "Adds a patch")
      .def("degreeElevate", &Class::degreeElevate, "Degree elevate all patches by the given amount", py::arg("i") = 1, py::arg("dir") = -1)
      .def("degreeIncrease", &Class::degreeIncrease, "Degree increase all patches by the given amount", py::arg("i") = 1, py::arg("dir") = -1)
      .def("degreeDecrease", &Class::degreeDecrease, "Degree decrease all patches by the given amount", py::arg("i") = 1, py::arg("dir") = -1)
      .def("degreeReduce", &Class::degreeReduce, "Degree reduce all patches by the given amount", py::arg("amount"))
      .def("setDegree", &Class::setDegree, "Set the degree of all patches to the given value", py::arg("degree"))
      .def("reduceContinuity", &Class::reduceContinuity, "Reduce the continuity of all patches by the given amount", py::arg("amount") = 1)
      .def("uniformRefine", &Class::uniformRefine, "Uniformly refine all patches by the given amount", py::arg("numKnots") = 1, py::arg("mul") = 1, py::arg("dir") = -1)
      ;
  }

#endif
}
