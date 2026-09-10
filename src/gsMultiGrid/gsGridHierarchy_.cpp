#include <gsMultiGrid/gsGridHierarchy.h>
#include <gsMultiGrid/gsGridHierarchy.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsGridHierarchy<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsGridHierarchy(py::module &m)
{
    using Class = gsGridHierarchy<real_t>;
    using MB = gsMultiBasis<real_t>;
    using BC = gsBoundaryConditions<real_t>;
    using Opts = gsOptionList;
    using Mat = gsSparseMatrix<real_t, RowMajor>;

    py::class_<Class>(m, "gsGridHierarchy")
        // Static factory methods for refinement
        .def_static("buildByRefinement",
                    static_cast<Class (*)(MB, const BC&, const Opts&, index_t, index_t, index_t, index_t)>
                    (&Class::buildByRefinement),
                    py::arg("multiBasis"), py::arg("boundaryConditions"), py::arg("options"),
                    py::arg("levels"), py::arg("numberOfKnots")=1, py::arg("multiplicityOfKnots")=1,
                    py::arg("unk")=0,
                    "Build hierarchy by uniform refinement")

        // Static factory methods for coarsening
        .def_static("buildByCoarsening",
                    static_cast<Class (*)(MB, index_t, const BC&, const Opts&, index_t, index_t, index_t)>
                    (&Class::buildByCoarsening),
                    py::arg("multiBasis"), py::arg("numComponents"), py::arg("boundaryConditions"),
                    py::arg("options"), py::arg("levels"), py::arg("degreesOfFreedom")=1, py::arg("unk")=0,
                    "Build hierarchy by uniform coarsening")

        // Static factory method for hierarchical coarsening (HSplines)
        .def_static("buildByHierarchicalCoarsening",
                    &Class::buildByHierarchicalCoarsening,
                    py::arg("multiBasis"), py::arg("numComponents"), py::arg("boundaryConditions"),
                    py::arg("options"), py::arg("unk")=0,
                    "Build hierarchy by hierarchical coarsening (for HSplines)")

        // Instance methods: accessors
        .def("nLevels",
             [](Class &self) -> index_t {
                 return static_cast<index_t>(self.getMultiBases().size());
             },
             "Get the number of levels in the hierarchy")

        .def("getMBasis",
             [](Class &self, index_t level) -> const MB& {
                 return self.getMultiBases()[level];
             },
             py::arg("level"), py::return_value_policy::reference_internal,
             "Get the multibasis at a given level")

        .def("getTransferMatrix",
             [](Class &self, index_t level) -> const Mat& {
                 return self.getTransferMatrices()[level];
             },
             py::arg("level"), py::return_value_policy::reference_internal,
             "Get the transfer matrix from level to level+1")
        ;
}

#endif

} // namespace gismo
