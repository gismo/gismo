#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsDofMapperCreator.h>
#include <gsAssembler/gsDofMapperCreator.hpp>

namespace gismo
{

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoxTopology & topology,
                                       const gsBoundaryConditions<real_t> & bc,
                                       index_t nComp, index_t unk,
                                       bool conforming, bool finalize);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       index_t nComp, bool conforming,
                                       bool finalize);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoxTopology & topology,
                                       index_t nComp, bool conforming,
                                       bool finalize);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoundaryConditions<real_t> & bc,
                                       index_t nComp, index_t unk, bool conforming,
                                       bool finalize);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoundaryConditions<real_t> & bc,
                                       dirichlet::strategy ds, iFace::strategy is,
                                       index_t nComp, index_t unk, bool finalize);

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsDofMapperCreator(py::module &m)
{
    py::enum_<dirichlet::strategy>(m, "DirichletStrategy")
        .value("elimination", dirichlet::strategy::elimination)
        .value("penalize",     dirichlet::strategy::penalize)
        .value("nitsche",      dirichlet::strategy::nitsche)
        .value("eliminatNormal", dirichlet::strategy::eliminatNormal)
        .value("none",         dirichlet::strategy::none)
        .export_values();

    py::enum_<iFace::strategy>(m, "InterfaceStrategy")
        .value("conforming", iFace::strategy::conforming)
        .value("glue",       iFace::strategy::glue)
        .value("dg",        iFace::strategy::dg)
        .value("smooth",    iFace::strategy::smooth)
        .value("none",      iFace::strategy::none)
        .export_values();

    m.def("createMapper",
          static_cast<gsDofMapper (*)(const gsFunctionSet<real_t>&,
                                      const gsBoxTopology&,
                                      const gsBoundaryConditions<real_t>&,
                                      index_t, index_t, bool, bool)>(&createMapper),
          "Create a gsDofMapper (full options)",
          py::arg("bases"),
          py::arg("topology"),
          py::arg("bc"),
          py::arg("nComp")=1,
          py::arg("unk")=-1,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          static_cast<gsDofMapper (*)(const gsFunctionSet<real_t>&,
                                      index_t, bool, bool)>(&createMapper),
          "Create a gsDofMapper (bases only)",
          py::arg("bases"),
          py::arg("nComp")=1,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          static_cast<gsDofMapper (*)(const gsFunctionSet<real_t>&,
                                      const gsBoxTopology&,
                                      index_t, bool, bool)>(&createMapper),
          "Create a gsDofMapper (bases + topology)",
          py::arg("bases"),
          py::arg("topology"),
          py::arg("nComp")=1,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          static_cast<gsDofMapper (*)(const gsFunctionSet<real_t>&,
                                      const gsBoundaryConditions<real_t>&,
                                      index_t, index_t, bool, bool)>(&createMapper),
          "Create a gsDofMapper (bases + boundary conditions)",
          py::arg("bases"),
          py::arg("bc"),
          py::arg("nComp")=1,
          py::arg("unk")=0,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          static_cast<gsDofMapper (*)(const gsFunctionSet<real_t>&,
                                      const gsBoundaryConditions<real_t>&,
                                      dirichlet::strategy, iFace::strategy,
                                      index_t, index_t, bool)>(&createMapper),
          "Create a gsDofMapper (bases + boundary conditions + strategies)",
          py::arg("bases"),
          py::arg("bc"),
          py::arg("ds"),
          py::arg("is"),
          py::arg("nComp")=1,
          py::arg("unk")=0,
          py::arg("finalize")=false);
}

#endif

} // namespace gismo
