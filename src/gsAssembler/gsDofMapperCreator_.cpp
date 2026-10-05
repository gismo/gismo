#include <gsCore/gsTemplateTools.h>

#include <gsAssembler/gsDofMapperCreator.h>
#include <gsAssembler/gsDofMapperCreator.hpp>

namespace gismo
{

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoxTopology & topology,
                                       const gsBoundaryConditions<real_t> & bc,
                                       index_t nComp, index_t unk,
                                       bool conforming, bool finalize,
                                       gsDofMapper::storage st);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       index_t nComp, bool conforming,
                                       bool finalize, gsDofMapper::storage st);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoxTopology & topology,
                                       index_t nComp, bool conforming,
                                       bool finalize, gsDofMapper::storage st);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoundaryConditions<real_t> & bc,
                                       index_t nComp, index_t unk, bool conforming,
                                       bool finalize, gsDofMapper::storage st);

TEMPLATE_INST gsDofMapper createMapper(const gsFunctionSet<real_t> & bases,
                                       const gsBoundaryConditions<real_t> & bc,
                                       dirichlet::strategy ds, iFace::strategy is,
                                       index_t nComp, index_t unk, bool finalize,
                                       gsDofMapper::storage st);

TEMPLATE_INST gsDofMapper createMapper(const std::vector<const gsFunctionSet<real_t>*> & basesPerComp,
                                       const gsBoxTopology & topology,
                                       const gsBoundaryConditions<real_t> & bc,
                                       index_t unk, bool conforming, bool finalize,
                                       gsDofMapper::storage st);

TEMPLATE_INST gsDofMapper createMapper(const std::vector<gsMultiBasis<real_t> > & basesPerComp,
                                       const gsBoundaryConditions<real_t> & bc,
                                       index_t unk, bool conforming, bool finalize,
                                       gsDofMapper::storage st);

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
          [](const gsFunctionSet<real_t> & bases, const gsBoxTopology & topology,
             const gsBoundaryConditions<real_t> & bc, index_t nComp, index_t unk,
             bool conforming, bool finalize)
          { return createMapper<real_t>(bases, topology, bc, nComp, unk, conforming, finalize); },
          "Create a gsDofMapper (full options)",
          py::arg("bases"),
          py::arg("topology"),
          py::arg("bc"),
          py::arg("nComp")=1,
          py::arg("unk")=-1,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          [](const gsFunctionSet<real_t> & bases, index_t nComp, bool conforming, bool finalize)
          { return createMapper<real_t>(bases, nComp, conforming, finalize); },
          "Create a gsDofMapper (bases only)",
          py::arg("bases"),
          py::arg("nComp")=1,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          [](const gsFunctionSet<real_t> & bases, const gsBoxTopology & topology,
             index_t nComp, bool conforming, bool finalize)
          { return createMapper<real_t>(bases, topology, nComp, conforming, finalize); },
          "Create a gsDofMapper (bases + topology)",
          py::arg("bases"),
          py::arg("topology"),
          py::arg("nComp")=1,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          [](const gsFunctionSet<real_t> & bases, const gsBoundaryConditions<real_t> & bc,
             index_t nComp, index_t unk, bool conforming, bool finalize)
          { return createMapper<real_t>(bases, bc, nComp, unk, conforming, finalize); },
          "Create a gsDofMapper (bases + boundary conditions)",
          py::arg("bases"),
          py::arg("bc"),
          py::arg("nComp")=1,
          py::arg("unk")=0,
          py::arg("conforming")=true,
          py::arg("finalize")=false);

    m.def("createMapper",
          [](const gsFunctionSet<real_t> & bases, const gsBoundaryConditions<real_t> & bc,
             dirichlet::strategy ds, iFace::strategy is,
             index_t nComp, index_t unk, bool finalize)
          { return createMapper<real_t>(bases, bc, ds, is, nComp, unk, finalize); },
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
