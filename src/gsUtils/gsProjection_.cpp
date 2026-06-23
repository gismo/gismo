#include <gsCore/gsTemplateTools.h>

//Prerequisits
#include <gsUtils/gsProjection.h>
#include <gsUtils/gsProjection.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsProjection<ProjectionNorm::L2,real_t>;
STRUCT_TEMPLATE_INST gsProjection<ProjectionNorm::H1,real_t>;
STRUCT_TEMPLATE_INST gsProjection<ProjectionNorm::H2,real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_enum_gsProjectionNorm(py::module &m)
{
    py::enum_<ProjectionNorm>(m, "ProjectionNorm")
        .value("L2",    ProjectionNorm::L2)
        .value("H1",    ProjectionNorm::H1)
        .value("H2",    ProjectionNorm::H2)
        .export_values();
}

// Helper struct to map ProjectionNorm to Python class names
template<ProjectionNorm Norm> struct projectionNormName;
template<> struct projectionNormName<L2>    { static constexpr const char* value = "gsL2Projection"; };
template<> struct projectionNormName<H1>    { static constexpr const char* value = "gsH1Projection"; };
template<> struct projectionNormName<H2>    { static constexpr const char* value = "gsH2Projection"; };

template<ProjectionNorm Norm> void pybind11_init_gsProjection(py::module &m)
{
    using Class = gsProjection<Norm, real_t>;
    py::class_<Class>(m, projectionNormName<Norm>::value)
    .def_static("project", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, geometryMap, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a geometry onto a basis (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, integrationBasis, geometryMap, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a geometry onto a basis (multi-patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, geometryMap, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a geometry onto a basis (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, integrationBasis, geometryMap, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a geometry onto a basis (single patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, geometryMap, sourceFunction, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a function onto a basis (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, integrationBasis, geometryMap, sourceFunction, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a function onto a basis (multi-patch) with custom integration basis and source function", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, geometryMap, sourceFunction, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a function onto a basis (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> C;
        real_t error;
        {
            py::gil_scoped_release release;
            error = Class::project(projectionBasis, integrationBasis, geometryMap, sourceFunction, C, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(error, C);
    }, "Project a function onto a basis (single patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    // system bindings
    .def_static("system", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, geometryMap, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for geometry projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, integrationBasis, geometryMap, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for geometry projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, geometryMap, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for geometry projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, integrationBasis, geometryMap, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for geometry projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, geometryMap, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for function projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, integrationBasis, geometryMap, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for function projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, geometryMap, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for function projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::system(projectionBasis, integrationBasis, geometryMap, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
        }
        return std::make_pair(systemMatrix, rhs);
    }, "Obtain system matrix and right-hand side for function projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    // matrix bindings
    .def_static("matrix", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, short_t targetDim, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        {
            py::gil_scoped_release release;
            Class::matrix(projectionBasis, geometryMap, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
        }
        return systemMatrix;
    }, "Obtain the mass matrix for geometry projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("targetDim") = -1, py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("matrix", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, short_t targetDim, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        {
            py::gil_scoped_release release;
            Class::matrix(projectionBasis, integrationBasis, geometryMap, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
        }
        return systemMatrix;
    }, "Obtain the mass matrix for geometry projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("targetDim") = -1, py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("matrix", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, short_t targetDim, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        {
            py::gil_scoped_release release;
            Class::matrix(projectionBasis, geometryMap, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
        }
        return systemMatrix;
    }, "Obtain the mass matrix for geometry projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("targetDim") = -1, py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("matrix", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, short_t targetDim, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsSparseMatrix<real_t> systemMatrix;
        {
            py::gil_scoped_release release;
            Class::matrix(projectionBasis, integrationBasis, geometryMap, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
        }
        return systemMatrix;
    }, "Obtain the mass matrix for geometry projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("targetDim") = -1, py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    // rhs bindings
    .def_static("rhs", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, geometryMap, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for geometry projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, integrationBasis, geometryMap, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for geometry projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, geometryMap, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for geometry projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, integrationBasis, geometryMap, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for geometry projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, geometryMap, sourceFunction, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for function projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, integrationBasis, geometryMap, sourceFunction, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for function projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, geometryMap, sourceFunction, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for function projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsBoundaryConditions<real_t> & bc, const gsOptionList & options, real_t alpha, real_t beta, real_t gamma)
    {
        gsMatrix<real_t> rhs;
        {
            py::gil_scoped_release release;
            Class::rhs(projectionBasis, integrationBasis, geometryMap, sourceFunction, rhs, bc, options, alpha, beta, gamma);
        }
        return rhs;
    }, "Obtain the right-hand side for function projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("bc") = gsBoundaryConditions<real_t>(), py::arg("options") = gsOptionList(), py::arg("alpha") = (real_t)1.0, py::arg("beta") = (real_t)1.0, py::arg("gamma") = (real_t)1.0)
    ;
}

// Explicit instantiations
template void pybind11_init_gsProjection<L2>(py::module &m);
template void pybind11_init_gsProjection<H1>(py::module &m);
template void pybind11_init_gsProjection<H2>(py::module &m);

#endif
}
