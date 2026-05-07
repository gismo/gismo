#include <gsCore/gsTemplateTools.h>

//Prerequisits
#include <gsUtils/gsProjection.h>
#include <gsUtils/gsProjection.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsProjection<ProjectionNorm::L2,real_t>;
STRUCT_TEMPLATE_INST gsProjection<ProjectionNorm::H1,real_t>;
STRUCT_TEMPLATE_INST gsProjection<ProjectionNorm::H2,real_t>;

// #ifdef GISMO_WITH_PYBIND11

// 	namespace py = pybind11;	

//     void pybind11_init_gsL2Projection(py::module &m)
//     {
//         using Class = gsL2Projection<real_t>;
//         py::class_<Class>(m, "gsL2Projection")
//         .def_static("project", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; // Release the GIL while performing the projection, to enable parallelism in Python if the projection is called from Python. The GIL will be re-acquired when the lambda returns.
//                 error = Class::project(projectionBasis, geometryMap, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a geometry onto a basis (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsMultiBasis<real_t> & integrationBasis, const gsFunctionSet<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; // Release the GIL while performing the projection, to enable parallelism in Python if the projection is called from Python. The GIL will be re-acquired when the lambda returns.
//                 error = Class::project(integrationBasis, projectionBasis, geometryMap, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a geometry onto a basis (multi-patch) with custom integration basis", py::arg("integrationBasis"), py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; 
//                 error = Class::project(projectionBasis, geometryMap, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a geometry onto a basis (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; 
//                 error = Class::project(projectionBasis, integrationBasis, geometryMap, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a geometry onto a basis (single patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; 
//                 error = Class::project(projectionBasis, geometryMap, sourceFunction, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a geometry onto a basis (multi-patch) with custom source function", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; 
//                 error = Class::project(projectionBasis, integrationBasis, geometryMap, sourceFunction, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a geometry onto a basis (multi-patch) with custom integration basis and source function", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; 
//                 error = Class::project(projectionBasis, geometryMap, sourceFunction, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a function onto a basis (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsMatrix<real_t> C;
//             real_t error;
//             { 
//                 py::gil_scoped_release release; 
//                 error = Class::project(projectionBasis, integrationBasis, geometryMap, sourceFunction, C, options); 
//             }
//             return std::make_pair(error, C);
//         }, "Project a function onto a basis (single patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, geometryMap, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for geometry projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsMultiBasis<real_t> & integrationBasis, const gsFunctionSet<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(integrationBasis, projectionBasis, geometryMap, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for geometry projection with custom integration basis (multi-patch)", py::arg("integrationBasis"), py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, geometryMap, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for geometry projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, integrationBasis, geometryMap, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for geometry projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, geometryMap, sourceFunction, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for function projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsFunctionSet<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, integrationBasis, geometryMap, sourceFunction, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for function projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, geometryMap, sourceFunction, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for function projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("system", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::system(projectionBasis, integrationBasis, geometryMap, sourceFunction, systemMatrix, rhs, options);
//             }
//             return std::make_pair(systemMatrix, rhs);
//         }, "Obtain system matrix and right-hand side for function projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         // Matrix bindings
//         .def_static("matrix", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, short_t targetDim, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             {
//                 py::gil_scoped_release release;
//                 Class::matrix(projectionBasis, geometryMap, systemMatrix, targetDim, options);
//             }
//             return systemMatrix;
//         }, "Obtain the mass matrix for geometry projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("targetDim") = 1, py::arg("options") = gsOptionList())
//         .def_static("matrix", [](const gsMultiBasis<real_t> & integrationBasis, const gsFunctionSet<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, short_t targetDim, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             {
//                 py::gil_scoped_release release;
//                 Class::matrix(integrationBasis, projectionBasis, geometryMap, systemMatrix, targetDim, options);
//             }
//             return systemMatrix;
//         }, "Obtain the mass matrix for geometry projection with custom integration basis (multi-patch)", py::arg("integrationBasis"), py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("targetDim") = 1, py::arg("options") = gsOptionList())
//         .def_static("matrix", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, short_t targetDim, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             {
//                 py::gil_scoped_release release;
//                 Class::matrix(projectionBasis, geometryMap, systemMatrix, targetDim, options);
//             }
//             return systemMatrix;
//         }, "Obtain the mass matrix for geometry projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("targetDim") = 1, py::arg("options") = gsOptionList())
//         .def_static("matrix", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, short_t targetDim, const gsOptionList & options)
//         {
//             gsSparseMatrix<real_t> systemMatrix;
//             {
//                 py::gil_scoped_release release;
//                 Class::matrix(projectionBasis, integrationBasis, geometryMap, systemMatrix, targetDim, options);
//             }
//             return systemMatrix;
//         }, "Obtain the mass matrix for geometry projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("targetDim") = 1, py::arg("options") = gsOptionList())
//         // Right-hand side bindings
//         .def_static("rhs", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, geometryMap, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for geometry projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsMultiBasis<real_t> & integrationBasis, const gsFunctionSet<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
//         {
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(integrationBasis, projectionBasis, geometryMap, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for geometry projection with custom integration basis (multi-patch)", py::arg("integrationBasis"), py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
//         {            
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, geometryMap, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for geometry projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
//         {            
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, integrationBasis, geometryMap, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for geometry projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, geometryMap, sourceFunction, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for function projection (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsFunctionSet<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
//         {
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, integrationBasis, geometryMap, sourceFunction, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for function projection with custom integration basis (multi-patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
//         {            
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, geometryMap, sourceFunction, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for function projection (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         .def_static("rhs", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
//         {            
//             gsMatrix<real_t> rhs;
//             {
//                 py::gil_scoped_release release;
//                 Class::rhs(projectionBasis, integrationBasis, geometryMap, sourceFunction, rhs, options);
//             }
//             return rhs;
//         }, "Obtain the right-hand side for function projection with custom integration basis (single patch)", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
//         // 
//         ;
//     }

// #endif

}
