#include <gsCore/gsTemplateTools.h>

//Prerequisits
#include <gsUtils/gsL2Projection.h>
#include <gsUtils/gsL2Projection.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsL2Projection<real_t>;

#ifdef GISMO_WITH_PYBIND11

	namespace py = pybind11;	

    void pybind11_init_gsL2Projection(py::module &m)
    {
        using Class = gsL2Projection<real_t>;
        py::class_<Class>(m, "gsL2Projection")
        .def_static("project", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; // Release the GIL while performing the projection, to enable parallelism in Python if the projection is called from Python. The GIL will be re-acquired when the lambda returns.
                error = Class::project(projectionBasis, geometryMap, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a geometry onto a basis (multi-patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsMultiBasis<real_t> & integrationBasis, const gsFunctionSet<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; // Release the GIL while performing the projection, to enable parallelism in Python if the projection is called from Python. The GIL will be re-acquired when the lambda returns.
                error = Class::project(integrationBasis, projectionBasis, geometryMap, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a geometry onto a basis (multi-patch) with custom integration basis", py::arg("integrationBasis"), py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; 
                error = Class::project(projectionBasis, geometryMap, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a geometry onto a basis (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; 
                error = Class::project(projectionBasis, integrationBasis, geometryMap, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a geometry onto a basis (single patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsMultiBasis<real_t> & projectionBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; 
                error = Class::project(projectionBasis, geometryMap, sourceFunction, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a geometry onto a basis (multi-patch) with custom source function", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsFunctionSet<real_t> & projectionBasis, const gsMultiBasis<real_t> & integrationBasis, const gsMultiPatch<real_t> & geometryMap, const gsFunctionSet<real_t> & sourceFunction, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; 
                error = Class::project(projectionBasis, integrationBasis, geometryMap, sourceFunction, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a geometry onto a basis (multi-patch) with custom integration basis and source function", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; 
                error = Class::project(projectionBasis, geometryMap, sourceFunction, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a function onto a basis (single patch)", py::arg("projectionBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
        .def_static("project", [](const gsBasis<real_t> & projectionBasis, const gsBasis<real_t> & integrationBasis, const gsGeometry<real_t> & geometryMap, const gsFunction<real_t> & sourceFunction, const gsOptionList & options)
        {
            gsMatrix<real_t> C;
            real_t error;
            { 
                py::gil_scoped_release release; 
                error = Class::project(projectionBasis, integrationBasis, geometryMap, sourceFunction, C, options); 
            }
            return std::make_pair(error, C);
        }, "Project a function onto a basis (single patch) with custom integration basis", py::arg("projectionBasis"), py::arg("integrationBasis"), py::arg("geometryMap"), py::arg("sourceFunction"), py::arg("options") = gsOptionList())
        ;
    }

#endif

}
