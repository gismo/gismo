#include <gsCore/gsTemplateTools.h>

//Prerequisits
#include<gsAssembler/gsQuadrature.h>
#include<gsHSplines/gsHTensorBasis.h>
#include<gsCore/gsLinearAlgebra.h>
#include<gsCore/gsFunction.h>
#include<gsNurbs/gsBSpline.h>
#include<gsUtils/gsCombinatorics.h>

#include <gsHSplines/gsQuasiInterpolate.h>
#include <gsHSplines/gsQuasiInterpolate.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsQuasiInterpolate<real_t>;


#ifdef GISMO_WITH_PYBIND11

	namespace py = pybind11;	

    void pybind11_init_gsQuasiInterpolate(py::module &m)
    {
        using Class = gsQuasiInterpolate<real_t>;
        py::class_<Class>(m, "gsQuasiInterpolate")
        .def_static("localIntpl", static_cast<gsMatrix<real_t> (*)(const gsBasis<real_t>&, const gsFunction<real_t>&, index_t, const gsMatrix<real_t>&)>(&Class::localIntpl), "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"), py::arg("index"), py::arg("ab"))
        .def_static("localIntpl", static_cast<gsMatrix<real_t> (*)(const gsBasis<real_t>&, const gsFunction<real_t>&, index_t)>(&Class::localIntpl), "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"), py::arg("index"))
        .def_static("localIntpl", static_cast<gsMatrix<real_t> (*)(const gsHTensorBasis<1,real_t>&, const gsFunction<real_t>&, index_t)>(&Class::localIntpl), "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"), py::arg("index"))
        .def_static("localIntpl", static_cast<gsMatrix<real_t> (*)(const gsHTensorBasis<2,real_t>&, const gsFunction<real_t>&, index_t)>(&Class::localIntpl), "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"), py::arg("index"))
        .def_static("localIntpl", static_cast<gsMatrix<real_t> (*)(const gsHTensorBasis<3,real_t>&, const gsFunction<real_t>&, index_t)>(&Class::localIntpl), "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"), py::arg("index"))
        .def_static("localIntpl", static_cast<gsMatrix<real_t> (*)(const gsHTensorBasis<4,real_t>&, const gsFunction<real_t>&, index_t)>(&Class::localIntpl), "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"), py::arg("index"))
        .def_static("localIntpl", [](const gsBasis<real_t>& bb, const gsFunction<real_t>& fun) 
                                    { 
                                        gsMatrix<real_t> result;
                                        Class::localIntpl(bb, fun, result);
                                        return result; 
                                    }, "Compute the local quasi-interpolation coefficients for a given basis function and a given function", py::arg("basis"), py::arg("function"))
        ;
    }

#endif


}
