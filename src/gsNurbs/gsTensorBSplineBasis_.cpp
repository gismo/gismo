/* Symbol export for G+Smo shared object */

//#define gsTensorBSplineBasis_EXPORT

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>
#include <gsCore/gsBasisFun.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.hpp>

namespace gismo
{

#define INST(D) CLASS_TEMPLATE_INST gsTensorBSplineBasis<D,real_t>;
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

#define INST(D) CLASS_TEMPLATE_INST internal::gsXml< gsTensorBSplineBasis<D,real_t> >;
GISMO_DIM_FOREACH(INST)
#undef INST

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

std::vector<gsKnotVector<real_t>> knot_vector_list_from_args(const py::args & knotVectors)
{
    std::vector<gsKnotVector<real_t>> knotVectorList;

    if (knotVectors.size() == 1 && py::isinstance<py::sequence>(knotVectors[0]))
    {
        const py::sequence sequence = py::reinterpret_borrow<py::sequence>(knotVectors[0]);
        knotVectorList.reserve(py::len(sequence));
        for (auto it : sequence)
            knotVectorList.push_back(py::cast<gsKnotVector<real_t>>(it));
    }
    else
    {
        knotVectorList.reserve(knotVectors.size());
        for (auto it : knotVectors)
            knotVectorList.push_back(py::cast<gsKnotVector<real_t>>(it));
    }

    return knotVectorList;
}

template <short_t D>
gsTensorBSplineBasis<D,real_t> tensor_basis_from_args(const py::args & knotVectors)
{
    std::vector<gsKnotVector<real_t>> knotVectorList = knot_vector_list_from_args(knotVectors);

    if (knotVectorList.size() != D)
        throw py::value_error("Expected either " + std::to_string(D) + " knot vectors as positional arguments or a single sequence containing exactly " + std::to_string(D) + " knot vectors.");

    return gsTensorBSplineBasis<D,real_t>(give(knotVectorList));
}

// C++11-compatible TMP dispatch: try dimension D, recurse down to 2
template<short_t D, typename std::enable_if<(D > 2), int>::type = 0>
py::object make_gsTensorBSplineBasis(size_t n, std::vector<gsKnotVector<real_t>> kvs)
{
    if ((short_t)n == D)
        return py::cast(new gsTensorBSplineBasis<D,real_t>(give(kvs)));
    return make_gsTensorBSplineBasis<D-1>(n, std::move(kvs));
}

template<short_t D, typename std::enable_if<(D == 2), int>::type = 0>
py::object make_gsTensorBSplineBasis(size_t n, std::vector<gsKnotVector<real_t>> kvs)
{
    if ((short_t)n == D)
        return py::cast(new gsTensorBSplineBasis<D,real_t>(give(kvs)));
    throw py::value_error("Expected 2 to " + std::to_string(GISMO_MAX_DIMENSION)
                          + " knot vectors as positional arguments or in a single sequence.");
}

void pybind11_init_gsTensorBSplineBasis_factory(py::module &m)
{
    m.def("gsTensorBSplineBasis", [](py::args knotVectors) -> py::object
    {
        std::vector<gsKnotVector<real_t>> knotVectorList = knot_vector_list_from_args(knotVectors);
        return make_gsTensorBSplineBasis<GISMO_MAX_DIMENSION>(knotVectorList.size(), std::move(knotVectorList));
    },
    "Factory constructor that dispatches to gsTensorBSplineBasis2..N based on the number of knot vectors");
}

template <short_t d>
void pybind11_init_gsTensorBSplineBasis(py::module &m)
{
        using Base = gsBasis<real_t>;
        using Class = gsTensorBSplineBasis<d,real_t>;

        py::class_<Class,Base> cls(m, ("gsTensorBSplineBasis" + std::to_string(d)).c_str());

        cls.def(py::init<>());

        cls.def(py::init([](py::args knotVectors)
            {
                return new Class(tensor_basis_from_args<d>(knotVectors));
            }),
            "Constructor from knot vectors passed either as varargs or as a single sequence");

        cls
        .def("knots", static_cast<      gsKnotVector<real_t>& (Class::*)(int)      > (&Class::knots), "Get the knot vector as a reference")
        .def("knots", static_cast<const gsKnotVector<real_t>& (Class::*)(int) const> (&Class::knots), "Get the knot vector as a const reference")
        .def("component", static_cast<      gsBSplineBasis<real_t>& (Class::*)(int)      > (&Class::component), "Returns the basis component as a reference")
        .def("component", static_cast<const gsBSplineBasis<real_t>& (Class::*)(int) const> (&Class::component), "Returns the basis component as a const reference")
        .def("size", static_cast<index_t (Class::*)() const> (&Class::size), "Returns the size")
        .def("dim", &Class::dim, "Returns the dimension")
        .def("eval", &Class::eval, "Evaluates points into a matrix")
        .def("function", &Class::function, "Returns the basis function i")
        .def("evalSingle", static_cast<gsMatrix<real_t> (Class::*)(index_t, const gsMatrix<real_t> &) const> (&Class::evalSingle), "Evaluates the basis function i")
        .def("degree", &Class::degree, "Returns the degree")
        ;
}

#define INST(D) template void pybind11_init_gsTensorBSplineBasis<D>(py::module &m);
GISMO_DIM_FOREACH_FROM2(INST)
#undef INST

#endif

}
