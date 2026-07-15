#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsFunction.h>
#include <gsCore/gsFunction.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsFunction<real_t> ;

#ifdef GISMO_WITH_PYBIND11  

namespace py = pybind11;

/* Define a trampoline class for gsFunction */
template <typename T>
class gsFunctionTrampoline : public gsFunction<T>
{
public:
    using Base = gsFunction<T>;
    using Base::Base; // Inherit constructors

    gsFunctionTrampoline() = default;

    // We must override the core virtual methods G+Smo uses.
    // For _into methods, pybind11 copies Eigen matrices by value when passing to Python,
    // so we cannot use PYBIND11_OVERRIDE directly for output-parameter methods.
    // Instead we pre-allocate result, create a non-owning numpy view into its memory,
    // and pass that view to Python so modifications are written back to the C++ matrix.

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        py::gil_scoped_acquire acquire;
        py::function overload = py::get_override(static_cast<const Base*>(this), "eval_into");
        if (overload) {
            result.resize(this->targetDim(), u.cols());
            // Use a capsule as base to prevent numpy from copying the data.
            // The capsule has a no-op destructor since result owns the memory.
            py::capsule base(result.data(), [](void*) {});
            py::array_t<T> result_view(
                {(py::ssize_t)result.rows(), (py::ssize_t)result.cols()},
                {(py::ssize_t)sizeof(T), (py::ssize_t)(result.rows() * sizeof(T))},
                result.data(), base);
            overload(u, result_view);
        } else {
            pybind11::pybind11_fail("Tried to call pure virtual function \"gsFunction::eval_into\"");
        }
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        py::gil_scoped_acquire acquire;
        py::function overload = py::get_override(static_cast<const Base*>(this), "deriv_into");
        if (overload) {
            result.resize(this->targetDim() * this->domainDim(), u.cols());
            py::capsule base(result.data(), [](void*) {});
            py::array_t<T> result_view(
                {(py::ssize_t)result.rows(), (py::ssize_t)result.cols()},
                {(py::ssize_t)sizeof(T), (py::ssize_t)(result.rows() * sizeof(T))},
                result.data(), base);
            overload(u, result_view);
        } else {
            Base::deriv_into(u, result);
        }
    }

    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        py::gil_scoped_acquire acquire;
        py::function overload = py::get_override(static_cast<const Base*>(this), "deriv2_into");
        if (overload) {
            const short_t ddim = this->domainDim();
            const index_t nders = ddim * (ddim + 1) / 2;
            result.resize(this->targetDim() * nders, u.cols());
            py::capsule base(result.data(), [](void*) {});
            py::array_t<T> result_view(
                {(py::ssize_t)result.rows(), (py::ssize_t)result.cols()},
                {(py::ssize_t)sizeof(T), (py::ssize_t)(result.rows() * sizeof(T))},
                result.data(), base);
            overload(u, result_view);
        } else {
            Base::deriv2_into(u, result);
        }
    }

    short_t domainDim() const override
    {
        PYBIND11_OVERRIDE_PURE(short_t, Base, domainDim, );
    }

    short_t targetDim() const override
    {
        PYBIND11_OVERRIDE_PURE(short_t, Base, targetDim, );
    }

    std::ostream &print(std::ostream &os) const override
    {
        py::gil_scoped_acquire acquire;
        py::function overload = py::get_overload(static_cast<const Base *>(this), "print");
        if (overload) {
            os << overload().cast<std::string>();
        } else {
            // Only error out if Python hasn't provided a print method
            GISMO_ERROR("print() is not implemented in the Python-derived class.");
        }
        return os;
    }

private:

    // The trampoline clones the base class, not itself!
    Base * clone_impl() const override
    {
        py::gil_scoped_acquire acquire;
        
        // 1. Check if the Python class implemented a 'clone' method
        py::function overload = py::get_overload(static_cast<const Base *>(this), "clone");
        
        if (overload) 
        {
            // 2. Call Python clone. It should return a NEW Python object.
            // pybind11 will wrap that new Python object in a NEW C++ Trampoline.
            auto new_py_obj = overload();
            
            // 3. Release ownership of the C++ pointer from Python to G+Smo
            return new_py_obj.release().template cast<Base*>();
        }

        // Fallback: If no Python clone exists, this is an error for a trampoline
        GISMO_ERROR("A Python-derived gsFunction must implement a clone() method "
                    "using 'return copy.copy(self)' to support parallel assembly.");
    }

};

void pybind11_init_gsFunction(py::module &m)
{
    using Base = gsFunctionSet<real_t>;
    using Class = gsFunction<real_t>;
    using Trampoline = gsFunctionTrampoline<real_t>;
    py::class_<Class, Trampoline, Base>(m, "gsFunction")
        .def(py::init<>())
        .def("jacobian",  &Class::jacobian, "Returns the Jacobian")
        .def("hessian",   &Class::hessian, "Returns the Hessian")
        .def("laplacian", &Class::laplacian, "Returns the Laplacian")
        .def("argMin", &Class::argMin, "Returns the position of the minimum",
             py::arg("accuracy") = 1e-6, py::arg("max_loop") = 100,
             py::arg("init") = gsMatrix<real_t>(),
             py::arg("damping_factor") = 1 )
  ;
}
#endif

}
