/** @file gsFunctionSet_.cpp

    @brief instantiation of gsFunctionSet

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsFunctionSet.hpp>

namespace gismo {

CLASS_TEMPLATE_INST gsFunctionSet<real_t>;

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

/* Define a trampoline class for gsFunctionSet */
template <typename T>
class gsFunctionSetTrampoline : public gsFunctionSet<T>
{
public:
    using Base = gsFunctionSet<T>;
    using Base::Base; // Inherit constructors

    const gsFunction<T>& piece(index_t p) const override
    {
        PYBIND11_OVERRIDE_PURE(const gsFunction<T>&, Base, piece, p);
    }

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
            py::capsule base(result.data(), [](void*) {});
            py::array_t<T> result_view(
                {(py::ssize_t)result.rows(), (py::ssize_t)result.cols()},
                {(py::ssize_t)sizeof(T), (py::ssize_t)(result.rows() * sizeof(T))},
                result.data(), base);
            overload(u, result_view);
        } else {
            Base::eval_into(u, result);
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
            result.resize(this->targetDim() * (ddim * (ddim + 1) / 2), u.cols());
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
        GISMO_ERROR("A Python-derived gsFunctionSet must implement a clone() method "
                    "using 'return copy.copy(self)' to support parallel assembly.");
    }
};

void pybind11_init_gsFunctionSet(py::module &m)
{
  using Class = gsFunctionSet<real_t>;
  using Trampoline = gsFunctionSetTrampoline<real_t>;
  py::class_<Class, Trampoline>(m, "gsFunctionSet")
    .def(py::init<>())
    .def("piece", &Class::piece, "Returns the piece with index p")
    .def("eval_into", &Class::eval_into, "Evaluates the function set into a matrix")
    .def("deriv_into", &Class::deriv_into, "Evaluates the first derivative into a matrix")
    .def("deriv2_into", &Class::deriv2_into, "Evaluates the second derivative into a matrix")
    .def("evalAllDers_into", &Class::evalAllDers_into, "Evaluates all derivatives upto certien order into a vector of matrices")
    .def("active", &Class::active, "Evaluates the actives and returns a matrix")
    .def("eval", &Class::eval, "Evaluates the function set and returns a matrix")
    .def("deriv", &Class::deriv, "Evaluates the first derivative and returns a matrix")
    .def("deriv2", &Class::deriv2, "Evaluates the second derivative and returns a matrix")
    .def("evalAllDers", &Class::evalAllDers, "Evaluates all derivatives upto certien order into a vector of matrices")
    .def("domainDim", &Class::domainDim, "Returns the domain dimension")
    .def("targetDim", &Class::targetDim, "Returns the target dimension")
    .def("__str__",
        [] (Class & self)
        {
            std::ostringstream os;
            self.print(os);
            return os.str();
        },
        "Returns a string with information about the object.")
  ;
}
#endif


}

namespace std {

TEMPLATE_INST void swap(gismo::gsFuncData<real_t> & f1, gismo::gsFuncData<real_t> & f2);


}
