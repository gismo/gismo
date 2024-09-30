
#include <gsCore/gsLinearAlgebra.h>
#include <gsPde/gsPointLoads.h>

namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsPointLoads
   */
  namespace py = pybind11;

  void pybind11_init_gsPointLoads(pybind11::module &m)
  {
    using Class = gsPointLoads<real_t>;
    py::class_<Class>(m, "gsPointLoads")

    // Constructors
    .def(py::init<>())

    // Member functions
    .def("clear", &Class::clear, "Clears the object")

    // .def("addLoad", static_cast<void (Class::*)(const gsVector<real_t> &, const gsVector<real_t> &           ) > (&Class::addLoad), "Adds a point load")
    // .def("addLoad", static_cast<void (Class::*)(const gsVector<real_t> &, const gsVector<real_t> &, int      ) > (&Class::addLoad), "Adds a point load")
    .def("addLoad", static_cast<void (Class::*)(const gsVector<real_t> &, const gsVector<real_t> &, int, bool) > (&Class::addLoad), "Adds a point load")

    // .def("addLoad", static_cast<void (Class::*)(const gsVector<real_t> &, const T                       ) > (&Class::addLoad), "Adds a point load")
    // .def("addLoad", static_cast<void (Class::*)(const gsVector<real_t> &, const T            , int      ) > (&Class::addLoad), "Adds a point load")
    // .def("addLoad", static_cast<void (Class::*)(const gsVector<real_t> &, const T            , int, bool) > (&Class::addLoad), "Adds a point load")

    .def("numLoads", &Class::numLoads, "Returns the number of point loads")
    ;
  }
#endif // GISMO_WITH_PYBIND11

}
