/** @file gsDofMapper_.cpp.in

    @brief instantiation of gsDofMapper

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan
**/

#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsDofMapper.h>
#include <gsCore/gsDofMapper.hpp>


namespace gismo {

    TEMPLATE_INST void gsDofMapper::init(
         const gsMultiBasis<real_t> & bases, index_t nComp);

    TEMPLATE_INST void gsDofMapper::init(
            std::vector<const gsMultiBasis<real_t> *> const & bases);

      TEMPLATE_INST void gsDofMapper::init(
        const gsMultiBasis<real_t>         &basis,
        const gsBoundaryConditions<real_t> &bc, int unk
        );

    TEMPLATE_INST void gsDofMapper::initSingle(
        const gsBasis<real_t> & bases, index_t nComp);

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsDofMapper(py::module &m)
{
    using Class = gsDofMapper;
    py::class_<Class>(m, "gsDofMapper")
    // Constructors
    .def(py::init<const gsMultiBasis<real_t> &>(),
        py::arg("nComp") = 1) //default arguments
    .def(py::init<const gsBasis<real_t> &>(),
        py::arg("nComp") = 1) //default arguments

    // Member functions
    .def("asVector",                    &Class::asVector,                   "Returns a vector taking flat local indices to global")
    .def("inverseAsVector",             &Class::inverseAsVector,            "Returns a vector taking global indices to flat local")
    // Following requires binding of gsBoxTopology
    // .def("setMatchingInterfaces",       &Class::setMatchingInterfaces,      "Called to initialize the gsDofMapper with matching interfaces after m_bases have already been set")
    .def("colapseDofs",                 &Class::colapseDofs,                "Calls matchDof() for all dofs on the given patch side i ps. Thus, the whole set of dofs collapses to a single global dof")
    .def("matchDof", &Class::matchDof, "Couples dof \a i of patch \a u with dof \a j of patch \a v such that they refer to the same global dof at component \a comp.")
    .def("matchDofs", &Class::matchDofs, "Couples dofs \a b1 of patch \a u with dofs \a b2 of patch \a v one by one such that they refer to the same global dof.")
    .def("markCoupled", &Class::markCoupled, "Mark the local dof \a i of patch \a k as coupled.")
    .def("markTagged", &Class::markTagged, "Mark a local dof \a i of patch \a k as tagged")
    .def("markCoupledAsTagged", &Class::markCoupledAsTagged, "Mark all coupled dofs as tagged")
    .def("markBoundary", &Class::markBoundary, "Mark the local dofs \a boundaryDofs of patch \a k as eliminated. ")
    .def("eliminateDof", &Class::eliminateDof, "Mark the local dof \a i of patch \a k as eliminated.")
    .def("finalize", &Class::finalize, "Must be called after all boundaries and interfaces have been marked to set up the dof numbering.")
    .def("isFinalized", &Class::isFinalized, "Checks whether finalize() has been called.")
    .def("isPermutation", &Class::isPermutation, "Returns true iff the mapper is a permuatation")

    .def("setIdentity", &Class::setIdentity, "Set this mapping to be the identity")
    .def("setShift", &Class::setShift, "Set the shift amount for the global numbering")
    .def("addShift", &Class::addShift, "Add a shift amount to the global numbering")

    .def("index", &Class::index, "Returns the global dof index associated to local dof \a i of patch \a k.")
    .def("bindex", &Class::bindex, "Returns the boundary index of local dof \a i of patch \a k.")
    .def("cindex", &Class::cindex, "Returns the coupled dof index")
    .def("tindex", &Class::tindex, "Returns the tagged dof index")
    .def("global_to_bindex", &Class::global_to_bindex, "Returns the boundary index of global dof \a gl.")
    .def("is_free_index", &Class::is_free_index, "Returns true if global dof \a gl is not eliminated.")
    .def("is_free", &Class::is_free, "Returns true if local dof \a i of patch \a k is not eliminated.")
    .def("is_boundary_index", &Class::is_boundary_index, "Returns true if global dof \a gl is eliminated")
    .def("is_boundary", &Class::is_boundary, "Returns true if local dof \a i of patch \a k is eliminated.")
    .def("is_coupled", &Class::is_coupled, "Returns true if local dof \a i of patch \a k is coupled.")
    .def("is_coupled_index", &Class::is_coupled_index, "Returns true if \a gl is a coupled dof.")
    .def("is_tagged", &Class::is_tagged, "Returns true if local dof \a i of patch \a k is tagged.")
    .def("is_tagged_index", &Class::is_tagged_index, "Returns true if \a gl is a tagged dof.")
    .def("numComponents", &Class::numComponents, "Returns the number of components present in the mapper")
    .def("size", static_cast<index_t (Class::*)()        const > (&Class::size), "Returns the total number of dofs (free and eliminated).")
    .def("size", static_cast<index_t (Class::*)(index_t) const > (&Class::size), "Returns the total number of dofs (free and eliminated).")
    .def("freeSize", static_cast<index_t (Class::*)()        const > (&Class::freeSize), "Returns the number of free (not eliminated) dofs.")
    .def("freeSize", static_cast<index_t (Class::*)(index_t) const > (&Class::freeSize), "Returns the number of free (not eliminated) dofs.")
    .def("coupledSize", &Class::coupledSize, "Returns the number of coupled (not eliminated) dofs.")
    .def("taggedSize", &Class::taggedSize, "Returns the number of tagged dofs.")
    .def("boundarySize", &Class::boundarySize, "Returns the number of eliminated dofs.")

    .def("offset", &Class::offset, "Returns the offset corresponding to patch \a k")
    .def("numPatches", &Class::numPatches, "Returns the number of patches present underneath the mapper")
    .def("mapSize", &Class::mapSize, "Returns the total number of patch-local degrees of freedom that are being mapped")
    .def("componentsSize", &Class::componentsSize, "Returns the components size")
    .def("patchSize", &Class::patchSize, "Returns the total number of patch-local DoFs that live on patch \a k for component \a c")
    .def("totalSize", &Class::totalSize, "Returns the total size of the mapper")
    .def("indexOnPatch", &Class::indexOnPatch, "For \a gl being a global index, this function returns true whenever \a gl corresponds to patch \a k")
    ;
}
#endif
}


