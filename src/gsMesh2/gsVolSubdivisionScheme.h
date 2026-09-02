/** @file gsVolSubdivisionScheme.h

    @brief Parent class for subdivision operations on a volumetric mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsIO/gsOptionList.h>
#include <gsMesh2/gsVolMesh.h>

namespace gismo
{

/** \brief The abstract base class for volumetric subdivision schemes.

    The volumetric counterpart of gsSubdivisionScheme, with the same shape: a
    scheme is assigned a target mesh, carries a gsOptionList, and mutates the
    target when subdivide() is called.

    The one behavioural difference is worth knowing about. A half-edge mesh can
    be refined by local surgery, so the surface schemes edit their mesh in
    place. A 3-map cannot: refining a cell replaces it by one cell per corner
    and rewires every dart, so an implementation of subdivide_impl() is expected
    to build the refined mesh and move-assign it over the target.

    Should not be instantiated directly.

    \tparam Scalar coefficient type of the vertex positions

    \sa gsSubdivisionScheme, gsVolCatmullClark, gsVolMesh
*/
template <class Scalar>
class gsVolSubdivisionScheme
{
protected: // Type definitions for mesh components.
    typedef typename gsVolMesh<Scalar>::Point    Point;
    typedef typename gsVolMesh<Scalar>::Vertex   Vertex;
    typedef typename gsVolMesh<Scalar>::Edge     Edge;
    typedef typename gsVolMesh<Scalar>::Face     Face;
    typedef typename gsVolMesh<Scalar>::Cell     Cell;
    typedef typename gsVolMesh<Scalar>::Corner   Corner;
    typedef typename gsVolMesh<Scalar>::Halfedge Halfedge;
    typedef typename gsVolMesh<Scalar>::Halfface Halfface;

public:
    /// Assign the mesh to be subdivided
    void assign(gsVolMesh<Scalar>* mesh) { m_mesh = mesh; }

protected: // Constructor
    /// \brief Default constructor without options.
    ///
    /// The member pointer `m_mesh` is not initialized, so some methods may fail
    /// if used without first assigning a mesh.
    gsVolSubdivisionScheme() : m_options(gsOptionList()), m_mesh(nullptr)
    { }

public: // Destructor
    virtual ~gsVolSubdivisionScheme() {}

protected: // Members
    /// Options used to customize the subdivision scheme.
    gsOptionList m_options;

    /// The mesh this subdivision scheme operates on.
    gsVolMesh<Scalar>* m_mesh;

public: // Options
    /// \brief Getter function for option manipulation.
    gsOptionList& options() { return m_options; }

    /// Returns the mesh
    gsVolMesh<Scalar> & mesh()
    { GISMO_ASSERT(nullptr!=m_mesh,"Invalid mesh"); return *m_mesh; }

    /// Returns the mesh in const context
    const gsVolMesh<Scalar> & mesh() const
    { GISMO_ASSERT(nullptr!=m_mesh,"Invalid mesh"); return *m_mesh; }

public: // Subdivision method

    /// \brief Repeated subdivision of the targeted mesh.
    ///
    /// \param repetitions how often to repeat the subdivision algorithm
    void subdivide(size_t repetitions = 1)
    {
        for (size_t i = 0; i < repetitions; ++i)
            this->subdivide_impl();
    }

protected:
    /// \brief Subdivides the targeted mesh based on the chosen algorithm.
    ///
    /// Has to be overridden by implementors.
    virtual void subdivide_impl() = 0;

public: // Validity
    /// \brief Validity of a subdivision algorithm for a given mesh.
    enum gsSubdivisionMeshValidity
    {
        /// Applying the algorithm will yield meaningless results or an error.
        INVALID = 0,
        /// The scheme can be applied to the mesh without problems.
        VALID = 1,
        /// No statement about the validity is made.
        UNDETERMINED = 2
    };

    /// \brief Checks if the given mesh is valid for this scheme, and targets it.
    gsSubdivisionMeshValidity check_mesh(gsVolMesh<Scalar>& mesh)
    {
        m_mesh = &mesh;
        return check_mesh();
    }

    /// \brief Checks if the targeted mesh is valid for this scheme.
    virtual gsSubdivisionMeshValidity check_mesh()
    {
        return gsSubdivisionMeshValidity::UNDETERMINED;
    }
};

} // namespace gismo
