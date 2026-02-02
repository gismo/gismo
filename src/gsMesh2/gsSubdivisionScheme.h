/** @file gsSubdivisionScheme.h

    @brief Parent class for subdivision operations on mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#pragma once

#include "gsCore/gsExport.h"
#include <gsIO/gsOptionList.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

/// The abstract base class for subdivision schemes.
/// Should not be instantiated directly.
class GISMO_EXPORT gsSubdivisionScheme
{
protected: // Type definitions for mesh components.
    using Point = gsSurfMesh::Point;
    using Vertex = gsSurfMesh::Vertex;
    using Face = gsSurfMesh::Face;
    using Halfedge = gsSurfMesh::Halfedge;
    using Edge = gsSurfMesh::Edge;

protected: // Constructors
    /// \brief Default constructor without options.
    gsSubdivisionScheme() : m_options(gsOptionList()) {}

    /// \brief Constructor with a set of options for child classes to use.
    explicit gsSubdivisionScheme(gsOptionList& options) : m_options(options) {}

public: // Destructor
    /// \brief Normal destructor.
    virtual ~gsSubdivisionScheme() {}

public: // Options
    /// GsOptions used to customize the subdivision scheme.
    gsOptionList m_options;

    /// \brief Getter function for option manipulation.
    ///
    /// Returns possible options for the chosen subdivision scheme.
    gsOptionList& options() { return m_options; }

public: // Subdivision method
    /// \brief Subdivides a mesh based on the chosen algorithm.
    /// The main subdivision method.
    /// Takes in a pointer to a mesh and mutates that mesh by sudividing it.
    /// Has to be overriden by implementors.
    ///
    /// \param mesh The mesh to be subdivided. Will be changed in the operation.
    virtual void subdivide(gsSurfMesh& mesh) = 0;

    /// \brief Repeated subdivision.
    ///
    /// Runs subdivision algorithm defined in `subdivide` multiple times.
    ///
    /// \param mesh The mesh to be subdivided. Will be changed in the operation.
    /// \param repetitions How often to repeath the subdivision algorithm.
    void subdivide(gsSurfMesh& mesh, size_t repetitions)
    {
        for (size_t i = 0; i < repetitions; ++i)
            this->subdivide(mesh);
    }

    /// \brief Subdivides a mesh in place.
    ///
    /// Returns a subdivided copy of the given mesh without mutating the
    /// original, using the algorithm defined in `subdivide`. Does not mutate
    /// the given mesh.
    ///
    /// \param mesh The mesh to be subdivided. Will not be changed in the
    /// operation.
    /// \param repetitions How often to repeath the subdivision algorithm.
    gsSurfMesh subdivided(const gsSurfMesh& mesh)
    {
        return this->subdivided(mesh, 1);
    }

    /// \brief Subdivides a mesh in place, multiple times.
    ///
    /// Returns a copy of the given mesh that has been subdivided multiple times
    /// without mutating the original, using the algorithm defined in
    /// `subdivide`. Does not mutate the given mesh.
    ///
    /// \param mesh The mesh to be subdivided. Will not be changed in the
    /// operation.
    /// \param repetitions How often to repeath the subdivision algorithm.
    gsSurfMesh subdivided(const gsSurfMesh& mesh, size_t repetitions)
    {
        auto mesh_copy = gsSurfMesh(mesh);
        this->subdivide(mesh_copy, repetitions);
        return mesh_copy;
    }

public: // Validity
    /// \brief Validity of a subdivision algorithm.
    ///
    /// An enumeration describing the validity of applying a given subdivision
    /// algorithm to a given mesh.
    enum gsSubdivisionMeshValidity
    {
        /// The scheme can be applied to the mesh without problems.
        VALID,
        /// Applying the algorithm to the mesh will yield meaningless results or
        /// cause an error.
        INVALID,
        /// No statement about the validity of applying the algorithm to the
        /// mesh is made.
        UNDETERMINED
    };

    /// \brief Checks if the given mesh is valid for this subdivision scheme.
    ///
    /// Must not necessarily be overwritten by children classes and will return
    /// an an indeterminate answer in that case. May also return an
    /// `UNDETERMINED` to signify that the scheme cannot decide (yet) if the
    /// given mesh is valid. If one of the other options is returned, the
    /// contained answer must be definitive.
    ///
    /// \param mesh The mesh we want to apply this algorithm to.
    virtual gsSubdivisionMeshValidity valid_mesh(const gsSurfMesh& mesh)
    {
        return gsSubdivisionMeshValidity::UNDETERMINED;
    }

}; // namespace internal

/// \brief Zero-Op Subdivision
///
/// An 'identity subdivision' that leaves the passed mesh untouched.
/// Potentially useful for chaining.
class GISMO_EXPORT gsIdentityScheme : gsSubdivisionScheme
{
public:
    gsIdentityScheme() : gsSubdivisionScheme() {}

public:
    void subdivide(gsSurfMesh& mesh) {}
};

} // namespace gismo
