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
    /// This constructor simply creates a default gsSubdivisionScheme.
    /// No options are set.
    /// The member pointers `m_mesh` is not initialized, so some methods may
    /// fail if used without first initializing it.
    gsSubdivisionScheme() : m_options(gsOptionList()), m_mesh(nullptr) {}

    /// \brief Constructor with a set of options for child classes to use.
    /// The member pointers `m_mesh` is not initialized, so some methods may
    /// fail if used without first initializing it.
    explicit gsSubdivisionScheme(gsOptionList& options)
        : m_options(options), m_mesh(nullptr)
    {
    }

    /// \brief Constructor with a mesh to operate on.
    explicit gsSubdivisionScheme(gsSurfMesh* mesh)
        : m_options(gsOptionList()), m_mesh(mesh)
    {
    }

    /// \brief Constructor with a mesh to operate on and a set of options for
    /// child classes to use.
    gsSubdivisionScheme(gsSurfMesh* mesh, gsOptionList& options)
        : m_options(options), m_mesh(mesh)
    {
    }

public: // Destructor
    /// \brief Normal destructor.
    virtual ~gsSubdivisionScheme() {}

protected: // Members
    /// GsOptions used to customize the subdivision scheme.
    gsOptionList m_options;

    /// The mesh this subdivision scheme operates on.
    /// Can be set with the constructor, or with certain methods.
    /// This mesh will be referred to as the 'targeted mesh' in further
    /// documentation.
    gsSurfMesh* m_mesh;

public: // Options
    /// \brief Getter function for option manipulation.
    ///
    /// Returns possible options for the chosen subdivision scheme.
    gsOptionList& options() { return m_options; }

public: // Subdivision method
    /// \brief Subdivides the targeted mesh based on the chosen algorithm.
    ///
    /// The main subdivision method. Mutates the mesh `m_mesh` managed by this
    /// by sudividing it. Has to be overriden by implementors.
    virtual void subdivide() = 0;

    /// \brief Subdivides a mesh based on the chosen algorithm.
    ///
    /// Takes in a pointer to a mesh and mutates that mesh by sudividing it.
    /// Sets the member mesh to the given mesh.
    ///
    /// \param mesh The mesh to be subdivided. Will be changed in the operation.
    virtual void subdivide(gsSurfMesh& mesh) { subdivide(mesh, 1); };

    /// \brief Repeated subdivision.
    ///
    /// Runs subdivision algorithm defined in `subdivide` multiple times.
    ///
    /// \param mesh The mesh to be subdivided. Will be changed in the operation.
    /// \param repetitions How often to repeath the subdivision algorithm.
    void subdivide(gsSurfMesh& mesh, size_t repetitions)
    {
        m_mesh = &mesh;
        for (size_t i = 0; i < repetitions; ++i)
            this->subdivide();
    }

    /// \brief Repeated subdivision on the targeted mesh.
    ///
    /// Runs subdivision algorithm defined in `subdivide` multiple times.
    /// Operates on the targeted mesh and mutates it.
    ///
    /// \param repetitions How often to repeath the subdivision algorithm.
    void subdivide(size_t repetitions)
    {
        for (size_t i = 0; i < repetitions; ++i)
            this->subdivide();
    }

    /// \brief Subdivides a mesh out-of-place.
    ///
    /// Returns a subdivided copy of the given mesh without mutating the
    /// original, using the algorithm defined in `subdivide`. Does not mutate
    /// the given mesh.
    ///
    /// \param mesh The mesh to be subdivided. Will not be changed in the
    /// operation.
    gsSurfMesh subdivided(gsSurfMesh& mesh)
    {
        return this->subdivided(mesh, 1);
    }

    /// \brief Subdivides the targeted mesh out-of-place.
    ///
    /// Returns a subdivided copy of the targeted mesh without mutating the
    /// original, using the algorithm defined in `subdivide`. Does not mutate
    /// the given mesh.
    gsSurfMesh subdivided() { return this->subdivided(1); }

    /// \brief Subdivides a mesh out-of-place, multiple times.
    ///
    /// Returns a copy of the given mesh that has been subdivided multiple times
    /// without mutating the original, using the algorithm defined in
    /// `subdivide`. Does not mutate the given mesh.
    /// Sets the targeted mesh to the given mesh.
    ///
    /// \param mesh The mesh to be subdivided. Will not be changed in the
    /// operation.
    /// \param repetitions How often to repeath the subdivision algorithm.
    gsSurfMesh subdivided(gsSurfMesh& mesh, size_t repetitions)
    {
        m_mesh = &mesh;
        auto mesh_copy = gsSurfMesh(mesh);
        this->subdivide(mesh_copy, repetitions);
        return mesh_copy;
    }

    /// \brief Subdivides the targeted mesh out-of-place, multiple times.
    ///
    /// Returns a copy of the targeted mesh that has been subdivided multiple
    /// times without mutating the original, using the algorithm defined in
    /// `subdivide`. Does not mutate the targeted mesh.
    ///
    /// \param repetitions How often to repeath the subdivision algorithm.
    gsSurfMesh subdivided(size_t repetitions)
    {
        auto mesh_copy = gsSurfMesh(*m_mesh);
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
        /// Applying the algorithm to the mesh will yield meaningless results or
        /// cause an error.
        INVALID = 0,
        /// The scheme can be applied to the mesh without problems.
        VALID = 1,
        /// No statement about the validity of applying the algorithm to the
        /// mesh is made.
        UNDETERMINED = 2,
    };

    /// \brief Checks if the given mesh is valid for this subdivision scheme.
    ///
    /// May return an `UNDETERMINED` to signify that this method was not
    /// overwritten by the implementor or the scheme cannot decide (yet) if the
    /// given mesh is valid. If one of the other options is returned, the
    /// contained answer must be definitive.
    ///
    /// This method also changes the targeted mesh to the given mesh.
    ///
    /// \return A \c gsSubdivisionMeshValidity value indicating the result.
    ///
    /// \param mesh The mesh we want to apply this algorithm to.
    gsSubdivisionMeshValidity check_mesh(gsSurfMesh& mesh)
    {
        m_mesh = &mesh;
        return check_mesh();
    }

    /// \brief Checks if the targeted mesh is valid for this subdivision scheme.
    ///
    /// Must not necessarily be overwritten by children classes and will return
    /// an an indeterminate answer in that case. May also return an
    /// `UNDETERMINED` to signify that the scheme cannot decide (yet) if the
    /// given mesh is valid. If one of the other options is returned, the
    /// contained answer must be definitive.
    ///
    /// \return A \c gsSubdivisionMeshValidity value indicating the result.
    ///
    /// \param mesh The mesh we want to apply this algorithm to.
    virtual gsSubdivisionMeshValidity check_mesh()
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
    void subdivide() {}
};

} // namespace gismo
