/** @file gsSubdivisionScheme.h

    @brief Parent class for subdivision operations on mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher, A. Mantzaflaris, D.Tolis
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsOptionList.h>
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

public:
    /// Assign the mesh to be subdivided
    void assign(gsSurfMesh* mesh) { m_mesh = mesh; }

protected: // Constructor
    /// \brief Default constructor without options.
    /// This constructor simply creates a default gsSubdivisionScheme.
    /// No options are set.
    /// The member pointers `m_mesh` is not initialized, so some methods may
    /// fail if used without first initializing it.
    gsSubdivisionScheme() : m_options(gsOptionList()), m_mesh(nullptr)
    {
        m_options.addSwitch("normalize", "Normalize limit normals and tangents", true);
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

    /// Returns the mesh
    gsSurfMesh & mesh() { GISMO_ASSERT(nullptr!=m_mesh,"Invalid mesh"); return *m_mesh; }

    /// Returns the mesh in const context
    const gsSurfMesh & mesh() const { GISMO_ASSERT(nullptr!=m_mesh,"Invalid mesh"); return *m_mesh; }

public: // Subdivision method

/// \brief Repeated subdivision on the targeted mesh.
    ///
    /// Runs subdivision algorithm defined in `subdivide` multiple times.
    /// Operates on the targeted mesh and mutates it.
    ///
    /// \param repetitions How often to repeath the subdivision algorithm.
    void subdivide(size_t repetitions)
    {
        for (size_t i = 0; i < repetitions; ++i)
            this->subdivide_impl();
    }

protected:
    /// \brief Subdivides the targeted mesh based on the chosen algorithm.
    ///
    /// The main subdivision method. Mutates the mesh `m_mesh` managed by this
    /// by sudividing it. Has to be overriden by implementors.
    virtual void subdivide_impl() = 0;

public:

    /// Compute vertex limit positions
    virtual gsSurfMesh::Vertex_property<Point> vertex_limits(std::string label = "v:limit")
    {GISMO_NO_IMPLEMENTATION}

    /// Compute vertex limit normals
    virtual gsSurfMesh::Vertex_property<Point> vertex_normal_limits(std::string label = "v:normal")
    {GISMO_NO_IMPLEMENTATION}

    /// Compute vertex limit tangent
    virtual gsSurfMesh::Vertex_property<Point> vertex_tangent_limits(std::string label = "v:tanvec")
    {GISMO_NO_IMPLEMENTATION}

    /// Compute face limit positions
    virtual gsSurfMesh::Face_property<Point> face_limits(std::string label = "f:limit")
    {GISMO_NO_IMPLEMENTATION}

    /// Compute face limit normals
    virtual gsSurfMesh::Face_property<Point> face_normal_limits(std::string label = "f:normal")
    {GISMO_NO_IMPLEMENTATION}

    /// Compute face limit tangent
    virtual gsSurfMesh::Face_property<Point> face_tangent_limits(std::string label = "f:tanvec")
    {GISMO_NO_IMPLEMENTATION}

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
    void subdivide_impl() GISMO_OVERRIDE {}
};

} // namespace gismo
