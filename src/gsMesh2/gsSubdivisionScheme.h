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
#include <gsMesh2/gsSurfMesh.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/// The base class for subdivision schemes.
/// Should not be instantiated directly.
class GISMO_EXPORT gsSubdivisionScheme
{

protected: // Types

    /// Type definitions for mesh components.
    typedef gsSurfMesh::Point Point;
    typedef gsSurfMesh::Vertex Vertex;
    typedef gsSurfMesh::Face Face;


protected: // Constructors

    /// Default constructor without options.
    gsSubdivisionScheme() : m_options(gsOptionList()) {}

    /// Constructor with a set of options for child classes to use.
    explicit gsSubdivisionScheme(gsOptionList& options) : m_options(options) {}

public: // Destructor
    virtual ~gsSubdivisionScheme(){} 

public: // Options

    /// GsOptions used to customize the subdivision scheme.
    gsOptionList m_options;
    /// Returns possible options for the chosen subdivision Scheme.
    gsOptionList& options() { return m_options; }

public: // Subdivision method

    /// The main subdivision method.
    /// Takes in a pointer to a mesh and mutates that mesh by sudividing it.
    /// Has to be overriden by implementors.
    virtual void subdivide(gsSurfMesh* mesh) = 0;

    /// Runs subdivision multiple times.
    void subdivide(gsSurfMesh* mesh, unsigned int repetitions)
    {
      for (unsigned int _i = 0; _i < repetitions; ++_i) 
        this->subdivide(mesh);
    }

    /// Returns a subdivided copy of the given mesh without mutating the original.
    gsSurfMesh subdivided(const gsSurfMesh& mesh)
    {
      return this->subdivided(mesh, 1);
    }

    /// Returns a copy of the original mesh subdivided multiple times, without mutating the original.
    gsSurfMesh subdivided(const gsSurfMesh& mesh, unsigned int repetitions)
    {
      auto mesh_copy = gsSurfMesh(mesh);
      this->subdivide(&mesh_copy, repetitions);
      return mesh_copy;
    }

public: // Validity

    enum gsSubdivisionMeshValidity
    {
      VALID,
      INVALID,
      UNDETERMINED
    };


    /// Checks if the given mesh is valid for this subdivision scheme.
    /// Must not necessarily be overwritten by children classes and will return an empty optional (i.e. no answer) in this case. May also return an empty optional to signify that the scheme cannot decide (yet) if the given mesh is valid.
    /// If a non-empty optional is returned, the contained answer must be definitive.
    virtual gsSubdivisionMeshValidity valid_mesh(const gsSurfMesh& mesh)
    {
      return gsSubdivisionMeshValidity::UNDETERMINED;
    }

};//namespace internal

/// An 'identity subdivision' that leaves the passed mesh untouched.
/// Potentially useful for chaining.
class GISMO_EXPORT gsIdentityScheme : gsSubdivisionScheme {
public:
    gsIdentityScheme() : gsSubdivisionScheme() {}
public:
    void subdivide(gsSurfMesh* mesh) {}
};

} // namespace gismo
