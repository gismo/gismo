#pragma once

//== INCLUDES =================================================================


#include <gsMesh2/Vector.h>
#include <gsCore/gsConfig.h>


//=============================================================================


namespace gismo {


//=============================================================================


/// Scalar type
typedef real_t Scalar;

/// Point type
typedef Vector<Scalar,3> Point;
//typedef Eigen::Vector<Scalar,3> Point;

/// 3D vector type
typedef Vector<Scalar,3> Vec3;

/// Normal type
typedef Vector<Scalar,3> Normal;

/// Color type
typedef Vector<Scalar,3> Color;

/// Texture coordinate type
typedef Vector<Scalar,3> Texture_coordinate;


//=============================================================================
} // namespace gismo
//=============================================================================

