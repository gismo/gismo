/** @file gsKdNode.h

    @brief Provides declaration of the tree node.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mokris
*/

# pragma once

namespace gismo {

/**
    @brief Struct of for an Axis-aligned bounding box

    Template parameters
    \param d is the dimension
    \param Z is the box-coordinate index type
    
    \ingroup HSplines
*/

template<short_t d, class Z = index_t>
struct gsAABB
{
public:
    typedef gsVector<Z, d> point;

    gsAABB(const point & low, const point & upp, index_t lvl)
    :
    first(low), second(upp), level(lvl)
    { }

    gsAABB(const point & low, const point & upp)
    :
    gsAABB(low, upp, -1)
    { }

    gsAABB(const point & upp)
    :
    second(upp), level(-1)
    {
        first.setZero();
    }

    gsAABB()
    :
    level(-1)
    {
        first.setZero();
        second.setZero();
    }

    /// Copy constructor (makes deep copy)
    gsAABB( const gsAABB<d,Z>& other )
    {
        operator=(other);
    }

    /// Move constructor
    gsAABB( gsAABB<d,Z>&& other )
    {
        operator=(give(other));
    }

    /// Assignment operator
    gsAABB<d,Z>& operator= ( const gsAABB<d,Z>& other )
    {
        if (this!=&other)
        {
            first  = other.first;
            second = other.second;
            level  = other.level;
        }
        return *this;
    }

    /// Move assignment operator
    gsAABB<d,Z>& operator= ( gsAABB<d,Z>&& other )
    {
        first  = give(other.first);
        second = give(other.second);
        level  = give(other.level);
        return *this;
    }

public:

    point first;
    point second;

    /// Level in which the box lives
    index_t level;

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
};

} // namespace gismo
