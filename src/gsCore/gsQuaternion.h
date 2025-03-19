/** @file gsPreCICE.h

    @brief Header file for using gsPreCICE extension

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Jingya Li
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsAffineFunction.h>

#include <gsUtils/gsCombinatorics.h>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#ifndef GISMO_DEG2RAD
#define GISMO_DEG2RAD (M_PI / 180.0)
#endif

#ifndef GISMO_RAD2DEG
#define GISMO_RAD2DEG (180.0 / M_PI)
#endif
// static constexpr double GISMO_DEG2RAD = M_PI / 180.0;
// static constexpr double GISMO_RAD2DEG = 180.0 / M_PI;

namespace gismo 
{
/**
 * @brief Class defining the quaternion objects, which is a four dimensional numbers,
 * unit quaternions can be used to represent rotations in 3d (Euler parameters)
 * 
 * @tparam T 
 */


template<class T>
class gsQuaternion
{
public:

    /**
     * @brief      Default constructor
     */
    gsQuaternion();

    /**
     * @brief      Constructor from four scalars, the first scalar is the real part
     * and the other three are the imaginary parts
     *
     * @param[in]  e0        The first scalar
     * @param[in]  e1        The second scalar
     * @param[in]  e2        The third scalar
     * @param[in]  e3        The fourth scalar
     */

    gsQuaternion(T e0, T e1, T e2, T e3);

    /**
     * @brief      Constructor from a vector and scalar on the real part
     *
     * @param[in]  vec       The vector of four scalars
     */

    gsQuaternion(T e1, const gsVector<T> & vec);

    /**
     * @brief Copy constructor
     * 
     */
    gsQuaternion(const gsQuaternion<T> & other);

    /**
     * @brief      Copy constructor from a quaternion of different type
     */

    template<class U>
    gsQuaternion(const gsQuaternion<U> & other);


    /**
     * @brief Access operator
     * 
     */

    T& e0() {return m_data[0];}
    T& e1() {return m_data[1];}
    T& e2() {return m_data[2];}
    T& e3() {return m_data[3];}

    const T& e0() const {return m_data[0];}
    const T& e1() const {return m_data[1];}
    const T& e2() const {return m_data[2];}
    const T& e3() const {return m_data[3];}

    /**
     * @brief Access operator for the data
     * 
     */
    T* data() {return m_data;}
    const T* data() const {return m_data;}

    /** 
     * @brief Setter and getter
     */
    void Set(T e0, T e1, T e2, T e3);

    /**
     * @brief Set from a scalar
     * 
     */

    void Set(T scalar);

    /**
     * @brief Set from another object
     * 
     */
    void Set(const gsQuaternion<T> & other);


    /**
     * @brief Set the quaternion from a vector 
     * 
     */
    void Set(T e1, const gsVector<T> & vec);

    /**
     * @brief Set to a null quaternion
     * 
     */

    void SetNull();

    /**
     * @brief Set to a unit quaternion
     * 
     */
    void SetUnit();

    /**
     * @brief Set the scalar part only
     * 
     */
    void SetScalar(T scalar);

    /**
     * @brief Set the vector part only
     * 
     */

    void SetVector(const gsVector<T> & vec);

    /**
     * @brief Reture true if the quaternion is null
     * 
     */
    bool IsNull() const;

    /**
     * @brief Return true if the quaternion is identity
     * 
     */
    bool IsIdentity() const;

    /**
     * @brief equals to another quaternion
     * 
     */
    bool Equals(const gsQuaternion<T> & other) const;

    /**
     * @brief Return true if the quaternion equals to another quaternion within a tolerance
     * 
     */
    bool Equals(const gsQuaternion<T> & other, T tol) const;

    /**
     * @brief Return true if the quaternion equals to another quaternion
     * 
     */
    bool operator==(const gsQuaternion<T> & other) const;

    /**
     * @brief Getter for the vector part
     * 
     */
    gsVector<T> GetVector() const;

    /**
     * @brief Get the x-axis  of a coordinate system defined by the quaternion
     * It assumes the quaternion is normalized, so it represents a rotation
     */
    gsVector<T> GetXAxis() const;

    /**
     * @brief Get the y-axis  of a coordinate system defined by the quaternion
     * It assumes the quaternion is normalized, so it represents a rotation
     */
    gsVector<T> GetYAxis() const;

    /**
     * @brief Get the z-axis  of a coordinate system defined by the quaternion
     * It assumes the quaternion is normalized, so it represents a rotation
     */
    gsVector<T> GetZAxis() const;

    //------------------------------norms-----------------------------------
    /**
     * @brief Compute the euclidean norm of the quaternion, that is its length or maginitude
     * 
     * @return T 
     */

    T norm() const;

    /**
     * @brief Compute the L2 euclidean norm of the quaternion
     * 
     * @return T 
     */
    T norm2() const;
    /**
     * @brief The infinity norm, which is the maximum absolute value of one of the components
     * 
     */
    T normInf() const;

    //------------------------------Operator Overloading-----------------------------------
    /**
     * @brief Operator overloading for the addition of two quaternions
     * 
     */
    gsQuaternion<T> operator+() const;
    /**
     * @brief   Operator overloading for the substraction of a quaternion
     * 
     * @return gsQuaternion<T> 
     */
    gsQuaternion<T> operator-() const;

    /**
     * @brief Operator overloading for the addition of two quaternions
     * 
     */
    gsQuaternion<T> operator+(const gsQuaternion<T> & other) const;
    gsQuaternion<T>& operator+=(const gsQuaternion<T> & other);
    /**
     * @brief Operator overloading for the subtraction of two quaternions
     * 
     */
    gsQuaternion<T> operator-(const gsQuaternion<T> & other) const;
    gsQuaternion<T>& operator-=(const gsQuaternion<T> & other);

    /**
     * @brief Quaternion product
     * @details The unit quaternions can represent rotations, the
     * product of two unit quaternions is the composition of the rotations
     * 
     * NOTE!!!: The quaternion product is not commutative
     */
    gsQuaternion<T> operator*(const gsQuaternion<T> & other) const;
    gsQuaternion<T>& operator*=(const gsQuaternion<T> & other);

    /**
     * @brief Subscript operator overloading
     */
    T& operator[](unsigned index);
    const T& operator[](unsigned index) const;

    /**
     * @brief Assignment operator, copy the quaternion
     * 
     */
    gsQuaternion<T>& operator=(const gsQuaternion<T> & other);

    /**
     * @brief Operator for sign change
     * 
     */
    // gsQuaternion<T> operator+() const;
    // gsQuaternion<T> operator-() const;

    /**
     * @brief Operator for making a conjugate quaternion
     * A conjugate quaternion has the vectorial part with changed sign.
     */
    gsQuaternion<T> operator!() const;

    /**
     * @brief Operator for quaternion sum
     * 
     */
    // gsQuaternion<T> operator+(const gsQuaternion<T> & other) const;
    // gsQuaternion<T>& operator+=(const gsQuaternion<T> & other);

    /**
     * @brief Operator for quaternion difference
     * 
     */
    // gsQuaternion<T> operator-(const gsQuaternion<T> & other) const;
    // gsQuaternion<T>& operator-=(const gsQuaternion<T> & other);

    /**
     * @brief Operator for specular quaternion product: A >> B = B * A
     * 
     */
    gsQuaternion<T> operator>>(const gsQuaternion<T> & other) const;
    gsQuaternion<T>& operator>>=(const gsQuaternion<T> & other);

    /**
     * @brief Scalar multiplication: operator for scaling the quaternion by a scalar value as
     * q * scalar
     */
    gsQuaternion<T> operator*(T scalar) const;
    gsQuaternion<T>& operator*=(T scalar);

    /**
     * @brief Operator for element-wise devision
     * 
     */
    gsQuaternion<T> operator/(const gsQuaternion<T> & other) const;
    gsQuaternion<T>& operator/=(const gsQuaternion<T> & other);

    /**
     * @brief Operator for element-wise division by a scalar
     * 
     */
    gsQuaternion<T> operator/(T scalar) const;
    gsQuaternion<T>& operator/=(T scalar);

    /**
     * @brief Operator for element-wise comparison
     * 
     */
    bool operator<=(const gsQuaternion<T> & other) const;
    bool operator>=(const gsQuaternion<T> & other) const;
    bool operator<(const gsQuaternion<T> & other) const;
    bool operator>(const gsQuaternion<T> & other) const;
    // bool operator==(const gsQuaternion<T> & other) const;
    bool operator!=(const gsQuaternion<T> & other) const;

    //------------------------------Functions-----------------------------------
    /**
     * @brief Set *this quaternion to the sum of A and B: this = A + B
     * 
     */
    void add(const gsQuaternion<T> & A, const gsQuaternion<T> & B);

    /**
     * @brief Set *this quaternion to the difference of A and B: this = A - B
     * 
     */
    void sub(const gsQuaternion<T> & A, const gsQuaternion<T> & B);

    /**
     * @brief Cross product of two quaternions
     * NOTE!!!: The quaternion product is not commutative
     */
    void cross(const gsQuaternion<T> & A, const gsQuaternion<T> & B);

    /**
     * @brief Scalar multiplication: this = q * scalar
     * 
     */
    void mul(const gsQuaternion<T> & A, T scalar);
    /**
     * @brief dot product of this quaternion with another quaternion
     * 
     */
    T dot(const gsQuaternion<T> & other) const;

    /**
     * @brief Scale this quaternion by a scalar: this *=s
     * 
     */
    void scale(T scalar);

    /**
     * @brief normalize this quaternion in place, so the euclidean length is 1
     * Return false if the original quaternion has zero length
     * Set it to [1,0,0,0] and return true otherwise
     */
    bool normalize();

    /**
     * @brief Return the normalized quaternion
     * 
     */
    gsQuaternion<T> GetNormalized() const;

    /**
     * @brief conjugate of the quaternion
     * 
     */

    void Conjugate(const gsQuaternion<T> & q);

    /**
     * @brief Calculate the conjugate in place
     * 
     */
    void Conjugate();

    /**
     * @brief Return a conjugate version of this quaternion
     * 
     */
    gsQuaternion<T> GetConjugate() const;

    /**
     * @brief Inverse of the quaternion
     * 
     */
    gsQuaternion<T> GetInverse() const;

    //----------------------------------Transformations-----------------------------------
    /**
     * @brief Rotate a vector on the basis of this quanternion: ressult = q * [0, A] * q^(-1)
     * 
     */
    gsVector<T> Rotate(const gsVector<T> & A) const;

    /**
     * @brief Rotate back a vector on the basis of this quanternion: ressult = q^(-1) * [0, A] * q
     * 
     */
    gsVector<T> RotateBack(const gsVector<T> & A) const;

    //------------------------------------Conversion-----------------------------------
    /**
     * @brief Set the quaternion from a rotation vector (ie. a 3D axis of rotation with length as angle of rotation)
     * defined in absolute coords.
     * 
     */
    void SetFromRotationVector(const gsVector<T> & vec);

    /**
     * @brief Get the rotation vector (ie. a 3D axis of rotation with length as angle of rotation)
     * 
     */

    gsVector<T> GetRotationVector() const;

    /**
     * @brief Set the quaternion from an angle of rotation and axis
     * 
     */
    void SetFromAngleAxis(T angle, const gsVector<T> & axis);

    /**
     * @brief Set the quaternion from an angle of rotation about x-axis
     * 
     */
    void SetFromAngleX(T angle) 
    { 
        gsVector<T> axis(3); // 创建一个三维向量
        axis(0) = 1;
        axis(1) = 0;
        axis(2) = 0;
        SetFromAngleAxis(angle, axis); 
    }

    /**
     * @brief Set the quaternion from an angle of rotation about y-axis
     * 
     */
    void SetFromAngleY(T angle) 
    { 
        
        gsVector<T> axis(3); // 创建一个三维向量
        axis(0) = 0;
        axis(1) = 1;
        axis(2) = 0;
        SetFromAngleAxis(angle, axis); 
    }

    /**
     * @brief Set the quaternion from an angle of rotation about z-axis
     * 
     */
    void SetFromAngleZ(T angle) 
    { 
        
        gsVector<T> axis(3); // 创建一个三维向量
        axis(0) = 0;
        axis(1) = 0;
        axis(2) = 1;
        SetFromAngleAxis(angle, axis); 
    }

    /**
     * @brief Get the angle of rotation and axis of the quaternion
     * Convert the quaternion to an angle of rotation and an axis, defined in absolute coords.
     */
    void GetAngleAxis(T & angle, gsVector<T> & axis) const;

    /**
     * @brief Set the quaternion from Cardan angles ZYX (Tait-Bryan sequence Z-Y'-X'', intrinsic rotations).
     * 
     */
    void SetCardanZyx(const gsVector<T> & q);

    /**
     * @brief Convert the quaternion to Cardan angles ZYX (Tait-Bryan sequence Z-Y'-X'', intrinsic rotations).
     * 
     */
    gsVector<T> GetCardanZyx() const;

    /**
     * @brief Set the quaternion from Cardan angles XYZ (Tait-Bryan sequence X-Y'-Z'', intrinsic rotations). 
     * 
     */
    void SetCardanXyz(const gsVector<T> & q);

    /**
     * @brief Convert the quaternion to Cardan angles XYZ (Tait-Bryan sequence X-Y'-Z'', intrinsic rotations).
     * 
     */
    gsVector<T> GetCardanXyz() const;

    /**
     * @brief Set the quaternion dq/dt. dt. Inputs: the vector of angular speed w specified in absolute coords,
     * and the rotation expressed as a quaternion q.
     */
    void SetDerivativeAbs(const gsVector<T> & w, const gsQuaternion<T> & q);

    /**
     * @brief Set the quaternion dq/dt. dt. Inputs: the vector of angular speed w specified in relative coords,
     * and the rotation expressed as a quaternion q.
     */
    void SetDerivativeRel(const gsVector<T> & w, const gsQuaternion<T> & q);

    
    /**
     * @brief Compute the vector of angular speed w specified in absolute coords, given the quaternion dq/dt. dt.
     * 
     */
    void GetAngularSpeedAbs(gsVector<T>& w, const gsQuaternion<T> & q);
    

    /**
     * @brief Compute the vector of angular speed w specified in relative coords, given the quaternion dq/dt. dt.
     * 
     */
    void GetAngularSpeedRel(gsVector<T>& w, const gsQuaternion<T> & q);


    /**
     * @brief Set the quaternion ddq/dtdt. Inputs: the vector of angular acceleration 'a' specified
     * in absolute coords, the rotation expressed as a quaternion q, the rotation speed
     * as a quaternion 'q_dt'.
     */

    void SetDt2FromAngAccAbs(const gsVector<T>& a, const gsQuaternion<T>& q, const gsQuaternion<T>& q_dt);

    /**
     * @brief Set the quaternion ddq/dtdt. Inputs: the vector of angular acceleration 'a' specified
     * in relative coords, the rotation expressed as a quaternion q, the rotation speed
     * as a quaternion 'q_dt'.
     */
    void SetDt2FromAngAccRel(const gsVector<T>& a, const gsQuaternion<T>& q, const gsQuaternion<T>& q_dt);

    /**
     * @brief Set the Dt From Angle Axis object
     * 
     * @param q 
     * @param angle_dt 
     * @param axis 
     */
    void SetDtFromAngleAxis(const gsQuaternion<T>& q, T angle_dt, const gsQuaternion<T>& axis);

    /**
     * @brief Set the Dt 2 From Angle Axis object
     * 
     * @param q 
     * @param q_dt 
     * @param angle_dtdt 
     * @param axis 
     */
    void SetDt2FromAngleAxis(const gsQuaternion<T>& q,
                             const gsQuaternion<T>& q_dt,
                             T angle_dtdt,
                             const gsQuaternion<T>& axis);

    

private:
    /// Stores all mesh names (might be useful later)
    gsVector<T> m_data = gsVector<T>(4);

    template<class U>
    friend class gsQuaternion;
};



//------------------------------Implementation-----------------------------------

//Default constructor
template<class T>
inline gsQuaternion<T>::gsQuaternion()
{
    // m_data.row(0) = 0;
    // m_data.row(1) = 0;
    // m_data.row(2) = 0;
    // m_data.row(3) = 0;
    m_data << 0, 0, 0, 0;
}

//Constructor from four scalars
template<class T>
inline gsQuaternion<T>::gsQuaternion(T e0, T e1, T e2, T e3)
{
    m_data << e0, e1, e2, e3;
}

//Constructor from a vector and scalar on the real part
template<class T>
inline gsQuaternion<T>::gsQuaternion(T scalar, const gsVector<T> & vec)
{
    m_data(0) =  scalar;
    m_data(1) = vec(0);
    m_data(2) = vec(1);
    m_data(3) = vec(2);
}
//Copy constructor
template<class T>
inline gsQuaternion<T>::gsQuaternion(const gsQuaternion<T> & other)
{
    m_data[0] = other.m_data[0];
    m_data[1] = other.m_data[1];
    m_data[2] = other.m_data[2];
    m_data[3] = other.m_data[3];
}

//Copy constructor from a quaternion of different type
template<class T>
template<class U>
inline gsQuaternion<T>::gsQuaternion(const gsQuaternion<U> & other)
{
    m_data[0] = static_cast<T>(other.m_data[0]);
    m_data[1] = static_cast<T>(other.m_data[1]);
    m_data[2] = static_cast<T>(other.m_data[2]);
    m_data[3] = static_cast<T>(other.m_data[3]);
}

//-----------------------------Setter and Getter-----------------------------------
template<class T>
inline T& gsQuaternion<T>::operator[](unsigned index)
{
    GISMO_ASSERT(index < 4, "Index out of bounds");
    return m_data[index];
}

template<class T>
inline const T& gsQuaternion<T>::operator[](unsigned index) const
{
    GISMO_ASSERT(index < 4, "Index out of bounds");
    return m_data[index];
}

//-----------------------Assignment Operator-----------------------------------
template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator=(const gsQuaternion<T> & other)
{
    m_data[0] = other.m_data[0];
    m_data[1] = other.m_data[1];
    m_data[2] = other.m_data[2];
    m_data[3] = other.m_data[3];
    return *this;
}

//-----------------------Operator Overloading-----------------------------------
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator+() const
{
    return *this;   
}

template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator-() const
{
    return gsQuaternion<T>(-m_data[0], -m_data[1], -m_data[2], -m_data[3]);
}

// Conjugate operator
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator!() const
{
    return gsQuaternion<T>(m_data[0], -m_data[1], -m_data[2], -m_data[3]);
}

// Quaternion sum
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator+(const gsQuaternion<T> & other) const
{
    return gsQuaternion<T>(m_data[0] + other.m_data[0], m_data[1] + other.m_data[1], m_data[2] + other.m_data[2], m_data[3] + other.m_data[3]);
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator+=(const gsQuaternion<T> & other)
{
    m_data[0] += other.m_data[0];
    m_data[1] += other.m_data[1];
    m_data[2] += other.m_data[2];
    m_data[3] += other.m_data[3];
    return *this;
}

// Quaternion difference
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator-(const gsQuaternion<T> & other) const
{
    return gsQuaternion<T>(m_data[0] - other.m_data[0], m_data[1] - other.m_data[1], m_data[2] - other.m_data[2], m_data[3] - other.m_data[3]); 
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator-=(const gsQuaternion<T> & other)
{
    m_data[0] -= other.m_data[0];
    m_data[1] -= other.m_data[1];
    m_data[2] -= other.m_data[2];
    m_data[3] -= other.m_data[3];
    return *this;
}

// Quaternion product
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator*(const gsQuaternion<T> & other) const
{
   gsQuaternion<T> result;
   result.cross(*this, other);
   return result;
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator*=(const gsQuaternion<T> & other)
{
    this->cross(*this, other);
    return *this;
}

// Operator >> for specular quaternion product
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator>>(const gsQuaternion<T> & other) const
{
    gsQuaternion<T> result;
    result.cross(other, *this);
    return result;
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator>>=(const gsQuaternion<T> & other)
{
    this->cross(other, *this);
    return *this;
}

// Scalar multiplication
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator*(T scalar) const
{
    return gsQuaternion<T>(m_data[0] * scalar, m_data[1] * scalar, m_data[2] * scalar, m_data[3] * scalar);
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator*=(T scalar)
{
    m_data[0] *= scalar;
    m_data[1] *= scalar;
    m_data[2] *= scalar;
    m_data[3] *= scalar;
    return *this;
}

// Element-wise division
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator/(const gsQuaternion<T> & other) const
{
    return gsQuaternion<T>(m_data[0] / other.m_data[0], m_data[1] / other.m_data[1], m_data[2] / other.m_data[2], m_data[3] / other.m_data[3]);
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator/=(const gsQuaternion<T> & other)
{
    m_data[0] /= other.m_data[0];
    m_data[1] /= other.m_data[1];
    m_data[2] /= other.m_data[2];
    m_data[3] /= other.m_data[3];
    return *this;
}

// Element-wise division by a scalar
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::operator/(T scalar) const
{
    return gsQuaternion<T>(m_data[0] / scalar, m_data[1] / scalar, m_data[2] / scalar, m_data[3] / scalar);
}

template<class T>
inline gsQuaternion<T>& gsQuaternion<T>::operator/=(T scalar)
{
    m_data[0] /= scalar;
    m_data[1] /= scalar;
    m_data[2] /= scalar;
    m_data[3] /= scalar;
    return *this;
}

//------------------------------Comparison-----------------------------------
template<class T>
inline bool gsQuaternion<T>::operator<=(const gsQuaternion<T> & other) const
{
    return m_data[0] <= other.m_data[0] && m_data[1] <= other.m_data[1] && m_data[2] <= other.m_data[2] && m_data[3] <= other.m_data[3];
}

template<class T>
inline bool gsQuaternion<T>::operator>=(const gsQuaternion<T> & other) const
{
    return m_data[0] >= other.m_data[0] && m_data[1] >= other.m_data[1] && m_data[2] >= other.m_data[2] && m_data[3] >= other.m_data[3];    
}

template<class T>
inline bool gsQuaternion<T>::operator<(const gsQuaternion<T> & other) const
{
    return m_data[0] < other.m_data[0] && m_data[1] < other.m_data[1] && m_data[2] < other.m_data[2] && m_data[3] < other.m_data[3];
}

template<class T>
inline bool gsQuaternion<T>::operator>(const gsQuaternion<T> & other) const
{
    return m_data[0] > other.m_data[0] && m_data[1] > other.m_data[1] && m_data[2] > other.m_data[2] && m_data[3] > other.m_data[3];
}

template<class T>
inline bool gsQuaternion<T>::operator==(const gsQuaternion<T> & other) const
{
    return m_data[0] == other.m_data[0] && m_data[1] == other.m_data[1] && m_data[2] == other.m_data[2] && m_data[3] == other.m_data[3];
}

template<class T>
inline bool gsQuaternion<T>::operator!=(const gsQuaternion<T> &other ) const
{
    return !(*this == other);
}

//-------------------------------------Getter and Setter-----------------------------------
template<class T>
inline void gsQuaternion<T>::Set(T e0, T e1, T e2, T e3)
{
    m_data[0] = e0;
    m_data[1] = e1;
    m_data[2] = e2;
    m_data[3] = e3;
}

// Set from another object
template<class T>
inline void gsQuaternion<T>::Set(const gsQuaternion<T>& q)
{
    m_data[0] = q.m_data[0];
    m_data[1] = q.m_data[1];
    m_data[2] = q.m_data[2];
    m_data[3] = q.m_data[3];
}

// Set from a scalar
template<class T>
inline void gsQuaternion<T>::Set(T scalar)
{
    m_data[0] = scalar;
    m_data[1] = scalar;
    m_data[2] = scalar;
    m_data[3] = scalar;
}

// Set null
template<class T>
inline void gsQuaternion<T>::SetNull()
{
    m_data[0] = 0;
    m_data[1] = 0;
    m_data[2] = 0;
    m_data[3] = 0;
}

// Set unit
template<class T>
inline void gsQuaternion<T>::SetUnit()
{
    m_data[0] = 1;
    m_data[1] = 0;
    m_data[2] = 0;
    m_data[3] = 0;
}

// Set the scalar part only
template<class T>
inline void gsQuaternion<T>::SetScalar(T scalar)
{
    m_data[0] = scalar;
}

// Set the vector part only
template<class T>
inline void gsQuaternion<T>::SetVector(const gsVector<T> & vec)
{
    m_data[1] = vec.row(0);
    m_data[2] = vec.row(1);
    m_data[3] = vec.row(2);
}

// Return true if the quaternion is null
template<class T>
inline bool gsQuaternion<T>::IsNull() const
{
    return m_data[0] == 0 && m_data[1] == 0 && m_data[2] == 0 && m_data[3] == 0;
}

// Return true if the quaternion is identity
template<class T>
inline bool gsQuaternion<T>::IsIdentity() const
{
    return m_data[0] == 1 && m_data[1] == 0 && m_data[2] == 0 && m_data[3] == 0;
}

// Return true if the quaternion equals to another quaternion
template<class T>
inline bool gsQuaternion<T>::Equals(const gsQuaternion<T> & other) const
{
    return m_data[0] == other.m_data[0] && m_data[1] == other.m_data[1] && m_data[2] == other.m_data[2] && m_data[3] == other.m_data[3];
}

template<class T>
inline bool gsQuaternion<T>::Equals(const gsQuaternion<T> & other, T tol) const
{
    return std::abs(m_data[0] - other.m_data[0]) <= tol && std::abs(m_data[1] - other.m_data[1]) <= tol && std::abs(m_data[2] - other.m_data[2]) <= tol && std::abs(m_data[3] - other.m_data[3]) <= tol;
}

template<class T>
inline gsVector<T> gsQuaternion<T>::GetVector() const
{
    gsVector<T> ret(3);  // 假设 gsVector<T>(3) 构造出一个 3 维向量
    ret(0) = m_data[1];
    ret(1) = m_data[2];
    ret(2) = m_data[3];
    return ret;
}

template<class T>
inline gsVector<T> gsQuaternion<T>::GetXAxis() const
{
    return gsVector<T>((m_data[0] * m_data[0] + m_data[1] * m_data[1]) * 2 - 1,
                       (m_data[1] * m_data[2] + m_data[0] * m_data[3]) * 2,
                       (m_data[1] * m_data[3] - m_data[0] * m_data[2]) * 2);
}

template<class T>
inline gsVector<T> gsQuaternion<T>::GetYAxis() const
{
    return gsVector<T>((m_data[1] * m_data[2] - m_data[0] * m_data[3]) * 2,
                       (m_data[0] * m_data[0] + m_data[2] * m_data[2]) * 2 - 1,
                       (m_data[2] * m_data[3] + m_data[0] * m_data[1]) * 2);
}

template<class T>
inline gsVector<T> gsQuaternion<T>::GetZAxis() const
{
    return gsVector<T>((m_data[1] * m_data[3] + m_data[0] * m_data[2]) * 2,
                       (m_data[2] * m_data[3] - m_data[0] * m_data[1]) * 2,
                       (m_data[0] * m_data[0] + m_data[3] * m_data[3]) * 2 - 1);
}


//------------------------------Norms-----------------------------------
template<class T>
inline T gsQuaternion<T>::norm2() const
{
    return this->dot(*this);
}

template<class T>
inline T gsQuaternion<T>::norm() const
{
    return std::sqrt(this->norm2());
}

template<class T>
inline T gsQuaternion<T>::normInf() const
{
    T e0e1 = std::max(std::abs(m_data[0]), std::abs(m_data[1]));
    T e0e1e2 = std::max(e0e1, std::abs(m_data[2]));
    return std::max(e0e1e2, std::abs(m_data[3]));
}

//------------------------------Functions-----------------------------------
template<class T>
inline void gsQuaternion<T>::add(const gsQuaternion<T> & A, const gsQuaternion<T> & B)
{
    m_data[0] = A.m_data[0] + B.m_data[0];
    m_data[1] = A.m_data[1] + B.m_data[1];
    m_data[2] = A.m_data[2] + B.m_data[2];
    m_data[3] = A.m_data[3] + B.m_data[3];
}

template<class T>
inline void gsQuaternion<T>::sub(const gsQuaternion<T> & A, const gsQuaternion<T> & B)
{
    m_data[0] = A.m_data[0] - B.m_data[0];
    m_data[1] = A.m_data[1] - B.m_data[1];
    m_data[2] = A.m_data[2] - B.m_data[2];
    m_data[3] = A.m_data[3] - B.m_data[3];
}

template <class T>
inline void gsQuaternion<T>::cross(const gsQuaternion<T>& A, const gsQuaternion<T>& B) 
{
    T w = A.m_data[0] * B.m_data[0] - A.m_data[1] * B.m_data[1] - A.m_data[2] * B.m_data[2] -
             A.m_data[3] * B.m_data[3];
    T x = A.m_data[0] * B.m_data[1] + A.m_data[1] * B.m_data[0] - A.m_data[3] * B.m_data[2] +
             A.m_data[2] * B.m_data[3];
    T y = A.m_data[0] * B.m_data[2] + A.m_data[2] * B.m_data[0] + A.m_data[3] * B.m_data[1] -
             A.m_data[1] * B.m_data[3];
    T z = A.m_data[0] * B.m_data[3] + A.m_data[3] * B.m_data[0] - A.m_data[2] * B.m_data[1] +
             A.m_data[1] * B.m_data[2];
    m_data[0] = w;
    m_data[1] = x;
    m_data[2] = y;
    m_data[3] = z;
}

template<class T>
inline T gsQuaternion<T>::dot(const gsQuaternion<T> & other) const
{
    return m_data[0] * other.m_data[0] + m_data[1] * other.m_data[1] + m_data[2] * other.m_data[2] + m_data[3] * other.m_data[3];
}

template<class T>
inline void gsQuaternion<T>::mul(const gsQuaternion<T> & A, T scalar)
{
    m_data[0] = A.m_data[0] * scalar;
    m_data[1] = A.m_data[1] * scalar;
    m_data[2] = A.m_data[2] * scalar;
    m_data[3] = A.m_data[3] * scalar;
}

template<class T>
inline void gsQuaternion<T>::scale(T scalar)
{
    m_data[0] *= scalar;
    m_data[1] *= scalar;
    m_data[2] *= scalar;
    m_data[3] *= scalar;
}

template<class T>
inline bool gsQuaternion<T>::normalize()
{
    T length = this->norm();
    if (length < std::numeric_limits<T>::min()) {
        m_data[0] = 1;
        m_data[1] = 0;
        m_data[2] = 0;
        m_data[3] = 0;
        return false;
    }
    this->scale(1 / length);
    return true;
}

template<class T>
inline gsQuaternion<T> gsQuaternion<T>::GetNormalized() const
{
    gsQuaternion<T> result(*this);
    result.normalize();
    return result;
}

template<class T>
inline void gsQuaternion<T>::Conjugate(const gsQuaternion<T> & q)
{
    m_data[0] = q.m_data[0];
    m_data[1] = -q.m_data[1];
    m_data[2] = -q.m_data[2];
    m_data[3] = -q.m_data[3];
}

template<class T>
inline void gsQuaternion<T>::Conjugate()
{
   Conjugate(*this);
}

template<class T>
inline gsQuaternion<T> gsQuaternion<T>::GetConjugate() const
{
    return gsQuaternion<T>(m_data[0], -m_data[1], -m_data[2], -m_data[3]);
}
template<class T>
inline gsQuaternion<T> gsQuaternion<T>::GetInverse() const
{
    gsQuaternion<T> conj = this->GetConjugate();
    conj.scale(1 / this->norm2());
    return conj;
}

template<class T>
inline gsVector<T> gsQuaternion<T>::Rotate(const gsVector<T> & A) const
{
    T r0r0 = m_data[0] * m_data[0];
    T r1r1 = m_data[1] * m_data[1];
    T r2r2 = m_data[2] * m_data[2];
    T r3r3 = m_data[3] * m_data[3];
    T r0r1 = m_data[0] * m_data[1];
    T r0r2 = m_data[0] * m_data[2];
    T r0r3 = m_data[0] * m_data[3];
    T r1r2 = m_data[1] * m_data[2];
    T r1r3 = m_data[1] * m_data[3];
    T r2r3 = m_data[2] * m_data[3];

    T x = (r0r0 + r1r1 - r2r2 - r3r3) * A(0) + 2 * (r1r2 - r0r3) * A(1) + 2 * (r1r3 + r0r2) * A(2);
    T y = 2 * (r0r3 + r1r2) * A(0) + (r0r0 - r1r1 + r2r2 - r3r3) * A(1) + 2 * (r2r3 - r0r1) * A(2);
    T z = 2 * (r1r3 - r0r2) * A(0) + 2 * (r2r3 + r0r1) * A(1) + (r0r0 - r1r1 - r2r2 + r3r3) * A(2);

    gsVector<T> res(3);
    res << x, y, z;
    return res;
    
}

template<class T>
inline gsVector<T> gsQuaternion<T>::RotateBack(const gsVector<T> & A) const
{
    return this->GetConjugate().Rotate(A);
}

template<class T>
inline void gsQuaternion<T>::SetFromRotationVector(const gsVector<T> & angle_axis)
{
    T theta_squared = angle_axis.norm();
    if (theta_squared > 1e-30) {
        T theta = std::sqrt(theta_squared);
        T half_theta = theta / 2;
        T k = std::sin(half_theta) / theta;
        m_data(0) = std::cos(half_theta);
        m_data(1) = angle_axis(0) * k;
        m_data(2) = angle_axis(1) * k;
        m_data(3) = angle_axis(2) * k;
    } else {
        // For almost zero rotation:
        T k(0.5);
        m_data(0) = T(1.0);
        m_data(1) = angle_axis(0) * k;
        m_data(2) = angle_axis(1) * k;
        m_data(3) = angle_axis(2) * k;
    }
}

template<class T>
inline gsVector<T> gsQuaternion<T>::GetRotationVector() const
{
    gsVector<T> angle_axis;
    T sin_squared = m_data(1) * m_data(1) + m_data(2) * m_data(2) + m_data(3) * m_data(3);
    // For non-zero rotation
    if (sin_squared > 1e-30) 
    {
        T sin_theta = std::sqrt(sin_squared);
        T k = 2 * std::atan2(sin_theta, m_data(0)) / sin_theta;
        angle_axis(0) = m_data(1) * k;
        angle_axis(1) = m_data(2) * k;
        angle_axis(2) = m_data(3) * k;
    } else {
        // For almost zero rotation
        T k(2.0);
        angle_axis(0) = m_data(1) * k;
        angle_axis(1) = m_data(2) * k;
        angle_axis(2) = m_data(3) * k;
    }
    return angle_axis;
}

template<class T>
inline void gsQuaternion<T>::SetFromAngleAxis(T angle, const gsVector<T> & axis)
{
    T half_angle = angle / 2;
    T sin_half_angle = std::sin(half_angle);
    m_data(0) = std::cos(half_angle);
    m_data(1) = axis(0) * sin_half_angle;
    m_data(2) = axis(1) * sin_half_angle;
    m_data(3) = axis(2) * sin_half_angle;
}

template<class T>
inline void gsQuaternion<T>::GetAngleAxis(T & angle, gsVector<T> & axis) const
{
    T sin_squared = m_data(1) * m_data(1) + m_data(2) * m_data(2) + m_data(3) * m_data(3);
    // For non-zero rotation
    if (sin_squared > 1e-30) 
    {
        T sin_theta = std::sqrt(sin_squared);
        angle = 2 * std::atan2(sin_theta, m_data(0));
        T k = 1 / sin_theta;
        axis(0) = m_data(1) * k;
        axis(1) = m_data(2) * k;
        axis(2) = m_data(3) * k;
        axis.normalize();
    } else {
        // For almost zero rotation
        angle = 0.0;
        axis(0) = 1; // m_data(1) * 2.0;
        axis(1) = 0;
        axis(2) = 0;
    }

    // Ensure that angle is always in [-PI, PI] range
    if (angle > M_PI) 
    {
        angle -= 2 * M_PI;
    } else if (angle < -M_PI) 
    {
        angle += 2 * M_PI;
    }
}

/**
 * @brief  Set the quaternion from Cardan angles XYZ (Tait-Bryan sequence X-Y'-Z'', intrinsic rotations).
 * 它认为先绕 X 轴旋转，再绕 Y 轴，最后绕 Z 轴
 * https://en.wikipedia.org/wiki/Euler_angles
 * @tparam T 
 * @param angles 
 */

template<class T>
inline void gsQuaternion<T>::SetCardanXyz(const gsVector<T> & angles)
{
    T t0 = std::cos(angles(2) * T(0.5));
    T t1 = std::sin(angles(2) * T(0.5));
    T t2 = std::cos(angles(0) * T(0.5));
    T t3 = std::sin(angles(0) * T(0.5));
    T t4 = std::cos(angles(1) * T(0.5));
    T t5 = std::sin(angles(1) * T(0.5));

    m_data(0) = t0 * t2 * t4 + t1 * t3 * t5;
    m_data(1) = t0 * t3 * t4 - t1 * t2 * t5;
    m_data(2) = t0 * t2 * t5 + t1 * t3 * t4;
    m_data(3) = t1 * t2 * t4 - t0 * t3 * t5;
}


template<class T>
inline gsVector<T> gsQuaternion<T>::GetCardanXyz() const
{
    gsVector<T> euler(3);
    T t0 = m_data(0) * m_data(0);
    T t1 = m_data(1) * m_data(1);
    T t2 = m_data(2) * m_data(2);
    T t3 = m_data(3) * m_data(3);

    // roll
    euler(0) = std::atan2(2 * (m_data(2) * m_data(3) + m_data(0) * m_data(1)), t3 - t2 - t1 + t0);
    // pitch
    euler(1) = -std::asin(2 * (m_data(1) * m_data(3) - m_data(0) * m_data(2)));
    // yaw
    euler(2) = std::atan2(2 * (m_data(1) * m_data(2) + m_data(3) * m_data(0)), t1 + t0 - t3 - t2);
    return euler;

}

template<class T>
inline void gsQuaternion<T>::SetCardanZyx(const gsVector<T> & angles)
{
    T c1 = std::cos(angles(2) / 2);
    T s1 = std::sin(angles(2) / 2);
    T c2 = std::cos(angles(0) / 2);
    T s2 = std::sin(angles(0) / 2);
    T c3 = std::cos(angles(1) / 2);
    T s3 = std::sin(angles(1) / 2);

    T c1c2 = c1 * c2;
    T s1s2 = s1 * s2;

    m_data(0) = c1c2 * c3 + s1s2 * s3;
    m_data(1) = c1c2 * s3 - s1s2 * c3;
    m_data(2) = c1 * s2 * c3 + s1 * c2 * s3;
    m_data(3) = s1 * c2 * c3 - c1 * s2 * s3;
}

template<class T>
inline gsVector<T> gsQuaternion<T>::GetCardanZyx() const
{
    gsVector<T> euler(3);
    T t0 = m_data(0) * m_data(0);
    T t1 = m_data(1) * m_data(1);
    T t2 = m_data(2) * m_data(2);
    T t3 = m_data(3) * m_data(3);


    euler(0) = std::asin(-2 * (m_data(1) * m_data(3) - m_data(2) * m_data(0)));

    euler(1) = std::atan2(2 * (m_data(1) * m_data(2) + m_data(3) * m_data(0)), t0 - t1 - t2 + t3);

    euler(2) = std::atan2(2 * (m_data(1) * m_data(2) + m_data(3) * m_data(0)), t0 + t1 - t2 - t3);
    return euler;
}

template<class T>
inline void gsQuaternion<T>::SetDerivativeAbs(const gsVector<T> & w, const gsQuaternion<T> & q)
{
    gsQuaternion<T> q_dt(0,w);
    *this = q_dt * q * T(0.5); // q_dt = 0.5 * {0, w} * q
}

template<class T>
inline void gsQuaternion<T>::SetDerivativeRel(const gsVector<T> & w, const gsQuaternion<T> & q)
{
    gsQuaternion<T> q_dt(0,w);
    *this = q * q_dt * T(0.5); // q_dt = 0.5 * q * {0, w}
}

template<class T>
inline void gsQuaternion<T>::GetAngularSpeedAbs(gsVector<T>& w, const gsQuaternion<T> & q)
{
    gsQuaternion<T> q_dt;
    q_dt.cross(*this, q.GetConjugate());
    w = q_dt.GetVector();
    w = w * T(2);
}

template<class T>
inline void gsQuaternion<T>::GetAngularSpeedRel(gsVector<T>& w, const gsQuaternion<T> & q)
{
    gsQuaternion<T> q_dt;
    q_dt.cross(q.GetConjugate(), *this);
    w = q_dt.GetVector();
    w = w * T(2);
}

/**
 * @brief 绝对版本（Abs）：适用于 a 在全局坐标系下的情况。
 * 
 * @tparam T 
 * @param a 
 * @param q 
 * @param q_dt 
 */

 template<class T>
 inline void gsQuaternion<T>::SetDt2FromAngAccAbs(const gsVector<T>& a,
                                                    const gsQuaternion<T>& q,
                                                    const gsQuaternion<T>& q_dt)
 {
     // 第一项：0.5 * {0, a} * q  
     // 其中 {0, a} 构造方式为 gsQuaternion(0, a)
     gsQuaternion<T> term1(0, a);
     
     // 求 q 的共轭（对于单位四元数，逆等于共轭）
     gsQuaternion<T> q_inv = q.GetConjugate();
     
     // 由 q_dt = 0.5 * {0, w} * q 可得 {0, w} = 2 * (q_dt * q_inv)
     gsQuaternion<T> w_quat = (q_dt * q_inv) * T(2);
     
     // 第二项：0.5 * {0, w} * q_dt
     gsQuaternion<T> term2 = (w_quat * q_dt) * T(0.5);
     
     // 最终，四元数的二阶导数为：
     *this = (term1 * q + term2) * T(0.5);
 }
/**
 * @brief 相对版本（Rel）：适用于 a 在局部坐标系（随物体旋转的参考系）下的情况。
 * 
 * @tparam T 
 * @param a 
 * @param q 
 * @param q_dt 
 */

 template<class T>
 inline void gsQuaternion<T>::SetDt2FromAngAccRel(const gsVector<T>& a,
                                                    const gsQuaternion<T>& q,
                                                    const gsQuaternion<T>& q_dt)
 {
     // 第一项：0.5 * q * {0, a}
     gsQuaternion<T> term1(0, a);
     
     // 由 q_dt = 0.5 * q * {0, w_rel} 可得 {0, w_rel} = 2 * (q.GetConjugate() * q_dt)
     gsQuaternion<T> q_inv = q.GetConjugate();
     gsQuaternion<T> w_rel_quat = (q_inv * q_dt) * T(2);
     
     // 第二项：0.5 * q_dt * {0, w_rel}
     gsQuaternion<T> term2 = (q_dt * w_rel_quat) * T(0.5);
     
     // 最后，q 的二阶导数为：
     *this = (q * term1 + term2) * T(0.5);
 }

template<class T>
inline void gsQuaternion<T>::SetDtFromAngleAxis(const gsQuaternion<T>& q, T angle_dt, const gsQuaternion<T>& axis)
{
    this->SetDerivativeAbs(angle_dt * axis, q);
}

template<class T>
inline void gsQuaternion<T>::SetDt2FromAngleAxis(const gsQuaternion<T>& q,
                                                 const gsQuaternion<T>& q_dt,
                                                 T angle_dtdt,
                                                 const gsQuaternion<T>& axis)

{
    this->SetDtFromAngleAxis(angle_dtdt * axis, q, q_dt);
}

} //namespace gismo
