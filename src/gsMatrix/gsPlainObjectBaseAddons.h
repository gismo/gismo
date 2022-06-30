/** @file gsPlainObjectBaseAddons.h

    @brief Provides extra member declarations to the Eigen PlainObjectBase class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

EIGEN_STRONG_INLINE void setPtr(Scalar * ptr)
{
    // Use with Caution for memory leaks! 
    m_storage.data() = ptr;
}
