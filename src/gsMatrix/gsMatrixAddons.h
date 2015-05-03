/** @file gsMatrixAddons.h

    @brief Provides extra member declarations to the Eigen MatrixBase class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

inline const internal::adjugate_impl<Derived> adjugate() const;

inline void adjugateInPlace();
