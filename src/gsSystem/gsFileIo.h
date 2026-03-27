/** @file gsModuleManager.h

    @brief The module manager for dynamically loaded modules for G+Smo.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsSystem/gsModule.h>

class gsFileIo : public gsModule
{
public:
    virtual std::string extension() = 0;

    virtual void read(std::string fname) {GISMO_NO_IMPLEMENTATION }

    virtual void write(std::string fname) {GISMO_NO_IMPLEMENTATION }
};

