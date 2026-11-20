/** @file gsFileIo.cpp

    @brief Implementation of gsFileIo.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
*/

#include <gsSystem/gsFileIo.h>

#include <algorithm>
#include <cctype>

namespace gismo
{

gsFileIo::~gsFileIo()
{ }

bool gsFileIo::canRead() const
{ return false; }

bool gsFileIo::canWrite() const
{ return false; }

bool gsFileIo::read(const std::string & fn, internal::gsXmlTree & data) const
{
    GISMO_UNUSED(fn);
    GISMO_UNUSED(data);
    return false;
}

bool gsFileIo::write(const std::string & fn, const internal::gsXmlTree & data) const
{
    GISMO_UNUSED(fn);
    GISMO_UNUSED(data);
    return false;
}

gsFileIo::uPtr gsFileIo::find(const std::string & ext)
{
    std::string key = ( !ext.empty() && '.' == ext[0] ) ? ext.substr(1) : ext;
    if ( key.empty() )
        return uPtr();
    std::for_each(key.begin(), key.end(),
                  [](char & a){ a = static_cast<char>(::tolower(a)); } );

    // The dynamic modules are loaded when the library is initialized
    // (gsDllMain.cpp); get() registers the built-in ones
    return gsRegistry<gsFileIo>::get().make(key);
}

} // namespace gismo
