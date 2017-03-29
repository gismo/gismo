/** @file gsFileManager.cpp

    @brief Utility class for finding files and handling paths

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gsIO/gsFileManager.h>
#include <iostream>
#include <fstream>
#include <gsCore/gsConfig.h>

namespace gismo
{

gsFileManager& gsFileManagerSingleton()
{
    static gsFileManager singleton;
    return singleton;
}

bool gsFileManager::fileExists(const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}

void gsFileManager::setSearchPaths(const std::string& paths)
{
    m_paths.clear();

    std::string::const_iterator a;
    std::string::const_iterator b = paths.begin();
    while (true)
    {
        a = b;
        while (b != paths.end() && (*b) != ';') { ++b; }

        std::string p(a,b);

        if (!p.empty())
        {
            if (*p.rbegin() != GISMO_PATH_SEPERATOR)
                p.push_back(GISMO_PATH_SEPERATOR);

            m_paths.push_back(p);
        }

        if ( b == paths.end() ) break;

        ++b;
    }
}

std::string gsFileManager::getSearchPaths() const
{
    std::string result;
    for (std::vector<std::string>::const_iterator it = m_paths.begin();
            it < m_paths.end(); ++it)
    {
        result += (*it) + ";";
    }
    return result;
}

bool gsFileManager::find( std::string& fn )
{
    if ( fileExists(fn) ) return true;

    if ( isFullyQualified(fn) || isRelative(fn) ) return false;

    for (std::vector<std::string>::const_iterator it = m_paths.begin();
            it < m_paths.end(); ++it)
    {
        const std::string tmp = (*it) + fn;
        if ( fileExists( tmp ) )
        {
            fn = tmp;
            return true;
        }
    }

    return false;
}

gsFileManager::gsFileManager()
{
#ifdef GISMO_SEARCH_PATHS
    setSearchPaths("" GISMO_SEARCH_PATHS);
#endif
}

} //namespace gismo


