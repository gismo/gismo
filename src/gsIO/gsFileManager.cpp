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
#include <gsUtils/gsUtils.h> 

#if defined _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

namespace gismo
{

class gsFileManagerData;

// Use a singleton to store data:
class GISMO_EXPORT gsFileManagerData
{
public:

    friend gsFileManagerData& gsFileManagerDataSingleton();
    friend class gsFileManager;

    void setSearchPaths(const std::string& paths);

private:
    gsFileManagerData()
    {
#ifdef GISMO_SEARCH_PATHS
        setSearchPaths("" GISMO_SEARCH_PATHS);
#endif
    }

    gsFileManagerData(const gsFileManagerData&);
    gsFileManagerData& operator= (const gsFileManagerData&);
    std::vector<std::string> m_paths;
};

gsFileManagerData& gsFileManagerDataSingleton()
{
    static gsFileManagerData singleton;
    return singleton;
}


bool gsFileManager::fileExists(const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}

char gsFileManager::getLocalPathSeperator()
{
#if defined _WIN32
    return '\\';
#else
    return '/';
#endif
}

bool gsFileManager::isFullyQualified(const std::string& fn)
{
#if defined _WIN32
    return util::starts_with(fn,"/")
        || util::starts_with(fn,"\\")
        || ( fn.size() > 2 && fn[1] == ':' && ( fn[2] == '/' || fn[2] == '\\' ) );
#else
    return util::starts_with(fn,"/");
#endif
}

bool gsFileManager::isRelative(const std::string& fn)
{
#if defined _WIN32
    return util::starts_with(fn,"./")
        || util::starts_with(fn,".\\")
        || util::starts_with(fn,"../")
        || util::starts_with(fn,"..\\");
#else
    return util::starts_with(fn,"./")
        || util::starts_with(fn,"../");
#endif
}

void _replace_slash_by_basckslash(std::string& str)
{
    for ( std::string::iterator it=str.begin(); it!=str.end(); it++ )
        if ( *it=='/' ) *it = '\\';
}

void gsFileManager::setSearchPaths(const std::string& paths)
{
    gsFileManagerDataSingleton().setSearchPaths(paths);
}

void gsFileManagerData::setSearchPaths(const std::string& paths)
{
    m_paths.clear();

    std::string::const_iterator a;
    std::string::const_iterator b = paths.begin();
    while (true)
    {
        a = b;
        while (b != paths.end() && (*b) != ';') { ++b; }

        std::string p(a,b);

#if defined _WIN32
        _replace_slash_by_basckslash(p);
#endif

        if (!p.empty())
        {
#if defined _WIN32
            if (*p.rbegin() != '\\')
                p.push_back('\\');
#else
            if (*p.rbegin() != '/')
                p.push_back('/');
#endif

            m_paths.push_back(p);
        }

        if ( b == paths.end() ) break;

        ++b;
    }
}

std::string gsFileManager::getSearchPaths()
{
    std::string result;
    gsFileManagerData& dat = gsFileManagerDataSingleton();
    for (std::vector<std::string>::const_iterator it = dat.m_paths.begin();
            it < dat.m_paths.end(); ++it)
    {
        result += (*it) + ";";
    }
    return result;
}

bool gsFileManager::find( std::string& fn )
{
#if defined _WIN32
    _replace_slash_by_basckslash(fn);
#endif

    if ( fileExists(fn) ) return true;

    if ( isFullyQualified(fn) || isRelative(fn) ) return false;

    gsFileManagerData& dat = gsFileManagerDataSingleton();

    for (std::vector<std::string>::const_iterator it = dat.m_paths.begin();
            it < dat.m_paths.end(); ++it)
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

bool gsFileManager::mkdir( std::string fn )
{
#if defined _WIN32
    _replace_slash_by_basckslash(fn); 
    return 0!=CreateDirectory(fn.c_str(),NULL);
#else
    return ::mkdir(fn.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}


bool gsFileManager::pathEqual( const std::string& p1, const std::string& p2 )
{
    const size_t sz = p1.size();

    if (sz != p2.size())
        return false;

    for (size_t i=0; i<sz; ++i)
    {
        if (!(
            p1[i] == p2[i]
            || ( p1[i] == '/' && p2[i] == '\\' )
            || ( p1[i] == '\\' && p2[i] == '/' )
        )) return false;

    }
    return true;

}


} //namespace gismo


