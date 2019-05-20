/** @file gsFileManager.cpp

    @brief Utility class for finding files and handling paths

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, A. Mantzaflaris, H. Weiner
*/

#include <gsIO/gsFileManager.h>
#include <iostream>
#include <fstream>
#include <gsCore/gsConfig.h>
#include <gsUtils/gsUtils.h>
#include <cstdlib>

#if defined _WIN32
#include <windows.h>
#include <direct.h>
#ifdef __MINGW32__
#include <sys/stat.h>
#endif
#else
#include <sys/stat.h>
#include <dlfcn.h>
#include <unistd.h>
#endif

namespace gismo
{

// Struct for storing data
struct GISMO_EXPORT gsFileManagerData {
    gsFileManagerData();
    std::vector<std::string> m_paths;
    const char* argv0;
};

gsFileManagerData& gsFileManagerDataSingleton()
{
    static gsFileManagerData singleton;
    return singleton;
}

void gsFileManager::setArgv0(const char* argv0)
{
    gsFileManagerDataSingleton().argv0 = argv0;
}

bool gsFileManager::fileExists(const std::string& name)
{
    return !find(name).empty();
}

bool gsFileManager::fileExistsInDataDir(const std::string& name)
{
    return !findInDataDir(name).empty();
}

namespace {
bool _dirExistsWithoutSearching(const std::string& path)
{
    struct stat info;
    return (0==stat(path.c_str(), &info)) && (info.st_mode & S_IFDIR);
}

bool _fileExistsWithoutSearching(const std::string& fn)
{
   // Note:
   //   std::ifstream s(fn.c_str()); return s.good() && ! s.eof();
   // is also possible; however that treats empty files as non-existing.
#if defined _WIN32
    DWORD dwAttrib = GetFileAttributes( fn.c_str() );
    return (dwAttrib != INVALID_FILE_ATTRIBUTES && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
#else
    struct stat buf;
    return ( (0==stat(fn.c_str(), &buf)) && (0!=S_ISREG(buf.st_mode)) );
#endif
}
} // anonymous namespace

const std::string& gsFileManager::getValidPathSeparators()
{
#if defined _WIN32 || defined __CYGWIN__
    static std::string ps("\\/");
#else
    static std::string ps("/");
#endif
    return ps;
}

char gsFileManager::getNativePathSeparator()
{
    return getValidPathSeparators()[0];
}


namespace {
inline bool _contains( const std::string& haystack, char needle )
{
    return haystack.find(needle) != std::string::npos;
}

} // anonymous namespace


bool gsFileManager::isFullyQualified(const std::string& fn)
{
    // TODO: valid
#if defined _WIN32
    return util::starts_with(fn,"/")
        || util::starts_with(fn,"\\")
        || ( fn.size() > 2 && fn[1] == ':' && ( fn[2] == '/' || fn[2] == '\\' ) );
#else
    return util::starts_with(fn,"/");
#endif
}

bool gsFileManager::isExplicitlyRelative(const std::string& fn)
{
    // TODO: valid
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

namespace {

inline void _replace_slash_by_basckslash(std::string& str)
{
    for ( std::string::iterator it=str.begin(); it!=str.end(); it++ )
        if ( *it=='/' ) *it = '\\';
}

inline bool _addSearchPaths(const std::string& in, std::vector<std::string>& out)
{
    bool ok = true;
    std::string p;
    std::string::const_iterator a;
    std::string::const_iterator b = in.begin();
    while (true)
    {
        a = b;
        while (b != in.end() && (*b) != ';') { ++b; }

        p.assign(a,b);

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
            if (_dirExistsWithoutSearching(p))
                out.push_back(p);
            else
                ok = false;
        }

        if ( b == in.end() ) break;

        ++b;
    }
    return ok;
}

} // anonymous namespace

// constructor; registers data from macro
gsFileManagerData::gsFileManagerData()
{
#ifdef GISMO_SEARCH_PATHS
    (void)_addSearchPaths("" GISMO_SEARCH_PATHS, m_paths);
#endif
}

bool gsFileManager::addSearchPaths(const std::string& paths)
{
    return _addSearchPaths( paths, gsFileManagerDataSingleton().m_paths );
}

bool gsFileManager::setSearchPaths(const std::string& paths)
{
    gsFileManagerDataSingleton().m_paths.clear();
    return _addSearchPaths( paths, gsFileManagerDataSingleton().m_paths );
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

std::string gsFileManager::find(std::string fn)
{
#if defined _WIN32
    _replace_slash_by_basckslash(fn);
#endif

    if ( _fileExistsWithoutSearching(fn) ) return fn;

    if ( isFullyQualified(fn) || isExplicitlyRelative(fn) ) return std::string();

    gsFileManagerData& dat = gsFileManagerDataSingleton();

    std::string tmp;
    for (std::vector<std::string>::const_iterator it = dat.m_paths.begin();
            it < dat.m_paths.end(); ++it)
    {
        tmp = (*it) + fn;
        if ( _fileExistsWithoutSearching(tmp) )
            return tmp;
    }

    return std::string();
}

std::string gsFileManager::findInDataDir(std::string fn)
{
#if defined _WIN32
    _replace_slash_by_basckslash(fn);
#endif

    // We know that GISMO_DATA_DIR ends with a path seperator, but
    // maybe the user does not know it.
    if ( fn[0] == '/' || fn[0] == '\\' ) fn.erase(0,1);

    std::string fn_out = GISMO_DATA_DIR + fn;

    if ( _fileExistsWithoutSearching(fn_out) ) return fn_out;

    return std::string();
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

std::string gsFileManager::getTempPath()
{
#if defined _WIN32
    TCHAR _temp[MAX_PATH];
    DWORD l = GetTempPath(/*length of buffer:*/MAX_PATH, _temp);
    GISMO_UNUSED(l);
    GISMO_ASSERT(l, "GetTempPath did return 0");
    return std::string(_temp);
#else

    // Typically, we should consider TMPDIR
    //   http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap08.html#tag_08_03
    //   https://en.wikipedia.org/wiki/TMPDIR&oldid=728654758
    char * _temp = getenv ("TMPDIR");
    // getenv returns NULL ptr if the variable is unknown (http://en.cppreference.com/w/cpp/utility/program/getenv).
    // If it is an empty string, we should also exclude it.
    if (_temp != NULL && _temp[0] != '\0')
    {
        // note: env variable needs no free
        return std::string(_temp);
    }

    // Okey, if first choice did not work, try this:
    _temp = getenv("TEMP");
    if (_temp != NULL && _temp[0] != '\0')
    {
        // note: env variable needs no free
        return std::string(_temp);
    }

    // And as third choice, use just current directory
    // http://man7.org/linux/man-pages/man2/getcwd.2.html
    _temp = getcwd(NULL, 0);
    GISMO_ASSERT(NULL!=_temp, "getcwd returned NULL.");
    std::string path(_temp);
    // The string is allocated using malloc, see the reference above
    std::free(_temp);
    return path;
#endif
}

std::string gsFileManager::getCurrentPath()
{
#if defined _WIN32
    TCHAR _temp[MAX_PATH];
    DWORD l = GetCurrentDirectory(/*length of buffer:*/MAX_PATH, _temp);
    GISMO_UNUSED(l);
    GISMO_ASSERT(l, "GetCurrentDirectory did return 0");
    return std::string(_temp);
#else
    // http://man7.org/linux/man-pages/man2/getcwd.2.html
    char* _temp = getcwd(NULL, 0);
    GISMO_ASSERT(NULL!=_temp, "getcwd returned NULL.");
    std::string path(_temp);
    // The string is allocated using malloc, see the reference above
    std::free(_temp);
    return path;
#endif
}

std::string gsFileManager::getExePath()
{
#if defined _WIN32
    TCHAR _temp[MAX_PATH];
    DWORD l = GetModuleFileName( NULL, _temp, MAX_PATH );
    GISMO_UNUSED(l);
    GISMO_ASSERT(l, "GetModuleFileName did return 0");
    GISMO_ASSERT(_fileExistsWithoutSearching(_temp),
        "The executable cannot be found where it is expected." );
    return getCanonicRepresentation( std::string(_temp) + "/../" );
#else
    const char* argv0 = gsFileManagerDataSingleton().argv0;
    if (!argv0)
    {
        gsWarn << "gsCmdLine::getValues has not been called. Therefore, "
             "the path is not available.\n";
        return getCurrentPath();
    }

    GISMO_ASSERT( _fileExistsWithoutSearching( isFullyQualified( argv0 )
        ? getCanonicRepresentation( std::string(argv0) )
        : getCanonicRepresentation( getCurrentPath() + "/" + argv0 ) ),
        "The executable cannot be found where it is expected." );

    return isFullyQualified( argv0 )
        ? getCanonicRepresentation( std::string(argv0) + "/../" )
        : getCanonicRepresentation( getCurrentPath() + "/" + argv0 + "/../" );
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

std::string gsFileManager::getExtension(std::string const & fn)
{
    if(fn.find_last_of(".") != std::string::npos)
    {
        std::string ext = fn.substr(fn.rfind(".")+1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        return ext;
    }
    return "";
}

std::string gsFileManager::getBasename(std::string const & fn)
{
    if(fn.find_last_of(".") != std::string::npos)
    {
        size_t pos1 = fn.find_last_of(getValidPathSeparators());
        size_t pos2 = fn.rfind(".");
        std::string name = fn.substr(pos1+1, pos2-pos1-1);
        return name;
    }
    return fn;
}

std::string gsFileManager::getFilename(std::string const & fn)
{
    size_t pos1 = fn.find_last_of(getValidPathSeparators());
    if(pos1 != std::string::npos)
    {
        std::string name = fn.substr(pos1+1);
        return name;
    }
    return fn;
}

namespace {
struct gsStringView {
    const char* m_begin;
    const char* m_end;

    gsStringView( const std::string& s, size_t b, size_t e )
        : m_begin(s.c_str()+b), m_end(s.c_str()+e) {}

    bool operator==(const char* c) const
    {
        for (const char* it=m_begin; it<m_end; ++it, ++c)
            if ( *it != *c ) return false;
        return *c == '\0';
    }

    const char* begin() const { return m_begin; }
    const char* end()   const { return m_end;   }

};
} // end anonymous namespace

std::string gsFileManager::getCanonicRepresentation(const std::string& s)
{
#if defined _WIN32
#else
    std::vector<gsStringView> parts;
    size_t last = 0;
    for (size_t i=0; i<s.size(); ++i)
        if (_contains(getValidPathSeparators(),s[i]))
        {
            parts.push_back(gsStringView(s,last,i));
            last = i + 1;
        }
    parts.push_back(gsStringView(s,last,s.size()));

    std::vector<gsStringView> result;
    size_t sz = parts.size();
    for (size_t i=0; i<sz; ++i)
    {
        if (parts[i] == "" && i > 0 && i < sz-1)
            ; // discard part
        else if (parts[i] == "." && i > 0)
            ; // discard
        else if (parts[i] == "..")
        {
            if (result.empty() || result.back() == "..")
                result.push_back(parts[i]);
            else if (result.back() == "")
                gsWarn << "Cannot go above root.\n";
            else if (result.back() == ".")
            {
                result.pop_back();
                result.push_back(parts[i]);
            }
            else
                result.pop_back();
        }
        else
            result.push_back(parts[i]);
    }

    std::string final_result;
    final_result.append( result[0].begin(), result[0].end() );
    for (size_t i=1; i<result.size(); ++i)
    {
        final_result.push_back( getNativePathSeparator() );
        final_result.append( result[i].begin(), result[i].end() );
    }
    return final_result;
#endif
}


// todo: return bool for success/failure
void gsFileManager::open(const std::string & fn)
{

#if defined __APPLE__
    const int ret = std::system( ("open " + fn + " &").c_str() );
#elif defined __unix__  //__linux__
    const int ret = std::system( ("xdg-open " + fn + " &").c_str() );
#elif defined _WIN32
    HINSTANCE hi = ShellExecute(GetDesktopWindow(), "open", fn.c_str(),
                                NULL, NULL, SW_SHOWNORMAL);
    const bool ret = !( (INT_PTR)hi>32);
#else
    GISMO_STATIC_ASSERT(0,"Platform not identified");
#endif
    //return ret;
    if (0!=ret)
        gsWarn<<"\nFailed to open file "<<fn<<
            " using OS preferred application.\n\n";
}

} //namespace gismo
