/** @file gsFileManager.h

    @brief Utility class for finding files and handling paths

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/// @brief This class checks if the given filename can be found
///  in one of the pre-defined search paths. It is possible to
///  register additional search paths.
///
/// @ingroup IO
class GISMO_EXPORT gsFileManager
{
public:

    /// Checks if the file exists (also in the search paths)
    static bool fileExists(const std::string& name);

    /// Checks if the file exists in GISMO_DATA_DIR
    static bool fileExistsInDataDir(const std::string& name);

    /// Get local path seperator
    static char getLocalPathSeperator();

    /// Checks if the path is fully qualified
    /// If a name starts with "/", it is considered fully qualified
    static bool isFullyQualified(const std::string& fn);

    /// Checks if the path is a relative path
    /// If a name starts with "./" or "../", it is considered fully qualified
    static bool isRelative(const std::string& fn);

    /// Set the search paths
    static void setSearchPaths(const std::string& paths);

    /// Get the defined search path
    static std::string getSearchPaths();

    /// \brief Find a file.
    ///
    /// \param fn The filename
    /// \returns  The full path or empty string
    ///
    /// If the file can be found, returns the full path.
    /// Otherwiese, returns empty string.
    ///
    /// If \a fn satisfied \a isFullyQualified or \a isRelative, it is kept unchanged
    static std::string find(std::string fn);

    /// \brief Find a file in GISMO_DATA_DIR
    ///
    /// \param fn The filename
    /// \returns  The full path or empty string
    ///
    /// If the file can be found, returns the full path.
    /// Otherwiese, returns empty string.
    static std::string findInDataDir(std::string fn);

    /// Make directory
    static bool mkdir( std::string fn );

    /// Checks paths for equality, ignoring slash vs. backslash
    static bool pathEqual( const std::string& p1, const std::string& p2 );

    /// Returns the extension of the filename \a fn
    static std::string getExtension(std::string const & fn)
    {
        if(fn.find_last_of(".") != std::string::npos)
        {
            std::string ext = fn.substr(fn.rfind(".")+1);
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
            return ext;
        }
        return "";
    }

    /// Returns the base name without extension of the filename \a fn
    static std::string getBasename(std::string const & fn)
    {
        if(fn.find_last_of(".") != std::string::npos)
        {
            std::size_t pos1 = fn.find_last_of("/\\");
            std::size_t pos2 = fn.rfind(".");
            std::string name = fn.substr(pos1+1, pos2-pos1-1);
            return name;
        }
        return fn;
    }

    /// Returns the filename without the path of \a fn
    static std::string getFilename(std::string const & fn)
    {
        std::size_t pos1 = fn.find_last_of("/\\");
        if(pos1 != std::string::npos)
        {
            std::string name = fn.substr(pos1+1);
            return name;
        }
        return fn;
    }
};


} // namespace gismo
