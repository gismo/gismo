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

    /// Get native path seperator
    static char getNativePathSeparator();

    /// Checks if the path is fully qualified
    /// If a name starts with "/", it is considered fully qualified
    static bool isFullyQualified(const std::string& fn);

    /// Checks if the path is a relative path
    /// If a name starts with "./" or "../", it is considered fully qualified
    static bool isRelative(const std::string& fn);

    /// Set the search paths
    static bool setSearchPaths(const std::string& paths);

    /// Add more search paths
    static bool addSearchPaths(const std::string& paths);

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
    /// If \a fn satisfied \a isFullyQualified or \a isRelative, it is kept unchanged.
    ///
    /// In any case, slashes are replaced by the native path separator.
    static std::string find(std::string fn);

    /// Checks if the file exists (also in the search paths)
    static bool fileExists(const std::string& name);

    /// Checks if the directory named \a path exists
    static bool dirExists(const std::string& path);

    /// \brief Find a file in GISMO_DATA_DIR
    ///
    /// \param fn The filename
    /// \returns  The full path or empty string
    ///
    /// If the file can be found, returns the full path.
    /// Otherwiese, returns empty string.
    ///
    /// In any case, slashes are replaced by the native path separator.
    static std::string findInDataDir(std::string fn);

    /// Checks if the file exists in GISMO_DATA_DIR
    static bool fileExistsInDataDir(const std::string& name);

    /// Auto-detect temp directory
    static std::string getTempPath();

    /// Get current directory
    static std::string getCurrentPath();

    /// Get path of executable
    static std::string getExePath();

    /// Make directory
    static bool mkdir( std::string fn );

    /// Checks paths for equality, ignoring slash vs. backslash
    static bool pathEqual( const std::string& p1, const std::string& p2 );

    /// Returns the extension of the filename \a fn
    static std::string getExtension(std::string const & fn);

    /// Returns the base name without extension of the filename \a fn
    static std::string getBasename(std::string const & fn);

    /// Returns the filename without the path of \a fn
    static std::string getFilename(std::string const & fn);

    /// \brief Returns the canonic representation of the path \a fn
    ///
    /// This reduces foo/baz/../bar or foo/./bar to foo/bar
    static std::string getCanonicRepresentation(const std::string & fn);

    /// Opens the file \a fn using the preferred application of the OS
    static void open(const std::string & fn);

private:

    // Return true iff file \a fn exists on the hard disk
    static bool fileNotPathExists(const std::string & fn);

    // The result of argv[0]
    // Since its static, it will be null by default
    // This is called by gsCmdLine and via the unittest runner
    // There is no need to call it otherwise
    static void setArgv0( const char * c );

    friend class gsCmdLine;
    friend class gsUnitTestSelector;
};


} // namespace gismo
