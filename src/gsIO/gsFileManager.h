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

/// @brief File-system related functionality.
///
/// Input paths to the fuctions of this class can be given using any
/// valid path seperator, this is in Unix only "/" and in Windows both
/// "/" and "\\".
///
/// Return values only contain the preferred native path seperator, which
/// is in Unix "/" and in Windows "\\".
///
/// @ingroup IO
class GISMO_EXPORT gsFileManager
{
public:

    /// Get preferred native path seperator
    static char getNativePathSeparator();

    /// Get valid path seperators
    static const std::string& getValidPathSeparators();
    
    /// @brief Checks if the path is fully qualified, also known as "absolute path"
    /// Under Unix, if a name starts with "/", it is considered fully qualified.
    /// Under Windows, it starts with the drive-letter followed by the path or
    /// with a "/" or a "\\" (refers to the current drive).
    static bool isFullyQualified(const std::string& fn);

    /// @brief Checks if the path is a relative path
    /// Under Unix, if a name starts with "./" or "../", it is considered relative.
    /// Under Windows, if a name starts with "./", ".\\", "../" or "..\\",
    /// it is considered relative.
    static bool isExplicitlyRelative(const std::string& fn);

    /// @brief Set the search paths
    /// Returns true iff all paths exist
    static bool setSearchPaths(const std::string& paths);

    /// @brief Add more search paths
    /// Returns true iff paths exist
    static bool addSearchPaths(const std::string& paths);

    /// @brief Get the defined search paths
    static std::string getSearchPaths();

    /// @brief Find a file.
    ///
    /// @param fn The filename
    /// @returns The full path or empty string
    ///
    /// If the fn \a isFullyQualified (like "/foo/bar.txt"), or
    /// if the fn \a isExplicitlyRelative (like "../foo/bar.txt"), the
    /// name is returned unchanged* if the file can be found. Otherwise,
    /// an empty string is returned.
    ///
    /// If the fn has the form "bar.txt" or "foo/bar.txt", the file
    /// is searched in the current directory and all search paths
    /// (cf. \a getSearchPaths). If the file can be found, the full
    /// path is returned. Otherwise, an empty string is returned.
    ///
    /// *: In any case, slashes are replaced by the native path separator.
    static std::string find(std::string fn);

    /// @brief Checks if the file exists
    ///
    /// If the fn \a isFullyQulaified (like "/foo/bar.txt"), or
    /// if the fn \a isExplicitlyRelative (like "../foo/bar.txt"), the
    /// only this path is considered.
    ///
    /// If the fn has the form "bar.txt" or "foo/bar.txt", the file
    /// is searched in the current directory and all search paths
    /// (cf. \a getSearchPaths).
    ///
    /// @see \a find
    static bool fileExists(const std::string& name);

    /// @brief Find a file in GISMO_DATA_DIR
    ///
    /// @param fn The filename
    /// @returns  The full path or empty string
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

    /// Get path of executable (without filename)
    static std::string getExePath();

    /// @brief Make directory
    ///
    /// Return true iff directory is available after calling
    /// this function. (This also holds if the directory has
    /// existed already.)
    static bool mkdir( std::string fn );

    /// @brief Checks paths for equality of paths
    ///
    /// If the path is not \a isFullyQualified, creates an absolute
    /// path using \a getCurrentPath. Then \a getCanonicRepresentation
    /// is called to get canonical representations which are compared.
    static bool pathEqual( const std::string& p1, const std::string& p2 );

    /// Returns the extension of the filename \a fn
    static std::string getExtension(std::string const & fn);

    /// Returns the base name without extension of the filename \a fn
    static std::string getBasename(std::string const & fn);

    /// Returns the filename without the path of \a fn
    static std::string getFilename(std::string const & fn);

    /// @brief Returns the canonic representation of the path \a fn
    ///
    /// This reduces foo/baz/../bar or foo/./bar to foo/bar. Moreover,
    /// the non-preferred path seperators are replaced by the preferred
    /// ones.
    ///
    /// This does not access the file system, the current directory, or
    /// the search paths.
    static std::string getCanonicRepresentation(const std::string & fn);

    /// Opens the file \a fn using the preferred application of the OS
    static void open(const std::string & fn);

private:
    // The result of argv[0]
    // Since its static, it will be null by default
    // This is called by gsCmdLine and via the unittest runner
    // There is no need to call it otherwise
    static void setArgv0( const char * c );

    friend class gsCmdLine;
    friend class gsUnitTestSelector;
};


} // namespace gismo
