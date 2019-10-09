/** @file gsFileManager.h

    @brief Utility class for finding files and handling paths

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, A. Mantzaflaris, J. Vogl
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
/// Return values only contain the preferred native path separator, which
/// is in Unix "/" and in Windows "\\".
///
/// Output Paths always end with path separator, files without.
///
/// Does not check for special system cases, like AUX under Windows.
/// (see https://en.wikipedia.org/wiki/Filename)
/// @ingroup IO
class GISMO_EXPORT gsFileManager
{
public:

    /// Get preferred native path seperator
    static char getNativePathSeparator();

    /// Get valid path seperators
    static const std::string& getValidPathSeparators();

    /// Get system-dependent invalid characters for paths and filenames
    static const std::string& getInvalidCharacters();
    
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
    /// if the fn \a isExplicitlyRelative (like "../foo/bar.txt"), then
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
    /// Otherwise, returns empty string.
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

    /// Get path to home directory.
    static std::string getHomePath();

    /// Calculates the relative path from absolute path \a from to path \a to.
    /// Example:
    /// from  : /a/b/cd/
    /// to    : /a/b/c/d
    /// result: ../c/d
    /// Return an empty string, if \a from isn't absolute.
    /// Expands \a to with getCurrentPath if not absolute.
    static std::string makeRelative(const std::string& from, const std::string& to);

    /// @brief Make directory
    ///
    /// Return true iff directory is available after calling
    /// this function. (This also holds if the directory has
    /// existed already.)
    /// If a relative path is given as \a fn, it will be expanded
    /// with getCurrentDirectory.
    static bool mkdir( std::string fn );

    /// @brief Checks paths for equality of paths
    ///
    /// If the path is not \a isFullyQualified, creates an absolute
    /// path using \a getCurrentPath. Then \a getCanonicRepresentation
    /// is called to get canonical representations which are compared.
    static bool pathEqual( const std::string& p1, const std::string& p2 );

    /// Returns the extension of the filename \a fn
    static std::string getExtension(std::string const & fn);

    /// Returns the base name without path and extension of the filename \a fn
    static std::string getBasename(std::string const & fn);

    /// Returns the filename without the path of \a fn
    static std::string getFilename(std::string const & fn);

    /// Returns the path without filename \a fn.
    /// If \a resolve is set, it will resolve relative paths to absolute ones
    /// with use of \a find.
    static std::string getPath(std::string const & fn, bool resolve = false);

    /// @brief Returns the canonic representation of the path \a fn
    ///
    /// This reduces foo/baz/../bar or foo/./bar to foo/bar. Moreover,
    /// the non-preferred path separators are replaced by the preferred
    /// ones.
    ///
    /// This does not access the file system, the current directory, or
    /// the search paths. Therefore, leading .. can't be replaced.
    /// Leading ./././ will be reduced to ./
    ///
    /// \param fn input path as std::string
    /// \param asPath If true, make sure the last character of the output is the native-path-separator.
    /// Else, let is as it was.
    /// \return canonical representation of \a fn as std::string
    static std::string getCanonicRepresentation(const std::string & fn, bool asPath = false);

    /// Opens the file \a fn using the preferred application of the OS
    static void open(const std::string & fn);

};


} // namespace gismo
