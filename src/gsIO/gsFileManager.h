/** @file gsFileManager.h

    @brief Utility class for finding files and handling paths

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <string>

namespace gismo 
{


/// @brief Returns the instance of \a gsFileManager
///
/// @ingroup IO
GISMO_EXPORT gsFileManager& gsFileManagerSingleton();

/// @brief This class checks if the given filename can be found
///  in one of the pre-defined search paths. It is possible to
///  register additional search paths.
///
/// This class is a singleton, which can only be instanciated with
/// \a gsFileManagerSingleton
///
/// @ingroup IO
class GISMO_EXPORT gsFileManager
{
public:

    /// Checks if the file exists
    static bool fileExists(const std::string& name);

    /// Checks if the path is fully qualified
    /// If a name starts with "/", it is considered fully qualfied
    static bool isFullyQualified(const std::string& fn) { return fn[0] == GISMO_PATH_SEPERATOR; }

    /// Checks if the path is a relative path
    /// If a name starts with "./" or "../", it is considered fully qualfied
    static bool isRelative(const std::string& fn)
    {
        return ( fn[0] == '.' && fn[1] == GISMO_PATH_SEPERATOR )
            || ( fn[0] == '.' && fn[1] == '.' && fn[2] == GISMO_PATH_SEPERATOR );
    }

    /// Set the search paths
    void setSearchPaths(const std::string& paths);

    /// Get the defined search path
    std::string getSearchPaths() const;

    /// \brief Find a file.
    ///
    /// \param fn[in|out]  The filename
    ///
    /// If the file can be found, returns true and replaces \a fn by the full path.
    /// Otherwiese, returns false and keeps the name unchanged.
    ///
    /// If \a fn satisfied \a isFullyQualified or \a isRelative, it is kept unchanged
    bool find( std::string& fn );

private:
    GISMO_EXPORT friend gsFileManager& gsFileManagerSingleton();
    gsFileManager();
    gsFileManager(const gsFileManager&);
    gsFileManager& operator= (const gsFileManager&);
    std::vector<std::string> m_paths;
};

} // namespace gismo
