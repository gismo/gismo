/** @file gsJSONValidator.h
  
    @brief Header file for JSON validation for preCICE

    https://spectralib.org/doc

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li (2023,..., TU Delft)
*/

#pragma once

#include <gsJSON/gsJSON.h>
#include <string>

namespace gismo 
{

class gsJSONValidator 
{
public:
    // Constructor with schema file
    gsJSONValidator(const std::string& schemaFile) 
    {
        // Load schema from file
        m_schema = gsJSON(schemaFile);
    }

    // Basic validation method
    // This is a simplified implementation that doesn't do full schema validation
    // but checks basic structure requirements instead
    bool validate(const gsJSON& document, std::string& error) 
    {
        // Basic structure validation
        try {
            // Check for required participant section
            if (!document.contains("participant")) {
                error = "Missing 'participant' section in configuration";
                return false;
            }
            
            // Check for required participant properties
            if (!document["participant"].contains("name")) {
                error = "Missing 'name' in participant section";
                return false;
            }
            
            if (!document["participant"].contains("preciceConfigFile")) {
                error = "Missing 'preciceConfigFile' in participant section";
                return false;
            }
            
            // Check for required meshes section
            if (!document.contains("meshes") || !document["meshes"].get<json>().is_array()) {
                error = "Missing or invalid 'meshes' array in configuration";
                return false;
            }
            
            // Validate each mesh
            for (const auto& mesh : document["meshes"]) 
            {
                if (!mesh.contains("name") || !mesh.contains("type")) 
                {
                    error = "Mesh missing required 'name' or 'type' property";
                    return false;
                }
                
                std::string type = mesh["type"].get<std::string>();
                if (type == "structured") 
                {
                    if (!mesh.contains("nx") || !mesh.contains("ny") || 
                        !mesh.contains("dim") || !mesh.contains("x-range") || 
                        !mesh.contains("y-range")) {
                        error = "Structured mesh missing required properties";
                        return false;
                    }
                } else if (type == "custom") 
                {
                    if (!mesh.contains("points")) 
                    {
                        error = "Custom mesh missing 'points' property";
                        return false;
                    }
                } else
                {
                    error = "Unknown mesh type: " + type;
                    return false;
                }
            }
            
            return true;
        } 
        catch (const std::exception& e) 
        {
            error = std::string("Validation error: ") + e.what();
            return false;
        }
    }

private:
    gsJSON m_schema;  // Not used in this simplified implementation
};

} // namespace gismo