/** @file gsOptimizerRegistry.h

    @brief Registry for optional optimizers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsOptimizer/gsOptimizer.h>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include <stdexcept>

namespace gismo
{

/**
 * @brief      Registry for runtime optimizer discovery
 *
 * This class provides a runtime registry for optimizer solvers.
 * Optional modules (like gsOptim, gsHLBFGS) register themselves
 * on load, and users can retrieve them by name.
 *
 * Usage:
 *   if (gsOptimizerRegistry::has("gsOptim-BFGS"))
 *       auto opt = gsOptimizerRegistry::get("gsOptim-BFGS");
 */
class gsOptimizerRegistry
{
public:
    /**
     * @brief      Register an optimizer factory
     *
     * @param      name     The name to register (e.g., "gsOptim-BFGS")
     * @param      factory  Factory function that creates the optimizer
     *
     * @note       If the name is already registered, this is a no-op (idempotent)
     */
    static void register_(const std::string& name,
                          std::function<gsOptimizer<real_t>::uPtr()> factory)
    {
        auto& instance = getMutableInstance();
        
        auto it = instance.factories.find(name);
        if (it != instance.factories.end())
        {
            return;
        }
        instance.factories[name] = std::move(factory);
    }

    /**
     * @brief      Get an optimizer by name
     *
     * @param      name     The optimizer name (e.g., "gsOptim-BFGS")
     *
     * @return     Unique pointer to the optimizer
     *
     * @throws     std::runtime_error if the optimizer is not found
     */
    static gsOptimizer<real_t>::uPtr get(const std::string& name)
    {
        auto& instance = getMutableInstance();
        auto it = instance.factories.find(name);
        if (it == instance.factories.end())
        {
            throw std::runtime_error("OptimizerRegistry: '" + name + "' not found. "
                                      "Available optimizers: " + listAsString());
        }
        return it->second();
    }

    /**
     * @brief      Check if an optimizer is available
     *
     * @param      name     The optimizer name
     *
     * @return     true if the optimizer is registered
     */
    static bool has(const std::string& name)
    {
        auto& instance = getMutableInstance();
        return instance.factories.find(name) != instance.factories.end();
    }

    /**
     * @brief      Get list of all available optimizers
     *
     * @return     Vector of optimizer names
     */
    static std::vector<std::string> list()
    {
        auto& instance = getMutableInstance();
        std::vector<std::string> names;
        names.reserve(instance.factories.size());
        for (const auto& kv : instance.factories)
            names.push_back(kv.first);
        return names;
    }

private:
    struct RegistryData
    {
        std::map<std::string, std::function<gsOptimizer<real_t>::uPtr()>> factories;
    };

    // Definition lives in gsOptimizerRegistry.cpp (compiled into libgismo.so).
    // GISMO_EXPORT ensures that MODULE .so files loaded via dlopen resolve this
    // symbol from the shared libgismo.so, giving one authoritative singleton.
    static GISMO_EXPORT RegistryData& getMutableInstance();

    static std::string listAsString()
    {
        auto& instance = getMutableInstance();
        if (instance.factories.empty())
            return "(none)";
        std::string result;
        for (const auto& kv : instance.factories)
        {
            if (!result.empty())
                result += ", ";
            result += kv.first;
        }
        return result;
    }
}; // end class

} // end namespace gismo
