/** @file knotVectorVar_example.cpp

    @brief Demonstrates using gsKnotVector with autodiff::var type for reverse-mode AD.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsNurbs/gsKnotVector.hpp>  // Include full implementation for header-only use
#include <autodiff/reverse/var.hpp>

using namespace gismo;
using autodiff::var;  // Use reverse-mode var type
// Note: We explicitly qualify val(), wrt(), derivatives() with full namespace to avoid ambiguity with forward-mode
// Note: This example uses header-only instantiation of gsKnotVector<var> since var cannot be 
// precompiled into the library due to Eigen compatibility issues with ternary operators.

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Demonstrates gsKnotVector with autodiff::var for reverse-mode AD.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "=== gsKnotVector<var> Example ===\n\n";

    // Example 1: Create a simple knot vector with var type
    gsInfo << "Example 1: Basic construction\n";
    {
        std::vector<var> knots = {var(0.0), var(0.0), var(0.5), var(1.0), var(1.0)};
        gsKnotVector<var> kv(knots, 1);
        
        gsInfo << "Created gsKnotVector<var> with " << kv.size() << " knots\n";
        gsInfo << "Degree: " << kv.degree() << "\n";
        gsInfo << "Knots: ";
        for (size_t i = 0; i < kv.size(); ++i) {
            gsInfo << val(kv[i]) << " ";
        }
        gsInfo << "\n\n";
    }

    // Example 2: Uniform knot vector with var
    gsInfo << "Example 2: Uniform knot vector\n";
    {
        var a(0.0), b(1.0);
        gsKnotVector<var> kv(a, b, 3, 2);  // [0,0,0.25,0.5,0.75,1,1]
        
        gsInfo << "Uniform knot vector: ";
        for (size_t i = 0; i < kv.size(); ++i) {
            gsInfo << val(kv[i]) << " ";
        }
        gsInfo << "\n";
        gsInfo << "Number of unique knots: " << kv.uSize() << "\n\n";
    }

    // Example 3: Knot insertion with var
    gsInfo << "Example 3: Knot insertion\n";
    {
        std::vector<var> initial = {var(0.0), var(0.0), var(1.0), var(1.0)};
        gsKnotVector<var> kv(initial, 1);
        
        gsInfo << "Initial knots: ";
        for (size_t i = 0; i < kv.size(); ++i) {
            gsInfo << val(kv[i]) << " ";
        }
        gsInfo << "\n";
        
        // Insert a new knot
        kv.insert(var(0.5));
        
        gsInfo << "After inserting 0.5: ";
        for (size_t i = 0; i < kv.size(); ++i) {
            gsInfo << val(kv[i]) << " ";
        }
        gsInfo << "\n\n";
    }

    // Example 4: Using var for derivative computation
    gsInfo << "Example 4: Derivatives with respect to knot positions\n";
    {
        // Create a knot vector where knot positions are variables
        var k1 = 0.0;  // Fixed start
        var k2 = 0.3;  // Variable knot
        var k3 = 0.7;  // Variable knot
        var k4 = 1.0;  // Fixed end
        
        std::vector<var> knots = {k1, k1, k2, k3, k4, k4};
        gsKnotVector<var> kv(knots, 1);
        
        // Compute some function of the knot spacing
        var spacing1 = k2 - k1;
        var spacing2 = k3 - k2;
        var spacing3 = k4 - k3;
        
        // Objective: minimize variance of spacings
        var mean = (spacing1 + spacing2 + spacing3) / 3.0;
        var variance = ((spacing1 - mean)*(spacing1 - mean) + 
                        (spacing2 - mean)*(spacing2 - mean) + 
                        (spacing3 - mean)*(spacing3 - mean)) / 3.0;
        
        gsInfo << "Knot positions: k2=" << autodiff::reverse::detail::val(k2) 
               << ", k3=" << autodiff::reverse::detail::val(k3) << "\n";
        gsInfo << "Spacings: [" << autodiff::reverse::detail::val(spacing1) 
               << ", " << autodiff::reverse::detail::val(spacing2) 
               << ", " << autodiff::reverse::detail::val(spacing3) << "]\n";
        gsInfo << "Variance of spacings: " << autodiff::reverse::detail::val(variance) << "\n";
        
        // Compute gradients
        auto [dvar_dk2, dvar_dk3] = autodiff::reverse::detail::derivatives(
            variance, autodiff::reverse::detail::wrt(k2, k3));
        
        gsInfo << "Gradient of variance w.r.t. k2: " << dvar_dk2 << "\n";
        gsInfo << "Gradient of variance w.r.t. k3: " << dvar_dk3 << "\n";
        gsInfo << "\n";
        gsInfo << "These gradients can be used for knot optimization!\n\n";
    }

    // Example 5: Demonstrating why var works with gsKnotVector
    gsInfo << "Example 5: Why autodiff::var works with gsKnotVector\n";
    {
        gsInfo << "autodiff::var provides:\n";
        gsInfo << "  - Copy constructor and assignment\n";
        gsInfo << "  - Comparison operators (<, >, ==, <=, >=)\n";
        gsInfo << "  - Standard arithmetic operators\n";
        gsInfo << "  - Compatible with std::vector and std::merge\n";
        gsInfo << "\nAll requirements for gsKnotVector<T> are satisfied!\n";
    }

    return 0;
}
