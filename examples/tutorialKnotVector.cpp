/** @file tutorialKnotVector

    @brief Tutorial on gsKnotVector class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <string>
#include <algorithm>
#include <gismo.h>

using namespace gismo;

// forward declaration of some utility functions
void printKnotVector(const gsKnotVector<>& kv, const std::string& name);
void printKnotVector(const gsKnotVector<>& kv);
std::vector<real_t> makeVectorOfKnots();
void print(const real_t& el);


int main(int argc, char* argv[])
{
    // ======================================================================
    // different construction of a knot vector
    // ======================================================================
    
    std::cout << "------------- Constructions -----------------------------\n";

    gsKnotVector<> kv0; 
    printKnotVector(kv0, "kv0");
    
    // only with degree
    gsKnotVector<> kv1(3);
    printKnotVector(kv1, "kv1");
    
    // degree and number of knots
    gsKnotVector<> kv2(3, 8);
    printKnotVector(kv2, "kv2");
    
    real_t a = 1.0; // starting knot
    real_t b = 3.0; // ending knot
    unsigned interior = 4; // number of interior knots
    unsigned multEnd = 3; // multiplicity at the two end knots
    gsKnotVector<> kv3(a, b, interior, multEnd);
    printKnotVector(kv3, "kv3");

    std::vector<real_t> knots = makeVectorOfKnots();
    gsKnotVector<> kv4(knots, 2, 1); // knots, degree, regularity
    printKnotVector(kv4, "kv4");

    gsKnotVector<> kv5(2, knots); // degree, knots
    printKnotVector(kv5, "kv5");

    gsKnotVector<> kv6;
    kv6.initUniform(5, 3); // number of knots, multiple ends
    printKnotVector(kv6, "kv6");
    
    gsKnotVector<> kv7;
    kv7.initClamped(a, b, 3, 5); // start, end, degree, number of interior knots
    printKnotVector(kv7, "kv7");

    
    // ======================================================================
    // looping over knots
    // ======================================================================

    // looping over all knots
    
    std::cout << "\n\n"
              << "------------- Looping over knots -----------------------\n"
              << "kv7: \n";
    for (gsKnotVector<>::iterator it = kv7.begin(); it != kv7.end(); it++)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n\n";
    
    // looping over unique knots
    for (gsKnotVector<>::uiterator it = kv7.ubegin(); it != kv7.uend(); it++)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n\n\n";

    
    // ======================================================================
    // some properties
    // ======================================================================


    std::cout << "------------- Some properties    -----------------------\n"
              << "kv7: \n\n";
    
    printKnotVector(kv7);

    std::cout << "kv7.size(): " << kv7.size() << "\n\n"
              << "kv7.findspan(1.5): " << kv7.findspan(1.5) << "\n\n"
              << "kv7.findspan(2): " << kv7.findspan(2) << "\n\n"
              << "kv7.has(2): " << kv7.has(2) << "\n\n"
              << "kv7.has(2.1): " << kv7.has(2.1) << "\n\n"
              << "kv7.isUniform(): " << kv7.isUniform() << "\n\n"
              << "kv7.numKnotSpans(): " << kv7.numKnotSpans() << "\n\n"
              << "kv7.isOpen(): " << kv7.isOpen() << "\n\n"
              << "kv7.multiplicity(2): " << kv7.multiplicity(2) << "\n\n"
              << "kv7.multiplicity(1): " << kv7.multiplicity(1) << "\n\n\n";

    
    // ======================================================================
    // some operations
    // ======================================================================

    std::cout << "------------- Some operations    -----------------------\n";
    printKnotVector(kv6, "kv6");
    
    
    std::vector<real_t> unique = kv6.unique();
    std::cout << "\nUnique knots: \n";
    std::for_each(unique.begin(), unique.end(), print);
    
    gsMatrix<>* greville = kv6.greville();
    std::cout << "\n\nGreville points: \n" << *greville << "\n\n";
    delete greville;
    
    
    std::vector<int> mult = kv6.multiplicities();
    std::cout << "Multiplicities: ";
    std::for_each(mult.begin(), mult.end(), print);
    std::cout << "\n\n";
    
    printKnotVector(kv6, "kv6");
    
    std::cout << "kv6.uniformRefine()\n";
    kv6.uniformRefine();
    printKnotVector(kv6);
    
    std::cout << "kv6.degreeElevate()\n";
    kv6.degreeElevate();
    printKnotVector(kv6);
    
    std::cout << "For other capabilites of gsKnotVector look at "
        "src/gsNurbs/gsKnotVector.h\n" << std::endl;

    return 0;
}

void print(const real_t& el)
{
    std::cout << el << " ";
}


void printKnotVector(const gsKnotVector<>& kv,
                     const std::string& name)
{
    std::cout << name << ":\n";
    kv.print(std::cout);
    std::cout << "\n" << std::endl;
}


void printKnotVector(const gsKnotVector<>& kv)
{
    for (gsKnotVector<>::const_iterator it = kv.begin(); it != kv.end(); it++)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n\n";
}



std::vector<real_t> makeVectorOfKnots()
{
    std::vector<real_t> knots;
    knots.push_back(0);
    knots.push_back(0.1);
    knots.push_back(0.5);
    knots.push_back(0.6);
    knots.push_back(0.9);
    knots.push_back(1);
    
    return knots;
}
