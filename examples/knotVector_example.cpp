/** @file knotVector_example.cpp

    @brief Tutorial on gsKnotVector class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh, S. Imperatore
*/

//![Include namespace]
#include <iostream>
#include <string>
#include <algorithm>
#include <gismo.h>

using namespace gismo;
//![Include namespace]

// forward declaration of some utility functions
void printKnotVector(const gsKnotVector<>& kv, const std::string& name);
void printKnotVector(const gsKnotVector<>& kv);
void print(const real_t& el);


int main(int argc, char* argv[])
{
    gsCmdLine cmd("Tutorial on gsKnotVector class.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    // ======================================================================
    // different construction of a knot vector
    // ======================================================================

    gsInfo << "------------- Constructions -------------\n";

    //! [empty initialization]
    // --- empty constructor
    gsKnotVector<> kv0;

    // --- only with degree
    gsKnotVector<> kv1(3);
    //! [empty initialization]
    printKnotVector(kv0, "kv0");
    printKnotVector(kv1, "kv1");

    //! [uniform initialization]
    // --- uniform initialization
    real_t a = 1; // starting knot
    real_t b = 3; // ending knot
    index_t interior = 4; // number of interior knots
    index_t multEnd = 3; // multiplicity at the two end knots

    gsKnotVector<> kv2(a, b, interior, multEnd);

    index_t multInt = 2; // multiplicity of the interior knots (default = 1)
    gsKnotVector<> kv3(a, b, interior, multEnd, multInt);
    
    index_t degree = 5; // degree of the spline space (default = -1);
    gsKnotVector<> kv4(a, b, interior, multEnd, multInt, degree);
    //! [uniform initialization]
    printKnotVector(kv2, "kv2");
    printKnotVector(kv3, "kv3");
    printKnotVector(kv4, "kv4");
    

    //! [knotContainer initialization]
    // --- construction from knot container
    std::vector<real_t> knotContainer;
    knotContainer.push_back(0);
    knotContainer.push_back(0.1);
    knotContainer.push_back(0.5);
    knotContainer.push_back(0.6);
    knotContainer.push_back(0.9);
    knotContainer.push_back(1);

    gsKnotVector<> kv5(knotContainer, 2); // knots, degree
    //! [knotContainer initialization]
    printKnotVector(kv5, "kv5");

    // --- more initializations
    //! [more uniform initialization]
    gsKnotVector<> kv6;
    kv6.initUniform(5, 3); // number of knots, multiple ends
    printKnotVector(kv6, "kv6");
    //! [more uniform initialization]

    //! [clamped initialization]
    gsKnotVector<> kv7;
    kv7.initClamped(a, b, 3, 5); // start, end, degree, number of interior knots
    printKnotVector(kv7, "kv7");
    //! [clamped initialization]


    // ======================================================================
    // some properties
    // ======================================================================

    gsInfo << "------------- Some properties -------------\n"
              << "kv7: \n\n";
    
    printKnotVector(kv7, "kv7");

    //! [kv properties]
    gsInfo  << "kv7.size(): "            << kv7.size()                   << "\n"
            << "kv7.findspan(1.5): "     << kv7.iFind(1.5) - kv7.begin() << "\n"
            << "kv7.findspan(2): "       << kv7.iFind(2) - kv7.begin()   << "\n"
            << "kv7.has(2): "            << kv7.has(2)                   << "\n"
            << "kv7.has(2.1): "          << kv7.has(2.1)                 << "\n"
            << "kv7.isUniform(): "       << kv7.isUniform()              << "\n"
            << "kv7.isOpen(): "          << kv7.isOpen()                 << "\n"
            << "kv7.multiplicity(4/3): " << kv7.multiplicity(4./3)       << "\n"
            << "kv7.numKnotSpans(): "    << kv7.uSize() - 1              << "\n\n";
    //! [kv properties]

    // ======================================================================
    // some operations
    // ======================================================================

    gsInfo << "------------- Some operations -------------\n";
    printKnotVector(kv6, "kv6");

    //! [kv operations]
    std::vector<real_t> unique = kv6.unique();
    gsInfo << "\nUnique knots: \n";
    std::for_each(unique.begin(), unique.end(), print);

    gsMatrix<> greville = kv6.greville();
    gsInfo << "\n\nGreville points: \n" << greville << "\n\n";

    std::vector<index_t> mult = kv6.multiplicities();
    gsInfo << "Multiplicities: ";
    std::for_each(mult.begin(), mult.end(), print);
    gsInfo << "\n\n";

    printKnotVector(kv6, "kv6");

    gsInfo << "kv6.uniformRefine()\n";
    kv6.uniformRefine();
    printKnotVector(kv6);

    gsInfo << "kv6.degreeElevate()\n";
    kv6.degreeElevate();
    printKnotVector(kv6);
    //! [kv operations]

    // ======================================================================
    // looping over knots
    // ======================================================================

    gsInfo << "\n"
              << "------------- Looping over knots -------------\n"
              << "kv4: \n";
    //![kv loop]
    for (gsKnotVector<>::iterator it = kv4.begin(); it != kv4.end(); it++)
    {
        gsInfo << *it << " ";
    }
    gsInfo << "\n\n";

    // looping over unique knots
    for (gsKnotVector<>::uiterator it = kv4.ubegin(); it != kv4.uend(); it++)
    {
        gsInfo << *it << " ";
    }
    gsInfo << "\n\n";
    //![kv loop]


    gsInfo << "For other capabilites of gsKnotVector look at "
        "src/gsNurbs/gsKnotVector.h\n" << "\n";

    return 0;
}

void print(const real_t& el)
{
    gsInfo << el << " ";
}


void printKnotVector(const gsKnotVector<>& kv,
                     const std::string& name)
{
    gsInfo << name << ":\n" << kv << "\n";
    gsInfo << "knot values:\n";
    for (gsKnotVector<>::const_iterator it = kv.begin(); it != kv.end(); it++)
    {
        gsInfo << *it << " ";
    }
    gsInfo << "\n\n";
}


void printKnotVector(const gsKnotVector<>& kv)
{
    for (gsKnotVector<>::const_iterator it = kv.begin(); it != kv.end(); it++)
    {
        gsInfo << *it << " ";
    }
    gsInfo << "\n\n";
}

