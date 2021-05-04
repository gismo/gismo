#!/usr/bin/python

""""
    @file KnotVector  example

    @brief Play with a Bspline KnotVector in Python

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
"""

import sys
sys.path.append('/Users/graziellacarboni/gismo/build/lib')

import pygismo as gs
import numpy as np

def printKnt(ku):
    for x in range(ku.size()):
        print(ku.knot(x))
    return ku

    
#Empty contructor
ke = gs.nurbs.gsKnotVector()

c = np.array([0.,0.,0.,0.,0.5,1.,1.,1.,1.])
ku = gs.nurbs.gsKnotVector(c,3)
print(f"ku contains {ku.size()} elements.")
print(f"ku contains {ku.uSize()} unique elements.")
printKnt(ku)

ku.insert(0.25,1)
print(f"Aftern insertion of {ku.knot(4)} ku contains {ku.size()} elements.")
print(f"Aftern insertion of {ku.knot(4)} ku contains {ku.uSize()} unique elements.")
printKnt(ku)

print(f"The value of the 2-nd unique knot is {ku.uValue(1)}.")
print(f"Number of knot intervals inside the domain: ", ku.numElements())
print("ku element multiplicities: ", ku.multiplicities())

#print("Find uPointer: ", ku.uFind(0)) to check
#print("Find iPointer: ", ku.iFind(0)) to check

print("ku is consistent: ", ku.check())
#print("Sanity check: ", ku.isConsistent([0.,0.,0.,0.,0.25,0.5,1.,1.,1.,1.], [4,1,1,4])) to check
print("0.2 is inside ku domain: ", ku.inDomain(0.2))
print("1.5 is inside ku domain: ", ku.inDomain(1.5))


