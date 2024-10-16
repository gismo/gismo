#!/usr/bin/python

""""
@file BSpline curve example

@brief Play with a B-spline curve in Python

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http:#mozilla.org/MPL/2.0/.

Author(s): S. Imperatore
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np

#################### STEP 1: read a multi-patch file #################
# Variables initialization with more descriptive names
paramMethod = 0 # Parameterization method
verboseMode = 0 # Verbose mode: 0 - no print, 1 - summary, 2 - full print
numRefine = 1 # Number of Uniform h-refinement
numElevate = 0 # Number of degree elevation steps to perform
AAPreconditionType = 0 # Preconditioner type for Anderson Acceleration
isInterfaceFree = False # Make interfaces free
isBRep = False # Input data is B-Rep

# Constants for input file and description
INPUT_FILE = "breps/other/TUDflame.xml"
DESCRIPTION = "Hi, give me a file (eg: .xml) containing multi-patch and I will try to parameterize it!"

# Load XML file - multi-patch computational domain
# ! [Read geometry]
# Check if the input file exists

# MultiPatch reader
mp = gs.core.gsMultiPatch();
print("Read file \"" + INPUT_FILE + "\"\n")
fd = gs.io.gsFileData(INPUT_FILE)
fd.getAnyFirst(mp);
# ! [Read geometry]

#! [construct a initial parameterization if the input data is B-Rep]
pardim = mp.domainDim()
geodim = mp.targetDim()
# If input data is B-Rep or dimension does not match, construct an initial parameterization
if (isBRep | pardim != geodim):
    springInitializer = gs.modelling.gsSpringPatch(mp)
    mp.clear()
    springInitializer.compute()
    mp.addPatch(springInitializer.result())
#! [construct a initial parameterization if the input data is B-Rep]

#! [align the orientations of a multi-patch parameterization]
mp.fixOrientation()
#! [align the orientations of a multi-patch parameterization]

#! [refine the geometry for better result]
# Elevate and p-refine the multi-patch geometry to order p + numElevate
if (numElevate != 0):
    mp.degreeElevate(numElevate)
    print(str(numElevate) + " degree elevation(s) performed!\n")


# h-refine each patch
if (numRefine != 0):
    for r in range(0,numRefine):
        mp.uniformRefine()
    print(str(numRefine) + " uniform h-refinement(s) performed!\n")

#! [perform analysis-suitable parameterization construction]
result = gs.core.gsMultiPatch();
if (geodim == 2):
    opt = gs.modelling.gsBarrierPatch2(mp,False)
    opt.options().setInt("Verbose", verboseMode)
    opt.options().setInt("ParamMethod", paramMethod)
    opt.options().setInt("AAPreconditionType", AAPreconditionType)
    opt.compute()
    result = opt.result()
elif (geodim == 3):
    opt = gs.modelling.gsBarrierPatch3(mp,not isInterfaceFree)
    opt.options().setInt("Verbose", verboseMode)
    opt.options().setInt("ParamMethod", paramMethod)
    opt.compute()
    result = opt.result()
else:
    gsInfo << "current version only support pardim = geodim = 2 or 3.\n"
    quit(1)
#! [perform analysis-suitable parameterization construction]

#! [output the parameterization result]
# outputFilename = INPUT_FILE
# outputFilename = gsFileManager::getBasename(outputFilename) # file name without a path
# outputFilename += "_result"
# outputFilename += "R"+std::to_string(numRefine)+"E"+std::to_string(numElevate)
# outputFilename += (!isInterfaceFree) ? "Fixed" : "Free"

# outputResult(result, outputFilename)
# gsInfo << "The results have been written into " << outputFilename << "\n"
#! [output the parameterization result]
