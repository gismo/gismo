#!/usr/bin/env python3

""""
    @file commandLineArg_example.py

    @brief Tutorial on how to use command line parser in G+Smo.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
"""

import sys
import pygismo as gs

# Variables that will take values from the command line
string='''none'''  # string variable default value
flNumber = 1.0     # flNumber variable default value
boolean = False    # boolean variable default value
plainString=''     # argument of reading plain string

intvec = []        # integer vector

"""
-----------------------------------------------------------------
First we Initialize the object that sets up and parses command line arguments

This defines by default 3 arguments that can be readily used:

 --,  --ignore_rest
Ignores the rest of the labeled arguments following this flag.

--version
Displays version information and exits.

-h,  --help
Displays usage information for all other arguments and exits.
"""

cmd = gs.io.gsCmdLine("Tutorial Command Line Arguments")

"""
-----------------------------------------------------------------
General syntax to add an argument:
cmd.addType("f", "flag", "Description", destination)
"f"    is the short flag: -f
"flag" is the long  flag: --flag (same effect as "-f")
"Description" describes what this argument is about
destination is the variable that will have the value of the input argument
"""

"""
-----------------------------------------------------------------
Adding a string argument, given by the "-s" (or "--stringArg") flag
If set, string is updated to the input value, otherwise string remains untouched
"""
#cmd.addString("s", "stringArg",
#              "Description of string command line argument.",
#              string)

"""
-----------------------------------------------------------------
Adding an int argument, given by the "-i" (or "--num") flag
If set, number is updated to the input value, otherwise number remains untouched
"""
cmd.addNewInt   ("i", "num",
                 "Description of int command line argument",
                 1)

"""
-----------------------------------------------------------------
Adding a float argument, given by the "-r" (or "--real") flag
If set, flNumber is updated to the input value, otherwise flNumber remains untouched
"""
cmd.addReal  ("r", "real",
              "Description of float command line argument",
              flNumber)

"""
-----------------------------------------------------------------
Adding a switch argument, given by the "--bool" flag
If set, boolean is updated to the input value, otherwise boolean remains untouched
"""
cmd.addSwitch("bool","Description of the switch argument.", boolean);

"""
-----------------------------------------------------------------
Extra plain argument (manually defined):
Plain arguments are given without a flag.
"""
#cmd.addPlainString("plain", "Description of the plain command line argument.", plainString)

"""
-----------------------------------------------------------------
Each flag can be only called once. The commands
cmd.addMultiString, cmd.addMultiInt and cmd.addMultiReal
allow one to register flags that can be used several times.
They store the data in a vector.
"""
#cmd.addMultiInt("m", "multiint", "Description of multiint command line argument.", intvec);

"""
-----------------------------------------------------------------
Read the arguments and update with the inputs, if given.
The program stopps here if there was an error in the command line
or the user invoked "--help" or "--version"
"""

try:
    cmd.getValues(sys.argv)
except:
    print("An error occured")

number= cmd.getInt("num")

print('''
Printing command line arguments:\n\n\n
Plain string:  ''' + plainString   + ''' \n\n
String:        ''' + string        + ''' \n\n
Float:         ''' + str(flNumber) + ''' \n\n
Integer:       ''' + str(number)   + ''' \n\n
Switch:        ''' + str(boolean)  + ''' \n\n
MultiInt:      ''' + str(intvec)   + ''' \n\n
''')
