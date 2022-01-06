#!/usr/bin/python

""""
    @file fileIO_example.py

    @brief File input/output using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs

