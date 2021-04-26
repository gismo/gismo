#!/usr/bin/env python3

""""
    @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
"""

import sys
import pygismo as gs
import numpy as np

c = np.array([[0.,1.],[1.,1.]])
b = gs.nurbs.gsBSpline(0.0,1.0,0,1,c,0,False)
b.degree(0)
b.sample(10)
