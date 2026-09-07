#!/usr/bin/python

""""
    @file function_trampoline.py

    @brief Illustrate how to overload a gsFunction with a python class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
"""

##################################################################
# This part is needed when pygismo is built inside your build 
# folder using `make pygismo`
import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.insert(0, gismo_path)
##################################################################


import numpy as np
import pygismo as gs


class MyFunction(gs.core.gsFunction):
    def __init__(self):
        super().__init__()

    def domainDim(self):
        return 2

    def targetDim(self):
        return 1

    def eval_into(self, u: np.ndarray, result: np.ndarray) -> None:
        x = u[0, :]
        y = u[1, :]
        result[:] = np.sin(np.pi * x) * np.sin(np.pi * y)
        print(result)

    def deriv_into(self, u: np.ndarray, result: np.ndarray) -> None:
        x = u[0, :]
        y = u[1, :]
        result[0, :] = np.pi * np.cos(np.pi * x) * np.sin(np.pi * y)  # df/dx
        result[1, :] = np.pi * np.sin(np.pi * x) * np.cos(np.pi * y)  # df/dy

    def __call__(self, u: np.ndarray) -> np.ndarray:
        return self.eval(u)


def main():
    # Create a function and a basis
    f = MyFunction()
    # Create evaluation points
    eval_points = np.array([[0.25, 0.25], [0.5, 0.5], [0.75, 0.75]]).T
    # Evaluate the function at the points
    f_values = f.eval(eval_points)
    result = np.empty((1, eval_points.shape[1]))
    f.eval_into(eval_points, result)
    print("Function values at evaluation points:")
    print(f_values)
    print(result)


if __name__ == "__main__":
    main()
