"""
    @file Fitting example

    @brief Demonstrates adaptive fitting of data samples in Python. Needs the gismo cppyy bindings.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
"""

import argparse
import time
from gismo_cppyy import gismo
import gismo_cppyy


def main():
    parser = argparse.ArgumentParser(
        prog="fitting_example",
        description="Adaptive fitting of data samples"
    )

    parser.add_argument("--save",
                        dest="save",
                        action="store_const",
                        help="Save result in XML format",
                        const=True, default=False)
    parser.add_argument("-c", "--pacor",
                        default=1,
                        help="Steps of parameter correction",
                        type=int,
                        dest="maxPcIter")
    parser.add_argument("-i", "--iter",
                        default=2,
                        help="number of iterations",
                        type=int,
                        dest="iter")
    parser.add_argument("-x", "--deg_x",
                        default=2,
                        help="degree in x direction",
                        type=int,
                        dest="deg_x")
    parser.add_argument("-y", "--deg_y",
                        default=2,
                        help="degree in y direction",
                        type=int,
                        dest="deg_y")
    parser.add_argument("-s", "--lambda",
                        default=1e-7,
                        help="smoothing coefficient",
                        type=float,
                        dest="smoothingcoef")
    parser.add_argument("-t", "--threshold",
                        default=1e-2,
                        help="error threshold (special valule -1)",
                        type=float,
                        dest="threshold")
    parser.add_argument("-p", "--refPercent",
                        default=0.1,
                        help="percentage of points to refine in each iteration",
                        type=float,
                        dest="refPercent")
    parser.add_argument("-q", "--extension",
                        default=2,
                        help="extension size",
                        type=int,
                        dest="extension")
    parser.add_argument("-r", "--urefine",
                        default=3,
                        help="initial uniform refinement steps",
                        type=int,
                        dest="numURef")
    parser.add_argument("-e", "--tolerance",
                        default=1e-2,
                        help="error tolerance (desired upper bound for pointwise error)",
                        type=float,
                        dest="tolerance")
    parser.add_argument("-d", "--data",
                        default="fitting/deepdrawingC.xml",
                        help="Input sample data",
                        type=str,
                        dest="fn")
    args = parser.parse_args()

    if args.deg_x < 1:
        raise ValueError("Degree x must be positive.")
    if args.deg_y < 1:
        raise ValueError("Degree y must be positive.")
    if args.extension < 0:
        raise ValueError("Extension must be non negative.")
    if args.tolerance < 0:
        raise ValueError("Error tolerance must be non negative.")
    if args.threshold > 0 and args.threshold > args.tolerance:
        raise ValueError("Refinement threshold is over tolerance.")

    # Read data
    # Surface fitting
    # Expected input is a file with matrices with:
    # id 0:  u,v   -- parametric coordinates, size 2 x N
    # id 1:  x,y,z -- corresponding mapped values, size 3 x N

    fd_in = gismo.gsFileData["real_t"](args.fn)
    uv = fd_in.getId[gismo.gsMatrix["real_t"]](0)
    xyz = fd_in.getId[gismo.gsMatrix["real_t"]](1)


    if not (uv.cols() == xyz.cols() and uv.rows() == 2 and xyz.rows() == 3):
        raise TypeError("Wrong input")


    # Determine the parameter domain by min/max of parameter values
    u_min = uv.row(0).minCoeff()
    u_max = uv.row(0).maxCoeff()
    v_min = uv.row(1).minCoeff()
    v_max = uv.row(1).maxCoeff()

    # Create knot vectors without interior knots
    u_knots = gismo.gsKnotVector["real_t"](u_min, u_max, 0, args.deg_x + 1)
    v_knots = gismo.gsKnotVector["real_t"](v_min, v_max, 0, args.deg_y + 1)

    # Create a tensor basis and apply initial uniform refinement
    T_tbasis = gismo.gsTensorBSplineBasis[2, "real_t"](u_knots, v_knots)
    for _ in range(args.numURef):
        T_tbasis.uniformRefine()

    # Create initial hierarchical basis
    THB = gismo.gsTHBSplineBasis[2, "real_t"](T_tbasis)


    # Specify extension size in u and v cells
    ext = [args.extension, args.extension]


    # Create hierarchical refinement object
    ref = gismo.gsHFitting[2, "real_t"](uv, xyz, THB, args.refPercent, ext, args.smoothingcoef)

    errors = ref.pointWiseErrors()

    # Print settings summary
    print(f"Fitting {xyz.cols()} samples.\n"
          f"----------------\n"
          f"Cell extension     : {ext}.")
    if args.threshold >= 0.0:
        print(f"Ref. threshold     : {args.threshold}.")
    else:
        print(f"Cell refinement    : {100 * args.refPercent}%.\n"
              f"Error tolerance    : {args.tolerance}\n"
              f"Smoothing parameter: {args.smoothingcoef}")

    for i in range(args.iter + 1):
        print(f"----------------\n"
              f"Iteration {i}")

        starttime = time.time()
        ref.nextIteration(args.tolerance, args.threshold, args.maxPcIter)
        runtime = time.time() - starttime

        print(f"Fitting time: {time}")
        print(f"Fitted with {ref.result().basis()}")
        print(f"Min distance : {ref.minPointError()} / Max distance: {ref.maxPointError()}")
        print(f"Points below tolerance {100. * ref.numPointsBelow(args.tolerance) / len(errors)}%")

        if ref.maxPointError() < args.tolerance:
            print(f"Error tolerance achieved after {i} iterations.")
            break

    print("----------------")

    if args.save:
        print("Done. Writing solution to file fitting_out.xml")
        fd = gismo.gsFileData["real_t"]()
        fd.add(ref.result())
        fd.dump("fitting_out")
    else:
        print("Done. No output created, re-run with --save to get a xml file containing the solution.")



if __name__ == "__main__":
    main()
