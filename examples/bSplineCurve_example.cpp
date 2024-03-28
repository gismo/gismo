/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, Ye Ji
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit
    bool trim = false; // If set to true, paraview file is generated and launched on exit

    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("trim", "Basic trim/merge operations", trim);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    // Make a BSpline curve
    gsKnotVector<> kv(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots1
    gsMatrix<> coefs(4, 3);
    coefs << 0, 0, 0,
             1, 2, 3,
             2, 1, 4,
             4, 4, 4;

    gsBSpline<> curve( kv, give(coefs));

    // Print the Bspline curve
    gsInfo << "I am a " << curve << "\n";

    if (plot)
    {
        // Output a paraview file
        gsWriteParaview( curve, "bsplinecurve", 100);
        gsFileManager::open("bsplinecurve.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    // Basic trim/merge operations on BSpline curves - @Ye
    if (trim)
    {
      gsInfo << "Original BSpline curve: " << curve << "\n";
      gsWriteParaview( curve, "originalCurve", 100);  // Output the original curve

      // Segment this BSpline curve between parameters 0.3 and 0.8
      gsBSpline<> segment = curve.segmentFromTo(0.3, 0.8);
      gsInfo << "Curve segment from u0 = 0.3 to u1 = 0.8: " << segment << "\n";
      gsWriteParaview(segment, "segment", 100);  // Output the curve segment

      // Split the curve at parameter 0.4 into two parts
      gsBSpline<> segmentLeft, segmentRight;
      curve.splitAt(0.4, segmentLeft, segmentRight);
      gsInfo << "Curve segment from u0 = 0.0 to u1 = 0.4: " << segmentLeft << "\n";
      gsInfo << "Curve segment from u0 = 0.4 to u1 = 1.0: " << segmentRight << "\n";
      gsWriteParaview( segmentLeft, "segmentLeft", 100);
      gsWriteParaview( segmentRight, "segmentRight", 100);

      // Merge the left and right segments back to the original curve
      // Note: Due to the segmentation, an inner knot value of 0.4 is introduced, while
      // the geometry remains exactly the same as the original one
      gsBSpline<> mergedCurve = segmentLeft;
      mergedCurve.merge(&segmentRight);
      gsInfo << "The merged curve: " << mergedCurve << "\n";
      gsWriteParaview( mergedCurve, "mergedCurve", 100);

      // convert it into bezier segments
      std::vector<gsBSpline<>> bezSegments = mergedCurve.toBezier();
      gsMultiPatch<> bezierContainer;
      for (const gsBSpline<>& bezSegment:bezSegments) {
        bezierContainer.addPatch(bezSegment);
      }
      gsWriteParaview( bezierContainer, "bezierContainer", 100);
    }
    else
      gsInfo << "Done. Re-run with --trim to learn basic trim/merge operations\n";

    return 0;
}
