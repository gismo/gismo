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
#include <gsIO/gsParaview.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit
    bool trim = false; // If set to true, trim/merge operations are displayed
    bool intersect = false; // If set to true, intersection example is displayed

    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("trim", "Basic trim/merge operations", trim);
    cmd.addSwitch("intersect", "Intersection operations", intersect);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [make curve]
    // Make a BSpline curve
    gsKnotVector<> kv(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots
    gsMatrix<> coefs(4, 3);
    coefs << 0, 0, 0,
             1, 2, 3,
             2, 1, 4,
             4, 4, 4;

    gsBSpline<> curve(kv, coefs);
    //! [make curve]

    // Print the Bspline curve
    gsInfo << "I am a " << curve << "\n";

    if (plot)
    {   
        // Output a paraview file
        coefs.transposeInPlace();
        gsParaview<real_t> pv;
        pv.options().setInt("numPoints", 100);
        pv.options().setSwitch("show", true);
        pv.write(curve, "bsplinecurve0");
        pv.options().setSwitch("plotElements", true);
        pv.options().setSwitch("plotControlNet", true);
        pv.write(curve, "bsplinecurve");
        pv.writePoints(coefs, "coefficients");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    // Basic trim/merge operations on BSpline curves - @Ye
    if (trim)
    {
      gsInfo << "Original BSpline curve: " << curve << "\n";
      gsParaview<real_t> pv;
      pv.options().setInt("numPoints", 100);
      pv.write(curve, "originalCurve");  // Output the original curve

      // Segment this BSpline curve between parameters 0.3 and 0.8
      gsBSpline<> segment = curve.segmentFromTo(0.3, 0.8);
      gsInfo << "Curve segment from u0 = 0.3 to u1 = 0.8: " << segment << "\n";
      pv.write(segment, "segment");  // Output the curve segment

      // Split the curve at parameter 0.4 into two parts
      gsBSpline<> segmentLeft, segmentRight;
      curve.splitAt(0.4, segmentLeft, segmentRight);
      gsInfo << "Curve segment from u0 = 0.0 to u1 = 0.4: " << segmentLeft << "\n";
      gsInfo << "Curve segment from u0 = 0.4 to u1 = 1.0: " << segmentRight << "\n";
      pv.write(segmentLeft, "segmentLeft");
      pv.write(segmentRight, "segmentRight");

      // Merge the left and right segments back to the original curve
      // Note: Due to the segmentation, an inner knot value of 0.4 is introduced, while
      // the geometry remains exactly the same as the original one
      gsBSpline<> mergedCurve = segmentLeft;
      mergedCurve.merge(&segmentRight);
      gsInfo << "The merged curve: " << mergedCurve << "\n";
      pv.write(mergedCurve, "mergedCurve");

      // convert it into bezier segments
      gsMultiPatch<> bezSegments = mergedCurve.toBezier();
      pv.write(bezSegments, "bezierContainer");
    }
    else
      gsInfo << "Done. Re-run with --trim to learn basic trim/merge operations\n";

    // Basic intersection operations between two BSpline curves - @Ye
    if (intersect)
    {
      gsMatrix<real_t> ctrPts1(4, 2);
      ctrPts1 << 0,0, 1,1, 2,1, 3,1;
      gsBSpline<real_t> bsp1(0, 1, 0, 3, ctrPts1);
      gsMatrix<real_t> ctrPts2(4, 2);
      ctrPts2 << 0,0, 1,2, 2,2, 3,0;
      gsBSpline<real_t> bsp2(0, 1, 0, 3, ctrPts2);

      auto intersectPts = bsp1.intersect(bsp2, 1e-5);

      gsInfo << intersectPts.size() << " intersections are found!" << "\n";
      gsMatrix<> iPts(bsp1.geoDim(), intersectPts.size());
      for (size_t j = 0; j < intersectPts.size(); ++j)
      {
        iPts.col(j) = intersectPts[j].getPoint();
      }
      if (!intersectPts.empty())
      {
        gsParaview<real_t> pv;
        pv.writePoints(iPts, "intersect");
      }

      gsParaview<real_t> pv2;
      pv2.options().setInt("numPoints", 2000);
      pv2.write(bsp1, "bsp1");
      pv2.write(bsp2, "bsp2");
    }
    else
      gsInfo << "Done. Re-run with --intersect to learn intersection operations\n";

    return 0;
}
