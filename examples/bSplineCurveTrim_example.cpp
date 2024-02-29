/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
  gsInfo << "This example demonstrates the basic trim/merge operations on BSpline curves.\n\n";

  // Create a BSpline curve with the same geometry as in bSplineCurve_example.cpp
  gsKnotVector<> kv(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots1
  gsMatrix<> coefs(4, 3);
  coefs << 0, 0, 0,
      1, 2, 3,
      2, 1, 4,
      4, 4, 4;
  gsBSpline<> curve( kv, give(coefs));

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

  return 0;
}
