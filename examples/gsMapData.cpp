/** @file gsMapData.cpp

    @brief A file where Dominik learns to use gsMapData.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#include <gismo.h>

using namespace gismo;

int main()
{
    gsTensorBSpline<2> tpBSpline;

    gsFileData<> fd("surfaces/teapot.xml");
    fd.getId<gsTensorBSpline<2>>(0, tpBSpline);

    gsMapData<> md(NEED_NORMAL);

    gsMatrix<real_t> points(2, 1);
    points(0, 0) = 0.2;
    points(1, 0) = 0.7;
    md.points = points;

    tpBSpline.computeMap(md);
    gsInfo << md.normals << std::endl;

    return 0;
}
