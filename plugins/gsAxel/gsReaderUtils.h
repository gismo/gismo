/** @file gsReaderUtils.h

    @brief This file provides declaration reader utilities for Axel

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <QtCore>
#include <QtXml>

int getDimension(QDomElement e);

void getKnots(QDomElement e, int i, gsKnotVector<> & result);

void getControlPoints(QDomElement e, int i, int dim, gsMatrix<> & result);

QList<int> getDegrees(QDomElement e);
