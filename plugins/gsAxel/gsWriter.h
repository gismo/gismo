/** @file gsWriter.h

    @brief This file provides declaration of the XML writer for Axel

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include "gsAxelPluginExport.h"

#include <dtkCoreSupport/dtkAbstractDataWriter.h>

//#include <axlCore/axlCoreExport.h>


class dtkAbstractData;
class gsWriterPrivate;


class GSAXELPLUGIN_EXPORT gsWriter : public dtkAbstractDataWriter
{
    Q_OBJECT

public :
    gsWriter(void);
    ~gsWriter(void);

public:
    QString identifier(void) const;
    QString description(void) const;
    QStringList handled(void) const;

    static bool registered(void);

public:
    bool canWrite(const QString& file);
    bool write(const QString& file);


private:
    gsWriterPrivate *d;

};

dtkAbstractDataWriter *creategsWriter(void);
