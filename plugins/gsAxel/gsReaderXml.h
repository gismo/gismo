/** @file gsReaderXml.h

    @brief This file provides declaration of the XML reader for Axel

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include "gsAxelPluginExport.h"

#include <dtkCoreSupport/dtkAbstractDataReader.h>

class GSAXELPLUGIN_EXPORT gsReaderXml : public dtkAbstractDataReader
{
    Q_OBJECT
    
public :
    gsReaderXml(void);
    ~gsReaderXml(void);
    
public:
    QString identifier(void) const;
    QString description(void) const;
    QStringList handled(void) const;
    
    static bool registered(void);
    
public:

    bool canRead(const QString& file );
    bool read   (const QString& file );

};

dtkAbstractDataReader *creategsReaderXml(void);
