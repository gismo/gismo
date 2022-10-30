/** @file gsReaderAxl.h

    @brief This file provides declaration of the XML reader for Axl

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include "gsAxelPluginExport.h"

#include <axlCore/axlAbstractDataReader.h>

class GSAXELPLUGIN_EXPORT gsReaderAxl : public axlAbstractDataReader
{
    Q_OBJECT
    
public :
    gsReaderAxl(void);
    ~gsReaderAxl(void);
    
public:
    QString identifier(void) const;
    QString description(void) const;
    QStringList handled(void) const;
    
    static bool registered(void);
    
public:
    bool accept(const QDomNode& node);
    bool reject(const QDomNode& node);

    axlAbstractData *read(const QDomNode& node);
};

dtkAbstractDataReader *creategsReaderAxl(void);
