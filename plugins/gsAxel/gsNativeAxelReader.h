/* goSurfaceBSplineReader.h ---
 *
 * Author: Julien Wintz
 * Copyright (C) 2008-2011 - Julien Wintz, Inria.
 * Created: Wed Sep 21 11:13:23 2011 (+0200)
 * Version: $Id$
 * Last-Updated: Wed Sep 21 11:49:25 2011 (+0200)
 *           By: Julien Wintz
 *     Update #: 21
 */

/* Commentary:
 *
 */

/* Change log:
 *
 */

#ifndef GOSURFACEBSPLINEREADER_H
#define GOSURFACEBSPLINEREADER_H

#include <axlCore/axlAbstractDataReader.h>

#include "gsAxelPluginExport.h"

class dtkAbstractData;

class GSAXELPLUGIN_EXPORT gsNativeAxelReader : public axlAbstractDataReader
{
    Q_OBJECT

public :
      gsNativeAxelReader(void);
     ~gsNativeAxelReader(void);

public:
    QString identifier(void) const;
    QString description(void) const;
    QStringList handled(void) const;

    static bool registered(void);

public:
    bool accept(const QDomNode& node);
    bool reject(const QDomNode& node);

    dtkAbstractData *read(const QDomNode& node);
};

dtkAbstractDataReader *creategsNativeAxelReader(void);

#endif //GOSURFACEBSPLINEREADER_H
