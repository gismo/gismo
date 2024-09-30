/** @file gsPdeSolver.h

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
#ifndef MYPLUGINDEFAULTPROCESS_H
#define MYPLUGINDEFAULTPROCESS_H

#include "MyPluginPluginExport.h"

#include <dtkCore/dtkAbstractProcess.h>

class axlAbstractData;

class MyPluginDefaultProcessPrivate;

class MYPLUGINPLUGIN_EXPORT MyPluginDefaultProcess: public dtkAbstractProcess
{
    Q_OBJECT
    
public:
    MyPluginDefaultProcess(void);
    virtual ~MyPluginDefaultProcess(void);
    
    virtual QString description(void) const;
    virtual QString identifier(void) const;
    
    
public:
    static bool registered(void);
    
public:
    dtkAbstractData *output(void);
    
public slots:
    void setInput(dtkAbstractData *data, int channel);
    
public slots:
    int update(void);
    
private:
    MyPluginDefaultProcessPrivate *d;
};

dtkAbstractProcess *createMyPluginDefaultProcess(void);

#endif //MYPLUGINDEFAULTPROCESS_H

