/** @file gsPdeSolverDialog.h

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#ifndef MYPLUGINDEFAULTPROCESSDIALOG_H
#define MYPLUGINDEFAULTPROCESSDIALOG_H

#include "MyPluginPluginExport.h"

#include <axlGui/axlInspectorToolFactory.h>


class axlAbstractData;

class MyPluginDefaultProcessDialogPrivate;

class MYPLUGINPLUGIN_EXPORT MyPluginDefaultProcessDialog : public axlInspectorToolInterface
{
    Q_OBJECT
    
public:
    MyPluginDefaultProcessDialog(QWidget *parent = 0);
    ~MyPluginDefaultProcessDialog(void);
    
    static bool registered(void);
    
signals:
    void dataInserted(axlAbstractData *data);
    
    
public slots:
    void run(void);
    
private:
    MyPluginDefaultProcessDialogPrivate *d;
};

axlInspectorToolInterface *createMyPluginDefaultProcessDialogPrivate(void);


#endif //MYPLUGINDEFAULTPROCESSDIALOG_H

