/** @file gsPdeSolverDialog.cpp

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include "MyPluginDefaultProcessDialog.h"

#include <axlCore/axlAbstractData.h>

#include <axlGui/axlInspectorObjectController.h>
#include <axlGui/axlInspectorObjectManagerTreeItem.h>

#include <dtkCore/dtkAbstractProcessFactory.h>
#include <dtkCore/dtkAbstractProcess.h>

#include <QtGui>

class MyPluginDefaultProcessDialogPrivate
{
    
};

MyPluginDefaultProcessDialog::MyPluginDefaultProcessDialog(QWidget *parent) : axlInspectorToolInterface(parent), d(new MyPluginDefaultProcessDialogPrivate)
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    
    QPushButton *button = new QPushButton("Run", this);
    connect(button, SIGNAL(clicked()), this, SLOT(run()));
    
    layout->addWidget(new QLabel("MyPluginDefaultProcessDialog", this));
    layout->addWidget(button);
}

MyPluginDefaultProcessDialog::~MyPluginDefaultProcessDialog(void)
{
    delete d;
    
    d = NULL;
}

bool MyPluginDefaultProcessDialog::registered(void)
{
    return axlInspectorToolFactory::instance()->registerInspectorTool("MyPluginDefaultProcess", createMyPluginDefaultProcessDialogPrivate);
}



void MyPluginDefaultProcessDialog::run(void)
{
    axlAbstractData *input1 = NULL;
    axlAbstractData *input2 = NULL;
    int numberOfDataSelected =0;
    
    foreach(axlInspectorObjectManagerTreeItem *item, m_controller->items()) {
        
        if(item->text(2) == "Selected")
        {
            if(numberOfDataSelected == 0)
                input1 = dynamic_cast<axlAbstractData *>(m_controller->data(item));
            
            if(numberOfDataSelected == 1)
                input2 = dynamic_cast<axlAbstractData *>(m_controller->data(item));
            
            numberOfDataSelected++;
        }
    }
    
    if(numberOfDataSelected == 2 && input1 && input2)
    {
        m_process->setInput(input1,0);
        m_process->setInput(input2,1);
        m_process->run();
        if(axlAbstractData *axldata = dynamic_cast<axlAbstractData *>(m_process->output()))
            emit dataInserted(axldata);
    } else {
        // In the case where not input is selected, we still run the process because in this
        // example the input data are not used.
        // But a error message should be sent.
        m_process->run();
        if(axlAbstractData *axldata = dynamic_cast<axlAbstractData *>(m_process->output()))
            emit dataInserted(axldata);
    }

}

axlInspectorToolInterface *createMyPluginDefaultProcessDialogPrivate(void)
{
    return new MyPluginDefaultProcessDialog;
}


