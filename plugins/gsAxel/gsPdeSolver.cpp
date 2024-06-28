/** @file gsPdeSolver.cpp

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include "MyPluginDefaultProcess.h"

#include "MyPluginPlugin.h"



#include <axlCore/axlAbstractData.h>

#include <axlCore/axlAbstractFieldParametricCurve.h>
#include <axlCore/axlAbstractFieldParametricSurface.h>
#include <axlCore/axlAbstractFieldParametricVolume.h>

#include <dtkCore/dtkAbstractProcessFactory.h>
#include <dtkCore/dtkAbstractDataFactory.h>

// /////////////////////////////////////////////////////////////////
// MyPluginDefaultProcessPrivate
// /////////////////////////////////////////////////////////////////

class MyPluginDefaultProcessPrivate
{
public:
    QPair<dtkAbstractData *, dtkAbstractData *> input;
    
    axlAbstractData *output;
};

// /////////////////////////////////////////////////////////////////
// MyPluginDefaultProcess
// /////////////////////////////////////////////////////////////////

MyPluginDefaultProcess::MyPluginDefaultProcess(void) : dtkAbstractProcess(), d(new MyPluginDefaultProcessPrivate)
{
    d->output = NULL;
}

MyPluginDefaultProcess::~MyPluginDefaultProcess(void)
{
    delete d;
    
    d = NULL;
}

void MyPluginDefaultProcess::setInput(dtkAbstractData *data, int channel)
{
    if(channel==0)
        d->input.first = data;
    else if(channel==1)
        d->input.second = data;
    else
        qDebug()<<"Only two channel available : 0 or 1";
}

dtkAbstractData *MyPluginDefaultProcess::output(void)
{
    return d->output;
}

int MyPluginDefaultProcess::update(void)
{
    // create a field.

    // field.input = // gsGeometryData( qui est this->d->input)

    // filed.fonction = result of solver

    // this->d->input->addField( field )

    // field.update() // Field Signal 
    
    // emit(    

    axlSphere* S = new axlSphere(0,0,0, 2);
    S->setColor(QColor(255,0,0));
    d->output=S;
    return 1;
}

bool MyPluginDefaultProcess::registered(void)
{
    return MyPluginPlugin::processFactSingleton->registerProcessType("MyPluginDefaultProcess", createMyPluginDefaultProcess, "axlAbstractDefault");
}

QString MyPluginDefaultProcess::description(void) const
{
    return "My MyPluginDefaultProcess description";
}

QString MyPluginDefaultProcess::identifier(void) const
{
    return "MyPluginDefaultProcess";
}

// /////////////////////////////////////////////////////////////////
// Type instanciation
// /////////////////////////////////////////////////////////////////

dtkAbstractProcess *createMyPluginDefaultProcess(void)
{
    return new MyPluginDefaultProcess;
}

