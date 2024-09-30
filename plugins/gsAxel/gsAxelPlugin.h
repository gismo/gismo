/** @file gsAxelPlugin.h

    @brief This file provides declaration of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

//#include <dtkCore/dtkPlugin.h>
// Transition files, to DTK1 (stable release version)
#include <dtkCoreSupport/dtkPlugin.h>

#include "gsAxelPluginExport.h"

#include <gismo.h>

class dtkAbstractDataFactory;
class dtkAbstractProcessFactory;

class gsAxelPluginPrivate;
class QMenu;

// the G+Smo menu

class GSAXELPLUGIN_EXPORT gsAxelPlugin : public dtkPlugin
{
    Q_OBJECT
    Q_INTERFACES(dtkPlugin)
    Q_PLUGIN_METADATA(IID "gismo.github.io" FILE "gsAxelPlugin.json")

public:
    gsAxelPlugin(QObject *parent = 0);
    ~gsAxelPlugin(void);
    
    virtual bool initialize(void);
    virtual bool uninitialize(void);
    
    virtual QString name(void) const;
    virtual QString description(void) const;
    
    virtual QStringList tags(void) const;
    virtual QStringList types(void) const;
    
public:
    static dtkAbstractDataFactory *dataFactSingleton;
    static dtkAbstractProcessFactory *processFactSingleton;
    
protected:
    void setupGismoMenu(void);
    
    
private slots:

    void aboutGismo();

private:
    gsAxelPluginPrivate *d;
};

