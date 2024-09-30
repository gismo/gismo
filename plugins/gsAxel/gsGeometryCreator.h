/** @file gsGeometryCreator.h

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
#pragma once

#include <axlGui/axlInspectorToolFactory.h>

#include "gsAxelPluginExport.h"

class gsGeometryCreatorPrivate;

class GSAXELPLUGIN_EXPORT gsGeometryCreator : public axlInspectorToolInterface {
    Q_OBJECT
    
public:
    gsGeometryCreator(QWidget* parent = 0);
    virtual ~gsGeometryCreator();
    
    static bool registered(void);
    
signals:
    void dataInserted(axlAbstractData* data);
    
public slots:
    void run(void);
    void loadBasis(void);
    void loadGeometry(void);
    void loadMultiPatch(void);

    void update_index(int i);
    void update_degree(int i);

private:
    gsGeometryCreatorPrivate* d;
    
};

dtkAbstractProcess* createProcessCreatorgsGeometryData(void);
axlInspectorToolInterface* creategsGeometryCreator(void);

