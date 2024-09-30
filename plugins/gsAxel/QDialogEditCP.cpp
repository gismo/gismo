/** @file QDialogEditCP.cpp

    @brief This file Provides implemetation of QDialogEditCP

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#include <QtWidgets>

#include "QDialogEditCP.h"
#include "gsGeometryDialog.h"

#include "gsGeometryData.h"


QDialogEditCP::QDialogEditCP(gsGeometryDialog * _dialog, QWidget *parent)
: QDialog(parent), dialog(_dialog) 
{
    label = new QLabel(tr("Coordinates:"));

    gismo::gsGeometry<> * g = getGeometryPointer(dialog->data()).get();
    const int cp = dialog->selectedCp();

    coordinatePoint_x = new QDoubleSpinBox(this);
    coordinatePoint_x->setRange(-1e6, 1e6);
    coordinatePoint_x->setDecimals(16);
    coordinatePoint_x->setValue(g->coef(cp,0));
    coordinatePoint_x->setSingleStep(0.1);

    const int gdim = g->geoDim();
    if ( gdim > 1)
    {
        coordinatePoint_y = new QDoubleSpinBox(this);
        coordinatePoint_y->setRange(-1e6, 1e6);
        coordinatePoint_y->setDecimals(16);
        coordinatePoint_y->setValue(g->coef(cp,1));
        coordinatePoint_y->setSingleStep(0.1);
        if ( gdim > 2)
        {
            coordinatePoint_z = new QDoubleSpinBox(this);
            coordinatePoint_z->setRange(-1e6, 1e6);
            coordinatePoint_z->setDecimals(16);
            coordinatePoint_z->setValue(g->coef(cp,2));
            coordinatePoint_z->setSingleStep(0.1);
        }
    }

    OkButton = new QPushButton(tr("&OK"),this);
    OkButton->setDefault(true);
    CancelButton = new QPushButton(tr("&Cancel"),this);
    //ApplyButton  = new QPushButton(tr("&Apply") ,this);

    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    mainLayout->addWidget(label);
    mainLayout->addWidget(coordinatePoint_x);
    mainLayout->addWidget(coordinatePoint_y);
    mainLayout->addWidget(coordinatePoint_z);
    mainLayout->addWidget(OkButton    );
    mainLayout->addWidget(CancelButton);
    
    setLayout(mainLayout);

    setWindowTitle(tr("Edit Control Point"));

    connect(OkButton    , SIGNAL(clicked()), this  , SLOT(accept()));
    connect(CancelButton, SIGNAL(clicked()), this  , SLOT(reject()));
}

double QDialogEditCP::xCoord()
{return coordinatePoint_x->value();}

double QDialogEditCP::yCoord()
{return coordinatePoint_y->value();}

double QDialogEditCP::zCoord()
{return coordinatePoint_z->value();}
