/** @file QDialogEditCP.h

    @brief This file provides declaration of the QDialogEditCP widget

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

#include <QDialog>

class QLabel;
class QDoubleSpinBox;
class QPushButton;

class gsGeometryDialog; // inspector interface

class QDialogEditCP : public QDialog
{
    Q_OBJECT

public:

    explicit QDialogEditCP(gsGeometryDialog * _data, QWidget *parent = 0);

public:
    
    double xCoord();
    double yCoord();
    double zCoord();
                   
public slots:
    //void onApply();

private:
    gsGeometryDialog * dialog;

    QLabel *label;
    QDoubleSpinBox *coordinatePoint_x;
    QDoubleSpinBox *coordinatePoint_y;
    QDoubleSpinBox *coordinatePoint_z;
    QPushButton * OkButton;
    QPushButton * CancelButton;
    QPushButton * ApplyButton;
};

