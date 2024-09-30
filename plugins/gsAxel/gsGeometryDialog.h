/** @file gsGeometryDialog.h

    @brief This file provides declaration of the gsGeometryDialog for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


#include <axlGui/axlInspectorObjectFactory.h>
#include <axlGui/axlInspectorToolFactory.h>

#include "gsAxelPluginExport.h"


class dtkAbstractData;
class axlAbstractData;
class QDialogEditCP;

class gsGeometryDialogPrivate;

class GSAXELPLUGIN_EXPORT gsGeometryDialog : public axlInspectorObjectInterface
{
    Q_OBJECT
    
public:
    gsGeometryDialog(QWidget *parent = 0);
    ~gsGeometryDialog(void);
    
    QSize sizeHint(void) const;
    
    static bool registered(void);

    axlAbstractData * data();

    int selectedCp();

signals:
    
    void colorChanged(QColor color,  dtkAbstractData *data);
    
    void dataChangedByShader(dtkAbstractData *data, QString isophoteShaderXml);
    void dataChangedByOpacity(dtkAbstractData *data, double opacity);
    
    void dataChangedByColor(dtkAbstractData *data, double red, double green, double blue);
    void dataChangedByGeometry(dtkAbstractData *data);
    
    void update(void);

    void samplingChanged(void);

    void dataInserted(axlAbstractData* data);
    
public slots:
    void setData(dtkAbstractData *data);
    
    void onColorChanged(QColor color);
    void onOpacityChanged(int opacity);
    
    void onSamplingDataChanged_u(int numSamples);
    void onSamplingDataChanged_v(int numSamples);
    void onSamplingDataChanged_w(int numSamples);

    void onIndexSelected(int i);

    void openShader(void);
    void onShaderStateChanged(bool isShader);
    void onShaderChanged(QString);
    void onLineEditShaderChanged(QString);
    
    void showBasis(void);

    void saveGeometry(void);

    void refineGeometry(void);

    void insertKnot(void);

    void onEditCP(void);


    /* Joker
    void jokerPlay(void);
    //*/

    // void onShowCurrentPoint(double u, double v);
    // void onMoveCurrentPoint(double u, double v);
    // void onHideCurrentPoint(double u, double v);

private:
    void updateText();

    void initComboBoxShaderValue(void);
    void initWidget(void);
    int initOpacityValue(void);
    QString initShaderValue(void);
    QColor initColorValue(void);
    
private:

    gsGeometryDialogPrivate *d;
};

axlInspectorObjectInterface *creategsGeometryDialog(void);

