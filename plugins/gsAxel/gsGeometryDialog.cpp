/** @file gsGeometryDialog.cpp

    @brief This file Provides implemetation of G+Smo geometry dialog for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include "gsGeometryDialog.h"

#include "gsGeometryData.h"
#include "gsBasisData.h"
#include "gsGeometryCreator.h"

#include "QDialogEditCP.h"

#include <axlCore/axlAbstractData.h>
#include <axlCore/axlMenuFactory.h>
#include <axlGui/axlInspectorUtils.h>
#include <dtkCoreSupport/dtkAbstractData.h>

#include <dtkGuiSupport/dtkColorButton.h>
#include <dtkGuiSupport/dtkSplitter.h>

#include <QtGui>

#include <gsIO/gsFileData.h>

#define DEFAULT_SAMPLES 20

// Process
#include <axlCore/axlAbstractFieldGenerator.h>



//
#define CALL_DATA_METHOD(method) \
if ( gsAxelCurve * obj = dynamic_cast<gsAxelCurve *>(d->data) ) \
        obj-> method ();  else \
    if ( gsAxelSurface * obj = dynamic_cast<gsAxelSurface *>(d->data) ) \
        obj-> method (); else \
    if ( gsAxelVolume * obj = dynamic_cast<gsAxelVolume *>(d->data) ) \
        obj-> method (); else \
    if ( gsAxelTrimSurf * obj = dynamic_cast<gsAxelTrimSurf *>(d->data) ) \
        obj-> method ();
#define CALL_DATA_METHOD1(method,arg1)  \
if ( gsAxelCurve * obj = dynamic_cast<gsAxelCurve *>(d->data) ) \
        obj-> method (arg1);  else \
    if ( gsAxelSurface * obj = dynamic_cast<gsAxelSurface *>(d->data) ) \
        obj-> method (arg1); else \
    if ( gsAxelVolume * obj = dynamic_cast<gsAxelVolume *>(d->data) ) \
        obj-> method (arg1); else \
    if ( gsAxelTrimSurf * obj = dynamic_cast<gsAxelTrimSurf *>(d->data) ) \
        obj-> method (arg1);
//

class gsGeometryDialogPrivate
{
public:
    axlAbstractData *data;

    dtkColorButton *colorButton;
    
    QComboBox *comboBoxShader;
    QCheckBox *checkBoxShader;
    QLineEdit *lineEditShader;
    QPushButton *buttonShader;    
    QSlider *sliderOpacity;

    QSpinBox *spinBoxSampling_u;
    QSpinBox *spinBoxSampling_v;
    QSpinBox *spinBoxSampling_w;

    // Control point box
    QSpinBox       *cpIndex;
    QPushButton    *cpEdit;

    QPushButton *basisButton;    

    QPushButton *saveButton;    

    QPushButton * refineButton;

    QPushButton * insertKnotButton;

    /* Joker
    QPushButton * jokerButton;
    //*/

    QTextEdit *info;

};

gsGeometryDialog::gsGeometryDialog(QWidget *parent) : axlInspectorObjectInterface(parent), d(new gsGeometryDialogPrivate)
{
    d->data = NULL;
    
    d->comboBoxShader = NULL;
    d->sliderOpacity = NULL;
    d->checkBoxShader = NULL;
    d->lineEditShader = NULL;
    d->buttonShader = NULL;
    d->colorButton = NULL;
    d->info = NULL;
}

void gsGeometryDialog::initWidget(void) {

    //OPACITY//
    d->sliderOpacity = new QSlider(Qt::Horizontal, this);
    
    QHBoxLayout *layoutOpacity = new QHBoxLayout;
    layoutOpacity->addWidget(new QLabel("Opacity",this));
    layoutOpacity->addWidget(d->sliderOpacity);
    d->sliderOpacity->setMaximum(100);
    d->sliderOpacity->setValue(initOpacityValue());
    
    //SHADER//
    d->comboBoxShader = new QComboBox(this);
    d->comboBoxShader->setInsertPolicy(QComboBox::InsertAlphabetically);
    d->checkBoxShader = new QCheckBox(this);
    d->lineEditShader = new QLineEdit(this);
    d->buttonShader = new QPushButton(this);
    d->buttonShader->setText("open");
    
    d->lineEditShader->setText(initShaderValue());
    this->initComboBoxShaderValue();
    
    if(d->lineEditShader->text().isEmpty())
    {
        d->lineEditShader->setEnabled(false);
        d->buttonShader->setEnabled(false);
        d->comboBoxShader->setEnabled(false);
    }
    else
        d->checkBoxShader->setChecked(true);
    
    //COLOR//
    d->colorButton = new dtkColorButton(this);
    
    QHBoxLayout *layoutColorButton = new QHBoxLayout;
    layoutColorButton->addWidget(new QLabel("Color",this));
    layoutColorButton->addWidget(d->colorButton);
    d->colorButton->setColor(this->initColorValue());
    
    //GUI//
    QVBoxLayout *layoutShader = new QVBoxLayout;
    QHBoxLayout *layoutShader1 = new QHBoxLayout;
    
    QLabel *labelShader  = new QLabel("Shader",this);
    layoutShader1->addWidget(labelShader);
    layoutShader1->addWidget(d->checkBoxShader);
    layoutShader1->addWidget(d->comboBoxShader);
    layoutShader1->addWidget(d->buttonShader);
    
    layoutShader1->setStretchFactor(labelShader, 2);
    layoutShader1->setStretchFactor(d->checkBoxShader, 1);
    layoutShader1->setStretchFactor(d->comboBoxShader, 4);
    layoutShader1->setStretchFactor(d->buttonShader, 3);
    
    layoutShader->addLayout(layoutShader1);
    layoutShader->addWidget(d->lineEditShader);
    QVBoxLayout *layoutTop = new QVBoxLayout;
    layoutTop->addWidget(new QLabel("Gismo Geometry Properties", this));
    layoutTop->addLayout(layoutColorButton);
    layoutTop->addLayout(layoutOpacity);
    layoutTop->addLayout(layoutShader);

    // Sampling
    QHBoxLayout *layoutSampling = new QHBoxLayout;
    layoutSampling->addWidget(new QLabel("Sampling",this));
    d->spinBoxSampling_u = new QSpinBox(this);
    d->spinBoxSampling_u->setMaximum(5000);
    d->spinBoxSampling_u->setMinimum(5);
    d->spinBoxSampling_u->setValue(DEFAULT_SAMPLES);
    layoutSampling->addWidget(d->spinBoxSampling_u);
    connect(d->spinBoxSampling_u, SIGNAL(valueChanged(int)), this, SLOT(onSamplingDataChanged_u(int)));
    int pdim = getGeometryPointer(d->data)->parDim();
    if (  pdim > 1)
	{
	    d->spinBoxSampling_v = new QSpinBox(this);
	    d->spinBoxSampling_v->setMaximum(5000);
	    d->spinBoxSampling_v->setMinimum(5);
	    d->spinBoxSampling_v->setValue(DEFAULT_SAMPLES);
	    layoutSampling->addWidget(d->spinBoxSampling_v);
	    connect(d->spinBoxSampling_v, SIGNAL(valueChanged(int)), this, SLOT(onSamplingDataChanged_v(int)));
	    if ( pdim > 2)
		{	    
		    d->spinBoxSampling_w = new QSpinBox(this);
		    d->spinBoxSampling_w->setMaximum(5000);
		    d->spinBoxSampling_w->setMinimum(5);
		    d->spinBoxSampling_w->setValue(DEFAULT_SAMPLES);
		    layoutSampling->addWidget(d->spinBoxSampling_w);
		    connect(d->spinBoxSampling_w, SIGNAL(valueChanged(int)), this, SLOT(onSamplingDataChanged_w(int)));
		}
	}
    layoutTop->addLayout(layoutSampling);

    // CONTROL POINT
    //gismo::gsGeometry<> * g = getGeometryPointer(d->data);
    d->cpIndex           = new QSpinBox(this);
    //d->cpIndex->setButtonSymbols(QAbstractSpinBox::NoButtons);
    d->cpIndex->setMinimum(0  );
    d->cpIndex->setMaximum(1e6);
    QHBoxLayout *layout_cpid = new QHBoxLayout;
    layout_cpid->addWidget(new QLabel("Control point",this));
    layout_cpid->addWidget(d->cpIndex);
    layoutTop->addLayout(layout_cpid);
    d->cpEdit = new QPushButton("Edit",this);
    layoutTop->addWidget(d->cpEdit);
    layout_cpid->addWidget(d->cpEdit);
    connect(d->cpEdit, SIGNAL(clicked()), this, SLOT(onEditCP()) );
/*
    d->coordinatePoint_x = new QDoubleSpinBox(this);
    d->coordinatePoint_x->setRange(-1e6, 1e6);
    d->coordinatePoint_x->setValue(g->coef(0,0));
    d->coordinatePoint_x->setSingleStep(0.1);
    d->coordinatePoint_y = new QDoubleSpinBox(this);
    d->coordinatePoint_y->setRange(-1e6, 1e6);
    d->coordinatePoint_y->setValue(g->coef(0,1));
    d->coordinatePoint_y->setSingleStep(0.1);
    d->coordinatePoint_z = new QDoubleSpinBox(this);
    d->coordinatePoint_z->setRange(-1e6, 1e6);
    d->coordinatePoint_z->setValue(g->coef(0,2));
    d->coordinatePoint_z->setSingleStep(0.1);
    QHBoxLayout *layout_coord = new QHBoxLayout;
    layout_coord->addWidget(d->coordinatePoint_x);
    //QHBoxLayout *layout_coord_y = new QHBoxLayout;
    //layout_cpid->addWidget(new QLabel("y",this));
    layout_cpid->addWidget(d->coordinatePoint_y);
    layout_cpid->addWidget(d->coordinatePoint_z);
    layoutTop->addLayout(layout_coord);
*/
    // when changing the value, grab the new index
    //connect(d->cpIndex, SIGNAL(valueChanged(int)), this, SLOT(onControlPointIndexChanged(int)));

    // When selecting a control point inside the view, update the coordinate boxes
    connect(d->data, SIGNAL(indexSelected(int)), this, SLOT(onIndexSelected(int)));

    //Basis button
    d->basisButton = new QPushButton("Show basis",this);
    layoutTop->addWidget(d->basisButton);
    connect(d->basisButton, SIGNAL(clicked()), this, SLOT(showBasis()));

    //Save button
    d->saveButton = new QPushButton("Save geometry",this);
    layoutTop->addWidget(d->saveButton);
    connect(d->saveButton, SIGNAL(clicked()), this, SLOT(saveGeometry()));

    //Refine button
    d->refineButton = new QPushButton("Refine here",this);
    layoutTop->addWidget(d->refineButton);
    connect(d->refineButton, SIGNAL(clicked()), this, SLOT(refineGeometry()));

    // Insert knot button
    d->insertKnotButton = new QPushButton("Insert knot",this);
    layoutTop->addWidget(d->insertKnotButton);
    connect(d->insertKnotButton, SIGNAL(clicked()), this, SLOT(insertKnot()));

    /* Joker
    d->jokerButton = new QPushButton("Joker",this);
    layoutTop->addWidget(d->jokerButton);
    connect(d->jokerButton, SIGNAL(clicked()), this, SLOT(jokerPlay()));
    //*/

    d->info = new QTextEdit("G+Smo");
    d->info->setTextInteractionFlags(Qt::NoTextInteraction);
    updateText();

    layoutTop->addWidget(d->info);
    // QSize myEditSize = d->info->document()->size().toSize();
    // myEditSize.setWidth(QWIDGETSIZE_MAX);
    // d->info->setMaximumSize(myEditSize);

    layoutTop->addStretch(1);    
    QWidget *top = new QWidget(this);
    top->setMaximumWidth(295);
    top->setLayout(layoutTop);

    connect(d->comboBoxShader, SIGNAL(currentIndexChanged(QString)), this, SLOT(onLineEditShaderChanged(QString)));
    connect(d->checkBoxShader, SIGNAL(clicked(bool)), this, SLOT(onShaderStateChanged(bool)));
    connect(d->buttonShader, SIGNAL(clicked()), this, SLOT(openShader()));
    connect(d->colorButton, SIGNAL(colorChanged(QColor)), this, SLOT(onColorChanged(QColor)));
    connect(d->sliderOpacity, SIGNAL(sliderMoved(int)), this, SLOT(onOpacityChanged(int)));
    connect(d->lineEditShader, SIGNAL(textChanged(QString)), this, SLOT(onShaderChanged(QString)));



    // G+Smo menu
    QMenu * m = axlMenuFactory::instance()->menus().at(0);
    QAction * a = m->actions().at(0);
    a->setEnabled(true);
    connect(a, SIGNAL(triggered()), this, SLOT(refineGeometry()));
    // note: function ptr can go there instead of SLOT
}


gsGeometryDialog::~gsGeometryDialog(void)
{
    delete d;
    
    d = NULL;

    QAction * a = axlMenuFactory::instance()->menus().at(0)->actions().at(0);
    a->setEnabled(false);
}

bool gsGeometryDialog::registered(void)
{
    return 
        axlInspectorObjectFactory::instance()->registerInspectorObject("SplineCurve", creategsGeometryDialog)
        &&
        axlInspectorObjectFactory::instance()->registerInspectorObject("SplineSurface", creategsGeometryDialog)
        &&
        axlInspectorObjectFactory::instance()->registerInspectorObject("TrimSurface", creategsGeometryDialog)
        &&
        axlInspectorObjectFactory::instance()->registerInspectorObject("SplineVolume", creategsGeometryDialog);
}

axlAbstractData * gsGeometryDialog::data(void)
{
    return d->data;
}

int gsGeometryDialog::selectedCp()
{
    return d->cpIndex->value();
}


axlInspectorObjectInterface *creategsGeometryDialog(void)
{
    return new gsGeometryDialog;
}

QSize gsGeometryDialog::sizeHint(void) const
{
    return QSize(600, 600);
}

void gsGeometryDialog::updateText()
{
    std::stringstream s;
    s << * getGeometryPointer(d->data) ;
    d->info->setText( s.str().c_str() );
}


void gsGeometryDialog::setData(dtkAbstractData *data)
{
    if(
          (d->data = dynamic_cast<gsAxelCurve   *>(data))  
       || (d->data = dynamic_cast<gsAxelSurface *>(data))
       || (d->data = dynamic_cast<gsAxelVolume  *>(data))
       || (d->data = dynamic_cast<gsAxelTrimSurf*>(data))
       )
    {
        initWidget();
    }
}

void gsGeometryDialog::onIndexSelected(int i)
{
    //--i;// 1-based numbering? --> no.
    d->cpIndex->setValue(i);
}

void gsGeometryDialog::onSamplingDataChanged_u(int numSamples)
{
    //CALL_DATA_METHOD(
    if ( gsAxelCurve * obj = dynamic_cast<gsAxelCurve *>(d->data) )
            obj->setNumSamples_u(numSamples);
    else 
    if ( gsAxelSurface * obj = dynamic_cast<gsAxelSurface *>(d->data) )
            obj->setNumSamples_u(numSamples);
    else
    if ( gsAxelVolume * obj = dynamic_cast<gsAxelVolume *>(d->data) )
            obj->setNumSamples_u(numSamples);
    else
    if ( gsAxelTrimSurf * obj = dynamic_cast<gsAxelTrimSurf *>(d->data) )
            obj->setNumSamples_u(numSamples);
    emit update();
}

void gsGeometryDialog::onSamplingDataChanged_v(int numSamples)
{
    if ( gsAxelCurve * obj = dynamic_cast<gsAxelCurve *>(d->data) )
        obj->setNumSamples_v(numSamples);
    else 
        if ( gsAxelSurface * obj = dynamic_cast<gsAxelSurface *>(d->data) )
            obj->setNumSamples_v(numSamples);
        else
            if ( gsAxelVolume * obj = dynamic_cast<gsAxelVolume *>(d->data) )
                obj->setNumSamples_v(numSamples);
            else
                if ( gsAxelTrimSurf * obj = dynamic_cast<gsAxelTrimSurf *>(d->data) )
                    obj->setNumSamples_v(numSamples);
    emit update();
}

void gsGeometryDialog::onSamplingDataChanged_w(int numSamples)
{
    if ( gsAxelCurve * obj = dynamic_cast<gsAxelCurve *>(d->data) )
            obj->setNumSamples_w(numSamples);
    else 
    if ( gsAxelSurface * obj = dynamic_cast<gsAxelSurface *>(d->data) )
            obj->setNumSamples_w(numSamples);
    else
    if ( gsAxelVolume * obj = dynamic_cast<gsAxelVolume *>(d->data) )
            obj->setNumSamples_w(numSamples);
    else
    if ( gsAxelTrimSurf * obj = dynamic_cast<gsAxelTrimSurf *>(d->data) )
            obj->setNumSamples_w(numSamples);
    emit update();
}


void gsGeometryDialog::onColorChanged(QColor color)
{
    d->data->setColor(color);
    
    emit dataChangedByColor(d->data, color.redF(), color.greenF(), color.blueF());
    
    emit update();
}

void gsGeometryDialog::initComboBoxShaderValue(void)
{
    if(d->comboBoxShader)
    {
        // First add item of axlShader.qrc, then find shader from shader path
        QDir dirShader( ":axlShader/shader/");
        dirShader.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
        
        QFileInfoList list = dirShader.entryInfoList();
        //        for (int i = 0; i < list.size(); ++i) {
        //            d->comboBoxShader->addItem(list.at(i).fileName());
        //        }
        
        QSettings settings("inria", "dtk");
        QString defaultPath;
        settings.beginGroup("shader");
        QString defaultPathShader = settings.value("path", defaultPath).toString();
        defaultPathShader.append("/");
        
        QDir defaultDirShader(defaultPathShader);
        QStringList filters;
        filters << "*.xml";
        defaultDirShader.setNameFilters(filters);
        QFileInfoList list2 = defaultDirShader.entryInfoList();
        
        list.append(list2);
        
        QStringList items;
        
        for (int i = 0; i < list.size(); ++i) {
            if(!items.contains(list.at(i).fileName()))
                items << list.at(i).fileName();
        }
        
        qSort(items.begin(), items.end(), caseInsensitiveLessThan);
        int indInitShader = -1;
        int indCurrentShader = -1;
        
        
        foreach(QString item, items)
        {
            indCurrentShader++;
            d->comboBoxShader->addItem(item);
            
            QFileInfo currentFileInfo(d->lineEditShader->text());
            
            if(currentFileInfo.exists())
            {
                if(item == currentFileInfo.fileName())
                    indInitShader =indCurrentShader;
            }
        }
        
        //init the value from the lineEditShader.
        if(indInitShader != -1)
            d->comboBoxShader->setCurrentIndex(indInitShader);
        
    }
}

void gsGeometryDialog::onLineEditShaderChanged(QString shader)
{
    // First add item of axlShader.qrc, then find shader from shader path
    QDir dirShader( ":axlShader/shader/");
    dirShader.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
    
    QFileInfo currentFile(dirShader, shader);
    if(!currentFile.exists())
    {
        QSettings settings("inria", "dtk");
        QString defaultPath;
        settings.beginGroup("shader");
        QString defaultPathShader = settings.value("path", defaultPath).toString();
        defaultPathShader.append("/");
        
        QDir defaultDirShader(defaultPathShader);
        currentFile = QFileInfo(defaultDirShader, shader);
        
    }
    
    d->lineEditShader->setText(currentFile.absoluteFilePath());
}

void gsGeometryDialog::openShader()
{
    if(d->lineEditShader->isEnabled())
    {
        QString fileToOpen;
        fileToOpen = QFileDialog::getOpenFileName(this, tr("Open shader"), "", tr("xml document (*.xml)"));
        d->lineEditShader->setText(fileToOpen);
    }
}

void gsGeometryDialog::onShaderChanged(QString shader)
{
    d->data->setShader(shader);
    
    //emit dataChangedByShader(d->data, d->lineEditShader->text());
    emit modifiedProperty(d->data, 2);

    emit update();
}

void gsGeometryDialog::onShaderStateChanged(bool isShader)
{
    if(isShader)
    {
        d->comboBoxShader->setEnabled(true);
        d->lineEditShader->setEnabled(true);
        d->buttonShader->setEnabled(true);
        emit dataChangedByShader(d->data, d->lineEditShader->text());
    }
    else
    {
        d->comboBoxShader->setEnabled(false);
        d->lineEditShader->setEnabled(false);
        d->buttonShader->setEnabled(false);
        emit dataChangedByShader(d->data, "");
    }
    emit update();
}

void gsGeometryDialog::onOpacityChanged(int opacity)
{
    double opacity_d = 1.0 - 0.01 * opacity; // range from 0.00 to 1.00
    d->data->setOpacity(opacity_d);
    
    //emit dataChangedByOpacity(d->data, opacity_d);
    emit modifiedProperty(d->data, 1);
    emit update();
}

int gsGeometryDialog::initOpacityValue(void)
{
    double initOpacity = 0.0;
    double opacity = d->data->opacity();
    if(opacity > initOpacity)
        initOpacity = opacity;
    
    return 100 * (1.0 - initOpacity);
}

QString gsGeometryDialog::initShaderValue(void)
{
    return  d->data->shader();
}


QColor gsGeometryDialog::initColorValue(void)
{
    return d->data->color();
}


void gsGeometryDialog::showBasis(void)
{
    gsBasisPointer myGismoData = 
	    getGeometryPointer(d->data)->basis().clone().release();
	std::cout << "Plotting basis "<< * myGismoData <<"\n";
    
	// Create basis object and add it in the object list
    gsBasisData * myData = new gsBasisData(myGismoData);
    // myData->setColor(QColor("#0080ff"));
    // double opacity = 1.0 - 0.01 * d->sliderOpacity->value();
    // myData->setOpacity(opacity);
    
	emit dataInserted(myData);
}

void gsGeometryDialog::saveGeometry(void)
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File to Gismo XML format"),"~/",tr("Gismo files (*.xml);;"));

    if ( fileName.size() )
    {
        gsGeometryPointer myGismoData = getGeometryPointer(d->data);
	if ( myGismoData )
	    {
		std::cout << "Saving to "<< fileName.toUtf8().constData() <<"\n";
                gismo::gsFileData<double> out;
                out<< *myGismoData ;
                out.dump(fileName.toUtf8().constData());
            }
    }	
}

void gsGeometryDialog::insertKnot(void)
{
    gsGeometryPointer myGismoData = getGeometryPointer(d->data);
    gsAxelSurface * surf = dynamic_cast<gsAxelSurface *>(d->data);
    
    if ( ! surf )
    {
        gsWarn<<"Not a surface.\n";
        return;
    }

    // Get index of the selected control point
    const int parameter = surf->getParameter();
    
    // knot coordinates
    real_t u, v;

    if ( gismo::gsTHBSpline<2> * hb = dynamic_cast<gismo::gsTHBSpline<2>*>(myGismoData.get()) )
    {
      // get level of basis function
      const int lvl = hb->basis().levelOf(parameter);

      gsDebugVar(lvl);

      const gismo::gsTensorBSplineBasis<2,real_t> & 
      tb = hb->basis().tensorLevel(lvl);

      const unsigned ind = hb->basis().flatTensorIndexOf(parameter);
      const gismo::gsVector<unsigned,2> tind = tb.tensorIndex(ind);

      const int degu = tb.degree(0);
      const int degv = tb.degree(1);

      const unsigned us = tind[0] + (degu + 1) / 2  ;
      const unsigned vs = tind[1] + (degv + 1) / 2  ;
      
      gsDebugVar(us);
      gsDebugVar(vs);

      u = tb.knot(0,us);
      v = tb.knot(1,vs);

      gsDebugVar(u);
      gsDebugVar(v);

      // Inserting knots u and v ..
      //hb->increaseMultiplicity( lvl,0,u,1);// knot, direction, incrAmount

      hb->basis().increaseMultiplicity(1,0,0.5,1);


      gsDebugVar( tb.knots(0).detail() );

    }
}


/* Joker
void gsGeometryDialog::jokerPlay(void)
{
    gsGeometryPointer myGismoData = getGeometryPointer(d->data);
	if ( myGismoData )
    {
        gsInfo<<"Joker!\n";

        emit update();
    }
}
//*/

void gsGeometryDialog::onEditCP(void)
{
    QDialogEditCP cpEdit(this);

    const int ret = cpEdit.exec();
    if ( ret == QDialog::Accepted)
    {
        gismo::gsGeometry<> * g = getGeometryPointer(d->data).get();
        const int cp  = d->cpIndex->value();
        const int gdim = g->geoDim();
        g->coef(cp,0) = cpEdit.xCoord();
        if ( gdim > 1)
        {
            g->coef(cp,1) = cpEdit.yCoord();
                if ( gdim > 2)
                    g->coef(cp,2) = cpEdit.zCoord();
        }

        d->data->touchStructure();// should be just onControlPointChange
        d->data->touchGeometry();
        
        emit update();

        //CALL_DATA_METHOD( samplingChanged );
    }
}

void gsGeometryDialog::refineGeometry(void)
{
    gsGeometryPointer myGismoData = getGeometryPointer(d->data);
    gsAxelSurface * surf = dynamic_cast<gsAxelSurface *>(d->data);
    
    if ( ! surf )
    {
        gsWarn<<"Not a surface.\n";
        return;
    }

    // Get index of the selected control point
    const int parameter = surf->getParameter(); // to do: replace by signal

    if ( gismo::gsTHBSpline<2> * hb = dynamic_cast<gismo::gsTHBSpline<2>*>(myGismoData.get()) )
    {
         gismo::gsMatrix<unsigned, 2, 2> elements;
         hb->basis().elementSupport_into(parameter, elements);
         gsInfo<<"element support: \n"<<  elements <<"\n";

         std::vector<unsigned> box;
         const int lvl = hb->basis().levelOf(parameter) + 1;// increase level by 1
         gsInfo<<"level: \n"<< lvl-1  <<"\n";
         box.push_back(lvl);
         box.insert(box.end(), elements.data(), elements.data()+4 );
         box[1]= elements(0,0) << 1; // Take indices to the next level
         box[2]= elements(1,0) << 1;
         box[3]= elements(0,1) << 1;
         box[4]= elements(1,1) << 1;

         gsInfo<<"Refining control point "<<parameter<<".\n";
         gsInfo<<"with support: "<<  box[1]<<", "<<box[2]<<", "<<box[3]<<", "<<box[4] <<"\n";

         hb->basis().refineElements_withCoefs(hb->coefs(), box);

         updateText();

         surf->updateControlGrid();

         emit update();
    }
    else
    {
        gsWarn<<"We cannot refine nothing else than a hier. spline yet.\n";
    }


    // m_process: member of the gsGeometryData
    //axlFieldSpatialCoordinatesCreator * m_process = new axlFieldSpatialCoordinatesCreator();
    // somewhere to delete it

    //process->setInput(d->data, 0); //+ PDE solution
    //process->update();
    
    // not needed, in general we add it to the 
    //emit dataInserted(process->output);
}
