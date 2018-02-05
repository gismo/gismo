/** @file gsGeometryCreator.cpp

    @brief This file Provides implemetation of G+Smo geometry dialog for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <axlCore/axlAbstractProcess.h>
#include <dtkCoreSupport/dtkAbstractProcessFactory.h>
#include <dtkCoreSupport/dtkAbstractDataFactory.h>
#include <dtkLog/dtkLog.h>
#include <dtkGuiSupport/dtkColorButton.h>

#include "gsAxelPlugin.h"
#include "gsGeometryData.h"
#include "gsBasisData.h"
#include "gsMultiPatchData.h"
#include "gsGeometryCreator.h"

#include <QtGui>

#include <gsIO/gsReadFile.h>
#include <gsNurbs/gsNurbsCreator.h>

/* // gsSmoothPatches
#include "gsSmoothPatches/gsCompositeBasis.h"
#include "gsSmoothPatches/gsCompositeBSplineBasis.h"
#include "gsSmoothPatches/gsCompositeGeom.h"
//*/

class gsGeometryCreatorPrivate {
public:
    QSlider *sliderOpacity;

    dtkColorButton *colorButton;

    int geometry;
    int degree;
};

gsGeometryCreator::gsGeometryCreator(QWidget* parent) :
    axlInspectorToolInterface(parent),
    d(new gsGeometryCreatorPrivate) {
    
    //int qLineEditWidth = 200;
    
    //OPACITY//
    d->sliderOpacity = new QSlider(Qt::Horizontal, this);
    
    QHBoxLayout *layoutOpacity = new QHBoxLayout;
    layoutOpacity->addWidget(new QLabel("Opacity",this));
    layoutOpacity->addWidget(d->sliderOpacity);
    d->sliderOpacity->setMaximum(100);
    
    //COLOR//
    d->colorButton = new dtkColorButton(this);
    QHBoxLayout *layoutColorButton = new QHBoxLayout;
    layoutColorButton->addWidget(new QLabel("Color",this));
    layoutColorButton->addWidget(d->colorButton);
    d->colorButton->setColor(QColor("#0080ff"));
    
    //Geometry selector//
    QComboBox * select_geo = new QComboBox;
    //QListView * listView = new QListView(select_geo);
    select_geo->addItem("Disc (NURBS)");
    select_geo->addItem("Sphere (NURBS)");
    select_geo->addItem("Square (B-spline)");
    select_geo->addItem("Cube (B-spline)");
    select_geo->addItem("Helf-cube (B-spline)");
    select_geo->addItem("1/4 Annulus(NURBS)");
    select_geo->addItem("Circle (NURBS)");
    select_geo->addItem("B-spline 1/4-annulus");
    connect(select_geo, SIGNAL(currentIndexChanged(int)), this, SLOT(update_index(int)));
    d->geometry= 0; // default value is the first item

    QSpinBox * degree_box = new QSpinBox;
    degree_box->setMaximum(100);
    degree_box->setMinimum(0);
    degree_box->setValue(2);
    connect(degree_box, SIGNAL(valueChanged(int)), this, SLOT(update_degree(int)));
    d->degree= 2; // default value

    //RUN BUTTON//
    QPushButton* createButton = new QPushButton("Create", this);
    connect(createButton, SIGNAL(clicked()), this, SLOT(run()));

    //LOAD GEO BUTTON//
    QPushButton* loadGButton = new QPushButton("Load Geometry", this);
    connect(loadGButton, SIGNAL(clicked()), this, SLOT(loadGeometry()));

    //LOAD MULTIPATCH BUTTON//
    QPushButton* loadMButton = new QPushButton("Load multipatch", this);
    connect(loadMButton, SIGNAL(clicked()), this, SLOT(loadMultiPatch()));

    //LOAD BASIS BUTTON//
    QPushButton* loadBButton = new QPushButton("Load Basis", this);
    connect(loadBButton, SIGNAL(clicked()), this, SLOT(loadBasis()));
    
    //LAYOUTS//
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    mainLayout->addWidget(new QLabel("gsGeometryCreator"));
    mainLayout->addLayout(layoutOpacity);
    mainLayout->addLayout(layoutColorButton);
    mainLayout->addWidget(select_geo);
    mainLayout->addWidget(degree_box);
    mainLayout->addWidget(createButton);

    mainLayout->addWidget(loadGButton);
    mainLayout->addWidget(loadBButton);
    mainLayout->addWidget(loadMButton);
}

gsGeometryCreator::~gsGeometryCreator() 
{
    delete d;
    d = NULL;
}

bool gsGeometryCreator::registered(void) 
{
    gsAxelPlugin::processFactSingleton->registerProcessType("gsGeometryCreator", createProcessCreatorgsGeometryData, "gismo"); // more than one ?
    return
	axlInspectorToolFactory::instance()->registerInspectorTool("gsGeometryCreator", creategsGeometryCreator);
}

void gsGeometryCreator::update_index(int i) 
{
  d->geometry= i;
  //std::cout << "gsGeometryCreator got index "<< i <<"\n";
}

void gsGeometryCreator::update_degree(int i) {
  d->degree= i;
  //std::cout << "gsGeometryCreator got index "<< i <<"\n";
}


void gsGeometryCreator::loadGeometry(void) {

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File with G+Smo"),"~/",tr("Gismo files (*.xml);;GeoPDE files (*.txt);;Axel files (*.axl);;All files (*.*)"));

    if ( fileName.size() )
    {
        gismo::gsGeometry<double>::Ptr myGismoData = gismo::gsReadFile<double>(fileName.toUtf8().constData() );
	if ( myGismoData )
	    {
		std::cout << "Loaded a "<< * myGismoData <<"\n";
		
		// Create Axel object
		axlAbstractData * myData = newGeometryData(myGismoData);	
		myData->setColor(QColor("#0080ff"));
		double opacity = 1.0 - 0.01 * d->sliderOpacity->value();
		myData->setOpacity(opacity);
		
		emit dataInserted(myData);
	    }
    }
}

void gsGeometryCreator::loadBasis(void) 
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File with G+Smo"),"~/",tr("Gismo files (*.xml);;GeoPDE files (*.txt);;Axel files (*.axl);;All files (*.*)"));

    if ( fileName.size() )
    {
        gismo::gsBasis<double>::uPtr myGismoData = gismo::gsReadFile<double>(fileName.toUtf8().constData() ) ;
        if ( myGismoData )
	    {
            
            std::cout << "Loaded a "<< * myGismoData <<"\n";
            
            // Create Axel object
            gsBasisData * myData = new gsBasisData(myGismoData.release());
		
            myData->setColor(QColor("#0080ff"));
            double opacity = 1.0 - 0.01 * d->sliderOpacity->value();
            myData->setOpacity(opacity);
            
            emit dataInserted(myData);
	    }
    }
}

void gsGeometryCreator::loadMultiPatch(void) 
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File with G+Smo"),"~/",tr("Gismo files (*.xml);;GeoPDE files (*.txt);;Axel files (*.axl);;All files (*.*)"));

/* // SmoothPatch checkbox

    QFileDialog fbox(this,
                     tr("Open File with G+Smo"),"~/",
                     tr("Gismo files (*.xml);;GeoPDE files (*.txt);;Axel files (*.axl);;All files (*.*)") );

    fbox.setOption(QFileDialog::DontUseNativeDialog);

    QCheckBox * chbox = new QCheckBox(&fbox);
    chbox->setText("Smoothen patches");

    // Add checkbox and show
    QHBoxLayout * hbl = new QHBoxLayout();
    hbl->addWidget(chbox);
    QGridLayout * layout = dynamic_cast<QGridLayout*>(fbox.layout());
    const int numRows = layout->rowCount();
    layout->addLayout( hbl, numRows,0,1,-1);//span all columns
    fbox.show();
    
    QString fileName;
    if (fbox.exec()) fileName = fbox.selectedFiles()[0];
//*/

    if ( fileName.size() )
    {
        gismo::gsMultiPatch<double>::uPtr myGismoData = gismo::gsReadFile<double>(fileName.toUtf8().constData() ) ;
        if (myGismoData)
	    {
            /* // SmoothPatches

            if ( chbox->isChecked() )
            {
                std::cout << "Loading a smooth "<< * myGismoData <<"\n";

                // Assume surface
                axlAbstractData * myData = new gsAxelSurface(
                    new gismo::gsCompositeGeom<2,real_t>(*myGismoData) );
                delete myGismoData;

                myData->setColor(QColor("#0080ff"));
                double opacity = 1.0 - 0.01 * d->sliderOpacity->value();
                myData->setOpacity(opacity);
                emit dataInserted(myData);
            }
            else
            //*/
            {
                std::cout << "Loading a "<< *myGismoData <<"\n";

                // Create Axel multipatch object
                gsMultiPatchData * myData = new gsMultiPatchData(myGismoData.release());

                myData->setColor(QColor("#0080ff"));
                double opacity = 1.0 - 0.01 * d->sliderOpacity->value();
                myData->setOpacity(opacity);
                emit dataInserted(myData);
            }
        }
        
    }
}

void gsGeometryCreator::run(void) 
{

    gismo::gsGeometry<double>::Ptr myGismoData;
   
    std::string filename;
    switch ( d->geometry )
    {
    case 1:
        filename = "sphere.xml"; //default example
        myGismoData = gismo::gsReadFile<double>(filename);
        //gismo::gsNurbsCreator<>::NurbsSphere();// some error in the coefs
        break;
    case 2:
        myGismoData = gismo::gsNurbsCreator<>::BSplineSquare(d->degree);
        break;
    case 3:
        myGismoData = gismo::gsNurbsCreator<>::BSplineCube(d->degree);
        break;
    case 4:
        myGismoData = gismo::gsNurbsCreator<>::BSplineHalfCube();
        myGismoData->coefs().array() += 2;
        break;
    case 5:
        myGismoData = gismo::gsNurbsCreator<>::NurbsQuarterAnnulus();
        myGismoData->coefs().array() -= 2.5;
        break;
    case 6:
        myGismoData = gismo::gsNurbsCreator<>::NurbsCircle();
        break;
    case 7:
        myGismoData = gismo::gsNurbsCreator<>::BSplineQuarterAnnulus(d->degree);
        break;
        
    default : // = 0
        myGismoData = gismo::gsNurbsCreator<>::NurbsDisk();
        break;
        
    };
    
    // Create Axel object
    axlAbstractData * myData = newGeometryData(myGismoData);
    std::cout << "gsGeometryCreator created a "<< * myGismoData <<"\n";
    
    myData->setColor(d->colorButton->color());
    double opacity = 1.0 - 0.01 * d->sliderOpacity->value();
    myData->setOpacity(opacity);
    
    emit dataInserted(myData);
}


dtkAbstractProcess* createProcessCreatorgsGeometryData(void) 
{
//dtkAbstractProcess* createProcessCreatorMyPluginDefaultData(void) {
    axlAbstractProcess* process = new axlAbstractProcess;
    process->setDescription("processCreator GismoGeometry create a new gsGeometryData");
    process->setIdentifier("gsGeometryCreator");

    return process; // to pass the factory
}

axlInspectorToolInterface* creategsGeometryCreator(void) 
{
    return new gsGeometryCreator;
}

