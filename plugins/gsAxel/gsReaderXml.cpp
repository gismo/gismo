/** @file gsReaderXml.cpp

    @brief This file provides implementation of the XML reader for Axel

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include "gsReaderXml.h"
#include "gsGeometryData.h"
#include "gsMultiPatchData.h"
#include "gsBasisData.h"

#include "gsAxelPlugin.h"

#include <dtkCoreSupport/dtkAbstractData.h>
#include <dtkCoreSupport/dtkAbstractDataFactory.h>

// G+SMO include
#include <gismo.h>

// /////////////////////////////////////////////////////////////////
// gsReaderXml
// /////////////////////////////////////////////////////////////////

gsReaderXml::gsReaderXml(void)
{
    this->setObjectName(this->identifier());   
}

gsReaderXml::~gsReaderXml(void)
{
    
}

QString gsReaderXml::identifier(void) const
{
    return "gsReaderXml";
}

QString gsReaderXml::description(void) const
{
    return "gsReaderXml";
}

QStringList gsReaderXml::handled(void) const
{
    return QStringList() << ".xml" << "G+Smo XML file (.xml)";
}

bool gsReaderXml::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataReaderType("gsReaderXml", QStringList(), creategsReaderXml);
}


bool gsReaderXml::canRead(const QString& file)
{
    std::string fn(file.toUtf8().constData() );
    // Identify XML extension
    std::string ext;
    ext = fn.substr(fn.find_last_of(".") + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    return (ext=="xml"); // || ext=="axl"
}


bool gsReaderXml::read(const QString& file)
{
    if ( file.size() )
    {
        gismo::gsFileData<double> filedata(file.toUtf8().constData());
      
        if ( filedata.has< gismo::gsGeometry<double> >() )
        {
            gsGeometryPointer myGismoData( filedata.getFirst<gismo::gsGeometry<double> >().release() );
            if ( myGismoData )
            {
                gsInfo << "Opened "<< * myGismoData <<"\n";
                // Create Axel object
                axlAbstractData * myData = newGeometryData(myGismoData);
                this->setData(myData);
                return true;
            }
        }

        if ( filedata.has< gismo::gsMultiPatch<double> >() )
	    {
            gsMultiPatchPointer myGismoData = 
                filedata.getFirst<gismo::gsMultiPatch<double> >().release();
            if ( myGismoData )
            {
                gsInfo << "Opened "<< * myGismoData <<"\n";
                // Create Axel object
                axlAbstractData * myData = new gsMultiPatchData(myGismoData);
                this->setData(myData);
                return true;
            }
	    }

        if ( filedata.has< gismo::gsBasis<double> >() )
        {
            gsBasisPointer myGismoData = 
                filedata.getFirst<gismo::gsBasis<double> >().release();
            if ( myGismoData )
            {
                gsInfo << "Opened "<< * myGismoData <<"\n";
                // Create Axel object
                axlAbstractData * myData = new gsBasisData(myGismoData);
                this->setData(myData);
                return true;
            }

        }
    }
    return false;
    // On linux
    //std::cout << file.toUtf8().constData() <<"\n";
    // or this if you on Windows 
    //std::cout << file.toLocal8Bit().constData() <<"\n";
}

dtkAbstractDataReader *creategsReaderXml(void)
{
    return new gsReaderXml;
}

