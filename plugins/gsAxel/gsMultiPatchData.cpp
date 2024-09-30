/** @file gsMultiPatchData.cpp

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
#include "gsAxelPluginExport.h"

#include <dtkCoreSupport/dtkAbstractDataFactory.h>

#include "gsMultiPatchData.h"
#include "gsGeometryData.h"

#include <gsCore/gsMultiPatch.h>

#include <QColor>

#define DEFAULT_SAMPLES 20

 
gsMultiPatchPointer getMultiPatchPointer( axlAbstractDataComposite * axlData)
{
    // To do: probably static_cast is suitable here
    if(gsMultiPatchData *gsData = dynamic_cast<gsMultiPatchData *>(axlData))
        return gsData->getGismoPointer() ;
    else
	std::cout <<"Problem, axlAbstractData does not downcast to gismo implementation.\n";
    return NULL;
}


// /////////////////////////////////////////////////////////////////
// gsGeometryData
// /////////////////////////////////////////////////////////////////


gsMultiPatchData::gsMultiPatchData(void) : axlAbstractDataComposite()
{
 m_mpatch = NULL;
}


gsMultiPatchData::gsMultiPatchData(gsMultiPatchPointer mpatch) : axlAbstractDataComposite(), m_mpatch(mpatch)
{
//sampling.resize( ( m_mpatch->dim() > 2 ? m_basis->dim() : 2) );
sampling.resize(2) ;

sampling.setConstant(DEFAULT_SAMPLES);

for ( gismo::gsMultiPatch<double>::const_iterator it = m_mpatch->begin(); 
      it != m_mpatch->end(); ++it )
{
    
    axlAbstractData * myData = newGeometryData( gismo::memory::make_shared_not_owned(*it) );
    myData->setColor(QColor("#0080ff"));   
    this->add( myData );
}

 this->setColor(QColor("#0080ff"));   

}



gsMultiPatchData::~gsMultiPatchData(void)
{
    delete m_mpatch;
}


bool gsMultiPatchData::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataType("gsMultiPatchData", creategsMultiPatchData);
}


QString gsMultiPatchData::description(void) const
{
    QString result = "gsMultiPatchData : A Gismo multipatch object";
    return result;
}


QString gsMultiPatchData::identifier(void) const
{
    return "gsMultiPatchData";
}


QString gsMultiPatchData::objectType(void) const
{
    // return "axlAbstractData";
    return "axlAbstractDataComposite";
}

gismo::gsVector<unsigned> gsMultiPatchData::numSamples(void)
{
    return sampling;
}

int gsMultiPatchData::numSamples_u(void)
{
    return sampling[0];
}

int gsMultiPatchData::numSamples_v(void)
{
    return sampling[1];
}

void gsMultiPatchData::setNumSamples_u(int smpl)
{
    sampling[0]=smpl;
//emit this->samplingChanged();
//    this->touch();
}

void gsMultiPatchData::setNumSamples_v(int smpl)
{
    sampling[1]=smpl;
//    emit this->samplingChanged();
//    this->touch();
}

#undef DEFAULT_SAMPLES

// /////////////////////////////////////////////////////////////////
// Type instanciation
// /////////////////////////////////////////////////////////////////

dtkAbstractData *creategsMultiPatchData(void)
{
    return new gsMultiPatchData;
}
