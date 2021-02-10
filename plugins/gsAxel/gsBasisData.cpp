/** @file gsBasisData.cpp

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
#include "gsAxelPluginExport.h"

#include <dtkCoreSupport/dtkAbstractDataFactory.h>

#include "gsBasisData.h"// templated generic data

#define DEFAULT_SAMPLES 100


gsBasisPointer getBasisPointer( axlAbstractData * axlData)
{
    // To do: probably static_cast is suitable here
    if(gsBasisData *gsData = dynamic_cast<gsBasisData *>(axlData))
        return gsData->getGismoPointer() ;
    else
	std::cout <<"Problem, axlAbstractData does not downcast to gismo implementation.\n";
    return NULL;
}


// /////////////////////////////////////////////////////////////////
// gsGeometryData
// /////////////////////////////////////////////////////////////////


gsBasisData::gsBasisData(void) : axlAbstractData()
{

}


gsBasisData::gsBasisData(gsBasisPointer basis) : axlAbstractData(), m_basis(basis)
{
sampling.resize( ( m_basis->dim() > 2 ? m_basis->dim() : 2) );
sampling.setConstant(DEFAULT_SAMPLES);
}



gsBasisData::~gsBasisData(void)
{
    delete m_basis;
}


bool gsBasisData::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataType("gsBasisData", creategsBasisData);
}


QString gsBasisData::description(void) const
{
    QString result = "gsBasisData : A Gismo basis object";
    return result;
}


QString gsBasisData::identifier(void) const
{
    return "gsBasisData";
}


QString gsBasisData::objectType(void) const
{
    return "axlAbstractData";
}

gismo::gsVector<unsigned> gsBasisData::numSamples(void)
{
    return sampling;
}

int gsBasisData::numSamples_u(void)
{
    return sampling[0];
}

int gsBasisData::numSamples_v(void)
{
    return sampling[1];
}

void gsBasisData::setNumSamples_u(int smpl)
{
    sampling[0]=smpl;
//emit this->samplingChanged();
//    this->touch();
}

void gsBasisData::setNumSamples_v(int smpl)
{
    sampling[1]=smpl;
//    emit this->samplingChanged();
//    this->touch();
}

#undef DEFAULT_SAMPLES

// /////////////////////////////////////////////////////////////////
// Type instanciation
// /////////////////////////////////////////////////////////////////

dtkAbstractData *creategsBasisData(void)
{
    return new gsBasisData;
}
