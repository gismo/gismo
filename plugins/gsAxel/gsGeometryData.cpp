/** @file gsGeometryData.cpp

    @brief This file Provides declaration of G+Smo geometry data for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include "gsAxelPluginExport.h"

#include "gsGeometryData.h"// templated geometry data


// Identifier for different data
const QString dataId<axlAbstractCurveBSpline>::id   = "SplineCurve";
const QString dataId<axlAbstractSurfaceBSpline>::id = "SplineSurface";
const QString dataId<axlAbstractSurfaceTrimmed>::id = "TrimSurface";
const QString dataId<axlAbstractVolumeBSpline>::id  = "SplineVolume";


gsGeometryPointer getGeometryPointer( axlAbstractData * axlData)
{
    if(gsAxelCurve *gsData = dynamic_cast<gsAxelCurve *>(axlData))
        return gsData->getGismoPointer() ;
    else if (gsAxelSurface *gsData = dynamic_cast<gsAxelSurface *>(axlData))
        return gsData->getGismoPointer() ;
    else if (gsAxelTrimSurf *gsData = dynamic_cast<gsAxelTrimSurf *>(axlData))
        return gsData->getGismoPointer() ;
    else if(gsAxelVolume *gsData = dynamic_cast<gsAxelVolume *>(axlData))
         return gsData->getGismoPointer() ;
    else
	std::cout <<"Problem, axlAbstractData does not downcast to gismo implementation.\n";
    return NULL;
}

axlAbstractData * newGeometryData ( gsGeometryPointer gsData )
{
    //std::cout << "newGeometryData of dimension = " << gsData->parDim() <<"\n";
    switch ( gsData->parDim() )
	{
	case 1:
	    return new gsAxelCurve(gsData);
	    break;
	case 2:
	    if ( dynamic_cast<gismo::gsTrimSurface<double>*>(gsData.get()) )
		return new gsAxelTrimSurf(gsData);
	    else
		return new gsAxelSurface(gsData);
	    break;
	case 3:
	    return new gsAxelVolume(gsData);
	    break;
	default:
	    std::cout << "newGeometryData: Problem, dimension = "
		      << gsData->parDim() <<"\n";
	    return NULL;
       };
}

// /////////////////////////////////////////////////////////////////
// Template & Type instanciation
// /////////////////////////////////////////////////////////////////

template class GSAXELPLUGIN_EXPORT gsGeometryData< axlAbstractCurveBSpline   > ;
template class GSAXELPLUGIN_EXPORT gsGeometryData< axlAbstractSurfaceBSpline > ;
template class GSAXELPLUGIN_EXPORT gsGeometryData< axlAbstractVolumeBSpline  > ;
template class GSAXELPLUGIN_EXPORT gsGeometryData< axlAbstractSurfaceTrimmed > ;

dtkAbstractData *creategsGeometryData1(void)
{
    return new gsAxelCurve;
}

dtkAbstractData *creategsGeometryData2(void)
{
    return new gsAxelSurface;
}

dtkAbstractData *creategsGeometryData3(void)
{
    return new gsAxelVolume;
}

dtkAbstractData *creategsTrimSurf(void)
{
    return new gsAxelTrimSurf;
}

