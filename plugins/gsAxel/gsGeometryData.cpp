
#include "gsAxelPluginExport.h"

#include "gsGeometryData.h"// templated geometry data



gsGeometryPointer getGeometryPointer( axlAbstractData * axlData)
{
    // To do: probably static_cast is suitable here
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
	    if ( dynamic_cast<gismo::gsTrimSurface<double>*>(gsData) )
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

