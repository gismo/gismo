/* (C) RICAM */

/* Put a short description of your plugin here */

/* Angelos.Mantzaflaris@oeaw.ac.at */

#include "gsReader.h"
#include "gsGeometryData.h"
#include "gsMultiPatchData.h"
#include "gsBasisData.h"

#include "gsAxelPlugin.h"

#include <dtkCore/dtkAbstractData.h>
#include <dtkCore/dtkAbstractDataFactory.h>

// G+SMO include
#include <gismo.h>

// /////////////////////////////////////////////////////////////////
// gsReader
// /////////////////////////////////////////////////////////////////

gsReader::gsReader(void)
{
    
}

gsReader::~gsReader(void)
{
    
}

QString gsReader::identifier(void) const
{
    return "gsReader";
}

QString gsReader::description(void) const
{
    return "gsReader";
}

QStringList gsReader::handled(void) const
{
    return QStringList() << "G+Smo XML file (.xml)";
}

bool gsReader::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataReaderType("gsReader", QStringList(), creategsReader);
}

bool gsReader::canRead(const QString& file)
{
  std::string fn(file.toUtf8().constData() );
  // Identify XML extension
  std::string ext;
  ext = fn.substr(fn.find_last_of(".") + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  return (ext== "xml");
}

bool gsReader::read(const QString& file)
{
  if ( file.size() )
    {
      gismo::gsFileData<double> filedata(file.toUtf8().constData());
      
      if ( filedata.has< gismo::gsGeometry<double> >() )
        {
	  gsGeometryPointer myGismoData = 
	    filedata.getFirst<gismo::gsGeometry<double> >();
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
		filedata.getFirst<gismo::gsMultiPatch<double> >();
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
	    filedata.getFirst<gismo::gsBasis<double> >();
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

dtkAbstractDataReader *creategsReader(void)
{
    return new gsReader;
}

