
#include "gsBasisData.h"
#include "gsBasisConverter.h"
#include "gsAxelPlugin.h"

#include <axlCore/axlAbstractData.h>

#include <axlCore/axlMesh.h>
#include <axlCore/axlPoint.h>

#include <dtkCoreSupport/dtkAbstractDataFactory.h>

#include <gsUtils/gsPointGrid.hpp>

class gsBasisConverterPrivate
{
public:
    gsBasisData * data;
};

gsBasisConverter::gsBasisConverter(void) : axlAbstractDataConverter(), d(new gsBasisConverterPrivate)
{
    d->data = NULL;
}

gsBasisConverter::~gsBasisConverter(void)
{
    delete d;

    d = NULL;
}

QString gsBasisConverter::description(void) const
{
    return "Converter from gsBasisData to axlMesh";
}

QStringList gsBasisConverter::fromTypes(void) const
{
    return QStringList() << "gsBasisData"; //"gsBasisData";
}

QString gsBasisConverter::toType (void) const
{
    return "axlMesh";
}

bool gsBasisConverter::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataConverterType("gsBasisConverter", QStringList(), "axlMesh", creategsBasisConverter);
}

void gsBasisConverter::setData(dtkAbstractData *data)
{
    if ( !( d->data = dynamic_cast<gsBasisData *>(data) ) )
        gsInfo <<" Problem: gsBasisConverter received the wrong object.\n";
}

axlMesh *gsBasisConverter::toMesh(void)
{
  using namespace gismo;
  if ( ! d->data ) 
  {
    gsWarn << "gsBasisConverter was called on a NULL pointer\n";
    return new axlMesh;
  }
    // Fetch pointer to gismo geometry
    gsBasisPointer basis = d->data->getGismoPointer() ;
    gsDebug << "gsBasisConverter is converting "<< * basis <<"\n";

    axlMesh *basismesh = new axlMesh;

    for ( int k = 0; k != basis->size(); ++k)
    {
	axlMesh *mesh = new axlMesh;

    // Note: Became QVector in recent Axel
    QVector<axlPoint *> pointSet;
    //QList<axlPoint *> pointSet;

	int pdim= basis->dim();
	
	gsMatrix<double> ab = basis->support() ;
	gsVector<double> a = ab.col(0);
	gsVector<double> b = ab.col(1);
	
	gsVector<unsigned> np = d->data->numSamples();
	gsMatrix<double> pts = gsPointGrid(a,b, d->data->numSamples() ) ;
	
	const gsMatrix<double> eval_geo = basis->evalSingle ( k, pts ) ;//pts
		
	if ( pdim > 2 )
    {
        gsWarn<<"Info: The dimension is to big, projecting into first 2 coordinatess..\n";
		pdim=2;
		pts.conservativeResize(2,eval_geo.cols());
    }
	
	// Vertices
	for ( index_t j=0; j<eval_geo.cols(); ++j)
    {
		pointSet.append( new axlPoint( pts(0,j),
                                       ( pdim>1 ? pts(1,j) : 0 ),
                                       eval_geo(0,j) ) );
    }
    
	mesh->setVertices( pointSet );

	// Curve
	if ( pdim == 1 )
	    {
		mesh->push_back_new_edge();
		for(unsigned i = 0; i < np[0]; i++) {
		    mesh->edgePushBack(0,i);
		}
	    }// end dim 1
	else   // Surface
	    if ( pdim == 2 )
		{
		    unsigned ind1, ind2;
		    for(unsigned j = 0; j < np[0] - 1; j++) {
			for(unsigned i= 0; i <np[1] - 1; i++) {
			    
			    ind1 =  j * np[0] + i;
			    ind2 = ind1 + np[1];
	  
			    QVector<int> firstTriangle;
			    QVector<int> secondTriamgle;
			    
			    firstTriangle.push_back(ind1);
			    firstTriangle.push_back(ind1 + 1);
			    firstTriangle.push_back(ind2);
			    
			    secondTriamgle.push_back(ind2);
			    secondTriamgle.push_back(ind1 + 1);
			    secondTriamgle.push_back(ind2 + 1);
	  
			    mesh->push_back_face(firstTriangle);
			    mesh->push_back_face(secondTriamgle);
			}
		    }
		    
		    // Add the loop of boundary edges
		    axlMesh::Edge e0,e1,e2,e3;
		    for(unsigned i = 0; i < np[0] ; i++) {
			e0<<i;
		    }
		    mesh->push_back_edge(e0);
		    
		    for(unsigned j = 0; j < np[1] ; j++) {
			e1<<(np[0]-1+j*np[0]);
		    }
		    mesh->push_back_edge(e1);
		    
		    for(unsigned i = 0; i < np[0] ; i++) {
			e2<<(np[1]-1)*np[0] + (np[0]-1-i);
		    }
		    mesh->push_back_edge(e2);
		    
		    for(unsigned j = 0; j < np[1] ; j++) {
			e3<<(np[1]-1-j)*np[0];
		    }
		    mesh->push_back_edge(e3);
		}// end dim 2
	
	basismesh->append(mesh);
    }

    return basismesh;
    dtkWarn() << "Trouble: could not convert gsBasisData into an axlMesh";
  return NULL;
}

dtkAbstractDataConverter *creategsBasisConverter(void)
{
    return new gsBasisConverter;
}
