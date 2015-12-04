
#include "gsAxelPlugin.h"

#include "gsGeometryData.h"
#include "gsGeometryConverter.h"

#include <axlCore/axlAbstractData.h>

#include <axlCore/axlMesh.h>
#include <axlCore/axlPoint.h>

#include <dtkCoreSupport/dtkAbstractDataFactory.h>

#include <gsUtils/gsPointGrid.hpp>

class gsGeometryConverterPrivate
{
public:
    //    gsGeometryData<> *data;
    axlAbstractData *data;
};

gsGeometryConverter::gsGeometryConverter(void) : axlAbstractDataConverter(), d(new gsGeometryConverterPrivate)
{
    d->data = NULL;
}

gsGeometryConverter::~gsGeometryConverter(void)
{
    delete d;

    d = NULL;
}

QString gsGeometryConverter::description(void) const
{
    return "Converter from gsGeometryConverter to axlMesh";
}

QStringList gsGeometryConverter::fromTypes(void) const
{
    return QStringList() << "gsGeometryConverter" << "axlAbstractGeometry";
}

QString gsGeometryConverter::toType (void) const
{
    return "axlMesh";
}

bool gsGeometryConverter::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataConverterType("gsGeometryConverter", QStringList(), "axlMesh", createGsGeometryConverter);
}

void gsGeometryConverter::setData(dtkAbstractData *data)
{
    d->data = dynamic_cast<axlAbstractData *>(data);
    // if(gsAxelCurve *gsData = dynamic_cast<gsAxelCurve *>(data))
    //d->data = gsData;
	//if(gsAxelSurface *gsData = dynamic_cast<gsAxelSurface *>(data))
        //d->data = gsData;
	//if(gsAxelVolume *gsData = dynamic_cast<gsAxelVolume *>(data))
        //d->data = gsData;

}

axlMesh *gsGeometryConverter::toMesh(void)
{
  using namespace gismo;

    // Fetch pointer to gismo geometry
  gsGeometryPointer Geo = getGeometryPointer( d->data );
    std::cout << "gsGeometryConverter is converting "<< * Geo  <<"\n";

    int ndim= Geo->geoDim();    
    unsigned pdim= Geo->parDim();
    int npts = 1000; // to do: take this from gsGeometryData

    gsMatrix<double> ab = Geo->parameterRange() ;
    gsVector<double> a = ab.col(0);
    gsVector<double> b = ab.col(1);
    gsVector<unsigned> np = uniformSampleCount(a,b, npts );// to do, cast to int
    gsMatrix<double> pts = gsPointGrid(a,b,np) ;

    gsMatrix<double>  eval_geo = Geo->eval  ( pts ) ;//pts

    if ( 3 - ndim > 0 )
    {
        eval_geo.conservativeResize(3,eval_geo.cols() );
        eval_geo.bottomRows(3-ndim).setZero();
    }

    axlMesh *mesh = new axlMesh;

    // Note: newer version of axel use QVector
    QVector<axlPoint *> pointSet;
    //QList<axlPoint *> pointSet;

    // Vertices
    for(index_t i = 0; i < eval_geo.cols(); i++) {
      pointSet.append( new axlPoint( eval_geo(0,i),
				     eval_geo(1,i),
				     eval_geo(2,i) ) );
    }
    mesh->setVertices( pointSet );

    // Curve
    if ( pdim == 1 )
    {
      mesh->push_back_new_edge();
      for(unsigned i = 0; i < np[0]; i++) {
        mesh->edgePushBack(0,i);
    }
    return mesh;
    }// end dim 1


    // Surface
    if ( pdim == 2 )
    {
      int ind1, ind2;
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

      return mesh;
    }// end dim 2

    if ( pdim == 3 )
    {
    int ind1, ind2, ind3, ind4;

    for(unsigned k= 0; k < np[2] - 1; k++)
    {
        for(unsigned j = 0; j < np[1] - 1; j++)
        {
            for(unsigned i= 0; i <np[0] - 1; i++)
            {

                ind1 = k* np[0] *np[1] + j * np[0] + i;
                ind2 = ind1 + np[0];
                ind3 = ind1 + np[0]*np[1] ;
                ind4 = ind2 + np[0]*np[1] ;


                QVector<int> Triangle1;
                QVector<int> Triangle2;
                QVector<int> Triangle3;
                QVector<int> Triangle4;
                QVector<int> Triangle5;
                QVector<int> Triangle6;
                QVector<int> Triangle7;
                QVector<int> Triangle8;
                QVector<int> Triangle9;
                QVector<int> Triangle10;
                QVector<int> Triangle11;
                QVector<int> Triangle12;
                QVector<int> Triangle13;
                QVector<int> Triangle14;
                QVector<int> Triangle15;
                QVector<int> Triangle16;
                QVector<int> Triangle17;
                QVector<int> Triangle18;

                                Triangle1.push_back(ind1);
                                Triangle1.push_back(ind2 + 1);
                                Triangle1.push_back(ind3);

                                Triangle2.push_back(ind3);
                                Triangle2.push_back(ind1 + 1);
                                Triangle2.push_back(ind2 + 1);

                                Triangle3.push_back(ind1);
                                Triangle3.push_back(ind1 + 1);
                                Triangle3.push_back(ind3);

                Triangle4.push_back(ind1);
                Triangle4.push_back(ind1 + 1);
                Triangle4.push_back(ind2 + 1);

                                Triangle5.push_back(ind3);
                                Triangle5.push_back(ind1 + 1);
                                Triangle5.push_back(ind3 +1);

                Triangle6.push_back(ind1+1);
                Triangle6.push_back(ind3 + 1);
                Triangle6.push_back(ind2 + 1);

                                Triangle7.push_back(ind3);
                                Triangle7.push_back(ind2+1);
                                Triangle7.push_back(ind4 +1);

                Triangle8.push_back(ind3);
                Triangle8.push_back(ind3 + 1);
                Triangle8.push_back(ind4 + 1);

                Triangle9.push_back(ind2+1);
                Triangle9.push_back(ind3 + 1);
                Triangle9.push_back(ind4+1);

                                Triangle10.push_back(ind3);
                                Triangle10.push_back(ind3 + 1);
                                Triangle10.push_back(ind2 + 1);

                Triangle11.push_back(ind2);
                Triangle11.push_back(ind2+1);
                Triangle11.push_back(ind4);

                                Triangle12.push_back(ind2);
                                Triangle12.push_back(ind3);
                                Triangle12.push_back(ind2+1);

                Triangle13.push_back(ind3);
                Triangle13.push_back(ind4);
                Triangle13.push_back(ind2);


                Triangle14.push_back(ind2);
                Triangle14.push_back(ind1);
                Triangle14.push_back(ind2+1);

                Triangle15.push_back(ind1);
                Triangle15.push_back(ind3 );
                Triangle15.push_back(ind2 );


                                Triangle16.push_back(ind4);
                                Triangle16.push_back(ind3);
                                Triangle16.push_back(ind2+1);

                Triangle17.push_back(ind4);
                Triangle17.push_back(ind3 );
                Triangle17.push_back(ind4+1 );

                Triangle18.push_back(ind2+1);
                Triangle18.push_back(ind4+1);
                Triangle18.push_back(ind4);

                mesh->push_back_face(Triangle1);
                mesh->push_back_face(Triangle2);
                 mesh->push_back_face(Triangle3);
                mesh->push_back_face(Triangle4);
                                mesh->push_back_face(Triangle5);
                mesh->push_back_face(Triangle6);
                                mesh->push_back_face(Triangle7);
                mesh->push_back_face(Triangle8);
                mesh->push_back_face(Triangle9);
                                mesh->push_back_face(Triangle10);
                mesh->push_back_face(Triangle11);
                                mesh->push_back_face(Triangle12);
                mesh->push_back_face(Triangle13);
                mesh->push_back_face(Triangle14);
                mesh->push_back_face(Triangle15);
                                mesh->push_back_face(Triangle16);
                mesh->push_back_face(Triangle17);
                mesh->push_back_face(Triangle18);


            }
        }
    }
    return mesh;
    }// end dim 3


    dtkWarn() << "Trouble: could not convert gsGeometryData into an axlMesh";
  return NULL;
}

dtkAbstractDataConverter *createGsGeometryConverter(void)
{
    return new gsGeometryConverter;
}
