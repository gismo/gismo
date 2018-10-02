/** @file gsGeometryData.hpp

    @brief This file Provides implemetation of G+Smo geometry data for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <dtkCoreSupport/dtkAbstractDataFactory.h>

#include <typeinfo>

// Gismo includes
#include <gsCore/gsGeometry.h>
//#include <gsThbs/gsHTensorBasis.h>
#include <gsHSplines/gsHTensorBasis.h>
//#include <gsNurbs/gsTensorBSpline.h>
//#include <gsNurbs/gsTensorNurbs.h>
#include <gsModeling/gsTrimSurface.h>


// /////////////////////////////////////////////////////////////////
// gsGeometryData
// /////////////////////////////////////////////////////////////////

template <class axlObj>
gsGeometryData<axlObj>::gsGeometryData(void) 
    : axlObj(), m_geometry(NULL), 
      sampling_u(DEFAULT_SAMPLES), sampling_v(DEFAULT_SAMPLES), sampling_w(DEFAULT_SAMPLES)
{

}

template <class axlObj>
gsGeometryData<axlObj>::gsGeometryData(gsGeometryPointer geo)
    : axlObj(), m_geometry(geo), 
      sampling_u(DEFAULT_SAMPLES), sampling_v(DEFAULT_SAMPLES), sampling_w(DEFAULT_SAMPLES)
{
    updateControlGrid();
}

/*
template <class axlObj>
gsGeometryData<axlObj>::~gsGeometryData(void)
{

}
*/

template <class axlObj>
bool gsGeometryData<axlObj>::registered(void)
{
    // is dynamic cast better ?
    if ( typeid(axlObj) == typeid(axlAbstractCurveBSpline) )
     	 return gsAxelPlugin::dataFactSingleton->registerDataType("SplineCurve", creategsGeometryData1);
    if ( typeid(axlObj) == typeid(axlAbstractSurfaceBSpline) )
	 return gsAxelPlugin::dataFactSingleton->registerDataType("SplineSurface", creategsGeometryData2);
    if ( typeid(axlObj) == typeid(axlAbstractVolumeBSpline) )
    	 return gsAxelPlugin::dataFactSingleton->registerDataType("SplineVolume", creategsGeometryData3);
    if ( typeid(axlObj) == typeid(axlAbstractSurfaceTrimmed) )
    	 return gsAxelPlugin::dataFactSingleton->registerDataType("TrimSurface", creategsTrimSurf);

    return false;
}

template <class axlObj>
QString gsGeometryData<axlObj>::description(void) const
{
    QString result = "A Gismo ";
    result.append(dataId<axlObj>::id);
    return result;
}

template <class axlObj>
QString gsGeometryData<axlObj>::identifier(void) const
{
    return dataId<axlObj>::id;
    //return "gsGeometryData";
}

template <class axlObj>
QString gsGeometryData<axlObj>::objectType(void) const
{
    return "axlAbstractData";
}

template <class axlObj>
int gsGeometryData<axlObj>::numCoefs(void) const
{
    return m_geometry->coefsSize();
}

template <class axlObj>
int gsGeometryData<axlObj>::countControlPoints(void) const
{
    return m_geometry->coefsSize();
}

template <class axlObj>
int gsGeometryData<axlObj>::countControlPoints_u(void) const
{
    return m_geometry->basis().component(0).size();
}
template <class axlObj>
int gsGeometryData<axlObj>::countControlPoints_v(void) const
{
    return m_geometry->basis().component(1).size();
}
template <class axlObj>
int gsGeometryData<axlObj>::countControlPoints_w(void) const
{
    return m_geometry->basis().component(2).size();
}

template <class axlObj>
int gsGeometryData<axlObj>::order(void) const
{
    return m_geometry->basis().degree(0)+1;
}
template <class axlObj>
int gsGeometryData<axlObj>::order_u(void) const
{
    return m_geometry->basis().degree(0)+1;
}
template <class axlObj>
int gsGeometryData<axlObj>::order_v(void) const
{
    return m_geometry->basis().degree(1)+1;
}
template <class axlObj>
int gsGeometryData<axlObj>::order_w(void) const
{
    return m_geometry->basis().degree(2)+1;
}

template <class axlObj>
void gsGeometryData<axlObj>::eval(axlPoint *point, double u)
{
    gismo::gsMatrix<double> ev;
    gismo::gsMatrix<double> uv(1,1);
    uv<< u ; 
    m_geometry->eval_into(uv,ev);
    point->setCoordinates(ev(0,0), ev(1,0), ( ev.rows()>2 ? ev(2,0) : 0 ) );
} 

template <class axlObj>
void gsGeometryData<axlObj>::eval(axlPoint *point, double u,double v)
{
    gismo::gsMatrix<double> ev;
    gismo::gsMatrix<double> uv(2,1);
    uv<< u, v ; 
    m_geometry->eval_into(uv,ev);
    point->setCoordinates(ev(0,0), ev(1,0), ( ev.rows()>2 ? ev(2,0) : 0 ) );
} 

template <class axlObj>
void gsGeometryData<axlObj>::eval(axlPoint *point, double u, double v, double w)
{
    gismo::gsMatrix<double> ev;
    gismo::gsMatrix<double> uv(3,1);
    uv<< u, v, w ; 
    m_geometry->eval_into(uv,ev);
    point->setCoordinates(ev(0,0), ev(1,0), ev(2,0) );
} 

template <class axlObj>
double gsGeometryData<axlObj>::startParam(void)
{
    return m_geometry->parameterRange()(0,0);
}

template <class axlObj>
double gsGeometryData<axlObj>::startParam_u(void)
{
    return m_geometry->parameterRange()(0,0);
}

template <class axlObj>
double gsGeometryData<axlObj>::startParam_v(void)
{
    return m_geometry->parameterRange()(1,0);
}

template <class axlObj>
double gsGeometryData<axlObj>::startParam_w(void)
{
    return m_geometry->parameterRange()(2,0);
}

template <class axlObj>
double gsGeometryData<axlObj>::endParam(void)
{
    return m_geometry->parameterRange()(0,1);
}

template <class axlObj>
double gsGeometryData<axlObj>::endParam_u(void)
{
    return m_geometry->parameterRange()(0,1);
}

template <class axlObj>
double gsGeometryData<axlObj>::endParam_v(void)
{
    return m_geometry->parameterRange()(1,1);
}

template <class axlObj>
double gsGeometryData<axlObj>::endParam_w(void)
{
    return m_geometry->parameterRange()(2,1);
}

template <class axlObj>
int gsGeometryData<axlObj>::numSamples(void)
{
    return sampling_u;
}

template <class axlObj>
int gsGeometryData<axlObj>::numSamples_u(void)
{
    return sampling_u;
}

template <class axlObj>
int gsGeometryData<axlObj>::numSamples_v(void)
{
    return sampling_v;
}

template <class axlObj>
int gsGeometryData<axlObj>::numSamples_w(void)
{
    return sampling_w;
}

template <class axlObj>
void gsGeometryData<axlObj>::setNumSamples(int sampling)
{
    sampling_u=sampling;
    emit this->samplingChanged();
    this->touchGeometry();
}
template <class axlObj>
void gsGeometryData<axlObj>::setNumSamples_u(int sampling)
{
    sampling_u=sampling;
    emit this->samplingChanged();
    this->touchGeometry();
}
template <class axlObj>
void gsGeometryData<axlObj>::setNumSamples_v(int sampling)
{
    sampling_v=sampling;
    emit this->samplingChanged();
    this->touchGeometry();
}
template <class axlObj>
void gsGeometryData<axlObj>::setNumSamples_w(int sampling)
{
    sampling_w=sampling;   
    emit this->samplingChanged();
    this->touchGeometry();
}

// template <class axlObj>
// int gsGeometryData<axlObj>::knotVectorSize_w(void) const
// {
//     return 0;
//     //    return m_geometry->basis().knots().size();
// }

template <class axlObj>
double gsGeometryData<axlObj>::getCoord(int n, int v) const 
{
    return m_geometry->coef( n-1, v);
}
template <class axlObj>
double gsGeometryData<axlObj>::getCoord(int n, int m, int v) const 
{
    return m_geometry->coef( n-1+(m-1)*this->countControlPoints_u()
			    , v);
}
template <class axlObj>
double gsGeometryData<axlObj>::getCoord(int n, int m, int k, int v) const 
{
    return m_geometry->coef( n-1+(m-1)*this->countControlPoints_u()
                             +(k-1)*this->countControlPoints_v()*this->countControlPoints_u()
			    , v);
}

template <class axlObj>
axlPoint gsGeometryData<axlObj>::getCoef(int n) const
{
    //std::cout<< "! getCoeff"<< n <<"\n";
    // Zero-based numbering
    const int pos = n-1;
    return axlPoint( m_geometry->coef(pos,0), m_geometry->coef(pos,1), 
                   ( m_geometry->geoDim()>2 ? m_geometry->coef(pos,2) : 0 ) );

}

template <class axlObj>
axlPoint gsGeometryData<axlObj>::getCoef(int n, int m) const 
{
    gsInfo <<"getCoeff"<< n <<" "<< m <<"\n";
    int pos = (n-1+(m-1)*this->countControlPoints_u()) ;
    return axlPoint( m_geometry->coef(pos,0), m_geometry->coef(pos,1), 
		    (m_geometry->geoDim()>2 ? m_geometry->coef(pos,2) : 0 ) );
}
template <class axlObj>
axlPoint gsGeometryData<axlObj>::getCoef(int n, int m, int k) const 
{
    int pos = (n-1+(m-1)*this->countControlPoints_u()+
              (k-1)*this->countControlPoints_v()*this->countControlPoints_u() ) ;
    return axlPoint( m_geometry->coef(pos,0), m_geometry->coef(pos,1), 
		    (m_geometry->geoDim()>2 ? m_geometry->coef(pos,2) : 0 ) );
}

template <class axlObj>
bool gsGeometryData<axlObj>::setCoef(int n, double *controlPoint)
{
    // Zero-based numbering
    n--;
    m_geometry->coef(n,0) = controlPoint[0];
    m_geometry->coef(n,1) = controlPoint[1];
    if (m_geometry->geoDim()>2)
        m_geometry->coef(n,2) = controlPoint[2];
    this->touchGeometry();
    return true;
} 
template <class axlObj>
bool gsGeometryData<axlObj>::setCoef(int i, int j, double *controlPoint)
{
    int pos = (i+j*this->countControlPoints_u()) ;
    m_geometry->coef(pos,0) = controlPoint[0];
    m_geometry->coef(pos,1) = controlPoint[1];
    if (m_geometry->geoDim()>2)
        m_geometry->coef(pos,2) = controlPoint[2];
    this->touchGeometry();
    return true;
} 
template <class axlObj>
bool gsGeometryData<axlObj>::setCoef(int i, int j, int k, double *controlPoint)
{
    int pos = i+j*this->countControlPoints_u()+
        k*this->countControlPoints_v()*this->countControlPoints_u() ;
    m_geometry->coef(pos,0) = controlPoint[0];
    m_geometry->coef(pos,1) = controlPoint[1];
    if (m_geometry->geoDim()>2)
	m_geometry->coef(pos,2) = controlPoint[2];
    this->touchGeometry();
    return true;
} 

template <class axlObj>
bool gsGeometryData<axlObj>::setCoef(int i, int j, int v, double c)
{
    m_geometry->coef(i+j*this->countControlPoints_u(),v) = c;
    this->touchGeometry();
    return true;
} 
template <class axlObj>
bool gsGeometryData<axlObj>::setCoef(int i, int j, int k, int v, double c)
{
    m_geometry->coef(i+j*this->countControlPoints_u()+
        k*this->countControlPoints_v()*this->countControlPoints_u(),v) = c;
    this->touchGeometry();
    return true;
} 


template <class axlObj>
void gsGeometryData<axlObj>::normal(axlPoint *normal, double u,double v)
{
    gismo::gsMatrix<double> tmp;
    gismo::gsMatrix<double> uv(2,1);
    uv<< u, v ; 
    // m_geometry->deriv_into(uv,ev);// this uses gsFunction as default..
    
    // to do: add on every gsGeom.
    // measure, measure_into, jac_into and maybe gradient(i)
    // deriv_into
    int n = m_geometry->geoDim();
    if ( n==2)
	{
	    //tmp = m_geometry->jacobian(uv);
 	    //gismo::gsMatrix<double> ev = tmp->topRows(2).cross(tmp->bottomRows(2) );
 	    //normal->setValues(ev(0,0), ev(1,0),  1 );
	    normal->setCoordinates(0, 0 , 1 );
	}
    else
	{
	    tmp = m_geometry->jacobian(uv);
	    //std::cout<< "J=\n"<< *tmp <<"\n";
	    gismo::gsMatrix<double> ev = 
		tmp.block<3,1>(0,0).cross( tmp.block<3,1>(0,1) );
	    // 	    gismo::gsMatrix<double> ev = tmp->topRows(2)* tmp->bottomRows(2);
	    normal->setCoordinates(ev(0,0), ev(1,0), ev(2,0) );
	}

}


template <class axlObj>
void gsGeometryData<axlObj>::setParameter(int parameter)
{
    m_parameter = parameter;
    gsInfo<< "Selected control point:"<< m_parameter <<"\n";
    //gsInfo<< "Coordinates           :"<< m_geometry->coef(m_parameter) <<"\n";
    
    /*
    const gismo::gsHTensorBasis<2,double> * hb = 
        dynamic_cast<const gismo::gsHTensorBasis<2,double>*>(&m_geometry->basis());
    if ( hb )
    {
        const int lvl =  hb->levelOf(m_parameter);
        gsInfo<< "Level                 :"<<  lvl <<"\n";
        const int ti = hb->flatTensorIndexOf(m_parameter);
        gsInfo<< "Flat tensor index          :"<<ti <<"\n";
        gsInfo<< "Tensor index          :"<< hb->tensorLevel(lvl).tensorIndex(ti).transpose() <<"\n";
    }
    //*/

}


template <class axlObj>
void gsGeometryData<axlObj>::updateControlGrid()
{
    if( dynamic_cast<gismo::gsHTensorBasis<2,double> *>( &m_geometry->basis() ) ||
        dynamic_cast<gismo::gsHTensorBasis<3,double> *>( &m_geometry->basis() ) ||
        dynamic_cast<gismo::gsHTensorBasis<1,double> *>( &m_geometry->basis() ) 
        )
    {
        // Remove previous control grid
        this->resetControlPointConnections();

        // Construct the control mesh
        gismo::gsMesh<double> msh;
        m_geometry->controlNet(msh);

        // Pass it to axel
        for (typename std::vector< gismo::gsEdge<double> >::const_iterator 
                 it=msh.edges().begin(); it!=msh.edges().end(); ++it)
        {    
            this->defineControlPointConnection( it->source->getId(), it->target->getId() );
        }

        // Update the structure in the vtkView
        this->touchStructure();
        
        // Note: we do not have access to the actor of the object. a
        // line like the one below must be done though the interface
        // instead, by connecting a slot of axlAbstractActor with a
        // signal of axlActorBSpline (and derived classes)
        // this->getActor()->onModifiedControlPoints();

        // this->touchGeometry();
        // this->touchGeometry();

        // Q: what does updated() do ?
        this->updated();

    }//end if
}



#undef DEFAULT_SAMPLES
