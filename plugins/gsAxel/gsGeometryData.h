/** @file gsGeometryData.h

    @brief This file Provides declaration of G+Smo geometry data for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include "gsAxelPlugin.h"

#include <axlCore/axlAbstractCurveBSpline.h>
#include <axlCore/axlAbstractSurfaceBSpline.h>
#include <axlCore/axlAbstractVolumeBSpline.h>
#include <axlCore/axlAbstractSurfaceTrimmed.h>

#include <axlCore/axlPoint.h>

#define DEFAULT_SAMPLES 20

// Forward Declarations
namespace gismo
{
    template <class T>      class gsGeometry; 
    template <class T>      class gsMesh; 
}

template <class axlObj> class gsGeometryData;

typedef gismo::gsGeometry<double>::Ptr gsGeometryPointer;

// Short names for different objects
typedef gsGeometryData<axlAbstractCurveBSpline>   gsAxelCurve;
typedef gsGeometryData<axlAbstractSurfaceBSpline> gsAxelSurface;
typedef gsGeometryData<axlAbstractVolumeBSpline>  gsAxelVolume;
typedef gsGeometryData<axlAbstractSurfaceTrimmed> gsAxelTrimSurf;

// Identifier for different data
template <class axlObj>
struct dataId
{
    using axlObj::ERROR_NO_IMPLEMENTATION;
};

template <>
struct dataId<axlAbstractCurveBSpline>
{
    static const QString id;
};

template <>
struct dataId<axlAbstractSurfaceBSpline>
{
    static const QString id;
};

template <>
struct dataId<axlAbstractSurfaceTrimmed>
{
    static const QString id;
};

template <>
struct dataId<axlAbstractVolumeBSpline>
{
    static const QString id;
};


// Creator functions (suffix is parametric dim of object)
dtkAbstractData *creategsGeometryData1(void);
dtkAbstractData *creategsGeometryData2(void);
dtkAbstractData *creategsGeometryData3(void);
dtkAbstractData *creategsTrimSurf(void);

// Helpers for commonly used routines
gsGeometryPointer getGeometryPointer( axlAbstractData * axlData);
axlAbstractData * newGeometryData   ( gsGeometryPointer gsData );


class gsGeometryBase
{
public:
    virtual axlAbstractData * toAxel() = 0;
};

template <class axlObj>
class gsGeometryData : public axlObj, public gsGeometryBase
{
    //As long as we don't need any custom signals or slots and we
    //don't need to access meta-data provided by QObject's methods, we
    //can omit Q_OBJECT macro. In fact, we let axlObj manage all that.
    //Q_OBJECT
    
public:
    typedef axlObj Base;

public:
    /// Creates an uninitialized gsGeometryData, which can only be
    /// assigned to some gsGeometry later.
    gsGeometryData(void);

    /// Creates axel geometry data and ties it to a gismo pointer
    gsGeometryData(gsGeometryPointer geo);
    
    /// Virtual destructor, enables safe inheritance.
    virtual ~gsGeometryData(void) { }
    
    axlObj * toAxel() {return static_cast<axlObj*>(this); }

    virtual QString description(void) const;
    virtual QString identifier(void) const;
    virtual QString objectType(void) const;
    
    static bool registered(void);

public:
    // countControlPoints
    int numCoefs() const;
    int countControlPoints(void)   const;
    int countControlPoints_u(void) const;
    int countControlPoints_v(void) const;
    int countControlPoints_w(void) const;

    // order
    int order  (void) const;
    int order_u(void) const;
    int order_v(void) const;
    int order_w(void) const;

    // eval
    void eval(axlPoint *point, double u);
    void eval(axlPoint *point, double u, double v);   
    void eval(axlPoint *point, double u, double v, double w);   

    // startParam
    double startParam(void);
    double startParam_u(void);
    double startParam_v(void);
    double startParam_w(void);

    // endParam
    double endParam(void);
    double endParam_u(void);
    double endParam_v(void);
    double endParam_w(void);

    // numSamples
    int numSamples(void);
    int numSamples_u(void);
    int numSamples_v(void);
    int numSamples_w(void);
    void setNumSamples  (int sampling);
    void setNumSamples_u(int sampling);
    void setNumSamples_v(int sampling);
    void setNumSamples_w(int sampling);

    // Get knot vector sizes
    // int knotVectorSize  (void) const;
    // int knotVectorSize_u(void) const;
    // int knotVectorSize_v(void) const;
    //int knotVectorSize_w(void) const;

    // Get coefficient
    double   getCoord(int n, int v)  const;// curve
    double   getCoord(int n, int m, int v) const;// surface
    double   getCoord(int n, int m, int k, int v) const;// volume

    axlPoint getCoef(int n) const  ;// curve
    axlPoint getCoef(int n, int m) const;// surface
    axlPoint getCoef(int n, int m, int k) const;// volume

    // Set coefficient
    bool setCoef(int n, double *controlPoint);
    bool setCoef(int i, int j, int k, double *controlPoint);
    bool setCoef(int i, int j, double *controlPoint);

    bool setCoef(int i, int j, int v, double c);    
    bool setCoef(int i, int j, int k, int v, double c);

    //bool rational(void) const;
    
    void normal(axlPoint *normal, double u,double v);

    //    bool setCoef(int n, double *controlPoint);

    // TrimSurf
    //void inDomain(axlPoint *normal, double u,double v);

public:

    //
    virtual void setParameter(int parameter);

    int getParameter() { return m_parameter;}

    void updateControlGrid();

public:
    gsGeometryPointer & getGismoPointer() { return m_geometry; };

private:
    gsGeometryPointer m_geometry;

    int sampling_u;
    int sampling_v;
    int sampling_w;

    int m_parameter;
};


#include "gsGeometryData.hpp"
