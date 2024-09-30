/** @file gsBasisData.h

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/#pragma once

#include "gsAxelPlugin.h"

#include <axlCore/axlAbstractData.h>
//#include <axlCore/axlAbstractDataComposite.h> // TO DO

#include <axlCore/axlPoint.h>

#include <gsCore/gsBasis.h>
#include <gsCore/gsLinearAlgebra.h>

// Forward Declarations
namespace gismo{
template <class T>      class gsBasis; 
}
class gsBasisData;

typedef gismo::gsBasis<double> * gsBasisPointer;

dtkAbstractDataConverter *creategsBasisConverter(void);

class GSAXELPLUGIN_EXPORT gsBasisData : public axlAbstractData // TO DO: inherit from axlAbstractDataComposite
{
    //As long as we don't need any custom signals or slots and we
    //don't need to access meta-data provided by QObject's methods, we
    //can omit Q_OBJECT macro. In fact, we let axlObj manage all that.
    Q_OBJECT
    
public:
    /// Creates an uninitialized SplineSurface, which can only be assigned to
    /// or read(...) into.
    gsBasisData(void);

    gsBasisData(gsBasisPointer basis);
    
    /// Virtual destructor, enables safe inheritance.
    virtual ~gsBasisData(void);
    
    virtual QString description(void) const;
    virtual QString identifier(void) const;
    virtual QString objectType(void) const;
    
    static bool registered(void);

public:
    // numSamples
    gismo::gsVector<unsigned> numSamples(void);
    int numSamples_u(void);
    int numSamples_v(void);

    void setNumSamples_u(int);
    void setNumSamples_v(int);

public:
    gsBasisPointer & getGismoPointer() { return m_basis; };

private:
    gsBasisPointer m_basis;

    gismo::gsVector<unsigned> sampling;
};

dtkAbstractData *creategsBasisData(void);

