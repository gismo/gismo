/** @file gsXmlUtils.hpp

    @brief Provides implementation of input/output XML utilities struct.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <fstream>

// The concrete-type includes are gone: since the serialization
// inversion (S3) this file only defines the abstract-base dispatch,
// resolved through gsXmlRegistry at runtime. Modules own their
// concrete gsXml specializations (gs*Xml.hpp / class .hpp files).
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBoundary.h>

#include <gsIO/gsXmlGenericUtils.hpp>

#define GSXML_PUT_DYNAMIC_CAST(TYPE) \
if (const TYPE * g = dynamic_cast<const TYPE *>(ptr)) \
    return gsXml<TYPE>::put(*g, data);

#define GSXML_GET_TYPE(TYPE) \
if (s == gsXml<TYPE>::type().c_str()) \
    return gsXml<TYPE>::get(node);

namespace gismo {

namespace internal {

/*
 * Getting Xml data
 */

/// Get a solid
/// Get a Mesh
/// Get a Matrix from XML data
/// Get a SparseMatrix from XML data
/*
 * Getting Bases from XML data
 */


/// Get a TensorNurbsBasis from XML data
/// Get a Nurbs from XML data
/// Get a Tensor Nurbs from XML data
/// Get a TrimSurface from XML data
/// Get a Geometry from XML data
template<class T>
class gsXml< gsGeometry<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsGeometry<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return ""; }

    static gsGeometry<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Geometry") ),
                      "Something went wrong, was waiting for a Geometry tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "Geometry without a type in the xml file\n";
            return NULL;
        }
        // resolved through the runtime registry: each module registers
        // the concrete types it owns (see gs*XmlRegistration.h)
        return gsXmlDispatch< gsGeometry<T> >::get(node);
    }


    static gsXmlNode * put (const gsGeometry<T> & obj,
                            gsXmlTree & data)
	{
        // resolved through the runtime registry put-chain, which
        // reproduces the historical dynamic_cast order via priorities
        return gsXmlDispatch< gsGeometry<T> >::put(obj, data);
	}

};

/// Get a Function from XML data
template<class T>
class gsXml< gsFunction<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsFunction<T>);
    static std::string tag () { return "Function"; }
    static std::string type () { return ""; }

    static gsFunction<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Function") ),
                      "Something went wrong, was waiting for a Function tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "Function without a type in the xml file\n";
            return NULL;
        }
        // resolved through the runtime registry (see gs*XmlRegistration.h)
        return gsXmlDispatch< gsFunction<T> >::get(node);
    }

    static gsXmlNode * put (const gsFunction<T> & obj,
                            gsXmlTree & data)
	{
        // resolved through the runtime registry put-chain
        return gsXmlDispatch< gsFunction<T> >::put(obj, data);
	}

};

/// Get a Curve from XML data
template<class T>
class gsXml< gsCurve<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsCurve<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return ""; }

    static gsCurve<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Geometry") ),
                      "Something went wrong, was waiting for a Geometry tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "Geometry without a type in the xml file\n";
            return NULL;
        }
        // resolved through the runtime registry: each module registers
        // the concrete types it owns (see gs*XmlRegistration.h)
        return gsXmlDispatch< gsCurve<T> >::get(node);
    }


    static gsXmlNode * put (const gsCurve<T> & obj,
                            gsXmlTree & data)
	{
        // resolved through the runtime registry put-chain, which
        // reproduces the historical dynamic_cast order via priorities
        return gsXmlDispatch< gsCurve<T> >::put(obj, data);
	}

};

/// Get a Surface from XML data
template<class T>
class gsXml< gsSurface<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsSurface<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return ""; }

    static gsSurface<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Geometry") ),
                      "Something went wrong, was waiting for a Geometry tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "Geometry without a type in the xml file\n";
            return NULL;
        }
        // resolved through the runtime registry: each module registers
        // the concrete types it owns (see gs*XmlRegistration.h)
        return gsXmlDispatch< gsSurface<T> >::get(node);
    }

    static gsXmlNode * put (const gsSurface<T> & obj,
                            gsXmlTree & data)
	{
        // resolved through the runtime registry put-chain, which
        // reproduces the historical dynamic_cast order via priorities
        return gsXmlDispatch< gsSurface<T> >::put(obj, data);
	}
};

/// Get a Basis from XML data
template<class T>
class gsXml< gsBasis<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBasis<T>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return ""; }

    static gsBasis<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Basis") ),
                      "Something went wrong, waiting for a basis tag.\n" );

        gsXmlAttribute * gtype = node->first_attribute("type");
        if ( ! gtype )
        {
            gsWarn<< "Basis without a type in the xml file.\n";
            return NULL;
        }
        // resolved through the runtime registry: each module registers
        // the concrete types it owns (see gs*XmlRegistration.h)
        return gsXmlDispatch< gsBasis<T> >::get(node);
    }

    static gsXmlNode * put (const gsBasis<T> & obj,
                            gsXmlTree & data )
    {
        // resolved through the runtime registry put-chain, which
        // reproduces the historical dynamic_cast order via priorities
        return gsXmlDispatch< gsBasis<T> >::put(obj, data);
	}
};

/// Get a Pde from XML data
template<class T>
class gsXml< gsPde<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsPde<T>);
    static std::string tag () { return "Pde"; }
    static std::string type () { return ""; }

    static gsPde<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"Pde"),
                      "Something went wrong. Expected Pde tag." );

        // resolved through the runtime registry (see gsPdeXmlRegistration.h)
        return gsXmlDispatch< gsPde<T> >::get(node);
    }

    static gsXmlNode * put (const gsPde<T> &,
                            gsXmlTree & )
    {
        return NULL;
    }
};


// Get a Multipatch
/// Get a MultiBasis from XML data
/// Get a PlanarDomain from XML data
/// Get a CurveLoop from XML data
/// Get a Curve fitting class data from XML
/// Get a Poisson Pde from XML data
/*
/// Get a Poisson Pde from XML data
*/

/*
/// Get a Boundary Value Problem from XML
*/

}// end namespace internal

}// end namespace gismo

//#undef GSXML_COMMON_FUNCTIONS
//#undef TMPLA2
//#undef TMPLA3
#undef GSXML_GET_TYPE
#undef GSXML_PUT_DYNAMIC_CAST
