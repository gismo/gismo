/** @file gsXmlUtils.hpp

    @brief Provides implementation of input/output XML utilities struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <fstream>

//#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsRationalBasis.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsBoundary.h>

#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorNurbs.h>

#include <gsHSplines/gsHBSpline.h>
#include <gsHSplines/gsTHBSpline.h>

#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsSolid.h>
#include <gsModeling/gsCurveFitting.h>

#include <gsUtils/gsMesh/gsMesh.h>

//#include <gsTrBezier/gsTriangularBezierBasis.h>
//#include <gsTrBezier/gsTriangularBezier.h>
#include <gsPde/gsPoissonPde.h>
#include <gsPde/gsSurfacePoissonPde.h>


#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo {

namespace internal {
 
/*
 * Getting Xml data
 */

/// Get a solid
template<class T>
class gsXml< gsSolid<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsSolid<T>);
    static std::string tag () { return "Solid"; }
    static std::string type () { return ""; }

    static gsSolid<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"Solid"), 
                      "Something went wrong. Expected Solid tag." );

        gsSolid<T> * m = new gsSolid<T>;

        int n  = atoi ( node->first_attribute("vertices")->value() ) ;
        int nVol  = atoi ( node->first_attribute("volumes")->value() ) ;
        T x,y, z;
        gsXmlNode * tmp = node->first_node("Vertex");
        std::istringstream str;
        str.str( tmp->value() );

        int nf;
        int vertID;
        int trimID;
        int ntest(0);
        std::vector< std::vector< gsSolidHeVertex<T>* > > vert;

        // get vertices
        for (int i=0; i<n; ++i)
        {
            ntest++;
            gsGetReal(str, x);
            gsGetReal(str, y);
            gsGetReal(str, z);
            m->addHeVertex(x,y,z);
        }
        GISMO_ASSERT( ntest==n, 
                      "Number of vertices does not match the Solid tag." );

        // get faces and surfaces
        gsXmlNode * toplevel = node->parent();// the geometry patches should be siblings of node
        n  = atoi ( node->first_attribute("faces")->value() ) ;
        tmp = node->first_node("Face");
        std::istringstream strf;
        strf.str( tmp->value() );
        ntest = 0;
        for (int iface=0; iface<n; iface++)
        {
            vert.clear();
            ntest++;
            size_t nLoops = 0;
            do
            {
                nLoops++;
                vert.push_back(std::vector< gsSolidHeVertex<T>* >());
                gsGetInt(strf, nf);
                // read num vertices on this loop
                for (int ivert=0; ivert<nf; ivert++) // read vertices
                {
                    gsGetInt(strf, vertID);
                    vert[nLoops-1].push_back(m->vertex[vertID]);
                }
                // next number is either:
                // * the trim surface id,
                // * -1 to indicate that the surface is automatically computed, or
                // * -2 to indicate that there are further internal loops
                gsGetInt(strf, trimID);
            } while(trimID <= -2); // -2 indicates that there are vertices remaining
            if (trimID>-1)
                m->addFace(vert, getById< gsTrimSurface<T> >( toplevel, trimID ) );
            else if (trimID==-1 && nLoops == 1)
                m->addFace(vert[0]);
            else if (trimID==0 && nLoops == 1)
                GISMO_ERROR("Faces must have unequal 0 as id (last value: increase from 1 or use -1 for all)");  // otherwise SEGFAULT happens
            else if (trimID==-1)
            {
                gsWarn<<"\nAutomatic creation of trimmed surfaces is only supported for a single loop\n";
            }
            else
            {
                gsWarn<<"\n ID of the trimmed surface trimID=" <<trimID<<" is invalid (must be >=-1)\n";
            }
        }
        m->setHeMate();
        // read in volumes. (optional)
        gsXmlNode * nodeVol = node->first_node("Volume");
        if(nodeVol == NULL)
        {
            GISMO_ASSERT(nVol == 1, "More than one volume but faces for volumes not specified");
            m->addVolume(m->face);
        }
        else
        {
            // set volumes if more than one
            std::istringstream strVol;
            strVol.str( nodeVol->value() );
            std::vector<gsSolidHalfFace<T> *> volFaces;
            for(int i = 0; i < nVol; i++)
            {
                volFaces.clear();
                int numFaces;
                gsGetInt(strVol, numFaces);
                for(int j = 0; j < numFaces; j++)
                {
                    int faceId;
                    gsGetInt(strVol, faceId);
                    volFaces.push_back(m->face[faceId]);
                }
                m->addVolume(volFaces);
            }
        }
        assert(ntest==n);// check if the number of surfaces in the input file is correct

        return m;
    }

    static gsXmlNode * put (const gsSolid<T> & obj,
                            gsXmlTree & data )
    {
        // Make Vertex node
        size_t nVert = obj.vertex.size();
        gsMatrix<T> vert(nVert, 3);
        for(size_t i = 0; i < nVert; i++)
        {
            vert.row(i) = obj.vertex[i]->coords;
        }
        gsXmlNode * nodeVertex = putMatrixToXml(vert, data, "Vertex");
        // Make Face node
        size_t nFace = obj.face.size();
        std::ostringstream strf;
        std::vector<int> faceVerts;
        for(size_t i = 0; i < nFace; i++)
        {
            size_t nLoops = obj.face[i]->loop.size();
            for(size_t loopIdx = 0; loopIdx < nLoops; loopIdx++)
            {
                faceVerts.clear();
                gsSolidHalfEdge<T> *eFirst = obj.face[i]->loop[loopIdx];
                gsSolidHalfEdge<T> *e = eFirst;
                do // loop over boundary of face, collecting vertex ids
                {
                    faceVerts.push_back(e->source->getId());
                    e = e->next;
                } while(e != eFirst);
                size_t nfv = faceVerts.size();
                // write # vertices, vertex ids
                strf << nfv << " ";
                for(size_t j = 0; j < nfv; j++)
                {
                    strf << faceVerts[j] << " ";
                }
                // if we aren't done, write a -2 to indicate that there are more loops
                if(loopIdx < nLoops - 1) strf << -2 << " ";
            }
            // write trim surf id (may as well make it the same as the face id)
            strf << i << "\n";
        }
        gsXmlNode* nodeFace = internal::makeNode("Face", strf.str(), data);
        // Make Volume node
        int nVol = obj.nVolumes();
        std::ostringstream strVol;
        for(int i = 0; i < nVol; i++)
        {
            size_t nVolF = obj.volume[i]->face.size();
            strVol << nVolF << " ";
            for(size_t j = 0; j < nVolF; j++)
            {
                strVol << obj.volume[i]->face[j]->getId();
                if(j < nVolF - 1) strVol << " ";
            }
            strVol << "\n";
        }
        gsXmlNode * nodeVolume = internal::makeNode("Volume", strVol.str(), data);
        // Create solid
        gsXmlNode * nodeSolid = internal::makeNode("Solid", data);
        // Set attributes for # vertices, edges, faces, volumes.
        nodeSolid->append_attribute(makeAttribute("vertices", nVert, data));
        nodeSolid->append_attribute(makeAttribute("faces", nFace, data));
        nodeSolid->append_attribute(makeAttribute("volumes", nVol, data));
        //if multiple solids are being written someone else will have to decide the id
        //nodeSolid->append_attribute(makeAttribute("id", 0, data));
        nodeSolid->append_node(nodeVertex);
        nodeSolid->append_node(nodeFace);
        nodeSolid->append_node(nodeVolume);

        // write trimmed surfaces to the root
        for(size_t i = 0; i < nFace; i++)
        {
            gsXmlNode* nodeTS = gsXml< gsTrimSurface<T> >::put(*(obj.face[i]->surf), data);
            data.appendToRoot(nodeTS);
        }

        return nodeSolid;
    }
};


/// Get a Mesh
template<class T>
class gsXml< gsMesh<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsMesh<T>);
    static std::string tag () { return "Mesh"; }
    static std::string type () { return "off"; }

    static gsMesh<T> * get (gsXmlNode * node)
    {
        assert( ( !strcmp( node->name(),"Mesh") )
                &&  ( !strcmp(node->first_attribute("type")->value(),"off") ) );
      
        gsMesh<T> * m = new gsMesh<T>;
        std::istringstream str;
        str.str( node->value() );
      
        unsigned n  = atoi ( node->first_attribute("vertices")->value() ) ;
        T x,y, z;
        for (unsigned i=0; i<n; ++i)
        {
            gsGetReal(str, x);
            gsGetReal(str, y);
            gsGetReal(str, z);
            m->addVertex(x,y,z);
        }
      
        n  = atoi ( node->first_attribute("faces")->value() ) ;
        unsigned c = 0;
        std::vector<int> face;
        for (unsigned i=0; i<n; ++i)
        {
            gsGetInt(str, c);
            face.resize(c);
            for (unsigned j=0; j<c; ++j)
                gsGetInt(str, face[j]);
            m->addFace(face);
        }
        m->cleanMesh();
        return m;
    }

    static gsXmlNode * put (const gsMesh<T> &,
                            gsXmlTree & )
    {
        return NULL;
    }
};


/// Get a Matrix from XML data
template<class T>
class gsXml< gsMatrix<T> >
{
private:
    gsXml() { }
    typedef gsMatrix<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "Matrix"; }
    static std::string type() { return ""; }
  
    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        GISMO_ASSERT( !strcmp( node->name(),"Matrix"), 
                      "Something went wrong. Expected Matrix tag." );
        
        unsigned rows  = atoi ( node->first_attribute("rows")->value() ) ;
        unsigned cols  = atoi ( node->first_attribute("cols")->value() ) ;
        getMatrixFromXml<T>(node, rows, cols, obj);
    }
    
    static gsXmlNode * put (const gsMatrix<T> & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * mat_data = putMatrixToXml(obj,data);
        // Record matrix dimensions
        mat_data->append_attribute( 
            makeAttribute("rows", obj.rows(), data) );
        mat_data->append_attribute( 
            makeAttribute("cols", obj.cols(), data) );
        
        return mat_data;
    }
};


/// Get a SparseMatrix from XML data
template<class T>
class gsXml< gsSparseMatrix<T> >
{
private:
    gsXml() { }
    typedef gsSparseMatrix<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "SparseMatrix"; }
    static std::string type() { return ""; }
  
    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        GISMO_ASSERT( !strcmp( node->name(),"SparseMatrix"), 
                      "Something went wrong. Expected SparseMatrix tag." );

        const index_t rows  = atoi ( node->first_attribute("rows")->value() ) ;
        const index_t cols  = atoi ( node->first_attribute("cols")->value() ) ;

        gsSparseEntries<T> entries;
        getSparseEntriesFromXml<T>(node, entries);

        obj.resize(rows,cols);
        obj.setFrom(entries);
    }
    
    static gsXmlNode * put (const gsSparseMatrix<T> & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * mat_data = putSparseMatrixToXml(obj,data);

        mat_data->append_attribute( 
            makeAttribute("rows", obj.rows(), data) );
        mat_data->append_attribute( 
            makeAttribute("cols", obj.cols(), data) );

        return mat_data;
    }
};

/*
 * Getting Bases from XML data
 */
    
/// Get a NurbsBasis from XML data
template<class T>
class gsXml< gsNurbsBasis<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsNurbsBasis<T>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "NurbsBasis"; }

    static gsNurbsBasis<T> * get (gsXmlNode * node)
    {
        return getRationalBasisFromXml<gsNurbsBasis<T> >(node);
    }

    static gsXmlNode * put (const gsNurbsBasis<T> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml(obj,data);
    }
};


/// Get a TensorNurbsBasis from XML data
template<short_t d, class T>
class gsXml< gsTensorNurbsBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorNurbsBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "TensorNurbsBasis"+to_string(d); }

    static gsTensorNurbsBasis<d,T> * get (gsXmlNode * node)
    {
        return getRationalBasisFromXml< gsTensorNurbsBasis<d,T> >(node);
    }
    
    static gsXmlNode * put (const gsTensorNurbsBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml< gsTensorNurbsBasis<d,T> >(obj,data);
    }
};


/// Get a Nurbs from XML data
template<class T>
class gsXml< gsNurbs<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsNurbs<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "Nurbs"; }

    static gsNurbs<T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsNurbs<T> >(node);
    }
    
    static gsXmlNode * put (const gsNurbs<T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj,data);
    }
};

/// Get a Tensor Nurbs from XML data
template<short_t d, class T>
class gsXml< gsTensorNurbs<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorNurbs<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "TensorNurbs"+to_string(d); }

    static gsTensorNurbs<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTensorNurbs<d,T> >( node );
    }
    
    static gsXmlNode * put (const gsTensorNurbs<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj,data);
    }
};

/// Get a TrimSurface from XML data
template<class T>
class gsXml< gsTrimSurface<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTrimSurface<T>);
    static std::string tag  () { return "TrimSurface"; }
    static std::string type () { return ""; }

    static gsTrimSurface<T> * get (gsXmlNode * node)
    {
        assert( !strcmp( node->name(),"TrimSurface") );
        
        gsXmlNode * tmp = node->first_node("Geometry");
        gsSurface<T> * geo =  gsXml<gsSurface<T> >::get (tmp) ;
        
        tmp = node->first_node("PlanarDomain");
        gsPlanarDomain<T> * pd  =  gsXml<gsPlanarDomain<T> >::get (tmp) ;
        
        return new gsTrimSurface<T>( geo, pd );
    }
    
    static gsXmlNode * put (const gsTrimSurface<T> & obj,
                            gsXmlTree & data )
    {
        gsXmlNode* nodeTS = internal::makeNode("TrimSurface", data);
        gsXmlNode* nodeGeom = gsXml< gsGeometry<T> >::put(*(obj.getTP()), data);
        gsXmlNode* nodeDom = gsXml< gsPlanarDomain<T> >::put(obj.domain(), data);

        nodeTS->append_node(nodeGeom);
        nodeTS->append_node(nodeDom);
        return nodeTS;
    }
};

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
        std::string s = gtype->value() ;
      
        if ( s == "BSpline"    ) 
            return gsXml< gsBSpline<T> >::get(node);
        if ( s == "Nurbs"      ) 
            return gsXml< gsNurbs<T> >::get(node);
        if ( s == "HBSpline2"  )  
            return gsXml< gsHBSpline<2,T> >::get(node); 
        if ( s == "HBSpline3"  )  
            return gsXml< gsHBSpline<3,T> >::get(node); 
        if ( s == "THBSpline1" )  
            return gsXml< gsTHBSpline<1,T> >::get(node);
        if ( s == "THBSpline2" )  
            return gsXml< gsTHBSpline<2,T> >::get(node);
        if ( s == "THBSpline3" )
            return gsXml< gsTHBSpline<3,T> >::get(node);


        if ( s == "TensorBSpline1" ) 
            return gsXml< gsTensorBSpline<1,T> >::get(node);      
        if ( s == "TensorBSpline2" ) 
            return gsXml< gsTensorBSpline<2,T> >::get(node);
        if ( s == "TensorBSpline3" ) 
            return gsXml< gsTensorBSpline<3,T> >::get(node);
        if ( s == "TensorBSpline4" ) 
            return gsXml< gsTensorBSpline<4,T> >::get(node);
        if ( s == "TensorNurbs2" ) 
            return gsXml< gsTensorNurbs<2,T> >::get(node);
        if ( s == "TensorNurbs3" ) 
            return gsXml< gsTensorNurbs<3,T> >::get(node);
        if ( s == "TensorNurbs4" ) 
            return gsXml< gsTensorNurbs<4,T> >::get(node);

        //if ( s == "TrimSurface" )
        //    return gsXml< gsTrimSurface<T> >::get(node);

        //if ( s == "TriangularBezier2" )
        //    return gsXml< gsTriangularBezier<2,T> >::get(node);

        gsWarn<<"gsXmlUtils: getGeometry: No known geometry \""<<s<<"\". Error.\n";
        return NULL;
    }


    static gsXmlNode * put (const gsGeometry<T> & obj,
                            gsXmlTree & data)
	{
	    const gsGeometry<T> * ptr = & obj;

	    if ( const gsBSpline<T> * g = 
             dynamic_cast<const gsBSpline<T> *>( ptr ) )
		    return gsXml< gsBSpline<T> >::put(*g,data);
        
	    if ( const gsNurbs<T> * g = 
             dynamic_cast<const gsNurbs<T> *>( ptr ) )
		    return gsXml< gsNurbs<T> >::put(*g,data);
        
	    if ( const gsTensorBSpline<2,T> * g = 
             dynamic_cast<const gsTensorBSpline<2,T> *>( ptr ) )
            return gsXml< gsTensorBSpline<2,T> >::put(*g,data);
        
	    if ( const gsTensorBSpline<3,T> * g = 
             dynamic_cast<const gsTensorBSpline<3,T> *>( ptr ) )
            return gsXml< gsTensorBSpline<3,T> >::put(*g,data);

	    if ( const gsTensorBSpline<4,T> * g = 
             dynamic_cast<const gsTensorBSpline<4,T> *>( ptr ) )
            return gsXml< gsTensorBSpline<4,T> >::put(*g,data);
        
	    if ( const gsTensorNurbs<2,T> * g = 
             dynamic_cast<const gsTensorNurbs<2,T> *>( ptr ) )
            return gsXml< gsTensorNurbs<2,T> >::put(*g,data);
        
	    if ( const gsTensorNurbs<3,T> * g = 
             dynamic_cast<const gsTensorNurbs<3,T> *>( ptr ) )
            return gsXml< gsTensorNurbs<3,T> >::put(*g,data);

	    if ( const gsTensorNurbs<4,T> * g = 
             dynamic_cast<const gsTensorNurbs<4,T> *>( ptr ) )
            return gsXml< gsTensorNurbs<4,T> >::put(*g,data);

	    if ( const gsTHBSpline<1,T> * g = 
             dynamic_cast<const gsTHBSpline<1,T> *>( ptr ) )
	        return gsXml< gsTHBSpline<1,T> >::put(*g,data);
        
	    if ( const gsTHBSpline<2,T> * g = 
             dynamic_cast<const gsTHBSpline<2,T> *>( ptr ) )
	        return gsXml< gsTHBSpline<2,T> >::put(*g,data);

	    if ( const gsTHBSpline<3,T> * g = 
             dynamic_cast<const gsTHBSpline<3,T> *>( ptr ) )
	        return gsXml< gsTHBSpline<3,T> >::put(*g,data);
	    
	    if ( const gsTrimSurface<T> * g = 
             dynamic_cast<const gsTrimSurface<T> *>( ptr ) )
		    return gsXml< gsTrimSurface<T> >::put(*g,data);
        
	    if ( const gsHBSpline<1,T> * g = 
	    	 dynamic_cast<const gsHBSpline<1,T> *>( ptr ) )
            return gsXml< gsHBSpline<1,T> >::put(*g,data);

	    if ( const gsHBSpline<2,T> * g = 
	    	 dynamic_cast<const gsHBSpline<2,T> *>( ptr ) )
            return gsXml< gsHBSpline<2,T> >::put(*g,data);

	    if ( const gsHBSpline<3,T> * g = 
	    	 dynamic_cast<const gsHBSpline<3,T> *>( ptr ) )
            return gsXml< gsHBSpline<3,T> >::put(*g,data);

        //if ( const gsTriangularBezier<2,T> * g =
        //     dynamic_cast<const gsTriangularBezier<2,T> *>( ptr ) )
        //    return gsXml< gsTriangularBezier<2,T> >::put(*g,data);
        
		gsWarn<<"gsXmlUtils: put Geometry: No known object "<< obj <<"Error.\n";
        return NULL;
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
        std::string s = gtype->value() ;
      
        if ( s == "BSpline"    ) 
            return gsXml< gsBSpline<T> >::get(node);
        if ( s == "Nurbs"      ) 
            return gsXml< gsNurbs<T> >::get(node);
      
        gsWarn<<"gsXmlUtils: getCurve: No known curve \""<<s<<"\". Error.\n";
        return NULL;
    }


    static gsXmlNode * put (const gsCurve<T> & obj,
                            gsXmlTree & data)
	{
	    const gsGeometry<T> * ptr = & obj;

	    if ( const gsBSpline<T> * g = 
             dynamic_cast<const gsBSpline<T> *>( ptr ) )
		    return gsXml< gsBSpline<T> >::put(*g,data);
        
	    if ( const gsNurbs<T> * g = 
             dynamic_cast<const gsNurbs<T> *>( ptr ) )
		    return gsXml< gsNurbs<T> >::put(*g,data);
                   
	    if ( const gsHBSpline<1,T> * g = 
	    	 dynamic_cast<const gsHBSpline<1,T> *>( ptr ) )
            return gsXml< gsHBSpline<1,T> >::put(*g,data);
        
		gsWarn<<"gsXmlUtils: put Curve: No known object "<< obj <<"Error.\n";
        return NULL;
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
        std::string s = gtype->value() ;
      
        if ( s == "HBSpline2"  )  
            return gsXml< gsHBSpline<2,T> >::get(node); 
        if ( s == "THBSpline2" )  
            return gsXml< gsTHBSpline<2,T> >::get(node);
      
        if ( s == "TensorBSpline2" ) 
            return gsXml< gsTensorBSpline<2,T> >::get(node);
        if ( s == "TensorNurbs2" ) 
            return gsXml< gsTensorNurbs<2,T> >::get(node);

        gsWarn<<"gsXmlUtils: getSurface: No known surface \""<<s<<"\". Error.\n";
        return NULL;
    }

    static gsXmlNode * put (const gsSurface<T> & obj,
                            gsXmlTree & data)
	{
	    const gsGeometry<T> * ptr = & obj;

	    if ( const gsTensorBSpline<2,T> * g = 
             dynamic_cast<const gsTensorBSpline<2,T> *>( ptr ) )
            return gsXml< gsTensorBSpline<2,T> >::put(*g,data);
        
	    if ( const gsTensorNurbs<2,T> * g = 
             dynamic_cast<const gsTensorNurbs<2,T> *>( ptr ) )
            return gsXml< gsTensorNurbs<2,T> >::put(*g,data);
        
	    if ( const gsTHBSpline<2,T> * g = 
             dynamic_cast<const gsTHBSpline<2,T> *>( ptr ) )
	        return gsXml< gsTHBSpline<2,T> >::put(*g,data);
	    
	    if ( const gsHBSpline<2,T> * g = 
	    	 dynamic_cast<const gsHBSpline<2,T> *>( ptr ) )
            return gsXml< gsHBSpline<2,T> >::put(*g,data);
        
		gsWarn<<"gsXmlUtils: put Geometry: No known object "<< obj <<"Error.\n";
        return NULL;
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
        GISMO_ASSERT( ( !strcmp( node->name(),"Basis") ), "Something went wrong, waiting for a basis." );
        
        gsXmlAttribute * btype = node->first_attribute("type");
        if ( ! btype )
        {
            gsWarn<< "Basis without a type in the xml file.\n";
            return NULL;
        }
        std::string s = btype->value() ;

        if ( s == "BSplineBasis" )      
            return gsXml< gsBSplineBasis<T> >::get(node);
        if ( s == "NurbsBasis"   )
            return gsXml< gsNurbsBasis<T>   >::get(node);

        if ( s == "HBSplineBasis" )
            return gsXml< gsHBSplineBasis<1,T> >::get(node);
        if ( s == "HBSplineBasis2" )
            return gsXml< gsHBSplineBasis<2,T> >::get(node);
        if ( s == "HBSplineBasis3" )
            return gsXml< gsHBSplineBasis<3,T> >::get(node);
        if ( s == "HBSplineBasis4" )
            return gsXml< gsHBSplineBasis<4,T> >::get(node);

        if ( s == "THBSplineBasis" )
            return gsXml< gsTHBSplineBasis<1,T> >::get(node);
        if ( s == "THBSplineBasis2" )
            return gsXml< gsTHBSplineBasis<2,T> >::get(node);
        if ( s == "THBSplineBasis3" )
            return gsXml< gsTHBSplineBasis<3,T> >::get(node);
        if ( s == "THBSplineBasis4" )
            return gsXml< gsTHBSplineBasis<4,T> >::get(node);

        //if ( s == "gsTriangularBezierBasis2" )
        //    return gsXml< gsTriangularBezierBasis<2,T> >::get(node);

        if ( s == "TensorBSplineBasis2" )
            return gsXml< gsTensorBSplineBasis<2, T> >::get(node);
        if ( s == "TensorBSplineBasis3" )
            return gsXml< gsTensorBSplineBasis<3, T> >::get(node);
        if ( s == "TensorBSplineBasis4" )
            return gsXml< gsTensorBSplineBasis<4, T> >::get(node);

        if ( s == "TensorNurbsBasis2" )
            return gsXml< gsTensorNurbsBasis<2, T> >::get(node);
        if ( s == "TensorNurbsBasis3" )
            return gsXml< gsTensorNurbsBasis<3, T> >::get(node);
        if ( s == "TensorNurbsBasis4" )
            return gsXml< gsTensorNurbsBasis<4, T> >::get(node);
               
        gsWarn<<"gsXmlUtils: getBasis: No known basis \""<<s<<"\". Error.\n";
        return NULL;
    }

    static gsXmlNode * put (const gsBasis<T> & obj,
                            gsXmlTree & data )
    {
        const gsBasis<T> * ptr = & obj;

        if ( const gsBSplineBasis<T> * g = 
             dynamic_cast<const gsBSplineBasis<T> *>( ptr ) )
            return gsXml< gsBSplineBasis<T> >::put(*g,data);

        if ( const gsNurbsBasis<T> * g = 
             dynamic_cast<const gsNurbsBasis<T> *>( ptr ) )
            return gsXml< gsNurbsBasis<T> >::put(*g,data);

        // Tensor B-spline
        if ( const gsTensorBSplineBasis<2, T> * g = 
             dynamic_cast<const gsTensorBSplineBasis<2, T> *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<2, T> >::put(*g,data);

        if ( const gsTensorBSplineBasis<3, T> * g = 
             dynamic_cast<const gsTensorBSplineBasis<3, T> *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<3, T> >::put(*g,data);

        if ( const gsTensorBSplineBasis<4, T> * g = 
             dynamic_cast<const gsTensorBSplineBasis<4, T> *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<4, T> >::put(*g,data);

        // Tensor Nurbs
        if ( const gsTensorNurbsBasis<2, T> * g = 
             dynamic_cast<const gsTensorNurbsBasis<2, T> *>( ptr ) )
            return gsXml< gsTensorNurbsBasis<2, T> >::put(*g,data);

        if ( const gsTensorNurbsBasis<3, T> * g = 
             dynamic_cast<const gsTensorNurbsBasis<3, T> *>( ptr ) )
            return gsXml< gsTensorNurbsBasis<3, T> >::put(*g,data);

        if ( const gsTensorNurbsBasis<4, T> * g = 
             dynamic_cast<const gsTensorNurbsBasis<4, T> *>( ptr ) )
            return gsXml< gsTensorNurbsBasis<4, T> >::put(*g,data);


        // Tensor-Hier. B-splines
        if ( const gsHTensorBasis<1,T>  * g =
             dynamic_cast<const gsHTensorBasis<1,T> *>( ptr ) )
            return gsXml< gsHTensorBasis<1,T> >::put(*g,data);

        if ( const gsHTensorBasis<2,T>  * g =
             dynamic_cast<const gsHTensorBasis<2,T> *>( ptr ) )
            return gsXml< gsHTensorBasis<2,T> >::put(*g,data);

        if ( const gsHTensorBasis<3,T>  * g =
             dynamic_cast<const gsHTensorBasis<3,T> *>( ptr ) )
            return gsXml< gsHTensorBasis<3,T> >::put(*g,data);

        if ( const gsHTensorBasis<4,T>  * g =
             dynamic_cast<const gsHTensorBasis<4,T> *>( ptr ) )
            return gsXml< gsHTensorBasis<4,T> >::put(*g,data);

        if ( const gsTHBSplineBasis<3,T>  * g =
             dynamic_cast<const gsTHBSplineBasis<3,T> *>( ptr ) )
            return gsXml< gsTHBSplineBasis<3,T> >::put(*g,data);

        //if ( const gsTriangularBezierBasis<2,T>  * g =
        //     dynamic_cast<const gsTriangularBezierBasis<2,T> *>( ptr ) )
        //    return gsXml< gsTriangularBezierBasis<2,T> >::put(*g,data);

        gsWarn<<"gsXmlUtils put: getBasis: No known basis \""<<obj<<"\". Error.\n";
        return NULL;
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
        
        std::string s = node->first_attribute("type")->value() ;
        if ( s == "PoissonPde" )
            return gsXml< gsPoissonPde<T> >::get(node);
        if ( s == "SurfacePoissonPde" )
            return gsXml< gsSurfacePoissonPde<T> >::get(node);
        
        gsWarn<<"gsXmlUtils: getPde: No known Pde \""<<s<<"\". Error.\n";
        return NULL;
    }
    
    static gsXmlNode * put (const gsPde<T> &,
                            gsXmlTree & )
    {
        return NULL;
    }
};


// Get a Multipatch
template<class T>
class gsXml< gsMultiPatch<T> >
{
private:
    gsXml() { }
    typedef gsMultiPatch<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "MultiPatch"; }
    static std::string type () { return ""; }
    
    GSXML_GET_POINTER(Object);
    
    static void get_into (gsXmlNode * node, Object & obj)
    {
        GISMO_ASSERT( !strcmp( node->name(),"MultiPatch"), 
                      "Something went wrong. Expected Multipatch tag." );
        
        // the geometry patches should be siblings of node
        gsXmlNode * toplevel = node->parent();
        
        const int d = atoi( node->first_attribute("parDim")->value() );
        
        gsXmlNode * tmp = node->first_node("patches");
        std::istringstream str ;
        str.str( tmp->value() );
        
        std::vector< gsGeometry<T> *> patches;
        std::map<int,int> ids;
        if ( ! strcmp( tmp->first_attribute("type")->value(),"id_range") )
        {
            int first, last;
            gsGetInt(str, first);
            gsGetInt(str, last);
            for ( int i = first; i<=last; ++i )
            {
                GISMO_ASSERT( searchId(i, toplevel) != NULL, 
                              "No Geometry with Id "<<i<<" found in the XML data.");
                patches.push_back( getById< gsGeometry<T> >( toplevel, i ) );
                patches.back()->setId(i);
                ids[i] = i - first;
            }
        }
        else if ( ! strcmp( tmp->first_attribute("type")->value(),"id_index") )
        {
            int c = 0;
            for (int pindex; gsGetInt(str, pindex);)
            {
                GISMO_ASSERT( searchId(pindex, toplevel) != NULL, 
                              "No Geometry with Id "<<pindex<<" found in the XML data.");
                patches.push_back( getById< gsGeometry<T> >( toplevel, pindex ) );
                patches.back()->setId(pindex);
                ids[pindex] = c++;
            }
        }
        else
        {
            gsWarn<<"Unknown tag in XML multipatch object.\n";
        }


        //patches: 2 0 1
        // before offset range: 5 3 4

        // Boundaries and interfaces are also 3,4,5 so we need to translate them t0 0,1,2
        
        // Read boundary
        std::vector< patchSide > boundaries;
        tmp = node->first_node("boundary");
        if (tmp)
            getBoundaries(tmp, ids, boundaries);
        
        // Read interfaces
        std::vector< boundaryInterface > interfaces;
        tmp = node->first_node("interfaces");
        if (tmp)
            getInterfaces(tmp, d, ids, interfaces);

        obj = gsMultiPatch<T>(patches, boundaries, interfaces);        
    }

    static gsXmlNode * put (const gsMultiPatch<T> & obj,
                            gsXmlTree & data)
    {
        // First insert all geometries
        int max_id = data.maxId();
        gsXmlNode * tmp;
        for ( typename gsMultiPatch<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            tmp = gsXml<gsGeometry<T> >::put(**it,data);
            data.appendToRoot(tmp);
        }
        
        std::ostringstream str;
        str<< max_id+1 <<" "<< data.maxId();
        tmp = internal::makeNode("patches" , str.str(), data);
        tmp->append_attribute( internal::makeAttribute("type", "id_range", data) );
        str.clear(); str.str("");
        
        // Make MultiPatch node
        gsXmlNode * mp_node = internal::makeNode("MultiPatch" , data);
        mp_node->append_attribute( internal::makeAttribute("parDim", obj.parDim() , data) );
        mp_node->append_node(tmp);
      
        appendBoxTopology(obj, mp_node, data);

        return mp_node;
    }
    
};

/// Get a MultiBasis from XML data
template <class T>
class gsXml< gsMultiBasis<T> >
{
private:
    gsXml() { }
    typedef gsMultiBasis<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    GSXML_GET_POINTER(Object);
    static std::string tag() { return "MultiBasis"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode* node, Object & result)
    {
        GISMO_ASSERT( !strcmp( node->name(), "MultiBasis" ),
                      "Something went wrong. Expected MultiBasis tag." );

        gsXmlNode* topLevel = node->parent();

        const int d = atoi( node->first_attribute("parDim")->value() );

        gsXmlNode* patchNode = node->first_node("patches");
        std::istringstream iss;
        iss.str( patchNode->value() );

        typename gsMultiBasis<T>::BasisContainer bases;
        std::map<int, int> ids;
        if ( !strcmp( patchNode->first_attribute("type")->value(), "id_range") )
        {
            int first, last;
            gsGetInt(iss, first);
            gsGetInt(iss, last);
            for (int i = first; i <= last; ++i)
            {
                bases.push_back( getById< gsBasis<T> >( topLevel, i ) );
                ids[i] = i - first;
            }
        }
        else if ( !strcmp( patchNode->first_attribute("type")->value(), "id_index") )
        {
            int c = 0;
            for ( int pindex; gsGetInt(iss, pindex); )
            {
                bases.push_back( getById< gsBasis<T> >( topLevel, pindex ) );
                ids[pindex] = c++;
            }
        }
        else
        {
            gsWarn << "unknown tag in XML multipatch object \n";
        }

        // Read boundary
        std::vector< patchSide > boundaries;
        gsXmlNode * tmp = node->first_node("boundary");
        if (tmp)
            getBoundaries(tmp, ids, boundaries);
        
        // Read interfaces
        std::vector< boundaryInterface > interfaces;
        tmp = node->first_node("interfaces");
        if (tmp)
            getInterfaces(tmp, d, ids, interfaces);

        gsBoxTopology topology( d, bases.size(), boundaries, interfaces);

        result = gsMultiBasis<T>(bases, topology);
        freeAll(bases);
    }

    static gsXmlNode* put(const gsMultiBasis<T>& obj,
                          gsXmlTree& data)
    {
        // Insert all the basis
        int max_id = data.maxId();
        for ( typename gsMultiBasis<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            gsXmlNode* basisXml = gsXml< gsBasis<T> >::put(**it, data);
            data.appendToRoot( basisXml );
        }

        std::ostringstream oss;
        oss<<  max_id+1 <<" "<< data.maxId();
        gsXmlNode* node = internal::makeNode("patches", oss.str(), data);
        node->append_attribute( internal::makeAttribute("type", "id_range", data) );
        oss.clear();
        oss.str("");

        gsXmlNode* mbNode = internal::makeNode(tag(), data);
        mbNode->append_attribute( internal::makeAttribute("parDim", obj.dim(), data) );
        mbNode->append_node(node);

        appendBoxTopology(obj.topology(), mbNode, data);

        return mbNode;
    }

};


/// Get a PlanarDomain from XML data
template<class T>
class gsXml< gsPlanarDomain<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsPlanarDomain<T>);
    static std::string tag  () { return "PlanarDomain"; }
    static std::string type () { return ""; }

    static gsPlanarDomain<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"PlanarDomain"), 
                      "Something went wrong. Expected PlanarDomain tag." );
      
        std::vector<gsCurveLoop<T>*> loops;
        for (gsXmlNode * tmp = node->first_node("CurveLoop"); 
             tmp; tmp = tmp->next_sibling("CurveLoop"))
            loops.push_back( gsXml<gsCurveLoop<T> >::get(tmp) ) ;

        return new gsPlanarDomain<T>( loops );
    }

    static gsXmlNode * put (const gsPlanarDomain<T> & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * pl = internal::makeNode("PlanarDomain", data);
		gsXmlNode* tmp;

        // get number of loops
        int nl = obj.numLoops();

        for (int i=0; i!=nl; ++i)
        {
            tmp = internal::gsXml< gsCurveLoop<T> >::put(obj.loop(i), data );
            tmp->append_attribute( makeAttribute("index", i, data) );
            pl->append_node(tmp);
        }    
        return pl;
    }
};

/// Get a CurveLoop from XML data
template<class T>
class gsXml< gsCurveLoop<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsCurveLoop<T>);
    static std::string tag  () { return "CurveLoop"; }
    static std::string type () { return ""; }

    static gsCurveLoop<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"CurveLoop"), 
                      "Something went wrong. Expected CurveLoop tag." );
      
        std::vector<gsCurve<T>* > curves;

        for (gsXmlNode * tmp = node->first_node("Geometry"); 
             tmp; tmp = tmp->next_sibling("Geometry"))
            curves.push_back( gsXml<gsCurve<T> >::get(tmp) ) ;
      
        return new gsCurveLoop<T>( curves );
    }

    static gsXmlNode * put (const gsCurveLoop<T> & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * cl = internal::makeNode("CurveLoop", data);
		gsXmlNode* tmp;

        // get number of curves
        int nc = obj.numCurves();

        for (int i=0; i!=nc; ++i)
        {
            tmp = internal::gsXml< gsGeometry<T> >::put(obj.curve(i), data );
            tmp->append_attribute( makeAttribute("index", i, data) );
            cl->append_node(tmp);
        }    
        return cl;
    }
};


/// Get a Curve fitting class data from XML
template<class T>
class gsXml< gsCurveFitting<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsCurveFitting<T>);
    static std::string tag () { return "CurveFitting"; }
    static std::string type() { return ""; }

    static gsCurveFitting<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"CurveFitting"), 
                      "Something went wrong. Expected CurveFitting tag." );
      
        bool closed = (atoi(node->first_attribute("closed")->value() ) != 0);

        // Read knot-vector
        gsXmlNode   * tmp = node->first_node("KnotVector");
        gsKnotVector<T> * kv = gsXml< gsKnotVector<T> >::get(tmp);

        // Read parameter values
        tmp = node->first_node("Matrix");
        gsMatrix<T> * parval = gsXml< gsMatrix<T> >::get(tmp);

        // Read points
        tmp = tmp->next_sibling("Matrix");
        gsMatrix<T> * pts =  gsXml< gsMatrix<T> >::get(tmp);

        gsCurveFitting<T> * cf = new gsCurveFitting<T>(*parval,*pts,*kv,closed);
        delete parval;
        delete pts;
        delete kv;
        return cf ;
    }

    static gsXmlNode * put (const gsCurveFitting<T> &,
                            gsXmlTree & )
    {
        return NULL;
    }
};


/// Get a Poisson Pde from XML data
template<class T>
class gsXml< gsPoissonPde<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsPoissonPde<T>);
    static std::string tag () { return "Pde"; }
    static std::string type () { return "PoissonPde"; }

    static gsPoissonPde<T> * get (gsXmlNode * node)
    {
        assert( ( !strcmp( node->name(),"Pde") ) && 
                ( !(
                    strcmp( node->first_attribute("type")->value(),"PoissonPde")
                    && strcmp( node->first_attribute("type")->value(),"SurfacePoissonPde")
                     )) );

        // Read the dimension
        GISMO_ASSERT( node->first_attribute("dim"), "xml reader: No dim found" ) ;
        short_t d = atoi( node->first_attribute("dim")->value() );

        
        unsigned tDim = 0;
        gsXmlAttribute * targetDim = node->first_attribute("targetDim");
        
        if ( targetDim )
            tDim = atoi( targetDim->value() );

        if ( tDim >= 1 )
        {
            gsXmlNode * tmp = node->first_node("rhs");
            gsFunctionExpr<T>  rhs_fnct;
            getFunctionFromXml(tmp, rhs_fnct);
            
            tmp = node->first_node("solution");
            if ( tmp )
            {
                gsFunctionExpr<T> msol;
                getFunctionFromXml(tmp, msol);
                
                return new gsPoissonPde<T>(rhs_fnct, d, msol );
            }
            else
            {
                return new gsPoissonPde<T>( rhs_fnct, d );
            }
        }

        // Read right hand side function
        gsXmlNode * tmp = node->first_node("rhs");	
        gsFunctionExpr<T> rhs(tmp->value(), d);

        // Read exact solution, if one exists in the file
        tmp = node->first_node("solution");	
        if ( tmp )
        {
            gsFunctionExpr<T> sol(tmp->value(), d);
            //gsDebugVar (*sol);
            return new gsPoissonPde<T>(rhs, d, sol );
        }
        else
        {
            return new gsPoissonPde<T>( rhs, d );
        }
    }
    
    static gsXmlNode * put (const gsPoissonPde<T> &,
                            gsXmlTree & )
    {
        return NULL;
    }

};

/*
/// Get a Poisson Pde from XML data
template<class T>
class gsXml< gsSurfacePoissonPde<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsSurfacePoissonPde<T>);
    static std::string tag () { return "Pde"; }
    static std::string type () { return "SurfacePoissonPde"; }

    static gsSurfacePoissonPde<T> * get (gsXmlNode * node)
    {
        assert( ( !strcmp( node->name(),"Pde") ) && 
                ( !strcmp( node->first_attribute("type")->value(),"SurfacePoissonPde") ) );

        // Read the dimension
        assert( node->first_attribute("dim") ) ;
        unsigned d = atoi( node->first_attribute("dim")->value() );

        // Read right hand side function
        gsXmlNode   * tmp = node->first_node("rhs");	
        gsFunctionExpr<T> * rhs = new gsFunctionExpr<T>(tmp->value());

        // Read exact solution, if one exists in the file
        tmp = node->first_node("solution");	
        if ( tmp )
        {
            gsFunctionExpr<T> * sol = new gsFunctionExpr<T>(tmp->value());
            //gsDebugVar (*sol);
            return new gsSurfacePoissonPde<T>(rhs, d, sol );
        }
        else
        {
            return new gsSurfacePoissonPde<T>( rhs, d );
        }
    }
    
    static gsXmlNode * put (const gsSurfacePoissonPde<T> & obj,
                            gsXmlTree & data )
    {
        return NULL;
    }

};
*/

/*
/// Get a Boundary Value Problem from XML
template<class T>
class gsXml< gsBVProblem<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBVProblem<T>);
    static std::string tag () { return "BVProblem"; }
    static std::string type() { return ""; }

    static gsBVProblem<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"BVProblem"), 
                      "Something went wrong. Expected BVProblem tag." );

        gsBVProblem<T> * bvp;

        // Read the Pde
        gsPde<T>  * pde = gsXml< gsPde<T> >::get( node->first_node("Pde") );
	
        // Read domain
        int domain = atoi(node->first_attribute("domain")->value() );
        gsXmlNode * toplevel = node->parent();// the geometry patches should be siblings of node
        std::string dtag; // = getTag(toplevel, domain );
        for (gsXmlNode * child = node->first_node(); 
             child; child = child->next_sibling() )
            if (  atoi(child->first_attribute("id")->value() ) == domain )
            {
                dtag = child->name();
                break;
            }
    
        if ( dtag == "Geometry" )
        {
            gsGeometry<T> * geo = getById< gsGeometry<T> >(toplevel,domain);
            bvp = new gsBVProblem<T>(geo, pde);
        }
        else if ( dtag == "MultiPatch" )
        {
            // to do: memory to delete?
            gsMultiPatch<T> * mp = getById< gsMultiPatch<T> >(toplevel,domain);
            bvp = new gsBVProblem<T>(*mp, pde);
        }
        else
        {
            GISMO_ERROR("Invalid tag");
        }
	
        // Read in boundary conditions
        for (gsXmlNode * child = node->first_node("bc"); 
             child; child = child->next_sibling("bc") )
        {
            gsFunctionExpr<T> * ff = 
                new gsFunctionExpr<T>(child->first_attribute("function")->value() );
            std::istringstream str;
            str.str( child->value() );
        
            if ( !strcmp(child->first_attribute("type")->value(), "dirichlet") )
            {		       
                for (int side; gsGetInt(str, side);) 
                    bvp->addCondition( static_cast<boxSide>(side),
                                       condition_type::dirichlet, ff);
            }
            else if ( !strcmp(child->first_attribute("type")->value(), "neumann") )
            {		       
                for (int side; gsGetInt(str, side);) 
                    bvp->addCondition( static_cast<boxSide>(side),
                                       condition_type::neumann, ff);
            }		
        }
    
        return bvp ;
    }

    static gsXmlNode * put (const gsBVProblem<T> & obj,
                            gsXmlTree & data )
    {
        return NULL;
    }
};
*/

}// end namespace internal

}// end namespace gismo

//#undef GSXML_COMMON_FUNCTIONS
//#undef TMPLA2
//#undef TMPLA3
