/** @file gsModelingXml.hpp

    @brief XML serialization for the gsModeling types (implementation
    header).

    Moved here from gsIO/gsXmlUtils.hpp: each module owns its gsXml
    specializations (modularization stream S3, step A6).

    ONLY include this file from gsXmlRegistration_.cpp (lib mode) or through
    gsModelingXmlRegistration.h in header-only mode - a second inclusion in another
    translation unit of the library would poison the symbol visibility
    of the explicit instantiations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>
#include <gsModeling/gsSolid.h>
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsCurveFitting.h>

namespace gismo {
namespace internal {

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

} // namespace internal
} // namespace gismo
