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
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsMFunctionExpr.h>
#include <gsCore/gsRationalBasis.h>
#include <gsCore/gsMultiPatch.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsSolid.h>
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsPde/gsBVProblem.h>
#include <gsModeling/gsCurveFitting.h>
#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsHBSpline.h>
//#include <gsTrBezier/gsTriangularBezierBasis.h>
//#include <gsTrBezier/gsTriangularBezier.h>
#include <gsPde/gsPoissonPde.h>
#include <gsPde/gsSurfacePoissonPde.h>

namespace gismo {

namespace internal {

/* // Note: no default impl.
/// Generic get XML class: specializations provide implementation
template<class Object>
std::string gsXml<Object>::tag ()
{
    //Next line will produce compile-time error
    //Object::Object_does_not_exist_ERROR;
    return "";
}
    
template<class Object>
std::string gsXml<Object>::type ()
{
    // Produce compile-time error
    //Object::Object_does_not_exist_ERROR;
    return "";
}
    
template<class Object>
Object * gsXml<Object>::get (gsXmlNode * node)
{
    //Next line will produce compile-time error
    //Object::Object_does_not_exist_ERROR;

    GISMO_ERROR("XML reading error. Is the object readable?");
    return NULL;
}

template<class Object>
gsXmlNode * gsXml<Object>::put (const Object & obj, gsXmlTree & data)
{
    //Next line will produce compile-time error
    //Object::Object_does_not_exist_ERROR;

    GISMO_ERROR("XML reading error. Is the object readable?");
    return NULL;
}
*/
 
////////////////////////////////////////////////////////
// Implementations of common XML get/put functions
////////////////////////////////////////////////////////
   
template<class T>
void getMatrixFromXml ( gsXmlNode * node, unsigned const & rows, 
                        unsigned const & cols, gsMatrix<T> & result ) 
{
    //gsWarn<<"Reading "<< node->name() <<" matrix of size "<<rows<<"x"<<cols<<"Geometry..\n";
    std::istringstream str;
    str.str( node->value() );
    result.resize(rows,cols);
 
    for (unsigned i=0; i<rows; ++i)
        for (unsigned j=0; j<cols; ++j)
            if ( !(str >> result(i,j) ) )
            {
                gsWarn<<"XML Warning: Reading matrix of size "<<rows<<"x"<<cols<<" failed.\n";
                gsWarn<<"Tag: "<< node->name() <<", Matrix entry: ("<<i<<", "<<j<<").\n";
                return;
            }
}

template<class T>
void getSparseEntriesFromXml ( gsXmlNode * node,
                              gsSparseEntries<T> & result ) 
{
    result.clear();

    std::istringstream str;
    str.str( node->value() );
    index_t r,c;
    T val;

    while( (str >> r) && (str >> c) && (str >> val) ) 
        result.add(r,c,val);
}


template<class T>
void getFunctionFromXml ( gsXmlNode * node, gsMFunctionExpr<T> & result ) 
{
    gsWarn<<"Reading "<< node->name() <<" function\n";

    GISMO_ASSERT( node->first_attribute("dim"), "xml reader: No dim found" ) ;
    int d = atoi( node->first_attribute("dim")->value() );
  
    std::vector< std::string > expr_strings;

    for (gsXmlNode * child = node->first_node("component"); 
         child; child = child->next_sibling() )
    {
        expr_strings.push_back(  child->value() );
    }

    result = gsMFunctionExpr<T>( expr_strings, d );
}


template<class T>
gsXmlNode * putMatrixToXml ( gsMatrix<T> const & mat, gsXmlTree & data, std::string name) 
{
    std::ostringstream str;
    str << std::setprecision(FILE_PRECISION);
    // Write the matrix entries
    for (index_t i=0; i< mat.rows(); ++i)
    {
        for (index_t j=0; j<mat.cols(); ++j)
            str << mat(i,j)<< " ";
        str << "\n";
    }

    // Create XML tree node
    gsXmlNode* new_node = internal::makeNode(name, str.str(), data);        
    return new_node;
};

template<class T>
gsXmlNode * putSparseMatrixToXml ( gsSparseMatrix<T> const & mat, 
                                   gsXmlTree & data, std::string name)
{
    typedef typename gsSparseMatrix<T>::InnerIterator cIter;

    std::ostringstream str;
    str << std::setprecision(FILE_PRECISION);
    const index_t nCol = mat.cols();

    for (index_t j=0; j != nCol; ++j) // for all columns
        for ( cIter it(mat,j); it; ++it ) // for all non-zeros in column
        {
            // Write the matrix entry
            str <<it.index() <<" "<<j<<" " << it.value() << "\n";
        }
    
    // Create XML tree node
    gsXmlNode* new_node = internal::makeNode(name, str.str(), data);        
    return new_node;
};

template<class T>
gsXmlNode * makeNode( const std::string & name, 
                      const gsMatrix<T> & value, gsXmlTree & data,
                      bool transposed)
{
    std::ostringstream oss;
    // Set precision
    oss << std::setprecision(FILE_PRECISION);
  
    if ( transposed )
        for ( index_t j = 0; j< value.rows(); ++j)
        {
            for ( index_t i = 0; i< value.cols(); ++i)
                oss << value(j,i) <<" ";
            //oss << "\n";
        }
    else
        for ( index_t j = 0; j< value.cols(); ++j)
        {
            for ( index_t i = 0; i< value.rows(); ++i)
                oss << value(i,j) <<" ";
            //oss << "\n";
        }
  
    return makeNode(name, oss.str(), data);
}

template<class Object>
Object * getTensorBasisFromXml ( gsXmlNode * node)
{
    GISMO_ASSERT( !strcmp( node->name(),"Basis"), "Wrong tag name." ) ;
    //&&  ( !strcmp(node->first_attribute("type")->value(),
    //              internal::gsXml<Object>::type().c_str() ) ),
    //"Reading tensor basis failed because of wrong tag name.");
    
	// Component container
	std::vector<typename Object::CoordinateBasis* > bb;

    // Special case of reading a 1D tensor basis
    if ( !strcmp(node->first_attribute("type")->value(),
                 gsXml<typename Object::CoordinateBasis>::type().c_str() ) )
    {
        bb.push_back( gsXml<typename Object::CoordinateBasis>::get(node) );
        return new Object( bb );
    }

    gsXmlNode * tmp = node->first_node("Basis");
    GISMO_ASSERT( tmp , "Wrong data in the xml file.");
    unsigned d = Object::Dim;
	for ( unsigned i = 0; i!=d; ++i)
    {
        bb.push_back( gsXml<typename Object::CoordinateBasis>::get(tmp) );
        tmp =  tmp->next_sibling("Basis");        
    }
    //gsDebugVar( bb.size() );
    return new Object( bb );
}

template<class Object>
gsXmlNode * putTensorBasisToXml ( Object const & obj, gsXmlTree & data)
{
    // Write the component bases
    static const unsigned d = Object::Dim;
    if (d==1)
        return gsXml<typename Object::CoordinateBasis>::put(obj.component(0), data);

    // Add a new node (without data)
    gsXmlNode* tp_node = internal::makeNode("Basis" , data);        
    tp_node->append_attribute( makeAttribute("type", 
                                             internal::gsXml<Object>::type().c_str(), data) );
    
	for ( unsigned i = 0; i!=d; ++i )
    {
		gsXmlNode* tmp = 
		    internal::gsXml< typename Object::CoordinateBasis >::put(obj.component(i), data );
		tmp->append_attribute( makeAttribute("index", i, data) );
		tp_node->append_node(tmp);
    }
    
    // All set, return the basis
    return tp_node;
}

template<class Object>
Object * getRationalBasisFromXml ( gsXmlNode * node)
{
    GISMO_ASSERT( ( !strcmp( node->name(),"Basis") )
                  &&  ( !strcmp(node->first_attribute("type")->value(),
                                internal::gsXml<Object>::type().c_str() ) ),
                  "Something is wrong with the XML data: There should be a node with a "<<internal::gsXml<Object>::type().c_str()<<" Basis.");

    // Read source basis
    gsXmlNode * tmp = node->first_node("Basis");
    typename Object::SourceBasis * src = gsXml<typename Object::SourceBasis>::get(tmp) ;
       
    // Read weights
    tmp = node->first_node("weights");
    gsMatrix<typename Object::Scalar_t>   weights;
    getMatrixFromXml<typename Object::Scalar_t>( tmp, src->size(), 1, weights);
    return new Object( src, give(weights) );
}


template<class Object>
Object * getHTensorBasisFromXml ( gsXmlNode * node)
{
    GISMO_ASSERT( ( !strcmp( node->name(),"Basis") )
                  &&  ( !strcmp(node->first_attribute("type")->value(),
                                internal::gsXml<Object>::type().c_str() ) ),
                  "Something is wrong with the XML data: There should be a node with a "<<
                  internal::gsXml<Object>::type().c_str()<<" Basis.");

    typedef typename Object::Scalar_t T;
    static const int d = Object::Dim;
    
    // Read max level
    //unsigned lvl = atoi( node->first_attribute("levels")->value() );
    gsXmlNode * tmp = node->first_node("Basis");
    GISMO_ASSERT( tmp , "Expected to find a basis node.");
    
    // Read the Tensor-product basis
    gsTensorBSplineBasis<d,T> * tp = 
        gsXml<gsTensorBSplineBasis<d,T> >::get(tmp);
    
    // Initialize the HBSplineBasis
    std::istringstream str;
    
    // Insert all boxes
    unsigned c;
    std::vector<unsigned int> all_boxes;
    for (tmp = node->first_node("box"); 
         tmp; tmp = tmp->next_sibling("box"))
    {
        all_boxes.push_back(atoi( tmp->first_attribute("level")->value() ));
        str.clear();
        str.str( tmp->value() );
        for( unsigned i = 0; i < 2*d; i++)
        {
            str>> c;
            all_boxes.push_back(c);
        }
    }
    Object * hbs = new Object(*tp, all_boxes);
    delete tp;
    return hbs;
}

template<class Object>
gsXmlNode * putHTensorBasisToXml ( Object const & obj, gsXmlTree & data)
{
    typedef typename Object::Scalar_t T;
    const int d = obj.dim();

    // Add a new node (without data)
    gsXmlNode* tp_node = internal::makeNode("Basis" , data);
    
    tp_node->append_attribute( makeAttribute("type",
                                             internal::gsXml<Object>::type().c_str(), data) );

    //tp_node->append_attribute( makeAttribute( "levels",2 ,data )); // deprecated
  
    // Write the component bases
    gsXmlNode * tmp = putTensorBasisToXml(obj.tensorLevel(0), data);
    tp_node->append_node(tmp);
    
    //Output boxes
    gsMatrix<unsigned> box(1,2*d);

    for( typename Object::hdomain_type::const_literator lIter = 
             obj.tree().beginLeafIterator(); lIter.good() ; lIter.next() )
    {
        if ( lIter->level > 0 )
        {
            box.leftCols(d)  = lIter.lowerCorner().transpose();
            box.rightCols(d) = lIter.upperCorner().transpose();
       
            tmp = putMatrixToXml( box, data, "box" );
           
            tmp->append_attribute( makeAttribute("level", to_string(lIter->level), data ) );
            tp_node->append_node(tmp);
        }
    }

/*
// Write box history (deprecated)
typename Object::boxHistory const & boxes = obj.get_inserted_boxes();

for(unsigned int i = 0; i < boxes.size(); i++)
{
box.leftCols(d)  = boxes[i].lower.transpose();
box.rightCols(d) = boxes[i].upper.transpose();

tmp = putMatrixToXml( box, data, "box" );
tmp->append_attribute( makeAttribute("level", to_string(boxes[i].level), data ) );
tp_node->append_node(tmp);
}
//*/

    // All set, return the basis
    return tp_node;
}


template<class Object>
gsXmlNode * putRationalBasisToXml ( Object const & obj, gsXmlTree & data)
{
    // Add a new node
    gsXmlNode* rat_node = internal::makeNode("Basis" , data);        
    rat_node->append_attribute( makeAttribute("type",
                                              internal::gsXml< Object >::type().c_str(), data) );
	
    // Write the source basis
	gsXmlNode* tmp = 
        internal::gsXml< typename Object::SourceBasis >::put(obj.source(), data );
    rat_node->append_node(tmp);
    
    // Write the weights
    tmp = putMatrixToXml( obj.weights(), data, "weights" );
    rat_node->append_node(tmp);
	
	// All done, return the node
	return rat_node;
}

/*
template<class Object>
Object * getById(gsXmlNode * node, const int & id)
{
    std::string tag = internal::gsXml<Object>::tag();
    for (gsXmlNode * child = node->first_node(tag.c_str());
         child; child = child->next_sibling(tag.c_str()))
    {
        if (  atoi(child->first_attribute("id")->value() ) == id )
            return internal::gsXml<Object>::get(child);
    }
    std::cerr<<"gsXmlUtils Warning: "<< internal::gsXml<Object>::tag() 
             <<" with id="<<id<<" not found.\n";
    return NULL;
}
*/


// Helper to get the tag of an id
// std::string getTag(gsXmlNode * node, const int & id)
// {
//     for (gsXmlNode * child = node->first_node(); 
// 	 child; child = child->next_sibling() )
// 	if (  atoi(child->first_attribute("id")->value() ) == id )
// 	    return child->name();
    
//     std::cerr<<"gsXmlUtils Warning: Tag with id="<<id<<" not found.\n";
//     return "";
// }
  

////////////////////////////////////////////////////////
// Getting Xml data
////////////////////////////////////////////////////////

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
            str >>std::ws>>x>>std::ws>>y>>std::ws>>z>>std::ws;        
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
                strf>>std::ws>> nf >>std::ws; // read num vertices on this loop
                for (int ivert=0; ivert<nf; ivert++) // read vertices
                {
                    strf>>std::ws>> vertID >>std::ws;
                    vert[nLoops-1].push_back(m->vertex[vertID]);
                }
                // next number is either:
                // * the trim surface id,
                // * -1 to indicate that the surface is automatically computed, or
                // * -2 to indicate that there are further internal loops
                strf>>std::ws>> trimID >>std::ws;
            } while(trimID <= -2); // -2 indicates that there are vertices remaining
            if (trimID>-1)
                m->addFace(vert, getById< gsTrimSurface<T> >( toplevel, trimID ) );
            else if (trimID==-1 && nLoops == 1)
                m->addFace(vert[0]);
            else if (trimID==-1)
            {
                gsWarn<<"\nAutomatic creation of trimmed surfaces is only supported for a single loop\n";
                exit(1);
            }
            else
            {
                gsWarn<<"\n ID of the trimmed surface trimID=" <<trimID<<" is invalid (must be >=-1)\n";
                exit(1);
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
                strVol >> std::ws >> numFaces;
                for(int j = 0; j < numFaces; j++)
                {
                    int faceId;
                    strVol >> std::ws >> faceId;
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
        gsXmlNode * root = data.first_node("xml");
        for(size_t i = 0; i < nFace; i++)
        {
            gsXmlNode* nodeTS = gsXml< gsTrimSurface<T> >::put(*(obj.face[i]->surf), data);
            nodeTS->append_attribute(makeAttribute("id", i, data) );
            root->append_node(nodeTS);
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
            str >>std::ws>>x>>std::ws>>y>>std::ws>>z>>std::ws;
            m->addVertex(x,y,z);
        }
      
        n  = atoi ( node->first_attribute("faces")->value() ) ;
        unsigned c;
        std::vector<int> face;
        for (unsigned i=0; i<n; ++i)
        {
            str >>std::ws>> c ;
            face.resize(c);
            for (unsigned j=0; j<c; ++j)
                str >>std::ws>> face[j] ;
            m->addFace(face);
        }
        return m;
    }

    static gsXmlNode * put (const gsMesh<T> & obj,
                            gsXmlTree & data )
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
public:
    GSXML_COMMON_FUNCTIONS(gsMatrix<T>);
    static std::string tag () { return "Matrix"; }
    static std::string type() { return ""; }
  
    static gsMatrix<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"Matrix"), 
                      "Something went wrong. Expected Matrix tag." );
        
        unsigned rows  = atoi ( node->first_attribute("rows")->value() ) ;
        unsigned cols  = atoi ( node->first_attribute("cols")->value() ) ;
        gsMatrix<T> * tmp = new gsMatrix<T>;
        getMatrixFromXml<T>(node,rows,cols, *tmp);
        return tmp;
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
public:
    GSXML_COMMON_FUNCTIONS(gsSparseMatrix<T>);
    static std::string tag () { return "SparseMatrix"; }
    static std::string type() { return ""; }
  
    static gsSparseMatrix<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"SparseMatrix"), 
                      "Something went wrong. Expected SparseMatrix tag." );

        const index_t rows  = atoi ( node->first_attribute("rows")->value() ) ;
        const index_t cols  = atoi ( node->first_attribute("cols")->value() ) ;
        gsSparseMatrix<T> * result = new gsSparseMatrix<T>(rows,cols);

        gsSparseEntries<T> entries;
        getSparseEntriesFromXml<T>(node, entries);

        result->setFrom(entries);
        //result->makeCompressed(); // needed ?
        return result;
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

////////////////////////////////////////////////////////
// Getting Bases from XML data
////////////////////////////////////////////////////////

    
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

        // (!) Deprecated 
/*      assert( ( !strcmp( node->name(),"Basis") )
        &&  ( !strcmp(node->first_attribute("type")->value(),"NurbsBasis") ) );
        
        gsXmlNode * tmp = node->first_node("KnotVector");
        // if type: == Plain, == Compact .. 
        gsKnotVector<T> kv = * safe( gsXml<gsKnotVector<T> >::get (tmp) );

        return new gsNurbsBasis<T>( kv, kv.degree() ); */
    }

    static gsXmlNode * put (const gsNurbsBasis<T> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml(obj,data);
    }
};

/// Get a HTensorBasis from XML data
template<unsigned d, class T>
class gsXml< gsHTensorBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsHTensorBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return ""; } // tag ?

    static gsHTensorBasis<d,T> * get (gsXmlNode * node)
    {
        gsXmlAttribute * btype = node->first_attribute("type");
        if ( ! btype )
        {
            gsWarn<< "Basis without a type in the xml file.\n";
            return NULL;
        }
        std::string s = btype->value() ;
        if ( s.compare(0, 9, "HBSplineB" , 9 ) == 0 ) // needs correct d as well
            return gsXml< gsHBSplineBasis<d,T> >::get(node);
        if ( s.compare(0, 10,"THBSplineB", 10) == 0 )
            return gsXml< gsTHBSplineBasis<d,T> >::get(node);

        gsWarn<<"gsXmlUtils: gsHTensorBasis: No known basis \""<<s<<"\". Error.\n";
        return NULL;
    }
    
    static gsXmlNode * put (const gsHTensorBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        const gsBasis<T> * ptr = & obj;

        // Hier. B-splines
        if ( const gsHBSplineBasis<d,T>  * g =
             dynamic_cast<const gsHBSplineBasis<d,T> *>( ptr ) )
            return gsXml< gsHBSplineBasis<d,T> >::put(*g,data);
        
        // Truncated hier. B-splines
        if ( const gsTHBSplineBasis<d,T>  * g =
             dynamic_cast<const gsTHBSplineBasis<d,T> *>( ptr ) )
            return gsXml< gsTHBSplineBasis<d,T> >::put(*g,data);
        
        gsWarn<<"gsXmlUtils put: getBasis: No known basis \""<<obj<<"\". Error.\n";
        return NULL;
    }
};


/// Get a TensorNurbsBasis from XML data
template<unsigned d, class T, class Basis_t>
class gsXml< gsTensorNurbsBasis<d,T,Basis_t> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorNurbsBasis<TMPLA3(d,T,Basis_t)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "TensorNurbsBasis"+to_string(d); }

    static gsTensorNurbsBasis<d,T,Basis_t> * get (gsXmlNode * node)
    {
        return getRationalBasisFromXml< gsTensorNurbsBasis<d,T,Basis_t> >(node);
    }
    
    static gsXmlNode * put (const gsTensorNurbsBasis<d,T,Basis_t> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml< gsTensorNurbsBasis<d,T,Basis_t> >(obj,data);
    }
};


/// Get a HBSpline from XML data
template<class T, unsigned d>
class gsXml< gsHBSpline<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsHBSpline<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "HBSpline"+to_string(d); }

    static gsHBSpline<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsHBSpline<d,T> >(node);
    }
    
    static gsXmlNode * put (const gsHBSpline<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsHBSpline<d,T> >(obj,data);
    }
};


/// Get a THBSpline from XML data
template<unsigned d, class T>
class gsXml< gsTHBSpline<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTHBSpline<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "THBSpline"+to_string(d); }

    static gsTHBSpline<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTHBSpline<d,T> >(node);
    }

    static gsXmlNode * put (const gsTHBSpline<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsTHBSpline<d,T> >(obj,data);
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
template<unsigned d, class T>
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

/// Get a Hierarchical B-spline basis from XML data
template<unsigned d, class T>
class gsXml< gsHBSplineBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsHBSplineBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "HBSplineBasis"+ (d>1 ? to_string(d):""); }
    
    static gsHBSplineBasis<d,T> * get (gsXmlNode * node)
    {
        return getHTensorBasisFromXml< gsHBSplineBasis<d,T> > (node);
    }
  
    static gsXmlNode * put (const gsHBSplineBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putHTensorBasisToXml< gsHBSplineBasis<d,T> > (obj, data);
    }
};


/// Get a Truncated Hierarchical B-spline basis from XML data
template<unsigned d, class T>
class gsXml< gsTHBSplineBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTHBSplineBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "THBSplineBasis"+ (d>1 ? to_string(d):""); }

    static gsTHBSplineBasis<d,T> * get (gsXmlNode * node)
    {
        return getHTensorBasisFromXml< gsTHBSplineBasis<d,T> > (node);
    }

    static gsXmlNode * put (const gsTHBSplineBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putHTensorBasisToXml< gsTHBSplineBasis<d,T> > (obj, data);
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
            return gsXml< gsBSpline<T, gsKnotVector<T> > >::get(node);
        if ( s == "Nurbs"      ) 
            return gsXml< gsNurbs<T,   gsKnotVector<T> > >::get(node);
        if ( s == "HBSpline2"  )  
            return gsXml< gsHBSpline<2,T> >::get(node); 
        if ( s == "HBSpline3"  )  
            return gsXml< gsHBSpline<3,T> >::get(node); 
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

	    if ( const gsBSpline<T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsBSpline<T,gsKnotVector<T> > *>( ptr ) )
		    return gsXml< gsBSpline<T, gsKnotVector<T> > >::put(*g,data);
        
	    if ( const gsNurbs<T> * g = 
             dynamic_cast<const gsNurbs<T> *>( ptr ) )
		    return gsXml< gsNurbs<T> >::put(*g,data);
        
	    if ( const gsTensorBSpline<2,T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorBSpline<2,T,gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorBSpline<2,T,gsKnotVector<T> > >::put(*g,data);
        
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

	    if ( const gsTensorBSpline<2,T,gsCompactKnotVector<T> > * g = 
             dynamic_cast<const gsTensorBSpline<2,T,gsCompactKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorBSpline<2,T,gsCompactKnotVector<T> > >::put(*g,data);
        
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
            return gsXml< gsBSpline<T, gsKnotVector<T> > >::get(node);
        if ( s == "Nurbs"      ) 
            return gsXml< gsNurbs<T,   gsKnotVector<T> > >::get(node);
      
        gsWarn<<"gsXmlUtils: getCurve: No known curve \""<<s<<"\". Error.\n";
        return NULL;
    }


    static gsXmlNode * put (const gsCurve<T> & obj,
                            gsXmlTree & data)
	{
	    const gsGeometry<T> * ptr = & obj;

	    if ( const gsBSpline<T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsBSpline<T,gsKnotVector<T> > *>( ptr ) )
		    return gsXml< gsBSpline<T, gsKnotVector<T> > >::put(*g,data);
        
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
            return gsXml< gsBSplineBasis<T, gsKnotVector<T> > >::get(node);
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
            return gsXml< gsTensorBSplineBasis<2, T, gsKnotVector<T> > >::get(node);
        if ( s == "TensorBSplineBasis3" )
            return gsXml< gsTensorBSplineBasis<3, T, gsKnotVector<T> > >::get(node);
        if ( s == "TensorBSplineBasis4" )
            return gsXml< gsTensorBSplineBasis<4, T, gsKnotVector<T> > >::get(node);

        if ( s == "TensorNurbsBasis2" )
            return gsXml< gsTensorNurbsBasis<2, T, gsKnotVector<T> > >::get(node);
        if ( s == "TensorNurbsBasis3" )
            return gsXml< gsTensorNurbsBasis<3, T, gsKnotVector<T> > >::get(node);
        if ( s == "TensorNurbsBasis4" )
            return gsXml< gsTensorNurbsBasis<4, T, gsKnotVector<T> > >::get(node);
               
        gsWarn<<"gsXmlUtils: getBasis: No known basis \""<<s<<"\". Error.\n";
        return NULL;
    }

    static gsXmlNode * put (const gsBasis<T> & obj,
                            gsXmlTree & data )
    {
        const gsBasis<T> * ptr = & obj;

        if ( const gsBSplineBasis<T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsBSplineBasis<T,gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsBSplineBasis<T, gsKnotVector<T> > >::put(*g,data);

        if ( const gsNurbsBasis<T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsNurbsBasis<T,gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsNurbsBasis<T, gsKnotVector<T> > >::put(*g,data);

        // Tensor B-spline
        if ( const gsTensorBSplineBasis<2, T, gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorBSplineBasis<2, T, gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<2, T, gsKnotVector<T> > >::put(*g,data);

        if ( const gsTensorBSplineBasis<3, T, gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorBSplineBasis<3, T, gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<3, T, gsKnotVector<T> > >::put(*g,data);

        if ( const gsTensorBSplineBasis<4, T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorBSplineBasis<4, T,gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<4, T,gsKnotVector<T> > >::put(*g,data);

        // Tensor Nurbs
        if ( const gsTensorNurbsBasis<2, T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorNurbsBasis<2, T, gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorNurbsBasis<2, T, gsKnotVector<T> > >::put(*g,data);

        if ( const gsTensorNurbsBasis<3, T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorNurbsBasis<3, T, gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorNurbsBasis<3, T, gsKnotVector<T> > >::put(*g,data);

        if ( const gsTensorNurbsBasis<4, T,gsKnotVector<T> > * g = 
             dynamic_cast<const gsTensorNurbsBasis<4, T, gsKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorNurbsBasis<4, T, gsKnotVector<T> > >::put(*g,data);


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

        // Tensor B-spline (compact knot-vector)
        if ( const gsTensorBSplineBasis<2, T,gsCompactKnotVector<T> > * g = 
             dynamic_cast<const gsTensorBSplineBasis<2, T,gsCompactKnotVector<T> > *>( ptr ) )
            return gsXml< gsTensorBSplineBasis<2, T,gsCompactKnotVector<T> > >::put(*g,data);
    
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
    
    static gsXmlNode * put (const gsPde<T> & obj,
                            gsXmlTree & data )
    {
        return NULL;
    }
};


/// Get a Multipatch
template<class T>
class gsXml< gsMultiPatch<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsMultiPatch<T>);
    static std::string tag () { return "MultiPatch"; }
    static std::string type () { return ""; }
    
    static gsMultiPatch<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"MultiPatch"), 
                      "Something went wrong. Expected Multipatch tag." );
        
        gsXmlNode * toplevel = node->parent();// the geometry patches should be siblings of node
        
        const int d = atoi( node->first_attribute("parDim")->value() );

        // temporaries for interface reading
        gsVector<index_t> dirMap(d);
        gsVector<bool>    dirOrient(d);

        gsXmlNode * tmp = node->first_node("patches");
        std::istringstream str ;
        str.str( tmp->value() );
        
        std::vector< gsGeometry<T> *> patches;
        std::map<int,int> ids;
        if ( ! strcmp( tmp->first_attribute("type")->value(),"id_range") )
        {
            int first, last;
            str >> std::ws >>  first >>  std::ws >> last >> std::ws ;
            for ( int i = first; i<=last; ++i )
            {
                patches.push_back( getById< gsGeometry<T> >( toplevel, i ) );
                ids[i] = i - first;
            }
        }
        else if ( ! strcmp( tmp->first_attribute("type")->value(),"id_index") )
        {
            int c = 0;
            for (int pindex; str >> pindex;)
            {
                patches.push_back( getById< gsGeometry<T> >( toplevel, pindex ) );
                ids[pindex] = c++;
            }
        }
        else
        {
            gsWarn<<"Unknown tag in XML multipatch object.\n";
        }
        
        std::string line;
        gsVector<int> p(4); // patch-side-patch-side
        // Read interfaces
        std::vector< boundaryInterface > interfaces;
        tmp = node->first_node("interfaces");
        if (tmp)
        {
            str.clear();
            str.str( tmp->value() );
            
            // Read interface (groups or size 4 + 2*d)
            
            while ( str>>std::ws >> p[0] ) // While there are more ints (groups or size 4+d-1)
            {
                for ( int i=1; i<4; ++i)
                    if ( ! (str >> std::ws >> p[i] >> std::ws) )
                        gsWarn<<"Error reading interface.\n";
                
                // Get ids
                p[0] = ids[ p[0] ];
                p[2] = ids[ p[2] ];

// /*
                // Read the matching direction permutation
                for ( int i=0; i!=d; ++i)
                    if ( !(str >> std::ws >> dirMap[i]) )
                        gsWarn<<"Error reading interface direction map.\n";

                // Read the interface orientation
                for ( int i=0; i!=d; ++i)
                    if ( !(str >> std::ws >> dirOrient[i]) )
                        gsWarn<<"Error reading interface orientation.\n";

                interfaces.push_back( boundaryInterface(p, dirMap, dirOrient) );
//*/

/*           // OLD format: read in Orientation flags
                gsVector<bool> orient(d-1);// orientation flags
                int k;
                for ( int i=0; i!=d-1; ++i)
                {
                    if ( !(str >> std::ws >> k >> std::ws) )
                        gsWarn<<"Error reading interface orientation.\n";
                    orient[i]= (k>0);
                }
                interfaces.push_back( boundaryInterface(p,orient) ) ;
//*/
            }
        }
        
        // Read boundary
        std::vector< patchSide > boundaries;
        tmp = node->first_node("boundary");
        if (tmp)
        {
            str.clear();
            str.str( tmp->value() );
            while ( str>>std::ws >> p[0]  )
            {
                p[0] = ids[ p[0] ];
                str >> std::ws >> p[1] ;
                boundaries.push_back( patchSide(p[0], p[1]) );
            }
        }
        
        return new gsMultiPatch<T>(patches, boundaries, interfaces);
    }
    
    static gsXmlNode * put (const gsMultiPatch<T> & obj,
                            gsXmlTree & data)
    {
        gsXmlNode * toplevel = data.first_node("xml");
        
        // First insert all geometries
        int id_start = 0;
        int id_end   = id_start;   
        gsXmlNode* tmp;
        for ( typename gsMultiPatch<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            tmp = gsXml<gsGeometry<T> >::put(**it,data);
            tmp->append_attribute( internal::makeAttribute("id", id_end++, data) );
            toplevel->append_node(tmp);
        }  
        
        //value : id_start id_end
        std::ostringstream str;
        str<<  id_start<<" "<< --id_end;
        tmp = internal::makeNode("patches" , str.str(), data);
        tmp->append_attribute( internal::makeAttribute("type", "id_range", data) );
        str.clear(); str.str("");
        
        // Make MultiPatch node
        gsXmlNode * mp_node = internal::makeNode("MultiPatch" , data);
        mp_node->append_attribute( internal::makeAttribute("parDim", obj.parDim() , data) );
        mp_node->append_node(tmp);
      
        if ( obj.nInterfaces() != 0 )
        {      
            for ( typename gsMultiPatch<T>::const_iiterator it = obj.iBegin();
                  it != obj.iEnd(); ++it )
            {
                str<< it->first().patch  <<" " << int(it->first().side())<<" "
                   << it->second().patch <<" " << int(it->second().side())<<" "
                   << it->dirMap().transpose()         <<" "
                   << it->dirOrientation().transpose() <<"\n";
            }
            tmp = internal::makeNode("interfaces", str.str(),  data);
            mp_node->append_node(tmp);
            str.clear(); str.str("");
        }
        
        if ( obj.nBoundary() )
        {
            for ( typename gsMultiPatch<T>::const_biterator it = obj.bBegin();
                  it != obj.bEnd(); ++it )
            {
                str<< it->patch <<" "<< int(it->side())<<"\n";
            }
            tmp = internal::makeNode("boundary", str.str(), data);
            mp_node->append_node(tmp);
            str.clear(); str.str("");
        }
        
        return mp_node;
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

    static gsXmlNode * put (const gsCurveFitting<T> & obj,
                            gsXmlTree & data )
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
        unsigned d = atoi( node->first_attribute("dim")->value() );

        
        unsigned tDim = 0;
        gsXmlAttribute * targetDim = node->first_attribute("targetDim");
        
        if ( targetDim )
            tDim = atoi( targetDim->value() );

        if ( tDim >= 1 )
        {
            gsXmlNode * tmp = node->first_node("rhs");
            gsMFunctionExpr<T>  rhs_fnct;
            getFunctionFromXml(tmp, rhs_fnct);
            
            tmp = node->first_node("solution");
            if ( tmp )
            {
                gsMFunctionExpr<T> msol;
                getFunctionFromXml(tmp, msol);
                
                return new gsPoissonPde<T>(rhs_fnct, d, msol );
            }
            else
            {
                return new gsPoissonPde<T>( rhs_fnct, d );
            }
        }

        // Read right hand side function
        gsXmlNode   * tmp = node->first_node("rhs");	
        gsFunctionExpr<T> rhs(tmp->value());

        // Read exact solution, if one exists in the file
        tmp = node->first_node("solution");	
        if ( tmp )
        {
            gsFunctionExpr<T> sol(tmp->value());
            //gsDebugVar (*sol);
            return new gsPoissonPde<T>(rhs, d, sol );
        }
        else
        {
            return new gsPoissonPde<T>( rhs, d );
        }
    }
    
    static gsXmlNode * put (const gsPoissonPde<T> & obj,
                            gsXmlTree & data )
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
                for (int side; str >> side;) 
                    bvp->addCondition( static_cast<boxSide>(side),
                                       condition_type::dirichlet, ff);
            }
            else if ( !strcmp(child->first_attribute("type")->value(), "neumann") )
            {		       
                for (int side; str >> side;) 
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

}// end namespace internal

}// end namespace gismo

//#undef GSXML_COMMON_FUNCTIONS
//#undef TMPLA2
//#undef TMPLA3
