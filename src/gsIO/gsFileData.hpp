/** @file gsFileData.hpp

    @brief Implementation of utility class which holds I/O XML data to
    read/write to/from files

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once


//#include <string>
//#include <iostream>
//#include <fstream>

#include <gsNurbs/gsKnotVector.h>

#include <rapidxml/rapidxml.hpp>       // External file
#include <rapidxml/rapidxml_print.hpp> // External file

#ifdef gsOpennurbs_ENABLED
#include <gsOpennurbs/gsReadOpenNurbs.h>
#endif

#ifdef gsOpenCascade_ENABLED
#include <gsOpenCascade/gsReadOcct.h>
#endif

#ifdef GISMO_WITH_PSOLID               // Extension files
#include <gsParasolid/gsReadParasolid.h>
#endif


#include <gzstream/gzstream.h>
#include <gsIO/gsFileManager.h>

namespace gismo {

template<class T>
gsFileData<T>::gsFileData()
{
    data = new FileData;
    data->makeRoot();
}

template<class T>
gsFileData<T>::gsFileData(String const & fn, bool recursive)
{
    data = new FileData;
    data->makeRoot();
    this->read(fn, recursive);
}

template<class T>
gsFileData<T>::~gsFileData()
{
    data->clear();
    delete data;
}


template<class T> void
gsFileData<T>::clear()
{
    data->clear();
    data->makeRoot(); // ready to re-use
}


template<class T>
std::ostream & gsFileData<T>::print(std::ostream &os) const
{
    //rapidxml::print_no_indenting
    os<< *data;
    return os;
}


template<class T> void
gsFileData<T>::dump(std::string const & fname)  const
{ save(fname); }


template<class T> void
gsFileData<T>::addComment(std::string const & message)
{
    gsXmlNode * comment = internal::makeComment(message, *data);
    data->appendToRoot(comment);
}

template<class T> void
gsFileData<T>::save(std::string const & fname, bool compress)  const
{
    gsXmlNode * comment = internal::makeComment("This file was created by G+Smo "
                                                GISMO_VERSION, *data);
    data->prepend_node(comment);

    gsXmlNode * declNode = data->allocate_node(rapidxml::node_type::node_declaration);
    declNode->append_attribute(data->allocate_attribute("version","1.0"));
    declNode->append_attribute(data->allocate_attribute("encoding","UTF-8"));
    data->prepend_node(declNode);

    if (compress)
    {
        saveCompressed(fname);
        return;
    }

    String tmp = gsFileManager::getExtension(fname);
    if (tmp != "xml" )
        tmp = fname + ".xml";
    else
        tmp = fname;

    m_lastPath = tmp;

    std::ofstream fn( tmp.c_str() );
    //rapidxml::print_no_indenting
    fn<< *data;
    fn.close();
    data->remove_node( data->first_node() );
}

template<class T> void
gsFileData<T>::saveCompressed(std::string const & fname)  const
{
    String tmp = gsFileManager::getExtension(fname);
    if (tmp != "gz" )
    {
        if (tmp != "xml" )
            tmp = fname + ".xml.gz";
        else
            tmp = fname + ".gz";
    }
    else
        tmp = fname;

    m_lastPath = tmp;

    ogzstream fn( tmp.c_str() );
    fn << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    //rapidxml::print_no_indenting
    fn<< *data;
    fn.close();
}

template<class T> void
gsFileData<T>::ioError(int lineNumber, const std::string& str)
{
    gsWarn<<"gsFileData: Problem with file "<<m_lastPath
          <<": IO error near line "<<lineNumber<<std::endl;
    throw std::runtime_error(str + " failed");
}

template<class T>
bool gsFileData<T>::read(String const & fn, bool recursive)
{
    m_lastPath = gsFileManager::find(fn);
    if ( m_lastPath.empty() )
    {
        gsWarn<<"gsFileData: Problem with file "<<fn<<": File not found.\n";
        gsWarn<<"search paths: "<< gsFileManager::getSearchPaths()<<"\n";
        return false;
    }

    // Identify filetype by extension
    String ext = gsFileManager::getExtension(fn);

    if (ext== "xml")
        return readXmlFile(m_lastPath, recursive);
    else if (ext== "gz" && util::ends_with(m_lastPath, ".xml.gz") )
        return readXmlGzFile(m_lastPath, recursive);
    else if (ext== "txt")
        return readGeompFile(m_lastPath);
    else if (ext== "g2")
        return readGoToolsFile(m_lastPath);
    else if (ext== "axl")
        return readAxelFile(m_lastPath);
    else if (ext== "off")
        return readOffFile(m_lastPath);
#ifdef gsOpennurbs_ENABLED
    else if (ext== "3dm")
        return read3dmFile(m_lastPath);
#endif
#ifdef gsOpenCascade_ENABLED
    else if (ext== "brep")
        return readBrepFile(m_lastPath);
    //else if (ext== "iges")
    //    return readIgesFile(m_lastPath);
    //else if (ext== "step")
    //    return readStepFile(m_lastPath);
#endif
#ifdef GISMO_WITH_PSOLID
    else if (ext== "xmt_txt")
        return readParasolidFile(m_lastPath);
    else if (ext== "x_t")
        return readParasolidFile(m_lastPath);
    else if (ext== "xmt_bin")
        return readParasolidFile(m_lastPath);
#endif
    else if (ext== "obj")
        return readObjFile(m_lastPath);
    else if (ext== "stl")
        return readStlFile(m_lastPath);
    else if (ext=="igs" || ext== "iges")
        return readIgesFile(m_lastPath);
//    else if (ext=="bv")
//        return readBezierView(m_lastPath);
    else if (ext=="x3d")
        return readX3dFile(m_lastPath);
    else if (ext=="csv")
        return readCsvFile(m_lastPath);
    else
    {
        gsWarn<<"gsFileData: Problem with file "<<fn<<": Unknown extension \"."<<ext<<"\".\n";
        return false;
    }
}

/*---------- Native Gismo format */

template<class T>
bool gsFileData<T>::readXmlFile( String const & fn, bool recursive)
{
    // Open file
    std::ifstream file(fn.c_str(), std::ios::in);
    if ( file.fail() )
    {gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    return readGismoXmlStream(file, recursive);
}

template<class T>
bool gsFileData<T>::readXmlGzFile( String const & fn, bool recursive)
{
    // Open file
    igzstream file(fn.c_str(), std::ios::in);
    if ( file.fail() )
    {gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    return readGismoXmlStream(file, recursive);
}


template<class T>
bool gsFileData<T>::readGismoXmlStream(std::istream & is, bool recursive)
{
    m_buffer.emplace_back();
    std::vector<char> & buffer = m_buffer.back();
    buffer.assign( std::istreambuf_iterator<char>(is.rdbuf() ),
                   std::istreambuf_iterator<char>() );
    buffer.push_back('\0');
    // Load file contents
    data->parse<0>(&buffer[0], true);

    gsXmlNode * ln = data->last_node("xml");
    if ( !ln )
    {
        gsWarn<< "gsFileData: Problem with file "<<m_lastPath
              <<": Invalid XML file, no root tag <xml> found.\n";
        assert( ln ) ;
    }

    if (recursive)
    {
        std::list<std::string> ifn;
        for (gsXmlNode * child = ln->first_node("xmlfile") ;
            child; child = child->next_sibling("xmlfile") )
        {
            ifn.push_back( child->value() );
            //ln->remove_node(child);
        }

        std::string cfn = gsFileManager::getPath(m_lastPath);
        for (auto & f : ifn)
        {
            readXmlFile(cfn + f);
        }
    }

    gsXmlNode * root = data->getRoot();
    root->merge_sibling(ln);
    data->remove_node(ln);

    // TO DO: Check if it contains unknown tags...
    return true;
}

template<class T>
void gsFileData<T>::getInclude(gsFileData<T> & res, index_t id, real_t time, std::string label)
{   
    // Ensures that only one argument is actually provided
    GISMO_ENSURE((id!=-1 ^  time!=-1. ^  label!="") &&
                !(id!=-1 && time!=-1. && label!=""),
                "gsFileData::getInclude("<<id<<","<<time<<","<<label<<"), too many arguments provided!");
    std::string attr_name, attr_string;

    // Attribute name and string representation, depending on pprovided input.
    if ( id!=-1 )
    {
        attr_name   = "id";
        attr_string = std::to_string(id); 
    }
    else if ( time!=-1)
    {
        attr_name   = "time";
        attr_string = std::to_string(time);
    }
    else if ( label!="")
    {
        attr_name   = "label";
        attr_string = label;
    } 

    bool found=false;
    gsXmlNode * root = getXmlRoot();
    const gsXmlAttribute * attribute;

    gsXmlNode * nd = internal::searchNode(root, attr_name, attr_string, "xmlfile");
    if (nd)
    {
        std::string filename = gsFileManager::getPath(m_lastPath) +  nd->value();
        res.clear();
        res.read(filename);
        return;
    }
    GISMO_ERROR("Include with " << attr_name << "=" << attr_string << " does not exist!");
}

/*---------- Axl file */

template<class T>
bool gsFileData<T>::readAxelFile( String const & fn )
{
    // Open file
    std::ifstream file(fn.c_str(), std::ios::in);
    if ( file.fail() )
    {gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    std::vector<char> buffer(
        std::istreambuf_iterator<char>(file.rdbuf() ),
        std::istreambuf_iterator<char>() );
    buffer.push_back('\0');

    // Read Axel Xml data
    FileData axldata;
    axldata.parse<0>(&buffer[0]);

    // Look for the root <axl>
    gsXmlNode * node = axldata.first_node("axl");
    String s;

    // Translate to Gismo XML
    for (gsXmlNode * child = node->first_node(); child; child = child->next_sibling())
    {
        s = child->name();
        if ( s == "curve" )
        { readAxelCurve(child);  }
        if ( s == "surface" )
        { readAxelSurface(child); }
        if ( s == "mesh" )
        { readAxelMesh(child); }
    }

    buffer.clear();
    axldata.clear();
    return true;
};

template<class T>
bool gsFileData<T>::readAxelCurve(gsXmlNode * node )
{
    std::stringstream str;

    //bool rational(true);

    gsXmlNode* g = internal::makeNode("Geometry", *data);
    g->append_attribute( internal::makeAttribute("type", "BSpline", *data) );
    data->appendToRoot(g);
    gsXmlNode * parent= g;

    gsXmlNode* tmp = node->first_node("dimension");
    String geoDim = tmp->value();

    gsXmlNode* b = internal::makeNode("Basis", *data);
    b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
    parent->append_node(b);

    unsigned d;
    tmp= node->first_node("order");
    str.str( tmp->value() );
    str >> d; d-=1;

    tmp = node->first_node("knots");
    g = internal::makeNode("KnotVector", tmp->value(), *data);
    g->append_attribute( internal::makeAttribute("degree", d, *data ) );
    b->append_node(g);

    // Coefficients
    tmp = node->first_node("points");
    g = internal::makeNode("coefs", tmp->value(), *data);
    g->append_attribute( internal::makeAttribute("geoDim", geoDim, *data ) );
    parent->append_node(g);

    return true;
};

template<class T>
bool gsFileData<T>::readAxelSurface(gsXmlNode * node )
{
    std::stringstream str;

    gsXmlNode* g = internal::makeNode("Geometry", *data);
    g->append_attribute( internal::makeAttribute("type", "TensorBSpline2", *data) );
    data->appendToRoot(g);
    gsXmlNode * parent = g;

    //gsXmlNode* tmp = node->first_node("dimension");// dimension is 3

    unsigned d[2];
    gsXmlNode * tmp = node->first_node("order");
    str.clear();
    str.str( tmp->value() );
    str >> d[0] >> d[1];
    d[0]-=1; d[1]-=1;

    // Tensor Basis
    gsXmlNode* tb = internal::makeNode("Basis", *data);
    tb->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis2", *data) );
    //tb->append_attribute( internal::makeAttribute("parDim", "2", *data ) );
    parent->append_node(tb);

    gsXmlNode* b = internal::makeNode("Basis", *data);
    b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
    b->append_attribute( internal::makeAttribute("index", "0", *data) );
    tb->append_node(b);
    tmp = node->first_node("knots");
    g = internal::makeNode("KnotVector", tmp->value(), *data);
    g->append_attribute( internal::makeAttribute("degree", d[0], *data ) ) ;
    b->append_node(g);

    b = internal::makeNode("Basis", *data);
    b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
    b->append_attribute( internal::makeAttribute("index", "1", *data) );
    tb->append_node(b);
    tmp = tmp->next_sibling("knots");
    g = internal::makeNode("KnotVector", tmp->value(), *data);
    g->append_attribute( internal::makeAttribute("degree", d[1], *data ) ) ;
    b->append_node(g);

    // Coefficients
    tmp = node->first_node("points");
    g = internal::makeNode("coefs", tmp->value(), *data);
    g->append_attribute( internal::makeAttribute("geoDim", "3", *data ) );
    parent->append_node(g);

    return true;
};

template<class T>
bool gsFileData<T>::readAxelMesh(gsXmlNode * node )
{
    std::ostringstream str;
    gsXmlNode * tmp = node->first_node("points");
    GISMO_ASSERT(tmp,"No points in mesh.");
    str<< tmp->value();
    tmp = node->first_node("faces");
    if (tmp)
        str<< tmp->value();
    internal::gsXmlNode* g = internal::makeNode("Mesh", str.str(), *data);
    g->append_attribute(internal::makeAttribute("type", "off", *data ) );

    tmp = node->first_node("count");
    GISMO_ASSERT(tmp,"No point count in mesh.");
    std::istringstream iss(tmp->value());
    int n;
    iss >> n;
    g->append_attribute(internal::makeAttribute("vertices", n, *data ) );
    iss >> n; iss >> n; //skips one deliberately
    g->append_attribute( internal::makeAttribute("faces", n, *data ) );
    internal::gsXmlNode* parent = data->first_node("xml") ;
    parent->append_node(g);

    tmp =  node->first_node("edges");
    if (tmp)
        gsWarn <<"Skipping edges in axl file.\n";
    return true;
}

// ******************************** //
// GoTools g2 file
// ******************************** //

template<class T>
bool gsFileData<T>::readGoToolsFile( String const & fn )
{
    //Input file
    std::ifstream file(fn.c_str(),std::ios::in);
    if ( !file.good() )
    {gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    std::istringstream lnstream;
    lnstream.unsetf(std::ios_base::skipws);

    String line;

    // Node for a Geometry object
    gsXmlNode * g;
    // Node for a basis object
    gsXmlNode * src;

    // Temporaries
    bool rational;
    int  ncp, deg, c, parDim, geoDim;

    // Structure:
    // type, version
    // geoDim, rational
    // n_coefs_u order_u
    // knots_u
    // n_coefs_v order_v
    // knots_v
    // n_coefs_w order_w
    // knots_w
    // coefficients

    while (!file.eof() && getline(file, line))
    {
        while (!file.eof() && line == "") getline(file, line);

        // Read entity type and version line
        lnstream.clear();
        lnstream.str(line);
        if ( !(lnstream >> std::ws >> ncp )          ) continue;
        if ( !(lnstream >> std::ws >> c   ) || c > 1 ) continue;
        if ( !(lnstream >> std::ws >> c   ) || c >9  ) continue;
        if ( !(lnstream >> std::ws >> c   ) || c >9  ) continue;

        switch (ncp)
        {
        case 100:  // Class_SplineCurve
            parDim = 1;
            break;
        case 200:  // Class_SplineSurface
            parDim = 2;
            break;
        case 700:  // Class_SplineVolume
            parDim = 3;
            break;
        case 210:  // Class_BoundedSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools trimmed surface (ClassType="<<ncp<<") not implemented.\n";

        case 110:  // Class_CurveOnSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools CurveOnSurface (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 120:  // Class_Line
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Line (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 130:  // Class_Circle
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Circle (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 140:  // Class_Ellipse
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Ellipse (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 150:  // Class_BoundedCurve
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools BoundedCurve (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 160:  // Class_Hyperbola
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Hyperbola (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 170:  // Class_Parabola
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Parabola (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 211:  // Class_SurfaceOnVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools SurfaceOnVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 220:  // Class_GoBaryPolSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoBaryPolSurface (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 230:  // Class_GoHBSplineParamSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoHBSplineParamSurface (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 240:  // Class_CompositeSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools CompositeSurface (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 250:  // Class_Plane
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Plane (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 260:  // Class_Cylinder
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Cylinder (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 270:  // Class_Sphere
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Sphere (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 280:  // Class_Sphere
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Sphere (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 290:  // Class_Torus
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Torus (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 291:  // Class_SurfaceOfRevolution
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools SurfaceOfRevolution (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 292:  // Class_Disc
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Disc (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 293:  // Class_LRSplineSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools LRSplineSurface (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 294:  // Class_TSplineSurface
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools TSplineSurface (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 300:  // Class_Go3dsObject
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Go3dsObject (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 310:  // Class_GoHeTriang
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoHeTriang (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 320:  // Class_GoSdTriang
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoSdTriang (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 330:  // Class_GoQuadMesh
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoQuadMesh (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 340:  // Class_GoHybridMesh
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoHybridMesh (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 350:  // Class_ParamTriang
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoHybridMesh (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 360:  // Class_GoVrmlGeometry
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoVrmlGeometry (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 400:  // Class_PointCloud
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools PointCloud (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 410:  // Class_LineCloud
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools LineCloud (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 500:  // Class_GoTriangleSets
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools GoTriangleSets (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 510:  // Class_RectGrid
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools RectGrid (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 710:  // Class_BoundedVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools BoundedVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 720:  // Class_Parallelepiped
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools Parallelepiped (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 721:  // Class_SphereVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools SphereVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 722:  // Class_CylinderVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools CylinderVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 723:  // Class_ConeVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools ConeVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 724:  // Class_TorusVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools TorusVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        case 793:  // Class_LRSplineVolume
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Reading GoTools LRSplineVolume (ClassType="<<ncp<<") not implemented.\n";
            continue;
            break;
        default:
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Unknown GoTools entity (ClassType="<<ncp<<").\n";
            continue;
            break;
        }

        // Read Geometry dimension and type
        getline(file, line);
        while (!file.eof() && line == "") getline(file, line);
        lnstream.clear();
        lnstream.str(line);
        lnstream >> std::ws >>  geoDim >>  std::ws >> rational >> std::ws ;

        g = internal::makeNode("Geometry", *data);
        if (parDim==1)
	    {
            g->append_attribute( internal::makeAttribute("type",
                                                         (rational ? "Nurbs" : "BSpline"), *data) );
        }
        else
        {
            g->append_attribute( internal::makeAttribute("type",
                                                         (rational ? "TensorNurbs" : "TensorBSpline")+internal::to_string(parDim), *data) );
        }

        data->appendToRoot(g);
        src = internal::makeNode("Basis", *data);
        if (parDim>1)
            src->append_attribute( internal::makeAttribute("type",
                                                           "TensorBSplineBasis"+internal::to_string(parDim), *data) );

        if ( rational )
        {
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": RATIONAL GoTools input is not supported/tested/working.\n";
            // Rational tensor basis
            gsXmlNode* rtb = internal::makeNode("Basis", *data);
            rtb->append_attribute( internal::makeAttribute("type",
                                                           (parDim==1 ? "NurbsBasis" : "TensorNurbsBasis"+internal::to_string(parDim)), *data) );
            rtb->append_node(src);
            g->append_node(rtb);
        }
        else
        {
            g->append_node(src);
        }

        ncp = 1;
        for (int i=0; i<parDim; ++i)
		{
            // Get numCoeffs_i, order_i
            getline(file, line);
            while (!file.eof() && line == "") getline(file, line);
            lnstream.clear();
            lnstream.str(line);
            // Reading the degree
            lnstream >> std::ws >>  c >>  std::ws >> deg >> std::ws ;
            deg--;
            ncp *= c;

            // Reading a coordinate-wise basis (knot-vector)
            getline(file, line);
            while (!file.eof() && line == "") getline(file, line);
            lnstream.clear();
            lnstream.str(line);
            gsXmlNode* b;
            if (parDim > 1)
            {
                b = internal::makeNode("Basis", *data);
                b->append_attribute( internal::makeAttribute("index", i, *data) );
            }
            else
                b = src;

            b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
            gsXmlNode* k = internal::makeNode("KnotVector", lnstream.str(), *data);
            k->append_attribute( internal::makeAttribute("degree", deg, *data) ) ;
            b->append_node(k);
            if (parDim > 1)
                src->append_node(b);
        }

        // w, w*cp_x, w*cp_y, w*cp_z: coordinates of the weighted control points (rational)
        // otherwise
        // cp_x, cp_y, cp_z: coordinates of the control points (non-rational)
        // The control points are numbered in a reverse lexicographic order
        std::ostringstream coefstream;

        for (int i = 0; i < ncp; i++) // Assumes each coefficient on a new line
        {
            if ( getline(file, line) )
            {
                while (!file.eof() && line == "") getline(file, line);
                coefstream << line.substr(0,line.size()) << std::endl;
            }
            else
            {
                gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                      <<": Failed to read coefficients.\n";
                return false;
            }
        }
        src = internal::makeNode("coefs", *data);
        src->value( internal::makeValue( coefstream.str(), *data) );
        src->append_attribute( internal::makeAttribute("geoDim", geoDim, *data ) );
        g->append_node(src);
    }

    return true;
}

//template<class T>
//bool gsFileData<T>::readGoToolsSpline(gsXmlNode * node )
//{ }
//template<class T>
//bool gsFileData<T>::readGoToolsTrimSurf(gsXmlNode * node )
//{ }

/*---------- GeoPdes txt file */

template<class T>
bool gsFileData<T>::readGeompFile( String const & fn )
{
    //Input file
    std::ifstream file(fn.c_str(),std::ios::in);
    if ( file.fail() )
    { gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    std::istringstream lnstream;
    lnstream.unsetf(std::ios_base::skipws);

    std::stringstream str;

    String line;
    //int patch_count(0);

    int N,Np,Ni(0),Ns(0);
    std::vector<gsKnotVector<T> *> knots;
    T tmp;

    //Parsing file
    while (!file.eof() && getline(file, line))
        // ind = line.find_first_not_of(' ');
        if (line[0] != '#') break;

    if (file.eof())
    {
        gsWarn<<"gsFileData: Problem with file "<<fn<<": Reached end of file.\n";
        return false;
    }

    //  N : dimension of the geometry
    //  Np: number of patches to construct the geometry
    //  Ni: total number of interfaces, each one connecting two patches
    //  Ns: total number of subdomains, formed by the union of patches
    lnstream.clear();
    lnstream.str(line);
    lnstream >> std::ws >>  N >>  std::ws >> Np >> std::ws ;
    if (Np>1)
        lnstream >> std::ws >> Ni  >> std::ws >> Ns >>  std::ws ;

    // Start id by 1, to match numbering in GeoPDEs file:
    //max_id=0;

    gsXmlNode* g;

    //gsDebug<<"Reading N="<<N<<" and Np="<< Np  <<"\n";
    gsVector<int> p(N);
    gsVector<int> ncp(N);
    bool patch(true);
    String bdr, ifc;

    while ( !file.eof() ) // Read next line
    {
        while (!file.eof() && getline(file, line))
            if ( line[0] != '#' ) break;

        // lnstream.clear();
        // lnstream.str(line);
        // if ( ! isdigit(line[0]) )
        // 	lnstream >> token;

        std::transform(line.begin(),line.end(),line.begin(),::tolower);
        //gsDebug<< "token=\""<<token<<"\"\n";

        if ( line == "" )// avoid empty lines
        {
            continue;
        }
        /* //Note: no need to read topology
           else if ( line.find("interface")!=String::npos )
           {
           while (!file.eof() && getline(file, line))
           if (line[0] != '#') break;
           ifc.append(line);
           while (!file.eof() && getline(file, line))
           if (line[0] != '#') break;
           ifc.append(" ");
           ifc.append(line);
           while (!file.eof() && getline(file, line))
           if (line[0] != '#') break;
           ifc.append(" ");
           ifc.append(line);
           ifc.append("\n");
           patch=false;
           }
           else if ( line.find("boundary")!=String::npos )
           {
           while (!file.eof() && getline(file, line))
           if (line[0] != '#') break;
           // getting number of sides
           int nb;
           lnstream.clear();
           lnstream.str(line) ;
           lnstream >> std::ws >> nb;
           for ( int i = 0; i<nb; ++i)
           {
           while (!file.eof() && getline(file, line))
           if (line[0] != '#') break;
           bdr.append(line);
           bdr.append(" ");
           patch=false;
           }
           }
        */
        else if ( ( line.find("patch")!=String::npos ) || patch==true )
        {
            // gsDebug <<"Patch "<<  line <<"\n";
            //GISMO_ASSERT( patch_count++ < Np, "Something went wrong while reading GeoPDEs file." );

            // p(i): the degree in each Cartesian direction (N integers)
            if ( ! isdigit(line[0]) )
                while (!file.eof() && getline(file, line))
                    if (line[0] != '#') break;
            lnstream.clear();
            lnstream.str(line) ;
            for (int i=0;i<N;++i)
                lnstream >> std::ws >> p[i] ;
            //gsDebug<<"Reading degrees OK "<< p.transpose() <<"\n";

            // ncp(i): the number of control points in each direction (N integers)
            while (!file.eof() && getline(file, line))
                if (line[0] != '#') break;
            lnstream.clear();
            lnstream.str(line);
            for (int i=0;i<N;++i)
                lnstream >> std::ws >> ncp[i] ;
            unsigned sz= ncp.prod() ;
            //gsDebug<<"Reading ncps OK "<< ncp.transpose() <<"\n";

            g = internal::makeNode("Geometry", *data);
            g->append_attribute( internal::makeAttribute("type", "TensorNurbs"+internal::to_string(N), *data) );
            data->appendToRoot(g);

            // Rational tensor basis
            gsXmlNode* rtb = internal::makeNode("Basis", *data);
            rtb->append_attribute( internal::makeAttribute("type", "TensorNurbsBasis"+internal::to_string(N), *data) );
            g->append_node(rtb);

            // Read source basis
            gsXmlNode* src = internal::makeNode("Basis", *data);
            rtb->append_node(src);

            if (N==1)
            {
                src->append_attribute( internal::makeAttribute("type", "NurbsBasis", *data) );
                while (!file.eof() && getline(file, line))
                    if (line[0] != '#') break;
                lnstream.clear();
                lnstream.str(line);
                gsXmlNode* k = internal::makeNode("KnotVector", lnstream.str(), *data);
                k->append_attribute( internal::makeAttribute("degree", p[0], *data ) ) ;
                src->append_node(k);
            }
            else // N>1
            {
                src->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis"+internal::to_string(N), *data) );
                //src->append_attribute( internal::makeAttribute("parDim", N, *data ) );
                for (int i=0;i<N;++i)
                {
                    while (!file.eof() && getline(file, line))
                        if (line[0] != '#') break;
                    lnstream.clear();
                    lnstream.str(line);
                    gsXmlNode* b = internal::makeNode("Basis", *data);
                    b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
                    b->append_attribute( internal::makeAttribute("index", i, *data) );
                    src->append_node(b);
                    gsXmlNode* k = internal::makeNode("KnotVector", lnstream.str(), *data);
                    k->append_attribute( internal::makeAttribute("degree", p[i], *data) ) ;
                    b->append_node(k);
                }
            }

            // cp_x, cp_y, cp_z: coordinates of the weighted control points
            //   (see Section 4.2 of The NURBS Book, L. Piegl & W. Tiller)
            //   (N rows, each one with prod_{i=1}^{N} ncp(i) float values)
            // The control points are numbered in a reverse lexicographic order
            gsMatrix<T> coefs(sz,N ) ;
            for (int i=0;i<N;++i)
            {
                while (!file.eof() && getline(file, line))
                    if (line[0] != '#') break;
                lnstream.clear();
                lnstream.str(line);
                for (unsigned k=0;k< sz;++k)
                {
                    lnstream >> std::ws >> tmp ;
                    coefs(k,i) = tmp ;
                }
            }
            //gsDebug<<"Reading coefs OK\n"<< *coefs <<"\n";

            // weights: weight associated to each basis function (or control point)
            //          (prod(ncp ) float values)
            gsMatrix<T>  weights(sz,1) ;
            while (!file.eof() && getline(file, line))
                if (line[0] != '#') break;
            lnstream.clear();
            lnstream.str(line);
            for (unsigned k=0;k< sz;++k)
            {
                lnstream >> std::ws >> tmp ;
                weights(k,0) = tmp ;
                coefs.row(k) /= tmp; //Divide weighted coeffient by the weight
            }
            //gsDebug<<"Reading weights OK\n"<< *weights <<"\n";

            // if ( weights == gsMatrix<T>::Ones(sz,1) )
            //      gsDebug<<"gsFileData: In fact weights are all equal to 1.\n";

            gsXmlNode* c = internal::makeNode("weights", weights, *data, true);
            rtb->append_node(c);

            c = internal::makeNode("coefs", coefs, *data, false);
            c->append_attribute( internal::makeAttribute("geoDim", N, *data ) );
            g->append_node(c);
        }
    }

    /*
    // Note: no need to read multipatch structure
    if ( Np > 1 )
    {
    g = internal::makeNode("MultiPatch", *data);
    g->append_attribute( internal::makeAttribute("parDim",N, *data) );
    parent->append_node(g);
    str.clear();
    str << 1 <<" "<< Np;
    gsXmlNode* c = internal::makeNode("patches", str.str(), *data);
    c->append_attribute( internal::makeAttribute("type","id_range", *data) );
    g->append_node(c);

    c = internal::makeNode("interfaces", ifc, *data);
    g->append_node(c);
    c = internal::makeNode("boundary",bdr, *data);
    g->append_node(c);
    }
    */
    return true;
};

/*---------- SurfLab/BezierView */

/*
  template<class T>
  bool gsFileData<T>::readBezierView( String const & fn )
  {
  //Input file
  std::ifstream file(fn.c_str(),std::ios::in);
  if ( !file.good() )
  {gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

  std::istringstream lnstream;
  lnstream.unsetf(std::ios_base::skipws);

  String line;

  //Note: patch kind/type definitions
  // #define POLY   1        // polyhedron
  // #define TRIANG 3        // triangular patch
  // #define TP_EQ  4        // tensorproduct with same degree in both dir.
  // #define TP     5        // general tensorproduct
  // #define TRIM_CURVE 6
  // #define TP_BSP 7        // general b-spline tensorproduct
  // #define RATIONAL 8      probably rational tensor product Bezier patch
  // #define PNTRI    9      // PN triangle patch, containing points and normals
  // #define PNTP    10      // PN quads patch, containing points and normals
  int kind;

  char string[255];
  int degu, degv, num_normals = 0, num_points;

  gsXmlNode* parent = data->first_node("xml") ;
  // Node for a Geometry object
  gsXmlNode * g;
  // Node for a basis object
  gsXmlNode * src;

  while (!file.eof() )
  {
  // skip group ids
  while (!file.eof() && getline(file, line))
  if ( line.find("group")!=String::npos ||
  line.find("Group")!=String::npos )
  break;

  // get kind/type
  lnstream.clear();
  lnstream.str(line) ;
  lnstream >> std::ws >>  kind >>  std::ws ;

  switch (kind) :
  {
  case 4 :
  case 8 :
  case 10:
  case 5:

  // read degrees
  lnstream >> std::ws >> degu >>  std::ws ;
  if(kind==5)
  degv = degu;
  else
  lnstream >> std::ws >> degv >>  std::ws ;

  if(kind==10)
  {
  int Ndegu, Ndegv;
  lnstream >> std::ws >> Ndegu >>  std::ws >> Ndegv >> std::ws ;
  num_normals = 3 * (Ndegu+1)*(Ndegv+1);
  }

  num_points = (degu+1)*(degv+1);
  points_dim = (kind==8 ? 4 : 3);

  // read in all control points
  gsMatrix<T> coefs(num_points, points_dim);
  for (i=0;i<num_points;i++)
  {
  lnstream.clear();
  lnstream.str(line);
  for (int j=0;i<points_dim;++k)
  lnstream >> std::ws >> coefs(i,k) ;
  }

  if(kind==10)
  for (i=0;i<num_normals;i++)
  sstr.ignore(128, std::ws);

  if(kind==8) // 4D control points
  {
  gsMatrix<T>  weights =  coefs.row(3);
  coefs.resize(gsEigen::NoChange,3);
  gsDebug<<"weights: "<< weights.transpose() <<"\n";
  }

  g = internal::makeNode("Geometry", *data);
  src = internal::makeNode("Basis", *data);
  g->append_node(src);
  gsXmlNode* b = internal::makeNode("coefs", coefs, *data, true);
  g->append_node(b);

  src->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis2", *data) );

  String kv(4*degv+3,' ');
  for (i=0;i<degv+1;i++)
  [
  kv[2*i] = '0';
  }

  b = internal::makeNode("Basis", *data);
  b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
  b->append_attribute( internal::makeAttribute("index", 0, *data) );
  b = internal::makeNode("KnotVector", lnstream.str(), *data);
  b->append_attribute( internal::makeAttribute("degree", degu, *data) ) ;

  src->append_node(b);
  b = internal::makeNode("Basis", *data);
  b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
  b->append_attribute( internal::makeAttribute("index", 1, *data) );
  src->append_node(b);

  break;

  case 3 : // Triangular patch
  case 10:

  default:
  break;
  }
  }
  }
*/

/*---------- OFF mesh from .off file */

template<class T>
bool gsFileData<T>::readOffFile( String const & fn )
{
    //https://stackoverflow.com/questions/47125387/stringstream-and-binary-data
    //std::istringstream buffer;
    //buffer << file.rdbuf();

    //https://stackoverflow.com/questions/38874200/trying-to-replace-scanf-with-sstream

    /* //verb-read
    std::ifstream buffer(fn);
    std::ostringstream bb; bb << buffer.rdbuf();    
    gsXmlNode* m = internal::makeNode("SurfMesh", *data);
    m->append_attribute( internal::makeAttribute("type", "off", *data) );        
    m->value( internal::makeValue( bb.str(), *data) );
    data->appendToRoot(m);
    return true;
    */

    //Input file
    std::ifstream file(fn.c_str(),std::ios::in);
    if ( !file.good() )
    { gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    gsXmlNode* g = internal::makeNode("SurfMesh", *data);
    g->append_attribute( internal::makeAttribute("type", "off", *data) );
    data->appendToRoot(g);

    String line;
    std::istringstream lnstream;
    lnstream.unsetf(std::ios_base::skipws);
    std::ostringstream tmp;

    getline(file, line);
    if ( line.compare(0,3,"OFF") != 0)
        return false;

    getline(file, line);
    int nverts, nfaces, nedges(0);
    lnstream.str(line);
    lnstream >> std::ws >>  nverts >>
        std::ws >> nfaces >>
        std::ws >> nedges ;

    g->append_attribute( internal::makeAttribute("vertices", nverts, *data) );
    g->append_attribute( internal::makeAttribute("faces"   , nfaces, *data) );
    g->append_attribute( internal::makeAttribute("edges"   , nedges, *data) );

    for (int i = 0; i < nverts; i++)
        if ( getline(file, line) )
            tmp << line.substr(0,line.size()) << std::endl;
        else
            return false;

    for (int i = 0; i < nfaces; i++)
        if ( getline(file, line) )
            tmp << line.substr(0,line.size()) << std::endl;
        else
            return false;

    g->value( internal::makeValue( tmp.str(), *data) );
    tmp.clear();

    return true;
}

/*---------- STL mesh file */

template<class T>
bool gsFileData<T>::readStlFile( String const & fn )
{
    bool solid(false),facet(false),loop(false);
    //Input file
    std::ifstream file(fn.c_str(),std::ios::in);
    if ( !file.good() )
    { gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

    gsXmlNode* g = internal::makeNode("Mesh", *data);
    g->append_attribute( internal::makeAttribute("type", "off", *data) );
    data->appendToRoot(g);

    std::ostringstream triangles;
    triangles.unsetf(std::ios_base::skipws);
    std::ostringstream vertices;
    vertices.unsetf(std::ios_base::skipws);
    unsigned nvert(0), nfaces(0), tmp(0);

    unsigned lineNumber(0);
    std::string str;

    while( !file.eof() && getline(file, str) )
    {
        std::transform(str.begin(),str.end(),str.begin(),::tolower);
        if(str.find("solid")!=String::npos && str.find("endsolid")==String::npos)
        {
            if(solid) ioError(lineNumber,"startSolid");
            solid=true;
        }
        else if(str.find("endsolid")!=String::npos)
        {
            if(!solid || facet || loop) ioError(lineNumber,"endSolid");
            solid=false;
        }
        else if(str.find("facet")!=String::npos && str.find("endfacet")==String::npos)
        {
            if(!solid || facet || loop) ioError(lineNumber,"startFacet");
            facet=true;
        }
        else if(str.find("endfacet")!=String::npos)
        {
            if(!solid || !facet || loop) ioError(lineNumber,"endFacet");
            facet=false;
        }
        else if(str.find("outer")!=String::npos)
        {
            if(!solid || !facet || loop) ioError(lineNumber,"startLoop");
            loop=true;
        }
        else if(str.find("endloop")!=String::npos)
        {
            if(!solid || !facet || !loop )
                ioError(lineNumber,"endLoop");
            triangles<< tmp;
            for (unsigned i= nvert-tmp; i!=nvert; ++i)
                triangles<<" "<< i;
            triangles<<"\n";
            nfaces++;
            loop=false;
            tmp = 0;
        }
        else if(str.find("vertex")!=String::npos)
        {
            if(!solid || !facet || !loop )
                ioError(lineNumber,"vertex");
            tmp++;
            nvert++;
            size_t pos=str.rfind("vertex")+7;
            assert(pos!=String::npos);
            vertices << str.substr(pos, str.size()-pos) <<"\n";
        }
    }

    g->append_attribute( internal::makeAttribute("vertices", nvert,  *data) );
    g->append_attribute( internal::makeAttribute("faces"   , nfaces, *data) );
    vertices << triangles.str() ;
    g->value( internal::makeValue( vertices.str(), *data) );

    return true;
}


template<class T>
bool gsFileData<T>::readObjFile( String const & fn )
{
    GISMO_UNUSED(fn);
    //gsWarn<<"Assuming Linux file, please convert dos2unix first.\n";

    //Input file
    std::ifstream file(fn.c_str(),std::ios::in);
    if ( !file.good() )
    { gsWarn<<"gsFileData: Problem with file "<<fn<<": Cannot open file stream.\n"; return false; }

/*
// -- mesh start
    std::string s;
    while(file>>s)
   {
         switch(*s.c_str())
         {
         case 'v':
              {
                    vertex v;
                    file>>v.x>>v.y>>v.z;
                    this->vertices.push_back(v);
              }
              break;
         case 'f':
              {
                    face f;
                    file>>f.v1>>f.v2>>f.v3;
                    faces.push_back(f);
              }
              break;
         default:
             break;
        }
   }
//-- mesh end
*/
#if FALSE

    std::istringstream lnstream;
    lnstream.unsetf(std::ios_base::skipws);

    String token(" "), bdr, ifc;

    //Parsing file
    while (!file.eof() && getline(file, line))
    {
        lnstream.clear();
        lnstream.str(line);
        lnstream >> std::ws;
        if (lnstream.eof())
            continue;
        else if (lnstream.peek() == '#')
            continue;
        else {
            lnstream >> token;

            //vertex (v)
            if (token == "v")
            {
                double x;
                Vertex v;
                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> x >>  std::ws ;
                    v.push_back(x);
                }
                points.push_back(v);
            }
            //vertex (vp)
            else if (token == "vp")
            {
                double x;
                Vertex v;
                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> x >>  std::ws ;
                    v.push_back(x);
                }
                points2d.push_back(v);
            }
            else if (token == "cstype")
            {
                lnstream >>std::ws >> cstype >> std::ws;
                std::string t;
                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> t >> std::ws ;
                    cstype+= "_"+ t;
                    //gsDebug<<"File has a "<<cstype<<std::endl;
                }
            }
            else if (token == "deg")
            {
                lnstream >> std::ws >> du >> std::ws >> dv>> std::ws;
                gsDebug<<"du, dv: "<<du<<", "<<dv <<std::endl;
            }
            else if (token == "curv")
            {
                type="";
                cstype="";
                gsWarn<<"gsFileData: Problem with file "<<fn<<": Ignore curv.\n";

                continue;
            }
            else if (token == "curv2")
            {
                c2 = new curv2();
                c2->deg= du;
                int cp;
                type=token;

                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> cp ;
                    //gsDebug<<"control point2d_"<< cp<< " is "<< (points2d[cp-1])[0]  <<std::endl;
                    c2->control_points.push_back(cp-1);
                }
                //cstype="";

                continue;
            }
            else if (token == "surf")
            {
                surf_count++;
                int cp;
                type=token;
                lnstream >> std::ws >> start_u >> std::ws >> end_u;
                lnstream >> std::ws >> start_v >> std::ws >> end_v;
                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> cp ;
                    //gsDebug<<"control point_"<< cp<< " is "<< (points[cp-1])[0]  <<std::endl;
                    control_points.push_back(cp-1);
                }
            }
            else if (token == "parm")
            {
                char c;
                double val;
                lnstream >> std::ws >> c >> std::ws;

                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> val >> std::ws;
                    //gsDebug<<"knot_"<<c<<": "<<val  <<std::endl;
                    if (type == "curv2")
                        c2->knots.push_back(val);
                    else
                        knots[var[c]].push_back(val);
                }

            }
            else if (token == "trim")// trim loop
            {
                double u,v;
                int cv;
                trim_loop t;

                while (lnstream && !lnstream.eof())
                {
                    lnstream >> std::ws >> u >> std::ws >> v >> std::ws >> cv >> std::ws;
                    t.push_back(cv);
                }

                gsWarn<<"gsFileData: Problem with file "<<fn<<": Ignore trim loop.\n";
            }
            else if (token == "end")//always in the end
            {
                gsDebug<<"End reading "<< cstype <<" "<<type<<std::endl;

                if (cstype=="bspline" && type=="surf")
                {
                    curves.push_back(c2);
                }


                if ( trims.size() )
                    axl<<"<surface type=\"trimmed\" number =\"1\">\n";

                //output spline
                if (cstype=="bspline" && type=="surf")
                {
                    axl<<"<surface type=\"bspline\" name=\"bs_"<<surf_count<<"\">\n";
                    axl<<"<dimension>3</dimension>\n";
                    axl<<"<number>"<<knots[0].size()-du-1 <<" "<<knots[1].size()-dv-1<<"</number>\n";
                    axl<<"<order>"<<du+1<<" "<<dv+1<<"</order>\n";
                    axl<<"<knots>";
                    for (int i=0; i<knots[0].size(); i++ )
                        axl<< knots[0][i]<< " ";
                    axl<<"</knots>\n";
                    axl<<"<knots>";
                    for (int i=0; i<knots[1].size(); i++ )
                        axl<< knots[1][i]<< " ";
                    axl<<"</knots>\n";
                    axl<<"<points>\n";
                    for (int i=0; i<control_points.size(); i++ )
                    {
                        for (int j=0; j<points[control_points[i]].size(); j++ )
                            axl<< points[control_points[i]][j] << " ";
                        axl << "\n";
                    }
                    axl<<"</points>\n";
                    axl<< "</surface>\n";
                    gsDebug<<"Got "<< cstype <<" "<<type<<" "<<surf_count<<std::endl;
                }
                else if (cstype=="rat_bspline" && type=="surf")
                {
                    axl<<"<surface type=\"bspline\" rational=\"1\" name=\"bs_"<<surf_count<<"\">\n";
                    axl<<"<dimension>3</dimension>\n";
                    axl<<"<number>"<<knots[0].size()-du-1 <<" "<<knots[1].size()-dv-1<<"</number>\n";
                    axl<<"<order>"<<du+1<<" "<<dv+1<<"</order>\n";
                    axl<<"<knots>";
                    for (int i=0; i<knots[0].size(); i++ )
                        axl<< knots[0][i]<< " ";
                    axl<<"</knots>\n";
                    axl<<"<knots>";
                    for (int i=0; i<knots[1].size(); i++ )
                        axl<< knots[1][i]<< " ";
                    axl<<"</knots>\n";
                    axl<<"<points>\n";
                    for (int i=0; i<control_points.size(); i++ )
                    {
                        // RATIONAL ?
                        for (int j=0; j<points[control_points[i]].size()-1; j++ )
                        {
                            //gsWarn<<"ok "<< j <<std::endl;
                            axl<< points[control_points[i]][j] << " ";
                        }
                        //if (j==4) axl<<1;
                        axl << "\n";
                    }
                    axl<<"</points>\n";
//


                    for (int i=0; i<trims.size(); i++ )
                    {
                        axl<< "<curveloop>\n";
                        for (int j=0; j<trims[i].size(); j++ )
                        {
//            write_curve( axl, trims[i][j] );
                        }
                        axl<< "</curveloop>\n";
                    }

//
                    axl<< "</surface>\n";

                    gsDebug<<"Got "<< cstype <<" "<<type<<" "<<surf_count<<std::endl;
                }
                else if (cstype=="bspline" && type=="curv2")
                {

                }
                else if (cstype!="")
                {
                    gsWarn<<"gsFileData: Problem with file "<<fn<<": Ignoring "<< cstype<<" "<<type <<std::endl;
                }
                //if (surf_count==5) break;
                //delete knots, degrees, control points
                knots.clear();
                knots.push_back(Vertex());
                knots.push_back(Vertex());
                control_points.clear();
                trims.clear();
                //control_points2d.clear();
                type="";
                cstype="";
            }
            else { // unknown token
                std::string message = "ignoring line ‘" + line + "’";
            }
        }
    }


#endif

    return true;
}


template<class T>
bool gsFileData<T>::readBrepFile( String const & fn )
{
#ifdef gsOpenCascade_ENABLED
    return extensions::gsReadBrep( fn.c_str(), *data);
#else
    GISMO_UNUSED(fn);
    return false;
#endif
}

namespace
{

#define MAXENTITY 1000    /* maximum number of entities recognised from IGES file */
#define FIELD_L 26  /*length of fixed length field in text file (sign, digits - double float) */

enum entity_j_name{
    E_TYPE= 0,
    PD_PTR= 1,
    PD_CNT= 2,
    ELAYER= 3
};

int read_iges_line(FILE *f, char *s)
{
    /* read filerow from iges*/
    int i;
    int c;
    for (i=0; i < 80; i++)
    {
        c=getc(f);
        if (c == EOF) return -1;
        s[i]=(char) c;
        if (i==0)
        { /* eat CR LF in line beginning and reset "i" to zero */
            switch (s[i])
            {
            case '\x0D': i--; break;
            case '\x0A': i--; break;
            }
        }
    }
    return 0;
}

void copy_field(char *s1, int poz, int x, char *s2)
{
    /* copy "x" chars from position "poz" from string "s1" to string "s2"*/
    int i=poz, j=0;
    while (j < x)
    {
        s2[j]=s1[i];
        i++;
        j++;
    }
}

void format_number(char *n_s)
{
    /* format number string */
    int n_len = strlen(n_s);
    /* regularize E,e,D,d -> e */
    for (int i=0; i < n_len; i++)
    {
        if(n_s[i] == 'E' || n_s[i] == 'e' || n_s[i] == 'D' || n_s[i] == 'd')
            n_s[i]='e';
    }
    /* eliminate e0 e+0 e-0 e+000 at the end
     * cycle from the end - on 'e', move the end, else jump out
     * */
    if(n_len > 2)
    {
        if(n_s[n_len - 1] == '0')
        {
            for (int i=n_len-1; i >= 0; i--)
            {
                switch (n_s[i])
                {
                case '0':
                case '+':
                case '-':
                    break;
                case 'e':
                    n_s[i]= '\0';
                    n_len= i;
                    goto endfor_fn;
                default : goto endfor_fn;
                }
            }
        endfor_fn:  /* jump out of FOR from inside of SWITCH */
            ;
        }
    }
}

void parse_d_entry(int entity[MAXENTITY][4], int & emark, int & entity_sum, char *ret1, char *ret2)
{
    /* Parse single entry (pair of lines) in D section. */
    char fld[9];                   /* for reading the value of array */
    int kod, pd_ptr, pd_count, layer, form;
    fld[8]='\0';
    /* Read:
     * entity type (A), form (B), pointer to PD line (C), PD count (D), layer (E)
     * and if the entity is recognised (406-3,106-2,116,126,128) load it into entity[][]
     * */
    copy_field(ret1,0,8,fld);  kod=     atoi(fld);    /* A */
    copy_field(ret2,32,8,fld); form=    atoi(fld);    /* B */
    copy_field(ret1,8,8,fld);  pd_ptr=  atoi(fld);    /* C */
    copy_field(ret2,24,8,fld); pd_count=atoi(fld);    /* D */
    copy_field(ret1,32,8,fld); layer=   atoi(fld);    /* E */
    GISMO_UNUSED(form);
    switch (kod)
    {
/*
  case 406:
  if (form != 3) break;                       // only layers (406-3)
  entity[emark][E_TYPE]=kod;
  entity[emark][PD_PTR]=pd_ptr;
  entity[emark][PD_CNT]=pd_count;
  entity[emark][ELAYER]=layer;
  entity_sum++;
  emark++;
  break;
  case 106:
  if (o_extract_pt == 0) break;
  if (form != 2) break;                       // only points (106-2)
  if (o_extract_unblanked == 1 && ret1[65] == '1') break; // only unblanked !
  entity[emark][E_TYPE]=kod;
  entity[emark][PD_PTR]=pd_ptr;
  entity[emark][PD_CNT]=pd_count;
  entity[emark][ELAYER]=layer;
  entity_sum++;
  emark++;
  break;
  case 116:
  if (o_extract_pt == 0) break;
  if (o_extract_unblanked == 1 && ret1[65] == '1') break;
  entity[emark][E_TYPE]=kod;
  entity[emark][PD_PTR]=pd_ptr;
  entity[emark][PD_CNT]=pd_count;
  entity[emark][ELAYER]=layer;
  entity_sum++;
  emark++;
  break;
*/
    case 126:
        entity[emark][E_TYPE]=kod;
        entity[emark][PD_PTR]=pd_ptr;
        entity[emark][PD_CNT]=pd_count;
        entity[emark][ELAYER]=layer;
        entity_sum++;
        emark++;
        break;
    case 128:
        //if (only_extract_unblanked == 1 && ret1[65] == '1') break;
        entity[emark][E_TYPE]=kod;
        entity[emark][PD_PTR]=pd_ptr;
        entity[emark][PD_CNT]=pd_count;
        entity[emark][ELAYER]=layer;
        entity_sum++;
        emark++;
        break;
    }
}

void read_iges_pd128(char *s, int begin, std::stringstream & ss)
{
    /* read 128 PD in IGES */
    static char t[FIELD_L+1];
    static int phase; /* 0-header, 1-cps, 2-knots, 3-weights */
    static int m;       /* marker in strings x, y, z */
    static int c, endknot;    /* counter of points */
    static int q[5],b[4]; /* b[0,1]-maxindex_of_ctrl_pts, b[2,3]-degree_of_Bspl */
    //static int sw[5]; /* breaks: 10, +b[0]+b[2]+2, +b[1]+b[3]+2, +(b[0]+1)*(b[1]+1), +(b[0]+1)*(b[1]+1) */
    int i;  /* i-order of the character */
    static int j;  /* j-order of head. */
    if (begin == 1)
    {
        phase=0;
        c=0;
        m=0;
        j=0;
        //sw[0]=10;
    }
    for (i=0; i <= 64; i++)
    {
        switch (phase)
        {
        case 0:  // PHASE header
            switch (j)
            {
            case 1:
            case 2:
            case 3:
            case 4:
                if (s[i] == ',')
                {
                    t[m]='\0';
                    b[j-1]=atoi(t);
                    j++; m=0;
                }
                else
                { t[m]=s[i]; m++; }
                break;
            default: //j=0, 5..9
                if (s[i] == ',')
                {
                    // j=0: 128
                    //(j=1..4: data in b[] ) 
                    // 5 integer parameters:
                    // j=5: closed in first direction
                    // j=6: closed in second direction
                    // j=7: rational
                    // j=8: periodic in first direction
                    // j=9: periodic in second direction

                    // b[0]+b[2]+2  u-knots
                    // b[1]+b[2]+2  v-knots
                    // (b[0]+1)*(b[1]+1) weights
                    // control points
                    t[m]= '\0';
                    m= 0;
                    //printf("%d: Got number: %s\n", j, t);
                    if (0!=j) q[j-5]= atoi(t);
                    if (7==j) ss<<q[0]<<",";
                    j++;
                }
                else
                {                                  /* digit + - E e D d */
                    t[m]= s[i];
                    m++;
                    break;
                }
                if ( 10 == j )
                {
                    endknot = b[0]+b[2]+2+b[1]+b[3]+2;   
                    phase=1; // knots start
                }
            }
            break;
        case 3:  // PHASE CP
            if (s[i] == ' ') return;
            else if (s[i] == ',' || s[i] == ';')  /* the last coordinate */
            {
                t[m]='\0';
                format_number(t);
                ss <<" "<< t;
                m= 0;
                c++;
                if (c == 3*(b[0]+1)*(b[1]+1))
                {
                    c=0;
                    return;
                }
            }
            else
            {                                  /* digit + - E e D d */
                t[m]= s[i];
                m++;
            }
            break;
        case 1:  // PHASE knots
            if (s[i] == ' ') return;
            else if (s[i] == ',')
            {                 /* next knot */
                t[m]= '\0';
                format_number(t);
                ss <<" "<< t;
                m= 0;
                c++;
                if (c == b[0]+b[2]+2) ss<<","<<b[2]<<",";
            }
            else {                                  /* digit + - E e D d */
                t[m]= s[i];
                m++;
                break;
            }
        
            if (c == endknot)
            {
                ss<<","<<b[3]<<",";
                c=0;
                phase=2; //weights start
            }
            break;
        case 2:  // PHASE weights
            if (s[i] == ' ') return;
            else if (s[i] == ',')
            {                 /* next knot */
                t[m]= '\0';
                format_number(t);
                ss <<" "<< t;
                m= 0;
                c++;
            }
            else {                                  /* digit + - E e D d */
                t[m]= s[i];
                m++;
                break;
            }

            if (c == (b[0]+1)*(b[1]+1))
            {
                ss<<",";
                c=0;
                phase=3; //control points start
            }
            break;
        }
    }
}

void read_iges_pd126(char *s, int begin, std::stringstream & ss)
{
    /* read 126 PD in IGES */
    static char t[FIELD_L+1];
    static int phase; /* 0-header, 1-cps, 2-knots, 3-weights */
    static int m;       /* marker in strings x, y, z */
    static int c, endknot;    /* counter of points */
    static int b[2]; /* b[0]-maxindex_of_ctrl_pts, b[1]-degree_of_Bspl */
    int i;  /* i-order of the character */
    static int j;  /* j-order of head. */
    if (begin == 1)
    {
        phase=0;
        c=0;
        m=0;
        j=0;
    }
    for (i=0; i <= 64; i++)
    {
        switch (phase)
        {
        case 0:  // PHASE header
            switch (j)
            {
            case 1:
            case 2:
                if (s[i] == ',')
                {
                    t[m]='\0';
                    b[j-1]=atoi(t);
                    j++; m=0;
                }
                else
                { t[m]=s[i]; m++; }
                break;
            default:  //j=0, 3..6
                if (s[i] == ',')
                {
                    // j=0: 126
                    //(j=1..2: data in b[] ) 
                    // j=3: planar
                    // j=4: open/closed curve
                    // j=5: rational
                    // j=6: periodic
                    // b[0]+b[1]+2  knots
                    // b[0]+1 weights
                    // control points
                    t[m]= '\0';
                    m= 0;
                    //printf("%d: Got number: %s\n", j, t);
                    j++;
                }
                else
                {                                  /* digit + - E e D d */
                    t[m]= s[i];
                    m++;
                    break;
                }
                if ( 7 == j )
                {
                    endknot = b[0]+b[1]+2;
                    phase=1; // knots start
                }
            }
            break;
        case 3:  // PHASE CP
            if (s[i] == ' ') return;
            else if (s[i] == ',' || s[i] == ';')  /* the last coordinate */
            {
                t[m]='\0';
                format_number(t);
                ss <<" "<< t;
                m= 0;
                c++;
                if (c == 3*(b[0]+1))
                {
                    c=0;
                    return;
                }
            }
            else
            {                                  /* digit + - E e D d */
                t[m]= s[i];
                m++;
            }
            break;
        case 1:  // PHASE knots
            if (s[i] == ' ') return;
            else if (s[i] == ',')
            {                 /* next knot */
                t[m]= '\0';
                format_number(t);
                ss <<" "<< t;
                m= 0;
                c++;
            }
            else {                                  /* digit + - E e D d */
                t[m]= s[i];
                m++;
                break;
            }

            if (c == endknot)
            {
                ss<<","<<b[1]<<",";
                c=0;
                phase=2; //weights start
            }
            break;
        case 2:  // PHASE weights
            if (s[i] == ' ') return;
            else if (s[i] == ',')
            {                 /* next knot */
                t[m]= '\0';
                format_number(t);
                ss <<" "<< t;
                m= 0;
                c++;
            }
            else {                                  /* digit + - E e D d */
                t[m]= s[i];
                m++;
                break;
            }

            if (c == (b[0]+1))
            {
                ss<<",";
                c=0;
                phase=3; //control points start
            }
            break;
        }
    }
}


// todo: Write file


}//namespace

template<class T>
bool gsFileData<T>::readIgesFile( String const & fn )
{
    //Input file
    FILE * fr = fopen(fn.c_str(), "rb");

    int entity[MAXENTITY][4];
    int entity_sum=0,         /* A number of entities in array entity[] */
    emark=0;              /* for cycling over entity[] */

    if ( NULL==fr )
    {
        gsWarn<<"gsFileData: Problem with input "<<fn<<": Cannot open file.\n";
        return false;
    }

    char line[81],          /* text content of current line from IGES */
        pairline[81];       /* pair line (in D section of IGES file) */
    char pd_seq_s[8]; pd_seq_s[7]='\0';
    int  pd_seq;           /* parameter data sequence number */
    int err_code;

    if (read_iges_line(fr, line) == -1)
    { gsWarn<<"IGES file empty.\n"; return false; }
    if (line[72] == 'F')
    { gsWarn<<"IGES file is not ASCII file (binary mode is not supported).\n"; return false; }

    /* Skip S and G sections and find the 1st line of the 1st D entry,
     * then load 2nd line of the first D entry as well (2).
     * Analyze the pair of lines and if the enity is recognized, load it
     * into array entity[][] (3).
     * */
    while (line[72] != 'D')
    {
        err_code=read_iges_line(fr, line);
        if (err_code == -1) { gsWarn<<"IGES file empty.\n"; return false; }
    }
    err_code=read_iges_line(fr, pairline);             /* (2) */
    if (err_code == -1) { gsWarn<<"IGES file empty.\n"; return false; }
    parse_d_entry(entity,emark,entity_sum,line, pairline);                    /* (3) */

    /* Repeat the same on subsequent entries in D section */
    while (line[72] == 'D')
    {
        err_code=read_iges_line(fr, line);
        if (err_code == -1) { gsWarn<<"IGES file empty.\n"; return false; }
        if (line[72] == 'D')
        {
            err_code=read_iges_line(fr, pairline);
            if (err_code == -1) { gsWarn<<"IGES file empty.\n"; return false; }
            parse_d_entry(entity,emark,entity_sum,line, pairline);
            if (entity_sum >= MAXENTITY) { gsWarn<<"Too many entities (more than "<<MAXENTITY<<" in the file.\n"; return false; }
        }
    }

    /* For every entity in array entity[] find PD line with the sequence nr
     * and process subsequent lines (find coordinates and write to text file).
     * */
    copy_field(line,73,7,pd_seq_s);
    pd_seq = atoi(pd_seq_s);
    for (emark=0; emark < entity_sum; emark++)
    {
        while (pd_seq != entity[emark][PD_PTR]) /* find line in PD */
        {
            read_iges_line(fr, line);
            copy_field(line,73,7,pd_seq_s);
            pd_seq=atoi(pd_seq_s);
        }

        std::stringstream ss;
        std::string token;

        switch (entity[emark][E_TYPE])
        {
            /*
              case 406:
              for (i=1; i <= entity[emark][PD_CNT]; i++) {
              read_pd406(line, i);
              read_iges_line(fr, line);
              copy_field(line,73,7,pd_seq_s); pd_seq=atoi(pd_seq_s);
              }
              break;
              case 106:
              if (o_write_layername == 1) {
              fprintf(fw, "%s\n", layer[entity[emark][ELAYER]]);
              }
              fprintf(fw, "\n");
              for (i=1; i <= entity[emark][PD_CNT]; i++) {
              read_pd106_writetxt(line, i);
              read_iges_line(fr, line);
              copy_field(line,73,7,pd_seq_s); pd_seq=atoi(pd_seq_s);
              }
              fprintf(fw, "\n");
              break;
              case 116:
              for (i=1; i <= entity[emark][PD_CNT]; i++) {
              read_pd116_writetxt(line, i);
              read_iges_line(fr, line);
              copy_field(line,73,7,pd_seq_s); pd_seq=atoi(pd_seq_s);
              }
              break;
            */
        case 126:
            // layer name: layer[entity[emark][ELAYER]];
            for (int i=1; i <= entity[emark][PD_CNT]; i++)
            {
                read_iges_pd126(line, i, ss);
                read_iges_line(fr, line);
                copy_field(line,73,7,pd_seq_s);
                pd_seq=atoi(pd_seq_s);
            }
            break;
        case 128:
            for (int i=1; i <= entity[emark][PD_CNT]; i++)
            {
                read_iges_pd128(line, i, ss);
                read_iges_line(fr, line);
                copy_field(line,73,7,pd_seq_s);
                pd_seq=atoi(pd_seq_s);
            }

            //Debug
            //gsInfo<<"\nREAD\n"; while(std::getline(ss, token, ',')) gsInfo << token << '\n'; break;

            std::getline(ss, token, ',');
            bool rat = (token!="0");
            gsXmlNode* g = internal::makeNode("Geometry", *data);
            g->append_attribute( internal::makeAttribute("type",
                               (rat?"TensorNurbs2":"TensorBSpline2"), *data) );
            data->appendToRoot(g);
            gsXmlNode* rtb(nullptr), * src;
            src = internal::makeNode("Basis", *data);
            if (rat)
            {
                rtb = internal::makeNode("Basis", *data);
                rtb->append_attribute( internal::makeAttribute("type", "TensorNurbsBasis2", *data) );
                g->append_node(rtb);
                rtb->append_node(src);
            }
            else
                g->append_node(src);
            src->append_attribute( internal::makeAttribute("type","TensorBSplineBasis2", *data) );
            gsXmlNode* b = internal::makeNode("Basis", *data);
            b->append_attribute( internal::makeAttribute("index", 0, *data) );
            b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
            src->append_node(b);
            std::getline(ss, token, ',');
            gsXmlNode* k = internal::makeNode("KnotVector", token, *data);
            std::getline(ss, token, ',');
            k->append_attribute( internal::makeAttribute("degree", token, *data) ) ;
            b->append_node(k);
            b = internal::makeNode("Basis", *data);
            b->append_attribute( internal::makeAttribute("index", 1, *data) );
            b->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
            src->append_node(b);
            std::getline(ss, token, ',');
            k = internal::makeNode("KnotVector", token, *data);
            std::getline(ss, token, ',');
            k->append_attribute( internal::makeAttribute("degree", token, *data) ) ;
            b->append_node(k);
            std::getline(ss, token, ',');
            if (rat)
            {
                k = internal::makeNode("weights", token, *data);
                rtb->append_node(k);
            }
            std::getline(ss, token, ',');
            k = internal::makeNode("coefs", token, *data);
            k->append_attribute( internal::makeAttribute("geoDim", 3, *data ) );
            g->append_node(k);
            break;
        }
    }

    const bool fc = ( 0==fclose(fr) );
    if (fc) gsWarn<< "File closing didn't succeeded!\n";
    return fc;
}

#undef MAXENTITY
#undef FIELD_L

template<class T>
void gsFileData<T>::addX3dShape(gsXmlNode * shape)
{
    // assert shape->name()==Shape

    gsXmlNode * patch;

    //node = node->first_node("NurbsTrimmedSurface");
    //node = node->first_node("NurbsCurve2D");
    //node = node->first_node("NurbsCurve");

    int p;
    char * ch = 0;
    std::istringstream str;

    for (gsXmlNode * node = shape->first_node("NurbsPatchSurface");
         node; node = node->next_sibling("NurbsPatchSurface") )
    {
        // Read TensorBSplineBasis
        gsXmlNode* tp_node = internal::makeNode("Basis" , *data);
        tp_node->append_attribute( internal::makeAttribute("type",
                                                           "TensorBSplineBasis2", *data) );

        // gsDebug<<"node "<< node <<"\n";
        // gsDebug<<"uOrder "<< node->first_attribute("uOrder") <<"\n";
        // gsDebug<<"uKnot "<< node->first_attribute("uKnot")<<"\n";

        p  = atoi(node->first_attribute("uOrder")->value()) -1 ;
        gsXmlAttribute * kv_attr = node->first_attribute("uKnot");

        if ( kv_attr )
            ch = node->first_attribute("uKnot")->value() ;
        else
        {
            // make 0..1 knots by default
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Setting knots to [0..1] by default not implemented";
        }
        gsXmlNode * kv_node = internal::makeNode("KnotVector", String(ch), *data);
        kv_node->append_attribute( internal::makeAttribute("degree", p, *data) );
        gsXmlNode* bs_node = internal::makeNode("Basis" , *data);
        bs_node->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
        bs_node->append_attribute( internal::makeAttribute("index", 0, *data) );
        bs_node->append_node(kv_node);
        tp_node->append_node(bs_node);

        p  = atoi(node->first_attribute("vOrder")->value()) -1 ;
        kv_attr = node->first_attribute("vKnot");
        if ( kv_attr )
            ch = kv_attr->value() ;
        else
        {
            // make 0..1 knots by default
            gsWarn<<"gsFileData: Problem with file "<<m_lastPath
                  <<": Setting knots to [0..1] by default not implemented";
        }
        kv_node = internal::makeNode("KnotVector", String(ch), *data);
        kv_node->append_attribute( internal::makeAttribute("degree", p, *data) );
        bs_node = internal::makeNode("Basis" , *data);
        bs_node->append_attribute( internal::makeAttribute("type", "BSplineBasis", *data) );
        bs_node->append_attribute( internal::makeAttribute("index", 1, *data) );
        bs_node->append_node(kv_node);
        tp_node->append_node(bs_node);

        patch =  internal::makeNode("Geometry" , *data);
        gsXmlAttribute * rat_weights = node->first_attribute("weight");
        if ( rat_weights )
        {
            // If nurbs, read weights etc
            ch = rat_weights->value();
            patch->append_attribute(internal::makeAttribute("type",
                                                            "TensorNurbs2", *data));

            gsXmlNode* nurbs_node = internal::makeNode("Basis" , *data);
            nurbs_node->append_attribute( internal::makeAttribute("type",
                                                                  "TensorNurbsBasis2", *data) );
            nurbs_node->append_node(tp_node);
            gsXmlNode * weights_node = internal::makeNode("weights",String(ch), *data);
            nurbs_node->append_node(weights_node);
            patch->append_node(nurbs_node);
        }
        else
        {
            // Attach basis to the patch
            patch->append_attribute(internal::makeAttribute("type",
                                                            "TensorBSpline2", *data));
            patch->append_node( tp_node);
        }

        // Attach the control points to the patch
        ch = node->first_node("Coordinate")->first_attribute("point")->value();
        gsXmlNode* cp_node = internal::makeNode("coefs",String(ch), *data);
        cp_node->append_attribute( internal::makeAttribute("geoDim", 3, *data) );
        patch->append_node( cp_node);

        // Attach patch to gismo xml tree
        data->appendToRoot(patch);
    }
}

template<class T>
void gsFileData<T>::addX3dTransform(gsXmlNode * trans)
{

    gsXmlAttribute * attr = trans->first_attribute("translation");
    if ( attr )
        gsDebug<<"Translate "<< attr->value() <<"\n";// (x,y,z)


    attr = trans->first_attribute("rotation");
    if ( attr )
        gsDebug<<"Rotate "<< attr->value() <<"\n";// (x,y,z,angle)


    attr = trans->first_attribute("scale");
    if ( attr )
        gsDebug<<"Scale "<< attr->value() <<"\n";// (x,y,z)

//<transform dim="3">
// all children optional
//<translation>x,y,z</translation>
//<rotation>x y z angle</rotation>
//<scale>x y z</scale>
//<Matrix>---</Matrix>
//</transform>

}


template<class T>
bool gsFileData<T>::readX3dFile( String const & fn )
{
    // http://www.web3d.org/x3d/content/examples/NURBS/
    // Open file
    std::ifstream file(fn.c_str(), std::ios::in);
    if ( file.fail() )
    {gsWarn<<"gsFileData: Problem with file "<<fn<<": cannot open file stream.\n"; return false;}

    std::vector<char> buffer(
        std::istreambuf_iterator<char>(file.rdbuf() ),
        std::istreambuf_iterator<char>() );
    buffer.push_back('\0');

    // Read X3D data
    FileData x3ddata;
    x3ddata.parse<0>(&buffer[0]);

    // Look for the root <X3D>
    gsXmlNode * x3d = x3ddata.first_node("X3D");

    // Looking for shapes
    for (gsXmlNode * scene = x3d->first_node("Scene");
         scene; scene = scene->next_sibling("Scene") )
    {
        for (gsXmlNode * shape = scene->first_node("Shape");
             shape; shape = shape->next_sibling("Shape") )
            addX3dShape( shape );

        for (gsXmlNode * trans = scene->first_node("Transform");
             trans; trans = trans->next_sibling("Transform") )
        {
            // Descent to transforms
            gsXmlNode * trans_rec = trans;
            while ( true )
            {
                addX3dTransform( trans_rec );
                gsXmlNode * tmp = trans_rec->first_node("Transform");
                if ( tmp )
                    trans_rec = tmp;
                else
                    break;
            }

            for (gsXmlNode * shape = trans_rec->first_node("Shape");
                 shape; shape = shape->next_sibling("Shape") )
            {
                addX3dShape( shape );
            }

            for (gsXmlNode * coll  = trans->first_node("Collision");
                 coll; coll = trans->next_sibling("Collision") )
            {
                gsDebug<<"Reach collision tag.\n";
                for (gsXmlNode * shape = coll->first_node("Shape");
                     shape; shape = shape->next_sibling("Shape") )
                {
                    gsDebug<<"Reach shape in tag.\n";
                    addX3dShape( shape );
                }
            }
        }
    }

    return true;
}

template<class T>
bool gsFileData<T>::read3dmFile( String const & fn )
{
#ifdef gsOpennurbs_ENABLED
    return extensions::gsReadOpenNurbs( fn.c_str(), *data);
#else
    GISMO_UNUSED(fn);
    return false;
#endif
}


template<class T>
bool gsFileData<T>::readParasolidFile( String const & fn )
{
    // Remove extension and pass to parasolid
    //int lastindex = fn.find_last_of(".");
    //return extensions::gsReadParasolid( fn.substr(0, lastindex).c_str(), *data);
#ifdef GISMO_WITH_PSOLID
    return extensions::gsReadParasolid( fn.c_str(), *data);
#else
    GISMO_UNUSED(fn);
    return false;
#endif
}

template<class T>
bool gsFileData<T>::readCsvFile( String const & fn )
{
    std::ifstream indata;
    indata.open(fn.c_str());
    std::string cell, line, mstr;
    index_t rows = 0, nv = 0;
    std::istringstream lnstream;
    lnstream.unsetf(std::ios_base::skipws);
    while (std::getline(indata, line))
    {
        lnstream.clear();
        lnstream.str(line);
        while (std::getline(lnstream, cell, ','))
        {
            ++nv;
            mstr += cell + " ";
        }
        ++rows;
    }
    gsXmlNode * nd =  internal::makeNode("Matrix", mstr, *data);
    nd->append_attribute( internal::makeAttribute("format","ascii",*data) );
    nd->append_attribute( internal::makeAttribute("rows",rows,*data) );
    nd->append_attribute( internal::makeAttribute("cols",nv/rows,*data) );

    data->appendToRoot(nd);
    return true;
}

template<class T>
std::string
gsFileData<T>::contents () const
{
    std::ostringstream os;
    os << "--- \n";
    int i(1);
    for (gsXmlNode * child = data->first_node("xml")->first_node();
         child; child = child->next_sibling() )
    {
        os << i++ <<". " << child->name() ;
        for (gsXmlAttribute * attr = child->first_attribute();
             attr; attr = attr->next_attribute() )
            os << ", "<< attr->name() << "="<< attr->value();
        os <<"\n";
    }
    os << "--- \n";
    return os.str();
};

template<class T> inline
int gsFileData<T>::numTags() const
{
    int i(0);
    for (gsXmlNode * child = data->first_node("xml")->first_node() ;
         child; child = child->next_sibling() )
        ++i;
    return i;
}

template<class T> inline
typename gsFileData<T>::gsXmlNode *
gsFileData<T>::getXmlRoot() const
{
    return data->getRoot();
}

template<class T> inline
void gsFileData<T>::deleteXmlSubtree(gsXmlNode * node)
{
    node->parent()->remove_node(node);
    // TO do: delete recursively ?
    delete node;
}

template<class T> inline
typename gsFileData<T>::gsXmlNode *
gsFileData<T>::getFirstNode(const std::string & name, const std::string & type) const
{
    gsXmlNode * root = data->first_node("xml");
    if ( ! root )
    {
        gsWarn<< "gsFileData: Problem with file "<<m_lastPath
              <<": Invalid XML file, no root tag <xml> found.\n";
        assert( root ) ;
    }

    if ( type == "" )
        return root->first_node( name.c_str() );
    else
    {
        for (gsXmlNode * child = root->first_node( name.c_str() ) ;
             child; child = child->next_sibling( name.c_str() ) )
            if ( !strcmp( child->first_attribute("type")->value(), type.c_str() ) )
                return child;
        return NULL;
    }
}

template<class T> inline
typename gsFileData<T>::gsXmlNode *
gsFileData<T>::getAnyFirstNode(const std::string & name, const std::string & type) const
{
    gsXmlNode * root = data->first_node("xml");
    assert( root ) ;
    if ( type == "" )
        // Searching up to third level of the XML tree
        for (gsXmlNode * child = root->first_node() ;
             child; child = child->next_sibling() )
        {
            if (!strcmp( child->name(), name.c_str() ) )
                return child;
            // Level 2
            for (gsXmlNode * child2 = child->first_node() ;
                 child2; child2 = child2->next_sibling() )
            {
                if ( !strcmp( child2->name(), name.c_str() ) )
                    return child2;
                // Level 3
                for (gsXmlNode * child3 = child2->first_node() ;
                     child3; child3 = child3->next_sibling() )
                    if ( !strcmp( child3->name(), name.c_str() ) )
                        return child3;
            }
        }
    else
        // Searching up to third level of the XML tree
        for (gsXmlNode * child = root->first_node() ;
             child; child = child->next_sibling() )
        {
            if (!strcmp( child->name(), name.c_str() ) &&
                !strcmp( child->first_attribute("type")->value(), type.c_str() ) )
                return child;
            // Level 2
            for (gsXmlNode * child2 = child->first_node() ;
                 child2; child2 = child2->next_sibling() )
            {
                if ( !strcmp( child2->name(), name.c_str() ) &&
                     !strcmp( child2->first_attribute("type")->value(), type.c_str() ))
                    return child2;
                // Level 3
                for (gsXmlNode * child3 = child2->first_node() ;
                     child3; child3 = child3->next_sibling() )
                    if ( !strcmp( child3->name(), name.c_str() ) &&
                         !strcmp( child3->first_attribute("type")->value(), type.c_str()))
                        return child3;
            }
        }
    return NULL;
}

template<class T> inline
typename gsFileData<T>::gsXmlNode *
gsFileData<T>::getNextSibling(gsXmlNode* const & node, const std::string & name,
                              const std::string & type)
{
    if ( type == "" )
        return node->next_sibling( name.c_str() );
    else
    {
        for (gsXmlNode * next = node->next_sibling( name.c_str() );
             next; next = next->next_sibling( name.c_str() ) )
            if ( !strcmp( next->first_attribute("type")->value(), type.c_str() ) )
                return next;
        return NULL;
    }
}


};// namespace gismo
