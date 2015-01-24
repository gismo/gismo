/** @file gsReadParasolid.cpp

    @brief Provides implementation of gsReadParasolid functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsParasolid/gsFrustrum.h>

#include <gsParasolid/gsReadParasolid.h>

#include <sstream>
 #include <string>
#include <fstream>
#include <iomanip>

#include <rapidxml/rapidxml.hpp> // External file


namespace gismo {

namespace extensions {

bool gsReadParasolid( const char * fname, internal::gsXmlTree & data )
{
    gsDebug<< "Reading parasolid..\n";

    PK_PART_receive_o_t receive_options;
    PK_ERROR_code_t err;
    PK_PART_t * part = 0;  // Pointer to input part
    int numParts = 0;

    // Set receive options
    PK_PART_receive_o_m(receive_options);
    receive_options.transmit_format = PK_transmit_format_text_c;

    // Start Parasolid session
    gsPKSession::start();

    // Receive part (body, solid, sheet, wireframe, assembly)
    err = PK_PART_receive(fname, &receive_options, &numParts, &part);
    PARASOLID_ERROR(PK_PART_receive, err);

    if( err == PK_ERROR_schema_access_error )  
    {  
        // Schema folder was not found
        gsWarn<< " Could not find schema files (sch_*.sch_txt)"
              << " that allow reading the input file.\n";
        return false;
    } 
    
    gsDebugVar(numParts);

    if ( numParts == 0 )
        gsDebug<< "No parts received.\n";

    // Read the parts
    for ( int j = 0; j!= numParts; ++j) // for all parts
    {
        readPK_PART_geoms(part[j], data);
        readPK_PART(part[j], data);
    }

    PK_MEMORY_free(part);
        
    // Stop Parasolid session
    gsPKSession::stop();
    
    return true;    
}

bool readPK_SURF( const PK_SURF_t & surf, internal::gsXmlTree & data  )
{
    PK_BSURF_t  bsurf;     // B-spline object

    // if ( surf is a BSURF)
    {
        PK_UVBOX_t  uvbox;     // Parameter interval
        PK_ERROR_code_t err = PK_SURF_ask_uvbox( surf, &uvbox);
        PARASOLID_ERROR(PK_SURF_ask_uvbox, err);    
       
        // /*
        // Args v1: surface, bounds, force_cubic, force_non_rational, tolerance, 
        // created_BSURF, b_surf_is_exact
        PK_LOGICAL_t exact;
        err = PK_SURF_make_bsurf( surf, uvbox, false, true, 1e-3, &bsurf, &exact );
        PARASOLID_ERROR(PK_SURF_make_bsurf, err);

        if ( ! exact )
            gsWarn<< "(!) Surface was approximated.\n";

        //*/

         /*
        // Args v2: surface, bounds, force_cubic, force_non_rational, tolerance, 
        // created_BSURF, b_surf_is_exact
        PK_SURF_make_bsurf_o_s cr_opts;
        PK_SURF_make_bsurf_t cr_stat;
        PK_achieved_cont_t cr_cont;
        //cr_opts.tolerance = 1e-3;
        double tol;
        err = PK_SURF_make_bsurf_2(surf, uvbox, &cr_opts, &cr_stat, &bsurf, &tol, &cr_cont );
        PARASOLID_ERROR(PK_SURF_make_bsurf_2, err);
        //*/
        
        // Read in surface to XML
        return readPK_BSURF(bsurf, data);
    }
}

//bool readPK_GEOM( const PK_GEOM_t & pgeo, internal::gsXmlTree & data  )

bool readPK_BSURF( const PK_BSURF_t & bsurf, internal::gsXmlTree & data  )
{
    gsDebug<< "Reading BSUFR..\n";

    internal::gsXmlNode* parent = data.first_node("xml") ;

    PK_BSURF_sf_t sf;      // B-spline data (standard form)

    // Get standard form
    PK_ERROR_code_t err = PK_BSURF_ask(bsurf, &sf);
    PARASOLID_ERROR(PK_BSURF_ask, err);
/*
    gsInfo <<" Reading a PK_BSURF, with: \n";
    gsInfo <<"Degree: "<<  sf.u_degree <<", "<<  sf.v_degree <<".\n";
    gsInfo <<"Num. of coefs: "<<  sf.n_u_vertices <<", "<< sf.n_v_vertices <<".\n";
    gsInfo <<"dim of coefs: "<< sf.vertex_dim <<"\n";
    gsInfo <<"is_rational: "<< static_cast<bool>(sf.is_rational) <<"\n";

    gsInfo <<"Num knots u : "<< sf.n_u_knots <<"\n";
    //gsInfo <<"Vals knots u: "<< sf.u_knot     << "\n";
    //gsInfo <<"Mult knots u: "<< sf.u_knot_mult <<"\n";

    gsInfo <<"Num knots v : "<< sf.n_v_knots <<"\n";
    //gsInfo <<"Vals knots v: "<< sf.v_knot     << "\n";
    //gsInfo <<"Mult knots v: "<< sf.v_knot_mult <<"\n";

    gsInfo <<"is_u_periodic: "<< static_cast<bool>(sf.is_u_periodic) <<"\n";
    gsInfo <<"is_v_periodic: "<< static_cast<bool>(sf.is_v_periodic) <<"\n";
    gsInfo <<"is_u_closed: "<< static_cast<bool>(sf.is_u_closed) <<"\n";
    gsInfo <<"is_v_closed: "<< static_cast<bool>(sf.is_v_closed) <<"\n";
    gsInfo <<"self_intersecting: "<< static_cast<bool>(sf.self_intersecting) <<"\n";
    gsInfo <<"convexity: "<< sf.convexity <<"\n";
//*/

    std::stringstream str;
    str << std::setprecision(16);
    internal::gsXmlNode* g = internal::makeNode("Geometry", data);
    g->append_attribute( internal::makeAttribute("type", "TensorBSpline2", data) );
    parent->append_node(g);
    parent= g;
    
    // Tensor Basis
    internal::gsXmlNode* tb = internal::makeNode("Basis", data);
    tb->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis2", data) );
    parent->append_node(tb);
    
    // v-direction
    internal::gsXmlNode* b = internal::makeNode("Basis", data);
    b->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
    b->append_attribute( internal::makeAttribute("index", "1", data) );
    tb->append_node(b);
    str.clear(); str.str("");
    for ( int i=0; i!=sf.n_v_knots; ++i )
        for ( int j=0; j!=sf.v_knot_mult[i]; ++j )
            str << sf.v_knot[i] <<" ";
    g = internal::makeNode("KnotVector", str.str(), data);
    g->append_attribute( internal::makeAttribute("degree", sf.v_degree, data ) ) ;
    b->append_node(g);

    // u-direction
    b = internal::makeNode("Basis", data);
    b->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
    b->append_attribute( internal::makeAttribute("index", "0", data) );
    tb->append_node(b);
    str.clear(); str.str("");
    for ( int i=0; i!=sf.n_u_knots; ++i )
        for ( int j=0; j!=sf.u_knot_mult[i]; ++j )
            str << sf.u_knot[i] <<" ";
    g = internal::makeNode("KnotVector", str.str(), data);
    g->append_attribute( internal::makeAttribute("degree", sf.u_degree, data ) ) ;
    b->append_node(g);
        
    // Coefficients
    str.clear(); str.str("");
    int c = 0;
    for ( int j=0; j!=sf.n_u_vertices; ++j )
        for ( int i=0; i!=sf.n_v_vertices; ++i )
            for ( int k=0; k!=sf.vertex_dim; ++k )
                str << sf.vertex[c++] <<" ";
    g = internal::makeNode("coefs", str.str(), data);
    g->append_attribute( internal::makeAttribute("geoDim", "3", data ) );
    parent->append_node(g);
    
    //Clean up
    free( sf.u_knot     );
    free( sf.u_knot_mult);
    free( sf.v_knot     );
    free( sf.v_knot_mult);
    free( sf.vertex     );
    return true;
}


/*
bool readPK_BCURVE( const PK_BCURVE_t& bcurv, internal::gsXmlTree & data  )
{
    internal::gsXmlNode* parent = data.first_node("xml") ;

    PK_BCURVE_sf_t sf;      // B-spline data (standard form)

    // Get standard form
    PK_ERROR_code_t err = PK_BCURVE_ask(bsurf, &sf);
    PARASOLID_ERROR(PK_BCURVE_ask, err);

}
*/


bool readPK_PART( const PK_PART_t & part, internal::gsXmlTree & data  )
{
    gsDebug<< "Reading part..\n";

    // The word part is used to mean assembly or body.
	// If the passed part is an assembly loop through each part in
	// it. Each part may also be an assembly.
	PK_CLASS_t part_class = 0;
	PK_ENTITY_ask_class(part, &part_class);

	if (part_class == PK_CLASS_body)
	{
        return readPK_BODY(part, data);
    }
	else if (part_class == PK_CLASS_assembly)
	{
		PK_PART_t * parts = 0;
		int       n_parts = 0;
        PK_ERROR_code_t err = PK_ASSEMBLY_ask_parts(part, &n_parts, &parts);
        PARASOLID_ERROR(PK_ASSEMBLY_ask_parts, err);
        
        gsDebugVar(n_parts);

        for(int i_part = 0 ; i_part < n_parts ; i_part++)
        {
            readPK_PART(parts[i_part], data);	// recursive
        }

        PK_MEMORY_free(parts);
        return true;
	}

    return false;
}


bool readPK_BODY( const PK_BODY_t & body, internal::gsXmlTree & data  )
{
    gsDebug<< "Reading body..\n";

    PK_ERROR_code_t err;
    PK_SURF_t   surf;      // Surface on face
    PK_FACE_t * faces = 0;
    int n_faces = 0;

    // ask for faces
    err = PK_BODY_ask_faces(body, &n_faces, &faces);
    PARASOLID_ERROR(PK_BODY_ask_faces, err);
    
    for(int i = 0 ; i < n_faces ; i++)
    {
        //gsDebug<< numFaces << " faces received from part "<<j<<".\n";
        
        // ask for surface of face
        err = PK_FACE_ask_surf( faces[i], &surf );
        PARASOLID_ERROR(PK_FACE_ask_surf, err);

        readPK_SURF(surf, data);
    }

    PK_MEMORY_free(faces);    
    return true;
}


bool readPK_PART_geoms( const PK_BODY_t & body, internal::gsXmlTree & data  )
{
    gsDebug<< "Reading part geometries..\n";

    PK_ERROR_code_t err;
    PK_LOGICAL_t b;
    int n_geo = 0;
    PK_GEOM_t * geoms = 0;

    // Surfaces
    err = PK_PART_ask_geoms(body, &n_geo, &geoms);
    gsDebugVar(n_geo);

    PARASOLID_ERROR(PK_BODY_ask_surfs, err);
    for(int i = 0 ; i < n_geo; i++)
    {
        PK_ENTITY_is_surf(geoms[i], &b);
        if (b)
            readPK_SURF(geoms[i], data);
        else
            gsWarn<<"Unknown entity..\n";

        PK_ENTITY_is_curve(geoms[i], &b);
        if (b)
            gsWarn<<"Curve entity..\n";
    }

    PK_MEMORY_free(geoms);

    return true;
}


}//extensions

}//gismo

#undef PARASOLID_ERROR
