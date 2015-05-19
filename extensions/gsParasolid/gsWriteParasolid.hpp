/** @file gsWriteParasolid.hpp

    @brief Provides implementation of gsWriteParasolid function.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsParasolid/gsFrustrum.h>

#include <gsParasolid/gsReadParasolid.h>
#include <gsParasolid/gsPKSession.h>

#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsSurface.h>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsBSpline.h>



namespace gismo {

namespace extensions {

/*
Notes:

// The sheet is a kind of part
// PK_PART_t part = sheet;

This function attaches surfaces to faces.
PK_FACE_attach_surfs

//This function creates a sheet body from a collection of faces.
PK_FACE_make_sheet_body
PK_FACE_make_sheet_bodies

//This function creates solid bodies from a collection of faces
PK_FACE_make_solid_bodies

    //PK_FACE_make_sheet_body
    //PK_FACE_make_solid_bodies

    //PK_PART_find_entity_by_ident


    // Create parasolid sheet-body out of the surface
    //PK_UVBOX_t  uvbox;      // Parameter domain intevals in u and v
    //PK_SURF_ask_uvbox(bsurf, &uvbox); // Get the parameter box
    //PK_BODY_t   sheet;      // Sheet is a kind of body
    //err = PK_SURF_make_sheet_body(bsurf, uvbox, &sheet);
    //PARASOLID_ERROR(PK_SURF_make_sheet_body, err);  

*/


template<class T>
bool gsWriteParasolid( const gsMultiPatch<T> & gssurfs, std::string const & filename )
{
    PK_ERROR_code_t err;
    PK_BODY_t  part;     // Empty part
    PK_GEOM_t   geo[ gssurfs.nPatches() ];      // Geometries

    // Start Parasolid session
    gsPKSession::start();

    // Disable continuity and self-intersection checks on geometrical
    // data of bodies
    PK_LOGICAL_t checks(0);
    PK_SESSION_set_check_continuity(checks);
    PK_SESSION_set_check_self_int(checks);

    // Create Parasolid geometries
    int count = 0;
    for (typename gsMultiPatch<T>::const_iterator 
             it = gssurfs.begin(); it != gssurfs.end(); ++it)
    {
        createPK_GEOM(**it, geo[count++] );
    }

/*
    // Push all geometries as orphans in an assembly part
    PK_ASSEMBLY_create_empty(&part);
    err = PK_PART_add_geoms(part, count, geo);
    PARASOLID_ERROR(PK_PART_add_geoms, err);  
//*/

    // Make a sheet body out of each geometry 
    PK_BODY_t  parts[count];
    for (int i = 0; i!= count; ++i )
    {
        PK_UVBOX_t  uvbox;     // Parameter interval
        PK_ERROR_code_t err = PK_SURF_ask_uvbox(geo[i], &uvbox);
        PARASOLID_ERROR(PK_SURF_ask_uvbox, err);    

        PK_SURF_make_sheet_body(geo[i], uvbox, &part);
        PARASOLID_ERROR(PK_SURF_make_sheet_body, err);    

        parts[i] = part;
    }

    // Write it out to the file
    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;
    err = PK_PART_transmit(count, parts, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);

    // Stop Parasolid session
    gsPKSession::stop();

    return err;
}


template<class T>
bool gsWriteParasolid( const gsGeometry<T> & ggeo, std::string const & filename )
{
    PK_ERROR_code_t err;
    PK_ASSEMBLY_t  part;
    PK_GEOM_t  pgeo; // Parasolid geometric entity

    // Start Parasolid session
    gsPKSession::start();

    // Disable continuity and self-intersectiion checks
    PK_LOGICAL_t checks(0);
    PK_SESSION_set_check_continuity(checks);
    PK_SESSION_set_check_self_int(checks);

    // Create Parasolid geometry
    createPK_GEOM(ggeo, pgeo);
    
    // Create parasolid part out of the geometry
    PK_ASSEMBLY_create_empty(&part);
    err = PK_PART_add_geoms(part, 1, &pgeo);
    PARASOLID_ERROR(PK_PART_add_geoms, err);  

    // Write it out to the file
    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;
    err = PK_PART_transmit(1, &part, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);

    // Stop Parasolid session
    gsPKSession::stop();

    return err;
}


template<class T> void
createPK_GEOM( const gsGeometry<T> & ggeo, 
             PK_GEOM_t & pgeo)
{
    // Identify input gismo geometry
    if ( const gsTensorBSpline<2,T> * tbsp = 
         dynamic_cast<const gsTensorBSpline<2,T> *>(&ggeo) )
    {
        createPK_BSURF(*tbsp, pgeo);
    }
// the following lines produce warnings, because writing a multipatch already assumes 
// that the geometries are surfaces
//     else if ( const gsBSpline<>* bspl = 
// 	      dynamic_cast< const gsBSpline<>* >(&ggeo) )
//     {
// 	createPK_BCURVE(*bspl, pgeo);
//     }
    else
    {
        gsInfo << "Cannot write "<<ggeo<<" to parasolid file.\n";
    }
}


template<class T> void
createPK_BSURF( const gsTensorBSpline<2,T> & bsp, 
             PK_BSURF_t & bsurf)
{
    // Translate to parasolid standard form, ie fill up parasolid
    // spline data record
    PK_BSURF_sf_t sform;   // B-spline data holder (standard form)

    // Degrees
    sform.u_degree      = bsp.basis().degree(1);
    sform.v_degree      = bsp.basis().degree(0);

    // Knots in u-direction
    std::vector<T> gknot0 = bsp.basis().knots(1).unique();
    std::vector<int> gmult0 = bsp.basis().knots(1).multiplicities();
    sform.n_u_knots     = gknot0.size();
    sform.u_knot        = gknot0.data();
    sform.u_knot_mult   = gmult0.data();

    // Knots in v-direction
    std::vector<T> gknot1 = bsp.basis().knots(0).unique();
    std::vector<int> gmult1 = bsp.basis().knots(0).multiplicities();
    sform.n_v_knots     = gknot1.size();
    sform.v_knot        = gknot1.data();
    sform.v_knot_mult   = gmult1.data();

    // Control points
    sform.n_u_vertices  = bsp.basis().size(1);
    sform.n_v_vertices  = bsp.basis().size(0);
    gsMatrix<T> coefs   = bsp.coefs();
    const int n = bsp.geoDim();
    if ( n < 3 )
    {
        coefs.conservativeResize(Eigen::NoChange, 3);
        coefs.rightCols(3-n).setZero();
    }
    coefs.transposeInPlace();
    coefs.resize(3*bsp.basis().size(), 1);
    sform.vertex_dim    = 3; // always 3 for surfaces
    sform.vertex        = coefs.data();

    // Attributes
    sform.is_rational   = PK_LOGICAL_false;
    sform.form          = PK_BSURF_form_unset_c;
    sform.u_knot_type   = PK_knot_unset_c;
    sform.v_knot_type   = PK_knot_unset_c;
    sform.is_u_periodic = PK_LOGICAL_false;
    sform.is_v_periodic = PK_LOGICAL_false;
    sform.is_u_closed   = PK_LOGICAL_false;
    sform.is_v_closed   = PK_LOGICAL_false;
    sform.self_intersecting = PK_self_intersect_unset_c;
    sform.convexity         = PK_convexity_unset_c;

    // Create parasolid surface with the previous spline data
    PK_ERROR_code_t err = PK_BSURF_create(&sform, &bsurf);
    PARASOLID_ERROR(PK_BSURF_create, err);
}

template<class T> void
createPK_BCURVE( const gsBSpline<T>& curve, 
		 PK_BCURVE_t& bcurve)
{
    PK_BCURVE_sf_t sform; // B-curve data holder (standard form)

    // Degree
    sform.degree = curve.degree();
    
    // Knots
    std::vector<T> knots = curve.basis().knots().unique();
    std::vector<int> mult = curve.basis().knots().multiplicities();
    sform.n_knots = knots.size();
    sform.knot = knots.data();
    sform.knot_mult = mult.data();


    // Control points
    sform.n_vertices = curve.basis().size();
    gsMatrix<T> coefs = curve.coefs();
    const int n = curve.geoDim();
    if (n < 3)
    {
	coefs.conservativeResize(Eigen::NoChange, 3);
	coefs.rightCols(3 - n).setZero();
    }
    coefs.transposeInPlace();
    coefs.resize(3 * curve.basis().size(), 1);
    sform.vertex_dim = 3;
    sform.vertex = coefs.data();
    
    // Attributes
    sform.is_rational = PK_LOGICAL_false;
    sform.form = PK_BCURVE_form_unset_c;
    sform.knot_type = PK_knot_unset_c;
    sform.is_periodic = PK_LOGICAL_false;
    sform.is_closed = PK_LOGICAL_false;
    sform.self_intersecting = PK_self_intersect_unset_c;
    
    PK_ERROR_code_t err = PK_BCURVE_create(&sform, &bcurve);
    PARASOLID_ERROR(PK_BCURVE_create, err);
}

}//extensions

}//gismo

#undef PARASOLID_ERROR
