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

#include <gsUtils/gsMesh/gsMesh.h>

#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsTHBSplineBasis.h>


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
    std::cout << "write parasolid mulitpatch" << std::endl;
    
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


template<class T>
bool gsWriteParasolid( const gsMesh<T>& mesh, const std::string & filename)
{
    gsPKSession::start();

    PK_LOGICAL_t checks(0);
    PK_SESSION_set_check_continuity(checks);
    PK_SESSION_set_check_self_int(checks);
    
    PK_ASSEMBLY_t assembly;
    exportMesh(mesh, assembly);

    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;

    PK_ERROR_code_t err = PK_PART_transmit(1, &assembly, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);
    
    gsPKSession::stop();
    
    return err;
}


template<class T>
bool gsWriteParasolid(const gsTHBSpline<2, T>& thb, const std::string& filename)
{
    gsPKSession::start();
    PK_LOGICAL_t checks(0);
    PK_SESSION_set_check_continuity(checks);
    PK_SESSION_set_check_self_int(checks);
    
    PK_ASSEMBLY_t assembly;
    exportTHBsurface<real_t>(thb, assembly);
    
    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;    

    PK_ERROR_code_t err = PK_PART_transmit(1, &assembly, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);
    
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
        createPK_BSURF< T, gsKnotVector<T> >(*tbsp, pgeo);
    }
// the following lines produce warnings if called from multipatch version of gsWriteParasolid, 
// because it already assumes that the geometries are surfaces
    else if ( const gsBSpline<>* bspl = 
	      dynamic_cast< const gsBSpline<>* >(&ggeo) )
    {
	createPK_BCURVE(*bspl, pgeo);
    }
    else
    {
        gsInfo << "Cannot write "<<ggeo<<" to parasolid file.\n";
    }
}


template<class T, class KnotVectorType> void
createPK_BSURF(const gsTensorBSpline< 2, T, KnotVectorType > & bsp, 
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

template<class T> void
exportMesh(const gsMesh<T>& mesh,
	   PK_ASSEMBLY_t& assembly)
{
    // tried to:
    //  - make one wire body out of all mesh edges with PK_CURVE_make_wire_body_2
    //    and it doesn't work, because edges cross each other
    //  - make a boolean union of all mesh edges with PK_BODY_boolean_2
    //    and it doesn't work function fails to make a union body (I don't know 
    //    the reason)

    PK_ERROR_t err = PK_ASSEMBLY_create_empty(&assembly);
    PARASOLID_ERROR(PK_ASSEMBLY_create_empty, err);

    gsKnotVector<T> kv(0, 1, 0, 2);
    gsMatrix<T> coefs(2, 3);
    gsBSpline<T> bspl(kv, coefs);

    gsMatrix<T> newCoefs(2, 3);
    for (int i = 0; i != mesh.numEdges; ++i)
    {
	newCoefs.row(0) = mesh.edge[i].source->coords.transpose();
	newCoefs.row(1) = mesh.edge[i].target->coords.transpose();
	
	if ((newCoefs.row(0) - newCoefs.row(1)).norm() < 1e-6)
	{
	    continue;
	}

	bspl.setCoefs(newCoefs);
	
	PK_BCURVE_t bcurve;
	createPK_BCURVE(bspl, bcurve);
	
	err = PK_PART_add_geoms(assembly, 1, &bcurve);
	PARASOLID_ERROR(PK_PART_add_geoms, err);
    }
}


// This function is here as workaround around Parasolid / MTU visualization bug.
// The bug is:
//     Let say that you make a trimmed sheet in a such way that you make a 
//     hole into geometry. So there is a region in geometry which will be trimmed 
//     away. If multiple coefficients, which define this region, are equal to 
//     (0, 0, 0) then Parasolid / MTU visualization complains because the surface 
//     self intersects. Parasolid / MTU visulation doesn't realize that this reigon
//     will be trimmed away and it is not important.
//
// The workaround:
//     We change the problematic coefficients, from the trimmed-away area, such that 
//     the bspline approximates the THB surface.
//
// TODO: refactor this function: split it to multiple parts
void
makeValidGeometry(const gsTHBSpline<2>& surface,
		  gsTensorBSpline<2, real_t, gsCompactKnotVector<> >& bspline)
{
    // check how many coefs are zero
    
    gsMatrix<>& coefs = bspline.coefs();
    gsVector<int> globalToLocal(coefs.rows());
    globalToLocal.setConstant(-1);
    
    int numZeroRows = 0;
    for (index_t row = 0; row != coefs.rows(); row++)
    {
	bool zeroRow = true;
	for (index_t col = 0; col != coefs.cols(); col++)
	{
	    if (coefs(row, col) != 0.0)
	    {
		zeroRow = false;
		break;
	    }
	}
	
	if (zeroRow)
	{
	    globalToLocal(row) = numZeroRows;
	    numZeroRows++;
	}
    }
 
    // if more than 2 coefs ontinue

    if (numZeroRows < 2)
    {
	return;
    }
    
    // evaluate basis functions with zero coefs on grevielle points

    gsMatrix<> anchors;
    bspline.basis().anchors_into(anchors);

    gsMatrix<> params(2, coefs.rows());
    gsVector<> localToGlobal(coefs.rows());
    index_t counter = 0;

    for (index_t fun = 0; fun != globalToLocal.rows(); fun++)
    {
	if (globalToLocal(fun) != -1)
	{
	    params.col(counter) = anchors.col(fun);
	    localToGlobal(counter) = fun;
	    counter++;
	}
    }
    
    params.conservativeResize(2, counter);
    localToGlobal.conservativeResize(counter);
    
    gsMatrix<> points;
    surface.eval_into(params, points);
    
    // do least square fit 

    gsSparseMatrix<> A(localToGlobal.rows(), localToGlobal.rows());
    A.setZero();
    gsMatrix<> B(localToGlobal.rows(), 3);
    B.setZero();
    
    gsMatrix<> value;
    gsMatrix<unsigned> actives;
    
    for (index_t k = 0; k != params.cols(); k++)
    {
	const gsMatrix<>& curr_param = params.col(k);
	
	bspline.basis().eval_into(curr_param, value);
	bspline.basis().active_into(curr_param, actives);
	const index_t numActive = actives.rows();
	
	for (index_t i = 0; i != numActive; i++)
	{
	    const int I = globalToLocal(actives(i, 0));
	    if (I != -1)
	    {
		B.row(I) += value(i, 0) * points.col(k);
		
		for (index_t j = 0; j != numActive; j++)
		{
		    const int J = globalToLocal(actives(j, 0));
		    if (J != -1)
		    {
			A(I, J) += value(i, 0) * value(j, 0);
		    }
		}
	    }
	}
    }

    A.makeCompressed();
    Eigen::BiCGSTAB< gsSparseMatrix<>,  Eigen::IncompleteLUT<real_t> > solver(A);
    if (solver.preconditioner().info() != Eigen::Success)
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        return;
    }

    // update coefficients

    gsMatrix<> x = solver.solve(B);
    
    for (index_t row = 0; row != x.rows(); row++)
    {
	coefs.row(localToGlobal(row)) = x.row(row);
    }
}



template <class T> void
exportTHBsurface(const gsTHBSpline<2, T>& surface,
		 PK_ASSEMBLY_t& assembly)
{
    PK_ERROR_t err = PK_ASSEMBLY_create_empty(&assembly);
    PARASOLID_ERROR(PK_ASSEMBLY_create_empty, err);

    gsTHBSplineBasis<2>::AxisAlignedBoundingBox boundaryAABB;
    gsTHBSplineBasis<2>::TrimmingCurves trimCurves;
   
    const gsTHBSplineBasis<2>& basis= surface.basis();
    
    basis.decomposeDomain(boundaryAABB, trimCurves);
    
    for (unsigned level = 0; level != boundaryAABB.size(); level++)
    {

	for (unsigned box = 0; box != boundaryAABB[level].size(); box++)
	{
	    gsTensorBSpline<2, T, gsCompactKnotVector<> > bspline =
		basis.getBSplinePatch(boundaryAABB[level][box], level, surface.coefs());
	    
	    makeValidGeometry(surface, bspline);

	    PK_BSURF_t bsurf;
	    createPK_BSURF<T, gsCompactKnotVector<> >(bspline, bsurf);
	    
	    std::vector<PK_CURVE_t> curves;
	    std::vector<PK_INTERVAL_t> intervals;
	    std::vector<int> trim_loop;
	    std::vector<int> trim_set;
	    
	    for (unsigned curve = 0; curve != trimCurves[level][box].size(); curve++)
	    {
		for (unsigned seg = 0; seg != trimCurves[level][box][curve].size(); seg++)
		{
		    real_t x1 = trimCurves[level][box][curve][seg][0];
		    real_t y1 = trimCurves[level][box][curve][seg][1];
		    real_t x2 = trimCurves[level][box][curve][seg][2];
		    real_t y2 = trimCurves[level][box][curve][seg][3];
		    
		    PK_CURVE_t line;
		    PK_INTERVAL_t intervalDummy;
		    PK_INTERVAL_t interval;
		    
		    PK_SURF_make_curve_isoparam_o_t options;
		    PK_SURF_make_curve_isoparam_o_m(options);
		    
		    if (x1 == x2)
		    {	
			err = PK_SURF_make_curve_isoparam(bsurf, x1, PK_PARAM_direction_v_c, 
							  &options, &line, &intervalDummy);
			PARASOLID_ERROR(PK_SURF_make_curve_isoparam, err);
			
			interval.value[0] = y1;
			interval.value[1] = y2;
		    }
		    else
		    {
			err = PK_SURF_make_curve_isoparam(bsurf, y1, PK_PARAM_direction_u_c, 
							  &options, &line, &intervalDummy);
			PARASOLID_ERROR(PK_SURF_make_curve_isoparam, err);
			
			interval.value[0] = x1;
			interval.value[1] = x2;
		    }
		    
		    curves.push_back(line);
		    intervals.push_back(interval);
		    trim_loop.push_back(curve);
		    trim_set.push_back(0);
		}
	    }
	    
	    PK_SURF_trim_data_t trim_data;
	    trim_data.n_spcurves = static_cast<int>(curves.size());
	    trim_data.spcurves = curves.data();
	    trim_data.intervals = intervals.data();
	    trim_data.trim_loop = trim_loop.data();
	    trim_data.trim_set = trim_set.data();
	    
	    PK_SURF_make_sheet_trimmed_o_t options;
	    PK_SURF_make_sheet_trimmed_o_m(options);
	    options.check_loops = PK_LOGICAL_true;
	    
	    PK_BODY_t trimSurface;
	    PK_check_state_t state;
	    err = PK_SURF_make_sheet_trimmed(bsurf, trim_data, 1e-8, &options, 
					     &trimSurface, &state);
	    PARASOLID_ERROR(PK_SURF_make_sheet_trimmed, err);	    
	    
	    if (state != PK_BODY_state_ok_c)
	    {
		std::cout << "Something went wrong." << std::endl;
	    }
	    
	    PK_INSTANCE_sf_t sform;
	    sform.assembly = assembly;
	    sform.transf = PK_ENTITY_null;
	    sform.part = trimSurface;

	    PK_INSTANCE_t instance;
	    err = PK_INSTANCE_create(&sform, &instance);
	    PARASOLID_ERROR(PK_INSTANCE_create, err);
	}
    }
}



}//extensions

}//gismo

#undef PARASOLID_ERROR
