/** @file gsWriteParasolid.hpp

    @brief Provides implementation of gsWriteParasolid function.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mokris
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

#include <gsIO/gsWriteParaview.h>


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
class gsTrimData
{
public:
    gsTrimData(unsigned level,std::vector<index_t>& AABBBox,std::vector<std::vector<std::vector<T> > >polylines)
        : m_level(level),m_AABBBox(AABBBox),m_Polylines(polylines)
    { }

    gsTrimData() : m_level(0),m_AABBBox(),m_Polylines()
    { }

public:
    unsigned m_level;
    std::vector<index_t> m_AABBBox;
    std::vector<std::vector<std::vector<T> > > m_Polylines;
};

// forward declarations
void makeValidGeometry(const gsTHBSpline<2>& surface,
                       gsTensorBSpline<2, real_t>& bspline);

void getInterval(const bool directionU,
                 const real_t param1,
                 const real_t param2,
                 const real_t paramConst,
                 const gsTensorBSpline<2, real_t>& bspline,
                 const PK_CURVE_t line,
                 PK_INTERVAL_t& result);

bool validMultiplicities(const std::vector<index_t>& mult,
                         const int deg);

template <class T>
bool exportCheck(const gsTHBSpline<2, T>& surface);

template <class T>
bool exportTHBsurface( const gsTHBSpline<2, T>& surface,
                       gsTHBSplineBasis<2>::TrimmingCurves trimCurves,
                       gsTHBSplineBasis<2>::AxisAlignedBoundingBox boundaryAABB,
                       PK_ASSEMBLY_t& body );

template <class T>
void getTrimCurvesAndBoundingBoxes(const gsTHBSpline<2, T>& surface,
                                   std::vector<index_t>& boxes,
                                   gsTHBSplineBasis<2>::TrimmingCurves& trimCurves,
                                   gsTHBSplineBasis<2>::AxisAlignedBoundingBox& boundaryAABB);

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
    exportTHBsurface<T>(thb, assembly);

    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;

    PK_ERROR_code_t err = PK_PART_transmit(1, &assembly, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);

    gsPKSession::stop();

    return err;
}

template<class T>
bool gsWriteParasolid( const gsTHBSpline<2, T>& thb,const std::vector<T>& par_boxes, const std::string& filename )
{
    gsPKSession::start();
    PK_LOGICAL_t checks(0);
    PK_SESSION_set_check_continuity(checks);
    PK_SESSION_set_check_self_int(checks);

    PK_ASSEMBLY_t assembly;
    exportTHBsurface<T>(thb,par_boxes,assembly);

    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;

    PK_ERROR_code_t err = PK_PART_transmit(1, &assembly, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);
    gsPKSession::stop();
    return err;
}

template<class T>
bool gsWritePK_SHEET(const gsTensorBSpline<2, T>& tp, const std::string& filename)
{
    gsPKSession::start();
    PK_LOGICAL_t checks(0);
    PK_SESSION_set_check_continuity(checks);
    PK_SESSION_set_check_self_int(checks);

    PK_ERROR_code_t err;

    PK_BSURF_t bsurf;
    createPK_BSURF<T>(tp, bsurf, false, false);

    PK_UVBOX_t uv_box;
    uv_box.param[0] = 0;
    uv_box.param[1] = 0;
    uv_box.param[2] = 1;
    uv_box.param[3] = 1;
    PK_BODY_t body;
    err = PK_SURF_make_sheet_body(bsurf, uv_box, &body);
    PARASOLID_ERROR(PK_SURF_make_sheet_body, err);

    // PK_ASSEMBLY_t assembly;
    // PK_ASSEMBLY_create_empty(&assembly);
    // err = PK_PART_add_geoms(assembly, 1, &bsurf);
    // PARASOLID_ERROR(PK_PART_add_geoms, err);

    PK_PART_transmit_o_t transmit_options;
    PK_PART_transmit_o_m(transmit_options);
    transmit_options.transmit_format = PK_transmit_format_text_c;

    err = PK_PART_transmit(1, &body, filename.c_str(), &transmit_options);
    PARASOLID_ERROR(PK_PART_transmit, err);

    gsPKSession::stop();

    return err;
}


template<class T>
bool createPK_GEOM( const gsGeometry<T> & ggeo,
                    PK_GEOM_t & pgeo)
{
    // Identify input gismo geometry
    if ( const gsTensorBSpline<2,T> * tbsp =
         dynamic_cast<const gsTensorBSpline<2,T> *>(&ggeo) )
    {
        return createPK_BSURF(*tbsp, pgeo);
    }
// the following lines produce warnings if called from multipatch version of gsWriteParasolid,
// because it already assumes that the geometries are surfaces
    else if ( const gsBSpline<>* bspl =
              dynamic_cast< const gsBSpline<>* >(&ggeo) )
    {
        return createPK_BCURVE(*bspl, pgeo);
    }
    else
    {
        gsInfo << "Cannot write "<<ggeo<<" to parasolid file.\n";
        return false;
    }
}


template<class T>
bool createPK_BSURF(const gsTensorBSpline< 2, T> & bsp,
                    PK_BSURF_t & bsurf,
                    bool closed_u,
                    bool closed_v)
{
    typedef typename gsKnotVector<T>::mult_t mult_t;
    std::vector<mult_t> mult;

    // Gismo and Parasolid store the coefficients in different order.
    // Therefore we create the BSURF with u and v swapped and then transpose it.
    for (index_t dim = 0; dim != 2; dim++)
    {
        const int deg = bsp.basis().degree(dim);
        mult = bsp.basis().knots(dim).multiplicities();

        if (!validMultiplicities(mult, deg))
        {
            return false;
        }
    }

    // Translate to parasolid standard form, ie fill up parasolid
    // spline data record
    PK_BSURF_sf_t sform;   // B-spline data holder (standard form)

    // Degrees
    sform.u_degree      = bsp.basis().degree(1);
    sform.v_degree      = bsp.basis().degree(0);

    // Knots in u-direction
    std::vector<T> gknot0 = bsp.basis().knots(1).unique();
    std::vector<mult_t> gmult0 = bsp.basis().knots(1).multiplicities();
    sform.n_u_knots     = gknot0.size();
    sform.u_knot        = gknot0.data();
    sform.u_knot_mult   = gmult0.data();

    // Knots in v-direction
    std::vector<T> gknot1 = bsp.basis().knots(0).unique();
    std::vector<mult_t> gmult1 = bsp.basis().knots(0).multiplicities();
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

    if (closed_u)
    {
        sform.is_v_closed = PK_LOGICAL_true;
    }

    if (closed_v)
    {
        sform.is_u_closed = PK_LOGICAL_true;
    }

    sform.self_intersecting = PK_self_intersect_unset_c;
    sform.convexity         = PK_convexity_unset_c;

    // Create parasolid surface with the previous spline data
    PK_ERROR_code_t err = PK_BSURF_create(&sform, &bsurf);
    PARASOLID_ERROR(PK_BSURF_create, err);

    // Transposition (new on 2019-02-26).
    PK_BSURF_reparameterise_o_t options;
    PK_BSURF_reparameterise_o_m(options);
    options.transpose = PK_LOGICAL_true;

    err = PK_BSURF_reparameterise(bsurf,&options);
    PARASOLID_ERROR(PK_BSURF_reparameterise, err);

    return true;
}

template<class T>
bool createPK_BCURVE( const gsBSpline<T>& curve,
                      PK_BCURVE_t& bcurve)
{
    typedef typename gsKnotVector<T>::mult_t mult_t;
    
    PK_BCURVE_sf_t sform; // B-curve data holder (standard form)

    // Degree
    sform.degree = curve.degree();

    // Knots
    std::vector<T> knots = curve.basis().knots().unique();
    std::vector<mult_t> mult = curve.basis().knots().multiplicities();
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

    return true;
}

template<class T>
bool exportMesh(const gsMesh<T>& mesh,
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
    for (size_t i = 0; i != mesh.numEdges(); ++i)
    {
        newCoefs.row(0) = mesh.edges()[i].source->transpose();
        newCoefs.row(1) = mesh.edges()[i].target->transpose();

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

    return true;
}


template <class T>
bool exportTHBsurface(const gsTHBSpline<2, T>& surface,
                      PK_ASSEMBLY_t& assembly)
{
    if(!exportCheck(surface))
        return false;

    gsTHBSplineBasis<2>::AxisAlignedBoundingBox boundaryAABB;
    gsTHBSplineBasis<2>::TrimmingCurves trimCurves;

    const gsTHBSplineBasis<2>& basis= surface.basis();

    basis.decomposeDomain(boundaryAABB, trimCurves);

    std::vector<gsTrimData<T> > trimData;

    bool success = getTrimDataFromBoundaryTrimCurves(boundaryAABB,trimCurves,trimData);

    if(!success)
        return false;

    return exportTHBsurface(surface,trimData,assembly);
}

template <class T>
bool exportTHBsurface(const gsTHBSpline<2, T>& surface,
                      const std::vector<T>& par_boxes,
                      PK_ASSEMBLY_t& assembly)
{
    if(par_boxes.size()==0)
        return exportTHBsurface(surface,assembly);

    if(!exportCheck(surface))
        return false;

    std::vector<gsTrimData<T> > trimData;

    bool success = getTrimCurvesAndBoundingBoxes<T>(surface,par_boxes,trimData);

    if(!success)
        return false;

    return exportTHBsurface(surface,trimData,assembly);
}


// ================================================================================
// utilities

// returns true la all multiplicities in muls are less (<) than deg + 1
// parasolid restriction
bool validMultiplicities(const std::vector<index_t>& mult,
                         const int deg)
{
    for (size_t i = 1; i != mult.size() - 1; i++)
    {
        if (mult[i] == deg + 1)
        {
            std::cout <<
                "Only B-Splines with (inner) multiplicity less (<) "
                "than degree + 1 are supported. \n"
                "Parasolid restriction."
                      << std::endl;
            return false;
        }
    }
    return true;
}


// computes the parameter interval of the iso-curve (line) on the surface (bspline)
// curve is between two points on the surface defined by input parameters
// - directionU tells if the curve is is-curve in parametric directionU
// - paramConst if the fixed parameter in parametric direction defined by directionU
// - param1 & param2 are non fixed parameters defined in the other direction than before
// - param1, param2, paramConst defines two points -- between these two points we would like
//                              to define parameter range for our curve
void getInterval(const bool directionU,
                 const real_t param1,
                 const real_t param2,
                 const real_t paramConst,
                 const gsTensorBSpline<2, real_t>& bspline,
                 const PK_CURVE_t line,
                 PK_INTERVAL_t& result)
{
    gsMatrix<> param(2, 2);
    if (directionU)
    {
        param(0, 0) = param1;
        param(1, 0) = paramConst;

        param(0, 1) = param2;
        param(1, 1) = paramConst;
    }
    else
    {
        param(0, 0) = paramConst;
        param(1, 0) = param1;

        param(0, 1) = paramConst;
        param(1, 1) = param2;
    }

    gsMatrix<> points;
    bspline.eval_into(param, points);


    for (index_t i = 0; i != 2; i++)
    {
        PK_VECTOR_t position;
        position.coord[0] = 0;
        position.coord[1] = 0;
        position.coord[2] = 0;

        for (index_t row = 0; row != points.rows(); row++)
        {
            position.coord[row] = points(row, i);
        }

        real_t p;

        PK_ERROR_t err = PK_CURVE_parameterise_vector(line, position, &p);
        PARASOLID_ERROR(PK_CURVE_parameterise_vector, err);

        result.value[i] = p;
    }

    // try to correct the problem of periodic surfaces
    if( param1!=param2 && result.value[0]==result.value[1] )
    {
        result.value[0]=param1;
        result.value[1]=param2;
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
                  gsTensorBSpline<2, real_t>& bspline)
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

    // if more than 2 coefs are zero, continue

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

    gsMatrix<> support = bspline.basis().support();

    for (index_t fun = 0; fun != globalToLocal.rows(); fun++)
    {
        if (globalToLocal(fun) != -1)
        {
            params.col(counter) = anchors.col(fun);
            localToGlobal(counter) = fun;
            for (index_t row = 0; row != support.rows(); row++)
            {
                if (params(row, counter) < support(row, 0))
                {
                    params(row, counter) = support(row, 0);
                }

                if (support(row, 1) < params(row, counter))
                {
                    params(row, counter) = support(row, 1);
                }
            }
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
    gsMatrix<> B(localToGlobal.rows(), surface.geoDim());
    B.setZero();

    gsMatrix<> value;
    gsMatrix<index_t> actives;

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
    gsSparseSolver<real_t>::BiCGSTABILUT solver(A);
    if (solver.preconditioner().info() != Eigen::Success)
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        return;
    }

    gsMatrix<> x = solver.solve(B);

    for (index_t row = 0; row != x.rows(); row++)
    {
        coefs.row(localToGlobal(row)) = x.row(row);
    }
}

template <class T>
bool exportCheck(const gsTHBSpline<2, T>& surface)
{
    typedef typename gsKnotVector<T>::mult_t mult_t;
    
    for (index_t dim = 0; dim != 2; dim++)
    {
        typedef std::vector< gsTensorBSplineBasis< 2, real_t>* > Bases;
        const Bases& bases = surface.basis().getBases();
        const int deg = (bases[0])->degree(dim);
        std::vector<mult_t> mult = (bases[0])->knots(dim).multiplicities();

        if (!validMultiplicities(mult, deg))
        {
            return false;
        }
    }

    return true;
}

template <class T>
bool exportTHBsurface( const gsTHBSpline<2, T>& surface,
                       std::vector<gsTrimData<T> >& trimData,
                       PK_ASSEMBLY_t& assembly )
{
    const gsTHBSplineBasis<2>& basis= surface.basis();

    PK_ERROR_t err = PK_SESSION_set_check_continuity(PK_LOGICAL_false);
    PARASOLID_ERROR(PK_SESSION_set_check_continuity, err);

    err = PK_ASSEMBLY_create_empty(&assembly);
    PARASOLID_ERROR(PK_ASSEMBLY_create_empty, err);

    for (size_t box = 0;box<trimData.size();++box)
    {
        unsigned level = trimData[box].m_level;
        std::vector<index_t> AABBBox = trimData[box].m_AABBBox;
        std::vector<std::vector<std::vector<T> > > polylines = trimData[box].m_Polylines;

        gsTensorBSpline<2, T> bspline =
            basis.getBSplinePatch(AABBBox, level, surface.coefs());

        //std::cout << "bspline size: " << bspline.basis().size() << std::endl;
        //std::cout << "bspline coef size: " << bspline.coefs().rows() << std::endl;
        makeValidGeometry(surface, bspline);

        //gsWriteParaview(bspline,"paraviewtest" + util::to_string(box),1000,true,true);

        PK_BSURF_t bsurf;
        createPK_BSURF<T>(bspline, bsurf); // swap

        std::vector<PK_CURVE_t> curves;
        std::vector<PK_INTERVAL_t> intervals;
        std::vector<int> trim_loop;
        std::vector<int> trim_set;

        for (unsigned loop = 0; loop != polylines.size(); loop++)
        {
            for (unsigned seg = 0; seg != polylines[loop].size(); seg++)
            {
                real_t y1 = polylines[loop][seg][0];
                real_t x1 = polylines[loop][seg][1];
                real_t y2 = polylines[loop][seg][2];
                real_t x2 = polylines[loop][seg][3];

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

                    getInterval(true, y1, y2, x1, bspline, line, interval);
                }
                else
                {
                    err = PK_SURF_make_curve_isoparam(bsurf, y1, PK_PARAM_direction_u_c,
                                                      &options, &line, &intervalDummy);
                    PARASOLID_ERROR(PK_SURF_make_curve_isoparam, err);

                    getInterval(false, x1, x2,  y1, bspline, line, interval);
                }

                curves.push_back(line);
                intervals.push_back(interval);
                trim_loop.push_back(loop);
                trim_set.push_back(0);

                // ------------------------------------------------------------
                // very usefull for debugging / logging
//            PK_PART_transmit_o_t transmit_options;
//            PK_PART_transmit_o_m(transmit_options);
//            transmit_options.transmit_format = PK_transmit_format_text_c;

//            std::string name = "level_" + internal::to_string(level) +
//            "_box_" + internal::to_string(box) +
//            "_curve_" + internal::to_string(curve) +
//            "_seg_" + internal::to_string(seg);

//            PK_GEOM_copy_o_t copyOptions;
//            PK_GEOM_copy_o_m(copyOptions);

//            PK_GEOM_copy_r_t result;
//            err = PK_GEOM_copy(1, &line, &copyOptions, &result);
//            PARASOLID_ERROR(PK_GEOM_copy, err);

//            PK_CURVE_t copyLine = *(result.copied_geoms);

//            PK_BODY_t lineBody;
//            PK_CURVE_make_wire_body_o_t option;
//            PK_CURVE_make_wire_body_o_m(option);
//            int line_new_edges;
//            PK_EDGE_t** edges = NULL;
//            int** edge_index = NULL;

//            err = PK_CURVE_make_wire_body_2(1, &copyLine, &interval, &option,
//                            &lineBody, &line_new_edges, edges, edge_index);
//            PARASOLID_ERROR(PK_CURVE_make_wire_body_2, err);



//            PK_ERROR_code_t err = PK_PART_transmit(1, &lineBody, name.c_str(), &transmit_options);
//            PARASOLID_ERROR(PK_PART_transmit, err);
                // ------------------------------------------------------------

            }
        }
        // ----------------------------------------------------------------------
        // very usefull for debugging / logging

    // copy
//        PK_GEOM_copy_o_t copyOptions;
//        PK_GEOM_copy_o_m(copyOptions);
//        PK_GEOM_copy_r_t result;
//        err = PK_GEOM_copy(1, &bsurf, &copyOptions, &result);
//        PARASOLID_ERROR(PK_GEOM_copy, err);
//        PK_BSURF_t copyBsurf = *(result.copied_geoms);


//        PK_ASSEMBLY_t part;
//        PK_ASSEMBLY_create_empty(&part);
//        err = PK_PART_add_geoms(part, 1, &copyBsurf);
//        PARASOLID_ERROR(PK_PART_add_geoms, err);

//        PK_PART_transmit_o_t transmit_options;
//        PK_PART_transmit_o_m(transmit_options);
//        transmit_options.transmit_format = PK_transmit_format_text_c;


//        std::string name = "level_" + internal::to_string(level) +
//                       "_box_" + internal::to_string(box);

//        PK_ERROR_code_t err = PK_PART_transmit(1, &part, name.c_str(), &transmit_options);
//        PARASOLID_ERROR(PK_PART_transmit, err);

        // ----------------------------------------------------------------------

	// Extra transposition (2019-11-14) to compensate for the transposition in createPK_BSURF.
	// PK_BSURF_reparameterise_o_t rep_options;
	// PK_BSURF_reparameterise_o_m(rep_options);
	// rep_options.transpose = PK_LOGICAL_true;

	// err = PK_BSURF_reparameterise(bsurf,&rep_options);
	// PARASOLID_ERROR(PK_BSURF_reparameterise, err);
	// End of the transposition.


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
        // set the name of the body to identify it
        std::string name = "BoxID_" + internal::to_string(box);
        PK_ATTDEF_t attdef;
        PK_ATTRIB_t attribut;
        PK_ATTDEF_find("SDL/TYSA_NAME",&attdef);
        if(attdef == PK_ENTITY_null)
            std::cout<<"entity null for tysa name" <<std::endl;
        PK_ERROR_code_t err2 = PK_ATTRIB_create_empty(trimSurface,attdef,&attribut);
        if(err2 == PK_ERROR_existing_attrib || err2 == PK_ERROR_wrong_entity)
            std::cout<<"attribute already exists" <<std::endl;
        PK_ATTRIB_set_string(attribut,0,name.c_str());
        if (state != PK_BODY_state_ok_c)
        {
            std::cout << "Something went wrong. state("  << state << ")" << std::endl;
        }

        PK_INSTANCE_sf_t sform;
        sform.assembly = assembly;
        sform.transf = PK_ENTITY_null;
        sform.part = trimSurface;

        PK_INSTANCE_t instance;
        err = PK_INSTANCE_create(&sform, &instance);
        PARASOLID_ERROR(PK_INSTANCE_create, err);

    }

    err = PK_SESSION_set_check_continuity(PK_LOGICAL_true);
    PARASOLID_ERROR(PK_SESSION_set_check_continuity, err);

    return true;
}

template <class T>
bool getParBoxAsIndexBoxInLevel(const gsTHBSplineBasis<2, T>& basis,unsigned lvl,const std::vector<real_t>& par_box,
                                std::vector<index_t>& index_box)
{
    T lowU=par_box[0];
    T lowV=par_box[1];
    T upU=par_box[2];
    T upV=par_box[3];
    unsigned lowerIndexU,upperIndexU,lowerIndexV,upperIndexV;
    const typename gsBSplineTraits<2,T>::Basis & tBasis = *(basis.getBases()[lvl]);
    const gsKnotVector<T> & tKvU = tBasis.component(0).knots();
    if( !( *(tKvU.begin())<=lowU && lowU<=upU &&
           upU<=tKvU.get().end()[-tKvU.degree()-1]) )
    {
        if( !(lowU<=upU) )
        {
            std::cout << "Error in the box." << std::endl;
            return false;
        }
        std::cout<<"Box in u is not in Range, it will be cut at domain boundary."<<std::endl;
        if( !( *(tKvU.begin())<=lowU) )
            lowU=*(tKvU.begin());
        if( !(upU<=tKvU.get().end()[-tKvU.degree()-1]) )
            upU=tKvU.get().end()[-tKvU.degree()-1];
    }
    lowerIndexU=(tKvU.uFind(lowU)-tKvU.uFind(0));
    upperIndexU=(tKvU.uFind(upU)-tKvU.uFind(0))+1;

    const gsKnotVector<T> & tKvV = tBasis.component(1).knots();
    if( !( *(tKvV.begin())<=lowV && lowV<=upV &&
           upV<=tKvV.get().end()[-tKvV.degree()-1]) )
    {
        if( !(lowV<=upV) )
        {
            std::cout << "Error in the box." << std::endl;
            return false;
        }
        std::cout<<"Box in v is not in Range, it will be cut at domain boundary."<<std::endl;
        if( !( *(tKvV.begin())<=lowV) )
            lowV=*(tKvV.begin());
        if( !(upV<=tKvV.get().end()[-tKvV.degree()-1]) )
            upV=tKvV.get().end()[-tKvV.degree()-1];
    }
    lowerIndexV=(tKvV.uFind(lowV)-tKvV.uFind(0));
    upperIndexV=(tKvV.uFind(upV)-tKvV.uFind(0))+1;
    index_box.clear();
    index_box.push_back(lvl);
    index_box.push_back(lowerIndexU);
    index_box.push_back(lowerIndexV);
    index_box.push_back(upperIndexU);
    index_box.push_back(upperIndexV);
    return true;
}

template <class T>
bool parBoxesIntersect(const std::vector<T>& par_boxes)
{
    T lowerUi,lowerVi,upperUi,upperVi,lowerUj,lowerVj,upperUj,upperVj;
    unsigned d=2;
    unsigned boxSize=2*d;
    for(unsigned i = 0;i<par_boxes.size();i+=boxSize)
    {
        lowerUi=par_boxes[i];
        lowerVi=par_boxes[i+1];
        upperUi=par_boxes[i+2];
        upperVi=par_boxes[i+3];
        for(unsigned j = 0;j<i;j+=boxSize)
        {
            lowerUj=par_boxes[j];
            lowerVj=par_boxes[j+1];
            upperUj=par_boxes[j+2];
            upperVj=par_boxes[j+3];

            if(upperUi>lowerUj&&lowerUi<upperUj&&
                    upperVi>lowerVj&&lowerVi<upperVj)
                return true;
        }
    }
    return false;
}

template <class T>
bool getTrimDataFromBoundaryTrimCurves(gsTHBSplineBasis<2>::AxisAlignedBoundingBox boundaryAABB,
gsTHBSplineBasis<2>::TrimmingCurves trimCurves,std::vector<gsTrimData<T> >& trimdata)
{
    trimdata.clear();
    for(unsigned level = 0;level<trimCurves.size();++level)
    {
        for(unsigned component=0;component<trimCurves[level].size();++component)
        {
            gsTrimData<T> td(level,boundaryAABB[level][component],trimCurves[level][component]);
            trimdata.push_back(td);
        }
    }
    return true;
}

template <class T>
bool getTrimCurvesAndBoundingBoxes(const gsTHBSpline<2, T>& surface,
                                   const std::vector<T>& par_boxes,
                                   std::vector<gsTrimData<T> >& trimdata)
{
    if(parBoxesIntersect(par_boxes))
    {
        std::cout<<"Given boxes have intersections."<<std::endl;
        return false;
    }

    const gsTHBSplineBasis<2, T>* basis = static_cast< const gsTHBSplineBasis<2,T>* > (&surface.basis());
    unsigned maxLevel = basis->tree().getMaxInsLevel();
    unsigned lvl;
    std::vector<std::vector<index_t> >aabbBoxesForCheck;
    std::vector<T> par_box;
    std::vector<index_t> index_box;
    unsigned d=2;
    unsigned boxSize=2*d;
    bool success;
    trimdata.clear();
    for(unsigned i = 0;i<par_boxes.size();i=i+boxSize)
    {
        par_box.clear();
        par_box.push_back(par_boxes[i]);
        par_box.push_back(par_boxes[i+1]);
        par_box.push_back(par_boxes[i+2]);
        par_box.push_back(par_boxes[i+3]);
        success=getParBoxAsIndexBoxInLevel(*basis,maxLevel,par_box,index_box);

        if(!success)
            return false;

        gsVector<index_t,2>lower;
        lower << index_box[1],index_box[2];
        gsVector<index_t,2>upper;
        upper << index_box[3],index_box[4];
        lvl=basis->tree().query4(lower, upper, index_box[0]);

        success=getParBoxAsIndexBoxInLevel(*basis,lvl,par_box,index_box);
        if(!success)
            return false;

        //std::cout << "par-box: "<< par_box[0] << " "<< par_box[1] << " " << par_box[2] << " " << par_box[3] << std::endl;
        //std::cout << "el-box: "<< index_box[1] << " "<< index_box[2] << " "
        //          << index_box[3] << " " << index_box[4] << " in level: " << lvl << std::endl;

        std::vector<index_t> aabb_box;
        aabb_box.push_back(index_box[1]<<(maxLevel-lvl));
        aabb_box.push_back(index_box[2]<<(maxLevel-lvl));
        aabb_box.push_back(index_box[3]<<(maxLevel-lvl));
        aabb_box.push_back(index_box[4]<<(maxLevel-lvl));

        std::vector<std::vector<std::vector<T> > > trimCurveComp;
        std::vector<std::vector<T> >trimCurveBox;
        std::vector<T> trimCurve;
        trimCurve.push_back(par_box[0]);
        trimCurve.push_back(par_box[1]);
        trimCurve.push_back(par_box[0]);
        trimCurve.push_back(par_box[3]);
        trimCurveBox.push_back(trimCurve);
        trimCurve.clear();
        trimCurve.push_back(par_box[0]);
        trimCurve.push_back(par_box[3]);
        trimCurve.push_back(par_box[2]);
        trimCurve.push_back(par_box[3]);
        trimCurveBox.push_back(trimCurve);
        trimCurve.clear();
        trimCurve.push_back(par_box[2]);
        trimCurve.push_back(par_box[1]);
        trimCurve.push_back(par_box[2]);
        trimCurve.push_back(par_box[3]);
        trimCurveBox.push_back(trimCurve);
        trimCurve.clear();
        trimCurve.push_back(par_box[0]);
        trimCurve.push_back(par_box[1]);
        trimCurve.push_back(par_box[2]);
        trimCurve.push_back(par_box[1]);
        trimCurveBox.push_back(trimCurve);
        trimCurveComp.push_back(trimCurveBox);

        gsTrimData<T> td(lvl,aabb_box,trimCurveComp);
        trimdata.push_back(td);

        aabbBoxesForCheck.push_back(aabb_box);
    }
    return true;
}



}//extensions

}//gismo

#undef PARASOLID_ERROR
