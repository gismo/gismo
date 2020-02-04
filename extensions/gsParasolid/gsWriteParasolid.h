/** @file gsWriteParasolid.h

    @brief Provides declaration of gsWriteParasolid functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

/*
  (BSURF_sf 
  :u_degree          
  :v_degree          
  :n_u_vertices      
  :n_v_vertices      
  :vertex_dim        
  :is_rational       nil
  :vertex            (vertex-sf cp0)
  :n_u_knots         
  :u_knot            
  :u_knot_mult       
  :n_v_knots         
  :v_knot            
  :v_knot_mult       
*/


#pragma once

#include <string>

#include <gsCore/gsForwardDeclarations.h>
#include <gsHSplines/gsTHBSplineBasis.h>

#include <gsParasolid/gsPKSession.h>

typedef int PK_GEOM_t;
typedef int PK_BSURF_t;
typedef int PK_BCURVE_t;
typedef int PK_BODY_t;
typedef int PK_ASSEMBLY_t;
struct PK_UVBOX_s;

namespace gismo {

namespace extensions {

    template<class T>
    class gsTrimData;
    
    /// Writes a gsSurface to a parasolid file
    /// \param gssurf a surface
    /// \param fname filename (without extension)
    template<class T>
    bool gsWriteParasolid( const gsGeometry<T> & ggeo, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsMultiPatch<T> & gssurfs, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsMesh<T>& mesh, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsTHBSpline<2, T>& thb, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsTHBSpline<2, T>& thb, const std::vector<T>&par_boxes, std::string const & filename );

    /// Converts \a tp into a PK_SHEET and writes it to filename.xmt_txt.
    template<class T>
    bool gsWritePK_SHEET(const gsTensorBSpline<2, T>& tp, const std::string& filename);


    /// Translates a gsTensorBSpline to a PK_BSURF_t
    /// \param[in] bsp B-spline surface
    /// \param[out] bsurf Parasolid spline surface
    template<class T>
    bool createPK_BSURF( const gsTensorBSpline< 2,T> & bsp, PK_BSURF_t & bsurf,
			 bool closed_u = false, bool closed_v = false );

    /// Translates a gsBSpline to a PK_BCURVE_t
    /// \param[in] curve B-Spline surve
    /// \param[out] bcurve Parasolid spline curve
    template<class T> 
    bool createPK_BCURVE( const gsBSpline<T>& curve, PK_BCURVE_t& bcurve );

    /// Translates a gsGeometry to a PK_GEOM_t
    /// \param[in] ggeo inpute G+SMO geometry
    /// \param[out] pgeo Parasolid geometric entity
    template<class T> 
    bool createPK_GEOM( const gsGeometry<T> & ggeo, PK_GEOM_t & pgeo );

    /// Translates a gsMesh to PK_BODY_t
    /// \param[in] mesh input G+Smo mesh
    /// \param[out] body Parasolid wire body
    template<class T> 
    bool exportMesh( const gsMesh<T>& mesh, PK_BODY_t& body );

    /// Translates a THB-Spline surface to PK_BODY_t
    /// \param[in] surface THB-Spline surface
    /// \param[out] body Parasolid body
    template<class T> 
    bool exportTHBsurface( const gsTHBSpline<2, T>& surface, PK_ASSEMBLY_t& body );

    /// Translates a THB-Spline surface to PK_BODY_t
    /// \param[in] surface THB-Spline surface
    /// \param[in] boxes which give the splitting of the domain
    /// \param[out] body Parasolid body
    template<class T>
    bool exportTHBsurface( const gsTHBSpline<2, T>& surface, const std::vector<T>& par_boxes, PK_ASSEMBLY_t& body );

    template <class T>
    bool getTrimCurvesAndBoundingBoxes(const gsTHBSpline<2, T>& surface,
                                       const std::vector<T>& par_boxes,
                                       std::vector<gsTrimData<T> >& trimdata);

    /// Translates a box in the parameter space to a box in index space of the given lvl, which
    /// contains the parameter box
    /// \param[in] thb spline basis
    /// \param[in] lvl
    /// \param[in] the parameter box of size 4 [lowU,lowV,upU,upV]
    /// \param[out] the index box of size 5 [lvl,lowIndexU,lowIndexV,upIndexU,upIndexV]
    template <class T>
    bool getParBoxAsIndexBoxInLevel(const gsTHBSplineBasis<2, T>& basis,unsigned lvl,
                                    const std::vector<real_t>& par_box,std::vector<unsigned>& index_box);

    /// Checks if there are intersections in a given vector of par_boxes
    /// \param[in] vector of real_t, always 4 elements describe a box: [lowU,lowV,upU,upV]
    template <class T>
    bool parBoxesIntersect(const std::vector<T>& par_boxes);


}//extensions

}//gismo

