#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsWriteParaview.hpp>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{
  
TEMPLATE_INST
void gsWriteParaview(const gsField<T> & field, std::string const & fn, 
                     unsigned npts, bool mesh);

TEMPLATE_INST
void gsWriteParaview(const gsGeometry<T> & Geo, std::string const & fn, 
                     unsigned npts, bool mesh, bool ctrlNet);

TEMPLATE_INST
void gsWriteParaview( std::vector<gsGeometry<T> *> const & Geo, std::string const & fn, 
                      unsigned npts, bool mesh, bool ctrlNet);

TEMPLATE_INST
void gsWriteParaview(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & domain,
                     std::string const & fn, unsigned npts);

TEMPLATE_INST
void gsWriteParaview_basisFnct(int i, gsBasis<T> const& basis, std::string const & fn, 
                               unsigned npts );

TEMPLATE_INST
void gsWriteParaview(gsGeometrySlice<T> const& Geo, std::string const & fn, unsigned npts );

TEMPLATE_INST
void gsWriteParaview(gsFunction<T> const& func, gsMatrix<T> const& supp, std::string const & fn, unsigned npts );

TEMPLATE_INST
void gsWriteParaview(gsBasis<T> const& basis, std::string const & fn, 
                     unsigned npts, bool mesh);

TEMPLATE_INST
void gsWriteParaviewPoints(gsMatrix<T> const& X, gsMatrix<T> const& Y, std::string const & fn);

TEMPLATE_INST
void gsWriteParaviewPoints(gsMatrix<T> const& X, gsMatrix<T> const& Y, gsMatrix<T> const& z, std::string const & fn);

TEMPLATE_INST
void gsWriteParaviewPoints(gsMatrix<T> const& X, gsMatrix<T> const& Y, gsMatrix<T> const& z, gsMatrix<T> const& v, std::string const & fn);

TEMPLATE_INST
void gsWriteParaviewPoints(gsMatrix<T> const& points, std::string const & fn);

TEMPLATE_INST
void gsWriteParaviewTPgrid(gsMatrix<T> const& points,
                           gsMatrix<T> const& data,
                           const gsVector<index_t> & np,
                           std::string const & fn);

TEMPLATE_INST
void gsWriteParaview(gsSolid<T> const& sl, std::string const & fn, unsigned numPoints_for_eachCurve, int vol_Num,
                     T edgeThick, gsVector3d<T> const & translate, int color_convex,
                     int color_nonconvex, int color_eloop, std::vector<unsigned> const & eloop);

TEMPLATE_INST
void gsWriteParaviewSolid(gsSolid<T> const  & sl, 
                     std::string const & fn, 
                     unsigned numSamples );

TEMPLATE_INST
void gsWriteParaview(gsMesh<T> const& sl, std::string const & fn, bool pvd);

TEMPLATE_INST
void gsWriteParaview(const std::vector<gsMesh<T> >& sl, std::string const & fn);

//TEMPLATE_INST
//void gsWriteParaview(gsHeMesh<T> const& sl, std::string const & fn);

TEMPLATE_INST
void gsWriteParaview(gsPlanarDomain<T> const & pdomain, 
                     std::string const & fn, unsigned npts);

TEMPLATE_INST
void gsWriteParaview(const gsTrimSurface<T> & ts, std::string const & fn, 
                     unsigned npts, bool trimCurves);

TEMPLATE_INST
void gsWriteParaview(const gsVolumeBlock<T>& volBlock,
                     std::string const& fn,
                     unsigned npts);

TEMPLATE_INST
void gsWriteParaviewTrimmedCurve(const gsTrimSurface<T>& surf,
                                 const unsigned idLoop,
                                 const unsigned idCurve,
                                 const std::string fn,
                                 unsigned npts);

TEMPLATE_INST
void writeSinglePatchField(const gsFunction<T> & geometry,
                           const gsFunction<T> & parField,
                           const bool isParam,
                           std::string const & fn, unsigned npts);


} // namespace gismo

#undef T
#undef uZ
#undef Z
