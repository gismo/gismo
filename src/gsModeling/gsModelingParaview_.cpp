#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsModelingParaview.h>
#include <gsModeling/gsModelingParaview.hpp>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{


TEMPLATE_INST
void gsWriteParaview(gsSolid<T> const& sl, std::string const & fn, unsigned numPoints_for_eachCurve, int vol_Num,
                     T edgeThick, gsVector3d<T> const & translate, int color_convex,
                     int color_nonconvex, int color_eloop, std::vector<unsigned> const & eloop);

TEMPLATE_INST
void gsWriteParaviewSolid(gsSolid<T> const  & sl,
                     std::string const & fn,
                     unsigned numSamples );

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

} // namespace gismo

#undef T
#undef uZ
#undef Z
