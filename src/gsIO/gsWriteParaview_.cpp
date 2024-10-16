#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsWriteParaview.hpp>
// #include <gsCore/gsMultiPatch.h>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{
  
TEMPLATE_INST
void gsWriteParaview(const gsField<T> & field, std::string const & fn, 
                     unsigned npts, bool mesh, const std::string pDelim);

TEMPLATE_INST
void gsWriteParaview(const gsGeometry<T> & Geo, std::string const & fn, 
                     unsigned npts, bool mesh, bool ctrlNet);

TEMPLATE_INST
void gsWriteParaview( std::vector<gsGeometry<T> *> const & Geo, std::string const & fn, 
                      unsigned npts, bool mesh, bool ctrlNet, const std::string pDelim);

TEMPLATE_INST
void gsWriteParaview(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & domain,
                     std::string const & fn, unsigned npts);

TEMPLATE_INST
void gsWriteParaview_basisFnct(int i, gsBasis<T> const& basis, std::string const & fn, 
                               unsigned npts );

TEMPLATE_INST
void gsWriteParaview(gsGeometrySlice<T> const& Geo, std::string const & fn, unsigned npts );

TEMPLATE_INST
void gsWriteParaview(gsFunctionSet<T> const& func, std::string const & fn, unsigned npts);

TEMPLATE_INST
void gsWriteParaview(gsFunctionSet<T> const& geo, gsFunctionSet<T> const& func,std::string const & fn, unsigned npts, const std::string pDelim);

TEMPLATE_INST
void gsWriteParaview(gsMappedSpline<2,T> const& mspline, std::string const & fn,unsigned npts);

TEMPLATE_INST
void gsWriteParaview(gsFunctionSet<T> const& geom, gsMappedBasis<2,T>  const& mbasis,std::string const & fn,unsigned npts,const bool fullsupport, const std::vector<index_t> indices);

TEMPLATE_INST
void gsWriteParaview(gsMultiPatch<T> const& mp, gsMultiBasis<T> const& mb,std::string const & fn, unsigned npts);

TEMPLATE_INST
void gsWriteParaview(gsFunction<T> const& func, gsMatrix<T> const& supp, std::string const & fn, unsigned npts, bool graph);

TEMPLATE_INST
void gsWriteParaview(gsBasis<T> const& basis, std::string const & fn, 
                     unsigned npts, bool mesh);

TEMPLATE_INST
void gsWriteParaview(gsHBox<2,T> & hbox, std::string const & fn);

TEMPLATE_INST
void gsWriteParaview(gsHBoxContainer<2,T> & hbox, std::string const & fn);

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
void gsWriteParaview(gsMesh<T> const& sl, std::string const & fn, const gsMatrix<T>& params);

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
void gsWriteParaviewBdr(gsMultiPatch<T> const & patches,
                     std::string const & fn,
                     unsigned npts, bool ctrlNet);

TEMPLATE_INST
void gsWriteParaviewIfc(gsMultiPatch<T> const & patches,
                     std::string const & fn,
                     unsigned npts, bool ctrlNet);

TEMPLATE_INST
void gsWriteParaview(gsMultiPatch<T> const & patches,
                     typename gsBoundaryConditions<T>::bcContainer const & bcs,
                     std::string const & fn,
                     unsigned npts, bool ctrlNet);

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

TEMPLATE_INST
void writeSingleCompMesh(const gsBasis<T> & basis, const gsGeometry<T> & Geo,
                         std::string const & fn, unsigned resolution);

TEMPLATE_INST
void writeSingleHBox(gsHBox<2,T> & box, std::string const & fn);

TEMPLATE_INST
void writeSingleControlNet(const gsGeometry<T> & Geo, std::string const & fn);

///////////////////////////////////////////////////////////////////////

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsWriteParaview(py::module &m)
{
    m.def("gsWriteParaview",
            static_cast<void (*)(const gsGeometry<real_t> &, std::string const &, unsigned, bool, bool)>(&gsWriteParaview),
            "Writes a geometry to Paraview",
            py::arg("Geo"),
            py::arg("fn"),
            py::arg("npts")=1000,
            py::arg("mesh")=false,
            py::arg("ctrlNet")=false);
    m.def("gsWriteParaview",
            static_cast<void (*)(const gsBasis<real_t> &, std::string const &, unsigned, bool)>(&gsWriteParaview),
            "Writes a basis to Paraview",
            py::arg("basis"),
            py::arg("fn"),
            py::arg("npts")=1000,
            py::arg("mesh")=false);
    m.def("gsWriteParaview",
            static_cast<void (*)(const gsFunctionSet<real_t> &, std::string const &, unsigned)>(&gsWriteParaview),
            "Writes a geometry to Paraview",
            py::arg("fun"),
            py::arg("fn"),
            py::arg("npts")=1000);
    m.def("gsWriteParaview",
            static_cast<void (*)(const gsFunctionSet<real_t> &, const gsFunctionSet<real_t> &, std::string const &, unsigned, const std::string)>(&gsWriteParaview),
            "Writes a geometry to Paraview",
            py::arg("geo"),
            py::arg("func"),
            py::arg("fn"),
            py::arg("npts")=1000,
            py::arg("pDelim")="");

    m.def("gsWriteParaviewPoints",
            static_cast<void (*)(const gsMatrix<real_t> &, const gsMatrix<real_t> &, std::string const &)>(&gsWriteParaviewPoints),
            "Writes points to Paraview");
    m.def("gsWriteParaviewPoints",
            static_cast<void (*)(const gsMatrix<real_t> &, const gsMatrix<real_t> &, const gsMatrix<real_t> &,std::string const &)>(&gsWriteParaviewPoints),
            "Writes points to Paraview");
    m.def("gsWriteParaviewPoints",
            static_cast<void (*)(const gsMatrix<real_t> &, const gsMatrix<real_t> &, const gsMatrix<real_t> &, const gsMatrix<real_t> &,std::string const &)>(&gsWriteParaviewPoints),
            "Writes points to Paraview");

}


#endif


} // namespace gismo

#undef T
#undef uZ
#undef Z
