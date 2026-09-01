#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsWriteParaview.hpp>
#include <gsMesh2/gsSurfMesh.h>
// #include <gsCore/gsMultiPatch.h>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{

TEMPLATE_INST
void gsWriteParaview(const gsField<T> & field, std::string const & fn,
                     unsigned npts, bool mesh, const std::string pDelim, bool skipPvd);

TEMPLATE_INST
void gsWriteParaview(const gsGeometry<T> & Geo, std::string const & fn,
                     unsigned npts, bool mesh, bool ctrlNet, bool skipPvd);

TEMPLATE_INST
void gsWriteParaview( std::vector<gsGeometry<T> *> const & Geo, std::string const & fn,
                      unsigned npts, bool mesh, bool ctrlNet, const std::string pDelim);

TEMPLATE_INST
void gsWriteParaviewBezier(const gsMultiPatch<real_t> & mPatch, std::string const & filename, bool ctrlNet);

TEMPLATE_INST
void gsWriteParaviewUnstructuredGrid(const gsMultiPatch<T> & mPatch, std::string const & fn,
                                     unsigned npts, bool export_base64, bool skipPvd);

TEMPLATE_INST
void gsWriteParaviewUnstructuredGrid(const gsField<T> & field, std::string const & fn,
                                     unsigned npts, bool export_base64, bool skipPvd);

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
void gsWriteParaview(gsFunctionSet<T> const& geo, gsFunctionSet<T> const& func,std::string const & fn, unsigned npts, const std::string pDelim, bool skipPvd);

TEMPLATE_INST
void gsWriteParaview(const gsMultiPatch<T> & Geo, std::string const & fn, 
                     unsigned npts, bool mesh, bool ctrlNet, const std::string pDelim);

TEMPLATE_INST
void gsWriteParaview(gsMappedSpline<2,T> const& mspline, std::string const & fn,unsigned npts, bool skipPvd);

TEMPLATE_INST
void gsWriteParaview(gsFunctionSet<T> const& geom, gsMappedBasis<2,T>  const& mbasis,std::string const & fn,unsigned npts,const bool fullsupport, const std::vector<index_t> indices, bool skipPvd);

TEMPLATE_INST
void gsWriteParaview(gsMultiPatch<T> const& mp, gsMultiBasis<T> const& mb,std::string const & fn, unsigned npts, bool skipPvd);

TEMPLATE_INST
void gsWriteParaview(gsFunction<T> const& func, gsMatrix<T> const& supp, std::string const & fn, unsigned npts, bool graph);

TEMPLATE_INST
void gsWriteParaview(gsBasis<T> const& basis, std::string const & fn,
                     unsigned npts, bool mesh);

TEMPLATE_INST
void gsWriteParaview(gsBasis<T> const& basis,
                     std::vector<index_t> const & indices,
                     std::string const & fn,
                     unsigned npts, bool mesh);

TEMPLATE_INST
void gsWriteParaview(const gsMatrix<T> & box, std::string const & fn, const gsVector<T> & values);
TEMPLATE_INST
void gsWriteParaview(const gsMatrix<T> & box, std::string const & fn, const std::vector<T> & values);
TEMPLATE_INST
void gsWriteParaview(const gsMatrix<T> & box, std::string const & fn, T value);

TEMPLATE_INST
void gsWriteParaview(const gsHBox<2,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void gsWriteParaview(const gsHBox<3,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void gsWriteParaview(const gsHBoxContainer<2,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void gsWriteParaview(const gsHBoxContainer<3,T> & hbox, std::string const & fn, short_t mode);

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
void gsWriteParaview(const gsSurfMesh<T> & sm,
                     std::string const & fn);

TEMPLATE_INST
void gsWriteParaview(const gsSurfMesh<T> & sm,
                     std::string const & fn,
                     std::vector<std::string> props);

TEMPLATE_INST
void gsWriteHalfedgesParaview(const gsSurfMesh<T>& sm, 
                                    const std::string& fn,
                                    real_t eps);                     

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
void writeSingleBasisMesh(const gsBasis<T> & basis, std::string const & fn);

TEMPLATE_INST
void writeSingleCompMesh(const gsBasis<T> & basis, const gsGeometry<T> & Geo,
                         std::string const & fn, unsigned resolution);

TEMPLATE_INST
void writeSingleBox(const gsMatrix<T> & box, std::string const & fn, T value);

TEMPLATE_INST
void writeSingleHBox(const gsHBox<2,T> & box, std::string const & fn);

TEMPLATE_INST
void writeSingleHBox(const gsHBox<3,T> & box, std::string const & fn);

TEMPLATE_INST
void writeSingleControlNet(const gsGeometry<T> & Geo, std::string const & fn);

///////////////////////////////////////////////////////////////////////

#ifdef GISMO_WITH_PYBIND11

namespace py = pybind11;

void pybind11_init_gsWriteParaview(py::module &m)
{
    m.def("gsWriteParaview",
            [](const gsField<real_t> &field, const std::string &filename, unsigned numPoints, bool mesh, const std::string &pDelim)
            { gsWriteParaview(field, filename, numPoints, mesh, pDelim); },
            "Writes a field to Paraview",
            py::arg("field"),
            py::arg("filename"),
            py::arg("numPoints")=1000,
            py::arg("mesh")=false,
            py::arg("pDelim")="");
    m.def("gsWriteParaview",
            [](const gsGeometry<real_t> &geometry, const std::string &filename, unsigned numPoints, bool mesh, bool ctrlNet, bool skipPvd)
            { gsWriteParaview(geometry, filename, numPoints, mesh, ctrlNet, skipPvd); },
            "Writes a geometry to Paraview",
            py::arg("geometry"),
            py::arg("filename"),
            py::arg("numPoints")=1000,
            py::arg("mesh")=false,
            py::arg("ctrlNet")=false,
            py::arg("skipPvd")=false);
    m.def("gsWriteParaview",
            [](const gsBasis<real_t> &basis, const std::string &filename, unsigned numPoints, bool mesh)
            { gsWriteParaview(basis, filename, numPoints, mesh); },
            "Writes a basis to Paraview",
            py::arg("basis"),
            py::arg("filename"),
            py::arg("numPoints")=1000,
            py::arg("mesh")=false);
    m.def("gsWriteParaview",
            [](const gsFunctionSet<real_t> &function, const std::string &filename, unsigned numPoints)
            { gsWriteParaview(function, filename, numPoints); },
            "Writes a geometry to Paraview",
            py::arg("function"),
            py::arg("filename"),
            py::arg("numPoints")=1000);
    m.def("gsWriteParaview",
            [](const gsFunctionSet<real_t> &geometry, const gsFunctionSet<real_t> &function, const std::string &filename, unsigned numPoints, const std::string &pDelim, bool skipPvd)
            { gsWriteParaview(geometry, function, filename, numPoints, pDelim, skipPvd); },
            "Writes a geometry to Paraview",
            py::arg("geometry"),
            py::arg("function"),
            py::arg("filename"),
            py::arg("numPoints")=1000,
            py::arg("pDelim")="",
            py::arg("skipPvd")=false);
    m.def("gsWriteParaview",
            [](const gsMultiPatch<real_t> &geometry, const std::string &filename, unsigned numPoints, bool mesh, bool ctrlNet, const std::string &pDelim)
            { gsWriteParaview(geometry, filename, numPoints, mesh, ctrlNet, pDelim); },
            "Writes a multi-patch to paraview",
            py::arg("geometry"),
            py::arg("filename"),
            py::arg("numPoints")=1000,
            py::arg("mesh")=false,
            py::arg("controlNet")=false,
            py::arg("delimiter")="_");

    m.def("gsWriteParaviewPoints",
            [](const gsMatrix<real_t> &X, const gsMatrix<real_t> &Y, const std::string &fn)
            { gsWriteParaviewPoints(X, Y, fn); },
            "Writes points to Paraview");
    m.def("gsWriteParaviewPoints",
            [](const gsMatrix<real_t> &X, const gsMatrix<real_t> &Y, const gsMatrix<real_t> &ZMat, const std::string &fn)
            { gsWriteParaviewPoints(X, Y, ZMat, fn); },
            "Writes points to Paraview");
    m.def("gsWriteParaviewPoints",
            [](const gsMatrix<real_t> &X, const gsMatrix<real_t> &Y, const gsMatrix<real_t> &ZMat, const gsMatrix<real_t> &V, const std::string &fn)
            { gsWriteParaviewPoints(X, Y, ZMat, V, fn); },
            "Writes points to Paraview");

}


#endif


} // namespace gismo

#undef T
#undef uZ
#undef Z
