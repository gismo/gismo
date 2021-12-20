/** @file gsPyBind11.cpp

    @brief PyBind11 main module file

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#include <gsCore/gsConfig.h>
#include <gsCore/gsExport.h>

#include <gismo.h>

#ifdef GISMO_BUILD_PYBIND11

namespace py = pybind11;

/**
   @brief Creates G+Smo Python module
*/
PYBIND11_MODULE(pygismo, m) {

  m.attr("__name__") = "pygismo";
  m.attr("__version__") = GISMO_VERSION;
  m.doc() = "G+Smo (Geometry + Simulation Modules)";

  py::add_ostream_redirect(m, "ostream_redirect");
  
  py::module assembler = m.def_submodule("assembler");

  assembler.attr("__name__") = "pygismo.assembler";
  assembler.attr("__version__") = GISMO_VERSION;
  assembler.doc() = "G+Smo (Geometry + Simulation Modules): Assembler module";

  py::module core = m.def_submodule("core");

  core.attr("__name__") = "pygismo.core";
  core.attr("__version__") = GISMO_VERSION;
  core.doc() = "G+Smo (Geometry + Simulation Modules): Core module";
  
  gismo::pybind11_init_gsFunction( core );
  gismo::pybind11_init_gsFunctionExpr( core );
  gismo::pybind11_init_gsMultiPatch( core );
  gismo::pybind11_init_gsMultiBasis( core );

  gismo::pybind11_enum_gsBoundary( core );

  py::module hsplines = m.def_submodule("hsplines");

  hsplines.attr("__name__") = "pygismo.hspline";
  hsplines.attr("__version__") = GISMO_VERSION;
  hsplines.doc() = "G+Smo (Geometry + Simulation Modules): HSplines module";

  gismo::pybind11_init_gsHBSplineBasis2( core );
  gismo::pybind11_init_gsHBSplineBasis3( core );
  gismo::pybind11_init_gsHBSplineBasis4( core );
  gismo::pybind11_init_gsHBSpline2( core );
  gismo::pybind11_init_gsHBSpline3( core );
  gismo::pybind11_init_gsHBSpline4( core );
  gismo::pybind11_init_gsTHBSplineBasis2( core );
  gismo::pybind11_init_gsTHBSplineBasis3( core );
  gismo::pybind11_init_gsTHBSplineBasis4( core );
  gismo::pybind11_init_gsTHBSpline2( core );
  gismo::pybind11_init_gsTHBSpline3( core );
  gismo::pybind11_init_gsTHBSpline4( core );
  
  py::module io = m.def_submodule("io");

  io.attr("__name__") = "pygismo.io";
  io.attr("__version__") = GISMO_VERSION;
  io.doc() = "G+Smo (Geometry + Simulation Modules): IO module";

  gismo::pybind11_init_gsCmdLine( io );
  gismo::pybind11_init_gsFileData( io );
  gismo::pybind11_init_gsOptionList (io );  

  py::module matrix = m.def_submodule("matrix");

  matrix.attr("__name__") = "pygismo.matrix";
  matrix.attr("__version__") = GISMO_VERSION;
  matrix.doc() = "G+Smo (Geometry + Simulation Modules): Matrix module";
  
  gismo::pybind11_init_gsMatrix<real_t>(matrix,"Real"); //gsMatrixReal
  gismo::pybind11_init_gsMatrix<index_t>(matrix,"Int"); //gsMatrixInt
  gismo::pybind11_init_gsSparseMatrix<real_t>(matrix,"Real"); //gsSparseMatrixReal
  gismo::pybind11_init_gsSparseMatrix<index_t>(matrix,"Int"); //gsSparseMatrixInt
  
  py::module modelling = m.def_submodule("modelling");

  modelling.attr("__name__") = "pygismo.modelling";
  modelling.attr("__version__") = GISMO_VERSION;
  modelling.doc() = "G+Smo (Geometry + Simulation Modules): Modelling module";

  py::module mpi = m.def_submodule("mpi");
  
  mpi.attr("__name__") = "pygismo.mpi";
  mpi.attr("__version__") = GISMO_VERSION;
  mpi.doc() = "G+Smo (Geometry + Simulation Modules): MPI module";

  py::module multigrid = m.def_submodule("multigrid");
  
  multigrid.attr("__name__") = "pygismo.multigrid";
  multigrid.attr("__version__") = GISMO_VERSION;
  multigrid.doc() = "G+Smo (Geometry + Simulation Modules): MultiGrid module";

  py::module nurbs = m.def_submodule("nurbs");

  nurbs.attr("__name__") = "pygismo.nurbs";
  nurbs.attr("__version__") = GISMO_VERSION;
  nurbs.doc() = "G+Smo (Geometry + Simulation Modules): NURBS module";

  gismo::pybind11_init_gsKnotVector( nurbs );
  gismo::pybind11_init_gsBSpline( nurbs );
  gismo::pybind11_init_gsBSplineBasis( nurbs );
  gismo::pybind11_init_gsTensorBSpline2( nurbs );
  gismo::pybind11_init_gsTensorBSpline3( nurbs );
  gismo::pybind11_init_gsTensorBSpline4( nurbs );
  gismo::pybind11_init_gsTensorBSplineBasis2( nurbs );
  gismo::pybind11_init_gsTensorBSplineBasis3( nurbs );
  gismo::pybind11_init_gsTensorBSplineBasis4( nurbs );

  
  py::module pde = m.def_submodule("pde");

  pde.attr("__name__") = "pygismo.pde";
  pde.attr("__version__") = GISMO_VERSION;
  pde.doc() = "G+Smo (Geometry + Simulation Modules): Pde module";

  gismo::pybind11_enum_gsBoundaryConditions( core );
  gismo::pybind11_init_gsBoundaryConditions( pde );

  py::module solver = m.def_submodule("solver");

  solver.attr("__name__") = "pygismo.solver";
  solver.attr("__version__") = GISMO_VERSION;
  solver.doc() = "G+Smo (Geometry + Simulation Modules): Solver module";

  py::module tensor = m.def_submodule("tensor");

  tensor.attr("__name__") = "pygismo.tensor";
  tensor.attr("__version__") = GISMO_VERSION;
  tensor.doc() = "G+Smo (Geometry + Simulation Modules): Tensor module";

  py::module utils = m.def_submodule("utils");

  utils.attr("__name__") = "pygismo.utils";
  utils.attr("__version__") = GISMO_VERSION;
  utils.doc() = "G+Smo (Geometry + Simulation Modules): Utils module";
}

#endif // GISMO_BUILD_PYBIND11
