/** @file gsScalePoints.cpp

    @brief Scale point cloud in [0,1]^D

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>


#ifdef gsParasolid_ENABLED
#include <gsParasolid/gsWriteParasolid.h>
#endif

using namespace gismo;

/* input : a parametrized point cloud with parameters and points
  output : interior parameters and corresponding points
           boundary points ordered anticlockwise south-east-north-west.
*/
template<class T>
void sortPointCloud(gsMatrix<T> & parameters,
                    gsMatrix<T> & points,
                    gsMatrix<T> & uv_interiors,
                    gsMatrix<T> & uv_south,
                    gsMatrix<T> & uv_east,
                    gsMatrix<T> & uv_north,
                    gsMatrix<T> & uv_west,
                    gsMatrix<T> & p_interiors,
                    gsMatrix<T> & p_south,
                    gsMatrix<T> & p_east,
                    gsMatrix<T> & p_north,
                    gsMatrix<T> & p_west,
                    std::vector<index_t> & corners)
{
  std::vector<index_t> interiors;
  std::vector<index_t> b_west;
  std::vector<index_t> b_east;
  std::vector<index_t> b_south;
  std::vector<index_t> b_north;

  // Determine the parameter domain by mi/max of parameter values
  T u_min = parameters.row(0).minCoeff(),
    u_max = parameters.row(0).maxCoeff(),
    v_min = parameters.row(1).minCoeff(),
    v_max = parameters.row(1).maxCoeff();

  gsVector<T> curr_point(2,1);
  for(index_t i=0; i < parameters.cols(); i++)
  {
    curr_point = parameters.col(i);
    if( (u_min < curr_point(0)) && (curr_point(0) < u_max) && (v_min < curr_point(1)) && (curr_point(1) < v_max) )
    {
      interiors.push_back(i);
    }
    else // not interior point
    {
      if( (math::abs(curr_point(0) - u_min) < 1e-15) && (curr_point(1) > v_min) )
      {// west edge
        b_west.push_back(i);
      }
      else if( (math::abs(curr_point(0) - u_max) < 1e-15) && curr_point(1) < v_max)
      {// east edge
        b_east.push_back(i);
      }
      else if( (math::abs(curr_point(1) - v_min) < 1e-15) && (curr_point(0) < u_max) )
      {// south edge
        b_south.push_back(i);
      }
      else
      {// north edge
        b_north.push_back(i);
      }
    }
  }

  gsInfo << "There are " << interiors.size() << " interior points.\n";
  corners.push_back(interiors.size()); // c1
  gsDebugVar(interiors.size());
  gsInfo << "There are " << b_south.size() << " south points.\n";
  corners.push_back(interiors.size() + b_south.size()); // c2
  gsDebugVar(interiors.size() + b_south.size());
  gsInfo << "There are " << b_east.size() << " east points.\n";
  corners.push_back(interiors.size() + b_south.size() + b_east.size()); // c3
  gsDebugVar(interiors.size() + b_south.size() + b_east.size());
  gsInfo << "There are " << b_north.size() << " north points.\n";
  corners.push_back(interiors.size() + b_south.size() + b_east.size() + b_north.size());
  gsDebugVar(interiors.size() + b_south.size() + b_east.size() + b_north.size());
  gsInfo << "There are " << b_west.size() << " west points.\n";


  uv_interiors.resize(2, interiors.size());
  p_interiors.resize(3, interiors.size());
  for( index_t i = 0; i < interiors.size(); i++ )
  {
    uv_interiors.col(i) = parameters.col(interiors[i]);
    p_interiors.col(i) = points.col(interiors[i]);
  }

  uv_west.resize(2, b_west.size());
  gsMatrix<T> tmp_west(3, b_west.size());
  for( index_t i = 0; i < b_west.size(); i++ )
  {
    uv_west.col(i) = parameters.col(b_west[i]);
    tmp_west.col(i) = points.col(b_west[i]);
  }

  uv_east.resize(2, b_east.size());
  gsMatrix<T> tmp_east(3, b_east.size());
  for( index_t i = 0; i < b_east.size(); i++ )
  {
    uv_east.col(i) = parameters.col(b_east[i]);
    tmp_east.col(i) = points.col(b_east[i]);
  }

  uv_south.resize(2, b_south.size());
  gsMatrix<T> tmp_south(3, b_south.size());
  for( index_t i = 0; i < b_south.size(); i++ )
  {
    uv_south.col(i) = parameters.col(b_south[i]);
    tmp_south.col(i) = points.col(b_south[i]);
  }

  uv_north.resize(2, b_north.size());
  gsMatrix<T> tmp_north(3, b_north.size());
  for( index_t i = 0; i < b_north.size(); i++ )
  {
    uv_north.col(i) = parameters.col(b_north[i]);
    tmp_north.col(i) = points.col(b_north[i]);
  }

  uv_south.transposeInPlace();
  uv_east.transposeInPlace();
  uv_north.transposeInPlace();
  uv_west.transposeInPlace();

  gsInfo << "------------------------------\n";
  gsInfo << "south:\n";
  std::vector<index_t> tmp = uv_south.idxByColumn(0);
  gsDebugVar(uv_south);
  p_south.resize(tmp_south.rows(), tmp_south.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_south.col(i) = tmp_south.col(tmp[i]);
  }
  uv_south.transposeInPlace();
  gsInfo << "------------------------------\n";

  tmp.clear();
  gsInfo << "east:\n";
  tmp = uv_east.idxByColumn(1);
  gsDebugVar(uv_east);
  p_east.resize(tmp_east.rows(), tmp_east.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_east.col(i) = tmp_east.col(tmp[i]);
  }
  uv_east.transposeInPlace();
  gsInfo << "------------------------------\n";

  tmp.clear();
  gsInfo << "north:\n";
  gsDebugVar(uv_north);
  gsDebugVar(tmp_north);
  gsInfo << "-------- idx by column: --------\n";
  tmp = uv_north.idxByColumn(0);

  for (std::vector<index_t>::iterator it = tmp.begin(); it != tmp.end(); ++it)
     gsInfo << *it <<"\n";
  gsInfo << "----------------\n";


  std::reverse(tmp.begin(),tmp.end());
  gsDebugVar(uv_north);
  gsVector<T> tcol = uv_north.col(0).reverse();
  uv_north.col(0) = tcol;
  tcol = uv_north.col(1).reverse();
  uv_north.col(1) = tcol;
  gsDebugVar(uv_north);
  // for (std::vector<index_t>::iterator it = tmp.begin(); it != tmp.end(); ++it)
  //   gsInfo << *it <<"\n";
  p_north.resize(tmp_north.rows(), tmp_north.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_north.col(i) = tmp_north.col(tmp[i]);
  }

  gsDebugVar(uv_north);
  gsDebugVar(p_north);

  uv_north.transposeInPlace();
  gsInfo << "------------------------------\n";



  tmp.clear();
  gsInfo << "west:\n";
  tmp = uv_west.idxByColumn(1);

  gsDebugVar(uv_west);
  tcol = uv_west.col(0).reverse();
  uv_west.col(0) = tcol;
  tcol = uv_west.col(1).reverse();
  uv_west.col(1) = tcol;
  gsDebugVar(uv_west);
  std::reverse(tmp.begin(),tmp.end());

  p_west.resize(tmp_west.rows(), tmp_west.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_west.col(i) = tmp_west.col(tmp[i]);
  }
  uv_west.transposeInPlace();
  gsInfo << "------------------------------\n";

  gsWriteParaviewPoints(uv_interiors, "uv_interiors");
  gsWriteParaviewPoints(p_interiors, "p_interiors");

  gsWriteParaviewPoints(uv_west, "uv_west");
  gsWriteParaviewPoints(tmp_west, "p_west");

  gsWriteParaviewPoints(uv_east, "uv_east");
  gsWriteParaviewPoints(tmp_east, "p_east");

  gsWriteParaviewPoints(uv_south, "uv_south");
  gsWriteParaviewPoints(tmp_south, "p_south");

  gsWriteParaviewPoints(uv_north, "uv_north");
  gsWriteParaviewPoints(tmp_north, "p_north");


  parameters.resize(uv_interiors.rows(), points.cols());
  parameters << uv_interiors.row(0), uv_south.row(0), uv_east.row(0), uv_north.row(0), uv_west.row(0),
                uv_interiors.row(1), uv_south.row(1), uv_east.row(1), uv_north.row(1), uv_west.row(1);

  points.resize(p_interiors.rows(), parameters.cols());
  points << p_interiors.row(0), p_south.row(0), p_east.row(0), p_north.row(0), p_west.row(0),
                p_interiors.row(1), p_south.row(1), p_east.row(1), p_north.row(1), p_west.row(1),
                p_interiors.row(2), p_south.row(2), p_east.row(2), p_north.row(2), p_west.row(2);

} // end sortPointCloud




int main(int argc, char *argv[])
{
    // Options with default values
    std::string fn = "fitting/deepdrawingC";
    bool bsort = false;

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addSwitch("b", "bsort", "Sort the input points as interiors-boundaries; boundariy points anticlockwise.", bsort);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read data]
    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    std::string ext = (gsFileManager::getExtension(gsFileManager::getFilename(fn))=="") ? ".xml" : "";
    gsFileData<> fd_in(fn+ext);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data

    // Check if matrix sizes are OK
    GISMO_ENSURE( uv.cols() == xyz.cols() && uv.rows() < xyz.rows(),
                  "Wrong input");

    gsWriteParaviewPoints(uv, "uv");
    gsWriteParaviewPoints(xyz, "xyz");

    gsMatrix<> uv_interiors, uv_south, uv_east, uv_north, uv_west, p_interiors, p_south, p_east, p_north, p_west;
    std::vector<index_t> cccc;
    sortPointCloud(uv, xyz, uv_interiors, uv_south, uv_east, uv_north, uv_west, p_interiors, p_south, p_east, p_north, p_west, cccc);

    real_t p_min = xyz.minCoeff(),
           p_max = xyz.maxCoeff();
    real_t den = p_max - p_min;

    gsMatrix<> points(xyz.rows(), xyz.cols());
    points = (1/den)*(xyz - p_min * gsMatrix<>::Ones(xyz.rows(), xyz.cols()));
    gsWriteParaviewPoints(points, "scaled_points");

    gsFileData<> fd;
    fd << uv ;
    fd << points ;

    gsMatrix<> mm(2, 5);
    mm << 1.2, 0.3, 0.1, -0.015, 0.05,
          0, 1, 2, 3, 4;
    gsInfo << mm << "\n";
    
    std::string scalename = fn + "_scalePoints.xml";
    fd.dump(scalename);

    return 0;
  }
