/** @file gsGoTools.h

    @brief Input and output Utilities concerning GoTools export

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Speh
*/

#pragma once

#include <fstream>

#include <gsUtils/gsMesh/gsMesh.h>
#include <gsUtils/gsCombinatorics.h>

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHDomainIterator.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsFileManager.h>

namespace gismo {


/// @brief
/// Writes body part of GoTools (.g2) format to a file.
///
/// \param bspl tensor B-Spline
/// \param out file stream
///
/// \ingroup IO
template <short_t d, typename T>
void gsWriteGoToolsBodySpline(const gsTensorBSpline<d, T>& bspl,
                              std::ofstream& out)
{
    // The format is (similar for d < 3)
    //
    // - dimension of geometry space, whether or not the volume is
    //   rational: 1=rational, 0=non-rational
    // - the number of coefficients in the first parameter direction,
    //   the polynomial order in this direction (i.e. degree+1)
    // - the knot vector in the first parameter direction, multiple
    //   knots are represented by giving the knot value several times
    // - the number of coefficients in the second parameter direction,
    //   the polynomial order in this direction (i.e. degree+1)
    // - the knot vector in the second parameter direction
    // - the number of coefficients in the third parameter direction,
    //   the polynomial order in this direction (i.e. degree+1)
    // - the knot vector in the third parameter direction
    // - the volume coefficients



    out << std::setprecision(15);

    out << bspl.geoDim() << " " << 0 << "\n";

    for (short_t dim = 0; dim < bspl.parDim(); dim++)
    {
        out << bspl.basis().size(dim) << " " <<
               bspl.basis().degree(dim) + 1 << "\n";

        const gsKnotVector<T> kv = bspl.basis().knots(dim);
        typename gsKnotVector<T>::const_iterator iter;
        for (iter = kv.begin(); iter != kv.end(); iter++)
        {
            out << *iter << " ";
        }
        out << "\n";
    }

    const gsMatrix<T>& coefs = bspl.coefs();

    for (int row = 0; row < coefs.rows(); row++)
    {
        out << coefs.row(row) << "\n";

    }

    out << std::endl;
}


/// Writes bspline in GoTools (.g2) format to a file.
///
/// \param bspl tensor B-Spline
/// \param out file stream
template <short_t d, typename T>
void gsWriteGoToolsSpline(const gsTensorBSpline<d, T>& bspl,
                          std::ofstream& out)
{

    // first we write header

    // format ID VERSION, where
    // ID : Class_SplineCurve = 100,
    //      Class_SplineSurface = 200,
    //      Class_SplineVolume = 700
    // VERSION 1 0 0

    // about the headers
    // https://github.com/SINTEF-Geometry/GoTools/blob/master/gotools-core/include/GoTools/geometry/ClassType.h
    // each element in the patch must have its own header

    if (d == 1)
    {
        out << 100;
    }
    else if (d == 2)
    {
        out << 200;
    }
    else if (d == 3)
    {
        out << 700;
    }
    else
    {
        gsWarn << "Dimension is too high: GoTools does not support dimension "
                  "higher than two."
                  "Aborting ...";
        return;
    }

    out << " 1 0 0\n";

    gsWriteGoToolsBodySpline<d, T>(bspl, out);

}


/// Writes geometry in GoTools (.g2) format to a file.
///
/// \param geom geometry
/// \param fileName output file name
template <typename T>
void gsWriteGoTools(const gsGeometry<T>& geom,
                    const std::string& fileName)
{
    std::string fn(fileName);

    // check the extension
    std::string ext = gsFileManager::getExtension(fileName);
    if (ext != "g2")
    {
         fn += ".g2";
    }

    // opening file
    std::ofstream file;
    file.open(fn.c_str());
    if (!file.is_open())
    {
        gsWarn << "Can not open file: " << fileName << "\n"
                  "Aborting ...";
        return;
    }

    gsWriteGoTools(geom, file);

    file.close();

}


/// Writes geometry in GoTools (.g2) format to a file.
///
/// \param geom geometry
/// \param out file stream
template <typename T>
void gsWriteGoTools(const gsGeometry<T>& geom,
                    std::ofstream& out)
{

    // check which object we have and call appropriate function
    const gsTensorBSpline<1, T>* bspline1 =
            dynamic_cast<const gsTensorBSpline<1, T>* >(&geom);

    if (bspline1 != NULL)
    {
        gsWriteGoToolsSpline<1, T>(*bspline1, out);
        return;
    }

    const gsTensorBSpline<2, T>* bspline2 =
            dynamic_cast<const gsTensorBSpline<2, T>* >(&geom);

    if (bspline2 != NULL)
    {
        gsWriteGoToolsSpline<2, T>(*bspline2, out);
        return;
    }

    const gsTensorBSpline<3, T>* bspline3 =
            dynamic_cast<const gsTensorBSpline<3, T>* >(&geom);

    if (bspline3 != NULL)
    {
        gsWriteGoToolsSpline<3, T>(*bspline3, out);
        return;
    }

    gsWarn << "Can not write your geometry, unsupported format.\n"
              "Currently supported formats: BSplines.\n"
              "Aborting ...\n";

    return;
}


/// Writes multi patch in GoTools (.g2) format to a file.
///
/// \param multiPatch multi patch
/// \param fileName output file name
template <typename T>
void gsWriteGoTools(const gsMultiPatch<T>& multiPatch,
                    const std::string& fileName)
{
    std::string fn(fileName);

    // check the extension
    std::string ext = gsFileManager::getExtension(fileName);
    if (ext != "g2")
    {
        fn += ".g2";
    }

    // opening file
    std::ofstream file;
    file.open(fn.c_str());
    if (!file.is_open())
    {
        gsWarn << "Can not open file: " << fileName << "\n"
                  "Aborting ...";
        return;
    }

    typename gsMultiPatch<T>::const_iterator iter;
    for (iter = multiPatch.begin(); iter != multiPatch.end(); ++iter)
    {
        gsWriteGoTools<T>(**iter, file);
    }

    file.close();

}


} // end namespace gismo
