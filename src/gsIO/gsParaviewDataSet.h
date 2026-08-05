/** @file gsParaviewDataSet.h

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsOptionList.h>
#include <gsIO/gsParaviewUtils.h>
#include <gsExpressions/gsExprHelper.h>

#include <fstream>

namespace gismo
{

template <class E, class T>
std::vector<std::string> toParaview(const expr::_expr<E>& expr,
                                    gsExprEvaluator<T>& evaltr,
                                    unsigned nPts, unsigned precision,
                                    std::string label,
                                    const bool& export_base64);

/**
    \brief This class represents a group of vtk (Paraview) files
    that refer to one multiPatch, for one timestep.

    This class is used by gsParaviewCollection to manage said files,
    but can be used by the user explicitly as well.

    \ingroup IO
*/
template <class T>
class GISMO_EXPORT gsParaviewDataSet
{
private:
    std::string m_basename;
    std::vector<std::string> m_filenames;
    const gsMultiPatch<T>* m_geometry;
    gsExprEvaluator<T>* m_evaltr;
    gsOptionList m_options;
    bool m_isSaved;

public:
    typedef memory::unique_ptr<gsParaviewDataSet<T> > uPtr;

    /// @brief Constructor without evaluator (field-only workflows)
    gsParaviewDataSet(std::string basename,
                      const gsMultiPatch<T>& geometry,
                      gsOptionList options = defaultOptions());

    /// @brief Constructor with evaluator (expression workflows)
    gsParaviewDataSet(std::string basename,
                      const gsMultiPatch<T>& geometry,
                      gsExprEvaluator<T>& eval,
                      gsOptionList options = defaultOptions());

    gsParaviewDataSet() = delete;

    /// @brief Evaluates an expression, and writes that data to the vtk files.
    template <class E>
    void addField(const expr::_expr<E>& expr, std::string label)
    {
        GISMO_ENSURE(m_evaltr, "expression fields need an evaluator");
        GISMO_ENSURE(!m_isSaved,
                     "You cannot add more fields if the gsParaviewDataSet has "
                     "been saved.");

        const unsigned nPts = getNumPoints(m_options);
        const unsigned precision = static_cast<unsigned>(m_options.askInt("precision", 5));
        const bool export_base64 = m_options.askSwitch("base64", false);

        const std::vector<std::string> tags =
            toParaview(expr, *m_evaltr, nPts, precision, label, export_base64);
        const std::vector<std::string> fnames = filenames();
        const gsMultiPatch<T>& geometry = *m_geometry;

        for (index_t k = 0; k != geometry.nPieces(); ++k)
        {
            std::ofstream file;
            file.open(fnames[k].c_str(), std::ios_base::app);
            file << tags[k];
            file.close();
        }
    }

    // Just here to stop the recursion
    void addFields(std::vector<std::string>) { }

    /// @brief Recursive form of addField() for expressions
    template <class E, typename... Rest>
    void addFields(std::vector<std::string> labels, const expr::_expr<E>& expr, Rest... rest)
    {
        GISMO_ENSURE(sizeof...(Rest) == labels.size() - 1,
                     "The length of labels must match the number of expressions provided");
        std::vector<std::string> newlabels(labels.cbegin() + 1, labels.cend());
        addField(expr, labels[0]);
        addFields(newlabels, rest...);
    }

    /// @brief Evaluates a gsField and writes that data to the vtk files.
    template <class U>
    void addField(const gsField<U> field, std::string label)
    {
        GISMO_ENSURE(!m_isSaved,
                     "You cannot add more fields if the gsParaviewDataSet has "
                     "been saved.");
        const gsMultiPatch<T>& geometry = *m_geometry;
        GISMO_ENSURE(
            (field.parDim() == geometry.domainDim() &&
             field.geoDim() == geometry.targetDim() &&
             field.nPieces() == geometry.nPieces() &&
             field.patches().coefsSize() == geometry.coefsSize()),
            "Provided gsField and stored geometry are not compatible!");

        const unsigned nPts = getNumPoints(m_options);
        const unsigned precision = static_cast<unsigned>(m_options.askInt("precision", 5));
        const bool export_base64 = m_options.askSwitch("base64", false);

        const std::vector<std::string> tags =
            toParaview(field, nPts, precision, label, export_base64);
        const std::vector<std::string> fnames = filenames();

        for (index_t k = 0; k != geometry.nPieces(); ++k)
        {
            std::ofstream file;
            file.open(fnames[k].c_str(), std::ios_base::app);
            file << tags[k];
            file.close();
        }
    }

    /// @brief Recursive form of addField() for gsField
    template <class U, typename... Rest>
    void addFields(std::vector<std::string> labels, const gsField<U> field, Rest... rest)
    {
        std::vector<std::string> newlabels(labels.cbegin() + 1, labels.cend());
        addField(field, labels[0]);
        addFields(newlabels, rest...);
    }

    /// @brief Returns the names of the files created by this gsParaviewDataSet.
    const std::vector<std::string> filenames();

    void save();

    bool isEmpty();

    bool isSaved();

    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt("plot.npts", "Number of points per-patch.", 1000);
        opt.addInt("precision", "Number of decimal digits.", 5);
        opt.addInt("plotElements.resolution", "Drawing resolution for element mesh.", -1);
        opt.addSwitch("makeSubfolder", "Export vtk files to subfolder ( below the .pvd file ).", true);
        opt.addSwitch("base64", "Export in base64 binary format", false);
        opt.addString("subfolder", "Name of subfolder where the vtk files will be stored.", "");
        opt.addSwitch("plot.elements", "Controls plotting of element mesh.", false);
        opt.addSwitch("plotControlNet", "Controls plotting of control point grid.", false);
        return opt;
    }

    gsOptionList& options() { return m_options; }

private:
    static unsigned getNumPoints(const gsOptionList& opts);
    static bool getPlotElements(const gsOptionList& opts);
    void initFilenames();
};

} // End namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParaviewDataSet.hpp)
#endif
