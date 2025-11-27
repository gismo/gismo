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
#include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
#include <gsCore/gsDofMapper.h>         // Only to make linker happy
#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsIO/gsParaviewUtils.h>

#include<fstream>

namespace gismo
{
/**
    \brief This class represents a group  of vtk (Paraview) files
    that refer to one multiPatch, for one timestep.

    This class is used by gsParaviewCollection to manage said files,
    but can be used by the user explicitly as well.

    \ingroup IO
*/
class GISMO_EXPORT gsParaviewDataSet // a collection of .vts files
{
private:
    std::string m_basename;
    std::vector<std::string> m_filenames;
    const gsMultiPatch<real_t> * m_geometry;
    gsExprEvaluator<real_t> * m_evaltr;
    gsOptionList m_options;
    bool m_isSaved;

public:
    /// @brief Basic constructor
    /// @param basename The basename that will be used to create all the individual filenames
    /// @param geometry A gsMultiPatch of the geometry that will be exported and where the fields are defined
    /// @param eval Optional. A gsExprEvaluator, necessary when working with gsExpressions for evaluation purposes
    /// @param options A set of options, if unspecified, defaultOptions() is called.
    gsParaviewDataSet(std::string basename,
                      gsMultiPatch<real_t> * const geometry,
                      gsExprEvaluator<real_t> * eval=nullptr,
                      gsOptionList options=defaultOptions());

    gsParaviewDataSet():m_basename(""),
                        m_geometry(nullptr),
                        m_evaltr(nullptr),
                        m_options(defaultOptions()),
                        m_isSaved(false)
                        {}

    /// @brief Evaluates an expression, and writes that data to the vtk files.
    /// @tparam E
    /// @param expr The gsExpression to be evaluated
    /// @param label The name that will be displayed in Paraview for this field.
    template <class E>
    void addField(const expr::_expr<E>& expr, std::string label) {
        GISMO_ENSURE(!m_isSaved,
                     "You cannot add more fields if the gsParaviewDataSet has "
                     "been saved.");
        // evaluates the expression and appends it to the vts files
        // for every patch
        const unsigned nPts = m_options.askInt("numPoints", 1000);
        const unsigned precision = m_options.askInt("precision", 5);
        const bool export_base64 = m_options.askSwitch("base64", false);

        // gsExprEvaluator<real_t> ev;
        // gsMultiBasis<real_t> mb(*m_geometry);
        // ev.setIntegrationElements(mb);
        const std::vector<std::string> tags =
            toVTK(expr, m_evaltr, nPts, precision, label, export_base64);
        std::vector<std::string> fnames = filenames();

        for (index_t k = 0; k != m_geometry->nPieces();
             k++)  // For every patch.
        {
            std::ofstream file;
            file.open(fnames[k].c_str(), std::ios_base::app);  // Append to file
            file << tags[k];
            file.close();
        }
    }

    // Just here to stop the recursion
    void addFields(std::vector<std::string>){ }


    /// @brief Recursive form of addField()
    /// @tparam E
    /// @tparam ...Rest
    /// @param labels Vector of strings, containing the names of the fields as they will be shown in ParaView.
    /// @param expr The expressions to be evaluated ( arbitrary number of them )
    /// @param ...rest
    template <class E, typename... Rest>
    void addFields(std::vector<std::string> labels, const expr::_expr<E> & expr, Rest... rest) {
        // keep all but first label
        GISMO_ENSURE( sizeof...(Rest) == labels.size() - 1, "The length of labels must match the number of expressions provided" );
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());


        addField(   expr, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
        addFields(   newlabels, rest...);   // Recursion
    }

    /// @brief Evaluates a gsField ( the function part ), and writes that data
    /// to the vtk files.
    /// @tparam T
    /// @param field The gsField to be evaluated
    /// @param label The name that will be displayed in Paraview for this field.
    template <class T>
    void addField(const gsField<T> field, std::string label) {
        GISMO_ENSURE(!m_isSaved,
                     "You cannot add more fields if the gsParaviewDataSet has "
                     "been saved.");
        GISMO_ENSURE(
            (field.parDim() == m_geometry->domainDim() &&
             field.geoDim() == m_geometry->targetDim() &&
             field.nPieces() == m_geometry->nPieces() &&
             field.patches().coefsSize() == m_geometry->coefsSize()),
            "Provided gsField and stored geometry are not compatible!");
        // evaluates the field  and appends it to the vts files
        // for every patch
        const unsigned nPts = m_options.askInt("numPoints", 1000);
        const unsigned precision = m_options.askInt("precision", 5);
        const bool export_base64 = m_options.askSwitch("base64", false);

        const std::vector<std::string> tags =
            toVTK(field, nPts, precision, label, export_base64);
        const std::vector<std::string> fnames = filenames();

        for (index_t k = 0; k != m_geometry->nPieces();
             k++)  // For every patch.
        {
            std::ofstream file;
            file.open(fnames[k].c_str(), std::ios_base::app);  // Append to file
            file << tags[k];
            file.close();
        }
    }

    /// @brief Recursive form of addField()
    /// @tparam T
    /// @tparam ...Rest
    /// @param labels Vector of strings, containing the names of the fields as they will be shown in ParaView.
    /// @param field The gsFields to be evaluated ( arbitrary number of them )
    /// @param ...rest
    template <class T, typename... Rest>
    void addFields(std::vector<std::string> labels, const gsField<T> field, Rest... rest) {
        // keep all but first label
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());

        addField(   field, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
        addFields(   newlabels, rest...);   // Recursion
    }

    /// @brief Returns the names of the files created by this gsParaviewDataSet.
    /// @return A vector of strings
    const std::vector<std::string> filenames();

    void save();

    bool isEmpty();

    bool isSaved();

    /// @brief Accessor to the current options.
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt("numPoints", "Number of points per-patch.", 1000);
        opt.addInt("precision", "Number of decimal digits.", 5);
        opt.addInt("plotElements.resolution", "Drawing resolution for element mesh.", -1);
        opt.addSwitch("makeSubfolder", "Export vtk files to subfolder ( below the .pvd file ).", true);
        opt.addSwitch("base64", "Export in base64 binary format", false);
        opt.addString("subfolder","Name of subfolder where the vtk files will be stored.", "");
        opt.addSwitch("plotElements", "Controls plotting of element mesh.", false);
        opt.addSwitch("plotControlNet", "Controls plotting of control point grid.", false);
        return opt;
    }

    gsOptionList & options() {return m_options;}

private:

 void initFilenames();
};
} // End namespace gismo
