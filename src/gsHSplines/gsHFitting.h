/** @file gsHFitting.h

    @brief Adaptive fitting using hierarchical splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Giannelli, G. Kiss, S. Imperatore, A. Mantzaflaris, D. Mokriš
*/

#pragma once

#include <gsModeling/gsFitting.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>
#include <gsHSplines/gsHBoxUtils.h>

namespace gismo {

/**
    \brief
    This class applies hierarchical fitting of parametrized point clouds.

    \tparam T coefficient type

    \ingroup HSplines
*/

template <short_t d, class T>
class gsHFitting : public gsFitting<T>
{
public:

    typedef typename gsBSplineTraits<d,T>::Basis tensorBasis;
    typedef typename gsFitting<T>::tdm_method tdm_method;

public:
    /// Default constructor
    gsHFitting();

    /**
    *\brief gsHFitting: Main constructor of the h-fitting class
    * \param param_values a matrix containing the parameter values that parametrize the \a points
    * \param points matrix containing the points to be fitted
    * \param basis  hiearchical basis to use for fitting
    * \param refin percentage of errors to refine (if this strategy is chosen)
    * \param extension extension to apply to marked cells
    * \param lambda smoothing weight
    */
    gsHFitting(gsMatrix<T> const & param_values,
               gsMatrix<T> const & points,
               gsHTensorBasis<d,T> & basis,
               T refin, const std::vector<unsigned> & extension,
               T lambda = 0)
    : gsFitting<T>(param_values, points, basis)
    {
        GISMO_ASSERT((refin >=0) && (refin <=1),
                     "Refinement percentage must be between 0 and 1." );
        GISMO_ASSERT(extension.size() == d, "Extension is not of the right dimension");
        GISMO_ASSERT( (gsAsConstVector<unsigned>(extension).array()>=0).all(),
                      "Extension must be a positive number.");

        m_ref    = refin;     //how many % to refine

        m_ext    = extension;

        m_lambda = lambda;    // Smoothing parameter

        m_max_error = m_min_error = 0;

        m_pointErrors.reserve(m_param_values.cols());
    }

public:

    /**
     * @brief iterativeRefine: iteratively refine the basis
     * @param iterations maximum number of iterations
     * @param tolerance (>=0) if the max error is below the tolerance the refinement stops
     * @param err_threshold if non negative all cells with errors
     * bigger than the threshold are refined /
     * If it is equal to -1 the m_ref percentage is used
     * 0 = global refinement
     */
    void iterativeRefine(int iterations, T tolerance, T err_threshold = -1);

    /**
     * @brief nextIteration: perform one iterazion of adaptive refinement with PDM fitting without boundary constraints for \ref iterativeRefine;
     * @param tolerance (>=0) if the maximum error is below the tolerance the refinement stops;
     * @param err_threshold the same as in \ref iterativeRefine.
     */
    bool nextIteration(T tolerance, T err_threshold, index_t maxPcIter = 0);

    /**
     * @brief Like \a nextIteration without \a fixedSides but keeping the values
     * on these sides unchanged throughout the fit.
     */
    bool nextIteration(T tolerance, T err_threshold,
                       const std::vector<boxSide>& fixedSides,
                       index_t maxPcIter = 0);


    /** @brief nextRefinement: One step of the refinement of \ref iterativeRefine;
    * @param tolerance (>=0) if the maximum error is below the tolerance the refinement stops;
    * @param err_threshold same as in \ref iterativeRefine;
    * @param maxPcIter number of parameter correction steps;
    * @param admissibleRef if true, the marking for refinement is admissible.
    */
    bool nextRefinement(T tolerance, T err_threshold, index_t maxPcIter = 0, bool admissibleRef = false);

    /**
     * @brief Like \a nextRefinement without \a fixedSides but keeping the values
     * on these sides unchanged throughout the fit.
     */
    bool nextRefinement(T tolerance, T err_threshold,
                       const std::vector<boxSide>& fixedSides,
                       index_t maxPcIter = 0,
                       bool admissibleRef = false);

    /**
     * @brief nextIteration_tdm: perform one iterazion of adaptive refinement with HDM fitting and with boundary constraints;
     * @param tolerance (>=0) if the maximum error is below the tolerance the refinement stops;
     * @param err_threshold the same as in \ref iterative_refine.
     * @param interpIdx is the index of the boundary points to compute with PDM;
     * @param admissibleRef if true, the refinement is admissible.
     */
    bool nextIteration_tdm(T tolerance, T err_threshold, index_t maxPcIter, T mu, T sigma, const std::vector<index_t> & interpIdx, tdm_method method, bool admissibleRef = false);

    /**
     * @brief Like \a nextIteration_tdm without \a fixedSides but keeping the values
     * on these sides unchanged throughout the fit.
     */
    bool nextIteration_tdm(T tolerance, T err_threshold,
                          const std::vector<boxSide>& fixedSides,
                          index_t maxPcIter,
                          T mu, T sigma,
                          const std::vector<index_t>& interpIdx,
                          tdm_method method,
                          bool admissibleRef);

    /**
     * @brief nextIteration_pdm: perform one iterazion of adaptive refinement with PDM fitting and with boundary constraints;
     * @param tolerance (>=0) if the maximum error is below the tolerance the refinement stops;
     * @param err_threshold the same as in \ref iterativeRefine.
     * @param interpIdx is the index of the boundary points to compute with PDM;
     * @param admissibleRef if true, the refinement is admissible.
     */
    bool nextIteration_pdm(T tolerance, T err_threshold, index_t maxPcIter, const std::vector<index_t> & interpIdx, bool admissibleRef = false);

    /**
     * @brief Like \a nextIteration_pdm without \a fixedSides but keeping the values
     * on these sides unchanged throughout the fit.
     */
    bool nextIteration_pdm(T tolerance, T err_threshold,
                          const std::vector<boxSide>& fixedSides,
                          index_t maxPcIter,
                          const std::vector<index_t>& interpIdx,
                          bool admissibleRef);

    /// Return the refinement percentage
    T getRefPercentage() const
    {
        return m_ref;
    }

    /// Returns the chosen cell extension
    const std::vector<unsigned> & get_extension() const
    {
        return m_ext;
    }

    /// Sets the refinement percentage
    void setRefPercentage(double refPercent)
    {
        GISMO_ASSERT((refPercent >=0) && (refPercent <=1), "Invalid percentage" );
        m_ref = refPercent;
    }

    /// Sets the cell extension
    void setExtension(std::vector<unsigned> const & extension)
    {
        GISMO_ASSERT( (gsAsConstVector<unsigned>(extension).array()>=0).all(),
                      "Extension must be a positive number.");
        GISMO_ASSERT(extension.size()== static_cast<size_t>(this->m_basis.dim()),
                     "Error in dimension");
        m_ext = extension;
    }

    /// Returns boxes which define refinment area.
    std::vector<index_t> getBoxes(const std::vector<T>& errors,
                                   const T threshold);

protected:
    /// Appends a box around \a parameter to the \a boxes only if the box is not already in boxes
    virtual void appendBox(std::vector<index_t>& boxes,
                   std::vector<index_t>& cells,
                   const gsVector<T>& parameter);

    /// Automatic set the refinement threshold
    T setRefineThreshold(const std::vector<T>& errors);

    /// Checks if \a a_cell is already inserted in container \a cells
    static bool isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                                      const std::vector<index_t>& cells);

    /// Appends a \a box to the end of \a boxes (This function also works for cells)
    static void append(std::vector<index_t>& boxes,
                       const gsVector<index_t>& box)
    {
        for (index_t col = 0; col != box.rows(); col++)
            boxes.push_back(box[col]);
    }

    /**
    * @brief getMarkedHBoxesFromBasis_max: returns the markd cells to refine admissibiliy.
    * @param basis: the hierarchical basis from which we extract the elements of the domain
    * @param error: the pointwise parameter error
    * @param parameters: the sites on which the point-wise error is computed
    * @param threshold: the threshold to mark for refinement.
    */
    gsHBoxContainer<d> getMarkedHBoxesFromBasis_max(const gsHTensorBasis<d,T>& basis,
                                                    const std::vector<T>& errors,
                                                    const gsMatrix<T>& parameters,
                                                    T threshold,
                                                    T extension);

protected:

    /// How many % to refine - 0-1 interval
    T m_ref;

    /// Smoothing parameter
    T m_lambda;

    /// Size of the extension
    std::vector<unsigned> m_ext;

    using gsFitting<T>::m_param_values;
    using gsFitting<T>::m_points;
    using gsFitting<T>::m_basis;
    using gsFitting<T>::m_result;

    using gsFitting<T>::m_pointErrors;
    using gsFitting<T>::m_max_error;
    using gsFitting<T>::m_min_error;
};

//perform one iterazion of adaptive refinement for PDM fitting, without boundary constraints
template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration(T tolerance, T err_threshold,
                                     index_t maxPcIter)
{
    std::vector<boxSide> dummy;
    return nextIteration(tolerance, err_threshold, dummy, maxPcIter);
}


//perform one iterazion of adaptive refinement for HDM fitting with boundary constraints
template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration_tdm(T tolerance, T err_threshold,
                                         index_t maxPcIter, T mu, T sigma,
                                         const std::vector<index_t>& interpIdx,
                                         tdm_method method,
                                         bool admissibleRef)
{
    std::vector<boxSide> dummy;
    return nextIteration_tdm(tolerance, err_threshold, dummy, maxPcIter, mu, sigma, interpIdx, method, admissibleRef);
}

//perform one iterazion of adaptive refinement for PDM fitting with boundary constratints
template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration_pdm(T tolerance, T err_threshold,
                                         index_t maxPcIter,
                                         const std::vector<index_t>& interpIdx,
                                         bool admissibleRef)
{
    std::vector<boxSide> dummy;
    return nextIteration_pdm(tolerance, err_threshold, dummy, maxPcIter, interpIdx, admissibleRef);
}

//One step of the refinement
template<short_t d, class T>
bool gsHFitting<d, T>::nextRefinement(T tolerance, T err_threshold, index_t maxPcIter, bool admissibleRef)
{
    std::vector<boxSide> dummy;
    return nextRefinement(tolerance, err_threshold, dummy, maxPcIter, admissibleRef);
}

//perform one iterazion of adaptive refinement for PDM fitting, without boundary constraints
template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration(T tolerance, T err_threshold,
                                     const std::vector<boxSide>& fixedSides,
                                     index_t maxPcIter)
{
    // INVARIANT
    // look at iterativeRefine
    if ( m_pointErrors.size() != 0 )
    {

        if ( m_max_error > tolerance )
        {
            // if err_treshold is -1 we refine the m_ref percent of the whole domain
            T threshold = (err_threshold >= 0) ? err_threshold : setRefineThreshold(m_pointErrors);
            std::vector<index_t> boxes = getBoxes(m_pointErrors, threshold);
            if(boxes.size()==0)
                return false;

            gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d,T> *> (this->m_basis);
            basis->refineElements(boxes);

            // If there are any fixed sides, prescribe the coefs in the finer basis.
            if(m_result != NULL && fixedSides.size() > 0)
            {
                m_result->refineElements(boxes);
                gsFitting<T>::setConstraints(fixedSides);
            }
        }
        else
        {
            return false;
        }
    }

    // SOLVE one fitting step and compute the errors
    this->compute(m_lambda);

    //parameter correction without boudary constraints
    this->parameterCorrection(1e-7, maxPcIter, 1e-4);//closestPoint accuracy, orthogonality tolerance

    // ESTIMATES the point-wise approximation error
    this->computeErrors();

    return true;
}

//perform one iterazion of adaptive refinement for PDM fitting with boundary constraints
template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration_pdm(T tolerance, T err_threshold,
                                        const std::vector<boxSide>& fixedSides,
                                        index_t maxPcIter,
                                        const std::vector<index_t>& interpIdx,
                                        bool admissibleRef)
{
    // look at iterativeRefine
    if ( m_pointErrors.size() != 0 )
    {
        if ( m_max_error > tolerance )
        {
            gsHBoxContainer<2> markedRef;
            std::vector<index_t> boxes;

            // if err_treshold is -1 we refine the m_ref percent of the whole domain
            T threshold = (err_threshold >= 0) ? err_threshold : setRefineThreshold(m_pointErrors);

            gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d,T> *> (this->m_basis);

            if (admissibleRef)
            {
                markedRef = getMarkedHBoxesFromBasis_max(*basis, m_pointErrors, m_param_values, threshold, 2.);
                boxes = markedRef.toRefBoxes();
            }
            else
            {
                boxes = getBoxes(m_pointErrors, threshold);
            }

            if(boxes.size()==0)
                return false;

            basis->refineElements(boxes);
            m_result->refineElements(boxes);

            // If there are any fixed sides, prescribe the coefs in the finer basis.
            if(m_result != NULL && fixedSides.size() > 0)
            {
                m_result->refineElements(boxes);
                gsFitting<T>::setConstraints(fixedSides);
            }
        }
        else
        {
            return false;
        }
    }

    // SOLVE one PDM fitting step
    this->compute(m_lambda);

    //spply maxPcIter parameter correction steps separating interior and boundary points
    this->parameterCorrectionSepBoundary_pdm(1e-6, maxPcIter,interpIdx);//closestPoint accuracy, orthogonality tolerance

    // ESTIMATE the point-wise approximation error
    this->computeErrors();

    return true;
}



//perform one iterazion of adaptive refinement for HDM fitting with boundary constraints
template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration_tdm(T tolerance, T err_threshold,
                                        const std::vector<boxSide>& fixedSides,
                                        index_t maxPcIter, T mu, T sigma,
                                        const std::vector<index_t>& interpIdx,
                                        tdm_method method,
                                        bool admissibleRef)
{
    if ( m_pointErrors.size() != 0 )
    {

    gsHBoxContainer<2> markedRef;
    std::vector<index_t> boxes;

    if ( m_max_error > tolerance )
    {
      // if err_treshold is -1 we refine the m_ref percent of the whole domain
      T threshold = (err_threshold >= 0) ? err_threshold : setRefineThreshold(m_pointErrors);

      // MARK
      gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d,T> *> (this->m_basis);

      // MARK-ADMISSIBLE
      if (admissibleRef)
      {
        markedRef = getMarkedHBoxesFromBasis_max(*basis, m_pointErrors, m_param_values, threshold, 2.);
        boxes = markedRef.toRefBoxes();
      }
      else
      {
        boxes = getBoxes(m_pointErrors, threshold);
      }

      if(boxes.size()==0)
          return false;

      basis->refineElements(boxes);
      m_result->refineElements(boxes);

      // If there are any fixed sides, prescribe the coefs in the finer basis.
      if(m_result != NULL && fixedSides.size() > 0)
        {
          m_result->refineElements(boxes);
          gsFitting<T>::setConstraints(fixedSides);
        }
    }
    else
    {
      return false;
    }
    }

    // SOLVE one HDM fitting step
    this->compute_tdm(m_lambda, mu, sigma, interpIdx, method);


    // apply maxPcIter parameter correction steps separating interior and boundary points
    gsInfo << "Parameter correction: parameterCorrectionSepBoundary_tdm\n";
    this->parameterCorrectionSepBoundary_tdm(1e-6, maxPcIter, mu, sigma, interpIdx, method); // refit

    // ESTIMATE the point-wise approximation error
    this->computeErrors();

    return true;
}

//perform one iterazion of adaptive refinement for PDM fitting without boundary constraints
template<short_t d, class T>
bool gsHFitting<d, T>::nextRefinement(T tolerance, T err_threshold,
                                      const std::vector<boxSide>& fixedSides,
                                      index_t maxPcIter,
                                      bool admissibleRef)
{
    // INVARIANT
    // look at iterativeRefine
    if ( m_pointErrors.size() != 0 )
    {

        if ( m_max_error > tolerance )
        {
            std::vector<index_t> boxes;
            gsHBoxContainer<2> markedRef;
            // if err_treshold is -1 we refine the m_ref percent of the whole domain
            T threshold = (err_threshold >= 0) ? err_threshold : setRefineThreshold(m_pointErrors);

            gsHTensorBasis<d, T>* basis = static_cast<gsHTensorBasis<d,T> *> (this->m_basis);

            // MARK
            if (admissibleRef)
            {
              markedRef = getMarkedHBoxesFromBasis_max(*basis, m_pointErrors, m_param_values, threshold, 2.);
              boxes = markedRef.toRefBoxes();
            }
            else
            {
              boxes = getBoxes(m_pointErrors, threshold);
            }

            if(boxes.size()==0)
                return false;

            basis->refineElements(boxes);
            m_result->refineElements(boxes);

            if(m_result != NULL && fixedSides.size() > 0)
            {
                m_result->refineElements(boxes);
                gsFitting<T>::setConstraints(fixedSides);
            }
        }
        else
        {
            return false;
        }
        // const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis);
        // m_result = bb->makeGeometry( give(this->m_result->coefs()) ).release();

    }
    else
    {
      this->compute(m_lambda);
    }

    this->parameterCorrection(1e-7, maxPcIter, 1e-4);
    this->computeErrors();

    return true;
}


// iteratively refine the basis
template<short_t d, class T>
void gsHFitting<d, T>::iterativeRefine(int numIterations, T tolerance, T err_threshold)
{
    // INVARIANT:
    // m_pointErrors contains the point-wise errors of the fitting
    // therefore: if the size of m_pointErrors is 0, there was no fitting up to this point

    if ( m_pointErrors.size() == 0 )
    {
        this->compute(m_lambda);
        this->computeErrors();
    }

    bool newIteration;
    for( int i = 0; i < numIterations; i++ )
    {
        newIteration = nextIteration( tolerance, err_threshold );
        if( m_max_error <= tolerance )
        {
            break;
        }
        if( !newIteration )
        {
            break;
        }
    }
}

// Returns boxes which define refinment area
template <short_t d, class T>
std::vector<index_t> gsHFitting<d, T>::getBoxes(const std::vector<T>& errors,
                                                 const T threshold)
{
    // cells contains lower corners of elements marked for refinment from maxLevel
    std::vector<index_t> cells;

    // boxes contains elements marked for refinement from differnet levels,
    // format: { level lower-corners  upper-corners ... }
    std::vector<index_t> boxes;

    for (size_t index = 0; index != errors.size(); index++)
    {
        if (threshold <= errors[index])
        {
            appendBox(boxes, cells, this->m_param_values.col(index));
        }
    }

    return boxes;
}

// Appends a box around parameter to the boxes only if the box is not already in boxes
template <short_t d, class T>
void gsHFitting<d, T>::appendBox(std::vector<index_t>& boxes,
                                  std::vector<index_t>& cells,
                                  const gsVector<T>& parameter)
{
    gsTHBSplineBasis<d, T>* basis = static_cast< gsTHBSplineBasis<d,T>* > (this->m_basis);
    const int maxLvl = basis->maxLevel();
    const tensorBasis & tBasis = *(basis->getBases()[maxLvl]);

    // get a cell
    gsVector<index_t, d> a_cell;

    for (short_t dim = 0; dim != d; dim++)
    {
        const gsKnotVector<T> & kv = tBasis.component(dim).knots();
        a_cell(dim) = kv.uFind(parameter(dim)).uIndex();
    }

    if (!isCellAlreadyInserted(a_cell, cells))
    {
        append(cells, a_cell);

        // get level of a cell
        gsVector<index_t, d> a_cell_upp = a_cell + gsVector<index_t, d>::Ones();
        const int cell_lvl = basis->tree().query3(a_cell, a_cell_upp, maxLvl) + 1;

        // get the box
        gsVector<index_t> box(2 * d + 1);
        box[0] = cell_lvl;
        for (short_t dim = 0; dim != d; dim++)
        {
            const unsigned numBreaks = basis->numBreaks(cell_lvl, dim) - 1 ;

            unsigned lowIndex = 0;
            if (cell_lvl < maxLvl)
            {
                const unsigned shift = maxLvl - cell_lvl;
                lowIndex = (a_cell(dim) >> shift);
            }
            else
            {
                const unsigned shift = cell_lvl - maxLvl;
                lowIndex = (a_cell(dim) << shift);
            }

            // apply extensions
            index_t low = ( (lowIndex > m_ext[dim]) ? (lowIndex - m_ext[dim]) : 0 );
            index_t upp = ( (lowIndex + m_ext[dim] + 1 < numBreaks) ?
                             (lowIndex + m_ext[dim] + 1) : numBreaks );

            box[1 + dim    ] = low;
            box[1 + d + dim] = upp;
        }

        append(boxes, box);
    }
}

// Check if a cell is already inserted in (the refinement) container of cells
template <short_t d, class T>
bool gsHFitting<d, T>::isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                                             const std::vector<index_t>& cells)
{

    for (size_t i = 0; i != cells.size(); i += a_cell.rows())
    {
        index_t commonEntries = 0;
        for (index_t col = 0; col != a_cell.rows(); col++)
        {
            if (cells[i + col] == a_cell[col])
            {
                commonEntries++;
            }
        }

        if (commonEntries == a_cell.rows())
        {
            return true;
        }
    }

    return false;
}

// Automatic set of the refinement threshold
template<short_t d, class T>
T gsHFitting<d, T>::setRefineThreshold(const std::vector<T>& errors )
{
    std::vector<T> errorsCopy = errors;
    const size_t i = cast<T,size_t>(errorsCopy.size() * (1.0 - m_ref));
    typename std::vector<T>::iterator pos = errorsCopy.begin() + i;
    std::nth_element(errorsCopy.begin(), pos, errorsCopy.end());
    return *pos;
}


// Check if a point is inside a cell
template <class T>
bool is_point_inside_cell(const gsMatrix<T>& parameter,
                          const gsMatrix<T>& element)
{
    const real_t x = parameter(0, 0);
    const real_t y = parameter(1, 0);

    return element(0, 0) <= x && x < element(0, 1) &&
           element(1, 0) <= y && y < element(1, 1);
}

// Check if a point is inside a cell
template <class T>
bool is_point_inside_cell(const T x,
                          const T y,
                          const gsMatrix<T>& element)
{
    bool condition = (element(0, 0) <= x && x <= element(0, 1) && element(1, 0) <= y && y <= element(1, 1));
    return condition;
}


// Returns the maximum error at the parameters inside the a cell
template<class T>
T getCellMaxError(const gsMatrix<T>& a_cell,
                  const std::vector<T>& errors,
                  const gsMatrix<T>& parameters){

    std::vector<T> a_cellErrs;
    T cell_max_err = 0;
    for(index_t it=0; it < parameters.cols(); it++){
        const T xx = parameters.col(it)(0);
        const T yy = parameters.col(it)(1);
            if (is_point_inside_cell(xx, yy, a_cell))
            {
                a_cellErrs.push_back(errors[it]);
            }
        }

    for(typename std::vector<T>::iterator errIt = a_cellErrs.begin(); errIt != a_cellErrs.end(); ++errIt){
      if (*errIt > cell_max_err){
        cell_max_err = *errIt;
      }
    }
    return cell_max_err;
}

// returns the markd cells to refine admissibiliy.
template <short_t d, class T>
gsHBoxContainer<d> gsHFitting<d, T>::getMarkedHBoxesFromBasis_max(const gsHTensorBasis<d,T>& basis,
                                                const std::vector<T>& errors,
                                                const gsMatrix<T>& parameters,
                                                T threshold,
                                                T extension)
{
    gsHBoxContainer<d> markedHBoxes;
    typename gsBasis<T>::domainIter domItEnd =  basis.domain()->endAll();
    for (auto domIt = basis.domain()->beginAll(); domIt<domItEnd; ++domIt )    // loop over all elements
    {
        gsMatrix<T> elMatrix(d,d);
        elMatrix.col(0)<< domIt.lowerCorner(); // first column  = lower corner
        elMatrix.col(1)<< domIt.upperCorner(); // second column = upper corner
        T cellMaxError = getCellMaxError(elMatrix, errors, parameters);
        if (cellMaxError >= threshold)
        {
            gsHDomainIterator<T,d> * domHIt = nullptr;
            domHIt = dynamic_cast<gsHDomainIterator<T,2> *>(domIt.get());
            gsHBox<d> a_box(domHIt);
            gsHBoxContainer<d> tmp(gsHBoxUtils<d,T>::markAdmissible(a_box,extension));
            markedHBoxes.add(tmp);
        }
    }
    return markedHBoxes;
}

}// namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsFitting.hpp)
// #endif
