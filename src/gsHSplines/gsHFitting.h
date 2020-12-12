/** @file gsHFitting.h

    @brief Adaptive fitting using hierarchical splines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Giannelli, G. Kiss
*/

#pragma once

#include <gsModeling/gsFitting.h>
#include <gsHSplines/gsHTensorBasis.h>

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

public:
    /// Default constructor
    gsHFitting();

    /**
        \brief
        Main constructor of the fitting class

        \param param_values a matrix containing the parameter values
        that parametrize the \a points

        \param points The points to be fitted

        \param basis  Hiearchical basis to use for fitting

        \param refin Percentage of errors to refine (if this strategy is chosen)

        \param extension Extension to apply to marked cells

        \param lambda Smoothing parameter
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
     * @brief iterative_refine iteratively refine the basis
     *
     *
     * @param iterations maximum number of iterations
     *
     * @param tolerance (>=0) if the max error is below the tolerance the refinement stops
     *
     * @param err_threshold if non negative all cells with errors
     * bigger than the threshold are refined /
     * If it is equal to -1 the m_ref percentage is used
     * 0 = global refinement
     */
    void iterativeRefine(int iterations, T tolerance, T err_threshold = -1);

    /**
     * @brief nextIteration One step of the refinement of iterative_refine(...);
     * @param tolerance (>=0) if the maximum error is below the tolerance the refinement stops;
     * @param err_threshold the same as in iterative_refine(...).
     */
    bool nextIteration(T tolerance, T err_threshold);

    /**
     * @brief Like \a nextIteration without \a fixedSides but keeping the values
     * on these sides unchanged throughout the fit.
     */
    bool nextIteration(T tolerance, T err_threshold,
		       const std::vector<boxSide>& fixedSides);

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

    /// Sets constraints in such a way that the previous values at \a
    /// fixedSides of the geometry remain intact.
    void setConstraints(const std::vector<boxSide>& fixedSides);

    /// Set constraints in such a way that the resulting geometry on
    /// each of \a fixedSides will coincide with the corresponding
    /// curve in \a fixedCurves.
    void setConstraints(const std::vector<boxSide>& fixedSides,
			const std::vector<gsBSpline<T> >& fixedCurves);

protected:
    /// Appends a box around parameter to the boxes only if the box is not
    /// already in boxes
    virtual void appendBox(std::vector<index_t>& boxes,
                   std::vector<index_t>& cells,
                   const gsVector<T>& parameter);

    /// Identifies the threshold from where we should refine
    T setRefineThreshold(const std::vector<T>& errors);

    /// Checks if a_cell is already inserted in container of cells
    static bool isCellAlreadyInserted(const gsVector<index_t, d>& a_cell,
                                      const std::vector<index_t>& cells);

    /// Appends a box to the end of boxes (This function also works for cells)
    static void append(std::vector<index_t>& boxes,
                       const gsVector<index_t>& box)
    {
        for (index_t col = 0; col != box.rows(); col++)
            boxes.push_back(box[col]);
    }

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

template<short_t d, class T>
void gsHFitting<d, T>::setConstraints(const std::vector<boxSide>& fixedSides)
{
    if(fixedSides.size() == 0)
	return;

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;

    for(std::vector<boxSide>::const_iterator it=fixedSides.begin(); it!=fixedSides.end(); ++it)
    {
	gsMatrix<index_t> ind = this->m_basis->boundary(*it);
	for(index_t r=0; r<ind.rows(); r++)
	{
	    index_t fix = ind(r,0);
	    // If it is a new constraint, add it.
	    if(std::find(indices.begin(), indices.end(), fix) == indices.end())
	    {
		indices.push_back(fix);
		coefs.push_back(this->m_result->coef(fix));
	    }
	}
    }

    gsFitting<T>::setConstraints(indices, coefs);
}

template<short_t d, class T>
void gsHFitting<d, T>::setConstraints(const std::vector<boxSide>& fixedSides,
				      const std::vector<gsBSpline<T> >& fixedCurves)
{
    if(fixedSides.size() == 0)
	return;

    GISMO_ASSERT(fixedCurves.size() == fixedSides.size(),
		 "fixedCurves and fixedSides are of different sizes.");

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;
    for(size_t s=0; s<fixedSides.size(); s++)
    {
	gsMatrix<T> coefsThisSide = fixedCurves[s].coefs();
	gsMatrix<index_t> indicesThisSide = m_basis->boundaryOffset(fixedSides[s],0);
	GISMO_ASSERT(coefsThisSide.rows() == indicesThisSide.rows(),
		     "Coef number mismatch between prescribed curve and basis side.");

	for(index_t r=0; r<indicesThisSide.rows(); r++)
	{
	    index_t fix = indicesThisSide(r,0);
	    // If it is a new constraint, add it.
	    if(std::find(indices.begin(), indices.end(), fix) == indices.end())
	    {
		indices.push_back(fix);
		coefs.push_back(coefsThisSide.row(r));
	    }
	}
    }

    gsFitting<T>::setConstraints(indices, coefs);
}

template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration(T tolerance, T err_threshold)
{
    std::vector<boxSide> dummy;
    return nextIteration(tolerance, err_threshold, dummy);
}

template<short_t d, class T>
bool gsHFitting<d, T>::nextIteration(T tolerance, T err_threshold,
				     const std::vector<boxSide>& fixedSides)
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
		setConstraints(fixedSides);
	    }

            gsDebug << "inserted " << boxes.size() / (2 * d + 1) << " boxes.\n";
        }
        else
        {
            gsDebug << "Tolerance reached.\n";
            return false;
        }
    }

    // We run one fitting step and compute the errors
    this->compute(m_lambda);
    this->computeErrors();

    return true;
}

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
            gsDebug << "Tolerance reached at iteration: " << i << "\n";
            break;
        }
        if( !newIteration )
        {
            gsDebug << "No more Boxes to insert at iteration: " << i << "\n";
            break;
        }
    }
}

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

template<short_t d, class T>
T gsHFitting<d, T>::setRefineThreshold(const std::vector<T>& errors )
{
    std::vector<T> errorsCopy = errors;
    const size_t i = cast<T,size_t>(errorsCopy.size() * (1.0 - m_ref));
    typename std::vector<T>::iterator pos = errorsCopy.begin() + i;
    std::nth_element(errorsCopy.begin(), pos, errorsCopy.end());
    return *pos;
}


}// namespace gismo
