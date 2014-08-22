
#pragma once

#include <gsModeling/gsFitting.h>
#include <gsThbs/gsHTensorBasis.h>

namespace gismo {

/** 
    \brief
    This class applies hierarchical fitting of parametrized point clouds.

    \tparam T coefficient type

*/

template <class T>
class gsHierarchical_refinement : public gsFitting<T>
{
public:
    /// Default constructor
    gsHierarchical_refinement();

    /** 
        \brief
        Main constructor of the fitting class
        
        \param param_values a matrix containing the parameter values that parametrize the \a points

        \param points The points to be fitted

        \param basis  Hiearchical basis to use for fitting

        \param refin Wercentage of errors to refine (if this strategy is chosen)

        \param extension Wxtension tp apply to marked cells
        
        \param lambda Smoothing parameter
    */
    gsHierarchical_refinement(gsMatrix<T> const & param_values, 
                              gsMatrix<T> const & points, 
                              gsHTensorBasis<2,T> & basis, 
                              T refin, std::vector<int> extension, 
                              T lambda = 0)
    : gsFitting<T>(param_values, points, basis)
    {
        GISMO_ASSERT((refin >=0) && (refin <=1), "Refinement percentage must be between 0 and 1." );
        for(std::size_t i = 0; i < extension.size();i++)
        {
            GISMO_ASSERT(extension[i]>=0, "Extension must be a positive number.");
        }

        GISMO_ASSERT(extension.size()== static_cast<std::size_t>(basis.dim()), 
                     "Extension is not of the right dimension");

        m_ref    = refin;     //how many % to refine

        m_ext.resize(2);
        m_ext[0] = extension[0]; // size of the extension
        m_ext[1] = extension[1];

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
     * @param lambda  smoothing parameter
     *
     * @param err_type  0 = minimization term, 1 = euclidean distance
     *
     * @param err_threshold if non negative all cells with errors
     * bigger the the threshold are refined / 
     * If it is equal to -1 the m_ref percentage is used 
     * 0 = global refinement
     */
    void iterative_refine(int iterations, T tolerance, 
                          int err_type = 1, T err_threshold = -1);

    /// Return the errors for each point
    const std::vector<T> & pointWiseErrors() const
    {
        return m_pointErrors;
    }


    /// Return the refinement percentage
    T get_refinement() const
    {
        return m_ref;
    }

    /// Returns the chosen cell extension
    const std::vector<unsigned> & get_extension() const 
    {
        return m_ext;
    }

    /// Sets the refinement percentage
    void set_refinement(double refin)
    {
        GISMO_ASSERT((refin >=0) && (refin <=1), "Invalid percentage" );
        m_ref = refin;
    }

    /// Sets the cell extension
    void set_extension(std::vector<unsigned> const & extension)
    {
        for(std::size_t i = 0; i < extension.size();i++)
        {
            GISMO_ASSERT(extension[i]>=0, "Extension must be non negative.");
        }
        assert(extension.size()== unsigned(this->m_basis.dim()));
        m_ext = extension;
    }

    /// Returns the minimum point-wise error from the pount cloud (or zero if not fitted)
    T minPointError() const { return m_min_error; }

    /// Returns the maximum point-wise error from the pount cloud (or zero if not fitted)
    T maxPointError() const { return m_max_error; }

    /// Computes the number of points below the error threshold (or zero if not fitted)
    std::size_t numPointsBelow(T threshold) const 
    { 
        const std::size_t result= 
            std::count_if(m_pointErrors.begin(), m_pointErrors.end(), 
                          std::bind2nd(std::less<T>(), threshold));
        return result; 
    }

private:

    /// Identifies the threshold from where we should refine
    T to_refine(const std::vector<T> & errors);

    /// Identifies the cells in the highest level where the error is bigger then max_error
    std::vector<int> select_cells(std::vector<T> & errors, double max_error);

    /// Return the levels in which the cells are included
    std::vector<int> get_levels(const std::vector<int> & cells);

    /// Combine cells and levels to get the refinement boxes
    std::vector<unsigned> get_boxes(const std::vector<int> & cells, 
                                    const std::vector<int> & levels);

    /// Computes the euclidean error for each point
    void computeErrors();

private:

    /// How many % to refine - 0-1 interval
    T m_ref;

    /// Smoothing parameter
    T m_lambda; 

    /// Size of the extension
    std::vector<unsigned> m_ext;

    // All point-wise errors
    std::vector<T> m_pointErrors;

    /// Maximum point-wise error
    T m_max_error;

    /// Minimum point-wise error
    T m_min_error;

    using gsFitting<T>::m_param_values;
    using gsFitting<T>::m_points;
    using gsFitting<T>::m_basis;
    using gsFitting<T>::m_result;
};



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


template<class T>
void gsHierarchical_refinement<T>::iterative_refine(int iterations, 
                                                    T tolerance, //
                                                    int err_type, 
                                                    T err_threshold)
{
    gsHTensorBasis<2,T>  * h_basis = static_cast<gsHTensorBasis<2,T> *> (this->m_basis);

    // are we computing for the first time ?
    if ( m_pointErrors.size() == 0 )
    {
        this->compute(m_lambda);

        computeErrors();
        
        return;
    }
    
    std::vector<int> cells ; //cells to be refined
    std::vector<int> levels; //cells to be refined
    std::vector<unsigned> boxes;

    for(int i = 0; i < iterations; i++)
    {
        if ( m_max_error > tolerance ) // are we done ?
        {
            if(err_threshold>=0) // refinement strategy
            {
                cells = select_cells(m_pointErrors, err_threshold);
            }
            else
            {
                cells = select_cells(m_pointErrors, to_refine(m_pointErrors));
            }

            levels = get_levels(cells); //get the  level of the cells
            boxes  = get_boxes(cells,levels);

            h_basis->refine(boxes);

            gsInfo <<"inserted "<< levels.size() <<" boxes.\n";

            this->compute(m_lambda);

            computeErrors();
        }
        else
        {
            std::cout<<"Tolerance reached at iteration: "<< i <<"\n";
            break;
        }
    }
}



////private functions
template<class T>
T gsHierarchical_refinement<T>::to_refine(const std::vector<T> & errors)
{
    std::vector<T> tmp = errors; // temp copy
    const std::size_t i = static_cast<std::size_t>(errors.size()*(1.0-m_ref));
    typename std::vector<T>::iterator pos = tmp.begin() + i;
    std::nth_element(tmp.begin(), pos, tmp.end());
    return *pos;
}

template<class T>
std::vector<int> gsHierarchical_refinement<T>::get_levels(const std::vector<int> & cells)
{
    gsHTensorBasis<2,T>  * h_basis = static_cast<gsHTensorBasis<2,T> *> (this->m_basis);

    const int level = h_basis->get_max_inserted_level();
    //const int level = h_basis->maxAllowdLevel();

    std::vector<int> result;
    gsVector<unsigned int> lower(2);
    gsVector<unsigned int> upper(2);
    for( std::size_t i = 0; i < cells.size(); i+=2 )
    {
        lower[0] = cells[i];
        lower[1] = cells[i+1];
        upper[0] = cells[i]+1;
        upper[1] = cells[i+1]+1;
        result.push_back( h_basis->m_tree.query3(lower, upper, level) );
    }
    return result;
}

template<class T>
std::vector<unsigned> gsHierarchical_refinement<T>::get_boxes(const std::vector<int> & cells, 
                                                              const std::vector<int> & levels)
{
    gsHTensorBasis<2,T>  * h_basis = static_cast<gsHTensorBasis<2,T> *> (this->m_basis);

    const int maxLevel = h_basis->get_max_inserted_level(); 
    //const int level = h_basis->maxAllowdLevel();

    std::vector<unsigned> result;

    unsigned a;
    for(std::size_t i = 0; i < cells.size(); i+=2)
    {
        const int newLevel   = levels[i/2]+1;

        const unsigned sz0 = h_basis->numBreaks(newLevel,0) - 1 ;
        // equiv: sz0 = h_basis->m_bases[newLevel]->component(0).knots().uSize() - 1 ;
        const unsigned sz1 = h_basis->numBreaks(newLevel,1) - 1 ;
        // equiv: sz1 = h_basis->m_bases[newLevel]->component(1).knots().uSize() - 1 ;
        const unsigned shift = (maxLevel<newLevel ? newLevel - maxLevel : maxLevel - newLevel);

        result.push_back(newLevel);

        // Lower - x
        a  = ( (cells[i]+1) << shift );
        //equiv: a = Knot0.Uniquefindspan( lKnot0.uValue(cells[i]+1) );
        result.push_back( ( a > m_ext[0] ? a - m_ext[0] : 0 ) );

        // Lower - y
        a = ( (cells[i+1]+1) << shift );
        //equiv: a = Knot1.Uniquefindspan( lKnot1.uValue(cells[i+1]+1) );
        result.push_back( ( a > m_ext[1] ? a - m_ext[1] : 0 ) );

        // Upper - x
        a = ( (cells[i]+1) << shift );
        //equiv: a = Knot0.Uniquefindspan( lKnot0.uValue(cells[i]+1) );
        result.push_back( ( a + m_ext[0] + 1 < sz0 ? 
                            a + m_ext[0] + 1 : sz0 ) );

        // Upper - y
        a = ( (cells[i+1]+1) << shift );        
        //equiv: a  = Knot1.Uniquefindspan(lKnot1.uValue(cells[i+1]+1) );
        result.push_back( ( a + m_ext[1] + 1 < sz1 ? 
                            a + m_ext[1] + 1 : sz1 ) );
    }

    return result;
}

template<class T>
std::vector<int> gsHierarchical_refinement<T>::select_cells(std::vector<T> & errors, 
                                                            double max_error)
{
    gsHTensorBasis<2,T>  * h_basis = static_cast<gsHTensorBasis<2,T> *> (this->m_basis);
    std::vector<int> result;
    std::vector<int> temp;

    const int level = h_basis->get_max_inserted_level();
    //const int level = h_basis->maxAllowdLevel();

    const gsCompactKnotVector<T> & lKnot0  = 
        h_basis->m_bases[level]->component(0).knots();
    const gsCompactKnotVector<T> & lKnot1  = 
        h_basis->m_bases[level]->component(1).knots();

    for(std::size_t i= 0; i < errors.size(); i++)
    {
        if(errors[i]>=max_error)
        {
            //find the cell in max level

            temp.push_back( lKnot0.Uniquefindspan( m_param_values(0,i) ) );
            temp.push_back( lKnot1.Uniquefindspan( m_param_values(1,i) ) );

            // TODO; get the boxes here at once.

            bool insert= true;
            // Compare temp with the boxes in result
            //std::equal(result.begin(), result.begin(), temp.begin() );
            for(std::size_t j = 0; j < result.size(); j += temp.size())
            {
                std::size_t ins = 0;
                for(std::size_t k = 0; k < temp.size();k++)
                {
                    if(result[j+k] == temp[k])
                    {
                        ins++;
                    }
                }

                if(ins == temp.size())
                {   
                    //if the box is already in result do not insert it again
                    insert = false;
                    break;
                }
            }
            if(insert == true)
            {
                //insert the box if it is not in reasult
                for(std::size_t j = 0; j < temp.size();j++)
                {
                    result.push_back(temp[j]);
                }
            }
        }
        temp.clear();
    }

    return result;
}


template<class T>
void gsHierarchical_refinement<T>::computeErrors()
{
    m_pointErrors.clear();

    gsMatrix<T> val_i;
    //m_result->eval_into(m_param_values.col(0), val_i);
    m_result->eval_into(m_param_values, val_i);
    m_pointErrors.push_back( (m_points.row(0) - val_i.col(0).transpose()).norm() );
    m_max_error = m_min_error = m_pointErrors.back();
    
    for (index_t i = 1; i < m_points.rows(); i++)
    {
        //m_result->eval_into(m_param_values.col(i), val_i);

        const T err = (m_points.row(i) - val_i.col(i).transpose()).norm() ;

        m_pointErrors.push_back(err);

        if ( err > m_max_error ) m_max_error = err;
        if ( err < m_min_error ) m_min_error = err;
    }
}


};// namespace gismo
