/** @file gsOptProblem.hpp

    @brief Provides implementation of an optimization problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#ifdef GISMO_WITH_IPOPT
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#endif

#include <gsIO/gsFileManager.h>

namespace gismo
{

class gsOptProblemPrivate
{
public:
#ifdef GISMO_WITH_IPOPT
    // Pointer to IpOpt interface
    Ipopt::SmartPtr<Ipopt::TNLP> tnlp;
#endif

};

#ifdef GISMO_WITH_IPOPT
    /** @brief Interface for IpOpt optimization problem
     *
     */
template<typename T>
class gsIpOptTNLP : public Ipopt::TNLP
{
    typedef Ipopt::Index                     Index;
    typedef Ipopt::Number                    Number;
    typedef Ipopt::SolverReturn              SolverReturn;
    typedef Ipopt::IpoptData                 IpoptData;
    typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;

    public:
        gsIpOptTNLP(gsOptProblem<T> & op) : m_op(op) { }

    public:
        /**@name Overloaded from TNLP */
        //@{
        /** Method to return some info about the nlp */
        bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                          Index& nnz_h_lag, IndexStyleEnum& index_style)
        {
            // gsDebug<<"Getting get_nlp_info.\n";

            n = m_op.m_numDesignVars;
            m = m_op.m_numConstraints;

            // Nonzeros in the constaint jacobian
            nnz_jac_g = m_op.m_numConJacNonZero;

            // hessian of the lagrangian not supported yet
            nnz_h_lag = 0;

            // index style for row/col entries:
            // C_STYLE: 0-based, FORTRAN_STYLE: 1-based
            index_style = C_STYLE;

            return true;
        }

        /** Method to return the bounds for my problem */
        bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                             Index m, Number* g_l, Number* g_u)
        {
            //gsDebug<<"Getting get_bounds_info.\n";

            // to do: { memcpy(target, start, numVals); }

            copy_n( m_op.m_desLowerBounds.data(), n, x_l );
            copy_n( m_op.m_desUpperBounds.data(), n, x_u );

            copy_n( m_op.m_conLowerBounds.data(), m, g_l );
            copy_n( m_op.m_conUpperBounds.data(), m, g_u );

            return true;
        }

        /** Method to return the starting point for the algorithm */
        bool get_starting_point(Index n, bool init_x, Number* x,
                                bool init_z, Number* z_L, Number* z_U,
                                Index m, bool init_lambda,
                                Number* lambda)
        {
            //gsDebug<<"Getting get_starting_point.\n";

            // Here, we assume we only have starting values for the design variables
            copy_n( m_op.m_curDesign.data(), n, x );
            return true;
        }


        /** Method to return the objective value */
        bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
        {
            //gsDebug<<"Getting eval_f.\n";

            gsAsConstVector<T> xx(x, n);
            obj_value = m_op.evalObj( xx );
            return true;
        }

        /** Method to return the gradient of the objective */
        bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
        {
            //gsDebug<<"Getting eval_grad_f.\n";

            gsAsConstVector<T> xx(x     , n);
            gsAsVector<T> result (grad_f, n);
            m_op.gradObj_into(xx, result);
            return true;
        }

        /** Method to return the constraint residuals */
        bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
        {
            //gsDebug<<"Getting eval_g.\n";

            gsAsConstVector<T> xx(x, n);
            gsAsVector<T> result(g, m);
            m_op.evalCon_into(xx, result);
            return true;
        }

        /** Method to return:
         *   1) The structure of the jacobian (if "values" is NULL)
         *   2) The values of the jacobian (if "values" is not NULL)
         */
        bool eval_jac_g(Index n, const Number* x, bool new_x,
                        Index m, Index nele_jac, Index* iRow, Index *jCol,
                        Number* values)
        {
            if (values == NULL)
            {
                // pass the structure of the jacobian of the constraints
                copy_n( m_op.m_conJacRows.data(), nele_jac, iRow );
                copy_n( m_op.m_conJacCols.data(), nele_jac, jCol );
            }
            else
            {
                gsAsConstVector<T> xx(x     , n      );
                gsAsVector<T>  result(values, nele_jac);
                m_op.jacobCon_into(xx, result);
            }

            return true;
        }


        /** Method to return:
         *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
         *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
         */
        bool eval_h(Index n, const Number* x, bool new_x,
                    Number obj_factor, Index m, const Number* lambda,
                    bool new_lambda, Index nele_hess, Index* iRow,
                    Index* jCol, Number* values)
        {
            GISMO_ERROR("IpOpt Hessian option not supported yet!");
        }

        //@}



    /** This method is called once per iteration, after the iteration
        summary output has been printed.  It provides the current
        information to the user to do with it anything she wants.  It
        also allows the user to ask for a premature termination of the
        optimization by returning false, in which case Ipopt will
        terminate with a corresponding return status.  The basic
        information provided in the argument list has the quantities
        values printed in the iteration summary line.  If more
        information is required, a user can obtain it from the IpData
        and IpCalculatedQuantities objects.  However, note that the
        provided quantities are all for the problem that Ipopt sees,
        i.e., the quantities might be scaled, fixed variables might be
        sorted out, etc.  The status indicates things like whether the
        algorithm is in the restoration phase...  In the restoration
        phase, the dual variables are probably not not changing.
    */
    virtual bool intermediate_callback(Ipopt::AlgorithmMode mode,
                                       Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
    {
        /*
        Ipopt::SmartPtr< const Ipopt::Vector >  curr  = ip_data->curr()->x();
        Ipopt::SmartPtr< Ipopt::DenseVector > dv = MakeNewDenseVector ();
        curr->Copy(dv);
        m_op.m_curDesign = gsAsConstVector<T>(dx->Values(),m_op.m_curDesign.rows());
        */
        // gsInfo << "\n === intermediateCallback is called === \n\n";
        return m_op.intermediateCallback();

        //SmartPtr< const IteratesVector >  trial = ip_data->trial();
        //int it = ip_data->iter_count();
        //Number 	tol = ip_data->tol();
    }


        /** @name Solution Methods */
        //@{
        /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
        void finalize_solution(SolverReturn status,
                               Index n, const Number* x, const Number* z_L, const Number* z_U,
                               Index m, const Number* g, const Number* lambda,
                               Number obj_value,
                               const IpoptData* ip_data,
                               IpoptCalculatedQuantities* ip_cq)
        {
            m_op.m_curDesign = gsAsConstVector<T>(x,n);
            m_op.m_lambda = gsAsConstVector<T>(lambda,m);

            //m_op.finalize();
        }

        //@}
    private:
        gsOptProblem<T> & m_op;
    private:
        gsIpOptTNLP(const gsIpOptTNLP & );
        gsIpOptTNLP& operator=(const gsIpOptTNLP & );
};
#endif


template <typename T>
gsOptProblem<T>::gsOptProblem()
{
    #ifdef GISMO_WITH_IPOPT

    m_data       =  new gsOptProblemPrivate();
    m_data->tnlp =  new gsIpOptTNLP<T>(*this);
    #else
    m_data = NULL;
    #endif
}

template <typename T>
gsOptProblem<T>::~gsOptProblem()
{
    delete m_data;
}

template <typename T>
void gsOptProblem<T>::gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    const index_t n = u.rows();
    //GISMO_ASSERT((index_t)m_numDesignVars == n*m, "Wrong design.");

    gsMatrix<T> uu = u;//copy
    gsAsVector<T> tmp(uu.data(), n);
    gsAsConstVector<T> ctmp(uu.data(), n);
    index_t c = 0;

    // for all partial derivatives (column-wise)
    for ( index_t i = 0; i!=n; i++ )
    {
        // to do: add m_desLowerBounds m_desUpperBounds check
        tmp[i]  += T(0.00001);
        const T e1 = this->evalObj(ctmp);
        tmp[i]   = u[i] + T(0.00002);
        const T e3 = this->evalObj(ctmp);
        tmp[i]   = u[i] - T(0.00001);
        const T e2 = this->evalObj(ctmp);
        tmp[i]   = u[i] - T(0.00002);
        const T e4 = this->evalObj(ctmp);
        tmp[i]   = u[i];
        result[c++]= ( 8 * (e1 - e2) + e4 - e3 ) / T(0.00012);
    }
}

template <typename T>
void gsOptProblem<T>::solve()
{
#ifdef GISMO_WITH_IPOPT

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(true);

    Ipopt::ApplicationReturnStatus status;
    std::string path = gsFileManager::findInDataDir( "options/ipopt.opt" );
    status = app->Initialize( path );

    if (status != Ipopt::Solve_Succeeded)
    {
        gsWarn << "\n\n*** Error during initialization!\n";
        return;
    }

    status = app->OptimizeTNLP(m_data->tnlp);
    //if (status != Ipopt::Solve_Succeeded)
    //   gsInfo << "Optimization did not succeed.\n";
    
    // Retrieve some statistics about the solve
    numIterations  = app->Statistics()->IterationCount();
    finalObjective = app->Statistics()->FinalObjective();
    //gsInfo << "\n*** The problem solved in " << numIterations << " iterations!\n";
    //gsInfo << "*** The final value of the objective function is " << finalObjective <<".\n";

#else

    GISMO_NO_IMPLEMENTATION
#endif

}


} // end namespace gismo
