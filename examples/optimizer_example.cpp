/** @file optmizer_example

    @brief Toy example for optimizer and optimization problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsGradientDescent.h>

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif

#ifdef gsIpOpt_ENABLED
#include <gsIpOpt/gsIpOpt.h>
#endif

using namespace gismo;

/** 
 * @brief 
 * Simple optimization example, to demonstrate the definition of an
 * optimization problem using the base class gsOptProblem.
 *
 *  This class implements the following NLP.
 *
 * min_x f(x) = -(x1-2)^2      (objective function)
 *  s.t.
 *       0 = x0^2 + x1 - 1     (constraint)
 *       -1 <= x0 <= 1         (variable bounds)
 *
 */


// To define an optimization problem we inherit from gsOptProblem class
// and implement the default constructor and few inherited virtual functions

//! [OptProblemExample Class]
template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    //! [OptProblemExample Constructor]
    // The constructor defines all properties of our optimization problem
    gsOptProblemExample()
    {
        // Number of design variables
        m_numDesignVars  = 2;
        // Number of constraints
        m_numConstraints = 1;
        // Number of non-zeros in the Jacobian of the constraints
        m_numConJacNonZero = 2;

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);
        
        // x0 has a lower bound of -1 and an upper bound of 1
        m_desLowerBounds[0] = -1.0;
        m_desUpperBounds[0] =  1.0;
        
        // x1 has no upper or lower bound, so we set them to
        // a large negative and a large positive number.
        // The value that is interpretted as -/+infinity can be
        // set in the options, but it defaults to -/+1e19
        m_desLowerBounds[1] = -1.0e19;
        m_desUpperBounds[1] =  1.0e19;

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // we have one equality constraint, so we set the bounds on
        // this constraint to be equal (and zero).
        m_conLowerBounds[0] = 
        m_conUpperBounds[0] = 0;
        
        // we initialize x in bounds, in the upper right quadrant
        m_curDesign.resize(2,1);
        m_curDesign(0,0) = 0.5;
        m_curDesign(1,0) = 1.5;        
        
        // 
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);

        // element at 0,0: grad_{x0} g_{1}(x)
        m_conJacRows[0] = 0;
        m_conJacCols[0] = 0;
        // element at 0,1: grad_{x0} g_{1}(x)
        m_conJacRows[1] = 0;
        m_conJacCols[1] = 1;
    }
    //! [OptProblemExample Constructor]

public:

    //! [OptProblemExample evalObj]
    // The evaluation of the objective function must be implemented
    T evalObj( const gsAsConstVector<T> & u ) const
    {
        // return the value of the objective function
        const T x1 = u(1,0);
        //return -(x1 - 2.0) * (x1 - 2.0);
        return (x1 - 2.0) * (x1 - 2.0);
    }
    //! [OptProblemExample evalObj]

    //! [OptProblemExample gradObj_into]
    // The gradient of the objective function (resorts to finite differences if left unimplemented)
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        result.resize(m_numDesignVars,1);

        // grad_{x0} f(x): x0 is not in the objective
        result(0,0)  = 0.0;

        // grad_{x1} f(x):
        const T x1 = u(1,0);
        //result(1,0)  =   -2.0*(x1 - 2.0);
        result(1,0)  =   2.0*(x1 - 2.0);
    }
    //! [OptProblemExample gradObj_into]

    //! [OptProblemExample evalCon_into]
    // The evaluation of the constraints must be implemented
    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        // return the value of the constraints: g(x)
        result.resize(m_numConstraints,1);

        const T x0 = u(0,0);
        const T x1 = u(1,0);
        result[0]  = -(x0*x0 + x1 - 1.0);
    }
    //! [OptProblemExample evalCon_into]

    //! [OptProblemExample jacobCon_into]
    // The Jacobian of the constraints (resorts to finite differences if left unimplemented)
    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        result.resize(m_numDesignVars,1);

        const T x0 = u(0,0);

        // element at 0,0: grad_{x0} g_{1}(x)
        result[0] = -2.0 * x0;
        
        // element at 0,1: grad_{x1} g_{1}(x)
        result[1] = -1.0;
    }
    //! [OptProblemExample jacobCon_into]

private:

    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};
//! [OptProblem]

int main(int argc, char* argv[])
{
    //! [Parse command line]
    index_t solver  = 0;

    gsCmdLine cmd("Demonstrates the use of optimizers.");
    cmd.addInt( "s", "solver", "Solver used. 0:gsGradientDescent, 1:gsHLBFGS, 2:IpOpt", solver);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]


    // Define an optimizer object
    gsInfo << "\nOptimization problem:";
    gsInfo << "\nmin_x f(x) = -(x1-2)^2      (objective function)";
    gsInfo << "\n  s.t.";
    gsInfo << "\n      0 = x0^2 + x1 - 1     (constraint)";
    gsInfo << "\n      -1 <= x0 <= 1         (variable bounds)\n\n";

    //! [Optimizer selection]
    gsOptProblemExample<real_t> problem;

    gsOptimizer<real_t> * optimizer;
    switch (solver)
    {
        case 0 :
        optimizer = new gsGradientDescent<>(&problem);
        break;

#ifdef gsHLBFGS_ENABLED
        case 1 :
        optimizer = new gsHLBFGS<real_t>(&problem);
        optimizer->options().setInt("Verbose",2);
        break;
#endif
#ifdef gsIpOpt_ENABLED
        case 2:
        optimizer = new gsIpOpt<real_t>(&problem);

        //Set the momentum rate used for the step calculation (default is 0.0).
        //Defines how much momentum is kept from previous iterations.
        // optimizer->setMomentum(4);

        // Turn verbosity on, so the optimizer prints status updates after each
        // iteration.
        // optimizer->options().setInt("Verbose",14);
        break;
#endif
        default:
        GISMO_ERROR("No optimizer defined for option "<<solver<<"\n");
    }
    //! [Optimizer selection]
                
    // NO EFFECT FOR gsIpOpt
    //! [Optimizer options]
    // Set number of iterations as stop criterion.
    // Set it to 0 or negative for infinite iterations (default is 0).
    optimizer->options().setInt("MaxIterations",200);

    // Set the minimum length of the gradient.
    // The optimizer stops minimizing if the gradient length falls below this
    // value (default is 1e-9).
    optimizer->options().setReal("MinGradientLength",1e-9);

    // Set the minimum length of the step.
    // The optimizer stops minimizing if the step length falls below this
    // value (default is 1e-9).
    optimizer->options().setReal("MinStepLength",1e-9);
    //! [Optimizer options]

    gsVector<> in(2);
    in << 0.5, 1.5;        

    // Start the optimization
    optimizer->solve(in);

    // Print final design info
    gsInfo << "\nNumber of iterations : " << optimizer->iterations() <<"\n";
    gsInfo << "Final objective value: " << optimizer->objective() <<"\n";
    gsInfo << "Final design: " << optimizer->currentDesign().transpose() <<"\n";

    delete optimizer;
    return EXIT_SUCCESS;
}
