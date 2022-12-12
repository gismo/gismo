/** @file gsTrilinosEigenProblems.cpp

    @brief Wrappers for Trilinos parallel eigenvalue solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsMpi/gsMpi.h>
#include <gsTrilinos/SparseMatrix.h>
#include <gsTrilinos/Vector.h>

#include <gsTrilinos/gsTrilinosEigenProblem.h>

// Include eigensolvers 
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

// Include header to define eigenproblem Ax = \lambda*x
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.
#include "AnasaziEpetraAdapter.hpp"
// Include header for Epetra sparse matrix and multivector.
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

// Include selected communicator class required by Epetra objects
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


namespace gismo
{
namespace trilinos
{
namespace solver
{


class EigenProblemPrivate
{
    friend class EigenProblem;
    
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;
    typedef Anasazi::BasicEigenproblem<double, MV, OP> BasicEigenProblem;
    typedef Anasazi::Eigenproblem<double, MV, OP> AnasaziProblem;
    typedef Anasazi::SolverManager<double, MV, OP> SolverManager;

    /// Solver manager
    Teuchos::RCP<SolverManager> manager;

    /// Problem data
    Teuchos::RCP<AnasaziProblem> problem;

    /// Parameter list
    Teuchos::ParameterList params;
};


EigenProblem::EigenProblem(const SparseMatrix & A,
                         const AnasaziMethod & method)
: my(new EigenProblemPrivate)
{
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Set eigensolver parameters.
    const double tol = 1.0e-8; // convergence tolerance
    const int nev = 1; // number of eigenvalues for which to solve
    const int blockSize = 1; // block size (number of eigenvectors processed at once)
    const int numBlocks = 5; // restart length
    const int maxRestarts = 100; // maximum number of restart cycles

    // Pass parameters list into the eigensolver.
    my->params.set ("Which", "SM"); //LM
    my->params.set ("Block Size",  blockSize);
    my->params.set ("Num Blocks", numBlocks);
    my->params.set ("Maximum Restarts", maxRestarts);
    my->params.set ("Convergence Tolerance", tol);
    my->params.set ("Verbosity", Anasazi::Errors + Anasazi::Warnings +
                    Anasazi::TimingDetails + Anasazi::FinalSummary);
    
    // Create a set of initial vectors to start the eigensolver.
    // This needs to have the same number of columns as the block size.
    const Epetra_Map & Map = A.get()->OperatorDomainMap();
    RCP<MV> ivec = rcp (new MV (Map, blockSize));
    ivec->Random();

    // Create the eigenproblem.  This object holds all the stuff about
    // your problem that Anasazi will see (matrix A and inital vectors)
    my->problem.reset( new EigenProblemPrivate::
                       BasicEigenProblem( A.getRCP(), ivec) );

    // Tell the eigenproblem that the operator A is symmetric.
    my->problem->setHermitian (true);

    // Set the number of eigenvalues requested
    my->problem->setNEV(nev);

    // Tell the eigenproblem that you are finishing passing it information.
    GISMO_ENSURE(my->problem->setProblem(),
                 "Anasazi::BasicEigenproblem::setProblem() returned an error.");

    // Create the eigensolver
    switch (method)
    {
    case LOBPCG:
        my->manager.reset(new Anasazi::
                          LOBPCGSolMgr<double, MV, OP>(my->problem, my->params) );
        break;
    case BlockDavidson:
        my->manager.reset(new Anasazi::
                          BlockDavidsonSolMgr<double, MV, OP>(my->problem, my->params) );

    case BlockKrylovSchur:
        my->manager.reset(new Anasazi::
                          BlockKrylovSchurSolMgr<double, MV, OP>(my->problem, my->params) );        
        break;
    default:
        GISMO_ERROR("Error, method choice "<< method);
    break;
    }
}

EigenProblem::~EigenProblem()
{
    delete my;
}

void EigenProblem::solve() const
{
    typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;
    using Teuchos::RCP;
    
#ifdef HAVE_MPI
        Epetra_MpiComm comm (gsMpi::init().worldComm() );
#else
        Epetra_SerialComm comm;
#endif
        
    // Solve the eigenvalue problem.
    //
    // Note that creating the eigensolver is separate from solving it.
    // After creating the eigensolver, you may call solve() multiple
    // times with different parameters or initial vectors.  This lets
    // you reuse intermediate state, like allocated basis vectors.
    Anasazi::ReturnType returnCode = my->manager->solve();
    if (returnCode != Anasazi::Converged && comm.MyPID() == 0)
        gsWarn << "Anasazi eigensolver did not converge.\n";
    
    // Get the eigenvalues and eigenvectors from the eigenproblem.
    Anasazi::Eigensolution<double,MV> sol = my->problem->getSolution();

    // Anasazi returns eigenvalues as Anasazi::Value, so that if
    // Anasazi's Scalar type is real-valued (as it is in this case), but
    // some eigenvalues are complex, you can still access the
    // eigenvalues correctly.  In this case, there are no complex
    // eigenvalues, since the matrix pencil is symmetric.
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    RCP<MV> evecs = sol.Evecs;

    // Compute residuals.
    std::vector<double> normR (sol.numVecs);
    if (sol.numVecs > 0)
    {
        Teuchos::SerialDenseMatrix<int,double> T (sol.numVecs, sol.numVecs);

        // getA(), getM()
        const Epetra_Map & Map = my->problem->getOperator()->OperatorDomainMap();
        MV tempAevec (Map, sol.numVecs);
        T.putScalar (0.0);
        for (int i=0; i<sol.numVecs; ++i)
            T(i,i) = evals[i].realpart;

        my->problem->getOperator()->Apply (*evecs, tempAevec);
        MVT::MvTimesMatAddMv (-1.0, *evecs, T, 1.0, tempAevec);
        MVT::MvNorm (tempAevec, normR);
    }

    // Print the results on MPI process 0.
    if (comm.MyPID() == 0)
    {
        gsInfo << "Solver manager returned "
               << (returnCode == Anasazi::Converged ? "converged." : "unconverged.")
               << "\n\n"
               << "------------------------------------------------------\n"
               << std::setw(16) << "Eigenvalue"
               << std::setw(18) << "Direct Residual"
               << "\n------------------------------------------------------\n";
        for (int i=0; i<sol.numVecs; ++i)
            gsInfo << std::setw(16) << evals[i].realpart
                   << std::setw(18) << normR[i] / evals[i].realpart <<"\n";
        gsInfo << "------------------------------------------------------\n";
    }

}

};// namespace solver
};// namespace trilinos
};// namespace gismo
