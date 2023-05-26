/** @file gsSparseSolver.h

    @brief abstract interfaces for solvers and wrapper around Eigen solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris
*/

#pragma once

namespace gismo {

// forward declarations
template<typename T> class gsEigenCGIdentity;
template<typename T> class gsEigenCGDiagonal;
template<typename T> class gsEigenBiCGSTABIdentity;
template<typename T> class gsEigenBiCGSTABDiagonal;
template<typename T> class gsEigenBiCGSTABILUT;
template<typename T> class gsEigenSparseLU;
template<typename T> class gsEigenSparseQR;
template<typename T> class gsEigenSimplicialLDLT;
template<typename T> class gsEigenSimplicialDLT;

template<typename T> class gsEigenSuperLU;
template<typename T> class gsEigenPardisoLDLT;
template<typename T> class gsEigenPardisoLLT;
template<typename T> class gsEigenPardisoLU;

template<typename T> class gsEigenMINRES;
template<typename T> class gsEigenGMRES;
template<typename T> class gsEigenDGMRES;

/** @brief Abstract class for solvers.
    The solver interface is base on 3 methods:
    -compute set the system matrix (possibly compute the factorization or preconditioners)
    -solve solves for a given right hand side
    -succeed returns true if solving succeded according to solver dependent criteria
    (usually tolerance based)
    So in order to solve \f$ A x = b \f$ with a solver \a s two functions must be called:
    s.compute(A) and s.solve(b). The calls can be chained as in  s.compute(A).solve(b).


    Moreover, a collection of available sparse solvers is given as typedefs
    Example of usage:
    \code
    gsSparseMatrix<real_t> M; // sparse system matrix
    gsMatrix<real_t> b; // right-hand side
    gsSparseSolver<real_t>::CGDiagonal solver;
    solver.compute(M);
    gsMatrix<> x = solver.solve(b);
    \endcode
    
    See also
    http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
    and
    http://eigen.tuxfamily.org/dox/classEigen_1_1IterativeSolverBase.html

    \ingroup Matrix
*/
template <typename T = real_t>
class gsSparseSolver
{
public:
    typedef memory::unique_ptr<gsSparseSolver> uPtr;

    typedef gsEigenCGIdentity<T>           CGIdentity ;
    typedef gsEigenCGDiagonal<T>           CGDiagonal;
    typedef gsEigenBiCGSTABDiagonal<T>     BiCGSTABDiagonal;
    typedef gsEigenBiCGSTABIdentity<T>     BiCGSTABIdentity;
    typedef gsEigenBiCGSTABILUT<T>         BiCGSTABILUT;
    typedef gsEigenSparseLU<T>             LU;
    typedef gsEigenSparseQR<T>             QR;
    typedef gsEigenSimplicialLDLT<T>       SimplicialLDLT;
    typedef gsEigenSimplicialLDLT<T>       SimplicialLLT;

    // optionals
    typedef gsEigenSuperLU<T>              SuperLU;
    typedef gsEigenPardisoLDLT<T>          PardisoLDLT;
    typedef gsEigenPardisoLLT<T>           PardisoLLT;
    typedef gsEigenPardisoLU<T>            PardisoLU;

    typedef gsEigenMINRES<T>               MINRES;
    typedef gsEigenGMRES<T>                GMRES;
    typedef gsEigenDGMRES<T>               DGMRES;
    
public:
    typedef gsSparseMatrix<T> MatrixT;
    typedef gsMatrix<T>       VectorT;

public:
    virtual ~gsSparseSolver(){}

    virtual gsSparseSolver& compute (const MatrixT &matrix) = 0;

    virtual VectorT   solve   (const VectorT &rhs)    const = 0;

    virtual bool      succeed ()                      const = 0;

    virtual int info() const = 0;

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsSparseSolver\n";
        return os;
    }

    /// Prints the object as a string with extended details.
    virtual std::string detail() const 
    {
        std::ostringstream os;
        print(os);
        return os.str();
    }

    static uPtr get(const std::string & slv)
    {
        if (slv=="CGDiagonal")       return uPtr(new CGDiagonal());
        if (slv=="SimplicialLDLT")   return uPtr(new SimplicialLDLT());
        if (slv=="SimplicialLLT")   return uPtr(new SimplicialLLT());
#       ifdef GISMO_WITH_PARDISO
        if (slv=="PardisoLU")        return uPtr(new PardisoLU());
        if (slv=="PardisoLDLT")      return uPtr(new PardisoLDLT());
        if (slv=="PardisoLLT")       return uPtr(new PardisoLLT());
#       endif
#       ifdef GISMO_WITH_SUPERLU
        if (slv=="SuperLU")          return uPtr(new SuperLU());
#       endif
        if (slv=="BiCGSTABILUT")     return uPtr(new BiCGSTABILUT());
        if (slv=="BiCGSTABDiagonal") return uPtr(new BiCGSTABDiagonal());
        if (slv=="QR")               return uPtr(new QR());
        if (slv=="LU")               return uPtr(new LU());
        if (slv=="CGIdentity")       return uPtr(new CGIdentity());
        if (slv=="BiCGSTABIdentity") return uPtr(new BiCGSTABIdentity());
        // if (slv=="MINRES") return uPtr(new MINRES());
        // if (slv=="GMRES")  return uPtr(new GMRES());
        // if (slv=="DGMRES") return uPtr(new DGMRES());
        GISMO_ERROR("Solver \'"<< slv << "\' not known to G+Smo");
    }
};

/// \brief Print (as string) operator for sparse solvers
template<class T>
std::ostream &operator<<(std::ostream &os, const gsSparseSolver<T>& b)
{return b.print(os); }

#define GISMO_EIGEN_SPARSE_SOLVER(gsname, eigenName)                    \
    template<typename T>                                                \
    class gsname : public gsSparseSolver<T>, public gsEigenAdaptor<T>::eigenName \
    {                                                                   \
        typedef typename gsSparseSolver<T>::MatrixT MatrixT;            \
        typedef typename gsSparseSolver<T>::VectorT VectorT;            \
    protected:                                                          \
        index_t m_rows;                                                 \
        index_t m_cols;                                                 \
    public:                                                             \
        gsname()                                                        \
            : m_rows(0),m_cols(0)                                       \
        {}                                                              \
        gsname(const MatrixT &matrix)                                   \
            : gsEigenAdaptor<T>::eigenName(matrix), m_rows(matrix.rows()),m_cols(matrix.cols()) \
        {}                                                              \
        gsname& compute   (const MatrixT &matrix)                       \
        {                                                               \
            m_rows=matrix.rows();                                       \
            m_cols=matrix.cols();                                       \
            gsEigenAdaptor<T>::eigenName::compute(matrix);              \
            return *this;                                               \
        }                                                               \
        VectorT solve  (const VectorT &rhs) const                       \
        {                                                               \
            return gsEigenAdaptor<T>::eigenName::solve(rhs);            \
        }                                                               \
        bool succeed() const                                            \
        { return gsEigenAdaptor<T>::eigenName::info()==Eigen::Success;} \
        int info() const                                                \
        { return gsEigenAdaptor<T>::eigenName::info();}                 \
        index_t rows() const {return m_rows;}                           \
        index_t cols() const {return m_cols;}                           \
        std::ostream &print(std::ostream &os) const                     \
        {                                                               \
            os <<STRINGIFY(eigenName);                                  \
            return os;                                                  \
        }                                                               \
    };

GISMO_EIGEN_SPARSE_SOLVER (gsEigenCGIdentity,     CGIdentity)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenCGDiagonal,     CGDiagonal)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenBiCGSTABIdentity, BiCGSTABIdentity)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenBiCGSTABDiagonal, BiCGSTABDiagonal)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenBiCGSTABILUT,     BiCGSTABILUT)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenSparseLU,       SparseLU)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenSparseQR,       SparseQR)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenSimplicialLDLT, SimplicialLDLT)
GISMO_EIGEN_SPARSE_SOLVER (gsEigenSimplicialLLT, SimplicialLLT)

#ifdef GISMO_WITH_SUPERLU
    GISMO_EIGEN_SPARSE_SOLVER (gsEigenSuperLU, SuperLU)
#endif

#ifdef GISMO_WITH_PARDISO
    GISMO_EIGEN_SPARSE_SOLVER (gsEigenPardisoLDLT, PardisoLDLT)
    GISMO_EIGEN_SPARSE_SOLVER (gsEigenPardisoLLT, PardisoLLT)
    GISMO_EIGEN_SPARSE_SOLVER (gsEigenPardisoLU, PardisoLU)
#endif

//GISMO_EIGEN_SPARSE_SOLVER (gsEigenMINRES, MINRES)
//GISMO_EIGEN_SPARSE_SOLVER (gsEigenGMRES,  GMRES)
//GISMO_EIGEN_SPARSE_SOLVER (gsEigenDGMRES, DGMRES)


#undef GISMO_EIGEN_SPARSE_SOLVER

}
