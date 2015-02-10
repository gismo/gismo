#include <gsSolver/gsConjugateGradient.h>

namespace gismo
{

void gsConjugateGradient::initIteration(const gsConjugateGradient::VectorType& rhs, gsConjugateGradient::VectorType& x0, const gsLinearOperator& precond)
{
    GISMO_ASSERT(rhs.cols()== 1, "Implemented only for single column right hand side matrix");

    int n = m_mat.cols();
    int m = 1; // == rhs.cols();
    z.resize(n,m);
    tmp.resize(n,m);
    tmp2.resize(n,m);
    p.resize(n,m);

    m_mat.apply(x0,tmp2);  //apply the system matrix
    residual = rhs - tmp2; //initial residual

    precond.apply(residual, p);      //initial search direction

    absNew = Eigen::numext::real(residual.col(0).dot(p.col(0)));  // the square of the absolute value of r scaled by invM
    rhsNorm2 = rhs.squaredNorm();
    residualNorm2 = 0;
    threshold = m_tol*m_tol*rhsNorm2;
    m_numIter = 0;

}


bool gsConjugateGradient::step( gsConjugateGradient::VectorType& x, const gsLinearOperator& precond )
{
    m_mat.apply(p,tmp); //apply system matrix

    real_t alpha = absNew / p.col(0).dot(tmp.col(0));   // the amount we travel on dir
    x += alpha * p;                       // update solution
    residual -= alpha * tmp;              // update residual

    residualNorm2 = residual.squaredNorm();
    if(residualNorm2 < threshold)
        return true;

    precond.apply(residual, z);          // approximately solve for "A z = residual"

    real_t absOld = absNew;

    absNew = Eigen::numext::real(residual.col(0).dot(z.col(0)));     // update the absolute value of r
    real_t beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
    p = z + beta * p;                             // update search direction
    return false;
}

}
