#pragma once

#include <gismo.h>
#include <Eigen/Dense>
#include <gsModeling/gsFitting.h>
#include <gsModeling/gsFitting.hpp>
#include <gsSolver/gsSolverUtils.h>

namespace gismo {

template <class T>
class gsStandardFitting : public gsFitting<T>
{
public:
    gsStandardFitting(const gsMatrix<T> & param_values,
                  const gsMatrix<T> & points,
                  gsBasis<T> & basis)
        : gsFitting<T>(param_values, points, basis)
    { }
    void assembleSystemStandard(gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B);
    void computeStandard();
};

template<class T>
void gsStandardFitting<T>::assembleSystemStandard(gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B)
{
    index_t N1 = this->m_basis->component(0).size();    // number of basis functions in u direction
    const int num_points = this->m_points.rows();       // number of data points

    //for computing the value of the basis function
    gsMatrix<T> uvalues, vvalues;
    gsMatrix<index_t> uactives, vactives;

    // 1D basis -- u direction
    gsBasis<T>* ubasis;
    ubasis = &(this->m_basis->component(0));
    // 1D basis -- v direction
    gsBasis<T>* vbasis;
    vbasis = &(this->m_basis->component(1));

    // Evaluations and active basis at parameters
    ubasis -> eval_into(this->m_param_values.row(0),uvalues);
    vbasis -> eval_into(this->m_param_values.row(1),vvalues);
    ubasis -> active_into(this->m_param_values.row(0), uactives);
    vbasis -> active_into(this->m_param_values.row(1), vactives);

    index_t flops = 0;
    for(index_t k = 0; k != num_points; ++k)
    {
        const index_t numActive = uactives.rows();
        for (index_t i1 = 0; i1 != numActive; ++i1)
        {
            const index_t ii1 = uactives(i1,k);
            for (index_t i2 = 0; i2 != numActive; i2++)
            {
                const index_t ii2 = vactives(i2,k);
                const index_t ii = ii2*N1+ii1;
                real_t val1 = uvalues(i1,k)*vvalues(i2,k);
                flops += 1;
                m_B.row(ii) += val1 * this->m_points.row(k);
                for (index_t j1 = 0; j1 != numActive; ++j1)
                {
                    const index_t jj1 = uactives(j1,k);
                    real_t val2 = uvalues(j1,k)*val1;
                    flops += 1;
                    for (index_t j2 = 0; j2 != numActive; j2++)
                    {
                        const index_t jj2 = vactives(j2,k);
                        const index_t jj = jj2*N1+jj1;
                        A_mat(ii, jj) += val2*vvalues(j2,k);
                        flops += 2;
                    }
                }
            }
        }
    }
    //gsFileData<T> fd;
    //gsMatrix<T> C=A_mat.toDense();
    //fd << C;
    //fd.dump("SlowFittingMatrix");
    //gsFileData<> fb;
    //fb << m_B;
    //fb.dump("SlowFittingb");
    gsInfo << "Counted flops slow fitting        : " << flops << std::endl;
}

template<class T>
void gsStandardFitting<T>::computeStandard()
{
    // Wipe out previous result
    if (this-> m_result )
        delete this->m_result;

    const int num_basis=this->m_basis->size();
    const short_t dimension=this->m_points.cols();

    //left side matrix
    gsSparseMatrix<T> A_mat(num_basis , num_basis );
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i) // to do: improve
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis , dimension);
    m_B.setZero(); // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    gsStopwatch time;
    time.restart();
    assembleSystemStandard(A_mat, m_B);
    time.stop();
    gsInfo<<"Assembly time                     : "<< time <<"\n";

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    A_mat.makeCompressed();

    typename gsSparseSolver<T>::BiCGSTABILUT solver( A_mat );

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        this-> m_result = NULL;
        return;
    }
    // Solves for many right hand side  columns
    gsMatrix<T> x;
    x = solver.solve(m_B); //toDense()

    // If there were constraints, we obtained too many coefficients.
    x.conservativeResize(num_basis, Eigen::NoChange);

    // finally generate the B-spline curve
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();
}

}
