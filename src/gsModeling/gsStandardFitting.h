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
    gsMatrix<T> uvalue, vvalue, curr_point;
    gsMatrix<index_t> uactives, vactives;

    // 1D basis -- u direction
    gsBasis<T>* ubasis_tmp;
    ubasis_tmp = &(this->m_basis->component(0));
    gsBSplineBasis<T> ubasis = *(static_cast<gsBSplineBasis<T>*>(ubasis_tmp));
    // 1D basis -- v direction
    gsBasis<T>* vbasis_tmp;
    vbasis_tmp = &(this->m_basis->component(1));
    gsBSplineBasis<T> vbasis = *(static_cast<gsBSplineBasis<T>*>(vbasis_tmp));

    index_t flops = 0;
    for(index_t k = 0; k != num_points; ++k)
    {
        curr_point = this->m_param_values.col(k);
        gsMatrix<T> curr_u (1,1);
        curr_u(0,0) = curr_point(0);
        gsMatrix<T> curr_v (1,1);
        curr_v(0,0) = curr_point(1);

        //computing the values of the u- and v-basis functions at the current point
        ubasis.eval_into(curr_u, uvalue);
        vbasis.eval_into(curr_v, vvalue);

        // which functions have been computed i.e. which are active
        ubasis.active_into(curr_u, uactives);
        vbasis.active_into(curr_v, vactives);

        const index_t numActive = uactives.rows();

        for (index_t i1 = 0; i1 != numActive; ++i1)
        {
            const index_t ii1 = uactives.at(i1);
            for (index_t i2 = 0; i2 != numActive; i2++)
            {
                const index_t ii2 = vactives.at(i2);
                const index_t ii = ii2*N1+ii1;
                real_t val1 = uvalue.at(i1)*vvalue.at(i2);
                flops += 1;
                m_B.row(ii) += val1 * this->m_points.row(k);
                for (index_t j1 = 0; j1 != numActive; ++j1)
                {
                    const index_t jj1 = uactives.at(j1);
                    real_t val2 = uvalue.at(j1)*val1;
                    flops += 1;
                    for (index_t j2 = 0; j2 != numActive; j2++)
                    {
                        const index_t jj2 = vactives.at(j2);
                        const index_t jj = jj2*N1+jj1;
                        A_mat(ii, jj) += val2*vvalue.at(j2);
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
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis , num_basis );
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    // TODO: improve by taking constraints nonzeros into account.
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

    // --- Smoothing matrix computation
    //test degree >=3
    /*if(lambda > 0)
      this->applySmoothing(lambda, A_mat);*/

    /*if(m_constraintsLHS.rows() > 0)
    extendSystem(A_mat, m_B);*/

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
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

    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();
}

}
