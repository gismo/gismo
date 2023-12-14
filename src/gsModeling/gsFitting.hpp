/** @file gsFitting.hpp

    @brief Provides implementation of data fitting algorithms by least
    squares approximation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl, G. Kiss, A. Mantzaflaris, D. Mokris

*/

#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsAssembler/gsExpressions.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsNurbs/gsBSpline.h>
#include <gsTensor/gsTensorDomainIterator.h>



namespace gismo
{

template<class T>
gsFitting<T>::~gsFitting()
{
    if ( m_result!=nullptr )
      delete m_result;
}

template<class T>
gsFitting<T>::  gsFitting(gsMatrix<T> const & param_values,
                          gsMatrix<T> const & points,
                          gsBasis<T>  & basis)
{
    GISMO_ASSERT(points.cols()==param_values.cols(), "Pointset dimensions problem "<< points.cols() << " != " <<param_values.cols() );
    GISMO_ASSERT(basis.domainDim()==param_values.rows(), "Parameter values are inconsistent: "<< basis.domainDim() << " != " <<param_values.rows() );

    m_param_values = param_values;
    m_points = points;
    m_result = nullptr;
    m_basis = &basis;
    m_points.transposeInPlace();

    m_offset.resize(2);
    m_offset[0] = 0;
    m_offset[1] = m_points.rows();
}

template<class T>
gsFitting<T>::gsFitting(gsMatrix<T> const& param_values,
                        gsMatrix<T> const& points,
                        gsVector<index_t>  offset,
                        gsMappedBasis<2,T> & mbasis)
{
    m_param_values = param_values;
    m_points = points;
    m_result = nullptr;
    m_basis = &mbasis;
    m_points.transposeInPlace();
    m_offset = give(offset);
}

template<class T>
void gsFitting<T>::compute(T lambda)
{
  gsInfo << "START ---------------------------------------------- PDM coefs.\n";
    m_last_lambda = lambda;

    // Wipe out previous result
    if ( m_result!=nullptr )
        delete m_result;

    const int num_basis = m_basis->size();
    const short_t dimension = m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->domainDim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * m_basis->basis(0).degree(i) + 1 ) * 4;
    // TODO: improve by taking constraints nonzeros into account.
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis + m_constraintsRHS.rows(), dimension);
    m_B.setZero();  // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    assembleSystem(A_mat, m_B);


    // --- Smoothing matrix computation
    //test degree >=3
    if(lambda > 0)
      applySmoothing(lambda, A_mat);

    if(m_constraintsLHS.rows() > 0)
	extendSystem(A_mat, m_B);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>::BiCGSTABILUT solver( A_mat );

    if ( solver.preconditioner().info() != gsEigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";

        return;
    }
    // Solves for many right hand side  columns
    gsMatrix<T> x;

    x = solver.solve(m_B); //toDense()

    // If there were constraints, we obtained too many coefficients.
    x.conservativeResize(num_basis, gsEigen::NoChange);

    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve

    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give(x) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(x));

    gsInfo << "END ---------------------------------------------- PDM coefs.\n";
}




void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
  std::ofstream file(name.c_str());
  for(int  i = 0; i < matrix.rows(); i++)
  {
    for(int j = 0; j < matrix.cols(); j++)
    {
       std::string str = std::to_string(matrix(i,j));
       if(j+1 == matrix.cols())
       {
           file<<std::setprecision(10)<<str;
       }
       else
       {
           file<<std::setprecision(10)<<str<<',';
       }
    }
    file<<'\n';
  }
}


template<class T>
void gsFitting<T>::updateGeometry(gsMatrix<T> coefficients,
                                  gsMatrix<T> parameters)
{
  if (!this->m_result)
  {
    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give( coefficients ) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(coefficients));
  }
  else
  {
    this->m_result->coefs() = coefficients;
  }
  this->m_param_values = parameters;
  this->computeErrors();
}


template<class T>
void gsFitting<T>::initializeGeometry(const gsMatrix<T> & coefficients,
                                      const gsMatrix<T> & parameters)
{
  if (!this->m_result)
  {
    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give( coefficients ) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(coefficients));
  }
  else
  {
    this->m_result->coefs() = coefficients;
  }
  this->m_param_values = parameters;
}



template<class T>
void gsFitting<T>::compute_tdm(T lambda, T mu, T sigma, const std::vector<index_t> & interpIdx)
{
  gsInfo << "---------------------------------------------------------------------------------------------------------\n";
  gsInfo << "---------------------------------------------------------------------------------------------------------\n";
  gsInfo << "START compute_tdm(...)\n";

  time_t now = time(0);

    m_last_lambda = lambda;
    if ( !m_result )
    {
      gsInfo << "No existing geometry...\n";
        //compute(m_last_lambda);
        compute(0);
      gsInfo << "... now it does, as Ordinary Least Squares model.\n";
      gsMatrix<T> refCoefs = m_result->coefs();
      // gsInfo << "current geometry:\n" << *m_result << "\n";
      // gsInfo << "with coefs:\n" << refCoefs << "\n";

      // Wipe out previous result
      if ( m_result )
          delete m_result;

      if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
          m_result = bb->makeGeometry( give(refCoefs) ).release();
      else
          m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(refCoefs));
    }
    else
    {
      const index_t num_basis = m_basis->size();
      const index_t num_points = m_points.rows();
      const index_t dimension = m_points.cols();

      if( interpIdx.size() == 0)
      {
        gsInfo << "Input point cloud needs to be ordered:\n"
               << "interior points, south boundary curve, east boundary curve, north boundary curve, west boundary curve.\n";
      return;
      }

      // sn: compute the normals

      gsMatrix<T> refCoefs = m_result->coefs();

      gsWriteParaview(*m_result, std::to_string(now) + "_geo_in_ctdm", 10000);

      // gsInfo << "current geometry:\n" << *m_result << "\n";
      gsInfo << "with coefs:\n" << refCoefs << "\n";
      gsMatrix<T> pltcoefs = refCoefs.transpose();
      gsWriteParaviewPoints(pltcoefs, std::to_string(now) + "_coefs_in_ctdm");
      gsInfo << "and basis:\n" << m_result->basis()<< "\n";
      //gsInfo << "and params:\n" << this->returnParamValues();
      gsWriteParaviewPoints(this->returnParamValues(), std::to_string(now) + "_params_in_ctdm");

      // Wipe out previous result
      // if ( m_result )
          // delete m_result;

      // if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
          // m_result = bb->makeGeometry( give(refCoefs) ).release();
      // else
          // m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(refCoefs));

       // gsInfo << "It should be the same...\n";
       // gsDebugVar(*m_result);
       // gsInfo << "with same coefs:\n" << m_result->coefs() << "\n";
      //gsInfo << "and same params:\n" << this->returnParamValues()<< "\n";



      gsExprEvaluator<T> ev;
      auto G = ev.getMap(*m_result);
      gsMatrix<T> normals(3, m_param_values.cols());
      normals.setZero();
      gsSparseMatrix<T> N_diag(m_points.rows() * 3, m_points.rows());
      N_diag.setZero();
      //
      // gsDebugVar(m_param_values.col(0));
      // gsDebugVar(m_points.rows());
      // gsDebugVar(m_points.cols());
      //for(index_t j=0; j < interpIdx[0]; j++)

      //for(index_t j=0; j < m_param_values.cols(); j++)
      for(index_t j=0; j < interpIdx[0]; j++)
      {
        normals.col(j) = ev.eval(sn(G).normalized(), m_param_values.col(j));
        N_diag(j,j) = normals(0,j);
        N_diag(m_points.rows()+j,j) = normals(1,j);
        N_diag(2*m_points.rows()+j,j) = normals(2,j);
      }


      // nv: outer normals for boundary curves.
      // normals.col(j) = (point - surface); to treat the boundary.
      //if(false)
      {
      for(index_t j=interpIdx[0]; j < m_param_values.cols(); j++)
      {
        // south, east, north, west
        patchSide ps;
        if ( (interpIdx[1] <= j) && (j < interpIdx[2]) ) // east
        {
          GISMO_ENSURE(m_param_values(0,j)==1, internal::to_string(j) + "-th parameter NOT on EAST boundary");
          ps = patchSide(0,boundary::east);
          normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
        }
        else if ( (interpIdx[2] <= j) && (j < interpIdx[3]) )
        {
          GISMO_ENSURE(m_param_values(1,j)==1, internal::to_string(j) + "-th parameter NOT on NORTH boundary");
          ps = patchSide(0,boundary::north);
          normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
        }
        else if ( interpIdx[3] <= j)
        {
          GISMO_ENSURE(m_param_values(0,j)==0, internal::to_string(j) + "-th parameter NOT on WEST boundary");
          ps = patchSide(0,boundary::west);
          normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
        }
        else
        {
          GISMO_ENSURE(m_param_values(1,j)==0, internal::to_string(j) + "-th parameter NOT on SOUTH boundary");
          ps = patchSide(0,boundary::south);
          normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
        }


        N_diag(j,j) = normals(0,j);
        N_diag(m_points.rows()+j,j) = normals(1,j);
        N_diag(2*m_points.rows()+j,j) = normals(2,j);
      }// normals on boundary
      } //fi

      writeToCSVfile(internal::to_string(num_basis)+"_" + std::to_string(now)+"_normals.csv", normals);

      gsSparseMatrix<T> sparseColloc(m_points.rows(), num_basis);

      // gsInfo << "--------------------------------------------------------\n";
      // gsInfo << "space dimension = " << num_basis << "\n";
      // gsDebugVar(m_result->basis());
      sparseColloc = m_result->basis().collocationMatrix( m_param_values ) ;

      gsMatrix<T> tmp = sparseColloc;
      gsMatrix<T> Bb(m_points.rows() * 3, m_basis->size() * 3);
      Bb.setZero();
      Bb.block(0,0,tmp.rows(),tmp.cols()) = tmp;
      Bb.block(tmp.rows(),tmp.cols(),tmp.rows(),tmp.cols()) = tmp;
      Bb.block(2*tmp.rows(),2*tmp.cols(),tmp.rows(),tmp.cols()) = tmp;

      // writeToCSVfile(internal::to_string(num_basis)+"Bb.csv", Bb);

      gsSparseMatrix<T> B_mat(m_points.rows() * 3, m_basis->size() * 3);
      B_mat = Bb.sparseView();

      gsMatrix<T> X_tilde(m_points.rows() * 3, 1);
      // gsDebugVar(m_points.col(0));
      X_tilde << m_points.col(0), m_points.col(1), m_points.col(2);

      gsSparseMatrix<T> A_tilde(B_mat.cols(), B_mat.cols());
      gsSparseMatrix<T> Im(3*m_points.rows(), 3*m_points.rows());
      Im.setIdentity();
      gsMatrix<T> rhs(B_mat.cols(),1);

      if ( mu > 0 && sigma > 0)
      {
        gsInfo << mu << "*PDM + " << sigma << "*TDM.\n";
        A_tilde = B_mat.transpose() * ( mu * Im + sigma * N_diag * N_diag.transpose()) * B_mat;
        rhs =  B_mat.transpose() * ( mu * Im + sigma * N_diag * N_diag.transpose() ) * X_tilde ;
      }
      else if (mu > 0 && sigma == 0)
      {
        gsInfo << "PDM.\n";
        A_tilde = B_mat.transpose() *  B_mat;
        rhs = B_mat.transpose() * X_tilde ;
      }
      else if (sigma > 0 && mu == 0)
      {
        gsInfo << "TDM.\n";
        A_tilde = B_mat.transpose() * (  N_diag * N_diag.transpose() ) * B_mat;
        rhs =  B_mat.transpose() * ( N_diag * N_diag.transpose() ) * X_tilde ;
      }
      else
      {
        gsWarn<<  "No available fitting method. Aborting.\n";
        return;
      }

      gsSparseMatrix<T> m_G;
      m_G.resize(m_basis->size(), m_basis->size());
      gsSparseMatrix<T> G_mat(A_tilde.rows(), A_tilde.cols());
      gsMatrix<T> Gg(A_tilde.rows(), A_tilde.cols());
      Gg.setZero();

      if(lambda > 0)
      {
        gsInfo << "Penalization.\n";
        applySmoothing(lambda, m_G);
        Gg.block(0,0,m_G.rows(),m_G.cols()) = m_G;
        Gg.block(m_G.rows(),m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
        Gg.block(2*m_G.rows(),2*m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
        G_mat = Gg.sparseView();
        A_tilde += lambda * G_mat;
      }


      gsSparseMatrix<T> LM(A_tilde.rows(), A_tilde.cols());
      LM.setIdentity();
      A_tilde = A_tilde + 1e-7 * LM;

      gsInfo << "Solving the linear system.\n";
      A_tilde.makeCompressed();



      // gsMatrix<T> A_csv(A_tilde);
      // writeToCSVfile(internal::to_string(num_basis)+"_" + std::to_string(now)+"_A_tilde.csv", A_csv);
      // writeToCSVfile("X_tilde.csv", X_tilde);
      // writeToCSVfile("rhs.csv", rhs);
      // gsInfo << "Written\n";

      // typename gsSparseSolver<T>::BiCGSTABILUT solver( A_tilde );
      //
      // if ( solver.preconditioner().info() != gsEigen::Success )
      // {
      //     gsWarn<<  "The preconditioner failed. Aborting.\n";
      //
      //     return;
      // }

      typename gsSparseSolver<T>::QR solver(A_tilde);
      gsMatrix<T> sol_tilde = solver.solve(rhs); //toDense()
      gsInfo << "rank = " << solver.rank() / 3 << " = " << m_basis->size() << "\n";
      gsInfo << "Solved.\n";

      if (solver.info() != gsEigen::Success)
      {
        gsInfo << "QR: " << solver.lastErrorMessage();
      }

      // If there were constraints, we obtained too many coefficients.
      sol_tilde.conservativeResize(num_basis * 3, gsEigen::NoChange);

      // gsDebugVar(sol_tilde.rows());
      // gsDebugVar(sol_tilde.cols());

      gsMatrix<T> coefs_tilde(m_basis->size(), 3);
      for(index_t j=0; j<m_basis->size(); j++)
      {
        coefs_tilde(j,0) = sol_tilde(j);
        coefs_tilde(j,1) = sol_tilde(m_basis->size()+j);
        coefs_tilde(j,2) = sol_tilde(2*m_basis->size()+j);
      }

      gsInfo << "New coefs computed with TDM.\n";
      gsDebugVar(coefs_tilde);
      //gsInfo << "and params:\n" << this->returnParamValues() << "\n";
      // Wipe out previous result
      if ( m_result )
          delete m_result;

      if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
          m_result = bb->makeGeometry( give(coefs_tilde) ).release();
      else
          m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(coefs_tilde));

      gsInfo<< "Updated result.\n";
      //gsDebugVar(*m_result);

    }// fi

    gsWriteParaview(*m_result, std::to_string(now) + "_geo_out_ctdm", 10000);

    // gsInfo << "current geometry:\n" << *m_result << "\n";

    gsMatrix<T> pltoutc = m_result->coefs().transpose();
    gsWriteParaviewPoints(pltoutc, std::to_string(now) + "_coefs_out_ctdm");
    gsInfo << "and basis:\n" << m_result->basis()<< "\n";
    //gsInfo << "and params:\n" << this->returnParamValues();
    gsWriteParaviewPoints(this->returnParamValues(), std::to_string(now) + "_params_out_ctdm");

    gsInfo << "END compute_tdm(...)\n";
    gsInfo << "---------------------------------------------------------------------------------------------------------\n";
    gsInfo << "---------------------------------------------------------------------------------------------------------\n";
}



template<class T>
void gsFitting<T>::compute_tdmlm(T lambda, T lm, const std::vector<index_t> & interpIdx)
{
    m_last_lambda = lambda;
    if ( !m_result )
    {
      gsInfo << "No existing geometry...\n";
      compute(m_last_lambda);
      gsInfo << "... now it does.\n";
    }

    const index_t num_basis = m_basis->size();
    const index_t num_points = m_points.rows();
    const index_t dimension = m_points.cols();

    if( interpIdx.size() == 0)
    {
      gsInfo << "Input point cloud needs to be ordered:\n"
             << "interior points, south boundary curve, east boundary curve, north boundary curve, west boundary curve.\n";
    return;
    }



    gsMatrix<T> refCoefs = m_result->coefs();

    if ( m_result )
        delete m_result;

    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give(refCoefs) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(refCoefs));

    // sn: compute the normals
    // nv: outer normals for boundary curves ???
    gsExprEvaluator<T> ev;
    auto G = ev.getMap(*m_result);
    gsMatrix<T> normals(3, m_param_values.cols());
    normals.setZero();
    gsSparseMatrix<T> N_diag(m_points.rows() * 3, m_points.rows());
    N_diag.setZero();

    for(index_t j=0; j < interpIdx[0]; j++)
    {
      normals.col(j) = ev.eval(sn(G).normalized(), m_param_values.col(j));
      N_diag(j,j) = normals(0,j);
      N_diag(m_points.rows()+j,j) = normals(1,j);
      N_diag(2*m_points.rows()+j,j) = normals(2,j);
    }
    // writeToCSVfile(internal::to_string(num_basis)+"normals.csv", N_diag);

    // nv: outer normals for boundary curves.
    // if(false)
    {
    for(index_t j=interpIdx[0]; j < m_param_values.cols(); j++)
    {
      // south, east, north, west
      patchSide ps;
      if ( (interpIdx[1] <= j) && (j < interpIdx[2]) ) // east
      {
        GISMO_ENSURE(m_param_values(0,j)==1, internal::to_string(j) + "-th parameter NOT on EAST boundary");
        ps = patchSide(0,boundary::east);
        normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
      }
      else if ( (interpIdx[2] <= j) && (j < interpIdx[3]) )
      {
        GISMO_ENSURE(m_param_values(1,j)==1, internal::to_string(j) + "-th parameter NOT on NORTH boundary");
        ps = patchSide(0,boundary::north);
        normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
      }
      else if ( interpIdx[3] <= j)
      {
        GISMO_ENSURE(m_param_values(0,j)==0, internal::to_string(j) + "-th parameter NOT on WEST boundary");
        ps = patchSide(0,boundary::west);
        normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
      }
      else
      {
        GISMO_ENSURE(m_param_values(1,j)==0, internal::to_string(j) + "-th parameter NOT on SOUTH boundary");
        ps = patchSide(0,boundary::south);
        normals.col(j) = ev.evalBdr(nv(G).normalized(), m_param_values.col(j), ps);
      }


      N_diag(j,j) = normals(0,j);
      N_diag(m_points.rows()+j,j) = normals(1,j);
      N_diag(2*m_points.rows()+j,j) = normals(2,j);
    }// normals on boundary
  } //fi


    gsSparseMatrix<T> sparseColloc(m_points.rows(), num_basis);

    sparseColloc = m_result->basis().collocationMatrix( m_param_values ) ;

    gsMatrix<T> tmp = sparseColloc;
    gsMatrix<T> Bb(m_points.rows() * 3, m_basis->size() * 3);
    Bb.setZero();
    Bb.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    Bb.block(tmp.rows(),tmp.cols(),tmp.rows(),tmp.cols()) = tmp;
    Bb.block(2*tmp.rows(),2*tmp.cols(),tmp.rows(),tmp.cols()) = tmp;

    gsSparseMatrix<T> B_mat(m_points.rows() * 3, m_basis->size() * 3);
    B_mat = Bb.sparseView();

    gsMatrix<T> X_tilde(m_points.rows() * 3, 1);
    X_tilde << m_points.col(0), m_points.col(1), m_points.col(2);


    gsSparseMatrix<T> A_tilde(B_mat.cols(), B_mat.cols());
    gsSparseMatrix<T> Im(B_mat.cols(), B_mat.cols());
    Im.setIdentity();
    gsMatrix<T> rhs(B_mat.cols(),1);


    gsInfo << "TDMLM.\n";
    gsInfo << "( " << B_mat.rows() << " x " << B_mat.cols() << " )^t * (" << N_diag.rows() << " x " << N_diag.cols() << ")\n";
    A_tilde = B_mat.transpose() * (  N_diag * N_diag.transpose() ) * B_mat  + lm * Im;
    gsInfo << "( " << B_mat.rows() << " x " << B_mat.cols() << " )^t * (" << N_diag.rows() << " x " << N_diag.cols() << ") * ( " << N_diag.cols() << " x " << N_diag.rows() << ")\n";
    rhs =  B_mat.transpose() * ( N_diag * N_diag.transpose() ) * X_tilde;

    gsSparseMatrix<T> m_G;
    m_G.resize(m_basis->size(), m_basis->size());
    gsSparseMatrix<T> G_mat(A_tilde.rows(), A_tilde.cols());
    gsMatrix<T> Gg(A_tilde.rows(), A_tilde.cols());
    Gg.setZero();

    if(lambda > 0)
    {
      gsInfo << "Penalization.\n";
      applySmoothing(lambda, m_G);
      Gg.block(0,0,m_G.rows(),m_G.cols()) = m_G;
      Gg.block(m_G.rows(),m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
      Gg.block(2*m_G.rows(),2*m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
      G_mat = Gg.sparseView();
      A_tilde += lambda * G_mat;
    }

    gsInfo << "Solving the linear system.\n";
    A_tilde.makeCompressed();

    gsInfo << "System dimensions: " << A_tilde.rows() << " x " << A_tilde.cols() << "\n";
    typename gsSparseSolver<T>::QR solver(A_tilde);
    //QR CGIdentity CGDiagonal;
    //typename gsSparseSolver<T>::LU solver(A_tilde);
    gsMatrix<T> sol_tilde = solver.solve(rhs); //toDense()
    // gsInfo << "rank = " << solver.rank() / 3 << " = " << m_basis->size() << "\n";
    gsInfo << "Solved.\n";

    if (solver.info() != gsEigen::Success)
    {
      gsInfo << "problem with linear system.\n";
    }

    // gsDebugVar(sol_tilde.rows());
    // gsDebugVar(sol_tilde.cols());

    gsMatrix<T> coefs_tilde(m_basis->size(), 3);
    for(index_t j=0; j<m_basis->size(); j++)
    {
      coefs_tilde(j,0) = sol_tilde(j);
      coefs_tilde(j,1) = sol_tilde(m_basis->size()+j);
      coefs_tilde(j,2) = sol_tilde(2*m_basis->size()+j);
    }

    // gsDebugVar(coefs_tilde);

    // Wipe out previous result
    if ( m_result )
        delete m_result;

    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give(coefs_tilde) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(coefs_tilde));

    // gsDebugVar(*m_result);
}




template <class T>
bool gsFitting<T>::is_corner(gsMatrix<T> & p_domain,
                             gsVector<T> & param)
{
  bool corner_check = false;
  if( (math::abs(param(0) - p_domain(0,0)) < 1e-15) && (math::abs(param(1) - p_domain(1,0)) < 1e-15) ){
    // gsInfo << "param:\n" << param << "\n";
    // gsInfo << param(0,0) << " == " << p_domain(0,0) << "\n";
    // gsInfo << param(0,1) << " == " << p_domain(1,0) << "\n";
    corner_check = true;
  }
  else if( (math::abs(param(0) - p_domain(0,1)) < 1e-15) && (math::abs(param(1) - p_domain(1,0)) < 1e-15) ){
    corner_check = true;
  }
  else if( (math::abs(param(0) - p_domain(0,1)) < 1e-15) && (math::abs(param(1) - p_domain(1,1)) < 1e-15) ){
    corner_check = true;
  }
  else if( (math::abs(param(0) - p_domain(0,0)) < 1e-15) && (math::abs(param(1) - p_domain(1,1)) < 1e-15) ){
    corner_check = true;
  }
  else{
    corner_check = false;
  }
  return corner_check;
}

template <class T>
void gsFitting<T>::parameterCorrectionFixedBoundary(T accuracy,
                                                    index_t maxIter,
                                                    T mu, T sigma,
                                                    const std::vector<index_t>& interpIdx)
{
    index_t sepIndex = interpIdx[0];

    if ( !m_result )
    {
      compute(m_last_lambda);
    }

    for (index_t it = 0; it!=maxIter; ++it)
    {

      gsVector<T> newParam;
      gsMatrix<T> geoSupport = m_result->support();
// #     pragma omp parallel for default(shared) private(newParam)
      for (index_t i = 0; i < sepIndex; ++i)
      {
          // gsInfo << "Interior:\n" << newParam << "\n";
          const auto & curr = m_points.row(i).transpose();
          newParam = m_param_values.col(i);
          m_result->closestPointTo(curr, newParam, accuracy, true);

          // Decide whether to accept the correction or to drop it
          // gsInfo << "Correction:\n" << newParam << "\n";
          // gsInfo << "On the boundary ? " << "\n";
          // gsInfo << geoSupport(0,0) << "<" << newParam(0) << "<" << geoSupport(0,1) << "\n";
          // gsInfo << geoSupport(1,0) << "<" << newParam(1) << "<" << geoSupport(1,1) << "\n";
          // if( ( geoSupport(0,0) < newParam(0) ) && ( newParam(0) < geoSupport(0,1) ) && ( geoSupport(1,0) < newParam(1) ) && ( newParam(1) < geoSupport(1,1) ) )
          // {
            // gsInfo << "And we get here.\n";
          if ((m_result->eval(newParam) - curr).norm() < (m_result->eval(m_param_values.col(i))- curr).norm()){
            m_param_values.col(i) = newParam;
          }

          // (!) There might be the same parameter for two points
          // or ordering constraints in the case of structured/grid data
      }

      }

      // refit
      compute_tdm(m_last_lambda, mu, sigma, interpIdx);
}

template <class T>
void gsFitting<T>::parameterCorrectionSepBoundary(T accuracy,
                                                  index_t maxIter,
                                                  T mu, T sigma,
                                                  const std::vector<index_t>& interpIdx)
{
  index_t sepIndex = interpIdx[0];

    if ( !m_result )
      compute(m_last_lambda);

    for (index_t it = 0; it!=maxIter; ++it)
    {

      gsVector<T> newParam;
      gsMatrix<T> geoSupport = m_result->support();
// #     pragma omp parallel for default(shared) private(newParam)
      for (index_t i = 0; i < sepIndex; ++i)
      //for (index_t i = 1; i<m_points.rows()-1; ++i) //(!curve) skip first last pt
      {
          // gsInfo << "Interior:\n" << newParam << "\n";
          const auto & curr = m_points.row(i).transpose();
          newParam = m_param_values.col(i);
          m_result->closestPointTo(curr, newParam, accuracy, true);

          // Decide whether to accept the correction or to drop it
          // gsInfo << "Correction:\n" << newParam << "\n";
          // gsInfo << "On the boundary ? " << "\n";
          // gsInfo << geoSupport(0,0) << "<" << newParam(0) << "<" << geoSupport(0,1) << "\n";
          // gsInfo << geoSupport(1,0) << "<" << newParam(1) << "<" << geoSupport(1,1) << "\n";
          // if( ( geoSupport(0,0) < newParam(0) ) && ( newParam(0) < geoSupport(0,1) ) && ( geoSupport(1,0) < newParam(1) ) && ( newParam(1) < geoSupport(1,1) ) )
          // {
            // gsInfo << "And we get here.\n";
          if ((m_result->eval(newParam) - curr).norm() < (m_result->eval(m_param_values.col(i))- curr).norm()){
            m_param_values.col(i) = newParam;
          }
          //}

          // (!) There might be the same parameter for two points
          // or ordering constraints in the case of structured/grid data
      }

      // Correct the parameters, but keep them on the boundary


      for (index_t i = sepIndex; i < m_points.rows(); ++i)
      {
          gsVector<T> newBoundaryParam(1);
          gsVector<T> oldBoundaryParam(1);
          const auto & curr = m_points.row(i).transpose();
          newParam = m_param_values.col(i);

          //if (!is_corner(geoSupport, newParam)){

            // gsInfo << "West  : " << newParam(0) << "==" << geoSupport(0,0) << "\n";
            // gsInfo << "East  : " << newParam(0) << "==" << geoSupport(0,1) << "\n";
            // gsInfo << "South : " << newParam(1) << "==" << geoSupport(1,0) << "\n";
            // gsInfo << "North : " << newParam(1) << "==" << geoSupport(1,1) << "\n";


            gsVector<T> leftParam = m_param_values.col(i);
            gsVector<T> rightParam = m_param_values.col(i);

            if ( (i==sepIndex) && !is_corner(geoSupport, newParam) ){
              leftParam = m_param_values.col(m_points.rows()-1);
              rightParam = m_param_values.col(i+1);
            }

            if ( (i==m_points.rows()-1) && !is_corner(geoSupport, newParam) ){
              leftParam = m_param_values.col(i-1);
              rightParam = m_param_values.col(sepIndex);
            }

            if ( (i>sepIndex) && (i < m_points.rows()-1) && !is_corner(geoSupport, newParam) ){

              leftParam = m_param_values.col(i-1);
              rightParam = m_param_values.col(i+1);
            }

            if( math::abs(newParam(0) - geoSupport(0,0)) < 1e-15 )
            {
              // gsInfo << newParam(0) << " == " << geoSupport(0,0) << "\n";
              // gsInfo << "1 = West-boundary:\n" << newParam << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(1);


              newBoundaryParam << newParam(1);
              oldBoundaryParam << newParam(1);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "West-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";

              // assumption: boundary parameters in anti-clockwise order
              if( (leftParam(1) > newBoundaryParam(0)) && (newBoundaryParam(0) > rightParam(1)) ){
                newParam(1) = newBoundaryParam(0);
              }
              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }

            }
            else if ( math::abs(newParam(0) - geoSupport(0,1)) < 1e-15 )
            {
              // gsInfo << newParam(0) << " == " << geoSupport(0,1) << "\n";
              // gsInfo << "2 = East-boundary:\n" << newParam << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(2);


              newBoundaryParam << newParam(1);
              oldBoundaryParam << newParam(1);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "East-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";
              if( (leftParam(1) < newBoundaryParam(0)) && (newBoundaryParam(0) < rightParam(1)) ){
                newParam(1) = newBoundaryParam(0);
              }

              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }

            }
            else if ( math::abs(newParam(1) - geoSupport(1,0)) < 1e-15 )
            {
              // gsInfo << newParam(1) << " == " << geoSupport(1,0) << "\n";
              // gsInfo << "3 = South-boundary:\n" << newParam << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(3);


              newBoundaryParam << newParam(0);
              oldBoundaryParam << newParam(0);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "South-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";
              if( (leftParam(0) < newBoundaryParam(0)) && (newBoundaryParam(0) < rightParam(0)) ){
                newParam(0) = newBoundaryParam(0);
              }

              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }
            }
            else{
              // gsInfo << newParam(1) << " == " << geoSupport(1,1) << "\n";
              // gsInfo << i << "-th point:\n" << "\n";
              // gsInfo << rightParam(0)<< " < " << newParam(0) << " < " << leftParam(0) << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(4);

              // gsInfo << newParam(1) << "==" << geoSupport(1,1) << "\n";
              newBoundaryParam << newParam(0);
              oldBoundaryParam << newParam(0);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "North-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";
              // assumption: parameters in anti-clockwise order
              // gsInfo << "Projection:\n";
              // gsInfo << rightParam(0)<< " < " << newBoundaryParam(0) << " < " << leftParam(0) << "\n";
              if( (leftParam(0) > newBoundaryParam(0)) && (newBoundaryParam(0) > rightParam(0)) ){
                // gsInfo << "Before:\n" << newParam << "\n";
                newParam(0) = newBoundaryParam(0);
                // gsInfo << "After:\n" << newParam << "\n";
              }
              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }
            }


          //} // corner: to be moved or not.
      }

      // refit
      compute_tdm(m_last_lambda, mu, sigma, interpIdx);
    }
}

template <class T>
void gsFitting<T>::parameterProjectionSepBoundary(T accuracy,const std::vector<index_t>& interpIdx)
{

  gsInfo << "parameterProjectionSepBoundary(...)\n";
  if ( !m_result )
  {
    compute(m_last_lambda);
  }
//#       pragma omp parallel for default(shared) private(newParam)
  gsInfo << "Parameter projection for interior points...\n";
  index_t trackCorrection = 0;
  for (index_t i = 0; i < interpIdx[0]; ++i)
  {
    gsVector<T> newParam;
    const auto & curr = m_points.row(i).transpose();
    newParam = m_param_values.col(i);
    m_result->closestPointTo(curr, newParam, accuracy, true); // true: use initial point

    // Decide whether to accept the correction or to drop it
    if ((m_result->eval(newParam) - curr).norm()
            < (m_result->eval(m_param_values.col(i)) - curr).norm())
    {
      m_param_values.col(i) = newParam;
      trackCorrection += 1;
    }
  }
  gsInfo << "... applied to " << trackCorrection << " interior points, out of " << interpIdx[0] <<".\n";
      // parameter correction on boundary curves
      // south: (u,0)
  gsInfo << "Parameter correction for south boundary points...\n";
  trackCorrection = 0;
  for (index_t i = interpIdx[0]+1; i < interpIdx[1]; ++i)
  {
    gsVector<> newParam(1,1);
    gsVector<> oldParam(1,1);
    newParam(0,0) = m_param_values(0,i);
    oldParam(0,0) = m_param_values(0,i);
    const auto & curr = m_points.row(i).transpose();
    typename gsGeometry<T>::uPtr b = m_result->boundary(3); // south edge
    b->closestPointTo(curr, newParam, accuracy, true);

    if ((b->eval(newParam) - curr).norm()
            < (b->eval(oldParam) - curr).norm())
    {
    m_param_values(0,i) = newParam(0,0);
    trackCorrection += 1;
    }
  }
  gsInfo << "... applied to " << trackCorrection << " south points, out of " << interpIdx[1]-interpIdx[0] <<".\n";
  // east (1,v)
  for (index_t i = interpIdx[1]+1; i < interpIdx[2]; ++i)
  {
    gsVector<> newParam(1,1);
    gsVector<> oldParam(1,1);
    newParam(0,0) = m_param_values(1,i); // we consider the v of the i-th parameter
    oldParam(0,0) = m_param_values(1,i); // we consider the v of the i-th parameter
    const auto & curr = m_points.row(i).transpose();
    typename gsGeometry<T>::uPtr b = m_result->boundary(2); // east edge
    b->closestPointTo(curr, newParam, accuracy, true);

    if ((b->eval(newParam) - curr).norm()
      < (b->eval(oldParam) - curr).norm())
        m_param_values(1,i) = newParam(0,0);
  }
  //north (u,1)
  for (index_t i = interpIdx[2]+1; i < interpIdx[3]; ++i)
  {
    gsVector<> newParam(1,1);
    gsVector<> oldParam(1,1);
    newParam(0,0) = m_param_values(0,i); // we consider the u of the i-th parameter
    oldParam(0,0) = m_param_values(0,i); // we consider the u of the i-th parameter
    const auto & curr = m_points.row(i).transpose();
    typename gsGeometry<T>::uPtr b = m_result->boundary(4); // north edge
    b->closestPointTo(curr, newParam, accuracy, true);

    if ((b->eval(newParam) - curr).norm()
      < (b->eval(oldParam) - curr).norm())
      m_param_values(0,i) = newParam(0,0);
  }
  //west (0,v)
  for (index_t i = interpIdx[3]+1; i < m_points.rows(); ++i)
  {
    gsVector<> newParam(1,1);
    gsVector<> oldParam(1,1);
    newParam(0,0) = m_param_values(1,i); // we consider the v of the i-th parameter
    oldParam(0,0) = m_param_values(1,i); // we consider the v of the i-th parameter
    const auto & curr = m_points.row(i).transpose();
    typename gsGeometry<T>::uPtr b = m_result->boundary(1); // west edge
    b->closestPointTo(curr, newParam, accuracy, true);

    if ((b->eval(newParam) - curr).norm()
        < (b->eval(oldParam) - curr).norm())
        m_param_values(1,i) = newParam(0,0);
  }
}


template <class T>
void gsFitting<T>::parameterCorrectionSepBoundary_tdm(T accuracy,
                                                index_t maxIter,
                                                T mu, T sigma,
                                                const std::vector<index_t>& interpIdx)
{
  gsInfo << "parameterCorrectionSepBoundary_tdm(...)\n";
    if ( !m_result )
    {
      compute(m_last_lambda);
    }

    const index_t d = m_param_values.rows();
    const index_t n = m_points.cols();
    for (index_t it = 0; it<maxIter; ++it)
    {
      time_t now = time(0);
//#       pragma omp parallel for default(shared) private(newParam)
      gsWriteParaviewPoints(this->returnParamValues(), std::to_string(now) + "param_in_correction");
      gsInfo << "(a.) Projections.\n";
      parameterProjectionSepBoundary(accuracy, interpIdx);
      gsWriteParaviewPoints(this->returnParamValues(), std::to_string(now) + "param_out_correction");
      gsMatrix<T> pltcoefs = m_result->coefs().transpose();
      gsWriteParaviewPoints(pltcoefs, std::to_string(now) + "coefs_in_correction");
      gsInfo << "(b.) compute T DM coefs again;";
      compute_tdm(m_last_lambda, mu, sigma, interpIdx);
      pltcoefs = m_result->coefs().transpose();
      now = time(0);
      gsWriteParaviewPoints(pltcoefs, std::to_string(now) + "coefs_out_correction");
      gsInfo << "+++++++++++++++++++++++++++++++++++++++++++\n";
    }// step of PC


}

template <class T>
void gsFitting<T>::parameterCorrectionSepBoundary_pdm(T accuracy,
                                                  index_t maxIter,
                                                  const std::vector<index_t>& interpIdx)
{
  gsInfo << "parameterCorrectionSepBoundary_pdm(...)\n";
    if ( !m_result )
    {
      compute(m_last_lambda);
    }

    const index_t d = m_param_values.rows();
    const index_t n = m_points.cols();
    for (index_t it = 0; it<maxIter; ++it)
    {
//#       pragma omp parallel for default(shared) private(newParam)
      gsInfo << "(a) projections.\n";
      parameterProjectionSepBoundary(accuracy, interpIdx);
      gsInfo << "(b) P DM coefficients computation.\n";
      compute(m_last_lambda);
    }// step of PC
}


template <class T>
void gsFitting<T>::parameterCorrectionSepBoundary_tdmlm(T accuracy,
                                                  index_t maxIter,
                                                  T lm,
                                                  const std::vector<index_t>& interpIdx)
{
  gsInfo << "parameterCorrectionSepBoundary_tdmlm(...)\n";
    if ( !m_result )
    {
      compute(m_last_lambda);
    }

    const index_t d = m_param_values.rows();
    const index_t n = m_points.cols();
    for (index_t it = 0; it<maxIter; ++it)
    {
      gsInfo << "(a) projections.\n";
//#       pragma omp parallel for default(shared) private(newParam)
      parameterProjectionSepBoundary(accuracy, interpIdx);
      gsInfo << "(b) coefficients TDM LM computation.\n";
      compute_tdmlm(m_last_lambda, lm, interpIdx);
    }// step of PC
}






/*

for (index_t i = 0; i < interpIdx[0]; ++i)
{
  gsVector<> newParam;
  original.closestPointTo(X.col(i), newParam, 1e-10, true);
  params.col(i) = newParam;
}
// parameter correction on boundary curves
// south: (u,0)
for (index_t i = interpIdx[0]+1; i < interpIdx[1]; ++i)
{
  gsVector<> newParam;
  original.closestPointTo(X.col(i), newParam, 1e-10, true);
  params(0,i) = newParam(0,0);
}
// east (1,v)
for (index_t i = interpIdx[1]+1; i < interpIdx[2]; ++i)
{
  gsVector<> newParam;
  original.closestPointTo(X.col(i), newParam, 1e-10, true);
  params(1,i) = newParam(1,0);
}
//north (u,1)
for (index_t i = interpIdx[2]+1; i < interpIdx[3]; ++i)
{
  gsVector<> newParam;
  original.closestPointTo(X.col(i), newParam, 1e-10, true);
  params(0,i) = newParam(0,0);
}
//west (0,v)
for (index_t i = interpIdx[3]+1; i < X.cols(); ++i)
{
  gsVector<> newParam;
  original.closestPointTo(X.col(i), newParam, 1e-10, true);
  params(1,i) = newParam(1,0);
}
*/






template <class T>
void gsFitting<T>::parameterCorrection(T accuracy,
                                       index_t maxIter,
                                       T tolOrth)
{
    if ( !m_result )
        compute(m_last_lambda);

    const index_t d = m_param_values.rows();
    const index_t n = m_points.cols();

    for (index_t it = 0; it<maxIter; ++it)
    {
//#       pragma omp parallel for
        for (index_t i = 0; i<m_points.rows(); ++i)
        //for (index_t i = 1; i<m_points.rows()-1; ++i) //(!curve) skip first last pt
        {
            const auto & curr = m_points.row(i).transpose();
            gsVector<T> newParam = m_param_values.col(i);
            m_result->closestPointTo(curr, newParam, accuracy, true);

            // Decide whether to accept the correction or to drop it
            if ((m_result->eval(newParam) - curr).norm()
                    < (m_result->eval(m_param_values.col(i))
                        - curr).norm())
                    m_param_values.col(i) = newParam;

            // (!) There might be the same parameter for two points
            // or ordering constraints in the case of structured/grid data
        }

        // refit
        compute(m_last_lambda);
    }
}




// template<class T> inline
// gsSparseMatrix<T> gsBasis<T>::collocationMatrix(const gsMatrix<T> & u) const
// {
//     gsSparseMatrix<T> result( u.cols(), this->size() );
//     gsMatrix<T> ev;
//     gsMatrix<index_t> act;
//
//     eval_into  (u.col(0), ev);
//     active_into(u.col(0), act);
//     result.reservePerColumn( act.rows() );
//     for (index_t i=0; i!=act.rows(); ++i)
//         result.insert(0, act.at(i) ) = ev.at(i);
//
//     for (index_t k=1; k!=u.cols(); ++k)
//     {
//         eval_into  (u.col(k), ev );
//         active_into(u.col(k), act);
//         for (index_t i=0; i!=act.rows(); ++i)
//             result.insert(k, act.at(i) ) = ev.at(i);
//     }
//
//     result.makeCompressed();
//     return result;
// }
//
// template <class T>
// void gsFitting<T>::assembleTDMSystem(gsSparseMatrix<T>& A_mat,
//                                      gsMatrix<T>& m_rhs)
// {
//   // A_mat = B_mat^T * N * N^T * B_mat
//   // N = diag(n_x, n_y, n_z)
//   const int num_patches ( m_basis->nPieces() ); //initialize
//
//   //for computing the value of the basis function
//   gsMatrix<T> value, curr_point;
//   gsMatrix<index_t> actives;
//
//   for (index_t h = 0; h < num_patches; h++ )
//   {
//       auto & basis = m_basis->basis(h);
//
// //#   pragma omp parallel for default(shared) private(curr_point,actives,value)
//       for (index_t k = m_offset[h]; k < m_offset[h+1]; ++k)
//       {
//           curr_point = m_param_values.col(k); // compute the normal value to the point
//
//           //computing the values of the basis functions at the current point
//           basis.eval_into(curr_point, value);
//
//           // which functions have been computed i.e. which are active
//           basis.active_into(curr_point, actives);
//
//           const index_t numActive = actives.rows();
//
//           for (index_t i = 0; i != numActive; ++i)
//           {
//               const index_t ii = actives.at(i);
// //#           pragma omp critical (acc_m_B)
//               m_rhs.row(ii) += value.at(i) * m_points.row(k);
//               for (index_t j = 0; j != numActive; ++j)
// //#               pragma omp critical (acc_A_mat)
//                   A_mat(ii, actives.at(j)) += value.at(i) * value.at(j);
//           }
//       }
//   }
// }


template <class T>
void gsFitting<T>::assembleSystem(gsSparseMatrix<T>& A_mat,
                                  gsMatrix<T>& m_B)
{
    const int num_patches ( m_basis->nPieces() ); //initialize

    //for computing the value of the basis function
    gsMatrix<T> value, curr_point;
    gsMatrix<index_t> actives;

    for (index_t h = 0; h < num_patches; h++ )
    {
        auto & basis = m_basis->basis(h);

//#   pragma omp parallel for default(shared) private(curr_point,actives,value)
        for (index_t k = m_offset[h]; k < m_offset[h+1]; ++k)
        {
            curr_point = m_param_values.col(k);

            //computing the values of the basis functions at the current point
            basis.eval_into(curr_point, value);

            // which functions have been computed i.e. which are active
            basis.active_into(curr_point, actives);

            const index_t numActive = actives.rows();

            for (index_t i = 0; i != numActive; ++i)
            {
                const index_t ii = actives.at(i);
//#           pragma omp critical (acc_m_B)
                m_B.row(ii) += value.at(i) * m_points.row(k);
                for (index_t j = 0; j != numActive; ++j)
//#               pragma omp critical (acc_A_mat)
                    A_mat(ii, actives.at(j)) += value.at(i) * value.at(j);
            }
        }
    }
}

template <class T>
void gsFitting<T>::extendSystem(gsSparseMatrix<T>& A_mat,
				gsMatrix<T>& m_B)
{
    index_t basisSize = m_basis->size();

    // Idea: maybe we can resize here instead of doing it in compute().

    // This way does not work, as these block operations are read-only for sparse matrices.
    /*A_mat.topRightCorner(m_constraintsLHS.cols(), m_constraintsLHS.rows()) = m_constraintsLHS.transpose();
      A_mat.bottomLeftCorner(m_constraintsLHS.rows(), m_constraintsLHS.cols()) = m_constraintsLHS;*/
    m_B.bottomRows(m_constraintsRHS.rows()) = m_constraintsRHS;

    for (index_t k=0; k<m_constraintsLHS.outerSize(); ++k)
    {
	for (typename gsSparseMatrix<T>::InnerIterator it(m_constraintsLHS,k); it; ++it)
	{
	    A_mat(basisSize + it.row(), it.col()) = it.value();
	    A_mat(it.col(), basisSize + it.row()) = it.value();
	}
    }
}

template<class T>
gsSparseMatrix<T> gsFitting<T>::smoothingMatrix(T lambda) const
{
    m_last_lambda = lambda;

    const int num_basis=m_basis->size();

    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->domainDim(); ++i) // to do: improve
        nonZerosPerCol *= ( 2 * m_basis->basis(0).degree(i) + 1 );
    A_mat.reservePerColumn( nonZerosPerCol );
    const_cast<gsFitting*>(this)->applySmoothing(lambda, A_mat);
    return A_mat;
}

template<class T>
void gsFitting<T>::applySmoothing(T lambda, gsSparseMatrix<T> & A_mat)
{
    m_last_lambda = lambda;
    const int num_patches(m_basis->nPieces()); //initialize
    const short_t dim(m_basis->domainDim());
    const short_t stride = dim * (dim + 1) / 2;

    gsVector<index_t> numNodes(dim);
    gsMatrix<T> quNodes, der2, localA;
    gsVector<T> quWeights;
    gsMatrix<index_t> actives;

#   ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#   endif

    for (index_t h = 0; h < num_patches; h++)
    {
        auto & basis = m_basis->basis(h);

        //gsDebugVar(dim);
        //gsDebugVar(stride);

        for (short_t i = 0; i != dim; ++i)
        {
            numNodes[i] = basis.degree(i);//+1;
        }

        gsGaussRule<T> QuRule(numNodes); // Reference Quadrature rule

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();


#       ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#       else
        for (; domIt->good(); domIt->next() )
#       endif
        {
            // Map the Quadrature rule to the element and compute basis derivatives
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            basis.deriv2_into(quNodes, der2);
            basis.active_into(domIt->center, actives);
            const index_t numActive = actives.rows();
            localA.setZero(numActive, numActive);

            // perform the quadrature
            for (index_t k = 0; k < quWeights.rows(); ++k)
            {
                const T weight = quWeights[k] * lambda;

                for (index_t i = 0; i != numActive; ++i)
                    for (index_t j = 0; j != numActive; ++j)
                    {
                        T localAij = 0; // temporary variable

                        for (short_t s = 0; s < stride; s++)
                        {
                            // The pure second derivatives
                            // d^2u N_i * d^2u N_j + ...
                            if (s < dim)
                            {
                                localAij += der2(i * stride + s, k) *
                                    der2(j * stride + s, k);
                            }
                            // Mixed derivatives 2 * dudv N_i * dudv N_j + ...
                            else
                            {
                                localAij += T(2) * der2(i * stride + s, k) *
                                    der2(j * stride + s, k);
                            }
                        }
                        localA(i, j) += weight * localAij;
                    }
            }

            for (index_t i = 0; i != numActive; ++i)
            {
                const int ii = actives(i, 0);
                for (index_t j = 0; j != numActive; ++j)
                    A_mat(ii, actives(j, 0)) += localA(i, j);
            }
        }

    }
}

template<class T>
void gsFitting<T>::computeErrors()
{
    m_pointErrors.clear();

    gsMatrix<T> val_i;
    //->eval_into(m_param_values.col(0), val_i);
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


template<class T>
void gsFitting<T>::computeMaxNormErrors()
{
    m_pointErrors.clear();

    gsMatrix<T> values;
    m_result->eval_into(m_param_values, values);

    for (index_t i = 0; i != m_points.rows(); i++)
    {
        const T err = (m_points.row(i) - values.col(i).transpose()).cwiseAbs().maxCoeff();

        m_pointErrors.push_back(err);

        if ( i == 0 || m_max_error < err ) m_max_error = err;
        if ( i == 0 || err < m_min_error ) m_min_error = err;
    }
}



template<class T>
void gsFitting<T>::computeApproxError(T& error, int type) const

{
    gsMatrix<T> curr_point, results;

    const int num_patches(m_basis->nPieces());

    error = 0;

    for (index_t h = 0; h < num_patches; h++)
    {

        for (index_t k = m_offset[h]; k < m_offset[h + 1]; ++k)
        {
            curr_point = m_param_values.col(k);

            if (m_result)
                m_result->eval_into(curr_point, results);
            else
            {
                m_mresult.eval_into(h, curr_point, results);
            }

                const T err = (m_points.row(k) - results.transpose()).squaredNorm();

                switch (type) {
                case 0:
                    error += err;
                    break;
                case 1:
                    error += math::sqrt(err);
                    break;
                default:
                    gsWarn << "Unknown type in computeApproxError(error, type)...\n";
                    break;
                }
        }
    }
}

template<class T>
std::vector<T> gsFitting<T>::computeErrors(const gsMatrix<> & parameters,const gsMatrix<> & points)
{
  std::vector<T> min_max_mse;
  gsMatrix<> eval;
  m_result->eval_into(parameters, eval);

  gsMatrix<> pointWiseErrors(1, eval.cols());

  for (index_t col = 0; col != eval.cols(); col++)
  {
      pointWiseErrors(0, col) = (eval.col(col) - points.col(col)).norm();
  }

  real_t min_error = 1e6;
  real_t max_error = 0;
  real_t mse_error = 0;

  for (index_t i = 1; i < pointWiseErrors.cols(); i++)
  {
    const real_t err = pointWiseErrors(0,i) ;
    mse_error += err * err ;
    if ( err > max_error ) max_error = err;
    if ( err < min_error ) min_error = err;
  }

  min_max_mse.push_back(min_error);
  min_max_mse.push_back(max_error);
  min_max_mse.push_back(mse_error/pointWiseErrors.cols());

  return min_max_mse;
}

template<class T>
void gsFitting<T>::get_Error(std::vector<T>& errors, int type) const
{
    errors.clear();
    errors.reserve(m_points.rows());

    gsMatrix<T> curr_point, results;

    T err = 0;

    const int num_patches(m_basis->nPieces());

    for (index_t h = 0; h < num_patches; h++)
    {
        for (index_t k = m_offset[h]; k < m_offset[h + 1]; ++k)
        {
            curr_point = m_param_values.col(k);

            if (m_result)
                m_result->eval_into(curr_point, results);
            else
            {
                m_mresult.eval_into(h, curr_point, results);
            }

            results.transposeInPlace();

            err = (m_points.row(k) - results).template lpNorm<gsEigen::Infinity>();

                    switch (type)
                    {
                    case 0:
                        errors.push_back(err);
                        break;
                    case 1:
                        errors.push_back(math::sqrt(err));
                        break;
                    default:
                        gsWarn << "Unknown type in get_Error(errors, type)...\n";
                        break;
                    }
        }
    }
}

template<class T>
void gsFitting<T>::setConstraints(const std::vector<index_t>& indices,
				  const std::vector<gsMatrix<T> >& coefs)
{
    if(coefs.size() == 0)
	return;

    GISMO_ASSERT(indices.size() == coefs.size(),
		 "Prescribed indices and coefs must have the same length.");

    gsSparseMatrix<T> lhs(indices.size(), m_basis->size());
    gsMatrix<T> rhs(indices.size(), coefs[0].cols());

    index_t duplicates = 0;
    for(size_t r=0; r<indices.size(); r++)
    {
        index_t fix = indices[r];
        lhs(r-duplicates, fix) = 1;
        rhs.row(r-duplicates) = coefs[r];
    }

    setConstraints(lhs, rhs);
}

template<class T>
void gsFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides)
{
    if(fixedSides.size() == 0)
    return;

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;

    for(std::vector<boxSide>::const_iterator it=fixedSides.begin(); it!=fixedSides.end(); ++it)
    {
        gsMatrix<index_t> ind = this->m_basis->basis(0).boundary(*it);
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
    setConstraints(indices, coefs);
}

template<class T>
void gsFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides,
                      const std::vector<gsBSpline<T> >& fixedCurves)
{
    std::vector<gsBSpline<T> > tmp = fixedCurves;
    std::vector<gsGeometry<T> *> fixedCurves_input(tmp.size());
    for (size_t k=0; k!=fixedCurves.size(); k++)
        fixedCurves_input[k] = dynamic_cast<gsGeometry<T> *>(&(tmp[k]));
    setConstraints(fixedSides, fixedCurves_input);
}

template<class T>
void gsFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides,
                      const std::vector<gsGeometry<T> * >& fixedCurves)
{
    if(fixedSides.size() == 0)
    return;

    GISMO_ASSERT(fixedCurves.size() == fixedSides.size(),
         "fixedCurves and fixedSides are of different sizes.");

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;
    for(size_t s=0; s<fixedSides.size(); s++)
    {
    gsMatrix<T> coefsThisSide = fixedCurves[s]->coefs();
    gsMatrix<index_t> indicesThisSide = m_basis->basis(0).boundaryOffset(fixedSides[s],0);
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

    setConstraints(indices, coefs);
}


} // namespace gismo
