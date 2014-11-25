
#pragma once

#include <map>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>// for BasisTraits -- to do: move it to the basis.
#include <gsNurbs/gsKnotVector.h>
#include <gsUtils/gsQuadrature.h>

#include <gsUtils/gsPointGrid.h>
#include <gsCore/gsDofMapper.h>

#include <gsAssembler/gsGaussAssembler.h>

namespace gismo {

template<class T>
class jacDetFunction : public gsFunction<T>
{
private:
    jacDetFunction( ){ } 
public:
    jacDetFunction(const gsGeometry<T> & geom)
        : m_geometry(geom)
    { } 
    
    int domainDim() const {return m_geometry.domainDim();}
    
    int targetDim() const { return 1;}
    
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    { 
        result.resize( 1, u.cols() );
        for ( index_t i = 0; i != u.cols(); ++i)
        {
            result(0,i) =  m_geometry.jac(u.col(i))->determinant();
        } 
    }
    
protected:
    const gsGeometry<T> & m_geometry;
};

template<class T>
class stiffnessFactor : public gsFunction<T>
{
private:
    stiffnessFactor( ){ } 
public:
    stiffnessFactor(const gsGeometry<T> & geom)
        : m_geometry(geom)
    { } 
    
    int domainDim() const                     
    {
        return m_geometry.domainDim();
    }
    
    int targetDim() const                     
    {
        return m_geometry.geoDim() * m_geometry.geoDim();
    }
    
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        int d  = this->m_geometry.parDim() ;
        result.resize(d*d, u.cols() );
        // Temporary storage
        gsMatrix<T> A;
        T detA ;
        
        for (index_t k=0; k!= u.cols(); ++k)
        {
            const gsMatrix<T> & uc = u.col(k);
            this->m_geometry.deriv_into( uc, A ); 
            detA= fabs( A.determinant() ) ;
            A = A.inverse();
            A = detA * ( A * A.transpose() );
            for (int i=0; i!= d; ++i)
                for (int j=0; j!= d ; ++j)
                    result( j*d +i,k) = A(i,j);
        }
    }
    
protected:
    const gsGeometry<T> & m_geometry;
};

template<class T> void
gsInterpolationAssembler<T>::initialize(const gsBasis<T>& m_basis)
{
    if ( ! m_geometryFactor )
    {
        //delete m_geometryFactor;
        m_geometryFactor = NULL;
    }
    if ( ! m_rhsApprox )
    {
        // delete m_rhsApprox;
        m_rhsApprox      = NULL;
    }

    if ( m_basis.dim()==2)
    {
        gsTensorBSplineBasis<2,T>  * basisptr = 
            dynamic_cast< gsTensorBSplineBasis<2,T>* >(// CHANGE
                const_cast<gsBasis<T>*>(&m_basis));
        if (! basisptr )
            gsWarn<<"Cannot assemble any basis other than Tensor B-spline basis."<<"\n";
    }
    else if ( m_basis.dim()==3)
    {
        gsTensorBSplineBasis<3,T>  * basisptr = 
            dynamic_cast< gsTensorBSplineBasis<3,T>* >(// CHANGE
                const_cast<gsBasis<T>*>(&m_basis));
        if (! basisptr )
            gsWarn<<"Cannot assemble any basis other than Tensor B-spline basis."<<"\n";
    }
    
    if ( (this->m_geometry) )
    {
        index_t d = this->m_geometry->parDim();
        m_p.resize(d);
        m_sz.resize(d);
        for (index_t i=0;i!=d;++i)
        {
            m_p[i] = m_basis.degree(i);
            m_sz[i]= m_basis.component(i).size() - 1;
            if ( int(m_basis.component(i).size()) < 3*m_p[i] )
                gsWarn<< "Basis too small in direction "<<i
                      <<" , boundary function will interact! case not covered yet.\n";
        }
        this->makeLut();
        
        m_initialized = true;
    }      
    else
        m_initialized = false;
}
 
template<class T>
gsSparseMatrix<T> * 
gsInterpolationAssembler<T>::stiffness(const gsBasis<T>& B )
{
    // Inititalize
    this->initialize(B);
        
    // Approximation step
    if ( ! m_geometryFactor )
        computeStiffnessFactor(B);

    unsigned bb= B.size();
    
    // Indices over lattice points
    gsVector<index_t> II(m_d), J0(m_d), J1(m_d), K0(m_d), K1(m_d), KK;
    gsMatrix<T> cw_int1(m_d*m_d, m_d);
    gsVector<T> comp(m_d*m_d), comp2(m_d*m_d);
    gsCombinat< gsVector<index_t> > k_iter, i_iter, j_iter;
    
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(bb,bb);
    K->setZero();
    
    unsigned nzRowsPerCol=1;
    m_h.resize(m_d);
    for (int i=0;i!=m_d;++i)
    {
        const int pp = B.component(i).degree();
        nzRowsPerCol *= 2 * pp  + 1;
        m_h[i] = B.component(i).domain()->breaks()[pp+1] - 
            B.component(i).domain()->breaks()[pp];
    }
    K->reserve( gsVector<T>::Constant(bb,nzRowsPerCol));
    
    T sum;
    const gsMatrix<T> & cp = this->m_geometryFactor->coefs();

    i_iter.first_lattice_point(m_sz, II);  
    do
    {
        const int i = tindex(II);
        overlapRange(II,J0,J1);
        // to do: combinat iterator to march on the lex> in the box
        j_iter.first_lattice_point(J0,J1,J0);
        do
        {  
            const int j = tindex(J0);
            // Construct only the lower diagonal part 
            if (j>i)
                continue;
            overlapRange(II,J0,K0,K1);
            
            KK = K0;
            // row - ith integral, column gradient product
            for ( int t = 0; t!=m_d; ++t)
                tri_integral(KK[t], II[t], J0[t], t, cw_int1.col(t) );
            sum= cw_int1.rowwise().prod().dot( cp.row( tindex(KK) ) );          
            do
            {                        
                int t(0);
                while ( KK[t] == K1[t]) { ++t; }
                KK[t]+=1;
                tri_integral(KK[t], II[t], J0[t], t, cw_int1.col(t) );
                
                --t;
                while ( t>=0 )
                {
                    KK[t]= K0[t];
                    tri_integral(KK[t], II[t], J0[t], t, cw_int1.col(t) );
                    --t;
                }

                sum += cw_int1.rowwise().prod().dot( cp.row( tindex(KK) ) );
            }
            while ( KK != K1 );
       
            K->coeffRef(i,j) = sum;
        }
        while ( j_iter.next_lattice_point(J0) );
    } 
    while ( i_iter.next_lattice_point(II) );
    
    K->makeCompressed();
    return K;
}
    
    
template<class T>
gsVector<T> * 
gsInterpolationAssembler<T>::moments( const gsBasis<T>& B, gsFunction<T> const& f )
{
    // Inititalize
    this->initialize(B);
    
    unsigned bb= B.size();
    
    // Approximation step
    this->computeRhsApprox(B,f);

    // Indices over lattice points
    gsVector<index_t> II, J0, J1, K0(m_d),K1(m_d);
    gsCombinat< gsVector<index_t> > k_iter, i_iter, j_iter;
    
    gsVector<T> * b = new gsVector<T>(bb);
    b->setZero();
    
    T sum, prod;
    int m,c;
    const gsMatrix<T> & rcp = this->m_rhsApprox->coefs();
    i_iter.first_lattice_point(m_sz, II);  
    do
    {
        const int i = tindex(II);
        overlapRange(II,J0,J1);
        j_iter.first_lattice_point(J0,J1,J0);
        sum=0;
        do
        {  
            const int j = tindex(J0);
            prod=T(1);  // * m_h[dir]
            for ( int k = 0; k!=m_d; ++k)
            {
                signature(J0[k],II[k],k,c,m);
                prod *=  (m_rhs_lookup[k])(c,m);
            }
            sum += rcp(j, 0) * prod;
        }
        while ( j_iter.next_lattice_point(J0) );
        (*b)[i] = sum;
    }
    while ( i_iter.next_lattice_point(II) );
    
    return b;
}
    

template<class T>
gsSparseSystem<T>
gsInterpolationAssembler<T>::assemblePoisson( const gsBasis<T>& m_basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsPoissonPde<T> & pde, int patchIndex)
{
    // Eigen provides special optimizations for statically
    // sized small matrices, so take advantage of that by
    // dispatching statically on dimensions 1-4
    switch (this->geometry().parDim())
    {
    //case 1:  return assemblePoisson_impl<1>      (m_basis, mapper, ddof, rhs, patchIndex);
    case 2:  return assemblePoisson_impl<2>      (m_basis, mapper, ddof, pde, patchIndex);
    case 3:  return assemblePoisson_impl<3>      (m_basis, mapper, ddof, pde, patchIndex);
    // case 4:  return assemblePoisson_impl<4>      (m_basis, mapper, ddof, rhs, patchIndex);
    // compiler trouble with dynamic dispatch and tensor basis.
    //default: return assemblePoisson_impl<Dynamic>(m_basis, mapper, ddof, rhs, patchIndex);
    default:
        return gsSparseSystem<T>();
    }
}
  
  
template<class T>
template <int D>
gsSparseSystem<T>
gsInterpolationAssembler<T>::assemblePoisson_impl( 
						  const gsBasis<T>& m_basis, 
						  const gsDofMapper& mapper, 
						  const gsMatrix<T> & ddof, 
						  const gsPoissonPde<T> & pde,
						  int patchIndex)
{
    // Inititalize
    this->initialize(m_basis);
    
    // Approximation step
    if ( ! m_geometryFactor )
        computeStiffnessFactor(m_basis, *pde.rhs() );

    gsVector<index_t, D> II(m_d), J0(m_d), J1(m_d), K0(m_d), K1(m_d);
    gsCombinat< gsVector<index_t, D> > k_iter, i_iter, j_iter;
  
    // initialize stiffness matrix and reserve space for entries
    int m_idofs= mapper.freeSize();
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(m_idofs,m_idofs) ;
    K->setZero();
    
    unsigned nzRowsPerCol=1;
    m_h.resize(m_d);
    for (int i=0;i!=m_d;++i)
    {
        const int pp = m_basis.component(i).degree();
        nzRowsPerCol *= 2 * pp  + 1;
        std::vector<T> brks = m_basis.component(i).domain()->breaks();
        m_h[i] = brks[pp+1] - brks[pp];
    }
    K->reserve( gsVector<T>::Constant(m_idofs, nzRowsPerCol) );
    const T h_prod = m_h.prod(); 

    // initialize right-hand side vector to zero
    gsVector<T> * b = new gsVector<T>(m_idofs) ; 
    b->setZero();
    
    T sum, prod, corr;
    int v,m,c;
    const gsMatrix<T> & cp  = this->m_geometryFactor->coefs();
    const gsMatrix<T> & rcp = this->m_rhsApprox->coefs();
    
    i_iter.first_lattice_point(m_sz, II);  
    do// For all II
    {
      const int ii = mapper.index( tindex(II), patchIndex );
      if ( mapper.is_free_index(ii) )
      {

          // Right hand side: buggy
          overlapRange(II,J0,J1);
          sum=T(0);
          j_iter.first_lattice_point(J0,J1,J0);
          do// For all JJ
          {  
              const int jj = mapper.index( tindex(J0), patchIndex );	      
              prod= h_prod;
              for ( int k = 0; k!=D; ++k)
              {
                  signature(J0[k],II[k],k,c,m);
                  prod *=  m_rhs_lookup[k](c,m);
              }
              sum += rcp(jj,0) *  prod ;
          }
          while ( j_iter.next_lattice_point(J0) );
          // Load vector
          (*b)[ii] +=  sum;

          overlapRange(II,J0,J1);
          j_iter.first_lattice_point(J0,J1,J0);
          do// For all JJ
          {  
              const int jj         = mapper.index( tindex(J0), patchIndex );	      
              const bool j_is_bd   = mapper.is_boundary_index(jj);
              
              // Exploit symmetry (construct only lower triangular part)
              if ( (jj<=ii) || j_is_bd )
              {
                  overlapRange(II,J0,K0,K1);	      
                  sum=T(0);
                  k_iter.first_lattice_point(K0,K1,K0);
                  do// For all KK
                  {
                      const int k0 = tindex(K0) ;
                      for ( int r = 0; r!=D; ++r)
                          for ( int s = 0; s!=D; ++s)
                          {
                              prod=T(1);
                              for ( int k = 0; k!=D; ++k)
                              {
                                  signature(K0[k],II[k],J0[k],k,r,s,v,c,m,corr);
                                  prod *=  corr * m_lookup[7*k+v](c,m);
                              }
                              sum += cp(k0, r*D+s) * prod;
                          }
                  }
                  while (k_iter.next_lattice_point(K0));
                  
                  if (j_is_bd)
                  {
                      // Dirichlet contribution to the load vector
                      (*b)[ii] -= ddof( mapper.global_to_bindex(jj) ) * sum;
                  }
                  else// jj is free and not bigger than ii
                  {
                      // Set the stiffness matrix entry
                      K->coeffRef(ii, jj) = sum;
                  }
		}// end if jj<=ii
      }
	  while ( j_iter.next_lattice_point(J0) );
      }// end ii is free
    }
    while ( i_iter.next_lattice_point(II) );
    
    K->makeCompressed();
    return gsSparseSystem<T>(K,b);
}
    
template<class T>
void gsInterpolationAssembler<T>::makeLut()
{   
    gsMatrix<T> ngrid;
    gsVector<T> wgrid;

    for ( index_t r = 0; r!=m_p.size(); ++r ) 
    {
        const int p = m_p[r];

        // Initialize tables
        for (int k=0; k!=7; ++k) 
        {
            m_lookup[7*r+k] = gsMatrix<T>::Zero((p+1)*(p+2)/2,p+1);
        }

        m_rhs_lookup[r] = gsMatrix<T>::Zero((p+1)*(p+2)/2,p+1);


        // Make univariate basis of degree p
        gsBSplineBasis<T> basis(gsKnotVector<T>(0,4*p+1,4*p,p+1));

        //gsDebug<<"lut knots= "<< basis.knots() <<"\n";
        
        // Get Gauss grid points and weights
        gsVector<int> int_nodes(1); 
        int_nodes << (int)ceil( (3.0*p+1.0) / 2.0 ) ;// exact integration
        
        std::vector<std::vector<T> > breaks;
        //Set integration breaks
        breaks.push_back( basis.domain()->breaks() ) ;
        tensorIteratedGaussRule(ngrid, wgrid, int_nodes, breaks);
        
        // Fill in LUTs by Gauss quadrature
        // evaluate basis functions and gradients and create blocks to refer to them
        gsMatrix<T> values;
        basis.evalAllDers_into(ngrid, 1, values);
        typename gsMatrix<T>::Block ev = values.topRows(p+1);
        typename gsMatrix<T>::Block dv = values.bottomRows(p+1);
        gsMatrix<unsigned> act;
        basis.active_into(ngrid,act); // indices

        for (index_t q=0; q!= ev.cols(); ++q)// for all quadrature points
        {
            for (index_t b=0; b!=p+1; ++b)
            {      
                if ( act(b,q) > unsigned(p) ) // base basis function
                    continue;
                const int m  = p - act(b,q) ; // multiplicity -1 of b	
                for (index_t i=0; i!=act.rows(); ++i)
                {
                    const int ii = act(i,q) - act(b,q);
                    if ( ii>=0 && ii<=p )
                    {
                        //m_lookup[7*r+7]->at(ii,m) += //111 (for rhs)
                        m_rhs_lookup[r](ii,m) += //111 (for rhs)
                            wgrid[q] * ev(b,q) * ev(i,q) ;

                        for (index_t j=0; j!=act.rows(); ++j)
                        {
                            const int jj = act(j,q) - act(b,q);      	    
                            if ( jj>=ii && jj<=p )
                            {
                                //c is the upper triangular index
                                const index_t c = (p+1)*ii+jj-(ii*(ii+1)/2);
                                m_lookup[7*r+0](c,m) += //000
                                    wgrid[q] * ev(b,q) * ev(i,q) * ev(j,q);
                                m_lookup[7*r+1](c,m) += //001
                                    wgrid[q] * ev(b,q) * ev(i,q) * dv(j,q);
                                m_lookup[7*r+2](c,m) += //010
                                    wgrid[q] * ev(b,q) * dv(i,q) * ev(j,q);
                                m_lookup[7*r+3](c,m) += //011
                                    wgrid[q] * ev(b,q) * dv(i,q) * dv(j,q);
                                m_lookup[7*r+4](c,m) += //100
                                    wgrid[q] * dv(b,q) * ev(i,q) * ev(j,q);
                                m_lookup[7*r+5](c,m) += //101
                                    wgrid[q] * dv(b,q) * ev(i,q) * dv(j,q);
                                m_lookup[7*r+6](c,m) += //110
                                    wgrid[q] * dv(b,q) * dv(i,q) * ev(j,q);
                            }
                        }
                    }
                }
            }
        }
    }
}
    


template<class T>
void gsInterpolationAssembler<T>::computeStiffnessFactor(const gsBasis<T>& m_basis, 
							 const gsFunction<T> & rhs)  const 
{
  if ( this->m_geometryFactor )
    delete this->m_geometryFactor;
  if ( this->m_rhsApprox )
    delete this->m_rhsApprox;
  
  int d  = this->geometry().parDim() ;
  
  // Temporary storage
  gsMatrix<T>  A;
  gsMatrix<T> uc;
  T detA ;
  
  std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(NEED_VALUE | NEED_JACOBIAN) );

  // Approximating the geometry factor
  gsMatrix<T> pts;
  m_basis.anchors_into(pts);
  gsMatrix<T> vals ( d*d ,pts.cols() );
  gsMatrix<T> rvals( 1,pts.cols() );
  
  for (index_t k=0; k!= pts.cols(); ++k)
  {
      geoEval->evaluateAt( pts.col(k) );
      A.noalias() = geoEval->jacobians();
      detA = fabs( A.determinant() ) ;
      A= A.inverse();
      A= detA * A * A.transpose();
      rhs.eval_into( geoEval->values(), uc );  
      rvals(0,k) = detA * uc(0,0);
      for (int i=0; i!= d; ++i)
          for (int j=0; j!= d ; ++j)
              vals( i*d +j,k) = A(i,j);
  }
  
  this->m_geometryFactor = m_basis.interpolate( vals );
  this->m_rhsApprox      = m_basis.interpolate( rvals );
}
  
template<class T>
void gsInterpolationAssembler<T>::computeStiffnessFactor(const gsBasis<T>& m_basis)  const 
{
    if ( this->m_geometryFactor )
        delete this->m_geometryFactor;
    
    gsMatrix<T> vals1;
    stiffnessFactor<T> Af(*this->m_geometry);
    Af.eval_into( *m_basis.anchors(), vals1);
    this->m_geometryFactor = m_basis.interpolate(vals1);
}

template<class T>
void gsInterpolationAssembler<T>::computeMassFactor(const gsBasis<T>& m_basis)  const 
{
    if ( this->m_geometryFactor )
        delete this->m_geometryFactor;
    
    gsMatrix<T> vals1;
    jacDetFunction<T> Af(*this->m_geometry);
    Af.eval_into( *m_basis.anchors(), vals1);
    this->m_geometryFactor = m_basis.interpolate( vals1 );
}

template<class T>
void gsInterpolationAssembler<T>::computeExactStiffnessFactor(const gsBasis<T>& m_basis)  const 
{
    if ( this->m_geometryFactor )
        delete this->m_geometryFactor;
    
    gsTensorBSplineBasis<2,T> b( dynamic_cast<const gsTensorBSplineBasis<2,T>&>(m_basis) );
    b.reduceContinuity(1);
    for ( int k = 0; k<2; k++)
        b.component(k).degreeElevate( b.degree(k) );

    stiffnessFactor<T> Af(*this->m_geometry);

    // this->m_geometryFactor = b.interpolate(Af); // equiv.
    gsMatrix<T> vals1;
    Af.eval_into( *m_basis.anchors(), vals1);
    this->m_geometryFactor = m_basis.interpolate( vals1 );
}

template<class T>
void gsInterpolationAssembler<T>::computeExactMassFactor(const gsBasis<T>& m_basis)  const 
{
    if ( this->m_geometryFactor )
        delete this->m_geometryFactor;
    
    gsTensorBSplineBasis<2,T> b( dynamic_cast<const gsTensorBSplineBasis<2,T>&>(m_basis) );
    b.reduceContinuity(1);
    for ( int k = 0; k<2; k++)
        b.component(k).degreeElevate( b.degree(k) );

    jacDetFunction<T> Af(*this->m_geometry);

    // this->m_geometryFactor = b.interpolate(Af); // equiv.
    gsMatrix<T> vals1;
    Af.eval_into( *m_basis.anchors(), vals1);
    this->m_geometryFactor = m_basis.interpolate( vals1 );
}

template<class T>
void gsInterpolationAssembler<T>::computeRhsApprox(const gsBasis<T>& m_basis,
						   const gsFunction<T> & rhs)  const 
{
    if ( this->m_rhsApprox )
        delete this->m_rhsApprox;
    
    // Temporary storage
    T tmp ;
  
    // Approximating the fixed factor
    gsMatrix<T> pts;
    m_basis.anchors_into(pts);
    gsMatrix<T> rvals( 1,pts.cols() );
    
    for (index_t k=0; k!= pts.cols(); ++k)
    {
        const gsMatrix<T> & uc = pts.col(k);
        tmp= fabs( this->m_geometry->jac(uc)->determinant() );
        rvals(0,k) = tmp * rhs.eval( this->geometry().eval( uc ) )->at(0,0);
    }
    
    this->m_rhsApprox = m_basis.interpolate( rvals );
}

template<class T>
T gsInterpolationAssembler<T>::interpolationError(int const & n_p) const
{
    // degree of the underlying Gauss rule to use
    const int Degree = 2;

    // compute the tensor Gauss rule
    gsMatrix<T> nodes;
    gsVector<T> weights;
    gsMatrix<T> range = this->m_geometry->parameterRange();
    uniformGaussRule<T>(nodes, weights, n_p, Degree, range.col(0), range.col(1));

    const int numPts = nodes.cols();

    // evaluate u and v
    gsMatrix<T> tmp, A, gval;
    T detA;

    // perform the quadrature
    T sum = 0.0;
    for (index_t k = 0; k < numPts; ++k)
    {
      const gsMatrix<T> & uc = nodes.col(k);
      this->m_geometryFactor->eval_into(uc,gval);
      A = this->geometry().jac( uc );
      detA= fabs( A.determinant() ) ;
      tmp = A.inverse();
      tmp = detA * (tmp * tmp.transpose() );

      sum += weights[k] * (gval.reshapeCol(0,m_d,m_d) - tmp).squaredNorm() ;
    }

    return std::sqrt(sum);
}

}; // namespace gismo
