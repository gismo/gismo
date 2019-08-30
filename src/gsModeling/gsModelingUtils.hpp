/** @file gsModelingUtils.hpp

    @brief Utility functions required by gsModeling classes

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley
*/

#pragma once
#include <iostream>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsTensorBSpline.h>

#include <gsAssembler/gsGaussRule.h>


namespace gismo
{


template<class T>
gsMatrix<T> * innerProduct( const gsBasis<T>& B1, const gsBasis<T>& B2)
{
    gsMatrix<T> * K = new gsMatrix<T>(B1.size(), B2.size() ) ;
    K->setZero();
    
    int nGauss = int( ceil( double(B1.degree(0) + B2.degree(0) + 1)/2 ) );
    if (nGauss<1) nGauss=1;
    
    gsGaussRule<T> QuRule(nGauss); // Reference Quadrature rule
    gsMatrix<T> ngrid;          // tensor Gauss nodes
    gsVector<T> wgrid;          // tensor Gauss weights
    gsMatrix<index_t> act1, act2;
    gsMatrix<T>        ev1 , ev2;

    typename gsBasis<T>::domainIter domIt = B1.makeDomainIterator();
    for (; domIt->good(); domIt->next())
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), ngrid, wgrid );
  
        B1.eval_into(ngrid,ev1);
        B2.eval_into(ngrid,ev2);
        B1.active_into(ngrid,act1);
        B2.active_into(ngrid,act2);

        for (index_t k=0; k!= ngrid.cols(); ++k)
            for (index_t i=0; i!=act1.rows(); ++i)
                for (index_t j=0; j!=act2.rows(); ++j)
                    (*K)( act1(i,k) , act2(j,k) ) +=  wgrid(k) * ev1(i,k) * ev2(j,k) ;      
    }

    return K;
}

template<class T>
gsMatrix<T> * innerProduct1( const gsBasis<T>& B1, const gsBasis<T>& B2)
{
    gsMatrix<T> * K = new gsMatrix<T>(B1.size(), B2.size() ) ;
    K->setZero();

    int nGauss = int( ceil( double(B1.degree(0)-1 + B2.degree(0)-1 + 1)/2 ) );
    if (nGauss<1) nGauss=1;
    
    gsGaussRule<T> QuRule(nGauss); // Reference Quadrature rule
    gsMatrix<T> ngrid;          // tensor Gauss nodes
    gsVector<T> wgrid;          // tensor Gauss weights
    gsMatrix<index_t> act1, act2;
    gsMatrix<T>        ev1 , ev2;
    
    typename gsBasis<T>::domainIter domIt = B1.makeDomainIterator();
    for (; domIt->good(); domIt->next())
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), ngrid, wgrid );
        
        B1.deriv_into(ngrid,ev1);
        B2.deriv_into(ngrid,ev2);
        B1.active_into(ngrid,act1);
        B2.active_into(ngrid,act2);
        
        for (index_t k=0; k!= ngrid.cols(); ++k)
            for (index_t i=0; i!=act1.rows(); ++i)
                for (index_t j=0; j!=act2.rows(); ++j)
                    (*K)( act1(i,k) , act2(j,k) ) +=  wgrid(k) * ev1(i,k) * ev2(j,k) ;
        
    }
    
    return K;
}

template<class T>
gsMatrix<T> * innerProduct2( const gsBasis<T>& B1, const gsBasis<T>& B2)
{
    gsMatrix<T> * K = new gsMatrix<T>(B1.size(), B2.size() ) ;
    K->setZero();

    int nGauss = int( ceil( double(B1.degree(0)-2 + B2.degree(0)-2 + 1)/2 ) );
    if (nGauss<1) nGauss=1;

    gsGaussRule<T> QuRule(nGauss); // Reference Quadrature rule
    gsMatrix<T> ngrid;          // tensor Gauss nodes
    gsVector<T> wgrid;          // tensor Gauss weights
    gsMatrix<index_t> act1, act2;
    gsMatrix<T>        ev1 , ev2;
    
    typename gsBasis<T>::domainIter domIt = B1.makeDomainIterator();
    for (; domIt->good(); domIt->next())
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), ngrid, wgrid );
        
        B1.deriv2_into(ngrid,ev1);
        B2.deriv2_into(ngrid,ev2);
        B1.active_into(ngrid,act1);
        B2.active_into(ngrid,act2);
        
        for (index_t k=0; k!= ngrid.cols(); ++k)
            for (index_t i=0; i!=act1.rows(); ++i)
                for (index_t j=0; j!=act2.rows(); ++j)
                    (*K)( act1(i,k) , act2(j,k) ) +=  wgrid(k) * ev1(i,k) * ev2(j,k) ;      
    }
    
    return K;
}


/// intersection of two vectors
template <class T>
gsVector<T> vectorIntersect(gsVector<T> const & tangent1, 
                            gsVector<T> const & tangent2, 
                            gsMatrix<T> const & Vert1,
                            gsMatrix<T> const & Vert2)
{
    gsVector3d<T> abc1;
    gsVector3d<T> abc2;
    abc1(0) = tangent1(1);
    abc1(1) = -tangent1(0);
    abc1(2) = -( abc1(0)*Vert1(0,0) + abc1(1)*Vert1(0,1) );
    abc2(0) = tangent2(1);
    abc2(1) = -tangent2(0);
    abc2(2) = -( abc2(0)*Vert2(0,0) + abc2(1)*Vert2(0,1) );      
    gsVector<T> unknown(2);
    T detMatrixab = abc1(0)*abc2(1)-abc1(1)*abc2(0);
    if (detMatrixab!=0)
    {
        unknown(0) = -(1/detMatrixab)*(  abc2(1)*abc1(2)-abc1(1)*abc2(2) );
        unknown(1) = -(1/detMatrixab)*( -abc2(0)*abc1(2)+abc1(0)*abc2(2) );
    }
    else
    {
        unknown(0) = .5*Vert1(0) + .5*Vert2(0);
        unknown(1) = .5*Vert1(1) + .5*Vert2(1);
    };
    return unknown;
}

/// Angle between two vector: 0 <= angle <= pi
template <class T>
T conditionedAngle(gsVector3d<T> vec1, gsVector3d<T> vec2)
{   T dotp = vec1.dot(vec2)/( vec1.norm()*vec2.norm() );
    if ( dotp<-1 ) {
        dotp = -1;
        gsWarn<<"gsModelingUtils: truncations done for std::acos \n";
        return EIGEN_PI;}
    if ( dotp>1 ) {dotp = 1;
        gsWarn<<"gsModelingUtils: truncations done for std::acos \n";
        return 0;}

    T ag = math::acos(dotp);
    return ag;
}

/// Angle between two vector when viewing from a given normal vector
/// (the angle can be more than pi)
template <class T>
T conditionedAngle(gsVector3d<T> vec1, gsVector3d<T> vec2, gsVector3d<T> normal)
{
    T ag = conditionedAngle<T>(vec1,vec2);
    T cag = ( normal.dot( vec1.cross( vec2 ) ) >= 0 ) ? ag : T(2*EIGEN_PI-ag);
    return cag;
}


/// Find a critical point of a quadratic X^T A X subject to a linear
/// constraint C X = d.
template<class T>
gsVector<T> criticalPointOfQuadratic(gsMatrix<T> & A, gsMatrix<T> & C, gsVector<T> & d)
{
    index_t n = A.rows();
    index_t m = d.rows();
  
    gsMatrix<T> combined(n + m, n + m);
    combined.block(0, 0, n, n) = 2*A;
    combined.block(0, n, n, m) = C.transpose();
    combined.block(n, 0, m, n) = C;
    combined.block(n, n, m, m).setZero();
  
    gsVector<T> longD(n + m);
    longD.head(n).setZero();
    longD.tail(m) = d;
  
    gsMatrix<T> result = combined.fullPivLu().solve(longD);
    return result.block(0, 0, n, 1);
}

/// Find a critical point of a quadratic X^T A X + bX + c subject to a
/// linear constraint C X = d.
template<class T>
gsMatrix<T> criticalPointOfQuadratic(gsMatrix<T> const & A, gsMatrix<T> const & b, 
                                     gsMatrix<T> const & C, gsMatrix<T> const & d)
{
    index_t n = A.rows(); // dimension of X
    index_t m = d.rows(); // number of exact constraints
    assert(m<=n); // if not, the problem is ill defined
    assert(A.cols()==n); // A must be a square matrix
    assert(C.cols()==n);
    assert(d.cols()==1 && C.rows()==d.rows());
    gsMatrix<T> bt;
    if (b.rows()!=1) bt=b.transpose(); else bt=b;
  
    gsMatrix<T> combined(n + m, n + m);
    combined.block(0, 0, n, n) = 2*A;
    combined.block(0, n, n, m) = C.transpose();
    combined.block(n, 0, m, n) = C;
    combined.block(n, n, m, m).setZero();
  
    gsVector<T> longD(n + m);  
    longD.head(n) = -bt.transpose();
    longD.tail(m) = d;
  
    gsMatrix<T> result = combined.fullPivLu().solve(longD);
    return result.block(0, 0, n, 1);
}

/// Find a critical point of a quadratic X^T A X subject to a linear
/// constraint C X = d.
template<class T>
gsMatrix<T> criticalPointOfQuadratic(gsMatrix<T> const & A, gsMatrix<T> const & C, 
                                     gsMatrix<T> const & d)
{ 
    gsMatrix<T> b(1,A.cols()); b.setZero();
    return criticalPointOfQuadratic<T>(A,b,C,d); 
}

/// Find X which solves: min (AX-b)^T (AX-b), s.t. CX=d
template<class T>
gsMatrix<T> optQuadratic(gsMatrix<T> const & A, gsMatrix<T> const & b, 
                         gsMatrix<T> const & C, gsMatrix<T> const & d)
{ 
    return criticalPointOfQuadratic<T>( (A.transpose())*A, (-2)*( A.transpose() )*b, C, d);
}

/// Find X which solves: min w_1 (A_1 X-b_1)^T (A_1 X-b_1) + w_2 (A_2 X-b_2)^T (A_2 X-b_2), s.t. CX=d
template<class T>
gsMatrix<T> optQuadratic(gsMatrix<T> const & A1, gsMatrix<T> const & b1, 
                         T const & w1, gsMatrix<T> const & A2, 
                         gsMatrix<T> const & b2, T const & w2, 
                         gsMatrix<T> const & C, gsMatrix<T> const & d)
{ 
    return criticalPointOfQuadratic<T>( w1*(A1.transpose())*A1 + w2*(A2.transpose())*A2, 
                                        (-2)*w1*( A1.transpose() )*b1 + (-2)*w2*( A2.transpose() )*b2, C, d);
}
				      
/// Find X which solves: min w_1 (A_1 X-b_1)^T (A_1 X-b_1) + w_2 (A_2
/// X-b_2)^T (A_2 X-b_2) + w3 X'QX, s.t. CX=d
template<class T>
gsMatrix<T> optQuadratic(gsMatrix<T> const & A1, gsMatrix<T> const & b1, 
                         T const & w1, gsMatrix<T> const & A2, 
                         gsMatrix<T> const & b2, T const & w2, 
                         gsMatrix<T> const & C, gsMatrix<T> const & d,
                         T const & w3, gsMatrix<T> const & Q)
{ 
    return criticalPointOfQuadratic<T>( w1*(A1.transpose())*A1 + w2*(A2.transpose())*A2 + w3*Q, 
                                        (-2)*w1*( A1.transpose() )*b1 + (-2)*w2*( A2.transpose() )*b2, C, d);
}

/// Flip columes from left to right and vice versa
template <class T>
gsMatrix<T> flipLR(const gsMatrix<T> & mat)
{
    return mat.rowwise().reverse();
        
    // size_t const ncol = mat.cols();
    // gsMatrix<T> nMat(mat.rows(),ncol);  
    // for (size_t i=0; i<= ncol-1;  i++)
    // {
    //   nMat.col(i) = mat.col(ncol-1-i);
    // };
    // return nMat;
}

/// cross product of each pair of columes of two matrices and
/// normalize each columes
template <class T>
gsMatrix<T> crossNorm2Mat(gsMatrix<T> const & mat1,gsMatrix<T> const & mat2)
{
    assert(mat1.rows()==3 && mat2.rows()==3);
    assert(mat1.cols()==mat2.cols());
    size_t const nr = mat1.rows();
    size_t const nc = mat1.cols();
    gsMatrix<T> rcross(nr,mat1.cols()); rcross.setZero();
    gsMatrix<T> tem(3,1);
    for (size_t i=0; i<= nc-1; i++)
    {
        tem = ( gsVector3d<T>( mat1.col(i) ) ).cross( gsVector3d<T>( mat2.col(i) ) );   
        rcross.col(i) = tem / tem.norm();   
    };
    return rcross;
}

/// addConstraints
template <class T>
void addConstraints(gsMatrix<T> const & C1, gsMatrix<T> const & d1, 
                    gsMatrix<T> const & C2, gsMatrix<T> const & d2, 
                    gsMatrix<T> & C, gsMatrix<T> & d)
{
    int nr1 = C1.rows();
    int nc1 = C1.cols();  
    int nr2 = C2.rows();

    assert(nc1== C2.cols() );
    assert(nr1==d1.rows());
    assert(nr2==d2.rows());
    C.resize(nr1+nr2,nc1);
    C.block(0,0,nr1,nc1) = C1;
    C.block(nr1,0,nr2,nc1) = C2;
    d.resize(nr1+nr2,1);
    d.block(0,0,nr1,1) = d1;
    d.block(nr1,0,nr2,1) = d2;  
}

/// convert a with abs(a) < eps=2.220446049250313e-16 into 0
template <class T>
gsMatrix<T> convert2Zero(gsMatrix<T> const & mat)
{
    T eps=2.220446049250313e-16;
    
    gsMatrix<T> matc = mat;
    int n1 = mat.rows(); 
    int n2 = mat.cols();
    for (int i=0;i!=n1;i++)
    {
        for (int j=0;j!=n2;j++)
        {
            if (math::abs(mat(i,j))< eps) matc(i,j)=0.;
        }
    }
    return matc;
}

/// remove columes 0, nPoints, 2*nPoints,.. of a given matrix
template <class T>
void removeCol(gsMatrix<T> & mat, int const & removeEnds, int const & nPoints)
{
    assert(removeEnds==1 || removeEnds==2); 
    int nPeriod = mat.cols()/nPoints;
    assert( nPeriod*nPoints == mat.cols() );
    int ind1,ind2;
    if (removeEnds==1)
    {
        for (int i=nPeriod-1;i>=0;i--)
        {
            ind2 = i*nPoints + nPoints-1; 
            mat.removeCol(ind2);
        };
    };
    if (removeEnds==2)
    {
        for (int i=nPeriod-1;i>=0;i--)
        {
            ind2 = i*nPoints + nPoints-1;
            ind1 = i*nPoints ; //+0
            mat.removeCol(ind2);
            mat.removeCol(ind1);
        };
    };     
}

/**
   Interpolation with standard smoothing.
   TODO1: make the output as gsGeometry, gsBSpline for now; also use gsBasis as input
   TODO2: there should a different weight for approximating normal: w_nor
   Size of input matrices: each colummn represents a geometry point.
*/
template <class T>
gsBSpline<T> gsInterpolate(gsKnotVector<T> & kv,const gsMatrix<T> & preImage,
                           const gsMatrix<T> & image,
                           const gsMatrix<T> & preNormal,const gsMatrix<T> & normal,
                           const gsMatrix<T> & preImageApp,const gsMatrix<T> & imageApp,
                           T const & w_reg,T const & w_app,
                           gsMatrix<T> &outPointResiduals, gsMatrix<T> &outNormalResiduals)
{
    const int ntcp = kv.size()-kv.degree()-1;
    gsMatrix<T> tcp (ntcp, 2);

    // Quadratic forms which (approximately) constitute the beam strain energy
    gsBSplineBasis<T> bs(kv);
    gsMatrix<T> *Q = innerProduct2(bs, bs);

    // Exact constraints: point interpolation
    short_t dimPI = 1; // dimension of space of preImage, TODO: put dimPI, dimI to template<dimPI,...
    short_t dimI = 2;  // dimension of space of image
    GISMO_UNUSED(dimPI);
    int nip = image.cols(); // number of interpolating points
    int nn=normal.cols(); // number of prescribed normals
    gsMatrix<T> Nu, dNu, dNu_nm, NuApp;
    //--
    GISMO_ASSERT(dimPI==1 && dimI==2," "); // can be easily extended for other dimensions
    Nu   = bs.eval(preImage.row(0)); // u-BSplines
    dNu  = bs.deriv(preImage.row(0)); // u-BSplines
    Nu.transposeInPlace();
    dNu.transposeInPlace();
    dNu_nm  = bs.deriv(preNormal.row(0)); // u-BSplines for normals
    dNu_nm.transposeInPlace();
    gsMatrix<T> AdN = normal.row(0).asDiagonal() * dNu_nm;
    gsMatrix<T> BdN = normal.row(1).asDiagonal() * dNu_nm;

    // Approximate constraints
    NuApp = bs.eval(preImageApp.row(0));
    NuApp.transposeInPlace();
    gsMatrix<T> X0 = imageApp.row(0);
    gsMatrix<T> Y0 = imageApp.row(1);

    //-- resulting Saddle point linear System
    int nss = dimI*ntcp + dimI*nip + nn;
    gsMatrix<T> Ass(nss,nss);
    gsMatrix<T> bss(nss,1);
    Ass.setZero(); bss.setZero();
    gsMatrix<T> Hess = w_reg*2*(*Q) + w_app*2*(NuApp.transpose())* NuApp;
    //--- row 0
    Ass.block(0,0,ntcp,ntcp) = Hess;
    Ass.block(0,2*ntcp,ntcp,nip) = Nu.transpose();
    Ass.block(0,2*ntcp+2*nip,ntcp,nn) = AdN.transpose();
    bss.block(0,0,ntcp,1) = w_app*2*NuApp.transpose()*X0.transpose();
    //--- row 1
    Ass.block(ntcp,ntcp,ntcp,ntcp) = Hess;
    Ass.block(ntcp,2*ntcp+nip,ntcp,nip) = Nu.transpose();
    Ass.block(ntcp,2*ntcp+2*nip,ntcp,nn) = BdN.transpose();
    bss.block(ntcp,0,ntcp,1) = w_app*2*NuApp.transpose()*Y0.transpose();
    //--- row 2
    Ass.block(2*ntcp,0,nip,ntcp) = Nu;
    bss.block(2*ntcp,0,nip,1) = (image.row(0)).transpose();
    //--- row 3
    Ass.block(2*ntcp+nip,ntcp,nip,ntcp) = Nu;
    bss.block(2*ntcp+nip,0,nip,1) = (image.row(1)).transpose();
    //--- row 4
    Ass.block(2*ntcp+2*nip,0,nn,ntcp) = AdN;
    Ass.block(2*ntcp+2*nip,ntcp,nn,ntcp) = BdN;

    gsMatrix<T> result = Ass.fullPivLu().solve(bss);
    tcp.col(0) = result.block(0   , 0, ntcp, 1);
    tcp.col(1) = result.block(ntcp, 0, ntcp, 1);

//    gsDebug<< " parameterRange: \n"<< *trimLoop[sourceID]->basis().parameterRange()<<"\n";
//    gsDebug<<" Ass: \n"<<Ass<<"\n";
//    gsDebug<<" bss: \n"<<bss<<"\n";
//    gsDebug<<" result: \n"<<result<<"\n";
//    gsDebug<<" Q: \n"<< *Q<<"\n";
//    gsDebug<<" preimage: \n"<< preImage<<"\n";
//    gsDebug<<" prenormal: \n"<< preNormal<<"\n";
//    gsDebug<<" image: \n"<< image<<"\n";
//    gsDebug<<" normal: \n"<< normal<<"\n";
//    gsDebug<<" Nu: \n"<< *Nu<<"\n";
//    gsDebug<<" dNu: \n"<< *dNu<<"\n";
//    gsDebug<<" AdN: \n"<< AdN<<"\n";
//    gsDebug<<" BdN: \n"<< BdN<<"\n";
//    gsDebug<<" tcp: \n"<< tcp<<"\n";
//    gsDebug<<" preimageApp: \n"<< preImageApp<<"\n";
//    gsDebug<<" imageApp: \n"<< imageApp<<"\n";
//    gsDebug<<" residual of app x constraints: \n"<< *NuApp*tcp.col(0)-imageApp.row(0).transpose()<<std::endl;
//    gsDebug<<" residual of app y constraints: \n"<< *NuApp*tcp.col(1)-imageApp.row(1).transpose()<<std::endl;
//    gsDebug<<" residual of normal constraints: \n"<< AdN*tcp.col(0)+BdN*tcp.col(1)<<std::endl;

    outPointResiduals = (NuApp * tcp).transpose() - imageApp;
    outNormalResiduals = AdN * tcp.col(0) + BdN * tcp.col(1);
    //gsDebug << std::flush;

    delete Q;

    gsBSpline<T> tcurve(kv, give(tcp));

    return tcurve;
}


/// Create a surface (as a tensor product B-spline) satisfying conditions:
/// The evaluation of the surface at the columns of \a exactPoints are
/// constrained to equal the columns of \a exactValues. The evaluation
/// at the columns of \a appxPointsEdges (resp \a appxPointsInt) should be
/// approximately equal to the columns of \a appxValuesEdges (resp
/// \a appxValuesInt) with weighting \a wEdge (resp \a wInt). The normals
/// to the surface, evaluated at the columns of \a appxNormalPoints, should
/// be approximately equal to the columns of \a appxNormals with weighting
/// \a wNormal. Finally you can add a weighting \a wReg for the regularity
/// of the surface. The parameter \a force_normal is for a special case;
/// typically it should be set to false.
template<class T>
typename gsTensorBSpline<2,T>::Ptr gsInterpolateSurface(
    const gsMatrix<T> &exactPoints, const gsMatrix<T> &exactValues,
    const gsMatrix<T> &appxPointsEdges, const gsMatrix<T> &appxValuesEdges,
    const gsMatrix<T> &appxPointsInt, const gsMatrix<T> &appxValuesInt,
    const gsMatrix<T> &appxNormalPoints, const gsMatrix<T> &appxNormals,
    T wEdge, T wInt, T wNormal, T wReg,
    const gsKnotVector<T> &kv1, const gsKnotVector<T> &kv2,
    bool force_normal
    )
{
    GISMO_ASSERT(exactPoints.rows() == 2 && exactValues.rows() == 3, "Matrix input has incorrect dimension");
    GISMO_ASSERT(exactPoints.cols() == exactValues.cols(), "Matrix input has incorrect dimension");
    GISMO_ASSERT(appxPointsEdges.rows() == 2 && appxValuesEdges.rows() == 3, "Matrix input has incorrect dimension");
    GISMO_ASSERT(appxPointsEdges.cols() == appxValuesEdges.cols(), "Matrix input has incorrect dimension");
    GISMO_ASSERT(appxPointsInt.rows() == 2 && appxValuesInt.rows() == 3, "Matrix input has incorrect dimension");
    GISMO_ASSERT(appxPointsInt.cols() == appxValuesInt.cols(), "Matrix input has incorrect dimension");
    GISMO_ASSERT(appxNormalPoints.rows() == 2 && appxNormals.rows() == 3, "Matrix input has incorrect dimension");
    GISMO_ASSERT(appxNormalPoints.cols() == appxNormals.cols(), "Matrix input has incorrect dimension");


    int const patchDeg1 = kv1.degree();
    int const patchDeg2 = kv2.degree();
    int const n1 = kv1.size() - patchDeg1 - 1;
    int const n2 = kv2.size() - patchDeg2 - 1;

    gsBSplineBasis<T> bs1(kv1);
    gsBSplineBasis<T> bs2(kv2);
    gsMatrix<T>  Nu, Nv, dNu, dNv;
    gsMatrix<T> R;
    int npts;

    //Assemble Exact constraints for corners
    gsMatrix<T> ident1(n1, n1);
    gsMatrix<T> ident2(n2, n2);
    ident1.setIdentity();
    ident2.setIdentity();
    Nu  = bs1.evalFunc(exactPoints.row(0), ident1); // u-BSplines
    Nv  = bs2.evalFunc(exactPoints.row(1), ident2); // v-BSplines
    npts = exactPoints.cols();
    gsMatrix<T> Ecor(3*npts,3*n1*n2);
    gsMatrix<T> ecor(3*npts,1);
    Ecor.setZero();ecor.setZero();
    R = Nv.khatriRao(Nu); // R = M tensors N
    Ecor.block(0,0,npts,n1*n2) = R.transpose();
    Ecor.block(npts,n1*n2,npts,n1*n2) = R.transpose();
    Ecor.block(2*npts,2*n1*n2,npts,n1*n2) = R.transpose();
    ecor.block(0,0,npts,1) = (exactValues.row(0)).transpose();
    ecor.block(npts,0,npts,1) = (exactValues.row(1)).transpose();
    ecor.block(2*npts,0,npts,1) = (exactValues.row(2)).transpose();

    //Assemble Approximative constraints for inner points on edges
    Nu  = bs1.evalFunc(appxPointsEdges.row(0), ident1); // u-BSplines
    Nv  = bs2.evalFunc(appxPointsEdges.row(1), ident2); // v-BSplines
    npts = appxPointsEdges.cols();
    gsMatrix<T> AappEdge(3*npts,3*n1*n2);
    gsMatrix<T> bappEdge(3*npts,1);
    AappEdge.setZero();bappEdge.setZero();
    R = Nv.khatriRao(Nu); // R = M tensors N
    AappEdge.block(0,0,npts,n1*n2) = R.transpose();
    AappEdge.block(npts,n1*n2,npts,n1*n2) = R.transpose();
    AappEdge.block(2*npts,2*n1*n2,npts,n1*n2) = R.transpose();
    bappEdge.block(0,0,npts,1) = (appxValuesEdges.row(0)).transpose();
    bappEdge.block(npts,0,npts,1) = (appxValuesEdges.row(1)).transpose();
    bappEdge.block(2*npts,0,npts,1) = (appxValuesEdges.row(2)).transpose();

    //Assemble Approximate constraints for interior
    Nu  = bs1.evalFunc(appxPointsInt.row(0), ident1); // u-BSplines
    Nv  = bs2.evalFunc(appxPointsInt.row(1), ident2); // v-BSplines
    npts = appxPointsInt.cols();
    gsMatrix<T> AappInt(3*npts,3*n1*n2);
    gsMatrix<T> bappInt(3*npts,1);
    AappInt.setZero();bappInt.setZero();
    R = Nv.khatriRao(Nu); // R = M tensors N
    AappInt.block(0,0,npts,n1*n2) = R.transpose();
    AappInt.block(npts,n1*n2,npts,n1*n2) = R.transpose();
    AappInt.block(2*npts,2*n1*n2,npts,n1*n2) = R.transpose();
    bappInt.block(0,0,npts,1) = (appxValuesInt.row(0)).transpose();
    bappInt.block(npts,0,npts,1) = (appxValuesInt.row(1)).transpose();
    bappInt.block(2*npts,0,npts,1) = (appxValuesInt.row(2)).transpose();

    //Assemble Approximative constraints for normals
    Nu  = bs1.evalFunc(appxNormalPoints.row(0), ident1); // u-BSplines
    dNu = bs1.derivFunc(appxNormalPoints.row(0), ident1);
    Nv  = bs2.evalFunc(appxNormalPoints.row(1), ident2); // v-BSplines
    dNv = bs2.derivFunc(appxNormalPoints.row(1), ident2);
    gsMatrix<T> Nx,Ny,Nz; // x, y, z components of normals
    gsMatrix<T> trNormals = appxNormals.transpose();
    Nx = trNormals.col(0);
    Ny = trNormals.col(1);
    Nz = trNormals.col(2);
    if (force_normal==true) { Nx.setZero(); Ny.setZero();Nz.setOnes();} //TODO: do this automatically

    gsMatrix<T> dRdu = Nv.khatriRao(dNu); // R = M tensors N
    gsMatrix<T> dRdv = dNv.khatriRao(Nu); // R = M tensors N
    dRdu.transposeInPlace();
    dRdv.transposeInPlace();
    int nnor = Nx.rows(); // number of normals
    gsMatrix<T> Anor(2*nnor,3*n1*n2),bnor(2*nnor,1);
    Anor.setZero();
    bnor.setZero();
    Anor.block(0,0,nnor,n1*n2) = Nx.asDiagonal()*dRdu;  // sigma_u . normal=0, x part
    Anor.block(0,n1*n2,nnor,n1*n2) = Ny.asDiagonal()*dRdu;  // sigma_u . normal=0, y part
    Anor.block(0,2*n1*n2,nnor,n1*n2) = Nz.asDiagonal()*dRdu;  // sigma_u . normal=0, z part
    Anor.block(nnor,0,nnor,n1*n2) = Nx.asDiagonal()*dRdv;  // sigma_v . normal=0, x part
    Anor.block(nnor,n1*n2,nnor,n1*n2) = Ny.asDiagonal()*dRdv;  // sigma_v . normal=0, y part
    Anor.block(nnor,2*n1*n2,nnor,n1*n2) = Nz.asDiagonal()*dRdv;  // sigma_v . normal=0, z part

    // Quadratic forms which constitute the plate bending energy
    gsMatrix<T> * M = innerProduct(bs1, bs1);
    gsMatrix<T> * M1 = innerProduct1(bs1, bs1);
    gsMatrix<T> * M2 = innerProduct2(bs1, bs1);
    gsMatrix<T> *N, *N1, *N2;
    if (kv1==kv2) { N = M; N1 = M1; N2 = M2; } else
    {
        N  = innerProduct(bs2, bs2);
        N1 = innerProduct1(bs2, bs2);
        N2 = innerProduct2(bs2, bs2);
    };
    gsMatrix<T> Q1D = N->kron(*M2) + 2*N1->kron(*M1)+N2->kron(*M);
    gsMatrix<T> Q(3*n1*n2,3*n1*n2);
    Q.setZero();
    Q.block(0,0,n1*n2,n1*n2) = Q1D;
    Q.block(n1*n2,n1*n2,n1*n2,n1*n2) = Q1D;
    Q.block(2*n1*n2,2*n1*n2,n1*n2,n1*n2) = Q1D;

    // now solve the cps
    gsMatrix<T> coefA = wEdge*(AappEdge.transpose())*AappEdge + wInt*(AappInt.transpose())*AappInt + wNormal*(Anor.transpose())*Anor + wReg*Q;
    gsMatrix<T> coefb = (-2)*wEdge*(AappEdge.transpose())*bappEdge + (-2)*wInt*(AappInt.transpose())*bappInt + (-2)*wNormal*(Anor.transpose())*bnor;

    gsMatrix<T> cp = criticalPointOfQuadratic(
        coefA,
        coefb,
        Ecor,
        ecor);

    cp.resize( cp.rows() / 3, 3);

    typename gsTensorBSpline<2,T>::Ptr master(new gsTensorBSpline<2,T>( kv1, kv2, give(cp) ));

    delete M; delete M1; delete M2;
    if (kv1!=kv2) {delete N; delete N1; delete N2;}

    // check that the spline actually satisfies the exact constraints
    for(index_t idxConstr = 0; idxConstr < exactPoints.cols(); idxConstr++)
    {
        gsMatrix<T> pt = exactPoints.col(idxConstr);
        gsMatrix<T> expectVal = exactValues.col(idxConstr);
        gsMatrix<T> surfVal = master->eval(pt);
        GISMO_ASSERT((expectVal - surfVal).norm() < 0.0001, "Fit surface did not satisfy exact constraints");
    }

    return master;
}


} // namespace gismo
