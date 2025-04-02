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
#include <gsIO/gsFileData.h>
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

    typename gsBasis<T>::domainIter domIt = B1.domain()->beginAll();
    typename gsBasis<T>::domainIter domItEnd = B1.domain()->endAll();

    for (; domIt<domItEnd; ++domIt )
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), ngrid, wgrid );

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

    typename gsBasis<T>::domainIter domIt = B1.domain()->beginAll();
    typename gsBasis<T>::domainIter domItEnd = B1.domain()->endAll();
    for (; domIt<domItEnd; ++domIt )
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), ngrid, wgrid );

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

    typename gsBasis<T>::domainIter domIt = B1.domain()->beginAll();
    typename gsBasis<T>::domainIter domItEnd = B1.domain()->endAll();
    for (; domIt<domItEnd; ++domIt )
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), ngrid, wgrid );

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
    T cag = ( normal.dot( vec1.cross( vec2 ) ) >= 0 ) ? ag : (T)(2*EIGEN_PI)-ag;
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
    gsMatrix<T> Hess = w_reg*(T)(2)*(*Q) + w_app*(T)(2)*(NuApp.transpose())* NuApp;
    //--- row 0
    Ass.block(0,0,ntcp,ntcp) = Hess;
    Ass.block(0,2*ntcp,ntcp,nip) = Nu.transpose();
    Ass.block(0,2*ntcp+2*nip,ntcp,nn) = AdN.transpose();
    bss.block(0,0,ntcp,1) = w_app*(T)(2)*NuApp.transpose()*X0.transpose();
    //--- row 1
    Ass.block(ntcp,ntcp,ntcp,ntcp) = Hess;
    Ass.block(ntcp,2*ntcp+nip,ntcp,nip) = Nu.transpose();
    Ass.block(ntcp,2*ntcp+2*nip,ntcp,nn) = BdN.transpose();
    bss.block(ntcp,0,ntcp,1) = w_app*(T)(2)*NuApp.transpose()*Y0.transpose();
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
    gsMatrix<T> coefb = (T)(-2)*wEdge*(AappEdge.transpose())*bappEdge + (T)(-2)*wInt*(AappInt.transpose())*bappInt + (T)(-2)*wNormal*(Anor.transpose())*bnor;

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

    /**
       Constructs a gsSparseMatrix<T> with with block \a block repeated three times along the diagonal and
       saves it to \a result.
       TODO: Make more general.
     */
    template<class T>
    void threeOnDiag(const gsSparseMatrix<T>& block,
                    gsSparseMatrix<T>& result) //const
    {

      index_t brows = block.rows();
      index_t bcols = block.cols();

      result = gsSparseMatrix<T>(3*brows,3*bcols);
      result.reservePerColumn(block.nonZeros() / bcols);

      for(index_t j = 0; j != bcols; j++)
      {
        for (auto it = block.begin(j); it; ++it)
        {
            result.insertTo(it.row(), j, it.value());
            result.insertTo(brows+it.row(), bcols+j, it.value());
            result.insertTo(2*brows+it.row(), 2*bcols+j, it.value());
        }
      }
      result.makeCompressed();
    }

    /// Function to scale the input points \a xyz in [0,1]^D and saves it to \a points.
    template<class T>
    void scalePoints(const gsMatrix<T> & xyz,
                           gsMatrix<T> & points)
    {
        T p_min = xyz.minCoeff(),
        p_max = xyz.maxCoeff();
        T den = p_max - p_min;

        points.resize(xyz.rows(), xyz.cols());
        points = (1/den)*(xyz - p_min * gsMatrix<T>::Ones(xyz.rows(), xyz.cols()));
    }


    /// Scale the interval [tMin, tMax] to [0, 1].
    template<class T>
    T scaleTo01(T tMin, T t, T tMax)
    {
	    return (t - tMin) / (tMax - tMin);
    }

    /// Scale the matrix \a mT entries to [0, 1].
    template <class T>
    void scaleTo01(real_t tMin, gsMatrix<T>& mT, real_t tMax)
    {
	    for(index_t i=0; i<mT.rows(); i++)
		    for(index_t j=0; j<mT.cols(); j++)
			    mT(i, j) = scaleTo01(tMin, mT(i, j), tMax);
    }

    /// Scale the matrix \a xyz entries to [0, 1]^D.
    template <class T>
    void scaleTo01(gsMatrix<T>& xyz, bool verbose)
    {
	    T xyzMin = xyz.minCoeff();
	    T xyzMax = xyz.maxCoeff();
	    scaleTo01(xyzMin, xyz, xyzMax);

	    if(verbose)
		    gsInfo << std::setprecision(15)
			       << "Scaling TO [0, 1]^3 as follows.\n"
			       << "min:   " << xyzMin << "\n"
			       << "max:   " << xyzMax << "\n"
			       << "scale: " << T(1.0) / (xyzMax - xyzMin) << std::endl;
    }


    /// Scale the inteval [0,1] to [tMin, tMax].
    template <class T>
    T scaleFrom01(T tMin, T t, T tMax)
    {
	    return (tMax - tMin) * t + tMin;
    }


    /// Scale the matrix \a mT entries from [0, 1] to [tMin, tMax].
    template <class T>
    void scaleFrom01(T tMin, gsMatrix<T>& mT, T tMax, bool verbose)
    {
	    for(index_t i=0; i<mT.rows(); i++)
		    for(index_t j=0; j<mT.cols(); j++)
			    mT(i, j) = scaleFrom01(tMin, mT(i, j), tMax);

	    if(verbose)
		    gsInfo << "Scaling FROM [0, 1]^3.\n"
			       << "inverted scale: " << T(1.0) / (tMax - tMin) << std::scientific << std::endl;
    }


    /// Scale the geometry \a geo from [0, 1]^D to [tMin, tMax]^D.
    template <class T>
    void scaleFrom01(T tMin, gsGeometry<T>& geo, T tMax, bool verbose)
    {
	    gsMatrix<T> coefs = geo.coefs();
	    scaleFrom01(tMin, coefs, tMax, verbose);
	    geo.coefs() = coefs;
    }


    /// Scale the geometry contained in file \a fin from [0, 1]^D to [tMin, tMax]^D and save it to \a fout.
    template <class T>
    void scaleGeo(const std::string& fin,
			      const std::string& fout,
			      T tMin,
			      T tMax,
			      bool verbose)
    {
	    gsFileData<T> fdIn(fin);
	    typename gsGeometry<T>::uPtr geo = fdIn.template getFirst< gsGeometry<T> >();
	    scaleFrom01(tMin, *geo.get(), tMax, verbose);
	    gsFileData<T> fdOut;
	    fdOut << *geo.get();
	    fdOut.dump(fout);
    }
    
    
    /// Scale the points contained in file \a fin from [0, 1]^D to [tMin, tMax]^D and save it to \a fout.
    template <class T>
    void scalePts(const std::string& fin,
                  const std::string& fout,
			      index_t uvIdIn,
			      index_t uvIdOut,
			      index_t xyzIdIn,
			      index_t xyzIdOut,
			      T tMin,
			      T tMax,
			      bool verbose)
    {
        // Reading the inputs.
        gsFileData<T> fdIn(fin);
        gsMatrix<T> xyz, uv;
        fdIn.template getId<gsMatrix<T>>(uvIdIn,  uv);
        fdIn.template getId<gsMatrix<T>>(xyzIdIn, xyz);

        // Scaling the data.
        if(tMax < tMin) // defaults, i.e., scaling down
            scaleTo01(xyz, verbose);
        else
            scaleFrom01(tMin, xyz, tMax, verbose);

        // Writing the outputs
        gsFileData<T> fdOut;

        if(uvIdOut < xyzIdOut)
        {
            fdOut << uv;
            fdOut << xyz;
        }
        else
        {
            fdOut << xyz;
            fdOut << uv;
        }
        fdOut.dump(fout);
    }

    
    /** \brief sortPointCloud: sorts the point cloud into interior and boundary points.
    * parameters and points ordered by : interior (parameters/points) and
      boundary (parameters/points) ordered anticlockwise south-east-north-west edges,
      plus the 4 corner domains stored in a vector [c1, c2, c3, c4].
    * @param parameters : matrix of parameters
    * @param points : matrix of points 
    * @param corners : vector of corner domains indeces [c1, c2, c3, c4]
    */
    template<class T>
    void sortPointCloud(gsMatrix<T> & parameters,
                        gsMatrix<T> & points,
                        std::vector<index_t> & corners)
    {
        // The following matrices and vectors store the parameters and points values and indeces.
        // There is no need to store these information, we could also use only one matrix and 1 std::vector and overwirte them each time.
        gsMatrix<T> uv_interiors, uv_south, uv_east, uv_north, uv_west;
        gsMatrix<T> p_interiors, p_south, p_east, p_north, p_west;
        std::vector<index_t> interiors, b_west, b_east, b_south, b_north;

        // Determine the parameter domain by mi/max of parameter values
        T u_min = parameters.row(0).minCoeff(),
        u_max = parameters.row(0).maxCoeff(),
        v_min = parameters.row(1).minCoeff(),
        v_max = parameters.row(1).maxCoeff();

        gsVector<T> curr_point(2,1);
        for(index_t i=0; i < parameters.cols(); i++)
        {
            curr_point = parameters.col(i);
            if( (u_min < curr_point(0)) && (curr_point(0) < u_max) && (v_min < curr_point(1)) && (curr_point(1) < v_max) )
                interiors.push_back(i);
            else // not interior point
            {
                if( (math::abs(curr_point(0) - u_min) < 1e-15) && (curr_point(1) > v_min) )
                    b_west.push_back(i);//west edge
                else if( (math::abs(curr_point(0) - u_max) < 1e-15) && curr_point(1) < v_max)
                    b_east.push_back(i);// east edge
                else if( (math::abs(curr_point(1) - v_min) < 1e-15) && (curr_point(0) < u_max) )
                    b_south.push_back(i);// south edge
                else
                    b_north.push_back(i);// north edge
            }
        }
    
        corners.push_back(interiors.size()); // c1
        corners.push_back(interiors.size() + b_south.size()); // c2
        corners.push_back(interiors.size() + b_south.size() + b_east.size()); // c3
        corners.push_back(interiors.size() + b_south.size() + b_east.size() + b_north.size()); // c4

        uv_interiors.resize(2, interiors.size());
        p_interiors.resize(3, interiors.size());
        for( size_t i = 0; i < interiors.size(); i++ )
        {
            uv_interiors.col(i) = parameters.col(interiors[i]);
            p_interiors.col(i) = points.col(interiors[i]);
        }

        uv_west.resize(2, b_west.size());
        gsMatrix<T> tmp_west(3, b_west.size());
        for( size_t i = 0; i < b_west.size(); i++ )
        {
            uv_west.col(i) = parameters.col(b_west[i]);
            tmp_west.col(i) = points.col(b_west[i]);
        }

        uv_east.resize(2, b_east.size());
        gsMatrix<T> tmp_east(3, b_east.size());
        for( size_t i = 0; i < b_east.size(); i++ )
        {
            uv_east.col(i) = parameters.col(b_east[i]);
            tmp_east.col(i) = points.col(b_east[i]);
        }

        uv_south.resize(2, b_south.size());
        gsMatrix<T> tmp_south(3, b_south.size());
        for( size_t i = 0; i < b_south.size(); i++ )
        {
            uv_south.col(i) = parameters.col(b_south[i]);
            tmp_south.col(i) = points.col(b_south[i]);
        }

        uv_north.resize(2, b_north.size());
        gsMatrix<T> tmp_north(3, b_north.size());
        for( size_t i = 0; i < b_north.size(); i++ )
        {
            uv_north.col(i) = parameters.col(b_north[i]);
            tmp_north.col(i) = points.col(b_north[i]);
        }

        uv_south.transposeInPlace();
        uv_east.transposeInPlace();
        uv_north.transposeInPlace();
        uv_west.transposeInPlace();


        std::vector<index_t> tmp = uv_south.idxByColumn(0);
        p_south.resize(tmp_south.rows(), tmp_south.cols());
        for(size_t i = 0; i<tmp.size(); i++)
        {
            p_south.col(i) = tmp_south.col(tmp[i]);
        }
        uv_south.transposeInPlace();


        tmp.clear();
        tmp = uv_east.idxByColumn(1);
        p_east.resize(tmp_east.rows(), tmp_east.cols());
        for(size_t i = 0; i<tmp.size(); i++)
        {
            p_east.col(i) = tmp_east.col(tmp[i]);
        }
        uv_east.transposeInPlace();


        tmp.clear();
        tmp = uv_north.idxByColumn(0);
        std::reverse(tmp.begin(),tmp.end());
        gsVector<T> tcol = uv_north.col(0).reverse();
        uv_north.col(0) = tcol;
        tcol = uv_north.col(1).reverse();
        uv_north.col(1) = tcol;
        
        p_north.resize(tmp_north.rows(), tmp_north.cols());
        for(size_t i = 0; i<tmp.size(); i++)
        {
            p_north.col(i) = tmp_north.col(tmp[i]);
        }
        uv_north.transposeInPlace();

        tmp.clear();
        tmp = uv_west.idxByColumn(1);
        tcol = uv_west.col(0).reverse();
        uv_west.col(0) = tcol;
        tcol = uv_west.col(1).reverse();
        uv_west.col(1) = tcol;
        std::reverse(tmp.begin(),tmp.end());

        p_west.resize(tmp_west.rows(), tmp_west.cols());
        for(size_t i = 0; i<tmp.size(); i++)
        {
            p_west.col(i) = tmp_west.col(tmp[i]);
        }
        uv_west.transposeInPlace();


        // reordering of the input point cloud (parameters and points)
        parameters.resize(uv_interiors.rows(), points.cols());
        parameters << uv_interiors.row(0), uv_south.row(0), uv_east.row(0), uv_north.row(0), uv_west.row(0),
                    uv_interiors.row(1), uv_south.row(1), uv_east.row(1), uv_north.row(1), uv_west.row(1);

        points.resize(p_interiors.rows(), parameters.cols());
        points << p_interiors.row(0), p_south.row(0), p_east.row(0), p_north.row(0), p_west.row(0),
                p_interiors.row(1), p_south.row(1), p_east.row(1), p_north.row(1), p_west.row(1),
                p_interiors.row(2), p_south.row(2), p_east.row(2), p_north.row(2), p_west.row(2);

    } // end sortPointCloud


    
    /** \brief sampleGridGeometry: samples a grid point cloud from a given geometry
    * @param mp : input multi-patch
    * @param numPatch : patch number
    * @param numSamples : number of samples in each direction
    * @param params : output parameters
    * @param points : output points 
    */
    template<class T>
    void sampleGridGeometry(const gsMultiPatch<T> & mp,
                            const index_t & numPatch,
                            const index_t & numSamples,
                            gsMatrix<T> & params,
                            gsMatrix<T> & points)
    {
        GISMO_ASSERT( numPatch <= mp.nPatches()-1 , "Patch number not found, quitting.");
        
        const gsGeometry<T> & geometry = mp.patch(numPatch);
        
        gsVector<unsigned> numPtsVec(2);
        numPtsVec<<numSamples,numSamples;
        gsVector<T> a = geometry.support().col(0);
        gsVector<T> b = geometry.support().col(1);
        
        params = gsPointGrid(a,b, numPtsVec);
        geometry.eval_into(params, points);
    }


    /** \brief sampleScatteredGeometry: samples a scattered point cloud from a given geometry
     * @param mp : input multi-patch
     * @param numPatch : patch number
     * @param numSamples : number of interior samples
     * @param numBdr : number of boundary samples
     * @param params : output parameters
     * @param points : output points 
     */
    template<class T>
    void sampleScatteredGeometry(const gsMultiPatch<T> & mp,
                                const index_t & numPatch,
                                const index_t & numSamples,
                                index_t & numBdr,
                                gsMatrix<T> & params,
                                gsMatrix<T> & points)
    {
        GISMO_ASSERT( numPatch <= mp.nPatches()-1 , "Patch number not found, quitting.");
        
        const gsGeometry<T> & geometry = mp.patch(numPatch);
        

        // Sample the interior parameters
        gsVector<unsigned> numPtsVec(2);
        numPtsVec<<numSamples,numSamples;
        gsVector<T> a = geometry.support().col(0);
        gsVector<T> b = geometry.support().col(1);
        
        T urange= b(0)-a(0); // umax - umin
        
        gsMatrix<T> mu = gsMatrix<T>::Random(1,numSamples); // 3x3 Matrix filled with random numbers between (-1,1)
        mu = (mu + gsMatrix<T>::Constant(1,numSamples,1))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        mu = (mu + gsMatrix<T>::Constant(1,numSamples,a(0))); //set LO as the lower bound (offset)

        T vrange= b(1)-a(1); // vmax - vmin
        gsMatrix<T> mv= gsMatrix<T>::Random(1,numSamples); // 3x3 Matrix filled with random numbers between (-1,1)
        mv = (mv + gsMatrix<T>::Constant(1,numSamples,1))*vrange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        mv = (mv + gsMatrix<T>::Constant(1,numSamples,a(1))); //set LO as the lower bound (offset)
        
        gsMatrix<T> uv_interiors(2, numSamples);
        uv_interiors << mu, mv; //interior parameters

        // Sample the boundary and the corner parameters
        if (numBdr < 2)
            numBdr = cast<T,index_t>(math::ceil(math::sqrt(numSamples))) + 2; // number of boundary points
        
        gsMatrix<T> uv_boundary(2, numBdr*4-4);
        gsMatrix<T> b_0(1, numBdr-1);
        gsMatrix<T> b_1(1, numBdr-1);
        gsMatrix<T> b_2(1, numBdr-1);
        gsMatrix<T> b_3(1, numBdr-1);

        for (index_t pace=0; pace < 4; pace++)
        {
            if (pace == 0){
                gsMatrix<T> mu = gsMatrix<T>::Random(numBdr-2, 1); // 1xnumBdr-1 Matrix filled with random numbers between (-1,1)
                mu = (mu + gsMatrix<T>::Constant(numBdr-2, 1, 1))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
                mu = (mu + gsMatrix<T>::Constant(numBdr-2,1,a(0))); //set LO as the lower bound (offset)
                mu.sortByColumn(0);
                mu = mu.reshape(1, numBdr-2);
                b_0 << a(0), mu;
            }
            if (pace == 1){
                gsMatrix<T> mu = gsMatrix<T>::Random(numBdr-2, 1); // 1xnumBdr-1 Matrix filled with random numbers between (-1,1)
                mu = (mu + gsMatrix<T>::Constant(numBdr-2, 1, 1))*vrange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
                mu = (mu + gsMatrix<T>::Constant(numBdr-2,1, a(1))); //set LO as the lower bound (offset)
                mu.sortByColumn(0);
                mu = mu.reshape(1, numBdr-2);
                b_1 << a(1), mu;
            }
            if (pace == 2){
                gsMatrix<T> mu = gsMatrix<T>::Random(numBdr-2, 1); // 1xnumBdr-1 Matrix filled with random numbers between (-1,1)
                mu = (mu + gsMatrix<T>::Constant(numBdr-2, 1, 1))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
                mu = (mu + gsMatrix<T>::Constant(numBdr-2,1,a(0))); //set LO as the lower bound (offset)
                mu.sortByColumn(0);
                mu = mu.reshape(1, numBdr-2);
                b_2 << b(0), mu.reverse();
            }
            if (pace == 3){
                gsMatrix<T> mu = gsMatrix<T>::Random(numBdr-2, 1); // 1xnumBdr-1 Matrix filled with random numbers between (-1,1)
                mu = (mu + gsMatrix<T>::Constant(numBdr-2, 1, 1))*vrange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
                mu = (mu + gsMatrix<T>::Constant(numBdr-2,1,a(1))); //set LO as the lower bound (offset)
                mu.sortByColumn(0);
                mu = mu.reshape(1, numBdr-2);
                b_3 << b(1), mu.reverse();
            }
        }

        gsMatrix<T> u_zeros = gsMatrix<T>::Constant(1, numBdr-1, a(1));
        gsMatrix<T> v_zeros = gsMatrix<T>::Constant(1, numBdr-1, a(0));
        gsMatrix<T> u_ones = gsMatrix<T>::Constant(1, numBdr-1, b(1));
        gsMatrix<T> v_ones = gsMatrix<T>::Constant(1, numBdr-1, b(0));
        gsMatrix<T> zeros = gsMatrix<T>::Zero(1, numBdr-1);
        gsMatrix<T> ones  = gsMatrix<T>::Ones(1, numBdr-1);

    
        uv_boundary << b_0,     v_ones, b_2,    v_zeros, u_zeros, b_1,    u_ones, b_3;

        params.resize(2, numSamples + numBdr*4-4);
        params << uv_interiors, uv_boundary;

        geometry.eval_into(params, points);

    }



} // namespace gismo
