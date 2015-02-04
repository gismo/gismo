/** @file gsGeometry.hpp

    @brief Provides implementation of Geometry default operatiions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsUtils/gsMesh/gsMesh.h>
//#include <gsCore/gsBoundary.h>

namespace gismo
{

template<class T>
gsGeometry<T> * 
gsGeometry<T>::boundary(boxSide const& s) const
{
    gsMatrix<unsigned> *ind = this->basis().boundary(s); // get indices of the boundary DOF
    gsMatrix<T> coeffs (ind->size(), geoDim()); // create matrix for boundary coefficients
    
    for (index_t i=0; i != ind->size(); i++ )
    {
        coeffs.row(i) = m_coefs.row( (*ind)(i,0) );
    }
    delete ind;
    gsBasis<T> *Bs = this->basis().boundaryBasis(s);  // Basis for boundary side s
    gsGeometry *bgeo = Bs->makeGeometry( give(coeffs) );
    delete Bs;
    return bgeo;
}

template<class T>
void gsGeometry<T>::evaluateMesh(gsMesh<T>& mesh) const
{
    const int pDim = parDim();
    const int gDim = geoDim();

    gsMatrix<T> tmp;

    // For all vertices of the mesh, push forward the value by the
    // geometry mapping
    for ( int i = 0; i!= mesh.numVertices; ++i)
    {
        eval_into( mesh.vertex[i]->coords.topRows(pDim), tmp );
        mesh.vertex[i]->coords.topRows( gDim ) = tmp;
    }
}


template<class T>
void gsGeometry<T>::merge( gsGeometry * other )
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsGeometryEvaluator<typename gsGeometry<T>::Scalar_t> * 
gsGeometry<T>:: evaluator(unsigned flags) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsGeometry<T>::toMesh(gsMesh<T> & msh, int npoints) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsGeometry<T>::outerNormal_into(const gsMatrix<T>& u, gsMatrix<T> & result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
std::vector<gsGeometry<T> *> gsGeometry<T>:: boundary() const
{
    // TO DO: get boundary curves, using basis().boundary();
    GISMO_NO_IMPLEMENTATION
}


template<class T>
void gsGeometry<T>::degreeElevate(int const i) 
{
    gsBasis<T> * b = m_basis->clone();
    b->degreeElevate(i);
    
    gsMatrix<T> iVals, iPts = b->anchors();
    this->eval_into(iPts, iVals);
    gsGeometry<T> * g = b->interpolate(iVals, iPts);

    std::swap(m_basis, g->m_basis);
    g->coefs().swap(this->coefs());

    delete g;
    delete b;
}

template<class T>
void gsGeometry<T>::degreeElevate(int const dir, int const i)
{
    gsBasis<T> * b = m_basis->clone();
    b->component(dir).degreeElevate(i);
    
    gsMatrix<T> iVals, iPts = b->anchors();
    this->eval_into(iPts, iVals);
    gsGeometry<T> * g = b->interpolate(iVals, iPts);

    std::swap(m_basis, g->m_basis);
    g->coefs().swap(this->coefs());

    delete g;
    delete b;
}

template<class T>
typename gsMatrix<T>::uPtr
gsGeometry<T>::hessian(const gsMatrix<T>& u, unsigned coord) const
{  
    static const unsigned d = this->m_basis->dim();

    gsMatrix<T> B, *DD = new gsMatrix<T>(d,d);
    gsMatrix<T> tmp(d,d);
    gsMatrix<unsigned> ind;

    // coefficient matrix row k = coef. of basis function k
    const gsMatrix<T>& C = this->m_coefs; 
    // col j = nonzero second derivatives at column point u(..,j)
    m_basis->deriv2_into(u, B) ; 
    // col j = indices of active functions at column point u(..,j)
    m_basis->active_into(u, ind);  
  
    DD->setZero();
    unsigned j=0;// just one column
    //for ( unsigned j=0; j< u.cols(); j++ ) // for all points (columns of u)
    for ( index_t i=0; i< ind.rows() ; i++ ) // for all non-zero basis functions)
    {
        unsigned m=i*d;
        unsigned r= ind.rows()*d + i*d*(d-1)/2;
        //construct the Hessian of basis function ind(i,0)
        for (unsigned k=0; k<d; ++k ) // for all rows
        {
            tmp(k,k) = B(m+k,j);
            for (unsigned l=k+1; l<d; ++l ) // for all cols
                tmp(k,l) = tmp(l,k) = B(r++,0);
        }
        *DD += C(ind(i,j), coord) * tmp;
    }
  
    return typename gsMatrix<T>::uPtr(DD); 
}

}; // namespace gismo
