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
gsGeometry<T>::boundary(boundary::side const& s) const
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


}; // namespace gismo
