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
#include <gsCore/gsFuncData.h>

#include <gsCore/gsGeometrySlice.h>

namespace gismo
{

template<class T>
typename gsGeometry<T>::uPtr
gsGeometry<T>::boundary(boxSide const& s) const
{
    gsMatrix<unsigned> ind = this->basis().boundary(s); // get indices of the boundary DOF
    gsMatrix<T> coeffs (ind.size(), geoDim()); // create matrix for boundary coefficients
    
    for (index_t i=0; i != ind.size(); i++ )
    {
        coeffs.row(i) = m_coefs.row( (ind)(i,0) );
    }

    gsBasis<T> *Bs = this->basis().boundaryBasis(s);  // Basis for boundary side s
    uPtr bgeo = Bs->makeGeometry( give(coeffs) );

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
std::vector<gsGeometry<T>* >  gsGeometry<T>::uniformSplit(index_t dir) const
{
    GISMO_NO_IMPLEMENTATION
}

template<class T>
gsGeometrySlice<T> gsGeometry<T>::getIsoParametricSlice(index_t dir_fixed, T par) const
{
    return gsGeometrySlice<T>(this,dir_fixed,par);
}

template<class T>
typename gsMatrix<T>::RowXpr
gsGeometry<T>::coefAtCorner(boxCorner const & c)
{
    return this->m_coefs.row(this->basis().functionAtCorner(c));
}

template<class T>
typename gsMatrix<T>::ConstRowXpr
gsGeometry<T>::coefAtCorner(boxCorner const & c) const
{
    return this->m_coefs.row(this->basis().functionAtCorner(c));
}

template<class T>
void gsGeometry<T>::invertPoints(const gsMatrix<T> & points, 
                                 gsMatrix<T> & result, 
                                 const T accuracy)
{
    result.resize(parDim(), points.cols() );
    gsVector<T> arg;
    for ( index_t i = 0; i!= points.cols(); ++i)
    {         
        arg = parameterCenter();
        // int iter = 
        this->newtonRaphson(points.col(i), arg, true, accuracy, 100);
        result.col(i) = arg;
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
void gsGeometry<T>::degreeElevate(int const i, int const dir)
{
    typename gsBasis<T>::uPtr b = m_basis->clone();
 
    if ( dir == -1 )
        b->degreeElevate(i);
    else if (dir < parDim() )
        b->degreeElevate(i, dir);
    else
        GISMO_ERROR("Invalid direction "<< dir <<" to elevate.");

    gsMatrix<T> iVals, iPts = b->anchors();
    this->eval_into(iPts, iVals);
    typename gsGeometry<T>::uPtr g = b->interpolateData(iVals, iPts);

    std::swap(m_basis, g->m_basis);
    g->coefs().swap(this->coefs());
}

template<class T>
void gsGeometry<T>::degreeReduce(int const i, int const dir)
{
    typename gsBasis<T>::uPtr b = m_basis->clone();
 
    if ( dir == -1 )
        b->degreeReduce(i);
    else if (dir < parDim() )
        b->component(dir).degreeReduce(i);        
    else
        GISMO_ERROR("Invalid direction "<< dir <<" to degree-reduce.");

    gsMatrix<T> iVals, iPts = b->anchors();
    this->eval_into(iPts, iVals);
    typename gsGeometry<T>::uPtr g = b->interpolateData(iVals, iPts);

    std::swap(m_basis, g->m_basis);
    g->coefs().swap(this->coefs());
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



template <typename T>
void extractRows( const gsMatrix<T> &in, typename gsMatrix<unsigned>::constColumn actives, gsMatrix<T> &out)
{
    out.resize(actives.rows(), in.cols());
    for (index_t r=0; r<actives.rows();++r)
        out.row(r)=in.row(actives(r,0));
}

template <class T>
void
gsGeometry<T>::compute(const gsMatrix<T> & in, gsFuncData<T> & out) const
{  

    const unsigned flags = out.flags | NEED_ACTIVE;
    const index_t  numPt = in.cols();
    const index_t  numCo = m_coefs.cols();

    gsFuncData<T> tmp(flags);
    this->basis().compute(in, tmp);
    
    out.values.resize(out.maxDeriv()+1);
    out.dim.first  = tmp.dim.first;
    out.dim.second = numCo;
    if ( flags & SAME_ELEMENT )
    {
        gsMatrix<T> coefM;
        extractRows(m_coefs,tmp.active(0),coefM);

        if (flags & NEED_VALUE)
            out.values[0]=coefM.transpose()*tmp.values[0];
        if (flags & NEED_DERIV)
        {
            const index_t derS = tmp.derivSize();
            out.values[1].resize(derS*numCo,numPt);
            for (index_t p=0; p< numPt; ++p)
                out.values[1].reshapeCol(p, derS, numCo) = tmp.deriv(p)*coefM;
        }
        if (flags & NEED_DERIV2)
        {
            const index_t derS = tmp.deriv2Size();
            out.values[2].resize(derS*numCo,numPt);
            for (index_t p=0; p< numPt; ++p)
                out.values[2].reshapeCol(p, derS, numCo) = tmp.deriv2(p)*coefM;
        }
    } else
    {
        gsMatrix<T> coefM;
        const index_t derS = tmp.derivSize();
        const index_t der2S = tmp.deriv2Size();

        if (flags & NEED_VALUE)  out.values[0].resize(numCo,numPt);
        if (flags & NEED_DERIV)  out.values[1].resize(numCo*derS,numPt);
        if (flags & NEED_DERIV2) out.values[2].resize(numCo*der2S,numPt);

        for (index_t p=0; p<numPt;++p)
        {
            extractRows(m_coefs,tmp.active(p),coefM);
            if (flags & NEED_VALUE)
                out.values[0].reshapeCol(p,1,numCo) = tmp.eval(p)*coefM;
            if (flags & NEED_DERIV)
                out.values[1].reshapeCol(p, derS, numCo) = tmp.deriv(p)*coefM;
            if (flags & NEED_DERIV2)
                out.values[2].reshapeCol(p, der2S, numCo) = tmp.deriv2(p)*coefM;
        }
    }
}


} // namespace gismo
