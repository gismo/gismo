/** @file gsComposedBasis.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsComposedGeometry.h>

namespace gismo
{

template <class T>
class gsComposedBasis : public gsBasis<T>
{

    typedef memory::shared_ptr< gsComposedBasis > Ptr;
    typedef memory::unique_ptr< gsComposedBasis > uPtr;

    // /// Geometry Type
    typedef gsComposedGeometry<T> GeometryType;

    GISMO_CLONE_FUNCTION(gsComposedBasis)
//
    GISMO_MAKE_GEOMETRY_NEW

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

public:
    gsComposedBasis(const gsFunction<T> & composition, const gsBasis<T> & basis)
    :
    m_composition(&composition),
    m_basis(&basis)
    {
        GISMO_ENSURE(m_basis->domainDim()==m_composition->targetDim(),
            "Domain dimension of the basis "<<
            " should be equal to the target dimension of the composition "<<
            ", but basis.domainDim() = "<<basis.domainDim()<<
            " and composition.targetDim() = )"<<composition.targetDim());
    }

    short_t domainDim() const override { return m_composition->domainDim(); }
    short_t targetDim() const override { return m_basis->targetDim(); }

    gsMatrix<T> support() const override
    {
        gsMatrix<T> supp = m_basis->support();
        gsGridIterator<T,CUBE> pt(supp,math::pow(2,this->domainDim()));
        supp = pt.toMatrix();
        gsMatrix<T> result = supp;

        m_composition->invertPoints(supp,result,1e-10,true);

        supp.conservativeResize(this->domainDim(),2);
        for (size_t d=0; d!=this->domainDim(); d++)
            supp.row(d)<<result.row(d).array().minCoeff(),result.row(d).array().maxCoeff();

        return supp;
    } // This should be the inverse map

    gsMatrix<T> support(const index_t & i) const override
    {
        gsMatrix<T> supp = m_basis->support(i);
        gsGridIterator<T,CUBE> pt(supp,math::pow(2,this->domainDim()));
        supp = pt.toMatrix();
        gsMatrix<T> result = supp;

        m_composition->invertPoints(supp,result,1e-10,true);

        supp.conservativeResize(this->domainDim(),2);
        for (size_t d=0; d!=this->domainDim(); d++)
            supp.row(d)<<result.row(d).array().minCoeff(),result.row(d).array().maxCoeff();

        return supp;
    } // This should be the inverse map

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const override
    {
        gsMatrix<T> coords = m_composition->eval(u);
        this->_applyBounds(coords);
        m_basis->active_into(coords,result);
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {

        gsMatrix<T> coords = m_composition->eval(u);
        this->_applyBounds(coords);
        m_basis->eval_into(coords,result);
    }

    void evalSingle_into(index_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> coords = m_composition->eval(u);
        this->_applyBounds(coords);
        m_basis->evalSingle_into(i,coords,result);
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        const size_t DIM = m_basis->domainDim();
        index_t domainDim, targetDim;
        gsMatrix<T> coord, deriv, tmp, compderiv;

        m_composition->deriv_into(u,compderiv);
        domainDim = m_composition->domainDim();
        targetDim = m_composition->targetDim();

        m_composition->eval_into(u,coord);
        m_basis->deriv_into(coord,deriv);
        const size_t numAct = deriv.rows() / DIM;

        result.resize(numAct*domainDim*m_basis->targetDim(),u.cols());
        for (index_t k = 0; k!=u.cols(); k++)
        {
            gsAsMatrix<T,Dynamic,Dynamic> resultMat = compderiv.reshapeCol(k,domainDim,targetDim);
            gsAsMatrix<T,Dynamic,Dynamic> derivMat = deriv.reshapeCol(k,m_basis->domainDim(),m_basis->targetDim());
            for (index_t act = 0; act!=numAct; act++)
            {
                gsDebugVar(resultMat*deriv.block(act*DIM,k,domainDim*m_basis->targetDim(),1).reshaped(domainDim,m_basis->targetDim()));
                result.block(act*DIM,k,domainDim*m_basis->targetDim(),1).reshaped(domainDim,m_basis->targetDim()) = resultMat*deriv.block(act*DIM,k,domainDim*m_basis->targetDim(),1).reshaped(domainDim,m_basis->targetDim());

            }
        }

        // m_basis->deriv_into(coord,deriv);

        // tmp.resize(m_basis->targetDim()*domainDim,u.cols());

        //     gsAsMatrix<T,Dynamic,Dynamic> resultMat = result.reshapeCol(k,domainDim,targetDim);
        //     gsAsMatrix<T,Dynamic,Dynamic> derivMat = deriv.reshapeCol(k,m_basis->domainDim(),m_basis->targetDim());
        //     // The product has size:
        //     // (domainDim x targetDim) x (m_basis->domainDim(),m_basis->targetDim())
        //     //  =
        //     // (domainDim x m_basis->targetDim())
        //     gsAsMatrix<T,Dynamic,Dynamic> tmpMat = tmp.reshapeCol(k,domainDim,m_basis->targetDim());
        //     tmpMat = resultMat*derivMat;

        // }
        // result = tmp;
    }

    void derivSingle_into(index_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        index_t domainDim, targetDim;
        gsMatrix<T> coord, deriv, tmp, tmpderiv;

        m_composition->deriv_into(u,result);
        domainDim = m_composition->domainDim();
        targetDim = m_composition->targetDim();

        m_composition->eval_into(u,coord);
        m_basis->derivSingle_into(i,coord,deriv);

        tmp.resize(m_basis->targetDim()*domainDim,u.cols());
        for (index_t k = 0; k!=u.cols(); k++)
        {
            gsAsMatrix<T,Dynamic,Dynamic> resultMat = result.reshapeCol(k,domainDim,targetDim);
            gsAsMatrix<T,Dynamic,Dynamic> derivMat = deriv.reshapeCol(k,m_basis->domainDim(),m_basis->targetDim());
            // The product has size:
            // (domainDim x targetDim) x (m_basis->domainDim(),m_basis->targetDim())
            //  =
            // (domainDim x m_basis->targetDim())
            gsAsMatrix<T,Dynamic,Dynamic> tmpMat = tmp.reshapeCol(k,domainDim,m_basis->targetDim());
            tmpMat = resultMat*derivMat;

        }
        result = tmp;
    }

    // memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs) const override
    // {
    //     return memory::unique_ptr<gsGeometry<T> >(new gsGeometry<T>(*this, give(coefs)));
    //     // GISMO_NO_IMPLEMENTATION;
    // }

    ////// Pass throughs of basis

    /// See \ref gsBasis for documentation
    domainIter makeDomainIterator() const override { return m_basis->makeDomainIterator(); }

    /// See \ref gsBasis for documentation
    std::string detail() const override { return m_basis->detail(); };

    /// See \ref gsBasis for documentation
    size_t numElements() const override { return m_basis->numElements(); }

    /// See \ref gsBasis for documentation
    size_t numElements(boxSide const & s) const override { return m_basis->numElements(s); }

    /// See \ref gsBasis for documentation
    index_t size() const override {return m_basis->size(); }

    /// See \ref gsBasis for documentation
    void anchors_into(gsMatrix<T> & result) const override { m_basis->anchors_into(result); }


    /// See \ref gsBasis for documentation
    void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const override { m_basis->connectivity(nodes,mesh); }


    /// See \ref gsBasis for documentation


    /// See \ref gsBasis for documentation
    //

    //// NEEDS TO EVALUATE THE INVERSE MAP, NOT m_composition!!
    void mapMesh(gsMesh<T> & mesh) const
    {
        const int pDim = this->domainDim();

        gsMatrix<T> tmp, point;

        for (size_t i = 0; i!= mesh.numVertices(); ++i)
        {
            point = tmp = mesh.vertex(i).topRows(pDim);
            m_composition->invertPoints(point,tmp,1e-6,true);
            mesh.vertex(i).topRows(pDim) = tmp.topRows(pDim);
        }
    }





    std::ostream &print(std::ostream &os) const
    {
        os <<"Composite basis:\n";
        os << "* Compositoon "
           << " ( R^" << m_composition->domainDim() << " --> R^" << m_composition->targetDim() << "):\n"
           << m_composition<<"\n";
        os << "* Basis "
           << " ( R^" << m_basis->domainDim() << " --> R^" << m_basis->targetDim() << "):\n"
           << m_basis<<"\n";
        return os;
    }

private:
    void _applyBounds(gsMatrix<T> & coords) const
    {
        for (index_t k=0; k!=coords.cols(); k++)
        {
            coords.col(k) = coords.col(k).cwiseMax(m_basis->support().col(0));
            coords.col(k) = coords.col(k).cwiseMin(m_basis->support().col(1));
        }
    }

protected:
    const gsFunction<T> * m_composition;
    const gsBasis<T> * m_basis;

};
}
