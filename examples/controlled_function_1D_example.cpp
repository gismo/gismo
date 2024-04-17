/** @file geometry_example.cpp

    @brief Tutorial on gsGeometry abstract class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <gismo.h>
#include <gsCore/gsComposedBasis.h>

using namespace gismo;

template <class T>
class gsSquareDomain1D : public gsFunction<T>
{
    using Base = gsFunction<T> ;

public:
    gsSquareDomain1D()
    {
        index_t nInterior = 2;
        index_t degree = 2;
        gsKnotVector<T> kv(0,1,nInterior,degree+1);

        m_basis = gsBSplineBasis<T>(kv);
        m_domain = gsBSpline<T>(m_basis,m_basis.anchors().transpose());

        m_mapper = gsDofMapper(m_domain.basis(),m_domain.targetDim());

        gsMatrix<index_t> boundary = m_domain.basis().allBoundary();
        for (index_t a = 0; a!=boundary.rows(); a++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
                m_mapper.eliminateDof(boundary(a,0),0,d);
        m_mapper.finalize();

        m_parameters.resize(m_mapper.freeSize());
        // std::vector<index_t> i(m_mapper.freeSize());
        // std::vector<index_t> j(m_mapper.freeSize());
        for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
                if (m_mapper.is_free(k,0,d))
                {
                    m_parameters[m_mapper.index(k,0,d)] = m_domain.coefs()(k,d);
                    // i[m_mapper.index(k,0,d)] = k; // i index of free entries
                    // j[m_mapper.index(k,0,d)] = d; // j index of free entries
                }

        // This is a way to cast only the free coefficients to a vector, and change an entry of that vector.
        // However, it cannot be used in ''gsVector<T> & controls() override { return m_parameters; };''
        //
        // gsDebugVar(m_domain.coefs()(i,j).diagonal()(0));
        // m_domain.coefs()(i,j).diagonal()(0) = 0.5;
        // gsDebugVar(m_domain.coefs()(i,j).diagonal()(0));
    }

    const gsBSpline<T> & domain() const
    {
        return m_domain;
    }

    gsMatrix<T> support() const override
    {
        return m_domain.support();
    }

    short_t domainDim() const override
    {
        return m_domain.domainDim();
    }

    short_t targetDim() const override
    {
        return m_domain.domainDim();
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        m_domain.eval_into(u,result);
    }

    void updateGeom()
    {
        for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
            {
                if (m_mapper.is_free(k,0,d))
                    m_domain.coefs()(k,d) = m_parameters[m_mapper.index(k,0,d)];
            }
    }

    /// Returns the controls of the function
    const gsVector<T> & controls() const override { return m_parameters; };
          gsVector<T> & controls()       override { return m_parameters; };

    /// Returns the number of controls of the function
    size_t nControls() const override
    {
        return m_mapper.freeSize();
    }

    /// Returns the control derivative
    virtual void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override
    {
        gsMatrix<> tmp;

        result.resize(targetDim()*nControls(),points.cols());
        result.setZero();
        for (index_t p = 0; p!=points.cols(); p++)
        {
            gsAsMatrix<> res = result.reshapeCol(p,nControls(),targetDim());
            for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
                for (index_t d = 0; d!=m_domain.targetDim(); d++)
                    if (m_mapper.is_free(k,0,d))
                    {
                        m_domain.basis().evalSingle_into(k,points.col(p),tmp); // evaluate basis function k
                        res(m_mapper.index(k,0,d),d) = tmp(0,0); // tmp is a single value (1 point, 1 basis function)
                    }
        }
    }

protected:
    gsBSplineBasis<T> m_basis;
    gsBSpline<T> m_domain;
    gsDofMapper m_mapper;
    gsVector<T> m_parameters;
};


int main(int argc, char* argv[])
{

    // std::string input("surfaces/simple.xml");
    std::string output("");

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    // cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the geometry
    // ======================================================================

    // Define a square domain: this is the intermediate map
    gsSquareDomain1D<real_t> domain;

    gsKnotVector<> kv(0,1,10,4);
    gsBSplineBasis<> basis(kv);

    // Define a composite basis and composite geometry
    // The basis is composed by the square domain
    gsComposedBasis<real_t> cbasis(domain,basis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
    // The geometry is defined using the composite basis and some coefficients

    gsMatrix<> pars = domain.controls();
    pars *= 0.75;
    domain.controls() = pars.col(0);
    domain.updateGeom();

    gsWriteParaview(cbasis,"cbasis");


    return 0;
}


