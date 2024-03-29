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
#include <gsCore/gsControlledFunction.h>
#include <gsCore/gsComposedBasis.h>

using namespace gismo;

template <class T>
class gsSquareDomain : public gsControlledFunction<T>
{
    using Base = gsControlledFunction<T> ;

public:
    gsSquareDomain()
    {
        m_domain = *gsNurbsCreator<T>::BSplineSquare();
        m_domain.degreeElevate();
        m_domain.uniformRefine();
        // Mapper storing control points
        m_mapper = gsDofMapper(m_domain.basis(),m_domain.targetDim());

        gsDebugVar(m_domain.coefs());


        gsMatrix<index_t> boundary = m_domain.basis().allBoundary();
        for (index_t a = 0; a!=boundary.rows(); a++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
                m_mapper.eliminateDof(boundary(a,0),0,d);
        m_mapper.finalize();

        m_parameters.resize(m_mapper.freeSize());
        for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
                if (m_mapper.is_free(k,0,d))
                    m_parameters[m_mapper.index(k,0,d)] = m_domain.coefs()(k,d);

        gsDebugVar(m_parameters);
    }

    const gsTensorBSpline<2,T> & domain() const
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

protected:
    gsTensorBSpline<2,T> m_domain;
    gsDofMapper m_mapper;
    using Base::m_parameters;
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
    gsSquareDomain<real_t> domain;

    // Read the geometry from a file
    gsMultiPatch<> mp;
    gsReadFile<>("surfaces/simple.xml",mp);
    mp.degreeElevate();
    mp.uniformRefine();
    // mp.uniformRefine();


    // Collect the geometry and the basis
    const gsBasis<> & tbasis = mp.basis(0);
    const gsGeometry<> & tgeom = mp.patch(0);

    // Degine a composite basis and composite geometry
    // The basis is composed by the square domain
    gsComposedBasis<real_t> cbasis(domain,tbasis);
    // The geometry is defined using the composite basis and some coefficients
    gsComposedGeometry<real_t> cgeom(cbasis,tgeom.coefs());
    // Define a mesh for later
    gsMesh<> mesh(cbasis,10);


    // Plot the original objects
    gsWriteParaview(domain.domain(),"domain_before",1000,true,true);
    gsWriteParaview(cgeom,"geom_before",1000,false);
    gsWriteParaview(cbasis,"basis_before",1000,true);
    gsWriteParaview(cgeom,"geom_before",1000,false);
    cbasis.mapMesh(mesh);
    gsWriteParaview(mesh,"mesh_before",true);

    // Change the parametric domain
    gsMatrix<> pars = domain.parameters();
    pars.resize(pars.rows()/domain.domainDim(),domain.domainDim());
    for (index_t r = 0; r!= pars.rows(); r++)
    {
        // Rotate the inner control points on a circle
        real_t angle = math::atan2(pars(r,1)-0.5,pars(r,0)-0.5);
        angle += 3.1415926535/2.;
        pars.row(r)<<0.25*math::cos(angle)+0.5,0.25*math::sin(angle)+0.5;
    }
    pars.resize(pars.rows()*pars.cols(),1);
    // Register the changed parameters back to the parametric domain and update
    domain.parameters() = pars.col(0);
    domain.updateGeom();

    // Plot everything
    gsWriteParaview(domain.domain(),"domain_after",1000,true,true);
    gsWriteParaview(cgeom,"geom_after",1000,false);
    gsWriteParaview(cbasis,"basis_after",1000,true);
    gsWriteParaview(cgeom,"geom_after",1000,false);
    // TEST @Ye
//  gsMatrix<> u(2,1), xi;
//  u << 0.3, 0.5;
//  cbasis.eval_into(u, xi);
//  gsDebugVar(xi);
/////////////////////////////////////////////////////////////////////////////////
    cbasis.mapMesh(mesh);
    gsWriteParaview(mesh,"mesh_after",true);

    return 0;
}


