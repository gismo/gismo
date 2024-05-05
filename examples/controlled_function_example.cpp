/** @file controlled_function_example.cpp

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
class gsSquareDomain : public gsFunction<T>
{
    using Base = gsFunction<T> ;

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
    gsTensorBSpline<2,T> m_domain;
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
    gsSquareDomain<real_t> domain;

    // Read the geometry from a file
    gsMultiPatch<> mp;
    gsReadFile<>("surfaces/simple.xml",mp);
    // mp.degreeElevate();
    // mp.degreeElevate();
    mp.uniformRefine();
    mp.uniformRefine();

    // gsTensorBSpline<2,real_t> & tbsgeom = dynamic_cast<gsTensorBSpline<2,real_t> &>(mp.patch(0));
    // gsTHBSpline<2,real_t> thbgeom(tbsgeom);
    // std::vector<index_t> box1{1,0,0,4,4};
    // thbgeom.refineElements(box1);

    // Collect the geometry and the basis
    // gsTHBSplineBasis<2,real_t> thbbasis = thbgeom.basis();
    const gsBasis<> & tbasis = mp.basis(0); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
    const gsGeometry<> & tgeom = mp.patch(0); //G(u,v) -> deriv will give dG/du, dG/dv



    // Define a composite basis and composite geometry
    // The basis is composed by the square domain
    gsComposedBasis<real_t> cbasis(domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
    // The geometry is defined using the composite basis and some coefficients
    gsComposedGeometry<real_t> cgeom(cbasis,tgeom.coefs()); // G(u,v) = G(sigma(xi,eta))  -> deriv will give dG/dxi, dG/deta
    // Define a mesh for later
    gsMesh<> mesh(cbasis,10);


    gsVector<> pt(2);
    gsMatrix<> cders;
    pt.setConstant(0.25);
    gsComposedGeometry<real_t> cgeom2(domain,tgeom);
    cgeom2.control_deriv_into(pt,cders);
    gsDebugVar(cders.reshape(domain.nControls(),cgeom2.targetDim()));


    // Plot the original objects
    gsWriteParaview(domain.domain(),"domain_before",1000,true,true);
    gsWriteParaview(cgeom,"geom_before",1000,false);
    gsWriteParaview(cbasis,"basis_before",1000,true);
    gsWriteParaview(cgeom,"geom_before",1000,false);
    cbasis.mapMesh(mesh);
    gsWriteParaview(mesh,"mesh_before",true);

    // Change the parametric domain
    gsMatrix<> pars = domain.controls();
        // pars.resize(pars.rows()/domain.domainDim(),domain.domainDim());
        // for (index_t r = 0; r!= pars.rows(); r++)
        // {
        //     // Rotate the inner control points on a circle
        //     real_t angle = math::atan2(pars(r,1)-0.5,pars(r,0)-0.5);
        //     angle += 3.1415926535/2.;
        //     pars.row(r)<<0.25*math::cos(angle)+0.5,0.25*math::sin(angle)+0.5;
        // }
        // pars.resize(pars.rows()*pars.cols(),1);
    // pars *= 0.75;
    pars(0,0) *= 0.75;
    // Register the changed parameters back to the parametric domain and update
    domain.controls() = pars.col(0);
    domain.updateGeom();

    gsMatrix<> pts = gsPointGrid<>(cbasis.support(),2000);
    gsMatrix<> result;
    cbasis.eval_into(pts,result);
    std::string test = ((result.colwise().sum().array()>1-1e-12).all() && (result.colwise().sum().array()<1+1e-12).all()) ? "yes" : "no";
    gsInfo<<"Partition of unity: "<<test<<"\n";
    // gsDebugVar(result.colwise().sum());
    domain.eval_into(pts,result);

    gsWriteParaviewPoints(result,"points");

    gsWriteParaview(domain,tbasis.function(16),"basis5",1000);


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

    gsMatrix<> point(2,2);
    point.col(0).setConstant(0.25);
    point.col(1).setConstant(0.50);
    gsMatrix<> der;
    gsMatrix<index_t> act;
    cbasis.active_into(point,act);
    gsDebugVar(act);

    cbasis.deriv_into(point,der);
    gsDebugVar(der);

    for (index_t p=0; p!=act.cols(); p++)
        for (index_t k=0; k!=act.rows(); k++)
        {
            cbasis.derivSingle_into(act(k,p),point,der);
            gsDebugVar(act(k,p));
            gsDebugVar(der);
        }

    testBasis<real_t> tb(cbasis,point.col(0));
    gsMatrix<> ev;
    tb.eval_into(point.col(0),ev);
    gsDebugVar(ev);
    cbasis.eval_into(point.col(0),ev);
    gsDebugVar(ev);

    tb.deriv_into(point.col(0),der);
    gsDebugVar(der);
    tb.deriv_into2(point.col(0),der);
    gsDebugVar(der);




    gsVector<> pp = point.col(0);
    gsMatrix<> ev1, ev2;

    real_t delta = 0.0001;
    gsVector<> pt1 = pp, pt2 = pp;
    pt1.at(1) += delta;
    pt2.at(1) -= delta;
    cbasis.evalSingle_into(7,pt1,ev1);
    cbasis.evalSingle_into(7,pt2,ev2);


    gsVector<> der2 = (ev1-ev2)/(2*delta);
    gsDebugVar(der2);

    cbasis.derivSingle_into(7,pp,der);
    gsDebugVar(der);

    return 0;
}


