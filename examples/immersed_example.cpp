/** @file immersed_example.cpp

    @brief Tutorial on gsBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include <iostream>
#include <gismo.h>

// why needed ?
#include <gsHSplines/gsHDomain.hpp>

namespace gismo {

template<class T> class gsTrimmedDomain
{
//todo
};

template<class T>
class gsImplTrimmedDomain : gsTrimmedDomain<T>
{
private:
    typename gsFunction<T>::Ptr m_implFunction; // implicit function
    
public:
    gsImplTrimmedDomain(const gsFunction<T> & fnc)
    : m_implFunction(fnc.clone()) { }

    gsImplTrimmedDomain(const std::string & expression_string, short_t ddim)
    : m_implFunction(new gsFunctionExpr<T>(expression_string, ddim))
    { } 

    inline gsVector<int> sign(const gsMatrix<T> & u)
    {
        gsMatrix<T> val = m_implFunction->eval(u);
        return val.array().sign();
    }

    inline std::vector<bool> inDomain(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsMatrix<T> val = m_implFunction->eval(u);
        bool * r = res.data();
        for (T * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a<0 ? false : true);
        return res;
    }

    inline std::vector<bool> onBoundary(const gsMatrix<T> & u)
    {return std::vector<bool>(); }//todo

};

template<class T>
class gsImmersedGeometry// : public gsFunction<T>
{
private:
    typename gsFunction<T>::Ptr m_bgGeo; // background geometry
    gsImplTrimmedDomain<T> m_trDomain;  // trimmed domain

public:
    gsImmersedGeometry(const gsFunction<T> & bgGeo,
                       const gsImplTrimmedDomain<T> & trDomain) :
    m_bgGeo(bgGeo.clone()), m_trDomain(trDomain)
    { }

public://function interface?

public:
    const gsFunction<T> & background() const { return *m_bgGeo; }

    const gsImplTrimmedDomain<T> & domain() const { return m_trDomain; }

    /// Return a triangulation of the boundaty of the immersed geometry
    memory::unique_ptr<gsMesh<T> > toBoundaryMesh(int npoints = 50) const
    {
        gsVector<index_t>  numNodes(2); numNodes.setConstant(3);
        gsLobattoRule<real_t> qurule(numNodes);

        auto supp = m_bgGeo.support();
        // gsKdTree
    }

    /// Return a triangulation of the immersed geometry
    memory::unique_ptr<gsMesh<T> > toVolumeMesh(int npoints = 50) const
    {

    }

};

template<class T>
class gsTrimmedDomainIterator
{
    // elements:
    // bool active (true/false)
    // int position (-1,0,1) (out, cut, in)
private:
    gsDofMapper m_mapper;

    // gsKdTree
    
public:
    gsTrimmedDomainIterator(const gsImplTrimmedDomain<T> & trdom,
                            const gsBasis<T> & basis, index_t cutlevel = 0)
    {
        init(trdom,basis,cutlevel);
    }

    void init(const gsImplTrimmedDomain<T> & trdom,
              const gsBasis<T> & basis, index_t cutlevel = 0)
    {
        gsVector<index_t>  numNodes(2); numNodes.setConstant(3);
        gsLobattoRule<real_t> qurule(numNodes);

        //First pass:
        // mark element position (-1,0,1)
        // create dofmapper (mark active basis functions)

        //Second pass:
        // assign active elements
        // subdivide into cut-cells upto cutlevel
        
    }

    const gsDofMapper & mapper() const { return m_mapper; }
};

template<class T>
class gsCutCellRule
{ }; // todo

}//namespace gismo

using namespace gismo;

int main(int argc, char* argv[])
{
    gsCmdLine cmd("Immersing the geometry and the basis.");
    // cmd.addPlainString("input", "G+Smo input basis file.", input);
    // cmd.addString("o", "output", "Name of the output file.", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // defines the background geometry
    gsTensorBSpline<2,real_t> bg = *gsNurbsCreator<>::BSplineSquare(5,0,0);
    bg.knots(0).transform(-1,1);
    bg.knots(1).transform(-1,1);
    // defines the inside or the outside of the (parametric) domain
    gsImplTrimmedDomain<real_t> inOut("1 - x^2 - y^2", 2);
    // uses the two objects above to get an immersed geometry
    gsImmersedGeometry<real_t> igo(bg, inOut);

    // Background basis
    gsKnotVector<> kv(-1,1,9,2);
    gsTensorBSplineBasis<2,real_t> tbs(kv,kv);

    gsTrimmedDomainIterator<real_t> domIt(igo.background(), tbs, 0);
    

    
    gsVector<unsigned,2> upp;
    upp.setConstant(kv.numElements());
    //gsHDomain<2,unsigned> tree;
    //tree.init(upp,10);
    //tree.size();
    //tree.construct(inOut);
    //tree.printLeaves();
    //tree.coordinates(nullptr);

    return EXIT_SUCCESS;
}



