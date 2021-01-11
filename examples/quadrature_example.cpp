/** @file quadrature_example.cpp

    @brief Playing with quadrature!

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

// Look also in tuturialBasis for more functionality of gsBSplineBasis.

#include <gismo.h>

using namespace gismo;

// template<class T>
// class gsMap : public std::map<index_t,T>
// {
// public:
//     gsMap() {};

//     gsMap(const gsVector<T> & values)
//     {
//         this->setMap(values);
//     }

//     void setMap(const gsVector<T> & values)
//     {
//         gsSortedVector<T> sortValues(values);
//         for (index_t k=0; k!=sortValues.size(); k++)
//             this->operator[](k) = sortValues.at(k);
//     };

//     std::iterator find(const T key, std::iterator it)
//     {

//     }

//     std::const_iterator find(const T key, std::iterator it)
//     {

//     }
// };

template<class T>
class gsOverIntegrated GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsOverIntegrated()
    :
    m_basis(nullptr),
    m_interior(nullptr),
    m_boundary(nullptr)
    {}

    /// Initialize a tensor-product Gauss quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsOverIntegrated(const  gsBasis<T> & basis,
                            gsQuadRule<T> & quadInterior,
                            gsQuadRule<T> & quadBoundary)
    :
    m_basis(&basis),
    m_interior(quadInterior),
    m_boundary(quadBoundary)
    {
        m_start = m_basis->support().col(0);
        m_end = m_basis->support().col(1);
    };


    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsOverIntegrated() { }

public:
    // see gsQuadRule.h for documentation
    void setNodes( gsVector<index_t> const & numNodes,
                   unsigned digits = 0 )
    {
        m_interior.setNodes(numNodes,digits);
        m_boundary.setNodes(numNodes,digits);
    };

    using gsQuadRule<T>::setNodes; // unhide base

    /// \brief Dimension of the rule
    index_t dim() const { return m_basis->dim(); }

    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const
    {
        if ((lower-m_start).prod()==0 || (upper-m_end).prod()==0)
            m_boundary.mapTo( lower, upper, nodes, weights);
        else
            m_interior.mapTo( lower, upper, nodes, weights);
    }

private:
    const gsBasis<T> * m_basis;
    gsQuadRule<T> m_interior, m_boundary;
    mutable gsVector<T> m_start,m_end;

}; // class gsOverIntegrated

template<class T>
class gsTensorPatchRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsTensorPatchRule()
    :
    m_basis(nullptr)
    {}

    /// Initialize a tensor-product Gauss quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsTensorPatchRule(  const gsVector<T> & pointsX,
                        const gsVector<T> & pointsY,
                        const gsVector<T> & weightsX,
                        const gsVector<T> & weightsY,
                        const  gsBasis<T> & basis)
    :
    m_nodesX(pointsX),
    m_nodesY(pointsY),
    m_weightsX(weightsX),
    m_weightsY(weightsY),
    m_basis(&basis)
    {
        for (index_t k=0; k!=m_nodesX.size(); k++)
            m_mapX[m_nodesX.at(k)] = m_weightsX.at(k);

        for (index_t k=0; k!=m_nodesY.size(); k++)
            m_mapY[m_nodesY.at(k)] = m_weightsY.at(k);
    };


    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsTensorPatchRule() { }

public:
    // see gsQuadRule.h for documentation
    // void setNodes( gsVector<index_t> const & numNodes,
    //                unsigned digits = 0 )
    // {

    // }

    using gsQuadRule<T>::setNodes; // unhide base

    /// \brief Dimension of the rule
    index_t dim() const { return m_basis->dim(); }

    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const
    {
        nodes.resize(2,m_nodesX.size()*m_nodesY.size());
        weights.resize(m_nodesX.size()*m_nodesY.size());

        index_t k=0;
        for (auto itX = m_mapX.lower_bound(lower[0]); itX!=m_mapX.upper_bound(upper[0]); itX++) // lower_bound = geq, upper_bound= greather than
            for (auto itY = m_mapY.lower_bound(lower[1]); itY!=m_mapY.upper_bound(upper[1]); itY++, k++) // lower_bound = geq, upper_bound= greather than
                {
                    nodes(0,k) = itX->first;
                    nodes(1,k) = itY->first;
                    weights.at(k) = itX->second*itY->second;
                }

        nodes.conservativeResize(2,k);
        weights.conservativeResize(k);
    }

private:
    const gsVector<T> m_nodesX,m_nodesY;
    const gsVector<T> m_weightsX,m_weightsY;
    const gsBasis<T> * m_basis;

    std::vector<T> m_X, m_Y;
    std::vector<T> m_wX,m_wY;

    std::map<T,T> m_mapX,m_mapY;

    gsSortedVector<T> m_sortedNodes;


}; // class gsTensorPatchRule


// template<class T>
// class gsPatchRule GISMO_FINAL : public gsQuadRule<T>
// {
// public:

//     /// Default empty constructor
//     gsPatchRule()
//     :
//     m_basis(nullptr)
//     {};

//     /// Initialize a tensor-product Gauss quadrature rule for \a basis
//     /// using quA *deg_i + quB nodes (direction-wise)
//     gsPatchRule(const gsMatrix<T> & points,
//                 const gsVector<T> & weights,
//                 const  gsBasis<T> & basis)
//     :
//     m_nodes(&points),
//     m_weights(&weights),
//     m_basis(&basis)
//     {
//         // m_nodes = points;
//         // m_weights = weights;
//         m_start = m_basis->support().col(0);
//         m_end = m_basis->support().col(1);
//     };


//     //const unsigned digits = std::numeric_limits<T>::digits10 );

//     ~gsPatchRule() { };

// public:
//     // see gsQuadRule.h for documentation
//     void setNodes( gsVector<index_t> const & numNodes,
//                    unsigned digits = 0 )
//     {


//     };

//     using gsQuadRule<T>::setNodes; // unhide base

//     /// \brief Dimension of the rule
//     index_t dim() const { return m_basis->dim(); }

//     void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
//                        gsMatrix<T> & nodes, gsVector<T> & weights ) const
//     {

//     };

// private:
//     mutable gsVector<T> m_start,m_end;
//     const gsMatrix<T> * m_nodes;
//     const gsVector<T> * m_weights;
//     const gsBasis<T> * m_basis;

// }; // class gsPatchRule



int main(int argc, char* argv[])
{
    // ======================================================================
    // different construction of a knot vector
    // ======================================================================


    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    index_t interior = 4; // number of interior knots
    index_t multEnd = 3; // multiplicity at the two end knots
    index_t multInt = 1; // multiplicity at the interior knots
    bool plot = false;

    gsCmdLine cmd("Quadratire rules in G+Smo.");
    cmd.addReal("","starting","Starting knot",a);
    cmd.addReal("","ending","Ending knot",b);
    cmd.addInt("n","interior","Number of interior knots",interior);
    cmd.addInt("m","multI","Multiplicity at the interior knots",multInt);
    cmd.addInt("M","multE","Multiplicity at the two end knots",multEnd);
    cmd.addSwitch("plot","Plot with paraview",plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------- Constructions -----------------------------\n";

    gsKnotVector<> kv1(a, b, interior, multEnd, multInt);
    gsKnotVector<> kv2 = kv1;

    gsBSplineBasis<> bsb0(kv1);
    gsTensorBSplineBasis<2,real_t> tbsb(kv1,kv2);
    gsInfo<<tbsb<<"\n";


    gsWriteParaview(tbsb,"basis",1000,true);

    // ======================================================================
    // some properties
    // ======================================================================


    gsInfo << "------------- Some properties    -----------------------\n\n";

    gsInfo << "bsb0.size(): " << bsb0.size() << "\n\n"
              << "bsb0.numElements(): " << bsb0.numElements() << "\n\n"
              << "bsb0.degree(): " << bsb0.degree() << "\n\n";

    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << tbsb.dim() << "\n"
         << "Number of basis functions: " << tbsb.size() << "\n"
         << "Number of elements: " << tbsb.numElements() << "\n"
         << "Max degree of the basis: " << tbsb.maxDegree() << "\n"
         << "Min degree of the basis: " << tbsb.minDegree() << "\n"
         << "\n";

    // ======================================================================
    // some operations
    // ======================================================================

    // quadrature
    gsQuadRule<real_t> QuadRule,QuadRule2;
    gsMatrix<> points;
    gsVector<> weights;

    gsOptionList options;
    options.addReal("quA", "Number of quadrature points: quA*deg + quB", 0  );
    options.addInt ("quB", "Number of quadrature points: quA*deg + quB", 2    );

    QuadRule = gsGaussRule<real_t>(tbsb, 0,1);
    QuadRule2 = gsGaussRule<real_t>(tbsb, 0,2);

    gsOverIntegrated<real_t> MixedQuadRule(tbsb,QuadRule,QuadRule2);
    gsDebugVar(MixedQuadRule.dim());

    gsVector<> low(2);
    low<<0,0;
    gsVector<> up(2);
    up<<1,1;
    MixedQuadRule.mapTo(low,up,points,weights);


    typename gsBasis<real_t>::domainIter domIt =  // add patchInd to domainiter ?
                tbsb.makeDomainIterator();

    // gsMatrix<> allPoints(tbsb.dim(),tbsb.numElements()*QuadRule.numNodes());
    gsMatrix<> allPoints(tbsb.dim(),0);
    for (; domIt->good(); domIt->next() )
    {

        gsInfo<<"Element corners:\n"<<domIt->lowerCorner().transpose()<<"\n"<<domIt->upperCorner().transpose()<<"\n";

        // Map the Quadrature rule to the element
        MixedQuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);

        index_t start = allPoints.cols();
        allPoints.conservativeResize(Eigen::NoChange,allPoints.cols()+points.cols());
        allPoints.block(0,start,allPoints.rows(),points.cols()) = points;

        gsDebugVar(domIt->isBoundaryElement());
    }

    gsWriteParaviewPoints(allPoints,"quadPoints");


    // Make a random matrix
    gsVector<> random;
    random.setRandom(20);
    random += gsVector<>::Ones(random.size());
    random *= 0.5*(b-a);
    random += a*gsVector<>::Ones(random.size());

    gsVector<> unitWeights(random.size());
    unitWeights.setOnes();

    gsTensorPatchRule<real_t> tensorPatchRule(random,random,unitWeights,unitWeights,tbsb);
    // gsMap<real_t> map(random.row(0));
    gsDebugVar(random);
    // gsDebugVar(map[0]);
    // gsDebugVar(random(0,0));

    domIt->reset();
    index_t el=0;
    for (; domIt->good(); domIt->next() )
    {
        tensorPatchRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);

        gsDebugVar(points);
        gsDebugVar(weights);

        if (points.cols()!=0)
            gsWriteParaviewPoints(points,"randPoints_el" + std::to_string(el));
        el++;
    }

    gsMatrix<> rndPoints(1,random.size());
    rndPoints.row(0) = random;
            gsWriteParaviewPoints(rndPoints,"randPoints");



    gsWriteParaviewPoints(allPoints,"quadPoints");


    return 0;


    // typename gsBasis<real_t>::domainIter domIt =  // add patchInd to domainiter ?
    //             bsb0.makeDomainIterator();

    // gsOptionList options;
    // options.addReal("quA", "Number of quadrature points: quA*deg + quB", 1  );
    // options.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );

    // QuadRule = gsQuadrature::get(tbsb, options);
    // QuadRule.mapTo( a, b,points, weights);

    // gsDebugVar(points);

    // gsVector<> exact(bsb0.size());
    // exact.setZero();
    // gsMatrix<index_t> actives;
    // gsMatrix<real_t> values;
    // for (; domIt->good(); domIt->next() )
    // {

    //     gsInfo<<"Element corners:\n"<<domIt->lowerCorner().transpose()<<"\n"<<domIt->upperCorner().transpose()<<"\n";

    //     gsInfo<<domIt->side()<<"\n";

    //     // Map the Quadrature rule to the element
    //     QuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
    //                     points, weights);

    //     bsb0.active_into(points,actives);
    //     bsb0.eval_into(points,values);

    //     for (index_t p = 0; p!=actives.cols(); p++)
    //     {
    //         for (index_t r=0; r!=actives.rows(); r++)
    //         {
    //             exact.at(actives(r,p)) += weights.at(p) * values(r,p);
    //         }
    //     }
    // }
    // gsDebugVar(exact);

    // index_t ndof = bsb0.size();
    // index_t nquad = std::ceil(ndof/2.);

    // gsMatrix<> allPoints;
    // for (; domIt->good(); domIt->next() )
    // {

    //     gsInfo<<"Element corners:\n"<<domIt->lowerCorner().transpose()<<"\n"<<domIt->upperCorner().transpose()<<"\n";

    //     // Map the Quadrature rule to the element
    //     QuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
    //                     points, weights);

    //     bsb0.active_into(points,actives);
    //     bsb0.eval_into(points,values);

    //     for (index_t p = 0; p!=actives.cols(); p++)
    //     {
    //         for (index_t r=0; r!=actives.rows(); r++)
    //         {
    //             exact.at(actives(r,p)) += weights.at(p) * values(r,p);
    //         }
    //     }
    // }


    // points.resize(nquad,1);
    // points.col(0).setLinSpaced(nquad,a,b);
    // gsDebugVar(points);

    // bsb0.active_into(points.transpose(),actives);
    // bsb0.eval_into(points.transpose(),values);

    // gsMatrix<> shape(ndof,nquad);
    // shape.setZero();
    // for (index_t p = 0; p!=actives.cols(); p++)
    // {
    //     for (index_t r=0; r!=actives.rows(); r++)
    //     {
    //         shape(actives(r,p),p) = values(r,p);
    //     }
    // }
    // gsDebugVar(shape);

    // // gsVector<> resVec(2*nquad);
    // // res.segment(0,ndof) = shape * wq



    return 0;


    /*
        Other way of integration


        gsExprAssembler<> A;
        gsMultiBasis<> basis;
        basis.addBasis(&bsb0);
        A.setIntegrationElements(basis);
        typedef gsExprAssembler<>::space       space;
        gsExprEvaluator<> ev(A);
        space u = A.getSpace(basis);
        A.initSystem(true);
        A.assemble( u );
        gsDebugVar(A.rhs());
    */

}