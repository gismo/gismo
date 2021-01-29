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
gsSparseMatrix<T> sparseDiagonal(const gsVector<T> & vec)
{
    gsSparseMatrix<T> spdiag(vec.size(),vec.size());
    for (index_t k=0; k!=vec.size(); k++)
        spdiag.insert(k,k) = vec.at(k);
    spdiag.makeCompressed();
    return spdiag;
}

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

        /*
            To do: overload std::map.lower_bound s.t. it starts searching from index n
        */
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


template<class T>
class gsPatchRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsPatchRule()
    :
    m_basis(nullptr)
    {};

    /// Initialize a tensor-product Gauss quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsPatchRule(const  gsBSplineBasis<T> & basis,index_t degree, index_t regularity, bool overintegrate)
    :
    m_basis(basis),
    m_deg(degree),
    m_reg(regularity),
    m_over(overintegrate)
    {
        GISMO_ENSURE(basis.dim()==1,"Dimension must be equal to 1, not "<<basis.dim());
        GISMO_ENSURE(m_reg<m_deg,"regularity cannot be greater or equal to the order!");

        // m_nodes = points;
        // m_weights = weights;
        m_size = m_basis.size();
        m_knots = m_basis.knots();
        this->init();

        m_knots.greville_into(m_greville);
    };


    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsPatchRule() { };

public:
    // see gsQuadRule.h for documentation
    void setNodes( gsVector<index_t> const & numNodes,
                   unsigned digits = 0 )
    {


    };

    using gsQuadRule<T>::setNodes; // unhide base

    /// \brief Dimension of the rule
    index_t dim() const { return m_basis->dim(); }

    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const
    {
        m_qRule.mapTo(lower,upper,nodes,weights);
    };

    gsMatrix<T> integrate() const
    {
        typename gsBasis<real_t>::domainIter domIt =  // add patchInd to domainiter ?
                    m_basis2.makeDomainIterator();

        gsMatrix<T> result(m_basis2.size(),1);
        result.setZero();
        gsMatrix<T> nodes, tmp;
        gsMatrix<index_t> actives;
        gsVector<T> weights;
        gsInfo<<m_basis2.size()<<"\n";
        for (; domIt->good(); domIt->next() )
        {
            m_qRule.mapTo(domIt->lowerCorner(),domIt->upperCorner(),nodes,weights);
            m_basis2.active_into(nodes,actives);
            m_basis2.eval_into(nodes,tmp);
            tmp *= weights;
            for (index_t k = 0; k!=weights.size(); k++)
                result(actives.at(k),0) += tmp(k,0);
        }
        return result;
    };

    void init(bool overInt = false)
    {
        index_t pdiff = m_deg - m_basis.knots().degree();
        std::vector<index_t> multiplicities = m_knots.multiplicities();
        index_t rmin = *std::min_element(multiplicities.begin(), multiplicities.end());
        index_t rdiff = (m_deg-m_reg)-rmin ;

        if (pdiff>0)
            m_knots.degreeIncrease(pdiff);
        else if (pdiff<0)
            m_knots.degreeDecrease(-pdiff);
        if (rdiff>0)
            m_knots.increaseMultiplicity(rdiff);
        else if (rdiff>0)
            m_knots.reduceMultiplicity(-rdiff);

        // m_size = (m_knots.uSize()-1)*(m_deg-m_reg)+m_reg+1;
        m_basis2 = gsBSplineBasis<T>(m_knots);
        m_size = m_basis2.size();

        if (m_size % 2 == 1)
        {
            typedef typename gsKnotVector<T>::const_iterator knotIterator;
            knotIterator prevKnot = m_knots.begin();
            std::vector<T> diff;
            for (knotIterator it = m_knots.begin()+1; it!=m_knots.end(); it++ )
            {
                diff.push_back(*it-*prevKnot);
                prevKnot = it;
            }

            T max = *std::max_element(diff.begin(),diff.end());
            std::vector<index_t> maxIdx;
            index_t k=0;
            for (typename std::vector<T>::iterator it = diff.begin(); it!=diff.end(); it++,k++)
            {
                if (std::abs(*it-max)/(max)<1e-15)
                    maxIdx.push_back(k);
            }

            index_t i = maxIdx.at(std::ceil(maxIdx.size()/2.)-1);
            m_knots.insert( (m_knots.at(i) + m_knots.at(i+1))/2. );

            m_size++;
        }

        if (m_over) //(overInt)
        {
            // index_t rmax = *std::max_element(multiplicities.begin(), multiplicities.end());
            index_t numOver = m_knots.degree();

            T lowerLength = (m_knots(1)-m_knots.first())/(numOver+1);
            T upperLength = (m_knots.last()-m_knots(m_knots.uSize()-2))/(numOver+1);
            for (index_t k=0; k!=numOver; k++)
            {
                m_knots.insert(m_knots.first()+(k+1)*lowerLength);
                m_knots.insert(m_knots.last()-(k+1)*upperLength);
            }
        }

        for (gsKnotVector<>::iterator it = m_knots.begin(); it!=m_knots.end(); it++)
            gsInfo<<*it<<"\n";

        m_basis2 = gsBSplineBasis<T>(m_knots);
        m_size = m_basis2.size();

        m_qRule = gsGaussRule<T>(m_basis2,1,1);

    };

    void compute(T tol = 1e-10)
    {
        // initialize iteration
        m_integral = this->integrate();

        gsInfo<<"integral = "<<m_integral<<"\n";

        GISMO_ENSURE((m_size % 2 == 0),"Number of points should be even!");
        GISMO_ENSURE((m_basis.dim() == 1),"Basis dimension should be 1 but is "<<m_basis.dim());
        m_nodes.resize(m_size/2);
        m_weights.resize(m_size/2);
        for (index_t k = 0; k!=m_size/2; k++)
        {
            m_nodes.at(k)  = 0.5 * (m_greville.at(2*k) + m_greville.at(2*k+1));
            m_weights.at(k) = m_integral(2*k,0) + m_integral(2*k+1,0);
        }

        gsInfo<<"nodes = "<<m_nodes<<"\n";
        gsInfo<<"weights = "<<m_weights<<"\n";

        index_t itMax = 15;
        gsMatrix<index_t> actives;
        gsMatrix<T> vals,dvals,vals_tmp,dvals_tmp,res,dres;
        gsVector<T> update(m_size);
        vals.resize(m_size,m_size/2);   vals.setZero();
        dvals.resize(m_size,m_size/2);  dvals.setZero();
        typename gsSparseSolver<T>::QR solver;
        index_t it;
        for (it = 0; it != itMax; it++)
        {
            m_basis2.active_into(m_nodes.transpose(),actives);
            m_basis2.eval_into(m_nodes.transpose(),vals_tmp);
            m_basis2.deriv_into(m_nodes.transpose(),dvals_tmp);
            for (index_t act=0; act!=actives.rows(); act++)
                for (index_t pt=0; pt!=actives.cols(); pt++)
                {
                    vals(actives(act,pt),pt)    = vals_tmp(act,pt);
                    dvals(actives(act,pt),pt)   = dvals_tmp(act,pt);
                }

            res = vals * m_weights - m_integral;
            dres = gsMatrix<T>::Zero(m_size,m_size);

            dres.block(0,0,m_size,m_size/2) = vals;
            dres.block(0,m_size/2,m_size,m_size/2) = dvals * m_weights.asDiagonal();

            solver.compute(dres.sparseView());
            update = solver.solve(-res);

            m_weights+= update.head(m_size/2);
            m_nodes.col(0)  += update.tail(m_size/2);

            GISMO_ENSURE(m_nodes.minCoeff() > m_knots.first(),"Construction failed: min(nodes) < min(knots) minCoef = "<<m_nodes.minCoeff()<<"; min(knots) = "<<m_knots.first());
            GISMO_ENSURE(m_nodes.maxCoeff() < m_knots.last(),"Construction failed: max(nodes) > max(knots) maxCoef = "<<m_nodes.maxCoeff()<<"; min(knots) = "<<m_knots.last());

            if (res.norm() < tol)
            {
                gsInfo<<"Converged in "<<it<<" iterations\n";
                break;
            }
        }
        GISMO_ENSURE(it+1!=itMax,"Maximum iterations reached");


    }


public:
    mutable gsMatrix<T> m_greville;
    mutable gsKnotVector<T> m_knots;
    mutable gsVector<T> m_nodes;
    mutable gsVector<T> m_weights;
    mutable gsVector<T> m_integral;
    mutable index_t m_size;

    const index_t m_deg,m_reg;
    const bool m_over;

private:
    mutable gsVector<T> m_start,m_end;
    // const gsBasis<T> * m_basis;
    const gsBSplineBasis<T> & m_basis;
    mutable gsBSplineBasis<T> m_basis2;
    gsQuadRule<T> m_qRule;

}; // class gsPatchRule



int main(int argc, char* argv[])
{
    // ======================================================================
    // different construction of a knot vector
    // ======================================================================


    index_t order = 2;
    index_t regularity = 1;
    bool plot = false;
    bool overInt = false;

    gsCmdLine cmd("Quadrature rules in G+Smo.");
    cmd.addInt("p","deg","order of target space",order);
    cmd.addInt("r","reg","regularity of target space",regularity);
    cmd.addSwitch("plot","Plot with paraview",plot);
    cmd.addSwitch("over","overintegrate",overInt);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------- Constructions -----------------------------\n";

    gsKnotVector<> kv1(0, 1.0, 3, 3, 1);
    gsKnotVector<> kv2 = kv1;

    gsBSplineBasis<> bsb0(kv1);
    gsTensorBSplineBasis<2,real_t> tbsb(kv1,kv2);
    gsInfo<<tbsb<<"\n";


    gsWriteParaview(tbsb,"basis",1000,true);
    // bsb0.degreeIncrease();
    gsWriteParaview(bsb0,"basis_1D",1000,true);

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
    random *= 0.5;

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
        if (points.cols()!=0)
            gsWriteParaviewPoints(points,"randPoints_el" + std::to_string(el));
        el++;
    }

    gsMatrix<> rndPoints(1,random.size());
    rndPoints.row(0) = random;
    gsWriteParaviewPoints(rndPoints,"randPoints");



    gsWriteParaviewPoints(allPoints,"quadPoints");



    gsPatchRule<real_t> patchRule(bsb0,order,regularity,overInt);
    gsInfo<<"size = "<<patchRule.m_size<<"\n";
    gsInfo<<"integral = "<<patchRule.m_integral<<"\n";
    patchRule.compute();
    gsInfo<<"nodes = "<<patchRule.m_nodes<<"\n";
    gsInfo<<"weights = "<<patchRule.m_weights<<"\n";
    // gsDebugVar(bsb0.totalDegree());
    gsWriteParaviewPoints(patchRule.m_greville,"greville");




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