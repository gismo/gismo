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

template<class T>
class gsTensorPatchRule2 GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsTensorPatchRule2()
    {}

    /// Initialize a tensor-product Gauss quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsTensorPatchRule2(  const gsVector<T> & pointsX,
                        const gsVector<T> & pointsY,
                        const gsVector<T> & weightsX,
                        const gsVector<T> & weightsY)
    :
    m_nodesX(pointsX),
    m_nodesY(pointsY),
    m_weightsX(weightsX),
    m_weightsY(weightsY)
    {
        for (index_t k=0; k!=m_nodesX.size(); k++)
            m_mapX[m_nodesX.at(k)] = m_weightsX.at(k);

        for (index_t k=0; k!=m_nodesY.size(); k++)
            m_mapY[m_nodesY.at(k)] = m_weightsY.at(k);
    };


    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsTensorPatchRule2() { }

public:
    // see gsQuadRule.h for documentation
    // void setNodes( gsVector<index_t> const & numNodes,
    //                unsigned digits = 0 )
    // {

    // }

    using gsQuadRule<T>::setNodes; // unhide base

    /// \brief Dimension of the rule
    index_t dim() const { return 2; }

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

    std::vector<T> m_X, m_Y;
    std::vector<T> m_wX,m_wY;

    std::map<T,T> m_mapX,m_mapY;

    gsSortedVector<T> m_sortedNodes;


}; // class gsTensorPatchRule2

template<class T>
class gsTensorPatchRule3 GISMO_FINAL : public gsQuadRule<T>
{
public:

    /// Default empty constructor
    gsTensorPatchRule3()
    {}

    /// Initialize a tensor-product Gauss quadrature rule for \a basis
    /// using quA *deg_i + quB nodes (direction-wise)
    gsTensorPatchRule3(  const gsVector<T> & pointsX,
                        const gsVector<T> & pointsY,
                        const gsVector<T> & pointsZ,
                        const gsVector<T> & weightsX,
                        const gsVector<T> & weightsY,
                        const gsVector<T> & weightsZ)
    :
    m_nodesX(pointsX),
    m_nodesY(pointsY),
    m_nodesZ(pointsZ),
    m_weightsX(weightsX),
    m_weightsY(weightsY),
    m_weightsZ(weightsZ)
    {
        for (index_t k=0; k!=m_nodesX.size(); k++)
            m_mapX[m_nodesX.at(k)] = m_weightsX.at(k);

        for (index_t k=0; k!=m_nodesY.size(); k++)
            m_mapY[m_nodesY.at(k)] = m_weightsY.at(k);

        for (index_t k=0; k!=m_nodesZ.size(); k++)
            m_mapZ[m_nodesZ.at(k)] = m_weightsZ.at(k);
    };


    //const unsigned digits = std::numeric_limits<T>::digits10 );

    ~gsTensorPatchRule3() { }

public:
    // see gsQuadRule.h for documentation
    // void setNodes( gsVector<index_t> const & numNodes,
    //                unsigned digits = 0 )
    // {

    // }

    using gsQuadRule<T>::setNodes; // unhide base

    /// \brief Dimension of the rule
    index_t dim() const { return 3; }

    void mapTo( const gsVector<T>& lower, const gsVector<T>& upper,
                       gsMatrix<T> & nodes, gsVector<T> & weights ) const
    {
        nodes.resize(3,m_nodesX.size()*m_nodesY.size()*m_nodesZ.size());
        weights.resize(m_nodesX.size()*m_nodesY.size()*m_nodesZ.size());

        /*
            To do: overload std::map.lower_bound s.t. it starts searching from index n
        */
        index_t k=0;
        for (auto itX = m_mapX.lower_bound(lower[0]); itX!=m_mapX.upper_bound(upper[0]); itX++) // lower_bound = geq, upper_bound= greather than
            for (auto itY = m_mapY.lower_bound(lower[1]); itY!=m_mapY.upper_bound(upper[1]); itY++) // lower_bound = geq, upper_bound= greather than
                for (auto itZ = m_mapZ.lower_bound(lower[2]); itZ!=m_mapZ.upper_bound(upper[2]); itZ++, k++) // lower_bound = geq, upper_bound= greather than
                    {
                        nodes(0,k) = itX->first;
                        nodes(1,k) = itY->first;
                        nodes(2,k) = itZ->first;
                        weights.at(k) = itX->second*itY->second*itZ->second;
                    }

        nodes.conservativeResize(3,k);
        weights.conservativeResize(k);
    }

private:
    const gsVector<T> m_nodesX,m_nodesY,m_nodesZ;
    const gsVector<T> m_weightsX,m_weightsY,m_weightsZ;

    std::vector<T> m_X, m_Y,m_Z;
    std::vector<T> m_wX,m_wY,m_wZ;

    std::map<T,T> m_mapX,m_mapY,m_mapZ;

    gsSortedVector<T> m_sortedNodes;


}; // class gsTensorPatchRule3




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
    gsKnotVector<> kv3 = kv1;

    gsBSplineBasis<> bsb1(kv1);
    gsBSplineBasis<> bsb2(kv2);
    gsBSplineBasis<> bsb3(kv3);
    gsTensorBSplineBasis<2,real_t> tbsb2(kv1,kv2);
    gsTensorBSplineBasis<3,real_t> tbsb3(kv1,kv2,kv3);
    gsInfo<<tbsb2<<"\n";
    gsWriteParaview(bsb1,"basis_1D",1000);
    gsWriteParaview(tbsb2,"basis",1000,true);

    // ======================================================================
    // some properties
    // ======================================================================


    gsInfo << "------------- Some properties    -----------------------\n\n";

    gsInfo << "bsb1.size(): " << bsb1.size() << "\n\n"
              << "bsb1.numElements(): " << bsb1.numElements() << "\n\n"
              << "bsb1.degree(): " << bsb1.degree() << "\n\n";

    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << tbsb2.dim() << "\n"
         << "Number of basis functions: " << tbsb2.size() << "\n"
         << "Number of elements: " << tbsb2.numElements() << "\n"
         << "Max degree of the basis: " << tbsb2.maxDegree() << "\n"
         << "Min degree of the basis: " << tbsb2.minDegree() << "\n"
         << "\n";

    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << tbsb3.dim() << "\n"
         << "Number of basis functions: " << tbsb3.size() << "\n"
         << "Number of elements: " << tbsb3.numElements() << "\n"
         << "Max degree of the basis: " << tbsb3.maxDegree() << "\n"
         << "Min degree of the basis: " << tbsb3.minDegree() << "\n"
         << "\n";

    // ======================================================================
    // some operations
    // ======================================================================

    // Mixed Quadrature
    gsQuadRule<real_t> QuadRule,QuadRule2;
    gsMatrix<> points;
    gsVector<> weights;

    QuadRule = gsGaussRule<real_t>(tbsb2, 0,2);
    QuadRule2 = gsGaussRule<real_t>(tbsb2, 1,1);

    gsOverIntegrateRule<real_t> MixedQuadRule(tbsb2,QuadRule,QuadRule2);

    // Tensor Random Rule
    gsVector<> random;
    random.setRandom(20);
    random += gsVector<>::Ones(random.size());
    random *= 0.5;

    gsVector<> unitWeights(random.size());
    unitWeights.setOnes();

    gsTensorPatchRule2<real_t> tensorRandomRule(random,random,unitWeights,unitWeights);

    // Tensor Patch Rule
    gsOptionList options;
    options.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Legendre; (3) Patch-Rule",gsQuadrature::rule::PatchRule);
    options.addReal  ("quA", "Parameter for the order (in case of PatchRule) of the target space", order  );
    options.addInt   ("quB", "Parameter for the regularity (in case of PatchRule) of the target space", regularity    );
    options.addSwitch("overInt","Apply over-integration or not?",overInt);
    gsQuadRule<real_t> tensorPatchRule2D = gsQuadrature::get(tbsb2, options);

    // --------------------------------------------------------------------------------------

    typename gsBasis<real_t>::domainIter domIt =  // add patchInd to domainiter ?
                tbsb2.makeDomainIterator();

    // gsMatrix<> allPoints(tbsb.dim(),tbsb.numElements()*QuadRule.numNodes());
    gsMatrix<> GaussRule(tbsb2.dim(),0);
    gsMatrix<> MixedRule(tbsb2.dim(),0);
    gsMatrix<> TensorRandom(tbsb2.dim(),0);
    gsMatrix<> TensorPatch(tbsb2.dim(),0);
    index_t start;
    for (; domIt->good(); domIt->next() )
    {
        // Map the Quadrature rule to the element
        QuadRule2.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        gsInfo<<"points1 = "<<points<<"\n";
        gsInfo<<"weights1 = "<<weights.transpose()<<"\n";
        start = GaussRule.cols();
        GaussRule.conservativeResize(Eigen::NoChange,GaussRule.cols()+points.cols());
        GaussRule.block(0,start,GaussRule.rows(),points.cols()) = points;


        // Map the Quadrature rule to the element
        MixedQuadRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        gsInfo<<"points2 = "<<points<<"\n";
        gsInfo<<"weights2 = "<<weights.transpose()<<"\n";
        start = MixedRule.cols();
        MixedRule.conservativeResize(Eigen::NoChange,MixedRule.cols()+points.cols());
        MixedRule.block(0,start,MixedRule.rows(),points.cols()) = points;

        tensorRandomRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        gsInfo<<"points3 = "<<points<<"\n";
        gsInfo<<"weights3 = "<<weights.transpose()<<"\n";
        start = TensorRandom.cols();
        TensorRandom.conservativeResize(Eigen::NoChange,TensorRandom.cols()+points.cols());
        TensorRandom.block(0,start,TensorRandom.rows(),points.cols()) = points;

        tensorPatchRule2D.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                        points, weights);
        gsInfo<<"points4 = "<<points<<"\n";
        gsInfo<<"weights4 = "<<weights.transpose()<<"\n";
        start = TensorPatch.cols();
        TensorPatch.conservativeResize(Eigen::NoChange,TensorPatch.cols()+points.cols());
        TensorPatch.block(0,start,TensorPatch.rows(),points.cols()) = points;

    }

    gsWriteParaviewPoints(GaussRule,"Points_Original");
    gsWriteParaviewPoints(MixedRule,"Points_Mixed");
    gsWriteParaviewPoints(TensorRandom,"Points_Random");
    gsWriteParaviewPoints(TensorPatch,"Points_Patch");


    gsQuadRule<real_t> patchRule3D = gsQuadrature::get(tbsb3, options);
    ///////////////////////////
    gsVector<> lower(3);
    lower.setConstant(0);
    gsVector<> upper(3);
    upper.setConstant(1);
    gsMatrix<> Nodes;
    gsVector<> Weights;
    patchRule3D.mapTo(lower,upper,Nodes,Weights);
    gsInfo<<"Nodes = "<<Nodes<<"\n";
    gsInfo<<"Weights = "<<Weights<<"\n";

    gsWriteParaviewPoints(Nodes,"Points_Patch3D");

    ///////////////////////////


    return 0;
}