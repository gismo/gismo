//
// Created by Jingya Li on 18/12/2023.
//

#include <gismo.h>
#include <iostream>
#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsDofMapper.h>
using namespace gismo;

//TODO: write a class for the cable
//namespace gismo{
//    template<class T>
//    class gsCable{
//    public:
//        //Constructor
//        gsCable(gsMultiPatch<T> cable, gsMatrix<T>)
//    protected:
//        gsMatrix<T> localCabelMat;
//        gsMatrix<T> globalMat;
//
//        gsMapData<T> md;
//    protected:
//        //Basis value
//        std::vector<gsMatrix<T> > basisData;
//    };
//}

template<typename T>
gsMatrix<T> evaluateCable(const gsGeometry<T> &geo,
                          const gsMatrix<T> &slaveSamples){
    gsMatrix<T> derivedSlaveSamples;
    geo.eval_into(slaveSamples, derivedSlaveSamples);
    return derivedSlaveSamples;
}

template<typename T>
gsMatrix<T> getCableParameterValue(gsMultiPatch<T> cable, real_t numSamples)
{
    //TODO: add couple directions
    // add multiple cables
    index_t nCables = cable.size();
    if (nCables == 1){
        gsInfo<< "One cable is embedded\n";
        gsMatrix<> slaveSamples(2,numSamples+1);
        gsMatrix<> knotsCable(1,numSamples+1);
        for(index_t i =0; i < numSamples + 1; ++i){
            //Compute the parameter value (uniformly spaced)
            double_t temp = i/(numSamples);
            knotsCable.at(i) =  temp;
        }
        cable.patch(0).eval_into(knotsCable, slaveSamples);
        return slaveSamples;
    }
    else{
        gsInfo<<"Haven't implemented multiple cables yet\n";
        exit(1);
    }
}

//assemble(gsMatrix<T> qpoints,gsMatrix<T> qweights, expression)
//slaveSamples: 2 x numSamples, here are quadrature points
template<typename T>
gsMatrix<T> evaluateCableMassLocal(const gsGeometry<T> &geo,
                   const gsMatrix<T> &slaveSamples,
                   const gsBasis<T> &basis)
{
    gsMapData<T> md;
    // Set Geometry evaluation flags
    md.flags = NEED_MEASURE;
    gsMatrix<index_t> actives;
    //Set the curve quadrature points
    md.points = slaveSamples;

    // Assumes actives are the same for all quadrature points on the elements
    basis.active_into(md.points.col(0), actives);
    gsDebugVar(actives);
    gsMatrix<T> basisData;
    basis.eval_into(md.points, basisData);
    index_t numActive;
    numActive = actives.rows();

//    std::vector<gsMatrix<T> > basisData;
    // Evaluate basis functions on element
    //basis.evalAllDers_into( md.points, 1, basisData);
    basis.eval_into(md.points, basisData);

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geo.computeMap(md);
    gsInfo << "Got here 1\n";

//    gsMatrix<T> & bVals  = basisData[0];
//    gsMatrix<T> & bGrads = basisData[1];
//    gsDebugVar(bGrads.size());
    // Set up the weights, since the cable only follows the movement of the master surface
    // the weights are all 1
    gsVector<T> quWeights(slaveSamples.cols());
    quWeights.setOnes();


    gsMatrix<T> localMat;
    gsMatrix<T>  physGrad;

    localMat.noalias() = basisData * quWeights.asDiagonal() *
            md.measures.asDiagonal() * basisData.transpose();
//    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
//    {
//        // Multiply weight by the geometry measure
//        const T weight = quWeights[k] * md.measure(k);
//        gsDebugVar(weight);
//        gsDebugVar(bGrads.cols());
//        gsDebugVar(bGrads.rows());
//        gsDebugVar(md.measure(k));
//        //JL: should bGrads.cols == bGrads.rows?
//        transformGradients(md, k, bGrads, physGrad);
//        gsInfo << "Got here 3\n";
//        localMat.noalias() += weight * (physGrad.transpose() * physGrad);
//    }
    return localMat;
}

template<typename T>
gsMatrix<T> evaluateCableMassGlobal(const std::vector<gsMatrix<T> > &,
                                    const gsDofMapper & mapper,
                                    const gsGeometry<T> &geo,
                                    const gsMatrix<T> &slaveSamples,
                                    const gsBasis<T> &basis){
    gsMapData<T> md;
    // Set Geometry evaluation flags
    md.flags = NEED_MEASURE;
    gsMatrix<index_t> globalActives;
    //Set the curve quadrature points
    md.points = slaveSamples;

    // Assumes actives are the same for all quadrature points on the elements
    basis.active_into(md.points.col(0), globalActives);
    gsMatrix<T> basisData;
    basis.eval_into(md.points, basisData);
    index_t numActive;
    numActive = globalActives.rows();

    //TODO: add multipatch later
    index_t patchIndex = 0;
    gsVector<T> quWeights(slaveSamples.cols());
    quWeights.setOnes();


    gsMatrix<T> localMat;
    gsMatrix<T>  physGrad;

    localMat.noalias() = basisData * quWeights.asDiagonal() *
                         md.measures.asDiagonal() * basisData.transpose();
    //Get the dof for local matrix

    const index_t localDofs = localMat.rows();
    // Map patch-local DoFs to global DoFs
//    system.mapColIndices(globalActives, patchIndex, globalActives);

}



int main(int argc, char *argv[]){
    bool plot = true;


    //------------------ Generate Geometry ------------------
    /// Master NURBS surface: degree 2, Bezier patch
    gsKnotVector<> Kv(0,2,1,4);
    gsMatrix<> C(25,3);
    C<< 0.124215, 0.0620467, -0.909323,
            0.990116, -0.112269, -0.0161917,
            1.99052, -0.213628, 1.12389,
            2.99999, -0.176554, -0.0225916,
            3.93311, -0.232422, -0.736062,
            -0.0212573, 1.14758, -0.383038,
            1, 1, 0,
            2, 1, 0,
            3, 1, 0,
            4.0773, 1.08316, -0.416151,
            -0.0554789, 2.0257, -0.138956,
            1, 2, 0,
            1.68928, 1.82642, 3.25086,
            3, 2, 0,
            3.94707, 2.11565, -0.357777,
            -0.112924, 2.99649, -0.0987988,
            1, 3, 0,
            2, 3, 0,
            3, 3, 0,
            4.0241, 3.11986, -0.280967,
            -0.209772, 4.06133, -0.303981,
            0.90333, 4.17869, -0.250397,
            1.98218, 4.13487, -0.119872,
            3.10997, 4.14081, -0.364305,
            4.09176, 4.18805, -0.43237;

    gsTensorBSpline<2> master = gsTensorBSpline<2, real_t>( Kv, Kv, give(C) );
    gsTensorBSplineBasis<2> master_basis = master.basis();
    gsInfo << "master_basis = " << master_basis << "\n";

    /// An embedding curve
    // Make a BSpline curve
    gsKnotVector<> slave_knots(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots

    gsMatrix<> coefs(4, 2);
    coefs << 0, 0,
            0.1, 0.2,
            0.5, 0.25,
            0.8, 0.8;

    gsBSpline<> embedding_curve( slave_knots, give(coefs));
    gsMultiPatch<> cables;
    cables.addPatch(embedding_curve);
    gsMatrix<> slaveSamples = getCableParameterValue(cables, 10);
    gsDebugVar(slaveSamples);
    gsMatrix<> localMat = evaluateCableMassLocal(master, slaveSamples, master_basis);
    gsDebugVar(localMat);

    gsMatrix<> derivedSlaveSamples = evaluateCable(master, slaveSamples);
    gsDebugVar(derivedSlaveSamples);
    gsDebugVar(derivedSlaveSamples.dim());
//    gsDebugVar(localMat.rows());
//    gsDebugVar(localMat.cols());
//
//    gsWriteParaview(localMat, "localMat");

//    gsBSplineBasis<> curve_basis = embedding_curve.basis();
//
//
//    gsInfo << "curve_basis = " << curve_basis << "\n";
//
//    // Get the xi extension of the embedded cable (start and end knots)
//    gsInfo<<embedding_curve.knots()<<"\n";
//    gsVector<> xiExtension(2), etaExtension(2);
//    real_t xi_min = embedding_curve.knots().first();
//    real_t xi_max = embedding_curve.knots().last();
//    xiExtension << xi_min,xi_max;
//    gsInfo << "xiExtension = " << xiExtension << "\n";
//    // Get the eta extension of the embedded cable
//    etaExtension<< 0,0;
//    gsMultiPatch<> cables;
//    cables.addPatch(embedding_curve);
//    cables.size();
//
//    gsMatrix<> parameterCable = embedding_curve.coefs(); // The parameterization of the cable
//    // Get the running and the fixed parameters on the patch where the cable is embedded
//    if(xi_min==xi_max)  // If True, coupling is in the eta direction
//    {
//        gsInfo << "Coupling is in the eta direction\n";
//        gsVector<> couplingRegion = etaExtension;
//    }
//    else // If False, coupling is in the xi direction
//    {
//        gsVector<> couplingRegion = xiExtension;
//        // Find the correct spans for the coupled region
//        int spanStart = Kv.iFind(couplingRegion.at(0)) - Kv.begin();
//        int spanEnd = Kv.iFind(couplingRegion.at(1)) - Kv.begin();
//        // Corresponding to the coupled region surface knot span
//        // Construct the basis function for the embedded geometry
//
//
//    }
}


