/** @file gsG1AuxiliaryEdgeMultiplePatches.h
 *
    @brief Reparametrize the Geometry for one Interface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

# include <gismo.h>

# include <gsG1Basis/gsG1AuxiliaryPatch.h>
# include <gsG1Basis/gsG1OptionList.h>

// Andrea
# include <gsG1Basis/gsG1ASBasisEdge.h>
# include <gsG1Basis/gsG1ASGluingData.h>

// Pascal
# include <gsG1Basis/ApproxG1Basis/gsApproxG1BasisEdge.h>

//# include <gsG1Basis/ApproxG1Basis/gsApproxBetaSAssembler.h>
//# include <gsG1Basis/gsApproxGluingData2.h>
//# include <gsG1Basis/gsApproxGluingData3.h>
//# include <gsG1Basis/gsApproxGluingData4.h>

namespace gismo
{


class gsG1AuxiliaryEdgeMultiplePatches
{

public:

    // Constructor for one patch and it's boundary
    gsG1AuxiliaryEdgeMultiplePatches(const gsMultiPatch<> & sp, const size_t patchInd){
        auxGeom.push_back(gsG1AuxiliaryPatch(sp.patch(patchInd), patchInd));

    }

    // Constructor for one patch and it's boundary
    gsG1AuxiliaryEdgeMultiplePatches(const gsMultiPatch<> & sp, const gsMultiBasis<> & mb, const size_t patchInd){
        auxGeom.push_back(gsG1AuxiliaryPatch(sp.patch(patchInd), mb.basis(patchInd), patchInd));

    }

    // Constructor for two patches along the common interface
    gsG1AuxiliaryEdgeMultiplePatches(const gsMultiPatch<> & mp, const size_t firstPatch, const size_t secondPatch){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(firstPatch), firstPatch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(secondPatch), secondPatch));

    }

    // Constructor for two patches along the common interface
    gsG1AuxiliaryEdgeMultiplePatches(const gsMultiPatch<> & mp, const gsMultiBasis<> & mb, const size_t firstPatch, const size_t secondPatch){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(firstPatch), mb.basis(firstPatch), firstPatch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(secondPatch), mb.basis(secondPatch), secondPatch));

    }

    real_t getError() { return m_error; }
    index_t getPlus() { return plus; }

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    gsMultiPatch<> computeAuxTopology(){
        gsMultiPatch<> auxTop;
        for(unsigned i = 0; i <  auxGeom.size(); i++)
        {
            if(auxGeom[i].getPatch().orientation() == -1)
            {
                auxGeom[i].swapAxis();
//                gsInfo << "Changed axis on patch: " << auxGeom[i].getGlobalPatchIndex() << "\n";
            }
            auxTop.addPatch(auxGeom[i].getPatch());
        }
        auxTop.computeTopology();
        return auxTop;
    }


    gsMultiPatch<> reparametrizeG1Interface(){
        gsMultiPatch<> repTop(this->computeAuxTopology());

        if(repTop.interfaces()[0].second().side().index() == 1 && repTop.interfaces()[0].first().side().index() == 3)
            return repTop;

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (repTop.interfaces()[0].second().side().index())
        {
            case 1:
//                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << "\tLocal patch: " << repTop.interfaces()[0].second().patch << " not rotated\n";
                break;
            case 4: auxGeom[0].rotateParamClock();
//                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].second().patch << " rotated clockwise\n";
                break;
            case 3: auxGeom[0].rotateParamAntiClock();
//                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].second().patch << " rotated anticlockwise\n";
                break;
            case 2: auxGeom[0].rotateParamAntiClockTwice();
//                gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].second().patch << " rotated twice anticlockwise\n";
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (repTop.interfaces()[0].first().side().index())
        {
            case 3:
//                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " not rotated\n";
                break;
            case 4: auxGeom[1].rotateParamAntiClockTwice();
//                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " rotated twice anticlockwise\n";
                break;
            case 2: auxGeom[1].rotateParamAntiClock();
//                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " rotated anticlockwise\n";
                break;
            case 1: auxGeom[1].rotateParamClock();
//                gsInfo << "Global patch: " << auxGeom[1].getGlobalPatchIndex() <<"\tLocal patch: " << repTop.interfaces()[0].first().patch << " rotated clockwise\n";
                break;
            default:
                break;
        }

       return this->computeAuxTopology();
    }


    gsMultiPatch<> reparametrizeG1Boundary(const int bInd){
        gsMultiPatch<> repTop(this->computeAuxTopology());
        if(auxGeom[0].getOrient())
        {
            switch (bInd)
            {
                case 3:
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " not rotated\n";
                    break;
                case 2:
                    auxGeom[0].rotateParamClock();
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated clockwise\n";
                    break;
                case 4:
                    auxGeom[0].rotateParamAntiClockTwice();
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated twice anticlockwise\n";
                    break;
                case 1:
                    auxGeom[0].rotateParamAntiClock();
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated anticlockwise\n";
                    break;
            }
        }
        else
        {
            switch (bInd)
            {
                case 1:
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " not rotated\n";
                    break;
                case 4:
                    auxGeom[0].rotateParamClock();
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated clockwise\n";
                    break;
                case 2:
                    auxGeom[0].rotateParamAntiClockTwice();
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated twice anticlockwise\n";
                    break;
                case 3:
                    auxGeom[0].rotateParamAntiClock();
//                    gsInfo << "Global patch: " << auxGeom[0].getGlobalPatchIndex() << " rotated anticlockwise\n";
                    break;
            }
        }
        return this->computeAuxTopology();
    }


    void computeG1InterfaceBasis(gsG1OptionList g1OptionList)
    {
        gsMultiPatch<> mp_init;
        mp_init.addPatch(auxGeom[0].getPatch());// Right -> 0 ====> v along the interface
        mp_init.addPatch(auxGeom[1].getPatch()); // Left -> 1 ====> u along the interface

        gsMultiPatch<> test_mp(this->reparametrizeG1Interface()); // auxGeom contains now the reparametrized geometry

        gsMultiBasis<> test_mb(test_mp); // AFTER reparametrizeG1Interface()
        gsMultiPatch<> g1Basis_0, g1Basis_1;




        if(g1OptionList.getInt("user") == user::pascal)
        {

            std::vector<gsBSplineBasis<>> basis_pm;
            gsBSplineBasis<> basis_1 = dynamic_cast<gsBSplineBasis<> &>(test_mb.basis(0).component(1)); // 0 -> v, 1 -> u
            gsBSplineBasis<> basis_2 = dynamic_cast<gsBSplineBasis<> &>(test_mb.basis(1).component(0)); // 0 -> v, 1 -> u

            if (g1OptionList.getSwitch("info"))
            {
                gsInfo << "Basis_1 " << auxGeom[0].getBasis().basis(0) << "\n";
                gsInfo << "Basis_2 " << auxGeom[1].getBasis().basis(0) << "\n";

                gsInfo << "Basis_1 " << basis_1 << "\n";
                gsInfo << "Basis_2 " << basis_2 << "\n";
            }

            index_t p_1 = basis_1.degree();
            index_t p_2 = basis_2.degree();

            index_t p = p_1 > p_2 ? p_2 : p_1; // Minimum degree at the interface

            index_t m_r = g1OptionList.getInt("regularity") > p - 2 ? (p-2) : g1OptionList.getInt("regularity"); // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY


            // first,last,interior,mult_ends,mult_interior
            gsKnotVector<real_t> kv_plus(0,1,0,p+1,p-1-m_r); // p,r+1 //-1 bc r+1
            gsBSplineBasis<> basis_plus(kv_plus);

            if (basis_1.numElements() <= basis_2.numElements()) //
                for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i += basis_1.knots().multiplicityIndex(i))
                    basis_plus.insertKnot(basis_1.knot(i),p-1-m_r);
            else
                for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i += basis_2.knots().multiplicityIndex(i))
                    basis_plus.insertKnot(basis_2.knot(i),p-1-m_r);

            index_t innerKnotMulti = g1OptionList.getInt("innerKnotMulti");
            if (innerKnotMulti > 0 && p-m_r < innerKnotMulti+1)
            {
                basis_plus.insertKnot(0.25, m_r-1);
                basis_plus.insertKnot(0.5, m_r-1);
                basis_plus.insertKnot(0.75, m_r-1);
            }


            plus = basis_plus.size();

            basis_pm.push_back(basis_plus);
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Basis plus : " << basis_plus.knots().asMatrix() << "\n";

            gsKnotVector<real_t> kv_minus(0,1,0,p+1-1,p-1-m_r); // p-1,r //-1 bc p-1
            gsBSplineBasis<> basis_minus(kv_minus);

            if (basis_1.numElements() <= basis_2.numElements()) //
                for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i += basis_1.knots().multiplicityIndex(i))
                    basis_minus.insertKnot(basis_1.knot(i),p-1-m_r);
            else
                for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i += basis_2.knots().multiplicityIndex(i))
                    basis_minus.insertKnot(basis_2.knot(i),p-1-m_r);

            if (innerKnotMulti > 0 && p-m_r < innerKnotMulti+1)
            {
                basis_minus.insertKnot(0.25, m_r-1);
                basis_minus.insertKnot(0.5, m_r-1);
                basis_minus.insertKnot(0.75, m_r-1);
            }

            basis_pm.push_back(basis_minus);
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Basis minus : " << basis_minus.knots().asMatrix() << "\n";




            gsApproxGluingData<real_t> gluingData(test_mp, test_mb, false, g1OptionList); // Need both patches and bases

            if (g1OptionList.getSwitch("info"))
                gluingDataCondition(gluingData.get_alpha_S_tilde(0), gluingData.get_alpha_S_tilde(1), gluingData.get_beta_S_tilde(0), gluingData.get_beta_S_tilde(1));


            gsApproxG1BasisEdge<real_t> g1BasisEdge_0(test_mp.patch(0), auxGeom[0].getBasis().basis(0), basis_pm, gluingData, 1, false, g1OptionList);

            gsApproxG1BasisEdge<real_t> g1BasisEdge_1(test_mp.patch(1), auxGeom[1].getBasis().basis(0), basis_pm, gluingData, 0, false, g1OptionList);


            g1BasisEdge_0.setG1BasisEdge(g1Basis_0);
            g1BasisEdge_1.setG1BasisEdge(g1Basis_1);

            m_error = g1BasisEdge_0.get_error();

            if (g1OptionList.getSwitch("info"))
            {
                for (index_t i = 0; i < 3; i++)
                {
                    gsMatrix<> points(2, 1);
                    points.setZero();

                    gsInfo << "coefs: "
                           << g1Basis_0.patch(i).deriv(points)(0, 0) + g1Basis_1.patch(i).deriv(points)(1, 0) << "\n";
                    gsInfo << "coefs 2: "
                           << g1Basis_0.patch(i).deriv(points)(1, 0) + g1Basis_1.patch(i).deriv(points)(0, 0) << "\n";
                }
/*
                for (index_t i = g1BasisEdge_1.get_plus(); i < g1BasisEdge_1.get_plus()+2; i++)
                {
                    gsMatrix<> points(2, 1);
                    points.setZero();

                    gsInfo << "coefs: "
                           << g1Basis_0.patch(i).deriv(points)(0, 0) + g1Basis_1.patch(i).deriv(points)(1, 0) << "\n";
                    gsInfo << "coefs 2: "
                           << g1Basis_0.patch(i).deriv(points)(1, 0) + g1Basis_1.patch(i).deriv(points)(0, 0) << "\n";
                }

                gsInfo << "coefs mat: "
                       << g1Basis_0.patch(0).coefs().reshape(basis_1.size(), basis_2.size()).block(0, 0, 2, 5) << "\n";
                gsInfo << "coefs mat: "
                       << g1Basis_0.patch(1).coefs().reshape(basis_1.size(), basis_2.size()).block(0, 0, 2, 5) << "\n";
                gsInfo << "coefs mat 2: "
                       << g1Basis_0.patch(basis_plus.size()).coefs().reshape(basis_1.size(), basis_2.size())
                                   .block(0, 0, 2, 5) << "\n";
                gsInfo << "coefs mat 2: "
                       << g1Basis_0.patch(basis_plus.size() + 1).coefs().reshape(basis_1.size(), basis_2.size())
                                   .block(0, 0, 2, 5) << "\n";
*/
            }

            //g1ConditionRep(gluingData.get_alpha_S_tilde(0), gluingData.get_alpha_S_tilde(1), g1Basis_0,  g1Basis_1);
            //gluingDataCondition(g1BasisEdge_0.get_alpha(), g1BasisEdge_1.get_alpha(), g1BasisEdge_0.get_beta(), g1BasisEdge_1.get_beta());
        }
        else
        if(g1OptionList.getInt("user") == user::andrea)
        {
            gsG1ASGluingData<real_t> g1BasisEdge(test_mp, test_mb);
            gsG1ASBasisEdge<real_t> g1BasisEdge_0(test_mp.patch(0), test_mb.basis(0), 1, false, g1OptionList, g1BasisEdge);
            gsG1ASBasisEdge<real_t> g1BasisEdge_1(test_mp.patch(1), test_mb.basis(1), 0, false, g1OptionList, g1BasisEdge);
            g1BasisEdge_0.setG1BasisEdge(g1Basis_0);
            g1BasisEdge_1.setG1BasisEdge(g1Basis_1);
        }


//      Patch 0 -> Right
        auxGeom[0].parametrizeBasisBack(g1Basis_0);
//      Patch 1 -> Left
        auxGeom[1].parametrizeBasisBack(g1Basis_1);

    }


    void computeG1BoundaryBasis(gsG1OptionList g1OptionList, const int boundaryInd)
    {
        gsMultiPatch<> test_mp(this->reparametrizeG1Boundary(boundaryInd));
        gsMultiBasis<> test_mb(test_mp);
        gsMultiPatch<> g1Basis_edge;

        if(g1OptionList.getInt("user") == user::pascal)
        {
            if (g1OptionList.getSwitch("twoPatch"))
            {
                gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getBasis().basis(0).component(1)); // 0 -> u, 1 -> v
                for (index_t j = 0; j < 2; j++) // u
                {
                    for (index_t i = 2; i < basis_edge.size()-2; i++) // v
                    {
                        gsMatrix<> coefs;
                        coefs.setZero(auxGeom[0].getBasis().basis(0).size(),1);

                        coefs(i*(auxGeom[0].getBasis().basis(0).size()/basis_edge.size()) + j,0) = 1;

                        g1Basis_edge.addPatch(auxGeom[0].getBasis().basis(0).makeGeometry(coefs));
                    }
                }
            }
            else
            {
                std::vector<gsBSplineBasis<>> basis_pm;

                gsBSplineBasis<> basis_1 = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getBasis().basis(0).component(1)); // 0 -> v, 1 -> u

                index_t p_1 = basis_1.degree();

                index_t m_r = g1OptionList.getInt("regularity") > p_1 - 2 ? (p_1-2) : g1OptionList.getInt("regularity"); // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY

                // first,last,interior,mult_ends,mult_interior
                gsKnotVector<real_t> kv_plus(0,1,0,p_1+1,p_1-1-m_r); // p,r+1 //-1 bc r+1
                gsBSplineBasis<> basis_plus(kv_plus);


                for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                        basis_plus.insertKnot(basis_1.knot(i),p_1-1-m_r);

                basis_pm.push_back(basis_plus);

                gsKnotVector<real_t> kv_minus(0,1,0,p_1+1-1,p_1-1-m_r); // p-1,r //-1 bc p-1
                gsBSplineBasis<> basis_minus(kv_minus);

                for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                        basis_minus.insertKnot(basis_1.knot(i),p_1-1-m_r);

                basis_pm.push_back(basis_minus);

                gsApproxGluingData<real_t> gluingData(test_mp, auxGeom[0].getBasis(), false, g1OptionList);
                if (g1OptionList.getInt("gluingData") == gluingData::local)
                    gluingData.setLocalGluingData(basis_plus, basis_minus, "edge");
                else if (g1OptionList.getInt("gluingData") == gluingData::global)
                    gluingData.setGlobalGluingData(0,1);


                gsApproxG1BasisEdge<real_t> g1BasisEdge(test_mp, auxGeom[0].getBasis().basis(0), basis_pm, gluingData, 1, true, g1OptionList);
                g1BasisEdge.setG1BasisEdge(g1Basis_edge);
            }
        }
        else
        if(g1OptionList.getInt("user") == user::andrea)
        {
            gsG1ASGluingData<real_t> bdyGD; // Empty constructor creates the sol and solBeta in a suitable way to manage the GD on the boundary
            gsG1ASBasisEdge<real_t> g1BasisEdge(test_mp, test_mb, 1, true, g1OptionList, bdyGD);
            g1BasisEdge.setG1BasisEdge(g1Basis_edge);
        }

        auxGeom[0].parametrizeBasisBack(g1Basis_edge);
    }


    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i)
    {
        return auxGeom[i];
    }


    void gluingDataCondition(gsBSpline<> alpha_0, gsBSpline<> alpha_1, gsBSpline<> beta_0, gsBSpline<> beta_1)
{
    // BETA
    // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getPatch().basis().component(1)); // 0 -> v, 1 -> u
    index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

    gsKnotVector<> kv(0, 1, basis_edge.numElements()-1, 2 * m_p  + 1, 2 * m_p - 1 );
    gsBSplineBasis<> bsp(kv);

    gsMatrix<> greville = bsp.anchors();
    gsMatrix<> uv1, uv0, ev1, ev0;

    const index_t d = 2;
    gsMatrix<> D0(d,d);

    gsGeometry<>::Ptr beta_temp;

    uv0.setZero(2,greville.cols());
    uv0.bottomRows(1) = greville;

    uv1.setZero(2,greville.cols());
    uv1.topRows(1) = greville;

    const gsGeometry<> & P0 = auxGeom[0].getPatch(); // iFace.first().patch = 1
    const gsGeometry<> & P1 = auxGeom[1].getPatch(); // iFace.second().patch = 0
    // ======================================

    // ======== Determine bar{beta} ========
    for(index_t i = 0; i < uv1.cols(); i++)
    {
        P0.jacobian_into(uv0.col(i),ev0);
        P1.jacobian_into(uv1.col(i),ev1);

        D0.col(1) = ev0.col(0); // (DuFL, *)
        D0.col(0) = ev1.col(1); // (*,DuFR)

        uv0(0,i) = D0.determinant();
    }

    beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
    gsBSpline<> beta = dynamic_cast<gsBSpline<> &> (*beta_temp);

    uv0.setZero(2,greville.cols());
    uv0.bottomRows(1) = greville;

    // ======== Determine bar{alpha^(L)} == Patch 0 ========
    for (index_t i = 0; i < uv0.cols(); i++)
    {
        P0.jacobian_into(uv0.col(i), ev0);
        uv0(0, i) = 1 * ev0.determinant();

    }

    beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
    gsBSpline<> alpha0 = dynamic_cast<gsBSpline<> &> (*beta_temp);

    uv0.setZero(2,greville.cols());
    uv0.bottomRows(1) = greville;

    // ======== Determine bar{alpha^(L)} == Patch 0 ========
    for(index_t i = 0; i < uv0.cols(); i++)
    {
        P0.jacobian_into(uv0.col(i),ev0);
        D0 = ev0.col(1);
        real_t D1 = 1/ D0.norm();
        uv0(0,i) = - 1.0 * D1 * D1 * ev0.col(1).transpose() * ev0.col(0);
        //uv(0,i) = - ev.col(1).transpose() * ev.col(0);
    }
    beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
    gsBSpline<> beta0 = dynamic_cast<gsBSpline<> &> (*beta_temp);

    gsWriteParaview(beta0, "beta_exact", 2000);

    index_t p_size = 8;
    gsMatrix<> points(1, p_size);
    points.setRandom();
    points = points.array().abs();

    gsVector<> vec;
    vec.setLinSpaced(p_size,0,1);
    points = vec.transpose();

    gsInfo << "Points 2: " << points << " \n";
    gsInfo << "Beta 2: " << beta.eval(points) << " \n";
    gsInfo << "alpha1 2: " << alpha_1.eval(points) << " \n";
    gsInfo << "alpha0 2: " << alpha_0.eval(points) << " \n";
    gsInfo << "alpha0 formula 2: " << alpha_0.eval(points) - alpha0.eval(points) << " \n";
    gsInfo << "beta0 formula 2: " << beta_0.eval(points) - beta0.eval(points) << " \n";
    gsInfo << "beta_0 2: " << beta_0.eval(points) << " \n";
    gsInfo << "beta_1 2: " << beta_1.eval(points) << " \n";

    gsMatrix<> temp;
    temp = alpha_1.eval(points).cwiseProduct(beta_0.eval(points))
        + alpha_0.eval(points).cwiseProduct(beta_1.eval(points))
        - beta.eval(points);


    gsInfo << "Conditiontest Gluing data: \n" << temp.array().abs().maxCoeff() << "\n\n";


}

    void g1ConditionRep(gsBSpline<> alpha_0, gsBSpline<> alpha_1, gsMultiPatch<> g1Basis_0,  gsMultiPatch<> g1Basis_1)
{
    // BETA
    // first,last,interior,mult_ends,mult_interior,degree
    gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getPatch().basis().component(1)); // 0 -> v, 1 -> u
    index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

    gsKnotVector<> kv(0, 1, basis_edge.numElements()-1, 2 * m_p  + 1, 2 * m_p - 1 );
    gsBSplineBasis<> bsp(kv);

    gsMatrix<> greville = bsp.anchors();
    gsMatrix<> uv1, uv0, ev1, ev0;

    const index_t d = 2;
    gsMatrix<> D0(d,d);

    gsGeometry<>::Ptr beta_temp;

    uv0.setZero(2,greville.cols());
    uv0.bottomRows(1) = greville;

    uv1.setZero(2,greville.cols());
    uv1.topRows(1) = greville;

    const gsGeometry<> & P0 = auxGeom[0].getPatch(); // iFace.first().patch = 1
    const gsGeometry<> & P1 = auxGeom[1].getPatch(); // iFace.second().patch = 0
    // ======================================

    // ======== Determine bar{beta} ========
    for(index_t i = 0; i < uv1.cols(); i++)
    {
        P0.jacobian_into(uv0.col(i),ev0);
        P1.jacobian_into(uv1.col(i),ev1);
        D0.col(1) = ev0.col(0); // (DuFL, *)
        D0.col(0) = ev1.col(1); // (*,DuFR)

        uv0(0,i) = D0.determinant();
    }

    beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
    gsBSpline<> beta = dynamic_cast<gsBSpline<> &> (*beta_temp);



    index_t p_size = 10000;
    gsMatrix<> points(1, p_size);
    points.setRandom();
    points = points.array().abs();

    gsVector<> vec;
    vec.setLinSpaced(p_size,0,1);
    points = vec.transpose();

    gsMatrix<> points2d_0(2, p_size);
    gsMatrix<> points2d_1(2, p_size);

    points2d_0.setZero();
    points2d_1.setZero();
    points2d_0.row(1) = points; // v
    points2d_1.row(0) = points; // u

    real_t g1Error = 0;

    for (size_t i = 0; i < 1; i++)
    {
        gsMatrix<> temp;
        temp = alpha_1.eval(points).cwiseProduct(g1Basis_0.patch(i).deriv(points2d_0).topRows(1))
            + alpha_0.eval(points).cwiseProduct(g1Basis_1.patch(i).deriv(points2d_1).bottomRows(1))
            + beta.eval(points).cwiseProduct(g1Basis_0.patch(i).deriv(points2d_0).bottomRows(1));

        if (temp.array().abs().maxCoeff() > g1Error)
            g1Error = temp.array().abs().maxCoeff();
    }

    gsInfo << "Conditiontest G1 continuity Rep: \n" << g1Error << "\n\n";
}

protected:
    std::vector<gsG1AuxiliaryPatch> auxGeom;
    real_t m_error;

    index_t plus;
};
}

