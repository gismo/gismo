/** @file Monge_Ampere_example.cpp

    @brief Tutorial on how to use expression assembler to solve a non-linear Monge-Ampere equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]


void ProjectionNormalCPoints(gsMultiPatch<>& Psi, gsMultiPatch<> mp){
    // Projection normal of control points (exact geometry)
    int boxMaxNumber = mp.nBoxes();
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        if(!mp.isInterface( patchSide(boxNumber,1) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
        }
        }

        if(!mp.isInterface( patchSide(boxNumber,2) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0] + 1.;
        }
        }

        if(!mp.isInterface( patchSide(boxNumber,3) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
        }
        }
        if(!mp.isInterface( patchSide(boxNumber,4) ) ){
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
        }
        }
        }
};

void CorrectCornersLshape(gsMultiPatch<>& Psi, gsMultiPatch<> mp){
    //....To determine whether a geometry is an L-shape and identify which patch is located at a specific corner then correct the corner to be exact
    for (size_t iIntIndex = 0; iIntIndex < mp.nInterfaces(); ++iIntIndex)
    {
    for (size_t jIntIndex = 0; jIntIndex < mp.nInterfaces(); ++jIntIndex)
    {
        if( (jIntIndex <= iIntIndex) || (iIntIndex == mp.nInterfaces() - 1 && jIntIndex != 0))
            continue;            

        // ... Coonfig 1 L 
        if ( mp.bInterface(iIntIndex).first().patchIndex() == mp.bInterface(jIntIndex).first().patchIndex()) // test if patch has two interfaces
        {
        // gsInfo << "tring to find L-shape test 1" <<"\n";
        int boxNumber  = mp.bInterface(iIntIndex).first().patchIndex(); // patch of the corner
        int boxNumber1 = mp.bInterface(iIntIndex).second().patchIndex();
        int boxNumber2 = mp.bInterface(jIntIndex).second().patchIndex();
        if((mp.bInterface(iIntIndex).first().index() == 2 && mp.bInterface(jIntIndex).first().index() == 4) || (mp.bInterface(iIntIndex).first().index() == 4 && mp.bInterface(jIntIndex).first().index() == 2))   // means inner interfaces are 4 and 2
        {
            // gsInfo << ":\n:\n";
            // gsInfo << ".....\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(2).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).second().index() == 1){
                i_x = Psi.patch(boxNumber1).basis().boundary(1).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(3).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(3).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(1).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if ((mp.bInterface(iIntIndex).first().index() == 2 && mp.bInterface(jIntIndex).first().index() == 3) ||(mp.bInterface(iIntIndex).first().index() == 3 && mp.bInterface(jIntIndex).first().index() == 2))   // means inner interfaces are 3 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << ":  \n:  \n";
            int i_x = Psi.patch(boxNumber).basis().boundary(3).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];

            if(mp.bInterface(iIntIndex).second().index() == 1){
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = Psi.patch(boxNumber2).basis().boundary(4).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(4).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            }
        else if((mp.bInterface(iIntIndex).first().index() == 1 && mp.bInterface(jIntIndex).first().index() == 4) ||(mp.bInterface(iIntIndex).first().index() == 4 && mp.bInterface(jIntIndex).first().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "  :\n  :\n";
            // gsInfo << "..." <<"\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(1).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).second().index() == 2){
                i_x = Psi.patch(boxNumber1).basis().boundary(2).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(2).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if((mp.bInterface(iIntIndex).first().index() == 1 && mp.bInterface(jIntIndex).first().index() == 3) ||(mp.bInterface(iIntIndex).first().index() == 3 && mp.bInterface(jIntIndex).first().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << "  :\n  :\n";
            int i_x = 0;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            if(mp.bInterface(iIntIndex).first().index() == 2){
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            else{
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];                    
            }
            }
        }

        // ... Coonfig 2 
        if (  mp.bInterface(iIntIndex).first().patchIndex() == mp.bInterface(jIntIndex).second().patchIndex()){
        // gsInfo << "tring to find L-shape test 2" <<"\n";
        int boxNumber  = mp.bInterface(iIntIndex).first().patchIndex(); // patch of the corner
        int boxNumber1 = mp.bInterface(iIntIndex).second().patchIndex();
        int boxNumber2 = mp.bInterface(jIntIndex).first().patchIndex();
        if((mp.bInterface(iIntIndex).first().index() == 2 && mp.bInterface(jIntIndex).second().index() == 4) ||(mp.bInterface(iIntIndex).first().index() == 4 && mp.bInterface(jIntIndex).second().index() == 2))   // means inner interfaces are 4 and 2
        {
            // gsInfo << ":\n:\n";
            // gsInfo << ".....\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(2).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).second().index() == 1){
                i_x = Psi.patch(boxNumber1).basis().boundary(1).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(3).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(3).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(1).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if ((mp.bInterface(iIntIndex).first().index() == 2 && mp.bInterface(jIntIndex).second().index() == 3) ||(mp.bInterface(iIntIndex).first().index() == 3 && mp.bInterface(jIntIndex).second().index() == 2))   // means inner interfaces are 3 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << ":  \n:  \n";
            int i_x = Psi.patch(boxNumber).basis().boundary(3).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];

            if(mp.bInterface(iIntIndex).second().index() == 1){
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = Psi.patch(boxNumber2).basis().boundary(4).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(4).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            }
        else if((mp.bInterface(iIntIndex).first().index() == 1 && mp.bInterface(jIntIndex).second().index() == 4) ||(mp.bInterface(iIntIndex).first().index() == 4 && mp.bInterface(jIntIndex).second().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "  :\n  :\n";
            // gsInfo << "..." <<"\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(1).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).second().index() == 2){
                i_x = Psi.patch(boxNumber1).basis().boundary(2).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(2).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if((mp.bInterface(iIntIndex).first().index() == 1 && mp.bInterface(jIntIndex).second().index() == 3) ||(mp.bInterface(iIntIndex).first().index() == 3 && mp.bInterface(jIntIndex).second().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << "  :\n  :\n";
            int i_x = 0;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            if(mp.bInterface(iIntIndex).first().index() == 2){
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            else{
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];                    
            }
            }
        }

        // ... Coonfig 3
        if (  mp.bInterface(iIntIndex).second().patchIndex() == mp.bInterface(jIntIndex).first().patchIndex()){
        // gsInfo << "tring to find L-shape test 3" <<"\n";

        int boxNumber  = mp.bInterface(iIntIndex).second().patchIndex(); // patch of the corner
        int boxNumber1 = mp.bInterface(iIntIndex).first().patchIndex();
        int boxNumber2 = mp.bInterface(jIntIndex).second().patchIndex();
        if((mp.bInterface(iIntIndex).second().index() == 2 && mp.bInterface(jIntIndex).first().index() == 4) ||(mp.bInterface(iIntIndex).second().index() == 4 && mp.bInterface(jIntIndex).first().index() == 2))   // means inner interfaces are 4 and 2
        {
            // gsInfo << ":\n:\n";
            // gsInfo << ".....\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(2).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).first().index() == 1){
                i_x = Psi.patch(boxNumber1).basis().boundary(1).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(3).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(3).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(1).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if ((mp.bInterface(iIntIndex).second().index() == 2 && mp.bInterface(jIntIndex).first().index() == 3) ||(mp.bInterface(iIntIndex).second().index() == 3 && mp.bInterface(jIntIndex).first().index() == 2))   // means inner interfaces are 3 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << ":  \n:  \n";
            int i_x = Psi.patch(boxNumber).basis().boundary(3).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];

            if(mp.bInterface(iIntIndex).first().index() == 1){
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = Psi.patch(boxNumber2).basis().boundary(4).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(4).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            }
        else if((mp.bInterface(iIntIndex).second().index() == 1 && mp.bInterface(jIntIndex).first().index() == 4) ||(mp.bInterface(iIntIndex).second().index() == 4 && mp.bInterface(jIntIndex).first().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "  :\n  :\n";
            // gsInfo << "..." <<"\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(1).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).first().index() == 2){
                i_x = Psi.patch(boxNumber1).basis().boundary(2).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(2).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if((mp.bInterface(iIntIndex).second().index() == 1 && mp.bInterface(jIntIndex).first().index() == 3) ||(mp.bInterface(iIntIndex).second().index() == 3 && mp.bInterface(jIntIndex).first().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << "  :\n  :\n";
            int i_x = 0;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            if(mp.bInterface(iIntIndex).first().index() == 2){
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            else{
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];                    
            }
            }                
        }

        // ... Coonfig 4
        if ( mp.bInterface(iIntIndex).second().patchIndex() == mp.bInterface(jIntIndex).second().patchIndex()){
        gsInfo << "tring to find L-shape test 4" <<"\n";

        int boxNumber  = mp.bInterface(iIntIndex).second().patchIndex(); // patch of the corner
        int boxNumber1 = mp.bInterface(iIntIndex).first().patchIndex();
        int boxNumber2 = mp.bInterface(jIntIndex).first().patchIndex();
        if((mp.bInterface(iIntIndex).second().index() == 2 && mp.bInterface(jIntIndex).second().index() == 4) ||(mp.bInterface(iIntIndex).second().index() == 4 && mp.bInterface(jIntIndex).second().index() == 2))   // means inner interfaces are 4 and 2
        {
            // gsInfo << ":\n:\n";
            // gsInfo << ".....\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(2).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).first().index() == 1){
                i_x = Psi.patch(boxNumber1).basis().boundary(1).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(3).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(3).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(1).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if ((mp.bInterface(iIntIndex).second().index() == 2 && mp.bInterface(jIntIndex).second().index() == 3) ||(mp.bInterface(iIntIndex).second().index() == 3 && mp.bInterface(jIntIndex).second().index() == 2))   // means inner interfaces are 3 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << ":  \n:  \n";
            int i_x = Psi.patch(boxNumber).basis().boundary(3).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];

            if(mp.bInterface(iIntIndex).first().index() == 1){
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = Psi.patch(boxNumber2).basis().boundary(4).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            } 
            else{
                i_x = Psi.patch(boxNumber1).basis().boundary(4).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]+1.;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            }
        else if((mp.bInterface(iIntIndex).second().index() == 1 && mp.bInterface(jIntIndex).second().index() == 4) ||(mp.bInterface(iIntIndex).second().index() == 4 && mp.bInterface(jIntIndex).second().index() == 1))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "  :\n  :\n";
            // gsInfo << "..." <<"\n";
            int i_x = Psi.patch(boxNumber).basis().boundary(1).size() -1;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;

            if(mp.bInterface(iIntIndex).first().index() == 2){
                i_x = Psi.patch(boxNumber1).basis().boundary(2).size() -1;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = 0;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            } 
            else{
                i_x = 0;
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(3).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
                i_x = Psi.patch(boxNumber2).basis().boundary(2).size() -1;
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
                Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1]+1.;
            }
            }
        else if((mp.bInterface(iIntIndex).second().index() == 1 && mp.bInterface(jIntIndex).second().index() == 3) ||(mp.bInterface(iIntIndex).second().index() == 1 && mp.bInterface(jIntIndex).second().index() == 3))   // means inner interfaces are 4 and 2
        {
            // gsInfo << "..." <<"\n";
            // gsInfo << "  :\n  :\n";
            int i_x = 0;
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            if(mp.bInterface(iIntIndex).first().index() == 2){
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            }
            else{
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber1).coef( Psi.patch(boxNumber1).basis().boundary(4).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[0] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[0];
            Psi.patch(boxNumber2).coef( Psi.patch(boxNumber2).basis().boundary(2).at(i_x) ).array()[1] = mp.patch(boxNumber).coef( mp.patch(boxNumber).basis().boundary(1).at(0) ).array()[1];                    
            }
            }                
        }
    }
    }
};

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 4;
    index_t numElevate = 0;
    index_t maxIter = 50;
    double eps{1e-5}; // pinalization coefficient
    double l2errRes{0.}, tolerancePicard{1e-8};
    bool last = false, export_b64{false}, adaptiveMesh{true};
    // ...PNormalCP: Correct the normal part of the mapping and CornersLshape: adjust the corners of the three patches that form L.
    bool PNormalCP{false}, CornersLshape{false};
    if(adaptiveMesh){
        PNormalCP = true; 
        CornersLshape =true;
    }

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    // .... one single patch
   gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,2,1, 0.0, 0.0);
   //... patch 2 (L-shape)
   mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,0.0));
   mp.addInterface(0,2,2,1);

//    // ... patch 0-1
//    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,2,1, -1.0, -1.0);
//    // ... patch 2
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1,-1.,-2.0));
//    mp.addInterface(0,3,2,4);
//    // ... patch 3
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1,0.,-2.0));
//    mp.addInterface(2,2,3,1);
//    // ... patch 4
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,-2.0));
//    mp.addInterface(3,2,4,1);
//    // ... patch 5
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1.,-1.0));
//    mp.addInterface(4,4,5,3);
//    // ... patch 6
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1., 0.0));
//    mp.addInterface(5,4,6,3);
//    // ... patch 7
//    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 0.0, 0.0));
//    mp.addInterface(6,1,7,2);
//    mp.addInterface(1,2,7,1);
   // Get all interfaces and boundaries:
   mp.computeTopology();

    //..... Test 2
    // // Manufactured solition
    // gsFunctionExpr<> s("exp(0.5*(x**2 + y**2))",2);
    // // Manufactured Grad solition
    // gsFunctionExpr<> sN("x*exp(0.5*(x**2 + y**2))","y*exp(0.5*(x**2 + y**2))",2);
    // // Right-hand side function
    // gsFunctionExpr<> f("(1.+x**2+y**2)*exp(x**2 + y**2)",2);

    //..... Test 2
    // convex function
    gsFunctionExpr<> s("0.5*(x**2 + y**2)",2);
    // // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    gsFunctionExpr<> f("(1./(2.+cos(8.*pi*sqrt((x-0.5-0.25*0.)**2+(y-0.5)**2))))",2);
    //
    //gsFunctionExpr<> f("(1.+ 9./(1.+(10.*sqrt((x-0.7-0.25*0.)**2+(y-0.5)**2)*cos(atan2(y-0.5,x-0.7-0.25*0.) -20.*((x-0.7-0.25*0.)**2+(y-0.5)**2)))**2) )",2);
    //gsFunctionExpr<> f("( 1.+ 5.*exp(-50.*abs((x-0.5-0.25*cos(2.*pi*0.25))**2-(y-0.5-0.5 *sin(2.*pi*0.25))**2- 0.01)))",2);
    //gsFunctionExpr<> f("(1. + 5./cosh( 5.*((x-sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2 + 5./cosh( 5.*((x+sqrt(3)/2)**2+(y-0.5)**2 - (pi/2)**2) )**2)",2);
    gsInfo<<"Source function "<< f << "\n";

    gsInfo<<"The domain is "<< mp.detail() << "\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
   // For simplicity, set Neumann boundary conditions
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc.addCondition( *bit, condition_type::neumann, &sN );
   }
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    dbasis.degreeElevate(2);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(s, G);

    //gsFunctionExpr<> sI("0.5*(x**2+y**2)+x*y",2);
    //auto u_I = ev.getVariable(sI, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<>  h1err(numRefine+1); //l2err(numRefine+1) : The solution exists up to an additive constant.
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;
    if(adaptiveMesh)
    {
        //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
        for (int r=0; r<=numRefine; ++r)
        {
            dbasis.uniformRefine();
            mp.uniformRefine();
            if( last && r != numRefine)
                continue;
            u.setup(bc, dirichlet::l2Projection, 0);
            // Compute the system matrix and right-hand side

            //... nromalisation of density function
            auto CoeffDensity{ev.integral(ff.val() * meas(G))};
            // Initialize the system : start Computing the conductivity coeffeicient ...
            // Compute the Neumann terms defined on physical space
            auto g_N = A.getBdrFunction(G);
            auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
            // ...
            auto CoeffConductivity{Neumann_Int/ev.integral(pow(2.+2. * CoeffDensity/ff.val(), 0.5) * meas(G))};
            //... end 

            // Initialize the system :  identity mapping as initial guess
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< A.numDofs() <<std::flush;

            timer.restart();

            A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
            ,
            u*  CoeffConductivity * (-1.)*pow(2.+2. * CoeffDensity/ff.val(), 0.5) * meas(G) //rhs vector
            );
            
            // Compute the Neumann terms defined on physical space
            //auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

            ma_time += timer.stop();

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );

            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solver.compute( A.matrix() );
            solVector = solver.solve(A.rhs());

            slv_time += timer.stop();

            gsInfo<< "." <<std::flush; // Linear solving done

            // Picard loop
            index_t NiterPicard{0};
            gsMatrix<> sv0; //
            for(int ip{0}; ip<maxIter; ++ip)
            {
                gsMultiPatch<> UU;
                u_sol.extract(UU);
                gsWrite(UU, "U_solution");
                auto u_s = A.getCoeff(UU);

                //gsMultiBasis<> gbasis(dbasis);
                //gbasis.reduceContinuity(1);
                space v = A.getSpace(dbasis);
                gsMatrix<> vsolVector;
                solution v_sol = A.getSolution(v, vsolVector);
                A.initSystem(2);

                //gsVector<> pt(2); pt.setConstant(0.5);
                //ev.testEval( v, pt );
                //ev.testEval( igrad(u_sol,G), pt );

                // Obtain control points for the gradient of Psi
                A.assemble( v * v.tr() , v * igrad(u_s,G) );
                vsolVector = solver.compute(A.matrix()).solve(A.rhs());
                
                gsMultiPatch<> Psi;
                v_sol.extract(Psi);

                // ... correct boundary
                //if (PNormalCP)
                //    ProjectionNormalCPoints(Psi, mp);
                //if (CornersLshape)
                //    CorrectCornersLshape(Psi, mp);
                
                geometryMap PP = A.getMap(Psi);
                auto fp = A.getCoeff(f,PP);

                // ...  0  dirichlet for boundaries
                sv0 = solVector;
                u.setup(bc, dirichlet::l2Projection, 0);
            
                solution u_sol = A.getSolution(u, solVector);

                // Initialize the system
                A.initSystem();
                setup_time += timer.stop();

                //gsInfo<< A.numDofs() <<std::flush;

                timer.restart();
                // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
                
                // .. update Coeffeicient of conductivity
                CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(CoeffDensity/fp.val() - ihess(u_sol,G).det()), 0.5) * meas(G));

                // MAE system
                A.assemble(
                igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()* meas(G) //matrix
                ,
                u * CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(CoeffDensity/fp.val() - ihess(u_sol,G).det()), 0.5) * meas(G) //rhs vector
                );

                // Compute the Neumann terms defined on physical space
                auto g_N = A.getBdrFunction(G);
                A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose()   );
                

                gsInfo<< " ." <<std::flush;// Assemblying done

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                
                slv_time += timer.stop();

                gsInfo<< "." <<std::flush; // Linear solving done

                // omp_set_dynamic(0);     // Explicitly disable dynamic teams
                // omp_set_num_threads(1); // Use these threads for later parallel regions

                ++NiterPicard;
                l2errRes = (solVector-sv0).norm();// TODO
                if ( l2errRes < tolerancePicard ) break; // TODO
            }//for loop
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << NiterPicard << " L2 residual : "<<std::scientific<<l2errRes<<"\n";
            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(1); // Use these threads for later parallel regions

            timer.restart();
            //l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );

            h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
            err_time += timer.stop();
            gsInfo<< ". " <<std::flush; // Error computations done

        } //for loop
        //! [Solver loop]

    }
    else
    {
        //::::::::::::::::::::     manufactured exact solution         :::::::::::::::::::::::::
        for (int r=0; r<=numRefine; ++r)
        {
            dbasis.uniformRefine();
            mp.uniformRefine();
            if( last && r != numRefine)
                continue;
            u.setup(bc, dirichlet::l2Projection, 0);
            // Compute the system matrix and right-hand side

            // Initialize the system : start Computing the conductivity coeffeicient ...
            // Compute the Neumann terms defined on physical space
            auto g_N = A.getBdrFunction(G);
            auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
            // ...
            auto CoeffConductivity{Neumann_Int/ev.integral(pow(2.+2. * ff.val(), 0.5) * meas(G))};
            //... end 

            // Initialize the system
            A.initSystem();
            setup_time += timer.stop();

            gsInfo<< A.numDofs() <<std::flush;

            timer.restart();

            A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
            ,
            u*  CoeffConductivity * (-1.)*pow(2.+2. * ff.val(), 0.5) * meas(G) //rhs vector
            );
            
            // Compute the Neumann terms defined on physical space
            //auto g_N = A.getBdrFunction(G);
            A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

            ma_time += timer.stop();

            // gsDebugVar(A.matrix().toDense());
            // gsDebugVar(A.rhs().transpose()   );

            gsInfo<< "." <<std::flush;// Assemblying done

            timer.restart();
            solver.compute( A.matrix() );
            solVector = solver.solve(A.rhs());

            slv_time += timer.stop();

            gsInfo<< "." <<std::flush; // Linear solving done

            // Picard loop
            index_t NiterPicard{0};
            gsMatrix<> sv0; //
            for(int ip{0}; ip<maxIter; ++ip)
            {
                sv0 = solVector;
        //        u.setup(bc, dirichlet::interpolation, 0);
                u.setup(bc, dirichlet::l2Projection, 0);
                
                solution u_sol = A.getSolution(u, solVector);

                // Initialize the system
                A.initSystem();
                setup_time += timer.stop();

                //gsInfo<< A.numDofs() <<std::flush;

                timer.restart();
                // Compute the system matrix and right-hand side

                // .. update Coeffeicient of conductivity
                CoeffConductivity = Neumann_Int/ev.integral(pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(ff.val() - ihess(u_sol,G).det()), 0.5) * meas(G));

                // MAE system
                A.assemble(
                igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()  * meas(G) //matrix
                ,
                u * CoeffConductivity * (-1.) * pow( (ilapl(u_sol,G)*ilapl(u_sol,G).tr()).val() + 2.*(ff.val() - ihess(u_sol,G).det()), 0.5) *meas(G) //rhs vector
                );

                // Compute the Neumann terms defined on physical space
                auto g_N = A.getBdrFunction(G);
                A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );

                ma_time += timer.stop();

                // gsDebugVar(A.matrix().toDense());
                // gsDebugVar(A.rhs().transpose()   );
                

                gsInfo<< " ." <<std::flush;// Assemblying done

                timer.restart();
                solver.compute( A.matrix() );
                solVector = solver.solve(A.rhs());
                slv_time += timer.stop();

                gsInfo<< "." <<std::flush; // Linear solving done

                // omp_set_dynamic(0);     // Explicitly disable dynamic teams
                // omp_set_num_threads(1); // Use these threads for later parallel regions

                ++NiterPicard;
                l2errRes = (solVector-sv0).norm();// TODO
                if ( l2errRes < tolerancePicard ) break; // TODO
            }//for loop
        // ! end Picard loop
        gsInfo<< "\n Niter in Picard : " << NiterPicard << " L2 residual : "<<std::scientific<<l2errRes<<"\n";
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        timer.restart();
        //l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );

        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
        } //for loop
        //! [Solver loop]

    } // end of solver
    


    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";

    //! [Error and convergence rates]
    //gsInfo<< "\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (!last && numRefine>0)
    {
        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(igrad(u_sol,G),"gradient_numerical solution");
        if(adaptiveMesh){
        collection.addField(ff, "density function");
        collection.addField(ihess(u_sol,G).det(), "Jacobian function");}
        else{
        collection.addField(u_ex, "exact solution");
        }
        collection.saveTimeStep();
        collection.save();


        gsFileManager::open("ParaviewOutput/solution.pvd");


        gsMultiPatch<> UU;
        u_sol.extract(UU);
        gsWrite(UU, "U_solution");
        auto u_s = A.getCoeff(UU);

        gsMultiBasis<> gbasis(dbasis);
        gbasis.reduceContinuity(1);
        space v = A.getSpace(gbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(2);

        //gsVector<> pt(2); pt.setConstant(0.5);
        //ev.testEval( v, pt );
        //ev.testEval( igrad(u_sol,G), pt );

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * igrad(u_s,G) );
        vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        gsMultiPatch<> Psi;
        v_sol.extract(Psi);

        //... correct the boundary
        if (PNormalCP)
            ProjectionNormalCPoints(Psi, mp);
        if (CornersLshape)
            CorrectCornersLshape(Psi, mp);

        //geometryMap PP = A.getMap(Psi);
        //auto fp = A.getCoeff(f,PP);

        gsWrite(Psi, "Psi_mapping");
        gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main
