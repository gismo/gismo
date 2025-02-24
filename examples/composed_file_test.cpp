/** @file composed_file_test.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsNurbs/gsMobiusDomain.h>
#include <iomanip>

using namespace gismo;
//! [Include namespace]

void print(const real_t& el);

template <short_t _DIM, class T>
gsTensorBSplineBasis<_DIM,T> integrationBasis(const gsTensorBSplineBasis<_DIM,T> & basis1,
                                           const gsTensorBSplineBasis<_DIM,T> & basis2)
{
    gsTensorBSplineBasis<_DIM,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
    index_t targetDegree;
    for (size_t d = 0; d!=_DIM; d++)
    {
        // 1. Insert interior knots of composition basis
        for (typename gsKnotVector<T>::uiterator it = std::next(basis2.knots(d).ubegin());
                                                    it!= std::prev(basis2.knots(d).uend());
                                                    ++it)
            {
                if (ibasis.knots(d).has(*it))
                    continue;
                ibasis.insertKnot(*it,d);
            }
        // 2. Increase the degree
        targetDegree = ibasis.degree(d) * basis2.degree(d);
        ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);

    }
    return ibasis;
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    //std::string fn = "../filedata/monitor_results/monitor_example_face_nn_r0_e0_R3_E1/composedGeometry.xml";
    real_t qa = 1;
    index_t qb = 1;

    

    std::string exmgeom = "face_nn";
    
    

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addString( "f", "file", "Input XML file", exmgeom );
    cmd.addReal( "a", "qa", "control the number of quadrature points", qa );
    cmd.addInt( "b", "qb", "control the number of quadrature points", qb );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // Column headers
    std::vector<std::string> columns = {"input", "c_e0r0_R0E0", "c_e0r0_R1E0", "c_e0r0_R2E0", "c_e0r0_R3E0", "c_e0r0_R0E1", "c_e0r0_R1E1", "c_e0r0_R2E1", "c_e0r0_R3E1", "c_e0r0_R0E2", "c_e0r0_R1E2", "c_e0r0_R2E2", "c_e0r0_R3E2", "c_e0r0_R0E3", "c_e0r0_R1E3", "c_e0r0_R2E3", "c_e0r0_R3E3"};
    // Row labels
    std::vector<std::string> rows = {exmgeom, "area", "int-detJ", "int-FrobJ", "area-dist", "ang-dist"};

    // Open the file for writing
    
    //index_t RR = 0;
    //index_t EE = 1;
    gsMatrix <real_t> results(rows.size()-1, columns.size());
    results.setZero();
    
    for (index_t EE = 0; EE < 3; EE++){
        for (index_t RR = 0; RR <4; RR++){
            std::string fn = "../filedata/monitor_results/monitor_example_" + exmgeom + "_r0_e0_R" + std::to_string(RR) + "_E" + std::to_string(EE) + "/composedGeometry.xml";    
            
            
            
            //! [Read input file]
            gsFileData<> fd(fn);
            if(!fd.hasAny<gsGeometry<real_t>>())
                continue;
            //GISMO_ASSERT(,"The input file must contain a geometry.");
            gsGeometry<>::uPtr geom = fd.getFirst<gsGeometry<>>();

            // Extract the composition
            GISMO_ASSERT((dynamic_cast<gsComposedGeometry<real_t>*>(geom.get())),"The geometry must be a composed geometry.");
            gsComposedGeometry<>::uPtr cgeom = memory::make_unique(static_cast<gsComposedGeometry<real_t>*>(geom.release()));
            gsComposedBasis<>::uPtr    cbasis= memory::make_unique(cgeom->basis().clone().release());

            GISMO_ASSERT((dynamic_cast<gsGeometry<real_t>*>(&cgeom->composition())),"The composition must be a geometry.");
            gsGeometry<>::uPtr composition = memory::make_unique(static_cast<gsGeometry<real_t>*>(&cgeom->composition())->clone().release());
            // gsGeometry<>::uPtr geometry    = memory::make_unique(cgeom->geometry().clone().release());
            gsBasis<>::uPtr    basis       = memory::make_unique(cbasis->basis().clone().release());
            gsGeometry<>::uPtr geometry    = basis->makeGeometry(cgeom->coefs());

            //gsInfo<<"Composition:\n"<<*composition<<"\n";
            //gsInfo<<"Geometry:\n"<<*geometry<<"\n";
            //gsInfo<<"Basis:\n"<<*basis<<"\n";

            //gsWriteParaview(*geometry,"geometry");
            //gsWriteParaview(*cgeom,"cgeometry");

            GISMO_ASSERT((dynamic_cast<gsTensorBSplineBasis<2>*>(basis.get())),"The basis must be a tensor basis.");
            gsTensorBSplineBasis<2> * tbasis = static_cast<gsTensorBSplineBasis<2>*>(basis.get());
            GISMO_ASSERT((dynamic_cast<gsTensorBSplineBasis<2>*>(&composition->basis())),"The composition must be a tensor basis.");
            gsTensorBSplineBasis<2> * composition_basis = static_cast<gsTensorBSplineBasis<2>*>(&composition->basis());

            //gsInfo<<"Basis: "<<tbasis->knots(0).asMatrix()<<"\n";

            /*
            std::vector<index_t> mult0 = tbasis->knots(0).multiplicities();
            std::vector<index_t> mult1 = tbasis->knots(1).multiplicities();
            gsInfo << "Multiplicities u-direction: ";
            std::for_each(mult0.begin(), mult0.end(), print);
            gsInfo << "\n";
            gsInfo << "Multiplicities v-direction: ";
            std::for_each(mult1.begin(), mult1.end(), print);
            gsInfo << "\n\n";
            */
            gsTensorBSplineBasis<2> ibasis = integrationBasis(*tbasis,*composition_basis);

            gsExprEvaluator<> ev;
            ev.options().setSwitch("SameElement",false);
            ev.options().setReal("quA",qa);
            ev.options().setInt("quB",qb);

            ev.options().setSwitch("SameElement",false);
            ev.options().setReal("quA",qa);
            ev.options().setInt("quB",qb);
            gsMultiBasis<> mb(ibasis);
            ev.setIntegrationElements(mb);
            auto G = ev.getMap(*geometry);
            auto cG = ev.getMap(*cgeom);

            auto GArea = ev.integral(meas(G));
            auto cGArea = ev.integral(meas(cG));

            std::cout << std::setprecision(9) << std::endl;
            gsInfo<<"Area = "<< GArea  <<"\n";
            gsInfo<<"Area = "<< cGArea <<"\n";


            //jacobian determinant for a surface, i.e. the measure
            auto fform = jac(G).tr()*jac(G);
            auto detG = pow(fform.det().val(),0.5);
            auto G_frob = jac(G).norm();
            auto ang_G = pow( (pow(G_frob,2)-2*detG)/(pow(G_frob,2)+2*detG), 2);

            auto cfform = jac(cG).tr()*jac(cG);
            auto detcG = pow(cfform.det().val(),0.5);
            auto cG_frob = jac(cG).norm();
            auto ang_cG = pow( (pow(cG_frob,2)-2*detcG)/(pow(cG_frob,2)+2*detcG), 2);

            gsInfo << " ----------------- G ---------- cG ---------- \n";
            gsInfo << "file: " << fn << "\n";
            gsInfo << "distorion \n";
            
            gsInfo << "Area    :     " << ev.integral(detG*meas(G))/GArea << " ----- "<< ev.integral(detcG*meas(cG))/cGArea << std::scientific << "\n";
            gsInfo << "Angular :     " << ev.integral(ang_G*meas(G))      << " ----- " << ev.integral(ang_cG*meas(cG)) << std::scientific << "\n";

            // Write the data for this column (from all rows)
            
            if( RR == 0 && EE == 0){
                results(0,0) = GArea;
                results(1,0) = ev.integral(detG*meas(G));
                results(2,0) = ev.integral(G_frob*meas(G));
                results(3,0) = ev.integral(detG*meas(G))/GArea;
                results(4,0) = ev.integral(ang_G*meas(G));
                
                results(0,RR+1) = cGArea;
                results(1,RR+1) = ev.integral(detcG*meas(cG));
                results(2,RR+1) = ev.integral(cG_frob*meas(cG));
                results(3,RR+1) = ev.integral(detcG*meas(cG))/cGArea;
                results(4,RR+1) = ev.integral(ang_cG*meas(cG));
            }
            else{
                gsInfo << "RR: " << RR << " EE: " << EE << "\n";
                results(0,3*EE + RR+1) = cGArea;
                results(1,3*EE + RR+1) = ev.integral(detcG*meas(cG));
                results(2,3*EE + RR+1) = ev.integral(cG_frob*meas(cG));
                results(3,3*EE + RR+1) = ev.integral(detcG*meas(cG))/cGArea;
                results(4,3*EE + RR+1) = ev.integral(ang_cG*meas(cG));
            }

        
        }
    }

    gsInfo << results << "\n";
    
    gsWriteCsv("data.csv" , results, columns);

    // Read all lines of the file into a vector of strings
    std::ifstream inFile("data.csv");
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();  // Close the input file after reading

    // Open the file again for writing (this time in overwrite mode)
    std::ofstream outFile;
    outFile.open("../filedata/monitor_results/" + exmgeom + ".csv"); 

    outFile << rows[0] << ","<< lines[0] << "\n";  // Write the header line
    // Write the data for each row
    for (size_t i = 1; i < lines.size(); ++i) {

        outFile << rows[i] << ",";
        std::istringstream lineStream(lines[i]);
        std::string value;
        
        index_t countcols = 0;
        while (std::getline(lineStream, value, ',')) {
            real_t numValue = std::stod(value);  // Convert the string to a double
            if(countcols < columns.size()-1)
            {
                outFile << std::scientific << std::setprecision(6) << numValue <<",";
                countcols = countcols+1;
            }
            else
                outFile << std::scientific << std::setprecision(6) << numValue <<"\n"; // Format in scientific notation
        }
        
    }
    outFile.close(); 

    return EXIT_SUCCESS;

}// end main



void print(const real_t& el)
{
    gsInfo << el << " ";
}
