/** @file gsC1ArgyrisIO.h

    @brief Some Input Output stuff for Argyris space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

namespace gismo
{



class gsC1ArgyrisIO {

    public:
    /// Empty constructor
    ~gsC1ArgyrisIO() { }

    gsC1ArgyrisIO()
    { }

    void checkInput(gsMultiPatch<> & mp,
                            gsOptionList optionList)
    {
        bool inputCheck = true;
        //! [Assumptions for Geometry]
        // Check if the geometry is C^0
        // TODO

        // Check if the geometry is C^2
        for (size_t np = 0; np < mp.nPatches(); ++np)
        {
            gsTensorBSplineBasis<2, real_t> basis_patch = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(mp.basis(np));
            for (index_t dir = 0; dir < mp.domainDim(); dir++)
            {
                std::vector<index_t> mult_0 = basis_patch.knots(dir).multiplicities();
                if (mult_0.size() > 2)
                {
                    mult_0.erase(mult_0.begin()); // First
                    mult_0.pop_back(); // Last

                    index_t deg_0 = basis_patch.degree(dir);

                    for (std::vector<index_t>::const_iterator i = mult_0.begin(); i != mult_0.end(); ++i)
                        if (*i > deg_0 - 2)
                        {
                            inputCheck = false;
                            gsWarn <<"The geometry might(!) not suitable for using the Argyris space! It is not C^2! \n";
                        }
                }
            }
        }

        // Some more checks?
        //! [Assumptions for Geometry]

        //! [Assumptions for spline space]
        index_t discrete_p = optionList.getInt("discreteDegree");
        for (size_t numInt = 0; numInt < mp.nInterfaces(); ++numInt)
        {
            const boundaryInterface & item = mp.interfaces()[numInt];

            gsTensorBSplineBasis<2, real_t> basis_patch_1 = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(mp.basis(item.first().patch));
            gsTensorBSplineBasis<2, real_t> basis_patch_2 = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(mp.basis(item.second().patch));

            const index_t dir_1 = item.first().side().index() > 2 ? 0 : 1;
            const index_t dir_2 = item.first().side().index() > 2 ? 0 : 1;

            index_t deg_1 = basis_patch_1.degree(dir_1);
            index_t deg_2 = basis_patch_2.degree(dir_2);

            if (deg_1 + discrete_p < 3 || deg_2 + discrete_p < 3)
            {
                gsWarn << "The polynomial degree of the discrete space should be > 2! Use \"-p 1\" or \"-p 2\"\n";
                inputCheck = false;
            }
        }

        index_t discrete_r = optionList.getInt("discreteRegularity");
        for (size_t np = 0; np < mp.nPatches(); ++np)
        {
            gsTensorBSplineBasis<2, real_t> basis_patch = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(mp.basis(np));
            for (index_t dir = 0; dir < mp.domainDim(); dir++)
            {
                index_t deg_0 = basis_patch.degree(dir);
                if (deg_0 + discrete_p - 1 < discrete_r || discrete_r < 1)
                {
                    gsWarn << "The regularity of the discrete space should be 1 <= r <= p-1!\n";
                    inputCheck = false;
                }
            }
        }
        //! [Assumptions for spline space]


        if (!inputCheck)
        {
            gsInfo << "\n=================================================================================================";
            gsInfo << "\n\n";
            gsWarn << "Something with the input data does not fulfill the assumption for constructing the Argyris space!\n\n";
            gsInfo << "=================================================================================================\n";
        }
    }

    void gsWriteCSV( const gsMatrix<> & mat, std::string const & filename )
    {
        std::string tmp = gsFileManager::getExtension(filename);
        if (tmp != "csv" )
            tmp = filename + ".csv";
        else
            tmp = filename;

        std::ofstream file(tmp.c_str());
        for ( index_t i = 0 ; i != mat.rows(); ++i )
        {
            file << mat(i,0);
            for ( index_t j = 1 ; j != mat.cols(); ++j )
                file << "," << mat(i,j) ;
            file<< "\n";
        }
        file.close();
    }

    void writeLineString(std::ofstream & file, std::string command, std::string name )
    {
        file<<"# " + command + "\n";
        file<<name;
        file<<"\n";
    }

    void writeBlockMatrix(std::ofstream & file, std::string command, gsMatrix<> matrix, std::vector<std::string> colName , bool rate = false)
    {

        file<<"# Start " + command + "\n";
        // Colname
        for (std::vector<std::string>::const_iterator it = colName.begin(); it != colName.end(); it++)
        {
            if (it != std::prev(colName.end()))
                file<<*it<<',';
            else
                file<<*it;
        }
        file<<"\n";

        // Results
        for(int  i = 0; i < matrix.rows(); i++){
            for(int j = 0; j < matrix.cols(); j++){
                if (rate) {
                    if (j + 1 == matrix.cols())
                        file << std::fixed << std::setprecision(2) << matrix(i, j); // last is rate
                    else if (j%2 == 1 && j > 2)
                        file << std::fixed << std::setprecision(2) << matrix(i, j) << ',';
                    else if (j == 0)
                        file << std::fixed << std::setprecision(5) << matrix(i, j) << ',';
                    else if (j == 1)
                        file << std::fixed << std::setprecision(0) << matrix(i, j) << ',';
                    else
                        file << std::scientific << std::setprecision(5) << matrix(i, j) << ',';
                }
                else {
                    if (j + 1 == matrix.cols())
                        file << std::scientific << std::setprecision(5) << matrix(i, j);
                    else
                        file << std::scientific << std::setprecision(5) << matrix(i, j) << ',';
                }
            }
            file<<'\n';
        }
        file<<"# End " + command + "\n";
    }

    void writeBlockMatrix(std::ofstream & file, std::string command, gsMatrix<> matrix, gsMatrix<> rate, std::vector<std::string> colName)
    {

        file<<"# Start " + command + "\n";
        // Colname
        for (std::vector<std::string>::const_iterator it = colName.begin(); it != colName.end(); it++)
        {
            if (it != std::prev(colName.end()))
                file<<*it<<',';
            else
                file<<*it;
        }
        file<<"\n";

        // Results
        for(int  i = 0; i < matrix.rows(); i++){
            for(int j = 0; j < matrix.cols(); j++){
                file << std::scientific << std::setprecision(5) << matrix(i, j) << ',';
                if (j + 1 == matrix.cols())
                    file << std::fixed << std::setprecision(2) << rate(i, j); // last is rate
                else
                    file << std::fixed << std::setprecision(2) << rate(i, j) << ',';
            }
            file<<'\n';
        }
        file<<"# End " + command + "\n";
    }

    void saveMesh(std::string & path, gsMultiPatch<> & mp, gsMultiBasis<> & mb, unsigned resolution = 100)
    {
        std::string name = "linedata";
        std::string name2 = "points";
        std::string name3 = "points_sorted_";

        std::ofstream file_linedata(path + "/" + name + ".csv");
        std::ofstream file_points(path + "/" + name2 + ".csv");


        size_t offset = 0;
        for (size_t np = 0; np < mp.nPatches(); np++)
        {
            std::ofstream file_points_sorted(path + "/" + name3 + util::to_string(np) + ".csv");

            gsBasis<real_t> &basis = mb.basis(np);
            gsGeometry<real_t> &Geo = mp.patch(np);

            //basis.uniformRefine();

            gsMesh<real_t> sl(basis, resolution);
            Geo.evaluateMesh(sl);

            for (typename std::vector<gsVertex<real_t> *>::const_iterator it = sl.vertices().begin();
                 it != sl.vertices().end(); ++it) {
                const gsVertex<real_t> &vertex = **it;
                file_points << vertex[0] << " ";
                file_points << vertex[1] << "\n";
                //gsInfo << vertex[0] << " "; // 3D
                //gsInfo << vertex[1] << "\n"; // 3D
            }

            for (typename std::vector<gsEdge<real_t> >::const_iterator it = sl.edges().begin();
                 it != sl.edges().end(); ++it) {
                file_linedata << std::fixed << std::setprecision(4) << it->source->getId() + offset << " " << it->target->getId() + offset << "\n";
            }

            offset += sl.numVertices();

            gsMatrix<> points;
            points.setZero(2,resolution);
            gsVector<> point_temp;
            point_temp.setLinSpaced(resolution, 0, 1);

            // v = 0
            points.row(0) = point_temp;
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            // u = 1
            points.setOnes();
            points.row(1) = point_temp;
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            // v = 1
            points.setOnes();
            points.row(0) = point_temp.reverse();
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            // u = 0
            points.setZero();
            points.row(1) = point_temp.reverse();
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            file_points_sorted.close();
        }

        file_linedata.close();
        file_points.close();
    }

    void saveSolution(std::string & path, gsMultiPatch<> & mp, gsMultiBasis<> & mb, gsFunctionExpr<> & solVal, unsigned resolution = 50)
    {
        //std::string name = "solution";
        //std::ofstream file_points(path + "/" + name + ".csv");

        for (size_t np = 0; np < mp.nPatches(); np++)
        {
            std::string name = "solution" + std::to_string(np);
            std::ofstream file_points(path + "/" + name + ".csv");

            gsGeometry<real_t> &Geo = mp.patch(np);

            gsMatrix<real_t> ab = Geo.support();
            gsVector<real_t> a = ab.col(0);
            gsVector<real_t> b = ab.col(1);

            gsVector<unsigned> numpoints = uniformSampleCount(a, b, resolution*resolution);
            gsMatrix<real_t> pts = gsPointGrid(a, b, numpoints);

            gsMatrix<real_t> eval_geo = Geo.eval(pts);//pts
            gsMatrix<real_t> eval_field = solVal.eval(eval_geo);

            for (index_t it = 0; it < eval_geo.cols(); ++it) {
                file_points << eval_geo(0, it) << " ";
                file_points << eval_geo(1, it) << " ";
                file_points << eval_field(0, it) << "\n";
            }

            file_points.close();
        }
        //file_points.close();
    }

}; // class gsC1ArgyrisIO

} // namespace gismo