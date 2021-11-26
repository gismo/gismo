#include <iostream>
#include <gismo.h>

using namespace gismo;

enum CSVFlags
{
    SOLUTION    = 1 << 0, // Print the exact solution
    ERROR       = 1 << 1, // Print the error
    GEOMETRY    = 1 << 2, // Print the geometry
    MESH        = 1 << 3, // Print the mesh
    // Add more [...]
};

inline CSVFlags operator|(CSVFlags a, CSVFlags b)
{
    return static_cast<CSVFlags>(static_cast<int>(a) | static_cast<int>(b));
}

template<short_t d, class T>
class gsCSVOutput {

public:
    gsCSVOutput();

    gsCSVOutput( gsMultiPatch<T> & mp,
                 gsMultiBasis<T> & mb,
                 std::string path)
    : m_mp(mp), m_mb(mb), m_path(path)
    {

    }

public:

    void errorMatrix( gsMatrix<T> & matrix ) { m_matrix = matrix; }
    void jumpMatrix( gsMatrix<T> & matrix ) { m_jump = matrix; }
    void jumpRateMatrix( gsMatrix<T> & matrix ) { m_jumpRate = matrix; }

    /// \brief Save file contents to an csv file
    void saveCSVFile(std::string name = "dump")
    {
        if ((CSVFlags::SOLUTION & flags) == CSVFlags::SOLUTION) {
            gsInfo << "Solution \n";
        }
        if ((CSVFlags::ERROR & flags) == CSVFlags::ERROR) {
            gsInfo << "Error \n";
        }
        if ((CSVFlags::GEOMETRY & flags) == CSVFlags::GEOMETRY) {
            gsInfo << "Geometry \n";
        }
        if ((CSVFlags::MESH & flags) == CSVFlags::MESH) {
            gsInfo << "Mesh \n";
        }
    }

public:
    mutable unsigned flags;

     gsMultiPatch<T> & m_mp;
     gsMultiBasis<T> & m_mb;
     std::string m_path;

    gsMatrix<T> m_matrix, m_jump, m_jumpRate;

private:

    void writeLineString(std::ofstream & file, std::string key, std::string name )
    {
        file<<"# Start " + key + "\n";
        file<<name;
        file<<"\n";
    }

    void writeBlockMatrix(std::ofstream & file, std::string key, gsMatrix<> matrix, int prec = 5)
    {
        file<<"# Start " + key + "\n";
        // Results
        for(int  i = 0; i < matrix.rows(); i++){
            for(int j = 0; j < matrix.cols(); j++){
                if (j + 1 == matrix.cols())
                    file << std::scientific << std::setprecision(prec) << matrix(i, j); // Last
                else
                    file << std::scientific << std::setprecision(prec) << matrix(i, j) << ',';
            }
            file<<'\n';
        }
        file<<"# End " + key + "\n";
    }

    void saveMesh(gsMultiPatch<> & mp, gsMultiBasis<> & mb, unsigned resolution = 100)
    {
        std::string name = "linedata";
        std::string name2 = "points";
        std::string name3 = "points_sorted_";

        std::ofstream file_linedata(m_path + "/" + name + ".csv");
        std::ofstream file_points(m_path + "/" + name2 + ".csv");


        size_t offset = 0;
        for (size_t np = 0; np < mp.nPatches(); np++)
        {
            std::ofstream file_points_sorted(m_path + "/" + name3 + util::to_string(np) + ".csv");

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

    void saveSolution(gsMultiPatch<> & mp, gsMultiBasis<> & mb, gsFunctionExpr<> & solVal, unsigned resolution = 50)
    {


        for (size_t np = 0; np < mp.nPatches(); np++)
        {
            std::string name = "solution" + std::to_string(np);
            std::ofstream file_points(m_path + "/" + name + ".csv");

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
    }

};