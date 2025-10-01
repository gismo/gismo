/** @file

    @brief

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>
#include <fstream>

#include <gsG1ACC/gsSurfScheme.h>

using namespace gismo;

std::vector<std::string> get_bezier_point(std::string line) {
    std::stringstream ss(line);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vec(begin, end);
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::string>(ss, "\n"));
    return vec;
}


int main(int argc, char** argv)
{
    index_t r(0);
    std::string fn("C:/Users/jimt1/Documents/Git_Repos/freeform_repo/codes/ffs11/temp/filled_hole.bez");
    std::string fn_patch("");
    bool plot = false, save = false;
    index_t numSamples(500);
    bool plot_mesh = false;
    bool plot_net = false;
    bool check_mesh = true;
    bool mixed_deg = true;
    bool quad = false;
    bool check_mesh_dual = true;

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a Bez file");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("plot", "Plot the results", plot);
    cmd.addSwitch("save", "Save the results in files.", save);
    cmd.addSwitch("element", "Plot the element mesh (when applicable)", plot_mesh);
    cmd.addSwitch("controlNet", "Plot the control net (when applicable)", plot_net);
    cmd.addInt("s", "samples", "Number of samples to use for viewing", numSamples);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    
    std::ifstream file{ fn };
    std::string str;
    std::vector<std::string> lines;
    int Nu, Nv, dim;
    int line_c{ 0 };
    index_t deg1;
    gsMatrix<real_t> coefs;
    int patch_global_index=0,patch_local_index=0;
    std::vector<std::string> vec;
    gsMultiPatch<> mp;

    while (std::getline(file, str)) {

        // Reading degree of bezier patches (consider the same degree in each
        // direction)
        if (line_c == 0) {
            
            Nu = (str[3] >= 'A') ? (str[3] - 'A' + 10) : (str[3] - '0');
            Nv = (str[4] >= 'A') ? (str[4] - 'A' + 10) : (str[4] - '0');
            dim = (str[5] >= 'A') ? (str[5] - 'A' + 10) : (str[5] - '0');
            deg1 = Nu + 1;

        }
        else {
            lines.push_back(str);
        }
        line_c++;
    }
    gsKnotVector<> kv(0, 1, 0, deg1);
    gsTensorBSplineBasis<2> bb(kv, kv);
    mp.resize((line_c - 1) / pow(deg1,2));
    line_c = 0;
    for (std::string str:lines) {
        
        // Reading degree of bezier patches (consider the same degree in each
        // direction)

            if (patch_local_index == 0) {
                coefs.setZero(deg1 * deg1, 3);
                vec = get_bezier_point(str);
                coefs(patch_local_index, 0) = std::stod(vec[0]);
                coefs(patch_local_index, 1) = std::stod(vec[1]);
                coefs(patch_local_index, 2) = std::stod(vec[2]);
                patch_local_index++;
            }
            else {
                vec = get_bezier_point(str);
                coefs(patch_local_index, 0) = std::stod(vec[0]);
                coefs(patch_local_index, 1) = std::stod(vec[1]);
                coefs(patch_local_index, 2) = std::stod(vec[2]);
                patch_local_index++;
                if (patch_local_index == deg1 * deg1) {
                    mp.setPatch(patch_global_index, bb.makeGeometry(coefs));
                    patch_global_index++;
                    patch_local_index = 0; // re-initialize
                }
            }

        
    }





    if (plot)
    {
        gsInfo << "Export Bezier patches to Paraview\n";
        gsWriteParaview(mp, "sf_out", numSamples, plot_mesh, plot_net);
        gsFileManager::open("sf_out.pvd");
    }

    return EXIT_SUCCESS;
}


//=============================================================================
