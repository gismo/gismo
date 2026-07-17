
#include <gsMesh2/IO.h>

#include <clocale>


namespace gismo {

bool write_mesh(const gsSurfMesh& mesh, const std::string& filename)
{
    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);


    // extension determines reader
    if (ext == "off")
    {
        return write_off(mesh, filename);
    }
    else if (ext=="obj")
    {
        return write_obj(mesh, filename);
    }
    else if (ext=="stl")
    {
        return write_stl(mesh, filename);
    }

    // we didn't find a writer module
    return false;
}


}
