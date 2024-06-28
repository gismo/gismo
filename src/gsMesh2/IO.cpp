
#include <gsMesh2/IO.h>

#include <clocale>


//== NAMESPACE ================================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


bool read_mesh(gsSurfMesh& mesh, const std::string& filename)
{
    std::setlocale(LC_NUMERIC, "C");

    // clear mesh before reading from file
    mesh.clear();

    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "off")
    {
        return read_off(mesh, filename);
    }
    else if (ext == "obj")
    {
        return read_obj(mesh, filename);
    }
    else if (ext == "stl")
    {
        return read_stl(mesh, filename);
    }
    else if (ext == "poly")
    {
        return read_poly(mesh, filename);
    }

    // we didn't find a reader module
    return false;
}


//-----------------------------------------------------------------------------


bool write_mesh(const gsSurfMesh& mesh, const std::string& filename)
{
    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);


    // extension determines reader
    if (ext == "off")
    {
        return write_off(mesh, filename);
    }
    else if (ext=="obj")
    {
        return write_obj(mesh, filename);
    }
    else if (ext=="poly")
    {
        return write_poly(mesh, filename);
    }
    else if (ext=="stl")
    {
        return write_stl(mesh, filename);
    }

    // we didn't find a writer module
    return false;
}


//=============================================================================
} // namespace gismo
//=============================================================================
