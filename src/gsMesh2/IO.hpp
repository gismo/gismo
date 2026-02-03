
#include <gsMesh2/IO.h>

#include <clocale>
#include <algorithm>


//== NAMESPACE ================================================================


namespace gismo {


//== IMPLEMENTATION ===========================================================


template <class Scalar>
bool read_mesh(gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    std::setlocale(LC_NUMERIC, "C");

    // clear mesh before reading from file
    mesh.clear();

    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    // extension determines reader
    if (ext == "off")
    {
        return read_off<Scalar>(mesh, filename);
    }
    else if (ext == "obj")
    {
        return read_obj<Scalar>(mesh, filename);
    }
    else if (ext == "stl")
    {
        return read_stl<Scalar>(mesh, filename);
    }
    else if (ext == "poly")
    {
        return read_poly<Scalar>(mesh, filename);
    }

    // we didn't find a reader module
    return false;
}


//-----------------------------------------------------------------------------

template <class Scalar>
bool write_mesh(const gsSurfMesh<Scalar>& mesh, const std::string& filename)
{
    // extract file extension
    std::string::size_type dot(filename.rfind("."));
    if (dot == std::string::npos) return false;
    std::string ext = filename.substr(dot+1, filename.length()-dot-1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);


    // extension determines reader
    if (ext == "off")
    {
        return write_off<Scalar>(mesh, filename);
    }
    else if (ext=="obj")
    {
        return write_obj<Scalar>(mesh, filename);
    }
    else if (ext=="poly")
    {
        return write_poly<Scalar>(mesh, filename);
    }
    else if (ext=="stl")
    {
        return write_stl<Scalar>(mesh, filename);
    }

    // we didn't find a writer module
    return false;
}


//=============================================================================
} // namespace gismo
//=============================================================================
