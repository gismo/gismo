/** @file gsReadFile.h

    @brief Utility class which holds I/O XML data to read/write to/from files

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

/*
  Notes:
  QEPCAD:
  http://opus.bath.ac.uk/29503/
  http://opus.bath.ac.uk/29503/9/QEPCADexamplebank.txt

*/

#pragma once

#include <string>

#include <gsCore/gsDebug.h>

namespace gismo
{

template< class T>  class gsMultiPatch;

/**
  \brief Reads an object from a data file, if such the requested
  object exists in the file

  \sa gsFileData

  \tparam T arithmentic type
 */
template<class T = real_t >
class gsReadFile
{
public:

    /**
        \brief Opens a file and reads an object into a smartpointer
        (uPtr).

        This is used in conjuction with cast operators and depends on
        the type of object.

        For gsBasis there exists a cast opterator to
        std::vector<gsBasis<T>::uPtr>.

        If read fails, value inside smartpointer is null.
        Example of usage:
        \code{.cpp}
        gsFunctionExpr<>::uPtr expr;
        if(expr = gsReadFile<real_t>("/path/to/file.xml"))
            gsInfo << *expr;
        \endcode
        @param fn filename string
    */
    gsReadFile(std::string const & fn)
    : m_id(-1)
    {
        m_data.read(fn);
    }

    /*! @copydoc gsReadFile::gsReadFile(std::string const &)
     *
     * @param id
     */
    gsReadFile(std::string const & fn, index_t id)
    : m_id(id)
    {
        m_data.read(fn);
    }

    /**
        \brief Opens a file and reads an object into \a result.

        Example of usage:
        \code{.cpp}
        std::string filename = "/path/to/file.xml";
        gsMultiPatch<> MP;
        gsReadFile<>(filenane, MP);
        \endcode
        \exception std::runtime_error Error reading file.
       \param[in] fn filename string
       \param[out] result object to read into
    */
    template<class Obj>
    gsReadFile(std::string const & fn, Obj & result)
    : m_id(-1)
    {
        GISMO_ENSURE(m_data.read(fn), "Error reading file.");
        m_data.getAnyFirst(result);
    }

    /**
        \brief Opens a file and reads a gsMultiPatch object into \a result.

        Example of usage:
        \code{.cpp}
        std::string filename = "/path/to/file.xml";
        gsMultiPatch<> MP;
        gsReadFile<>(filenane, MP);
        \endcode
        \exception std::runtime_error Error reading file.
       \param[in] fn filename string
       \param[out] result gsMultiPatch object to read into
    */
    gsReadFile(std::string const & fn, gsMultiPatch<T> & result)
    : m_id(-1)
    {
        GISMO_ENSURE(m_data.read(fn), "Error reading file.");
        result.clear();
        memory::unique_ptr< gsMultiPatch<T> > mp =
            this->operator memory::unique_ptr< gsMultiPatch<T> >();
        if(mp)
            result = give(*mp);
    }

    ~gsReadFile() { m_data.clear(); }

private:

    /// File data as a Gismo xml tree
    gsFileData<T> m_data;

    index_t m_id;

public:

    /// Allows to read an Object from a file
    template<class Obj>
    operator memory::unique_ptr<Obj> ()
    {
        // Get the first object in the file
        if ( this->m_data.template hasAny< Obj >() )
            return  this->m_data.template getAnyFirst< Obj >();

        gsWarn<< "Failed to read object from file (not found).\n";
        return memory::unique_ptr<Obj>();
    }

    /// Allows to read a file into a gsGeometry
    operator memory::unique_ptr< gsGeometry<T> > ()
    {
        // Get the first geometry in the file
        if ( this->m_data.template hasAny< gsGeometry<T>  >() )
            return  this->m_data.template getAnyFirst< gsGeometry<T>  >();

        gsWarn<< "Failed to read gsGeometry from file (not found).\n";
        return memory::unique_ptr< gsGeometry<T> >();
    }

    /// Allows to read a file into a gsCurve
    operator memory::unique_ptr< gsCurve<T> > ()
    {
        // Get the first curve in the file
        if ( this->m_data.template hasAny< gsCurve<T>  >() )
            return  this->m_data.template getAnyFirst< gsCurve<T>  >();

        gsWarn<< "Failed to read gsCurve from file (not found).\n";
        return memory::unique_ptr< gsCurve<T> >();
    }

    /// Allows to read a file into a gsBasis
    operator memory::unique_ptr< gsBasis<T> > ()
    {
        // Get the first basis in the file
        if ( this->m_data.template hasAny< gsBasis<T>  >() )
            return  this->m_data.template getAnyFirst< gsBasis<T> >();

        gsWarn<< "Failed to read gsBasis from file (not found).\n";
        return memory::unique_ptr< gsBasis<T> >();
    }

    /// Allows to read a function expression
    operator memory::unique_ptr< gsFunctionExpr<T> > () const
    {

        if ( this->m_data.template hasAny< gsFunctionExpr<T> >() )
        {
            if ( -1==m_id )        // Get the first one in the file
                return  this->m_data.template getAnyFirst< gsFunctionExpr<T> >();
            else
                return  this->m_data.template getId<gsFunctionExpr<T> >(m_id);
        }

        gsWarn<< "Failed to read gsFunctionExpr from file (not found).\n";
        return memory::unique_ptr< gsFunctionExpr<T> >();
    }

    /// Allows to read a file into a gsBasis
    operator memory::unique_ptr< gsPlanarDomain<T> > ()
    {
        // Get the first basis in the file
        if ( this->m_data.template hasAny< gsPlanarDomain<T>  >() )
            return  this->m_data.template getAnyFirst< gsPlanarDomain<T> >();

        gsWarn<< "Failed to read gsPlanarDomain from file (not found).\n";
        return memory::unique_ptr< gsPlanarDomain<T> >();
    }

    /// Allows to convert a gsReadFile to a gsMultipatch
    operator memory::unique_ptr< gsMultiPatch<T> > ()
    {
        // Get the first MultiPatch tag, if one exists
        if ( this->m_data.template has< gsMultiPatch<T> >() )
            return  this->m_data.template getFirst< gsMultiPatch<T> >();

        // Else get all geometries and make a multipatch out of that
        if ( this->m_data.template has< gsGeometry<T> >() )
        {
            std::vector< memory::unique_ptr<gsGeometry<T> > > patches =
                this->m_data.template getAll< gsGeometry<T> >();
            std::vector< gsGeometry<T>* > releasedPatches = memory::release(patches);
            return memory::make_unique(new gsMultiPatch<T>( releasedPatches ));
        }

        gsWarn<< "Failed to read gsMultiPatch from file (not found).\n";
        return memory::unique_ptr< gsMultiPatch<T> >();
    }

    /// Allows to read a gsMesh
    operator memory::unique_ptr< gsMesh<T> > ()
    {
        // Get the first Mesh, if one exists
        if ( this->m_data.template has< gsMesh<T>  >() )
            return  this->m_data.template getFirst< gsMesh<T>  >();

        gsWarn<< "Failed to read gsMesh from file (not found).\n";
        return memory::unique_ptr< gsMesh<T> >();
    }

    /// Allows to read a file into a vector of gsBasis
    operator std::vector< memory::unique_ptr< gsBasis<T> > > ()
    {
        // Get all bases
        return  this->m_data.template getAll< gsBasis<T> >();
    }

    /// Allows to read a PDE
    operator memory::unique_ptr< gsPde<T> > ()
    {
        if ( this->m_data.template has< gsPde<T>  >() )
            return  this->m_data.template getFirst< gsPde<T>  >();

        gsWarn<< "Failed to read gsPde from file (not found).\n";
        return memory::unique_ptr< gsPde<T> >();
    }

    /// Read a poisson PDE
    operator memory::unique_ptr< gsPoissonPde<T> > ()
    {
        if ( this->m_data.template has< gsPoissonPde<T>  >() )
            return  this->m_data.template getFirst< gsPoissonPde<T>  >();

        gsWarn<< "Failed to read gsPoissonPde from file (not found).\n";
        return memory::unique_ptr< gsPoissonPde<T> >();
    }

    /// Allows to convert a gsReadFile to a sharead pointer
    template<class Obj>
    operator memory::shared_ptr<Obj> ()
    {
        return (memory::unique_ptr<Obj>)(*this);
    }

};  // class gsReadFile

/// \brief Write an arbitrary Gismo object to an XML file with the given filename.
template <typename Object>
void gsWrite(const Object& obj, const std::string& fname)
{
    gsFileData<real_t> fd;
    fd << obj;
    fd.dump(fname);
}

} // namespace gismo
