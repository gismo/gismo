/**
 * @file   gsReadFile.h
 * @author A. Mantzaflaris <Angelos.Mantzaflaris@oeaw.ac.at>
 * @date   March 2013
 * 
 * @brief  Reads a file
 * 
 *  Part of GISMO library, developed at RICAM-Linz, Austria
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
#include <gsIO/gsFileData.h>

namespace gismo {

template< class T>  class gsMultiPatch;

template<class T = real_t >
class gsReadFile
{
public:

    /** 
     * Reads a file into a gsReadFile object
     * 
     * @param fn filename string
     */
    gsReadFile(std::string const & fn)
        { 
            m_data.read(fn);
        };

    ~gsReadFile() { m_data.clear(); };
   
private:
    /// File data as a Gismo xml tree
    gsFileData<T> m_data;

public:

    std::ostream &print(std::ostream &os) const
        { os << "gsReadFile. .\n"; return os; };

    /// Alows to read an Object from a file
    template<class Obj>
    operator Obj * () 
    {
        // Get the first object in the file
        if ( this->m_data.template hasAny< Obj >() )
            return  this->m_data.template getAnyFirst< Obj >();

        gsWarn<< "Failed to read object from file (not found).\n";
        return NULL;
    }

    /// Alows to read a file into a gsGeometry
    operator gsGeometry<T> * () 
    {
        // Get the first geometry in the file
        if ( this->m_data.template hasAny< gsGeometry<T>  >() )
            return  this->m_data.template getAnyFirst< gsGeometry<T>  >();

        gsWarn<< "Failed to read gsGeometry from file (not found).\n";
        return NULL;
    }
    
    /// Alows to read a file into a gsCurve
    operator gsCurve<T> * () 
    {
        // Get the first curve in the file
        if ( this->m_data.template hasAny< gsCurve<T>  >() )
            return  this->m_data.template getAnyFirst< gsCurve<T>  >();

        gsWarn<< "Failed to read gsCurve from file (not found).\n";
        return NULL;
    }
    
    /// Alows to read a file into a gsBasis
    operator gsBasis<T> * () 
    {
        // Get the first basis in the file
        if ( this->m_data.template hasAny< gsBasis<T>  >() )
            return  this->m_data.template getAnyFirst< gsBasis<T> >();

        gsWarn<< "Failed to read gsBasis from file (not found).\n";
        return NULL;
    }
    
    /// Alows to read a file into a gsBasis
    operator gsPlanarDomain<T> * () 
    {
        // Get the first basis in the file
        if ( this->m_data.template hasAny< gsPlanarDomain<T>  >() )
            return  this->m_data.template getAnyFirst< gsPlanarDomain<T> >();

        gsWarn<< "Failed to read gsPlanarDomain from file (not found).\n";
        return NULL;
    }
    
    /// Alows to convert a gsReadFile to a gsMultipatch
    operator gsMultiPatch<T> * () 
    {
        // Get the first MultiPatch tag, if one exists -- TO DO
        if ( this->m_data.template has< gsMultiPatch<T>  >() )
            return  this->m_data.template getFirst< gsMultiPatch<T>  >();
        
        // Else get all geometries and make a multipatch out of that
        if ( this->m_data.template has< gsGeometry<T>  >() )
        {
            std::vector<gsGeometry<T>* > patches = 
                this->m_data.template getAll< gsGeometry<T>  >();
            return new gsMultiPatch<T>( patches );
        }
        
        gsWarn<< "Failed to read gsMultiPatch from file (not found).\n";
        return NULL;
    }
    
    /// Alows to read a gsMesh
    operator gsMesh<T> * () 
    {
        // Get the first Mesh, if one exists
        if ( this->m_data.template has< gsMesh<T>  >() )
            return  this->m_data.template getFirst< gsMesh<T>  >();
        
        gsWarn<< "Failed to read gsMesh from file (not found).\n";
        return NULL;
    }
    
    /// Alows to read a file into a vector of gsBasis
    operator std::vector< gsBasis<T> * > () 
    {
        // Get all bases
        return  this->m_data.template getAll< gsBasis<T> >();
    }
    
    /// Alows to read a PDE
    operator gsPde<T> * () 
    {
        if ( this->m_data.template has< gsPde<T>  >() )
            return  this->m_data.template getFirst< gsPde<T>  >();

        gsWarn<< "Failed to read gsPde from file (not found).\n";
        return NULL;
    }
    
    /// Read a poisson PDE
    operator gsPoissonPde<T> * () 
    {
        if ( this->m_data.template has< gsPoissonPde<T>  >() )
            return  this->m_data.template getFirst< gsPoissonPde<T>  >();

        gsWarn<< "Failed to read gsPoissonPde from file (not found).\n";
        return NULL;
    }
    
   
};  // class gsReadFile
 
};  // namespace gismo
 
