/** @file gsFileData.h

    @brief Utility class which holds I/O XML data to read/write to/from files

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <iostream>
#include <string>

#include <gsIO/gsXml.h>

namespace gismo
{

/**
   \brief This class represents an XML data tree which can be read
   from or written to a (file) stream

   \ingroup IO
 */
template<class T>
class gsFileData
{
public:
    typedef internal::gsXmlTree         FileData;
    typedef internal::gsXmlNode         gsXmlNode;
    typedef internal::gsXmlAttribute    gsXmlAttribute;
    typedef std::string                 String;

    typedef gsVector3d<T>               Point_3;

public:

    gsFileData();

    /**
     * Initializes a gsFileData object with the contents of a file
     *
     * @param fn filename string
     */
    explicit gsFileData(String const & fn);

    /**
     * Loads the contents of a file into a gsFileData object
     *
     * @param fn filename string
     *
     * Returns true on success, false on failure.
     */
    bool read(String const & fn) ;

    ~gsFileData();

    /// \brief Clear all data
    void clear();

    /// \brief Reports the number of objects which are held in the file data
    int numData() const { return data->numNodes();}

    /// \brief Save file contents to an xml file
    void save(String const & fname = "dump", bool compress = false) const;

    /// \brief Save file contents to compressed xml file
    void saveCompressed(String const & fname = "dump") const;

    /// \brief Dump file contents to an xml file
    void dump(String const & fname = "dump") const;

    void addComment(String const & message);

    // Returns the path of the last file where data was read from
    //
    // If the corresponding file did not exist, the return value is an
    // empty string.
    String lastPath() const { return m_lastPath; }

    /// Set the precision (number of decimals) used for writing floats
    /// to output files. A 32-bit float has a precision of about 8 digits.
    /// A 64-bit double has a precision of about 16.
    void setFloatPrecision(const unsigned k) { data->setFloatPrecision(k); }

    /// Returns the precision (number of decimals) used for writing floats
    /// to output files. 8 digits reflects to a 32-bit float, while 16 reflects
    /// to a 64-bit double.
    unsigned getFloatPrecision() const { return data->getFloatPrecision(); }

private:
    /// File data as an xml tree
    FileData * data;

    // Used to hold parsed data of native gismo XML files
    std::vector<char> m_buffer;

    // Holds the last path that was used in an I/O operation
    mutable String m_lastPath;

protected:

/*
 * File readers
 */

    /// Reads a file with xml extension
    bool readXmlFile( String const & fn );

    /// Reads a file with xml.gz extension
    bool readXmlGzFile( String const & fn );

    /// Reads Gismo's native XML file
    bool readGismoXmlStream(std::istream & is);

    /// Reads Axel file
    bool readAxelFile(String const & fn);
    bool readAxelSurface( gsXmlNode * node );
    bool readAxelCurve  ( gsXmlNode * node );

    /// Reads GeoPDEs txt file
    bool readGeompFile( String const & fn );

    /// Reads GoTools file
    bool readGoToolsFile(String const & fn);

    /// Reads Off mesh file
    bool readOffFile(String const & fn);

    /// Reads STL mesh file
    bool readStlFile(String const & fn);

    /// Reads Wavefront OBJ file
    bool readObjFile(String const & fn);

    /// Reads OpenCascade brep file
    bool readBrepFile(String const & fn);

    /// Reads Iges file
    bool readIgesFile(String const & fn);

    /// Reads X3D file
    bool readX3dFile(String const & fn);

    /// Reads 3DM file
    bool read3dmFile(String const & fn);

    /// Reads parasolid files
    bool readParasolidFile(String const & fn);

    bool readCsvFile(String const & fn);

    // Show the line number where something went wrong
    void ioError(int lineNumber,const String& str);

public:

    // Generic functions to fetch Gismo object
    // template<class Object>
    // inline Object * get( gsXmlNode * node)
    // {
    //     return internal::gsXml<Object>::get(node);// Using gsXmlUtils
    // }

    /// Searches and fetches the Gismo object with a given id
    template<class Object>
    inline memory::unique_ptr<Object> getId( const int & id)  const
    {
        return memory::make_unique( internal::gsXml<Object>::getId( getXmlRoot(), id ) );
    }

    /// Searches and fetches the Gismo object with a given id
    template<class Object>
    inline void getId( const int & id, Object& result)  const
    {
        memory::unique_ptr<Object> obj = getId<Object>(id);
        result = give(*obj);
    }

    /// Prints the XML tag of a Gismo object
    template<class Object>
    inline String tag() const
    { return internal::gsXml<Object>::tag(); }

    /// Prints the XML tag type of a Gismo object
    template<class Object>
    inline String type() const
    { return internal::gsXml<Object>::type(); }

    /// Returns true if an Object exists in the filedata
    template<class Object>
    inline bool has() const
    {
        return getFirstNode( internal::gsXml<Object>::tag(),
                             internal::gsXml<Object>::type() ) != 0 ;
    }

    /// Returns true if an Object exists in the filedata, even nested
    /// inside other objects
    template<class Object>
    inline bool hasAny() const
    {
        return getAnyFirstNode( internal::gsXml<Object>::tag(),
                                internal::gsXml<Object>::type() ) != 0 ;
    }

    /// Counts the number of Objects in the filedata
    template<class Object>
    inline int count() const
    {
        int i(0);
        for (gsXmlNode * child = getFirstNode( internal::gsXml<Object>::tag(),
                                               internal::gsXml<Object>::type() ) ;
             child; child = getNextSibling(child, internal::gsXml<Object>::tag(),
                                           internal::gsXml<Object>::type() ))
            ++i;
        return i;
    }


    /// Inserts an object to the XML tree
    template<class Object>
    void operator<<(const Object & obj)
    {
        this->add<Object>(obj);
    }

    /// Add the object to the Xml tree, same as <<
    template<class Object>
    void add (const Object & obj)
    {
        gsXmlNode* node =
            internal::gsXml<Object>::put(obj, *data);
        if ( ! node )
        {
            gsInfo<<"gsFileData: Trouble inserting "<<internal::gsXml<Object>::tag()
                         <<" to the XML tree. is \"put\" implemented ??\n";
        }
        else
        {
            data->appendToRoot(node);
        }
    }

    /// Returns the size of the data
    size_t bufferSize() const { return m_buffer.size(); };

    /// Prints the XML data as a string
    std::ostream &print(std::ostream &os) const;

/*
    /// Constructs the first Object found in the XML tree and assigns
    /// it to the pointer obj and then deletes it from the data tree.
    /// Returns true if there was something assigned, false if object
    /// did not exist.
    /// WARNING: Use getFirst<Object>() instead. This is buggy due to
    /// template resolution.
    template<class Object>
    bool operator>>(Object * obj)
    {
        gsWarn<< "getting "<< typeid(Object).name() <<"\n";
        gsXmlNode* node = getFirstNode(internal::gsXml<Object>::tag(),
                                       internal::gsXml<Object>::type() );
        if ( !node )
        {
            gsWarn<<"gsFileData: false!\n";
            return false;
        }
        else
        {
            obj = internal::gsXml<Object>::get(node);
            this->deleteXmlSubtree( node );
            return true;
        }
    }
*/

    /**
     * Returns the first Object found in the XML data as uPtr.
     * Doesn't look for nested objects.
     * Value of uPtr is null if nothing was found.
     * \code{.cpp}
     * gsFunctionExpr<>::uPtr expr;
     * if (expr = getFirst< gsFunctionExpr<> >())
     *     gsInfo << expr;
     * \endcode
     * @tparam Object Type of object.
     * @return An uPtr with the object inside, or null inside if no object was found.
     */
    template<class Object>
    inline memory::unique_ptr<Object> getFirst() const
    {
        gsXmlNode* node = getFirstNode(internal::gsXml<Object>::tag(),
                                       internal::gsXml<Object>::type() );
        if ( !node )
        {
            gsWarn<<"gsFileData: getFirst: Didn't find any "<<
                internal::gsXml<Object>::type()<<" "<<
                internal::gsXml<Object>::tag() <<". Error.\n";
            return memory::unique_ptr<Object>();
        }
        return memory::make_unique( internal::gsXml<Object>::get(node) );// Using gsXmlUtils
    }

    /**
     * Returns the first object of this type found in the XML data.
     * Doesn't look for nested objects.
     * Writes it into the parameter.
     * \code{.cpp}
     * gsMultiPatch<> mp;
     * if(getFirst(mp))
     *     gsInfo << mp;
     * \endcode
     * @tparam Object Type of object.
     * @param result Object read into.
     * @return True if result has been found, false if result was not found.
     */
    template<class Object>
    bool getFirst(Object & result) const
    {
        gsXmlNode* node = getFirstNode(internal::gsXml<Object>::tag(),
                                       internal::gsXml<Object>::type() );
        if ( !node )
        {
            gsWarn<<"gsFileData: getFirst: Didn't find any "<<
                internal::gsXml<Object>::type()<<" "<<
                internal::gsXml<Object>::tag() <<". Error.\n";
            return false;
        }
        internal::gsXml<Object>::get_into(node, result);// Using gsXmlUtils
        return true;
    }

    /// Returns a vector with all Objects found in the XML data
    template<class Object>
    inline std::vector< memory::unique_ptr<Object> > getAll()  const
    {
        std::vector< memory::unique_ptr<Object> > result;

        for (gsXmlNode * child = getFirstNode( internal::gsXml<Object>::tag(),
                                               internal::gsXml<Object>::type() ) ;
             child; child = getNextSibling(child, internal::gsXml<Object>::tag(),
                                           internal::gsXml<Object>::type() ))
        {
            result.push_back( memory::make_unique(internal::gsXml<Object>::get(child)) );
        }
        return result;
    }

    /**
     * Returns the first Object found in the XML data as uPtr.
     * Look also for nested objects.
     * Value of uPtr is null if nothing was found.
     * \code{.cpp}
     * gsFunctionExpr<>::uPtr expr;
     * if (expr = getFirst< gsFunctionExpr<> >())
     *     gsInfo << expr;
     * \endcode
     * @tparam Object Type of object.
     * @return An uPtr with the object inside, or null inside if no object was found.
     */
    template<class Object>
    inline memory::unique_ptr<Object> getAnyFirst() const
    {
        gsXmlNode* node = getAnyFirstNode(internal::gsXml<Object>::tag(),
                                          internal::gsXml<Object>::type() );
        if ( !node )
        {
            gsWarn <<"gsFileData: getAnyFirst: Didn't find any "<<
                internal::gsXml<Object>::type()<<" "<<
                internal::gsXml<Object>::tag() <<". Error.\n";
            return memory::unique_ptr<Object>();
      }
        return memory::make_unique( internal::gsXml<Object>::get(node) );// Using gsXmlUtils
    }

    /**
     * Returns the first object of this type found in the XML data.
     * Look also for nested objects.
     * Writes it into the parameter.
     * \code{.cpp}
     * gsMultiPatch<> mp;
     * if(getFirst(mp))
     *     gsInfo << mp;
     * \endcode
     * @tparam Object Type of object.
     * @param result Object read into.
     * @return True if result has been found, false if result was not found.
     */
    template<class Object>
    bool getAnyFirst(Object & result) const
    {
        gsXmlNode* node = getAnyFirstNode(internal::gsXml<Object>::tag(),
                                          internal::gsXml<Object>::type() );
        if ( !node )
        {
            gsWarn <<"gsFileData: getAnyFirst: Didn't find any "<<
                internal::gsXml<Object>::type()<<" "<<
                internal::gsXml<Object>::tag() <<". Error.\n";
            return false;
        }
        internal::gsXml<Object>::get_into(node, result);// Using gsXmlUtils
        return true;
    }

    /// Lists the contents of the filedata
    String contents () const;

    /// Counts the number of Objects/tags in the filedata
    int numTags () const;

private:

    gsXmlNode * getXmlRoot() const;
    static void deleteXmlSubtree (gsXmlNode* node);

    // getFirst ? (tag and or type)
    gsXmlNode * getFirstNode  ( const String & name = "",
                                const String & type = "" ) const;

    // getAny
    gsXmlNode * getAnyFirstNode( const String & name = "",
                                 const String & type = "" ) const;

    // getNext
    static gsXmlNode * getNextSibling( gsXmlNode* const & node,
                                       const String & name = "",
                                       const String & type = "" );

    // Helpers for X3D files
    void addX3dShape(gsXmlNode * shape);
    void addX3dTransform(gsXmlNode * shape);

}; // class gsFileData

// Print out operator
template<class T>
std::ostream &operator<<(std::ostream &os, const gsFileData<T> & fd)
{return fd.print(os); }

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFileData.hpp)
#endif
