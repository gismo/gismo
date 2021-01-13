/** @file gsXml.h

    @brief Provides declaration of input/output XML utilities struct.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsExport.h>

// Default memory sizes
// #define RAPIDXML_STATIC_POOL_SIZE  ( 64*1024 )
// #define RAPIDXML_DYNAMIC_POOL_SIZE ( 64*1024 )

#include <rapidxml/rapidxml.hpp>             // External file
#include <rapidxml/rapidxml_print.hpp>       // External file
//#include <rapidxml/rapidxml_utils.hpp>     // External file
//#include <rapidxml/rapidxml_iterators.hpp> // External file


/*
// Forward declare rapidxml structures
namespace rapidxml
{
    template<class Ch> class xml_node;
    template<class Ch> class xml_attribute;
    template<class Ch> class xml_document;
}
*/

#define GSXML_COMMON_FUNCTIONS(obj)             \
    static bool has(gsXmlNode * node)           \
    { return firstByTag(tag(), node) != 0;}     \
    static bool hasAny(gsXmlNode * node)        \
    { return anyByTag(tag(), node) != 0;}       \
    static bool count(gsXmlNode * node)         \
    { return countByTag(tag(), node) != 0; }    \
    static obj * getFirst (gsXmlNode * node)    \
    { return get(firstByTag(tag(), node)); }    \
    static obj * getAny (gsXmlNode * node)      \
    { return get(anyByTag(tag(), node)); }      \
    static  obj * getId (gsXmlNode * node, int id) \
    { return getById< obj >(node, id); }

#define GSXML_GET_POINTER(obj)          \
    static obj * get (gsXmlNode * node) \
    {   obj * result = new obj;         \
        get_into(node, *result);        \
        return result; }

#define TMPLA2(t1,t2)    t1,t2
#define TMPLA3(t1,t2,t3) t1,t2,t3

#ifdef GISMO_WITH_GMP
// Specialize file I/O to floating point format
#include<sstream>
inline std::istringstream &
operator>>(std::istringstream & is, mpq_class & var)
{
    // read as decimal
    std::string dn;
    if ( !(is >> dn) ) return is;
    const std::string::size_type comma( dn.find(".") );
    if( comma != std::string::npos )
    {
        const std::string::size_type exp = dn.size() - comma - 1;
        const mpz_class num( dn.erase(comma,1), 10);
        mpz_class den;
        mpz_ui_pow_ui(den.get_mpz_t(),10,exp);
        var = mpq_class(num, den);
    }
    else // integer or rational
        var.set_str(dn,10);

    //read as machine float
    //double tmp;
    //is >> tmp;
    //var = tmp;

    var.canonicalize();// remove common factors
    return is;
}

#include <fstream>// for paraview
template <class U> inline std::ofstream & operator<<
(std::ofstream &fs, __gmp_expr<U,U> & var)
{
    fs<<var.get_d();
    // write as rational
    //os << var.get_str(10);
    return fs;
}

#endif

namespace gismo {

template<class T>
inline bool gsGetReal(std::istream & is, T & var)
{
    GISMO_STATIC_ASSERT(!std::numeric_limits<T>::is_integer,
        "The second parameter needs to be an integer type.");
    std::string dn;
    if ( !(is >> dn) ) return false;
    const std::string::size_type slh( dn.find("/") );
    if( slh != std::string::npos )
    {
        var = strtod(dn.substr(0,slh).c_str(), NULL) /
              strtod(dn.substr(slh+1).c_str(), NULL) ;
    }
    else // integer or decimal
        var = strtod(dn.c_str(), NULL);

    return true;
}

template<class Z>
inline bool gsGetInt(std::istream & is, Z & var)
{
  GISMO_STATIC_ASSERT(std::numeric_limits<Z>::is_integer,
        "The second parameter needs to be an integer type.");
  //return static_cast<bool>(is >> var); //C++11
  return !(is >> var).fail();
}

#ifdef GISMO_WITH_GMP
template<>
inline bool gsGetReal(std::istream & is, mpq_class & var)
{
    // read as decimal
    bool ok = true;
    std::string dn;
    if ( !(is >> dn) ) return false;
    const std::string::size_type comma( dn.find(".") );
    if( comma != std::string::npos )
    {
        const std::string::size_type exp = dn.size() - comma - 1;
        const mpz_class num( dn.erase(comma,1), 10);// will throw on error
        mpz_class den;
        mpz_ui_pow_ui(den.get_mpz_t(),10,exp);
        var = mpq_class(num, den);
    }
    else // integer or rational
    {
        if ('+'==dn[0]) dn.erase(0, 1);
        ok = (0==var.set_str(dn,10));
    }

    // read as machine float
    //double tmp;
    //is >> tmp;
    //var = tmp;

    var.canonicalize();// remove common factors
    return ok;
}
#endif

//Note: automatic deduction of number traits, however using gsGetReal,
//gsGetInt can reveal type mistakes, so they are preferable
template <typename Z>
typename util::enable_if<std::numeric_limits<Z>::is_integer, bool>::type
gsGetValue(std::istream & is, Z & var)
{ return gsGetInt<Z>(is,var); }

template <typename T>
typename util::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
gsGetValue(std::istream & is, T & var)
{ return gsGetReal<T>(is,var); }

namespace internal {

typedef rapidxml::xml_node<char>        gsXmlNode;
typedef rapidxml::xml_attribute<char>   gsXmlAttribute;
typedef rapidxml::xml_document<char>    gsXmlTree;


/// Generic get XML class: specializations provide implementation
template<class Object>
class gsXml
{
private:
    gsXml() { }// Disallow instantization
public:

    static std::string tag ();
/*    {   // Next line will produce compile-time error
        // when class is not specialized for Object
        Object::Object_does_not_exist_ERROR;
        return "";
    }
//*/
    static std::string type ();
    static Object * get      (gsXmlNode * node);
    static void     get_into (gsXmlNode * node, Object & result);
    static gsXmlNode * put   (const Object & obj, gsXmlTree & data);

    // Common operations
    static bool     has      (gsXmlNode * node);
    static bool     count    (gsXmlNode * node);
    static Object * getFirst (gsXmlNode * node);
    //static void     getFirst_into (gsXmlNode * node);
    static Object * getAny   (gsXmlNode * node);
    //static void     getAny_into   (gsXmlNode * node);
    static Object * getId    (gsXmlNode * node, int id);
    //static void     getId_into   (gsXmlNode * node, int id, Object & result);
};

/// Helper to read an object by a given \em id value:
/// \param node parent node, we check his children to get the given \em id
/// \param id
template<class Object>
Object * getById(gsXmlNode * node, const int & id)
{
    std::string tag = internal::gsXml<Object>::tag();
    for (gsXmlNode * child = node->first_node(tag.c_str()); //note: gsXmlNode object in use
         child; child = child->next_sibling(tag.c_str()))
    {
        const gsXmlAttribute * id_at = child->first_attribute("id");
        if (id_at && atoi(id_at->value()) == id )
            return internal::gsXml<Object>::get(child);
    }
    std::cerr<<"gsXmlUtils Warning: "<< internal::gsXml<Object>::tag()
             <<" with id="<<id<<" not found.\n";
    return NULL;
}

/// Helper to fetch a node with a certain \em id value.
/// \param root parent node, we check his children for the given \em id
/// \param id the ID number which is seeked for
inline gsXmlNode * searchId(const int id, gsXmlNode * root)
{
    for (gsXmlNode * child = root->first_node();
         child; child = child->next_sibling())
    {
        const gsXmlAttribute * id_at = child->first_attribute("id");
        if ( id_at &&  atoi(id_at->value()) == id )
            return child;
    }
    gsWarn <<"gsXmlUtils: No object with id = "<<id<<" found.\n";
    return NULL;
}

/// Helper to allocate XML value
GISMO_EXPORT char * makeValue( const std::string & value, gsXmlTree & data);

/// Helper to allocate matrix in XML pool
template<class T>
char * makeValue(const gsMatrix<T> & value, gsXmlTree & data,
                 bool transposed);

/// Helper to allocate XML attribute
GISMO_EXPORT gsXmlAttribute *  makeAttribute( const std::string & name,
				              const std::string & value, gsXmlTree & data);

/// Helper to allocate XML attribute with unsigned int value
GISMO_EXPORT gsXmlAttribute *  makeAttribute( const std::string & name,
					           const unsigned & value, gsXmlTree & data);

/// Helper to allocate XML node
GISMO_EXPORT gsXmlNode *  makeNode( const std::string & name, gsXmlTree & data);

/// Helper to allocate XML node with value
GISMO_EXPORT gsXmlNode * makeNode( const std::string & name,
			             const std::string & value, gsXmlTree & data);

/// Helper to create an XML comment node
GISMO_EXPORT gsXmlNode *  makeComment(const std::string &, gsXmlTree & data);

/// Helper to convert small unsigned to string
GISMO_EXPORT std::string to_string(const unsigned & i);

/// Helper to count the number of Objects (by tag) that exist in the
/// XML tree
GISMO_EXPORT int countByTag(const std::string & tag, gsXmlNode * root );

/// Helper to count the number of Objects (by name and type) that
/// exist in the XML tree
GISMO_EXPORT int  countByTagType(const std::string & tag,
                                 const std::string & type,
                                 gsXmlNode * root );

/// Helper to get the first object (by tag) if one exists in
/// the XML tree
GISMO_EXPORT gsXmlNode * firstByTag(const std::string & tag,
                                    gsXmlNode * root );

/// Helper to get the first object (by tag and type) if one exists in
/// the XML tree
GISMO_EXPORT gsXmlNode * firstByTagType(const std::string & tag,
                                         const std::string & type,
                                         gsXmlNode * root );

// Helper which finds a node matching \a tag and \a type in the XML
// tree
//GISMO_EXPORT gsXmlNode * anyByTagType(const std::string & tag,
//                                       const std::string & type,
//                                       gsXmlNode * root );

/// Helper to get any object (by tag) if one exists in the XML tree
GISMO_EXPORT gsXmlNode * anyByTag(const std::string & tag,
                                  gsXmlNode * root );

GISMO_EXPORT void getBoundaries(gsXmlNode                * node,
                                std::map<int, int>       & ids,
                                std::vector< patchSide > & result);

GISMO_EXPORT void getInterfaces(gsXmlNode* node,
                                const int d,
                                std::map<int, int>& ids,
                                std::vector< boundaryInterface > & result);

GISMO_EXPORT void appendBoxTopology(const gsBoxTopology& topology,
                                    gsXmlNode* node,
                                    gsXmlTree& data);

/// Helper to allocate XML node with gsMatrix value
template<class T>
gsXmlNode * makeNode( const std::string & name,
                      const gsMatrix<T> & value, gsXmlTree & data,
                      bool transposed = false );

/// Helper to fetch functions
///\todo read gsFunction instead
template<class T>
void getFunctionFromXml ( gsXmlNode * node, gsFunctionExpr<T> & result );

/// Helper to fetch matrices
template<class T>
void getMatrixFromXml ( gsXmlNode * node,
                        unsigned const & rows,
                        unsigned const & cols,
                        gsMatrix<T> & result );

/// Helper to insert matrices into XML
template<class T>
gsXmlNode * putMatrixToXml ( gsMatrix<T> const & mat,
                             gsXmlTree & data, std::string name = "Matrix");

/// Helper to fetch sparse entries
template<class T>
void getSparseEntriesFromXml ( gsXmlNode * node,
                               gsSparseEntries<T> & result );

/// Helper to insert sparse matrices into XML
template<class T>
gsXmlNode * putSparseMatrixToXml ( gsSparseMatrix<T> const & mat,
                                   gsXmlTree & data, std::string name = "SparseMatrix");

}// end namespace internal

}// end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXml.hpp)
#endif
