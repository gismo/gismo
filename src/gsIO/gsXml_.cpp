#include <gsCore/gsTemplateTools.h>

#include <gsIO/gsXml.h>
#include <gsIO/gsXml.hpp>


#define T real_t


namespace gismo
{

namespace internal
{

TEMPLATE_INST
gsXmlNode * makeNode( const std::string & name, 
                      const gsMatrix<T> & value, gsXmlTree & data,
                      bool transposed);

TEMPLATE_INST
void getMatrixFromXml ( gsXmlNode * node, 
                        unsigned const & rows, 
                        unsigned const & cols, 
                        gsMatrix<T> & result );

TEMPLATE_INST
void getFunctionFromXml ( gsXmlNode * node, gsFunctionExpr<T> & result);

// TEMPLATE_INST
// void getMatrixFromXml ( gsXmlNode * node, 
//                         unsigned const & rows, 
//                         unsigned const & cols, 
//                         gsMatrix<unsigned> & result );

TEMPLATE_INST
gsXmlNode * putMatrixToXml ( gsMatrix<T> const & mat, 
                             gsXmlTree & data, std::string name);

TEMPLATE_INST // used in gsXmlGenericUtils.hpp
gsXmlNode * putMatrixToXml ( gsMatrix<unsigned> const & mat, 
                             gsXmlTree & data, std::string name);

TEMPLATE_INST
void getSparseEntriesFromXml ( gsXmlNode * node, 
                               gsSparseEntries<T> & result );

TEMPLATE_INST
gsXmlNode * putSparseMatrixToXml ( gsSparseMatrix<T> const & mat, 
                                   gsXmlTree & data, std::string name);


} // namespace internal

} // namespace gismo

#undef T

