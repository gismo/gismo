#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsPointIterator.h>

#include <gsUtils/gsMultiIndexIterators.hpp>

#define T real_t
namespace gismo
{

CLASS_TEMPLATE_INST
gsTensorPointGridIterator<T,2,index_t>;
CLASS_TEMPLATE_INST
gsTensorPointGridIterator<T,3,index_t>;
CLASS_TEMPLATE_INST
gsTensorPointGridIterator<T,-1,index_t>;
CLASS_TEMPLATE_INST
gsSimplexPointGridIterator<T,2,index_t>;
CLASS_TEMPLATE_INST
gsSimplexPointGridIterator<T,3,index_t>;
CLASS_TEMPLATE_INST
gsSimplexPointGridIterator<T,-1,index_t>;

}

#undef T
