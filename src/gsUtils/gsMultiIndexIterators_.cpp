#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsMultiIndexIterators.h>
#include <gsUtils/gsMultiIndexIterators.hpp>

namespace gismo
{
// Quote from the C++ Standard:
// 'An explicit instantiation that names a class template
// specialization is also an explicit instantiation of the same kind
// (declaration or deﬁnition) of each of its members (not including
// members inherited from base classes) that has not been previously
// explicitly specialized in the translation unit containing the
// explicit instantiation, except as described below.'
// This means we need to manually instantize the templated base class
// gsMultiIndexIterator
CLASS_TEMPLATE_INST gsMultiIndexIterator<index_t,-1>;
CLASS_TEMPLATE_INST gsMultiIndexIterator<index_t,1>;
CLASS_TEMPLATE_INST gsMultiIndexIterator<index_t,2>;
CLASS_TEMPLATE_INST gsMultiIndexIterator<index_t,3>;
CLASS_TEMPLATE_INST gsMultiIndexIterator<index_t,4>;

CLASS_TEMPLATE_INST gsTensorGridIterator<index_t,-1>;
CLASS_TEMPLATE_INST gsTensorGridIterator<index_t,1>;
CLASS_TEMPLATE_INST gsTensorGridIterator<index_t,2>;
CLASS_TEMPLATE_INST gsTensorGridIterator<index_t,3>;
CLASS_TEMPLATE_INST gsTensorGridIterator<index_t,4>;

CLASS_TEMPLATE_INST gsTensorGridVertexIterator<index_t,-1>;
CLASS_TEMPLATE_INST gsTensorGridVertexIterator<index_t,1>;
CLASS_TEMPLATE_INST gsTensorGridVertexIterator<index_t,2>;
CLASS_TEMPLATE_INST gsTensorGridVertexIterator<index_t,3>;
CLASS_TEMPLATE_INST gsTensorGridVertexIterator<index_t,4>;

CLASS_TEMPLATE_INST gsTensorGridBoundaryIterator<index_t,-1>;
CLASS_TEMPLATE_INST gsTensorGridBoundaryIterator<index_t,1>;
CLASS_TEMPLATE_INST gsTensorGridBoundaryIterator<index_t,2>;
CLASS_TEMPLATE_INST gsTensorGridBoundaryIterator<index_t,3>;
CLASS_TEMPLATE_INST gsTensorGridBoundaryIterator<index_t,4>;

CLASS_TEMPLATE_INST gsSimplexIterator<index_t,-1>;
CLASS_TEMPLATE_INST gsSimplexIterator<index_t,1>;
CLASS_TEMPLATE_INST gsSimplexIterator<index_t,2>;
CLASS_TEMPLATE_INST gsSimplexIterator<index_t,3>;
//CLASS_TEMPLATE_INST gsSimplexIterator<index_t,3>;

CLASS_TEMPLATE_INST gsCompositionIterator<index_t,-1>;
CLASS_TEMPLATE_INST gsCompositionIterator<index_t,1>;
CLASS_TEMPLATE_INST gsCompositionIterator<index_t,2>;
CLASS_TEMPLATE_INST gsCompositionIterator<index_t,3>;
//CLASS_TEMPLATE_INST gsCompositionIterator<index_t,4>;

}
