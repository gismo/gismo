#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsMultiIndexIterators.h>
#include <gsUtils/gsMultiIndexIterators.hpp>

#define T real_t

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
TEMPLATE_INST
class gsMultiIndexIterator<index_t,-1>;
TEMPLATE_INST
class gsMultiIndexIterator<index_t,2>;
TEMPLATE_INST
class gsMultiIndexIterator<index_t,3>;

TEMPLATE_INST
class gsTensorGridIterator<index_t,-1>;
TEMPLATE_INST
class gsTensorGridIterator<index_t,2>;
TEMPLATE_INST
class gsTensorGridIterator<index_t,3>;

TEMPLATE_INST
class gsTensorGridVertexIterator<index_t,-1>;
TEMPLATE_INST
class gsTensorGridVertexIterator<index_t,2>;
TEMPLATE_INST
class gsTensorGridVertexIterator<index_t,3>;

TEMPLATE_INST
class gsTensorGridBoundaryIterator<index_t,-1>;
TEMPLATE_INST
class gsTensorGridBoundaryIterator<index_t,2>;
TEMPLATE_INST
class gsTensorGridBoundaryIterator<index_t,3>;

TEMPLATE_INST
class gsSimplexIterator<index_t,-1>;
TEMPLATE_INST
class gsSimplexIterator<index_t,2>;
TEMPLATE_INST
class gsSimplexIterator<index_t,3>;

TEMPLATE_INST
class gsCompositionIterator<index_t,-1>;
TEMPLATE_INST
class gsCompositionIterator<index_t,2>;
TEMPLATE_INST
class gsCompositionIterator<index_t,3>;
}

#undef T
