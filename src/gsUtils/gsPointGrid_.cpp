#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsPointGrid.h>
#include <gsUtils/gsPointGrid.hpp>

#define T real_t

namespace gismo
{

TEMPLATE_INST
gsVector<unsigned> uniformSampleCount (const gsVector<T>& lower, 
                                       const gsVector<T>& upper, 
                                       int numPoints);

TEMPLATE_INST
gsMatrix<T> uniformPointGrid(const gsVector<T>& lower, 
                             const gsVector<T>& upper, 
                             int numPoints);

TEMPLATE_INST
gsMatrix<T> gsPointGrid( gsVector<T> const & a, gsVector<T> const & b, 
                         gsVector<unsigned> const & np );

TEMPLATE_INST
void uniformIntervals(const gsVector<T>& lower, 
                      const gsVector<T>& upper, 
                      std::vector< std::vector<T> >& intervals, 
                      int numIntervals);

}

#undef T
