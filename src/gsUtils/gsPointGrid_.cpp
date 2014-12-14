#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsPointGrid.h>
#include <gsUtils/gsPointGrid.hpp>

#define T real_t

namespace gismo
{

TEMPLATE_INST
gsMatrix<T>::uPtr gsPointGrid( gsVector<T> const & a, gsVector<T> const & b, 
                                        gsVector<unsigned> const & np );

TEMPLATE_INST
void gsPointGrid( std::vector< gsVector<T>* > const & cwise, gsMatrix<T>& res);

TEMPLATE_INST
void gsPointGrid( std::vector< gsVector<T> > const & cwise, gsMatrix<T>& res);

TEMPLATE_INST
void tensorProduct( std::vector< gsVector<T>* > const & cwise, gsVector<T>& res);

TEMPLATE_INST
gsMatrix<T>::uPtr uniformPointGrid(const gsVector<T>& lower, 
                                            const gsVector<T>& upper, 
                                            int numPoints);

TEMPLATE_INST
gsVector<unsigned> uniformSampleCount (const gsVector<T>& lower, 
                                       const gsVector<T>& upper, 
                                       int numPoints);

TEMPLATE_INST
void uniformIntervals(const gsVector<T>& lower, 
                      const gsVector<T>& upper, 
                      std::vector< std::vector<T> >& intervals, 
                      int numIntervals);

}

#undef T
