/** @file 

    @brief 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>
#include <gsUtils/gsKNNTree.h>

using namespace gismo;

namespace gismo {

template<>
struct gsKNNTreeTraits< gsVector<real_t> >
{
    typedef gsVector<real_t,3> T;
    static inline std::size_t size() { return 3; }
  
  static inline bool islhalf(const T& lhs, const T& rhs, std::size_t axis) { return lhs[axis] < rhs[axis]; }
  
  static inline double fabs(const T& lhs, const T& rhs, std::size_t axis) { return std::abs(lhs[axis] - rhs[axis]); }
  
  static inline double distance(const T& lhs, const T& rhs)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < size(); ++i) {
      result += fabs(lhs, rhs, i);
    }
    return result;
  }
};

}


int main(int argc, char *argv[])
{
    std::string fn = "fitting/deepdrawingC.xml";
    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    std::vector<std::pair<gsVector<>, int> > data;
    for(index_t i = 0; i<xyz.cols();++i)
        data.push_back( std::make_pair(xyz.col(i),  i) );

    gsInfo <<"Number of points: "<< data.size() <<"\n";
        
    gsKNNTree<gsVector<>, int > tree(data);

    gsInfo <<"k-d tree size      "<< tree.size() <<"\n";
    gsInfo <<"k-d tree dimension "<< tree.dimension() <<"\n";
    gsInfo <<"k-d tree root      "<< tree.root() <<"\n";

    return EXIT_SUCCESS;
}
