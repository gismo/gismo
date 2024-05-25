/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/



#include <gismo.h>
#include <gsOptimizer/gsOptProblem.h>
#include <gsOptim/gsOptim.h>

using namespace gismo;

template <class T>
class ackley : public gsOptProblem<T>
{
public:
    ackley() {}

    T evalObj(const gsAsConstVector<T> & u) const
    {
        const T x = u[0];
        const T y = u[1];

        const T obj_val = 20 + math::exp(1) - 20*math::exp( -0.2*math::sqrt(0.5*(x*x + y*y)) ) - math::exp( 0.5*(math::cos(2 * EIGEN_PI * x) + math::cos(2 * EIGEN_PI * y)) );

        return obj_val;
    }
};

int main()
{
    gsVector<real_t> x = 2.0 * gsVector<real_t>::Ones(2); // initial values: (2,2)

    ackley<real_t> problem;

    std::vector<std::string> methods{"BFGS","LBFGS","CG","GD","NM","DE","DEPRMM","PSO"};
    std::vector<std::string> methods_excluded{"PSODV","SUMT"};
    for (typename std::vector<std::string>::const_iterator m = methods.begin(); m!=methods.end(); m++)
    {
        gsOptim<real_t>::uPtr solver = gsOptim<real_t>::get(*m,&problem);
        solver->solve(x);
        x = solver->currentDesign();
        gsInfo<<*m <<" solver: solution to Ackley test:\n" << x << "\n";
    }
    return 0;
}