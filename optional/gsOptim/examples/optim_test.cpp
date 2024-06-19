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

template <class T>
class rastrigin : public gsOptProblem<T>
{
public:
    rastrigin() {}

    T evalObj(const gsAsConstVector<T> & u) const
    {
        const index_t n = u.size();

        const real_t A = 10;

        real_t obj_val = A*n + u.array().pow(2).sum() - A * (2 * EIGEN_PI * u).array().cos().sum();

        return obj_val;
    }
};

template <class T>
class sphere : public gsOptProblem<T>
{
public:
    sphere() {}

    T evalObj(const gsAsConstVector<T> & u) const
    {
        real_t obj_val = u.dot(u);
        return obj_val;
    }

    void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {
        result = 2.0*u;
    }

};

template <class T>
class booth : public gsOptProblem<T>
{
public:
    booth() {}

    T evalObj(const gsAsConstVector<T> & u) const
    {
        real_t x_1 = u(0);
        real_t x_2 = u(1);

        real_t obj_val = std::pow(x_1 + 2*x_2 - 7.0,2) + std::pow(2*x_1 + x_2 - 5.0,2);
        return obj_val;;
    }

    void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {
        real_t x_1 = u(0);
        real_t x_2 = u(1);

        result(0) = 2*(x_1 + 2*x_2 - 7.0) + 2*(2*x_1 + x_2 - 5.0)*2;
        result(1) = 2*(x_1 + 2*x_2 - 7.0)*2 + 2*(2*x_1 + x_2 - 5.0);
    }
};

int main()
{
    gsVector<real_t> x;


    std::vector<std::string> methods;
    std::vector<std::string> methods_excluded{"PSODV","SUMT"};

    methods = {"DE","DEPRMM","PSO"};
    for (typename std::vector<std::string>::const_iterator m = methods.begin(); m!=methods.end(); m++)
    {
        ackley<real_t> problem;
        gsOptim<real_t>::uPtr solver = gsOptim<real_t>::get(*m,&problem);
        x = 2.0 * gsVector<real_t>::Ones(2); // initial values: (2,2)
        solver->solve(x);
        x = solver->currentDesign();
        gsInfo<<*m <<" solver: solution to Ackley test:\n" << x << "\n";
    }

    for (typename std::vector<std::string>::const_iterator m = methods.begin(); m!=methods.end(); m++)
    {
        rastrigin<real_t> problem;
        gsOptim<real_t>::uPtr solver = gsOptim<real_t>::get(*m,&problem);
        x = 2.0 * gsVector<real_t>::Ones(2); // initial values: (2,2)
        solver->solve(x);
        x = solver->currentDesign();
        gsInfo<<*m <<" solver: solution to Rastrigin test:\n" << x << "\n";
    }

    methods = {"BFGS","LBFGS","CG","GD","NM"};
    for (typename std::vector<std::string>::const_iterator m = methods.begin(); m!=methods.end(); m++)
    {
        sphere<real_t> problem;
        gsOptim<real_t>::uPtr solver = gsOptim<real_t>::get(*m,&problem);
        x = 2.0 * gsVector<real_t>::Ones(2); // initial values: (2,2)
        solver->solve(x);
        x = solver->currentDesign();
        gsInfo<<*m <<" solver: solution to Sphere test:\n" << x << "\n";
    }

    for (typename std::vector<std::string>::const_iterator m = methods.begin(); m!=methods.end(); m++)
    {
        booth<real_t> problem;
        gsOptim<real_t>::uPtr solver = gsOptim<real_t>::get(*m,&problem);
        x = 2.0 * gsVector<real_t>::Ones(2); // initial values: (2,2)
        solver->solve(x);
        x = solver->currentDesign();
        gsInfo<<*m <<" solver: solution to Booth's test:\n" << x << "\n";
    }
    return 0;
}