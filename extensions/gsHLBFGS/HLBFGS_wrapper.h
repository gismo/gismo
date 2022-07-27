#ifndef __HLBFGS_WRAPPER_H__
#define __HLBFGS_WRAPPER_H__

#include <functional>
#include <vector>
#include <cassert>

#include "gsHLBFGS.h"
//#include "stlbfgs/stlbfgs.h"

namespace UM
{
    struct LBFGS_Optimizer
    {
        typedef std::function<void(const std::vector<double> &x, double &f, std::vector<double> &g)> func_grad_eval;

        LBFGS_Optimizer(func_grad_eval func) : func_grad(func) {}
        void run(std::vector<double> &sol);

        const func_grad_eval func_grad;
        int maxiter = 10000; // Maximum number of quasi-Newton updates
        double gtol = 1e-10; // The iteration will stop when ||g||/max(1,||x||) <= gtol
        bool verbose = true;
    };
}
#endif //__HLBFGS_WRAPPER_H__
