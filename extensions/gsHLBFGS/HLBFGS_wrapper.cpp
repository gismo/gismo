#include "HLBFGS_wrapper.h"
#include <iostream>
#include <iomanip> // Header file needed to use setw
#if WIN32
// disable int to size_t warning
#pragma warning(disable : 4267)
#endif

namespace UM
{
    static const std::function<void(int N, double *x, double *prev_x, double *f, double *g)> *local_func_grad = nullptr;
    static void static_func_grad(int N, double *x, double *prev_x, double *f, double *g)
    {
        (*local_func_grad)(N, x, prev_x, f, g);
    }
    // static void static_newiter_callback(int, int, double*, double*, double*, double*) {}

    static void static_newiter_callback(int iter, int call_iter, double *x, double *f, double *g, double *gnorm)
    {
        std::cout << "# iter " << iter << ": #func eval. " << call_iter << ", f = " << *f << ", ||g|| = " << *gnorm << std::endl;
        /*std::cout << "# iter "<< std::setw(4) << iter << ": #func eval. " << std::setw(4) << call_iter <<
            ", f = " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::right) << std::setprecision(6) << *f <<
            ", ||g|| = " << *gnorm << std::endl;*/
    }

    void LBFGS_Optimizer::run(std::vector<double> &sol)
    {
        int hlbfgs_info[20] = {0};
        double parameter[20] = {0};
        INIT_HLBFGS(parameter, hlbfgs_info);
        hlbfgs_info[3] = 1;       // b_m1qn3_ ? 1 : 0; // determines whether we use m1qn3
        hlbfgs_info[4] = maxiter; // max iterations
        hlbfgs_info[5] = verbose; // verbose

        parameter[5] = gtol;
        parameter[6] = gtol;

        std::function<void(int N, double *x, double *prev_x, double *f, double *g)> wrapfunc = [&](int N, double *x, double *, double *f, double *g)
        {
            std::vector<double> array_x(N), array_g(N);
            for (int i = 0; i < N; i++)
                array_x[i] = x[i];

            func_grad(array_x, *f, array_g);

            for (int i = 0; i < N; i++)
                g[i] = array_g[i];
        };

        local_func_grad = &wrapfunc;
        HLBFGS((int)sol.size(),
               (int)5,
               sol.data(),
               static_func_grad,
               nullptr,
               HLBFGS_UPDATE_Hessian,
               static_newiter_callback,
               parameter,
               hlbfgs_info);
    }
}
