
#include <gsCore/gsConfig.h>
#include <gsCore/gsExport.h>
#include <gsCore/gsBasis.h>

namespace gismo{
template<class T> inline
std::vector<gsSparseMatrix<T>> collocationMatrix1(const gsBasis<T> & b, const gsMatrix<T> & u)
{
    int dim = b.domainDim();
    std::vector<gsSparseMatrix<T>> result(dim+1, gsSparseMatrix<T>( u.cols(), b.size() ));
    std::vector<gsMatrix<T>> ev;
    gsMatrix<index_t> act;

    b.evalAllDers_into  (u.col(0), 1, ev);
    b.active_into(u.col(0), act);
    result[0].reservePerColumn( act.rows() );
    result[1].reservePerColumn( act.rows() );
    if (dim==2)
        result[2].reservePerColumn( act.rows() );
    for (index_t i=0; i!=act.rows(); ++i)
    {
        result[0].insert(0, act.at(i) ) = ev[0].at(i);
        result[1].insert(0, act.at(i) ) = ev[1].at(dim*i);
        if (dim == 2)
            result[2].insert(0, act.at(i) ) = ev[1].at(dim*i+1);
    }
    for (index_t k=1; k!=u.cols(); ++k)
    {
        b.evalAllDers_into  (u.col(k), 1, ev );
        b.active_into(u.col(k), act);
        for (index_t i=0; i!=act.rows(); ++i)
        {
            result[0].insert(k, act.at(i) ) = ev[0].at(i);
            result[1].insert(k, act.at(i) ) = ev[1].at(dim*i);
            if (dim == 2)
                result[2].insert(k, act.at(i) ) = ev[1].at(dim*i +1);
        }
    }

    result[0].makeCompressed();
    result[1].makeCompressed();
    if (dim == 2)
        result[2].makeCompressed();
    return result;
}

TEMPLATE_INST
std::vector<gsSparseMatrix<real_t>> collocationMatrix1(const gsBasis<real_t> & b, const gsMatrix<real_t> & u);

}

void pybind11_init_PPN(pybind11::module &m)
{
    pybind11::module ppn = m.def_submodule("ppn");
    // .def("dim", &Class::dim, "Returns the dimension of the basis")
    ppn.def("collocationMatrix1", &gismo::collocationMatrix1<real_t>, "returns the collocation matrix and its derivatives.");
}
