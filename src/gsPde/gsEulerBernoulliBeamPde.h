
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{


/** @brief
    The differential equation describing the linear Euler-Bernoulli beam.

    This is a 1D model and can thus only be assembled on 1D geometries.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T=real_t>
class gsEulerBernoulliBeamPde : public gsPde<T>
{
private:
    gsEulerBernoulliBeamPde();
    using gsPde<T>::m_unknownDim;
    using gsPde<T>::m_solution;

public:
    gsEulerBernoulliBeamPde(
        const gsMultiPatch<T>         &domain,
        const gsBoundaryConditions<T> &bc,
        gsFunction<T> * rhs )
        : gsPde<T>(domain,bc), m_rhs(rhs)
    {
        m_unknownDim.push_back(1);
    }

    gsEulerBernoulliBeamPde(
        const gsMultiPatch<T>         &domain,
        const gsBoundaryConditions<T> &bc,
        const gsFunction<T> & rhs )
        : gsPde<T>(domain,bc), m_rhs(rhs.clone())
    {
        m_unknownDim.push_back(1);
    }
    // COMPATIBILITY CONSTRUCTOR, DO NOT USE
    gsEulerBernoulliBeamPde(gsFunction<T>   *rhs)
        : m_rhs(rhs)
    {
            m_unknownDim.push_back(1);
    }


    gsFunction<T> * rhs() const         { return m_rhs; }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "Euler-Bernoulli beam equation\n";
        return os;
    }

private:

    gsFunction<T> *m_rhs;

}; // class gsEulerBernoulliBeamPde

} // namespace gismo
