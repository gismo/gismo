/** @file gsTensorFunction.h

    @brief Provides declaration of TensorFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/** 
    \brief This class provides a tensor function represented by
    component functions with certain rank per component.

*/
template<class T>
class gsTensorFunction : public gsFunction<T>
{
private:
    gsTensorFunction( ){ } 

public:

    gsTensorFunction(const gsFunction<T> & f1, 
                     const gsFunction<T> & f2,
                     const gsVector<int> & ranks)
    : m_ranks(ranks)
    { 
        m_components.push_back(f1.clone());
        m_components.push_back(f2.clone());
    } 
    
    // to do, 3D
    //gsTensorFunction(const gsFunction<T> & f1, 

    ~gsTensorFunction( )
    { 
        freeAll(m_components);
    } 

public:

    int domainDim() const                     
    {
        return m_components[0]->domainDim() + m_components[1]->domainDim();
    }
    
    int targetDim() const                     
    {
        return m_ranks.size();
    }
    
    int rank( int i) const { return m_ranks[i]; }

    const gsVector<int> & ranks() const { return m_ranks; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT( u.rows() == 2, "wrong dimension" );
        const index_t n = m_ranks.size();
        result.resize(n, u.cols() );

        gsMatrix<T> ev1, ev2;
        for (index_t j=0; j!= u.cols(); ++j)
        {
            m_components[0]->eval_into(u.col(j).row(0), ev1);
            m_components[1]->eval_into(u.col(j).row(1), ev2);            

            int pos = 0;
            for (index_t k=0; k!=n; ++k)
            {
                result(k,j) = (ev1.middleRows(pos, m_ranks[k]).transpose() * 
                               ev2.middleRows(pos, m_ranks[k])
                              ).sum();
                pos += m_ranks[k];
            }
        }
    }

    // Part in variable \a s of Function coordinate \a k
    void evalComponent_into(int k, int s, const gsMatrix<T>& u, gsMatrix<T>& result) const
    { 
        result.resize(m_ranks[k], u.cols() );
        GISMO_ASSERT( u.rows() ==1, "expect 1d pts");

        gsMatrix<T> ev;
        const int pos = m_ranks.segment(1,k).sum();// coord. position

        for (index_t j=0; j!= u.cols(); ++j)
        {
            m_components[s]->eval_into(u.col(j), ev);
            result.col(j) = ev.middleRows(pos, m_ranks[k]);
        }
    }

    // Part in variable \a k of Function coordinate \a i
    const gsFunction<T> & component(int k) const
    {
        return *m_components[k];
    }
    
protected:
    std::vector< const gsFunction<T>*> m_components;
    gsVector<int> m_ranks;
};


}// namespace gismo

